######## 
# Get study design matrix and R packages
######## 
library(data.table)
# Create design matrix for which we'll pull in data
design <- expand.grid(time = gl(3, 1, labels = c("pre", "post", "during")), pop = gl(2, 1, labels = c("denmark", "london")))
design <- design[-3,]
                      #type = gl(3, 1, labels = c("exons", "gwas", "neutral"))
design2 <- expand.grid(time1 = gl(3, 1, labels = c("pre", "post", "during")), time2 = gl(3, 1, labels = c("pre", "post", "during")), pop = gl(2, 1, labels = c("denmark", "london")))
design2<-subset(design2, !(design2$time1 == design2$time2)); design2 <- design2[c(1,7:8,10),]


min_n=10; min_maf=0.05
drop_samples_missing=.5

method='sliding'
path = "" ##insert path to DATA directory inluding the permutation files generated in the bash script under /SampleNames/permute/permute_{trial}
# for method=bins
if (method=='bins') {
  bin_breaks <- c(0.07, 0.1, 0.25, 0.5)
  bins <- matrix(ncol=2, nrow=length(bin_breaks)); bins[,1] <- c(0,bin_breaks[-length(bin_breaks)]); bins[,2] <- c(bin_breaks)
}
# for method=sliding: sliding windows centered on each percentage with at least 200 nearby variants
if (method == 'sliding') {
  min_neutral_nsites <- 200
  if (min_neutral_nsites == '200') {
    mafs <- seq(0.005,0.40, 0.01); window_size <- 0.05*rep(1, length(mafs))
    window_size[mafs>=0.135] <- 0.06; window_size[mafs>=0.145] <- 0.07; window_size[mafs>=0.155] <- 0.08
    window_size[mafs>=0.165] <- 0.09; window_size[mafs>=0.175] <- 0.11; window_size[mafs>=0.185] <- 0.12
    window_size[mafs>=0.205] <- 0.13; window_size[mafs>=0.225] <- 0.15; window_size[mafs>=0.255] <- 0.16
    window_size[mafs>=0.275] <- 0.17; window_size[mafs>=0.305] <- 0.18; window_size[mafs>=0.335] <- 0.19
    window_size[mafs>=0.345] <- 0.2; window_size[mafs>=0.365] <- 0.21; window_size[mafs>=0.375] <- 0.22
    window_size[mafs>=0.395] <- 0.24
    bins <- cbind(mafs-0.005, mafs+0.005, mafs-window_size/2, mafs+window_size/2)
    bins[nrow(bins),2] <- 0.5
  }
  # we used the following line to check this once we had the set of neutral sites
  # for (i in 1:length(mafs)) { print(paste("maf: ", mafs[i], "; nearby sites: ", sum(info_neut$maf > mafs[i] - window_size[i]/2 & info_neut$maf <= mafs[i] + window_size[i]/2), sep=""))}
}

options(warn=-1)
######## 
# Get trimmed aDNA frequencies
######## 
# Info files with the number and alternate allele freq in each population
# mac >= 3; biallelic; minQ >= 30
# exons filtered for sites within bed file; rest mapped to the bed file
######## 
exon_table=list()
neutral_table=list()
gwas_table=list()
permute_range <- 1:100
for(i in permute_range){
    
    info_gwas <- fread(paste(path,"/genoliks/immune.frq",sep=""))[,-c(3:4)]
    info_neut <- fread(paste(path,"/genoliks/neutral.frq",sep=""))[,-c(3:4)]
    info_exon <- fread(paste(path,"/genoliks/exon.frq",sep=""))[,-c(3:4)]

    colnames(info_gwas) <- c("chr", "pos", "ref", "alt")
    colnames(info_neut) <- c("chr", "pos", "ref", "alt")
    colnames(info_exon) <- c("chr", "pos", "ref", "alt")

    info_gwas$ref <- gsub(":.*", "",info_gwas$ref); info_gwas$alt <- gsub(":.*", "",info_gwas$alt)
    info_neut$ref <- gsub(":.*", "",info_neut$ref); info_neut$alt <- gsub(":.*", "",info_neut$alt)
    info_exon$ref <- gsub(":.*", "",info_exon$ref); info_exon$alt <- gsub(":.*", "",info_exon$alt)

    paste(info_gwas$chr,info_gwas$pos,sep="_") -> info_gwas$site
    paste(info_neut$chr,info_neut$pos,sep="_") -> info_neut$site
    paste(info_exon$chr,info_exon$pos,sep="_") -> info_exon$site
    
    exon_table[[i]] <- info_exon
    gwas_table[[i]] <- info_gwas
    neutral_table[[i]] <- info_neut

    #assign(paste("info_exon_trial",i,sep="_"),info_exon)
    #assign(paste("info_gwas_trial",i,sep="_"),info_gwas)
    #assign(paste("info_neut_trial",i,sep="_"),info_neut)
    
    #exon_table <- append(exon_table,paste("info_exon_trial",i,sep="_"),after=i-1)
    #gwas_table <- append(gwas_table,paste("info_gwas_trial",i,sep="_"),after=i-1)
    #neutral_table <- append(neutral_table,paste("info_neut_trial",i,sep="_"),after=i-1)
}


##run these prior to caluclating the permuted AFs

##run initial AF calculations with the original keep files provided on Github
for (x in permute_range){
    for (i in c(3:nrow(design),1:2)) {
      # read in exon data
      #exon_table <- paste("info_exon_trial",x,sep="_")
      #assign(info_exon,paste("info_exon_trial",x,sep="_"))
      info_exon = exon_table[[x]]
      name=paste(path,"/genoliks/genolik.exons_",design$pop[i],"_", design$time[i],".genolik",sep="")
      n2=paste(design$pop[i],design$time[i],sep="_"); d <- as.data.frame(fread(name)); rm(name)
      snames=read.delim(paste(path,"/SampleNames/", design$pop[i], "_", design$time[i],"_exons.012.indv",sep=""), header=F)
      ncol(d)-2 == 3*nrow(snames)
      if (i == 3) {paste(d[,1],d[,2],sep="_") -> info_exon$site_check}
      # remove site. change missing data to NA. get number of samples. 
      d <- d[,-(1:2)] ; d[d == -1] <- NA; n <- ncol(d)/3
      # report data on missingness per sample
      assign(paste("missing",n2, "exon", sep="_"), value = (colSums(is.na(d))/nrow(d))[seq(1,3*n,3)])
      # remove samples missing too much data
      ids_directory <-  paste(path,"/SampleNames/keep",sep="")
      keep_ids <- read.csv(paste(ids_directory,n2,"exon",sep="_"),sep="",header=FALSE)
      keep <-  which(snames[,'V1'] %in% keep_ids$V1) 
      #keep <- which(get(paste("missing",n2, "exon", sep="_")) <= drop_samples_missing)
      assign(paste("keep",n2,"exon",sep="_"), value=snames[keep,])
      k2 <- as.data.frame(matrix(nrow=nrow(d), ncol=1)); for (k in keep) {k2 <- cbind(k2, d[,(k*3-2):(k*3)])}
      d <- k2[,-1]; rm(keep, k2, k); n <- ncol(d)/3
      # get the allele frequency, and the number of individuals it was called from
      info_exon <- exon_table[[x]]
      info_exon[[paste(n2,"alternate",sep=".")]] <- rowMeans(d[,seq(3,n*3,3)],na.rm=T) + rowMeans(d[,seq(2,n*3,3)],na.rm=T)/2
      info_exon[[paste(n2,"called",sep=".")]] <- (n*3-rowSums(is.na(d)))/3
      # Report number of samples included 
      #assign(paste("samples",n2, "exon", sep="_"), n)
      # cleanup
      #assign(paste("info_exon_trial",x,sep="_"),info_exon)
      exon_table[[x]] = info_exon
      rm(d, n)

      # read in gwas data
      info_gwas = gwas_table[[x]]
      name=paste(path,"/genoliks/genolik.gwas_",design$pop[i],"_", design$time[i],".genolik",sep="")
      n2=paste(design$pop[i],design$time[i],sep="_"); d <- as.data.frame(fread(name)); rm(name)
      snames=read.delim(paste(path,"/SampleNames/", design$pop[i], "_", design$time[i],"_immune.012.indv",sep=""), header=F)
      ncol(d)-2 == 3*nrow(snames)
      if (i == 3) {paste(d[,1],d[,2],sep="_") -> info_gwas$site_check}
      # remove site. change missing data to NA. get number of samples. 
      d <- d[,-(1:2)] ; d[d == -1] <- NA; n <- ncol(d)/3
      # report data on missingness per sample
      assign(paste("missing",n2, "gwas", sep="_"), value = (colSums(is.na(d))/nrow(d))[seq(1,3*n,3)])
      # remove samples missing too much data
      #keep <- which(get(paste("missing",n2, "gwas", sep="_")) <= drop_samples_missing)
      ids_directory <-  paste(path,"/SampleNames/keep",sep="")
      keep_ids <- read.csv(paste(ids_directory,n2,"gwas",sep="_"),sep="",header=FALSE)
      keep <-  which(snames[,'V1'] %in% keep_ids$V1) 
      assign(paste("keep",n2,"gwas",sep="_"), value=snames[keep,])
      k2 <- as.data.frame(matrix(nrow=nrow(d), ncol=1)); for (k in keep) {k2 <- cbind(k2, d[,(k*3-2):(k*3)])}
      d <- k2[,-1]; rm(keep, k2, k); n <- ncol(d)/3
      # get the allele frequency, and the number of individuals it was called from
      info_gwas[[paste(n2,"alternate",sep=".")]] <- rowMeans(d[,seq(3,ncol(d),3)],na.rm=T) + rowMeans(d[,seq(2,ncol(d),3)],na.rm=T)/2
      info_gwas[[paste(n2,"called",sep=".")]] <- (ncol(d)-rowSums(is.na(d)))/3
      # Report number of samples included 
      assign(paste("samples",n2, "gwas", sep="_"), n)
      gwas_table[[x]] = info_gwas
      # cleanup
      rm(d, n)

      # read in neutral data
      info_neut = neutral_table[[x]]
      name=paste(path,"/genoliks/genolik.neutral_",design$pop[i],"_", design$time[i],".genolik",sep="")
      n2=paste(design$pop[i],design$time[i],sep="_"); d <- as.data.frame(fread(name)); rm(name)
      snames=read.delim(paste(path,"/SampleNames/", design$pop[i], "_", design$time[i],"_neutral.012.indv",sep=""), header=F)
      ncol(d)-2 == 3*nrow(snames)
      if (i == 3) {paste(d[,1],d[,2],sep="_") -> info_neut$site_check}
      # remove site. change missing data to NA. get number of samples.
      d <- d[,-(1:2)] ; d[d == -1] <- NA; n <- ncol(d)/3
      # report data on missingness per sample
      assign(paste("missing",n2, "neutral", sep="_"), value = (colSums(is.na(d))/nrow(d))[seq(1,3*n,3)])
      # remove samples missing too much data
      #keep <- which(get(paste("missing",n2, "neutral", sep="_")) <= drop_samples_missing)
      ids_directory <-  paste(path,"/SampleNames/keep",sep="")
      keep_ids <- read.csv(paste(ids_directory,n2,"neutral",sep="_"),sep="",header=FALSE)
      keep <-  which(snames[,'V1'] %in% keep_ids$V1) 
      assign(paste("keep",n2,"neutral",sep="_"), value=snames[keep,])
      k2 <- as.data.frame(matrix(nrow=nrow(d), ncol=1)); for (k in keep) {k2 <- cbind(k2, d[,(k*3-2):(k*3)])}
      d <- k2[,-1]; rm(keep, k2, k); n <- ncol(d)/3
      # get the allele frequency, and the number of individuals it was called from
      info_neut[[paste(n2,"alternate",sep=".")]] <- rowMeans(d[,seq(3,ncol(d),3)],na.rm=T) + rowMeans(d[,seq(2,ncol(d),3)],na.rm=T)/2
      info_neut[[paste(n2,"called",sep=".")]] <- (ncol(d)-rowSums(is.na(d)))/3
      # Report number of samples included 
      assign(paste("samples",n2, "neutral", sep="_"), n)
      neutral_table[[x]] = info_neut
      # cleanup
      rm(d, n)
    }; rm(i,n2)
}

##update exonic AFs with permuted samples


for (x in permute_range){
    info_exon = exon_table[[x]]
    name_pre = fread(paste(path,"/genoliks/genolik.exons_london_pre.genolik",sep=""),sep=" ")
    name_post = fread(paste(path,"/genoliks/genolik.exons_london_post.genolik",sep=""),sep=" ")
    name_post <- name_post[,3:191]
    d <- cbind(name_pre,name_post)
    pre_names = fread(paste(path,"/SampleNames/london_pre_exons.012.indv",sep=""),sep=" ",header=FALSE)
    post_names = fread(paste(path,"/SampleNames/london_post_exons.012.indv",sep=""),sep=" ",header=FALSE)
    snames = rbind(pre_names,post_names)
    d <- d[,-(1:2)] ; d[d == -1] <- NA; n <- ncol(d)/3
    #assign(paste("missing_london_pre_exon", sep="_"), value = (colSums(is.na(d))/nrow(d))[seq(1,3*n,3)])
    keep_ids_pre <- read.csv(paste(path,"/SampleNames/permute/permute_",x,"/keep_london_pre_exon",sep=""),sep="",header=FALSE)
    keep_ids_post <- read.csv(paste(path,"/SampleNames/permute/permute_",x,"/keep_london_post_exon",sep=""),sep="",header=FALSE)
    keep_pre <-  which(snames$V1 %in% keep_ids_pre$V1) 
    keep_post <-  which(snames$V1 %in% keep_ids_post$V1)
    assign(paste("keep_london_pre_exon",sep="_"), value=snames[keep_pre,])
    assign(paste("keep_london_post_exon",sep="_"), value=snames[keep_post,])
    k2_pre <- as.data.frame(matrix(nrow=nrow(d), ncol=1)); for (k in keep_pre) {k2_pre <- cbind(k2_pre, d[,(k*3-2):(k*3)])}
    d_pre <- k2_pre[,-1]; rm(keep_pre, k2_pre, k); n_pre <- ncol(d_pre)/3
    k2_post <- as.data.frame(matrix(nrow=nrow(d), ncol=1)); for (k in keep_post) {k2_post <- cbind(k2_post, d[,(k*3-2):(k*3)])}
    d_post <- k2_post[,-1]; rm(keep_post, k2_post, k); n_post <- ncol(d_post)/3
    info_exon[[paste("london_pre","alternate",sep=".")]] <- rowMeans(d_pre[,seq(3,n_pre*3,3)],na.rm=T) + rowMeans(d_pre[,seq(2,n_pre*3,3)],na.rm=T)/2
    info_exon[[paste("london_pre","called",sep=".")]] <- (n_pre*3-rowSums(is.na(d_pre)))/3
    assign(paste("samples","london_pre", "exon", sep="_"), n_pre)
    info_exon[[paste("london_post","alternate",sep=".")]] <- rowMeans(d_post[,seq(3,n_post*3,3)],na.rm=T) + rowMeans(d_post[,seq(2,n_post*3,3)],na.rm=T)/2
    info_exon[[paste("london_post","called",sep=".")]] <- (n_post*3-rowSums(is.na(d_post)))/3
    assign(paste("samples","london_post", "exon", sep="_"), n_post)
    exon_table[[x]] = info_exon
    #print(paste("New Exon Pre Cohort - original pre:", length(keep_ids_pre[keep_ids_pre$V1 %in% pre_names$V1,]),"original post:",length(keep_ids_pre[keep_ids_pre$V1 %in% post_names$V1,]),sep=" "))
    #print(paste("New Exon Post Cohort - original pre:", length(keep_ids_post[keep_ids_post$V1 %in% pre_names$V1,]),"original post:",length(keep_ids_post[keep_ids_post$V1 %in% post_names$V1,]),sep=" "))
    rm(d, n)
        
    }



##update gwas AFs with permuted samples

for (x in permute_range){
    info_gwas = gwas_table[[x]]
    name_pre = fread(paste(path,"/genoliks/genolik.gwas_london_pre.genolik",sep=""),sep=" ")
    name_post = fread(paste(path,"/genoliks/genolik.gwas_london_post.genolik",sep=""),sep=" ")
    name_post <- name_post[,3:296]
    d <- cbind(name_pre,name_post)
    pre_names = fread(paste(path,"/SampleNames/london_pre_immune.012.indv",sep=""),sep=" ",header=FALSE)
    post_names = fread(paste(path,"/SampleNames/london_post_immune.012.indv",sep=""),sep=" ",header=FALSE)
    snames = rbind(pre_names,post_names)
    d <- d[,-(1:2)] ; d[d == -1] <- NA; n <- ncol(d)/3
    #assign(paste("missing_london_pre_gwas", sep="_"), value = (colSums(is.na(d))/nrow(d))[seq(1,3*n,3)])
    keep_ids_pre <- read.csv(paste(path,"/SampleNames/permute/permute_",x,"/keep_london_pre_gwas",sep=""),sep="",header=FALSE)
    keep_ids_post <- read.csv(paste(path,"/SampleNames/permute/permute_",x,"/keep_london_post_gwas",sep=""),sep="",header=FALSE)
    keep_pre <-  which(snames$V1 %in% keep_ids_pre$V1) 
    keep_post <-  which(snames$V1 %in% keep_ids_post$V1)
    assign(paste("keep_london_pre_gwas",sep="_"), value=snames[keep_pre,])
    assign(paste("keep_london_post_gwas",sep="_"), value=snames[keep_post,])
    k2_pre <- as.data.frame(matrix(nrow=nrow(d), ncol=1)); for (k in keep_pre) {k2_pre <- cbind(k2_pre, d[,(k*3-2):(k*3)])}
    d_pre <- k2_pre[,-1]; rm(keep_pre, k2_pre, k); n_pre <- ncol(d_pre)/3
    k2_post <- as.data.frame(matrix(nrow=nrow(d), ncol=1)); for (k in keep_post) {k2_post <- cbind(k2_post, d[,(k*3-2):(k*3)])}
    d_post <- k2_post[,-1]; rm(keep_post, k2_post, k); n_post <- ncol(d_post)/3
    info_gwas[[paste("london_pre","alternate",sep=".")]] <- rowMeans(d_pre[,seq(3,n_pre*3,3)],na.rm=T) + rowMeans(d_pre[,seq(2,n_pre*3,3)],na.rm=T)/2
    info_gwas[[paste("london_pre","called",sep=".")]] <- (n_pre*3-rowSums(is.na(d_pre)))/3
    assign(paste("samples","london_pre", "gwas", sep="_"), n_pre)
    info_gwas[[paste("london_post","alternate",sep=".")]] <- rowMeans(d_post[,seq(3,n_post*3,3)],na.rm=T) + rowMeans(d_post[,seq(2,n_post*3,3)],na.rm=T)/2
    info_gwas[[paste("london_post","called",sep=".")]] <- (n_post*3-rowSums(is.na(d_post)))/3
    assign(paste("samples","london_post", "gwas", sep="_"), n_post)
    #print(paste("New GWAS Pre Cohort - original pre:", length(keep_ids_pre[keep_ids_pre$V1 %in% pre_names$V1,]),"original post:",length(keep_ids_pre[keep_ids_pre$V1 %in% post_names$V1,]),sep=" "))
    #print(paste("New GWAS Post Cohort - original pre:", length(keep_ids_post[keep_ids_post$V1 %in% pre_names$V1,]),"original post:",length(keep_ids_post[keep_ids_post$V1 %in% post_names$V1,]),sep=" "))
    gwas_table[[x]] = info_gwas 
    rm(d, n)
    }
                        


##update neutral AFs with permuted samples

for (x in permute_range){
    info_neut = neutral_table[[x]]
    name_pre = fread(paste(path,"/genoliks/genolik.neutral_london_pre.genolik",sep=""),sep=" ")
    name_post = fread(paste(path,"/genoliks/genolik.neutral_london_post.genolik",sep=""),sep=" ")
    name_post <- name_post[,3:296]
    d <- cbind(name_pre,name_post)
    pre_names = fread(paste(path,"/SampleNames/london_pre_neutral.012.indv",sep=""),sep=" ",header=FALSE)
    post_names = fread(paste(path,"/SampleNames/london_post_neutral.012.indv",sep=""),sep=" ",header=FALSE)
    snames = rbind(pre_names,post_names)
    d <- d[,-(1:2)] ; d[d == -1] <- NA; n <- ncol(d)/3
    #assign(paste("missing_london_pre_neutral", sep="_"), value = (colSums(is.na(d))/nrow(d))[seq(1,3*n,3)])
    keep_ids_pre <- read.csv(paste(path,"/SampleNames/permute/permute_",x,"/keep_london_pre_neutral",sep=""),sep="",header=FALSE)
    keep_ids_post <- read.csv(paste(path,"/SampleNames/permute/permute_",x,"/keep_london_post_neutral",sep=""),sep="",header=FALSE)
    keep_pre <-  which(snames$V1 %in% keep_ids_pre$V1) 
    keep_post <-  which(snames$V1 %in% keep_ids_post$V1)
    assign(paste("keep_london_pre_neutral",sep="_"), value=snames[keep_pre,])
    assign(paste("keep_london_post_neutral",sep="_"), value=snames[keep_post,])
    k2_pre <- as.data.frame(matrix(nrow=nrow(d), ncol=1)); for (k in keep_pre) {k2_pre <- cbind(k2_pre, d[,(k*3-2):(k*3)])}
    d_pre <- k2_pre[,-1]; rm(keep_pre, k2_pre, k); n_pre <- ncol(d_pre)/3
    k2_post <- as.data.frame(matrix(nrow=nrow(d), ncol=1)); for (k in keep_post) {k2_post <- cbind(k2_post, d[,(k*3-2):(k*3)])}
    d_post <- k2_post[,-1]; rm(keep_post, k2_post, k); n_post <- ncol(d_post)/3
    info_neut[[paste("london_pre","alternate",sep=".")]] <- rowMeans(d_pre[,seq(3,n_pre*3,3)],na.rm=T) + rowMeans(d_pre[,seq(2,n_pre*3,3)],na.rm=T)/2
    info_neut[[paste("london_pre","called",sep=".")]] <- (n_pre*3-rowSums(is.na(d_pre)))/3
    assign(paste("samples","london_pre", "neutral", sep="_"), n_pre)
    info_neut[[paste("london_post","alternate",sep=".")]] <- rowMeans(d_post[,seq(3,n_post*3,3)],na.rm=T) + rowMeans(d_post[,seq(2,n_post*3,3)],na.rm=T)/2
    info_neut[[paste("london_post","called",sep=".")]] <- (n_post*3-rowSums(is.na(d_post)))/3
    assign(paste("samples","london_post", "neutral", sep="_"), n_post)
    #print(paste("New Neutral Pre Cohort - original pre:", length(keep_ids_pre[keep_ids_pre$V1 %in% pre_names$V1,]),"original post:",length(keep_ids_pre[keep_ids_pre$V1 %in% post_names$V1,]),sep=" "))
    #print(paste("New Neutral Post Cohort - original pre:", length(keep_ids_post[keep_ids_post$V1 %in% pre_names$V1,]),"original post:",length(keep_ids_post[keep_ids_post$V1 %in% post_names$V1,]),sep=" "))
    neutral_table[[x]] = info_neut
    rm(d, n)
    }




#sum(info_gwas$site == info_gwas$site_check) == nrow(info_gwas)
#sum(info_exon$site == info_exon$site_check) == nrow(info_exon)
#sum(info_neut$site == info_neut$site_check) == nrow(info_neut)

#for (f in ls(pattern="missing")) {print(paste(f,sum(get(f) < 0.5), sep=": "))}; rm(f)


#######

########
# calculate alternate and minor allele frequency
# calculate a "London" frequency and a "Denmark" frequency as the mean between time points within that population
# NOTE: THIS TAKES US UP TO THE FIRST 20 COLUMNS
########
# let's do gwas first
for (x in permute_range){
    
    info_gwas = gwas_table[[x]]
    info_exon = exon_table[[x]]
    info_neut= neutral_table[[x]]
    
    info_gwas$london.alternate <- rowMeans(cbind(info_gwas$london_pre.alternate, info_gwas$london_during.alternate, info_gwas$london_post.alternate))
    info_gwas$denmark.alternate <- rowMeans(cbind(info_gwas$denmark_pre.alternate, info_gwas$denmark_post.alternate))
    info_gwas$alternate <- rowMeans(cbind(info_gwas$london.alternate, info_gwas$denmark.alternate))

    # let's do exons
    info_exon$london.alternate <- rowMeans(cbind(info_exon$london_pre.alternate, info_exon$london_during.alternate, info_exon$london_post.alternate))
    info_exon$denmark.alternate <- rowMeans(cbind(info_exon$denmark_pre.alternate, info_exon$denmark_post.alternate))
    info_exon$alternate <- rowMeans(cbind(info_exon$london.alternate, info_exon$denmark.alternate))

    # let's do neutral
    info_neut$london.alternate <- rowMeans(cbind(info_neut$london_pre.alternate, info_neut$london_during.alternate, info_neut$london_post.alternate))
    info_neut$denmark.alternate <- rowMeans(cbind(info_neut$denmark_pre.alternate, info_neut$denmark_post.alternate))
    info_neut$alternate <- rowMeans(cbind(info_neut$london.alternate, info_neut$denmark.alternate))

    # add in maf as well as the alternate allele frequency
    info_gwas$maf <- apply(cbind(info_gwas$alternate, 1-info_gwas$alternate), 1, min)
    info_exon$maf <- apply(cbind(info_exon$alternate, 1-info_exon$alternate), 1, min)
    info_neut$maf <- apply(cbind(info_neut$alternate, 1-info_neut$alternate), 1, min)
    
    gwas_table[[x]] =  info_gwas
    exon_table[[x]] = info_exon 
    neutral_table[[x]] = info_neut
    
    }




########
# Calculate Fst
########

for (x in permute_range){
    
    info_gwas = gwas_table[[x]]
    info_exon = exon_table[[x]]
    info_neut= neutral_table[[x]]

    for (i in 1:nrow(design)) {
      n2=paste(design$pop[i],design$time[i],sep="_")
      # expected heterozygosity
      info_gwas[[paste(n2,"ehet",sep=".")]] <- 2*info_gwas[[paste(n2,"alternate",sep=".")]]*(1-info_gwas[[paste(n2,"alternate",sep=".")]])
      info_exon[[paste(n2,"ehet",sep=".")]] <- 2*info_exon[[paste(n2,"alternate",sep=".")]]*(1-info_exon[[paste(n2,"alternate",sep=".")]])
      info_neut[[paste(n2,"ehet",sep=".")]] <- 2*info_neut[[paste(n2,"alternate",sep=".")]]*(1-info_neut[[paste(n2,"alternate",sep=".")]])
    }; rm(i, n2)
    # pairwise means, pairwise expectations, and Fst
    for (i in 1:nrow(design2)) {
      n2=paste(design2$pop[i],design2$time1[i],design2$time2[i],sep="_")
      info_gwas[[paste(n2,"mean",sep=".")]] <- rowMeans(cbind(info_gwas[[paste(design2$pop[i],"_",design2$time1[i],".alternate",sep="")]], info_gwas[[paste(design2$pop[i],"_",design2$time2[i],".alternate",sep="")]]))
      info_gwas[[paste(n2,"ehet",sep=".")]] <- 2*info_gwas[[paste(n2,"mean",sep=".")]]*(1-info_gwas[[paste(n2,"mean",sep=".")]])
      info_gwas[[paste(n2,"fst",sep=".")]] <- (info_gwas[[paste(n2,"ehet",sep=".")]] - (info_gwas[[paste(design2$pop[i],"_",design2$time1[i],".ehet",sep="")]] + info_gwas[[paste(design2$pop[i],"_",design2$time2[i],".ehet",sep="")]])/2)/info_gwas[[paste(n2,"ehet",sep=".")]]

      info_exon[[paste(n2,"mean",sep=".")]] <- rowMeans(cbind(info_exon[[paste(design2$pop[i],"_",design2$time1[i],".alternate",sep="")]], info_exon[[paste(design2$pop[i],"_",design2$time2[i],".alternate",sep="")]]))
      info_exon[[paste(n2,"ehet",sep=".")]] <- 2*info_exon[[paste(n2,"mean",sep=".")]]*(1-info_exon[[paste(n2,"mean",sep=".")]])
      info_exon[[paste(n2,"fst",sep=".")]] <- (info_exon[[paste(n2,"ehet",sep=".")]] - (info_exon[[paste(design2$pop[i],"_",design2$time1[i],".ehet",sep="")]] + info_exon[[paste(design2$pop[i],"_",design2$time2[i],".ehet",sep="")]])/2)/info_exon[[paste(n2,"ehet",sep=".")]]

      info_neut[[paste(n2,"mean",sep=".")]] <- rowMeans(cbind(info_neut[[paste(design2$pop[i],"_",design2$time1[i],".alternate",sep="")]], info_neut[[paste(design2$pop[i],"_",design2$time2[i],".alternate",sep="")]]))
      info_neut[[paste(n2,"ehet",sep=".")]] <- 2*info_neut[[paste(n2,"mean",sep=".")]]*(1-info_neut[[paste(n2,"mean",sep=".")]])
      info_neut[[paste(n2,"fst",sep=".")]] <- (info_neut[[paste(n2,"ehet",sep=".")]] - (info_neut[[paste(design2$pop[i],"_",design2$time1[i],".ehet",sep="")]] + info_neut[[paste(design2$pop[i],"_",design2$time2[i],".ehet",sep="")]])/2)/info_neut[[paste(n2,"ehet",sep=".")]]

    }; rm(i,n2)
    info_gwas -> info_gwas_trim; rm(info_gwas)
    info_neut -> info_neut_trim; rm(info_neut)
    info_exon -> info_exon_trim; rm(info_exon)
    
    gwas_table[[x]] =  info_gwas_trim
    exon_table[[x]] = info_exon_trim 
    neutral_table[[x]] = info_neut_trim
    
    }

######## 
# Filter for 10 individual per time point
######## 
# Report number of sites before filtering

for (x in permute_range){
    info_gwas_trim = gwas_table[[x]]
    info_exon_trim = exon_table[[x]]
    info_neut_trim= neutral_table[[x]]
    
    nrow(info_exon_trim); nrow(info_gwas_trim); nrow(info_neut_trim)
    min_n <- 10
    ######## 
    attach(info_exon_trim)
    info_exon_trim <- subset(info_exon_trim, apply(X = cbind(london_during.called, london_pre.called, london_post.called, denmark_pre.called, denmark_post.called), 1, min) >= min_n)
    detach(info_exon_trim)

    attach(info_neut_trim)
    info_neut_trim <- subset(info_neut_trim, apply(X = cbind(london_during.called, london_pre.called, london_post.called, denmark_pre.called, denmark_post.called), 1, min) >= min_n)
    detach(info_neut_trim)

    attach(info_gwas_trim)
    info_gwas_trim <- subset(info_gwas_trim, apply(X = cbind(london_during.called, london_pre.called, london_post.called, denmark_pre.called, denmark_post.called), 1, min) >= min_n)
    detach(info_gwas_trim)
    ######## 
    # Report number of sites after filtering
    nrow(info_exon_trim); nrow(info_gwas_trim); nrow(info_neut_trim)
    info_gwas <- info_gwas_trim; rm(info_gwas_trim)
    info_exon <- info_exon_trim; rm(info_exon_trim)
    info_neut <- info_neut_trim; rm(info_neut_trim)
    ######## 
    gwas_table[[x]] =  info_gwas
    exon_table[[x]] = info_exon 
    neutral_table[[x]] = info_neut
    
    }

########
# Get difference between time points, then p-values and candidate loci
########
# differences between time points
for (x in permute_range){
    info_gwas = gwas_table[[x]]
    info_exon = exon_table[[x]]
    info_neut= neutral_table[[x]]

    info_exon$delta_L13 <- info_exon$london_post.alternate-info_exon$london_pre.alternate
    info_exon$delta_L12 <- info_exon$london_during.alternate-info_exon$london_pre.alternate
    info_exon$delta_L23 <- info_exon$london_post.alternate-info_exon$london_during.alternate
    info_exon$delta_D13 <- info_exon$denmark_post.alternate-info_exon$denmark_pre.alternate

    info_neut$delta_L13 <- info_neut$london_post.alternate-info_neut$london_pre.alternate
    info_neut$delta_L12 <- info_neut$london_during.alternate-info_neut$london_pre.alternate
    info_neut$delta_L23 <- info_neut$london_post.alternate-info_neut$london_during.alternate
    info_neut$delta_D13 <- info_neut$denmark_post.alternate-info_neut$denmark_pre.alternate

    info_gwas$delta_L13 <- info_gwas$london_post.alternate-info_gwas$london_pre.alternate
    info_gwas$delta_L12 <- info_gwas$london_during.alternate-info_gwas$london_pre.alternate
    info_gwas$delta_L23 <- info_gwas$london_post.alternate-info_gwas$london_during.alternate
    info_gwas$delta_D13 <- info_gwas$denmark_post.alternate-info_gwas$denmark_pre.alternate
    
    gwas_table[[x]] =  info_gwas
    exon_table[[x]] = info_exon 
    neutral_table[[x]] = info_neut
    
    }

########
# trim out some of the columns: leaving just Fst, site, maf
# site information: 1:5
# alternate and minor allele frequency: 7, 9,11,13,15, 19, 20
# fst: 28, 31, 34, 37
# delta_AF: 38:41
########
for (x in permute_range){
    info_gwas = gwas_table[[x]]
    info_exon = exon_table[[x]]
    info_neut= neutral_table[[x]]

    #info_exon <- info_exon[,c(1:5, 7, 9, 11, 13, 15, 19, 20, 28, 31, 34, 37, 38:41)] missing site check column for some reason
    info_exon <- info_exon[,c(1:5, 6, 8, 10, 12, 14, 18, 19, 27, 30, 33, 36, 37:40)]
    info_neut <- info_neut[,c(1:5, 7, 9, 11, 13, 15, 19, 20, 28, 31, 34, 37, 38:41)]
    info_gwas <- info_gwas[,c(1:5, 7, 9, 11, 13, 15, 19, 20, 28, 31, 34, 37, 38:41)]
    
    gwas_table[[x]] =  info_gwas
    exon_table[[x]] = info_exon 
    neutral_table[[x]] = info_neut
    
    }
########

########
# Setup to get p-values
########
# group functional sites together
########
info_table <- list()
pvals_table <- list()
pvals_neut_table <- list()

for (x in permute_range){
    info_gwas = gwas_table[[x]]
    info_exon = exon_table[[x]]
    info_neut= neutral_table[[x]]
    
    info_neut$type <- "neutral"
    info_exon$type <- "exon"
    info_gwas$type <- "gwas"

    rbind(info_exon, info_gwas) -> info #; rm(info_exon, info_gwas)
    ########

    ########
    # Filter based on minor allele frequency
    ########
    info_neut <- subset(info_neut, info_neut$maf > min_maf)
    info <- subset(info, info$maf > min_maf)
    ########

    pvals <- info[,1:5]
    neut_pvals <- info_neut[,1:5]
    
    info_table[[x]] = info
    neutral_table[[x]] = info_neut
    pvals_table[[x]] = pvals
    pvals_neut_table[[x]] = neut_pvals
    
    }

########
# calculate p-values using method
########
#library(parallel)

calc_pval_L13 <- function(site) {
  input_maf <- info$maf[site]
  input_fst <- info$london_post_pre.fst[site]
  neut_maf <- info_neut$maf
  neut_fst <- info_neut$london_post_pre.fst
  
  if (method == "bins") {
    tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
    tmp_fst <- neut_fst[neut_maf > tmp_bin[1] & neut_maf <= tmp_bin[2]]
  }
  if (method == "sliding") {
    tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
    tmp_fst <- neut_fst[neut_maf > tmp_bin[3] & neut_maf <= tmp_bin[4]]
  }

  tmp_percentiles <- ecdf(tmp_fst)
  tmp_p <- 1-tmp_percentiles(input_fst)
  return(tmp_p)
}
#pvals$L13.pval_split <- do.call("c",lapply(1:nrow(info), FUN = calc_pval_L13))
calc_pval_L12 <- function(site) {
  input_maf <- info$maf[site]
  input_fst <- info$london_during_pre.fst[site]
  neut_maf <- info_neut$maf
  neut_fst <- info_neut$london_during_pre.fst
  
  if (method == "bins") {
    tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
    tmp_fst <- neut_fst[neut_maf > tmp_bin[1] & neut_maf <= tmp_bin[2]]
  }
  if (method == "sliding") {
    tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
    tmp_fst <- neut_fst[neut_maf > tmp_bin[3] & neut_maf <= tmp_bin[4]]
  }
  
  tmp_percentiles <- ecdf(tmp_fst)
  tmp_p <- 1-tmp_percentiles(input_fst)
  return(tmp_p)
}
#pvals$L12.pval_split <- do.call("c",lapply(1:nrow(info), FUN = calc_pval_L12))
calc_pval_L23 <- function(site) {
  input_maf <- info$maf[site]
  input_fst <- info$london_during_post.fst[site]
  neut_maf <- info_neut$maf
  neut_fst <- info_neut$london_during_post.fst
  
  if (method == "bins") {
    tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
    tmp_fst <- neut_fst[neut_maf > tmp_bin[1] & neut_maf <= tmp_bin[2]]
  }
  if (method == "sliding") {
    tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
    tmp_fst <- neut_fst[neut_maf > tmp_bin[3] & neut_maf <= tmp_bin[4]]
  }
  
  tmp_percentiles <- ecdf(tmp_fst)
  tmp_p <- 1-tmp_percentiles(input_fst)
  return(tmp_p)
}
#pvals$L23.pval_split <- do.call("c",lapply(1:nrow(info), FUN = calc_pval_L23))

calc_pval_D13 <- function(site) {
  input_maf <- info$maf[site]
  input_fst <- info$denmark_post_pre.fst[site]
  neut_maf <- info_neut$maf
  neut_fst <- info_neut$denmark_post_pre.fst
  
  if (method == "bins") {
    tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
    tmp_fst <- neut_fst[neut_maf > tmp_bin[1] & neut_maf <= tmp_bin[2]]
  }
  if (method == "sliding") {
    tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
    tmp_fst <- neut_fst[neut_maf > tmp_bin[3] & neut_maf <= tmp_bin[4]]
  }
  
  tmp_percentiles <- ecdf(tmp_fst)
  tmp_p <- 1-tmp_percentiles(input_fst)
  return(tmp_p)
}
#pvals$D13.pval_split <- do.call("c",lapply(1:nrow(info), FUN = calc_pval_D13))

{
  calc_neut_pval_L13 <- function(site) {
    input_maf <- info_neut$maf[site]
    input_fst <- info_neut$london_post_pre.fst[site]
    neut_maf <- info_neut$maf
    neut_fst <- info_neut$london_post_pre.fst
    
    if (method == "bins") {
      tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
      tmp_fst <- neut_fst[neut_maf > tmp_bin[1] & neut_maf <= tmp_bin[2]]
    }
    if (method == "sliding") {
      tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
      tmp_fst <- neut_fst[neut_maf > tmp_bin[3] & neut_maf <= tmp_bin[4]]
    }
    tmp_percentiles <- ecdf(tmp_fst)
    tmp_p <- 1-tmp_percentiles(input_fst)
    return(tmp_p)
  }
  #neut_pvals$L13.pval_split <- do.call("c",lapply(1:nrow(info_neut), FUN = calc_neut_pval_L13))
  calc_neut_pval_L12 <- function(site) {
    input_maf <- info_neut$maf[site]
    input_fst <- info_neut$london_during_pre.fst[site]
    neut_maf <- info_neut$maf
    neut_fst <- info_neut$london_during_pre.fst
    
    if (method == "bins") {
      tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
      tmp_fst <- neut_fst[neut_maf > tmp_bin[1] & neut_maf <= tmp_bin[2]]
    }
    if (method == "sliding") {
      tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
      tmp_fst <- neut_fst[neut_maf > tmp_bin[3] & neut_maf <= tmp_bin[4]]
    }
    
    tmp_percentiles <- ecdf(tmp_fst)
    tmp_p <- 1-tmp_percentiles(input_fst)
    return(tmp_p)
  }
  #neut_pvals$L12.pval_split <- do.call("c",lapply(1:nrow(info_neut), FUN = calc_neut_pval_L12))
  calc_neut_pval_L23 <- function(site) {
    input_maf <- info_neut$maf[site]
    input_fst <- info_neut$london_during_post.fst[site]
    neut_maf <- info_neut$maf
    neut_fst <- info_neut$london_during_post.fst
    
    if (method == "bins") {
      tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
      tmp_fst <- neut_fst[neut_maf > tmp_bin[1] & neut_maf <= tmp_bin[2]]
    }
    if (method == "sliding") {
      tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
      tmp_fst <- neut_fst[neut_maf > tmp_bin[3] & neut_maf <= tmp_bin[4]]
    }
    
    tmp_percentiles <- ecdf(tmp_fst)
    tmp_p <- 1-tmp_percentiles(input_fst)
    return(tmp_p)
  }
  #neut_pvals$L23.pval_split <- do.call("c",lapply(1:nrow(info_neut), FUN = calc_neut_pval_L23))
  
  calc_neut_pval_D13 <- function(site) {
    input_maf <- info_neut$maf[site]
    input_fst <- info_neut$denmark_post_pre.fst[site]
    neut_maf <- info_neut$maf
    neut_fst <- info_neut$denmark_post_pre.fst
    
    if (method == "bins") {
      tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
      tmp_fst <- neut_fst[neut_maf > tmp_bin[1] & neut_maf <= tmp_bin[2]]
    }
    if (method == "sliding") {
      tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
      tmp_fst <- neut_fst[neut_maf > tmp_bin[3] & neut_maf <= tmp_bin[4]]
    }
    
    tmp_percentiles <- ecdf(tmp_fst)
    tmp_p <- 1-tmp_percentiles(input_fst)
    return(tmp_p)
  }
  #neut_pvals$D13.pval_split <- do.call("c",lapply(1:nrow(info_neut), FUN = calc_neut_pval_D13))
  
}
#rm(list=ls(pattern="calc_neut_pval")); rm(list=ls(pattern="calc_pval_"))


########
#rm(list=ls(pattern="calc_gwas_pval")); rm(list=ls(pattern="calc_neut_pval")); rm(list=ls(pattern="calc_exon_pval"))


library(parallel)

for (x in permute_range){
    pvals = pvals_table[[x]]
    info=info_table[[x]]
    neut_pvals = pvals_neut_table[[x]]
    info_neut=neutral_table[[x]]
    
    pvals$L13.pval_split <- do.call("c",mclapply(1:nrow(info), FUN = calc_pval_L13))
    pvals$L12.pval_split <- do.call("c",mclapply(1:nrow(info), FUN = calc_pval_L12))
    pvals$L23.pval_split <- do.call("c",mclapply(1:nrow(info), FUN = calc_pval_L23))
    pvals$D13.pval_split <- do.call("c",mclapply(1:nrow(info), FUN = calc_pval_D13))

    neut_pvals$L13.pval_split <- do.call("c",mclapply(1:nrow(info_neut), FUN = calc_neut_pval_L13))
    neut_pvals$L12.pval_split <- do.call("c",mclapply(1:nrow(info_neut), FUN = calc_neut_pval_L12))
    neut_pvals$L23.pval_split <- do.call("c",mclapply(1:nrow(info_neut), FUN = calc_neut_pval_L23))
    neut_pvals$D13.pval_split <- do.call("c",mclapply(1:nrow(info_neut), FUN = calc_neut_pval_D13))
    
    pvals_table[[x]] = pvals
    pvals_neut_table[[x]] = neut_pvals
    }

library(ggplot2); library(patchwork)#; library(metap)
#load("./DATA/DATA_part1.RData")
## Print some metadata
print(paste("method: ", method, sep=""))
if (method == 'sliding') {print(paste("minimum window: 5%, minimum neutral sites: ", min_neutral_nsites, sep=""))}
print(paste("drop individuals missing more than X proportion of sites of data: ", drop_samples_missing, sep=""))
print(paste("minimum number of samples per population per gt calls: ", min_n, sep=""))
print(paste("minimum minor allele frequency: ", min_maf, sep=""))
print(paste("numer of candidate loci: ", nrow(info), sep=""))
print(paste("numer of neutral loci: ", nrow(info_neut), sep=""))

#**Minor allele frequency spectrum** of analyzed variants (after 5% maf and 10+ of samples per time point). 
#{r, echo=F}


#{r, echo=F}
library(ggplot2); library(patchwork)

for (x in permute_range){
    pvals = pvals_table[[x]]
    info=info_table[[x]]
    neut_pvals = pvals_neut_table[[x]]
    info_neut=neutral_table[[x]]
    
    info <- info[,1:21]
    info_neut <- info_neut[,1:21]

    info <- cbind(info, pvals[,6:9])
    colnames(info)[22:25] <- c("L13.pval", "L12.pval", "L23.pval", "D13.pval")
    info_neut <- cbind(info_neut, neut_pvals[,6:9])
    colnames(info_neut)[22:25] <- c("L13.pval", "L12.pval", "L23.pval", "D13.pval")
    
    info_table[[x]] = info
    neutral_table[[x]] = info_neut
    }

library(ggplot2); library(patchwork)
## we need to get a matrix that has the maf window, the percentile, and the degree of enrichment
tmp_bin_breaks <- c(0.1,.2,.3,.4,0.5)
tmp_bins <-  matrix(ncol=2, nrow=length(tmp_bin_breaks)); tmp_bins[,1] <- c(0,tmp_bin_breaks[-length(tmp_bin_breaks)]); tmp_bins[,2] <- c(tmp_bin_breaks); rm(tmp_bin_breaks)
print("number of sites in each maf bin")
for (i in 1:nrow(tmp_bins)) {
  print(paste(100*tmp_bins[i,1], "% to ", 100*tmp_bins[i,2], "%: ", 
              nrow(info[info$maf > tmp_bins[i,1] & info$maf <= tmp_bins[i,2],]), " sites", sep=""))
}
res <- NULL

London_p_maf10_e05 <- list()
London_p_maf10_e01 <- list()
London_p_maf30_e01 <- list()


for (x in permute_range){
    pvals = pvals_table[[x]]
    info=info_table[[x]]
    neut_pvals = pvals_neut_table[[x]]
    info_neut=neutral_table[[x]]
    
    for (tmp_enrich in c(0.001, 0.005, seq(0.01,0.1,0.01), seq(0.1,0.2,0.05))) {
      for (tmp_bin in 1:nrow(tmp_bins)) {
        for (tmp_pop in c("London", "Denmark", "L12", "L23")) {
          tmp_data <- info[info$maf > tmp_bins[tmp_bin,1] & info$maf <= tmp_bins[tmp_bin,2],]
        tmp_res <- as.data.frame(matrix(ncol=6, nrow=1)); colnames(tmp_res) <- c('maf_low', 'maf_high', 'population', 'enrichment', 'observed', 'expected')
        tmp_res$population <- tmp_pop
        tmp_res$maf_low <- tmp_bins[tmp_bin,1]; tmp_res$maf_high <- tmp_bins[tmp_bin,2]
        tmp_res$enrichment <- tmp_enrich  
        if (tmp_pop == "London") {tmp_res$observed <- sum(tmp_data$L13.pval < tmp_enrich)}
        if (tmp_pop == "Denmark") {tmp_res$observed <- sum(tmp_data$D13.pval < tmp_enrich)}
        if (tmp_pop == "L12") {tmp_res$observed <- sum(tmp_data$L12.pval < tmp_enrich)}
        if (tmp_pop == "L23") {tmp_res$observed <- sum(tmp_data$L23.pval < tmp_enrich)}
        tmp_res$expected <- nrow(tmp_data)*tmp_enrich
        rbind(res, tmp_res) -> res; rm(tmp_res, tmp_data)
        }; rm(tmp_pop)
      }; rm(tmp_bin)
    }; rm(tmp_enrich)
    res$maf_bin <- paste(100*res$maf_low, "% to ", 100*res$maf_high, "%", sep="")
    res$fc <- log2(res$observed/res$expected)
    #print("binomial test: 95%")
    tmp <- res[res$population == "London" & res$maf_high > 0.1 & res$enrichment == 0.05,]; tmp
    london_maf10_e05 <- binom.test(sum(tmp$observed), n = (1/.05)*sum(tmp$expected), p=0.05); log2(sum(tmp$observed)/sum(tmp$expected))
    London_p_maf10_e05[[x]] = london_maf10_e05$p.value
    #print("binomial test: 99%")
    tmp <- res[res$population == "London" & res$maf_high > 0.1 & res$enrichment == 0.01,]; tmp
    london_maf10_e01 <- binom.test(sum(tmp$observed), n = (1/.01)*sum(tmp$expected), p=0.01); log2(sum(tmp$observed)/sum(tmp$expected)); sum(tmp$observed)/sum(tmp$expected)
    London_p_maf10_e01[[x]] = london_maf10_e01$p.value
    #print("binomial test: 99%, 30%+")
    tmp <- res[res$population == "London" & res$maf_high >= 0.3 & res$enrichment == 0.01,]; tmp
    london_maf30_e01 <- binom.test(sum(tmp$observed), n = (1/.01)*sum(tmp$expected), p=0.01); log2(sum(tmp$observed)/sum(tmp$expected)); sum(tmp$observed)/sum(tmp$expected); rm(tmp)
    London_p_maf30_e01[[x]] = london_maf30_e01$p.value
}
write.csv(London_p_maf10_e05,"London_p_maf10_e05.txt",row.names=FALSE)
write.csv(London_p_maf10_e01,"London_p_maf10_e01.txt",row.names=FALSE)
write.csv(London_p_maf30_e01,"London_p_maf30_e01.txt",row.names=FALSE)


