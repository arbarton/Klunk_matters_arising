#######
# Define ML estimator
#######

#Estimate the allele frequency
likelihood <- function(p, data){
  gt.freq <- c((1-p)^2, 2*p*(1-p), p^2)
  ll <- sum(log(rowSums(t(t(data)*gt.freq))))
  return(ll)
}

estimate_af_ml<-function(d){
  n<-NROW(d)
  af<-rep(NA, n)
  for(i in 1:n){
    gl<-d[i,,drop=FALSE]
    gl <- gl[!is.na(gl)]
    gl <- matrix(gl, ncol=3, byrow=T)
    gl <- gl[!is.na(gl[,1]),,drop=FALSE]
    opt <- optimize(likelihood, interval=c(0,1), maximum=TRUE, data=gl)
    af[i] <- opt$maximum
  }
  return(af)
}

########
# Following code adapted from original Klunk et al. pipeline
########

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

########

min_n=10; min_maf=0.05
drop_samples_missing=.5

method='sliding'

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
######## 

######## 
# Get trimmed aDNA frequencies
######## 
# Info files with the number and alternate allele freq in each population
# mac >= 3; biallelic; minQ >= 30
# exons filtered for sites within bed file; rest mapped to the bed file
######## 
info_gwas <- fread("./DATA/genoliks/immune.frq")[,-c(3:4)]
info_neut <- fread("./DATA/genoliks/neutral.frq")[,-c(3:4)]
info_exon <- fread("./DATA/genoliks/exon.frq")[,-c(3:4)]

colnames(info_gwas) <- c("chr", "pos", "ref", "alt")
colnames(info_neut) <- c("chr", "pos", "ref", "alt")
colnames(info_exon) <- c("chr", "pos", "ref", "alt")

info_gwas$ref <- gsub(":.*", "",info_gwas$ref); info_gwas$alt <- gsub(":.*", "",info_gwas$alt)
info_neut$ref <- gsub(":.*", "",info_neut$ref); info_neut$alt <- gsub(":.*", "",info_neut$alt)
info_exon$ref <- gsub(":.*", "",info_exon$ref); info_exon$alt <- gsub(":.*", "",info_exon$alt)

paste(info_gwas$chr,info_gwas$pos,sep="_") -> info_gwas$site
paste(info_neut$chr,info_neut$pos,sep="_") -> info_neut$site
paste(info_exon$chr,info_exon$pos,sep="_") -> info_exon$site


for (i in c(3:nrow(design),1:2)) {
  # read in exon data
  name=paste("./DATA/genoliks/genolik.exons_",design$pop[i],"_", design$time[i],".genolik",sep="")
  n2=paste(design$pop[i],design$time[i],sep="_"); d <- as.data.frame(fread(name)); rm(name)
  snames=read.delim(paste("./DATA/genoliks/SampleNames/", design$pop[i], "_", design$time[i],"_exons.012.indv",sep=""), header=F)
  ncol(d)-2 == 3*nrow(snames)
  if (i == 3) {paste(d[,1],d[,2],sep="_") -> info_exon$site_check}
  # remove site. change missing data to NA. get number of samples. 
  d <- d[,-(1:2)] ; d[d == -1] <- NA; n <- ncol(d)/3
  # report data on missingness per sample
  assign(paste("missing",n2, "exon", sep="_"), value = (colSums(is.na(d))/nrow(d))[seq(1,3*n,3)])
  # remove samples missing too much data
  keep <- which(get(paste("missing",n2, "exon", sep="_")) <= drop_samples_missing)
  assign(paste("keep",n2,"exon",sep="_"), value=snames[keep,]); write.table(snames[keep,], paste("./DATA/genoliks/SampleNames/keep",n2,"exon",sep="_"), row.names=F, col.names=F, sep="\t", quote=F)
  k2 <- as.data.frame(matrix(nrow=nrow(d), ncol=1)); for (k in keep) {k2 <- cbind(k2, d[,(k*3-2):(k*3)])}
  d <- k2[,-1]; rm(keep, k2, k); n <- ncol(d)/3
  # get the allele frequency, and the number of individuals it was called from
  info_exon[[paste(n2,"alternate",sep=".")]] <- estimate_af_ml(d)
  info_exon[[paste(n2,"called",sep=".")]] <- (n*3-rowSums(is.na(d)))/3
  # Report number of samples included 
  assign(paste("samples",n2, "exon", sep="_"), n)
  # cleanup
  rm(d, n)
  
  # read in gwas data
  name=paste("./DATA/genoliks/genolik.gwas_",design$pop[i],"_", design$time[i],".genolik",sep="")
  n2=paste(design$pop[i],design$time[i],sep="_"); d <- as.data.frame(fread(name)); rm(name)
  snames=read.delim(paste("./DATA/genoliks/SampleNames/", design$pop[i], "_", design$time[i],"_immune.012.indv",sep=""), header=F)
  ncol(d)-2 == 3*nrow(snames)
  if (i == 3) {paste(d[,1],d[,2],sep="_") -> info_gwas$site_check}
  # remove site. change missing data to NA. get number of samples. 
  d <- d[,-(1:2)] ; d[d == -1] <- NA; n <- ncol(d)/3
  # report data on missingness per sample
  assign(paste("missing",n2, "gwas", sep="_"), value = (colSums(is.na(d))/nrow(d))[seq(1,3*n,3)])
  # remove samples missing too much data
  keep <- which(get(paste("missing",n2, "gwas", sep="_")) <= drop_samples_missing)
  assign(paste("keep",n2,"gwas",sep="_"), value=snames[keep,]); write.table(snames[keep,], paste("./DATA/genoliks/SampleNames/keep",n2,"gwas",sep="_"), row.names=F, col.names=F, sep="\t", quote=F)
  k2 <- as.data.frame(matrix(nrow=nrow(d), ncol=1)); for (k in keep) {k2 <- cbind(k2, d[,(k*3-2):(k*3)])}
  d <- k2[,-1]; rm(keep, k2, k); n <- ncol(d)/3
  # get the allele frequency, and the number of individuals it was called from
  info_gwas[[paste(n2,"alternate",sep=".")]] <- estimate_af_ml(d)
  info_gwas[[paste(n2,"called",sep=".")]] <- (ncol(d)-rowSums(is.na(d)))/3
  # Report number of samples included 
  assign(paste("samples",n2, "gwas", sep="_"), n)
  # cleanup
  rm(d, n)
  
  # read in neutral data
  name=paste("./DATA/genoliks/genolik.neutral_",design$pop[i],"_", design$time[i],".genolik",sep="")
  n2=paste(design$pop[i],design$time[i],sep="_"); d <- as.data.frame(fread(name)); rm(name)
  snames=read.delim(paste("./DATA/genoliks/SampleNames/", design$pop[i], "_", design$time[i],"_neutral.012.indv",sep=""), header=F)
  ncol(d)-2 == 3*nrow(snames)
  if (i == 3) {paste(d[,1],d[,2],sep="_") -> info_neut$site_check}
  # remove site. change missing data to NA. get number of samples.
  d <- d[,-(1:2)] ; d[d == -1] <- NA; n <- ncol(d)/3
  # report data on missingness per sample
  assign(paste("missing",n2, "neutral", sep="_"), value = (colSums(is.na(d))/nrow(d))[seq(1,3*n,3)])
  # remove samples missing too much data
  keep <- which(get(paste("missing",n2, "neutral", sep="_")) <= drop_samples_missing)
  assign(paste("keep",n2,"neutral",sep="_"), value=snames[keep,]); write.table(snames[keep,], paste("./DATA/genoliks/SampleNames/keep",n2,"neutral",sep="_"), row.names=F, col.names=F, sep="\t", quote=F)
  k2 <- as.data.frame(matrix(nrow=nrow(d), ncol=1)); for (k in keep) {k2 <- cbind(k2, d[,(k*3-2):(k*3)])}
  d <- k2[,-1]; rm(keep, k2, k); n <- ncol(d)/3
  # get the allele frequency, and the number of individuals it was called from
  info_neut[[paste(n2,"alternate",sep=".")]] <- estimate_af_ml(d)
  info_neut[[paste(n2,"called",sep=".")]] <- (ncol(d)-rowSums(is.na(d)))/3
  # Report number of samples included 
  assign(paste("samples",n2, "neutral", sep="_"), n)
  # cleanup
  rm(d, n)
}; rm(i,n2)
#save.image("./r01.summarized_allele_frequencies.RData")
## saved for John & Matthias

sum(info_gwas$site == info_gwas$site_check) == nrow(info_gwas)
sum(info_exon$site == info_exon$site_check) == nrow(info_exon)
sum(info_neut$site == info_neut$site_check) == nrow(info_neut)

for (f in ls(pattern="missing")) {print(paste(f,sum(get(f) < 0.5), sep=": "))}; rm(f)


########

########
# calculate alternate and minor allele frequency
# calculate a "London" frequency and a "Denmark" frequency as the mean between time points within that population
# NOTE: THIS TAKES US UP TO THE FIRST 20 COLUMNS
########
# let's do gwas first
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


########

########
# Calculate Fst
########
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
######## 

######## 
# Filter for 10 individual per time point
######## 
# Report number of sites before filtering
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

########
# Get difference between time points, then p-values and candidate loci
########
# differences between time points
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

########
# trim out some of the columns: leaving just Fst, site, maf
# site information: 1:5
# alternate and minor allele frequency: 7, 9,11,13,15, 19, 20
# fst: 28, 31, 34, 37
# delta_AF: 38:41
########
info_exon <- info_exon[,c(1:5, 7, 9, 11, 13, 15, 19, 20, 28, 31, 34, 37, 38:41)]
info_neut <- info_neut[,c(1:5, 7, 9, 11, 13, 15, 19, 20, 28, 31, 34, 37, 38:41)]
info_gwas <- info_gwas[,c(1:5, 7, 9, 11, 13, 15, 19, 20, 28, 31, 34, 37, 38:41)]
########

########
# Setup to get p-values
########
# group functional sites together
########
info_neut$type <- "neutral"
info_exon$type <- "exon"
info_gwas$type <- "gwas"

rbind(info_exon, info_gwas) -> info; rm(info_exon, info_gwas)
ML_info <- info
ML_neut <- info_neut
########

########
# Filter based on minor allele frequency
########
info_neut <- subset(info_neut, info_neut$maf > min_maf)
info <- subset(info, info$maf > min_maf)
########

pvals <- info[,1:5]
neut_pvals <- info_neut[,1:5]
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
pvals$L13.pval_split <- do.call("c",lapply(1:nrow(info), FUN = calc_pval_L13))
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
pvals$L12.pval_split <- do.call("c",lapply(1:nrow(info), FUN = calc_pval_L12))
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
pvals$L23.pval_split <- do.call("c",lapply(1:nrow(info), FUN = calc_pval_L23))

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
pvals$D13.pval_split <- do.call("c",lapply(1:nrow(info), FUN = calc_pval_D13))

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
  neut_pvals$L13.pval_split <- do.call("c",lapply(1:nrow(info_neut), FUN = calc_neut_pval_L13))
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
  neut_pvals$L12.pval_split <- do.call("c",lapply(1:nrow(info_neut), FUN = calc_neut_pval_L12))
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
  neut_pvals$L23.pval_split <- do.call("c",lapply(1:nrow(info_neut), FUN = calc_neut_pval_L23))
  
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
  neut_pvals$D13.pval_split <- do.call("c",lapply(1:nrow(info_neut), FUN = calc_neut_pval_D13))
  
}
rm(list=ls(pattern="calc_neut_pval")); rm(list=ls(pattern="calc_pval_"))


########
rm(list=ls(pattern="calc_gwas_pval")); rm(list=ls(pattern="calc_neut_pval")); rm(list=ls(pattern="calc_exon_pval"))
########

library(ggplot2); library(patchwork)#; library(metap)
## Print some metadata
print(paste("method: ", method, sep=""))
if (method == 'sliding') {print(paste("minimum window: 5%, minimum neutral sites: ", min_neutral_nsites, sep=""))}
print(paste("drop individuals missing more than X proportion of sites of data: ", drop_samples_missing, sep=""))
print(paste("minimum number of samples per population per gt calls: ", min_n, sep=""))
print(paste("minimum minor allele frequency: ", min_maf, sep=""))
print(paste("numer of candidate loci: ", nrow(info), sep=""))
print(paste("numer of neutral loci: ", nrow(info_neut), sep=""))


par(mfrow=c(2,3))
tmp <- subset(info, info$type == "exon")
plot(density(apply(cbind(tmp$london_pre.alternate,1-tmp$london_pre.alternate),1,min, na.rm=T)), col='red', main="London before", ylim=c(0,20))
tmp <- subset(info, info$type == "gwas")
lines(density(apply(cbind(tmp$london_pre.alternate,1-tmp$london_pre.alternate),1,min, na.rm=T)), col='green')
tmp <- info_neut
lines(density(apply(cbind(tmp$london_pre.alternate,1-tmp$london_pre.alternate),1,min, na.rm=T)), col='blue')
tmp <- subset(info, info$type == "exon")
plot(density(apply(cbind(tmp$london_post.alternate,1-tmp$london_post.alternate),1,min, na.rm=T)), col='red', main="London after", ylim=c(0,20))
tmp <- subset(info, info$type == "gwas")
lines(density(apply(cbind(tmp$london_post.alternate,1-tmp$london_post.alternate),1,min, na.rm=T)), col='green')
tmp <- info_neut
lines(density(apply(cbind(tmp$london_post.alternate,1-tmp$london_post.alternate),1,min, na.rm=T)), col='blue')
tmp <- subset(info, info$type == "exon")
plot(density(apply(cbind(tmp$london_during.alternate,1-tmp$london_during.alternate),1,min, na.rm=T)), col='red', main="London during", ylim=c(0,20))
tmp <- subset(info, info$type == "gwas")
lines(density(apply(cbind(tmp$london_during.alternate,1-tmp$london_during.alternate),1,min, na.rm=T)), col='green')
tmp <- info_neut
lines(density(apply(cbind(tmp$london_during.alternate,1-tmp$london_during.alternate),1,min, na.rm=T)), col='blue')
tmp <- subset(info, info$type == "exon")
plot(density(apply(cbind(tmp$denmark_pre.alternate,1-tmp$denmark_pre.alternate),1,min, na.rm=T)), col='red', main="denmark before", ylim=c(0,20))
tmp <- subset(info, info$type == "gwas")
lines(density(apply(cbind(tmp$denmark_pre.alternate,1-tmp$denmark_pre.alternate),1,min, na.rm=T)), col='green')
tmp <- info_neut
lines(density(apply(cbind(tmp$denmark_pre.alternate,1-tmp$denmark_pre.alternate),1,min, na.rm=T)), col='blue')
tmp <- subset(info, info$type == "exon")
plot(density(apply(cbind(tmp$denmark_post.alternate,1-tmp$denmark_post.alternate),1,min, na.rm=T)), col='red', main="denmark after", ylim=c(0,20))
tmp <- subset(info, info$type == "gwas")
lines(density(apply(cbind(tmp$denmark_post.alternate,1-tmp$denmark_post.alternate),1,min, na.rm=T)), col='green')
tmp <- info_neut
lines(density(apply(cbind(tmp$denmark_post.alternate,1-tmp$denmark_post.alternate),1,min, na.rm=T)), col='blue')
rm(tmp)
########

### Not controlling for maf appears biases results, so we'll use the pvalues controlling for maf.  

#no maf
# info <- cbind(info, pvals[,10:13])
# colnames(info)[22:25] <- c("L13.pval", "L12.pval", "L23.pval", "D13.pval")
# info_neut <- cbind(info_neut, neut_pvals[,10:13])
# colnames(info_neut)[22:25] <- c("L13.pval", "L12.pval", "L23.pval", "D13.pval")
#with maf
info <- info[,1:21]
info_neut <- info_neut[,1:21]

info <- cbind(info, pvals[,6:9])
colnames(info)[22:25] <- c("L13.pval", "L12.pval", "L23.pval", "D13.pval")
info_neut <- cbind(info_neut, neut_pvals[,6:9])
colnames(info_neut)[22:25] <- c("L13.pval", "L12.pval", "L23.pval", "D13.pval")
#rm(pvals, neut_pvals)


### Enrichment of highly differentiated sites split by maf bracket

## we need to get a matrix that has the maf window, the percentile, and the degree of enrichment
tmp_bin_breaks <- c(0.1,.2,.3,.4,0.5)
tmp_bins <-  matrix(ncol=2, nrow=length(tmp_bin_breaks)); tmp_bins[,1] <- c(0,tmp_bin_breaks[-length(tmp_bin_breaks)]); tmp_bins[,2] <- c(tmp_bin_breaks); rm(tmp_bin_breaks)
print("number of sites in each maf bin")
for (i in 1:nrow(tmp_bins)) {
  print(paste(100*tmp_bins[i,1], "% to ", 100*tmp_bins[i,2], "%: ", 
              nrow(info[info$maf > tmp_bins[i,1] & info$maf <= tmp_bins[i,2],]), " sites", sep=""))
}
res <- NULL
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
}; rm(tmp_enrich, tmp_bins)
res$maf_bin <- paste(100*res$maf_low, "% to ", 100*res$maf_high, "%", sep="")
res$fc <- log2(res$observed/res$expected)
print("binomial test: 95%")
tmp <- res[res$population == "London" & res$maf_high > 0.1 & res$enrichment == 0.05,]; tmp
binom.test(sum(tmp$observed), n = (1/.05)*sum(tmp$expected), p=0.05); log2(sum(tmp$observed)/sum(tmp$expected))
print("binomial test: 99%")
tmp <- res[res$population == "London" & res$maf_high > 0.1 & res$enrichment == 0.01,]; tmp
binom.test(sum(tmp$observed), n = (1/.01)*sum(tmp$expected), p=0.01); log2(sum(tmp$observed)/sum(tmp$expected)); sum(tmp$observed)/sum(tmp$expected)
print("binomial test: 99%, 30%+")
tmp <- res[res$population == "London" & res$maf_high >= 0.3 & res$enrichment == 0.01,]; tmp
binom.test(sum(tmp$observed), n = (1/.01)*sum(tmp$expected), p=0.01); log2(sum(tmp$observed)/sum(tmp$expected)); sum(tmp$observed)/sum(tmp$expected); rm(tmp)
print("binomial test: 95%")
tmp <- res[res$population == "Denmark" & res$maf_high > 0.1 & res$enrichment == 0.05,]; tmp
binom.test(sum(tmp$observed), n = (1/.05)*sum(tmp$expected), p=0.05); log2(sum(tmp$observed)/sum(tmp$expected))
print("binomial test: 99%")
tmp <- res[res$population == "Denmark" & res$maf_high > 0.1 & res$enrichment == 0.01,]; tmp
binom.test(sum(tmp$observed), n = (1/.01)*sum(tmp$expected), p=0.01); log2(sum(tmp$observed)/sum(tmp$expected)); sum(tmp$observed)/sum(tmp$expected)
print("binomial test: 99%, 30%+")
tmp <- res[res$population == "Denmark" & res$maf_high >= 0.3 & res$enrichment == 0.01,]; tmp
binom.test(sum(tmp$observed), n = (1/.01)*sum(tmp$expected), p=0.01); log2(sum(tmp$observed)/sum(tmp$expected)); sum(tmp$observed)/sum(tmp$expected); rm(tmp)
p1 <- ggplot(res[res$population == "London",], aes(x=1-enrichment , y=fc, color=maf_bin)) +
  geom_point(size=2) + geom_line() +
  theme_classic() + ylab("log2(enrichment)") +
  ggtitle("London: pre vs post") + #coord_cartesian(xlim=c(0,0.1), ylim=c(-0.25,3.75)) +
  geom_abline(slope=0, intercept=0, color='dark gray') + theme(legend.position = "bottom") + 
  guides(fill = guide_legend(label.position = "bottom")) + xlab("percentile of neutral variants")
p2 <- ggplot(res[res$population == "Denmark",], aes(x=1-enrichment , y=fc, color=maf_bin)) +
  geom_point(size=2) + geom_line() +
  theme_classic() + ylab("log2(enrichment)") +
  ggtitle("Denmark: pre vs post") + #coord_cartesian(xlim=c(0,0.1), ylim=c(-0.25,3.75)) +
  geom_abline(slope=0, intercept=0, color='dark gray') + theme(legend.position = "bottom") + 
  guides(fill = guide_legend(label.position = "bottom")) + xlab("percentile of neutral variants")
p3 <- ggplot(res[res$population == "L12" & res$maf_high > 0.1,], aes(x=1-enrichment , y=fc, color=maf_bin)) +
  geom_point(size=2) + geom_line() +
  theme_classic() + ylab("log2(enrichment)") +
  ggtitle("London: pre vs during") + #coord_cartesian(xlim=c(0,0.1), ylim=c(-0.25,3.75)) +
  geom_abline(slope=0, intercept=0, color='dark gray') + theme(legend.position = "bottom") + 
  guides(fill = guide_legend(label.position = "bottom")) + xlab("percentile of neutral variants")
p4 <- ggplot(res[res$population == "L23" & res$maf_high > 0.1,], aes(x=1-enrichment , y=fc, color=maf_bin)) +
  geom_point(size=2) + geom_line() +
  theme_classic() + ylab("log2(enrichment)") +
  ggtitle("London: during vs post") + #coord_cartesian(xlim=c(0,0.1), ylim=c(-0.25,3.75)) +
  geom_abline(slope=0, intercept=0, color='dark gray') + theme(legend.position = "bottom") + 
  guides(fill = guide_legend(label.position = "bottom")) + xlab("percentile of neutral variants")
panel_enrichment <- p1 & theme(legend.position = 'bottom')
panel_enrichment
panel_D_enrichment <- p2 & theme(legend.position = 'bottom')
splot_enrichment <- p2 + p3 + p4 + plot_layout(guides="collect") & theme(legend.position = 'bottom'); rm(p1,p2,p3,p4)
splot_enrichment
rm(res, tmp)


### Correlation in allele frequency betwen time points

par(mfrow=c(1,2))
plot(info_neut$london_pre.alternate, info_neut$london_post.alternate, main="London", xlab="before plague", ylab="after plague", frame=F)
# summary(lm(info_neut$london_pre.alternate ~ info_neut$london_post.alternate))
plot(info_neut$denmark_pre.alternate, info_neut$denmark_post.alternate, main="Denmark", xlab="before plague", ylab="after plague", frame=F)
# summary(lm(info_neut$denmark_pre.alternate ~ info_neut$denmark_post.alternate))
# 
# summary(lm(c(info_neut$london_post.alternate,info_neut$london_pre.alternate) ~ 
#              c(info_neut$denmark_post.alternate,info_neut$denmark_pre.alternate)))
# summary(lm(info_neut$london_pre.alternate ~ info_neut$denmark_pre.alternate))


### Because the only evidence of enrichment is for minor allele frequencies greater than 10%, let's only consider those sites when detecting candidate loci. 

par(mfrow=c(1,2))
info <- info[info$maf >= 0.1,]
info_neut <- info_neut[info_neut$maf >= 0.1,]
## replace p=0 with 1e-4 for metanalysis; metap doesn't work withi p=0
info$L13.pval[info$L13.pval == 0] <- 1e-5
info$L12.pval[info$L12.pval == 0] <- 1e-5
info$L23.pval[info$L23.pval == 0] <- 1e-5
info$D13.pval[info$D13.pval == 0] <- 1e-5
info_neut$L12.pval[info_neut$L12.pval == 0] <- 1e-5
info_neut$L23.pval[info_neut$L23.pval == 0] <- 1e-5
info_neut$D13.pval[info_neut$D13.pval == 0] <- 1e-5
# # part 1: get metap for london (jointp_london)
# for (i in 1:nrow(info)) {sumlog(p = cbind(info$L12.pval, info$L23.pval)[i,])$p -> info$jointp_london[i]}
# for (i in 1:nrow(info_neut)) {sumlog(p = cbind(info_neut$L12.pval, info_neut$L23.pval)[i,])$p -> info_neut$jointp_london[i]}
# # t1 <- -log10(info_neut$jointp_london); t2 <- -log10(info$jointp_london); 
# # qqplot(t1, t2, xlab="-log10(p neutral)", ylab="-log10(p)", frame=F, main="London metap"); abline(a=0,b=1,col='red')
# # rm(t1, t2)
# # 
# # part 2: get metap for london + denmark (jointp)
# for (i in 1:nrow(info)) {sumlog(p = cbind(info$jointp_london, info$D13.pval)[i,])$p -> info$jointp[i]}
# for (i in 1:nrow(info_neut)) {sumlog(p = cbind(info_neut$jointp_london, info_neut$D13.pval)[i,])$p -> info_neut$jointp[i]}
# # t1 <- -log10(info_neut$jointp); t2 <- -log10(info$jointp); 
# # qqplot(t1, t2, xlab="-log10(p neutral)", ylab="-log10(p)", frame=F, main="metap: replication cohorts"); abline(a=0,b=1,col='red')
# # rm(t1, t2)
# rm(i)
# 
# # jointp_13
# for (i in 1:nrow(info)) {sumlog(p = cbind(info$L13.pval, info$D13.pval)[i,])$p -> info$jointp_13[i]}
# info$meanfst_13 <- rowMeans(cbind(info$london_post_pre.fst, info$denmark_post_pre.fst))
# info$meanp_13 <- rowMeans(cbind(info$L13.pval, info$D13.pval))



## Filter based on concordant patterns
# significant in L13
candidate <- subset(info, info$L13.pval < 0.05)
candidate_neut <- subset(info_neut, info_neut$L13.pval < 0.05)
print(paste("First filter, highly differentiated in London (95th percentile):", nrow(candidate)))
# qualitative, opposite direction in individuals who died of the plague
candidate <- subset(candidate, !(sign(candidate$delta_L12) == sign(candidate$delta_L13)))
candidate_neut <- subset(candidate_neut, !(sign(candidate_neut$delta_L12) == sign(candidate_neut$delta_L13)))
print(paste("Second filter, opposite for individuals who died of the plague (no p-value filter):", nrow(candidate)))
# qualitative, opposite direction in individuals who died of the plague
candidate2 <- subset(candidate, sign(candidate$delta_L13) == sign(candidate$delta_D13) & candidate$D13.pval < 0.1)
candidate_neut2 <- subset(candidate_neut, sign(candidate_neut$delta_L13) == sign(candidate_neut$delta_D13) & candidate_neut$D13.pval < 0.1)
print(paste("Third filter, same direction for individuals from Denmark (90th percentile):", nrow(candidate2)))
# shared sign between L13 and D13
candidate <- subset(info, sign(info$delta_D13) == sign(info$delta_L13))
# opposite sign between L12 and L13 (meaning the frequency protective allele goes down then up, rather than down and partly back up)
candidate <- subset(candidate, !(sign(candidate$delta_L12) == sign(candidate$delta_L13)))
# at least one of the directions is highly differentiated relative to the distribution of null locations
# candidate2 <- subset(candidate, candidate$L13.pval < 0.05 & candidate$jointp < 0.05 )
# # shared sign between L13 and D13
# candidate_neut <- subset(info_neut, sign(info_neut$delta_D13) == sign(info_neut$delta_L13))
# # opposite sign between L12 and L23
# candidate_neut <- subset(candidate_neut, !(sign(candidate_neut$delta_L12) == sign(candidate_neut$delta_L13)))
# # at least one of the directions is highly differentiated relative to the distribution of null locations
# candidate_neut2 <- subset(candidate_neut, candidate_neut$L13.pval < 0.05 & candidate_neut$jointp < 0.05 )
# candidate2[order(candidate2$site),]
# nrow(candidate2)
# nrow(candidate_neut2)


### Manhattan plot: all candidate sites

#Filter for only candidate (exon/immune) sites which follow the expected patterns, removing those that don't. 

#This is 535 sites with no filtering based on significance. 
library(ggnewscale)
library(qqman); library(data.table); library(lattice)
for (i in 1:22) {candidate$chr[candidate$chr == paste("chr",i,sep="")] <- i}; rm(i)
candidate$chr <- as.numeric(candidate$chr)
library(dplyr)
don <- candidate %>% 
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(pos)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(candidate, ., by=c("chr"="chr")) %>%
  # Add a cumulative position of each SNP
  arrange(chr, pos) %>%
  mutate( BPcum=pos+tot)
axisdf = don %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
don$is.highlight <- 0; don$is.highlight[don$site %in% candidate2$site] <- 1
don$is.diff <- 0; don$is.diff[don$L13.pval < 0.05] <- 1
don$end <- c(1,don$BPcum[-length(don$BPcum)])
don$col = ifelse(as.numeric(don$chr) %% 2 == 0, 0, 1)
panel_manhattan <- ggplot(don, aes(x=BPcum, y=london_post_pre.fst)) + 
  # facet_grid(. ~ chr) + 
  # chromosome rectangles
  geom_rect(aes(xmin=BPcum, 
                xmax =end, 
                ymin = -Inf, 
                ymax = Inf, 
                fill = factor(col)))  + 
  scale_fill_manual(values=c('grey95', 'white')) +
  
  # points
  geom_point(data = don[don$chr %in% seq(1,22,2),], aes(color = 1.3*10^(.5-don[don$chr %in% seq(1,22,2),]$L13.pval)), size=1.3*10^(.5-don[don$chr %in% seq(1,22,2),]$D13.pval)) + 
  scale_color_gradient(low = "grey80", high = "#366B7D") + 
  
  new_scale_color() + 
  geom_point(data = don[don$chr %in% seq(2,22,2),], aes(color = 1.3*10^(.5-don[don$chr %in% seq(2,22,2),]$L13.pval)), size=1.3*10^(.5-don[don$chr %in% seq(2,22,2),]$D13.pval)) + 
  scale_color_gradient(low = "lightgray", high = "#5C438A") + 
  
  # custom axes:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  coord_cartesian(ylim=c(0,0.09), xlim=c(2487000, 2648153000)) +
  
  # Add highlighted points
  geom_point(data=subset(don, don$is.highlight==1), aes(x=BPcum, y=london_post_pre.fst), color="orange", size=2.2*10^(.5-subset(don, don$is.highlight==1)$D13.pval)) +
  
  # Add label for highlighted points
  geom_text(aes(label=subset(don, don$is.highlight==1)$site), data=subset(don, don$is.highlight==1), nudge_x=-150000000) +
  
  # legends and labels
  theme_classic() + xlab(label = "chromosome") + ylab(label="Fst: London before vs after plague") +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(), 
    axis.text = element_text(color='black'), text=element_text(size=16)
  )
panel_manhattan  
  

rm(design, design2, bins, axisdf, min_neutral_nsites, min_n, min_maf, method, mafs, window_size)
rm(list=ls(pattern="missing_"))
rm(list=ls(pattern="samples_"))
#candidate2$rsid <- c(
#  "rs1052025", "rs17473484", "rs11571319", "rs2549794"
#)
candidate2
min(info$maf)
# library(Rgb);library(data.table)
# gtf <- read.gtf("~/my_genomes/hg37/Homo_sapiens.GRCh37.87.chr.gtf.gz")
# 
# head(gtf)
# gtf <- gtf[gtf$feature == "gene" & gtf$gene_biotype == "protein_coding",]
# 
# info$gene <- NA
# for (i in 1:nrow(info)) {
#   tmp <- subset(gtf, gtf$start - 100000 <= info$pos[i] & gtf$end + 100000 >= info$pos[i])
#   if (nrow(tmp) >=1 ) {
#     tmp$dist <- abs(tmp$start - info$pos[i])
#     info$gene[i] <- tmp$gene_name[tmp$dist == min(tmp$dist)]
#   }
#   print(i)
# }
# write.table(unique(info$gene), "~/set_of_all_genes_maf10.txt", quote=F, row.names=F, col.names=F)
# 
# write.table(unique(info[info$L13.pval < 0.05 & !(sign(info$delta_L12) == sign(info$delta_L13)),]$gene), "~/genes_london_concordant_maf10.txt", quote=F, row.names=F, col.names=F)
# rm(gtf, tmp, i)
#rm(candidate, candidate_neut)
D_replicated_in_Denmark <- candidate2
A_all_loci <- info
B_L13_95th <- subset(info, info$L13.pval < 0.05)
C_replicated_in_plague_burials <- subset(B_L13_95th, !(sign(B_L13_95th$delta_L13) == sign(B_L13_95th$delta_L12)))
#rm(info, candidate2)

save('~/ML_afs.RData')