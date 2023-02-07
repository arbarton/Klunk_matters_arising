---
title: "Matters Arising simulations"
author: "Iain Mathieson"
date: "9/1/2022"
output:
  html_document: default
  pdf_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
set.seed(12345)
```

```{r}
#Pontus FST code
FST_WC=function(NpopA, NpopB, popAcount,popBcount){

	  npops= 2.0
	  nsamples = (NpopA + NpopB)
	  n_bar= (NpopA / npops) + (NpopB / npops)
	  samplefreq = ( (popAcount+popBcount) / (NpopA + NpopB) )
	  pop1freq = popAcount / (NpopA )
	  pop2freq = popBcount / (NpopB )
	  Npop1 = NpopA
	  Npop2 = NpopB
	  S2A= (1/ ( (npops-1.0) * n_bar) ) * ( ( (Npop1)* ((pop1freq-samplefreq)**2) ) + ( (Npop2)*((pop2freq-samplefreq)**2) ) )
	  nc = 1.0/(npops-1.0) * ( (Npop1+Npop2) - (((Npop1**2)+(Npop2**2)) / (Npop1+Npop2)) )
	  T_1 = S2A -( ( 1/(n_bar-1) ) * ( (samplefreq * (1-samplefreq)) -  ((npops-1)/npops)* S2A ) )
	  T_2 = (( (nc-1) / (n_bar-1) ) * samplefreq *(1-samplefreq) )   +  (1.0 +   (((npops-1)*(n_bar-nc))  / (n_bar-1)))       * (S2A/npops)

	  return (T_1 / T_2)
}

```

Start by running simulations to see what proportion of loci pass the FST and direction filters. Conservatively assume that the allele frequencies in every population are identical. 

```{r}

Beforesamplesize = 38*2
Duringsamplesize = 42*2
Aftersamplesize = 63*2
DaneBeforesamplesize = 29*2
DaneAftersamplesize = 34*2
DaneDuringsamplesize = 11*2 # Only for input format purposes; ignored by selection coefficient inference script.

freqs <- seq(0.1, 0.5, 0.1)
probs <- 0*freqs
all_df<-data.frame()

for(i in 1:length(freqs)){

frequency <- freqs[i]
N.sim <- 1000000

#Sample sizes in London and Denmark
Before_f_est = rbinom(N.sim,Beforesamplesize,frequency)/Beforesamplesize
During_f_est = rbinom(N.sim,Duringsamplesize,frequency)/Duringsamplesize    
After_f_est  = rbinom(N.sim,Aftersamplesize,frequency)/Aftersamplesize
Before_DK_est = rbinom(N.sim,DaneBeforesamplesize,frequency)/DaneBeforesamplesize
After_DK_est = rbinom(N.sim,DaneAftersamplesize,frequency)/DaneAftersamplesize
During_DK_est = rbinom(N.sim,11,frequency)/11 # Only for input format purposes; ignored by selection coefficient inference script.

#Direction conditions
cond1 <- (During_f_est < Before_f_est) & (After_f_est > Before_f_est) & (After_DK_est > Before_DK_est)
cond2 <- (During_f_est > Before_f_est) & (After_f_est < Before_f_est) & (After_DK_est < Before_DK_est)

#FST conditions
fst1 <- FST_WC(Beforesamplesize,Aftersamplesize, Beforesamplesize*Before_f_est, Aftersamplesize*After_f_est)
fst2 <- FST_WC(DaneBeforesamplesize,DaneAftersamplesize, DaneBeforesamplesize*Before_DK_est, DaneAftersamplesize*After_DK_est)
cond3 <- fst1>quantile(fst1, 0.95, na.rm=T)
cond4 <- fst2>quantile(fst2, 0.90, na.rm=T)

prob2 <- mean(cond3 & cond4 & (cond1 | cond2))
probs[i] <- prob2


#df = df[cond3 & cond4 & (cond1 | cond2),]
}

mean_prob=mean(probs)
plot(freqs, probs, type="b", ylim=c(0,0.0002), xlab="Allele frequency", ylab="Proportion passing filter")
```

A proportion P=`r mean_prob` of simulations pass these filters under the null. Does not depend on frequency. Next, sample the distribution of Fst values assuming no change in frequency. Produce historgrams of Fst and estimated selection coefficient for f=0.438. 

```{r fig.height = 5, fig.width = 5}
frequency=0.438 
N.sim=10000000

Before_f_est = rbinom(N.sim,Beforesamplesize,frequency)/Beforesamplesize
During_f_est = rbinom(N.sim,Duringsamplesize,frequency)/Duringsamplesize    
After_f_est  = rbinom(N.sim,Aftersamplesize,frequency)/Aftersamplesize
Before_DK_est = rbinom(N.sim,DaneBeforesamplesize,frequency)/DaneBeforesamplesize
After_DK_est = rbinom(N.sim,DaneAftersamplesize,frequency)/DaneAftersamplesize
During_DK_est = rbinom(N.sim,DaneDuringsamplesize,frequency)/DaneDuringsamplesize # Only for input format purposes; ignored by selection coefficient inference script.

fst1 <- FST_WC(Beforesamplesize,Aftersamplesize, Beforesamplesize*Before_f_est, Aftersamplesize*After_f_est)
fst2 <- FST_WC(DaneBeforesamplesize,DaneAftersamplesize, DaneBeforesamplesize*Before_DK_est, DaneAftersamplesize*After_DK_est)
cond1 <- (During_f_est < Before_f_est) & (After_f_est > Before_f_est) & (After_DK_est > Before_DK_est)
cond2 <- (During_f_est > Before_f_est) & (After_f_est < Before_f_est) & (After_DK_est < Before_DK_est)
cond3 <- fst1>quantile(fst1, 0.95, na.rm=T)
cond4 <- fst2>quantile(fst2, 0.90, na.rm=T)

cand <- cond3 & cond4 & (cond1 | cond2)

hist(fst1,main="",xlab="Fst",ylab="Density", probability=TRUE, breaks=100,  xlim=c(-0.02, 0.08), ylim=c(0,100))
abline(v=0.0247,lty=2,col="red",lwd=2)
text(0.05,9,paste("ERAP2, FST=0.0247\n p=",mean(fst1>0.0247),sep=" "),col="red")
```

Produce input files for selection coefficient inference.

```{r}
df = cbind(Before_f_est,During_f_est,After_f_est,Before_DK_est,During_DK_est,After_DK_est)
chr = c(rep(paste0("sim"),nrow(df)))
pos = c(1:nrow(df))

londonAll=cbind(chr,pos,(1-df[,1])*Beforesamplesize,Beforesamplesize,(1-df[,2])*Duringsamplesize,Duringsamplesize,(1-df[,3])*Aftersamplesize,Aftersamplesize)
colnames(londonAll)=c("chr","position","x1","n1","x2","n2","x3","n3")
write.table(file=paste0("london_all_loci_derived_counts_fA",frequency,".txt"),x = londonAll,quote = FALSE,sep = "\t",row.names = FALSE)
write.table(file=paste0("london_candidate_loci_derived_counts_fA",frequency,".txt"),x = londonAll[cand,],quote = FALSE,sep = "\t",row.names = FALSE)

daneAll=cbind(chr,pos,(1-df[,4])*DaneBeforesamplesize,DaneBeforesamplesize,(1-df[,5])*DaneDuringsamplesize,DaneDuringsamplesize,(1-df[,6])*DaneAftersamplesize,DaneAftersamplesize)
colnames(daneAll)=c("chr","position","x1","n1","x2","n2","x3","n3")
write.table(file=paste0("denmark_all_loci_derived_counts_fA",frequency,".txt"),x = daneAll, quote = FALSE,sep = "\t",row.names = FALSE)
write.table(file=paste0("denmark_candidate_loci_derived_counts_fA",frequency,".txt"),x = daneAll[cand,], quote = FALSE,sep = "\t",row.names = FALSE)

```
At this point run the selection inference command and plot the distribution of estimated selection coefficients conditional on passing the filter. 
```{r fig.height = 5, fig.width = 5}
#system("/anaconda/bin/python ./infer_selection_coefficients/Compute_LLs_est_shat.py _candidate_loci_derived_counts_fA0.438.txt Candidate --Ne 5e3,5e3")
shat<-read.table("Candidate_all_MaxLLRs.txt", header=TRUE)
hist(abs(shat$s2hat_all), xlab="Maximum likelihood estimate of selection coefficient",ylab="Density", probability=TRUE, main="")
abline(v=0.39,lty=2,col="red",lwd=2)
text(0.45,3,expression(paste("ERAP2, ", hat(s), "=0.39")),col="red")
```

Finally, estimate frequency of the ERAP2 SNP and plot. 

```{r fig.height = 5, fig.width = 5}
library(ggplot2)
SNP <- "rs2549794"
CHR <- "chr5"
POS <- 96244549
SRC <- "gwas"
N.BOOT <- 10000

likelihood <- function(p, data){
    gt.freq <- c((1-p)^2, 2*p*(1-p), p^2)
    ll <- sum(log(rowSums(t(t(data)*gt.freq))))
    return(ll)
}

ERAP2<-data.frame(Time=NA, Location=NA, Freq=NA, SD=NA)

#Load data and compute frequencies.
i=1
for(loc in c("London", "Denmark")){
    for(time in c("pre", "during", "post")){
      gl <- read.table(paste0("./Data/genoliks/genolik.", SRC, "_", loc, "_", time, ".genolik"), as.is=TRUE)
      #Which row are we looking at?
      site.id <- which(gl[,2]==POS & gl[,1]==CHR)
      #drop samples missing more than 50% of sites
      gl <- gl[,-c(1,2)]
      gl[gl==-1] <- NA
      keep <- (colSums(is.na(gl))/nrow(gl))<=0.5 #assume missingess is the same for each gt
      gl <- gl[site.id,keep]
      gl <- gl[!is.na(gl)]
      gl <- matrix(gl, ncol=3, byrow=T)
      gl <- gl[!is.na(gl[,1]),]

      opt <- optimize(likelihood, interval=c(0,1), maximum=TRUE, data=gl)
      freq <- opt$maximum
      freq <- 1-freq #ref is derived
      boots <- rep(NA, N.BOOT)
      for(k in 1:N.BOOT){
        boot.data <- gl[sample(NROW(gl), replace=TRUE),]
        boots[k] <-  optimize(likelihood, interval=c(0,1), maximum=TRUE, data=boot.data)$maximum
      }
      ERAP2[i,"Time"]<-time
      ERAP2[i,"Location"]<-loc
      ERAP2[i,"Freq"]<-freq
      ERAP2[i,"SD"]<-sd(boots)
      i=i+1
    }
}

rename.time<-c("Before BD", "BD", "After BD")
names(rename.time)<-c("pre", "during", "post")
ERAP2$Time<-rename.time[ERAP2$Time]
ERAP2$Location[ERAP2$Location=="London"]<-"London/England"

#Extra data (as decribed in paper)
extra<-data.frame(
  Time=c("1-2K BP", "Present", "1-2K BP", "Present"),
  Location=c("Denmark", "Denmark", "London/England", "London/England"),
  Freq=c(0.368, 0.429, 0.466, 0.438),
  SD=c(0.055331708, 0.001561046, 0.041283381, 0.000965098)
  )
ERAP2<-rbind(ERAP2, extra)
ERAP2$Time <- factor(ERAP2$Time, levels=c("1-2K BP", "Before BD", "BD", "After BD", "Present"))
ERAP2$Location <- factor(ERAP2$Location, levels=c("London/England", "Denmark"))

pd <- position_dodge(width = 0.2)
p <- ggplot(ERAP2, aes(x=Time, y=Freq, col=Location))
p <- p+geom_point(size=4, position=pd)
p <- p+theme_classic()+geom_errorbar(aes(ymin=Freq-1.96*SD, ymax=Freq+1.96*SD, width=0), lwd=2.5, position=pd)
p <- p+geom_hline(yintercept=0.438, linetype="dashed", color = "#F8766D")
p <- p+geom_hline(yintercept=0.429, linetype="dashed", color = "#00BFC4")
p <- p+ theme(legend.position = c(0.2, 0.8))
print(p)
```
