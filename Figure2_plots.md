# Matters arising: Generate Fig 2 plots


```R
knitr::opts_chunk$set(echo = TRUE)
```

## R functions for frequency estimation


```R
# ML estimator 
likelihood <- function(p, data){
  gt.freq <- c((1-p)^2, 2*p*(1-p), p^2)
  ll <- sum(log(rowSums(t(t(data)*gt.freq))))
  return(ll)
}

estimate_af_ml<-function(d){
  gl <- t(d)  
  gl <- gl[!is.na(gl[,1]),,drop=FALSE]
  opt <- optimize(likelihood, interval=c(0,1), maximum=TRUE, data=gl)
  af <- opt$maximum
  return(af)
}

# Estimator from Klunk et al.
estimate_af_klunk <- function(like){
  sum(rowMeans(like,na.rm=T)*0:2/2)
}

# Function that applies both estimators given GLs
get_all_af_estimates <- function(like){
  
  af_ml <- estimate_af_ml(like)
  af_klunk <- estimate_af_klunk(like)
  
  c(af_ml,af_klunk)
}
```

## Estimation of frequencies from simulated data


```R
# Function for simulating GLs given a mean depth 
sim_gls<-function(x,d=5,e=0.01,norm=FALSE){
  n<-length(x)
  dep<-rpois(n,d)
  nA<-rbinom(n,dep,c(e,0.5,1-e)[x+1])
  res<-rbind(dbinom(nA,dep,e),dbinom(nA,dep,0.5),dbinom(nA,dep,1-e))
  if(norm)
    res<-t(t(res)/colSums(res))
  res
}

# Simulate frequencies for M sites
set.seed(10)
true_af <- runif(M <- 30000)

# Specify number of individuals 
N <- 300

# Specify mean sequencing depth for these (for now 5 on average)

# - The mean depth varies between individuals but on average it is approx. 5
d <- runif(N)*10

# - All have mean depth of 5
d2 <- rep(5,N)

# For each of the M sites we simulate data and estimate frequencies
est_al_all<-matrix(NA,nrow=M,ncol=2)
est_al_all2<-matrix(NA,nrow=M,ncol=2)

print("Simulating GLs for the following mean depths averaged across individuals:")
print(c(mean(d),mean(d2)))

for(i in 1:M){
  # Simulate genotypes based on frequencies   
  g <- rbinom(N,2,true_af[i])
  # Simulate GLs from genotypes and depths   
  like <- sim_gls(g,d,norm=T)
  like2 <- sim_gls(g,d2,norm=T)
  # Estimate frequencies
  est_al_all[i,] <- get_all_af_estimates(like)
  est_al_all2[i,] <- get_all_af_estimates(like2)
}

```

    [1] "Simulating GLs for the following mean depths averaged across individuals:"
    [1] 4.785905 5.000000


## Generate plots for simulations


```R
library('ggplot2')
```


```R
sim_data <- data.frame(cbind(true_af,est_al_all))
colnames(sim_data) <- c('true_af','like','like2')

```


```R
true_ML_plot <- ggplot(data=sim_data,aes(x=true_af,y=like))+ 
    geom_point(pch=20,col=cols[1]) +
    xlab('True frequency') + 
    ylab('ML frequency estimate') + 
    xlim(0,1) + ylim(0,1) +
    theme_bw(base_size = 24) +
    theme(
        axis.ticks = element_line(size = 0.5),
        panel.grid = element_blank()
    )+ geom_abline(lwd=1)
true_ML_plot

true_Klunk_plot <- ggplot(data=sim_data,aes(x=true_af,y=like2))+ 
    geom_point(pch=20,col=cols[1]) +
    xlab('True frequency') + 
    ylab('Klunk et al. frequency estimate') + 
    xlim(0,1) + ylim(0,1) +
    theme_bw(base_size = 24) +
    theme(
        axis.ticks = element_line(size = 0.5),
        panel.grid = element_blank()
    )+ geom_abline(lwd=1)
true_Klunk_plot 
```


    
![png](output_9_0.png)
    



    
![png](output_9_1.png)
    


## Use real data to compare frequency estimators


```R
# Read in ML allele frequencies produced by ML_afs.r
load('~/ML_afs.RData')

# Read in allele frequencies from original pipeline before dropping MAF < 5%
load('~/Original_afs.RData')
```


```R
# Create unique site identifiers for each variant 
ML_info$site_type <- paste(ML_info$site,ML_info$type,sep="_")
ML_info$ML_maf <- ML_info$alternate
ML_neut$site_type <- paste(ML_neut$site,ML_neut$type,sep="_")
ML_neut$ML_maf <- ML_neut$alternate

info$site_type <- paste(info$site,info$type,sep="_")
info_neut$site_type <- paste(info_neut$site,info_neut$type,sep="_")
```


```R
# Merge frames 
info_all <- merge(info,ML_info[,c('site_type','ML_maf')],by='site_type')
neut_all <- merge(info_neut,ML_neut[,c('site_type','ML_maf')],by='site_type')

# bind candidate and neutral sites
all_sites <- rbind(info_all,neut_all)
```


```R
ML_Klunk_plot <- ggplot(data=all_sites,aes(x=ML_maf,y=alternate))+ 
    geom_point(pch=20,col=cols[1]) +
    xlab('ML frequency estimate') + 
    ylab('Klunk et al. frequency estimate') + 
    xlim(0,1) + ylim(0,1) +
    theme_bw(base_size = 24) +
    theme(
        axis.ticks = element_line(size = 0.5),
        panel.grid = element_blank()
    )+ geom_abline(lwd=1)
ML_Klunk_plot
```


    
![png](output_14_0.png)
    


## Create final figure


```R
load('~/ML_afs.RData')
```


```R
# Use Manhattan plot generated from ML_afs.r
panel_manhattan  
```


    
![png](output_17_0.png)
    



```R
library(patchwork)
```

    Warning message:
    “package ‘patchwork’ was built under R version 4.2.1”



```R
Figure_2 <- (true_ML_plot + true_Klunk_plot + ML_Klunk_plot)/ panel_manhattan + plot_annotation(tag_levels = 'A')
Figure_2
```


    
![png](output_19_0.png)
    



```R

```
