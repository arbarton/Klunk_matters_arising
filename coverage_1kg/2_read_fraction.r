#Rscript to create fraction of read depth to downsample with samtools

Klunk <- read.csv('samtools_depth_trimmed_incl.GRCH38.sort.match.tsv',sep="",header=FALSE)
my_range <- (1:100)
for(i in my_range){
    filename <- paste('samtools_depth_trial_',i,'.sort.tsv',sep="")
    trial <- read.csv(filename,sep="",header=FALSE)
    frac <- cbind(Klunk[1:2],round(Klunk[-1:-2]/trial[-1:-2],2))
    outfile <- paste('trial_',i,'_fraction.tsv',sep="")
    write.table(frac,outfile,row.names=FALSE,quote=FALSE,sep="\t")
}