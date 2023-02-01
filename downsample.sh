##Script for generating down-sampled bams matching Klunk et al coverage

#Scro[t to run for each trial to generate list of randomized samples
#!/bin/bash
trial=$1

UK10K_list="uk10k.tsv" #file with sample IDs and paths to bams for files to downsample

grep -v "IID" ${UK10K_list} | awk '{print $1}' | shuf -n 361 > sample_list.tmp.trial_${trial}.txt #generate sample list for each iteration

awk 'ARGIND==1 {UKids[$1]=1} ARGIND==2 {if (UKids[$1]==1) print $2}' sample_list.tmp.trial_${trial}.txt  ${UK10K_list} > bam_list.tmp.trial_${trial}.txt #generate bam list for each iteration

#for each iteration create a matrix of the read depth at each site for each sample in UK10k
~/reich_home/software/samtools-1.16.1/samtools depth -Q 20 -q 20 -H -b /home/arb38/scratch_dir/Klunk_bams/BD/all.variants.nochr.ukbb.bed -f bam_list.tmp.trial_${trial}.txt  > samtools_depth_trial_${trial}.tsv


##Rscript to make the table of the fraction of reads to grab
R
my_range <- (1:100)
Klunk <- read.csv('Klunk_dpeths.tsv',sep="",header=TRUE) #read samtools depth file for all original Klunk samples at all of the tested sites 
for(i in my_range){
    filename <- paste('./samtools_depth_trial_',i,'.tsv',sep="")
    trial <- read.csv(filename,sep="",header=TRUE)
    frac <- cbind(Klunk[1:2],round(Klunk[-1:-2]/trial[-1:-2],2))
    outfile <- paste('./trial_',i,'_fraction.tsv',sep="")
    write.table(frac,outfile,row.names=FALSE,quote=FALSE,sep="\t")
}
quit()


##make directories for output and files with site information and fractions
for Y in {1..100}; do for X in {1..361}; do mkdir ./run_${Y}/sample_${X}; done; done

for Y in {1..100}; do
for X in {1..361}; do
awk -v num=${X} '{print $1,$2,$(num+2)}' trial_${Y}_fraction.tsv | awk -v X=${X} -v Y=${Y} '{if ($3 >=1) print>"./run_"Y"/sample_"X"/sample_"X"_frac1.txt"; else print>"./run_"Y"/sample_"X"/sample_"X"_frac"$3".txt"}' 
awk -v num=${X} '{print $1,$2,$(num+2)}' trial_${Y}_fraction.tsv | awk -v X=${X} -v Y=${Y} '{if ($3 >=1) print $1,$2-1,$2 >"./run_"Y"/sample_"X"/sample"X".frac1.bam.bed"; else print $1,$2-1,$2 >"./run_"Y"/sample_"X"/sample"X".frac"$3".bam.bed"}' 
done
done

##make bams for each fraction set 
##make sure sites with the same or greater coverage preserves all reads
##run this for each trial as Y argument
Y=$1
for X in {1..361}
do

bamname=$(awk -v num=${X} 'NR==num {print $0}' bam_list.tmp.trial_${Y}.txt)
echo ${bamname}
for filename in ./run_${Y}/sample_${X}/sample_${X}_frac*.txt
        do while read -r a b c
            do fraction=$c ; echo $a $b $b >> ${filename}.sites
            done < $filename
        samtools view ${bamname} -h -q 20 -L ${filename}.sites -s ${fraction} -o ./run_${Y}/sample_${X}/sample${X}.frac${fraction}.bam
done
sed -i '/CHR/d' ./run_${Y}/sample_${X}/sample_${X}_frac1.txt.sites
samtools view ${bamname} -h -q 20 -L ./run_${Y}/sample_${X}/sample_${X}_frac1.txt.sites -o ./run_${Y}/sample_${X}/sample${X}.frac1.bam
rm ./run_${Y}/sample_${X}/sample_${X}_frac*.txt


##generate gvcfs for each bam and merge gvcfs per sample
##generating gvcfs separately prevents bam merging from changing the read depths

Y=$1
for X in {1..361}
do
rm ./run_${Y}/merge.list
rm ./run_${Y}/sample_${X}/sample_${X}.L.g.vcf.gz
rm ./run_${Y}/sample_${X}/sample_${X}.L.g.vcf.gz.tbi

sed -i '/X.CHROM/d' ./run_${Y}/sample_${X}/sample${X}.frac1.bam.bed
for bam in ./run_${Y}/sample_${X}/*bam; do
samtools index ${bam}
gatk HaplotypeCaller -ERC GVCF -I ${bam}  -R ../genome_auto.nochr.fa -L ${bam}.bed -O ${bam}.L.g.vcf.gz -mbq 20
done
ls ./run_${Y}/sample_${X}/*.L.g.vcf.gz > ./run_${Y}/sample_${X}/gvcf.list

gatk CombineGVCFs  -R ../genome_auto.nochr.fa  --variant ./run_${Y}/sample_${X}/gvcf.list -O  ./run_${Y}/sample_${X}/sample_${X}.L.g.vcf.gz

echo ./run_${Y}/sample_${X}/sample_${X}.L.g.vcf.gz >> ./run_${Y}/merge.list

#clean up extra files if low on space
rm ./run_${Y}/sample_${X}/sample${X}.frac*.L.g.vcf.gz
rm ./run_${Y}/sample_${X}/sample${X}.frac*.L.g.vcf.gz.tbi
rm ./run_${Y}/sample_${X}/sample${X}.*bam.bed
rm ./run_${Y}/sample_${X}/sample${X}.*bam.bai

done

##Prepare vcf to enter the rest of Klunk et al pipeline
##again run for each trial Y
Y=$1
rm ./run_${Y}/merge.list
for X in {1..361}
do
echo ./run_${Y}/sample_${X}/sample_${X}.L.g.vcf.gz >> ./run_${Y}/merge.list
done
gatk CombineGVCFs -R ../genome_auto.nochr.fa -O ./run_${Y}/run_${Y}_merged.g.vcf.gz -V ./run_${Y}/merge.list
gatk GenotypeGVCFs -V  ./run_${Y}/run_${Y}_merged.g.vcf.gz -R ../genome_auto.nochr.fa -O ./run_${Y}/run_${Y}_merged.geno.vcf.gz
bcftools reheader -s ../trimmed_samples.txt -o ./run_${Y}/run_${Y}_merged.reheader.geno.vcf.gz ./run_${Y}/run_${Y}_merged.geno.vcf.gz
rm ./run_${Y}/run_${Y}_merged.geno.vcf.gz
rm ./run_${Y}/run_${Y}_merged.geno.vcf.gz.tbi

#./run_${Y}/run_${Y}_merged.reheader.geno.vcf.gz can be used as input for next step in Klunk pipeline
#only subsequent change to pipeline should retain all samples retained in original run, rather than bsed on <50% missingness