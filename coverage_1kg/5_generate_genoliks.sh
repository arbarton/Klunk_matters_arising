#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-5:00:00
#SBATCH -n 2
#SBATCH --mem=6G
#SBATCH --array=1-100

Y=${SLURM_ARRAY_TASK_ID}
module load gatk

rm ./run_${Y}/merge.list
for X in {1..206}
do
echo ./run_${Y}/sample_${X}/sample_${X}.L.g.vcf.gz >> ./run_${Y}/merge.list
done
gatk CombineGVCFs -R GRCh38_full_analysis_set_plus_decoy_hla.fa -O ./run_${Y}/run_${Y}_merged.g.vcf.gz -V ./run_${Y}/merge.list
gatk GenotypeGVCFs -V  ./run_${Y}/run_${Y}_merged.g.vcf.gz -R igenomes/GRCh38_full_analysis_set_plus_decoy_hla.fa -O ./run_${Y}/run_${Y}_merged.geno.vcf.gz
bcftools reheader -s ./samples_list.txt -o ./run_${Y}/run_${Y}_merged.reheader.geno.vcf.gz ./run_${Y}/run_${Y}_merged.geno.vcf.gz
rm ./run_${Y}/run_${Y}_merged.geno.vcf.gz
rm ./run_${Y}/run_${Y}_merged.geno.vcf.gz.tbi


module load vcftools
vcftools --gzvcf ./run_${Y}/run_${Y}_merged.reheader.geno.vcf.gz  --mac 1 --max-alleles 2 --remove-indels --minQ 30 --freq --out ./run_${Y}/run_${Y}_merged
vcftools --gzvcf ./run_${Y}/run_${Y}_merged.reheader.geno.vcf.gz  --mac 1 --max-alleles 2 --remove-indels --minQ 30 --singletons --out ./run_${Y}/run_${Y}_merged

vcftools --gzvcf ./run_${Y}/run_${Y}_merged.reheader.geno.vcf.gz  --mac 1 --max-alleles 2 --remove-indels --bed neutral_grch38.bed --minQ 30 --recode --out ./run_${Y}/run_${Y}_neutral_v1
vcftools --vcf ./run_${Y}/run_${Y}_neutral_v1.recode.vcf --exclude-positions ./run_${Y}/run_${Y}_merged.singletons --recode --out ./run_${Y}/run_${Y}_neutral

vcftools --gzvcf ./run_${Y}/run_${Y}_merged.reheader.geno.vcf.gz --mac 1 --max-alleles 2 --remove-indels --bed exon_grch38.bed --minQ 30 --recode --out ./run_${Y}/run_${Y}_exon_v1
vcftools --vcf ./run_${Y}/run_${Y}_exon_v1.recode.vcf --exclude-positions ./run_${Y}/run_${Y}_merged.singletons --recode --out ./run_${Y}/run_${Y}_exon

vcftools --gzvcf ./run_${Y}/run_${Y}_merged.reheader.geno.vcf.gz  --mac 1 --max-alleles 2 --remove-indels --bed immune_grch38.bed --minQ 30 --recode --out ./run_${Y}/run_${Y}_immune_v1
vcftools --vcf ./run_${Y}/run_${Y}_immune_v1.recode.vcf --exclude-positions ./run_${Y}/run_${Y}_merged.singletons --recode --out ./run_${Y}/run_${Y}_immune


bcftools annotate -x ^INFO/DP,^FORMAT/GT,^FORMAT/AD,^FORMAT/DP,^FORMAT/GQ,^FORMAT/PL ./run_${Y}/run_${Y}_immune.recode.vcf | bcftools +setGT -- -ta -nu > ./run_${Y}/run_${Y}_joint.immune.forgenolik.vcf
bcftools annotate -x ^INFO/DP,^FORMAT/GT,^FORMAT/AD,^FORMAT/DP,^FORMAT/GQ,^FORMAT/PL ./run_${Y}/run_${Y}_exon.recode.vcf | bcftools +setGT -- -ta -nu > ./run_${Y}/run_${Y}_joint.exon.forgenolik.vcf
bcftools annotate -x ^INFO/DP,^FORMAT/GT,^FORMAT/AD,^FORMAT/DP,^FORMAT/GQ,^FORMAT/PL ./run_${Y}/run_${Y}_neutral.recode.vcf | bcftools +setGT -- -ta -nu > ./run_${Y}/run_${Y}_joint.neutral.forgenolik.vcf

## Exons
##new indv files are required that just have the individuals that were kept in the original analysis, I put these in a path called Sample_Names_1000 in my original Data folder
vcftools --vcf ./run_${Y}/run_${Y}_joint.exon.forgenolik.vcf --keep ./DATA/genoliks/SampleNames_1000/london_pre_exons.012.indv --recode --out ./run_${Y}/temp; vcftools --vcf ./run_${Y}/temp.recode.vcf --012 --out ./run_${Y}/london_pre_exons
sed '/^#/d' ./run_${Y}/temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ./LCLAE/filtbaboon1b 35 > ./run_${Y}/genolik.exons_london_pre.genolik
vcftools --vcf ./run_${Y}/run_${Y}_joint.exon.forgenolik.vcf --keep ./DATA/genoliks/SampleNames_1000/london_post_exons.012.indv  --recode --out ./run_${Y}/temp; vcftools --vcf ./run_${Y}/temp.recode.vcf --012 --out ./run_${Y}/london_post_exons
sed '/^#/d' ./run_${Y}/temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ./LCLAE/filtbaboon1b 51 > ./run_${Y}/genolik.exons_london_post.genolik
vcftools --vcf ./run_${Y}/run_${Y}_joint.exon.forgenolik.vcf --keep ./DATA/genoliks/SampleNames_1000/london_during_exons.012.indv --recode --out ./run_${Y}/temp; vcftools --vcf ./run_${Y}/temp.recode.vcf --012 --out ./run_${Y}/london_during_exons
sed '/^#/d' ./run_${Y}/temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ./LCLAE/filtbaboon1b 36 > ./run_${Y}/genolik.exons_london_during.genolik

vcftools --vcf ./run_${Y}/run_${Y}_joint.exon.forgenolik.vcf --keep ./DATA/genoliks/SampleNames_1000/denmark_pre_exons.012.indv --recode --out ./run_${Y}/temp; vcftools --vcf ./run_${Y}/temp.recode.vcf --012 --out ./run_${Y}/denmark_pre_exons
sed '/^#/d' ./run_${Y}/temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ./LCLAE/filtbaboon1b 25 > ./run_${Y}/genolik.exons_denmark_pre.genolik
vcftools --vcf ./run_${Y}/run_${Y}_joint.exon.forgenolik.vcf --keep ./DATA/genoliks/SampleNames_1000/denmark_post_exons.012.indv --recode --out ./run_${Y}/temp; vcftools --vcf ./run_${Y}/temp.recode.vcf --012 --out ./run_${Y}/denmark_post_exons
sed '/^#/d' ./run_${Y}/temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ./LCLAE/filtbaboon1b 27 > ./run_${Y}/genolik.exons_denmark_post.genolik

## Immune
vcftools --vcf ./run_${Y}/run_${Y}_joint.immune.forgenolik.vcf --keep ./DATA/genoliks/SampleNames_1000/london_pre_immune.012.indv --recode --out ./run_${Y}/temp; vcftools --vcf ./run_${Y}/temp.recode.vcf --012 --out ./run_${Y}/london_pre_immune; sed '/^#/d' ./run_${Y}/temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ./LCLAE/filtbaboon1b 37 > ./run_${Y}/genolik.gwas_london_pre.genolik
vcftools --vcf ./run_${Y}/run_${Y}_joint.immune.forgenolik.vcf --keep ./DATA/genoliks/SampleNames_1000/london_post_immune.012.indv --recode --out ./run_${Y}/temp; vcftools --vcf ./run_${Y}/temp.recode.vcf --012 --out ./run_${Y}/london_post_immune; sed '/^#/d' ./run_${Y}/temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ./LCLAE/filtbaboon1b 59 > ./run_${Y}/genolik.gwas_london_post.genolik
vcftools --vcf ./run_${Y}/run_${Y}_joint.immune.forgenolik.vcf --keep ./DATA/genoliks/SampleNames_1000/london_during_immune.012.indv --recode --out ./run_${Y}/temp; vcftools --vcf ./run_${Y}/temp.recode.vcf --012 --out ./run_${Y}/london_during_immune; sed '/^#/d' ./run_${Y}/temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ./LCLAE/filtbaboon1b 40 > ./run_${Y}/genolik.gwas_london_during.genolik

vcftools --vcf ./run_${Y}/run_${Y}_joint.immune.forgenolik.vcf --keep ./DATA/genoliks/SampleNames_1000/denmark_pre_immune.012.indv --recode --out ./run_${Y}/temp; vcftools --vcf ./run_${Y}/temp.recode.vcf --012 --out ./run_${Y}/denmark_pre_immune
sed '/^#/d' ./run_${Y}/temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ./LCLAE/filtbaboon1b 24 > ./run_${Y}/genolik.gwas_denmark_pre.genolik

vcftools --vcf ./run_${Y}/run_${Y}_joint.immune.forgenolik.vcf --keep ./DATA/genoliks/SampleNames_1000/denmark_post_immune.012.indv --recode --out ./run_${Y}/temp; vcftools --vcf ./run_${Y}/temp.recode.vcf --012 --out ./run_${Y}/denmark_post_immune
sed '/^#/d' ./run_${Y}/temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ./LCLAE/filtbaboon1b 33 > ./run_${Y}/genolik.gwas_denmark_post.genolik

## neutral
vcftools --vcf ./run_${Y}/run_${Y}_joint.neutral.forgenolik.vcf --keep ./DATA/genoliks/SampleNames_1000/london_pre_neutral.012.indv  --recode --out ./run_${Y}/temp; vcftools --vcf ./run_${Y}/temp.recode.vcf --012 --out ./run_${Y}/london_pre_neutral; sed '/^#/d' ./run_${Y}/temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ./LCLAE/filtbaboon1b 37 > ./run_${Y}/genolik.neutral_london_pre.genolik
vcftools --vcf ./run_${Y}/run_${Y}_joint.neutral.forgenolik.vcf --keep ./DATA/genoliks/SampleNames_1000/london_post_neutral.012.indv  --recode --out ./run_${Y}/temp; vcftools --vcf ./run_${Y}/temp.recode.vcf --012 --out ./run_${Y}/london_post_neutral; sed '/^#/d' ./run_${Y}/temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ./LCLAE/filtbaboon1b 63 > ./run_${Y}/genolik.neutral_london_post.genolik
vcftools --vcf ./run_${Y}/run_${Y}_joint.neutral.forgenolik.vcf --keep ./DATA/genoliks/SampleNames_1000/london_during_neutral.012.indv --recode --out ./run_${Y}/temp; vcftools --vcf ./run_${Y}/temp.recode.vcf --012 --out ./run_${Y}/london_during_neutral; sed '/^#/d' ./run_${Y}/temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ./LCLAE/filtbaboon1b 42 > ./run_${Y}/genolik.neutral_london_during.genolik

vcftools --vcf ./run_${Y}/run_${Y}_joint.neutral.forgenolik.vcf --keep ./DATA/genoliks/SampleNames_1000/denmark_pre_neutral.012.indv --recode --out ./run_${Y}/temp; vcftools --vcf ./run_${Y}/temp.recode.vcf --012 --out ./run_${Y}/denmark_pre_neutral
sed '/^#/d' ./run_${Y}/temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ./LCLAE/filtbaboon1b 27 > ./run_${Y}/genolik.neutral_denmark_pre.genolik

vcftools --vcf ./run_${Y}/run_${Y}_joint.neutral.forgenolik.vcf --keep ./DATA/genoliks/SampleNames_1000/denmark_post_neutral.012.indv  --recode --out ./run_${Y}/temp; vcftools --vcf ./run_${Y}/temp.recode.vcf --012 --out ./run_${Y}/denmark_post_neutral
sed '/^#/d' ./run_${Y}/temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ./LCLAE/filtbaboon1b 31 > ./run_${Y}/genolik.neutral_denmark_post.genolik

## get site information and alternate variants for table
vcftools --vcf ./run_${Y}/run_${Y}_joint.neutral.forgenolik.vcf --freq --out ./run_${Y}/neutral
vcftools --vcf ./run_${Y}/run_${Y}_joint.exon.forgenolik.vcf --freq --out ./run_${Y}/exon
vcftools --vcf ./run_${Y}/run_${Y}_joint.immune.forgenolik.vcf --freq --out ./run_${Y}/immune

mkdir ./run_${Y}/genoliks
mkdir ./run_${Y}/SampleNames

mv ./run_${Y}/*genolik ./run_${Y}/genoliks/
mv ./run_${Y}/*frq ./run_${Y}/genoliks/
mv ./run_${Y}/*indv ./run_${Y}/SampleNames/

rm  ./run_${Y}/*pos
rm  ./run_${Y}/*.012
rm  ./run_${Y}/*.log
rm ./run_${Y}/temp.recode.vc