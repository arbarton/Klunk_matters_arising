#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-12:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH -c 2
#SBATCH --array=1-206

Y=$1
X=${SLURM_ARRAY_TASK_ID}

module load gatk

#remove any old files from previous runs
rm ./run_${Y}/merge.list
rm ./run_${Y}/sample_${X}/sample_${X}.L.g.vcf.gz
rm ./run_${Y}/sample_${X}/sample_${X}.L.g.vcf.gz.tbi

#index bams and run haplotype caller on per-fraction bam
sed -i '/V1/d' ./run_${Y}/sample_${X}/sample${X}.frac1.bam.bed
for bam in ./run_${Y}/sample_${X}/*bam; do
samtools index ${bam}
gatk HaplotypeCaller -ERC GVCF -I ${bam}  -R GRCh38_full_analysis_set_plus_decoy_hla.fa -L ${bam}.bed -O ${bam}.L.g.vcf.gz -mbq 20
done
ls ./run_${Y}/sample_${X}/*.L.g.vcf.gz > ./run_${Y}/sample_${X}/gvcf.list

#combine by-fraction gvcfs for each sample
gatk CombineGVCFs  -R GRCh38_full_analysis_set_plus_decoy_hla.fa  --variant ./run_${Y}/sample_${X}/gvcf.list -O  ./run_${Y}/sample_${X}/sample_${X}.L.g.vcf.gz

#create list of gvcfs to combine for genotyping
echo ./run_${Y}/sample_${X}/sample_${X}.L.g.vcf.gz >> ./run_${Y}/merge.list

#remove larger intermediate files
rm ./run_${Y}/sample_${X}/sample${X}.frac*.L.g.vcf.gz
rm ./run_${Y}/sample_${X}/sample${X}.frac*.L.g.vcf.gz.tbi
rm ./run_${Y}/sample_${X}/sample${X}.*bam.bed
rm ./run_${Y}/sample_${X}/sample${X}.*bam.bai
