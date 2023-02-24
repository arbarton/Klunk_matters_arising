
#!/bin/bash

#SBATCH -p short
#SBATCH -t 0-02:00:00
#SBATCH --mem=2G


TRIAL=${SLURM_ARRAY_TASK_ID}

##Randomly generate list of 1kg samples to use

KG_list="1000g_samples.tsv"

awk '{print $1}' ${KG_list} | shuf -n 206 > sample_list.tmp.trial_${TRIAL}.txt 

awk 'ARGIND==1 {UKids[$1]=1} ARGIND==2 {if (UKids[$1]==1) print $2}' sample_list.tmp.trial_${TRIAL}.txt  ${KG_list} > bam_list.tmp.trial_${TRIAL}.txt 

##calculate read depth for each of the runs
./samtools-1.16.1/samtools depth -Q 20 -q 20 -H -b hglft_Klunk_variants_Grch38.bed \
    -f bam_list.tmp.trial_${TRIAL}.txt  > samtools_depth_trial_${TRIAL}.tsv

sort -k1,1 -k2,2n samtools_depth_trial_${TRIAL}.tsv> samtools_depth_trial_${TRIAL}.sort.tsv

mkdir run_${TRIAL}
for X in {1..206}; do mkdir run_${TRIAL}/sample_${X}; done


#Need to run Rscript to create fraction of read depth to downsample with samtools only needs to be run once but should be run at this point
#Rscript read_fraction.r
