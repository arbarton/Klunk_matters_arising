                                                                                                                        1,1           All#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-12:00:00
#SBATCH -J make_bams
#SBATCH --array=1-100

#X=$1
Y=${SLURM_ARRAY_TASK_ID}

for X in {1..206}; do
awk -v num=${X} '{print $1,$2,$(num+2)}' trial_${Y}_fraction.tsv | awk -v X=${X} -v Y=${Y} '{if ($3 >=1) print>"./run_"Y"/sample_"X"/sample_"X"_frac1.txt"; else print>"./run_"Y"/sample_"X"/sample_"X"_frac"$3".txt"}' 
awk -v num=${X} '{print $1,$2,$(num+2)}' trial_${Y}_fraction.tsv | awk -v X=${X} -v Y=${Y} '{if ($3 >=1) print $1,$2-1,$2 >"./run_"Y"/sample_"X"/sample"X".frac1.bam.bed"; else print $1,$2-1,$2 >"./run_"Y"/sample_"X"/sample"X".frac"$3".bam.bed"}' 
done

for X in {1..206}; do

awk -v num=${X} '{print $1,$2,$(num+2)}' trial_${Y}_fraction.tsv | awk -v X=${X} -v Y=${Y} '{if ($3 >=1) print $1,$2-1,$2 >"./run_"Y"/sample_"X"/sample"X".frac1.bam.bed"; else print $1,$2-1,$2 >"./run_"Y"/sample_"X"/sample"X".frac"$3".bam.bed"}' 
done

for X in {1..206}
do
bamname=$(awk -v num=${X} 'NR==num {print $0}' bam_list.tmp.trial_${Y}.txt)
echo ${bamname}
for filename in ./run_${Y}/sample_${X}/sample_${X}_frac*.txt
        do while read -r a b c
            do fraction=$c ; echo $a $b $b >> ${filename}.sites
            done < $filename
        samtools view ${bamname} -h -q 20 -L ${filename}.sites -s ${fraction} -o ./run_${Y}/sample_${X}/sample${X}.frac${fraction}.bam
done
sed -i '/V1/d' ./run_${Y}/sample_${X}/sample_${X}_frac1.txt.sites
samtools view ${bamname} -h -q 20 -L ./run_${Y}/sample_${X}/sample_${X}_frac1.txt.sites -o ./run_${Y}/sample_${X}/sample${X}.frac1.bam
rm ./run_${Y}/sample_${X}/sample_${X}_frac*.txt
done
