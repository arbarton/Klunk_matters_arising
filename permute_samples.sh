#!/bin/bash
##the lists of samples are obtained from supplementary table 1 from Klunk et al. 2022
## for example london_exon_all is a list of IDs for all pre- and post-BD samples that were indicated as being included in the exon analysis in column 10
## or london_gwas_exon_neutral in the matched case is the list of all pre- and post-BD samples typed on all three panels GWAS,exon,neutral, london_exon_neutral is just those with Exon,neutral in column 10


##perform permutations with no fixed total of individuals
#for X in {1..100}
#do mkdir permute_${X}
#shuf london_exon_all > london_exon_all_trial_${X}
#head -35 london_exon_all_trial_${X} > permute_${X}/keep_london_pre_exon
#tail -51 london_exon_all_trial_${X} > permute_${X}/keep_london_post_exon
#rm london_exon_all_trial_${X}
#shuf london_gwas_all > london_gwas_all_trial_${X}
#head -37 london_gwas_all_trial_${X} > permute_${X}/keep_london_pre_gwas
#tail -59 london_gwas_all_trial_${X} > permute_${X}/keep_london_post_gwas
#rm london_gwas_all_trial_${X}
#shuf london_neutral_all > london_neutral_all_trial_${X}
#head -37 london_neutral_all_trial_${X} > permute_${X}/keep_london_pre_neutral
#tail -63 london_neutral_all_trial_${X} > permute_${X}/keep_london_post_neutral
#rm london_neutral_all_trial_${X}
#done

##perform permutations with fixed total of individuals in each target panel
for X in {1..100}
do mkdir permute_${X}
shuf london_gwas_exon_neutral > london_gwas_exon_neutral_trial_${X}
head -34  london_gwas_exon_neutral_trial_${X} >  london_gwas_exon_neutral_trial_${X}_pre 
tail -48  london_gwas_exon_neutral_trial_${X} >  london_gwas_exon_neutral_trial_${X}_post 
shuf london_exon_neutral > london_exon_neutral_trial_${X}
head -1  london_exon_neutral_trial_${X} >  london_exon_neutral_trial_${X}_pre 
tail -3  london_exon_neutral_trial_${X} >  london_exon_neutral_trial_${X}_post 
shuf london_GWAS_neutral > london_GWAS_neutral_trial_${X}
head -2  london_GWAS_neutral_trial_${X} >  london_GWAS_neutral_trial_${X}_pre 
tail -11  london_GWAS_neutral_trial_${X} >  london_GWAS_neutral_trial_${X}_post 
cp london_GWAS  london_GWAS_trial_${X}_pre
cp london_neutral  london_neutral_trial_${X}_pre
cat london_gwas_exon_neutral_trial_${X}_pre  london_GWAS_neutral_trial_${X}_pre london_GWAS_trial_${X}_pre > permute_${X}/keep_london_pre_gwas
cat london_gwas_exon_neutral_trial_${X}_pre  london_exon_neutral_trial_${X}_pre  > permute_${X}/keep_london_pre_exon
cat london_gwas_exon_neutral_trial_${X}_pre london_exon_neutral_trial_${X}_pre london_GWAS_neutral_trial_${X}_pre > permute_${X}/keep_london_pre_neutral
cat london_gwas_exon_neutral_trial_${X}_post london_GWAS_neutral_trial_${X}_post >  permute_${X}/keep_london_post_gwas
cat london_gwas_exon_neutral_trial_${X}_post  london_exon_neutral_trial_${X}_post  > permute_${X}/keep_london_post_exon
cat london_gwas_exon_neutral_trial_${X}_post  london_exon_neutral_trial_${X}_post  london_GWAS_neutral_trial_${X}_post  london_neutral_trial_${X}_pre >  permute_${X}/keep_london_post_neutral 

rm *trial_${X}*
done 
