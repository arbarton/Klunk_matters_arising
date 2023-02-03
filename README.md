# Matters arising

This repository contains scripts for running the analyses in *Insufficient evidence for natural selection associated with the Black Death*, Barton et al. (2023).

**Figure1_plots.md** contians the scripts to generate the density and q-q plots in figure 1. It reads in the p-values generated from the downsample.sh, coverage_100.r permute_samples.sh, and run_permutations.r. 

**Figure2_plots.md** contains the scripts to generate the plost in figure 2. The simulations need for 2a and 2b are included here. The data for figures 2c and 2d are read in from Klunk_afs.r and ML_afs.r.

**Klunk_afs.r** runs the original Klunk et al. pipeline and saves the data for use in Figure2_plots.md.

**ML_afs.r** calculates the allele frequencies and runs the rest of the Klunk et al. pipeline using a maximum likelihood estimator for use in Figure2_plots.md.

**downsample.sh** contains a number of scripts for generating the down-sampled coverage profiles created from UK10k samples and then running those through the rest of the pipeline.

**coverage_100.r** runs the genotype likelihood files generated by downsample.sh through the rest of the original Klunk et al. pipeline. It retains the same samples as were used in the original analysis.

**permute_samples.sh** shuffles the sample labels for the pre- and post-BD samples while retaining the original number of samples used for each target panel.

**run_permutations.r** runs the original Klunk et al. pipeline but reads in the suffled labels from permute_samples.sh.

