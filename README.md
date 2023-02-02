# Matters arising

This repository contains scripts for rerunning the analysis in: 

**downsample.sh** contains a number of scripts for generating the down-sampled coverage profiles created from UK10k samples and then running those through the rest of the pipeline

**Figure1_plots.md** contians the scripts to generate the density and q-q plots in figure 1. It reads in the p-values generated from the downsample.sh, coverage_100.r permute_samples.sh, and run_permutations.r. 

**Figure2_plots.md** contains the scripts to generate the plost in figure 2. The simulations need for 2a and 2b are included here. The data for figures 2c and 2d are read in from Klunk_afs.r and ML_afs.r.

**Klunk_afs.r** runs the original Klunk et al. pipeline and saves the data for use in Figure2_plots.md.

**ML_afs.r** calculates the allele frequencies and runs the rest of the Klunk et al. pipeline using a maximum likelihood estimator for use in Figure2_plots.md.

