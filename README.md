# Matters arising

This repository contains scripts for rerunning the analysis in: 

**downsample.sh** contains a number of scripts for generating the down-sampled coverage profiles created from UK10k samples and then running those through the rest of the pipeline

**Figure1_plots.md** contians the scripts to generate the density and q-q plots in figure 1. It reads in the p-values generated from the downsample.sh, permute_samples.sh, and run_permutations.r. 

**Figure2_plots.md** contains the scripts to generate the plost in figure 2. The simulations need for 2a and 2b are included here. The data for figures 2c and 2d are read in from Klunk_afs.r and ML_afs.r.
