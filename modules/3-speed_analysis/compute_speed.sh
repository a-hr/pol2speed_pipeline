#!/usr/bin/env bash

env_path=/scratch/heral/envs/pol2speed
cov_path=/scratch/heral/pol2speed/data

module load Python
conda activate $env_path

# ---- Compute speed ----
# Rscript compute_speed.R $cov_path 

# TODO:
#   - compute speed
#   - save speed table (per sample per intron)
#   - save test results
#   - save plots
