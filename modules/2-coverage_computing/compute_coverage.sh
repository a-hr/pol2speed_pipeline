#!/usr/bin/env zsh

eval "$(conda shell.bash hook)"
conda activate pol2speed

bams_dir=/Users/varo/Desktop/frailty/bams
bed_path=/Users/varo/Desktop/frailty/intron_ref.bed
output_path=/Users/varo/Desktop/frailty/out

# ---- Compute coverage ----
which python

python coverage.py \
    $bams_dir \
    $bed_path \
    $output_path \
    2> $output_path/compute_coverage.log