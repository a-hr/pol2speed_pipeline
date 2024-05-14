#!/usr/bin/env bash

env_path=/scratch/heral/envs/pol2speed
bam_path=/scratch/heral/pol2speed/data

module load Python
conda activate $env_path

bams_dir=
bed_path=
output_path=

# ---- Compute coverage ----
python coverage.py \
    $bams_dir \
    $bed_path \
    $output_path \
    2> $output_path/compute_coverage.log