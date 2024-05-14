#!/usr/bin/env bash

env_path=/scratch/heral/envs/pol2speed

output_path=/scratch/heral/pol2speed/outputs
tmp_path=/scratch/heral/pol2speed/temp

# ---- Generate reference table ----
mkdir -p $output_path/reference $tmp_path

bash modules/1-reference_generation/generate_reference.sh $output_path/reference human_intron_reference.csv $tmp_path $env_path

# ---- Compute coverage ----
bash modules/2-coverage_computing/compute_coverage.sh $bams_dir 

# ---- Compute speed ----
bash modules/3-speed_analysis/compute_speed.sh

echo "$(date +"%T") - Pipeline finished successfully"