#!/usr/bin/env bash

#SBATCH --qos=serial
#SBATCH --job-name=generate_reference
#SBATCH --mem=2G
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=1

output_path=/scratch/heral/ad-hoc_pipelines/fragility/1-reference_generation/output
out_name=intron_ref.bed
tmp_path=/scratch/heral/ad-hoc_pipelines/fragility/1-reference_generation/tmp
reference_path=/scratch/heral/ad-hoc_pipelines/fragility/1-reference_generation/assets

env_path=/scratch/heral/envs/pol2speed

module load Python
conda activate $env_path

bash generate_reference.sh $output_path $out_name $tmp_path cache $reference_path "keep"
