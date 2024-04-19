#!/usr/bin/env bash

output_path=$1
out_name=$2
tmp_path=$3

env_path=$4

module load Python
conda activate $env_path

base_ref_name=intron_base_ref.csv
no_pseudo_name=intron_ref_no_pseudogenes.csv
length_filt_name=intron_ref_length_filt.bed

# ---- Generate reference table ----
echo "$(date +"%T") - Generating base reference table"

mkdir -p $tmp_path $output_path

Rscript scripts/gen_base_annotations.R \
    $tmp_path \
    $output_path \
    $base_ref_name \
    2> $tmp_path/gen_base_annotations.log

# ---- operations to filter pseudogene regions ----
echo "$(date +"%T") - Filtering pseudogenes"

Rscript scripts/filter_pseudogenes.R \
    $output_path/$base_ref_name \
    $tmp_path/$no_pseudo_name \
    2> $tmp_path/filter_pseudogenes.log

# ---- operations to filter overlapping regions ----
echo "$(date +"%T") - Removing duplicated regions, filtering introns and removing shared exonic regions"

# convert csv table to bed
python scripts/csv_to_bed.py \
    $tmp_path/$no_pseudo_name

python scripts/mart_to_bed.py \
    $tmp_path/start_elocs.csv \
    $tmp_path/exons.bed

# filter introns by length
Rscript scripts/filter_length.R \
    $tmp_path/$(basename $no_pseudo_name .csv).bed \
    $tmp_path/$length_filt_name \
    2> $tmp_path/filter_length.log

# remove regions overlapping with exons
bash scripts/remove_overlaps.sh \
    $tmp_path/$length_filt_name \
    $tmp_path/start_elocs.csv \
    $tmp_path/collapsed_introns.bed \
    2> $tmp_path/remove_overlaps.log

# generate final reference table
Rscript scripts/gen_final_annotations.R \
    $tmp_path/collapsed_introns.bed \
    $tmp_path/$length_filt_name \
    $output_path/$out_name \
    $output_path/intron_metadata.tsv \
    2> $tmp_path/gen_final_annotations.log

echo "$(date +"%T") - Reference table generated! Cleaning up temporary files"
echo "Check $output_path for the final reference table"

# remove temporary files
rm $tmp_path/*.csv $tmp_path/*.bed