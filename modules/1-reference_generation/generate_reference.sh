#!/usr/bin/env bash

output_path=$1
out_name=$2
tmp_path=$3
ensembl_reference=$4  # either "cache" or "download"
ensembl_reference_path=$5  # path to cache or where to download
clean_up=$6  # "clean" or "keep"

gene_length_threshold=1e5
intron_length_threshold=5e4

base_ref_name=base_intron_reference.tab
no_pseudo_name=no_pseudo_intron_reference.tab
length_filt_name=length_filtered_intron_reference.bed

# ---- Get Ensembl reference ----
mkdir -p $output_path $tmp_path

if [ $ensembl_reference == "cache" ]; then
    echo "$(date +"%T") - Using cached Ensembl reference"
    cp $ensembl_reference_path/*.tab $tmp_path
elif [ $ensembl_reference == "download" ]; then
    echo "$(date +"%T") - Downloading Ensembl reference"
    Rscript scripts/download_reference.R \
        $ensembl_reference_path \
        2> $tmp_path/01_download_reference.log
    cp $ensembl_reference_path/*.tab $tmp_path
else
    echo "Invalid argument for ensembl_reference. Use 'cache' or 'download'"
    exit 1
fi

# ---- Generate reference table ----
echo "$(date +"%T") - Generating base reference table"

# check if base reference table already exists
if [ -f $output_path/$base_ref_name ]; then
    echo "Base reference table already exists. Skipping..."
else
    Rscript scripts/gen_base_annotations.R \
        $tmp_path \
        $output_path \
        $base_ref_name \
        $tmp_path \
        2> $tmp_path/02_gen_base_annotations.log
fi

# ---- operations to filter pseudogene regions ----
echo "$(date +"%T") - Filtering pseudogenes"

Rscript scripts/filter_pseudogenes.R \
    $output_path/$base_ref_name \
    $tmp_path/$no_pseudo_name \
    $tmp_path \
    2> $tmp_path/03_filter_pseudogenes.log

# ---- operations to filter overlapping regions ----
echo "$(date +"%T") - Removing duplicated regions, filtering introns and removing shared exonic regions"

# convert csv tables to bed
python scripts/tab_to_bed.py \
    $tmp_path/$no_pseudo_name \
    $tmp_path/$(basename $no_pseudo_name .tab).bed

python scripts/mart_to_bed.py \
    $tmp_path/start_elocs.csv \
    $tmp_path/start_elocs.tmp

cat $tmp_path/start_elocs.tmp \
    | sortBed > $tmp_path/start_elocs.bed
rm $tmp_path/start_elocs.tmp

# filter introns by length
Rscript scripts/filter_length.R \
    $tmp_path/$(basename $no_pseudo_name .tab).bed \
    $tmp_path/$length_filt_name \
    $tmp_path \
    $gene_length_threshold \
    $intron_length_threshold \
    2> $tmp_path/04_filter_length.log

# remove regions overlapping with exons
bash scripts/remove_overlaps.sh \
    $tmp_path/$length_filt_name \
    $tmp_path/start_elocs.bed \
    $tmp_path/collapsed_introns.bed \
    2> $tmp_path/05_remove_overlaps.log

# generate final reference table
Rscript scripts/gen_final_annotations.R \
    $tmp_path/collapsed_introns.bed \
    $tmp_path/$length_filt_name \
    $tmp_path/transcript_reference.tab \
    $output_path/$out_name \
    $output_path/intron_metadata.tsv \
    2> $tmp_path/06_gen_final_annotations.log

echo "$(date +"%T") - Reference table generated!"

# remove temporary files
if [ $clean_up == "clean" ]; then
    echo "Cleaning up temporary files..."
    rm -rf $tmp_path
else
    echo "Keeping temporary files..."
fi

echo "Check $output_path for the final reference table."
echo "Goodbye!"
