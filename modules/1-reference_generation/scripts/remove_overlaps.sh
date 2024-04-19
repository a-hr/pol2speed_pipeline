#!/usr/bin/env bash

bed_file=$1
human_exons=$2
outfile=$3

cat $bed_file \
    | bedtools sort \
    | bedtools subtract -s -a stdin -b $human_exons \
    | bedtools sort \
    | bedtools merge -i stdin -s -c 4,5,6 -o distinct \
    > $outfile
