library(biomaRt)
library(tidyverse)

options(scipen = 20)

# ---- INPUTS ----
args <- commandArgs(trailingOnly = TRUE)

out.dir <- args[1]

# ---- QUERIES ----
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

transcript.table <- getBM(
    attributes = c(
        "ensembl_transcript_id",
        "ensembl_gene_id",
        "transcription_start_site",
        "transcript_length",
        "start_position",
        "end_position",
        "chromosome_name",
        "strand",
        "external_gene_name",
        "gene_biotype"
    ),
    mart = ensembl
)

exon.table <- getBM(
    attributes = c(
        "ensembl_transcript_id",
        "ensembl_gene_id",
        "exon_chrom_start",
        "exon_chrom_end",
        "chromosome_name",
        "strand",
        "rank"
    ),
    mart = ensembl
)

# ---- SAVE FILES ----

write_delim(
    transcript.table,
    file = paste0(out.dir, "/", "transcript_reference.tab"),
    delim = "\t"
)

write_delim(
    exon.table,
    file = paste0(out.dir, "/", "exon_reference.tab"),
    delim = "\t"
)
