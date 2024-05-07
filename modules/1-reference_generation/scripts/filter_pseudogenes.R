library(tidyverse)

options(scipen = 20)

# ---- INPUTS ----
args <- commandArgs(trailingOnly = TRUE)

intron.base.table <- args[1]
out.name <- args[2]
ref.dir <- args[3]

# ----- QUERY ANNOTATIONS (biomaRt) -----
intron.df <- read_delim(intron.base.table)

transcript.list <- unique(intron.df$ensembl_transcript_id)

# launch the request
annotations <- read_delim(
    file = paste0(ref.dir, "/", "transcript_reference.tab"),
    delim = "\t"
    ) %>%
    dplyr::select(all_of(
        c(
            "ensembl_gene_id",
            "ensembl_transcript_id",
            "external_gene_name",
            "gene_biotype"
        )
    )) %>%
    dplyr::filter(ensembl_transcript_id %in% transcript.list)

# filter out the transcripts whose biotype is "pseudogene"
# unique(annotations$gene_biotype)
blacklist <-
    c(
        "processed_pseudogene",
        "unprocessed_pseudogene",
        "transcribed_processed_pseudogene",
        "transcribed_unprocessed_pseudogene",
        "IG_V_pseudogene",
        "IG_C_pseudogene",
        "transcribed_unitary_pseudogene",
        "TR_V_pseudogene",
        "unitary_pseudogene",
        "artifact"
    )

transcripts.remove <- annotations %>%
    filter(gene_biotype %in% blacklist) %>%
    pull(ensembl_transcript_id)

intron.df2 <- intron.df %>%
    filter(!ensembl_transcript_id %in% transcripts.remove) %>%
    write_delim(
        file = out.name,
        delim = "\t"
    )
