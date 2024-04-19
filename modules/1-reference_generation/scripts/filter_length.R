library(tidyverse)
library(patchwork)
library(biomaRt)

# ---- INPUTS ----
args <- commandArgs(trailingOnly = TRUE)

no.pseudo.bed <- args[1]
filtered.introns.bed <- args[2]

gene.length.threshold <- 1e5
intron.length.threshold <- 5e4

# Read in the data
introns.df <- read.table(
    no.pseudo.bed,
    sep = "\t",
    header = FALSE
) %>%
    dplyr::select(-V5) %>%
    rename(
        "chr" = "V1",
        "int.start" = "V2",
        "int.end" = "V3",
        "id.rank" = "V4",
        "strand" = "V6"
    ) %>%
    separate_wider_delim("id.rank", ": ", names = c("ensembl_transcript_id", "rank"))

# Get annotations
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

annotations <- getBM(
    filters = "ensembl_transcript_id",
    values = unique(introns.df$ensembl_transcript_id),
    attributes = c(
        "ensembl_transcript_id",
        "ensembl_gene_id",
        "start_position",
        "end_position"
    ),
    mart = ensembl
)

introns.ann <- introns.df %>%
    left_join(y = annotations, by = "ensembl_transcript_id") %>%
    mutate(
        gene_length = end_position - start_position,
        intron_length = int.end - int.start
    ) %>%
    dplyr::select(-any_of(c(
        "length",
        "start_position",
        "end_position"
    )))

# Filter genes shorter than 100kb
introns.filt <- introns.ann %>%
    dplyr::filter(
        gene_length >= gene.length.threshold,
        intron_length >= intron.length.threshold
    ) # 36k introns

introns.filt %>%
    group_by(int.start, int.end) %>%
    summarise(group = paste(ensembl_transcript_id, collapse = ",")) %>%
    ungroup() %>%
    nrow() # 15,414 collapsed introns

# Output BED
introns.filt %>%
    unite(id.rank, c("ensembl_transcript_id", "rank")) %>%
    dplyr::select(c(chr, int.start, int.end, id.rank, intron_length, strand)) %>%
    write.table(
        filtered.introns.bed,
        row.names = F,
        quote = F,
        sep = "\t"
    )

