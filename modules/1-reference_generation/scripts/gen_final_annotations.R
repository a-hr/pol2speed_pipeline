library(dplyr)
library(tidyverse)
library(biomaRt)


# ---- INPUTS ----
args <- commandArgs(trailingOnly = TRUE)

collapse.introns.bed <- args[1]
allintrons.100k.bed <- args[2]

measurable.introns.bed <- args[3]
intron.metadata.tsv <- args[4]

# ---- Load data ----
filtered.introns <- read.delim(
    collapse.introns.bed,
    header = F,
    col.names = c("Chr", "Start", "End", ".TranscriptID", "Length", "Strand")
) %>%
    mutate(IntID = paste("Intron", rownames(.), sep = "."))

all.introns <- read.delim(
    allintrons.100k.bed,
    header = F,
    col.names = c("Chr", "Start", "End", ".TranscriptID", "Length", "Strand")
)

# ---- Demultiplex ENSTIDs and assign intron ID ----
meta.introns <- filtered.introns %>%
    separate_longer_delim(cols = .TranscriptID, delim = ",") %>%
    separate_wider_delim(cols = .TranscriptID, delim = "_", names = c("TranscriptID", "Rank")) %>%
    dplyr::select(-c(Start, End, Length, Chr, Strand)) %>%
    mutate(.TranscriptID = paste(TranscriptID, Rank, sep = "_")) %>%
    left_join(y = all.introns, by = ".TranscriptID") %>%
    dplyr::select(-".TranscriptID")

# ---- Enrich with metadata ----
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

annotations <- getBM(
    filters = "ensembl_transcript_id",
    values = unique(meta.introns$TranscriptID),
    attributes = c(
        "ensembl_transcript_id",
        "ensembl_gene_id",
        "start_position",
        "end_position",
        "external_gene_name"
    ),
    mart = ensembl
)

annotations <- annotations %>%
    mutate(GeneLength = end_position - start_position) %>%
    dplyr::rename(
        "GeneID" = "ensembl_gene_id",
        "TranscriptID" = "ensembl_transcript_id",
        "GeneName" = "external_gene_name"
        ) %>% 
    dplyr::select(-c(start_position, end_position))

meta.introns <- meta.introns %>%
    dplyr::rename("IntronLength" = "Length") %>%
    left_join(y = annotations, by = "TranscriptID")

# ---- Output tables ----
filtered.introns %>%
    dplyr::select(c(Chr, Start, End, IntID, Length, Strand)) %>%
    rowwise() %>%
    mutate(Length = End - Start) %>%
    write.table(
        measurable.introns.bed,
        sep = "\t",
        quote = F,
        col.names = F,
        row.names = F
    )

meta.introns %>%
    write.table(
        intron.metadata.tsv,
        sep = "\t",
        quote = F,
        row.names = F
    )
