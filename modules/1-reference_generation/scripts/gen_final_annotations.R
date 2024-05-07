library(tidyverse)

options(scipen = 20)

# ---- INPUTS ----
args <- commandArgs(trailingOnly = TRUE)

collapse.introns.bed <- args[1]
length.filt.bed <- args[2]
transcript.ref <- args[3]
measurable.introns.bed <- args[4]
intron.metadata.tsv <- args[5]

# ---- Load data ----
filtered.introns <- read_delim(
    collapse.introns.bed,
    col_names = c("Chr", "Start", "End", ".TranscriptID", "Length", "Strand")
) %>%
    mutate(IntID = paste("Intron", rownames(.), sep = "."))

all.introns <- read_delim(
    length.filt.bed,
    col_names = c("Chr", "Start", "End", ".TranscriptID", "Length", "Strand")
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
annotations <- read_delim(
    file.path(transcript.ref),
    delim = "\t",
    col_select = c(
        "ensembl_transcript_id",
        "ensembl_gene_id",
        "start_position",
        "end_position",
        "external_gene_name"
    )
) %>%
    dplyr::filter(ensembl_transcript_id %in% unique(meta.introns$TranscriptID))

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
    write_delim(
        paste0(measurable.introns.bed, "2"),
        delim = "\t",
        quote = "none",
        col_names = F
    )

meta.introns %>%
    write_delim(
        paste0(intron.metadata.tsv, "2"),
        delim = "\t",
        quote = "none"
    )
