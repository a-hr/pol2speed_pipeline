library(tidyverse)

options(scipen = 20)

# ---- INPUTS ----
args <- commandArgs(trailingOnly = TRUE)

no.pseudo.bed <- args[1]
out.name <- args[2]
ref.dir <- args[3]
gene.length.threshold <- as.numeric(args[4])
intron.length.threshold <- as.numeric(args[5])

# Read in the data
introns.df <- read_delim(
        no.pseudo.bed,
        col_names = c("chr", "int.start", "int.end", "id.rank", "V5", "strand")
    ) %>%
    dplyr::select(-V5) %>%
    separate_wider_delim("id.rank", ": ", names = c("ensembl_transcript_id", "rank"))

# Get annotations
annotations <- read_delim(
        file = paste0(ref.dir, "/", "transcript_reference.tab")
    ) %>%
    dplyr::select(all_of(
        c(
            "ensembl_gene_id",
            "ensembl_transcript_id",
            "start_position",
            "end_position"
        )
    )) %>%
    dplyr::filter(ensembl_transcript_id %in% unique(introns.df$ensembl_transcript_id))

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
    write_delim(
        out.name,
        col_names = FALSE,
        delim = "\t"
    )

