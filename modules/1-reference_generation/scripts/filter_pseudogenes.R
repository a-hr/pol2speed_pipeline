library(biomaRt)
library(dplyr)

# ---- INPUTS ----
args <- commandArgs(trailingOnly = TRUE)

base.reference <- args[1]
out.dir <- args[2]
out.name <- args[3]

# ----- QUERY ANNOTATIONS (biomaRt) -----
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

intron.df <- read.csv2(base.reference)

transcript.list <- unique(intron.df$ensembl_transcript_id)

# launch the request
annotations <- getBM(
    filters = "ensembl_transcript_id",
    values = transcript.list,
    attributes = c(
        "ensembl_gene_id",
        "ensembl_transcript_id",
        "external_gene_name",
        "gene_biotype"
    ),
    mart = ensembl
)

# filter out the transcripts whose biotype is "pseudogene"
unique(annotations$gene_biotype)
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
    filter(!ensembl_transcript_id %in% transcripts.remove)

write.csv2(
    intron.df2,
    file = paste0(out.dir, "/", out.name),
    row.names = F
)
