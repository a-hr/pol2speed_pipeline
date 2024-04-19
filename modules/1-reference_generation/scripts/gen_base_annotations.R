require("biomaRt")
require("tidyverse")
library("dplyr")
require("GenomicRanges")
library("progress")

# ---- INPUTS ----
args <- commandArgs(trailingOnly = TRUE)

tmp.dir <- args[1]
out.dir <- args[2]
out.name <- args[3]

ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
chrs <- c(1:22, "MT", "X", "Y")
intron_thr <- 1

# ----- QUERY ANNOTATIONS (biomaRt) -----

# get transcription start site for each transcript
tss <-
    getBM(
        attributes = c(
            "transcription_start_site",
            "transcript_length",
            "chromosome_name",
            "strand",
            "ensembl_gene_id",
            "ensembl_transcript_id"
        ),
        mart = ensembl,
    ) %>%
    dplyr::rename(tss_start = transcription_start_site) %>%
    dplyr::filter(chromosome_name %in% chrs)

# get exon data for all the transcripts
elocs <- getBM(
    attributes = c(
        "exon_chrom_start",
        "exon_chrom_end",
        "rank",
        "chromosome_name",
        "strand",
        "ensembl_gene_id",
        "ensembl_transcript_id"
    ),
    mart = ensembl,
    filters = "chromosome_name",
    values = chrs
)

# get a list of gene ids with at least intron_thr + 1 exons
multi.exon.tc <- elocs %>%
    dplyr::filter(rank >= intron_thr + 1) %>%
    distinct(ensembl_transcript_id)

# sort by rank ascending and extract start and end coords as vecs to operate
valid_transcripts <- tss %>%
    dplyr::filter(
        ensembl_transcript_id %in% multi.exon.tc$ensembl_transcript_id &
        strand == 1
    ) %>%
    pull(ensembl_transcript_id)

# SAVE TABLES
write.table(
    valid_transcripts,
    file = paste0(tmp.dir, "/", "start_valid_transcripts.csv"),
    sep = ";",
    quote = F,
    row.names = F
)

write.table(
    elocs,
    file = paste0(tmp.dir, "/", "start_elocs.csv"),
    sep = ";",
    quote = F,
    row.names = F
)

write.table(
    tss,
    file = paste0(tmp.dir, "/", "start_tss.csv"),
    sep = ";",
    quote = F,
    row.names = F
)

# ----- COMPUTE FORWARD TABLE -----

introns <- data.frame(
    intron_start = c(),
    intron_end = c(),
    intron_len = c(),
    intron_rank = c(),
    intron_interval = c(),
    ensembl_transcript_id = c()
)

pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                       total = length(valid_transcripts),
                       width = 100)      # Width of the progress bar

for (sample in valid_transcripts) {
    pb$tick()
    
    # get the exons contained in the sample transcript
    sample.exons <- elocs %>%
        dplyr::filter(ensembl_transcript_id == sample) %>%
        dplyr::arrange(rank)
    
    # create vectors containing the start and end coordinates of exons
    ex.start <- sample.exons %>%
        dplyr::pull(exon_chrom_start)
    
    ex.end <- sample.exons %>%
        dplyr::pull(exon_chrom_end)
    
    # align both vectors to be able to subtract element wise
    in.end <- ex.start[-1]  # first start coord not needed
    in.start <- ex.end[-length(ex.end)]  # last start coord not needed
    in.len <- in.end - in.start
    
    # append to df
    ranks <- 1:length(in.len)
    introns <- rbind(
        introns,
        data.frame(
            intron_start = in.start,
            intron_end = in.end,
            intron_len = in.len,
            intron_rank = ranks,
            intron_interval = sapply(ranks, function (n) paste0(n, "~", n + 1)),
            ensembl_transcript_id = sample
        )
    )
}

# add chr to transcripts
introns <- merge(
    x = introns,
    y = tss[, c("ensembl_transcript_id", "chromosome_name")],
    by.x = "ensembl_transcript_id",
    by.y = "ensembl_transcript_id"
)

write.table(
    introns,
    file = paste0(tmp.dir, "/", "introns_multi-intron_fw.csv"),
    sep = ";",
    quote = F,
    row.names = F
)

# ----- COMPUTE REVERSE TABLE -----

# load annotations
# tss <- read.csv2("start_tss.csv")
# elocs <- read.csv2("start_elocs.csv")

multi.exon.tc <- elocs %>%
    dplyr::filter(rank >= intron_thr + 1) %>%
    distinct(ensembl_transcript_id)

valid_transcripts <- tss %>%
    dplyr::filter(
        ensembl_transcript_id %in% multi.exon.tc$ensembl_transcript_id &
            strand == -1
    ) %>%
    pull(ensembl_transcript_id)


introns <- data.frame(
    intron_start = c(),
    intron_end = c(),
    intron_len = c(),
    intron_rank = c(),
    intron_interval = c(),
    ensembl_transcript_id = c()
)

pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                       total = length(valid_transcripts),
                       complete = "=",   # Completion bar character
                       incomplete = "-", # Incomplete bar character
                       current = ">",    # Current bar character
                       clear = FALSE,    # If TRUE, clears the bar when finish
                       width = 100)      # Width of the progress bar


for (sample in valid_transcripts) {
    pb$tick()
    
    # get the exons contained in the sample transcript
    sample.exons <- elocs %>%
        dplyr::filter(ensembl_transcript_id == sample) %>%
        dplyr::arrange(rank)
    
    # create vectors containing the start and end coordinates of exons
    ex.start <- sample.exons %>%
        dplyr::pull(exon_chrom_start)
    
    ex.end <- sample.exons %>%
        dplyr::pull(exon_chrom_end)
    
    # align both vectors to be able to subtract element wise
    in.end <- ex.start[-length(ex.start)]
    in.start <- ex.end[-1]
    in.len <- in.end - in.start
    
    # append to df
    ranks <- 1:length(in.len)
    introns <- rbind(
        introns,
        data.frame(
            intron_start = in.start,
            intron_end = in.end,
            intron_len = in.len,
            intron_rank = ranks,
            intron_interval = sapply(ranks, function (n) paste0(n + 1, "~", n)),
            ensembl_transcript_id = sample
        )
    )
}

introns <- merge(
    x = introns,
    y = tss[, c("ensembl_transcript_id", "chromosome_name")],
    by = "ensembl_transcript_id"
)

write.table(
    introns,
    file = paste0(tmp.dir, "/", "introns_multi-intron_rv.csv"),
    sep = ";",
    quote = F,
    row.names = F
)

# ----- MERGE FW-RV -----
introns <- read.csv2("introns_multi-intron_fw.csv") %>%
    dplyr::bind_rows(read.csv2("introns_multi-intron_rv.csv"))

write.csv2(
    introns,
    paste0(out.dir, "/", out.name),
    row.names = F
)
