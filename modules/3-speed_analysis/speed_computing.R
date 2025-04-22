library("tidyverse")
library("furrr")

MY_THEME <- theme(
    text = element_text(family = "Roboto"),
    # axis.title.x = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "grey50"),
    panel.grid = element_line(color = "#b4aea9"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linetype = "dashed"),
    panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
    plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
    legend.background = element_rect(fill = "#fbf9f4"),
    plot.title = element_text(
        family = "Roboto",
        size = 16,
        face = "bold",
        color = "#2a475e",
        margin = margin(b = 20)
    )
)

mean.threshold <- 25
R2.threshold <- 0.65
min_datapoints <- 10

base_path <- "/Users/varo/Desktop/frailty/out/coverages_validation"
output_dir <- "/Users/varo/Desktop/frailty/contrast"

# ---- Load Metadata ----
print("Welcome to the load_coverages script.")
print("Loading experimental design...")

# Read the processed samples and the experimental design

bam_files <- read.table(
    file.path(base_path, "bam_files.txt"), 
    header = FALSE
) %>%
    # process V2 to remove everything but "pool" + \d\d
    mutate(Pool = gsub(".*/bams/(SRR\\d\\d\\d\\d\\d\\d\\d).*", "\\1", V2)) %>%
    rename(c("Index" = "V1", "File" = "V2"))

design <- read.table(file.path(base_path, "design.csv"), header = FALSE) %>%
    mutate(Pool = gsub(".*/bams/(SRR\\d\\d\\d\\d\\d\\d\\d).*", "\\1", V1)) %>%
    left_join(bam_files, by = c("Pool" = "Pool"))

# Get the paths of the csv in the base path
csv_files <- list.files(base_path, pattern = "csv$", full.names = TRUE)

# Remove "design.csv" from the list of csv files
csv_files <- csv_files[!grepl("design.csv", csv_files)]
print("Done.")

# ---- Load Coverage ----
print("Loading coverage data...")

# Use read_csv without prompting for column types
read_data <- function(file) {
    readr::read_delim(
        file,
        skip = 1,
        col_names = FALSE,
        col_types = cols(.default = col_double())
        ) %>%
    as.matrix()
}

t1 <- Sys.time()
plan(multisession)
data <- future_map(csv_files, read_data)
names(data) <- csv_files
plan(sequential)
t2 <- Sys.time()

print(paste("Done.", "Elapsed time:", round(t2 - t1, 2), "seconds."))

print("Initial filtering and processing...")
# Filter out the files with less than min_datapoints rows
data <- data[sapply(data, function(x) nrow(x) >= min_datapoints)]

# Reverse matrix rows in "-" strand samples
minus_strand_idx <- which(endsWith(names(data), "-.csv"))
data[minus_strand_idx] <- lapply(data[minus_strand_idx], function(A) A[nrow(A):1, ])

# Filter low coverage introns
keep <- sapply(data, function(x) {
    mean(colMeans(x)) >= mean.threshold
})

data <- data[which(keep)]
print("Done.")
print(paste("Keeping", paste0(length(keep), "/", length(csv_files)), "samples."))

# ---- Compute Slope ----

# Function to compute a linear model for a column
run_lm <- function(column_data) {
    .data <- data.frame(y = column_data, x = 1:length(column_data))
    lm(y ~ x, data = .data)
}

# Apply to each matrix in the list
print("Computing linear models...")
plan(multisession) # Use multiple cores
lm_results <- future_map(data, ~ apply(.x, 2, run_lm))
plan(sequential) # Reset to single-threaded processing
print("Done.")

print("Extracting attributes...")
r2s <- lapply(lm_results, function(x) {
    # take the mean of the R2 values of the lm list
    mean(sapply(x, function(y) summary(y)$adj.r.squared))
})

slopes <- lapply(lm_results, function(x) {
    # get the slope of the lm list
    sapply(x, function(y) coef(y)[2])
})
print("Done.")

# ---- Filter ----
print("Filtering...")
keep <- which(r2s >= R2.threshold)
r2s <- r2s[keep]
slopes <- slopes[keep]
print(paste("Keeping", paste0(length(keep), "/", length(data)), "samples."))

# ---- Compute Speed ----

speeds <- lapply(slopes, function(x) {
    -1 / x
})

# take the stem from the file name
names(speeds) <- gsub(".*\\/(.*)\\.csv", "\\1", names(speeds))

# ---- Compare Speed ----

# bind the speeds into a data frame (name to column "intron", speed to column "speed", sublist name to column "sample")
speeds_df <- data.frame(
        intron = rep(names(speeds), sapply(speeds, length)),
        speed = unlist(speeds),
        sample = rep(bam_files$Pool, length(speeds))
    ) %>% 
    arrange(intron, sample) %>%
    mutate(
        type = ifelse(sample %in% design$Pool[design$V2 == "Old"], "Old", "Young"),
        int.id = rep(1:5, length(unique(names(speeds))) * 2)
    )

rownames(speeds_df) <- NULL

# Save the data frame
write_csv(speeds_df, file.path(output_dir, "speeds_VAL.csv"))
