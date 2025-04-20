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

# ---- Load Metadata ----
print("Welcome to the load_coverages script.")
print("Loading experimental design...")

# Read the processed samples and the experimental design
base_path <- "/Users/varo/Desktop/frailty/out/coverages_validation"
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
write_csv(speeds_df, file.path("/Users/varo/Desktop/frailty/contrast", "speeds_VAL.csv"))

speeds_df <- read.csv(file.path("/Users/varo/Desktop/frailty/contrast", "speeds_BIOD.csv"))
# Perform the Wilcox rank test

wt <- wilcox.test(
    speed ~ type,
    data = speeds_df,
    paired = T,
    alternative = "greater"
)

print(wt)


# ---- Prepare the data for plotting----

# pool-level data
comparison_df <- speeds_df %>%
    group_by(intron, int.id) %>%
    summarise(
        speed.F = first(speed),
        speed.R = last(speed),
        sex = first(sex)
    ) %>%
    mutate(
        ratio = speed.F / speed.R,
        log2_ratio = log2(ratio),
        chr = str_extract(intron, "(?<=chr)(\\d{1,2}|[XYM])"),
        strand = substr(intron, nchar(intron), nchar(intron)),
        chr = factor(chr, levels = c(1:22, "X", "Y", "M"))
    )

# intron-level data
comparison_df2 <- speeds_df %>%
    group_by(intron, sex) %>%
    summarise(
        # take the median of the type "F" and "R" speeds
        speed.F = median(speed[type == "Frail"]),
        speed.R = median(speed[type == "Robust"]),
        sex = first(sex)
    ) %>%
    mutate(
        ratio = speed.F / speed.R,
        log2_ratio = log2(ratio),
        chr = str_extract(intron, "(?<=chr)(\\d{1,2}|[XYM])"),
        strand = substr(intron, nchar(intron), nchar(intron)),
        chr = factor(chr, levels = c(1:22, "X", "Y", "M"))
    )


# ---- Plots by chromosome ----

# global F/R speed ratio (histogram), colored by chromosome
comparison_df %>%
    ggplot(aes(x = log2_ratio, fill = chr)) +
    geom_histogram(bins = 20, alpha = .6) +
    theme_minimal() +
    xlab("Frail / Robust Speed Ratio") +
    ylab("Count") +
    theme(legend.position = "none")

comparison_df2 %>%
    ggplot(aes(x = log2_ratio, fill = chr)) +
    geom_histogram(bins = 20, alpha = .6) +
    theme_minimal() +
    xlab("Frail / Robust Speed Ratio") +
    ylab("Count") +
    theme(legend.position = "none")

# F/R speed ratio mean by chromosome
comparison_df %>%
    group_by(chr) %>%
    summarize(
        mean_ratio = mean(log2_ratio),
        nsamples = n()
    ) %>%
    filter(nsamples > 20) %>%  # a nivel pool (1 intron)
    ggplot(aes(x = chr, y = mean_ratio, fill = chr)) +
    geom_bar(alpha = .6, position = "dodge", stat = "identity") +
    # annotate the number of introns in each chromosome
    geom_text(aes(label = nsamples), vjust = -1) +
    theme_minimal() +
    xlab("Chromosome") +
    ylab("Median Frail / Robust Speed Ratio") +
    theme(legend.position = "none")

comparison_df2 %>%
    group_by(chr) %>%
    summarize(
        mean_ratio = mean(log2_ratio),
        nsamples = n()
    ) %>%
    filter(nsamples > 2) %>%  # a nivel intron/sex (1 intron)
    ggplot(aes(x = chr, y = mean_ratio, fill = chr)) +
    geom_bar(alpha = .6, position = "dodge", stat = "identity") +
    # annotate the number of introns in each chromosome
    geom_text(aes(label = nsamples), vjust = -1) +
    theme_minimal() +
    xlab("Chromosome") +
    ylab("Median Frail / Robust Speed Ratio") +
    theme(legend.position = "none")

# F/R speed ratio mean by chromosome + sex
comparison_df %>%
    group_by(chr, sex) %>%  # Include 'sex' in grouping
    summarize(
        mean_ratio = mean(log2_ratio),
        nsamples = n(),
        .groups = "drop"
    ) %>%
    ggplot(aes(x = chr, y = mean_ratio, fill = sex)) +  # Map 'sex' to fill
    geom_bar(alpha = .6, position = position_dodge(width = 0.8), stat = "identity") +  # Dodge bars by 'sex'
    # Annotate the number of introns in each chromosome, adjusted for dodge
    geom_text(aes(x = chr, y = max(mean_ratio) + 0.05, label = nsamples),  # Position annotations above bars
              inherit.aes = FALSE,
              alpha = 0.5,
              vjust = -0.5) +
    theme_minimal() +
    xlab("Chromosome") +
    ylab("Mean log2(Frail / Robust) Speed Ratio across chromosomes") +
    theme(legend.position = "top") +
    scale_fill_brewer(palette = "Set2")  # Optional: nice color palette

comparison_df2 %>%
    group_by(chr, sex) %>%  # Include 'sex' in grouping
    summarize(
        mean_ratio = mean(log2_ratio),
        nsamples = n(),
        .groups = "drop"
    ) %>%
    ggplot(aes(x = chr, y = mean_ratio, fill = sex)) +  # Map 'sex' to fill
    geom_bar(alpha = .6, position = position_dodge(width = 0.8), stat = "identity") +  # Dodge bars by 'sex'
    # Annotate the number of introns in each chromosome, adjusted for dodge
    geom_text(aes(x = chr, y = max(mean_ratio) + 0.05, label = nsamples),  # Position annotations above bars
              inherit.aes = FALSE,
              alpha = 0.5,
              vjust = -0.5) +
    theme_minimal() +
    xlab("Chromosome") +
    ylab("Mean log2(Frail / Robust) Speed Ratio across chromosomes") +
    theme(legend.position = "top") +
    scale_fill_brewer(palette = "Set2")  # Optional: nice color palette

# ---- Plots by sex ----

# ratio histogram

medians <- comparison_df %>%
    group_by(sex) %>%
    summarize(median_ratio = median(log2(ratio)))

comparison_df %>%
    ggplot(aes(x = log2(ratio), fill = sex)) +
    geom_histogram(bins = 20, alpha = .6, position = "identity") +
    theme_minimal() +
    # show median for each group
    geom_vline(data = medians, aes(xintercept = median_ratio, col = sex), linetype = "dashed", linewidth = 1) +
    xlab("log2(Frail / Robust) Speed Ratio") +
    ylab("Count") +
    theme(legend.position = "top")

medians <- comparison_df2 %>%
    group_by(sex) %>%
    summarize(median_ratio = median(log2(ratio)))

comparison_df2 %>%
    ggplot(aes(x = log2(ratio), fill = sex)) +
    geom_histogram(bins = 20, alpha = .6, position = "identity") +
    theme_minimal() +
    # show median for each group
    geom_vline(data = medians, aes(xintercept = median_ratio, col = sex), linetype = "dashed", linewidth = 1) +
    xlab("log2(Frail / Robust) Speed Ratio") +
    ylab("Count") +
    theme(legend.position = "top")

wt.sex <- wilcox.test(
    speed ~ sex,
    data = speeds_df,
    paired = T,
    alternative = "two.sided"
)
wt.sex

wt.sex2 <- wilcox.test(
    speed ~ sex,
    data = data2,
    paired = T,
    alternative = "two.sided"
)
wt.sex2

# ---- Plot by genomic position ----
df <- comparison_df %>%
    rowwise() %>%
    mutate(
        start.pos = str_extract(intron, "(?<=chr(?:[1-9]|1[0-9]|2[0-2]|[XYM])_)\\d+"),
        end.pos = str_extract(intron, "(?<=chr(?:[1-9]|1[0-9]|2[0-2]|[XYM])_\\d{7,9}_)\\d+"),
        int.len = as.numeric(end.pos) - as.numeric(start.pos)
    )

# scatter plot of speed vs start position (global speeds) => speeds are not directly comparable, requires ratio
df %>%
    ggplot(aes(x = as.numeric(start.pos), y = ratio, col = type)) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    xlab("Start Position") +
    ylab("Pol2 Speed") +
    theme(legend.position = "top")

# ---- Pair plot ----
p.val <- wt$p.value

plot.pair <- speeds_df %>%
    mutate(
        pairing = paste(intron, int.id, sep = "-"),
        group = factor(type, levels = c("Robust", "Frail"))
    ) %>% 
    ggplot(aes(x = group, y = log10(speed + 1))) +
    geom_boxplot(
        width = 0.3,
        size = 1.5,
        fatten = 1.5,
        colour = "gray70",
        outlier.shape = NA
    ) +
    geom_point(colour = "red",
               size = .2,
               alpha = 0.5) +
    geom_line(
        aes(group = pairing),
        colour = "gray30",
        linetype = "11",
        linewidth = .1
    ) +
    theme_minimal() +
    xlab("Condition") +
    ylab("log10(Pol2 Speed)") +
    annotate("text", x = .7, y = .8, label = paste("Wilcoxon p.val:", round(p.val, 5)), vjust = -1)
plot.pair

data.grouped <- speeds_df %>%
    mutate(
        pairing = paste(intron, int.id, sep = "-"),
        group = factor(type, levels = c("Robust", "Frail"))
    ) %>% 
    dplyr::select(-c(sample, type, int.id)) %>% 
    pivot_wider(names_from = group, values_from = speed, names_prefix = "speed.") %>%
    mutate(speed.dif = (speed.Frail - speed.Robust)/speed.Robust, sense = if_else(speed.dif >= 0, "Up", "Down"))

plot.sense <- data.grouped %>%
    group_by(intron) %>%
    summarise(
        net.speed.change = mean(speed.dif),
        sense = if_else(net.speed.change >= 0, "Up", "Down"),
    ) %>%
    ggplot(aes(x = sense, fill = sense)) +
    geom_bar(stat = "count", width = .1) +
    theme_minimal() +
    xlab("Speed Change Direction (Frail - Robust)") +
    ylab("Intron Count") +
    theme(
        legend.position = "none",
    )

avg.speed.change <- data.grouped %>%
    select(-c(speed.Frail, speed.Robust)) %>%
    group_by(intron) %>%
    summarise(
        net.speed.change = mean(speed.dif),
    ) %>% 
    summarise(
        avg.speed.change = mean(net.speed.change),
        sd.speed.change = sd(net.speed.change)
    )

layout <- "
##AAA##
#BBBBB#
#BBBBB#
#BBBBB#
#BBBBB#
"

plot.sense + plot.pair + plot_layout(design = layout, heights = c(1, 2)) +
    plot_annotation(
        title = "Pol2 Elongation Speed Change in Frail vs Robust individuals",
        theme = MY_THEME
    ) &
    MY_THEME


ggsave("final_plot.png", plot = final.plot, device = "png", height = 800, width = 400, units = "px")

# ---- Pair plot grouped by intron ----
p.val <- wt$p.value

grouped.df <- speeds_df %>%
    group_by(intron, type, sex) %>%
    summarise(
        speed = median(speed),
        .groups = "drop"
    ) %>% 
    mutate(
        group = factor(type, levels = c("Robust", "Frail")),
    )

intron.data <- grouped.df %>%
    group_by(intron, sex) %>%
    summarise(
        speed_mean = mean(speed),
        speed_max = max(abs(speed)),
        .groups = "drop"
    ) %>%
    mutate(
        key = paste0(intron, sex)
    ) %>%
    select(c(key, speed_mean, speed_max))

grouped.df <- grouped.df %>%
    mutate(
        key = paste0(intron, sex)
    ) %>%
    left_join(intron.data, by = "key") %>%
    mutate(
        speed = speed / speed_max,
        group = factor(type, levels = c("Robust", "Frail")),
        pairing = paste0(intron, sex)
    ) %>%
    filter(
        speed > 0
    )

plot.pair <- grouped.df %>%
    ggplot(aes(x = group, y = speed)) +
    geom_boxplot(
        width = 0.3,
        size = 1.5,
        fatten = 1.5,
        colour = "gray70",
        outlier.shape = NA
    ) +
    geom_point(
        aes(col = sex),
        size = 2,
    ) +
    geom_line(
        aes(group = pairing, col = sex),
        linetype = "11",
        linewidth = .5
    ) +
    theme_minimal() +
    xlab("Condition") +
    ylab("RNA-Pol II Speed (scaled)") +
    annotate("text", x = .7, y = .8, label = paste("Wilcoxon p.val:", round(p.val, 5)), vjust = -1)

plot.sense <- grouped.df %>% 
    dplyr::select(c(key, sex, group, speed)) %>% 
    pivot_wider(names_from = group, values_from = speed, names_prefix = "speed.") %>% 
    mutate(
        net.speed.change = speed.Frail - speed.Robust,,
        sense = if_else(net.speed.change >= 0, "Up", "Down"),
    ) %>% 
    ggplot(aes(x = sense, fill = sex)) +
    geom_bar(stat = "count", width = .3, position = position_dodge(width = 1)) +
    theme_minimal() +
    xlab("Speed Change Direction (Frail - Robust)") +
    ylab("Intron Count") +
    theme(
        legend.position = "none",
    )

avg.speed.change <- data.grouped %>%
    select(-c(speed.Frail, speed.Robust)) %>%
    group_by(intron) %>%
    summarise(
        net.speed.change = mean(speed.dif),
    ) %>% 
    summarise(
        avg.speed.change = mean(net.speed.change),
        sd.speed.change = sd(net.speed.change)
    )

layout <- "
##AAA##
#BBBBB#
#BBBBB#
#BBBBB#
#BBBBB#
"

plot.sense + plot.pair + plot_layout(design = layout, heights = c(1, 2)) +
    plot_annotation(
        title = "RNA-Pol II Elongation Speed Change in Frail vs Robust individuals",
        theme = MY_THEME
    ) &
    MY_THEME


ggsave("paired-plot BIOD.png", plot = final.plot, device = "png", height = 800, width = 400, units = "px")
