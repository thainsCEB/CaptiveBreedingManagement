#!/usr/bin/env Rscript

# Inbreeding landscape from BCFtools ROH output (with optional metadata and PLINK comparison)

library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(gghalves)
library(patchwork)
library(wesanderson)
library(forcats)
library(purrr)
library(ggnewscale)
library(scales)

# Default theme
theme_emily <- function(base_size = 12, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = base_size, color = "black"),
      axis.title = element_text(size = base_size + 1, face = "bold", color = "black"),
      axis.line = element_line(linewidth = 0.4, color = "black"),
      axis.ticks = element_line(linewidth = 0.4, color = "black"),
      legend.text = element_text(size = base_size),
      legend.title = element_text(size = base_size + 1, face = "bold"),
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
      strip.text = element_text(size = base_size, face = "bold"),
      strip.background = element_rect(fill = "white", color = NA)
    )
}


# Arguments
args <- commandArgs(trailingOnly = TRUE)
if ("-h" %in% args || "--help" %in% args || length(args) < 4) {
  cat("\nUsage: Rscript inbreeding_landscape.R <roh_file> <chrom_size_file> <chr_list_file> <prefix> [sample_info_file] [plink_file] [recom_rate]\n\n")
  quit(save = "no")
}

roh_file <- args[1]
chrom_size_file <- args[2]
chr_list_file <- args[3]
prefix <- args[4]
sample_info_file <- ifelse(length(args) >= 5, args[5], NA)
plink_file <- ifelse(length(args) >= 6, args[6], NA)

# Chromosomes
chr_include <- scan(chr_list_file, what = character(), quiet = TRUE)
chr_autosomes <- chr_include[!grepl("(^[XYWZ]$|^chr[XYWZ]$|W|Z)", chr_include, ignore.case = TRUE)]
chr_lengths <- fread(chrom_size_file, header = FALSE) %>% filter(V1 %in% chr_autosomes)
auto_length <- sum(chr_lengths$V2) / 1000

# Sample metadata
sample_info <- if (!is.na(sample_info_file)) {
  s <- fread(sample_info_file, header = FALSE)
  if (ncol(s) < 2) stop("Sample info file must contain at least 2 columns: FID and Population")
  colnames(s)[1:2] <- c("FID", "Population")
  s
} else {
  NULL
}

# ROH
roh <- fread(roh_file)
expected_cols <- c("RG", "FID", "CHR", "Start", "End", "BP", "nSNPS", "Q")
if (!all(expected_cols %in% colnames(roh))) {
  colnames(roh) <- expected_cols  # fallback if no header
}
roh <- roh %>%
  mutate(across(c(Start, End, BP, nSNPS, Q), as.numeric))

roh <- roh %>%
  mutate(CHR = as.character(CHR)) %>%
  filter(CHR %in% chr_autosomes) %>%
  mutate(KB = BP / 1000) %>%
  filter(KB > 500)

roh <- if (!is.null(sample_info)) {
  left_join(roh, sample_info, by = "FID")
} else {
  roh %>% mutate(Population = "Unknown")
}

# Color palette
Populations <- unique(roh$Population)
col_palette <- setNames(wes_palette("FantasticFox1", max(3, length(Populations)), type = "continuous")[seq_along(Populations)], Populations)

# FROH
froh <- roh %>%
  group_by(FID, Population) %>%
  summarise(KBAVG = mean(KB), KBSUM = sum(KB), .groups = "drop") %>%
  mutate(percent_genome = (KBSUM / auto_length) * 100, FROH = KBSUM / auto_length)

# Summary
roh_summary <- roh %>%
  group_by(FID) %>%
  summarise(
    ROH_count = n(),
    ROH_sum_KB = sum(KB),
    ROH_mean_KB = mean(KB),
    ROH_median_KB = median(KB),
    .groups = "drop"
  ) %>%
  left_join(froh, by = "FID") %>%
  arrange(desc(FROH))

write.csv(roh_summary, paste0(prefix, "_froh_summary.csv"), row.names = FALSE, quote = FALSE)

# Violin plot
png(paste0(prefix, "_froh_violin.png"), width = 1800, height = 1200, res = 200)
print(
  ggplot(froh, aes(x = Population, y = percent_genome, fill = Population)) +
    geom_half_point(side = "l", shape = 21, stroke = 0.1, alpha = 0.5, size = 4) +
    geom_half_boxplot(side = "r", outlier.color = NA, width = 0.6, linewidth = 0.3, color = "black", alpha = 0.8) +
    scale_fill_manual(values = col_palette, name = "Population") +
    scale_y_continuous(labels = label_number(accuracy = 0.1)) +
    theme_emily() +
    xlab("Population") + ylab("% genome in ROH")
)
dev.off()

# FROH histogram
png(paste0(prefix, "_froh_histogram.png"), width = 1800, height = 1200, res = 200)
print(
  ggplot(froh, aes(FROH, fill = Population)) +
    geom_histogram(position = "identity", alpha = 0.5, bins = 15) +
    scale_fill_manual(values = col_palette) +
    theme_emily() +
    ylab("Individuals") + xlab(expression(italic(F["ROH"]))) +
    geom_vline(aes(xintercept = median(FROH)), linetype = "dashed", color = "grey40", linewidth = 0.7)
)
dev.off()

# ROH length histogram
png(paste0(prefix, "_roh_length_histogram.png"), width = 1800, height = 1200, res = 200)
print(
  roh %>%
    group_by(FID) %>%
    summarise(mean = mean(KB), .groups = "drop") %>%
    ggplot(aes(mean)) +
    geom_histogram(bins = 20, color = "white", fill = "grey70") +
    theme_emily() +
    ylab("Individuals") + xlab("Average ROH length per individual (Kb)")
)
dev.off()

# ROH class
roh_class <- roh %>%
  mutate(length_Mb = KB / 1000) %>%
  mutate(class = case_when(
    length_Mb >= 25.0 ~ 7,
    length_Mb >= 12.5 ~ 6,
    length_Mb >= 6.25 ~ 5,
    length_Mb >= 3.125 ~ 4,
    length_Mb >= 1.5625 ~ 3,
    length_Mb >= 0.78125 ~ 2,
    length_Mb >= 0.390625 ~ 1,
    TRUE ~ NA_real_
  )) %>%
  mutate(length_class = case_when(
    class == 1 ~ "0.4–0.8",
    class == 2 ~ "0.8–1.6",
    class == 3 ~ "1.6–3.1",
    class == 4 ~ "3.1–6.2",
    class == 5 ~ "6.2–12.5",
    class == 6 ~ "12.5–25",
    class == 7 ~ ">25",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(length_class)) %>%
  group_by(FID, length_class, Population) %>%
  summarise(prop_IBD = (sum(length_Mb) / (auto_length/1000)) * 100, .groups = "drop") %>%
  mutate(length_class = factor(length_class, levels = c("0.4–0.8", "0.8–1.6", "1.6–3.1", "3.1–6.2", "6.2–12.5", "12.5–25", ">25")))

png(paste0(prefix, "_roh_class_boxplot.png"), width = 1800, height = 1200, res = 200)
print(
  ggplot(roh_class, aes(length_class, prop_IBD, fill = Population)) +
    geom_half_point(side = "l", shape = 21, stroke = 0.1, alpha = 0.5, size = 2) +
    geom_half_boxplot(side = "r", outlier.color = NA, width = 0.8, linewidth = 0.3, color = "black", alpha = 0.8) +
    scale_fill_manual(values = col_palette) +
    theme_emily() +
    xlab("ROH length (Mb)") + ylab("% genome in ROH")
)
dev.off()

# ROH count vs sum
roh_nsum <- roh %>%
  group_by(FID, Population) %>%
  summarise(n = n(), sum = sum(KB), .groups = "drop")

model <- lm(n ~ sum, data = roh_nsum)
cat("ROH count ~ sum model:\n")
print(coef(model))

png(paste0(prefix, "_roh_count_vs_sum.png"), width = 1800, height = 1200, res = 200)
print(
  ggplot(roh_nsum, aes(sum, n, color = Population)) +
    geom_point(size = 3) +
    geom_abline(intercept = coef(model)[1], slope = coef(model)[2]) +
    scale_color_manual(values = col_palette) +
    theme_emily() +
    scale_x_continuous(labels = label_number(accuracy = 1)) +
    xlab("Sum of ROH lengths (Kb)") + ylab("Total number of ROH")
)
dev.off()

# PLINK comparison
if (!is.na(plink_file)) {
  plink <- fread(plink_file)
  plink <- plink %>%
    filter(CHR %in% chr_autosomes) %>%
    mutate(KB = KB) %>%
    group_by(FID) %>%
    summarise(KBSUM_plink = sum(KB), .groups = "drop") %>%
    mutate(FROH_plink = KBSUM_plink / auto_length) %>%
    left_join(froh, by = "FID")
  
  r2 <- summary(lm(FROH ~ FROH_plink, data = plink))$r.squared
  
  p_corr <- ggplot(plink, aes(FROH, FROH_plink)) +
    geom_point(aes(color = Population), size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", color = "grey40", se = FALSE) +
    scale_color_manual(values = col_palette) +
    theme_emily() +
    xlab(expression("BCFtools " * F[ROH])) +
    ylab(expression("PLINK " * F[ROH])) +
    annotate("text", x = 0.1, y = max(plink$FROH_plink, na.rm = TRUE) * 0.9,
             label = paste("R² =", round(r2, 3)), size = 4)
  
  ggsave(paste0(prefix, "_plink_vs_bcftools.png"), p_corr, width = 6, height = 5)
}

# --- Timing of Inbreeding Section ---
# Pull recombination rate (cM/Mb) and generation time (years) from args
recomb_rate <- as.numeric(ifelse(length(args) >= 7, args[7], NA))
generation_time_years <- as.numeric(ifelse(length(args) >= 8, args[8], NA))

if (is.na(recomb_rate) || is.na(generation_time_years)) {
  stop("Error: Please provide recombination rate (cM/Mb) and generation time (years) as the 7th and 8th arguments.")
}

# Calculate timing for each ROH
roh_class <- roh %>%
  mutate(length_Mb = KB / 1000) %>%
  mutate(length_class = case_when(
    length_Mb >= 25.0 ~ ">25",
    length_Mb >= 12.5 ~ "12.5-25",
    length_Mb >= 6.25 ~ "6.2-12.5",
    length_Mb >= 3.125 ~ "3.1-6.2",
    length_Mb >= 1.5625 ~ "1.6-3.1",
    length_Mb >= 0.78125 ~ "0.8-1.6",
    length_Mb >= 0.390625 ~ "0.4-0.8",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(length_class)) %>%
  mutate(
    length_class = factor(length_class,
                          levels = c("0.4-0.8", "0.8-1.6", "1.6-3.1",
                                     "3.1-6.2", "6.2-12.5", "12.5-25", ">25"),
                          ordered = TRUE
    ),
    generations_since_inbreeding = 1 / (2 * (recomb_rate / 100) * length_Mb),
    years_since_inbreeding = generations_since_inbreeding * generation_time_years
  )

# Summary stats per class
class_stats <- roh_class %>%
  group_by(length_class) %>%
  summarise(
    mean_gen = mean(generations_since_inbreeding, na.rm = TRUE),
    mean_years = mean(years_since_inbreeding, na.rm = TRUE),
    min_gen = min(generations_since_inbreeding, na.rm = TRUE),
    max_gen = max(generations_since_inbreeding, na.rm = TRUE),
    min_years = min(years_since_inbreeding, na.rm = TRUE),
    max_years = max(years_since_inbreeding, na.rm = TRUE),
    .groups = "drop"
  )

# Max y for annotations
y_max <- max(roh_class$generations_since_inbreeding, na.rm = TRUE)

# Build plot with mean generations on top and mean years below, with background for readability
p_timing <- ggplot(roh_class, aes(x = length_class, y = generations_since_inbreeding, fill = Population)) +
  geom_boxplot(outlier.color = NA, width = 0.6, color = "black", alpha = 0.8) +
  
  # Mean annotation (two lines, white background for visibility)
  geom_label(
    data = class_stats,
    aes(
      x = length_class,
      y = mean_gen,
      label = paste0(round(mean_gen, 1), "g\n", round(mean_years, 0), "y")
    ),
    vjust = -0.6,
    fontface = "bold",
    size = 3,
    fill = "white",
    label.size = 0.2,
    label.r = unit(0.15, "lines"),
    inherit.aes = FALSE
  ) +
  
  # Top annotations for recombination rate & generation time
  annotate("text",
           x = length(class_stats$length_class) - 0.2,
           y = y_max * 1.12,
           label = paste("Mean autosomal recombination rate =", recomb_rate, "cM/Mb"),
           hjust = 1, size = 5, fontface = "italic") +
  annotate("text",
           x = length(class_stats$length_class) - 0.2,
           y = y_max * 1.06,
           label = paste("Generation time =", generation_time_years, "years"),
           hjust = 1, size = 5, fontface = "italic") +
  
  scale_fill_manual(values = col_palette) +
  theme_emily() +
  xlab("ROH length class (Mb)") +
  ylab("Generations since inbreeding") +
  scale_y_continuous(
    sec.axis = sec_axis(~ . * generation_time_years,
                        name = "Years since inbreeding"),
    labels = label_number(accuracy = 0.1)
  )

# Save high-res output
ggsave(
  filename = paste0(prefix, "_roh_class_timing_boxplot.png"),
  plot = p_timing,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
)
