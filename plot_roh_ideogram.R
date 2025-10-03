#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(data.table)
  library(gtools)
  library(ggideogram)
})

# === OPTIONS ===
option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Input BED with ROH and nonROH"),
  make_option(c("-c", "--chromsizes"), type = "character", help = "Chromosome sizes file (chr\\tsize)"),
  make_option(c("-o", "--outfile"), type = "character", default = "roh_ideogram.png", help = "Output plot [PNG]"),
  make_option(c("-p", "--prefix"), type = "character", default = NULL, help = "Plot title prefix")
)
opt <- parse_args(OptionParser(option_list = option_list))

# === COLOR PALETTE ===
roh_palette <- c(
  "2" = "#000081",  # 1.0–1.6 Mb
  "3" = "#00b3ff",  # 1.6–3.1 Mb
  "4" = "#7bff7b",  # 3.1–6.2 Mb
  "5" = "#ffc600",  # 6.2–12.5 Mb
  "6" = "#f30900",  # 12.5–25 Mb
  "7" = "#800000"   # >25 Mb
)

roh_labels <- c(
  "2" = "1.0–1.6 Mb", "3" = "1.6–3.1 Mb", "4" = "3.1–6.2 Mb",
  "5" = "6.2–12.5 Mb", "6" = "12.5–25 Mb", "7" = ">25 Mb"
)

# === LOAD CHROM SIZES ===
chromsizes <- fread(opt$chromsizes, header = FALSE, col.names = c("chrom", "end"))
chromsizes$start <- 0
chromsizes <- chromsizes[
  grepl("^(SUPER_|chr|[0-9]+$)", chrom) & 
    !grepl("^scaffold_", chrom)
]
chromsizes$Chromosome <- gsub("^chr|^SUPER_", "", chromsizes$chrom)
chromsizes$Chromosome <- factor(chromsizes$Chromosome, levels = rev(mixedsort(unique(chromsizes$Chromosome))))

# === LOAD BED FILE ===
roh <- fread(opt$input, header = FALSE, fill = TRUE)
setnames(roh, c("chrom", "start", "end", "status", "Class"))

roh <- roh[
  grepl("^(SUPER_|chr|[0-9]+$)", chrom) & 
    !grepl("^scaffold_", chrom)
]

# === Filter for ROH class ≥ 2 ===
roh <- roh[!(status == "ROH" & as.integer(Class) < 2)]

# === Prep coordinates and chromosome factor ===
roh$Chromosome <- gsub("^chr|^SUPER_", "", roh$chrom)
roh$Chromosome <- factor(roh$Chromosome, levels = levels(chromsizes$Chromosome))
roh$ymin <- roh$start
roh$ymax <- roh$end

# === Clean Class column ===
roh$Class <- ifelse(roh$status == "ROH", roh$Class, NA)
roh$Class <- factor(roh$Class, levels = names(roh_palette))
roh$Class <- droplevels(roh$Class)

# === PLOT ===
p <- ggplot() +
  geom_ideogram(
    data = chromsizes,
    aes(x = Chromosome, ymin = start, ymax = end, chrom = Chromosome),
    radius = unit(4, 'pt'), width = 0.5,
    linewidth = 1,
    chrom.col = "#252525", chrom.lwd = .5
  ) +
  geom_ideogram(
    data = roh[status == "nonROH"],
    aes(x = Chromosome, ymin = ymin, ymax = ymax, chrom = Chromosome),
    width = 0.4
  ) +
  geom_ideogram(
    data = roh[status == "ROH" & !is.na(Class)],
    aes(x = Chromosome, ymin = ymin, ymax = ymax, chrom = Chromosome, fill = Class),
    width = 0.4
  ) +
  scale_fill_manual(
    name = "ROH Length",
    values = roh_palette,
    labels = roh_labels,
    drop = TRUE,
    guide = guide_legend(
      title.position = "top",
      ncol = 1,
      override.aes = list(color = NA)
    )
  ) +
  coord_flip() +
  labs(
    title = if (!is.null(opt$prefix)) paste0("ROH Landscape: ", opt$prefix) else "ROH Landscape"
  ) +
  theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    legend.position = 'right')

ggsave(opt$outfile, plot = p, width = 15, height = 9, dpi = 700)
