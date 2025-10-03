#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(tools)
  library(scales)
  library(ggrepel) # for non-overlapping labels
})

# === OPTIONS ===
option_list <- list(
  make_option(c("-a", "--angsd"), type = "character", help = "ANGSD heterozygosity file (Sample<tab>Heterozygosity)"),
  make_option(c("-c", "--corrected"), type = "character", help = "Corrected heterozygosity file from PLINK + callable genome size"),
  make_option(c("-o", "--outprefix"), type = "character", default = "het_correlation", help = "Output prefix")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$angsd) || is.null(opt$corrected)) {
  stop("Both --angsd and --corrected files must be provided.\n")
}

# === LOAD ANGSD ===
angsd <- fread(opt$angsd, header = TRUE)
if (!all(c("Sample", "Heterozygosity") %in% names(angsd))) {
  stop("ANGSD file must have columns: Sample, Heterozygosity")
}
angsd[, ind := gsub("\\.het$", "", Sample)]
setnames(angsd, "Heterozygosity", "angsd_het")
angsd <- angsd[, .(ind, angsd_het)]

# === LOAD CORRECTED ===
corrected <- fread(opt$corrected, header = TRUE)
if (!all(c("IID", "Corrected_Het") %in% names(corrected))) {
  stop("Corrected file must have columns: IID, Corrected_Het")
}
setnames(corrected, "IID", "ind")
setnames(corrected, "Corrected_Het", "corr_het")
corrected <- corrected[, .(ind, corr_het)]

# === MERGE ===
merged <- merge(angsd, corrected, by = "ind")
if (nrow(merged) == 0) {
  stop("No matching individuals found between ANGSD and corrected files.")
}

# === CORRELATION ===
pearson_cor <- cor.test(merged$angsd_het, merged$corr_het, method = "pearson")
spearman_cor <- cor.test(merged$angsd_het, merged$corr_het, method = "spearman")
r_value <- as.numeric(pearson_cor$estimate) # full precision

cat("=== Correlation Results ===\n")
cat(sprintf("Pearson correlation: R = %.15f, p-value = %.4g\n", r_value, pearson_cor$p.value))
cat(sprintf("Spearman correlation: rho = %.15f, p-value = %.4g\n", spearman_cor$estimate, spearman_cor$p.value))

# Save correlation table
fwrite(data.table(
  method = c("Pearson", "Spearman"),
  estimate = c(r_value, spearman_cor$estimate),
  p_value = c(pearson_cor$p.value, spearman_cor$p.value)
), paste0(opt$outprefix, "_correlation.tsv"), sep = "\t")

# === PLOT ===
p <- ggplot(merged, aes(x = angsd_het, y = corr_het)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", linetype = "dashed") +
  geom_text_repel(aes(label = ind), size = 3, max.overlaps = Inf) +
  theme_bw() +
  labs(
    title = "Correlation between ANGSD and Corrected Heterozygosity",
    x = "ANGSD heterozygosity",
    y = "Corrected heterozygosity"
  ) +
  scale_x_continuous(labels = number_format(accuracy = 0.0000001), expand = expansion(mult = 0.05)) +
  scale_y_continuous(labels = number_format(accuracy = 0.0000001), expand = expansion(mult = 0.05)) +
  annotate(
    "text",
    x = min(merged$angsd_het),
    y = max(merged$corr_het),
    label = paste0("R = ", format(r_value, scientific = FALSE, digits = 15)),
    hjust = 0, vjust = 1, size = 4
  )


ggsave(paste0(opt$outprefix, "_scatter.png"), p, width = 6, height = 5, dpi = 300)

cat(sprintf("\nResults written to:\n- %s_correlation.tsv\n- %s_scatter.png\n", opt$outprefix, opt$outprefix))
