#!/usr/bin/env Rscript

## Procrustes PCA plotter for fastNGSadmix-style covariances
## - Uniform point sizes across shapes
## - Classic_10/20 palette via paletteer
## - Projected groups use shapes: triangle, diamond, square, star, cross (cycled)
## - Auto-install required packages
## - Optional labels for projected points via --label-projected (5th arg)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  message("Usage:\n",
          "  Rscript procrustes_pca.R <covar_dir> <ref.fam> <sample_group.tsv> <out_prefix> [--label-projected]\n",
          "Where <sample_group.tsv> has two columns: <sample_prefix> <group_label>\n",
          "Optional 5th arg: --label-projected   # add text labels to projected points")
  quit(save = "no", status = 2)
}

covar_dir  <- args[1]
fam_file   <- args[2]
group_tbl  <- args[3]
out_prefix <- args[4]

## Optional 5th argument: --label-projected / label-projected / TRUE / T / 1
do_label <- FALSE
if (length(args) >= 5) {
  flag <- tolower(gsub("^--", "", args[5]))
  do_label <- flag %in% c("label-projected", "true", "t", "1", "yes", "y")
}

## -------------------- Auto-install missing packages --------------------
ensure_pkgs <- function(pkgs) {
  userlib <- Sys.getenv("R_LIBS_USER")
  if (is.null(userlib) || userlib == "") {
    userlib <- file.path(path.expand("~"), "R",
                         paste0(R.version$platform, "-library"),
                         paste0(R.version$major, ".", strsplit(R.version$minor, "[.]")[[1]][1]))
  }
  if (!dir.exists(userlib)) dir.create(userlib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(unique(c(userlib, .libPaths())))
  if (is.null(getOption("repos")) || getOption("repos")["CRAN"] %in% c("@CRAN@", NULL, "")) {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
  }
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      message("Installing missing package: ", p)
      install.packages(p, lib = userlib, dependencies = TRUE)
    }
  }
}

base_pkgs <- c("vegan", "ggplot2", "paletteer")
if (do_label) base_pkgs <- c(base_pkgs, "ggrepel")
ensure_pkgs(base_pkgs)

suppressPackageStartupMessages({
  library(parallel)
  library(vegan)
  library(ggplot2)
  library(paletteer)
  if (do_label) library(ggrepel)
})

## -------------------- Inputs & reference setup --------------------
covar_files_names <- list.files(covar_dir, pattern = "_covar.txt$", full.names = FALSE)
if (length(covar_files_names) == 0L) stop("No *_covar.txt files found in: ", covar_dir)
covar_files_paths <- file.path(covar_dir, covar_files_names)

covar0 <- read.table(covar_files_paths[1], header = TRUE, check.names = FALSE)
if (nrow(covar0) < 2) stop("Covariance matrix appears too small: ", covar_files_paths[1])

fam <- read.table(fam_file, stringsAsFactors = FALSE)
if (!is.null(rownames(covar0))) {
  ref_ids <- rownames(covar0)[1:(nrow(covar0) - 1)]
  fam <- fam[fam$V2 %in% ref_ids, , drop = FALSE]
  fam <- fam[order(match(fam$V2, ref_ids)), , drop = FALSE]
}

ev <- eigen(as.matrix(covar0))
explained_var <- 100 * ev$values / sum(ev$values)

## -------------------- Sample -> group table --------------------
grp <- read.table(group_tbl, header = FALSE, stringsAsFactors = FALSE)
if (ncol(grp) < 2) stop("Group table must have 2 columns: <sample_prefix> <group_label>")
colnames(grp)[1:2] <- c("sample_prefix", "group_label")

file_prefixes <- sub("_covar.txt$", "", basename(covar_files_paths))
files_df <- data.frame(filesName = file_prefixes, files = covar_files_paths, stringsAsFactors = FALSE)

keep <- intersect(grp$sample_prefix, files_df$filesName)
if (length(keep) == 0) stop("No overlap between group table and *_covar.txt prefixes.")
grp <- grp[grp$sample_prefix %in% keep, , drop = FALSE]
files_df <- files_df[files_df$filesName %in% keep, , drop = FALSE]
grp <- merge(files_df, grp, by.x = "filesName", by.y = "sample_prefix", sort = FALSE)

## -------------------- Compute per-sample PC1/PC2 --------------------
message("Creating temporary directory: procustesPCs")
mk_status <- system("mkdir procustesPCs")  # 0 if created, 1 if exists

pc_to_file <- function(f) {
  r  <- read.table(f, header = TRUE, check.names = FALSE)
  ee <- eigen(as.matrix(r))
  outname <- sub("_covar.txt$", "", basename(f))
  write.table(ee$vectors[, 1:2],
              file = file.path("procustesPCs", outname),
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  invisible(NULL)
}

run_pc <- function(files, cores) {
  if (.Platform$OS.type == "windows") {
    lapply(files, pc_to_file)
  } else {
    parallel::mclapply(files, pc_to_file, mc.cores = min(max(1, floor(cores / 2)), 20))
  }
  invisible(NULL)
}
num_cores <- parallel::detectCores()
invisible(run_pc(grp$files, num_cores))

pcs_paths <- file.path("procustesPCs", grp$filesName)
if (any(!file.exists(pcs_paths))) stop("Missing PC files in procustesPCs/; check permissions and inputs.")

message("Procrustes-aligning ", length(pcs_paths), " projected samples")

m0 <- read.table(pcs_paths[1], header = FALSE)
if (nrow(m0) < 2) stop("PC matrix too small in: ", pcs_paths[1])
m0 <- m0[1:(nrow(m0) - 1), , drop = FALSE]

predict_proj <- function(pc_file) {
  m2 <- read.table(pc_file, header = FALSE)
  r  <- vegan::procrustes(m0, m2[1:(nrow(m2) - 1), , drop = FALSE])
  as.numeric(predict(r, m2[nrow(m2), , drop = FALSE]))
}
res_mat <- t(vapply(pcs_paths, predict_proj, numeric(2)))
colnames(res_mat) <- c("PC1", "PC2")

xlim <- range(c(m0[, 1], res_mat[, 1]))
ylim <- range(c(m0[, 2], res_mat[, 2]))

## -------------------- Build plotting tables --------------------
ref_df <- data.frame(
  Species = fam$V1,
  Sample  = fam$V2,
  PC1     = m0[, 1],
  PC2     = m0[, 2],
  Group   = "Reference",
  stringsAsFactors = FALSE
)

proj_df <- data.frame(
  Species = grp$group_label,
  Sample  = grp$filesName,
  PC1     = res_mat[, 1],
  PC2     = res_mat[, 2],
  Group   = grp$group_label,
  stringsAsFactors = FALSE
)

df <- rbind(ref_df, proj_df)

## -------------------- Shapes (uniform sizes later in ggplot) --------------------
proj_groups <- setdiff(unique(df$Group), "Reference")
shape_seq <- c(17, 18, 15, 8, 4)  # triangle, diamond, square, star, cross
shape_map <- c("Reference" = 16)
if (length(proj_groups) > 0) {
  shape_map[proj_groups] <- shape_seq[((seq_along(proj_groups) - 1) %% length(shape_seq)) + 1]
}

## -------------------- Colors: Classic_10 / Classic_20 --------------------
species_levels <- unique(df$Species)
if (length(proj_groups) >= 1) {
  target <- proj_groups[1]
  species_levels <- c(target, setdiff(species_levels, target))
}
df$Species <- factor(df$Species, levels = species_levels)

n_cols <- length(species_levels)
if (n_cols <= 10) {
  pal_vec <- as.vector(paletteer::paletteer_d("ggthemes::Classic_10", n = n_cols))
} else {
  base_pal <- as.vector(paletteer::paletteer_d("ggthemes::Classic_20", n = min(20, n_cols)))
  pal_vec  <- base_pal[((seq_len(n_cols) - 1) %% length(base_pal)) + 1]
}
names(pal_vec) <- species_levels

## -------------------- Save data tables --------------------
write.table(df[, c("Species", "Sample", "PC1", "PC2", "Group")],
            paste0(out_prefix, ".txt"),
            col.names = TRUE, row.names = FALSE, quote = FALSE)

legacy_points <- rbind(
  cbind(ref_df$Species, ref_df$Sample, ref_df$PC1, ref_df$PC2),
  cbind(proj_df$Species, proj_df$Sample, proj_df$PC1, proj_df$PC2)
)
write.table(legacy_points, paste0(out_prefix, ".legacy.txt"),
            col.names = FALSE, row.names = FALSE, quote = FALSE)

## -------------------- Plot (uniform sizes + optional labels) --------------------
set.seed(1)  # for stable ggrepel jittering/placement

p <- ggplot(df, aes(x = PC1, y = PC2)) +
  geom_point(
    data = subset(df, Group == "Reference"),
    aes(color = Species, shape = Group),
    size = 3.0, stroke = 0.9, alpha = 0.95
  ) +
  geom_point(
    data = subset(df, Group != "Reference"),
    aes(color = Species, shape = Group),
    size = 3.0, stroke = 0.9, alpha = 0.98
  ) +
  scale_colour_manual(values = pal_vec) +
  scale_shape_manual(values = shape_map) +
  guides(
    shape  = guide_legend(override.aes = list(size = 3.0, stroke = 0.9, alpha = 1)),
    colour = guide_legend(override.aes = list(size = 3.0, alpha = 1))
  ) +
  coord_cartesian(xlim = xlim, ylim = ylim, expand = TRUE) +
  xlab(paste0("PC1 (", sprintf("%.2f", explained_var[1]), "%)")) +
  ylab(paste0("PC2 (", sprintf("%.2f", explained_var[2]), "%)")) +
  ggtitle("Procrustes PCA") +
  theme_minimal(base_size = 12) +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

if (do_label) {
  p <- p +
    ggrepel::geom_text_repel(
      data = subset(df, Group != "Reference"),
      aes(label = Sample, color = Species),
      size = 3.0, max.overlaps = 50,
      box.padding = 0.35, point.padding = 0.25,
      min.segment.length = 0,
      segment.alpha = 0.6, segment.size = 0.3,
      seed = 1,
      show.legend = FALSE
    )
}

ggsave(filename = paste0(out_prefix, ".png"),
       plot = p, width = 7.5, height = 5.8, dpi = 300, bg = "white")

## -------------------- Cleanup --------------------
if (mk_status == 0) {
  system("rm -r procustesPCs")
}
