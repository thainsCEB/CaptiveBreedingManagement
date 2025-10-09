#!/usr/bin/env Rscript

## Procrustes PCA plotter for fastNGSadmix-style covariances
## - Shapes from group (col 2, legend "Group")
## - Colors from taxa/ref (col 3, legend "Taxa")
## - Optional ellipses per subgroup (group × taxa):
##     --ellipses (default OFF)
##     --exclude-ellipse=<csv>  # exact col2 or "group × taxa" OR substring in col3
##     --marker[=<color>]       # draw ALL ellipses in one color (e.g., --marker, --marker=#FF00FF)
##     --ellipse-style <solid|dashed|dotted>
##     --ellipse-lwd <num>      # line width (aliases: --ellipse-width, --ellipse-size)
## - Point sizing:
##     --point-size <num>       # ggplot2 size for both ref & projected points
##     --point-stroke <num>     # outline stroke for both point layers
## - Title via --title (auto-wrapped), position via --title-align=left|center|right|0..1
## - Text controls:
##     --font / --font-family "<family>"
##     --title-size <pt> --axis-title-size <pt> --axis-text-size <pt>
##     --legend-title-size <pt> --legend-text-size <pt> --label-size <pt>
## - Output sizing via --width, --height (inches) and --dpi
## - Auto-install packages; only PCA figure + data tables (no barplot)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  message("Usage:\n",
          "  Rscript procrustes_pca.R <covar_dir> <ref.fam> <sample_group.tsv> <out_prefix> [options]\n",
          "Where <sample_group.tsv> (no header) has 3 columns: <sample_prefix> <group_label> <ref_taxon_for_color>\n",
          "Options:\n",
          "  --label-projected                  Add text labels to projected points\n",
          "  --ellipses                         Draw ellipses for each (group × ref_taxon) subgroup\n",
          "  --exclude-ellipse=<csv>            Exclude: exact <group> (col2), exact '<group> × <taxa>', or substring in taxa (col3)\n",
          "  --marker[=<color>]                 Color ALL ellipses with one color (default black)\n",
          "  --ellipse-style=<solid|dashed|dotted>\n",
          "  --ellipse-lwd=<num>                Ellipse line width (aliases: --ellipse-width, --ellipse-size)\n",
          "  --point-size=<num>                 Point size for both layers (default 3.0)\n",
          "  --point-stroke=<num>               Point stroke width (default 0.9)\n",
          "  --title=<text>                     Set plot title (or: --title <text>)\n",
          "  --title-align=<left|center|right|0..1>\n",
          "  --font=<family> | --font-family=<family>\n",
          "  --title-size <pt>\n",
          "  --axis-title-size <pt> --axis-text-size <pt>\n",
          "  --legend-title-size <pt> --legend-text-size <pt>\n",
          "  --label-size <pt>\n",
          "  --width=<in> --height=<in>         Output size in inches (defaults 7.5 × 5.8)\n",
          "  --dpi=<n>                          Output DPI (default 300)\n")
  quit(save = "no", status = 2)
}

covar_dir  <- args[1]
fam_file   <- args[2]
group_tbl  <- args[3]
out_prefix <- args[4]

## -------------------- Parse optional flags --------------------
do_label <- FALSE
draw_ellipses <- FALSE
exclude_ellipse_patterns <- character(0)
ellipse_force_color <- NULL
ellipse_linetype <- "solid"
ellipse_line_size <- 1.0

plot_title <- "Procrustes PCA"
title_align_hjust <- 0.5
font_family <- NULL

out_width  <- 7.5
out_height <- 5.8
out_dpi    <- 300

point_size   <- 3.0
point_stroke <- 0.9

title_size <- NA_real_
axis_title_size <- NA_real_
axis_text_size  <- NA_real_
legend_title_size <- NA_real_
legend_text_size  <- NA_real_
label_text_size   <- NA_real_  # default later to 3.0

parse_exclude_value <- function(val) {
  if (is.null(val) || length(val) == 0) return(character(0))
  toks <- unlist(strsplit(val, ","))
  toks <- trimws(toks)
  toks[toks != ""]
}
parse_numeric_arg <- function(a, next_val) {
  if (grepl("=", a)) as.numeric(sub("^[^=]*=", "", a))
  else if (!is.null(next_val) && !startsWith(next_val, "--")) as.numeric(next_val)
  else NA_real_
}
parse_string_arg <- function(a, next_val) {
  if (grepl("=", a)) sub("^[^=]*=", "", a)
  else if (!is.null(next_val) && !startsWith(next_val, "--")) next_val
  else NA_character_
}
map_linetype <- function(x) {
  v <- tolower(x)
  if (v %in% c("solid","s")) return("solid")
  if (v %in% c("dashed","dash")) return("dashed")
  if (v %in% c("dotted","dot")) return("dotted")
  message("Warning: unknown --ellipse-style '", x, "'; using 'solid'.")
  "solid"
}

i <- 5
while (i <= length(args)) {
  a <- args[i]; a_next <- if (i + 1 <= length(args)) args[i+1] else NULL
  al <- tolower(gsub("^--", "", a))

  if (al %in% c("label-projected", "label_projected", "label")) {
    do_label <- TRUE

  } else if (al %in% c("ellipses", "ellipse", "draw-ellipses", "draw_ellipses")) {
    draw_ellipses <- TRUE

  } else if (grepl("^exclude(-ellipse)?=", al)) {
    exclude_ellipse_patterns <- c(exclude_ellipse_patterns, parse_exclude_value(sub("^[^=]*=", "", a)))

  } else if (al %in% c("exclude", "exclude-ellipse", "exclude_ellipse", "excludeellipse")) {
    if (!is.null(a_next) && !startsWith(a_next, "--")) { exclude_ellipse_patterns <- c(exclude_ellipse_patterns, parse_exclude_value(a_next)); i <- i + 1 }
    else message("Warning: ", a, " provided without a value; ignoring.")

  } else if (grepl("^marker(=|$)", al)) {
    if (grepl("^--?marker=", a, ignore.case = TRUE)) {
      ellipse_force_color <- sub("^--?marker=", "", a, ignore.case = TRUE); if (ellipse_force_color == "") ellipse_force_color <- "black"
    } else ellipse_force_color <- "black"

  } else if (grepl("^(ellipse-style|ellipse-linetype|ellipse_type)(=|$)", al)) {
    val <- parse_string_arg(a, a_next); if (!is.na(val)) { ellipse_linetype <- map_linetype(val); if (!grepl("=", a) && !startsWith(a_next, "--")) i <- i + 1 }

  } else if (grepl("^(ellipse-lwd|ellipse-width|ellipse-size)(=|$)", al)) {
    val <- parse_numeric_arg(a, a_next); if (!is.na(val)) { ellipse_line_size <- val; if (!grepl("=", a) && !startsWith(a_next, "--")) i <- i + 1 }

  } else if (grepl("^point-size(=|$)", al)) {
    val <- parse_numeric_arg(a, a_next); if (!is.na(val)) { point_size <- val; if (!grepl("=", a) && !startsWith(a_next, "--")) i <- i + 1 }

  } else if (grepl("^point-stroke(=|$)", al)) {
    val <- parse_numeric_arg(a, a_next); if (!is.na(val)) { point_stroke <- val; if (!grepl("=", a) && !startsWith(a_next, "--")) i <- i + 1 }

  } else if (grepl("^title(=|$)", al)) {
    val <- parse_string_arg(a, a_next); if (!is.na(val)) { plot_title <- val; if (!grepl("=", a) && !startsWith(a_next, "--")) i <- i + 1 }

  } else if (grepl("^(title-align|title_pos|title-position|titlepos)(=|$)", al)) {
    val <- parse_string_arg(a, a_next)
    if (!is.na(val)) {
      v <- tolower(val)
      if (v %in% c("left","centre","center","right")) title_align_hjust <- if (v == "left") 0 else if (v == "right") 1 else 0.5
      else { num <- suppressWarnings(as.numeric(v)); if (is.finite(num)) title_align_hjust <- max(0, min(1, num)) }
      if (!grepl("=", a) && !startsWith(a_next, "--")) i <- i + 1
    }

  } else if (grepl("^(font|font-family)(=|$)", al)) {
    val <- parse_string_arg(a, a_next); if (!is.na(val)) { font_family <- val; if (!grepl("=", a) && !startsWith(a_next, "--")) i <- i + 1 }

  } else if (grepl("^title-size(=|$)", al)) {
    val <- parse_numeric_arg(a, a_next); if (!is.na(val)) { title_size <- val; if (!grepl("=", a) && !startsWith(a_next, "--")) i <- i + 1 }
  } else if (grepl("^axis-title-size(=|$)", al)) {
    val <- parse_numeric_arg(a, a_next); if (!is.na(val)) { axis_title_size <- val; if (!grepl("=", a) && !startsWith(a_next, "--")) i <- i + 1 }
  } else if (grepl("^axis-text-size(=|$)", al)) {
    val <- parse_numeric_arg(a, a_next); if (!is.na(val)) { axis_text_size <- val; if (!grepl("=", a) && !startsWith(a_next, "--")) i <- i + 1 }
  } else if (grepl("^legend-title-size(=|$)", al)) {
    val <- parse_numeric_arg(a, a_next); if (!is.na(val)) { legend_title_size <- val; if (!grepl("=", a) && !startsWith(a_next, "--")) i <- i + 1 }
  } else if (grepl("^legend-text-size(=|$)", al)) {
    val <- parse_numeric_arg(a, a_next); if (!is.na(val)) { legend_text_size <- val; if (!grepl("=", a) && !startsWith(a_next, "--")) i <- i + 1 }
  } else if (grepl("^label-size(=|$)", al)) {
    val <- parse_numeric_arg(a, a_next); if (!is.na(val)) { label_text_size <- val; if (!grepl("=", a) && !startsWith(a_next, "--")) i <- i + 1 }

  } else if (grepl("^width(=|$)", al)) {
    val <- parse_numeric_arg(a, a_next); if (!is.na(val)) { out_width <- val; if (!grepl("=", a) && !startsWith(a_next, "--")) i <- i + 1 }
  } else if (grepl("^height(=|$)", al)) {
    val <- parse_numeric_arg(a, a_next); if (!is.na(val)) { out_height <- val; if (!grepl("=", a) && !startsWith(a_next, "--")) i <- i + 1 }
  } else if (grepl("^dpi(=|$)", al)) {
    val <- parse_numeric_arg(a, a_next); if (!is.na(val)) { out_dpi <- as.integer(val); if (!grepl("=", a) && !startsWith(a_next, "--")) i <- i + 1 }

  } else {
    message("Note: Unrecognized option '", a, "' (ignored).")
  }
  i <- i + 1
}
exclude_ellipse_patterns <- unique(tolower(exclude_ellipse_patterns))

## Title wrapping to fit figure width (approx chars-per-inch ~= 10)
wrap_text <- function(s, width_chars = 60) paste(strwrap(s, width = width_chars), collapse = "\n")
wrap_chars <- max(24, floor(out_width * 10))
plot_title_wrapped <- wrap_text(plot_title, width_chars = wrap_chars)

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
  if (is.null(getOption("repos")) || is.na(getOption("repos")["CRAN"]) ||
      getOption("repos")["CRAN"] %in% c("@CRAN@", NULL, "")) {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
  }
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) {
    message("Installing missing package: ", p)
    install.packages(p, lib = userlib, dependencies = TRUE)
  }
}
base_pkgs <- c("vegan", "ggplot2", "paletteer")
extra_pkgs <- if (do_label) c(base_pkgs, "ggrepel") else base_pkgs
ensure_pkgs(extra_pkgs)

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

## -------------------- Sample -> group (3 columns) --------------------
grp <- read.table(group_tbl, header = FALSE, stringsAsFactors = FALSE)
if (ncol(grp) < 3) stop("Group table must have 3 columns: <sample_prefix> <group_label> <ref_taxon_for_color>")
colnames(grp)[1:3] <- c("sample_prefix", "group_label", "ref_taxon")

file_prefixes <- sub("_covar.txt$", "", basename(covar_files_paths))
files_df <- data.frame(filesName = file_prefixes, files = covar_files_paths, stringsAsFactors = FALSE)

keep <- intersect(grp$sample_prefix, files_df$filesName)
if (length(keep) == 0) stop("No overlap between group table and *_covar.txt prefixes.")
grp <- grp[grp$sample_prefix %in% keep, , drop = FALSE]
files_df <- files_df[files_df$filesName %in% keep, , drop = FALSE]
grp <- merge(files_df, grp, by.x = "filesName", by.y = "sample_prefix", sort = FALSE)

## -------------------- Compute per-sample PC1/PC2 --------------------
tmpdir <- "procrustesPCs"
message("Creating temporary directory: ", tmpdir)
mk_status <- if (dir.exists(tmpdir)) 1L else { dir.create(tmpdir, recursive = TRUE); 0L }

pc_to_file <- function(f) {
  r  <- read.table(f, header = TRUE, check.names = FALSE)
  ee <- eigen(as.matrix(r))
  outname <- sub("_covar.txt$", "", basename(f))
  write.table(ee$vectors[, 1:2],
              file = file.path(tmpdir, outname),
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

pcs_paths <- file.path(tmpdir, grp$filesName)
if (any(!file.exists(pcs_paths))) stop("Missing PC files in ", tmpdir, "/; check permissions and inputs.")

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
  Species  = fam$V1,
  Sample   = fam$V2,
  PC1      = m0[, 1],
  PC2      = m0[, 2],
  Group    = "Reference",
  ColorPop = fam$V1,
  stringsAsFactors = FALSE
)

proj_df <- data.frame(
  Sample   = grp$filesName,
  Group    = grp$group_label,
  ColorPop = grp$ref_taxon,
  PC1      = res_mat[, 1],
  PC2      = res_mat[, 2],
  stringsAsFactors = FALSE
)

df <- rbind(
  ref_df[,  c("Sample","PC1","PC2","Group","ColorPop")],
  proj_df[, c("Sample","PC1","PC2","Group","ColorPop")]
)

## -------------------- Subgroups & exclusions for ellipses --------------------
proj_only <- subset(df, Group != "Reference")
proj_only$GroupKey <- paste0(proj_only$Group, " × ", proj_only$ColorPop)

excluded_keys <- character(0)
if (length(exclude_ellipse_patterns)) {
  gl <- tolower(proj_only$Group)
  rk <- tolower(proj_only$GroupKey)
  rt <- tolower(proj_only$ColorPop)
  mask <- rep(FALSE, nrow(proj_only))
  for (pat in exclude_ellipse_patterns) {
    if (pat == "") next
    mask <- mask |
      (gl == pat) |                 # exact match on group (col2)
      (rk == pat) |                 # exact match on combined key
      grepl(pat, rt, fixed = TRUE)  # substring within taxa (col3)
  }
  excluded_keys <- unique(proj_only$GroupKey[mask])
  if (length(excluded_keys))
    message("Excluding ellipses for: ", paste(excluded_keys, collapse = "; "))
}

## -------------------- Ellipses (robust, less line-like) --------------------
ellipse_level <- 0.95
ellipse_min_axis_frac <- 0.05   # >= 5% of min axis span
global_span_x <- diff(range(c(df$PC1), na.rm = TRUE))
global_span_y <- diff(range(c(df$PC2), na.rm = TRUE))
min_axis <- ellipse_min_axis_frac * min(global_span_x, global_span_y)
default_r <- 0.02 * min(global_span_x, global_span_y)

make_cov_ellipse <- function(x, y, level = 0.95, n = 360, min_axis = 0.01) {
  cx <- mean(x); cy <- mean(y)
  Sigma <- stats::cov(cbind(x, y))
  if (!all(is.finite(Sigma))) return(NULL)
  eig <- eigen(Sigma, symmetric = TRUE)
  lam <- pmax(eig$values, .Machine$double.eps)   # stability
  chi <- sqrt(stats::qchisq(level, df = 2))
  axes <- pmax(sqrt(lam) * chi, min_axis)        # minimum axis lengths
  R <- eig$vectors %*% diag(axes, nrow = 2)
  th <- seq(0, 2*pi, length.out = n)
  pts <- t(R %*% rbind(cos(th), sin(th)))
  data.frame(x = pts[,1] + cx, y = pts[,2] + cy)
}

ellipse_df <- data.frame(x = numeric(0), y = numeric(0), GroupKey = character(0), EllipseColor = character(0))
if (draw_ellipses && nrow(proj_only)) {
  ellipse_list <- list()
  for (key in sort(unique(proj_only$GroupKey))) {
    if (key %in% excluded_keys) next
    sub <- proj_only[proj_only$GroupKey == key, , drop = FALSE]
    ell <- NULL
    if (nrow(sub) >= 2) {
      ell <- try(make_cov_ellipse(sub$PC1, sub$PC2, level = ellipse_level, n = 360, min_axis = min_axis), silent = TRUE)
      if (inherits(ell, "try-error") || is.null(ell)) ell <- NULL
    }
    if (is.null(ell)) {
      ## Fallback small circle if single point/degenerate
      cx <- mean(sub$PC1); cy <- mean(sub$PC2)
      th <- seq(0, 2*pi, length.out = 360)
      ell <- data.frame(x = cx + default_r * cos(th), y = cy + default_r * sin(th))
    }
    ucol <- unique(sub$ColorPop)
    ellipse_col <- if (!is.null(ellipse_force_color)) ellipse_force_color else if (length(ucol) == 1) ucol[1] else "MixedGroup"
    ell$GroupKey <- key
    ell$EllipseColor <- ellipse_col
    ellipse_list[[key]] <- ell
  }
  if (length(ellipse_list)) ellipse_df <- do.call(rbind, ellipse_list)
}

## -------------------- Shapes (uniform sizes) --------------------
proj_groups <- setdiff(unique(df$Group), "Reference")
shape_seq <- c(17, 18, 15, 8, 4)  # triangle, diamond, square, star, cross
shape_map <- c("Reference" = 16)
if (length(proj_groups) > 0) {
  shape_map[proj_groups] <- shape_seq[((seq_along(proj_groups) - 1) %% length(shape_seq)) + 1]
}

## -------------------- Colors palette --------------------
color_levels <- unique(c(ref_df$ColorPop, proj_df$ColorPop, if (nrow(ellipse_df)) as.character(ellipse_df$EllipseColor)))
n_cols <- length(color_levels)
if (n_cols <= 10) {
  pal_vec <- as.vector(paletteer::paletteer_d("ggthemes::Classic_10", n = n_cols))
} else {
  base_pal <- as.vector(paletteer::paletteer_d("ggthemes::Classic_20", n = min(20, n_cols)))
  pal_vec  <- base_pal[((seq_len(n_cols) - 1) %% length(base_pal)) + 1]
}
names(pal_vec) <- color_levels
if ("MixedGroup" %in% names(pal_vec)) pal_vec["MixedGroup"] <- "#4D4D4D"

df$ColorPop <- factor(df$ColorPop, levels = color_levels)
if (nrow(ellipse_df)) ellipse_df$EllipseColor <- factor(ellipse_df$EllipseColor, levels = color_levels)

## -------------------- Save data tables --------------------
write.table(df[, c("Sample", "PC1", "PC2", "Group", "ColorPop")],
            paste0(out_prefix, ".txt"),
            col.names = TRUE, row.names = FALSE, quote = FALSE)

legacy_points <- rbind(
  data.frame(ColorPop = ref_df$ColorPop,
             Sample   = ref_df$Sample,
             PC1      = ref_df$PC1,
             PC2      = ref_df$PC2,
             stringsAsFactors = FALSE),
  data.frame(ColorPop = proj_df$ColorPop,
             Sample   = proj_df$Sample,
             PC1      = proj_df$PC1,
             PC2      = proj_df$PC2,
             stringsAsFactors = FALSE)
)
write.table(legacy_points, paste0(out_prefix, ".legacy.txt"),
            col.names = FALSE, row.names = FALSE, quote = FALSE)

## -------------------- Plot --------------------
set.seed(1)
base_family <- ifelse(is.null(font_family), "", font_family)

p <- ggplot()

if (nrow(ellipse_df) > 0) {
  if (!is.null(ellipse_force_color)) {
    p <- p + geom_path(
      data = ellipse_df,
      aes(x = x, y = y, group = GroupKey),
      color = ellipse_force_color, size = ellipse_line_size, alpha = 0.95,
      linetype = ellipse_linetype
    )
  } else {
    p <- p + geom_path(
      data = ellipse_df,
      aes(x = x, y = y, group = GroupKey, color = EllipseColor),
      size = ellipse_line_size, alpha = 0.95,
      linetype = ellipse_linetype
    )
  }
}

p <- p +
  geom_point(
    data = subset(df, Group == "Reference"),
    aes(x = PC1, y = PC2, color = ColorPop, shape = Group),
    size = point_size, stroke = point_stroke, alpha = 0.95
  ) +
  geom_point(
    data = subset(df, Group != "Reference"),
    aes(x = PC1, y = PC2, color = ColorPop, shape = Group),
    size = point_size, stroke = point_stroke, alpha = 0.98
  ) +
  scale_colour_manual(values = pal_vec, name = "Taxa") +
  scale_shape_manual(values = shape_map,  name = "Group") +
  guides(
    shape  = guide_legend(override.aes = list(size = point_size, stroke = point_stroke, alpha = 1)),
    colour = guide_legend(override.aes = list(size = point_size, alpha = 1))
  ) +
  coord_equal(xlim = xlim, ylim = ylim, expand = TRUE) +
  xlab(paste0("PC1 (", sprintf("%.2f", explained_var[1]), "%)")) +
  ylab(paste0("PC2 (", sprintf("%.2f", explained_var[2]), "%)")) +
  ggtitle(plot_title_wrapped) +
  theme_minimal(base_size = 12, base_family = base_family)

## ---- Apply detailed theme (font, sizes, title alignment) ----
th <- theme(
  legend.position = "right",
  panel.grid.minor = element_blank(),
  plot.title = element_text(face = "bold", hjust = title_align_hjust)
)
if (!is.null(font_family)) th <- th + theme(text = element_text(family = font_family))
if (is.finite(title_size))        th <- th + theme(plot.title   = element_text(face = "bold", hjust = title_align_hjust, size = title_size, family = font_family))
if (is.finite(axis_title_size))   th <- th + theme(axis.title   = element_text(size = axis_title_size))
if (is.finite(axis_text_size))    th <- th + theme(axis.text    = element_text(size = axis_text_size))
if (is.finite(legend_title_size)) th <- th + theme(legend.title = element_text(size = legend_title_size))
if (is.finite(legend_text_size))  th <- th + theme(legend.text  = element_text(size = legend_text_size))
p <- p + th

## ---- Optional labels ----
label_size_final <- if (is.finite(label_text_size)) label_text_size else 3.0
if (do_label) {
  p <- p +
    ggrepel::geom_text_repel(
      data = subset(df, Group != "Reference"),
      aes(x = PC1, y = PC2, label = Sample, color = ColorPop),
      size = label_size_final, max.overlaps = 50,
      box.padding = 0.35, point.padding = 0.25,
      min.segment.length = 0,
      segment.alpha = 0.6, segment.size = 0.3,
      seed = 1,
      show.legend = FALSE
    )
}

ggsave(filename = paste0(out_prefix, ".png"),
       plot = p, width = out_width, height = out_height, dpi = out_dpi, bg = "white")

## -------------------- Cleanup --------------------
if (mk_status == 0) unlink(tmpdir, recursive = TRUE, force = TRUE)
