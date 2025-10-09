#!/usr/bin/env Rscript

## Final: visible kinship + contrasting inbreeding palette, top-half only, per-group
## Now with --font, --text-size, and combined panel output via gridExtra

quiet_install_and_load <- function(pkgs){
  ip <- rownames(installed.packages())
  miss <- pkgs[!pkgs %in% ip]
  if (length(miss))
    suppressMessages(suppressWarnings(
      install.packages(miss, repos="https://cloud.r-project.org", quiet=TRUE)
    ))
  invisible(lapply(pkgs, function(p)
    suppressPackageStartupMessages(suppressMessages(library(p, character.only=TRUE)))))
}
quiet_install_and_load(
  c("optparse","readr","dplyr","tidyr","stringr","ggplot2","ggnewscale","tools","scales","RColorBrewer","gridExtra")
)

suppressPackageStartupMessages({
  library(optparse); library(readr); library(dplyr); library(tidyr)
  library(stringr); library(ggplot2); library(ggnewscale); library(tools); library(scales)
  library(RColorBrewer); library(gridExtra)
})

opt <- OptionParser(option_list=list(
  make_option(c("-r","--ngsrelate"), type="character", help="ngsRelate file"),
  make_option(c("-m","--metadata"), type="character", help="Metadata TSV with Sample,Group"),
  make_option(c("-o","--outdir"), type="character", default="kinship_by_group", help="Output directory"),
  make_option(c("-n","--named"), action="store_true", default=FALSE, help="Use named-column layout"),
  make_option(c("-t","--metric"), type="character", default="king", help="Kinship metric: king|rab|theta"),
  make_option(c("-O","--order"), type="character", default="metadata", help="Order within group: metadata|cluster"),
  make_option(c("-w","--width"),  type="double", default=5, help="Per-plot width (inches)"),
  make_option(c("-H","--height"), type="double", default=5, help="Per-plot height (inches)"),
  make_option(c("-d","--dpi"), type="integer", default=300, help="DPI for output"),

  # NEW: typography and combined panel controls
  make_option(c("-f","--font"), type="character", default="sans", help="Font family for all text"),
  make_option(c("-s","--text-size"), type="double", default=10, help="Base text size"),
  make_option(c("--combine"), action="store_true", default=TRUE, help="Also write a combined panel image [default]"),
  make_option(c("--no-combine"), action="store_false", dest="combine", help="Do not write the combined panel"),
  make_option(c("--ncol"), type="integer", default=2, help="Columns in combined panel"),
  make_option(c("--nrow"), type="integer", default=NA, help="Rows in combined panel (auto if NA)"),
  make_option(c("--outfile"), type="character", default="combined.png", help="Filename for combined panel (inside outdir)")
))
args <- parse_args(opt)

# safety helpers
sanitize_font <- function(x) if (is.null(x) || length(x)==0 || is.na(x) || x=="") "sans" else as.character(x)
sanitize_num  <- function(x, default) { sx <- suppressWarnings(as.numeric(x)); if (length(sx)==0 || is.na(sx) || sx <= 0) default else sx }
args$font <- sanitize_font(args$font)
args$text_size <- sanitize_num(args$text_size, 10)

stop_if <- function(cond,msg){ if (cond) stop(msg, call.=FALSE) }
stop_if(is.null(args$ngsrelate), "Provide --ngsrelate (-r)")
stop_if(is.null(args$metadata),  "Provide --metadata (-m)")

meta <- suppressMessages(readr::read_tsv(args$metadata, col_types=cols(.default=col_character())))
stop_if(!all(c("Sample","Group") %in% names(meta)), "Metadata must have Sample,Group columns")
meta <- meta %>% group_by(Group) %>% filter(n() >= 2) %>% ungroup()

dir.create(args$outdir, showWarnings=FALSE, recursive=TRUE)

colmap <- function(named){
  if (named) list(ida=3, idb=4, king=33, rab=15, theta=17)
  else       list(ida=1, idb=2, king=31, rab=13, theta=15)
}
metric_col <- function(m, cm){
  m <- tolower(m)
  if (m=="king")  return(cm$king)
  if (m=="rab")   return(cm$rab)
  if (m=="theta") return(cm$theta)
  stop("Unknown --metric. Use king|rab|theta")
}
guess_inb_cols <- function(named) if (named) c(21,22) else c(19,20)

read_ngs_file <- function(path, named, metric){
  df <- suppressWarnings(read.table(path, header=FALSE, stringsAsFactors=FALSE, check.names=FALSE))
  cm <- colmap(named); kin_col <- metric_col(metric, cm)
  inb_cols <- guess_inb_cols(named)

  kin_raw <- as.character(df[[kin_col]])
  kin_num <- suppressWarnings(as.numeric(kin_raw))
  kin_num[is.na(kin_num) & kin_raw == "-"] <- 0   # map "-" to zero

  Fa <- Fb <- NULL
  if (ncol(df) >= max(inb_cols)) {
    Fa <- suppressWarnings(as.numeric(df[[inb_cols[1]]]))
    Fb <- suppressWarnings(as.numeric(df[[inb_cols[2]]]))
  }

  tibble(
    ida = as.character(df[[cm$ida]]),
    idb = as.character(df[[cm$idb]]),
    kin = kin_num,
    Fa  = if (!is.null(Fa)) Fa else NA_real_,
    Fb  = if (!is.null(Fb)) Fb else NA_real_
  )
}

order_group_samples <- function(ids, method, kin_tbl, meta){
  ids <- unique(ids)
  if (method=="metadata") {
    return(meta %>% filter(Sample %in% ids) %>%
             mutate(idx=row_number()) %>% arrange(idx) %>% pull(Sample))
  } else if (method=="cluster") {
    sub <- kin_tbl %>% filter(ida %in% ids & idb %in% ids)
    M <- matrix(0, length(ids), length(ids), dimnames=list(ids,ids))
    if (nrow(sub)>0) {
      for (i in seq_len(nrow(sub))) {
        a <- sub$ida[i]; b <- sub$idb[i]; v <- sub$kin[i]
        M[a,b] <- v; M[b,a] <- v
      }
      hc <- hclust(as.dist(1 - M / max(1e-12, max(M, na.rm=TRUE))))
      return(ids[hc$order])
    } else return(sort(ids))
  } else sort(ids)
}

pairs_all <- read_ngs_file(args$ngsrelate, args$named, args$metric)
pairs_all <- pairs_all %>% filter(ida %in% meta$Sample & idb %in% meta$Sample)
meta2 <- meta %>% distinct(Sample, Group)
pairs_all <- pairs_all %>%
  inner_join(meta2, by=c("ida"="Sample")) %>% rename(GroupA=Group) %>%
  inner_join(meta2, by=c("idb"="Sample")) %>% rename(GroupB=Group) %>%
  filter(GroupA == GroupB) %>% rename(Group=GroupA) %>%
  select(Group, ida, idb, kin, Fa, Fb)
stop_if(nrow(pairs_all)==0, "No within-group pairs to plot.")

inb <- bind_rows(
  pairs_all %>% transmute(Sample=ida, inb=Fa),
  pairs_all %>% transmute(Sample=idb, inb=Fb)
) %>% group_by(Sample) %>%
  summarize(inb = suppressWarnings(median(inb, na.rm=TRUE)), .groups="drop")

plots <- list()

for (g in sort(unique(pairs_all$Group))) {
  sub <- pairs_all %>% filter(Group==g)
  ids <- unique(c(sub$ida, sub$idb))
  ids <- order_group_samples(ids, args$order, sub, meta)

  diag_df <- tibble(Sample=ids) %>%
    left_join(inb, by="Sample") %>%
    mutate(
      i=factor(Sample, levels=ids),
      j=factor(Sample, levels=ids),
      type="inbreeding", val=inb
    )

  a1 <- sub %>% select(ida,idb,kin)
  a2 <- sub %>% transmute(ida=idb, idb=ida, kin=kin)
  both <- bind_rows(a1,a2)

  # top half including diagonal
  top_half <- expand.grid(i=ids, j=ids, stringsAsFactors=FALSE) %>%
    as_tibble() %>% filter(match(i,ids) <= match(j,ids)) %>%
    left_join(both, by=c("i"="ida","j"="idb"))

  df_plot <- bind_rows(
    diag_df,
    top_half %>% transmute(
      i=factor(i,levels=ids),
      j=factor(j,levels=rev(ids)),
      type="kinship",
      val=kin
    )
  )

  p <- ggplot() +
    # Kinship scale (cool blues)
    geom_tile(data=df_plot %>% filter(type=="kinship"),
              aes(x=i, y=j, fill=val), color=NA) +
    scale_fill_gradient2(
      name="Kinship\ncoefficient",
      limits=c(0,0.5),
      oob=scales::squish,
      low="#f7f7f7", mid="#9ecae1", high="#08306b",
      midpoint=0.25,
      guide=guide_colorbar(reverse=FALSE)
    ) +
    ggnewscale::new_scale_fill() +
    # Inbreeding scale (continuous YlOrBr)
    geom_tile(data=df_plot %>% filter(type=="inbreeding"),
              aes(x=i, y=j, fill=val), color=NA) +
    scale_fill_gradientn(
      name = "Inbreeding\ncoefficient",
      limits = c(0, 0.5),
      oob = scales::squish,
      colours = RColorBrewer::brewer.pal(9, "YlOrBr"),
      na.value = "grey90",
      guide = guide_colorbar(reverse = FALSE)
    ) +
    labs(title=g, x=NULL, y=NULL) +
    coord_equal() +
    theme_minimal(base_size=args$text_size) +
    theme(
      text = element_text(family=args$font, size=args$text_size),
      panel.grid=element_blank(),
      axis.text.x=element_text(angle=45,hjust=1,vjust=1),
      plot.title=element_text(hjust=0.5)
    )

  # save individual
  outpath <- file.path(args$outdir, paste0(g,".png"))
  ggsave(outpath, p, width=args$width, height=args$height, dpi=args$dpi, bg="white")
  message("Saved: ", outpath)

  plots[[length(plots)+1]] <- p
}

# Combined panel (optional, default TRUE)
if (isTRUE(args$combine) && length(plots) > 0){
  # auto rows if NA
  if (is.na(args$nrow)) {
    args$nrow <- ceiling(length(plots) / max(1, args$ncol))
  }
  g <- do.call(gridExtra::arrangeGrob, c(plots, ncol=args$ncol, nrow=args$nrow))
  comb_path <- file.path(args$outdir, args$outfile)

  # set canvas size as layout * per-plot size
  width_total  <- args$width  * args$ncol
  height_total <- args$height * args$nrow

  ggsave(comb_path, g, width=width_total, height=height_total, dpi=args$dpi, bg="white")
  message("Saved combined panel: ", comb_path)
}

message("All done: outputs in ", args$outdir)
