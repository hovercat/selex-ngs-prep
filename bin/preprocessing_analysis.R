#!/usr/bin/env Rscript
library(argparse)
library(tidyr)
library(dplyr)

parser <- argparse::ArgumentParser(description="")
parser$add_argument('-i', type="character", help="Input CSV with preprocessing losses")
parser$add_argument('-o', type="character", help="Plot name", required=FALSE)
parser$add_argument('-p', type="character", help="Percentage plot name ", required=FALSE)
parser$add_argument('-f', action="store_true", help="Overwrite existing output file?")

#args <- parser$parse_args(c("-i", "output/preprocessing_analysis/preprocessing_analysis_fwd.csv"))
args <- parser$parse_args()

if (is.null(args$o) || args$o == "") {
  args$o = paste0(args$i, ".png")
}
if (is.null(args$p) || args$p == "") {
  args$p = paste0(args$i, ".perc.png")
}

if (file.exists(args$o)) {
  if (args$f) {
    file.remove(args$o)
  } else {
    stop(paste0("Output file ", args$o, " already exists. Add -f to force overwrite."))
  }
}

prep <- read.table(args$i, header = TRUE, sep = '\t') %>%
  mutate(
    `-> Passed` = restrict_size,
    `4. Size Filter` = merge - restrict_size,
    `3. Merging Loss` = filter - merge,
    `2. Quality Filter` = trim - filter,
    `1. No Primers Found` = raw - trim,
    
    raw = NULL,
    trim = NULL,
    filter = NULL,
    merge = NULL,
    trim = NULL,
    restrict_size = NULL
    
  )
col_order = colnames(prep)[-1:-2]
prep_long <- prep %>%
  pivot_longer(c(-round_id, -round_name), names_to="Discarded", values_to="read_count")

prep_long$round_name <- factor(prep_long$round_name, levels = levels(factor(prep_long$round_name))[order(unique(prep_long$round_id), decreasing=TRUE)])
prep_long$Discarded <- factor(prep_long$Discarded, levels = col_order)

library(ggplot2)
library(viridis)

gg <- ggplot(prep_long, aes(x=round_name, y=read_count)) +
  geom_col(aes(fill=Discarded)) +
  coord_flip() +
  scale_fill_viridis_d(begin=0.35, end=0.9) +
  xlab("Round Name") +
  ylab("Filter Performance")


gg_perc <- ggplot(prep_long, aes(x=round_name, y=read_count)) +
  geom_col(position="fill", aes(fill=Discarded)) +
  coord_flip() +
  scale_fill_viridis_d(begin=0.35, end=0.9) +
  xlab("Round Name") +
  ylab("Filter Performance Percentage") +
  scale_y_continuous(labels = scales::percent)

ggsave(args$o, gg)
ggsave(args$p, gg_perc)

