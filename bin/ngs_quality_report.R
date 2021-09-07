#!/usr/bin/env Rscript
#packages <- c("argparse", "here", "BiocManager", "rmarkdown", "knitr", "ggplot2")
#if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
#  install.packages(setdiff(packages, rownames(installed.packages())), repos="https://cran.wu.ac.at/")  
#}

if(!"dada2" %in% rownames(installed.packages())) {
  library(BiocManager)
  BiocManager::install("dada2")
}

library(argparse)
library(here)
library(ggplot2)
library(dada2)

parser <- argparse::ArgumentParser(description="")
parser$add_argument('--title', type="character", default="FASTQ Sequence Quality Plots")
parser$add_argument('--fwd_reads', nargs='+', required=TRUE, help="relative path to forward reads")
parser$add_argument('--rev_reads', nargs='+', required=FALSE, help="realtive path to reverse reads")
parser$add_argument('--out_html', type="character", required=FALSE, default="./fastq_quality_report.html")
parser$add_argument('--out_dir', type="character", required=FALSE, default="./fastq_quality/")
parser$add_argument('-f', action="store_true", default=FALSE, help="Overwrite existing output file.")

#args <- parser$parse_args(c("--fwd_reads", "input_data/R0_S1_L001_R1_001.fastq", "input_data/R2_S2_L001_R1_001.fastq", 
#                            "--rev_reads", "input_data/R0_S1_L001_R2_001.fastq", "input_data/R2_S2_L001_R2_001.fastq",
#                            "-f", "--out_html", "fastq_quality.html", "--out_dir", "qualls"))
args <- parser$parse_args() 

# Check if user provided equal length forward and reverse read list or if reverse list is empty
if (length(parser$fwd_reads) != length(parser$rev_reads) && length(c(parser$rev_reads)) > 0) {
  stop("Lists of Forward and Reverse FASTQ files have to be of equal length.")
}

# Check if output files exist
if (file.exists(args$out_html)) {
  if (args$f) {
    file.remove(args$out_html)
  } else {
    stop(paste0("Output file ", args$out_html, " already exists. Add -f to force overwrite."))
  }
}
if (file.exists(args$out_dir)) {
  if (args$f) {
    unlink(args$out_dir, recursive = TRUE)
  } else {
    stop(paste0("Output directory ", args$out_dir, " already exists. Add -f to force overwrite."))
  }
}
dir.create(args$out_dir)

# Create forwad and reverse plot for all files provided together
gg_fwd_general <- dada2::plotQualityProfile(args$fwd_reads, aggregate = TRUE)
gg_rev_general <- NULL
ggsave("ngs_quality_forward.png", path=args$out_dir, plot=gg_fwd_general)
if (!is.null(args$rev_reads) && length(args$rev_reads) > 0) {
  gg_rev_general <- plotQualityProfile(args$rev_reads, aggregate = TRUE)
  ggsave("ngs_quality_reverse.png", plot=gg_rev_general, path=args$out_dir)
}

# Create a plot for every forward read
gg_fwd <- list()
for (i in 1:length(args$fwd_reads)) {
  fwd = args$fwd_reads[i]
  gg_fwd[[fwd]] <- plotQualityProfile(fwd)
  ggsave(paste0(basename(fwd), ".ngs_quality.png"), gg_fwd[[fwd]], path=args$out_dir)
}

# Create a plot for every reverse read
gg_rev <- list()
if (length(args$rev_reads) > 0) {
  for (i in 1:length(args$rev_reads)) {
    rev = args$rev_reads[i]
    gg_rev[[rev]] <- plotQualityProfile(rev)
    ggsave(paste0(basename(rev), ".ngs_quality.png"), gg_rev[[rev]], path=args$out_dir)
  }
}

# Prepare params list for RMarkdown
gg_params <- list()
gg_params$title = args$title
gg_params$date = args$date
gg_params$gg_fwd_general <- gg_fwd_general
gg_params$gg_rev_general <- gg_rev_general
gg_params$gg_fwd = gg_fwd
gg_params$gg_rev = gg_rev

# Start rendering
rmarkdown::render(
  here('bin', 'ngs_quality_report_plot.Rmd'), 
  params=gg_params, 
  output_file=args$out_html,
  output_dir="."
)

