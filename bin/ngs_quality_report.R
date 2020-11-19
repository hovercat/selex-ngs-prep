#!/usr/bin/env Rscript
packages <- c("argparse", "here", "BiocManager", "rmarkdown", "knitr")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

if(!"dada2" %in% rownames(installed.packages())) {
  BiocManager::install("dada2")
}

library(argparse)
library(here)

parser <- argparse::ArgumentParser(description="")
parser$add_argument('--title', type="character", default="FASTQ Sequence Quality Plots")
parser$add_argument('--author', type="character", default="")
parser$add_argument('--ngs_run_date', type="character", default=paste0(Sys.Date()))
parser$add_argument('--date', type="character", default=paste0(Sys.Date()))
parser$add_argument('--fwd_reads', nargs='+', required=TRUE, help="relative path to forward reads")
parser$add_argument('--rev_reads', nargs='+', required=FALSE, help="realtive path to reverse reads")

# temporary arguments (will be removed before rendering)
parser$add_argument('-o', type="character", required=FALSE, default="./fastq_quality_report.html")
parser$add_argument('-f', action="store_true", default=FALSE, help="Overwrite existing output file.")

#args <- parser$parse_args(c("--fwd_reads", "input_data/R0_S1_L001_R1_001.fastq", "input_data/R0_S1_L001_R1_001.fastq", 
#                            #"--rev_reads", "input_data/R0_S1_L001_R1_001.fastq", "input_data/R0_S1_L001_R1_001.fastq",
#                            "-f"))
args <- parser$parse_args()

if (length(parser$fwd_reads) != length(parser$rev_reads) && length(c(parser$rev_reads)) > 0) {
  stop("Lists of Forward and Reverse FASTQ files have to be of equal length.")
}

if (file.exists(args$o)) {
  if (args$f) {
    file.remove(args$o)
  } else {
    stop(paste0("Output file ", args$o, " already exists. Add -f to force overwrite."))
  }
}

output_file = args$o
args$o <- NULL
args$f <- NULL
args$wd = getwd()

rmarkdown::render(
  here('bin', 'ngs_quality_report.Rmd'), 
  params=args, 
  output_file=output_file,
  output_dir="."
)

