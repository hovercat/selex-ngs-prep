#!/usr/bin/env Rscript
#packages <- c("argparse", "here", "stringi", "stringr", "knitr", "rmarkdown", "dplyr", "tidyr", "ggplot2")
#if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
#  install.packages(setdiff(packages, rownames(installed.packages())), repos="https://cran.wu.ac.at/")  
#}

library(argparse)
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

parser <- argparse::ArgumentParser(description="")
parser$add_argument('--title', type="character", default="")
parser$add_argument('--date', type="character", default=paste0(Sys.Date()))
parser$add_argument('-i', "--input", nargs='+', required=FALSE, help="realtive path to csv files with nt data")
parser$add_argument('-o', type="character", required=FALSE, default="./nt_composition.html")
parser$add_argument('-f', action="store_true", default=FALSE, help="Overwrite existing output file.")

args <- parser$parse_args()

if (file.exists(args$o)) {
  if (args$f) {
    file.remove(args$o)
  } else {
    stop(paste0("Output file ", args$o, " already exists. Add -f to force overwrite."))
  }
}

# Read all CSV files and collect them into one data frame (adding round_name and sorting index along the way)
collected_csv <- NULL
for (j in 1:length(args$input)) {
  file_name = basename(args$input[j])
  round_name = substr(file_name, 0, stringr::str_locate(file_name, "\\.")[1]-1)
  round_csv <- read.table(args$input[j], header=TRUE, sep='\t')
  round_csv <- round_csv %>%
    mutate(
      round_index = j,
      round_name = round_name
    )
  collected_csv <- rbind(collected_csv, round_csv)
}

# Elongating the CSV, so every round, nt_position and nt has their own row.
collected_csv_long <- collected_csv %>%
  pivot_longer(c(-round_index, -round_name, -nt_position), names_to="nt", values_to="count")

# Calulation the mean nt-frequencies for every SELEX round
collected_csv_mean <- collected_csv_long %>%
  group_by(round_name, round_index, nt) %>%
  summarize(
    count = mean(count)
  ) %>%
  group_by(round_name, round_index) %>%
  mutate(
    perc = count/sum(count)
  )

# Plotting nt distribution averaged over SELEX rounds
gg_mean <- ggplot(collected_csv_mean, aes(x=round_name, y=perc, fill=factor(nt))) +
  geom_col(position="fill") +
  geom_text(
    aes(label = format(perc*100, digit=2, nsmall = 1, group = factor(nt))),
    position=position_stack(vjust = 0.5),
    size=4.1,
    fontface="bold",
    colour="white"
  ) +
  labs(fill = "", colour="") +
  xlab("SELEX Round") +
  ylab("Nucleotide proportions") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_viridis_d(direction=-1, begin=0.1, end=0.9) +
  theme_bw() +
  ggtitle("Nucleotide Distribution") +
  theme(axis.text = element_text(size=11,face="bold"),
        axis.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        title = element_text(size=15, face="bold"))

# Saving the plot to file
ggsave("nt_distribution.png", gg_mean)

# Plotting nt distribution for every single SELEX round
gg_nt_content_rounds <- list()
for (round_index in levels(factor(collected_csv_long$round_index))) {
  round_csv <- collected_csv_long[collected_csv_long$round_index == round_index,]
  round_name <- levels(factor(round_csv$round_name))
  
  gg_nt_content_round <- ggplot(round_csv, aes(x=factor(nt_position, ordered=FALSE), y=count, fill=factor(nt))) +
    geom_col(position = "fill") + 
    labs(fill = "", colour="") +
    xlab("Nucleotide position") +
    ylab("Nucleotide proportions") +
    scale_x_discrete(breaks=c(seq(from=0, to=max(round_csv$nt_position), by=5 ))) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_fill_viridis_d(direction=-1, begin=0.1, end=0.9) +
    theme_bw() +
    ggtitle(paste0("Nucleotide Distribution for ", round_name)) +
    theme(axis.text = element_text(size=11,face="bold"),
          axis.title = element_text(size=12, face="bold"),
          legend.text = element_text(size=12, face="bold"),
          title = element_text(size=15, face="bold"))
  gg_nt_content_rounds[[round_index]] = gg_nt_content_round # Use = instead of <- for pointer purposes
  
  # Saving plot to file
  ggsave(paste0(round_name, ".nt_distribution.png"), gg_nt_content_round)
}

# Preparing HTML rendering
args_rmd = list()
args_rmd$author = args$author
args_rmd$gg_mean <- gg_mean
args_rmd$gg_rounds = gg_nt_content_rounds
args_rmd$round_names <- unique(factor(collected_csv_long$round_name))
args_rmd$round_index <- unique(factor(collected_csv_long$round_index))

# Rendering HTML
rmarkdown::render(
  here('bin', 'selex_nt_composition_plot.Rmd'), 
  params=args_rmd, 
  output_file=args$o,
  output_dir="."
)
