#!/usr/bin/env Rscript
library("argparse")
library("dplyr")
library("tidyr")

# Parameter validation
args_parser = ArgumentParser(
  prog="SELEXrpm",
  description=""
)
args_parser$add_argument("--in-derep-csv", "-i", help="Input file must be csv as created by SELEXderep.", required=TRUE)
args_parser$add_argument("--out-file", "-o", default="rpm.csv", help="Output file. Default: rpm.csv")
#args <- args_parser$parse_args(c("-i", "output/fasta.prepped/out.csv"))
args <- args_parser$parse_args()

# Data input
round_names <- c()
df_selex <-  read.csv(args$in_derep_csv, sep="\t", header=TRUE)
rownames(df_selex) <- df_selex$sequence
round_names <- colnames(df_selex)[-1]

# Set every sequence count which is 0 to NA, so the long format is more sparse.
df_selex[df_selex == 0] = NA
df_selex_long <- df_selex %>%
  pivot_longer(-sequence, names_to="round", values_to="count", values_drop_na=TRUE)
#df_selex_long <- df_selex %>% gather("round", "count", -"id", -"seq", na.rm = TRUE)
df_selex <- NULL # clear wide format from ram

df_rpm <- df_selex_long %>%
  group_by(round) %>% # calculate rpm, exponent_rpm and exponent
  mutate(
    rpm = (count/sum(count) * (10^6))
  ) %>%
  ungroup() %>%
  pivot_wider(sequence, names_from=round, values_from=rpm, values_fill=0)


# write rpm to csv file
write.table(
  df_rpm,
  paste(args$out_file),
  append=FALSE,
  quote=FALSE,
  row.names = FALSE,
  sep="\t",
  na="0"
)
