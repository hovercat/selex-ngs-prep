#!/usr/bin/env Rscript

library("argparse")
library("dplyr")
library("tidyr")
library("latex2exp")
library("ggplot2")

# Parameter validation
args_parser = ArgumentParser(
  prog="SELEXanalyse_log_duplicates",
  description=""
)
args_parser$add_argument("--in-csv", "-i", help="Input file must be csv as created by SELEX_dereplicate SELEX_rpm", required=TRUE)
args_parser$add_argument("--log-base", "-l", type="integer", default=10, help="Bins are spaced logarithmically, e.g. for log-base of 10: [1, 10, 100, 1000, ...]. Default: 10")
args_parser$add_argument("--min-abundancy", "-m", type="integer", default=1, help="Determines how often the sequence must be present to be considered. Default: 1")
args_parser$add_argument("--out-unique-csv", default="selex_success.unqiue.csv", help="")
args_parser$add_argument("--out-csv", default="selex_success.csv", help="")
#args <- args_parser$parse_args(c("-i", "work/b4/5cc0b1e226b06da72db1b079a33c17/selex.aptamers.rpm.csv", "-l", "2"))
args <- args_parser$parse_args()


# Data input
df <-  read.csv(args$in_csv, sep="\t", header=TRUE) %>%
  pivot_longer(-sequence, names_to="round", values_to="count") %>%
  filter(count > 0) %>%
  group_by(round) %>% # calculate exponent
  mutate(
    exponent = floor(log(count, base=args$log_base))
  ) %>%
  ungroup()

round_names <- levels(factor(df$round))
crossed <- round_names %>%
  crossing(1:max(df$exponent))
colnames(crossed) <- c("round", "exponent")

df_binning <- df %>%
  group_by(round, exponent) %>% 
  summarize(bin_size_total=sum(count), bin_size_unique=n()) %>%
  right_join(crossed)
df_binning[is.na(df_binning)] <- 0


selex_success <- df_binning %>% 
  pivot_wider(exponent, names_from=round, values_from=bin_size_total) %>%
  arrange(exponent)

write.table(
  selex_success,
  args$out_csv,
  append=FALSE,
  quote=FALSE,
  row.names = FALSE,
  sep="\t",
  na="0"
)

gg <- ggplot(data=df_binning, mapping = aes(x=round, y=bin_size_total, fill=factor(exponent))) +
  geom_col(position="fill") +
  labs(fill="") +
  xlab("SELEX Round") +
  ylab("% contribution of frequency category to total reads") +
  scale_fill_viridis_d(direction=-1, begin=0, end=1, option="magma", labels=unname(TeX(paste0(args$log_base, "^{", levels(factor(df_binning$exponent)), "}")))) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(axis.text = element_text(size=11,face="bold"),
        axis.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=11, face="bold"))

ggsave(paste0(args$out_csv, ".png"), gg)

selex_success_unique <- df_binning %>% 
  pivot_wider(exponent, names_from=round, values_from=bin_size_unique)

write.table(
  selex_success_unique,
  args$out_unique_csv,
  append=FALSE,
  quote=FALSE,
  row.names = FALSE,
  sep="\t",
  na="0"
)

gg_unique <- ggplot(data=df_binning, mapping = aes(x=round, y=bin_size_unique, fill=factor(exponent))) +
  geom_col(position="fill") +
  labs(fill="") +
  xlab("SELEX Round") +
  ylab("% contribution of frequency category to total reads") +
  scale_fill_viridis_d(direction=-1, begin=0, end=1, option="magma", labels=unname(TeX(paste0(args$log_base, "^{", levels(factor(df_binning$exponent)), "}")))) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(axis.text = element_text(size=11,face="bold"),
        axis.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=11, face="bold"))

ggsave(paste0(args$out_unique_csv, ".png"), gg_unique)
