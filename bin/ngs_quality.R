#!/usr/bin/env Rscript
library(argparse)
# library(Biostrings)
# library(ShortRead)
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(hrbrthemes)
library(dada2)


cmd_parser <- argparse::ArgumentParser(description="")
cmd_parser$add_argument('-i', type="character", required=TRUE)
#cmd_parser$add_argument('-t', type="character", required=TRUE, help="fwd or rev")
#cmd_parser$add_argument('-r', type="character", required=TRUE, help="round name")
cmd_parser$add_argument('-o', type="character", required=FALSE, default="./quality.png")
args <- cmd_parser$parse_args(c("-i", "../input_data/R0_S1_L001_R1_001.fastq"))

gg <- plotQualityProfile(args$i)
ggsave(args$o, gg)

# 
# fastq_seqs <- ShortRead::readFastq(args$i)
# fastq_quality <- as(quality(fastq_seqs), "matrix") %>%
#   reshape2::melt(varnames=c("aptamer_id", "nt_position"), value.name="quality")
# 
# # Adding dummy values for NA
# #nt_positions <- 1:max(fastq_quality$nt_position)
# #qualities <- 0:40
# #df_nt <- data.frame(nt_position = nt_positions)
# #df_quality <- data.frame(quality = qualities)
# #df_dummy <- df_nt %>% full_join(df_quality, by=character()) # cross join
# #df_dummy$aptamer_id <- NA
# #fastq_quality <- rbind(fastq_quality, df_dummy)
# 
# fastq_quality_raster <- fastq_quality %>%
#   group_by(nt_position, quality) %>%
#   summarize(
#     #n = sum(!is.na(aptamer_id))
#     n = n()
#   ) %>%
#   group_by(nt_position) %>%
#   mutate(
#     n_perc = n/sum(n)
#   )
# 
# quality_mean <- fastq_quality %>% group_by(nt_position) %>%
#   summarize(
#     mean_quality = mean(quality, na.rm = TRUE),
#     sd_quality = sd(quality, na.rm = TRUE),
#     n = n()
#   )
# 
# gg <- ggplot(fastq_quality_raster, aes(x=nt_position, y=quality, fill=n_perc)) +
#   geom_tile(na.rm=TRUE) +
#   scale_fill_viridis_c(option="D") +
#   theme_classic()
#   #scale_fill_gradientn(colors=c("yellow", "black"), values=c(0, 1))
# gg
# 
# gg <- ggplot(fastq_quality, aes(x=nt_position, y=quality, group=nt_position)) +
#   geom_boxplot(outlier.shape = NA)
# 
# gg
