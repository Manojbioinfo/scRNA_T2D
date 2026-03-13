# -------------------------------
# Load packages
# -------------------------------
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(zoo)
library(tidyverse)
atac_alpha=read_tsv("senic/INS_250kb_Alpha.bed", col_names = F)
atac_alpha$Cell_Type="Alpha"
atac_beta=read_tsv("senic/INS_250kb_Beta.bed", col_names = F)
atac_beta$Cell_Type="Beta"

atac=rbind(atac_alpha, atac_beta)
head(atac)


atac <- atac %>%
  rename(
    chr = X1,
    start = X2,
    end = X3,
    name = X4,
    score = X5,
    strand = X6
    
  )
head(atac)

loop_alpha=read_tsv("AnchorHIC/alpha_loops_INS_250kb_region.txt", col_names = T)
loop_alpha$Cell_Type="Alpha"
loop_beta=read_tsv("AnchorHIC/beta_loops_INS_250kb_region.txt", col_names = T)
loop_beta$Cell_Type="Beta"

loop=rbind(loop_alpha, loop_beta)
head(loop)
colnames(loop) <- c("anchor1", "anchor2", "count", "score", "pvalue", "Cell_Type")
loop <- loop %>%
  select(-score, -pvalue)
head(loop)

library(openxlsx)
write.xlsx(atac,"Atac_human.xlsx")
write.xlsx(loop,"Loop_human.xlsx")



##########chipseq
#Histone Modification Enrichment Regions Around INS1 and INS2 Genes in Knockout (KO) Samples
list.files("Mouse_bed_merged/")
library(readr)
library(dplyr)
library(purrr)

# list files
files <- list.files("Mouse_bed_merged/", full.names = TRUE)

# read and combine
final_df <- files %>%
  map_dfr(~ read_tsv(.x, col_names = FALSE) %>%
            mutate(file_name = basename(.x)))

# view result
head(final_df)
library(dplyr)
library(stringr)

final_df <- final_df %>%
  mutate(file_name = str_remove(file_name, "_250kb_merged\\.bed")) %>%
  separate(file_name, into = c("Histone", "Condition", "INS"), sep = "_")
library(dplyr)

final_df <- final_df %>%
  rename(
    chrom = X1,
    start = X2,
    end = X3,
    score = X4
  ) %>%
  mutate(gene = INS)%>%select((-INS))
head(final_df)
table(final_df$gene)
write.xlsx(final_df,"Chip_mouse.xlsx")
