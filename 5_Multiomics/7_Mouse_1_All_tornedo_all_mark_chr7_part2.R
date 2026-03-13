# ===========================
# R Script: ChIP-seq full data + Ins2 ±250kb
# ===========================
#Epichaos env
#library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)

full_signal_mat=readRDS("Ins2_results_mouse/full_genome_signal_chr7_1kb.rds")


chr_length <-145441459 # mouse chr7 length (GRCm38) #https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&chromInfoPage=&token=0.tX-l5aBhpjzAcgi0NhvEaOm4woq0chzAvNKxZdi2cJ4tcSBaUKfEeC1XMF_XM_-QJ3dYS74V1tXfpMacGi5ABUZAhRMtvBmVCc4YkWRvrWrhEOLb3pVsmMWm2f_C8MC0LDjGIUw6_Jn1Bh4KBUjdURAOCi_gwn3AQAauifqa67qgv2z62o2IqdV8ubTGRo9Bjcn6Kl_v9tgjPN2hE96q9Wnt-i1aZ7WU_PEGZeK2MrPt2y8_yBu7-t88IYCWJLXvGLdeQG6LQtZSUT4KIK1YBXUAaJRPVfSJBqHBv0k26WoD-m4xiJ1fYmbGMGC9g3rDmY1DgS400PE1I9zTfV5B4GCJJM2isDe7ybJKIIR7uCZBd_5SPV_7V3EyQcwv1AbRUIbMu7vfoolxV8uYD34_YGhkkeJQ0WbRwU7Q8SFBxa3nMOiNjrTZ_9jyOYnIv0agXlrZpaFgzfrNQxa7HmrGZvmWdLuk_INxN9TKsAn2aHqtGwrhuyhTMUpUYjMDgbsEE1L1Ku83bwetpKVc-FyD7a387jCQ7jEGa-O0GGawWvtCmEW7snh6JtFN6qXZscWcKFPNZkwQkPJjR6FvKlQtMWeePcke2Mzv44mKdwnZX0u0e5bMn20XdXZTM2t4XA0bTRF9wRzY4sH2wI6tE-Rc_HWA97Niwdcpyvac-rNTd8PiPNgWVkDEuhaYqTFv8Nxhpyd1NqSGDSrR692Q9Zn45RYoaNpHWKu-eTsv_8F1z_5Yt2C4lcKfhg2wH7BXx2hZfnNPqj0MU7fPgGSkutZA46hGIljE8xvkq1-KzFirIMc5tpGxgxiVcjecCJ5scby8Ic9kzcW3pO7fcdm9FIyh47ykEGQZN9Il8J0CUfoyK0d9X9xj1i1czTOwO7IhojqmIrA1qCstyTtf3uTIIy2uCKr5hSsQKH9-GqdntrK_fv4.KB5T3wed60LJYnvtkRUoxg.a70a477bfa3428af8aee453ada06f8359f87ae660194ce78878f7089f4cd2cb6&pix=1900
bin_size <- 1000

genome_bins <- GRanges(
  seqnames = "chr7",
  ranges = IRanges(
    start = seq(1, chr_length, by = bin_size),
    end   = pmin(seq(1, chr_length, by = bin_size) + bin_size - 1, chr_length)
  )
)
genome_bins$bin_id <- seq_along(genome_bins)
genome_bins$position_kb <- start(genome_bins) / 1000



# ===========================
# 6️⃣ Subset Ins2 ±250kb from genome-wide bins
# ===========================
ins_region <- GRanges(
  seqnames = "chr7",
  ranges = IRanges(
    start = 142678656 - 250000,
    end   = 142679726 +250000
  )
)



ins_bins <- genome_bins[overlapsAny(genome_bins, ins_region), ]
ins_signal <- full_signal_mat[full_signal_mat$bin_id %in% ins_bins$bin_id, ]





saveRDS(ins_signal, "Ins2_results_mouse/Ins2_signal_1kb_bins.rds")

# ===========================
# 7️⃣ Prepare tornado plot (WT, KO, WT−KO)
# ===========================
signal_long <- ins_signal %>%
  pivot_longer(
    cols = -c(bin_id, position_kb),
    names_to = "sample",
    values_to = "signal"
  ) %>%
  mutate(
    mark = sub("_.*", "", sample),
    condition = ifelse(grepl("_WT", sample), "WT", "KO")
  )




diff_long <- signal_long %>%
  group_by(position_kb, mark, condition) %>%
  summarise(signal = mean(signal), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = signal) %>%
  mutate(WT_minus_KO = WT - KO)

plot_df <- diff_long %>%
  pivot_longer(cols = c("WT", "KO", "WT_minus_KO"),
               names_to = "condition",
               values_to = "signal")

plot_df$condition <- factor(plot_df$condition, levels = c("WT", "KO", "WT_minus_KO"))

saveRDS(plot_df, "Ins2_results_mouse/Ins2_plot_df_1kb.rds")

# ===========================
# 8️⃣ Save tornado plot PDF
# ===========================
library(dplyr)
library(ggplot2)

# 1. Remove WT_minus_KO
plot_df2 <- plot_df %>%
  filter(condition != "WT_minus_KO")

# 2. Compute a tornado ranking per bin (based on mean WT+KO signal)
rank_df <- plot_df2 %>%
  group_by(position_kb) %>%
  summarise(
    tornado_score = mean(signal, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(tornado_score) %>%
  mutate(tornado_rank = row_number())

# 3. Join ranks back
plot_df2 <- plot_df2 %>%
  left_join(rank_df, by = "position_kb")

# 4. Plot tornado heatmap (y = rank, NOT position)
pdf("Ins2_results_mouse/ChIPseq_Ins2_tornado_10kb.pdf", width = 12, height = 8)

ggplot(plot_df2, aes(x = condition, y = tornado_rank, fill = signal)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "navy",
    mid = "#FFEA99",
    high = "darkred",
    midpoint = median(plot_df2$signal, na.rm = TRUE)
  ) +
  facet_wrap(~ mark, ncol = 2) +
  labs(
    title = "Ins2 tornado plot (WT vs KO)",
    x = "Condition",
    y = "Ranked bins (signal-sorted)"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+scale_y_reverse()


dev.off()
