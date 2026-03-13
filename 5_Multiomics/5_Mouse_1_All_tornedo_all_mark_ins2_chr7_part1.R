# ===========================
# R Script: ChIP-seq full data + Ins2 ±250kb
# ===========================
#Epichaos env
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)

# ===========================
# 1️⃣ Set BigWig files
# ===========================
bw_files <- list(
  H3K4me3_WT1 = "Data/GSM5514712_307.bw",
  H3K4me3_WT2 = "Data/GSM5514713_407.bw",
  H3K4me3_KO1 = "Data/GSM5514714_308.bw",
  H3K4me3_KO2 = "Data/GSM5514715_408.bw",
  
  H3K4me1_WT1 = "Data/GSM5514716_305.bw",
  H3K4me1_WT2 = "Data/GSM5514717_405.bw",
  H3K4me1_KO1 = "Data/GSM5514718_306.bw",
  H3K4me1_KO2 = "Data/GSM5514719_406.bw",
  
  H3K27ac_WT1 = "Data/GSM5514720_401.bw",
  H3K27ac_WT2 = "Data/GSM5514721_411.bw",
  H3K27ac_KO1 = "Data/GSM5514722_402.bw",
  H3K27ac_KO2 = "Data/GSM5514723_412.bw",
  
  H3K27me3_WT1 = "Data/GSM5514724_303.bw",
  H3K27me3_WT2 = "Data/GSM5514725_403.bw",
  H3K27me3_KO1 = "Data/GSM5514726_304.bw",
  H3K27me3_KO2 = "Data/GSM5514727_404.bw"
)

# ===========================
# 2️⃣ Import full BigWig signal
# ===========================
import_bw_full <- function(bw_file) import(BigWigFile(bw_file))
bw_signal_list <- lapply(names(bw_files), function(name) {
  message("Importing BigWig: ", name, " -> ", bw_files[[name]])
  import(BigWigFile(bw_files[[name]]))
})
names(bw_signal_list) <- names(bw_files)


# ===========================
# 3️⃣ Tile genome-wide data (for simplicity, use 10kb bins across chromosomes)
# NOTE: For memory efficiency, you can subset chromosomes of interest
# ===========================
# Example: here we only create bins for chr7 (for demo, extend to all chromosomes if needed)
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
# 4️⃣ Extract genome-wide signal per bin
# ===========================
extract_signal_from_list <- function(track_list, bins) {
  res <- list()
  for (track_name in names(track_list)) {
    message("Processing track: ", track_name)
    track <- track_list[[track_name]]
    
    ov <- findOverlaps(bins, track)
    signal_vec <- sapply(seq_along(bins), function(i) {
      idx <- subjectHits(ov)[queryHits(ov) == i]
      if(length(idx) == 0) return(0)
      mean(track$score[idx])
    })
    
    res[[track_name]] <- signal_vec
  }
  df <- as.data.frame(res)
  df$bin_id <- bins$bin_id
  df$position_kb <- bins$position_kb
  return(df)
}



full_signal_mat <- extract_signal_from_list(bw_signal_list, genome_bins)
### done uptill above
### working froim down


full_signal_mat <- as.data.frame(full_signal_mat)
full_signal_mat$bin_id <- genome_bins$bin_id
full_signal_mat$position_kb <- genome_bins$position_kb

head(full_signal_mat)
### try to apply this score in this full_signal_mat
#Define the normalization function
norm_q <- function(x) {
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

full_signal_mat_norm <- full_signal_mat

cols_to_norm <- 1:(ncol(full_signal_mat) - 2)

full_signal_mat_norm[, cols_to_norm] <-
  apply(full_signal_mat[, cols_to_norm], 2, norm_q)

head(full_signal_mat_norm)

cols <- 1:(ncol(full_signal_mat_norm) - 2)

min_max_df <- data.frame(
  column = colnames(full_signal_mat_norm)[cols],
  min = apply(full_signal_mat_norm[, cols], 2, min, na.rm = TRUE),
  max = apply(full_signal_mat_norm[, cols], 2, max, na.rm = TRUE)
)



# ===========================
# 5️⃣ Save genome-wide signal as RDS
# ===========================
dir.create("Ins2_results_mouse", showWarnings = FALSE)
saveRDS(full_signal_mat, "Ins2_results_mouse/full_genome_signal_chr7_1kb.rds")

saveRDS(full_signal_mat_norm, "Ins2_results_mouse/full_genome_signal_chr7_1kb_normalised.rds")



