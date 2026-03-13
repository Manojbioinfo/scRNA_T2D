# ========================================
# ATAC signal density plotting (Alpha & Beta) - COMPLETE
# ========================================
rm(list = ls())

library(tidyverse)
library(GenomicRanges)
library(GenomeInfoDb)
library(zoo)
library(rtracklayer)

# -------------------------------
# 1. Define INS Region ±50 kb
# -------------------------------
center_start <- 2181009
center_end   <- 2182439
extension    <- 250000    # 250kb
min_score_threshold <- 40
smooth_window <- 500

start_pos <- max(1, center_start - extension)
end_pos   <- center_end + extension

target_region <- GRanges(
  seqnames = "chr11",
  ranges   = IRanges(start = start_pos, end = end_pos)
)
seqlevelsStyle(target_region) <- "UCSC"

message("Region: chr11:", start_pos, "-", end_pos, " (", end_pos - start_pos + 1, " bp)")

# -------------------------------
# 2. Helper Function (Fixed)
# -------------------------------
plot_atac_coverage <- function(bed_file, label, color = "firebrick") {
  
  message("Processing: ", label)
  
  if (!file.exists(bed_file)) {
    stop("File not found: ", bed_file)
  }
  
  # Read BED
  bed <- read_tsv(bed_file, col_names = FALSE, show_col_types = FALSE)
  colnames(bed)[1:5] <- c("chr", "start", "end", "name", "score")
  
  # BED → GRanges
  gr <- GRanges(
    seqnames = bed$chr,
    ranges   = IRanges(start = bed$start, end = bed$end),
    score    = bed$score
  )
  seqlevelsStyle(gr) <- "UCSC"
  
  # Subset region + filter
  gr_subset <- subsetByOverlaps(gr, target_region)
  gr_clean <- gr_subset[gr_subset$score >= min_score_threshold]
  
  message("  Total peaks: ", length(gr_subset), " → Clean: ", length(gr_clean))
  
  if (length(gr_clean) == 0) {
    warning("No peaks found for ", label)
    return(NULL)
  }
  
  # Coverage with full chromosome length
  hg19_chr11_len <- 135006516 
  cov_rle <- coverage(
    gr_clean, 
    weight = "score", 
    width = setNames(hg19_chr11_len, "chr11")
  )[["chr11"]]
  
  # Safe window extraction
  safe_end <- min(end_pos, length(cov_rle))
  window_signal <- as.numeric(cov_rle[start_pos:safe_end])
  
  # Pad if needed
  needed_length <- end_pos - start_pos + 1
  if (length(window_signal) < needed_length) {
    window_signal <- c(window_signal, rep(0, needed_length - length(window_signal)))
  }
  
  # Data frame
  df_cov <- tibble(
    position = seq(start_pos, end_pos),
    signal = window_signal
  ) %>%
    mutate(
      signal_smooth = zoo::rollmean(signal, k = smooth_window, fill = 0)
    )
  
  # -------------------------------
  # Raw Plot
  # -------------------------------
  p_raw <- ggplot(df_cov, aes(position, signal)) +
    geom_area(fill = alpha(color, 0.8)) +
    labs(
      title = paste(label, "- Raw (>", min_score_threshold, ")"),
      subtitle = paste("chr11:", start_pos, "-", end_pos),
      x = "Genomic position (bp)",
      y = "Signal"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(size = 12))
  
  ggsave(paste0("senic/",label, "_ATAC_raw_chr11.pdf"), p_raw, width = 10, height = 3)
  
  # -------------------------------
  # Smoothed Plot (FIXED geom_vline)
  # -------------------------------
  p_smooth <- ggplot(df_cov, aes(position, signal_smooth)) +
    geom_area(fill = alpha(color, 0.3)) +
    geom_line(color = color, linewidth = 1.2) +
    
    # INS TSS MARKER (NOW CORRECTLY PLACED)
    geom_vline(xintercept = center_start, linetype = "dashed", 
               color = "red", alpha = 0.8, linewidth = 0.8) +
    
    labs(
      title = paste(label, "- Smoothed (k =", smooth_window, ")"),
      subtitle = paste("chr11:", start_pos, "-", end_pos, "| Red = INS TSS"),
      x = "Genomic position (bp)",
      y = "Smoothed Signal",
      caption = paste("Peaks >=", min_score_threshold, "|", length(gr_clean), "peaks")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12),
      plot.caption = element_text(size = 8)
    )
  
  ggsave(paste0("senic/",label, "_ATAC_smoothed_chr11.pdf"), p_smooth, width = 10, height = 4)
  
  # Save data
  saveRDS(df_cov, paste0("senic/",label, "_ATAC_coverage_chr11.rds"))
  
  return(df_cov)
}

# -------------------------------
# 3. Run Analysis
# -------------------------------
alpha_df <- plot_atac_coverage("senic/INS_250kb_Alpha.bed", "Alpha", "#E69F00")
beta_df  <- plot_atac_coverage("senic/INS_250kb_Beta.bed",  "Beta",  "#56B4E9")

message("🎉 All plots generated successfully!")

# -------------------------------
# 4. BONUS: Overlay Alpha + Beta
# -------------------------------
if (!is.null(alpha_df) && !is.null(beta_df)) {
  df_combined <- bind_rows(
    alpha_df %>% mutate(cell_type = "Alpha"),
    beta_df  %>% mutate(cell_type = "Beta")
  )
  
  p_overlay <- ggplot(df_combined, aes(position, signal_smooth, color = cell_type)) +
    geom_line(linewidth = 1.5) +
    geom_vline(xintercept = center_start, linetype = "dashed", 
               color = "red", alpha = 0.8) +
    labs(
      title = "Alpha vs Beta: INS Locus (±50kb)",
      subtitle = "Smoothed ATAC signal (score > 40)",
      x = "Genomic position (bp)",
      y = "Smoothed Signal",
      color = "Cell Type"
    ) +
    scale_color_manual(values = c("Alpha" = "#E69F00", "Beta" = "#56B4E9")) +
    theme_minimal()
  
  ggsave("senic/Alpha_vs_Beta_INS_overlay.pdf", p_overlay, width = 12, height = 5)
  message("✓ Overlay plot saved!")
}
