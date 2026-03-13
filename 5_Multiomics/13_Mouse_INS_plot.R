# -------------------------------
# Load packages
# -------------------------------
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(zoo)

# -------------------------------
# Parameters
# -------------------------------
cis_ext <- 250000     # ±250 kb around region
smooth_window <- 50   # number of bins for smoothing

# -------------------------------
# Load data
# -------------------------------
# Big signal GRanges (must have 'score' column)
data <- readRDS("Mouse/GSM5514712_307.bw.score.rds")

# Gene coordinates (data.frame with chr, start, end, symbol)
gene_coords_df <- readRDS("mm10_all_gene_coordinates.rds")

# Target region (INS1)
ins1_region <- GRanges(
  seqnames = "chr19",
  ranges = IRanges(
    start = 52264297,
    end   = 52265015
  )
)

# -------------------------------
# Create cis window around INS1
# -------------------------------
cis_window <- GRanges(
  seqnames = seqnames(ins1_region),
  ranges = IRanges(
    start = start(ins1_region) - cis_ext,
    end   = end(ins1_region) + cis_ext
  )
)

# -------------------------------
# Subset signal within cis window
# -------------------------------
subset_data <- subsetByOverlaps(data, cis_window)

# Convert to data.frame for plotting
df <- as.data.frame(subset_data) %>%
  mutate(
    mid = start + floor((end - start)/2)
  )

# -------------------------------
# Smooth signal using rolling mean
# -------------------------------
df$score_smooth <- zoo::rollapply(df$score,
                                  width = smooth_window,
                                  FUN = mean,
                                  fill = NA,
                                  align = "center")

# -------------------------------
# Filter genes in this window
# -------------------------------
genes_in_window <- gene_coords_df %>%
  filter(chr == as.character(seqnames(ins1_region))) %>%
  filter(start <= (end(ins1_region) + cis_ext) &
           end >= (start(ins1_region) - cis_ext))

# -------------------------------
# Plot signal + genes
# -------------------------------
library(ggplot2)
library(dplyr)
library(zoo)

# df = smoothed signal data.frame
# genes_in_window = genes in the same window

# Add a layer y-offset for gene track
genes_in_window <- genes_in_window %>%
  mutate(y_gene = max(df$score_smooth, na.rm = TRUE) + 5)  # position above signal

ggplot() +
  # 1️⃣ Signal layer
  geom_line(data = df, aes(x = mid, y = score_smooth), color = "blue", size = 1) +
  geom_point(data = df, aes(x = mid, y = score), alpha = 0.3, size = 0.5) +
  
  # 2️⃣ Gene segments layer (horizontal bars)
  geom_segment(data = genes_in_window,
               aes(x = start, xend = end, y = y_gene, yend = y_gene),
               color = "red", size = 2) +
  
  # 3️⃣ Gene names on top of bars
  geom_text(data = genes_in_window,
            aes(x = (start + end)/2, y = y_gene + 2, label = symbol),
            angle = 45, hjust = 0, size = 3) +
  
  # 4️⃣ Theme and labels
  theme_minimal() +
  labs(
    x = paste0("Genomic coordinate (", as.character(seqnames(ins1_region)), ")"),
    y = "Score",
    title = "Signal around INS1 with gene annotations"
  ) +
  ylim(0, max(df$score_smooth, na.rm = TRUE) + 15)
