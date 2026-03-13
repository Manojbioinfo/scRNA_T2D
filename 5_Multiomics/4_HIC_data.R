rm(list=ls())


##############
# ================================
# Load libraries
# ================================
library(data.table)

# ================================
# Read input files
# ================================
alpha <- fread("Data/GSE195522_HiC_alpha.anchor_to_anchor.loops.txt.gz")
beta  <- fread("Data/GSE195522_HiC_beta.anchor_to_anchor.loops.txt.gz")

# ================================
# Function to parse "chr:start-end"
# ================================
parse_anchor <- function(x) {
  chr   <- sub(":.*", "", x)
  start <- as.integer(sub(".*:(\\d+)-.*", "\\1", x))
  end   <- as.integer(sub(".*-(\\d+)", "\\1", x))
  list(chr = chr, start = start, end = end)
}

# ================================
# Parse anchors for ALPHA
# ================================
a1 <- parse_anchor(alpha$V1)
a2 <- parse_anchor(alpha$V2)

alpha[, `:=`(
  chr1   = a1$chr,
  start1 = a1$start,
  end1   = a1$end,
  chr2   = a2$chr,
  start2 = a2$start,
  end2   = a2$end
)]

# ================================
# Parse anchors for BETA
# ================================
b1 <- parse_anchor(beta$V1)
b2 <- parse_anchor(beta$V2)

beta[, `:=`(
  chr1   = b1$chr,
  start1 = b1$start,
  end1   = b1$end,
  chr2   = b2$chr,
  start2 = b2$start,
  end2   = b2$end
)]

# ================================
# Define target region
# ================================
center_start <- 2181009
center_end   <- 2182439
extension    <- 250000

start_pos <- max(1, center_start - extension)
end_pos   <- center_end + extension
chrom     <- "chr11"

# ================================
# Filter loops overlapping region
# (either anchor overlaps)
# ================================
alpha_region <- alpha[
  (chr1 == chrom & end1 >= start_pos & start1 <= end_pos) |
    (chr2 == chrom & end2 >= start_pos & start2 <= end_pos)
]

beta_region <- beta[
  (chr1 == chrom & end1 >= start_pos & start1 <= end_pos) |
    (chr2 == chrom & end2 >= start_pos & start2 <= end_pos)
]

# ================================
# Keep only original columns (V1–V5)
# ================================
alpha_out <- alpha_region[, .(V1, V2, V3, V4, V5)]
beta_out  <- beta_region[,  .(V1, V2, V3, V4, V5)]


dir.create("AnchorHIC")

# ================================
# Write output files
# ================================
fwrite(
  alpha_out,
  file = "AnchorHIC/alpha_loops_INS_250kb_region.txt",
  sep = "\t"
)

fwrite(
  beta_out,
  file = "AnchorHIC/beta_loops_INS_250kb_region.txt",
  sep = "\t"
)

# ================================
# Done
# ================================
cat("Finished writing files:\n")
cat(" - alpha_loops_chr11_250kb_region.txt\n")
cat(" - beta_loops_chr11_250kb_region.txt\n")

