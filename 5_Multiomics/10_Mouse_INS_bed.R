## conda activate epiCHAOS
library(GenomicRanges)
library(GenomeInfoDb)
library(rtracklayer)

# ---------------------------------------------------------
# 1. BigWig files with sample metadata
# ---------------------------------------------------------
bw_files <- list(
  H3K4me3_WT1_GSM5514712 = "Data/GSM5514712_307.bw",
  H3K4me3_WT2_GSM5514713 = "Data/GSM5514713_407.bw",
  H3K4me3_KO1_GSM5514714 = "Data/GSM5514714_308.bw",
  H3K4me3_KO2_GSM5514715 = "Data/GSM5514715_408.bw",
  
  H3K4me1_WT1_GSM5514716 = "Data/GSM5514716_305.bw",
  H3K4me1_WT2_GSM5514717 = "Data/GSM5514717_405.bw",
  H3K4me1_KO1_GSM5514718 = "Data/GSM5514718_306.bw",
  H3K4me1_KO2_GSM5514719 = "Data/GSM5514719_406.bw",
  
  H3K27ac_WT1_GSM5514720 = "Data/GSM5514720_401.bw",
  H3K27ac_WT2_GSM5514721 = "Data/GSM5514721_411.bw",
  H3K27ac_KO1_GSM5514722 = "Data/GSM5514722_402.bw",
  H3K27ac_KO2_GSM5514723 = "Data/GSM5514723_412.bw",
  
  H3K27me3_WT1_GSM5514724 = "Data/GSM5514724_303.bw",
  H3K27me3_WT2_GSM5514725 = "Data/GSM5514725_403.bw",
  H3K27me3_KO1_GSM5514726 = "Data/GSM5514726_304.bw",
  H3K27me3_KO2_GSM5514727 = "Data/GSM5514727_404.bw"
)

# ---------------------------------------------------------
# 2. Output directory
# ---------------------------------------------------------
out_dir <- "Mouse_bed"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------
# 3. INS regions (±250 kb)
# ---------------------------------------------------------
ins1_region <- GRanges(
  seqnames = "chr19",
  ranges = IRanges(
    start = 52264297 - 250000,
    end   = 52265015 + 250000
  )
)

ins2_region <- GRanges(
  seqnames = "chr7",
  ranges = IRanges(
    start = 142678656 - 250000,
    end   = 142679726 + 250000
  )
)

seqlevelsStyle(ins1_region) <- "UCSC"
seqlevelsStyle(ins2_region) <- "UCSC"

# ---------------------------------------------------------
# 4. Process all BigWigs
# ---------------------------------------------------------
for (sample_id in names(bw_files)) {
  
  bw_path <- bw_files[[sample_id]]
  message("Processing: ", sample_id)
  
  # Import BigWig
  gr <- import(bw_path)
  genome(gr) <- "mm10"
  seqlevelsStyle(gr) <- "UCSC"
  
  # --------------------
  # INS1
  # --------------------
  ins1_hits <- subsetByOverlaps(gr, ins1_region)
  
  out_ins1 <- file.path(
    out_dir,
    paste0(sample_id, "_INS1_250kb.bed")
  )
  
  export(ins1_hits, out_ins1, format = "BED")
  
  # --------------------
  # INS2
  # --------------------
  ins2_hits <- subsetByOverlaps(gr, ins2_region)
  
  out_ins2 <- file.path(
    out_dir,
    paste0(sample_id, "_INS2_250kb.bed")
  )
  
  export(ins2_hits, out_ins2, format = "BED")
}

message("✅ All samples processed with GSM-traceable filenames")
