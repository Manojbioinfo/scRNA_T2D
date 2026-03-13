# =============================================
# Full RLM eQTM pipeline: Alpha & Beta RNA vs ATAC
# =============================================
#conda activate epichoas
#https://stuartlab.org/signac/articles/pbmc_multiomic

#conda activate epiCHAOS
# Function to add ATAC assay from BigWig
# --------------------------
library(rtracklayer)

add_atac_from_bwfile <- function(bw_file, out_dir = "Mouse") {
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # Load BigWig
  bw <- import(bw_file)
  genome(bw) <- "mm10"
  
  # Construct output filename
  out1 <- gsub("^Data", out_dir, bw_file)
  out1 <- gsub("\\.bw$", ".bw.score.rds", out1)
  
  saveRDS(bw, out1)
  
  return(out1)
}




bw_files <- list.files(
  "Data",
  pattern = "GSM.*\\.bw$",
  full.names = TRUE
)

out_files <- lapply(bw_files, add_atac_from_bwfile)
