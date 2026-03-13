# =============================================
# Full RLM eQTM pipeline: Alpha & Beta RNA vs ATAC
# =============================================
#conda activate epichoas
#https://stuartlab.org/signac/articles/pbmc_multiomic

#conda activate epiCHAOS


rm(list = ls())
gc()
Sys.setenv(http_proxy = "http://proxy.mh-hannover.de:8080")
Sys.setenv(https_proxy = "http://proxy.mh-hannover.de:8080")
options(download.file.method = "curl")
options(download.file.extra = "-L --proxy http://proxy.mh-hannover.de:8080")





# --------------------------
# Required packages
# --------------------------
required_packages <- c(
  "Signac", "Seurat", "GenomicRanges", "future", "plyranges", 
  "rtracklayer", "GenomicAlignments", "dplyr", "Matrix", 
  "TxDb.Hsapiens.UCSC.hg19.knownGene", "BSgenome.Hsapiens.UCSC.hg19", "org.Hs.eg.db"
)

missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(missing_packages) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(missing_packages)
}
lapply(required_packages, library, character.only = TRUE)

# --------------------------
# Load scRNA-seq Seurat object
# --------------------------
exp <- readRDS("../2_GSE195986/2_DEG_Analysis/Anno1.rds")

# Subset Alpha and Beta
alpha_obj <- subset(exp, subset = celltype == "Alpha")
beta_obj  <- subset(exp, subset = celltype == "Beta")

# --------------------------
# Function to add ATAC assay from BigWig
# --------------------------
if(!dir.exists("senic")) dir.create("senic")



add_atac_from_bwfile <- function(seurat_obj, bw_file) {
  
  # Load BigWig
  bw <- import(bw_file)
  genome(bw) <- "hg19"
  
  out1=gsub("Data","senic",bw_file)
  out1=gsub(".bw",".bw.score.rds",out1)
  
  saveRDS(bw,out1)

}



add_atac_from_bw <- function(seurat_obj, bw_file) {
  
  # Load BigWig
  bw <- import(bw_file)
  genome(bw) <- "hg19"
  
  
  # Promoter regions ±2kb
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  promoters_gr <- promoters(genes(txdb), upstream = 2000, downstream = 500)
  promoters_gr <- promoters_gr[!is.na(start(promoters_gr)) & !is.na(end(promoters_gr))]
  promoters_gr_reduced <- reduce(promoters_gr)
  
  # Map BigWig signal
  hits <- findOverlaps(promoters_gr_reduced, bw)
  gene_scores_vec <- tapply(bw$score[subjectHits(hits)], queryHits(hits), sum)
  
  # Full vector
  full_scores <- numeric(length(promoters_gr_reduced))
  full_scores[as.numeric(names(gene_scores_vec))] <- as.numeric(gene_scores_vec)
  
  # Build counts matrix (one column per cell)
  num_cells <- ncol(seurat_obj)
  cell_names <- colnames(seurat_obj)
  
  counts_matrix <- Matrix(rep(full_scores, times = num_cells),
                          nrow = length(full_scores),
                          ncol = num_cells,
                          sparse = TRUE)
  rownames(counts_matrix) <- paste0("promoter-", seq_along(promoters_gr_reduced))
  colnames(counts_matrix) <- cell_names
  
  # Create ChromatinAssay
  hg19_info <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)
  chrom_assay <- CreateChromatinAssay(
    counts = counts_matrix,
    ranges = promoters_gr_reduced,
    genome = hg19_info,
    min.cells = -1,
    min.features = -1,
    validate.fragments = FALSE
  )
  
  

  
  
  seurat_obj[["ATAC"]] <- chrom_assay
  DefaultAssay(seurat_obj) <- "ATAC"
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  
  # Run basic Signac LSI workflow
  seurat_obj <- RunTFIDF(seurat_obj)
  seurat_obj <- FindTopFeatures(seurat_obj, min.cutoff = 20)
  seurat_obj <- RunSVD(seurat_obj)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:15, reduction = 'lsi') #same as DEGs analysis
  
  return(seurat_obj)
}

# --------------------------
# Add ATAC assays for Alpha and Beta
# --------------------------
alpha_bw <- "Data/GSE194401_Alpha.bw"
beta_bw  <- "Data/GSE194401_Beta.bw"


add_atac_from_bwfile(alpha_obj, beta_bw)
add_atac_from_bwfile(beta_obj, beta_bw)



alpha_obj <- add_atac_from_bw(alpha_obj, alpha_bw)

if(!dir.exists("senic")) dir.create("senic")
saveRDS(alpha_obj, file = "senic/Alpha_ATAC_RNA_Seurat.rds")

beta_obj  <- add_atac_from_bw(beta_obj, beta_bw)
saveRDS(beta_obj, file = "senic/Beta_ATAC_RNA_Seurat.rds")
# --------------------------
# Merge Alpha + Beta into one Seurat object
# --------------------------
combined_obj <- merge(alpha_obj, y = beta_obj,add.cell.ids = c("Alpha", "Beta"))
saveRDS(combined_obj, file = "senic/AlphaBeta_ATAC_RNA_Seurat.rds")
# --------------------------
# Final check
# --------------------------

