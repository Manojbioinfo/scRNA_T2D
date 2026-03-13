
##conda activate epiCHAOS
library(GenomicRanges)
library(GenomeInfoDb)
library(rtracklayer)

# ---------------------------------------------------------
# 1. Load Alpha and Beta data
# ---------------------------------------------------------
alpha_data <- readRDS("senic/GSE194401_Alpha.bw.score.rds")
beta_data  <- readRDS("senic/GSE194401_Beta.bw.score.rds")

# ---------------------------------------------------------
# 2. Define INS gene region ±250 kb
#    INS (RefSeq: NC_000011.9)
# ---------------------------------------------------------
refseq_name  <- "NC_000011.9"
center_start <- 2181009
center_end   <- 2182439
extension    <- 250000  # 250 kb

safe_start <- max(1, center_start - extension)
safe_end   <- center_end + extension

ins_region <- GRanges(
  seqnames = refseq_name,
  ranges   = IRanges(start = safe_start, end = safe_end)
)

# ---------------------------------------------------------
# 3. Match chromosome naming style (UCSC)
# ---------------------------------------------------------
seqlevels(ins_region) <- "chr11"
seqlevelsStyle(ins_region) <- "UCSC"

seqlevelsStyle(alpha_data) <- "UCSC"
seqlevelsStyle(beta_data)  <- "UCSC"

# ---------------------------------------------------------
# 4. Extract INS region for Alpha and Beta
# ---------------------------------------------------------
alpha_ins <- subsetByOverlaps(alpha_data, ins_region)
beta_ins  <- subsetByOverlaps(beta_data,  ins_region)

cat("Alpha peaks in INS ±250kb:", length(alpha_ins), "\n")
cat("Beta peaks in INS ±250kb:",  length(beta_ins),  "\n")

# ---------------------------------------------------------
# 5. Export results
# ---------------------------------------------------------
export.bed(alpha_ins, "senic/INS_250kb_Alpha.bed")
export.bed(beta_ins,  "senic/INS_250kb_Beta.bed")


###############
