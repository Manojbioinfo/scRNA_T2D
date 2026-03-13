
library(GenomicRanges)
library(rtracklayer)
library(GenomicRanges)
library(biomaRt)




library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(GenomicRanges)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

all_genes <- genes(txdb)

all_genes$symbol <- mapIds(
  org.Mm.eg.db,
  keys = all_genes$gene_id,
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
)

saveRDS(all_genes, "mm10_all_genes_granges.rds")




all_genes <- readRDS("mm10_all_genes_granges.rds")
gene_coords_df <- data.frame(
  chr    = as.character(seqnames(all_genes)),
  start  = start(all_genes),
  end    = end(all_genes),
  strand = as.character(strand(all_genes)),
  symbol = all_genes$symbol,
  entrez = all_genes$gene_id
)

saveRDS(gene_coords_df, "mm10_all_gene_coordinates.rds")







#######

hits <- findOverlaps(new_regions, all_genes)

annotated_genes <- unique(all_genes[subjectHits(hits)])



hits <- findOverlaps(regions, genes)

gene_ids <- unique(genes$gene_id[subjectHits(hits)])

symbols <- mapIds(
  org.Mm.eg.db,
  keys = gene_ids,
  keytype = "ENTREZID",
  column = "SYMBOL"
)

result <- data.frame(
  entrez_id = gene_ids,
  symbol = symbols
)



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
