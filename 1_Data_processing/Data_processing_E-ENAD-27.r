

#https://biocellgen-public.svi.edu.au/mig_2019_scrnaseq-workshop/trajectory-inference.html
#install.packages(c("GSEABase","AUCell"))
#install.packages("SeuratObject")

library(Seurat)
library(Matrix)
#https://bioconductor.org/books/3.12/OSCA/cell-type-annotation.html#assigning-cell-labels-from-gene-sets
#https://igordot.github.io/clustermole/articles/example-bm-seurat.html
library(ggsci)



matrix_dir = "Data//"
barcode.path <- paste0(matrix_dir, "E-ENAD-27.aggregated_filtered_counts.mtx_cols")
features.path <- paste0(matrix_dir, "E-ENAD-27.aggregated_filtered_counts.mtx_rows")
matrix.path <- paste0(matrix_dir, "E-ENAD-27.aggregated_filtered_counts.mtx")
mat <- readMM(file = matrix.path)


feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1
mat=as.data.frame(mat)




counts_df <- as.data.frame(mat)
colnames(counts_df ) = barcode.names$V1
rownames(counts_df ) = feature.names$V1
library("AnnotationDbi")
library("org.Hs.eg.db")
library("clusterProfiler")
counts_df$gene_name <- mapIds(org.Hs.eg.db,keys=rownames(counts_df),column="SYMBOL",keytype="ENSEMBL",multiVals="first")
counts_df$gene_id <- rownames(counts_df)
dup_genes <- counts_df[duplicated(counts_df$gene_name),]
counts_df$feature <- ifelse((is.na(counts_df$gene_name) | counts_df$gene_name %in% dup_genes), counts_df$gene_id, counts_df$gene_name)
counts_df <- counts_df[complete.cases(rownames(counts_df)), ]
rownames(mat) <- make.unique(counts_df$feature)
head(mat [,1:3]) 



seurat_obj<- CreateSeuratObject(counts = mat, project = "E-ENAD-27", min.cells = 3, min.features = 100)

seurat_ids <- colnames(seurat_obj)

pheno_file=read.csv("ExpDesign-E-ENAD-27.tsv",sep = "\t")




# Extract IDs from the phenotype file
pheno_ids <- pheno_file$Assay  # Replace 'ID' with the actual column name in your phenotype file

# Match the IDs between Seurat and phenotype
matched_ids <- intersect(seurat_ids, pheno_ids)


# Arrange Seurat object and phenotype file based on the matched IDs
seurat_matched <- seurat_obj[, match(matched_ids, seurat_ids)]
pheno_matched <- pheno_file[match(matched_ids, pheno_ids), ]

# Now, seurat_matched and pheno_matched have the same order of IDs

matches <- sum(colnames(seurat_matched) %in% pheno_matched$Assay)

# Print the number of matches
print(matches)



##################

srat=seurat_matched
pheno2=pheno_matched
srat@meta.data$Age=pheno2$Sample.Characteristic.age.
srat@meta.data$Gender=pheno2$Sample.Characteristic.sex.

table(pheno2$Sample.Characteristic.disease.)


pheno2$Sample.Characteristic.disease.=gsub("normal","Healthy",pheno2$Sample.Characteristic.disease.)
pheno2$Sample.Characteristic.disease.=gsub("type II diabetes mellitus","T2D",pheno2$Sample.Characteristic.disease.)



## For simplicity again we will create a new column with name Disease and assign Sample.Characteristic.disease. 
## value to it
srat@meta.data$Disease=pheno2$Sample.Characteristic.disease.

table(srat$Disease)

saveRDS(srat,"27.seurat.org.rds",compress="xz")



data1=readRDS("27.seurat.org.rds")
dim(data1)
data1$nCount_RNA
