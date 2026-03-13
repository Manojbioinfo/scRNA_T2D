library(Seurat)
library(data.table)

library(GEOquery)
# Load necessary libraries
library(data.table)




# Extract the required file from the tar archive
system("tar -xvf GSE195986_RAW.tar")



# List of files to process
files <- c("GSM5857297_HT1.dropseq.dge.txt.gz", "GSM5857298_HT2.dropseq.dge.txt.gz", 
           "GSM5857299_HT3.dropseq.dge.txt.gz", "GSM5857300_HT4.dropseq.dge.txt.gz", 
           "GSM5857301_HT5.dropseq.dge.txt.gz", "GSM5857302_HT6.dropseq.dge.txt.gz", 
           "GSM5857303_HT7.dropseq.dge.txt.gz", "GSM5857304_T2D1.dropseq.dge.txt.gz", 
           "GSM5857305_T2D2.dropseq.dge.txt.gz", "GSM5857306_T2D3.dropseq.dge.txt.gz", 
           "GSM5857307_T2D4.dropseq.dge.txt.gz")

# Loop through each file, extract and create Seurat objects
seurat_objects <- list()

for (file in files) {
  # Unzip each file
  system(paste("gunzip", file))
  
  # Load the count matrix into R
  expression_matrix <- fread(sub(".gz", "", file), sep = "\t", header = TRUE, data.table = FALSE)
  
  # Convert to matrix format
  gene_names <- expression_matrix[[1]]  # First column contains gene names
  counts <- as.matrix(expression_matrix[, -1])  # Remove gene column
  rownames(counts) <- gene_names
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = counts, project = sub(".dropseq.dge.txt", "", file))
  
  # Store the Seurat object in a list
  seurat_objects[[sub(".dropseq.dge.txt.gz", "", file)]] <- seurat_obj
}

# Now seurat_objects is a list where each entry is a Seurat object for each file
# For example, to access the Seurat object for HT1, use:
# seurat_objects$GSM5857297_HT1


library(scCustomize)
library(Seurat)
library(GEOquery)

val_id=gsub(".dropseq.dge.txt.gz","",files)

merge_data=seurat_objects
merge_data <- Merge_Seurat_List(list_seurat = merge_data,add.cell.ids = val_id)
merge_data 
# Join the layers after merging
merge_data[["RNA"]] <- JoinLayers(merge_data[["RNA"]])
merge_data 

pheno=merge_data@meta.data
# Assuming the 'pheno' data frame is already loaded

# Create a new column 'Disease' and classify based on the 'orig.ident' column
pheno$Disease <- ifelse(grepl("HT", pheno$orig.ident), "Healthy", "T2D")

# View the updated pheno data frame
head(pheno)

head(pheno)


table(pheno$Disease)

merge_data@meta.data=pheno


saveRDS(merge_data,"GSE195986.scRNA.rds")
