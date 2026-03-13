#https://biocellgen-public.svi.edu.au/mig_2019_scrnaseq-workshop/trajectory-inference.html
#install.packages(c("GSEABase","AUCell"))
#install.packages("SeuratObject")

# Sys.setenv(CONDA_DEFAULT_ENV = "r4.4_env")
# current_env <- Sys.getenv("CONDA_DEFAULT_ENV")
# print(current_env)
# 
# Sys.setenv(http_proxy = "http://proxy.mh-hannover.de:8080")
# Sys.setenv(https_proxy = "http://proxy.mh-hannover.de:8080")
# options(download.file.method = "curl")
# options(download.file.extra = "-L --proxy http://proxy.mh-hannover.de:8080")
# # 
# # Sys.setenv(LD_LIBRARY_PATH = "/hpc/leinehome/guptaman/miniconda3/envs/r4.4_env/lib")
# # 

# Load necessary libraries
library(Seurat)
library(ggplot2)
library(ggsci)


library(Seurat)
library(Matrix)
library(ggsci)
library(ggplot2)
set.seed(1234)

rm(list=ls())


srat=readRDS("../1_Data/81608.seurat.org.rds")
srat@meta.data$Disease=gsub( "type 2 diabetes mellitus","T2D",srat@meta.data$Disease)
unique(srat@meta.data$Disease)
  
  
srat
# srat
# An object of class Seurat 
# 28616 features across 1600 samples within 1 assay 
# Active assay: RNA (28616 features, 0 variable features)
# 1 layer present: counts
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")

#srat$percent.mt




pdf("Pic.Manoj1.pdf", width =10, height =3)
VlnPlot(srat, features = c("nFeature_RNA","nCount_RNA"),ncol =2,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10),
        axis.text.x = element_text(angle=90))

dev.off()


jpeg("Pic.Manoj1.1.jpeg", width =5, height =5, units = "in", res = 600)
# Generate violin plot with ggsci color
VlnPlot(srat, features = c("nFeature_RNA","nCount_RNA"),ncol =2,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10),
        axis.text.x = element_text(angle=90))
dev.off()

jpeg("Pic.Manoj1.2.jpeg", width =10, height =5, units = "in", res = 600)
# Generate violin plot with ggsci color
VlnPlot(srat, features = c("nFeature_RNA","nCount_RNA" ),ncol =2,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10),
        axis.text.x = element_text(angle=90))
dev.off()





# Remove rows with missing values in either nFeature_RNA or nCount_RNA

data=subset(srat,nFeature_RNA>200 & nFeature_RNA<10000  & percent.mt < 25 )  ### This depends on the data quality
data
# 
# An object of class Seurat 
# 28616 features across 1270 samples within 1 assay 
# Active assay: RNA (28616 features, 0 variable features)
# 1 layer present: counts


pdf("Pic.Manoj2.pdf", width =10, height =3)
VlnPlot(data, features = c("nFeature_RNA","nCount_RNA"),ncol = 2,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10),
        axis.text.x = element_text(angle=90))
dev.off()




jpeg("Pic.Manoj2.1.jpeg", width =5, height =5, units = "in", res = 600)
# Generate violin plot with ggsci color
VlnPlot(data, features = c("nFeature_RNA","nCount_RNA"),ncol = 2,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10),
        axis.text.x = element_text(angle=90))
dev.off()

jpeg("Pic.Manoj2.2.jpeg", width =10, height =5, units = "in", res = 600)
# Generate violin plot with ggsci color
VlnPlot(data, features = c("nFeature_RNA","nCount_RNA"),ncol = 2,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10),
        axis.text.x = element_text(angle=90))
dev.off()






data <- FindVariableFeatures(data)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(data), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf("Pic.Manoj3.pdf", width =10, height =5)
plot1 + plot2
dev.off()


data <- NormalizeData(data)
data <- ScaleData(data)

data <- RunPCA(data, verbose = T)
print(data[["pca"]], dims = 1:5, nfeatures = 5)
pdf("Pic.Manoj4.pdf", width =10, height =5)
VizDimLoadings(data, dims = 1:2, reduction = "pca")
dev.off()

pdf("Pic.Manoj5.pdf", width =5, height =5)
DimPlot(data, reduction = "pca")
dev.off()

pdf("Pic.Manoj6.pdf", width =10, height =5)
DimHeatmap(data, dims = 1, cells = 500, balanced = TRUE)
dev.off()

min_pcs=15
pdf("Pic.Manoj7.pdf", width =10, height =8)
DimHeatmap(data, dims = 1:min_pcs, cells = 500, balanced = TRUE)
dev.off()

pdf("Pic.Manoj8.pdf", width =10, height =5)
ElbowPlot(data)
dev.off()



################ 


min_pcs=15
data <- FindNeighbors(data, dims = 1:min_pcs) # 10 PCS
data <- FindClusters(data, resolution =  0.5, verbose = T)
head(Idents(data), 5)
data <- RunUMAP(data, dims = 1:min_pcs)

pdf("Pic.Manoj9.pdf", width =9, height =7)
DimPlot(data, reduction = "umap", label = T)  
dev.off()

jpeg("Pic.Manoj9.jpeg", width =10, height =5, units = "in", res = 600)
DimPlot(data, reduction = "umap", label = T)  
dev.off()



############### cell type annotation
#https://igordot.github.io/clustermole/articles/example-bm-seurat.html
library(clustermole)
so=data
avg_exp_mat <- AverageExpression(so)
avg_exp_mat <- as.matrix(avg_exp_mat$RNA)
avg_exp_mat <- log1p(avg_exp_mat)

enrich_tbl <- clustermole_enrichment(expr_mat = avg_exp_mat, species = "hs")
library(openxlsx)
table(enrich_tbl$cluster)


write.xlsx(enrich_tbl,"Celltype label.xlsx")


#####
library(dplyr)

# Function to get the best matching cell type for each cluster
get_best_celltype <- function(df) {
  # First, try to get Pancreas-related cell type
  best_match <- df %>%
    arrange(score_rank) %>%
    filter(organ %in% c("Pancreas", "Pancreatic islet")) %>%
    slice_head(n = 1)
  
  # If still no match, assign "Unknown"
  if (nrow(best_match) == 0) {
    best_match <- df %>%
      arrange(score_rank) %>%
      slice_head(n = 1)  # Take best overall match
    best_match$celltype <- "Unknown"
  }
  
  return(best_match)
}

# Apply function to each cluster and arrange in increasing order of cluster numbers
best_celltypes <- enrich_tbl %>%
  group_by(cluster) %>%
  group_modify(~ get_best_celltype(.x)) %>%
  ungroup() %>%
  mutate(cluster = as.numeric(gsub("g", "", cluster))) %>%  # Convert "g0", "g1", etc. to numeric
  arrange(cluster) %>%  # Arrange clusters in increasing order
  mutate(cluster = paste0("g", cluster))  # Convert back to "g0", "g1", etc.

# View the sorted best cell types
print(best_celltypes)
best_celltypes$celltype <- gsub(" cell", "", best_celltypes$celltype)
best_celltypes$celltype <- gsub(" ", "", best_celltypes$celltype)
best_celltypes$celltype <- gsub("Pancreaticstellates", "Stellates", best_celltypes$celltype)
best_celltypes$celltype <- gsub("Gamma\\(PP\\)s", "Gamma", best_celltypes$celltype)
best_celltypes$celltype <- gsub("Alphas", "Alpha", best_celltypes$celltype)
best_celltypes$celltype <- gsub("Deltas", "Delta", best_celltypes$celltype)
best_celltypes$celltype <- gsub( "Ductals",  "Ductal", best_celltypes$celltype)
best_celltypes$celltype <- gsub( "Stellates" ,  "Stellate" , best_celltypes$celltype)
best_celltypes$celltype <- gsub("Acinars", "Acinar", best_celltypes$celltype)
best_celltypes$celltype <- gsub("Cancerstem\\(PancreaticCancer\\)", "Unknown", best_celltypes$celltype)
best_celltypes$celltype <- gsub("Peri-isletSchwanns", "Unknown", best_celltypes$celltype)
best_celltypes$celltype <- gsub("Pancreaticpolypeptide", "Unknown", best_celltypes$celltype)

sort(unique(best_celltypes$celltype))


manoj_p1=data@meta.data


# Ensure cluster column in best_celltypes matches the format in manoj_p1
manoj_p1$cluster <- paste0("g", manoj_p1$seurat_clusters)

# Use match to align clusters without changing order
manoj_p1$celltype <- best_celltypes$celltype[match(manoj_p1$cluster, best_celltypes$cluster)]
data@meta.data=manoj_p1

library(ggplot2)

library(ggsci)
library(scales)

# Get NEJM color palette (default 8 colors)
nejm_colors <- pal_lancet(alpha=0.6)(9)

# Show colors
#show_col(nejm_colors)








##############
#PISC=Peri-islet Schwann cells
#Gama=Gama delta T
new.cluster.ids <-best_celltypes$celltype

names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)


# Define unique colors for each cell type
celltype_colors <- c(
  "Alpha" = nejm_colors[1],
  "Beta" =  nejm_colors[2],
  "Delta" =  nejm_colors[3],
  "Gamma" =  nejm_colors[4],
  "Acinar" =  nejm_colors[5],
  "Ductal" =  nejm_colors[6],
  "Stellate" =  nejm_colors[7],
  "Unknown" = nejm_colors[8]
)





pdf("Pic.Manoj10.pdf", width = 5, height = 5)
p <- DimPlot(data, reduction = "umap", label = TRUE, alpha = 0.3) + 
  NoLegend() + 
  scale_color_manual(values = celltype_colors)

print(p)
dev.off()




jpeg("Pic.Manoj10.jpeg", width =5, height =5, units = "in", res = 600)
p <- DimPlot(data, reduction = "umap", label = TRUE, alpha = 0.3) + 
  NoLegend() + scale_color_nejm()
  #scale_color_manual(values = celltype_colors)

print(p)
dev.off()


jpeg("Pic.Manoj10.jpeg", width =5, height =5, units = "in", res = 600)
p <- DimPlot(data, reduction = "umap", label = TRUE, alpha = 0.3) + 
  NoLegend() + 
scale_color_manual(values = celltype_colors)

print(p)
dev.off()




jpeg("Pic.Manoj11.jpeg", width =10, height =5, units = "in", res = 600)
DimPlot(data, label = T , repel = T, label.size = 3,reduction = "umap", split.by = "Disease" ,alpha = 0.3) + NoLegend() + 
  scale_color_manual(values = celltype_colors)
dev.off()






pdf("Pic.Manoj11.pdf", width =10, height =5)
DimPlot(data, label = T , repel = T, label.size = 3,reduction = "umap", split.by = "Disease") + NoLegend() + 
  scale_color_manual(values = celltype_colors)
dev.off()



saveRDS(data,"Anno1.rds")






`