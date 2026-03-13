# Load Libraries
library(Seurat)
library(hdWGCNA)
library(tidyverse)
library(patchwork)
library(openxlsx)
library(WGCNA)
set.seed(134)
rm(list=ls())

# Save plot function
save_plots <- function(plot_obj, filename_base) {
  pdf(paste0(filename_base, ".pdf"), width = 8, height = 6)
  print(plot_obj)
  dev.off()
  
  jpeg(paste0(filename_base, ".jpeg"), width = 8 * 400, height = 6 * 400, res = 400)
  print(plot_obj)
  dev.off()
}


# Save plot function
save_plot_kme <- function(plot_obj, filename_base) {
  pdf(paste0(filename_base, ".pdf"), width =12, height = 6)
  print(plot_obj)
  dev.off()
  
  jpeg(paste0(filename_base, ".jpeg"), width = 12 * 400, height = 6 * 400, res = 400)
  print(plot_obj)
  dev.off()
}

# Save plot function
save_plot_feat<- function(plot_obj, filename_base) {
  pdf(paste0(filename_base, ".pdf"), width = 12, height = 5)
  print(plot_obj)
  dev.off()
  
  jpeg(paste0(filename_base, ".jpeg"), width = 12 * 400, height = 5 * 400, res = 400)
  print(plot_obj)
  dev.off()
}


# Load annotation and Seurat object
data <- readRDS("../2_DEG_Analysis/Anno1.rds")
data=data[,data$celltype!="Unknown"]
table(data$celltype)

# Initialize for WGCNA
seurat_obj <- SetupForWGCNA(
  seurat_obj = data,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "tutorial"
)

# Construct metacells
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  #group.by = c("cell_type", "Sample"),
  reduction = 'umap',
  k = 25,
  max_shared = 10,
  #ident.group = 'cell_type'
)

# Normalize
seurat_obj <- NormalizeMetacells(seurat_obj)

# Backup full object
seurat_obj_clean <- seurat_obj



##################
seurat_obj <- seurat_obj_clean

# Set expression data for current cluster
seurat_obj <- SetDatExpr(
  seurat_obj,
 # group_name = cluster_id,
  #group.by = 'cell_type',
  assay = 'RNA',
  layer = 'data'
)


#######################
cluster_id="Complete"

# Test soft powers
seurat_obj <- TestSoftPowers(seurat_obj, networkType = 'unsigned') ### WCGNA original packages suggested unsigned. IN unsigned they look for both positive and negative genes that are clustered togther . In signed only positive genes are clustered togther and. negative cluster speartely togther

plot_list <- PlotSoftPowers(seurat_obj)
save_plots(wrap_plots(plot_list, ncol = 2), "Pic_1_SoftPowerPlot")

# Construct co-expression network
seurat_obj <- ConstructNetwork(seurat_obj, tom_name = cluster_id)
save_plots(PlotDendrogram(seurat_obj, main = paste0(cluster_id, " Dendrogram")), "Pic_2_Dendrogram")

############################

# Compute module eigengenes
#seurat_obj <- ModuleEigengenes(seurat_obj, group.by.vars = "Sample")
seurat_obj <- ModuleEigengenes(seurat_obj)
# Connectivity & KMEs
#seurat_obj <- ModuleConnectivity(seurat_obj, group.by = "Sample", group_name = "Dataset")

seurat_obj <- ModuleConnectivity(seurat_obj)


p_kme <- PlotKMEs(seurat_obj, ncol = 6)
save_plot_kme(p_kme, "Pic_3_KME_Plot")



# Get modules and hub genes
modules <- GetModules(seurat_obj) %>% filter(module != "grey")
write.xlsx(modules, "Table1_modules.xlsx", rowNames = FALSE)

hub_genes <- GetHubGenes(seurat_obj, n_hubs = 50)
write.xlsx(hub_genes, "Table2_hub_genes.xlsx", rowNames = FALSE)

# Save feature plots
feature_plots <- ModuleFeaturePlot(seurat_obj, plot_ratio = 2, features = "hMEs", order = TRUE)
save_plot_feat(wrap_plots(feature_plots, ncol = 6), "Pic_4_Module_hMEs_FeaturePlot")

# Save final Seurat object
saveRDS(seurat_obj, file = paste0(cluster_id, "_hdWGCNA.rds"))
############################
############################
############################


