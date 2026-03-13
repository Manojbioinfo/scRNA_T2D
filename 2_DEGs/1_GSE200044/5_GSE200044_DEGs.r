# Assuming you have a Seurat object 'immune.combined'
library(Seurat)
library(Matrix)
library(ggsci)
library(ggplot2)
library(openxlsx)
# Set up your Seurat object and define clusters

all_files=list.files("sva_scrna/", full.names = T)
dir.create("DEGS")
i=1
for(i in 1:length(all_files))
{
  
  tryCatch(
    # This is what I want to do...
    {
      final_names=all_files[i]
      final_names2=gsub("sva_scrna","DEGS",final_names)
      
      
      seuratObj <- readRDS(final_names)
      
      cluster_names=gsub("sva_scrna//","",final_names)
      cluster_names=gsub(".withsva.rds","",cluster_names)
      
      
      Idents(seuratObj)=paste0(seuratObj$celltype,"_",seuratObj$Disease)
      
      
      cluster_cells_1=paste0(cluster_names,"_T2D")
      cluster_cells_2=paste0(cluster_names,"_Healthy")
      # Calculate your test statistic here
      # e.g., log-fold change, mean difference, etc.
      # Ensure you're using the same features for both clusters
      
      Idents(seuratObj)  # Returns the list above
      
      
      
      res=FindMarkers(seuratObj,latent.vars = c("age","gender",'SV1','SV2','SV3','SV4','SV5'),
                      test.use = "MAST", ident.1 = cluster_cells_1, ident.2 = cluster_cells_2, verbose = FALSE)
      
      
      res$p_val=ifelse(res$p_val==0,0.00001,res$p_val)
      
      ## change p value with 0 as 0.000000001 to avoid inf Z score as Zscore will be use for the META analysis
      
      res$Gene=rownames(res)
      
      
      
      # Calculate z-scores from p-values
      res$zscore <- qnorm(res$p_val / 2, lower.tail = FALSE) * sign(res$avg_log2FC)
      res$se <- res$avg_log2FC /res$zscore
      res$p_val_adj=p.adjust(res$p_val,method="fdr")
      
      
      
      Output=final_names2
      
      write.xlsx(res, Output)
    },
    # ... but if an error occurs, tell me what happened: 
    error=function(error_message) {
      message("This is my custom message.")
      message("And below is the error message from R:")
      message(error_message)
      
    }
  )
  
  
  
  
  
}

