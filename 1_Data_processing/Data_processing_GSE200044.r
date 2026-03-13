
library(scCustomize)
library(Signac)
library(Seurat)
library(GEOquery)
# Note that GSEMatrix=TRUE is the default

pheno=readRDS("GSE200044.pheno.rds")

pheno=dplyr::filter(pheno,`disease state:ch1`!="Pre-T2D")
pheno=pheno[-c(1:12),]


#saveRDS(pheno,"GSE200044.pheno.rds")


id=pheno$`sample id:ch1`
pancreas.list <-list()
i=1
for(i in 1:length(id))
{
  print(i)
  data=Read10X(paste0("10x_filtered_matrix/",id[i],"/"))
  metadata <-dplyr::filter(pheno,`sample id:ch1`==id[i])
  metadata <-dplyr::filter(metadata,`molecule_ch1`=="total RNA") #
  
 
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = data$`Gene Expression`, project = 'GSE200044')
  seurat_obj <- NormalizeData(seurat_obj, verbose = T)
  seurat_obj <- FindVariableFeatures(
    seurat_obj, selection.method = "vst",
    nfeatures = 2000, verbose =T
  )


  
  pbmc= seurat_obj
  pbmc@meta.data$Age=as.numeric(metadata$`age:ch1`)
  pbmc@meta.data$Gender=as.numeric(metadata$`gender:ch1`)
  pbmc@meta.data$Disease=as.numeric(metadata$`disease state:ch1`)
 
  outp=paste0("RNAseq/GSE200044_",id[i],"_2024.rds")
  pancreas.list[[i]]=pbmc
  
  saveRDS(pbmc , outp)
  
  
  rm(pbmc)
  
}
saveRDS(pancreas.list ,"GSE200044_1.rds")



merge_data=readRDS("GSE200044_1.rds")
merge_data <- Merge_Seurat_List(list_seurat = merge_data,add.cell.ids = id)
merge_data 
# Join the layers after merging
merge_data[["RNA"]] <- JoinLayers(merge_data[["RNA"]])
merge_data 
# Now you can access the combined counts matrix
pp <- merge_data[["RNA"]]$counts
merge_data
head(pp)
pheno=merge_data@meta.data
pheno$ID=rownames(pheno)
library(stringr)

# Extract the first word and create a new column
pheno$sample_ID <- str_extract(pheno$ID, "^[^_]+")

# View the updated data frame
head(pheno)

pheno2=readRDS("GSE200044.pheno.rds")

pheno2=dplyr::filter(pheno2,`disease state:ch1`!="Pre-T2D")
pheno2=pheno2[-c(1:12),]
pheno2$sample_ID=pheno2$`sample id:ch1`
pheno3=merge(pheno,pheno2, by="sample_ID")
pheno3$Gender=pheno3$`gender:ch1`
pheno3$Disease=pheno3$`disease state:ch1`
table(pheno3$Disease)
pheno3$Disease=gsub("Non-diabetic","Healthy",pheno3$Disease)
table(pheno3$Disease)
merge_data@meta.data=pheno3
saveRDS(merge_data,"GSE200044_RNA_2.rds")
saveRDS(pheno,"GSE200044_RNA_2.pheno.rds")







#library(Signac)
library(Seurat)
library(GEOquery)

pheno=readRDS("GSE200044.pheno.rds")

pheno$`sample id:ch1`
pheno$`disease state:ch1`
pheno$`gender:ch1`
pheno$`age:ch1`
library(dplyr)
pheno1=pheno%>%select(`sample id:ch1`,`disease state:ch1`,`gender:ch1`,`age:ch1`)

fil=list.files("RNAseq/", full.names = T)

# Add a new column to pheno1 with matching file names
pheno1$file_name <- sapply(pheno1$`sample id:ch1`, function(id) {
  fil[grep(paste0("_", id, "_"), fil)]
})

# Display the updated pheno1 dataframe
print(pheno1)






fil=c("RNAseq/GSE200044_C0026_2024.rds",
      "RNAseq/GSE200044_A0033_2024.rds",
      "RNAseq/GSE200044_A0031_2024.rds",
      "RNAseq/GSE200044_C0024_2024.rds",
      "RNAseq/GSE200044_A0024_2024.rds")




d1=readRDS(fil[1])
d2=readRDS(fil[2])
d3=readRDS(fil[3])
d4=readRDS(fil[4])
d5=readRDS(fil[5])



d1$Gender="M"
d1$Disease="Healthy"
d2$Gender="F"
d2$Disease="Healthy"
d3$Gender="M"
d3$Disease="T2D"
d4$Gender="F"
d4$Disease="T2D"
d5$Gender="M"
d5$Disease="T2D"




##############
comb1=merge(d1,d2)
comb1=merge(comb1,d3)
comb1=merge(comb1,d4)
comb1=merge(comb1,d5)


table(comb1@meta.data$Age)
table(comb1@meta.data$Disease)
table(comb1@meta.data$Gender)
saveRDS(comb1 ,"RNAseq_GSE200044.rds")
################

#hg19






# saveRDS(integrated ,"RNAseq_2024.rds")



