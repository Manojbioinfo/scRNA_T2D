

library(reticulate)
# will only work from terminal in  py3.8 env
#



# 5. Now you can use Python functions or reticulate integration safely


# Load required libraries
library(Seurat)
library(tidyverse)
#install.packages("DescTools", dependencies = T)
library(irlba) 
library(iasva, lib.loc = "../../../../../../Rpackage/myRlib")
library(sva)
library(Rtsne)
library(pheatmap, lib.loc = "../../../../../../Rpackage/myRlib")
library(corrplot)
library(DescTools, lib.loc = "../../../../../../Rpackage/myRlib")
library(RColorBrewer)
library(SummarizedExperiment)



## we have to run in base in local machine
set.seed(100)
color.vec <- brewer.pal(3, "Set1")

rm(list = ls())
#https://www.bioconductor.org/packages/devel/bioc/vignettes/iasva/inst/doc/detecting_hidden_heterogeneity_iasvaV0.95.html

counts=readRDS("counts.sva.rds")
pheno=readRDS("pheno.sva.rds")

#pheno$patient_id is cell identifier
pheno$patient_id <- sapply(strsplit(pheno$patient_id, "-"), `[`, 1)

mod <- model.matrix(~ disease +patient_id+ gender + geo_lib_size + age, data = pheno)

print(head(mod))

counts=as.matrix(counts)
print(class(counts))

# create a summarizedexperiment class
summ_exp <- SummarizedExperiment(assays = counts)
iasva.res<- iasva(summ_exp, mod[, -1],verbose = FALSE, 
                  permute = FALSE, num.sv = 5)


saveRDS(iasva.res,"iasva.res.allsva.rds")

