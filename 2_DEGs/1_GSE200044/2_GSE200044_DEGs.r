library(SummarizedExperiment)
library(matrixStats)
library(Seurat)
set.seed(100)
#color.vec <- brewer.pal(3, "Set1")

rm(list = ls())
#https://www.bioconductor.org/packages/devel/bioc/vignettes/iasva/inst/doc/detecting_hidden_heterogeneity_iasvaV0.95.html
seurat_obj_all <- readRDS("../2_DEG_Analysis/Anno1.rds")
pheno=seurat_obj_all@meta.data


norm_matrix <- seurat_obj_all[["RNA"]]$data




# 2. Compute variance for each gene (row)
gene_variances <- matrixStats::rowVars(as.matrix(norm_matrix))  # requires matrixStats package

# 3. Calculate number of top 10% genes
top_n <- ceiling(nrow(norm_matrix) * 0.10)

# 4. Get indices of top 10% most variable genes
top_var_genes_idx <- order(gene_variances, decreasing = TRUE)[1:top_n]

# 5. Get gene names (rownames) for the top variable genes
top_var_gene_names <- rownames(norm_matrix)[top_var_genes_idx]

# Optional: Subset the expression matrix
top_var_matrix <- norm_matrix[top_var_gene_names, ]



#norm_matrix=top_var_matrix 

counts <- top_var_matrix  # or norm_matrix / scaled_matrix

geo_lib_size <- colSums(log(counts + 1))
#barplot(geo_lib_size, xlab = "Cell", ylab = "Geometric Lib Size", las = 2)


lcounts <- log(counts + 1)

# PC1 and Geometric library size correlation
#pc1 <- irlba(lcounts - rowMeans(lcounts), 1)$v[, 1]
#cor(geo_lib_size, pc1)

set.seed(100)
pheno$patient_id=rownames(pheno)
# Keep only the first word before the first underscore
pheno$patient_id <- sub("_.*", "", pheno$patient_id)
pheno$patient_id[1:10]


print(table(pheno$Disease))

disease <- pheno$Disease
age=pheno$Age
age=gsub("year","",age)
age=as.numeric(age)
table(age)
gender=pheno$Gender
gender=ifelse(gender=="male","M","F")


length(disease)
length(age)
length(gender)
length(geo_lib_size)

gender=as.factor(gender)
#mod <- model.matrix(~disease+ geo_lib_size+age)


pheno$geo_lib_size <- colSums(log(counts + 1))

# Clean age
pheno$age <- as.numeric(pheno$age)

# Force factors and drop unused levels
pheno$disease <- factor(pheno$disease)
pheno$disease <- droplevels(pheno$disease)

pheno$gender <- factor(pheno$gender)
pheno$gender <- droplevels(pheno$gender)

# Double-check if factors have at least 2 levels
if (nlevels(pheno$disease) < 2) stop("`disease` has only one level.")
if (nlevels(pheno$gender) < 2) stop("`gender` has only one level.")


print(dim(pheno))
print(dim(counts))


pheno$disease=ifelse(pheno$disease=="Non-diabetic","Healthy",pheno$disease)
pheno$disease=ifelse(pheno$disease=="2","T2D",pheno$disease)
table(pheno$disease)
saveRDS(pheno,"pheno.sva.rds")
saveRDS(counts,"counts.sva.rds")


pheno=readRDS("pheno.sva.rds")

table(pheno$disease)
