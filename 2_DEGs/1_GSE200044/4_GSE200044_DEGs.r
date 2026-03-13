
sva=as.data.frame(readRDS("iasva.res.allsva.rds")$sv)
head(sva)
scrna=readRDS("../2_DEG_Analysis/Anno1.rds")
pheno=scrna@meta.data
pheno1=cbind(pheno,sva)
scrna@meta.data=pheno1


celltype=unique(scrna$celltype)

#https://www.reddit.com/r/bioinformatics/comments/txed6i/correcting_for_batch_effects_in_differential/


i=1

dir.create("sva_scrna")

for(i in 1:length(celltype))
{
  scrna1=scrna[,scrna$celltype==celltype[i]]
  scrna1
  
  out=paste0("sva_scrna/",celltype[i],".withsva.rds")
  print(out)
  saveRDS(scrna1,out)
}


table(scrna1$celltype)
