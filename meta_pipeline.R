rm(list=ls())
ptm <- proc.time()
source(file.path("/stockage/homes/bachanp/OpenDataset2.R"))
gselist <- OpenDataset2(config.fil="GSE_Curation.csv")
source(file.path("/stockage/homes/bachanp/SubtypeClassification.R"))
STL <- SubtypeClassification(gselist=gselist, config.file="GSE_Curation.csv")
source(file.path("/stockage/homes/bachanp/Merging.R"))
master.eset <- Merging(gselist=gselist, STL=STL, duplication.checker=TRUE)

#This shall be replace when the survivalscore function will be operational
##########################################################################
survival.score.list <- list()
gene.id <- sub("geneid_","",rownames(master.matrix))
for (i in 1:4){
  #ranking <- cbind(gene.id[floor(runif(190, 1,101))], runif(190,-10,10))  
  ranking <- cbind(gene.id, runif(nrow(master.matrix),-10,10))   
  colnames(ranking) <- c("EntrezGene.ID","Survival score")  
  survival.score.list[[i]] <- ranking
}
##########################################################################
source(file.path("/stockage/homes/bachanp/GseaAnalysis.R"))
gsea.output <- GseaAnalysis(survival.score.list=survival.score.list,gmt.file="c2.all.v3.1.entrez.gmt",matrix.subtype=matrix.subtype)

proc.time() - ptm

