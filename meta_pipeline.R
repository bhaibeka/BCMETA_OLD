rm(list=ls())
ptm <- proc.time()
source(file.path("/stockage/homes/bachanp/Function/OpenDataset2.R"))
gselist <- OpenDataset2(config.fil="GSE_Curation.csv",rescale.checker=TRUE)
source(file.path("/stockage/homes/bachanp/Function/SubtypeClassification.R"))
STL <- SubtypeClassification(gselist=gselist, config.file="GSE_Curation.csv")
source(file.path("/stockage/homes/bachanp/Function/Merging.R"))
master.eset <- Merging(gselist=gselist, STL=STL, duplication.checker=TRUE)
proc.time() - ptm


##########################################################################
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
source(file.path("/stockage/homes/bachanp/Function/GseaAnalysis.R"))
matrix.subtype <- as.matrix(pData(master.eset))[,20:length(pData(master.eset))]
gsea.output <- GseaAnalysis(survival.score.list=survival.score.list,gmt.file="c2.all.v3.1.entrez.gmt",matrix.subtype=matrix.subtype)

proc.time() - ptm

