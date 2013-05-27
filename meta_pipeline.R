rm(list=ls())
ptm <- proc.time()
source(file.path("/stockage/homes/bachanp/OpenDataset.R"))
uberlist <- OpenDataset(config.fil="GSE_Curation.csv",duplication.checker=TRUE)
gselist <- uberlist[[1]]
master.matrix <- uberlist[[2]]
source(file.path("/stockage/homes/bachanp/SubtypeClassification2.R"))
matrix.subtype <- SubtypeClassification2(matrix.exprs=master.matrix, config.file="GSE_Curation.csv")

#This shall be replace when the survivalscore function will be operational
##########################################################################
survival.score.list <- list()
gene.id <- sub("geneid_","",rownames(master.matrix))
for (i in 1:4){
  ranking <- cbind(gene.id[floor(runif(190, 1,101))], runif(190,-10,10))  
  colnames(ranking) <- c("EntrezGene.ID","Survival score")  
  survival.score.list[[i]] <- ranking
}
##########################################################################
source(file.path("/stockage/homes/bachanp/GseaAnalysis.R"))
gsea.output <- GseaAnalysis(survival.score.list=survival.score.list,gmt.file="c2.all.v3.1.entrez.gmt",matrix.subtype=matrix.subtype)

proc.time() - ptm

