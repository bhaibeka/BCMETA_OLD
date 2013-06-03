rm(list=ls())
ptm <- proc.time()
source(file.path("/stockage/homes/bachanp/Function/OpenDataset2.R"))
gselist <- OpenDataset2(config.fil="GSE_Curation.csv",rescale.checker=FALSE)
source(file.path("/stockage/homes/bachanp/Function/SubtypeClassification.R"))
STL <- SubtypeClassification(gselist=gselist, config.file="GSE_Curation.csv")
source(file.path("/stockage/homes/bachanp/Function/Merging.R"))
master.eset <- Merging(gselist=gselist, STL=STL, duplication.checker=TRUE,checker=TRUE)
source(file.path("/stockage/homes/bachanp/Function/SurvivalAnalysis2.R"))
survival.score.list <- SurvivalAnalysis2(master.eset=master.eset)
source(file.path("/stockage/homes/bachanp/Function/GseaAnalysis.R"))
#matrix.subtype <- as.matrix(pData(master.eset))[,22:length(pData(master.eset))]
gsea.output <- GseaAnalysis(survival.score.list=survival.score.list,gmt.file="c2.all.v3.1.entrez.gmt",matrix.subtype=matrix.subtype)

proc.time() - ptm

