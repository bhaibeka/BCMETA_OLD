rm(list=ls())
ptm <- proc.time()
#Opening the data from InsilicoDB
source(file.path("/stockage/homes/bachanp/BCMETA/OpenDataset2.R"))
gselist <- OpenDataset2(config.fil="GSE_Curation.csv",rescale=FALSE)

#Predicting Subtype
source(file.path("/stockage/homes/bachanp/Function/SubtypeClassification.R"))
STL <- SubtypeClassification(gselist=gselist, config.file="GSE_Curation.csv")

#Merging every dataset into one master data call master.eset
source(file.path("/stockage/homes/bachanp/BCMETA/Merging.R"))
master.eset <- Merging(gselist=gselist, STL=STL, duplication.checker=TRUE,survdata="dfs", time.cens=10)

#Make the survival analysis for every data set and every subtype
source(file.path("/stockage/homes/bachanp/Function/SurvivalAnalysis2.R"))
survival.score.list <- SurvivalAnalysis2(master.eset=master.eset)

#Create the GMT file for the GSEA analysis
source(file.path("/stockage/homes/bachanp/Function/WriteGMT.R"))
WriteGMT(input.file="signatures_human_orthologs.csv",output.file="Katrina.gmt")

#Make the GSEA analysis
source(file.path("/stockage/homes/bachanp/Function/GseaAnalysis.R"))
gsea.output <- GseaAnalysis(survival.score.list=survival.score.list,gmt.file="Katrina.gmt",matrix.subtype=matrix.subtype)

proc.time() - ptm

#Extracting all prediction value base on subtype
index <- which(colnames(pData(master.eset))=="subtype")
for (i in 1:length(gsea.output)){
  assign(paste("prediction",colnames(pData(master.eset))[index+i],sep="_"), gsea.output[[i]][[1]])
}

K.probe <- as.character(as.vector(read.m.file("signatures_human_orthologs.csv")[[1]][,2]))
all.probe <- as.vector(Biobase::featureData(master.eset)@data[,1])
index <- match(K.probe, all.probe)
index <- index[!is.na(index)]
K.score <- rep(NA,length(all.probe))
K.gene.score <- as.vector(read.m.file("signatures_human_orthologs.csv")[[1]][,3])
for (i in 1:length(index)){
  if(!is.na(index[i])){
    K.score[index[i]] <- K.gene.score[i]
  }
}

gene.score <- -log10(survival.score.list[[1]][,3])*sign(survival.score.list[[1]][,1]-0.5)
gene.score <- gene.score/max(abs(gene.score))
gene.score <- sort(gene.score, decreasing=TRUE)
gene.score <- as.matrix(gene.score)
K.score <- K.score/max(abs(K.gene.score))
#K.score <- sort(K.score, decreasing=TRUE)
K.score <- as.matrix(K.score)
heatmap.matrix <- cbind(gene.score,K.score)

library(gplots)

blueyellow <- function(N) {
  x <- (1:N)-1
  rgb(x, x, rev(x), maxColorValue=N)
}

#heatmap.2(heatmap.matrix, col=blueyellow(length(gene.score)),Rowv=FALSE,Colv=FALSE)
heatmap(heatmap.matrix, col=blueyellow(length(gene.score)),Rowv=FALSE,Colv=FALSE)
image(heatmap.matrix, col=blueyellow(length(gene.score)), zlim=c(-1,1))

  