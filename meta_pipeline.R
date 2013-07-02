rm(list=ls())
#################################################################
#####     Creation for preparing the different analysis     #####
#################################################################

#Opening the data from InsilicoDB
source(file.path("/stockage/homes/bachanp/Function/OpenDataset.R"))
gselist <- OpenDataset(config.file="GSE_Curation.csv")

#Predicting Subtype
source(file.path("/stockage/homes/bachanp/Function/SubtypeClassification.R"))
STL <- SubtypeClassification(gselist=gselist, config.file="GSE_Curation.csv")

#Merging every dataset into one master data call master.eset 
source(file.path("/stockage/homes/bachanp/Function/Merging.R"))
uberlist <- Merging(gselist=gselist, STL=STL, duplication.checker=FALSE,survdata="dfs", time.cens=10, method="unique")
master.eset <- uberlist[[1]]
GSCM <- uberlist[[2]]

#Make the survival analysis for every data set and every subtype
source(file.path("/stockage/homes/bachanp/Function/SurvivalAnalysis.R"))
survival.score.list <- SurvivalAnalysis(master.eset=master.eset)
#save(survival.score.list,file="All_SSL2")

#Create the GMT file for the GSEA analysis
source(file.path("/stockage/homes/bachanp/Function/WriteGMT.R"))
WriteGMT(input.file="EMT signature.csv",output.file="EMT.gmt")

##################################
######     GSEA analysis     #####
##################################
source(file.path("/stockage/homes/bachanp/Function/GseaAnalysis.R"))
gsea.output <- GseaAnalysis(survival.score.list=survival.score.list,gmt.file="Katrina.gmt",matrix.subtype=matrix.subtype)


#Extracting all prediction value base on subtype
index <- which(colnames(pData(master.eset))=="subtype")
for (i in 1:length(gsea.output)){
  assign(paste("prediction",colnames(pData(master.eset))[index+i],sep="_"), gsea.output[[i]][[1]])
}

################################################
#####      Random permutation analysis     ##### 
################################################

#rm(list=setdiff(ls(), "master.eset"))
n.perm <- 1000
source(file.path("/stockage/homes/bachanp/Function/RandomPermAnalysis.R"))
uberlist <- RandomPermAnalysis(master.eset=master.eset, inputfile="signatures_human_orthologs.csv", method="weighted.average", n.perm=n.perm)
metagene.ranking.list <- uberlist[[1]]
metagene.mat <- uberlist[[2]]
prognostic.list <- uberlist[[3]]

for (i in 1:length(prognostic.list)){
  assign(paste("prog.mat.",names(prognostic.list)[i],sep=""), prognostic.list[[i]])
}

####    Gaussian Type Random analysis graph    ####
pdf("LCC_signature.pdf")
z <- 1
for (i in 1:length(metagene.ranking.list[[z]])){
  g.output   <- qnorm(metagene.ranking.list[[z]][[i]][,3],mean=mean(metagene.ranking.list[[z]][[i]][,3]), sd=sd(metagene.ranking.list[[z]][[i]][,3]))
  histrv <- hist(g.output,20, main=sprintf("Signature score of metagene compare to a %g random permutations \nfor %s subtype",n.perm, names(metagene.ranking.list[[z]])[i]), xlab="qnorm of the signature score")
  abline(v=g.output[1], col = "red2", lwd=3)    
  op <- par(cex=.65)
  legend("topright", inset=c(0.5,0), legend=sprintf("%s \nqnorm=%g, metagene p.value=%g, \nsignificance p.value=%g", names(prognostic.list)[z], g.output[1], metagene.ranking.list[[z]][[i]][1,3], as.numeric(prognostic.list[[z]][i,2])), box.col="red2",text.font=2, y.intersp=2)
  par(op)
}
dev.off()

####    Boxplot type Random analysis graph    ####

pdf("RandomBoxplot.pdf")
for (i in 1:length(metagene.ranking.list)){
  Global <- metagene.ranking.list[[i]][[1]][,3]
  Lums <- metagene.ranking.list[[i]][[2]][,3]
  Basal <- metagene.ranking.list[[i]][[3]][,3]
  Her2 <- metagene.ranking.list[[i]][[4]][,3]
  LumB <- metagene.ranking.list[[i]][[5]][,3]
  LumA <- metagene.ranking.list[[i]][[6]][,3]
  Global.name <- sprintf("Global \np.value=%g, \nsignificance p.value=%g", metagene.ranking.list[[i]][[1]][1,3], as.numeric(prognostic.list[[i]][1,2]))
  Lums.name <- sprintf("Lums \np.value=%g, \nsignificance p.value=%g", metagene.ranking.list[[i]][[2]][1,3], as.numeric(prognostic.list[[i]][2,2]))
  Basal.name <- sprintf("Basal \np.value=%g, \nsignificance p.value=%g", metagene.ranking.list[[i]][[3]][1,3], as.numeric(prognostic.list[[i]][3,2]))
  Her2.name <- sprintf("Her2 \np.value=%g, \nsignificance p.value=%g", metagene.ranking.list[[i]][[4]][1,3], as.numeric(prognostic.list[[i]][4,2]))
  LumB.name <- sprintf("LumB \np.value=%g, \nsignificance p.value=%g", metagene.ranking.list[[i]][[5]][1,3], as.numeric(prognostic.list[[i]][5,2]))
  LumA.name <- sprintf("LumA \np.value=%g, \nsignificance p.value=%g", metagene.ranking.list[[i]][[6]][1,3], as.numeric(prognostic.list[[i]][6,2]))
  op <- par(cex.axis=0.65)
  boxplot(Global, Lums, Basal, Her2, LumB, LumA, names=c(Global.name, Lums.name, Basal.name, Her2.name, LumB.name, LumA.name),xlab="Subtype", ylab="p.value score",main=sprintf("Random analysis of %g for \n%s gene",n.perm, names(metagene.ranking.list)[i]))
  par(op)
  segments(x0=0.55,y0=metagene.ranking.list[[i]][[1]][1,3],x1=1.45,y1=metagene.ranking.list[[i]][[1]][1,3],col="red")
  segments(x0=1.55,y0=metagene.ranking.list[[i]][[2]][1,3],x1=2.45,y1=metagene.ranking.list[[i]][[2]][1,3],col="red")
  segments(x0=2.55,y0=metagene.ranking.list[[i]][[3]][1,3],x1=3.45,y1=metagene.ranking.list[[i]][[3]][1,3],col="red")
  segments(x0=3.55,y0=metagene.ranking.list[[i]][[4]][1,3],x1=4.45,y1=metagene.ranking.list[[i]][[4]][1,3],col="red")
  segments(x0=4.55,y0=metagene.ranking.list[[i]][[5]][1,3],x1=5.45,y1=metagene.ranking.list[[i]][[5]][1,3],col="red")
  segments(x0=5.55,y0=metagene.ranking.list[[i]][[6]][1,3],x1=6.45,y1=metagene.ranking.list[[i]][[6]][1,3],col="red")
}
dev.off()

################################################
#####      Subtype dependence analysis     ##### 
################################################

#rm(list=setdiff(ls(), "master.eset"))
source(file.path("/stockage/homes/bachanp/Function/SubtypeDependanceAnalysis.R"))
source(file.path("/stockage/homes/bachanp/Function/SubtypeBoxPlot.R"))
 
uberlist <- SubtypeDependanceAnalysis(master.eset=master.eset,inputfile="signatures_human_orthologs.csv", method="weighted.average")
dependence <- uberlist[[1]]
modified.master.eset <- uberlist[[2]]
Kruskal.test <- dependence[[1]]
Wilcox.test <- dependence[[2]]
median.matrix <- dependence[[3]]
SubtypeBoxPlot(master.eset=master.eset,Wilcox.test=Wilcox.test,modified.master.eset=modified.master.eset)





