RandomPermAnalysis <- function(master.eset,inputfile,method,n.perm){
  #########################################
  #####   Creation of the metagene   ######
  #########################################
  
  source(file.path("/stockage/homes/bachanp/Function/MetaGene.R"))
  #Fixe variable
  data <- t(exprs(master.eset))
  annot <- Biobase::featureData(master.eset)  
  n.perm <- n.perm
  
  #Geneset sample gene first iteration
  gene.mat <- read.m.file(inputfile)[[1]]
  method <- method
  sample.genes <- paste("geneid.", gene.mat$EntrezGene.ID, sep="")
  sample.coef <- gene.mat$coefficient
  
  sample.metagene <- MetaGene(genes=sample.genes, data=data, annot=annot, coefficients=sample.coef, method=method)[[1]]
  
  #Random sample gene
  randgene.length <- length(intersect(sample.genes, colnames(data)))
  rand.genes <- sample.genes[1:randgene.length]
  
  library(parallel)
  
  #Setting every parallelizing parameter
  nbcore <- 8
  availcore <- detectCores()
  if(nbcore > availcore) { nbcore <- availcore }
  options("mc.cores"=nbcore)
  splitix <- splitIndices(nx=n.perm, ncl=nbcore)
  
  #Metagene core function
  rand.metagene <- mclapply(splitix, function(splitix2,...){  
    sub.rand.metagene <- lapply(splitix2, function(x){
      rand.coef <- runif(randgene.length,min(sample.coef),max(sample.coef))
      subsub.rand.metagene <- MetaGene(genes=rand.genes, data=data, annot=annot, coefficients=rand.coef, method=method)[[1]]
      return(subsub.rand.metagene)  
    })
    sub.rand.metagene <- do.call(rbind, sub.rand.metagene)
  })
  rand.metagene <- do.call(rbind, rand.metagene)
  colnames(rand.metagene) <- colnames(exprs(master.eset))
  
  metagene.mat <- rbind(sample.metagene, rand.metagene)
  ##############################################
  #####   Calculating Concordance Index   ######
  ##############################################
  
  index <- which(colnames(pData(master.eset))=="subtype")
  subtype.name <- colnames(pData(master.eset))[(index+1):ncol(pData(master.eset))]
  pheno.data <- pData(master.eset)
  exprs.matrix <- exprs(master.eset)           
  stime <- as.numeric(as.vector(pheno.data$surv.time))
  sevent <- as.numeric(as.vector(pheno.data$surv.event))
  index1 <- which(colnames(pheno.data)=="subtype")
  strat <- as.vector(pData(master.eset)$dataset)
  metagene.ranking.list <- list()
  
  #Setting every parallelizing parameter
  nbcore <- 8
  availcore <- detectCores()
  if(nbcore > availcore) { nbcore <- availcore }
  options("mc.cores"=nbcore)
  splitix <- splitIndices(nx=nrow(metagene.mat), ncl=nbcore)
  source(file.path("/stockage/homes/bachanp/Function/ConcordanceIndex3.R"))
  
  #Concordance index core function
  for (j in 1:length(subtype.name)){  
    weights <- as.numeric(as.vector(pheno.data[,index1+j]))
    metagene.ranking <- mclapply(splitix, function(x){
      sub.metagene.ranking <- ConcordanceIndex3(exprs.matrix=metagene.mat[x, , drop=FALSE], stime=stime, sevent=sevent, strat=strat, weights=weights)  
      return(sub.metagene.ranking)
    })
    metagene.ranking <- do.call(rbind, metagene.ranking)
    metagene.ranking.list[[j]] <- metagene.ranking
  }
  return(list("metagena.ranking.list"=metagene.ranking.list, "metagene.mat"=metagene.mat))
}
  
  


