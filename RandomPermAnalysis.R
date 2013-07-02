RandomPermAnalysis <- function(master.eset, inputfile, method, n.perm){
  
  #########################################
  #####   Creation of the metagene   ######
  #########################################
  source(file.path("/stockage/homes/bachanp/Function/MetaGene.R"))
  metagene.rank.masterlist <- metagene.mat.list <- prognostic.list <- list()  
  #Fixe variable
  data <- t(exprs(master.eset))
  annot <- Biobase::featureData(master.eset)  
  #Geneset sample gene first iteration
  gene.mat.all <- read.m.file(inputfile)
  #gene.mat <- read.m.file(file.path("/home/bachanp/GMT and Signature file/signatures_human_orthologs.csv"))[[1]]
  #gene.mat.all <- read.m.file("modules_benchmark.csv")
  for (k in 1:length(gene.mat.all)){     
    gene.mat <- gene.mat.all[[k]]
    sample.genes <- paste("geneid.", gene.mat$EntrezGene.ID, sep="")
    sample.coef <- gene.mat$coefficient
    
    if (length(sample.genes) !=1){
      sample.metagene <- MetaGene(genes=sample.genes, data=data, annot=annot, coefficients=sample.coef, method=method)[[1]]
    }else{
      sample.metagene <- MetaGene(genes=sample.genes, data=data, annot=annot, coefficients=sample.coef, method=method)      
    }
      
    #Random sample gene
    randgene.length <- length(intersect(sample.genes, colnames(data)))    
    
    library(parallel)
    
    #Setting every parallelizing parameter
    nbcore <- 8
    availcore <- detectCores()
    if(nbcore > availcore) { nbcore <- availcore }
    options("mc.cores"=nbcore)
    splitix <- splitIndices(nx=n.perm, ncl=nbcore)
    
    #Metagene core function
    if (randgene.length !=1){
      rand.metagene <- mclapply(splitix, function(splitix2,...){  
        sub.rand.metagene <- lapply(splitix2, function(x){
            rand.genes <- rownames(exprs(master.eset))[sample(1:nrow(exprs(master.eset)),randgene.length,replace=F)]
            subsub.rand.metagene <- MetaGene(genes=rand.genes, data=data, annot=annot, coefficients=sample.coef, method=method)[[1]]
          return(subsub.rand.metagene)  
        })
        sub.rand.metagene <- do.call(rbind, sub.rand.metagene)
      })
      rand.metagene <- do.call(rbind, rand.metagene)
      colnames(rand.metagene) <- colnames(exprs(master.eset))
    
    }else{  
      
      index <- sample(1:nrow(exprs(master.eset)),n.perm,replace=F)
      if (any(index==which(rownames(exprs(master.eset))==sample.genes))){
        index<- index[-which(rownames(exprs(master.eset))==sample.genes)]
        index <- c(index,1)
        }    
      rand.metagene <- exprs(master.eset)[index,]
      }  
    
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
    names(metagene.ranking.list) <- c("Global", "Lums", "Basal", "Her2", "LumB", "LumA")
    
    prog.mat <- NULL
    for (i in 1:length(metagene.ranking.list)){
      prog.mat[i] <- sum(metagene.ranking.list[[i]][-1, "p"] <= metagene.ranking.list[[i]][1, "p"]) / n.perm
    }  
    prog.mat <- cbind(names(metagene.ranking.list), prog.mat)
    metagene.mat.list[[k]] <- metagene.mat
    prognostic.list[[k]]<- prog.mat    
    metagene.rank.masterlist[[k]] <- metagene.ranking.list
  }
  names(prognostic.list) <- names(metagene.rank.masterlist) <- names(metagene.mat.list) <- names(gene.mat.all)[1:k]
  return(list("metagene.rank.masterlist"=metagene.rank.masterlist, "metagene.mat.list"=metagene.mat.list,"prognostic.list"=prognostic.list))
}
  

