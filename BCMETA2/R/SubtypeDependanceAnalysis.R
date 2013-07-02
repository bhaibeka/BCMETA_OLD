SubtypeDependanceAnalysis <-
function(master.eset,gene.list,inputfile,method){
  # Classify GSE subtype and give out the probability of being that subtype
  #
  # Args:
  #   master.eset: A master.eset containing everything
  #   gene.list : a list of gene which are to be analyze, gene.list= "ALL" if all gene in the master.eset
  #               are to be consider
  #   method: that the user want to create the metagene
  #   inputfile: The inputfile to get the gene list and make the gene mapping
  #
  # Returns:
  #   The dependance of a specific gene over the subtype
  
  source(file.path("/stockage/homes/bachanp/Function/SubtypeDependance.R"))
  source(file.path("/stockage/homes/bachanp/Function/MetaGene.R"))
  
  metagene.matrix <- NULL
  gene.mat.all <- read.m.file(inputfile)
  for (i in 1:length(gene.mat.all)){
    gene.mat <- gene.mat.all[[i]]
    #total.gene.list <- paste("geneid.",read.m.file("signatures_human_orthologs.csv")[[1]]$EntrezGene.ID,sep="")
    data <- t(exprs(master.eset))
    annot <- Biobase::featureData(master.eset)  
    sample.genes <- paste("geneid.", gene.mat$EntrezGene.ID, sep="")
    sample.coef <- gene.mat$coefficient
    method <- "weighted.average" 
    temp <- master.eset
    if (length(sample.genes) !=1){
      sample.metagene <- MetaGene(genes=sample.genes, data=data, annot=annot, coefficients=sample.coef, method=method)[[1]]
    }else{
      sample.metagene <- MetaGene(genes=sample.genes, data=data, annot=annot, coefficients=sample.coef, method=method)      
    }
    metagene.matrix <- rbind(metagene.matrix, sample.metagene)
  }
  rownames(metagene.matrix) <- names(gene.mat.all)
  exprs(temp) <- rbind(metagene.matrix, exprs(temp))
  gene.list <- names(gene.mat.all)
  dep.output <- SubtypeDependance(master.eset=temp,gene.list=gene.list)
  return(list("dep.output"=dep.output,"temp"=temp))
}
