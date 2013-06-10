OpenDataset2 <- function(config.file, rescale=FALSE){
  #Open every GSE, specify in the configuration CSV, from InsilicoDB and create a list of Eset
  #structure call gselist
  #
  # Args:
  #   config.file: The string name of the configuration CSV file.
  #   rescale.checker: A marker, either TRUE or FALSE if you you want to rescale your dataset to be all
  #                        on the same bimodality
  #   
  # Returns:     
  #     The GSE list of the specify GSE number in the configuraiton CSV file 
  
  
  require(inSilicoDb2)
  require(genefu)
  #login
  InSilicoLogin(login="bhaibeka@gmail.com", password="747779bec8a754b91076d6cc1f700831")
  #File.To.Open (FTO)
  FTO <- read.csv(config.file)   
  gselist <- list() 
  load("~/NKIORIGINALPROBE.RData")
  load("~/METABRICORIGINALPROBE.RData")
  gselist[[1]] <- METABRICORIGINALPROBE
  phenoData(gselist[[1]])@data <- cbind(phenoData(gselist[[1]])@data, rep("GPL6947",nrow(phenoData(gselist[[1]])@data)))
  colnames(phenoData(gselist[[1]])@data)[length(colnames(phenoData(gselist[[1]])@data))] <- "platform"
  gselist[[2]]<- NKIORIGINALPROBE
  colnames(phenoData(gselist[[2]])@data)[length(colnames(phenoData(gselist[[2]])@data))] <- "platform"
  phenoData(gselist[[2]])@data[length(phenoData(gselist[[2]])@data)] <- "AGIR"
  counter <- 3
  
  #Opening all dataset from InsilicoDB
  for (i in 16:nrow(FTO)){
    if (as.character(FTO[i,4])==TRUE){    
      GPL <- getPlatforms(dataset=as.character(FTO[i,2]))[1]
      gselist[[counter]] <- getDataset(dataset=as.character(FTO[i,2]), curation=FTO[i,3], platform=GPL)        
      counter <- counter+1
    }
  }
  #Making sure that every rownames for matrix expression is well geneid.#####
  #and that every column of the clinical data matrix is well aligned for every GSE dataset
  for (i in 1:length(gselist)){
    rownames(exprs(gselist[[i]])) <- paste("geneid.",as.vector(Biobase::featureData(gselist[[i]])$ENTREZID),sep="")
    phenoData(gselist[[i]])@data<- phenoData(gselist[[i]])@data[,match(colnames(phenoData(gselist[[1]])@data),colnames(phenoData(gselist[[i]])@data))]    
  }
  #rescale each dataset individually     
  if (rescale) {
    for (i in 1:length(gselist)) {
      exprs(gselist[[i]]) <- t(((apply(exprs(gselist[[i]]), 1, function(x) { return(genefu::rescale(x, q=0.05, na.rm=TRUE)) })) - 0.5) * 2)
    }
  }
  #logout
  InSilicoLogout() 
  return(gselist)
}