OpenDataset <-
function(config.file){
  #Open every GSE, specify in the configuration CSV, from InsilicoDB and create a list of Eset
  #structure call gselist
  #
  # Args:
  #   config.file: The string name of the configuration CSV file.
  #   
  # Returns:     
  #     The GSE list of the specify GSE number in the configuraiton CSV file 
  
  
  require(inSilicoDb2)
  require(genefu)
  require(jetset)  
  source(file.path("/stockage/homes/bachanp/Function/ProbeGeneMapping.R"))
  source(file.path("/stockage/homes/bachanp/Function/PlatformMerging.R"))
  
  #login
  InSilicoLogin(login="bhaibeka@gmail.com", password="747779bec8a754b91076d6cc1f700831")
  #File.To.Open (FTO)
  FTO <- read.csv(config.file)   
  #FTO <- read.csv("GSE_Curation.csv")
  gselist <- list() 
  counter <- 1
  routine.hack1 <- NULL
  
  #Opening all dataset from InsilicoDB
  for (i in 16:nrow(FTO)){
    if (as.character(FTO[i,4])==TRUE){  
      if (as.character(FTO[i,2])=="GSE2034" || as.character(FTO[i,2])=="GSE5327"){
        routine.hack1 <- rbind(routine.hack1, c(as.character(FTO[i,2]), counter))
      }
      GPL <- getPlatforms(dataset=as.character(FTO[i,2]))
      temp <- getDatasets(dataset=as.character(FTO[i,2]), norm=as.character(FTO[i,6]), curation=FTO[i,3], features="PROBE")
      for (j in 1:length(GPL)){
        gselist[[counter]] <- temp[[j]]
        #making probe gene mapping for Affy structure
        if (as.character(FTO[i,5])=="Affy"){gselist[[counter]] <- ProbeGeneMapping(gselist[[counter]])}        
        else{
          rownames(exprs(gselist[[counter]])) <- paste("geneid.",Biobase::featureData(gselist[[counter]])$ENTREZID, sep="")
          index <- match(pData(gselist[[counter]])$id, colnames(exprs(gselist[[counter]])))
          exprs(gselist[[counter]]) <- exprs(gselist[[counter]])[,index]
          rownames(pData(gselist[[counter]])) <- as.character(pData(gselist[[counter]])$id)
        }        
        counter <- counter +1        
      }
      #Merging dataset if they're on different platforms
      if (length(GPL) > 1){
        index <- length(gselist)-(length(GPL)-1)
        gselist[[index]] <- PlatformMerging(gselist=gselist, GPL.length=length(GPL))
        for (k in 1:length(gselist)){gselist[[index+k]] <- NULL}
        counter <- counter - (length(GPL)-1)
      }
    } 
  }
  
  if (!is.null(routine.hack1)){
    if (nrow(routine.hack1)==2){
      exprs(gselist[[as.numeric(routine.hack1[1,2])]]) <- cbind(exprs(gselist[[as.numeric(routine.hack1[1,2])]]), exprs(gselist[[as.numeric(routine.hack1[2,2])]]))
      phenoData(gselist[[as.numeric(routine.hack1[1,2])]])@data <- rbind(phenoData(gselist[[as.numeric(routine.hack1[1,2])]])@data, phenoData(gselist[[as.numeric(routine.hack1[2,2])]])@data)
      gselist[[as.numeric(routine.hack1[2,2])]] <- NULL
    }
  }
  
  #Making sure that every column of the clinical data matrix is well aligned for every GSE dataset
  for (i in 1:length(gselist)){ 
    index=match(colnames(exprs(gselist[[i]])),rownames(pData(gselist[[i]])))    
    phenoData(gselist[[i]])@data <- phenoData(gselist[[i]])@data[index, ]
    
    #take only the breast cancer data if GSE containt more than 1 type of cancer
    if (any(colnames(pData(gselist[[i]]))=="Anatomical site")){
      index <- which(pData(gselist[[i]])$Anatomical=="breast")
      phenoData(gselist[[i]])@data <- phenoData(gselist[[i]])@data[index, ]
      exprs(gselist[[i]]) <- exprs(gselist[[i]])[ ,index]                     
    }
    
    #Elimate NILL patient
    temp <- pData(gselist[[i]])
    temp$platform <- NULL
    for (j in 1:nrow(temp)){        
      if (all(temp[j,]=="NILL")){        
        exprs(gselist[[i]]) <- exprs(gselist[[i]])[,-j]
        phenoData(gselist[[i]])@data <- phenoData(gselist[[i]])@data[-j,]
      }      
    }
    phenoData(gselist[[i]])@data<- phenoData(gselist[[i]])@data[,match(colnames(phenoData(gselist[[1]])@data),colnames(phenoData(gselist[[i]])@data))]    
    colnames(Biobase::featureData(gselist[[i]])@data) <- c("ENTREZID", "SYMBOL", "PROBE")
  }
  
  #logout
  InSilicoLogout() 
  return(gselist)
}
