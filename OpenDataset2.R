OpenDataset2 <- function(config.file){
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
  
  #Opening all dataset from InsilicoDB
  for (i in 16:nrow(FTO)){
    if (as.character(FTO[i,4])==TRUE){    
      GPL <- getPlatforms(dataset=as.character(FTO[i,2]))
      for (j in 1:length(GPL)){
        gselist[[counter]] <- getDatasets(dataset=as.character(FTO[i,2]), norm=as.character(FTO[i,6]), curation=FTO[i,3], features="PROBE")[[j]]       
        #making probe gene mapping for Affy structure
        if (as.character(FTO[i,5])=="Affy"){            
          gselist[[counter]] <- ProbeGeneMapping(gselist[[counter]])
        } else{        
          rownames(exprs(gselist[[counter]])) <- paste("geneid.",Biobase::featureData(gselist[[counter]])$ENTREZID, sep="")
          }
        if (length(GPL) > 1){
          index <- length(gselist)-GPL.length+1
          gselist[[index]] <- PlatformMerging(gselist=gselist, GPL.length=length(GPL))
          for (k in 1:length(gselist)){
            gselist[[index+k]] <- NULL
          }          
        }
        counter <- counter +1
      } 
    }
  }    
  
  #Making sure that every column of the clinical data matrix is well aligned for every GSE dataset
  for (i in 1:length(gselist)){    
    phenoData(gselist[[i]])@data<- phenoData(gselist[[i]])@data[,match(colnames(phenoData(gselist[[1]])@data),colnames(phenoData(gselist[[i]])@data))]    
  }

  #logout
  InSilicoLogout() 
  return(gselist)
}