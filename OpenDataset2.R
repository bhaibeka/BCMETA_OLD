OpenDataset2 <- function(config.file){
  #Open every GSE, specify in the configuration CSV, from InsilicoDB and create a list of Eset
  #structure call gselist
  #
  # Args:
  #   config.file: The string name of the configuration CSV file.
  #   duplication.checker: A marker, either TRUE or FALSE if you you want to verify
  #                        wheter or not you have duplicate sample into your master 
  #                        gene expression matrix.
  #   gene.treshold: If duplication.checker is TRUE, then the treshold is the number of most variable
  #                  you want to consider for your correlation statistical analysis.
  # Returns:  #   
  #     The GSE list of the specify GSE number in the configuraiton CSV file 
  
  
  library(inSilicoDb2)  
  InSilicoLogin(login="bhaibeka@gmail.com", password="747779bec8a754b91076d6cc1f700831")
  #File.To.Open (FTO)
  FTO <- read.csv(config.file) 
  gselist <- list() 
  counter <- 1
  #Opening all dataset from InsilicoDB
  for (i in 16:nrow(FTO)){
    if (as.character(FTO[i,4])==TRUE){    
      GPL <- getPlatforms(dataset=as.character(FTO[i,2]))[1]
      gselist[[counter]] <- getDataset(dataset=as.character(FTO[i,2]), curation=FTO[i,3], platform=GPL)      
      counter <- counter+1
    }
  } 
  InSilicoLogout() 
  return(gselist)
}