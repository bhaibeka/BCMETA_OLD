SubtypeClassification2 <- function(matrix.exprs, config.file){
  # Classify GSE subtype and give out the probability of being that subtype
  #
  # Args:
  #   matrix.exprs: Is a master matrix containing all the sample from the Opendataset function
  #   config.fil: he string name of the standard configuration CSV file.
  #
  # Returns:
  #   A list containing all of the classification with the probability matrix for each dataset
  
  library(genefu)
  # Error handling
  FTO <- read.csv(config.file)
  model <- as.character(FTO[3,2])  
  accepted.model <- c("ssp2003.robust","ssp2006.robust","pam50","scmgene.robust","scmod1.robust","scmod2.robust")
  if (!any(model==accepted.model)){
    stop("The model is invalid")
  }  
    
  #Formating all the matrix to use the clustering function.
  data <- t(matrix.exprs)  
  mapping <- matrix(cbind(colnames(data), sub("geneid_","",colnames(data))),ncol(data),2)  
  colnames(mapping) <- c("probe","EntrezGene.ID")
  rownames(mapping) <- mapping[ , "probe"]  
    
  #Depending on the model we use either prediction function
  accepted.intrinsic.model <- c("ssp2003.robust","ssp2006.robust","pam50")
  if (any(model==accepted.intrinsic.model)){
    subtype <- intrinsic.cluster.predict(sbt.model=eval(parse(text=model)), data=data, annot=mapping, do.mapping=TRUE, do.prediction.strength=FALSE)
    matrix_sub <- matrix(cbind(subtype$subtype,subtype$subtype.proba),length(subtype$subtype),ncol(cbind(subtype$subtype,subtype$subtype.proba)))
  } else {
    subtype <- subtype.cluster.predict(sbt.model=eval(parse(text=model)), data=data, annot=mapping, do.mapping=TRUE, do.prediction.strength=FALSE, do.BIC=FALSE, plot=FALSE)
    matrix_sub <- matrix(cbind(subtype$subtype2,subtype$subtype.proba2),length(subtype$subtype),ncol(cbind(subtype$subtype2,subtype$subtype.proba2)))
    colnames(matrix_sub) <- c("subtype2", "ER-|HER2-", "HER2+","ER+|HER2- High Prolif", "ER+|HER2- Low Prolif")      
  }
  rownames(matrix_sub) <- colnames(matrix.exprs)
  
  return(matrix_sub)
}