SubtypeClassification <- function(gselist, config.file){
  # Classify GSE subtype and give out the probability of being that subtype
  #
  # Args:
  #   gselist: It's a list containing all GSE to be classify
  #   config.fil: he string name of the standard configuration CSV file.
  #
  # Returns:
  #   A list containing all of the classification with the probability matrix for each dataset
  
  require(genefu)
  # Error handling
  FTO <- read.csv(config.file)
  model <- as.character(FTO[3,2])  
  accepted.model <- c("ssp2003.robust","ssp2006.robust","pam50","scmgene.robust","scmod1.robust","scmod2.robust")
  if (!any(model == accepted.model)){
    stop("The model is invalid")
  }
  
  #SubType.List (STL)
  STL<- list()  
  for (i in 1:length(gselist)){  
    #Formating all the matrix to use the clustering function.
    data <- t(exprs(gselist[[i]]))
    colnames(data) <- as.vector(paste("geneid",featureData(gselist[[i]])@data[,1],sep="_"))
    mapping <- matrix(cbind(colnames(data),as.vector(featureData(gselist[[i]])@data[,1])),length(rownames(exprs(gselist[[i]]))),2)  
    colnames(mapping) <- c("probe","EntrezGene.ID")
    rownames(mapping) <- mapping[ , "probe"]
    annot <- as.matrix(featureData(gselist[[i]])@data)
    colnames(annot)[1] <- "EntrezGene.ID"
    
    #Depending on the model we use either prediction function
    accepted.intrinsic.model <- c("ssp2003.robust","ssp2006.robust","pam50")
    if (any(model == accepted.intrinsic.model)){
      subtype <- intrinsic.cluster.predict(sbt.model=eval(parse(text=model)), data=data, annot=mapping, do.mapping=TRUE, do.prediction.strength=FALSE)
      matrix_sub <- matrix(cbind(subtype$subtype,subtype$subtype.proba),length(subtype$subtype),ncol(cbind(subtype$subtype,subtype$subtype.proba)))
      } else {
      subtype <- subtype.cluster.predict(sbt.model=eval(parse(text=model)), data=data, annot=mapping, do.mapping=TRUE, do.prediction.strength=FALSE, do.BIC=FALSE, plot=FALSE)
      matrix_sub <- matrix(cbind(subtype$subtype2,subtype$subtype.proba2),length(subtype$subtype),ncol(cbind(subtype$subtype2,subtype$subtype.proba2)))
      colnames(matrix_sub) <- c("subtype", "p1.ER-/HER2-", "p1.HER2+", "p1.ER+/HER2-", "subtype2", "p2.ER-/HER2-", "p2.HER2+","p2.ER+/HER2- High Prolif", "p2.ER+/HER2- Low Prolif")      
      }
    rownames(matrix_sub) <- colnames(exprs(gselist[[i]]))
    STL[[i]] <- matrix_sub
  }
    return(STL)
}