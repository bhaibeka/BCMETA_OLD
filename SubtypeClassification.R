SubtypeClassification <- function(gselist, config.file){
  # Classify GSE subtype and give out the probability of being that subtype
  #
  # Args:
  #   gselist: It's a list containing all GSE to be classify
  #   config.fil: he string name of the standard configuration CSV file.
  #
  # Returns:
  #   A list containing all of the classification with the probability matrix for each dataset
  
  library(genefu)
  library(Biobase)
  # Error handling
  FTO <- read.csv(config.file)
  model <- as.character(FTO[3,2])  
  accepted.model <- c("ssp2003.robust","ssp2006.robust","pam50","scmgene.robust","scmod1.robust","scmod2.robust")
  if (!any(model==accepted.model)){
    stop("The model is invalid")
  }
  
  #SubType.List (STL)
  STL<- list()  
  ##############################
  ###  Creating STL matrix  ####
  ##############################
  
  for (i in 1:length(gselist)){     
    #Formating all the matrix to use the clustering function.
    data <- t(exprs(gselist[[i]]))   
    colnames(data)<- rownames(exprs(gselist[[i]]))
    mapping <- matrix(cbind(colnames(data),as.vector(Biobase::featureData(gselist[[i]])@data[,1])),length(rownames(exprs(gselist[[i]]))),2)  
    colnames(mapping) <- c("probe","EntrezGene.ID")
    rownames(mapping) <- mapping[ , "probe"]
    annot <- as.matrix(Biobase::featureData(gselist[[i]])@data)
    colnames(annot)[1] <- "EntrezGene.ID"
        
    #Depending on the model we use either prediction function
    accepted.intrinsic.model <- c("ssp2003.robust","ssp2006.robust","pam50")
    if (any(model==accepted.intrinsic.model)){
      subtype <- intrinsic.cluster.predict(sbt.model=eval(parse(text=model)), data=data, annot=mapping, do.mapping=TRUE, do.prediction.strength=FALSE)
      matrix_sub <- cbind(subtype$subtype, rep(1,nrow(data)), subtype$subtype.proba[ , 3] + subtype$subtype.proba[ , 4], subtype$subtype.proba)
      colnames(matrix_sub) <- c("subtype", "Global population", "Lums", "Basal", "Her2", "LumB", "LumA") 
      } else {
      subtype <- subtype.cluster.predict(sbt.model=eval(parse(text=model)), data=data, annot=mapping, do.mapping=TRUE, do.prediction.strength=FALSE, do.BIC=FALSE, plot=FALSE)
      matrix_sub <- cbind(subtype$subtype2, rep(1,nrow(data)), subtype$subtype.proba2[ , 3] + subtype$subtype.proba2[ , 4],subtype$subtype.proba2)
      colnames(matrix_sub) <- c("subtype", "Global population", "Lums", "Basal", "Her2", "LumB", "LumA")      
      }
    rownames(matrix_sub) <- colnames(exprs(gselist[[i]]))
    STL[[i]] <- matrix_sub
  }
    return(STL)
}