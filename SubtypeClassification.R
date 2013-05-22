SubtypeClassification <- function(gselist, model){
  # Classify GSE subtype and give out the probability of being that subtype
  #
  # Args:
  #   gselist: It's a list containing all GSE to be classify
  #   mode: Is a string expliciting the statistic model of prediction choosen between:
  #   ssp2003.robust, ssp2006.robust, pam50, scmgene.robust, scmod1.robust, scmod2.robust       
  #
  # Returns:
  #   The classification with the probability matrix
  
  library(genefu)
  # Error handling
  if (model==ssp2003.robust | model==ssp2006.robust | model==pam50 | model==scmgene.robust | model==scmod1.robust | model==scmod2.robust){
    stop("The model is invalid")
  }
  #SubType List (STL)
  STL<- list()
  #Formating all the matrix to use the clustering function.
  for (i in 1:length(gselist)){     
    data <- t(exprs(gselist[[i]]))
    colnames(data) <- as.vector(paste(featureData(gselist[[i]])@data[,1],c(1:length(featureData(gselist[[i]])@data[,1])),sep="."))
    mapping <- matrix(cbind(rownames(exprs(gselist[[i]])),as.vector(featureData(gselist[[i]])@data[,1])),length(rownames(exprs(gselist[[i]]))),2)  
    colnames(mapping) <- c("probe","EntrezGene.ID")
    annot <- as.matrix(featureData(gselist[[i]])@data)
    colnames(annot)[1] <- "EntrezGene.ID"
    
    #Depend on the model
    if (model==ssp2003.robust | model==ssp2006.robust | model==pam50 | model==scmgene.robust){
      subtype <- intrinsic.cluster.predict(sbt.model=model, data=data, annot=mapping, do.mapping=TRUE, do.prediction.strength=FALSE, do.BIC=FALSE, plot=FALSE)
    } else {
      subtype <- subtype.cluster.predict(sbt.model=model, data=data, annot=mapping, do.mapping=TRUE, do.prediction.strength=FALSE, do.BIC=FALSE, plot=FALSE)
    }
    matrix_sub <- matrix(cbind(subtype$subtype,subtype$subtype.proba,subtype$subtype2,subtype$subtype.proba2),length(subtype$subtype),9)
    #colnames(matrix_sub) <- c("subtype", "p1.ER-/HER2-", "p1.HER2+", "p1.ER+/HER2-", "subtype2", "p2.ER-/HER2-", "p2.HER2+","p2.ER+/HER2- High Prolif", "p2.ER+/HER2- Low Prolif")
    rownames(matrix_sub) <- colnames(exprs(gselist[[i]]))
    STL[[i]] <- matrix_sub
  }
    return(STL)
}