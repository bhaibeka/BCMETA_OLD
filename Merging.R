Merging <- function(gselist, STL, duplication.checker=TRUE, survdata=c("rfs", "dmfs", "dfs", "os"), time.cens=10){
  #Open every GSE, specify in the configuration CSV, from InsilicoDB and create a list of Eset
  #structure call gselist
  #
  # Args:
  #   gselist: The list containing all GSE file that need to be merge.
  #   STL: List containing all the subtype information for each GSE dataset.
  #                  you want to consider for your correlation statistical analysis.
  #   duplication.checker: A marker, either TRUE or FALSE if you you want to verify
  #                        wheter or not you have duplicate sample into your master 
  #                        gene expression matrix.
  #   survdata: For t.fs and e.fs
  #   time.cens: maximum follow up (years)
  # Returns:     
  #     The merging eset
  
  require(Biobase)
  require(survcomp)
  survdata <- match.arg(survdata)
  expr.rowname <- NULL
  gene.id <- NULL
  for (i in 1:length(gselist)) {
    expr.rowname <- c(expr.rowname, rownames(exprs(gselist[[i]])))
    gene.id <- c(gene.id, Biobase::featureData(gselist[[i]])$ENTREZID)
  }
  expr.rowname <- unique(expr.rowname)
  gene.id <- unique(gene.id)
  matrix.featureData <- cbind(gene.id, expr.rowname)
  colnames(matrix.featureData) <- c("ENTREZID", "SYMBOL")
  
  #Rearranging the expression matrix
  matrix.exprs <- NULL
  for (i in 1:length(gselist)) {
    temp <- matrix(NA,length(expr.rowname),ncol(exprs(gselist[[i]])))
    rownames(temp) <- as.vector(paste("geneid",Biobase::featureData(gselist[[i]])@data[,1],sep="."))
    colnames(temp) <- colnames(exprs(gselist[[i]]))    
    matcher <- match(rownames(exprs(gselist[[i]])),expr.rowname)  
    temp[matcher,1:ncol(exprs(gselist[[i]]))]=exprs(gselist[[i]])
    matrix.exprs <- cbind(matrix.exprs, temp)
  }

  # Find duplicate option
  if (duplication.checker){
    gene.treshold <- 1000
    gene.var <- apply(matrix.exprs, 1, var)
    #Most.Varian.Gene (MVG)
    index <- order(gene.var, decreasing=TRUE)[1:gene.treshold]
    #Matrix.of.Most.Variable.Gene (MMVG)
    MMVG <- matrix.exprs[index, , drop=FALSE]
    cor.matrix <- cor(MMVG, method="spearman", use="pairwise.complete.obs")
    ending <- ncol(cor.matrix)
    temp1 <- matrix.exprs
    temp2 <- cor.matrix
    #Double.Erase.Checker (DEC)
    DEC <- NULL
    GSM.erase <- NULL
    for (i in 1:ending){
      if (!any(DEC == colnames(cor.matrix)[i])){
        #Duplicate.Checker (DC)
        DC <- which(cor.matrix[ , i] >= 0.95)
        if (length(DC) != 0){
          message(sprintf("%s is a duplicate of %s", colnames(cor.matrix)[i], rownames(cor.matrix)[DC]))
          DEC <- c(DEC, rownames(cor.matrix)[DC])        
          GSM.erase <- c(GSM.erase, (colnames(cor.matrix)[i]))
          #Index.to.Delete (ID)
          ID <- which(colnames(temp1) == colnames(cor.matrix)[i])
          temp1 <- temp1[ , -ID, drop=FALSE]
          temp2 <- temp2[-i, -i, drop=FALSE]        
        }
      }
    }
    matrix.exprs <- temp1
    cor.matrix <- temp2
  }
  
  #Creating master phenoData matrix and master subtype matrix
  if (duplication.checker){ GSM.erase <- NULL }
  matrix.phenoData <- NULL 
  matrix.subtype <- NULL
  for (i in 1:length(gselist)){  
    matrix.phenoData <- rbind(matrix.phenoData, as.matrix(pData(gselist[[i]])))    
    matrix.subtype <- rbind(matrix.subtype, STL[[i]])
  }   
  
  #Collecting the survival data
  surv.time <- surv.event <- rep(NA, nrow(matrix.phenoData))
  names(surv.time) <- names(surv.event) <- rownames(matrix.phenoData)
  switch (survdata,
    "rfs"={
      surv.time <- matrix.phenoData[ , "t.rfs"]
      surv.event <- matrix.phenoData[ , "e.rfs"]
    },
    "dmfs"={
      surv.time <- matrix.phenoData[ , "t.dmfs"]
      surv.event <- matrix.phenoData[ , "e.dmfs"]
    },
    "os"={
      surv.time <- matrix.phenoData[ , "t.os"]
      surv.event <- matrix.phenoData[ , "e.os"]
    },
    "dfs"={
      #hybrid between rfs and dmfs
      #use rfs when available, dmfs otherwise
      surv.time <- matrix.phenoData[ , "t.rfs"]
      surv.time[is.na(matrix.phenoData[ , "t.rfs"])] <- matrix.phenoData[is.na(matrix.phenoData[ , "t.rfs"]), "t.dmfs"]
      surv.event <- matrix.phenoData[ , "e.rfs"]
      surv.event[is.na(matrix.phenoData[ , "e.rfs"])] <- matrix.phenoData[is.na(matrix.phenoData[ , "e.rfs"]), "e.dmfs"]
    }
  )
  ss <- survcomp::censor.time(surv.time=surv.time / 365, surv.event=surv.event, time.cens=10)
  matrix.phenoData <- cbind(matrix.phenoData, "surv.time"=ss[[1]], "surv.event"=ss[[2]])
  
  #Rearranging the master phenoData matrix
  index <- match(GSM.erase, rownames(matrix.phenoData))
  matrix.phenoData <- matrix.phenoData[-index, , drop=FALSE]  
  #Rearranging the master subtype matrix
  index <- match(GSM.erase, rownames(matrix.subtype))
  matrix.subtype <- matrix.subtype[-index, , drop=FALSE]
  
  #Creating all the AnnotadeDataFrame to implement in the master eSet
  matrix.phenoData <- cbind(matrix.phenoData,matrix.subtype)
  colnames(matrix.exprs) <- rownames(matrix.phenoData)     
  phenoData <- new("AnnotatedDataFrame", data=data.frame(matrix.phenoData), varMetadata=data.frame(labelDescription=colnames(matrix.phenoData)))
  featureData <- new("AnnotatedDataFrame", data=data.frame(matrix.featureData), varMetadata=data.frame(labelDescription=c("EntrezID","Symbol")))
  master.eset <- new("ExpressionSet", phenoData = phenoData, exprs = matrix.exprs)
  Biobase::featureData(master.eset) <- featureData
  
  return(master.eset)
}
         
  
  
  