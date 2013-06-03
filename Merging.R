Merging <- function(gselist, STL, duplication.checker, checker){
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
  #   checker: For t.fs and e.fs
  # Returns:     
  #     The merging eset
  
  
  expr.rowname <- NULL
  gene.id <- NULL
  for (i in 1:length(gselist)){
    expr.rowname <- c(expr.rowname, rownames(exprs(gselist[[i]])))
    gene.id <- c(gene.id, Biobase::featureData(gselist[[i]])$ENTREZID)
  }
  expr.rowname <- unique(expr.rowname)
  gene.id <- unique(gene.id)
  matrix.featureData <- cbind(gene.id, expr.rowname)
  colnames(matrix.featureData) <- c("ENTREZID", "SYMBOL")
  
  #Rearranging the expression matrix
  for (i in 1:length(gselist)){
    temp <- matrix(NA,length(expr.rowname),ncol(exprs(gselist[[i]])))
    rownames(temp) <- as.vector(paste("geneid",Biobase::featureData(gselist[[i]])@data[,1],sep="_"))
    colnames(temp) <- colnames(exprs(gselist[[i]]))    
    matcher <- match(rownames(exprs(gselist[[i]])),expr.rowname)  
    temp[matcher,1:ncol(exprs(gselist[[i]]))]=exprs(gselist[[i]])
    if (i==1){
      matrix.exprs <- temp
    } else{
      matrix.exprs <- cbind(matrix.exprs,temp)
    }
  }

  # Find duplicate option
  if (duplication.checker==TRUE){
    gene.treshold <- 3000
    gene.var <- apply(matrix.exprs,1,var)
    #Most.Varian.Gene (MVG)
    MVG <- names(sort(gene.var, decreasing=TRUE)[1:gene.treshold])
    index <- match(MVG,rownames(matrix.exprs))
    #Matrix.of.Most.Variable.Gene (MMVG)
    MMVG <- matrix.exprs[index,]
    cor.matrix <- cor(MMVG) 
    ending <- ncol(cor.matrix)
    temp1 <- matrix.exprs
    temp2 <- cor.matrix
    #Double.Erase.Checker (DEC)
    DEC <- NULL
    GSM.erase <- NULL
    for (i in 1:ending){
      if (!any(DEC == colnames(cor.matrix)[i])){
        #Duplicate.Checker (DC)
        DC <- which(cor.matrix[,i]>0.95 & cor.matrix[,i]<0.999)
        if (length(DC)!=0){
          print(sprintf("%s is a duplicate of %s",colnames(cor.matrix)[i],rownames(cor.matrix)[DC]))
          DEC <- c(DEC, rownames(cor.matrix)[DC])        
          GSM.erase <- c(GSM.erase, (colnames(cor.matrix)[i]))
          #Index.to.Delete (ID)
          ID <- which(colnames(temp1)==colnames(cor.matrix)[i])
          temp1 <- temp1[,-ID]        
          temp2 <- temp2[-i, -i]                  
        }
      }
    }
    matrix.exprs <- temp1
    cor.matrix <- temp2
  }
  
  #Creating master phenoData matrix and master subtype matrix
  if (duplication.checker==FALSE){GSM.erase <- NULL}
  matrix.phenoData <- NULL 
  matrix.subtype <- NULL
  for (i in 1:length(gselist)){  
    matrix.phenoData <- rbind(matrix.phenoData, as.matrix(pData(gselist[[i]])))    
    matrix.subtype <- rbind(matrix.subtype, STL[[i]])
  }   
  
  #Creating the t.fs hybrid  
  if (checker==TRUE){    
    t.fs <- NULL
    e.fs <- NULL
    index.rfs <- which(colnames(matrix.phenoData)=="t.rfs")
    index.dmfs <- which(colnames(matrix.phenoData)=="t.dmfs") 
    index.erfs <- which(colnames(matrix.phenoData)=="e.rfs")
    index.edmfs <- which(colnames(matrix.phenoData)=="e.dmfs")
    for (i in 1:nrow(matrix.phenoData)){
      if (matrix.phenoData[i,index.dmfs]=="NA"){
        t.fs <- c(t.fs, matrix.phenoData[i,index.rfs])
        e.fs <- c(e.fs, matrix.phenoData[i,index.erfs])
      }else{
        t.fs <- c(t.fs, matrix.phenoData[i,index.dmfs])
        e.fs <- c(e.fs, matrix.phenoData[i,index.edmfs])
      }
    }
    matrix.phenoData <- cbind(matrix.phenoData, t.fs, e.fs)
  }
  # Rearranging the master phenoData matrix
  index <- match(GSM.erase, rownames(matrix.phenoData))
  for (i in index){
    matrix.phenoData <- matrix.phenoData[-i,]
  }  
  # Rearranging the master subtype matrix
  index <- match(GSM.erase, rownames(matrix.subtype))
  for (i in index){
    matrix.subtype <- matrix.subtype[-i,]
  }  
  
  #Creating all the AnnotadeDataFrame to implement in the master eSet
  matrix.phenoData <- cbind(matrix.phenoData,matrix.subtype)
  colnames(matrix.exprs) <- rownames(matrix.phenoData)     
  phenoData <- new("AnnotatedDataFrame", data=data.frame(matrix.phenoData), varMetadata=data.frame(labelDescription=colnames(matrix.phenoData)))
  featureData <- new("AnnotatedDataFrame", data=data.frame(matrix.featureData), varMetadata=data.frame(labelDescription=c("EntrezID","Symbol")))
  master.eset <- new("ExpressionSet", phenoData = phenoData, exprs = matrix.exprs)
  Biobase::featureData(master.eset) <- featureData
  
  return(master.eset)
}
         
  
  
  