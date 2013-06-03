OpenDataset <- function(config.file,duplication.checker){
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
  # Returns:
  #   -if duplication.checker == TRUE, then the function returns:
  #     The GSE list of the specify GSE number in the configuraiton CSV file, the master
  #     gene expression matrix, the correlation matrix and the erased sample vector
  #   -if duplication.checker == FALSE, then the function returns:
  #     The GSE list of the specify GSE number in the configuraiton CSV file and the master
  #     gene expression matrix
  
    
  require(inSilicoDb2)  
  InSilicoLogin(login="bhaibeka@gmail.com", password="747779bec8a754b91076d6cc1f700831")
  #File.To.Open (FTO)
  FTO <- read.csv(config.file) 
  gselist <- list()
  expr.rowname <- NULL
  counter <- 1
  #Opening all dataset from InsilicoDB
  for (i in 16:nrow(FTO)){
    if (as.character(FTO[i,4])==TRUE){    
      GPL <- getPlatforms(dataset=as.character(FTO[i,2]))[1]
      gselist[[counter]] <- getDataset(dataset=as.character(FTO[i,2]), curation=FTO[i,3], platform=GPL)
      expr.rowname <- c(expr.rowname, rownames(exprs(gselist[[counter]])))
      counter <- counter+1
    }
  } 
  InSilicoLogout()
  #Eliminating all doubles to create a master gene vector
  expr.rowname <- unique(expr.rowname)
  
  #Merging all the genes expressions of all GSE for all samples into
  #one master matrix whom each row are matching with the master gene vector
  for (i in 1:length(gselist)){
    temp <- matrix(NA,length(expr.rowname),ncol(exprs(gselist[[i]])))
    rownames(temp) <- expr.rowname
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
  
  
  #Clinical Info Standard (CIS)
  CIS <- c("tissue", "age", "node", "treatment", "pgr", "grade", "platform", 
        "size", "er", "t.rfs", "e.rfs", "id", "series", "her2", "t.dmfs", 
        "e.dmfs")
  # Error handling
  for (i in 1:length(gselist)){    
    if (any(varLabels(gselist[[i]]) != CIS) == TRUE){
      warning(sprintf("the %s dataset as been incorrectly curated", FTO[i,2]))
    }
  }
  
  if (duplication.checker == TRUE){
    return(list("gselist"=gselist, "matrix.exprs"=matrix.exprs, "cor.matrix"=cor.matrix, "GSM.erase"=GSM.erase))
  } else{
    return(list("gselist"=gselist, "matrix.exprs"=matrix.exprs))
  }
}