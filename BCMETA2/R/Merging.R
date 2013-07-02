Merging <-
function(gselist, STL, duplication.checker=TRUE, survdata=c("rfs", "dmfs", "dfs", "os"), time.cens=10,method){
  #Open every GSE, specify in the configuration CSV, from InsilicoDB and create a list of Eset
  #structure call gselist
  #
  # Args:
  #   gselist: The list containing all GSE file that need to be merge.
  #   STL: List containing all the subtype information for each GSE dataset.        
  #   duplication.checker: A marker, either TRUE or FALSE if you you want to verify
  #                        wheter or not you have duplicate sample into your master 
  #                        gene expression matrix.
  #   survdata: For t.fs and e.fs
  #   time.cens: maximum follow up (years)
  #   Method: either "unique" or "intersect" is use to for selecting geneid
  #
  # Returns:     
  #     The merging eset
  
  require(Biobase)
  require(survcomp)
  survdata <- match.arg(survdata)  
  
  #Creating gene array containing all ENTREZID depending on the chosen method
  gene.id <- NULL
  symbol <- NULL
  #gene.symbol.comparison.matrix
  GSCM <- NULL  
  for (i in 1:length(gselist)) {
    if (method == "intersect"){
      if (i ==1){
        gene.id <- rownames(exprs(gselist[[i]]))
        gene.id2 <- sapply(gene.id, function(x){strsplit(x, "\\.")[[1]][2]})
        symbol <- Biobase::featureData(gselist[[i]])@data$SYMBOL[match(gene.id,Biobase::featureData(gselist[[i]])@data$ENTREZID)]
      }
      else{
        gene.id <- intersect(rownames(exprs(gselist[[i]])), gene.id)
        gene.id2 <- sapply(gene.id, function(x){strsplit(x, "\\.")[[1]][2]})
        symbol <- intersect(Biobase::featureData(gselist[[i]])@data$SYMBOL[match(gene.id2,Biobase::featureData(gselist[[i]])@data$ENTREZID)], symbol)
      } 
    } else{
      gene.id <- c(gene.id, rownames(exprs(gselist[[i]])))
      gene.id2 <- sapply(rownames(exprs(gselist[[i]])), function(x){strsplit(x, "\\.")[[1]][2]})
      symbol <- c(symbol, Biobase::featureData(gselist[[i]])@data$SYMBOL[match(gene.id2,Biobase::featureData(gselist[[i]])@data$ENTREZID)])
    }        
  }
  
  GSCM <- cbind(gene.id,symbol)
  if (method == "unique"){
    gene.id <- unique(gene.id)
    index <- match(gene.id, GSCM[,1])
    GSCM <- GSCM[index,]
  }  
  
  matrix.featureData <- as.matrix(gene.id)
  colnames(matrix.featureData) <- "ENTREZID"
  
  #Rearranging the expression matrix
  matrix.exprs <- NULL
  for (i in 1:length(gselist)) {
    temp <- matrix(NA,length(gene.id),ncol(exprs(gselist[[i]])))
    colnames(temp) <- colnames(exprs(gselist[[i]]))   
    if (method == "unique"){
      matcher <- match(rownames(exprs(gselist[[i]])),gene.id)  
      matcher <- na.exclude(matcher)
      temp[matcher,]=exprs(gselist[[i]])
    }
    if (method == "intersect"){
      matcher <- match(gene.id,rownames(exprs(gselist[[i]])))  
      matcher <- na.exclude(matcher)
      temp=exprs(gselist[[i]])[matcher,]
    }
    matrix.exprs <- cbind(matrix.exprs, temp)
  }
  matrix.exprs.name <- NULL
  if (method == "unique"){
    for (i in 1:length(gselist)) {
      matrix.exprs.name <- c(matrix.exprs.name, rownames(exprs(gselist[[i]])))    
    }  
  } else {matrix.exprs.name <- gene.id}
  matrix.exprs.name <- unique(matrix.exprs.name)
  rownames(matrix.exprs) <- matrix.exprs.name
  
  ###### Scaling ########  
  message("#######        rescaling in process       #######")
  temp <- colnames(matrix.exprs)
  matrix.exprs <- sapply(1:ncol(matrix.exprs),function(x){((genefu::rescale(matrix.exprs[,x],na.rm=TRUE,q=0.05))-0.5)*2})
  matrix.exprs <- normalizeBetweenArrays(matrix.exprs) 
  colnames(matrix.exprs) <- temp
  message("#######       rescaling terminated        #######")
  
  #######################################
  ######   Find duplicate option   ######
  #######################################
  if (duplication.checker){
    gene.treshold <- 1000
    gene.var <- lapply(1:nrow(matrix.exprs), function(x){var(matrix.exprs[x,], na.rm=TRUE)})
    gene.var <- do.call(cbind, gene.var)
    #Most.Varian.Gene (MVG)
    index <- order(gene.var, decreasing=TRUE)[1:gene.treshold]
    #Matrix.of.Most.Variable.Gene (MMVG)
    MMVG <- matrix.exprs[index, , drop=FALSE]    
    
    #Setting every parallelizing parameter
    nbcore <- 16
    availcore <- detectCores()
    if(nbcore > availcore) { nbcore <- availcore }
    options("mc.cores"=nbcore)
    splitix <- splitIndices(nx=ncol(MMVG), ncl=nbcore)
    
    #Creating correlation matrix
    cor.matrix <- mclapply(splitix, function(splitix2, ...){      
      cor.matrix.temp <- sapply(1:ncol(MMVG), function(x){cor(MMVG[,splitix2, drop=FALSE], MMVG[,x, drop=FALSE], method="spearman", use="pairwise.complete.obs")})
      return(cor.matrix.temp)      
      },MMVG=MMVG)    
    
    #Retrieving parallel output
    cor.matrix <- do.call(rbind, cor.matrix)
    colnames(cor.matrix) <- rownames(cor.matrix) <- colnames(MMVG)
    #cor.matrix <- cor(MMVG, method="spearman", use="pairwise.complete.obs") 
    
    #Starting duplicated analysis      
    temp1 <- matrix.exprs
    temp2 <- cor.matrix
    #Double.Erase.Checker (DEC)
    DEC <- NULL
    GSM.erase <- NULL
    diag(cor.matrix) <- NA
    for (i in 1:ncol(cor.matrix)){
      if (!any(DEC == colnames(cor.matrix)[i])){
        #Duplicate.Checker (DC)        
        DC <- which(!is.na(cor.matrix[ , i]) & round(cor.matrix[ , i],digits=3) >= 0.96)
        if (length(DC) != 0){
          message(sprintf("%s is a duplicate of %s ", colnames(cor.matrix)[i], rownames(cor.matrix)[DC]))
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
  if (duplication.checker==FALSE){ GSM.erase <- NULL }
  matrix.phenoData <- NULL 
  matrix.subtype <- NULL
  for (i in 1:length(gselist)){     
    matrix.phenoData <- rbind(matrix.phenoData, as.matrix(pData(gselist[[i]])))    
    matrix.subtype <- rbind(matrix.subtype, STL[[i]])
  }   
  index1 <- c("age","size","node","er","pgr","her2","grade","t.rfs","e.rfs","t.dmfs","e.dmfs","t.os","e.os")
  matrix.phenoData[, index1] <- as.numeric(matrix.phenoData[, index1])
  
  #Collecting the survival data  
  surv.time <- surv.event <- rep(NA, nrow(matrix.phenoData))
  names(surv.time) <- names(surv.event) <- rownames(matrix.phenoData)
  switch (survdata,
    "rfs"={
      surv.time <- as.numeric(matrix.phenoData[ , "t.rfs"])
      surv.event <- as.numeric(matrix.phenoData[ , "e.rfs"])
    },
    "dmfs"={
      surv.time <- as.numeric(matrix.phenoData[ , "t.dmfs"])
      surv.event <- as.numeric(matrix.phenoData[ , "e.dmfs"])
    },
    "os"={
      surv.time <- as.numeric(matrix.phenoData[ , "t.os"])
      surv.event <- as.numeric(matrix.phenoData[ , "e.os"])
    },
    "dfs"={
      #hybrid between rfs and dmfs
      #use rfs when available, dmfs otherwise
      surv.time <- as.numeric(matrix.phenoData[ , "t.rfs"])
      surv.time[is.na(as.numeric(matrix.phenoData[ , "t.rfs"]))] <- as.numeric(matrix.phenoData[is.na(as.numeric(matrix.phenoData[ , "t.rfs"])), "t.dmfs"])
      surv.event <- as.numeric(matrix.phenoData[ , "e.rfs"])
      surv.event[is.na(as.numeric(matrix.phenoData[ , "e.rfs"]))] <- as.numeric(matrix.phenoData[is.na(as.numeric(matrix.phenoData[ , "e.rfs"])), "e.dmfs"])
    }
  )  
  ss <- survcomp::censor.time(surv.time=surv.time / 365, surv.event=surv.event, time.cens=time.cens)
  matrix.phenoData <- cbind(matrix.phenoData, "surv.time"=ss[[1]], "surv.event"=ss[[2]])
  
  
  if (!is.null(GSM.erase)){
    #Rearranging the master phenoData matrix
    index <- match(GSM.erase, rownames(matrix.phenoData))
    matrix.phenoData <- matrix.phenoData[-index, , drop=FALSE]
    #Rearranging the master subtype matrix
    index <- na.exclude(match(GSM.erase, rownames(matrix.subtype)))
    matrix.subtype <- matrix.subtype[-index, , drop=FALSE]
  }
  
  #Creating all the AnnotadeDataFrame to implement in the master eSet
  matrix.phenoData <- cbind(matrix.phenoData,matrix.subtype)
  colnames(matrix.exprs) <- rownames(matrix.phenoData)     
  phenoData <- new("AnnotatedDataFrame", data=data.frame(matrix.phenoData), varMetadata=data.frame(labelDescription=colnames(matrix.phenoData)))
  featureData <- new("AnnotatedDataFrame", data=data.frame(matrix.featureData), varMetadata=data.frame(labelDescription="EntrezID"))
  master.eset <- new("ExpressionSet", phenoData = phenoData, exprs = matrix.exprs)
  Biobase::featureData(master.eset) <- featureData
  
  return(list("master.eset"=master.eset, "GSCM"=GSCM))
}
