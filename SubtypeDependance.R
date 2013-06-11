SubtypeDependance <- function(master.eset,gene.list){
  # Classify GSE subtype and give out the probability of being that subtype
  #
  # Args:
  #   master.eset: A master.eset containing everything
  #   gene.list : a list of gene which are to be analyze
  #
  # Returns:
  #   The dependance of a specific gene over the subtype
  
  ###############################################################################################
  ###########################        Kluskal.Wallis test idea         ########################### 
  ###############################################################################################
  gene.id <- rownames(exprs(master.eset))
  library(parallel)
  
  #Setting every parallelizing parameter
  nbcore <- 8
  availcore <- detectCores()
  if(nbcore > availcore) { nbcore <- availcore }
  options("mc.cores"=nbcore)
  splitix <- splitIndices(nx=length(gene.id), ncl=nbcore)
  
  #######  Kruskal core function   #########
  
  Kruskal <- mclapply(splitix, function(splitix2,...){
    sub.kruskal <- lapply(splitix2, function(x){
      gene.exprs <- exprs(master.eset)[gene.id[x], ]
      subtype.factor <- factor(pData(master.eset)$subtype)
      subsub.kruskal <- kruskal.test(gene.exprs, subtype.factor) 
      gene.med <- median(gene.exprs, na.rm=TRUE)
      res <- cbind(gene.id[x], gene.med, as.numeric(subsub.kruskal[[1]]), as.numeric(subsub.kruskal[[3]]))
      return(res)       
    })
    sub.kruskal <- do.call(rbind, sub.kruskal) 
    return(sub.kruskal)
  }, gene.id=gene.id, master.eset=master.eset)
  Kruskal <- do.call(rbind, Kruskal)  
  colnames(Kruskal) <- c("Gene ID", "median", "KW.score", "p.value")
  rownames(Kruskal) <- gene.id
  
  
  #############################################################################################
  #########################           Wilcoxson test idea            ##########################
  #############################################################################################
  
  #######   Indexation of everything   #########
  
  # Indexing everything for proper subtype classification
  index2 <- which(colnames(pData(master.eset))=="subtype")
  subtype.name <- colnames(pData(master.eset))[(index2+1):ncol(pData(master.eset))]
  subtype.comb <- combinations(n=length(subtype.name), r=2)
  
  #Changing clinical data subtype to a more suitable subtype name
  list.index<- list() 
  temp <- as.matrix(pData(master.eset))
  index3 <- which(colnames(temp)=="subtype")
  temp[which(temp[, index3]=="ER+/HER2- High Prolif"),index3] <- "LumB"
  temp[which(temp[, index3]=="ER+/HER2- Low Prolif"),index3] <- "LumA"
  temp[which(temp[, index3]=="ER-/HER2-"),index3] <- "Basal"
  temp[which(temp[, index3]=="HER2+"),index3] <- "Her2"
 
  #Generating List of index for classification
  for (i in 1:length(subtype.name)){    
    if (subtype.name[[i]]=="Global.population"){list.index[["Global"]] <- 1:nrow(pData(master.eset))}
    if (subtype.name[[i]]=="Lums"){list.index[[subtype.name[[i]]]] <- NA}  
    list.index[[subtype.name[[i]]]] <- which(temp[,index3]==subtype.name[i])
  }
  
  #Rearranging List of index
  if (!is.null(list.index$Lums)){list.index$Lums <- c(list.index$LumA , list.index$LumB)} 
  list.index <- list(list.index$Basal,list.index$Her2,list.index$LumB,list.index$LumA,list.index$Lums, list.index$Global)
  
  #######   Preparation for test core function   #########  
  
  comb <- combinations(length(list.index),2)  
  list.coxson <- list()
  name.cox <- c("Basal","Her2","Lums","LumB","LumA","Global")
  gene.id <- rownames(exprs(master.eset))[1:80]
  #gene.id <- rownames(exprs(master.eset))
  ###########################################
  #######    Main core function     #########
  ###########################################  
  library(parallel)
  
  #Setting every parallelizing parameter
  nbcore <- 8
  availcore <- detectCores()
  if(nbcore > availcore) { nbcore <- availcore }
  options("mc.cores"=nbcore)
  splitix <- splitIndices(nx=length(gene.id), ncl=nbcore)
  
  #######  Median main core function   #########
  
  median.matrix <- mclapply(splitix, function(splitix2,...){
    sub.median.matrix <- lapply(splitix2, function(x){
      subsub.median.matrix <- lapply(1:length(list.index), function(z){median(exprs(master.eset)[gene.id[x],list.index[[z]]],na.rm=TRUE)})
      subsub.median.matrix <- do.call(cbind, subsub.median.matrix) 
      subsub.median.matrix <- cbind(gene.id[x], subsub.median.matrix)
      return(subsub.median.matrix)
      })
      sub.median.matrix <- do.call(rbind, sub.median.matrix)      
  })
  median.matrix <- do.call(rbind, median.matrix) 
  colnames(median.matrix) <- c("Gene ID", name.cox)
  rownames(median.matrix) <- gene.id
  
  #######  Wilcoxson matrix core function   #########
  
  Wilcox <- mclapply(splitix, function(splitix2,...){
    sub.wilcox <- lapply(splitix2, function(x){      
      subsub.wilcox <- lapply(1:nrow(comb), function(z){
        gene.exprs1 <- exprs(master.eset)[gene.id[x], list.index[[comb[z,1]]]]
        gene.exprs2 <- exprs(master.eset)[gene.id[x], list.index[[comb[z,2]]]]
        subsubsub.wilcox <- wilcox.test(gene.exprs1, gene.exprs2, conf.int = TRUE)
        res <- cbind(paste(name.cox[comb[z,1]], name.cox[comb[z,2]], sep="."), as.numeric(subsubsub.wilcox[[1]]), as.numeric(subsubsub.wilcox[[3]]), as.numeric(subsubsub.wilcox[[8]][1]), as.numeric(subsubsub.wilcox[[8]][2]), as.numeric(subsubsub.wilcox[[9]]))
        return(res)
      })  
      subsub.wilcox <- do.call(cbind, subsub.wilcox) 
      subsub.wilcox <- cbind(gene.id[x], subsub.wilcox)
      return(subsub.wilcox)
    })
    sub.wilcox <- do.call(rbind, sub.wilcox)
    return(sub.wilcox)
  })  
  Wilcox <- do.call(rbind, Wilcox)   
  matrix.cox.name <- "Gene ID"
  for (i in 1:nrow(comb)){
    matrix.cox.name <- c(matrix.cox.name, "Comparison", "W.score", "p.value", "Lower Marge", "Higher Marge", "distance")
  }
  colnames(Wilcox) <- matrix.cox.name
  rownames(Wilcox) <- gene.id
  
 
  
  
  