PlatformMerging <-
function(gselist, GPL.length){
  #Open every GSE, specify in the configuration CSV, from InsilicoDB and create a list of Eset
  #structure call gselist
  #
  # Args:
  #     gselist: Containt the gselist with every GSE. 
  #     GPL.length: Is the length of the GPL vector, it specify how many platform that need to be merge
  #   
  # Returns:     
  #     The GSE which each platform has been merge
  
  
  index <- length(gselist)-GPL.length+1
  for (i in 1:GPL.length-1){
    missing.genes <- setdiff(rownames(exprs(gselist[[index+i]])),rownames(exprs(gselist[[index]])))
    exprs(gselist[[index]]) <- rbind(exprs(gselist[[index]]), exprs(gselist[[index+i]])[match(missing.genes, rownames(exprs(gselist[[index+i]]))), ,drop=FALSE])
    gene.id <- sapply(missing.genes, function(x){strsplit(x, "\\.")[[1]][2]})
    index2 <- match(gene.id, Biobase::featureData(gselist[[index+i]])@data[ ,1])
    gene.annotation.matrix <- Biobase::featureData(gselist[[index+i]])@data[index2, ]
    colnames(gene.annotation.matrix) <- colnames(Biobase::featureData(gselist[[index]])@data)
    Biobase::featureData(gselist[[index]])@data <- rbind(Biobase::featureData(gselist[[index]])@data, gene.annotation.matrix)    
  }
  return(gselist[[index]])
}
