ProbeGeneMapping2 <- function(gse_input){
  #Open every GSE, specify in the configuration CSV, from InsilicoDB and create a list of Eset
  #structure call gselist
  #
  # Args:
  #     gse_input: the GSE eset that need to be map.  
  #   
  # Returns:     
  #     The GSE with the best probe gene mapping 
  
  
  gse <- gse_input
  gid <- as.character(Biobase::featureData(gse)@data$ENTREZID)
  names(gid) <- rownames(exprs(gse))
  ugid <- unique(gid)
  names(ugid) <- paste("geneid", ugid, sep=".")
  exprs.matrix <- t(exprs(gse))
  rr <- genefu::geneid.map(geneid1=gid, data1=exprs.matrix, geneid2=ugid)
  exprs.matrix <- exprs.matrix[ , match(names(rr$geneid1),colnames(exprs.matrix)), drop=FALSE]
  index <- match(names(rr$geneid1),names(gid))
  gene.id <- gid[index]
  gene.id2 <- sapply(gene.id, function(x){strsplit(x, " ///")[[1]][1]})
  gene.id3 <- paste("geneid",gene.id2,sep=".")
  colnames(exprs.matrix) <- gene.id3
  exprs(gse) <- t(exprs.matrix)
  gse_output <- gse  
  
  return(gse_output)
}