ProbeGeneMapping <- function(gse_input){
  #Open every GSE, specify in the configuration CSV, from InsilicoDB and create a list of Eset
  #structure call gselist
  #
  # Args:
  #     gse_input: the GSE eset that need to be map.  
  #   
  # Returns:     
  #     The GSE with the best probe gene mapping 
  
  gse <- gse_input
  js <- jetset.bhk::jscores(chip=annotation(gse), probeset=rownames(exprs(gse)))
  js <- js[rownames(exprs(gse)), , drop=FALSE]
  rownames(js) <- rownames(exprs(gse))
  ## identify the best probeset for each entrez gene id
  geneid1 <- as.character(js[ ,"EntrezID"])
  names(geneid1) <- rownames(js)
  geneid2 <- sort(unique(geneid1))
  names(geneid2) <- paste("geneid", geneid2, sep=".")
  gix1 <- !is.na(geneid1)
  gix2 <- !is.na(geneid2)
  geneid.common <- intersect(geneid1[gix1], geneid2[gix2])
  ## probes corresponding to common gene ids
  gg <- names(geneid1)[is.element(geneid1, geneid.common)]
  gid <- geneid1[is.element(geneid1, geneid.common)]
  ## duplicated gene ids
  gid.dupl <- unique(gid[duplicated(gid)])
  gg.dupl <- names(geneid1)[is.element(geneid1, gid.dupl)]
  ## unique gene ids
  gid.uniq <- gid[!is.element(gid, gid.dupl)]
  gg.uniq <- names(geneid1)[is.element(geneid1, gid.uniq)]
  ## which are the best probe for each gene
  js <- data.frame(js, "best"=FALSE)
  js[gg.uniq, "best"] <- TRUE
  ## data for duplicated gene ids
  if(length(gid.dupl) > 0) {      
    ## use jetset oevrall score to select the best probeset
    myscore <- js[gg.dupl,"overall"]
    myscore <- cbind("probe"=gg.dupl, "gid"=geneid1[gg.dupl], "score"=myscore)
    myscore <- myscore[order(as.numeric(myscore[ , "score"]), decreasing=TRUE, na.last=TRUE), , drop=FALSE]
    myscore <- myscore[!duplicated(myscore[ , "gid"]), , drop=FALSE]
    js[myscore[ ,"probe"], "best"] <- TRUE        
  }           
  
  #Reinscribing the gselist annotated frame with proper probe-gene mapping
  index <- which(js$best==TRUE) 
  probe <- rownames(js)[index]
  gene.id <- js$EntrezID[index]
  gene.symbol <- js$symbol[index]        
  mapping <- cbind(gene.id,gene.symbol,probe)
  featureData <- new("AnnotatedDataFrame", data=data.frame(mapping), varMetadata=data.frame(labelDescription=c("ENTREZID" , "SYMBOL", "PROBE")))
  Biobase::featureData(gse) <- featureData
  index <- match(probe,rownames(exprs(gse)))
  matrix.exprs <- exprs(gse)[index,]
  rownames(matrix.exprs) <- paste("geneid.", gene.id, sep="")
  colnames(matrix.exprs) <- colnames(exprs(gse))
  exprs(gse) <- matrix.exprs 
  
  gse_output <- gse
  return(gse_output)
}