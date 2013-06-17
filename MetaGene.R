MetaGene <- function(genes, data, annot, coefficients, method=c("principal.component", "weighted.average"), na.rm=FALSE) {
  #Take a specify gene set and create a single metagene base it's the selected method
  #
  # Args:
  #   genes: A vector of gene EntrezID that will create the metagene
  #   data: A matrix of gene expression        
  #   coefficients: A vector containing the coefficient that caracterize the gene from the geneset
  #
  # Returns:     
  #     the metagene
  
  
  method <- match.arg(method)
  n0 <- length(genes)
  genes2 <- intersect(genes, colnames(data))
  coefficients <- coefficients[match(genes2,genes)]
  genes <- genes2
  #if(length(genes) < n0) { warning(sprintf("%i genes from geneset are not present in the gene expression dataset", n0 - length(genes))) }
  ccix <- complete.cases(data[ , genes])
  nn <- sum(ccix)
  #if(nn > 0 && !na.rm) { stop("missing values are present") }  
  switch(method,
         "principal.component"= {
           metagene <- rep(NA, nrow(data))
           names(metagene) <- rownames(data)
           if(n0 > 0 && nn >= 3) {
             ## compute the first principal component
             prm <- prcomp(x=data[ccix, genes, drop = FALSE],na.action = na.omit)
             for (i in 1:nrow(prm[[2]])){
               if ((prm[[2]][i,"PC1"] < 0 && coefficient[i] > 0) || (prm[[2]][i,"PC1"] > 0 && coefficient[i] < 0)){
                 prm[[2]][i,"PC1"] <- prm[[2]][i,"PC1"]*(-1)
               }
             }
             temp <- prm$x[ , 1]
             metagene[names(temp)] <- temp             
           }
         },
         "weighted.average"={
           probe <- genes
           gene.id <- sapply(genes, function(x){strsplit(x, "\\.")[[1]][2]})
           gene.mat <- cbind(probe, gene.id, coefficients)
           colnames(gene.mat) <- c("probe","EntrezGene.ID","coefficient")
           metagene <- sig.score(gene.mat, data=data, annot=annot, do.mapping = FALSE, mapping, size = 0, cutoff = NA, signed = TRUE, verbose = FALSE)
           
         })
  return(metagene)
}








