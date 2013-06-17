## function computing gene-drug associations
geneset.score <- function(genes, data, coeffcients, model=FALSE, method=c("principal.component", "weighted.average"), na.rm=FALSE) {
  method <- match.arg(method)
  n0 <- length(genes)
  genes <- intersect(genes, colnames(data))
  if(length(genes) < n0) { warning(sprintf("%i genes from geneset are not present in the gene expression dataset", n0 - length(genes))) }
  ccix <- complete.cases(data[ , genes])
  nn <- sum(ccix)
  if(nn > 0 && !na.rm) { stop("missing values are present") }
  cl <- rep(NA, nrow(data))
  names(cl) <- rownames(data)
  switch(method,
    "principal.component"= {
      if(n0 > 0 && nn >= 3) {
        ## compute the first principal component
        prm <- prcomp(x=data[ccix, genes, drop = FALSE])
        tt <- prm$x[ , 1]
        cl[names(tt)] <- tt
        if(model) { cl <- list("pc1"=cl, "model"=prm) }
      }
    },
    "weighted.average"={
      ## signature.score
  })
  return(cl)
}

genes <- paste("geneid.", read.m.file("signatures_human_orthologs.csv")[[1]]$EntrezGene.ID, sep="")
data <- t(exprs(master.eset))
metagene <- geneset.score(genes=genes, data=data, model=TRUE, na.rm=TRUE)
PC1 <- princomp(na.exclude(metagene))[[1]]

gene.mat <- as.matrix(read.m.file("signatures_human_orthologs.csv")[[1]])
gene.mat[,1] <- sapply(1:nrow(gene.mat), function(x){paste("geneid.", as.numeric(gene.mat[x,2]),sep="")})
annot <- Biobase::featureData(master.eset)
data <- t(exprs(master.eset))
metagene2 <- sig.score(gene.mat, data=data, annot=annot, do.mapping = FALSE, mapping, size = 0, cutoff = NA, signed = TRUE, verbose = FALSE)






