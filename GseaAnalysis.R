GseaAnalysis <- function(survival.score.list, gmt.file, matrix.subtype){
  #Take a survival score matrix containing all the ranking score for each gene of the geneset 
  #you want to make the GSEA analysis
  #
  # Note: In order for this function to work properly, you need to have the following in your
  #       working directory:- gsea2-2.0.12.jar
  #                         - foo.R
  #                         - your specify gmt file or c2.all.v3.1.entrez.gmt (by default)
  #
  # Args:
  #   survival.score.list: A list containing many survival score matrix. Where each matrix is created
  #                        such a way that the first column of each matrix containt the gene's name
  #                        the second column containt the survival score and the last column containt 
  #                        the p value.
  #   gmt.file: The GMT file that you which to compare your genes with
  #   matrix.subtype: To know which geneset is associated to which substype
  #
  # Returns:
  #   In return the GSEA analysis containing all the useful informations in the specify pathway      
  if(missing(gmt.file)) {gmt.file <- "c2.all.v3.1.entrez.gmt"}
  
  gmt.file <- "c2.all.v3.1.entrez.gmt"  
  
  #In order to create the pathway for the survival score rnk file.
  working.directory <- getwd()
  dir.create(file.path(working.directory,"saveres2", "ranking"), recursive=TRUE, showWarnings=FALSE)
  rank.files <- NULL
  for (i in 1:length(survival.score.list)){
    ranking.matrix <- survival.score.list[[i]]
    gene.id <- sapply(rownames(survival.score.list[[i]]), function(x){strsplit(x, "\\_")[[1]][2]})
    ranking <- cbind(gene.id, ranking.matrix[,1])
    colnames(ranking) <- c("EntrezGene.ID","Survival score")
    #temp <- colnames(matrix.subtype)[i+1]
    temp <- i
    filepath=file.path(working.directory,"saveres2","ranking",paste("ranking_", sprintf("%s.rnk",temp),sep=""))
    write.table(ranking, file=filepath, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
    rank.files <- c(rank.files, filepath)
  }  
  
  #Initializing gsea.prerank function path and parameters
  dir.create(file.path(working.directory, "saveres2", "GSEA", "reports"), recursive=TRUE, showWarnings=FALSE)
  source(file.path(working.directory, "foo.R"))
  gsea.exec <- file.path(working.directory, "gsea2-2.0.12.jar")
  gsea.nperm <- 10
  min.geneset.size <- 15
  max.geneset.size <- 500
  genesets.filen <- file.path(working.directory,gmt.file)
  gsea.out <- file.path(working.directory)    
  
  #Applying GSEA function on all the file 
  #Setting every parallelizing parameter
  nbcore <- 4
  availcore <- detectCores()
  if(nbcore > availcore) { nbcore <- availcore }
  options("mc.cores"=nbcore)
  splitix <- splitIndices(nx=length(survival.score.list), ncl=nbcore)
  
  gsea.res <- mclapply(splitix, function(splitix2, ...) { 
    gseatt <- gsea.prerank(exe.path=gsea.exec, gmt.path=genesets.filen, rank.path=rank.files[splitix2], gsea.collapse=FALSE, nperm=gsea.nperm, scoring.scheme="weighted", make.sets=TRUE, include.only.symbols=FALSE, plot.top.x=length(nrow(ranking.matrix)), set.max=max.geneset.size, set.min=min.geneset.size, zip.report=FALSE, gsea.out=gsea.out, replace.res=FALSE, gsea.seed=987654321)
    ## results
    tt <- gseatt[[1]]  
    ## leading edge genes
    tt2 <- gseatt[[2]][!sapply(gseatt[[2]], is.null)]
    tt2 <- lapply(tt2, function(x) { return(paste("geneid", x, sep="_")) })
    return(list("geneset.res"=tt, "geneset.core"=tt2))
  }, rank.files=rank.files, gsea.exec=gsea.exec, genesets.filen=genesets.filen, gsea.nperm=gsea.nperm, max.geneset.size=max.geneset.size, min.geneset.size=min.geneset.size, gsea.out=gsea.out)
  names(gsea.res) <- names(rank.files)
  return(gsea.res)
}

