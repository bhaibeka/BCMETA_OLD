SurvivalAnalysis <-
function(master.eset){
  #Open every GSE, specify in the configuration CSV, from InsilicoDB and create a list of Eset
  #structure call gselist
  #
  # Args:
  #   master.eset: The huge eset that containt everything
  #  
  # Returns:     
  #     The list containing all ranking file 
  
  ##################################################################
  ##### Precomputing every possible data to enhance time performance
  ##################################################################
  rm(list=setdiff(ls(), "master.eset"))
  index <- which(colnames(pData(master.eset))=="subtype")
  subtype.name <- colnames(pData(master.eset))[(index+1):ncol(pData(master.eset))]
  #"subtype", "Global population", "Lums", "Basal", "Her2", "LumB", "LumA"
  survival.score.list <- list()
  for (j in 1:length(subtype.name)){ 
  
    pheno.data <- pData(master.eset)
    exprs.matrix <- exprs(master.eset)           
    stime <- as.numeric(as.vector(pheno.data$surv.time))
    sevent <- as.numeric(as.vector(pheno.data$surv.event))
    index <- which(colnames(pheno.data)=="subtype")
    weights <- as.numeric(as.vector(pheno.data[,index+j]))   
    strat <- as.vector(pData(master.eset)$dataset)
         
    #######################################
    ##### Parallelizing survival function
    #######################################
    library(parallel)
      
    #Setting every parallelizing parameter
    nbcore <- 16
    availcore <- detectCores()
    if(nbcore > availcore) { nbcore <- availcore }
    options("mc.cores"=nbcore)
    splitix <- splitIndices(nx=nrow(exprs.matrix), ncl=nbcore)
    source(file.path("/stockage/homes/bachanp/Function/ConcordanceIndex3.R"))
      
    #Main survival core function
    survival.ranking <- mclapply(splitix, function(splitix2, ...){      
      gene.ranking <- ConcordanceIndex3(exprs.matrix=exprs.matrix[splitix2, , drop=FALSE], stime=stime, sevent=sevent, strat=strat, weights=weights)      
      return(gene.ranking)
    }, stime=stime, exprs.matrix=exprs.matrix, sevent=sevent, strat=strat, weights=weights) 
      
    #Retrieving parallel output
    survival.ranking <- do.call(rbind, survival.ranking)
    #survival.ranking <- survival.ranking[, , drop=FALSE]  
    survival.score.list[[j]] <- survival.ranking  
  }
  return(survival.score.list)
}
