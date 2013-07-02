ConcordanceIndex3 <-
function(exprs.matrix, stime, sevent, strat, weights){
  # Compute the concordance index for each individual gene
  # based on the survcomp package
  #
  # Args:
  #  gene.id: The gene of importance to extract C index
  #  exprs.matrix: Matrix containing all the gene expression for the dataset
  #  stime: survival time
  #  sevent: survival event
  #  strat: strata factors (usually datasets)
  #  weights: sample weights (subtype probabilities for instance)
  #
  # Returns:
  #  a matrix containing the concordance index, its standard error and corresponding two-sided p-value 
  
  require(survcomp)
  if(missing(strat)) {
    strat <- rep(1, ncol(exprs.matrix))
    names(strat) <- colnames(exprs.matrix)
  }
  res <- apply(exprs.matrix, 1, function(x, stime, sevent, strat, weights) {
    cc <- survcomp::concordance.index(x=x, surv.time=stime, surv.event=sevent, weights=weights, strat=strat, outx=TRUE, method="noether", alternative="two.sided", na.rm=TRUE)
    return(c("cindex"=cc$c.index, "se"=cc$se, "p"=cc$p.value))
  }, stime=stime, sevent=sevent, strat=strat, weights=weights) 
  return(t(res))
}
