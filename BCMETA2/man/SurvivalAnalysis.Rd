\name{SurvivalAnalysis}
\alias{SurvivalAnalysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Survival Analysis
}
\description{
The function calculates the corcondance index, standard deviation and hand out prognostic value for every gene in the specified eSet's matrix expression.
}
\usage{
SurvivalAnalysis(master.eset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{master.eset}{
the merged eSet structure (or it could be a normal eSet).}
}
\details{
This function is parallelize to enhance computational performance.
}
\value{
It returns the survival.score.list, which is a list containing genome wide Concordance Index, Standard Deviation and P value for every molecular subtype.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Pierre-Olivier Bachant-Winner
Benjamin Haibe-Kains
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (master.eset) 
{
    rm(list = setdiff(ls(), "master.eset"))
    index <- which(colnames(pData(master.eset)) == "subtype")
    subtype.name <- colnames(pData(master.eset))[(index + 1):ncol(pData(master.eset))]
    survival.score.list <- list()
    for (j in 1:length(subtype.name)) {
        pheno.data <- pData(master.eset)
        exprs.matrix <- exprs(master.eset)
        stime <- as.numeric(as.vector(pheno.data$surv.time))
        sevent <- as.numeric(as.vector(pheno.data$surv.event))
        index <- which(colnames(pheno.data) == "subtype")
        weights <- as.numeric(as.vector(pheno.data[, index + 
            j]))
        strat <- as.vector(pData(master.eset)$dataset)
        library(parallel)
        nbcore <- 16
        availcore <- detectCores()
        if (nbcore > availcore) {
            nbcore <- availcore
        }
        options(mc.cores = nbcore)
        splitix <- splitIndices(nx = nrow(exprs.matrix), ncl = nbcore)
        source(file.path("/stockage/homes/bachanp/Function/ConcordanceIndex3.R"))
        survival.ranking <- mclapply(splitix, function(splitix2, 
            ...) {
            gene.ranking <- ConcordanceIndex3(exprs.matrix = exprs.matrix[splitix2, 
                , drop = FALSE], stime = stime, sevent = sevent, 
                strat = strat, weights = weights)
            return(gene.ranking)
        }, stime = stime, exprs.matrix = exprs.matrix, sevent = sevent, 
            strat = strat, weights = weights)
        survival.ranking <- do.call(rbind, survival.ranking)
        survival.score.list[[j]] <- survival.ranking
    }
    return(survival.score.list)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
