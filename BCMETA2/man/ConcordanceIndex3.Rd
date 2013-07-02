\name{ConcordanceIndex3}
\alias{ConcordanceIndex3}
\title{
Concordance index.
}
\description{
Takes expression gene matrix as input along with an survival time and survival event vector and calculates the survival corcondance index for every gene in the expression gene matrix for each molecular subtype. Note that strat is an optionnal vector that enable the user to stratify for a subgroup of survival data.
}
\usage{
ConcordanceIndex3(exprs.matrix, stime, sevent, strat, weights)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exprs.matrix}{
Is the expression gene matrix that contains every gene expression. The row names should be the geneid.#### and the column names should be the patient ID.
}
  \item{stime}{
A vector of survival time (either dmfs, rfs, os or the hybrid dfs)
}
  \item{sevent}{
A survival event vector containing only 0 or 1 for censoring
}
  \item{strat}{
A stratification subset. Usely we stratified by dataset 
}
  \item{weights}{
A vector containing patient molecular subtype weight. Predicted using genefu::subtype.cluster.predict or genefu::intrisic.cluster.predict
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Return a matrix containing 3 column (Corcondance Index, Standard Deviation, p.value) for every gene in the expression matrix.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Pierre-Olivier Bachant-Winner
Benjamin Haibe-Kains
}
\note{
This is a parallelize function that use multiple core to accelerated computation.
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
function (exprs.matrix, stime, sevent, strat, weights) 
{
    require(survcomp)
    if (missing(strat)) {
        strat <- rep(1, ncol(exprs.matrix))
        names(strat) <- colnames(exprs.matrix)
    }
    res <- apply(exprs.matrix, 1, function(x, stime, sevent, 
        strat, weights) {
        cc <- survcomp::concordance.index(x = x, surv.time = stime, 
            surv.event = sevent, weights = weights, strat = strat, 
            outx = TRUE, method = "noether", alternative = "two.sided", 
            na.rm = TRUE)
        return(c(cindex = cc$c.index, se = cc$se, p = cc$p.value))
    }, stime = stime, sevent = sevent, strat = strat, weights = weights)
    return(t(res))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Corcondance index }
\keyword{ Survival }% __ONLY ONE__ keyword per line
