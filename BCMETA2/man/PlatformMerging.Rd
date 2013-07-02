\name{PlatformMerging}
\alias{PlatformMerging}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Platform merging
}
\description{
Some GSE expression matrix were measure for two differents microchips array. In order to account for these cases, we merge the expression matrix of the concern GSE dataset together, thus preserving all available information in the same dataset.
}
\usage{
PlatformMerging(gselist, GPL.length)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gselist}{
The gselist that is being constructed within the OpenDataset routine
}
  \item{GPL.length}{
The length of the annotation vector.
}
}
\details{
Note this function is not meant to be use outside of the OpenDataset routine. It's simply a subroutine of the latter.
}
\value{
return the merge GSE into the gselist.
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
function (gselist, GPL.length) 
{
    index <- length(gselist) - GPL.length + 1
    for (i in 1:GPL.length - 1) {
        missing.genes <- setdiff(rownames(exprs(gselist[[index + 
            i]])), rownames(exprs(gselist[[index]])))
        exprs(gselist[[index]]) <- rbind(exprs(gselist[[index]]), 
            exprs(gselist[[index + i]])[match(missing.genes, 
                rownames(exprs(gselist[[index + i]]))), , drop = FALSE])
        gene.id <- sapply(missing.genes, function(x) {
            strsplit(x, "\\.")[[1]][2]
        })
        index2 <- match(gene.id, Biobase::featureData(gselist[[index + 
            i]])@data[, 1])
        gene.annotation.matrix <- Biobase::featureData(gselist[[index + 
            i]])@data[index2, ]
        colnames(gene.annotation.matrix) <- colnames(Biobase::featureData(gselist[[index]])@data)
        Biobase::featureData(gselist[[index]])@data <- rbind(Biobase::featureData(gselist[[index]])@data, 
            gene.annotation.matrix)
    }
    return(gselist[[index]])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ merging }
\keyword{ platform }% __ONLY ONE__ keyword per line
