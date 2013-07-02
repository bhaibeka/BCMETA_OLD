\name{SubtypeClassification}
\alias{SubtypeClassification}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SubtypeClassification(gselist, config.file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gselist}{
%%     ~~Describe \code{gselist} here~~
}
  \item{config.file}{
%%     ~~Describe \code{config.file} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
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
function (gselist, config.file) 
{
    library(genefu)
    library(Biobase)
    FTO <- read.csv(config.file)
    model <- as.character(FTO[3, 2])
    accepted.model <- c("ssp2003.robust", "ssp2006.robust", "pam50", 
        "scmgene.robust", "scmod1.robust", "scmod2.robust")
    if (!any(model == accepted.model)) {
        stop("The model is invalid")
    }
    STL <- list()
    for (i in 1:length(gselist)) {
        data <- t(exprs(gselist[[i]]))
        colnames(data) <- rownames(exprs(gselist[[i]]))
        mapping <- matrix(cbind(colnames(data), as.vector(Biobase::featureData(gselist[[i]])@data[, 
            1])), length(rownames(exprs(gselist[[i]]))), 2)
        colnames(mapping) <- c("probe", "EntrezGene.ID")
        rownames(mapping) <- mapping[, "probe"]
        annot <- as.matrix(Biobase::featureData(gselist[[i]])@data)
        colnames(annot)[1] <- "EntrezGene.ID"
        accepted.intrinsic.model <- c("ssp2003.robust", "ssp2006.robust", 
            "pam50")
        if (any(model == accepted.intrinsic.model)) {
            subtype <- intrinsic.cluster.predict(sbt.model = eval(parse(text = model)), 
                data = data, annot = mapping, do.mapping = TRUE, 
                do.prediction.strength = FALSE)
            matrix_sub <- cbind(subtype$subtype, rep(1, nrow(data)), 
                subtype$subtype.proba[, 3] + subtype$subtype.proba[, 
                  4], subtype$subtype.proba)
            colnames(matrix_sub) <- c("subtype", "Global population", 
                "Lums", "Basal", "Her2", "LumB", "LumA")
        }
        else {
            subtype <- subtype.cluster.predict(sbt.model = eval(parse(text = model)), 
                data = data, annot = mapping, do.mapping = TRUE, 
                do.prediction.strength = FALSE, do.BIC = FALSE, 
                plot = FALSE)
            matrix_sub <- cbind(subtype$subtype2, rep(1, nrow(data)), 
                subtype$subtype.proba2[, 3] + subtype$subtype.proba2[, 
                  4], subtype$subtype.proba2)
            colnames(matrix_sub) <- c("subtype", "Global population", 
                "Lums", "Basal", "Her2", "LumB", "LumA")
        }
        rownames(matrix_sub) <- colnames(exprs(gselist[[i]]))
        STL[[i]] <- matrix_sub
    }
    return(STL)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
