\name{SubtypeDependanceAnalysis}
\alias{SubtypeDependanceAnalysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subtype Dependence Analysis
}
\description{
This function verify if a gene set is subtype dependent. In order to do so, it firstly generates a metagene from the user's gene set using a given method. Then, it applies the Kruskal-Wallis test and the Wilcoxson test on the metagene to verify whether it is subtype dependent or not. 
}
\usage{
SubtypeDependanceAnalysis(master.eset, gene.list, inputfile, method)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{master.eset}{
the merged eSet structure (or it could be a normal eSet).
}
  \item{inputfile}{
 a CSV file containing the name of the gene set, the probe, the EntrezGeneID and the signature coefficient.
}
  \item{method}{
  It's a method use for generating the metagene. Either weigthed average or principal component
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

\item{temp }{The modified gene expression matrix containing the metagene as the first row}
\item{dep.out}{the dependence list containing every information on the subtype dependency}
%% ...
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
function (master.eset, gene.list, inputfile, method) 
{
    source(file.path("/stockage/homes/bachanp/Function/SubtypeDependance.R"))
    source(file.path("/stockage/homes/bachanp/Function/MetaGene.R"))
    metagene.matrix <- NULL
    gene.mat.all <- read.m.file(inputfile)
    for (i in 1:length(gene.mat.all)) {
        gene.mat <- gene.mat.all[[i]]
        data <- t(exprs(master.eset))
        annot <- Biobase::featureData(master.eset)
        sample.genes <- paste("geneid.", gene.mat$EntrezGene.ID, 
            sep = "")
        sample.coef <- gene.mat$coefficient
        method <- "weighted.average"
        temp <- master.eset
        if (length(sample.genes) != 1) {
            sample.metagene <- MetaGene(genes = sample.genes, 
                data = data, annot = annot, coefficients = sample.coef, 
                method = method)[[1]]
        }
        else {
            sample.metagene <- MetaGene(genes = sample.genes, 
                data = data, annot = annot, coefficients = sample.coef, 
                method = method)
        }
        metagene.matrix <- rbind(metagene.matrix, sample.metagene)
    }
    rownames(metagene.matrix) <- names(gene.mat.all)
    exprs(temp) <- rbind(metagene.matrix, exprs(temp))
    gene.list <- names(gene.mat.all)
    dep.output <- SubtypeDependance(master.eset = temp, gene.list = gene.list)
    return(list(dep.output = dep.output, temp = temp))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Subtype}
\keyword{ Dependence}% __ONLY ONE__ keyword per line
