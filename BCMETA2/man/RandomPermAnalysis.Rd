\name{RandomPermAnalysis}
\alias{RandomPermAnalysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Random Gene Permutation Analysis
}
\description{
The function seek out to calculate the random P value. It will generated a metagene out of the user's gene set and compare its prognostic P value against a 1000 (or n.perm) random gene set. It will then calculate how much random gene set have a better P value than the user's gene set and hand out what we call the significance P value or random P value.
}
\usage{
RandomPermAnalysis(master.eset, inputfile, method, n.perm)
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
  \item{n.perm}{
number of random gene set to test against the user's gene set.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

\item{metagene.ranking.list }{a list containing, for every molecular subtype, a matrix which itself contains 3 columns: metagene Corcondance Index, metagene Standard Deviation and metagene P value}
\item{metagene.mat }{a list containing, for every molecular subtype, a matrix of gene expression including the metagene's expression}
\item{prognostic.list}{A summary of how well the user metagene perform against the 1000 (or n.perm) random gene sets}
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
function (master.eset, inputfile, method, n.perm) 
{
    source(file.path("/stockage/homes/bachanp/Function/MetaGene.R"))
    metagene.rank.masterlist <- metagene.mat.list <- prognostic.list <- list()
    data <- t(exprs(master.eset))
    annot <- Biobase::featureData(master.eset)
    gene.mat.all <- read.m.file(inputfile)
    for (k in 1:length(gene.mat.all)) {
        gene.mat <- gene.mat.all[[k]]
        sample.genes <- paste("geneid.", gene.mat$EntrezGene.ID, 
            sep = "")
        sample.coef <- gene.mat$coefficient
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
        randgene.length <- length(intersect(sample.genes, colnames(data)))
        library(parallel)
        nbcore <- 8
        availcore <- detectCores()
        if (nbcore > availcore) {
            nbcore <- availcore
        }
        options(mc.cores = nbcore)
        splitix <- splitIndices(nx = n.perm, ncl = nbcore)
        if (randgene.length != 1) {
            rand.metagene <- mclapply(splitix, function(splitix2, 
                ...) {
                sub.rand.metagene <- lapply(splitix2, function(x) {
                  rand.genes <- rownames(exprs(master.eset))[sample(1:nrow(exprs(master.eset)), 
                    randgene.length, replace = F)]
                  subsub.rand.metagene <- MetaGene(genes = rand.genes, 
                    data = data, annot = annot, coefficients = sample.coef, 
                    method = method)[[1]]
                  return(subsub.rand.metagene)
                })
                sub.rand.metagene <- do.call(rbind, sub.rand.metagene)
            })
            rand.metagene <- do.call(rbind, rand.metagene)
            colnames(rand.metagene) <- colnames(exprs(master.eset))
        }
        else {
            index <- sample(1:nrow(exprs(master.eset)), n.perm, 
                replace = F)
            if (any(index == which(rownames(exprs(master.eset)) == 
                sample.genes))) {
                index <- index[-which(rownames(exprs(master.eset)) == 
                  sample.genes)]
                index <- c(index, 1)
            }
            rand.metagene <- exprs(master.eset)[index, ]
        }
        metagene.mat <- rbind(sample.metagene, rand.metagene)
        index <- which(colnames(pData(master.eset)) == "subtype")
        subtype.name <- colnames(pData(master.eset))[(index + 
            1):ncol(pData(master.eset))]
        pheno.data <- pData(master.eset)
        exprs.matrix <- exprs(master.eset)
        stime <- as.numeric(as.vector(pheno.data$surv.time))
        sevent <- as.numeric(as.vector(pheno.data$surv.event))
        index1 <- which(colnames(pheno.data) == "subtype")
        strat <- as.vector(pData(master.eset)$dataset)
        metagene.ranking.list <- list()
        nbcore <- 8
        availcore <- detectCores()
        if (nbcore > availcore) {
            nbcore <- availcore
        }
        options(mc.cores = nbcore)
        splitix <- splitIndices(nx = nrow(metagene.mat), ncl = nbcore)
        source(file.path("/stockage/homes/bachanp/Function/ConcordanceIndex3.R"))
        for (j in 1:length(subtype.name)) {
            weights <- as.numeric(as.vector(pheno.data[, index1 + 
                j]))
            metagene.ranking <- mclapply(splitix, function(x) {
                sub.metagene.ranking <- ConcordanceIndex3(exprs.matrix = metagene.mat[x, 
                  , drop = FALSE], stime = stime, sevent = sevent, 
                  strat = strat, weights = weights)
                return(sub.metagene.ranking)
            })
            metagene.ranking <- do.call(rbind, metagene.ranking)
            metagene.ranking.list[[j]] <- metagene.ranking
        }
        names(metagene.ranking.list) <- c("Global population", 
            "Lums", "Basal", "Her2", "LumB", "LumA")
        prog.mat <- NULL
        for (i in 1:length(metagene.ranking.list)) {
            prog.mat[i] <- sum(metagene.ranking.list[[i]][-1, 
                "p"] <= metagene.ranking.list[[i]][1, "p"])/n.perm
        }
        prog.mat <- cbind(names(metagene.ranking.list), prog.mat)
        metagene.mat.list[[k]] <- metagene.mat
        prognostic.list[[k]] <- prog.mat
        metagene.rank.masterlist[[k]] <- metagene.ranking.list
    }
    names(prognostic.list) <- names(metagene.rank.masterlist) <- names(metagene.mat.list) <- names(gene.mat.all)[1:k]
    return(list(metagene.rank.masterlist = metagene.rank.masterlist, 
        metagene.mat.list = metagene.mat.list, prognostic.list = prognostic.list))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ prognostic}
\keyword{  random }% __ONLY ONE__ keyword per line
