\name{GseaAnalysis}
\alias{GseaAnalysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GSEA analysis
}
\description{
A function that calls the java GSEA tool. It required a gmt input file containing all the gene that the user which to analyse. The function then call the java tool and run the enrichement analysis for every molecular subtype and test the user's gene set against a collection of 1000 random gene set enrichement scores for significance handing out the enrichement p value.
}
\usage{
GseaAnalysis(survival.score.list, gmt.file, matrix.subtype)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{survival.score.list}{
A survival list where each list element is a matrix that must contain 3 columns (Corcondance index, Standard deviation and p.value) for every gene. Note that each list element as a different matrix associated with molecular subtype.
}
  \item{gmt.file}{
The gmt file containing the user gene set.
}
  \item{matrix.subtype}{
A simple matrix containing every molecular subtype weight information for every patient.
}
}
\details{
This function required many pathway information in order to work properly. Must importanly the user must have in is disposition the GSEA java tool (gsea2-2.0.12.jar) in his working directory. The function will generated many folder (one for each molecular subtype) containing the GSEA analysis. If the user which to rerun the analysis he must absolutely deleted the previously generated folder. Otherwise, the java tool will return and error saying that the file already exist and will not remake the analysis.
}
\value{
All the information should be in the generated folder.
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
function (survival.score.list, gmt.file, matrix.subtype) 
{
    if (missing(gmt.file)) {
        gmt.file <- "c2.all.v3.1.entrez.gmt"
    }
    working.directory <- getwd()
    dir.create(file.path(working.directory, "saveres2", "ranking"), 
        recursive = TRUE, showWarnings = FALSE)
    rank.files <- NULL
    for (i in 1:length(survival.score.list)) {
        ranking.matrix <- survival.score.list[[i]]
        gene.id <- sapply(rownames(survival.score.list[[i]]), 
            function(x) {
                strsplit(x, "\\.")[[1]][2]
            })
        score <- -log10(survival.score.list[[i]][, 3]) * sign(survival.score.list[[i]][, 
            1] - 0.5)
        ranking <- cbind(gene.id, score)
        colnames(ranking) <- c("EntrezGene.ID", "Survival score")
        temp <- i
        filepath = file.path(working.directory, "saveres2", "ranking", 
            paste("ranking_", sprintf("\%s.rnk", temp), sep = ""))
        write.table(ranking, file = filepath, col.names = FALSE, 
            row.names = FALSE, sep = "\t", quote = FALSE)
        rank.files <- c(rank.files, filepath)
    }
    dir.create(file.path(working.directory, "saveres2", "GSEA", 
        "reports"), recursive = TRUE, showWarnings = FALSE)
    source(file.path(working.directory, "foo.R"))
    gsea.exec <- file.path(working.directory, "gsea2-2.0.12.jar")
    gsea.nperm <- 1000
    min.geneset.size <- 15
    max.geneset.size <- 500
    genesets.filen <- file.path(working.directory, gmt.file)
    gsea.out <- file.path(working.directory)
    nbcore <- length(survival.score.list)
    availcore <- detectCores()
    if (nbcore > availcore) {
        nbcore <- availcore
    }
    options(mc.cores = nbcore)
    splitix <- splitIndices(nx = length(survival.score.list), 
        ncl = nbcore)
    gsea.res <- mclapply(splitix, function(splitix2, ...) {
        gseatt <- gsea.prerank(exe.path = gsea.exec, gmt.path = genesets.filen, 
            rank.path = rank.files[splitix2], gsea.collapse = FALSE, 
            nperm = gsea.nperm, scoring.scheme = "weighted", 
            make.sets = TRUE, include.only.symbols = FALSE, plot.top.x = 20, 
            set.max = max.geneset.size, set.min = min.geneset.size, 
            zip.report = FALSE, gsea.out = gsea.out, replace.res = FALSE, 
            gsea.seed = 987654321)
        tt <- gseatt[[1]]
        tt2 <- gseatt[[2]][!sapply(gseatt[[2]], is.null)]
        tt2 <- lapply(tt2, function(x) {
            return(paste("geneid", x, sep = "."))
        })
        return(list(geneset.res = tt, geneset.core = tt2))
    }, rank.files = rank.files, gsea.exec = gsea.exec, genesets.filen = genesets.filen, 
        gsea.nperm = gsea.nperm, max.geneset.size = max.geneset.size, 
        min.geneset.size = min.geneset.size, gsea.out = gsea.out)
    names(gsea.res) <- names(rank.files)
    return(gsea.res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ GSEA }
\keyword{ enrichement }% __ONLY ONE__ keyword per line
