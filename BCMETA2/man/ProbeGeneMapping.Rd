\name{ProbeGeneMapping}
\alias{ProbeGeneMapping}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Probe Gene Mapping
}
\description{
Take a GSE (eSet input) that contains an expression matrix with probe as row names (probe expression matrix), and seek out the must variant one and apply the corresponding gene id. 
}
\usage{
ProbeGeneMapping(gse_input)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gse_input}{
An input eSet structure that contains an expression matrix with probes in row names.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Hand out the GSE (or eSet structure) with the gene expression matrix.
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
function (gse_input) 
{
    gse <- gse_input
    js <- jetset.bhk::jscores(chip = annotation(gse), probeset = rownames(exprs(gse)))
    js <- js[rownames(exprs(gse)), , drop = FALSE]
    rownames(js) <- rownames(exprs(gse))
    geneid1 <- as.character(js[, "EntrezID"])
    names(geneid1) <- rownames(js)
    geneid2 <- sort(unique(geneid1))
    names(geneid2) <- paste("geneid", geneid2, sep = ".")
    gix1 <- !is.na(geneid1)
    gix2 <- !is.na(geneid2)
    geneid.common <- intersect(geneid1[gix1], geneid2[gix2])
    gg <- names(geneid1)[is.element(geneid1, geneid.common)]
    gid <- geneid1[is.element(geneid1, geneid.common)]
    gid.dupl <- unique(gid[duplicated(gid)])
    gg.dupl <- names(geneid1)[is.element(geneid1, gid.dupl)]
    gid.uniq <- gid[!is.element(gid, gid.dupl)]
    gg.uniq <- names(geneid1)[is.element(geneid1, gid.uniq)]
    js <- data.frame(js, best = FALSE)
    js[gg.uniq, "best"] <- TRUE
    if (length(gid.dupl) > 0) {
        myscore <- js[gg.dupl, "overall"]
        myscore <- cbind(probe = gg.dupl, gid = geneid1[gg.dupl], 
            score = myscore)
        myscore <- myscore[order(as.numeric(myscore[, "score"]), 
            decreasing = TRUE, na.last = TRUE), , drop = FALSE]
        myscore <- myscore[!duplicated(myscore[, "gid"]), , drop = FALSE]
        js[myscore[, "probe"], "best"] <- TRUE
    }
    index <- which(js$best == TRUE)
    probe <- rownames(js)[index]
    gene.id <- js$EntrezID[index]
    gene.symbol <- js$symbol[index]
    mapping <- cbind(gene.id, gene.symbol, probe)
    featureData <- new("AnnotatedDataFrame", data = data.frame(mapping), 
        varMetadata = data.frame(labelDescription = c("ENTREZID", 
            "SYMBOL", "PROBE")))
    Biobase::featureData(gse) <- featureData
    index <- match(probe, rownames(exprs(gse)))
    matrix.exprs <- exprs(gse)[index, ]
    rownames(matrix.exprs) <- paste("geneid.", gene.id, sep = "")
    colnames(matrix.exprs) <- colnames(exprs(gse))
    exprs(gse) <- matrix.exprs
    gse_output <- gse
    return(gse_output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ mapping }
\keyword{ probe }% __ONLY ONE__ keyword per line
