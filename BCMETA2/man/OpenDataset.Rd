\name{OpenDataset}
\alias{OpenDataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Open Dataset
}
\description{
This function allow user to interactively download dataset from the InsilicoDB database. Every database is downloaded as an eSet R structure. The function then create a list containing every eSet structure, while making sure that every eSet is properly curated. 
}
\usage{
OpenDataset(config.file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{config.file}{
The configuration file. It's allowed to user to selected which GSE he wants to include is his study
}
}
\details{
-All dataset in InsilicoDB have been previously curated using a very specific curation required for this pipeline to work.
-This function calls 2 subroutines: 1 for the probe-gene-mapping and the other one for merging together dataset that would have been made on different micro-chip platform.
}
\value{
Return the eSet list.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Pierre-Olivier Bachant-Winner
Benjamin Haibe-Kains
}
\note{
For some combination of GSE we have to call some routine hack to properly merge them. User discretion is advised.
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
function (config.file) 
{
    require(inSilicoDb2)
    require(genefu)
    require(jetset)
    source(file.path("/stockage/homes/bachanp/Function/ProbeGeneMapping.R"))
    source(file.path("/stockage/homes/bachanp/Function/PlatformMerging.R"))
    InSilicoLogin(login = "bhaibeka@gmail.com", password = "747779bec8a754b91076d6cc1f700831")
    FTO <- read.csv(config.file)
    gselist <- list()
    counter <- 1
    routine.hack1 <- NULL
    for (i in 16:nrow(FTO)) {
        if (as.character(FTO[i, 4]) == TRUE) {
            if (as.character(FTO[i, 2]) == "GSE2034" || as.character(FTO[i, 
                2]) == "GSE5327") {
                routine.hack1 <- rbind(routine.hack1, c(as.character(FTO[i, 
                  2]), counter))
            }
            GPL <- getPlatforms(dataset = as.character(FTO[i, 
                2]))
            temp <- getDatasets(dataset = as.character(FTO[i, 
                2]), norm = as.character(FTO[i, 6]), curation = FTO[i, 
                3], features = "PROBE")
            for (j in 1:length(GPL)) {
                gselist[[counter]] <- temp[[j]]
                if (as.character(FTO[i, 5]) == "Affy") {
                  gselist[[counter]] <- ProbeGeneMapping(gselist[[counter]])
                }
                else {
                  rownames(exprs(gselist[[counter]])) <- paste("geneid.", 
                    Biobase::featureData(gselist[[counter]])$ENTREZID, 
                    sep = "")
                  index <- match(pData(gselist[[counter]])$id, 
                    colnames(exprs(gselist[[counter]])))
                  exprs(gselist[[counter]]) <- exprs(gselist[[counter]])[, 
                    index]
                  rownames(pData(gselist[[counter]])) <- as.character(pData(gselist[[counter]])$id)
                }
                counter <- counter + 1
            }
            if (length(GPL) > 1) {
                index <- length(gselist) - (length(GPL) - 1)
                gselist[[index]] <- PlatformMerging(gselist = gselist, 
                  GPL.length = length(GPL))
                for (k in 1:length(gselist)) {
                  gselist[[index + k]] <- NULL
                }
                counter <- counter - (length(GPL) - 1)
            }
        }
    }
    if (!is.null(routine.hack1)) {
        if (nrow(routine.hack1) == 2) {
            exprs(gselist[[as.numeric(routine.hack1[1, 2])]]) <- cbind(exprs(gselist[[as.numeric(routine.hack1[1, 
                2])]]), exprs(gselist[[as.numeric(routine.hack1[2, 
                2])]]))
            phenoData(gselist[[as.numeric(routine.hack1[1, 2])]])@data <- rbind(phenoData(gselist[[as.numeric(routine.hack1[1, 
                2])]])@data, phenoData(gselist[[as.numeric(routine.hack1[2, 
                2])]])@data)
            gselist[[as.numeric(routine.hack1[2, 2])]] <- NULL
        }
    }
    for (i in 1:length(gselist)) {
        index = match(colnames(exprs(gselist[[i]])), rownames(pData(gselist[[i]])))
        phenoData(gselist[[i]])@data <- phenoData(gselist[[i]])@data[index, 
            ]
        if (any(colnames(pData(gselist[[i]])) == "Anatomical site")) {
            index <- which(pData(gselist[[i]])$Anatomical == 
                "breast")
            phenoData(gselist[[i]])@data <- phenoData(gselist[[i]])@data[index, 
                ]
            exprs(gselist[[i]]) <- exprs(gselist[[i]])[, index]
        }
        temp <- pData(gselist[[i]])
        temp$platform <- NULL
        for (j in 1:nrow(temp)) {
            if (all(temp[j, ] == "NILL")) {
                exprs(gselist[[i]]) <- exprs(gselist[[i]])[, 
                  -j]
                phenoData(gselist[[i]])@data <- phenoData(gselist[[i]])@data[-j, 
                  ]
            }
        }
        phenoData(gselist[[i]])@data <- phenoData(gselist[[i]])@data[, 
            match(colnames(phenoData(gselist[[1]])@data), colnames(phenoData(gselist[[i]])@data))]
        colnames(Biobase::featureData(gselist[[i]])@data) <- c("ENTREZID", 
            "SYMBOL", "PROBE")
    }
    InSilicoLogout()
    return(gselist)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ GSE }
\keyword{ InsilicoDB }% __ONLY ONE__ keyword per line
