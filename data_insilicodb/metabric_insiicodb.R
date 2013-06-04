rm(list=ls())

library(survcomp)
library(genefu)

badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

saveres <- file.path("METABRIC")
if(!file.exists(saveres)) { dir.create(saveres, recursive=TRUE, showWarnings=FALSE) }
saveres2 <- file.path(saveres, "esets")
if(!file.exists(saveres2)) { dir.create(saveres2, recursive=TRUE, showWarnings=FALSE) }
      
load("METABRIC_GE_ENTREZ.RData")

## reformat data
# data <- rbind(tumor.ge, normal.ge)
data <- tumor.ge
rownames(data) <- gsub(badchars, "", rownames(data))

## reformat demo
clinn <- c("id", "dataset" , "series", "age", "size", "node", "er", "pgr", "her2", "grade", "t.rfs", "e.rfs", "t.dmfs", "e.dmfs", "t.os", "e.os", "treatment", "tissue")
# demo <- data.frame(matrix(NA, nrow=nrow(tumor.info) + nrow(normal.info), ncol=length(clinn), dimnames=list(c(rownames(tumor.info), rownames(normal.info)), clinn)))
demo <- data.frame(matrix(NA, nrow=nrow(tumor.info), ncol=length(clinn), dimnames=list(c(rownames(tumor.info)), clinn)))
demo[rownames(tumor.info), intersect(clinn, colnames(tumor.info))] <- tumor.info[ , intersect(clinn, colnames(tumor.info))]
demo[rownames(tumor.info), "t.rfs"] <- tumor.info[ , "t.dfs"]
demo[rownames(tumor.info), "e.rfs"] <- tumor.info[ , "e.dfs"]
demo[rownames(tumor.info), "treatment"] <- tumor.info[ , "Treatment"]
demo[rownames(tumor.info), "tissue"] <- "TUMOR"
demo[rownames(tumor.info), "series"] <- tumor.info[ , "Site"]
demo[ , "dataset"] <- "METABRIC"
# demo[rownames(normal.info), intersect(clinn, colnames(normal.info))] <- normal.info[ , intersect(clinn, colnames(normal.info))]
# demo[rownames(normal.info), "tissue"] <- "NORMAL"
rownames(demo) <- gsub(badchars, "", rownames(demo))
demo[ , "id"] <- rownames(demo)

## refomrat annot
annot <- data.frame("ENTREZID"=as.character(annot.ge[ , "EntrezID"]), "SYMBOL"=as.character(annot.ge[ , "Symbol"]), "GENENAME"=as.character(annot.ge[ , "nuID"]))


write.csv(annot, file=file.path(saveres, "annot_metabric.csv"), row.names=FALSE)
write.csv(demo, file=file.path(saveres, "demo_metabric.csv"), row.names=FALSE)
write.csv(t(data), file=file.path(saveres, "data_metabric.csv"))

## create individual eSets
ddn <- "METABRIC"

message(sprintf("Read %s data", ddn))
## read clinical info
sampleinfo <- read.csv(file.path(toupper(ddn), sprintf("demo_%s.csv", tolower(ddn))), stringsAsFactors=FALSE)
rownames(sampleinfo) <- gsub(badchars, "", sampleinfo[ , "id"])
## read annotations
annotge <- read.csv(file.path(toupper(ddn), sprintf("annot_%s.csv", tolower(ddn))), stringsAsFactors=FALSE)
rownames(annotge) <- annotge[ , "GENENAME"]
## read gene expression data
datage <- read.csv(file.path(toupper(ddn), sprintf("data_%s.csv", tolower(ddn))), stringsAsFactors=FALSE, row.names=1)
rownames(datage) <- rownames(annotge)
colnames(datage) <- gsub(badchars, "", colnames(datage))
datage <- data.matrix(datage)

## create individual eSets
for(i in 1:nrow(sampleinfo)) {
  message(sprintf("Creation of an eSet for sample %s", rownames(sampleinfo)[i]))
  eseti <- ExpressionSet(assayData=datage[ , i, drop=FALSE], phenoData=AnnotatedDataFrame(data=sampleinfo[i, , drop=FALSE]), featureData=AnnotatedDataFrame(annotge))
  experimentData(eseti)@preprocessing <- list("normalization"="ORIGINAL", package="illuminaHT12v3", version="0")
  annotation(eseti) <- "illuminaHT12v3"
  ## object name
  objn <- sprintf("%sORIGINALPROBE", rownames(sampleinfo)[i])
  assign(objn, eseti)
  save(list=objn, compress=TRUE, file=file.path(toupper(ddn), "esets", sprintf("%s.RData", objn)))
}
