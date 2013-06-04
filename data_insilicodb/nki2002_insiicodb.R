
library(survcomp)
library(genefu)

badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

saveres <- file.path("NKI")
if(!file.exists(saveres)) { dir.create(saveres, recursive=TRUE, showWarnings=FALSE) }
saveres2 <- file.path(saveres, "esets")
if(!file.exists(saveres2)) { dir.create(saveres2, recursive=TRUE, showWarnings=FALSE) }
      
load("nki2002.RData")
rownames(data) <- gsub(badchars, "", rownames(data))
rownames(demo) <- gsub(badchars, "", rownames(demo))

gid1 <- as.character(annot[ , "EntrezGene.ID"])
names(gid1) <- rownames(annot)
gid2 <- sort(unique(annot[ , "EntrezGene.ID"]))
names(gid2) <- paste("geneid", gid2, sep=".")
rr <- genefu::geneid.map(geneid1=gid1, data1=data, geneid2=gid2)

data <- rr$data1
annot <- annot[colnames(data), , drop=FALSE]
colnames(data) <- rownames(annot) <- names(rr$geneid2)

## refomrat annot
annot <- data.frame("ENTREZID"=as.character(annot[ , "EntrezGene.ID"]), "SYMBOL"=as.character(annot[ , "HUGO.gene.symbol"]), "GENENAME"=as.character(annot[ , "probe"]))

## reformat demo
demo <- demo[ , c("id", "dataset" , "series", "age", "size", "node", "er", "pgr", "her2", "grade", "t.rfs", "e.rfs", "t.dmfs", "e.dmfs", "t.os", "e.os", "treatment", "tissue"), drop=FALSE]
demo[ , "id"] <- rownames(demo)
demo[ , "tissue"] <- "TUMOR"
tt <- demo[ , "treatment"]
demo[!is.na(tt) & tt == 0, "treatment"] <- "NONE"
demo[!is.na(tt) & tt == 1, "treatment"] <- "CT"
demo[!is.na(tt) & tt == 2, "treatment"] <- "HT"

write.csv(t(data), file=file.path(saveres, "data_nki.csv"))
write.csv(annot, file=file.path(saveres, "annot_nki.csv"), row.names=FALSE)
write.csv(demo, file=file.path(saveres, "demo_nki.csv"), row.names=FALSE)

## create individual eSets
ddn <- "NKI"

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
  experimentData(eseti)@preprocessing <- list("normalization"="ORIGINAL", package="agilentRosetta", version="0")
  annotation(eseti) <- "agilentRosetta"
  ## object name
  objn <- sprintf("%sORIGINALPROBE", rownames(sampleinfo)[i])
  assign(objn, eseti)
  save(list=objn, compress=TRUE, file=file.path(toupper(ddn), "esets", sprintf("%s.RData", objn)))
}
