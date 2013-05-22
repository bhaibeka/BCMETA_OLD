###################################################
### Benjamin Haibe-Kains
## 27/05/2009
## All rights reserved
###################################################

rm(list=ls())

## type of survival analysis
meta <- FALSE

##censoring time for urvival data
tc <- 10 * 365

##cutoffs for km survival curves
myprobs <- c(0.33, 0.66)
#myprobs <- c(0.5)

##years for survivap proportions in km survival curves
yrs <- c(3,5,10)

##filtering
myfilt <- list("node"=c("all"), "er"="all", "age"="all", "treatment"="all", "tissue"=1)
#myfilt <- list("node"=0, "er"="all", "age"="all", "treatment"=0, "tissue"=1)

## survival data (dfs or os)
nsurv <- "dfs"

## multiple testing correction
## mycorr <- c("none", "bonferroni", "fdr")
#mtcorr <- c("fdr")
#mtcorr <- c("bonferroni")
mtcorr <- c("none")

## subtype identification model
## sbtclust <- c("ssp2003", "ssp2006", "pam50", "scmgene", "scmod1", "scmod2")
#sbtclust <- "scmod1"
sbtclust <- "scmgene"

## bimodality of ESR1 and ERBB2
## bimodn <- c("gene", "mod1", "mod2")
bimodn <- "gene"

##datasets
## all datasets
#ddn <- c("NKI", "STNO2", "NCI", "KOO", "MSK", "UPP", "STK", "VDX", "UNT", "MAINZ", "UNC2", "DUKE", "DUKE2", "CAL", "TRANSBIG", "EMC2", "NCH", "LUND", "LUND2", "FNCLCC", "MDA", "MDA4", "MDA6", "NCCS", "MGH", "TAM", "TAM2", "VDX3", "IRB", "TOP", "IGR2", "LH", "MDA3", "DFHCC", "DFHCC2", "EORTC10994", "HLP", "MAQC2", "MCCC")
#ddn <- c("VDX", "TRANSBIG", "UNT", "MAINZ")
## only VDX
#ddn <- c("VDX")
ddn <- NULL

#location <- c("laptop", "hydra")
location <- "laptop"

## create directory
saveres <- "saveres"
if(!file.exists(saveres)) { system(sprintf("mkdir %s", saveres)) }

## database
rdp <- "~/MyWork/microarray_data"
ddn.db <- ddn.db <- cbind(
	"dataset"=c("NKI", "STNO2", "NCI", "KOO", "MSK", "UPP", "STK", "VDX", "UNT", "MAINZ", "UNC2", "DUKE", "DUKE2", "CAL", "TRANSBIG", "EMC2", "NCH", "LUND", "LUND2", "FNCLCC", "MDA", "MDA4", "MDA6", "NCCS", "MGH", "TAM", "TAM2", "VDX3", "IRB", "TOP", "IGR2", "LH", "MDA3", "DFHCC", "DFHCC2", "EORTC10994", "HLP", "MAQC2", "MCCC", "MUG", "DFHCC3", "PNC", "EXPO", "LUH", "UCSF", "RH3", "UNC4", "MDA5"),
	"location"=c("nki2002", "sorlie2003", "sotiriou2003", "huang2003", "minn2005", "miller2005", "pawitan2005", "minn2007", "sotiriou2006", "schmidt2008", "hoadley2007", "bild2006", "bonnefoi2007", "chin2006", "transbig2006affy", "bos2009", "naderi2007", "nimeusmalmstrom2008", "saal2007", "campone2007", "hess2006", "liedtke2008", "andre2009", "yu2008", "ma2004", "loi2008", "chanrion2008", "zhang2008", "lu2008", "desmedt2009", "andre2007", "harris2007", "pusztai2008", "li2010", "silver2010", "farmer2005", "natrajan2009", "shi2010", "waddell2009", "calabro2008", "richardson2006", "dedeurwaerder2010", "expO_breast", "gruvberger2001", "korkola2003", "cizkova2010", "prat2010", "symmans2010"),
	"platform"=c("agilent", "cdna.stanford", "cdna.nci", "affy.u95", "affy", "affy", "affy", "affy", "affy", "affy", "agilent2", "affy.u95", "affy.3x", "affy", "affy", "affy", "agilent3", "swegene", "swegene", "umgc.ircna", "affy", "affy", "affy", "affy", "agilent4", "affy", "aminolink", "affy", "affy", "affy", "affy", "affy", "affy", "affy", "affy", "affy", "illumina", "affy", "illumina", "operon","affy", "affy", "affy", "cdna.luh", "cdna.ucsf", "affy", "agilent99", "affy"),
	"reference"=c("citep{veer2002gene,vandevijver2002gene}", "citep{sorlie2003repeated}", "citep{sotiriou2003breast}", "citep{huang2003gene}", "citep{minn2005genes}", "citep{miller2005expression}", "citep{pawitan2005gene}", "citep{wang2005geneexpression,minn2007lung}", "citep{sotiriou2006gene}", "citep{schmidt2008humoral}", "citep{hoadley2007egfr}", "citep{bild2006oncogenic}", "citep{bonnefoi2007validation}", "citep{chin2006genomic}", "citep{desmedt2007strong}", "citep{bos2009genes}", "citep{naderi2007geneexpression}", "citep{nimeusmalmstrom2008gene}", "citep{saal2007poor}", "citep{campone2007prediction}", "citep{hess2006pharmacogenomic}", "citep{liedtke2008pik3caactivating}", "citep{andre2009molecular}", "citep{yu2008precisely}", "citep{ma2004twogene}", "citep{loi2008predicting}", "citep{chanrion2008gene}", "citep{zhang200876gene}", "citep{lu2008predicting}", "citep{desmedt2009multifactor}", "citep{andre2007igr2}", "citep{harris2007predictors}", "citep{pusztai2008mda3}", "citep{li2010amplification}", "citep{silver2010efficacy}", "citep{farmer2005identification}", "citep{natrajan2009integrative}", "citep{shi2010maqc2}", "citep{waddell2009subtypes}", "citep{calabro2008effects}", "citep{richardson2006x}", "citep{dedeurwaerder2010epigenetic}", "citep{bittner2005expo}", "citep{gruvberger2001estrogen}", "citep{korkola2003differentiation}", "citep{cizkova2010gene}", "citep{prat2010phenotypic}", "citep{symmans2010genomic}"),
	"dfs"=c("dmfs", "rfs", "rfs", NA, "dmfs", "rfs", "rfs", "dmfs", "dmfs", "dmfs", "rfs", NA, NA, "dmfs", "dmfs", "dmfs", "dmfs", NA, NA, NA, NA, NA, NA, NA, "dmfs", "dmfs", "dmfs", "dmfs", NA, NA, NA, NA, NA, "dmfs", NA, NA, NA, NA, NA, NA, NA, "rfs", NA, NA, "dmfs", NA, "rfs", "dmfs"), 
	"os"=c("os", "os", NA, NA, NA, NA, NA, NA, NA, NA, "os", "os", NA, "os", "os", NA, "os", NA, NA, NA, NA, NA, NA, NA, NA, NA, "os", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "os", NA, NA, "os", NA, "os", NA), 
	"auto"=c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE)
)
dimnames(ddn.db)[[1]] <- as.character(ddn.db[ ,1])

## subtype identification
if(is.element(sbtclust, c("ssp2003", "ssp2006", "pam50"))) {
	sbtclustn <- sbtclustnn <- c("Basal", "Her2", "LumB", "LumA", "Normal")
	sbtclustn <- factor(x=sbtclustn, levels=sbtclustn)
	sbtclustn.lum <- c("LumA", "LumB")
	sbtclustn.lumn <- "Lum"
	sbtclustn2 <- c("Basal", "Her2", "Lum", "Normal")
	sbtclustn2 <- factor(x=sbtclustn2, levels=sbtclustn2)
	sbtclustn.col <- c("darkred", "darkgreen", "darkorange", "darkviolet", "black")
	names(sbtclustn.col) <- sbtclustn
	sbtclustn.lum.col <- c("darkred", "darkgreen", "darkblue")
	names(sbtclustn.lum.col) <- c("Basal", "Her2", sbtclustn.lumn, "Normal")
}
if(is.element(sbtclust, c("scmgene", "scmod1", "scmod2"))) {
	sbtclustn <- c("ER-/HER2-", "HER2+", "ER+/HER2- High Prolif", "ER+/HER2- Low Prolif")
	sbtclustn <- factor(x=sbtclustn, levels=sbtclustn)
	sbtclustnn <- c("Basal", "Her2", "LumB", "LumA")
	sbtclustn.lum <- c("ER+/HER2- High Prolif", "ER+/HER2- Low Prolif")
	sbtclustn.lumn <- "ER+/HER2-"
	sbtclustn2 <- c("ER-/HER2-", "HER2+", "ER+/HER2-")
	sbtclustn2 <- factor(x=sbtclustn2, levels=sbtclustn2)
	sbtclustn.col <- c("darkred", "darkgreen", "darkorange", "darkviolet")
	names(sbtclustn.col) <- sbtclustn
	sbtclustn.lum.col <- c("darkred", "darkgreen", "darkblue")
	names(sbtclustn.lum.col) <- c("ER-/HER2-", "HER2+", sbtclustn.lumn)
}

## bimodality
switch(bimodn,
"gene"={ bimodn <- c("bimod.ESR1.gene", "bimod.ERBB2.gene") },
"mod1"={ bimodn <- c("bimod.ESR1.mod1", "bimod.ERBB2.mod1") },
"mod2"={ bimodn <- c("bimod.ESR1.mod2", "bimod.ERBB2.mod2") })

###################################################
###  loadlib
###################################################
library(gplots)
library(survival)
library(Hmisc)
library(coin)
library(sma)
library(genefu)
library(survcomp)

###################################################
###  functions
###################################################
source("meta_surv_gene_foo.R")

###################################################
###  alldata
###################################################
##demo for all patients
demo.all <- read.csv("meta_demo.csv", stringsAsFactors=FALSE)
dimnames(demo.all)[[1]] <- as.character(demo.all[ ,1])

## remove duplicates
idtoremove <- read.csv("find_duplicates_microarray_data.csv", stringsAsFactors=FALSE)[ , 1]
## actually remove duplicated ids
demo.all <- demo.all[!is.element(dimnames(demo.all)[[1]], idtoremove), , drop=FALSE]

##standard name of columns except genes/signatures
snc <- c("samplename", "dataset", "series", "id", "age", "size", "node", "er", "pgr", "her2", "grade", "t.rfs", "e.rfs", "t.dmfs", "e.dmfs", "t.os", "e.os", "treatment", "tissue", c(paste("bimod.ESR1", c("gene", "mod1", "mod2"), sep="."), paste("bimod.ERBB2", c("gene", "mod1", "mod2"), sep=".")), c("ssp2003", "ssp2006", "pam50", "scmgene", "scmod1", "scmod2"))

## clinical info
cinfo <- c("samplename", "dataset", "series", "id", "age", "size", "node", "er", "pgr", "her2", "grade", paste("t", c("rfs", "dmfs", "os"), sep="."), paste("e", c("rfs", "dmfs", "os"), sep="."), "treatment", "tissue")

## clinical variables to consider in the survival analysis
#nc <- c("age", "size", "node", "er", "grade", "treatment") #names of clinical variables
nc <- c("age", "size", "node", "er", "grade") #names of clinical variables

## levels for each clinical variable
ncl <- list("age"=c("0", "1"), "size"=c("0", "1"), "node"=c("0", "1"), "er"=c("0", "1"), "grade"=c("0", "0.5", "1"), "treatment"=c("0","1"))

## genes or signatures to consider in the survival analysis
nm <- dimnames(demo.all)[[2]][!is.element(dimnames(demo.all)[[2]], snc)]

## clinical info
clinic.all <- demo.all[ ,cinfo,drop=FALSE]
dimnames(clinic.all)[[1]] <- dimnames(demo.all)[[1]]
## modify some clinical variables
clinic.all[ ,"age"] <- ifelse(clinic.all[ ,"age"] > 50, 1, 0)
clinic.all[ ,"size"] <- ifelse(clinic.all[ ,"size"] > 2, 1, 0)
clinic.all[ ,"grade"] <- (clinic.all[ ,"grade"] - 1) / 2
if(is.element("treatment", nc)) { clinic.all[ ,"treatment"] <- ifelse(clinic.all[ ,"treatment"] > 0, 1, 0) }
## censoring survival data
stt <- set <- rep(NA, nrow(clinic.all))
names(stt) <- names(set) <- dimnames(clinic.all)[[1]]
udn <- sort(unique(clinic.all[ , "dataset"]))
if(is.null(ddn)) { ddn <- udn }
for(i in 1:length(udn)) {
	myx <- clinic.all[ , "dataset"] == udn[i]
	if(!is.na(ddn.db[udn[i], nsurv])) {
		stt[myx] <- clinic.all[myx, paste("t", ddn.db[udn[i], nsurv], sep=".")]
		set[myx] <- clinic.all[myx, paste("e", ddn.db[udn[i], nsurv], sep=".")]
	} else { stt[myx] <- set[myx] <- NA }
}
tt <- censor.time(surv.time=stt, surv.event=set, time.cens=tc)
clinic.all <- cbind(clinic.all, "surv.time"=tt[[1]], "surv.event"=tt[[2]])

##gene/signature score
mscore.all <- demo.all[ ,nm,drop=FALSE]
dimnames(mscore.all)[[1]] <- dimnames(demo.all)[[1]]

##subtype
subtype.all <- factor(x=demo.all[ , sbtclust], levels=sbtclustn)
names(subtype.all) <- dimnames(demo.all)[[1]]

##bimodality of ESR1 and ERBB2
bimod.all <- demo.all[ ,bimodn,drop=FALSE]
dimnames(bimod.all)[[1]] <- dimnames(demo.all)[[1]]

###################################################
###  datams
###################################################
ddn.db2 <- ddn.db[ddn, ,drop=FALSE]
data.all <- dn <- NULL
for(i in 1:nrow(ddn.db2)) {
	rest <- sbt.probat <- NULL
	cat(sprintf("%s\t", ddn.db2[i,"dataset"]))
	dn <- c(dn, ddn.db2[i,"dataset"])
	
	## load data
	#load(sprintf("%s/%s/%s.RData", rdp, ddn.db2[i,"location"], ddn.db2[i,"location"]))
	
	#myx <- clinic.all[ ,"dataset"] == ddn.db2[i,"dataset"] & clinic.all[ ,"tissue"] == 1 & clinic.all[ ,"treatment"] == 0
	myx <- clinic.all[ ,"dataset"] == ddn.db2[i,"dataset"]
	myx[is.na(myx)] <- FALSE
	clinic <- clinic.all[myx, ,drop=FALSE]
	mscore <- mscore.all[myx, ,drop=FALSE]
	subtype <- subtype.all[myx]
	bimod <- bimod.all[myx, ,drop=FALSE]
	data.all <- c(data.all, list(list("clinic"=clinic, "mscore"=mscore, "subtype"=subtype, "bimod"=bimod)))
	names(data.all)[length(data.all)] <- ddn.db2[i,"dataset"]
}

###################################################
### saveres
###################################################
myf <- function(x) {
	names(x$clinic)[is.element(names(x$clinic), "surv.time")] <- "time"
	names(x$clinic)[is.element(names(x$clinic), "surv.event")] <- "event"
	return(x)
}
nt <- sum(unlist(lapply(data.all, function(x) { return(length(x$subtype)) })))

data.rmdupl <- NULL
for(i in 1:length(data.all)) {
data.rmdupl <- rbind(data.rmdupl, cbind(data.all[[i]]$clinic, data.all[[i]]$mscore, "subtype"=as.character(data.all[[i]]$subtype)))
}
write.csv(data.rmdupl, row.names=FALSE, file=sprintf("%s/meta_demo_rmdupl.csv", saveres))
data.all <- lapply(data.all, myf)

###################################################
### datafiltering
###################################################
data.all2 <- data.all
torm <- rep(FALSE, length(data.all))
for(j in 1:length(data.all)) {
	myx <- rep(TRUE, nrow(data.all[[j]]$clinic))
	for(i in 1:length(myfilt)) {	
		if(myfilt[[i]][1] != "all") { myx <- myx & is.element(data.all[[j]]$clinic[ ,names(myfilt)[i]], myfilt[[i]]) }
	}
	if(all(!myx)) { torm[j] <- TRUE }
	data.all2[[j]]$clinic <- data.all2[[j]]$clinic[myx, , drop=FALSE]
	data.all2[[j]]$mscore <- data.all2[[j]]$mscore[myx, , drop=FALSE]
	data.all2[[j]]$gcscore <- data.all2[[j]]$gcscore[myx, , drop=FALSE]
	data.all2[[j]]$subtype <- data.all2[[j]]$subtype[myx]
	data.all2[[j]]$bimod <- data.all2[[j]]$bimod[myx, , drop=FALSE]	
}
data.all2 <- data.all2[!torm]

data.all <- data.all2

myf.demod <- function(x, subtype) {
nc.info <- c(nc, "event")
ncl.info <- c(ncl, list("event"=c("0", "1")))

demod <- matrix(NA, nrow=length(x), ncol=length(nc.info))
dimnames(demod) <- list(names(x), nc.info)
for(i in 1:length(x)) {
	if(subtype[1] != "all") { myx <- is.element(x[[i]]$subtype, subtype) & !is.na(x[[i]]$subtype) }
	else { myx <- rep(TRUE, length(x[[i]]$subtype)) }
	dd <- x[[i]]$clinic[myx, ,drop=FALSE]
	for(j in 1:length(nc.info)) {
		nc.infoc <- rep(0, length(ncl.info[[j]]))
		names(nc.infoc) <- ncl.info[[j]]
		temp <- table(dd[ ,nc.info[j]])
		nc.infoc[names(temp)] <- temp
		demod[i,j] <- paste(nc.infoc, collapse="/")
	}
}
myf.demod2 <- function(x) {
	y <- strsplit(x, "/")
	temp <- NULL
	for(i in 1:length(y[[1]])) {
		temp <- c(temp, sum(as.numeric(unlist(lapply(y, function(x, i) { return(x[i]) }, i=i)))))
	}
	return(temp)
}
demod <- rbind(demod, "total"=unlist(lapply(apply(demod, 2, myf.demod2), paste, collapse="/")))
return(demod)
}

###################################################
### demod
###################################################
nnn <- c(as.list(c("all", as.character(sbtclustn))), list(sbtclustn.lum))
names(nnn) <- c("all", as.character(sbtclustnn), "Lum")
for(h in 1:length(nnn)) {
	demod <- myf.demod(x=data.all, subtype=nnn[[h]])
	write.csv(demod, file=sprintf("%s/demo_untreated_%s.csv", saveres, names(nnn)[h]))
}

###################################################
### kmsubtype
###################################################
sbt <- strat <- st <- se <- NULL
for(i in 1:length(data.all)) {
	st <- c(st, data.all[[i]]$clinic[ ,"time"])
	se <- c(se, data.all[[i]]$clinic[ ,"event"])
	strat <- c(strat, as.character(data.all[[i]]$clinic[ ,"dataset"]))
	sbt <- c(sbt, as.character(data.all[[i]]$subtype))
}
sbt2 <- sbt
sbt2[is.element(sbt, sbtclustn.lum)] <- sbtclustn.lumn
usbt2 <- unique(sbt2)
usbt2 <- usbt2[!is.na(usbt2)]
usbt2 <- usbt2[match(names(sbtclustn.lum.col), usbt2)]

dd <- data.frame("time"=st / 365, "event"=se, "strata"=strat, "subtype"=factor(x=sbt2, levels=sbtclustn2))
dd <- dd[complete.cases(dd), ]

pdf(sprintf("%s/allkmsubtypelum.pdf", saveres))
km.coxph.plot(formula.s=formula(Surv(time, event) ~ subtype), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of survival", main.title="", sub.title=NULL, leg.text=paste(names(sbtclustn.lum.col), "   "), leg.pos="bottomright", leg.inset=0.0, v.line=NULL, h.line=NULL, .col=sbtclustn.lum.col, .lty=c(1,1,1), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, bty="n", leg.bty="n", verbose=FALSE)
dev.off()

dd <- data.frame("time"=st / 365, "event"=se, "strata"=strat, "subtype"=factor(x=sbt, levels=sbtclustn))
pdf(sprintf("%s/allkmsubtype.pdf", saveres))
km.coxph.plot(formula.s=formula(Surv(time, event) ~ subtype), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of survival", main.title="", sub.title=NULL, leg.text=paste(names(sbtclustn.col), "   "), leg.pos="bottomright", leg.inset=0.0, v.line=NULL, h.line=NULL, .col=sbtclustn.col, .lty=c(1,1,1), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, bty="n", leg.bty="n", verbose=FALSE)
dev.off()

###################################################
### mscorepairwisecor
###################################################
if(ncol(data.all[[1]]$mscore) > 1) {
nnn <- c(as.list(c("all", as.character(sbtclustn))), list(sbtclustn.lum))
names(nnn) <- c("all", as.character(sbtclustnn), "Lum")
for(h in 1:length(nnn)) {
	dd <- myf.subtype(x=data.all, subtype=nnn[[h]])
	
	if(meta) { ## meta-estimation
		mycorz <- array(NA, c(length(nm), length(nm), length(dd)))
		dimnames(mycorz) <- list(nm, nm, names(dd))
		mycorz.se <- array(NA, c(length(nm), length(nm), length(dd)))
		dimnames(mycorz.se) <- list(nm, nm, names(dd))
		for(k in 1:length(dd)) {
			if(nrow(dd[[k]]$mscore) >= 3) {
				rr <- cor(dd[[k]]$mscore, use="pairwise.complete.obs", method="pearson")
				for(i in 2: nrow(rr)) {
					for(j in 1:(i-1)) {
						mycorz[i,j,k] <- fisherz(rr[i,j], inv=FALSE)
						nn <- sum(complete.cases(dd[[k]]$mscore[ ,c(i,j)]))
						if(nn >3) { mycorz.se[i,j,k] <- 1 / sqrt(nn - 3) }
						else { mycorz.se[i,j,k] <- NA }
					}
				}
			}
		}

		mycor.meta <- matrix(NA, nrow=length(nm), ncol=length(nm))
		dimnames(mycor.meta) <- list(nm, nm)
		for(i in 2: length(nm)) {
			for(j in 1:(i-1)) {
				xx <- NULL
				xx.se <- NULL
				for(k in 1:length(dd)) {
					xx <- c(xx, mycorz[i,j,k])
					xx.se <- c(xx.se, mycorz.se[i,j,k])
				}
				#print(cbind(xx, xx.se))
				mycor.meta[i,j] <- mycor.meta[j,i] <- fisherz(combine.est(x=xx, x.se=xx.se, na.rm=TRUE,hetero=FALSE)$estimate, inv=TRUE)
			}
		}
	} else { ## pool data
		ddt <- NULL
		for(k in 1:length(dd)) {
			ddt <- rbind(ddt, dd[[k]]$mscore)
		}
	 	mycor.meta <- cor(ddt, use="pairwise.complete.obs", method="pearson")
	}
	diag(mycor.meta) <- 1

	## pairwise correlations heatmap
	dd <- mycor.meta
	nix <- apply(dd, 1, function(x) { return(sum(is.na(x))) })
	dd <- dd[nix < (nrow(dd) - 1), nix < (nrow(dd) - 1)]

	pdf(sprintf("%s/hkb-mpairwcorheatmap_%s.pdf", saveres, names(nnn)[h]), width=14, height=14,onefile=FALSE)
	plot.cor(dd, labels=dimnames(dd)[[1]], zlim=c(-1,1))
	dev.off()

	dd <- mycor.meta
	write.csv(dd, sprintf("%s/pairw_cor_%s.csv", saveres, names(nnn)[h]))
}
}

###################################################
### survival meta analysis
###################################################
nnn <- c(as.list(c("all", as.character(sbtclustn))), list(sbtclustn.lum))
names(nnn) <- c("all", as.character(sbtclustnn), "Lum")
for(h in 1:length(nnn)) {
	#format data
	data.all.2 <- lapply(data.all, myf.subtype2, subtype=nnn[[h]])
	names(data.all.2) <- names(data.all)
	#merge datasets
	data.all.3 <- NULL
	for(i in 1:length(data.all.2)) {
		data.all.3$input <- rbind(data.all.3$input, data.all.2[[i]]$input)
		data.all.3$strata <- c(data.all.3$strata, as.character(data.all.2[[i]]$strata))
		data.all.3$output <- rbind(data.all.3$output, data.all.2[[i]]$output)
		data.all.3$subtype <- c(data.all.3$subtype, data.all.2[[i]]$subtype)
		data.all.3$bimod <- rbind(data.all.3$bimod, data.all.2[[i]]$bimod)
	}
	#remove non complete cases
	data.all.4 <- NULL
	cc.ix <- complete.cases(data.all.3$input, data.all.3$strata, data.all.3$output, data.all.3$subtype)
	data.all.4$input <- data.all.3$input[cc.ix, ,drop=FALSE]
	data.all.4$strata <- as.character(data.all.3$strata[cc.ix])
	data.all.4$output <- data.all.3$output[cc.ix, ,drop=FALSE]
	data.all.4$subtype <- as.character(data.all.3$subtype[cc.ix])
	data.all.4$bimod <- data.all.3$bimod[cc.ix, ,drop=FALSE]

	mymetasurv <- NULL
	if(meta) { ## meta-estimation
		rr <- surv.feature.rel(data=data.all.2, prior=NULL, method="hr")
		tt <- NULL
		for(i in 1:ncol(data.all.2[[1]]$input)) { tt <- c(tt,  list(t(sapply(rr[[1]], function(x,y) { return(x$estimate[y,]) }, y=i)))) }
		names(tt) <- dimnames(data.all.2[[1]]$input)[[2]]
		tt2 <- lapply(tt, function(x) { res <- t(apply(x, 1, function(x) { return(c("hr"=exp(x["x"]), "lower.95"=exp(x["x"] - qnorm(p=0.975) * x["x.se"]), "upper.95"=exp(x["x"] + qnorm(p=0.975) * x["x.se"]), "p"=pchisq((x["x"] / x["x.se"])^2, df=1, lower.tail=FALSE), "n"=x["n"])) })); res <- cbind("dataset"=dimnames(x)[[1]], res); return(res) })
		write.m.file(tt2, file=sprintf("%s/univ_meta_surv2_%s.csv", saveres, names(nnn)[h]))
		rr <- rr$meta
		
	} else { ## pool data
		rr <- surv.feature.rel(data=list("pool"=data.all.3), prior=NULL, method="hr")$meta
	}
	univ.meta.res <- NULL
	for(i in c(nc, nm)) {
		mme2 <- rr$model[ , ,is.element(dimnames(rr$model)[[3]], i)]
		univ.meta.res <- rbind(univ.meta.res, cbind("hr"=exp(mme2["coef"]), "lower.95"=exp(mme2["coef"] - qnorm(p=0.975) * mme2["coef.se"]), "upper.95"=exp(mme2["coef"] + qnorm(p=0.975) * mme2["coef.se"]), "p"=pchisq((mme2["coef"] / mme2["coef.se"])^2, df=1, lower.tail=FALSE), "n"=mme2["n"]))
	}
	dimnames(univ.meta.res)[[1]] <- c(nc, nm)
	if(mtcorr != "none") { ## multiple testing correction
		univ.meta.res <- cbind(univ.meta.res, p.adjust(univ.meta.res[ , "p"], method=mtcorr))
		dimnames(univ.meta.res)[[2]][ncol(univ.meta.res)] <- mtcorr
	}
	write.csv(univ.meta.res, file=sprintf("%s/univ_meta_surv_%s.csv", saveres, names(nnn)[h]))

	## univforestplot
	nnc <- length(nc)
	nnm <- length(nm)
	myspace <- "        "
	rn <- c(dimnames(univ.meta.res)[[1]][1:nnc], NA,  dimnames(univ.meta.res)[[1]][(nnc + 1):(nnc + nnm)])
	rn <- cbind(rn, rep(myspace, length(rn)))
	r.mean <- c(univ.meta.res[1:nnc,"hr"], NA, univ.meta.res[(nnc + 1):(nnc + nnm),"hr"])
	r.lower <- c(univ.meta.res[1:nnc,"lower.95"], NA, univ.meta.res[(nnc + 1):(nnc + nnm),"lower.95"])
	r.upper <- c(univ.meta.res[1:nnc,"upper.95"], NA, univ.meta.res[(nnc + 1):(nnc + nnm),"upper.95"])
	bs <- c(univ.meta.res[1:nnc,"n"], NA, univ.meta.res[(nnc + 1):(nnc + nnm),"n"]) /  max(univ.meta.res[ ,"n"], na.rm=TRUE)
	#bs <- c(rep(0.3, nnc), NA, rep(0.3, nnm))
	#names(bs) <- rn[ ,1]

	pdf(sprintf("%s/univforestplot_%s.pdf", saveres, names(nnn)[h]), height=20, width=4)
	myforestplot(labeltext=as.matrix(rn), mean=log2(r.mean), lower=log2(r.lower), upper=log2(r.upper), zero=0, align=c("l"), graphwidth=unit(1, "inches"), x.ticks=seq(-3,3,1), xlab=paste("log2 hazard ratio", myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue"), clip=c(-3,3), box.size=bs/2)
	dev.off()

	dimnames(univ.meta.res)[[2]] <- paste("ALL", dimnames(univ.meta.res)[[2]], sep=".")
	mymetasurv <- cbind(mymetasurv, univ.meta.res)

	## discretize data
	if(meta) { ## meta-estimation
		data.all.2.disc <- lapply(data.all.2, myf.subv, y=myprobs, subv=nm)
		#merge datasets
		data.all.3.disc <- NULL
		for(i in 1:length(data.all.2.disc)) {
			data.all.3.disc$input <- rbind(data.all.3.disc$input, data.all.2.disc[[i]]$input)
			data.all.3.disc$strata <- c(data.all.3.disc$strata, as.character(data.all.2.disc[[i]]$strata))
			data.all.3.disc$output <- rbind(data.all.3.disc$output, data.all.2.disc[[i]]$output)
		}
	} else { ## pool data
		data.all.3.disc <- data.all.3
		temp <- data.all.3$input
		mm <- apply(temp, 2, quantile, probs=myprobs, na.rm=TRUE)
		if(is.null(ncol(mm))) { mm <- rbind(mm) }
		for(i in 1:ncol(temp)) {
			if(!all(is.na(temp[ ,i]))) {
				ccc <- unique(mm[ , i])
				ccc <- ccc[!is.na(ccc)]
				if(length(ccc) == length(myprobs)) {
					temp[ ,i] <- as.numeric(cut2(x=temp[ ,i], cuts=mm[ ,i], levels.mean=TRUE))
				}
			}
		}
		data.all.3.disc$input <- temp
	}
	pdfn <- sprintf("km_%s_", names(nnn)[[h]])
	lll <- NULL
	for(i in 1:length(nm)) {
		ggg <- data.all.3.disc$input[ ,nm[i]]
		if(!all(is.na(ggg))) { ggg <- (ggg - min(ggg, na.rm=TRUE)) / (max(ggg, na.rm=TRUE) - min(ggg, na.rm=TRUE)) }
		dd <- data.frame("time"=data.all.3.disc$output[ ,"surv.time"] / 365, "event"=data.all.3.disc$output[ ,"surv.event"], "strata"=data.all.3.disc$strata, "group"=ggg)
		if(!all(!complete.cases(dd))) {
			lograp <- 1 - pchisq(survdiff(Surv(time, event) ~ group, data=dd, rho=0)$chisq, length(table(dd$group))-1)
			rr <- coxph(Surv(time, event) ~ strata(strata) + group, data=dd)
			#ppp <- pchisq(q=(rr$coefficients / sqrt(drop(rr$var)))^2, df=1, lower.tail=FALSE)
			ppp <- summary(rr)$sctest["pvalue"]
			rr <- cbind("hr"=exp(rr$coefficients), "lower.95"=exp((rr$coefficients - qnorm(p=0.975) * sqrt(drop(rr$var)))), "upper.95"=exp((rr$coefficients + qnorm(p=0.975) * sqrt(drop(rr$var)))), "p"=ppp, "n"=rr$n)
			dimnames(rr)[[1]] <- nm[i]
			myfile <- sprintf("%s/%s%s.pdf", saveres, pdfn, nm[i])
			tt <- list(myfile, dd, drop(rr), lograp)
			lll <- c(lll, list(tt))
		}
	}
	if(mtcorr != "none") {
		pcorr <- p.adjust(unlist(lapply(lll, function(x) { return(x[[3]]["p"]) })), method=mtcorr)
		pcorr2 <- p.adjust(unlist(lapply(lll, function(x) { return(x[[4]]) })), method=mtcorr)
		lll <- mapply(function(x, y, z) {
			x[[3]]["p"] <- y
			x[[4]] <- z
			return(x)
		}, x=lll, y=pcorr, z=pcorr2, SIMPLIFY=FALSE)
	}
	## actually plotting the survival curves
	lapply(lll, function(x) {
		pdf(x[[1]], height=7, width=7)
		km.coxph.plot(formula.s=formula(Surv(time, event) ~ group), data.s=x[[2]], sub.s="all", x.label="Time (years)", y.label="Probability of survival", main.title="", sub.title=NULL, leg.text=paste(c("Low", "Intermediate", "High"), "     "), leg.pos="bottomright", leg.inset=0.0, v.line=NULL, h.line=NULL, .col=c("darkblue", "darkgreen", "darkred"), .lty=c(1,1,1), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, o.text=ifelse(x[[4]] < 0.10, sprintf("Logrank p=%.1E\nHR=%.2g, 95%%CI[%.2g,%.2g], p=%.1E", x[[4]], x[[3]]["hr"], x[[3]]["lower.95"], x[[3]]["upper.95"], x[[3]]["p"]), sprintf("Logrank p=%.1E\n", x[[4]])), bty="n", leg.bty="n", verbose=FALSE)
		dev.off()
	})
	## save the discretized data
	data.rmdupl.temp <- data.frame(data.rmdupl[rownames(data.all.3.disc$input),!is.element(colnames(data.rmdupl), "subtype")], data.all.3.disc$input[ , nm, drop=FALSE], "subtypes"=data.rmdupl[rownames(data.all.3.disc$input), "subtype"])
	write.csv(data.rmdupl.temp, row.names=FALSE, file=sprintf("%s/meta_demo_discr_pool_rmdupl_%s.csv", saveres, names(nnn)[[h]]))
}
cat("\n")
