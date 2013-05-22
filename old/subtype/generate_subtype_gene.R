
###################################################
### code chunk number 2: loadlib
###################################################
library(xtable)
library(gplots)

library(genefu)
## version should be higher than 1.5.0
vv <- as.numeric(unlist(strsplit(packageDescription("genefu")$Version, "[.]")))
vv <- sum(rev(vv) * 10^(0:(length(vv)-1)))
if(vv < 150) { stop("version of genefu should be >= 1.5.0")}


###################################################
### code chunk number 3: selectedsets
###################################################
ddn <- c("NKI", "UNC4", "EXPO", "STNO2", "NCI", "KOO", "MSK", "UPP", "STK", "VDX", "UNT", "MAINZ", "DUKE", "DUKE2", "CAL", "TRANSBIG", "EMC2", "NCH", "LUND", "LUND2", "FNCLCC", "MDA", "MDA4", "MDA6", "NCCS", "MGH", "TAM", "TAM2", "VDX3", "IRB", "TOP", "IGR2", "MDA3", "DFHCC", "DFHCC2", "EORTC10994", "HLP", "MAQC2", "MCCC", "MUG", "DFHCC3", "PNC", "LUH", "UCSF") ## all datasets
#ddn <- c("NKI", "UNC4", "TOP", "MGH", "TAM", "TAM2", "MDA3", "VDX3", "UCSF")
#ddn <- c("VDX", "TRANSBIG", "UNT", "MAINZ")
#ddn <- c("VDX", "TRANSBIG", "MAINZ", "UNT", "UPP", "EMC2", "NKI")


###################################################
### code chunk number 4: metagenelist
###################################################
glist2 <- read.m.file(file="meta_gene_list.csv")


###################################################
### code chunk number 5: roptions
###################################################
#location <- c("laptop", "hydra")
location <- "laptop"

## create directory
saveres <- "saveres"
if(!file.exists(saveres)) { system(sprintf("mkdir %s", saveres)) }
#fig <- "fig"
#if(!file.exists(fig)) { system(sprintf("mkdir %s", fig)) }
boxres <- "boxplot"
if(!file.exists(sprintf("%s/%s", saveres, boxres))) { system(sprintf("mkdir %s/%s", saveres, boxres)) }
sbtplotres <- "subtypes"
if(!file.exists(sprintf("%s/%s", saveres, sbtplotres))) { system(sprintf("mkdir %s/%s", saveres, sbtplotres)) }

## avoid mapping if you are sure that the signatures use the same probe ids that in the datasets
domap.sigs <- TRUE

## quantile for the rescaling
resq <- 0.05

## datasets
rdp <- "../../microarray_data"
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

## colours and symbols
sbtnn <- c("ER-/HER2-", "HER2+", "ER+/HER2-") 
sbtn <- c("ER-/HER2-", "HER2+", "ER+/HER2- High Prolif", "ER+/HER2- Low Prolif")
sbtn.intrinsic <- c("Basal", "Her2", "LumB", "LumA", "Normal")
colg <- c("darkred", "darkgreen", "darkorange", "darkviolet")
colg.intrinsic <- c("darkred", "darkgreen", "darkorange", "darkviolet", "black")
pchg <- c(17, 0, 10, 10)
pchg.intrinsic <- c(17, 0, 10, 10, 5)
names(pchg) <- names(colg) <- sbtn
names(pchg.intrinsic) <- names(colg.intrinsic) <- sbtn.intrinsic


###################################################
### code chunk number 6: clinicinfo
###################################################
cinfo <- c("samplename", "dataset", "series", "id", "age", "size", "node", "er", "pgr", "her2", "grade", "t.rfs", "e.rfs", "t.dmfs", "e.dmfs", "t.os", "e.os", "treatment", "tissue")


###################################################
### code chunk number 7: autodatasets
###################################################
ddn.db2 <- ddn.db[ddn, ,drop=FALSE]
map.list.all <- t(unlist(lapply(glist2, nrow)))
dimnames(map.list.all) <- list("total", names(glist2))
res <- dn <- demo.all <- sbt.proba <- NULL
for(i in 1:nrow(ddn.db2)) {
	rest <- sbt.probat <- NULL
	if(ddn.db2[i,"auto"]) {
		cat(sprintf("%s\t", ddn.db2[i,"dataset"]))
		dn <- c(dn, ddn.db2[i,"dataset"])
		## the dataset does not require specific code
		
		## load data
		load(sprintf("%s/%s/%s.RData", rdp, ddn.db2[i,"location"], ddn.db2[i,"location"]))
		annot[ ,"EntrezGene.ID"] <- as.character(annot[ ,"EntrezGene.ID"])
		
		## collect clinical information
		demo.all <- rbind(demo.all, apply(X=demo[ ,cinfo,drop=FALSE], MARGIN=2, FUN=as.character))
		
		## bimodality of ESR1 and ERBB2
		cat(sprintf(" -> bimodality"))
		## bimodality for gene alone
		rest <- cbind(rest, "bimod.ESR1.gene"=bimod(x=mod1$ESR1[1, ,drop=FALSE], data=data, annot=annot, do.mapping=ifelse(ddn.db2[i,"platform"] == "affy", FALSE, TRUE), model="V")$status)
		rest <- cbind(rest, "bimod.ERBB2.gene"=bimod(x=mod1$ERBB2[1, ,drop=FALSE], data=data, annot=annot, do.mapping=ifelse(ddn.db2[i,"platform"] == "affy", FALSE, TRUE), model="V")$status)
		## bimodality for ESR1 and ERBB2 mod1
		rest <- cbind(rest, "bimod.ESR1.mod1"=bimod(x=mod1$ESR1, data=data, annot=annot, do.mapping=ifelse(ddn.db2[i,"platform"] == "affy", FALSE, TRUE), model="V")$status)
		rest <- cbind(rest, "bimod.ERBB2.mod1"=bimod(x=mod1$ERBB2, data=data, annot=annot, do.mapping=ifelse(ddn.db2[i,"platform"] == "affy", FALSE, TRUE), model="V")$status)
		## bimodality for ESR1 and ERBB2 mod2
		rest <- cbind(rest, "bimod.ESR1.mod2"=bimod(x=mod1$ESR1, data=data, annot=annot, do.mapping=TRUE, model="V")$status)
		rest <- cbind(rest, "bimod.ERBB2.mod2"=bimod(x=mod2$ERBB2, data=data, annot=annot, do.mapping=TRUE, model="V")$status)	
		
		## subtypes identification	
		cat(sprintf(" -> subtypes"))
		pdf(sprintf("%s/%s/subtypes_id_plot_%s.pdf", saveres, sbtplotres, ddn.db2[i,"dataset"]), height=7, width=7)
		## SSP2003
		tt <- intrinsic.cluster.predict(sbt.model=ssp2003.robust, data=data, annot=annot, do.mapping=TRUE, do.prediction.strength=FALSE)
		rest <- cbind(rest, "ssp2003"=tt$subtype)
		dimnames(tt$subtype.proba)[[2]] <- paste("ssp2003", dimnames(tt$subtype.proba)[[2]], sep=".")
		sbt.probat <- cbind(sbt.probat, "ssp2003"=tt$subtype.proba)
		## SSP2006
		tt <- intrinsic.cluster.predict(sbt.model=ssp2006.robust, data=data, annot=annot, do.mapping=TRUE, do.prediction.strength=FALSE)
		rest <- cbind(rest, "ssp2006"=tt$subtype)
		dimnames(tt$subtype.proba)[[2]] <- paste("ssp2006", dimnames(tt$subtype.proba)[[2]], sep=".")
		sbt.probat <- cbind(sbt.probat, "ssp2006"=tt$subtype.proba)
		## PAM50
		tt <- intrinsic.cluster.predict(sbt.model=pam50.robust, data=data, annot=annot, do.mapping=TRUE, do.prediction.strength=FALSE)
		rest <- cbind(rest, "pam50"=tt$subtype)
		dimnames(tt$subtype.proba)[[2]] <- paste("pam50", dimnames(tt$subtype.proba)[[2]], sep=".")
		sbt.probat <- cbind(sbt.probat, "pam50"=tt$subtype.proba)
		## SCMGENE
		tt <- subtype.cluster.predict(sbt.model=scmgene.robust, data=data, annot=annot, do.mapping=ifelse(ddn.db2[i,"platform"] == "affy", FALSE, TRUE), do.prediction.strength=FALSE, do.BIC=FALSE, plot=TRUE)
		title(main=sprintf("SCMGENE\n%s", ddn.db2[i,"dataset"]))
		rest <- cbind(rest, "scmgene"=tt$subtype2)
		dimnames(tt$subtype.proba2)[[2]] <- paste("scmgene", dimnames(tt$subtype.proba2)[[2]], sep=".")
		sbt.probat <- cbind(sbt.probat, "scmgene"=tt$subtype.proba2)
		## SCMOD1
		tt <- subtype.cluster.predict(sbt.model=scmod1.robust, data=data, annot=annot, do.mapping=ifelse(ddn.db2[i,"platform"] == "affy", FALSE, TRUE), do.prediction.strength=FALSE, do.BIC=FALSE, plot=TRUE)
		title(main=sprintf("SCMOD1\n%s", ddn.db2[i,"dataset"]))
		rest <- cbind(rest, "scmod1"=tt$subtype2)
		dimnames(tt$subtype.proba2)[[2]] <- paste("scmod1", dimnames(tt$subtype.proba2)[[2]], sep=".")
		sbt.probat <- cbind(sbt.probat, "scmod1"=tt$subtype.proba2)
		## SCMOD2
		tt <- subtype.cluster.predict(sbt.model=scmod2.robust, data=data, annot=annot, do.mapping=TRUE, do.prediction.strength=FALSE, do.BIC=FALSE, plot=TRUE)
		title(main=sprintf("SCMOD2\n%s", ddn.db2[i,"dataset"]))
		rest <- cbind(rest, "scmod2"=tt$subtype2)
		dimnames(tt$subtype.proba2)[[2]] <- paste("scmod2", dimnames(tt$subtype.proba2)[[2]], sep=".")
		sbt.probat <- cbind(sbt.probat, "scmod2"=tt$subtype.proba2)
		dev.off()
		
		## signature scores computations
		cat(sprintf(" -> signatures"))
		mysig <- mymap <- NULL
		for(j in 1:length(glist2)) {
			tt <- sig.score(x=glist2[[j]], data=data, annot=annot, do.mapping=domap.sigs, signed=TRUE)
			mymap <- c(mymap, tt$mappin[1])
			mysig <- cbind(mysig, (rescale(tt$score, q=resq, na.rm=TRUE) - 0.5) * 2)
		}
		dimnames(mysig)[[2]] <- names(glist2)
		rest <- cbind(rest, mysig)
		
		map.list.all <- rbind(mymap, map.list.all)
		res <- rbind(res, rest)
		sbt.proba <- rbind(sbt.proba, sbt.probat)
				
		## boxplots without points
		cat(sprintf(" -> boxplots"))
		sbt.models <- c("ssp2003", "ssp2006", "pam50", "scmgene", "scmod1", "scmod2")
		for(j in 1:ncol(mysig)) {
			## for each signature
			if(!all(!complete.cases(mysig[ ,j]))) {
				pdf(sprintf("%s/%s/boxplot_%s_%s.pdf", saveres, boxres, dimnames(mysig)[[2]][j], dn[length(dn)]), width=9, height=9)
				for(k in 1:length(sbt.models)) {
					## for each subtypes identifications
					dd2 <- NULL
					## create a list with the signature scores for each subtype
					oo <- sort(unique(as.character(rest[ ,sbt.models[k]])))
					for(m in 1:length(oo)) { dd2 <- c(dd2, list(mysig[rest[ ,sbt.models[k]] == oo[m], j])) }
					names(dd2) <- oo
					mycol <- mypch <- NULL
					if(is.element(sbt.models[k], c("ssp2003", "ssp2006", "pam50"))) {
						oo <- oo[match(sbtn.intrinsic, oo)]
						mycolg <- colg.intrinsic[oo]
						mypchg <- pchg.intrinsic[oo]
						dd2 <- dd2[sbtn.intrinsic]
					} else {
						oo <- oo[match(sbtn, oo)]
						mycolg <- colg[oo]
						mypchg <- pchg[oo]
						dd2 <- dd2[sbtn]
					}
					for(m in 1:length(oo)) {
					mycol <- c(mycol, rep(mycolg[oo[m]], length(dd2[[m]])))
					mypch <- c(mypch,  rep(mypchg[oo[m]], length(dd2[[m]])))
				}
				kwt <- NULL
				uu <- unique(as.character(rest[ ,sbt.models[k]]))
				uu <- uu[!is.na(uu)]
				if(length(uu) > 1) { kwt <- kruskal.test(mysig[ ,j] ~ as.factor(rest[ ,sbt.models[k]])) } else { kwt$p.value <- NA }
				dd3 <- c(list(" "=NA), list("  "=NA), dd2)
				par(mar=c(10.1,4.1,4.1,2.1))
				rr <- graphics::boxplot(x=dd3, las=3, outline=F, ylim=c(-2,2), main=sprintf("%s in %s\n%s", dimnames(mysig)[[2]][j], dn[length(dn)], sbt.models[k]), col="lightgrey")
				dd3 <- dd3[!is.na(dd3)]
				smartlegend(x="left", y="top", legend=c(paste(names(dd3), sapply(dd3, length), sep=": "), "", sprintf("K-W p = %.1E", kwt$p.value)), bg="white")
				}
				dev.off()
			}
		}
		
		## boxplots with points
		cat(sprintf(" -> boxplots2"))
		sbt.models <- c("ssp2003", "ssp2006", "pam50", "scmgene", "scmod1", "scmod2")
		for(j in 1:ncol(mysig)) {
			## for each signature
			if(!all(!complete.cases(mysig[ ,j]))) {
				pdf(sprintf("%s/%s/boxplotplus2_%s_%s.pdf", saveres, boxres, dimnames(mysig)[[2]][j], dn[length(dn)]), width=9, height=9)
				for(k in 1:length(sbt.models)) {
					## for each subtypes identifications
					dd2 <- NULL
					## create a list with the signature scores for each subtype
					oo <- sort(unique(as.character(rest[ ,sbt.models[k]])))
					for(m in 1:length(oo)) { dd2 <- c(dd2, list(mysig[rest[ ,sbt.models[k]] == oo[m], j])) }
					names(dd2) <- oo
					mycol <- mypch <- NULL
					if(is.element(sbt.models[k], c("ssp2003", "ssp2006", "pam50"))) {
						oo <- oo[match(sbtn.intrinsic, oo)]
						mycolg <- colg.intrinsic[oo]
						mypchg <- pchg.intrinsic[oo]
						dd2 <- dd2[sbtn.intrinsic]
					} else {
						oo <- oo[match(sbtn, oo)]
						mycolg <- colg[oo]
						mypchg <- pchg[oo]
						dd2 <- dd2[sbtn]
					}
					for(m in 1:length(oo)) {
					mycol <- c(mycol, rep(mycolg[oo[m]], length(dd2[[m]])))
					mypch <- c(mypch,  rep(mypchg[oo[m]], length(dd2[[m]])))
				}
				kwt <- NULL
				uu <- unique(as.character(rest[ ,sbt.models[k]]))
				uu <- uu[!is.na(uu)]
				if(length(uu) > 1) { kwt <- kruskal.test(mysig[ ,j] ~ as.factor(rest[ ,sbt.models[k]])) } else { kwt$p.value <- NA }
				dd3 <- c(list(" "=NA), list("  "=NA), dd2)
				par(mar=c(10.1,4.1,4.1,2.1))
				rr <- boxplotplus2(x=dd3, .las=3, outline=F, .jit=0.75, .ylim=c(-2,2), pt.cex=0.75, pt.col=c(NA, NA, mycol), pt.pch=c(NA, NA, mypch), main=sprintf("%s in %s\n%s", dimnames(mysig)[[2]][j], dn[length(dn)], sbt.models[k]))
				rr <- rr[rr != 0]
				smartlegend(x="left", y="top", legend=c(paste(names(rr), rr, sep=": "), "", sprintf("K-W p = %.1E", kwt$p.value)), bg="white", pch=c(mypchg, NA, NA), col=c(mycolg, "white", "white"))
				}
				dev.off()
			}
		}
	
		cat(sprintf("\n"))
	}
}


###################################################
### code chunk number 8: mandatasettop
###################################################
dnt <- "TOP"
if(is.element(dnt, ddn)) {
	rest <- sbt.probat <- NULL
	cat(sprintf("%s\t", dnt))
	dn <- c(dn, dnt)

	## load data
	load(sprintf("%s/%s/%s.RData", rdp, ddn.db[dnt,"location"], ddn.db[dnt,"location"]))

	## collect clinical information
	demo.all <- rbind(demo.all, apply(X=demo[ ,cinfo,drop=FALSE], MARGIN=2, FUN=as.character))

	## bimodality of ESR1 and ERBB2
	cat(sprintf(" -> bimodality"))
	## bimodality for gene alone
	rest <- cbind(rest, "bimod.ESR1.gene"=bimod(x=mod1$ESR1[1, ,drop=FALSE], data=data, annot=annot, do.mapping=ifelse(ddn.db[dnt,"platform"] == "affy", FALSE, TRUE), model="V")$status)
	rest <- cbind(rest, "bimod.ERBB2.gene"=bimod(x=mod1$ERBB2[1, ,drop=FALSE], data=data, annot=annot, do.mapping=ifelse(ddn.db[dnt,"platform"] == "affy", FALSE, TRUE), model="V")$status)
	## bimodality for ESR1 and ERBB2 mod1
	rest <- cbind(rest, "bimod.ESR1.mod1"=NA)
	rest <- cbind(rest, "bimod.ERBB2.mod1"=bimod(x=mod1$ERBB2, data=data, annot=annot, do.mapping=ifelse(ddn.db[dnt,"platform"] == "affy", FALSE, TRUE), model="V")$status)
	## bimodality for ESR1 and ERBB2 mod2
	rest <- cbind(rest, "bimod.ESR1.mod2"=NA)
	rest <- cbind(rest, "bimod.ERBB2.mod2"=bimod(x=mod2$ERBB2, data=data, annot=annot, do.mapping=TRUE, model="V")$status)	

	## subtypes identification	
	cat(sprintf(" -> subtypes"))
	## SSP2003
	tt <- intrinsic.cluster.predict(sbt.model=ssp2003.robust, data=data, annot=annot, do.mapping=TRUE, do.prediction.strength=FALSE)
	rest <- cbind(rest, "ssp2003"=tt$subtype)
	dimnames(tt$subtype.proba)[[2]] <- paste("ssp2003", dimnames(tt$subtype.proba)[[2]], sep=".")
	sbt.probat <- cbind(sbt.probat, "ssp2003"=tt$subtype.proba)
	## SSP2006
	tt <- intrinsic.cluster.predict(sbt.model=ssp2006.robust, data=data, annot=annot, do.mapping=TRUE, do.prediction.strength=FALSE)
	rest <- cbind(rest, "ssp2006"=tt$subtype)
	dimnames(tt$subtype.proba)[[2]] <- paste("ssp2006", dimnames(tt$subtype.proba)[[2]], sep=".")
	sbt.probat <- cbind(sbt.probat, "ssp2006"=tt$subtype.proba)
	## PAM50
	tt <- intrinsic.cluster.predict(sbt.model=pam50.robust, data=data, annot=annot, do.mapping=TRUE, do.prediction.strength=FALSE)
	rest <- cbind(rest, "pam50"=tt$subtype)
	dimnames(tt$subtype.proba)[[2]] <- paste("pam50", dimnames(tt$subtype.proba)[[2]], sep=".")
	sbt.probat <- cbind(sbt.probat, "pam50"=tt$subtype.proba)
	## SCMGENE
	rest <- cbind(rest, "scmgene"=NA)
	sbt.probat2 <- matrix(NA, nrow=nrow(data), ncol=length(sbtn), dimnames=list(dimnames(data)[[1]], paste("scmgene", sbtn, sep=".")))
	sbt.probat <- cbind(sbt.probat, "scmgene"=sbt.probat2)
	## SCMOD1
	rest <- cbind(rest, "sbtgmod1"=NA)
	sbt.probat2 <- matrix(NA, nrow=nrow(data), ncol=length(sbtn), dimnames=list(dimnames(data)[[1]], paste("scmgene", sbtn, sep=".")))
	sbt.probat <- cbind(sbt.probat, "scmod1"=sbt.probat2)
	## SCMOD2
	rest <- cbind(rest, "scmgene"=NA)
	sbt.probat2 <- matrix(NA, nrow=nrow(data), ncol=length(sbtn), dimnames=list(dimnames(data)[[1]], paste("scmod2", sbtn, sep=".")))
	sbt.probat <- cbind(sbt.probat, "scmod2"=sbt.probat2)

	## signature scores computations
	cat(sprintf(" -> signatures"))
	mysig <- mymap <- NULL
	for(j in 1:length(glist2)) {
		tt <- sig.score(x=glist2[[j]], data=data, annot=annot, do.mapping=domap.sigs, signed=TRUE)
		mymap <- c(mymap, tt$mappin[1])
		mysig <- cbind(mysig, (rescale(tt$score, q=resq, na.rm=TRUE) - 0.5) * 2)
	}
	dimnames(mysig)[[2]] <- names(glist2)
	rest <- cbind(rest, mysig)
	
	map.list.all <- rbind(mymap, map.list.all)
	res <- rbind(res, rest)
	sbt.proba <- rbind(sbt.proba, sbt.probat)

	## boxplots
	cat(sprintf(" -> boxplots"))
	sbt.models <- c("ssp2003", "ssp2006", "pam50")
	for(j in 1:ncol(mysig)) {
		## for each signature
		if(!all(!complete.cases(mysig[ ,j]))) {
			pdf(sprintf("%s/%s/boxplotplus2_%s_%s.pdf", saveres, boxres, dimnames(mysig)[[2]][j], dnt), width=9, height=9)
			for(k in 1:length(sbt.models)) {
				## for each subtypes identifications
				dd2 <- NULL
				## create a list with the signature scores for each subtype
				oo <- sort(unique(as.character(rest[ ,sbt.models[k]])))
				for(m in 1:length(oo)) { dd2 <- c(dd2, list(mysig[rest[ ,sbt.models[k]] == oo[m], j])) }
				names(dd2) <- oo
				mycol <- mypch <- NULL
				if(is.element(sbt.models[k], c("ssp2003", "ssp2006", "pam50"))) {
					oo <- oo[match(sbtn.intrinsic, oo)]
					mycolg <- colg.intrinsic[oo]
					mypchg <- pchg.intrinsic[oo]
					dd2 <- dd2[sbtn.intrinsic]
				} else {
					oo <- oo[match(sbtn, oo)]
					mycolg <- colg[oo]
					mypchg <- pchg[oo]
					dd2 <- dd2[sbtn]
				}
				for(m in 1:length(oo)) {
				mycol <- c(mycol, rep(mycolg[oo[m]], length(dd2[[m]])))
				mypch <- c(mypch,  rep(mypchg[oo[m]], length(dd2[[m]])))
			}
			kwt <- NULL
			uu <- unique(as.character(rest[ ,sbt.models[k]]))
			uu <- uu[!is.na(uu)]
			if(length(uu) > 1) { kwt <- kruskal.test(mysig[ ,j] ~ as.factor(rest[ ,sbt.models[k]])) } else { kwt$p.value <- NA }
			dd3 <- c(list(" "=NA), list("  "=NA), dd2)
			par(mar=c(10.1,4.1,4.1,2.1))
			rr <- boxplotplus2(x=dd3, .las=3, outline=F, .jit=0.75, .ylim=c(-2,2), pt.cex=0.75, pt.col=c(NA, NA, mycol), pt.pch=c(NA, NA, mypch), main=sprintf("%s in %s\n%s", dimnames(mysig)[[2]][j], dnt, sbt.models[k]))
			rr <- rr[rr != 0]
			smartlegend(x="left", y="top", legend=c(paste(names(rr), rr, sep=": "), "", sprintf("K-W p = %.1E", kwt$p.value)), bg="white", pch=c(mypchg, NA, NA), col=c(mycolg, "white", "white"))
			}
			dev.off()
		}
	}		
	cat(sprintf("\n"))
}


###################################################
### code chunk number 9: mandatasetmgh
###################################################
dnt <- "MGH"
if(is.element(dnt, ddn)) {
	rest <- sbt.probat <- NULL
	cat(sprintf("%s\t", dnt))
	dn <- c(dn, dnt)

	## load data
	load(sprintf("%s/%s/%s.RData", rdp, ddn.db[dnt,"location"], ddn.db[dnt,"location"]))

	## collect clinical information
	demo.all <- rbind(demo.all, apply(X=demo[ ,cinfo,drop=FALSE], MARGIN=2, FUN=as.character))

	## bimodality of ESR1 and ERBB2
	cat(sprintf(" -> bimodality"))
	## bimodality for gene alone
	rest <- cbind(rest, "bimod.ESR1.gene"=bimod(x=mod1$ESR1[1, ,drop=FALSE], data=data, annot=annot, do.mapping=ifelse(ddn.db[dnt,"platform"] == "affy", FALSE, TRUE), model="V")$status)
	rest <- cbind(rest, "bimod.ERBB2.gene"=bimod(x=mod1$ERBB2[1, ,drop=FALSE], data=data, annot=annot, do.mapping=ifelse(ddn.db[dnt,"platform"] == "affy", FALSE, TRUE), model="V")$status)
	## bimodality for ESR1 and ERBB2 mod1
	rest <- cbind(rest, "bimod.ESR1.mod1"=NA)
	rest <- cbind(rest, "bimod.ERBB2.mod1"=NA)
	## bimodality for ESR1 and ERBB2 mod2
	rest <- cbind(rest, "bimod.ESR1.mod2"=NA)
	rest <- cbind(rest, "bimod.ERBB2.mod2"=bimod(x=mod2$ERBB2, data=data, annot=annot, do.mapping=TRUE, model="V")$status)	

	## subtypes identification	
	cat(sprintf(" -> subtypes"))
	pdf(sprintf("%s/%s/subtypes_id_plot_%s.pdf", saveres, sbtplotres, ddn.db[dnt,"dataset"]), height=7, width=7)
	## SSP2003
	tt <- intrinsic.cluster.predict(sbt.model=ssp2003.robust, data=data, annot=annot, do.mapping=TRUE, do.prediction.strength=FALSE)
	rest <- cbind(rest, "ssp2003"=tt$subtype)
	dimnames(tt$subtype.proba)[[2]] <- paste("ssp2003", dimnames(tt$subtype.proba)[[2]], sep=".")
	sbt.probat <- cbind(sbt.probat, "ssp2003"=tt$subtype.proba)
	## SSP2006
	tt <- intrinsic.cluster.predict(sbt.model=ssp2006.robust, data=data, annot=annot, do.mapping=TRUE, do.prediction.strength=FALSE)
	rest <- cbind(rest, "ssp2006"=tt$subtype)
	dimnames(tt$subtype.proba)[[2]] <- paste("ssp2006", dimnames(tt$subtype.proba)[[2]], sep=".")
	sbt.probat <- cbind(sbt.probat, "ssp2006"=tt$subtype.proba)
	## PAM50
	tt <- intrinsic.cluster.predict(sbt.model=pam50.robust, data=data, annot=annot, do.mapping=TRUE, do.prediction.strength=FALSE)
	rest <- cbind(rest, "pam50"=tt$subtype)
	dimnames(tt$subtype.proba)[[2]] <- paste("pam50", dimnames(tt$subtype.proba)[[2]], sep=".")
	sbt.probat <- cbind(sbt.probat, "pam50"=tt$subtype.proba)
	## SCMGENE
	tt <- subtype.cluster.predict(sbt.model=scmgene.robust, data=data, annot=annot, do.mapping=ifelse(ddn.db[dnt,"platform"] == "affy", FALSE, TRUE), do.prediction.strength=FALSE, do.BIC=FALSE, plot=TRUE)
	title(main=sprintf("SCMGENE\n%s", ddn.db[dnt,"dataset"]))
	rest <- cbind(rest, "scmgene"=tt$subtype2)
	dimnames(tt$subtype.proba2)[[2]] <- paste("scmgene", dimnames(tt$subtype.proba2)[[2]], sep=".")
	sbt.probat <- cbind(sbt.probat, "scmgene"=tt$subtype.proba2)
	## SCMOD1
	tt <- subtype.cluster.predict(sbt.model=scmod1.robust, data=data, annot=annot, do.mapping=ifelse(ddn.db[dnt,"platform"] == "affy", FALSE, TRUE), do.prediction.strength=FALSE, do.BIC=FALSE, plot=TRUE)
	title(main=sprintf("SCMOD1\n%s", ddn.db[dnt,"dataset"]))
	rest <- cbind(rest, "scmod1"=tt$subtype2)
	dimnames(tt$subtype.proba2)[[2]] <- paste("scmod1", dimnames(tt$subtype.proba2)[[2]], sep=".")
	sbt.probat <- cbind(sbt.probat, "scmod1"=tt$subtype.proba2)
	## SCMOD2
	tt <- subtype.cluster.predict(sbt.model=scmod2.robust, data=data, annot=annot, do.mapping=TRUE, do.prediction.strength=FALSE, do.BIC=FALSE, plot=TRUE)
	title(main=sprintf("SCMOD2\n%s", ddn.db[dnt,"dataset"]))
	rest <- cbind(rest, "scmod2"=tt$subtype2)
	dimnames(tt$subtype.proba2)[[2]] <- paste("scmod2", dimnames(tt$subtype.proba2)[[2]], sep=".")
	sbt.probat <- cbind(sbt.probat, "scmod2"=tt$subtype.proba2)

	## signature scores computations
	cat(sprintf(" -> signatures"))
	mysig <- mymap <- NULL
	for(j in 1:length(glist2)) {
		tt <- sig.score(x=glist2[[j]], data=data, annot=annot, do.mapping=domap.sigs, signed=TRUE)
		mymap <- c(mymap, tt$mappin[1])
		mysig <- cbind(mysig, (rescale(tt$score, q=resq, na.rm=TRUE) - 0.5) * 2)
	}
	dimnames(mysig)[[2]] <- names(glist2)
	rest <- cbind(rest, mysig)
	map.list.all <- rbind(mymap, map.list.all)
	res <- rbind(res, rest)
	sbt.proba <- rbind(sbt.proba, sbt.probat)

	## boxplots
	cat(sprintf(" -> boxplots"))
	sbt.models <- c("ssp2003", "ssp2006", "pam50")
	for(j in 1:ncol(mysig)) {
		## for each signature
		if(!all(!complete.cases(mysig[ ,j]))) {
			pdf(sprintf("%s/%s/boxplotplus2_%s_%s.pdf", saveres, boxres, dimnames(mysig)[[2]][j], dnt), width=9, height=9)
			for(k in 1:length(sbt.models)) {
				## for each subtypes identifications
				dd2 <- NULL
				## create a list with the signature scores for each subtype
				oo <- sort(unique(as.character(rest[ ,sbt.models[k]])))
				for(m in 1:length(oo)) { dd2 <- c(dd2, list(mysig[rest[ ,sbt.models[k]] == oo[m], j])) }
				names(dd2) <- oo
				mycol <- mypch <- NULL
				if(is.element(sbt.models[k], c("ssp2003", "ssp2006", "pam50"))) {
					oo <- oo[match(sbtn.intrinsic, oo)]
					mycolg <- colg.intrinsic[oo]
					mypchg <- pchg.intrinsic[oo]
					dd2 <- dd2[sbtn.intrinsic]
				} else {
					oo <- oo[match(sbtn, oo)]
					mycolg <- colg[oo]
					mypchg <- pchg[oo]
					dd2 <- dd2[sbtn]
				}
				for(m in 1:length(oo)) {
				mycol <- c(mycol, rep(mycolg[oo[m]], length(dd2[[m]])))
				mypch <- c(mypch,  rep(mypchg[oo[m]], length(dd2[[m]])))
			}
			kwt <- NULL
			uu <- unique(as.character(rest[ ,sbt.models[k]]))
			uu <- uu[!is.na(uu)]
			if(length(uu) > 1) { kwt <- kruskal.test(mysig[ ,j] ~ as.factor(rest[ ,sbt.models[k]])) } else { kwt$p.value <- NA }
			dd3 <- c(list(" "=NA), list("  "=NA), dd2)
			par(mar=c(10.1,4.1,4.1,2.1))
			rr <- boxplotplus2(x=dd3, .las=3, outline=F, .jit=0.75, .ylim=c(-2,2), pt.cex=0.75, pt.col=c(NA, NA, mycol), pt.pch=c(NA, NA, mypch), main=sprintf("%s in %s\n%s", dimnames(mysig)[[2]][j], dnt, sbt.models[k]))
			rr <- rr[rr != 0]
			smartlegend(x="left", y="top", legend=c(paste(names(rr), rr, sep=": "), "", sprintf("K-W p = %.1E", kwt$p.value)), bg="white", pch=c(mypchg, NA, NA), col=c(mycolg, "white", "white"))
			}
			dev.off()
		}		
	}		
	cat(sprintf("\n"))
}


###################################################
### code chunk number 10: mandatasettam
###################################################
dnt <- "TAM"
if(is.element(dnt, ddn)) {
	rest <- sbt.probat <- NULL
	cat(sprintf("%s\t", dnt))
	dn <- c(dn, dnt)

	## load data
	load(sprintf("%s/%s/%s.RData", rdp, ddn.db[dnt,"location"], ddn.db[dnt,"location"]))

	## collect clinical information
	demo.all <- rbind(demo.all, apply(X=demo[ ,cinfo,drop=FALSE], MARGIN=2, FUN=as.character))

	## bimodality of ESR1 and ERBB2
	cat(sprintf(" -> bimodality"))
	## bimodality for gene alone
	#rest <- cbind(rest, "bimod.ESR1.gene"=NA)
	rest <- matrix(NA, nrow=nrow(data), ncol=1, dimnames=list(dimnames(data)[[1]], "bimod.ESR1.gene"))
	rest <- cbind(rest, "bimod.ERBB2.gene"=bimod(x=mod1$ERBB2[1, ,drop=FALSE], data=data, annot=annot, do.mapping=ifelse(ddn.db[dnt,"platform"] == "affy", FALSE, TRUE), model="V")$status)
	## bimodality for ESR1 and ERBB2 mod1
	rest <- cbind(rest, "bimod.ESR1.mod1"=NA)
	rest <- cbind(rest, "bimod.ERBB2.mod1"=bimod(x=mod1$ERBB2, data=data, annot=annot, do.mapping=ifelse(ddn.db[dnt,"platform"] == "affy", FALSE, TRUE), model="V")$status)
	## bimodality for ESR1 and ERBB2 mod2
	rest <- cbind(rest, "bimod.ESR1.mod2"=NA)
	rest <- cbind(rest, "bimod.ERBB2.mod2"=bimod(x=mod2$ERBB2, data=data, annot=annot, do.mapping=TRUE, model="V")$status)	

	## subtypes identification	
	cat(sprintf(" -> subtypes"))
	pdf(sprintf("%s/%s/subtypes_id_plot_%s.pdf", saveres, sbtplotres, ddn.db[dnt,"dataset"]), height=7, width=7)
	## SSP2003
	tt <- intrinsic.cluster.predict(sbt.model=ssp2003.robust, data=data, annot=annot, do.mapping=TRUE, do.prediction.strength=FALSE)
	rest <- cbind(rest, "ssp2003"=tt$subtype)
	dimnames(tt$subtype.proba)[[2]] <- paste("ssp2003", dimnames(tt$subtype.proba)[[2]], sep=".")
	sbt.probat <- cbind(sbt.probat, "ssp2003"=tt$subtype.proba)
	## SSP2006
	tt <- intrinsic.cluster.predict(sbt.model=ssp2006.robust, data=data, annot=annot, do.mapping=TRUE, do.prediction.strength=FALSE)
	rest <- cbind(rest, "ssp2006"=tt$subtype)
	dimnames(tt$subtype.proba)[[2]] <- paste("ssp2006", dimnames(tt$subtype.proba)[[2]], sep=".")
	sbt.probat <- cbind(sbt.probat, "ssp2006"=tt$subtype.proba)
	## PAM50
	tt <- intrinsic.cluster.predict(sbt.model=pam50.robust, data=data, annot=annot, do.mapping=TRUE, do.prediction.strength=FALSE)
	rest <- cbind(rest, "pam50"=tt$subtype)
	dimnames(tt$subtype.proba)[[2]] <- paste("pam50", dimnames(tt$subtype.proba)[[2]], sep=".")
	sbt.probat <- cbind(sbt.probat, "pam50"=tt$subtype.proba)
	## SCMGENE
	tt <- subtype.cluster.predict(sbt.model=scmgene.robust, data=data, annot=annot, do.mapping=ifelse(ddn.db[dnt,"platform"] == "affy", FALSE, TRUE), do.prediction.strength=FALSE, do.BIC=FALSE, plot=TRUE)
	title(main=sprintf("SCMGENE\n%s", ddn.db[dnt,"dataset"]))
	rest <- cbind(rest, "scmgene"=tt$subtype2)
	dimnames(tt$subtype.proba2)[[2]] <- paste("scmgene", dimnames(tt$subtype.proba2)[[2]], sep=".")
	sbt.probat <- cbind(sbt.probat, "scmgene"=tt$subtype.proba2)
	## SCMOD1
	tt <- subtype.cluster.predict(sbt.model=scmod1.robust, data=data, annot=annot, do.mapping=ifelse(ddn.db[dnt,"platform"] == "affy", FALSE, TRUE), do.prediction.strength=FALSE, do.BIC=FALSE, plot=TRUE)
	title(main=sprintf("SCMOD1\n%s", ddn.db[dnt,"dataset"]))
	rest <- cbind(rest, "scmod1"=tt$subtype2)
	dimnames(tt$subtype.proba2)[[2]] <- paste("scmod1", dimnames(tt$subtype.proba2)[[2]], sep=".")
	sbt.probat <- cbind(sbt.probat, "scmod1"=tt$subtype.proba2)
	## SCMOD2
	tt <- subtype.cluster.predict(sbt.model=scmod2.robust, data=data, annot=annot, do.mapping=TRUE, do.prediction.strength=FALSE, do.BIC=FALSE, plot=TRUE)
	title(main=sprintf("SCMOD2\n%s", ddn.db[dnt,"dataset"]))
	rest <- cbind(rest, "scmod2"=tt$subtype2)
	dimnames(tt$subtype.proba2)[[2]] <- paste("scmod2", dimnames(tt$subtype.proba2)[[2]], sep=".")
	sbt.probat <- cbind(sbt.probat, "scmod2"=tt$subtype.proba2)

	## signature scores computations
	cat(sprintf(" -> signatures"))
	mysig <- mymap <- NULL
	for(j in 1:length(glist2)) {
		tt <- sig.score(x=glist2[[j]], data=data, annot=annot, do.mapping=domap.sigs, signed=TRUE)
		mymap <- c(mymap, tt$mappin[1])
		mysig <- cbind(mysig, (rescale(tt$score, q=resq, na.rm=TRUE) - 0.5) * 2)
	}
	dimnames(mysig)[[2]] <- names(glist2)
	rest <- cbind(rest, mysig)
	map.list.all <- rbind(mymap, map.list.all)
	res <- rbind(res, rest)
	sbt.proba <- rbind(sbt.proba, sbt.probat)

	## boxplots
	cat(sprintf(" -> boxplots"))
	sbt.models <- c("ssp2003", "ssp2006", "pam50")
	for(j in 1:ncol(mysig)) {
		## for each signature
		if(!all(!complete.cases(mysig[ ,j]))) {
			pdf(sprintf("%s/%s/boxplotplus2_%s_%s.pdf", saveres, boxres, dimnames(mysig)[[2]][j], dnt), width=9, height=9)
			for(k in 1:length(sbt.models)) {
				## for each subtypes identifications
				dd2 <- NULL
				## create a list with the signature scores for each subtype
				oo <- sort(unique(as.character(rest[ ,sbt.models[k]])))
				for(m in 1:length(oo)) { dd2 <- c(dd2, list(mysig[rest[ ,sbt.models[k]] == oo[m], j])) }
				names(dd2) <- oo
				mycol <- mypch <- NULL
				if(is.element(sbt.models[k], c("ssp2003", "ssp2006", "pam50"))) {
					oo <- oo[match(sbtn.intrinsic, oo)]
					mycolg <- colg.intrinsic[oo]
					mypchg <- pchg.intrinsic[oo]
					dd2 <- dd2[sbtn.intrinsic]
				} else {
					oo <- oo[match(sbtn, oo)]
					mycolg <- colg[oo]
					mypchg <- pchg[oo]
					dd2 <- dd2[sbtn]
				}
				for(m in 1:length(oo)) {
				mycol <- c(mycol, rep(mycolg[oo[m]], length(dd2[[m]])))
				mypch <- c(mypch,  rep(mypchg[oo[m]], length(dd2[[m]])))
			}
			kwt <- NULL
			uu <- unique(as.character(rest[ ,sbt.models[k]]))
			uu <- uu[!is.na(uu)]
			if(length(uu) > 1) { kwt <- kruskal.test(mysig[ ,j] ~ as.factor(rest[ ,sbt.models[k]])) } else { kwt$p.value <- NA }
			dd3 <- c(list(" "=NA), list("  "=NA), dd2)
			par(mar=c(10.1,4.1,4.1,2.1))
			rr <- boxplotplus2(x=dd3, .las=3, outline=F, .jit=0.75, .ylim=c(-2,2), pt.cex=0.75, pt.col=c(NA, NA, mycol), pt.pch=c(NA, NA, mypch), main=sprintf("%s in %s\n%s", dimnames(mysig)[[2]][j], dnt, sbt.models[k]))
			rr <- rr[rr != 0]
			smartlegend(x="left", y="top", legend=c(paste(names(rr), rr, sep=": "), "", sprintf("K-W p = %.1E", kwt$p.value)), bg="white", pch=c(mypchg, NA, NA), col=c(mycolg, "white", "white"))
			}
			dev.off()
		}		
	}		
	cat(sprintf("\n"))
}


###################################################
### code chunk number 11: mandatasettam2
###################################################
dnt <- "TAM2"
if(is.element(dnt, ddn)) {
	rest <- sbt.probat <- NULL
	cat(sprintf("%s\t", dnt))
	dn <- c(dn, dnt)

	## load data
	load(sprintf("%s/%s/%s.RData", rdp, ddn.db[dnt,"location"], ddn.db[dnt,"location"]))

	## collect clinical information
	demo.all <- rbind(demo.all, apply(X=demo[ ,cinfo,drop=FALSE], MARGIN=2, FUN=as.character))

	## bimodality of ESR1 and ERBB2
	cat(sprintf(" -> bimodality"))
	## bimodality for gene alone
	rest <- cbind(rest, "bimod.ESR1.gene"=bimod(x=mod1$ESR1[1, ,drop=FALSE], data=data, annot=annot, do.mapping=ifelse(ddn.db[dnt,"platform"] == "affy", FALSE, TRUE), model="V")$status)
	rest <- cbind(rest, "bimod.ERBB2.gene"=bimod(x=mod1$ERBB2[1, ,drop=FALSE], data=data, annot=annot, do.mapping=ifelse(ddn.db[dnt,"platform"] == "affy", FALSE, TRUE), model="V")$status)
	## bimodality for ESR1 and ERBB2 mod1
	rest <- cbind(rest, "bimod.ESR1.mod1"=NA)
	rest <- cbind(rest, "bimod.ERBB2.mod1"=bimod(x=mod1$ERBB2, data=data, annot=annot, do.mapping=ifelse(ddn.db[dnt,"platform"] == "affy", FALSE, TRUE), model="V")$status)
	## bimodality for ESR1 and ERBB2 mod2
	rest <- cbind(rest, "bimod.ESR1.mod2"=NA)
	rest <- cbind(rest, "bimod.ERBB2.mod2"=bimod(x=mod2$ERBB2, data=data, annot=annot, do.mapping=TRUE, model="V")$status)	

	## subtypes identification	
	cat(sprintf(" -> subtypes"))
	pdf(sprintf("%s/%s/subtypes_id_plot_%s.pdf", saveres, sbtplotres, ddn.db[dnt,"dataset"]), height=7, width=7)
	## SSP2003
	tt <- intrinsic.cluster.predict(sbt.model=ssp2003.robust, data=data, annot=annot, do.mapping=TRUE, do.prediction.strength=FALSE)
	rest <- cbind(rest, "ssp2003"=tt$subtype)
	dimnames(tt$subtype.proba)[[2]] <- paste("ssp2003", dimnames(tt$subtype.proba)[[2]], sep=".")
	sbt.probat <- cbind(sbt.probat, "ssp2003"=tt$subtype.proba)
	## SSP2006
	tt <- intrinsic.cluster.predict(sbt.model=ssp2006.robust, data=data, annot=annot, do.mapping=TRUE, do.prediction.strength=FALSE)
	rest <- cbind(rest, "ssp2006"=tt$subtype)
	dimnames(tt$subtype.proba)[[2]] <- paste("ssp2006", dimnames(tt$subtype.proba)[[2]], sep=".")
	sbt.probat <- cbind(sbt.probat, "ssp2006"=tt$subtype.proba)
	## PAM50
	tt <- intrinsic.cluster.predict(sbt.model=pam50.robust, data=data, annot=annot, do.mapping=TRUE, do.prediction.strength=FALSE)
	rest <- cbind(rest, "pam50"=tt$subtype)
	dimnames(tt$subtype.proba)[[2]] <- paste("pam50", dimnames(tt$subtype.proba)[[2]], sep=".")
	sbt.probat <- cbind(sbt.probat, "pam50"=tt$subtype.proba)
	## SCMGENE
	tt <- subtype.cluster.predict(sbt.model=scmgene.robust, data=data, annot=annot, do.mapping=ifelse(ddn.db[dnt,"platform"] == "affy", FALSE, TRUE), do.prediction.strength=FALSE, do.BIC=FALSE, plot=TRUE)
	title(main=sprintf("SCMGENE\n%s", ddn.db[dnt,"dataset"]))
	rest <- cbind(rest, "scmgene"=tt$subtype2)
	dimnames(tt$subtype.proba2)[[2]] <- paste("scmgene", dimnames(tt$subtype.proba2)[[2]], sep=".")
	sbt.probat <- cbind(sbt.probat, "scmgene"=tt$subtype.proba2)
	## SCMOD1
	tt <- subtype.cluster.predict(sbt.model=scmod1.robust, data=data, annot=annot, do.mapping=ifelse(ddn.db[dnt,"platform"] == "affy", FALSE, TRUE), do.prediction.strength=FALSE, do.BIC=FALSE, plot=TRUE)
	title(main=sprintf("SCMOD1\n%s", ddn.db[dnt,"dataset"]))
	rest <- cbind(rest, "scmod1"=tt$subtype2)
	dimnames(tt$subtype.proba2)[[2]] <- paste("scmod1", dimnames(tt$subtype.proba2)[[2]], sep=".")
	sbt.probat <- cbind(sbt.probat, "scmod1"=tt$subtype.proba2)
	## SCMOD2
	tt <- subtype.cluster.predict(sbt.model=scmod2.robust, data=data, annot=annot, do.mapping=TRUE, do.prediction.strength=FALSE, do.BIC=FALSE, plot=TRUE)
	title(main=sprintf("SCMOD2\n%s", ddn.db[dnt,"dataset"]))
	rest <- cbind(rest, "scmod2"=tt$subtype2)
	dimnames(tt$subtype.proba2)[[2]] <- paste("scmod2", dimnames(tt$subtype.proba2)[[2]], sep=".")
	sbt.probat <- cbind(sbt.probat, "scmod2"=tt$subtype.proba2)

	## signature scores computations
	cat(sprintf(" -> signatures"))
	mysig <- mymap <- NULL
	for(j in 1:length(glist2)) {
		tt <- sig.score(x=glist2[[j]], data=data, annot=annot, do.mapping=domap.sigs, signed=TRUE)
		mymap <- c(mymap, tt$mappin[1])
		mysig <- cbind(mysig, (rescale(tt$score, q=resq, na.rm=TRUE) - 0.5) * 2)
	}
	dimnames(mysig)[[2]] <- names(glist2)
	rest <- cbind(rest, mysig)
	
	map.list.all <- rbind(mymap, map.list.all)
	res <- rbind(res, rest)
	sbt.proba <- rbind(sbt.proba, sbt.probat)

	## boxplots
	cat(sprintf(" -> boxplots"))
	sbt.models <- c("ssp2003", "ssp2006", "pam50")
	for(j in 1:ncol(mysig)) {		
		## for each signature
		if(!all(!complete.cases(mysig[ ,j]))) {
			pdf(sprintf("%s/%s/boxplotplus2_%s_%s.pdf", saveres, boxres, dimnames(mysig)[[2]][j], dnt), width=9, height=9)
			for(k in 1:length(sbt.models)) {
				## for each subtypes identifications
				dd2 <- NULL
				## create a list with the signature scores for each subtype
				oo <- sort(unique(as.character(rest[ ,sbt.models[k]])))
				for(m in 1:length(oo)) { dd2 <- c(dd2, list(mysig[rest[ ,sbt.models[k]] == oo[m], j])) }
				names(dd2) <- oo
				mycol <- mypch <- NULL
				if(is.element(sbt.models[k], c("ssp2003", "ssp2006", "pam50"))) {
					oo <- oo[match(sbtn.intrinsic, oo)]
					mycolg <- colg.intrinsic[oo]
					mypchg <- pchg.intrinsic[oo]
					dd2 <- dd2[sbtn.intrinsic]
				} else {
					oo <- oo[match(sbtn, oo)]
					mycolg <- colg[oo]
					mypchg <- pchg[oo]
					dd2 <- dd2[sbtn]
				}
				for(m in 1:length(oo)) {
				mycol <- c(mycol, rep(mycolg[oo[m]], length(dd2[[m]])))
				mypch <- c(mypch,  rep(mypchg[oo[m]], length(dd2[[m]])))
			}
			kwt <- NULL
			uu <- unique(as.character(rest[ ,sbt.models[k]]))
			uu <- uu[!is.na(uu)]
			if(length(uu) > 1) { kwt <- kruskal.test(mysig[ ,j] ~ as.factor(rest[ ,sbt.models[k]])) } else { kwt$p.value <- NA }
			dd3 <- c(list(" "=NA), list("  "=NA), dd2)
			par(mar=c(10.1,4.1,4.1,2.1))
			rr <- boxplotplus2(x=dd3, .las=3, outline=F, .jit=0.75, .ylim=c(-2,2), pt.cex=0.75, pt.col=c(NA, NA, mycol), pt.pch=c(NA, NA, mypch), main=sprintf("%s in %s\n%s", dimnames(mysig)[[2]][j], dnt, sbt.models[k]))
			rr <- rr[rr != 0]
			smartlegend(x="left", y="top", legend=c(paste(names(rr), rr, sep=": "), "", sprintf("K-W p = %.1E", kwt$p.value)), bg="white", pch=c(mypchg, NA, NA), col=c(mycolg, "white", "white"))
			}
			dev.off()
		}		
	}		
	cat(sprintf("\n"))
}


###################################################
### code chunk number 12: mandatasetvdx3
###################################################
dnt <- "VDX3"
if(is.element(dnt, ddn)) {
	rest <- sbt.probat <- NULL
	cat(sprintf("%s\t", dnt))
	dn <- c(dn, dnt)

	## load data
	load(sprintf("%s/%s/%s.RData", rdp, ddn.db[dnt,"location"], ddn.db[dnt,"location"]))

	## collect clinical information
	demo.all <- rbind(demo.all, apply(X=demo[ ,cinfo,drop=FALSE], MARGIN=2, FUN=as.character))

	## bimodality of ESR1 and ERBB2
	cat(sprintf(" -> bimodality"))
	## bimodality for gene alone
	#rest <- cbind(rest, "bimod.ESR1.gene"=NA)
	rest <- matrix(NA, nrow=nrow(data), ncol=1, dimnames=list(dimnames(data)[[1]], "bimod.ESR1.gene"))
	rest <- cbind(rest, "bimod.ERBB2.gene"=NA)
	## bimodality for ESR1 and ERBB2 mod1
	rest <- cbind(rest, "bimod.ESR1.mod1"=NA)
	rest <- cbind(rest, "bimod.ERBB2.mod1"=NA)
	## bimodality for ESR1 and ERBB2 mod2
	rest <- cbind(rest, "bimod.ESR1.mod2"=NA)
	rest <- cbind(rest, "bimod.ERBB2.mod2"=NA)	

	## subtypes identification	
	cat(sprintf(" -> subtypes"))
	## SSP2003
	tt <- intrinsic.cluster.predict(sbt.model=ssp2003.robust, data=data, annot=annot, do.mapping=TRUE, do.prediction.strength=FALSE)
	rest <- cbind(rest, "ssp2003"=tt$subtype)
	dimnames(tt$subtype.proba)[[2]] <- paste("ssp2003", dimnames(tt$subtype.proba)[[2]], sep=".")
	sbt.probat <- cbind(sbt.probat, "ssp2003"=tt$subtype.proba)
	## SSP2006
	tt <- intrinsic.cluster.predict(sbt.model=ssp2006.robust, data=data, annot=annot, do.mapping=TRUE, do.prediction.strength=FALSE)
	rest <- cbind(rest, "ssp2006"=tt$subtype)
	dimnames(tt$subtype.proba)[[2]] <- paste("ssp2006", dimnames(tt$subtype.proba)[[2]], sep=".")
	sbt.probat <- cbind(sbt.probat, "ssp2006"=tt$subtype.proba)
	## PAM50
	tt <- intrinsic.cluster.predict(sbt.model=pam50.robust, data=data, annot=annot, do.mapping=TRUE, do.prediction.strength=FALSE)
	rest <- cbind(rest, "pam50"=tt$subtype)
	dimnames(tt$subtype.proba)[[2]] <- paste("pam50", dimnames(tt$subtype.proba)[[2]], sep=".")
	sbt.probat <- cbind(sbt.probat, "pam50"=tt$subtype.proba)
	## SCMGENE
	rest <- cbind(rest, "scmgene"=NA)
	sbt.probat2 <- matrix(NA, nrow=nrow(data), ncol=length(sbtn), dimnames=list(dimnames(data)[[1]], paste("scmgene", sbtn, sep=".")))
	sbt.probat <- cbind(sbt.probat, "scmgene"=sbt.probat2)
	## SCMOD1
	rest <- cbind(rest, "sbtgmod1"=NA)
	sbt.probat2 <- matrix(NA, nrow=nrow(data), ncol=length(sbtn), dimnames=list(dimnames(data)[[1]], paste("scmgene", sbtn, sep=".")))
	sbt.probat <- cbind(sbt.probat, "scmod1"=sbt.probat2)
	## SCMOD2
	rest <- cbind(rest, "scmgene"=NA)
	sbt.probat2 <- matrix(NA, nrow=nrow(data), ncol=length(sbtn), dimnames=list(dimnames(data)[[1]], paste("scmod2", sbtn, sep=".")))
	sbt.probat <- cbind(sbt.probat, "scmod2"=sbt.probat2)

	## signature scores computations
	cat(sprintf(" -> signatures"))
	mysig <- mymap <- NULL
	for(j in 1:length(glist2)) {
		tt <- sig.score(x=glist2[[j]], data=data, annot=annot, do.mapping=domap.sigs, signed=TRUE)
		mymap <- c(mymap, tt$mappin[1])
		mysig <- cbind(mysig, (rescale(tt$score, q=resq, na.rm=TRUE) - 0.5) * 2)
	}
	dimnames(mysig)[[2]] <- names(glist2)
	rest <- cbind(rest, mysig)
	
	map.list.all <- rbind(mymap, map.list.all)
	res <- rbind(res, rest)
	sbt.proba <- rbind(sbt.proba, sbt.probat)

	## boxplots
	cat(sprintf(" -> boxplots"))
	sbt.models <- c("ssp2003", "ssp2006", "pam50")
	for(j in 1:ncol(mysig)) {		
		## for each signature
		if(!all(!complete.cases(mysig[ ,j]))) {
			pdf(sprintf("%s/%s/boxplotplus2_%s_%s.pdf", saveres, boxres, dimnames(mysig)[[2]][j], dnt), width=9, height=9)
			for(k in 1:length(sbt.models)) {
				## for each subtypes identifications
				dd2 <- NULL
				## create a list with the signature scores for each subtype
				oo <- sort(unique(as.character(rest[ ,sbt.models[k]])))
				for(m in 1:length(oo)) { dd2 <- c(dd2, list(mysig[rest[ ,sbt.models[k]] == oo[m], j])) }
				names(dd2) <- oo
				mycol <- mypch <- NULL
				if(is.element(sbt.models[k], c("ssp2003", "ssp2006", "pam50"))) {
					oo <- oo[match(sbtn.intrinsic, oo)]
					mycolg <- colg.intrinsic[oo]
					mypchg <- pchg.intrinsic[oo]
					dd2 <- dd2[sbtn.intrinsic]
				} else {
					oo <- oo[match(sbtn, oo)]
					mycolg <- colg[oo]
					mypchg <- pchg[oo]
					dd2 <- dd2[sbtn]
				}
				for(m in 1:length(oo)) {
				mycol <- c(mycol, rep(mycolg[oo[m]], length(dd2[[m]])))
				mypch <- c(mypch,  rep(mypchg[oo[m]], length(dd2[[m]])))
			}
			kwt <- NULL
			uu <- unique(as.character(rest[ ,sbt.models[k]]))
			uu <- uu[!is.na(uu)]
			if(length(uu) > 1) { kwt <- kruskal.test(mysig[ ,j] ~ as.factor(rest[ ,sbt.models[k]])) } else { kwt$p.value <- NA }
			dd3 <- c(list(" "=NA), list("  "=NA), dd2)
			par(mar=c(10.1,4.1,4.1,2.1))
			rr <- boxplotplus2(x=dd3, .las=3, outline=F, .jit=0.75, .ylim=c(-2,2), pt.cex=0.75, pt.col=c(NA, NA, mycol), pt.pch=c(NA, NA, mypch), main=sprintf("%s in %s\n%s", dimnames(mysig)[[2]][j], dnt, sbt.models[k]))
			rr <- rr[rr != 0]
			smartlegend(x="left", y="top", legend=c(paste(names(rr), rr, sep=": "), "", sprintf("K-W p = %.1E", kwt$p.value)), bg="white", pch=c(mypchg, NA, NA), col=c(mycolg, "white", "white"))
			}
			dev.off()
		}		
	}		
	cat(sprintf("\n"))
}


###################################################
### code chunk number 13: mandatasetlh
###################################################
dnt <- "LH"
if(is.element(dnt, ddn)) {
	cat(sprintf("%s\t", ddn.db[dnt,"dataset"]))
	rest <- sbt.probat <- NULL
	dn <- c(dn, ddn.db[dnt,"dataset"])
	## the dataset does not require specific code
	
	## load data
	load(sprintf("%s/%s/%s.RData", rdp, ddn.db[dnt,"location"], ddn.db[dnt,"location"]))
	
	## collect clinical information
	demo.all <- rbind(demo.all, apply(X=demo[ ,cinfo,drop=FALSE], MARGIN=2, FUN=as.character))
	
	## bimodality of ESR1 and ERBB2
	cat(sprintf(" -> bimodality"))
	## bimodality for gene alone
	rest <- cbind(rest, "bimod.ESR1.gene"=bimod(x=mod1$ESR1[1, ,drop=FALSE], data=data, annot=annot, do.mapping=ifelse(ddn.db[dnt,"platform"] == "affy", FALSE, TRUE), model="V")$status)
	rest <- cbind(rest, "bimod.ERBB2.gene"=NA)
	## bimodality for ESR1 and ERBB2 mod1
	rest <- cbind(rest, "bimod.ESR1.mod1"=bimod(x=mod1$ESR1, data=data, annot=annot, do.mapping=ifelse(ddn.db[dnt,"platform"] == "affy", FALSE, TRUE), model="V")$status)
	rest <- cbind(rest, "bimod.ERBB2.mod1"=NA)
	## bimodality for ESR1 and ERBB2 mod2
	rest <- cbind(rest, "bimod.ESR1.mod2"=bimod(x=mod1$ESR1, data=data, annot=annot, do.mapping=TRUE, model="V")$status)
	rest <- cbind(rest, "bimod.ERBB2.mod2"=NA)	
	
	## subtypes identification	
	cat(sprintf(" -> subtypes"))
	## SSP2003
	tt <- intrinsic.cluster.predict(sbt.model=ssp2003.robust, data=data, annot=annot, do.mapping=TRUE, do.prediction.strength=FALSE)
	rest <- cbind(rest, "ssp2003"=tt$subtype)
	dimnames(tt$subtype.proba)[[2]] <- paste("ssp2003", dimnames(tt$subtype.proba)[[2]], sep=".")
	sbt.probat <- cbind(sbt.probat, "ssp2003"=tt$subtype.proba)
	## SSP2006
	tt <- intrinsic.cluster.predict(sbt.model=ssp2006.robust, data=data, annot=annot, do.mapping=TRUE, do.prediction.strength=FALSE)
	rest <- cbind(rest, "ssp2006"=tt$subtype)
	dimnames(tt$subtype.proba)[[2]] <- paste("ssp2006", dimnames(tt$subtype.proba)[[2]], sep=".")
	sbt.probat <- cbind(sbt.probat, "ssp2006"=tt$subtype.proba)
	## PAM50
	tt <- intrinsic.cluster.predict(sbt.model=pam50.robust, data=data, annot=annot, do.mapping=TRUE, do.prediction.strength=FALSE)
	rest <- cbind(rest, "pam50"=tt$subtype)
	dimnames(tt$subtype.proba)[[2]] <- paste("pam50", dimnames(tt$subtype.proba)[[2]], sep=".")
	sbt.probat <- cbind(sbt.probat, "pam50"=tt$subtype.proba)
	## SCMGENE
	rest <- cbind(rest, "scmgene"=NA)
	sbt.probat2 <- matrix(NA, nrow=nrow(data), ncol=length(sbtn), dimnames=list(dimnames(data)[[1]], paste("scmgene", sbtn, sep=".")))
	sbt.probat <- cbind(sbt.probat, "scmgene"=sbt.probat2)
	## SCMOD1
	rest <- cbind(rest, "sbtgmod1"=NA)
	sbt.probat2 <- matrix(NA, nrow=nrow(data), ncol=length(sbtn), dimnames=list(dimnames(data)[[1]], paste("scmgene", sbtn, sep=".")))
	sbt.probat <- cbind(sbt.probat, "scmod1"=sbt.probat2)
	## SCMOD2
	rest <- cbind(rest, "scmgene"=NA)
	sbt.probat2 <- matrix(NA, nrow=nrow(data), ncol=length(sbtn), dimnames=list(dimnames(data)[[1]], paste("scmod2", sbtn, sep=".")))
	sbt.probat <- cbind(sbt.probat, "scmod2"=sbt.probat2)
	
	## signature scores computations
	cat(sprintf(" -> signatures"))
	mysig <- mymap <- NULL
	for(j in 1:length(glist2)) {
		tt <- sig.score(x=glist2[[j]], data=data, annot=annot, do.mapping=domap.sigs, signed=TRUE)
		mymap <- c(mymap, tt$mappin[1])
		mysig <- cbind(mysig, (rescale(tt$score, q=resq, na.rm=TRUE) - 0.5) * 2)
	}
	dimnames(mysig)[[2]] <- names(glist2)
	rest <- cbind(rest, mysig)
	
	map.list.all <- rbind(mymap, map.list.all)
	res <- rbind(res, rest)
	sbt.proba <- rbind(sbt.proba, sbt.probat)
	
	## boxplots
	cat(sprintf(" -> boxplots"))
	sbt.models <- c("ssp2003", "ssp2006", "pam50")
	for(j in 1:ncol(mysig)) {		
		## for each signature
		if(!all(!complete.cases(mysig[ ,j]))) {
			pdf(sprintf("%s/%s/boxplotplus2_%s_%s.pdf", saveres, boxres, dimnames(mysig)[[2]][j], dnt), width=9, height=9)
			for(k in 1:length(sbt.models)) {
				## for each subtypes identifications
				dd2 <- NULL
				## create a list with the signature scores for each subtype
				oo <- sort(unique(as.character(rest[ ,sbt.models[k]])))
				for(m in 1:length(oo)) { dd2 <- c(dd2, list(mysig[rest[ ,sbt.models[k]] == oo[m], j])) }
				names(dd2) <- oo
				mycol <- mypch <- NULL
				if(is.element(sbt.models[k], c("ssp2003", "ssp2006", "pam50"))) {
					oo <- oo[match(sbtn.intrinsic, oo)]
					mycolg <- colg.intrinsic[oo]
					mypchg <- pchg.intrinsic[oo]
					dd2 <- dd2[sbtn.intrinsic]
				} else {
					oo <- oo[match(sbtn, oo)]
					mycolg <- colg[oo]
					mypchg <- pchg[oo]
					dd2 <- dd2[sbtn]
				}
				for(m in 1:length(oo)) {
				mycol <- c(mycol, rep(mycolg[oo[m]], length(dd2[[m]])))
				mypch <- c(mypch,  rep(mypchg[oo[m]], length(dd2[[m]])))
			}
			kwt <- NULL
			uu <- unique(as.character(rest[ ,sbt.models[k]]))
			uu <- uu[!is.na(uu)]
			if(length(uu) > 1) { kwt <- kruskal.test(mysig[ ,j] ~ as.factor(rest[ ,sbt.models[k]])) } else { kwt$p.value <- NA }
			dd3 <- c(list(" "=NA), list("  "=NA), dd2)
			par(mar=c(10.1,4.1,4.1,2.1))
			rr <- boxplotplus2(x=dd3, .las=3, outline=F, .jit=0.75, .ylim=c(-2,2), pt.cex=0.75, pt.col=c(NA, NA, mycol), pt.pch=c(NA, NA, mypch), main=sprintf("%s in %s\n%s", dimnames(mysig)[[2]][j], dnt, sbt.models[k]))
			rr <- rr[rr != 0]
			smartlegend(x="left", y="top", legend=c(paste(names(rr), rr, sep=": "), "", sprintf("K-W p = %.1E", kwt$p.value)), bg="white", pch=c(mypchg, NA, NA), col=c(mycolg, "white", "white"))
			}
			dev.off()
		}		
	}		
	cat(sprintf("\n"))
}


###################################################
### code chunk number 14: mandatasetmda3
###################################################
dnt <- "MDA3"
if(is.element(dnt, ddn)) {
	cat(sprintf("%s\t", ddn.db[dnt,"dataset"]))
	rest <- sbt.probat <- NULL
	dn <- c(dn, ddn.db[dnt,"dataset"])
	## the dataset does not require specific code
	
	## load data
	load(sprintf("%s/%s/%s.RData", rdp, ddn.db[dnt,"location"], ddn.db[dnt,"location"]))
	
	## collect clinical information
	demo.all <- rbind(demo.all, apply(X=demo[ ,cinfo,drop=FALSE], MARGIN=2, FUN=as.character))
	
	## bimodality of ESR1 and ERBB2
	cat(sprintf(" -> bimodality"))
	## bimodality for gene alone
	rest <- cbind(rest, "bimod.ESR1.gene"=bimod(x=mod1$ESR1[1, ,drop=FALSE], data=data, annot=annot, do.mapping=ifelse(ddn.db[dnt,"platform"] == "affy", FALSE, TRUE), model="V")$status)
	rest <- cbind(rest, "bimod.ERBB2.gene"=NA)
	## bimodality for ESR1 and ERBB2 mod1
	rest <- cbind(rest, "bimod.ESR1.mod1"=bimod(x=mod1$ESR1, data=data, annot=annot, do.mapping=ifelse(ddn.db[dnt,"platform"] == "affy", FALSE, TRUE), model="V")$status)
	rest <- cbind(rest, "bimod.ERBB2.mod1"=bimod(x=mod1$ERBB2, data=data, annot=annot, do.mapping=ifelse(ddn.db[dnt,"platform"] == "affy", FALSE, TRUE), model="V")$status)
	## bimodality for ESR1 and ERBB2 mod2
	rest <- cbind(rest, "bimod.ESR1.mod2"=bimod(x=mod1$ESR1, data=data, annot=annot, do.mapping=TRUE, model="V")$status)
	rest <- cbind(rest, "bimod.ERBB2.mod2"=NA)	
	
	## subtypes identification	
	cat(sprintf(" -> subtypes"))
	## SSP2003
	tt <- intrinsic.cluster.predict(sbt.model=ssp2003.robust, data=data, annot=annot, do.mapping=TRUE, do.prediction.strength=FALSE)
	rest <- cbind(rest, "ssp2003"=tt$subtype)
	dimnames(tt$subtype.proba)[[2]] <- paste("ssp2003", dimnames(tt$subtype.proba)[[2]], sep=".")
	sbt.probat <- cbind(sbt.probat, "ssp2003"=tt$subtype.proba)
	## SSP2006
	tt <- intrinsic.cluster.predict(sbt.model=ssp2006.robust, data=data, annot=annot, do.mapping=TRUE, do.prediction.strength=FALSE)
	rest <- cbind(rest, "ssp2006"=tt$subtype)
	dimnames(tt$subtype.proba)[[2]] <- paste("ssp2006", dimnames(tt$subtype.proba)[[2]], sep=".")
	sbt.probat <- cbind(sbt.probat, "ssp2006"=tt$subtype.proba)
	## PAM50
	tt <- intrinsic.cluster.predict(sbt.model=pam50.robust, data=data, annot=annot, do.mapping=TRUE, do.prediction.strength=FALSE)
	rest <- cbind(rest, "pam50"=tt$subtype)
	dimnames(tt$subtype.proba)[[2]] <- paste("pam50", dimnames(tt$subtype.proba)[[2]], sep=".")
	sbt.probat <- cbind(sbt.probat, "pam50"=tt$subtype.proba)
	## SCMGENE
	rest <- cbind(rest, "scmgene"=NA)
	sbt.probat2 <- matrix(NA, nrow=nrow(data), ncol=length(sbtn), dimnames=list(dimnames(data)[[1]], paste("scmgene", sbtn, sep=".")))
	sbt.probat <- cbind(sbt.probat, "scmgene"=sbt.probat2)
	## SCMOD1
	rest <- cbind(rest, "sbtgmod1"=NA)
	sbt.probat2 <- matrix(NA, nrow=nrow(data), ncol=length(sbtn), dimnames=list(dimnames(data)[[1]], paste("scmgene", sbtn, sep=".")))
	sbt.probat <- cbind(sbt.probat, "scmod1"=sbt.probat2)
	## SCMOD2
	rest <- cbind(rest, "scmgene"=NA)
	sbt.probat2 <- matrix(NA, nrow=nrow(data), ncol=length(sbtn), dimnames=list(dimnames(data)[[1]], paste("scmod2", sbtn, sep=".")))
	sbt.probat <- cbind(sbt.probat, "scmod2"=sbt.probat2)
	
	## signature scores computations
	cat(sprintf(" -> signatures"))
	mysig <- mymap <- NULL
	for(j in 1:length(glist2)) {
		tt <- sig.score(x=glist2[[j]], data=data, annot=annot, do.mapping=domap.sigs, signed=TRUE)
		mymap <- c(mymap, tt$mappin[1])
		mysig <- cbind(mysig, (rescale(tt$score, q=resq, na.rm=TRUE) - 0.5) * 2)
	}
	dimnames(mysig)[[2]] <- names(glist2)
	rest <- cbind(rest, mysig)
	
	map.list.all <- rbind(mymap, map.list.all)
	res <- rbind(res, rest)
	sbt.proba <- rbind(sbt.proba, sbt.probat)
	
	## boxplots
	cat(sprintf(" -> boxplots"))
	sbt.models <- c("ssp2003", "ssp2006", "pam50")
	for(j in 1:ncol(mysig)) {		
		## for each signature
		if(!all(!complete.cases(mysig[ ,j]))) {
			pdf(sprintf("%s/%s/boxplotplus2_%s_%s.pdf", saveres, boxres, dimnames(mysig)[[2]][j], dnt), width=9, height=9)
			for(k in 1:length(sbt.models)) {
				## for each subtypes identifications
				dd2 <- NULL
				## create a list with the signature scores for each subtype
				oo <- sort(unique(as.character(rest[ ,sbt.models[k]])))
				for(m in 1:length(oo)) { dd2 <- c(dd2, list(mysig[rest[ ,sbt.models[k]] == oo[m], j])) }
				names(dd2) <- oo
				mycol <- mypch <- NULL
				if(is.element(sbt.models[k], c("ssp2003", "ssp2006", "pam50"))) {
					oo <- oo[match(sbtn.intrinsic, oo)]
					mycolg <- colg.intrinsic[oo]
					mypchg <- pchg.intrinsic[oo]
					dd2 <- dd2[sbtn.intrinsic]
				} else {
					oo <- oo[match(sbtn, oo)]
					mycolg <- colg[oo]
					mypchg <- pchg[oo]
					dd2 <- dd2[sbtn]
				}
				for(m in 1:length(oo)) {
				mycol <- c(mycol, rep(mycolg[oo[m]], length(dd2[[m]])))
				mypch <- c(mypch,  rep(mypchg[oo[m]], length(dd2[[m]])))
			}
			kwt <- NULL
			uu <- unique(as.character(rest[ ,sbt.models[k]]))
			uu <- uu[!is.na(uu)]
			if(length(uu) > 1) { kwt <- kruskal.test(mysig[ ,j] ~ as.factor(rest[ ,sbt.models[k]])) } else { kwt$p.value <- NA }
			dd3 <- c(list(" "=NA), list("  "=NA), dd2)
			par(mar=c(10.1,4.1,4.1,2.1))
			rr <- boxplotplus2(x=dd3, .las=3, outline=F, .jit=0.75, .ylim=c(-2,2), pt.cex=0.75, pt.col=c(NA, NA, mycol), pt.pch=c(NA, NA, mypch), main=sprintf("%s in %s\n%s", dimnames(mysig)[[2]][j], dnt, sbt.models[k]))
			rr <- rr[rr != 0]
			smartlegend(x="left", y="top", legend=c(paste(names(rr), rr, sep=": "), "", sprintf("K-W p = %.1E", kwt$p.value)), bg="white", pch=c(mypchg, NA, NA), col=c(mycolg, "white", "white"))
			}
			dev.off()
		}		
	}		
	cat(sprintf("\n"))
}


###################################################
### code chunk number 15: poolboxplot
###################################################
## retrieve signature scores for all datasets
cat(sprintf("Pooled boxplots"))
sbt.models <- c("ssp2003", "ssp2006", "pam50", "scmgene", "scmod1", "scmod2")
mysig <- res[ , names(glist2), drop=FALSE]
mysbt <- res[ , sbt.models, drop=FALSE]
for(j in 1:ncol(mysig)) {
	## for each signature
	if(!all(!complete.cases(mysig[ ,j]))) {
		pdf(sprintf("%s/%s/boxplot_pool_%s.pdf", saveres, boxres, dimnames(mysig)[[2]][j]), width=9, height=9)
		for(k in 1:length(sbt.models)) {
			## for each subtypes identifications
			dd2 <- NULL
			## create a list with the signature scores for each subtype
			
			
			oo <- sort(unique(as.character(mysbt[ ,sbt.models[k]])))
			
			for(m in 1:length(oo)) { dd2 <- c(dd2, list(as.numeric(mysig[mysbt[ ,sbt.models[k]] == oo[m], j]))) }
			names(dd2) <- oo
			mycol <- mypch <- NULL
			if(is.element(sbt.models[k], c("ssp2003", "ssp2006", "pam50"))) {
				oo <- oo[match(sbtn.intrinsic, oo)]
				mycolg <- colg.intrinsic[oo]
				mypchg <- pchg.intrinsic[oo]
				dd2 <- dd2[sbtn.intrinsic]
			} else {
				oo <- oo[match(sbtn, oo)]
				mycolg <- colg[oo]
				mypchg <- pchg[oo]
				dd2 <- dd2[sbtn]
			}
			for(m in 1:length(oo)) {
			mycol <- c(mycol, rep(mycolg[oo[m]], length(dd2[[m]])))
			mypch <- c(mypch,  rep(mypchg[oo[m]], length(dd2[[m]])))
		}
		kwt <- NULL
		uu <- unique(as.character(mysbt[ ,sbt.models[k]]))
		uu <- uu[!is.na(uu)]
		if(length(uu) > 1) { kwt <- kruskal.test(mysig[ ,j] ~ as.factor(mysbt[ ,sbt.models[k]])) } else { kwt$p.value <- NA }
		dd3 <- c(list(" "=NA), list("  "=NA), dd2)
		par(mar=c(10.1,4.1,4.1,2.1))
		rr <- graphics::boxplot(x=dd3, las=3, outline=F, ylim=c(-2,2), main=sprintf("%s (all datasets)\n%s", dimnames(mysig)[[2]][j], sbt.models[k]), col="lightgrey")
		dd3 <- dd3[!is.na(dd3)]
		smartlegend(x="left", y="top", legend=c(paste(names(dd3), sapply(dd3, length), sep=": "), "", sprintf("K-W p = %.1E", kwt$p.value)), bg="white")
		}
		dev.off()
	}
}


###################################################
### code chunk number 16: poolboxplotplus2
###################################################
## retrieve signature scores for all datasets
cat(sprintf("Pooled boxplots"))
sbt.models <- c("ssp2003", "ssp2006", "pam50", "scmgene", "scmod1", "scmod2")
mysig <- res[ , names(glist2), drop=FALSE]
mysbt <- res[ , sbt.models, drop=FALSE]
for(j in 1:ncol(mysig)) {
	## for each signature
	if(!all(!complete.cases(mysig[ ,j]))) {
		pdf(sprintf("%s/%s/boxplotplus2_pool_%s.pdf", saveres, boxres, dimnames(mysig)[[2]][j]), width=9, height=9)
		for(k in 1:length(sbt.models)) {
			## for each subtypes identifications
			dd2 <- NULL
			## create a list with the signature scores for each subtype
			
			
			oo <- sort(unique(as.character(mysbt[ ,sbt.models[k]])))
			
			for(m in 1:length(oo)) { dd2 <- c(dd2, list(as.numeric(mysig[mysbt[ ,sbt.models[k]] == oo[m], j]))) }
			names(dd2) <- oo
			mycol <- mypch <- NULL
			if(is.element(sbt.models[k], c("ssp2003", "ssp2006", "pam50"))) {
				oo <- oo[match(sbtn.intrinsic, oo)]
				mycolg <- colg.intrinsic[oo]
				mypchg <- pchg.intrinsic[oo]
				dd2 <- dd2[sbtn.intrinsic]
			} else {
				oo <- oo[match(sbtn, oo)]
				mycolg <- colg[oo]
				mypchg <- pchg[oo]
				dd2 <- dd2[sbtn]
			}
			for(m in 1:length(oo)) {
			mycol <- c(mycol, rep(mycolg[oo[m]], length(dd2[[m]])))
			mypch <- c(mypch,  rep(mypchg[oo[m]], length(dd2[[m]])))
		}
		kwt <- NULL
		uu <- unique(as.character(mysbt[ ,sbt.models[k]]))
		uu <- uu[!is.na(uu)]
		if(length(uu) > 1) { kwt <- kruskal.test(mysig[ ,j] ~ as.factor(mysbt[ ,sbt.models[k]])) } else { kwt$p.value <- NA }
		dd3 <- c(list(" "=NA), list("  "=NA), dd2)
		par(mar=c(10.1,4.1,4.1,2.1))
		rr <- boxplotplus2(x=dd3, .las=3, outline=F, .jit=0.75, .ylim=c(-2,2), pt.cex=0.75, pt.col=c(NA, NA, mycol), pt.pch=c(NA, NA, mypch), main=sprintf("%s (all datasets)\n%s", dimnames(mysig)[[2]][j], sbt.models[k]))
		rr <- rr[rr != 0]
		smartlegend(x="left", y="top", legend=c(paste(names(rr), rr, sep=": "), "", sprintf("K-W p = %.1E", kwt$p.value)), bg="white", pch=c(mypchg, NA, NA), col=c(mycolg, "white", "white"))
		}
		dev.off()
	}
}


###################################################
### code chunk number 17: saveres
###################################################
## map.list.all contains the mapping for all signatures in each dataset
dimnames(map.list.all) <- list(c(dn[length(dn):1], "total"), names(glist2))
write.csv(t(map.list.all), file=sprintf("%s/mapping.csv", saveres))
## res contains all the information for each dataset and demo.all contains all the clinical information
write.csv(cbind(demo.all, res), row.names=FALSE, file=sprintf("%s/meta_demo.csv", saveres))
## sbt.proba contains the probabilities to belong to each subtype as computed by the subtype clustering models
write.csv(sbt.proba, file=sprintf("%s/subtype_proba.csv", saveres))


