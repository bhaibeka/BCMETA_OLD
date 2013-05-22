###################################################
###  myf
###################################################
##pairwcor subtype
myf.subtype <- function(x, subtype) {
dd <- x
for(i in 1:length(x)) {
	if(subtype[1] != "all") { myx <- is.element(x[[i]]$subtype, subtype) & !is.na(x[[i]]$subtype) }
	else { myx <- rep(TRUE, length(x[[i]]$subtype)) }
	dd[[i]]$clinic <- dd[[i]]$clinic[myx, , drop=FALSE]
	dd[[i]]$mscore <- dd[[i]]$mscore[myx, , drop=FALSE]
	dd[[i]]$subtype <- dd[[i]]$subtype[myx]
	dd[[i]]$bimod <- dd[[i]]$bimod[myx, , drop=FALSE]
}
return(dd)
}

##subtype info
myf.subtype.info <- function(x) {
	x <- summary(x$subtype)
	if(all(names(x) != "NA's")) { x <- c(x, "NA's"=0)}
	return(x)
}

##format data
myf.subtype2 <- function(x, subtype) {
	#myx <- complete.cases(x$clinic[ ,c(nc, "dataset", "series", "time", "event", x$mscore[ ,nm])
	myx <- rep(TRUE, nrow(x$clinic))
	if(length(subtype) > 1 || subtype != "all") {
		myx <- myx & is.element(x$subtype, subtype) & !is.na(x$subtype)
	}
	nn <- dimnames(x$clinic)[[1]][myx]
	mysubtype <- as.character(x$subtype[myx])
	names(mysubtype) <- nn
	mybimod <- x$bimod[myx, , drop=FALSE]
	dimnames(mybimod)[[1]] <- nn
	stra <- x$clinic[myx,"dataset"]
	names(stra) <- nn
	myout <- censor.time(surv.time=x$clinic[myx,"time"], surv.event=x$clinic[myx,"event"], time.cens=tc)
	myout <- cbind(myout[[1]], myout[[2]])
	dimnames(myout) <- list(nn, c("surv.time", "surv.event"))
	myin <- cbind(x$clinic[myx,nc], x$mscore[myx,nm])
	dimnames(myin)[[2]] <- c(nc, nm)
	return(list("input"=myin, "strata"=stra, "output"=myout, "subtype"=mysubtype, "bimod"=mybimod))
}

##discretized data
myf.subv <- function(x, y, subv) {
	y <- sort(y)
	temp <- x$input[ ,subv,drop=FALSE]
	mm <- apply(temp, 2, quantile, probs=y, na.rm=TRUE)
	for(i in 1:ncol(temp)) {
		if(!all(is.na(temp[ ,i]))) {
			ccc <- unique(mm[ , i])
			ccc <- ccc[!is.na(ccc)]
			if(length(ccc) == length(y)) {
				temp[ ,i] <- as.numeric(cut2(x=temp[ ,i], cuts=mm[ ,i], levels.mean=TRUE))
			} else { temp[ ,i] <- rep(NA, nrow(temp)) }
		}
	}
	x$input <- temp
	return(x)
}

##################################################
#function to estimate the relevance of each feature
#when added to an existing cox model (prior)
#
##data: list of datasets (for meta analysis). Each dataset is a list itself that contains an "input" item (matrix of input variables), a "strata" item (matrix used for stratification) and a "output" item (matrix with two columns, "surv.time" and "surv.event")
##prior: names of input variables that should be included in the cox models by default
##method: method used for estimating the relevance of a variable. "likelihood.ratio" method uses the likelihood ratio test (difference in loglik between the model before including the variable and after the inclusion). "cvpl" method uses the cross-validated partial likelihood.  
##################################################
surv.feature.rel <- function(data, prior=NULL, method=c("likelihood.ratio", "wald", "hr", "cvpl"), ...) {

	##############
	#internal functions
	##############
	
	internal.likratio <- function(priord, newd, surv.data, stratad) {
		res <- c(NA, NA)
		mm <- matrix(NA, nrow=ncol(priord) + ncol(newd), ncol=3)
		dimnames(mm) <- list(c(dimnames(priord)[[2]], dimnames(newd)[[2]]), c("coef", "coef.se", "n"))
		options(warn=2)
		if(ncol(priord) > 0) { rr1 <- try(coxph(Surv(surv.data[[1]], surv.data[[2]]) ~ strata(stratad) + as.matrix(priord))) }
		else { rr1 <- try(coxph(Surv(surv.data[[1]], surv.data[[2]]) ~ strata(stratad) + NULL)) }
		rr2 <- try(coxph(Surv(surv.data[[1]], surv.data[[2]]) ~ strata(stratad) + as.matrix(cbind(priord, newd))))
		options(warn=0)
		if(class(rr1) != "try-error" && class(rr2) != "try-error") {
			res <- c(2 * (rr2$loglik[length(rr2$loglik)] - rr1$loglik[length(rr1$loglik)]), rr2$n)
			mm[ ,1] <- rr2$coefficients
			mm[ ,2] <- sqrt(diag(rr2$var))
			mm[ ,3] <- rr2$n
		}
		names(res) <- c("stat", "n")
		
		return(list("x"=res, "model"=mm))
	}
	
	internal.wald <- function(priord, newd, surv.data, stratad) {
		require(MASS) || stop("library MASS is not available!")
		res <- c(NA, NA)
		mm <- matrix(NA, nrow=ncol(priord) + ncol(newd), ncol=3)
		dimnames(mm) <- list(c(dimnames(priord)[[2]], dimnames(newd)[[2]]), c("coef", "coef.se", "n"))
		options(warn=2)
		rr2 <- try(coxph(Surv(surv.data[[1]], surv.data[[2]]) ~ strata(stratad) + as.matrix(cbind(priord, newd))))
		options(warn=0)
		if(class(rr2) != "try-error") {
			res <- c(drop(t(rr2$coefficients) %*% ginv(rr2$var) %*% rr2$coefficients), rr2$n)
			mm[ ,1] <- rr2$coefficients
			mm[ ,2] <- sqrt(diag(rr2$var))
			mm[ ,3] <- rr2$n
		}
		names(res) <- c("stat", "n")
		
		return(list("x"=res, "model"=mm))
	}
	
		internal.hr <- function(priord, newd, surv.data, stratad) {
		require(MASS) || stop("library MASS is not available!")
		res <- c(NA, NA, NA)
		nv <- ncol(priord) + ncol(newd)
		mm <- matrix(NA, nrow=nv, ncol=3)
		dimnames(mm) <- list(c(dimnames(priord)[[2]], dimnames(newd)[[2]]), c("coef", "coef.se", "n"))
		options(warn=2)
		rr2 <- try(coxph(Surv(surv.data[[1]], surv.data[[2]]) ~ strata(stratad) + as.matrix(cbind(priord, newd))))
		options(warn=0)
		if(class(rr2) != "try-error") {
			res <- c(rr2$coefficients[nv], sqrt(diag(rr2$var)[nv]), rr2$n)
			mm[ ,1] <- rr2$coefficients
			mm[ ,2] <- sqrt(diag(rr2$var))
			mm[ ,3] <- rr2$n
		}
		names(res) <- c("x", "x.se", "n")
		
		return(list("x"=res, "model"=mm))
	}
	
	internal.cvpl <- function(priord, newd, surv.data, stratad, ...) {
		res <- c(NA, NA, NA)
		mm <- matrix(NA, nrow=ncol(priord) + ncol(newd), ncol=3)
		dimnames(mm) <- list(c(dimnames(priord)[[2]], dimnames(newd)[[2]]), c("coef", "coef.se", "n"))
		#if(ncol(priord) == 0) { priord <- NULL }
		rr <- cvpl(data=newd, surv.time=surv.data[[1]], surv.event=surv.data[[2]], prior=priord, strata.cox=stratad, ...)
		if(all(res$convergence)) {
			res <- c(mean(rr$pl, na.rm=TRUE), sd(rr$pl, na.rm=TRUE), rr$n)
			options(warn=2)
			rr2 <- try(coxph(Surv(surv.data[[1]], surv.data[[2]]) ~ strata(stratad) + as.matrix(cbind(priord, newd))))
			options(warn=0)
			if(class(rr2) != "try-error") {
				mm[ ,1] <- rr2$coefficients
				mm[ ,2] <- sqrt(diag(rr2$var))
				mm[ ,2] <- rr2$n
			}
		}
		names(res) <- c("x", "x.se", "n")
		
		return(list("x"=res, "model"=mm))
	}
	
	##############
	
	method <- match.arg(method)
	
	xxx <- NULL 
	#list with 2 items:
	# matrix with error estimates or statistics and standard errors (if error estimates)
	# 3d array with coefficients and their standard errors for each model
	xx <- NULL #meta estimators of xxx
	
	nnv <- lapply(data, function(x) { return(dimnames(x$input)[[2]]) })
	if(!all(nnv[[1]] == unique(unlist(nnv)))) { stop("the number/names of input variables are not the same over datasets!") }
	nnv <- nnv[[1]]
	if(!is.null(prior)) {
		if(!all(is.element(prior, nnv))) { stop("input variable names in prior are not valid!") }
		nnv <- nnv[!is.element(nnv, prior)]	 
	}
	nv <- length(nnv)
	if(nv == 0) { nnv <- NULL }
	nnd <- names(data)
	nd <- length(nnd)
	
	for(i in 1:nd) {
		xxd <- NULL
		xx.sed <- NULL
		xx.nd <- NULL
		xx2d <- array(NA, dim=c(length(prior) + 1, 3, nv))
		#3d array containing the coefficients 
		#and their standard errors for each model
		ifelse(nv == 0, nv2 <- 1, nv2 <- nv) #if no new data
		for(j in 1:nv2) {
			ee <- FALSE
			prior.idata <- data[[i]]$input[ ,prior,drop=FALSE]
			dimnames(prior.idata)[[2]] <- as.character(prior)
			new.idata <-  data[[i]]$input[ ,nnv[j],drop=FALSE]
			dimnames(new.idata) <- list(dimnames(data[[i]]$input)[[1]], nnv[j])
			surv.data <- list("surv.time"=data[[i]]$output[ ,"surv.time"], "surv.event"=data[[i]]$output[ ,"surv.event"])
			stra <- data[[i]]$strata
			if(is.null(stra)) {
				stra <- rep(1, nrow(new.idata))
				names(stra) <- dimnames(new.idata)[[1]]
			}
			cc.ix <- complete.cases(surv.data[[1]], surv.data[[2]], cbind(prior.idata, new.idata), stra)
			surv.data[[1]] <- surv.data[[1]][cc.ix]
			surv.data[[2]] <- surv.data[[2]][cc.ix]
			stra <- stra[cc.ix]
			new.idata <- new.idata[cc.ix, ,drop=FALSE]
			prior.idata <- prior.idata[cc.ix, ,drop=FALSE]
			if(all(!cc.ix) || nv == 0 || length(unique(new.idata[ ,1])) == 1) {
				xxd <- c(xxd, NA)
				xx.sed <- c(xx.sed, NA)
				xx.nd <- c(xx.nd, NA)
			}
			else {
				switch(method,
				"likelihood.ratio"={
					rr <- internal.likratio(priord=prior.idata, newd=new.idata, surv.data=surv.data, stratad=stra)
					xxd <- c(xxd, pchisq(rr$x[1], df=1, lower.tail=FALSE))
					xx.sed <- c(xx.sed, NA)
					xx.nd <- c(xx.nd, rr$x[2])
					xx2d[ , ,j] <- rr$model
				},
				"wald"={
					rr <- internal.wald(priord=prior.idata, newd=new.idata, surv.data=surv.data, stratad=stra)
					xxd <- c(xxd, pchisq(rr$x[1], df=length(prior) + 1, lower.tail=FALSE))
					xx.sed <- c(xx.sed, NA)
					xx.nd <- c(xx.nd, rr$x[2])
					xx2d[ , ,j] <- rr$model
				},
				"hr"={
					rr <- internal.hr(priord=prior.idata, newd=new.idata, surv.data=surv.data, stratad=stra, ...)
					xxd <- c(xxd, rr$x[1])
					xx.sed <- c(xx.sed, rr$x[2])
					xx.nd <- c(xx.nd, rr$x[3])
					xx2d[ , ,j] <- rr$model
				},
				"cvpl"={
					rr <- internal.cvpl(priord=prior.idata, newd=new.idata, surv.data=surv.data, stratad=stra, ...)
					xxd <- c(xxd, rr$x[1])
					xx.sed <- c(xx.sed, rr$x[2])
					xx.nd <- c(xx.nd, rr$x[3])
					xx2d[ , ,j] <- rr$model
				})
			}
		}
		names(xxd) <- names(xx.sed) <- names(xx.nd) <- nnv
		dimnames(xx2d) <- list(c(prior, "newvar"), c("coef", "coef.se", "n"), nnv)
		xxx <- c(xxx, list(list("estimate"=cbind("x"=xxd, "x.se"=xx.sed, "n"=xx.nd), "model"=xx2d)))
	}
	names(xxx) <- nnd
	
	#meta analysis
	if(nd == 1) {
		xx <- xxx[[1]]
		switch(method,
		"likelihood.ratio"={
			#nothing to do
		},
		"wald"={
			#nothing to do
		},
		"hr"={
			xx[[1]] <- cbind("x"=pchisq((xx[[1]][ ,1] / xx[[1]][ ,2])^2, df=1, lower.tail=FALSE), "x.se"=NA, "n"=xx[[1]][ ,3]) 
		},
		"cvpl"={
			#nothing to do
		})
	}
	else {
		if(nv == 0) {
			xx <- list(matrix(NA, nrow=1, ncol=3), NULL)
		}
		else {
			xx <- list(NULL, array(NA, dim=c(length(prior) + 1, 3, nv)))
			for(i in 1:nv) {
				my.x <- unlist(lapply(xxx, function(x, ix) { return(x[[1]][ix,1]) }, ix=i))
				my.x.se <- unlist(lapply(xxx, function(x, ix) { return(x[[1]][ix,2]) }, ix=i))
				my.n <- unlist(lapply(xxx, function(x, ix) { return(x[[1]][ix,3]) }, ix=i))
				my.model <-lapply(xxx, function(x, ix) { return(x[[2]][ , ,ix]) }, ix=i)
				#list of matrices describing the model in several datasets
				temp1 <- matrix(NA, nrow=length(prior) + 1, ncol=3)
				dimnames(temp1) <- list(c(prior, nnv[i]), c("coef", "coef.se", "n"))
				for(j in 1:(length(prior) + 1)) {
					my.x2 <- unlist(lapply(my.model, function(x, ix) { ifelse(is.matrix(x), return(x[ix,1]), return(x[1])) }, ix=j))
					my.x.se2 <- unlist(lapply(my.model, function(x, ix) { ifelse(is.matrix(x), return(x[ix,2]), return(x[2])) }, ix=j))
					my.n2 <- unlist(lapply(my.model, function(x, ix) { ifelse(is.matrix(x), return(x[ix,3]), return(x[3])) }, ix=j))
					rr <- combine.est(x=my.x2, x.se=my.x.se2, hetero=FALSE, na.rm=TRUE)
					temp1[j, ] <- c(unlist(rr), sum(my.n2, na.rm=TRUE))
				}
				xx[[2]][ , ,i] <- temp1
				dimnames(xx[[2]]) <- list(c(prior, "newvar"), c("coef", "coef.se", "n"), nnv)
				
				switch(method,
				"likelihood.ratio"={
					xx[[1]] <- rbind(xx[[1]], c("x"=combine.test(p=my.x, method="z.transform", na.rm=TRUE), "x.se"=NA, "n"=sum(my.n, na.rm=TRUE)))
				},
				"wald"={
					xx[[1]] <- rbind(xx[[1]], c("x"=combine.test(p=my.x, method="z.transform", na.rm=TRUE), "x.se"=NA, "n"=sum(my.n, na.rm=TRUE)))
				},
				"hr"={
					rr <- combine.est(x=my.x, x.se=my.x.se, hetero=FALSE, na.rm=TRUE)
					xx[[1]] <- rbind(xx[[1]], c("x"=pchisq((rr[[1]] / rr[[2]])^2, df=1, lower.tail=FALSE), "x.se"=NA, "n"=sum(my.n, na.rm=TRUE)))
				},
				"cvpl"={
					rr <- combine.est(x=my.x, x.se=my.x.se, hetero=FALSE, na.rm=TRUE)
					xx[[1]] <- rbind(xx[[1]], c("x"=rr[[1]], "x.se"=rr[[2]], "n"=sum(my.n, na.rm=TRUE)))
				})
			}
		}
		dimnames(xx[[1]]) <- list(nnv, c("x", "x.se", "n"))
		names(xx) <- c("estimate", "model")
	}
	
	return(list("all"=xxx, "meta"=xx))	
}

#################################################
#function inspired from forestplot (rmeta package)
##################################################
myforestplot<-function(labeltext,mean,lower,upper,align=NULL, is.summary=FALSE, clip=c(-Inf,Inf), xlab="", zero= 0,graphwidth=unit(2,"inches"),col=meta.colors(),xlog=FALSE, box.size=NULL, x.ticks=NULL, ...){

  require("grid") || stop("`grid' package not found")
  require("rmeta") || stop("`rmeta' package not found")

  ## Function to draw a non-summary rect-plus-CI
  drawNormalCI <- function(LL, OR, UL, size, bcol, lcol) {
    size=0.75*size

    clipupper<-convertX(unit(UL, "native"), "npc", valueOnly=TRUE) > 1
    cliplower<-convertX(unit(LL, "native"), "npc", valueOnly=TRUE) < 0
    box<- convertX(unit(OR, "native"), "npc", valueOnly=TRUE)
    clipbox <- box<0 || box>1
    
    
    ## Draw arrow if exceed col range
    ## convertX() used to convert between coordinate systems
    if (clipupper || cliplower){
      ends<-"both"
      lims<-unit(c(0, 1), c("npc", "npc"))
      if (!clipupper) {
        ends<-"first"
        lims<-unit(c(0, UL), c("npc","native"))
      }
      if (!cliplower) {
        ends<-"last"
        lims<-unit(c(LL, 1), c("native", "npc"))
      }
      grid.lines(x=lims, y=0.5,arrow=arrow(ends=ends,length=unit(0.05, "inches")),
                 gp=gpar(col=lcol))

      if (!clipbox)
          grid.rect(x=unit(OR, "native"),
                    width=unit(size, "snpc"), height=unit(size, "snpc"),
                    gp=gpar(fill=bcol,col=bcol))
      
      } else   {
      ## Draw line white if totally inside rect
      grid.lines(x=unit(c(LL, UL), "native"), y=0.5,
                 gp=gpar(col=lcol))
      grid.rect(x=unit(OR, "native"),
                width=unit(size, "snpc"), height=unit(size, "snpc"),
                gp=gpar(fill=bcol,col=bcol))
      if ((convertX(unit(OR, "native") + unit(0.5*size, "lines"), "native", valueOnly=TRUE) > UL) &&
          (convertX(unit(OR, "native") - unit(0.5*size, "lines"), "native", valueOnly=TRUE) < LL))
        grid.lines(x=unit(c(LL, UL), "native"), y=0.5, gp=gpar(col=lcol))
    }
  }
  
  ## Function to draw a summary "diamond"
  drawSummaryCI <- function(LL, OR, UL, size, scol) {
    grid.polygon(x=unit(c(LL, OR, UL, OR), "native"),
                 y=unit(0.5 + c(0, 0.5*size, 0, -0.5*size), "npc"),gp=gpar(fill=scol,col=scol))
  }
  
  plot.new()
  ## calculate width based on labels with something in every column
  widthcolumn<-!apply(is.na(labeltext),1,any)
  
  nc<-NCOL(labeltext)
  nr<-NROW(labeltext)
  labels<-vector("list",nc)
	if(length(col$lines) < nr) { col$lines <- rep(col$lines[1], nr) }
  	if(length(col$box) < nr) { col$box <- rep(col$box[1], nr) }
  	if(length(col$summary) < nr) { col$summary <- rep(col$summary[1], nr) }
  
  if (is.null(align))
    align<-c("l",rep("r",nc-1))
  else
    align<-rep(align,length=nc)
  
  is.summary<-rep(is.summary,length=nr)
  
  for(j in 1:nc){
    labels[[j]]<-vector("list", nr)
    for(i in 1:nr){
      if (is.na(labeltext[i,j]))
        next
      x<-switch(align[j],l=0,r=1,c=0.5)
      just<-switch(align[j],l="left",r="right",c="center")
      labels[[j]][[i]]<-textGrob(labeltext[i,j], x=x,just=just,
                                 gp=gpar(fontface=if(is.summary[i]) "bold" else "plain",
                                                                    col=rep(col$text,length=nr)[i]) )
    }
  }
  
  colgap<-unit(3,"mm")
  
  colwidths<-unit.c(max(unit(rep(1,sum(widthcolumn)),"grobwidth",labels[[1]][widthcolumn])),colgap)
  if (nc>1){
    for(i in 2:nc)
      colwidths<-unit.c(colwidths, max(unit(rep(1,sum(widthcolumn)),"grobwidth",labels[[i]][widthcolumn])),colgap)
    
  }
  colwidths<-unit.c(colwidths,graphwidth)
  
  pushViewport(viewport(layout=grid.layout(nr+1,nc*2+1,
                          widths=colwidths,
                          heights=unit(c(rep(1, nr),0.5), "lines"))))
  
  
  
  cwidth<-(upper-lower)

  xrange<-c(max(min(lower,na.rm=TRUE),clip[1]), min(max(upper,na.rm=TRUE),clip[2]))
  if(is.null(box.size)) {
  	info<-1/cwidth
  	info<-info/max(info[!is.summary], na.rm=TRUE)
  	info[is.summary]<-1
  }
  else { info <- box.size }
  for(j in 1:nc){
    for(i in 1:nr){
      if (!is.null(labels[[j]][[i]])){
        pushViewport(viewport(layout.pos.row=i,layout.pos.col=2*j-1))
        grid.draw(labels[[j]][[i]])
          popViewport()
      }
    }
  }
  
  pushViewport(viewport(layout.pos.col=2*nc+1, xscale=xrange))
  grid.lines(x=unit(zero, "native"), y=0:1,gp=gpar(col=col$zero))
  if (xlog){
    ticks<-pretty(exp(xrange))
    ticks<-ticks[ticks>0]
    if (min(lower,na.rm=TRUE)<clip[1]) ticks<-c(exp(clip[1]),ticks)
    if (max(upper,na.rm=TRUE)>clip[2]) ticks<-c(ticks,exp(clip[2]))
    xax<-xaxisGrob(gp=gpar(cex=0.6,col=col$axes),at=log(ticks),name="xax")
    xax1<-editGrob(xax, gPath("labels"), label=format(ticks,digits=2))
    grid.draw(xax1)
  } else {
    grid.xaxis(at=x.ticks, gp=gpar(cex=1,col=col$axes))
  }
  grid.text(xlab, y=unit(-3, "lines"),gp=gpar(col=col$axes))
  popViewport()
  for (i in 1:nr) {
    if (is.na(mean[i])) next
    pushViewport(viewport(layout.pos.row=i, layout.pos.col=2*nc+1,
                          xscale=xrange))
    if (is.summary[i])
      drawSummaryCI(lower[i], mean[i], upper[i], info[i], scol=col$summary[i])
    else
      drawNormalCI(lower[i], mean[i], upper[i], info[i], bcol=col$box[i], lcol=col$lines[i]) 
    popViewport()
  }
  popViewport()
}
