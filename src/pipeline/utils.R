## utils script
#======================================================================================================
#' Method to calculate reference interval (percentile) curves for objects of class 'RWDRICurve'
#' 
#' @param x				(object) of class 'RWDRICurve'
#' @param covarValue	(numeric) vector of covariate values for which prediction of percentiles shall be performed
#' @param RIperc		(numeric) value specifying the percentiles, which define the reference interval curves
#' @param Scale			(character) specifying if percentiles are calculated on the original scale ("Or") or the transformed scale ("Tr") (only applied if method is iterative refineR)
#' 
#' @return				(list) with elements RIMat (matrix with RI curves), covarValue (vector with covariate values) and RIperc (vector with percentiles). 
#' 
#' @author Christopher Rank \email{christopher.rank@@roche.com}, Tatjana Ammer \email{tatjana.ammer@@roche.com}

getRICurve <- function(x, covarValue=NULL, RIperc=c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975), Scale=c("original", "transformed"), ...) {
	
	args = list(...)
	stopifnot(class(x) == "RWDRICurve")	
	stopifnot((is.null(covarValue)&x$Params$method!="GAMLSS") | is.numeric(covarValue))
	stopifnot(is.numeric(RIperc) & min(RIperc)>=0 & max(RIperc)<=1)
	Scale    <- match.arg(Scale[1], choices = c("original", "transformed"))
	
	RIperc 	  <- sort(RIperc)
	covarNameGamlss <- x$Params$covarNameGamlss	
	families  <- x$Models$family
	fam 	  <- 1

	assign("families", families, envir = .GlobalEnv)
	assign("fam", fam, envir = .GlobalEnv)
	
	withData <- args$withData
	stopifnot(is.null(withData) | is.logical(withData))
	
	
	if(withData)
		RIMat <- centiles.pred(x$Models, xname = covarNameGamlss, xvalues = covarValue, cent = RIperc*100, data = x$Models$.user$data) 
	else 
		RIMat <- centiles.pred(x$Models, xname = covarNameGamlss, xvalues = covarValue, cent = RIperc*100)
		
		
	RIMat <- t(RIMat[,c(-1)])
			
	rownames(RIMat) <- paste0("Perc_", RIperc)
	colnames(RIMat) <- paste0(covarName, "_", signif(covarValue, 3))		
	
	return(list(RIMat=RIMat, covarValue=covarValue, RIperc=RIperc))	
}




#' Method to calculate confidence intervals
#' 
#' @param pointEst		(object) of class 'RWDRICurve' with point estimate
#' @param estBS			(list) with objects of class 'RWDRICurve' with bootstrap estimates
#' @param covarValue	(numeric) vector of covariate values for which prediction of percentiles shall be performed
#' @param withData		(logical) if GAMLSS percentile estimation needs the input data, default: TRUE
#' @param RIperc		(numeric) value specifying the percentiles, which define the reference interval curves
#' @param CIprop		(numeric) value specifying the central region for estimation of confidence intervals
#' @param onlyCI		(logical) indicating if only confidence intervals should be computed, i.e. if estimation is carried out in a new R session and files have to be read beforehand
#' @param RICurve		(list) with elements RIMat, covarValue, and RIperc for point estimate (output of \code{\link{getRICurve}} function)
#' @param RIMats		(list) with elements RIMat, covarValue, and RIperc for each bootstrap estimate (output of \code{\link{getRICurve}} function)
#' @param ...			additional parameters
#' 
#' @return				(list) with elements RICurve (list with point estimates), RIMats (list with boostrap estimates), CILow (matrix with estimates for lower CI), CIHigh (matrix with estimates for upper CI)
#' 
#' @author Tatjana Ammer \email{tatjana.ammer@@roche.com}

estimateCIs <- function(pointEst, estBS, covarValue = seq(0,7000,1), withData =TRUE, RIperc = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975), CIprop = 0.95, onlyCI = FALSE,
		RICurve = NULL, RIMats = NULL, ...){
	
	args = list(...)
	
	extend = TRUE
	if(!is.null(args$extend))
		extend = args$extends
	
	if(!onlyCI & is.null(RICurve) & is.null(RIMats)){
		
		families <- pointEst$Models$family
		fam 	 <- 1
		invisible(capture.output(RICurve <- getRICurve(pointEst, covarValue = covarValue, withData = withData)))
		
		invisible(capture.output(
						RIMats <- lapply(estBS, function(x){ 
									getRICurve(x, covarValue = covarValue, withData =TRUE)
									
								})
				))
		
	}
	
	RIperc <- sort(RIperc)
	
	CILow = list()
	CIHigh = list()
	
	for(i in 1:length(RIperc)) {
		
		RIBS <- NULL 
		
		for(l in 1:length(RIMats)) {
			RIBS <- rbind(RIBS, RIMats[[l]]$RIMat[i,])
		}
		
		RIBS <- as.data.frame(RIBS)
		
		CILow[[i]] = apply(RIBS, MARGIN = 2, FUN = function(x){ 
					as.numeric(quantile(x,(1-CIprop)/2, na.rm = TRUE))}
		)
		
		CIHigh[[i]] = apply(RIBS, MARGIN = 2, FUN = function(x){ 
					as.numeric(quantile(x,1-(1-CIprop)/2, na.rm = TRUE))}
		)
		
		if(extend){
			# extend CI to include point estimate
			CILow[[i]][RICurve$RIMat[i,] < CILow[[i]]]  <- RICurve$RIMat[i,][RICurve$RIMat[i,] < CILow[[i]]]
			CIHigh[[i]][RICurve$RIMat[i,] > CIHigh[[i]]]  <- RICurve$RIMat[i,][RICurve$RIMat[i,] > CIHigh[[i]]]
		}
			
	}
	
	return(list(RICurve = RICurve, RIMats = RIMats, CILow = CILow, CIHigh = CIHigh))
}


#' Util function to apply a gaussian smoothing function to a vector 
#' @params RIVec		(numeric) vector specifying the vector to be smoothed
#' @params spanX		(numeric) indicating the smoothing strength
#' 
#' @return (numeric) smoothed vector 
#' 
#' @author Christopher Rank \email{christopher.rank@@roche.com}, Tatjana Ammer \email{tatjana.ammer@@roche.com}
smoothVec <- function(RIVec,  spanX=0.3) {
	
	stopifnot(is.numeric(spanX) & ((spanX>0.0 & spanX<1.0) | spanX>1.0))
	
	if(spanX < 1.0) {
		
		xShift <- max(1, ceiling(spanX*length(RIVec)/2))		
		sdDist <- spanX*length(RIVec)/6
		spanXPrior <- spanX/15
		
	} else {
		
		xShift <- max(1, ceiling(spanX/2))
		sdDist <- spanX/6
		spanXPrior <- spanX/ncol(RIMat)/15	
		
	}	
	
	# resize x	
	RIVecLarge <- vector(mode ='numeric', length = length(RIVec) + 2*xShift)
	
	RIVecLarge[       						1:          xShift] 	<- RIVec[1]
	RIVecLarge[   (length(RIVecLarge)-xShift+1):length(RIVecLarge)] <- RIVec[length(RIVec)]
	RIVecLarge[   (1:length(RIVec))+xShift] 						<- RIVec
	
	# apply gaussian filtering along age direction
	for(icol in 1:length(RIVec)) {
		
		# extract sliding window
		temp <- RIVecLarge[icol:(icol+2*xShift)]
		
		weightsDist <- dnorm(x=(-xShift):(xShift), mean=0, sd=sdDist)
		weightsConc <- 1
		
		weights <- weightsDist*weightsConc						
		weights[is.na(temp)] <- 0
		
		# calculate weighted average
		tempSmooth <- sum(temp*weights/sum(weights), na.rm=TRUE)
		
		RIVec[ icol] <- tempSmooth			
	}	
	
	return(RIVec)	
}



str2lang <- function(s){
	parse(text = s, keep.source=FALSE)[[1]]
} 



## util function for transforming x value before applying gamlss 
findPower <- function(y, x, data = NULL,  lim.trans = c(0, 1.5), prof=FALSE, k=2,  c.crit = 0.01, step=0.1)  
{
	cat("*** Checking for transformation for x ***", "\n") 
	ptrans	<- function(x, p) if (abs(p)<=0.0001) log(x) else I(x^p)
	fn 		<- function(p) GAIC(gamlss(y~pb(ptrans(x,p)), c.crit = c.crit, trace=FALSE), k=k)
	
	if (prof) # profile dev
	{
		pp <- seq(lim.trans[1],lim.trans[2], step) 
		pdev <- rep(0, length(pp)) 
		for (i in 1:length(pp)) 
		{
			pdev[i] <- fn(pp[i])  
			#   cat(pp[i], pdev[i], "\n")
		}
		plot(pdev~pp, type="l")
		points(pdev~pp,col="blue")
		par <- pp[which.min(pdev)]
		cat('*** power parameters ', par,"***"," \n") 
	} else
	{
		fn 	<- function(p) GAIC(gamlss(y~pb(ptrans(x,p)), c.crit = c.crit, trace=FALSE), k=k)
		par <- optimise(fn, lower=lim.trans[1], upper=lim.trans[2])$minimum
		cat('*** power parameters ', par,"***"," \n") 
	}  
	return(par)
}

# util function for applying a transformation to x values
ptrans<- function(x, p) if (p==0) log(x) else I(x^p)

