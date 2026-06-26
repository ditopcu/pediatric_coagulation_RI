#' Plot method for reference interval (percentile) curves with confidence intervals

#' @param allRes			(list) with elements RICurve (list with point estimates), RIMats (list with boostrap estimates), CILow (matrix with estimates for lower CI), CIHigh (matrix with estimates for upper CI)
#' @param withCI			(logical) specifying if confidence intervals should be plotted
#' @param RIperc			(numeric) value specifying the percentiles, which define the reference interval curves
#' @param xlim				(numeric) vector specifying the limits in x-direction	
#' @param ylim				(numeric) vector specifying the limits in y-direction	
#' @param xlab				(character) specifying the x-axis label	
#' @param ylab				(character) specifying the y-axis label	
#' @param title				(character) specifying plot title
#' @param cols				(character) specifying the colors used for the percentiles curves
#' @param ltys				(numeric) specifying the line types 
#' 
#' @return				No return value. Instead, a plot is generated.
#' 
#' @author Tatjana Ammer \email{tatjana.ammer@@roche.com}

plotRICurveEstimates <- function(allRes = NULL,  withCI =TRUE, RIperc =c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975), 
		xlim = NULL,  ylim = NULL, xlab ="Age", ylab = "Concentration [Units]",
		title ="", cols = NULL, ltys = c(4,3,2,1,2,3,4)){
	
	
	if(is.null(xlim))
		xlim <- c(min(allRes$RICurve$covarValue), max(allRes$RICurve$covarValue))
	
	if(is.null(ylim))
		ylim <- c(min(allRes$RICurve$RIMat), max(allRes$RICurve$RIMat))
	
	
	if(is.null(cols) & all(RIperc == c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)))
		cols <- c("2.5%" = "red", "10%" = "orange", "25%" = "green3", "50%" = "blue", "75%" = "green3", "90%" = "orange", "97.5%" = "red")
	else if(is.null(cols))
		cols <- rep("black",length(RIperc))				
		
	
	if(is.null(ltys))
		ltys = rep(1, nrow(allRes$RICurve$RIMat))
	
	
	
	title =""
	
	plot(1, type ="n", xaxt ="n", yaxt ="n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, font.lab = 1, cex.lab = 1.2, bty = "l")
	
	refineR:::addGrid(pretty(xlim), pretty(ylim), col = "grey90")
	axis(side = 1, at = pretty(xlim), labels =pretty(xlim), lwd.ticks = 1)
	axis(side = 2, at = pretty(ylim), las = 1)
	box()
	
	mtext(title, side = 3, line = 1.5, cex = 1.3, font = 2)

	for(i in 1:nrow(allRes$RICurve$RIMat)){
		lines(allRes$RICurve$covarValue, allRes$RICurve$RIMat[i,], col = cols[i], lwd = 2)
		if(withCI){
			polygon(x = c(allRes$RICurve$covarValue, rev(allRes$RICurve$covarValue)), y = c(allRes$CILow[[i]], rev(allRes$CIHigh[[i]])), 
				col = refineR:::as.rgb(cols[i],0.3), border =NA)
		}
	
	}
	
}


#' Plot method for reference interval (percentile) curves with confidence intervals with x-axis split like shown in Figure 2.

#' @param allRes			(list) with elements RICurve (list with point estimates), RIMats (list with boostrap estimates), CILow (matrix with estimates for lower CI), CIHigh (matrix with estimates for upper CI)
#' @param withCI			(logical) specifying if confidence intervals should be plotted
#' @param RIperc			(numeric) value specifying the percentiles, which define the reference interval curves
#' @param xlim				(numeric) vector specifying the limits in x-direction	
#' @param ylim				(numeric) vector specifying the limits in y-direction	
#' @param xlab				(character) specifying the x-axis label	
#' @param ylab				(character) specifying the y-axis label	
#' @param title				(character) specifying plot title
#' @param cols				(character) specifying the colors used for the percentiles curves
#' @param ltys				(numeric) specifying the line types 
#' 
#' @return				No return value. Instead, a plot is generated.
#' 
#' @author Tatjana Ammer \email{tatjana.ammer@@roche.com}
plotRICurveEstimatesSplitAxis <- function(allRes = NULL, withCI =TRUE,  RIperc =c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975), 
		xlim = NULL,  ylim = NULL, xlab ="Age", ylab = "Concentration [Units]",
		title ="", cols = NULL, ltys = NULL){
	
		
	if(is.null(ylim))
		ylim <- c(min(allRes$RICurve$RIMat), max(allRes$RICurve$RIMat))
	
	if(is.null(cols) & all(RIperc == c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)))
		cols <- c("2.5%" = "red", "10%" = "orange", "25%" = "green3", "50%" = "blue", "75%" = "green3", "90%" = "orange", "97.5%" = "red")
	else if(is.null(cols))
		cols <- rep("black",length(RIperc))				
	
	
	if(is.null(ltys))
		ltys = c(4,3,2,1,2,3,4)

	covarValue = allRes$RICurve$covarValue
	
	dayScale <- c(0,30)
	monthScale <- c(1,12)
	yearScale <- c(1,18)
	
	# compute index of subsets for day/month/year plot 
	daysIndex 	<- which(covarValue >= dayScale[1] & covarValue <= dayScale[2])
	monthsIndex <- which(covarValue > monthScale[1]*30.4375 & covarValue <= monthScale[2]*30.4375)
	yearsIndex 	<- which(covarValue > yearScale[1]*365.25 & covarValue <= yearScale[2]*365.25)
	
	
	RICurvePointEst <- allRes$RICurvePointEst
	RIMats 			<- allRes$RIMats
	CILow 			<- allRes$CILow
	CIHigh 			<- allRes$CIHigh
	
	# plot point estimate with confidence intervals 
	layout(mat = matrix(c(1,2,3),nrow = 1, ncol = 3), widths = c(1.4,1,2))
	plotSubplot(RICurvePointEst = allRes$RICurve, RIMatsCI = allRes$RIMats, index = daysIndex, convFactor = 1, cols = cols, xlim = dayScale, xlab ="Age [Days]", ylim = ylim, 
			ylab =ylab, withCI = withCI, CILow = CILow, CIHigh = CIHigh, ltys = ltys, xticks = c(0, 10, 20, 30))
	plotSubplot(RICurvePointEst = allRes$RICurvePointEst, RIMatsCI = allRes$RIMats, index = monthsIndex, convFactor = 30.4375, cols = cols, xlim = monthScale, xlab ="Age [Months]", ylim = ylim, 
			withCI =TRUE, CILow = allRes$CILow, CIHigh =  allRes$CIHigh, ltys = ltys, xticks = c(1, 3, 6, 9, 12 ))
	plotSubplot(RICurvePointEst = allRes$RICurvePointEst, RIMatsCI = allRes$RIMats, index = yearsIndex, convFactor = 365.25, cols = cols, xlim = yearScale, xlab ="Age [Years]", ylim = ylim,
			withCI =TRUE, CILow =  allRes$CILow, CIHigh =  allRes$CIHigh, ltys = ltys, xticks = c(1,3,6, 9, 12, 15, 18))
	
}




#' Helper function to produce sub plots for plots with x-axis breaks with confidence intervals 
#' 
#' @param RICurvePointEst	(list) with elements RIMat (matrix with RI curves), covarValue (vector with covariate values) and RIperc (vector with percentiles).
#' @param RIMatsCI	
#' @param index				(numeric) vector specifying the indices part of the current subplot 
#' @param convFactor		(numeric) specifying the conversion factor applied to convert age from days to months/years
#' @param cols				(character) specifying the colors used for the percentiles curves
#' @param xlim				(numeric) vector specifying the limits in x-direction	
#' @param xlab				(character) specifying the x-axis label	
#' @param ylim				(numeric) vector specifying the limits in y-direction	
#' @param ylab				(character) specifying the y-axis label
#' @param withCI 			(logical) specifying if confidence intervals should be plotted
#' @param CILow				
#' @param CIHigh 
#' @param ltys 				(numeric) specifying the line types 
#' @param percIndex
#' @param xticks 	
#' 
#' @return				No return value. Instead, a plot is generated.
#' 
#' @author Tatjana Ammer \email{tatjana.ammer@@roche.com}
plotSubplot <- function(RICurvePointEst = NULL,RIMatsCI = NULL, index, convFactor, cols, xlim, xlab, ylim, ylab = NULL, withCI =TRUE, CILow = NULL, CIHigh = NULL, ltys = c(4,3,2,1,2,3,4), 
		percIndex = NULL, xticks = NULL){
	
	if(!is.null(ylab)){
		par(mar = c(5.1, 6.2, 1, 0.8))
		plot(1, type ="n",xlab = xlab, ylab =paste0(ylab,"\n"), ylim = ylim, xlim = xlim,
				yaxt ="n", xaxt ="n", font.lab = 1, cex.lab = 1.7, bty ="l")
		axis(side = 2, at = pretty(ylim), las = 1, cex.axis = 1.7)
		
	}else {
		par(mar = c(5.1, 1, 1, 0.8))
		plot(1, type ="n", xlab = xlab, ylab = "", ylim = ylim, xlim = xlim,
				yaxt ="n", xaxt ="n", font.lab = 1, cex.lab = 1.7, bty = "n")
	}
	
	
	if(is.null(xticks)){
		xticks <- pretty(xlim)
		
		if(min(xticks) < xlim[1])
			xticks <- xticks[-1]
		
		if(max(xticks) > xlim[2])
			xticks <- xticks[-length(xticks)]
	}
	
	
	if(is.null(ltys)){
		ltys = rep(1, nrow(RICurvePointEst$RIMat))
	}
	
	
	refineR:::addGrid(x = xticks, y = pretty(ylim))
	axis(side = 1, at = xlim, labels =c("",""), lwd.ticks = 0)
	axis(side = 1, at = xticks, las = 1, cex.axis = 1.7)
	
	
	if(withCI & !is.null(CILow) & !is.null(CIHigh)){
		if(is.null(RICurvePointEst))
			RICurvePointEst = RIMatsCI[[1]]
		
		if(!is.null(percIndex))
			iterInd <- percIndex
		else 
			iterInd <- 1:nrow(RICurvePointEst$RIMat)
		
		for(i in iterInd){
			lines(RICurvePointEst$covarValue[index]/convFactor, RICurvePointEst$RIMat[i,index], col = cols[i], lwd = 2, lty = ltys[i])
			polygon(x = c(RICurvePointEst$covarValue[index]/convFactor, rev(RICurvePointEst$covarValue[index]/convFactor)),
					y = c(CILow[[i]][index], rev( CIHigh[[i]][index])), 
					col = refineR:::as.rgb(cols[i],0.3), border =NA)
		}
	}else{
		
		for(i in 1:nrow(RICurvePointEst$RIMat)){
			lines(RICurvePointEst$covarValue[index]/convFactor, RICurvePointEst$RIMat[i,index], col = cols[i], lwd = 2, lty = ltys[i])
		}
	}
}

