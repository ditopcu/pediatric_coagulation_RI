#============================================================================================
# Methods for GAMLSS Approach
#============================================================================================

#' Function to run the whole pipeline for the point estimate
#' 
#' @param inputData				(data.frame) with columns containing the covariate, value and an identifier, column names should correspond to 'covarName', 'colnameID', 'colnameValue'
#' @param covarName				(character) specifying the name of column containing the covariate
#' @param colnameID				(character) specifying the name of column containing the unique identifier
#' @param colnameValue			(character) specifying the name of column containing the value
#' @param outputDir				(character) specifying the directory where the results should be written to
#' @param filename				(character) specifying the root filename
#' @param NCores				(integer) specifying the number of cores used for parallelization, if NULL, number of available cores are used
#' @param ... 					additional parameters
#' 
#' @return	(object) of class 'RWDRICurve' containing the point estimate 
#' 
#' @author Tatjana Ammer \email{tatjana.ammer@@roche.com}

runPipelinePointEst <- function(inputData = NULL, covarName ="Age", colnameID = "PID", colnameValue ="Value", outputDir = getwd(), filename = NULL, NCores = NULL, ...){
	
	args = list(...)
	
	plan(multisession, workers = max(1, min(detectCores(logical = TRUE), NCores)))
	
	# check of input parameters 
	stopifnot(!is.null(inputData))
	
		
	if(is.null(filename))
		filename <- "Pipeline"
	
	# set all default parameters 
	families <- args$families
	if(is.null(families))
		families <- c("BCCG", "BCCGo", "BCT", "BCTo", "BCPEo", "LOGNO")
	
	tol <- args$tol 
	if(is.null(tol))
		tol <- 0.01
	
	minN <- args$minN
	if(is.null(minN))
		minN  <- 1000
	
	unique <- args$unique
	if(is.null(unique))
		unique <- TRUE
	
	Data <- inputData
	
	if(file.exists(file.path(outputDir, paste0(filename, "_refineR_Est_PointEst", ".RData")))){
		res <- get(load(file = file.path(outputDir, paste0(filename, "_refineR_Est_PointEst", ".RData"))))
		
	}else {
		
		aaa <- defineAgeGroupsWithTol(Data = Data, covarName = covarName, colnameID = colnameID, colnameValue = colnameValue, minN = minN, unique = unique, tol = tol)
		res <- runRefineRUniqueID(x = aaa$ageGroupsRIs, colnameValue = colnameValue, colnameID = colnameID)
		
		save(res, file = file.path(outputDir, paste0(filename, "_refineR_Est_PointEst", ".RData")))
	}
	
	
	if(file.exists(file.path(outputDir, paste0(filename, "_gamlss_PointEst", ".RData")))){
		
		pointEstGamlss<- get(load(file = file.path(outputDir, paste0(filename, "_gamlss_PointEst",  ".RData"))))
		
	}else{
		pointEstGamlss <- estimateModelGAMLSS(x = res$algoResults, Data= Data, colnameValue = colnameValue ,  
				families =families, pp = NULL, checkAllCombinations = TRUE)
		
		save(pointEstGamlss, file = file.path(outputDir, paste0(filename, "_gamlss_PointEst",   ".RData")))
	}
	
	
	plan(sequential) 
	
	return(pointEstGamlss)
}


#' Function to bootstrap and run the whole pipeline
#' 
#' @param inputData				(data.frame) with columns containing the covariate, value and an identifier, column names should correspond to 'covarName', 'colnameID', 'colnameValue'
#' @param pointEstGamlss		(object) of class 'RWDRICurve' containing the previously estimated point estimate
#' @param NBootstrap			(integer) specifying the number of bootstrap repetitions
#' @param covarName				(character) specifying the name of column containing the covariate
#' @param colnameID				(character) specifying the name of column containing the unique identifier
#' @param colnameValue			(character) specifying the name of column containing the value
#' @param outputDir				(character) specifying the directory where the results should be written to
#' @param filename				(character) specifying the root filename
#' @param NCores				(integer) specifying the number of cores used for parallelization, if NULL, minimum of number of available cores and NBoostrap are used
#' @param ... 					additional parameters
#' 
#' @return	(list) with results from bootstrapping the pipeline
#' 
#' @author Tatjana Ammer \email{tatjana.ammer@@roche.com}

runPipelineCIs <- function(inputData = NULL, pointEstGamlss = NULL, NBootstrap = 5, covarName ="Age", colnameID = "PID", colnameValue ="Value", outputDir = getwd(), filename =NULL, NCores = NULL,
		 ...){
	 
	 args = list(...)
	 
	 # check of input parameters 
	 stopifnot(!is.null(inputData))
	 stopifnot(!is.null(pointEstGamlss))
	 
	 if(is.null(filename))
		 filename <- "Pipeline"
	 
	
	 # check if all files are on disk, if not induce only sequential mode to load pre-computed datasets
	 filesOnDisk <- NULL
	 
	 for(i in 1:NBootstrap) {
		 
		 filesOnDisk <- c(filesOnDisk, file.exists(file.path(outputDir, paste0(filename, "_refineR_Est_",i, ".RData"))))
		 filesOnDisk <- c(filesOnDisk, file.exists(file.path(outputDir, paste0(filename, "_gamlssModel_Est_",i, ".RData"))))
	 }
	 
	 if(!all(filesOnDisk)) {
		 plan(multisession, workers = max(1, min(NBootstrap, detectCores(logical = TRUE), NCores)))
		 
	 } else {
		 plan(sequential) 
	 }
	
	 
	# set all default parameters 
	tol <- args$tol 
	if(is.null(tol))
		tol <- 0.01
	
	minN <- args$minN
	if(is.null(minN))
		minN  <- 1000
	
	unique <- args$unique
	if(is.null(unique))
		unique <- TRUE
	

	seed <- 123
	
	resBS <- NULL 
	
	finalFamily <- pointEstGamlss$Models$family[1]
	pp  		<- pointEstGamlss$Params$transformP
	
	if(length(pointEstGamlss$Models$sigma.coefficients) <= 1 & length(pointEstGamlss$Models$nu.coefficients) <= 1 & length(pointEstGamlss$Models$tau.coefficients) <= 1){
		depParamsOrig <- 1
	}else if(length(pointEstGamlss$Models$sigma.coefficients) == 2 & length(pointEstGamlss$Models$nu.coefficients) <= 1 & length(pointEstGamlss$Models$tau.coefficients) <= 1){
		depParamsOrig <- 2
	}else if(length(pointEstGamlss$Models$sigma.coefficients) == 2 & length(pointEstGamlss$Models$nu.coefficients) == 2 & length(pointEstGamlss$Models$tau.coefficients) <= 1){
		depParamsOrig <- 3
	}else if(length(pointEstGamlss$Models$sigma.coefficients) == 2 & length(pointEstGamlss$Models$nu.coefficients) == 2 & length(pointEstGamlss$Models$tau.coefficients) == 2){
		depParamsOrig <- 4
	}
	
	
	iterations <- 1:(NBootstrap)
		
	resBS <- future_lapply(iterations, FUN = function(i){
				
				
				Data <- inputData[sample(nrow(inputData), size = nrow(inputData), replace =TRUE),]		
				
				if(file.exists(file.path(outputDir, paste0(filename, "_refineR_Est_",i, ".RData")))){
					
					res <- get(load(file = file.path(outputDir, paste0(filename, "_refineR_Est_",i, ".RData"))))
					
				}else {
					aaa <- defineAgeGroupsWithTolBS(Data = Data, covarName = covarName, colnameID = colnameID, colnameValue = colnameValue, minN = minN, unique = unique, tol = tol)
					res <- runRefineRUniqueID(x = aaa$ageGroupsRIs, colnameValue = colnameValue, colnameID = colnameID)
					
					save(res, file = file.path(outputDir, paste0(filename, "_refineR_Est_",i, ".RData")))
					
				}
				
				
				if(file.exists(file.path(outputDir, paste0(filename, "_gamlssModel_Est_",i,  ".RData")))){
					
					gamlssModel <- get(load(file = file.path(outputDir, paste0(filename, "_gamlssModel_Est_",i,  ".RData"))))
					
				}else{
					
					gamlssModel <- tryCatch({
								estimateModelGAMLSS(x = res$algoResults, Data = Data, colnameValue = colnameValue ,
										families = finalFamily, pp = pp, depOrig = depParamsOrig)
							}, 
							error = function(e){
								message("gamlss failed")
								message(e)
								obj <- list()
								class(obj) <- "RWDRICurve"
								obj$.user <- list(
										data = Data,
										call = match.call()
								)
								obj$Params <- list(pp = pp)
								obj$Params <- list(method = "GAMLSS", covarName = res$algoResults$Params$covarName, transformP = pp, family = finalFamily)
								obj$Models <- NULL
								obj$refineRModels <- res$algoResults$Models		
								return(obj)
								
							}
					
					)
					save(gamlssModel, file = file.path(outputDir, paste0(filename, "_gamlssModel_Est_",i,  ".RData")))
										
				}
				return(gamlssModel)
			}
			,
			future.packages=c("refineR", "gamlss"), 
			future.seed=future_lapply(1:length(iterations), FUN = function(x) .Random.seed, future.chunk.size = Inf, future.seed = seed)
	)
	
	
	plan(sequential) 
	
	return(resBS)
}



#' Method to run refineR  for objects of class 'RWDRICurve'
#' 
#' @param x				(object) of class 'RWDRICurve'
#' @param colnameValue	(character) specifying the column name for which values the estimation of reference intervals should be performed
#' @param colNameID		(character) specifying the name of the covariate that is used in the filter out multiple samples from one patient
#' @param RIperc		(numeric) value specifying the percentiles, which define the reference interval curves
#' @param ...			NCores 
#' 
#' @return				(list) with elements RIMat (matrix with RI curves), algoResults (RWDRI Curve object with estimated refineR models). 
#' 
#' @author Tatjana Ammer \email{tatjana.ammer@@roche.com}

runRefineRUniqueID <- function(x = NULL, colnameValue ="Value", colnameID = "PID",  RIperc = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975), model = "BoxCox",seed = 123, ...){
	
	args = list(...)
	
	RIMat <- NULL
	
	indModels <- x$Models
	
	# run refineR in parallel 
	indModels <- future_lapply(indModels, function(m){
				findRI(m$Data, model = model)
			},future.packages=c("refineR"), future.seed=future_lapply(1:length(indModels), FUN = function(x) .Random.seed, future.chunk.size = Inf, future.seed = seed))
	

	for (i in 1:length(x$Models)){
		
		x$Models[[i]]$N				<- length(indModels[[i]]$Data)
		x$Models[[i]]$Model			<- indModels[[i]]$Model
		
		x$Models[[i]]$Lambda 		<- indModels[[i]]$Lambda 
		x$Models[[i]]$Mu 			<- indModels[[i]]$Mu 
		x$Models[[i]]$Sigma 		<- indModels[[i]]$Sigma
		x$Models[[i]]$P 			<- indModels[[i]]$P 
		x$Models[[i]]$Shift 		<- indModels[[i]]$Shift 
		x$Models[[i]]$Cost 			<- indModels[[i]]$Cost
		x$Models[[i]]$abOr 			<- indModels[[i]]$abOr
		x$Models[[i]]$roundingBase 	<- indModels[[i]]$roundingBase
		x$Models[[i]]$rf 			<- indModels[[i]]$rf
		
		
		ri			<- getRI(x$Models[[i]], RIperc = RIperc)
		
		
		tmp <- rbind(c(i, x$Models[[i]]$Age, ri$PointEst, x$Models[[i]]$N))

		RIMat <- rbind(RIMat, tmp)
		
	}
	
	return(list(RIMat = RIMat, algoResults = x))
}



#' Method to define age groups with initial width increasing linearly with age
#' 
#' @param Data			(data.frame) with columns indicating the covariate, the values used for estimation, and an identifier used for filtering multiple samples per subject
#' @param covarName 	(character) specifying the column name of the covariate, default "Age"
#' @param colNameID		(character) specifying the name of the column that is used in the filter out multiple samples from one patient, default "PID"
#' @param colnameValue	(character) specifying the column name for which values the estimation of reference intervals should be performed, default "Value"
#' @param minN 			(integer) specifying the minimum number of samples that needs to be reached for each age group, i.e. the groups are filled with 
#' 							samples from left and right of the central age(s) until minN is reached
#' @param unique		(logical) indicating if each subject should only be represented once in each age group
#' @param tol 			(numeric) specifying the rate in which the width increases depending on the covariate
#' 
#' @return				(list) with elements ageGroupsRIs (RWDRI Curve object with defined age groups), ageGroupsInfo (matrix with info about defined age groups). 
#' 
#' @author Tatjana Ammer \email{tatjana.ammer@@roche.com}

defineAgeGroupsWithTol <- function(Data = NULL, covarName = "Age", colnameID ="PID", colnameValue = "Value", minN = 1000, unique =TRUE, tol = 0.01){
	
	x <- list()
	class(x) <- "RWDRICurve"
	
	x$Params <- list(method="GAMLSS_refineR", covarName = covarName)	
	
	minAge <- min(Data[,covarName])
	maxAge <- max(Data[,covarName])
	
	# compute size of age groups with tolerance	
	ageGroups <- data.frame(ind = 1, min = 0, max = 0)
	i <- 1
	
	while (ageGroups$max[i] < maxAge) {
		
		minA <- ageGroups$max[i] + 1 
		maxA <- minA 
		
		while(((maxA - minA + 1)/(minA+1)) < tol){
			maxA <- maxA + 1
		}
		i <- i+1
		ageGroups[i,] <- list(ind = 1, min = minA,  max = maxA)
	}
	
	
	# enumerate pids for deciding whether to take samples from left or right 
	pids   <- unique(Data[,colnameID]) 
	pidsNr <- data.frame(IDS = pids, pidsNr = 1:length(pids))  
	colnames(pidsNr)[1] <- colnameID 
	Data <- merge(Data, pidsNr, by = colnameID)
	
	models <- list()
	
	ageGroups$lowerBorder <- NA
	ageGroups$upperBorder <- NA
	ageGroups$N 		  <- NA 
	
	# traverse through age groups and fill until minN is reached
	for(r in 1:nrow(ageGroups)){
		
		minA <- ageGroups$min[r]
		maxA <- ageGroups$max[r]
		
		# get data for central age group
		tmp <- Data[Data[,covarName] >= minA & Data[,covarName] <= maxA,] 
		tmp <- tmp[order(tmp[,covarName]),]
		
		#get unique IDs
		if(unique){
			if(nrow(tmp) != 0){
				tmp$Include <- FALSE
				tmp[match(unique(tmp[,colnameID]), tmp[,colnameID]),]$Include <- TRUE
				n   <- nrow(tmp[tmp$Include ==TRUE,])
				ids <- unique(tmp[,colnameID])
			} else {
				n <- 0 
			}
		}else {
			n <- nrow(tmp)
		}
		
		lowerB <- minA 
		upperB <- maxA 
		
		lowerBPrev <- minA
		upperBPrev <- maxA
		
		# expand group until minN is reached
		while(n < minN){
			
			if(lowerB > minAge & upperB < maxAge){
				lowerB <- lowerB - 1 
				upperB <- upperB + 1 
				
				left <- Data[Data[,covarName] >= lowerB & Data[,covarName] < lowerBPrev,]
				left <- left[order(left[,covarName]),]
				
				right <- Data[Data[,covarName] > upperBPrev & Data[,covarName] <= upperB,]
				right <- right[order(right[,covarName]),]
				
				lowerBPrev <- lowerB
				upperBPrev <- upperB
				# get unique IDs
				
				if(unique){
					
					left  <- removeDuplIDs(left, ids, colnameID)		
					right <- removeDuplIDs(right, ids, colnameID)
					
					
					idsLeft  <- left[left$Include ==TRUE, colnameID]
					idsRight <- right[right$Include == TRUE, colnameID] 
					
					# if pids are both in left and right extension, take samples from left if pidsNr are uneven, and from right if pidsNr are even
					if(nrow(left[left$Include ==TRUE & left[,colnameID] %in% idsRight & left$pidsNr%%2 ==0,])> 0)
						left[left$Include ==TRUE & left[,colnameID] %in% idsRight & left$pidsNr%%2 ==0,]$Include <- FALSE
					
					if(nrow(right[right$Include ==TRUE & right[,colnameID] %in% idsLeft & right$pidsNr%%2 ==1,])> 0)
						right[right$Include ==TRUE & right[,colnameID] %in% idsLeft & right$pidsNr%%2 ==1,]$Include <- FALSE
					
					ids <- c(ids, left[left$Include ==TRUE, colnameID], right[right$Include == TRUE,colnameID]) 
					
					tmp <- rbind(left, tmp, right) 
					n   <- nrow(tmp[tmp$Include ==TRUE,])
					
				}else {
					n <- nrow(tmp)
				}
				
			}else if(lowerB > minAge){ # right end
				
				lowerB <- lowerB - 1
				
				left <- Data[Data[,covarName] >= lowerB & Data[,covarName] < lowerBPrev,]
				left <- left[order(left[,covarName]),]
				
				lowerBPrev <- lowerB
				
				# get unique IDs
				if(unique){
					left <- removeDuplIDs(left, ids, colnameID)
					
					ids <- c(ids, left[left$Include == TRUE,colnameID])
					
					tmp <- rbind(tmp, left) 
					n   <- nrow(tmp[tmp$Include ==TRUE,])
					
				}else {
					n <- nrow(tmp)
				}
				
			}else if(upperB < maxAge){ # left end
				upperB <- upperB + 1
				
				right <- Data[Data[,covarName] > upperBPrev & Data[,covarName] <= upperB,]
				right <- right[order(right[,covarName]),]
				
				upperBPrev <- upperB
				
				# get unique IDs
				if(unique){
					right <- removeDuplIDs(right, ids, colnameID)
					
					ids <- c(ids, right[right$Include == TRUE, colnameID]) 
					
					tmp <- rbind(tmp, right) 
					n   <- nrow(tmp[tmp$Include ==TRUE,])
					
				}else {
					n <- nrow(tmp)
				}
			}			
		}
		
		ageGroups[r, c("lowerBorder", "upperBorder", "N")] <- c(lowerB, upperB, n)
		
		# init RWDRI object
		obj 		<- list()
		obj$DataSet	<- tmp
		
		if(unique)
			obj$Data 	<- tmp[tmp$Include ==TRUE, colnameValue]
		else 
			obj$Data	<- tmp[,colnameValue]
		
		
		obj$lowerB   <- lowerB
		obj$upperB	 <- upperB
		obj$startAge <- minA
		obj$stopAge  <- maxA
		obj$Age		 <- mean(minA, maxA)
		obj$Method 	 <- NA 
		obj$Lambda 	 <- NA 
		obj$Mu 		 <- NA
		obj$Sigma 	 <- NA 
		obj$P 		 <- NA
		obj$Cost	 <- NA
		obj$Shift	 <- 0 
		obj$roundingBase <- NA
		
		class(obj) <- "RWDRI"
		
		models[[r]] <-obj
	}
	x$Models <- models
	
	return(list(ageGroupsRIs = x, ageGroupsInfo = ageGroups))
}



#' Method to define age groups with initial width increasing linearly with age for bootstrapping (i.e. keeping identical samples)
#' 
#' @param Data			(data.frame) (Resampled) dataset with columns indicating the covariate, the values used for estimation, and an identifier used for filtering multiple samples per subject
#' @param covarName 	(character) specifying the column name of the covariate, default "Age"
#' @param colNameID		(character) specifying the name of the column that is used in the filter out multiple samples from one patient, default "PID"
#' @param colnameValue	(character) specifying the column name for which values the estimation of reference intervals should be performed, default "Value"
#' @param minN 			(integer) specifying the minimum number of samples that needs to be reached for each age group, i.e. the groups are filled with 
#' 							samples from left and right of the central age(s) until minN is reached
#' @param unique		(logical) indicating if each subject should only be represented once in each age group
#' @param tol 			(numeric) specifying the rate in which the width increases depending on the covariate
#' 
#' @return				(list) with elements ageGroupsRIs (RWDRI Curve object with defined age groups), ageGroupsInfo (matrix with info about defined age groups). 
#' 
#' @author Tatjana Ammer \email{tatjana.ammer@@roche.com}

defineAgeGroupsWithTolBS <- function(Data = NULL, covarName = "Age", colnameID ="PID", colnameValue = "Value", minN = 1000, unique =TRUE, tol = 0.01){
	
	x <- list()
	class(x) <- "RWDRICurve"
	
	x$Params <- list(method="GAMLSS_refineR", covarName = covarName)	
	
	minAge <- min(Data[,covarName])
	maxAge <- max(Data[,covarName])
	
	# compute size of age groups wilt tolerance	
	ageGroups <- data.frame(ind = 1, min = 0, max = 0)
	i <- 1
	
	while (ageGroups$max[i] < maxAge) {
		
		minA <- ageGroups$max[i] + 1 
		maxA <- minA 
		
		while(((maxA - minA + 1)/(minA+1)) < tol){
			maxA <- maxA + 1
		}
		i <- i+1
		ageGroups[i,] <- list(ind = 1, min = minA,  max = maxA)
	}
	
	
	# enumerate pids for deciding whether to take samples from left or right 
	pids   <- unique(Data[,colnameID]) 
#	pids   <- sample(pids)  # shuffle  PIDS?  
	pidsNr <- data.frame(IDS = pids, pidsNr = 1:length(pids))  
	colnames(pidsNr)[1] <- colnameID 
	Data <- merge(Data, pidsNr, by = colnameID)
	
	Data$Key <- paste0(Data$PID, "_", Data$Value, "_", Data$Age)
	
	models <- list()
	
	ageGroups$lowerBorder <- NA
	ageGroups$upperBorder <- NA
	ageGroups$N 		  <- NA 
	
	# traverse through age groups and fill until minN is reached
	for(r in 1:nrow(ageGroups)){
		
		minA <- ageGroups$min[r]
		maxA <- ageGroups$max[r]
		
		# get data for central age group
		tmp <- Data[Data[,covarName] >= minA & Data[,covarName] <= maxA,] 
		tmp <- tmp[order(tmp[,covarName]),]
		
		#get unique IDs
		if(unique){
			if(nrow(tmp) != 0){
				tmp$Include <- FALSE
				tmp[match(unique(tmp[,colnameID]), tmp[,colnameID]),]$Include <- TRUE
				
				combis <- unique(tmp[tmp$Include == TRUE,]$Key)
				
				if(nrow(tmp[tmp$Include == FALSE & tmp$Key %in% combis,]) > 0 )
					tmp[tmp$Include == FALSE & tmp$Key %in% combis,]$Include <- TRUE 
				
				n   <- nrow(tmp[tmp$Include ==TRUE,])
				ids <- unique(tmp[,colnameID])
				
			} else {
				n <- 0 
			}
			
		}else {
			n <- nrow(tmp)
		}
		
		lowerB <- minA 
		upperB <- maxA 
		
		lowerBPrev <- minA
		upperBPrev <- maxA
		
		# expand group until minN is reached
		while(n < minN){
			
			if(lowerB > minAge & upperB < maxAge){
				lowerB <- lowerB - 1 
				upperB <- upperB + 1 
				
				left <- Data[Data[,covarName] >= lowerB & Data[,covarName] < lowerBPrev,]
				left <- left[order(left[,covarName]),]
				
				right <- Data[Data[,covarName] > upperBPrev & Data[,covarName] <= upperB,]
				right <- right[order(right[,covarName]),]
				
				lowerBPrev <- lowerB
				upperBPrev <- upperB
				# get unique IDs
				
				if(unique){
					
					left  <- removeDuplIDs(left, ids, colnameID)		
					right <- removeDuplIDs(right, ids, colnameID)
					
					
					leftCombis <- unique(left[left$Include == TRUE,]$Key)
					
					if(nrow(left[left$Include == FALSE & left$Key %in% leftCombis,])>  0)
						left[left$Include == FALSE & left$Key %in% leftCombis,]$Include <- TRUE 
					
					rightCombis <- unique(right[right$Include == TRUE,]$Key)
					
					if(nrow(right[right$Include == FALSE & right$Key %in% rightCombis,])>  0)
						right[right$Include == FALSE & right$Key %in% rightCombis,]$Include <- TRUE 
					
					
					idsLeft  <- left[left$Include ==TRUE, colnameID]
					idsRight <- right[right$Include == TRUE, colnameID] 
					
					# if pids are both in left and right extension, take samples from left if pidsNr are uneven, and from right if pidsNr are even
					if(nrow(left[left$Include ==TRUE & left[,colnameID] %in% idsRight & left$pidsNr%%2 ==0,])> 0)
						left[left$Include ==TRUE & left[,colnameID] %in% idsRight & left$pidsNr%%2 ==0,]$Include <- FALSE
					
					if(nrow(right[right$Include ==TRUE & right[,colnameID] %in% idsLeft & right$pidsNr%%2 ==1,])> 0)
						right[right$Include ==TRUE & right[,colnameID] %in% idsLeft & right$pidsNr%%2 ==1,]$Include <- FALSE
					
					ids <- c(ids, left[left$Include ==TRUE, colnameID], right[right$Include == TRUE,colnameID]) 
					
					tmp <- rbind(left, tmp, right) 
					n   <- nrow(tmp[tmp$Include ==TRUE,])
					
				}else {
					n <- nrow(tmp)
				}
				
			}else if(lowerB > minAge){ # right end
				
				lowerB <- lowerB - 1
				
				left <- Data[Data[,covarName] >= lowerB & Data[,covarName] < lowerBPrev,]
				left <- left[order(left[,covarName]),]
				
				lowerBPrev <- lowerB
				
				# get unique IDs
				if(unique){
					left <- removeDuplIDs(left, ids, colnameID)
					
					leftCombis <- unique(left[left$Include == TRUE,]$Key)
					
					if(nrow(left[left$Include == FALSE & left$Key %in% leftCombis,])>  0)
						left[left$Include == FALSE & left$Key %in% leftCombis,]$Include <- TRUE 
					
					ids <- c(ids, left[left$Include == TRUE,colnameID])
					
					tmp <- rbind(tmp, left) 
					n   <- nrow(tmp[tmp$Include ==TRUE,])
					
				}else {
					n <- nrow(tmp)
				}
				
			}else if(upperB < maxAge){ # left end
				upperB <- upperB + 1
				
				right <- Data[Data[,covarName] > upperBPrev & Data[,covarName] <= upperB,]
				right <- right[order(right[,covarName]),]
				
				upperBPrev <- upperB
				
				# get unique IDs
				if(unique){
					right <- removeDuplIDs(right, ids, colnameID)
					
					rightCombis <- unique(right[right$Include == TRUE,]$Key)
					
					if(nrow(right[right$Include == FALSE & right$Key %in% rightCombis,])>  0)
						right[right$Include == FALSE & right$Key %in% rightCombis,]$Include <- TRUE 
					
					ids <- c(ids, right[right$Include == TRUE, colnameID]) 
					
					tmp <- rbind(tmp, right) 
					n   <- nrow(tmp[tmp$Include ==TRUE,])
					
				}else {
					n <- nrow(tmp)
				}
			}			
		}
		
		ageGroups[r, c("lowerBorder", "upperBorder", "N")] <- c(lowerB, upperB, n)
		
		# init RWDRI object
		obj 		<- list()
		obj$DataSet	<- tmp
		
		if(unique)
			obj$Data 	<- tmp[tmp$Include ==TRUE, colnameValue]
		else 
			obj$Data	<- tmp[,colnameValue]
		
		
		obj$lowerB   <- lowerB
		obj$upperB	 <- upperB
		obj$startAge <- minA
		obj$stopAge  <- maxA
		obj$Age		 <- mean(minA, maxA)
		obj$Method 	 <- NA 
		obj$Lambda 	 <- NA 
		obj$Mu 		 <- NA
		obj$Sigma 	 <- NA 
		obj$P 		 <- NA
		obj$Cost	 <- NA
		obj$Shift	 <- 0 
		obj$roundingBase <- NA
		
		class(obj) <- "RWDRI"
		
		models[[r]] <-obj
	}
	x$Models <- models
	
	return(list(ageGroupsRIs = x, ageGroupsInfo = ageGroups))
}




#' Method to remove multiple samples from a subject for refineR estimation
#' 
#' @param subDf			(data.frame) subset of dataset including additional columns with occurence of ids, and if sample is included or not 
#' @param ids 			(numeric) vector indicating all ids currently in age group 
#' @param colNameID		(character) specifying the name of the column that is used to filter out multiple samples from one patient, default "PID"
#' 
#' @return				(data.frame) with updated column "Include" 
#' 
#' @author Tatjana Ammer \email{tatjana.ammer@@roche.com}

removeDuplIDs <- function(subDf, ids, colnameID = "PID"){
	
	if(nrow(subDf)){
		subDf$Include <- FALSE
		subDf[match(unique(subDf[,colnameID]), subDf[,colnameID]),]$Include <- TRUE
		
		# check if pids already included, if yes, set include to False (i.e. remove for refineR est)
		if(nrow(subDf[subDf$Include ==TRUE & (subDf[,colnameID] %in% ids),]) > 0 )
			subDf[subDf$Include ==TRUE & (subDf[,colnameID] %in% ids),]$Include <- FALSE
	}
	
	return(subDf)
}



#' Method to compute probability of being non-pathological using estimated refineR model
#' 
#' @param x				(object) of class 'RWDRICurve'
#' @param Data			(data.frame) containing the whole dataset
#' @param colNameValue	(character) specifying the column name for which values the estimation of reference intervals should be performed
#' @param colNameID		(character) specifying the name of the column that is used to filter out multiple samples from one subject
#' @param Nhist			(integer) number of bins in the histogram (derived automatically if not set)
#' 
#' @return				(list) with estimated non-pathological probabilities for each group (probsNonPathol) and an updated data.frame (weightData)
#' 							with an additional column indicating the computed probabilites   
#' 
#' @author Tatjana Ammer \email{tatjana.ammer@@roche.com}

computeProbNP <- function(x, Data, colnameValue, colnameID, Nhist = 60){
	
	models 	  <- x$Models
	covarName <- x$Params$covarName
	
	probsNonPathol   <- list()
	
	weightData <- NULL
	
	# compute prob of being non-pathol
	for (irow in 1:length(models)){
		m <- models[[irow]]
		
		start <- m$startAge
		stop  <- m$stopAge	
		
		# compute probability of being non-pathological 
		probNP 		 <- getWeights(model = m, Nhist = Nhist)
		# smooth estimated probabilies along concentration axes
		probNPSmooth <- smoothVec(RIVec = probNP$countsRatio, spanX = 0.15)
						
		# set values to be bound between 0 and 1
		probNPSmooth[probNPSmooth < 0 ] <- 0
		probNPSmooth[probNPSmooth > 1 ] <- 1
			
		probNP$countsRatioSmooth <- probNPSmooth
		
		probsNonPathol[[irow]] <- probNP		

		
		# approx function here and save probNP for each sample in m$DataSet and use that in for loop 		
		subData <- Data[Data[,covarName] >= start & Data[,covarName] <= stop,]
		pNP <- approx(x = probNP$mids, y = probNP$countsRatioSmooth,xout = unique(sort(subData[,colnameValue])), rule = 1)
		
		if(any(is.na(pNP$y)))
			pNP$y[is.na(pNP$y)] <- 0 
		
		
		pNP <- as.data.frame(pNP)
		colnames(pNP) <- c(colnameValue, "ProbNP")
		
		subData <- merge(subData, pNP, by = "Value", all.x = TRUE)
		
		weightData <- rbind(weightData, subData)
	
	}
	
	# adjust probability to occurence of subject in whole dataset
	if(!is.null(colnameID) ){
		weightData$ProbNPFinal <- weightData$ProbNP*(1/weightData$Freq)
		weightData[weightData[,colnameValue] == 0,colnameValue] <- 1e-20
	}
	
	return(list(probsNonPathol = probsNonPathol, Data = weightData))
}




#' Helper function to compute probability of being non-pathological using estimated refineR model
#' 
#' @param model			(object) of class 'RWDRI'
#' @param Nhist			(integer) number of bins in the histogram (derived automatically if not set)
#' 
#' @return				(list) with estimated non-pathological probabilities (for certain concentration regions (mids).  
#' 
#' @author Tatjana Ammer \email{tatjana.ammer@@roche.com}
getWeights <- function(model, Nhist = 60){
	
	# compute reference intervals for defining xlim
	RI 		<- getRI(x = model, RIperc = c(0.005, 0.995))
	
	ab 		<- model$abOr
	xlim 	<- range(c(ab, 0.98*min(RI$PointEst),  1.1*max(RI$PointEst)), na.rm = TRUE)
	
	if(is.na(model$roundingBase))
	{
		# define histogram breaks
		increment  <- diff(xlim)/Nhist	
		breaks1    <- seq(from = xlim[1] - Nhist*increment, to = xlim[2] + Nhist*increment, by = increment)		
		breaks1    <- breaks1[breaks1 > 1e-20]
		
	} else
	{		
		xlimDiff <- diff(xlim)
		binSize <- model$roundingBase*max(1, round(xlimDiff/model$roundingBase/Nhist))
		
		# adapt xlim	
		xlim[1] <- max(0.5*model$roundingBase, round(xlim[1]/ model$roundingBase)*model$roundingBase - 0.5*model$roundingBase)		
		xlim[2] <- xlim[1] + ceiling(xlimDiff/binSize)*binSize
		
		# define histogram breaks
		breaks1 <- seq(from=xlim[1], to=xlim[2], by=binSize)
	}
	
	# generate histogram of data
	hist1  	   <- hist(model$Data[model$Data >= min(breaks1) & model$Data <= max(breaks1)], breaks = breaks1, plot = FALSE)
	
	# sort vectors in increasing order
	sortIndex  <- 1:length(hist1$mids)
	mids	   <- hist1$mids 
	countsData <- hist1$counts
	
	breakL <- breaks1[1:(length(breaks1)-1)]
	breakR <- breaks1[2:length(breaks1)]			
	
	# Box Cox transformation of histogram breaks and histogram range
	breakL 	  <- suppressWarnings(BoxCox(breakL-model$Shift, model$Lambda))		
	breakR 	  <- suppressWarnings(BoxCox(breakR-model$Shift, model$Lambda))	
	
	maxPred <- NA 
	
	if(!is.na(model$Method) && (model$Method == "refineR" & model$PkgVersion >= "1.6.0")){
		
		pCorr <- BoxCox(c(max(min(x$Data-x$Shift), 1e-20), min(max(x$Data-x$Shift), 1e20)), lambda=x$Lambda)				
		pCorr <- pnorm(q=pCorr, mean=x$Mu, sd=x$Sigma)
		pCorr <- 1/(pCorr[2]-pCorr[1])
		countsPred <- pCorr*length(Data)*x$P*(pnorm(q = breakR, mean = x$Mu, sd = x$Sigma) - pnorm(q = breakL, mean = x$Mu, sd = x$Sigma))			
	
	}else {
			
		# theoretical prediction of bin counts
		countsPred <- length(model$Data)*model$P*(pnorm(q = breakR, mean = model$Mu, sd = model$Sigma) - pnorm(q = breakL, mean = model$Mu, sd = model$Sigma))			
	
	}
	
	countsPred[countsPred < 0] <- 0
	
	countsPred 	<- countsPred[sortIndex]
	maxPred 	<- max(countsPred)
	
	# calculate ratio of prediction and data
	countsRatio <- countsPred/countsData
	countsRatio[is.infinite(countsRatio)] <- 0 # division by 0  
	
	return(data.frame(mids, countsRatio))
}



# Helper function for GAMLSS estimation
str2lang <- function(s){
	parse(text = s, keep.source=FALSE)[[1]]
} 



#' Wrapper function to compute gamlss estimations 
#' 
#' @param x						(object) of class 'RWDRICurve'
#' @param Data					(data.frame) containing the whole dataset
#' @param families				(character) specifying which distribution families should be evaluated
#' @param colNameValue			(character) specifying the column name for which values the estimation of reference intervals should be performed
#' @param colNameID				(character) specifying the name of the column that is used to filter out multiple samples from one subject
#' @param pp					(numeric) specifying pre-computed transformation parameter, if NULL (default), parameter will be estimated
#' @param checkAllCombinations	(logical) specifying if all parameter combinations should be tested
#' @param depOrig				(integer) indicating number of parameters depending on covariate
#' 
#' 
#' @return				(list) with estimated non-pathological probabilities (for certain concentration regions (mids).  
#' 
#' @author Tatjana Ammer \email{tatjana.ammer@@roche.com}

estimateModelGAMLSS <- function(x, Data, families = c("BCCG", "BCCGo"), colnameID = "PID",colnameValue = "Value",  
		pp = NULL, checkAllCombinations =FALSE, depOrig =NULL){

	# get number of each pid 
	pidList <- sort(Data[,colnameID])
	pidList	<- as.data.frame(table(pidList), stringsAsFactors = FALSE)
	colnames(pidList) <- c(colnameID, "Freq")
	
	Data <- merge(Data, pidList, by = colnameID)

	colnames(Data)[colnames(Data) == colnameValue] <- "Value"
	
	
	weightData = NULL 
	# compute weights 
	probNPData <-  computeProbNP(x = x, Data = Data, colnameValue = "Value", colnameID ="PID", Nhist = 60)
		
	weightData <- probNPData$Data
	
	if(nrow(weightData[weightData$Value <= 0,]) > 0)
		weightData[weightData$Value <= 0,]$Value <- 1e-10
	
	
	finalModel <- runGAMLSS(weightData = weightData, covarName =x$Params$covarName, families = families, pp = pp, checkAllCombinations = checkAllCombinations, depOrig = depOrig)
	
	
	tt <- list()
	class(tt) <- "RWDRICurve"
	tt$Params <- list(method="GAMLSS", covarNameGamlss = "Age", transformP = finalModel$pp, family = finalModel$finalM$family[1])
	tt$Models <- finalModel$finalM
	
	tt$refineRModels <- x$Models
		
	return(tt)
	
}




#' Method to compute gamlss estimations.
#' 
#' @param weightData			(data.frame) containing the dataset with additional column specifying the computed non-pathological probabilities/weights 
#' @param covarName 			(character) specifying the column name of the covariate, default "Age"
#' @param colnameProbNP			(character) specifying the column name containing the computed probabilities of being non-pathological
#' @param families				(character) specifying which distribution families should be evaluated
#' @param pp					(numeric) specifying pre-computed transformation parameter, if NULL (default), parameter will be estimated
#' @param checkAllCombinations	(logical) specifying if all parameter combinations should be tested
#' @param depOrig				(integer) indicating number of parameters depending on covariate
#' 
#' @return				(list) with final model (finalM) and estimate transformation parameter (pp)  
#' 
#' @author Tatjana Ammer \email{tatjana.ammer@@roche.com}
runGAMLSS <- function(weightData, covarName, colnameProbNP = "ProbNPFinal", families = c( "BCCG", "BCCGo"),  pp = NULL, checkAllCombinations =FALSE, depOrig = NULL){
	
	if(is.null(pp)){
		ppData  <- weightData[weightData[,colnameProbNP] > 0.5,]
		pp 		<- findPower(y = ppData$Value, x = ppData[, covarName], data = ppData)
	}
	
	weightData$CovPP <- ptrans(weightData[,covarName], pp)
	
	
	weightData$Age <- weightData[,covarName]
	
	if(length(families) > 1){
		gamlssModels <- list()
		if(!checkAllCombinations){
			gamlssModels <- list()
			for(fam in 1:length(families)){
				gamlssModels[[fam]] <- tryCatch({
							gamlss(formula = Value ~pb(ptrans(Age, pp)), sigma.formula = ~pb(ptrans(Age, pp)), 
									nu.formula = ~pb(ptrans(Age, pp)), family = families[fam], data = weightData, weights = weightData[,colnameProbNP], na.rm =TRUE)
							
						},
						error = function(cond){
							message(paste("GAMLSS failed for family: ", families[fam]))
							message(cond)
							obj     <- list()
							obj$sbc <- Inf  
							return(obj)
						}
				)
			}	
			
			scores <- lapply(gamlssModels, function(x) x$sbc)
			
			finalM <- gamlssModels[[which.min(scores)]]

			
		}else {
			gamlssModels <- list()
			counter <- 1
			for(fam in 1:length(families)){
				if(families[fam] == "BCCG" | families[fam] =="BCCGo"){
					for(dep in 1:3){
						if(dep ==1){
							
							gamlssModels[[counter]] <- tryCatch({
										gamlss(formula = Value ~pb(ptrans(Age, pp)), sigma.formula = ~1, 
												nu.formula = ~1, family = families[fam], data = weightData, weights = weightData[,colnameProbNP], na.rm =TRUE)
									}, error = function(cond){
										message(paste("GAMLSS failed for family: ", families[fam]))
										message(cond)
										obj     <- list()
										obj$sbc <- Inf  
										return(obj)
									}
									)
						}else if(dep ==2){
							gamlssModels[[counter]] <- tryCatch({
										gamlss(formula = Value ~pb(ptrans(Age, pp)), sigma.formula = ~pb(ptrans(Age, pp)), 
										nu.formula = ~1, family = families[fam], data = weightData, weights = weightData[,colnameProbNP], na.rm =TRUE)
									}, error = function(cond){
										message(paste("GAMLSS failed for family: ", families[fam]))
										message(cond)
										obj     <- list()
										obj$sbc <- Inf  
										return(obj)
								}
								)
						}else if(dep ==3){
							gamlssModels[[counter]] <- tryCatch({
										gamlss(formula = Value ~pb(ptrans(Age, pp)), sigma.formula = ~pb(ptrans(Age, pp)), 
										nu.formula = ~pb(ptrans(Age, pp)), family = families[fam], data = weightData, weights = weightData[,colnameProbNP], na.rm =TRUE)
									}, error = function(cond){
										message(paste("GAMLSS failed for family: ", families[fam]))
										message(cond)
										obj     <- list()
										obj$sbc <- Inf  
										return(obj)
									}
									)
						}
						counter <- counter +1 
					}
				}else if(families[fam] =="BCT" | families[fam] =="BCTo"){
					for(dep in 1:4){
						if(dep ==1){
							gamlssModels[[counter]] <- tryCatch({
										gamlss(formula = Value ~1, sigma.formula = ~1, 
												nu.formula = ~1, tau.formula = ~1, family = families[fam], data = weightData, weights = weightData[,colnameProbNP], na.rm =TRUE)
									}, error = function(cond){
										message(paste("GAMLSS failed for family: ", families[fam]))
										message(cond)
										obj     <- list()
										obj$sbc <- Inf  
										return(obj)
									}
							)
							
						}else if(dep ==2){
							gamlssModels[[counter]] <- tryCatch({
										gamlss(formula = Value ~pb(ptrans(Age, pp)), sigma.formula = ~pb(ptrans(Age, pp)), 
												nu.formula = ~1, tau.formula = ~1, family = families[fam], data = weightData, weights = weightData[,colnameProbNP], na.rm =TRUE)
									}, error = function(cond){
										message(paste("GAMLSS failed for family: ", families[fam]))
										message(cond)
										obj     <- list()
										obj$sbc <- Inf  
										return(obj)
									}
							)
							
						}else if(dep ==3){
							gamlssModels[[counter]] <- tryCatch({
										gamlss(formula = Value ~pb(ptrans(Age, pp)), sigma.formula = ~pb(ptrans(Age, pp)), 
												nu.formula = ~pb(ptrans(Age, pp)), tau.formula = ~1, family = families[fam], data = weightData, weights = weightData[,colnameProbNP], na.rm =TRUE)
									}, error = function(cond){
										message(paste("GAMLSS failed for family: ", families[fam]))
										message(cond)
										obj     <- list()
										obj$sbc <- Inf  
										return(obj)
									}
							)
							
						}else if(dep ==4){
							gamlssModels[[counter]] <- tryCatch({
										gamlss(formula = Value ~pb(ptrans(Age, pp)), sigma.formula = ~pb(ptrans(Age, pp)), 
												nu.formula = ~pb(ptrans(Age, pp)), tau.formula = ~pb(ptrans(Age, pp)), family = families[fam], data = weightData, weights = weightData[,colnameProbNP], na.rm =TRUE)
									}, error = function(cond){
										message(paste("GAMLSS failed for family: ", families[fam]))
										message(cond)
										obj     <- list()
										obj$sbc <- Inf  
										return(obj)
									}
							)
							
						}
						counter <- counter +1 
					}
				}else if(families[fam] =="BCPEo"){
					for (dep in 1:4){
						if(dep ==1){
							gamlssModels[[counter]] <- tryCatch({
										gamlss(formula = Value ~1, sigma.formula = ~1, 
												nu.formula = ~1, tau.formula = ~1, family = families[fam], data = weightData, weights = weightData[,colnameProbNP], na.rm =TRUE)
									}, error = function(cond){
										message(paste("GAMLSS failed for family: ", families[fam]))
										message(cond)
										obj     <- list()
										obj$sbc <- Inf  
										return(obj)
									}
							)
							
						}else if(dep ==2){
							gamlssModels[[counter]] <- tryCatch({
										gamlss(formula = Value ~pb(ptrans(Age, pp)), sigma.formula = ~pb(ptrans(Age, pp)), 
												nu.formula = ~1, tau.formula = ~1, family = families[fam], data = weightData, weights = weightData[,colnameProbNP], na.rm =TRUE)
									}, error = function(cond){
										message(paste("GAMLSS failed for family: ", families[fam]))
										message(cond)
										obj     <- list()
										obj$sbc <- Inf  
										return(obj)
									}
							)
							
						}else if(dep ==3){
							gamlssModels[[counter]] <- tryCatch({
										gamlss(formula = Value ~pb(ptrans(Age, pp)), sigma.formula = ~pb(ptrans(Age, pp)), 
												nu.formula = ~pb(ptrans(Age, pp)), tau.formula = ~1, family = families[fam], data = weightData, weights = weightData[,colnameProbNP], na.rm =TRUE)
									}, error = function(cond){
										message(paste("GAMLSS failed for family: ", families[fam]))
										message(cond)
										obj     <- list()
										obj$sbc <- Inf  
										return(obj)
									}
							)
							
						}else if(dep ==4){
							gamlssModels[[counter]] <- tryCatch({
										gamlss(formula = Value ~pb(ptrans(Age, pp)), sigma.formula = ~pb(ptrans(Age, pp)), 
												nu.formula = ~pb(ptrans(Age, pp)), tau.formula = ~pb(ptrans(Age, pp)), family = families[fam], data = weightData, weights = weightData[,colnameProbNP], na.rm =TRUE)
									}, error = function(cond){
										message(paste("GAMLSS failed for family: ", families[fam]))
										message(cond)
										obj     <- list()
										obj$sbc <- Inf  
										return(obj)
									}
							)
							
						}
						counter <- counter +1 
					}
				}else if(families[fam] =="LOGNO"){
					for (dep in 1:2){
							if(dep ==1){
								gamlssModels[[counter]] <- tryCatch({
											gamlss(formula = Value ~1, sigma.formula = ~1, 
													 family = families[fam], data = weightData, weights = weightData[,colnameProbNP], na.rm =TRUE)
										}, error = function(cond){
											message(paste("GAMLSS failed for family: ", families[fam]))
											message(cond)
											obj     <- list()
											obj$sbc <- Inf  
											return(obj)
										}
								)
								
							}else if(dep ==2){
								gamlssModels[[counter]] <- tryCatch({
											gamlss(formula = Value ~pb(ptrans(Age, pp)), sigma.formula = ~pb(ptrans(Age, pp)), 
													family = families[fam], data = weightData, weights = weightData[,colnameProbNP], na.rm =TRUE)
										}, error = function(cond){
											message(paste("GAMLSS failed for family: ", families[fam]))
											message(cond)
											obj     <- list()
											obj$sbc <- Inf  
											return(obj)
										}
								)
						}
					}
					counter <- counter +1 
				}
			
			}	
		}
		
		
		scores <- lapply(gamlssModels, function(x) x$sbc)
		
		finalM <- gamlssModels[[which.min(scores)]]
		
		
	}else{
		if(depOrig == 1){ # only mu
			finalM <- tryCatch({
						gamlss(formula = Value ~pb(ptrans(Age, pp)), sigma.formula = ~1, 
								nu.formula = ~1, tau.formula = ~1, family = families, data = weightData, weights = weightData[,colnameProbNP], na.rm =TRUE)
					},
					error = function(cond){
						message(paste("GAMLSS failed for family: ", families))
						message(cond)
						obj     <- list()
						obj$sbc <- Inf  
						return(obj)
					})
		}else if(depOrig == 2){  # mu and sigma
			finalM <- tryCatch({
						gamlss(formula = Value ~pb(ptrans(Age, pp)), sigma.formula = ~pb(ptrans(Age, pp)), 
								nu.formula = ~1, tau.formula = ~1, family = families, data = weightData, weights = weightData[,colnameProbNP], na.rm =TRUE)
					},
					error = function(cond){
						message(paste("GAMLSS failed for family: ", families))
						message(cond)
						obj     <- list()
						obj$sbc <- Inf  
						return(obj)
					})
		}else if(depOrig ==3){ # mu, sigma, and nu
			finalM <- tryCatch({
						gamlss(formula = Value ~pb(ptrans(Age, pp)), sigma.formula = ~pb(ptrans(Age, pp)), 
								nu.formula = ~pb(ptrans(Age, pp)), tau.formula = ~1, family = families, data = weightData, weights = weightData[,colnameProbNP], na.rm =TRUE)
					},
					error = function(cond){
						message(paste("GAMLSS failed for family: ", families))
						message(cond)
						obj     <- list()
						obj$sbc <- Inf  
						return(obj)
					})
		}else if(depOrig ==4){ # mu, sigma, and tau
			finalM <- tryCatch({
						gamlss(formula = Value ~pb(ptrans(Age, pp)), sigma.formula = ~pb(ptrans(Age, pp)), 
								nu.formula = ~pb(ptrans(Age, pp)), tau.formula = ~pb(ptrans(Age, pp)), family = families, data = weightData, weights = weightData[,colnameProbNP], na.rm =TRUE)
					},
					error = function(cond){
						message(paste("GAMLSS failed for family: ", families))
						message(cond)
						obj     <- list()
						obj$sbc <- Inf  
						return(obj)
					})
		}else {
			finalM <- tryCatch({
						gamlss(formula = Value ~pb(ptrans(Age, pp)), sigma.formula = ~pb(ptrans(Age, pp)), 
								nu.formula = ~pb(ptrans(Age, pp)), family = families, data = weightData, weights = weightData[,colnameProbNP], na.rm =TRUE)
					},
					error = function(cond){
						message(paste("GAMLSS failed for family: ", families))
						message(cond)
						obj     <- list()
						obj$sbc <- Inf  
						return(obj)
					})
		} # only mu

	
	}
	
	# traverse models and find best one using bic criterion 
	
	finalM$.user <- list(
			data = eval(finalM$call[["data"]]),
			call = match.call()
	)
	
	return(list(finalM = finalM, pp = pp))
	
}

