deniz <- function(allRes = NULL,  withCI =TRUE, RIperc =c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975), 
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
  axis(side = 1, at = pretty(xlim,n = 17), labels =pretty(xlim, n = 17), lwd.ticks = 1)
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








defineAgeGroupsWithTolBS <- function(Data = NULL, covarName = "Age", colnameID ="PID", colnameValue = "Value", minN = 1000, unique =TRUE, tol = 0.01){
  
  x <- list()
  class(x) <- "RWDRICurve"
  
  x$Params <- list(method="GAMLSS_refineR", covarName = covarName)	
  
  minAge <- min(Data[,covarName])
  maxAge <- max(Data[,covarName])
  
  # compute size of age groups wilt tolerance	
  ageGroups <- data.frame(ind = 1, min = 1, max = 1)
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




defineAgeGroupsWithTol <- function(Data = NULL, covarName = "Age", colnameID ="PID", colnameValue = "Value", minN = 1000, unique =TRUE, tol = 0.01){
  
  
  
  Data <- inputData
  
  print("sdsdsds") 
  covarName = "Age"
  colnameID ="PID"
  colnameValue = "Value"
  tol = 0.01
  unique =TRUE
  minN = 1000
  
  x <- list()
  class(x) <- "RWDRICurve"
  
  x$Params <- list(method="GAMLSS_refineR", covarName = covarName)	
  
  minAge <- min(Data[,covarName])
  maxAge <- max(Data[,covarName])
  
  # compute size of age groups with tolerance	
  ageGroups <- data.frame(ind = 1, min = 1, max = 1)
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
  
  # traverse through age groups and fill until minN is reached r <-1
  for(r in 1:nrow(ageGroups)){
    
  print(r)
    
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
