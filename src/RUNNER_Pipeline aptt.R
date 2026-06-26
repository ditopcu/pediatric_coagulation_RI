# load or install required packages
library(refineR)
library(gamlss)
library(parallel)
library(future)
library(future.apply)

# increase memory limit for parallelization
options('future.globals.maxSize' = 1500*1024^2)



# load required scripts
source("src/pipeline/utils.R")
source("src/pipeline/algoRICurves.R")
source("src/pipeline/plotRICurve.R")


# set output directory 
outputDir <- "cont_out/"

### read input data, adapt path to file and possibly function for reading the data
# the input data should have at least the following columns: 
# age (or continuous covariate) 
# value (concentration/ activity value for biomarker)
# patient id (or other column to identify if one sample originated from the same subject) 



m <- function() {

    
  inputData <- read.table(file = "data/coa_results/CSV cont aPTT 2025.10.17.csv", header = T, sep =";",
                                 stringsAsFactors =F)
    
    # inputData <- unfilteredData[unfilteredData$Age >= 1, ]
    # 
    #  
    # inputData  |> tibble::tibble()
    
    
    
    
    # specify the names as they occur in your input dataset
    covarName 		<- "Age"		# column with continuous covariate, e.g. Age
    colnameID 		<- "PID"		# column with unique subject identifier, e.g. PID
    colnameValue 	<- "Value"		# column with concentration/ activity values, e.g. Value
    
    # set core filename that will be extended for individual files that are written throughout the pipeline, e.g. name of the biomarker 
    # NOTE: if results with the same core filename already exist in the specified output directory, the existing results will be used and if some results are missin
    # 		these will be calculated. 
    filename <- "aPTT"
    
    # set number of bootstrap iterations, default = 5. Note, that increasing the number of bootstrap iterations, 
    # also increase the computation time. 
    NBootstrap 		<- 5
    
    # if you only want to use a subset of available cores set NCores to specified number
    NCores 			<- NULL
    
    # if you want to get the estimates for different percentiles, adapt RIperc here 
    RIperc <- c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)
    # RIperc <- c(0.025,  0.5, 0.975)
    
    #==========================================================================================================================================
    #==========================================================================================================================================
    #==========================================================================================================================================
    # start pipeline
    
    cat(paste0(Sys.time(), ": Pipeline started. \n"))
    
    seed <- 123
    set.seed(seed)
    
    cat(paste0(Sys.time(), ": Start of estimation for point estimate. \n"))
    
    pointEstGamlss <- runPipelinePointEst(inputData = inputData, covarName = covarName, colnameID = colnameID,  
    		colnameValue = colnameValue, outputDir = outputDir, filename = filename, NCores = NCores) 
    
    cat(paste0(Sys.time(), ": Point Estimation completed. \n"))
    cat(paste0(Sys.time(), ": Start of bootstrap iterations. \n"))
    
    bootstrapRes <- runPipelineCIs(inputData = inputData, pointEstGamlss = pointEstGamlss, NBootstrap = 5, 
    		covarName = covarName, colnameID = colnameID, colnameValue = colnameValue, outputDir = outputDir, filename = filename, NCores = NCores)
    
    cat(paste0(Sys.time(), ": Bootstrap iterations completed. \n"))
    
    # pipeline estimation finished
    #===========================================================================================================================================
    # evaluation and plotting 
    
    covarValue = c(seq(min(inputData[,covarName]), max(inputData[,covarName]),1))
    
    cat(paste0(Sys.time(), ": Computing continuous reference intervals and percentile curves with confidence intervals. \n"))
    
    allRes <- estimateCIs(pointEst = pointEstGamlss, estBS = bootstrapRes, covarValue = covarValue, withData = TRUE, 
                          RIperc = RIperc, CIprop = 0.90, onlyCI =FALSE, 
    		RICurve =NULL, RIMats = NULL)
    
    
    cat(paste0(Sys.time(), ": Plotting results \n"))
    
    
    
    plotRICurveEstimates(allRes = allRes, withCI = TRUE, RIperc = RIperc)
    
    deniz(allRes = allRes, withCI = TRUE, RIperc = RIperc, xlim = c(1,17),xlab = "Yaş",ylab = "aPTT, saniye",
                         cols =  c( "red",   "white",  "white",  "blue",
                                  "white",  "white",  "red")
                         
                         )
 
 
    
    # plot estimates on linear scale with confidence intervals (withCI = TRUE)
    png(filename = file.path(outputDir, paste0(filename,  "_", "CIs_linearScale.png")), width = 1200, height = 700,res = 100)
    # plotRICurveEstimates(allRes = allRes, withCI = TRUE, RIperc = RIperc)
    deniz(allRes = allRes, withCI = TRUE, RIperc = RIperc, xlim = c(1,17),xlab = "Yaş",ylab = "aPTT, saniye",
          cols =  c( "red",   "white",  "white",  "blue",
                     "white",  "white",  "red")
          
    )
    dev.off()
    
     
    
    tb <- allRes$RICurve$RIMat  |> 
      data.frame() |> 
      tibble::rownames_to_column(var = "Percentile") |> 
      tibble::tibble() 
    
    tb |> clipr::write_clip()
        

}
