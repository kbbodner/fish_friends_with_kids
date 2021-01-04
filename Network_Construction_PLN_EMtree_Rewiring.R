rm(list=ls())
library(dplyr,quietly=T) #need this to load data
library(reshape2)
library(tidygraph) #need to load data
library(parallel) #need for ResampleEMtree2
library(tidyverse)
#devtools::install_github("Rmomal/EMtree") #if using newer version of R
library(EMtree)
library(PLNmodels) #need for get_data
library(Metrics) #need RMSE calculation
library(igraph) #analyzing and plotting networks
library(qgraph) #need fo graphing layout

##############################################################################

#STEP 1: Organize data

#############################################################################

#FUNCTION NAME: get.fish.data

#PURPOSE: sets-up data for running adult-only and various stage-structured PLN models:
#no stage-structure, piscivore species,largest 25% species, largest 50% species, largest 75% species

#DESCRIPTION: reads in appropriate dataset (for adult-only; for stage-structured),
#selects species with obserations > 4, creates validation dataset by subsetting last season per site,
#and combines appropriate species/stage counts for output dataset (see combinations in INPUT)

#INPUT: one of the following strings:
        # "adultOnly" if only adults are needed (adult-only)
        # "noStage" if all adults and juveniles are combined
        # "piscivore" if stage classes are only for piscivore species (non-pisc combined)
        # "25Largest" if stage classes are only for 25 largest species by adult length (smaller species combined)
        # "50Largest" if stage classes are only for 50 largest species by adult length (smaller species combined)
        # "75Largest" if stage classes are only for 75 largest species by adult length (smaller species combined)

#OUTPUT: a list with first index as dataset for building PLN models and
#        the second index as the validation dataset (to predict on)

get.fish.data <- function(ageClass) {


  if (ageClass == "adultOnly") {
    
    #reading in the data
    allData <- read.csv("data/NEON_aggregate_adultOnly.csv")
    allData <- allData[,2:ncol(allData)]
    
    #obtaining the covariance matrix (environmental variables)
    Covar <- allData[,1]
    Covar <- cbind(Covar,allData[,2:9])
    names(Covar)[1] <- "Date"
    
    #Getting count data
    Counts <- allData[,10:ncol(allData)]
    Counts[] <- lapply(Counts, function(x) as.integer(x))
    
    #number of non zero entries for each column
    nonZero <- apply(Counts, 2, function(c) sum(c!=0))
    
    #keeping those only with more than 4 entries (observed at least 5 times)
    test1 <- subset(nonZero, nonZero>4)
    #obtianing the names of those species with greater than 4 entries
    test2 <- names(test1)
    #subsetting the coutns with only those species with greater than 4 entries
    Counts <- Counts[, (names(Counts) %in% test2)]
    
    #assingning species names
    for (i in 1:ncol(Counts)) {
      
      test <- strsplit(names(Counts)[i],".",fixed=TRUE)
      newName <- paste(substring(test[[1]][1],1,1),".",test[[1]][2],sep="")
      names(Counts)[i] <- newName
      
    }
    
    
  }
  
  if (ageClass == "noStage" || ageClass == "piscivore" || ageClass == "10Largest" ||
      ageClass == "25Largest" || ageClass == "50Largest" || ageClass == "75Largest") {
  
    #reading in the data
    allData <- read.csv("data/NEON_aggregate_stageStructured.csv")
    allData <- allData[,2:ncol(allData)]
    
    #obtaining the covariance matrix (environmental variables)
    Covar <- allData[,1]
    Covar <- cbind(Covar,allData[,2:9])
    names(Covar)[1] <- "Date"
    
    #Getting count data
    Counts <- allData[,10:ncol(allData)]
    Counts[] <- lapply(Counts, function(x) as.integer(x))
    
    #number of non zero entries for each column
    nonZero <- apply(Counts, 2, function(c) sum(c!=0))
    
    #keeping those only with more than 4 entries (observed at least 5 times)
    test1 <- subset(nonZero, nonZero>4)
    #obtianing the names of those species with greater than 4 entries
    test2 <- names(test1)
    #subsetting the coutns with only those species with greater than 4 entries
    Counts <- Counts[, (names(Counts) %in% test2)]
    
    #assingning species names
    for (i in 1:ncol(Counts)) {
      
      test <- strsplit(names(Counts)[i],".",fixed=TRUE)
      newName <- paste(substring(test[[1]][1],1,1),".",test[[1]][2],sep="")
      names(Counts)[i] <- newName
      
    }
  
  }

  ageClassFunction <- function(ageClass) {
    
    if (ageClass == "noStage") {
      
      #names of the fish species with age class
      columnNamesOfFish <- colnames(Counts)
      
      for (i in 1:length(columnNamesOfFish)) {
        
        #getting before the "_"
        nameHere <- str_split(columnNamesOfFish[i],"_")[[1]][1]
        columnNamesOfFish[i] <- nameHere
        
      }
      
      #changing the names
      colnames(Counts) <- columnNamesOfFish
      Counts <- t(rowsum(t(Counts), group = colnames(Counts)))
      
    } else if (ageClass == "piscivore") {
      
      #piscivorous fish names
      piscFish = c("A.natalis","C.carolinae","L.cyanellus","L.megalotis","S.trutta",
                   "S.fontinalis","S.atromaculatus","M.dolomieu")
      
      #names of the fish species with age class
      columnNamesOfFish <- colnames(Counts)
      
      for (i in 1:length(columnNamesOfFish)) {
        
        #variable to tell if piscFish name is in column name
        identifier = FALSE
        
        for (j in 1:length(piscFish)) {
          
          #if piscivorous fish name is the current name, then set to TRUE
          if (grepl(piscFish[j], names(Counts)[i], fixed=TRUE)) {
            
            identifier = TRUE
            
          }
          
        }
        
        if (identifier == FALSE) {
          
          #getting before the "_"
          nameHere <- str_split(columnNamesOfFish[i],"_")[[1]][1]
          columnNamesOfFish[i] <- nameHere
          
        }
        
      }
      
      #changing the names
      colnames(Counts) <- columnNamesOfFish
      Counts <- t(rowsum(t(Counts), group = colnames(Counts)))
      
    } else if (ageClass == "10Largest") {
      
      #10% largest species names
      largestFish = c("S.trutta","S.fontinalis","H.etowanum","M.carinatum")
      
      #names of the fish species with stage class
      columnNamesOfFish <- colnames(Counts)
      
      for (i in 1:length(columnNamesOfFish)) {
        
        #variable to tell if 10 largest name is in column name
        identifier = FALSE
        
        for (j in 1:length(largestFish)) {
          
          #if piscivorous fish name is the current name, then set to TRUE
          if (grepl(largestFish[j], names(Counts)[i], fixed=TRUE)) {
            
            identifier = TRUE
            
          }
          
        }
        
        if (identifier == FALSE) {
          
          #getting before the "_"
          nameHere <- str_split(columnNamesOfFish[i],"_")[[1]][1]
          columnNamesOfFish[i] <- nameHere
          
        }
        
      }
      
      #changing the names
      colnames(Counts) <- columnNamesOfFish
      Counts <- t(rowsum(t(Counts), group = colnames(Counts)))
      
    } else if (ageClass == "25Largest") {
      
      #25% largest species names
      largestFish = c("A.melas","A.natalis","H.etowanum","L.aepyptera",
                      "M.carinatum","N.leptacanthus","S.fontinalis","S.trutta")
      
      #names of the fish species with age class
      columnNamesOfFish <- colnames(Counts)
      
      for (i in 1:length(columnNamesOfFish)) {
        
        #variable to tell if 25 largest name is in column name
        identifier = FALSE
        
        for (j in 1:length(largestFish)) {
          
          #if piscivorous fish name is the current name, then set to TRUE
          if (grepl(largestFish[j], names(Counts)[i], fixed=TRUE)) {
            
            identifier = TRUE
            
          }
          
        }
        
        if (identifier == FALSE) {
          
          #getting before the "_"
          nameHere <- str_split(columnNamesOfFish[i],"_")[[1]][1]
          columnNamesOfFish[i] <- nameHere
          
        }
        
      }
      
      #changing the names
      colnames(Counts) <- columnNamesOfFish
      Counts <- t(rowsum(t(Counts), group = colnames(Counts)))
      
    } else if (ageClass == "50Largest") {
      
      #50% largest species names
      largestFish = c("A.melas","A.natalis","H.etowanum","I.gagei","L.aepyptera",
                      "L.chrysocephalus","L.cyanellus","L.macrochirus","L.megalotis",
                      "M.carinatum","N.funebris","N.leptacanthus",
                      "N.leptocephalus","S.atromaculatus","S.fontinalis","S.trutta")
      
      #names of the fish species with age class
      columnNamesOfFish <- colnames(Counts)
      
      for (i in 1:length(columnNamesOfFish)) {
        
        #variable to tell if 50 largest name is in column name
        identifier = FALSE
        
        for (j in 1:length(largestFish)) {
          
          #if piscivorous fish name is the current name, then set to TRUE
          if (grepl(largestFish[j], names(Counts)[i], fixed=TRUE)) {
            
            identifier = TRUE
            
          }
          
        }
        
        if (identifier == FALSE) {
          
          #getting before the "_"
          nameHere <- str_split(columnNamesOfFish[i],"_")[[1]][1]
          columnNamesOfFish[i] <- nameHere
          
        }
        
      }
      
      #changing the names
      colnames(Counts) <- columnNamesOfFish
      Counts <- t(rowsum(t(Counts), group = colnames(Counts)))
      
    } else if (ageClass == "75Largest") {
      
      #75% largest species names
      largestFish = c("A.melas","A.natalis","C.anomalum","C.carolinae","C.erythrogaster","F.olivaceus",
                      "H.etowanum","I.gagei","L.aepyptera","L.chrysocephalus","L.cyanellus","L.macrochirus","L.megalotis","M.carinatum","N.exilis",
                      "N.funebris","N.leptacanthus","N.leptocephalus","N.stramineus","R.atratulus","S.atromaculatus","S.fontinalis","S.trutta")
      
      #names of the fish species with age class
      columnNamesOfFish <- colnames(Counts)
      
      for (i in 1:length(columnNamesOfFish)) {
        
        #variable to tell if 70 largest name is in column name
        identifier = FALSE
        
        for (j in 1:length(largestFish)) {
          
          #if piscivorous fish name is the current name, then set to TRUE
          if (grepl(largestFish[j], names(Counts)[i], fixed=TRUE)) {
            
            identifier = TRUE
            
          }
          
        }
        
        if (identifier == FALSE) {
          
          #getting before the "_"
          nameHere <- str_split(columnNamesOfFish[i],"_")[[1]][1]
          columnNamesOfFish[i] <- nameHere
          
        }
        
      }
      
      #changing the names
      colnames(Counts) <- columnNamesOfFish
      Counts <- t(rowsum(t(Counts), group = colnames(Counts)))
      
    }
    
    return(Counts)
  
  }
  
  #Evaluating only the age class of interest
  Counts <- ageClassFunction(ageClass)
  
  Covar$Date <- as.Date(Covar$Date)
  Covar$SiteName <- as.factor(Covar$SiteName)
  
  #getting Covar and Counts together so that I can subset them
  test <- cbind(Covar,Counts)
  
  ###########################################
  #Isolating for specific sites to predict on
  ###########################################
  
  #only for Prin sites
  prinData <- subset(test,SiteName=="Prin")
  
  predictOnPrin<- prinData %>% 
    mutate(Date=as.Date(Date, format= "%y-%m-%d"))%>% 
    group_by(SiteName) %>%  
    arrange(desc(Date)) %>%
    slice(c(1))
  
  #build upon
  buildOnPrin<- prinData %>% 
    mutate(Date=as.Date(Date, format= "%y-%m-%d"))%>% 
    group_by(SiteName) %>%  
    arrange(desc(Date)) %>%
    slice(-c(1))
  
  #only for mcdi sites 
  mcdiData <- subset(test,SiteName=="Mcdi")
  
  predictOnMcdi<- mcdiData %>% 
    mutate(Date=as.Date(Date, format= "%y-%m-%d"))%>% 
    group_by(SiteName) %>%  
    arrange(desc(Date)) %>%
    slice(c(1,2,3))
  
  buildOnMcdi<- mcdiData %>% 
    mutate(Date=as.Date(Date, format= "%y-%m-%d"))%>% 
    group_by(SiteName) %>%  
    arrange(desc(Date)) %>%
    slice(-c(1,2,3))
  
  #only for mayf sites
  mayfData <- subset(test,SiteName=="Mayf")
  
  predictOnMayf<- mayfData %>% 
    mutate(Date=as.Date(Date, format= "%y-%m-%d"))%>% 
    group_by(SiteName) %>%  
    arrange(desc(Date)) %>%
    slice(c(1,2,3,4))
  
  buildOnMayf<- mayfData %>% 
    mutate(Date=as.Date(Date, format= "%y-%m-%d"))%>% 
    group_by(SiteName) %>%  
    arrange(desc(Date)) %>%
    slice(-c(1,2,3,4))
  
  #only for walk sites 
  walkData <- subset(test,SiteName=="Walk")
  
  predictOnWalk<- walkData %>% 
    mutate(Date=as.Date(Date, format= "%y-%m-%d"))%>% 
    group_by(SiteName) %>%  
    arrange(desc(Date)) %>%
    slice(c(1,2,3))
  
  buildOnWalk<- walkData %>% 
    mutate(Date=as.Date(Date, format= "%y-%m-%d"))%>% 
    group_by(SiteName) %>%  
    arrange(desc(Date)) %>%
    slice(-c(1,2,3))
  
  #only for king sites 
  kingData <- subset(test,SiteName=="King")
  
  predictOnKing<- kingData %>% 
    mutate(Date=as.Date(Date, format= "%y-%m-%d"))%>% 
    group_by(SiteName) %>%  
    arrange(desc(Date)) %>%
    slice(c(1,2,3,4))
  
  buildOnKing<- kingData %>% 
    mutate(Date=as.Date(Date, format= "%y-%m-%d"))%>% 
    group_by(SiteName) %>%  
    arrange(desc(Date)) %>%
    slice(-c(1,2,3,4))
  
  #only for hopb sites
  hopbData <- subset(test,SiteName=="Hopb")
  
  predictOnHopb<- hopbData %>% 
    mutate(Date=as.Date(Date, format= "%y-%m-%d"))%>% 
    group_by(SiteName) %>%  
    arrange(desc(Date)) %>%
    slice(c(1,2,3,4))
  
  buildOnHopb<- hopbData %>% 
    mutate(Date=as.Date(Date, format= "%y-%m-%d"))%>% 
    group_by(SiteName) %>%  
    arrange(desc(Date)) %>%
    slice(-c(1,2,3,4))
  ########################
  
  #The set of data to construct the PLN network on
  buildOn <- rbind(buildOnHopb,buildOnKing,buildOnMayf,buildOnMcdi,buildOnPrin,buildOnWalk)
  
  #The set of data to predict the PLN network on
  predictOn <- rbind(predictOnHopb,predictOnKing,predictOnMayf,predictOnMcdi,predictOnPrin,predictOnWalk)
  
  #obtaining the covariance matrix
  CovarBuild <- buildOn[,1:9]
  
  #Getting count data
  CountsBuild <- buildOn[,10:ncol(buildOn)] 
  
  #obtaining the covariance matrix to predict
  CovarPredict <- predictOn[,1:9]
  
  #Getting count data to predict
  CountsPredict <- predictOn[,10:ncol(predictOn)]
  
  CovarBuild$Date <- as.factor(CovarBuild$Date)
  CovarPredict$Date <- as.factor(CovarPredict$Date)
  
  fishSpecies <- prepare_data(CountsBuild,CovarBuild, offset = "TSS")
  fishSpeciesPredicted <- prepare_data(CountsPredict,CovarPredict, offset = "TSS")
  
  return(list(fishSpecies,fishSpeciesPredicted))
  
}


#########################################################################
#Running get.fish.data function using NEON fish data

#adult-only
df_adult.only<-get.fish.data("adultOnly")

#stage-structured
df_no.stage<-get.fish.data("noStage")
df_pisc<-get.fish.data("piscivore")
df_25.Largest<-get.fish.data("25Largest")
df_50.Largest<-get.fish.data("50Largest")
df_75.Largest<-get.fish.data("75Largest")

#######################################################################

#STEP 2: BUILDING PLN MODELS

###########################################################################
#FUNCTION NAME: modelMetrics

#PURPOSE: Used in PLNBuilding Function to obtain model metrics for each model

#DESCRIPTION: calculates RMSE and RMSE on counts>0 and combines with PLNmodel in-package evaluations

#INPUT: A PLNmodel, data for predicted model fishSpeciesPredicted

#OUTPUT: dataframe of a model's metrics including RMSE, RMSE>0, BIC, R-squared, log-likelihood,
#number of parameters, Integrated Classification Likelihood (ICL)

modelMetrics <- function(PLNMODEL,fishSpeciesPredicted) {
  
  #predicted values for all
  predictedValues <- predict(PLNMODEL, newdata = fishSpeciesPredicted, type = "response")
  predictedValuesPerSpecies <- predictedValues
  predictedValues <- as.vector(predictedValues)
  observed <- as.vector(fishSpeciesPredicted$Abundance)
  
  #predicted values of only those greater than 0
  isolatedPredictedValues <- predictedValues[fishSpeciesPredicted$Abundance > 0]
  isolatedObservedValues <- fishSpeciesPredicted$Abundance[fishSpeciesPredicted$Abundance > 0]
  
  #dataframe to return with testing values
  test <- PLNMODEL$criteria
  test$RMSE_all <- rmse(observed,predictedValues)
  test$RMSE_obs_greater_0 <- rmse(isolatedObservedValues,isolatedPredictedValues)
  
  return(test)
  
} #used to be called sample testing


#FUNCTION NAME:PLNBuilding

#PURPOSE: create a summary of model metrics for site, site + temperature, site + dissolved oxygen

#DESCRIPTION: creates three PLN models specified above using adult-only or stage-structured data
# and environmental data  and calculates model performance; calls modelMetrics()

#INPUT: dataframe of either adult-only counts or the selected stage-structured counts 
#data should be output of get.fish.data() function

#OUTPUT: dataframe of model metrics of site, site + temperature, and site + oxygen's including 
#RMSE, RMSE>0, BIC, R-squared, log-likelihood, number of parameters, Integrated Classification Likelihood (ICL)

PLNBuilding<-function(df){
  
  #The data that the model is built on
  fishSpecies <- df[[1]]
  
  #The data that the model is used to predict on
  fishSpeciesPredicted <- df[[2]]
  
  #PLNmodels that you would like to evaluate
  #Here: covariates site name, water temperature and dissolved oxygen
  modelSiteName <-PLN(Abundance ~ SiteName + offset(log(Offset)),data = fishSpecies, control = list(covariance = "full", trace = 0, algorithm="LBFGS"))
  modelSiteNameTemp <-PLN(Abundance ~ waterTemp + SiteName + offset(log(Offset)),data = fishSpecies, control = list(covariance = "full", trace = 0, algorithm="LBFGS"))
  modelOxySiteName<-PLN(Abundance ~ dissolvedOxygen + SiteName + offset(log(Offset)), data = fishSpecies, control = list(covariance = "full", trace = 0, algorithm="LBFGS"))
  
  #Printing the in an out of sample testing
  rbind(
    
    modelMetrics(modelSiteName, fishSpeciesPredicted),
    modelMetrics(modelSiteNameTemp, fishSpeciesPredicted),
    modelMetrics(modelOxySiteName, fishSpeciesPredicted)
    
  ) %>%
    as.data.frame(row.names=c("SiteName","SiteName+Temp","SiteName+Oxy")) %>%
    knitr::kable()
}

######################################################
#Build + Test PLN models our models

#adult PLN models
PLNBuilding(df_adult.only)

#stage-structured PLN models
PLNBuilding(df_no.stage)
PLNBuilding(df_pisc)
PLNBuilding(df_50.Largest)
PLNBuilding(df_25.Largest)
PLNBuilding(df_75.Largest)

####################################################################

#STEP 3: EMtree

####################################################################

#FUNCTION NAME: ResampleEMtree2

#PURPOSE: EMtree method modified to use the LBFGS method

#DESCRIPTION: Runs EMtree method (resampling procedure for edges probability) using the LBFGS method
#(see ??ResampleEMtree from the EMTree package description for more information)

#INPUTS: 
#counts - Data of observed counts with dimensions n x p, either a matrix, data.frame or tibble
#covar_matrix - matrix of covariates, should have the same number of rows as the count matrix.
#S - Total number of wanted sub-samples.
#maxIter - Maximum number of EMtree iterations at each sub-sampling.
#cond.tol - Tolerance for the psi matrix.
#cores - 	Number of cores, can be greater than 1 if data involves less than about 32 species.
#v - The proportion of observed data to be taken in each sub-sample. It is the ratio (sub-sample size)/n

#OUPUTS: a list which contains a dataframe, vector of EMtree maximum iteractions,and 
#vector of running times in each resampling

ResampleEMtree2 <- function (counts, covar_matrix = NULL, O = NULL, v = 0.8, S = 100, 
                             maxIter = 30, cond.tol = 1e-14, cores = 3) {
  counts = as.matrix(counts)
  n = nrow(counts)
  p = ncol(counts)
  P = p * (p - 1)/2
  V = round(v * n)
  Pmat = matrix(0, S, P)
  if (is.null(O)) {
    O = matrix(1, n, p)
  }
  if (is.null(covar_matrix)) {
    X = matrix(1, nrow = n, ncol = 1)
  }
  else {
    X = as.matrix(covar_matrix)
  }
  obj <- mclapply(1:S, function(b) {
    cat("\nS=", b, " ")
    set.seed(b)
    sample = sample(1:n, V, replace = F)
    counts.sample = counts[sample, ]
    X.sample = data.frame(X[sample, ])
    O.sample = O[sample, ]
    try({
      suppressWarnings(PLN.sample <- PLN(counts.sample ~ 
                                           -1 + offset(log(O.sample)) + ., data = X.sample, 
                                         control = list(trace = 0,algorithm="LBFGS")))
      inf1 <- EMtree(PLN.sample, maxIter = maxIter, cond.tol = cond.tol, 
                     verbatim = FALSE, plot = FALSE)[c("edges_prob", 
                                                       "maxIter", "timeEM", "alpha")]
      cat(" ", inf1$alpha)
      inf <- inf1[c("edges_prob", "maxIter", "timeEM")]
    }, silent = TRUE)
    if (!exists("inf")) 
      inf = list(edges_prob = matrix(NA, 1, P), maxIter = 0, 
                 timeEM = 0)
    return(inf)
  }, mc.cores = cores)
  lengths <- sapply(obj, function(x) {
    length(x)
  })
  if (mean(lengths) != 3) {
    indices <- which(lengths != 3)
    lapply(indices, function(x) {
      set.seed(x)
      sample = sample(1:n, V, replace = F)
      counts.sample = counts[sample, ]
      X.sample = data.frame(X[sample, ])
      O.sample = O[sample, ]
      suppressWarnings(PLN.sample <- PLN(counts.sample ~ 
                                           -1 + offset(log(O.sample)) + ., data = X.sample, 
                                         control = list(trace = 0,algorithm="LBFGS")))
      cat("need 2nd PLN")
      obj[[x]] <<- EMtree(PLN.sample, maxIter = maxIter, 
                          cond.tol = cond.tol, verbatim = FALSE, plot = FALSE)[c("edges_prob", 
                                                                                 "maxIter", "timeEM")]
    })
  }
  Pmat <- do.call(rbind, lapply(obj, function(x) {
    F_Sym2Vec(x$edges_prob)
  }))
  summaryiter = do.call(c, lapply(obj, function(x) {
    x$maxIter
  }))
  times <- do.call(c, lapply(obj, function(x) {
    x$timeEM
  }))
  return(list(Pmat = Pmat, maxIter = summaryiter, times = times))
}

#FUNCTION NAME: F_Sym2Vec
#PURPOSE:  Required to perform ResampleEMtree2
#INPUT: matrix
#DESCRIPTION: subsets and returns the lower triangular matrix values
#OUTPUT: vector
F_Sym2Vec <- function (A.mat) 
{
  return(A.mat[lower.tri(A.mat)])
}

########################################
#Run EMTree with adult-only and best-performing stage-structured model, larger species (50%)

#The abundance of the fish that PLN model builds on
counts_adultOnly = as.matrix(df_adult.only[[1]]$Abundance)
counts_50.largest=as.matrix(df_50.Largest[[1]]$Abundance)

#The covariates that were measured at time of sampling
covar_adultOnly = as_tibble(df_adult.only[[1]])
covar_50.largest = as_tibble(df_50.Largest[[1]])


#siteName and dissolvedOxygen were used as they performed best in PLN model building and selection
#If user decides to return to native ResampleEMtree, use "ResampleEMtree"
resample_output_adultOnly<-ResampleEMtree2(counts=counts_adultOnly, covar_matrix = c(covar_adultOnly$SiteName,covar_adultOnly$dissolvedOxygen), S=5, maxIter=50, cond.tol=1e-8, cores=1, v=0.70)

resample_output_50.largest<-ResampleEMtree2(counts=counts_50.largest, covar_matrix = c(covar_50.largest$SiteName,covar_50.largest$dissolvedOxygen), S=5, maxIter=50, cond.tol=1e-8, cores=1, v=0.70)


#The threshold (between 0-1) for assuming an interaction
#we selected 0.7 as a conservative guideline (before our network became disconnected)
thresholdValue = 0.7

#Obtaining the network with a given threshold
adult_network<-freq_selec(resample_output_adultOnly$Pmat,Pt=thresholdValue)
largest.fifty_network<-freq_selec(resample_output_50.largest$Pmat,Pt=thresholdValue) #50 largest


#################################################################
#FUNCTION: networkPlot

#PURPOSE: Create network visualization

#DESCRIPTION: creates adjacency matrix, and edge list, to visualize network, via
# the Fruchterman-Reingold layout approach in igraph.

#INPUT: matrix of connections, and species counts with species names/stages as colnames

#OUTPUT: network figure (igraph plot object) where edge thickness corresponds to higher 
# probabilities, and node size corresponds to node degree

networkPlot<-function(network, counts){
  g <- graph.adjacency(network,mode=c("undirected"),weighted=TRUE)
  e <- get.edgelist(g)
  l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g),
                                         area=8*(vcount(g)^2),repulse.rad=(vcount(g)^3.1))
  V(g)$vertex_degree<-degree(g)
  scale_factor <- 1.0
  network.fig<-plot(g, vertex.label.dist=1.5, vertex.label.color="black",vertex.color=c("gold"),
       vertex.label.cex = 1.2, edge.width = E(g)$weight,
       vertex.size = V(g)$vertex_degree * scale_factor, asp=0,
       vertex.size=degree(g)*0.25,edge.color="red",vertex.label.dist=0.0,
       layout=l,vertex.label.degree=0,edge.curved=0.0,
       vertex.frame.color="gold",vertex.label.font=1,
       vertex.label=colnames(counts))
  return(network.fig)
}

########################################################################
#Create figure and saving leda files

adult.network.fig<-networkPlot(adult_network, counts_adultOnly)
largest.fifty_network.network.fig<-networkPlot(largest.fifty_network, counts_50.largest)

#create network of 0s for no connection and 1 for connection
adult_network[adult_network>0]<-1
largest.fifty_network[largest.fifty_network>0]<-1

#turn matrix into graph object (igraph)
adult_graph<-graph_from_adjacency_matrix(adult_network, mode="undirected", weighted=NULL)
largest.fifty_graph<-graph_from_adjacency_matrix(largest.fifty_network, mode="undirected", weighted=NULL)

#write as leda files to run GCD-11 and GCM-11 (uncomment below to write leda files)
#write_graph(adult_graph,"data/adult_network.gw","leda")
#write_graph(largerSpecies_graph,"data/fifty_network.gw","leda")

##################################################################

#STEP :Random Networks through Rewiring

##################################################################

##FUNCTION NAME:rewiring

#PURPOSE: creates n # of random adult-only or largest 50 networks

#DESCRIPTION: rewires either the adult or stage-structured network using
#sample_degseq() in igraph (ensures no loops or multiple edges) iterating n times
#creating a list of n networks and writing n csv leda files.The random networks maintain 
#the same number of nodes and degree distribution.

#INPUT:a graph object from igraph package

#OUTPUT: a list of n networks (and n csv files)

rewiring<-function(graph,stage=T,n=100){
  n.edges<-gsize(graph) #number of edges
  n.nodes<-gorder(graph) #number of nodes
  degree<-degree(graph) #degree
  
  networkList<-list() #empty list to put networks
  
  if(stage==T){
    name="rewire_fifty"}else {name="rewire_adult"}
  
  for (i in 1:100){
    rewire_network=sample_degseq(degree, method="vl")
    simple<-is.simple(rewire_network)
    while(simple==F){
      rewire_network=sample_degseq(degree, method="vl")
      simple<-is.simple(rewire_network)
    }
    
    #create names and write leda files
    nameOfFile<- paste(paste(name,i,sep="_"),".gw",sep="")
    #write_graph(rewire_network,nameOfFile,"leda") #uncomment to write leda files
    
    networkList[[i]]<-rewire_network
    names(networkList)[i] <- nameOfFile
    
  }
  return(networkList)
}

#################################################################################
#Use rewiring function to create 100 rewired adult and stage-structured networks
rewiring(adult_graph,stage=F,n=100)
rewiring(largest.fifty_graph,stage=T,n=100)






