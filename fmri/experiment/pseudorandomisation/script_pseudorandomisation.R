# PSEUDORANDOMISATION SCRIPT #

# clear workspace and set working directory
rm(list = ls())

ifelse(!dir.exists("/storage/shared/research/cinn/2018/MAGMOT/pseudorandomisation/output_v2_maxRep-3"), dir.create("/storage/shared/research/cinn/2018/MAGMOT/pseudorandomisation/output_v2_maxRep-3"), FALSE)
setwd("/storage/shared/research/cinn/2018/MAGMOT/pseudorandomisation/output_v2_maxRep-4")
setwd("/storage/shared/research/cinn/2018/MAGMOT/pseudorandomisation/output_v2_maxRep-4_noSeed")
# setwd("/Users/stefaniemeliss/cinn/2018/MAGMOT/pseudorandomisation/output_v2_maxRep-3/")
# dir.create("/Users/stefaniemeliss/Desktop/output_v2/")
dir.create("/storage/shared/research/cinn/2018/MAGMOT/pseudorandomisation/output_v2_test/")
setwd("/storage/shared/research/cinn/2018/MAGMOT/pseudorandomisation/output_v2_test")

# define key variables
numStimuli <- 36
numBlocks <- 3
numberOfrandOrdermisationsNeeded <- 25

baseOrder <-c(1:numStimuli) # creates 1 top 36

set.seed(001) # this allows to replicate the results
#set.seed(001) # this allows to replicate the results

thresholds <- c(seq(1, 0.7, by = -0.02))
#thresholds <- c(seq(0.74, 0.7, by = -0.02))

maxRep <- c(3,4) # number of maximal sequential repetitions of a factor level within a trial order
#maxRep <- c(3) # number of maximal sequential repetitions of a factor level within a trial order
#maxRep <- c(4) # number of maximal sequential repetitions of a factor level within a trial order

maxDur <- 1800 # maximum duration of iteration
maxDurLocal <- 1800 # maximum duration of local iteration

goalGlobalRange <- .70 # define aim of iteration

qualityMeasurements <- c("threshold", "it_count", "duration", "globalRange", "maxRep")

iterationSummary <- data.frame(matrix(ncol = length(qualityMeasurements), nrow = length(thresholds))) # empty data frame
names(iterationSummary) <- qualityMeasurements

qualityMeasurementsLocally <- c("maxRep", "it_count_local", "duration_local", "globalRangeOptimised", "bestGlobalRange")
iterationSummaryLocally <- data.frame(matrix(ncol = length(qualityMeasurementsLocally), nrow = length(maxRep))) # empty data frame
names(iterationSummaryLocally) <- qualityMeasurementsLocally


#LOOP over the length of maximal repetitions of a curiosity class
for (r in seq_along(maxRep)){

  # LOOP over thresholds
  for (u in seq_along(thresholds)){
    globalRange <- 1
    ptm <- proc.time() # get timimg
    it_count <- 0

    while(globalRange >= thresholds[u]){

      trial_orders <- data.frame(baseOrder) # create data frame
      it_count <-  it_count+1 # interation counter

      # loop from 1 to number of randomisations needed (defined above)
      for (i in 1:numberOfrandOrdermisationsNeeded) { # do this 25 times so that each trial order satisfies the criteria of not having to many repetitive curiosity values
        max_consec <- maxRep[r] + 1

        while (max_consec > maxRep[r]){ # as long as the number of repetitions is larger than maximum allowed
          print(it_count <- it_count+1)
          high_curiosity <- sample(1:(numStimuli/2)) # randomise the order of the high and low curiosity
          low_curiosity <- sample(((numStimuli/2)+1):numStimuli)
          trial_order <- c()

          for (block in 1:numBlocks) { # to make sure we have an equal distribution of curiosity in the blocks
            start <- (block-1) * length(high_curiosity)/numBlocks     # this gives you the start of the range of indices -1 i.e. 0,6,12
            end <- start + (length(high_curiosity)/numBlocks)         # this gives you the end of the range of indices i.e. 6, 12, 18
            block_stim <- sample(c(high_curiosity[(start+1):end], low_curiosity[(start+1):end])) # now for each block, take the 1st, 2nd, 3rd set of randomised high and low vids, then shuffle then up within the block
            trial_order <- c(trial_order, block_stim) # then add them to the trial list
          }

          max_high <- max(rle(trial_order %in% high_curiosity)$lengths) # computes the lengths and values of runs of equal values in a vector
          max_low <- max(rle(trial_order %in% low_curiosity)$lengths)
          max_consec <- max(max_high, max_low) # checks which of them is higher, used as update for the while loop
        }

        trial_orders[,i] <- data.frame(trial_order) #save each of the trial orders as a colum in the data frame

        print("saving trial order")

        colnames(trial_orders)[i] <- paste0("trialOrder", i)  # name the column accourdingly
      }

      # compute the correlation for the whole data frame
      corOfTrialOrder <- data.frame(cor(trial_orders, method = "spearman"))
      corOfTrialOrder[corOfTrialOrder == 1] <- NA # replacing all correlation of an order with itself with NA

      minOfMin <- min(corOfTrialOrder, na.rm = T)
      maxOfMax <- max(corOfTrialOrder, na.rm = T)

      globalRange <- maxOfMax - minOfMin

      print(paste("global range:", globalRange))

      dur_ptm <- proc.time() - ptm # timing update
      if (dur_ptm [3]> maxDur){ # duration in seconds
        break # if iteration takes to much time, stop
      }
    }

    write.csv(trial_orders, file = paste0("TrialOrders_threshold-", thresholds[u],"maxRep-", maxRep[r], ".csv"), row.names = F) # save what you have simulated
    write.csv(corOfTrialOrder, file = paste0("interOrderCorrelations_threshold-", thresholds[u], "maxRep-", maxRep[r],".csv"), row.names = F)

    print(paste("writing files for threshold", thresholds[u], "and maxRep", maxRep[r]))

    randomisations <- names(corOfTrialOrder)

    for (o in seq_along(randomisations)) {
      assign(paste0("trialOrder_", o), corOfTrialOrder[, o])
    }
    allCorrelationCoefficients <- reshape2::melt(corOfTrialOrder) # produces warning: No id variables; using all as measure variables
    allOrdersPooled <- allCorrelationCoefficients$value

    # create a box plot to visualise the results of the estimation
    boxplot <- boxplot(allOrdersPooled, trialOrder_1, trialOrder_2, trialOrder_3, trialOrder_4, trialOrder_5, trialOrder_6, trialOrder_7, trialOrder_8, trialOrder_9, trialOrder_10,
                       trialOrder_11, trialOrder_12, trialOrder_13, trialOrder_14, trialOrder_15, trialOrder_16, trialOrder_17, trialOrder_18, trialOrder_19, trialOrder_20,
                       trialOrder_21, trialOrder_22, trialOrder_23, trialOrder_24, trialOrder_25,
                       main = paste("boxplot threshold", thresholds[u], "+ maxRep", maxRep[r]),
                       xlab = "pseudorandomised order",
                       ylab = "Spearman rank correlation",
                       ylim = c(-.5,.5),
                       xlim = c(1,26),
                       # names = c("overall",paste0("order_", rep(1:25)))
                       names = c("overall",paste0(rep(1:25)))
    )
    abline(h = mean(allOrdersPooled, na.rm = T), col = "red")
    abline(h = sd(allOrdersPooled, na.rm = T), col = "red", lty = 3)
    abline(h = -sd(allOrdersPooled, na.rm = T), col = "red", lty = 3)

    print(boxplot)
    dev.copy(jpeg,paste0("boxplotAll_threshold-", thresholds[u], "_maxRep-", maxRep[r]), width = 1200, height = 540)
    dev.off()

    iterationSummary[u,"threshold"] <- thresholds[u]
    iterationSummary[u,"it_count"] <- it_count
    iterationSummary[u,"duration"] <- dur_ptm [3]
    iterationSummary[u,"globalRange"] <- globalRange
    iterationSummary[u,"maxRep"] <- maxRep[r]

    # write a file containing information about the simulation itself
    write.csv(iterationSummary, file = paste0("maxRep-", maxRep[r], "_iterationSummary.csv"), row.names = F)

  }

  ###########################################################################################################################################################################
  ### start LOCAL iteration
  ###########################################################################################################################################################################

  iterationSummary <- read.csv(file = paste0("maxRep-", maxRep[r], "_iterationSummary.csv")) # using the information of each of the iterations

  bestIteration <- which.min(iterationSummary$globalRange) # look which of them has the lowest range
  bestOrder <- iterationSummary[bestIteration,] # save information about the best order

  bestTrialOrder <- read.csv(file = paste0("TrialOrders_threshold-", bestOrder$threshold,"maxRep-", maxRep[r], ".csv")) # read in the data frame of the best order
  bestTrialOrderOptimised <- bestTrialOrder
  corOfBestTrialOrder <- read.csv(file = paste0("interOrderCorrelations_threshold-", bestOrder$threshold, "maxRep-", maxRep[r],".csv")) # read in the correlation matrix of the best order

  rownames(corOfBestTrialOrder) <- names(corOfBestTrialOrder)

  bestGlobalRange <- bestOrder[,"globalRange"] # which is the best global range we have obtained thus far? information from iteration summary
  print(paste("start value of global range",bestGlobalRange))

  randomisations <- names(corOfBestTrialOrder)

  ptm_local <- proc.time() # set another timer
  it_count_local <- 0 # create another iteration counter

  while ( bestGlobalRange > goalGlobalRange){

    dur_ptm_local <- proc.time() - ptm_local
    if (dur_ptm_local [3]> maxDurLocal){ # duration in seconds
      break # quit iteration if it takes too long
    }

    it_count_local <-  it_count_local+1 # update iteration counter
    set.seed(it_count_local)
    #set.seed(it_count_local)
    min <- data.frame(apply(corOfBestTrialOrder, 1, min,  na.rm = T), randomisations) # define miminum correlation for each trial order
    max <- data.frame(apply(corOfBestTrialOrder, 1, max,  na.rm = T), randomisations) # define maximum correlation for each trial order
    descriptives <- merge(min, max, by = "randomisations")
    names(descriptives) <- c("randomisations", "min", "max")
    descriptives$range <- descriptives$max - descriptives$min

    localMax <- which.max(descriptives$max) # identify the worst fitting trial order
    localMin <- which.min(descriptives$min)
    localMaxRange <- which.max(descriptives$range)

    # extrema <- c("localMax", "localMin")
    extrema <- c("localMaxRange")

    for (l in seq_along(extrema)){ # for each of the defined extremas

      e <- get(extrema[l])
      max_consec <- maxRep[r] + 1

      while (max_consec > maxRep[r]){ # iterate a new trial order using the same criteria as previously
        print(paste("iteration", it_count_local <- it_count_local+1)) # update iteration

        high_curiosity <- sample(1:(numStimuli/2)) # randomise the order of the high and low curiosity
        low_curiosity <- sample(((numStimuli/2)+1):numStimuli)
        trial_order <- c()

        for (block in 1:numBlocks) { # this loop takes four high curiosity
          start <- (block-1) * length(high_curiosity)/numBlocks     # this gives you the start of the range of indices -1 i.e. 0,6,12
          end <- start + (length(high_curiosity)/numBlocks)         # this gives you the end of the range of indices i.e. 6, 12, 18
          block_stim <- sample(c(high_curiosity[(start+1):end], low_curiosity[(start+1):end])) # now for each block, take the 1st, 2nd, 3rd set of randomised high and low vids, then shuffle then up within the block
          trial_order <- c(trial_order, block_stim) # then add them to the trial list
        }

        max_high <- max(rle(trial_order %in% high_curiosity)$lengths)
        max_low <- max(rle(trial_order %in% low_curiosity)$lengths)
        max_consec <- max(max_high, max_low)

      }

      print(paste("new trial order to replace", extrema[l])) # print if you have found a new trial order
      print(trial_order)

      bestTrialOrderOptimised[,e] <- data.frame(trial_order) # replace trial order that was too extreme with the new one
      print("saving trial order")

      colnames(bestTrialOrderOptimised)[e] <- paste0("trialOrder", e)  # name the column accourdingly

      # compute the correlation for the whole data frame to check whether the new trial order is better then the old one
      corOfBestTrialOrderOptimised <- data.frame(cor(bestTrialOrderOptimised, method = "spearman"))
      corOfBestTrialOrderOptimised[corOfBestTrialOrderOptimised == 1] <- NA # replacing all correlation of an order with itself with NA

      minOfMinOpt <- min(corOfBestTrialOrderOptimised, na.rm = T)
      maxOfMaxOpt <- max(corOfBestTrialOrderOptimised, na.rm = T)

      bestGlobalRangeOptimised <- maxOfMaxOpt - minOfMinOpt # update bestGlobalRange variable for while loop

      print(paste("updated best global range", bestGlobalRangeOptimised))

      if (bestGlobalRangeOptimised < bestGlobalRange){ # only use the new trial order if it leeds to better values for the range than before
        bestTrialOrder <- bestTrialOrderOptimised
        bestGlobalRange <- bestGlobalRangeOptimised
        corOfBestTrialOrder <- corOfBestTrialOrderOptimised
      }

      print(paste("difference between actual global range and specified value", bestGlobalRange - goalGlobalRange))
    }

  }

  write.csv(bestTrialOrder, file = paste0("TrialOrders_locallyOptimised_maxRep-", maxRep[r], ".csv"), row.names = F)
  write.csv(corOfBestTrialOrder, file = paste0("interOrderCorrelations_locallyOptimised_maxRep-", maxRep[r],".csv"), row.names = F)

  print(paste("writing files for",  maxRep[r], "locally optimised"))

  randomisations <- names(corOfBestTrialOrder)

  for (o in seq_along(randomisations)) {
    assign(paste0("trialOrder_", o), corOfBestTrialOrder[, o])
  }
  allCorrelationCoefficients <- reshape2::melt(corOfBestTrialOrder) # produces warning: No id variables; using all as measure variables
  allOrdersPooled <- allCorrelationCoefficients$value

  # create a box plot to visualise the results of the estimation
  boxplot <- boxplot(allOrdersPooled, trialOrder_1, trialOrder_2, trialOrder_3, trialOrder_4, trialOrder_5, trialOrder_6, trialOrder_7, trialOrder_8, trialOrder_9, trialOrder_10,
                     trialOrder_11, trialOrder_12, trialOrder_13, trialOrder_14, trialOrder_15, trialOrder_16, trialOrder_17, trialOrder_18, trialOrder_19, trialOrder_20,
                     trialOrder_21, trialOrder_22, trialOrder_23, trialOrder_24, trialOrder_25,
                     main = paste("boxplot maxRep", maxRep[r], "locally optimised"),
                     xlab = "pseudorandomised order",
                     ylab = "Spearman rank correlation",
                     ylim = c(-.5,.5),
                     xlim = c(1,26),
                     # names = c("overall",paste0("order_", rep(1:25)))
                     names = c("overall",paste0(rep(1:25)))
  )
  abline(h = mean(allOrdersPooled, na.rm = T), col = "red")
  abline(h = sd(allOrdersPooled, na.rm = T), col = "red", lty = 3)
  abline(h = -sd(allOrdersPooled, na.rm = T), col = "red", lty = 3)

  print(boxplot)
  dev.copy(jpeg,paste0("boxplotAll_locallyOptimised_maxRep-", maxRep[r]), width = 1200, height = 540)
  dev.off()

  iterationSummaryLocally[r,"maxRep"] <- maxRep[r]
  iterationSummaryLocally[r,"it_count_local"] <- it_count_local
  iterationSummaryLocally[r,"duration_local"] <- dur_ptm_local [3]
  iterationSummaryLocally[r,"globalRangeOptimised"] <- bestGlobalRangeOptimised
  iterationSummaryLocally[r,"bestGlobalRange"] <- bestGlobalRange

  # write a file containing information about the simulation itself
  write.csv(iterationSummaryLocally, file = paste0("iterationSummaryLocally.csv"), row.names = F)


}
