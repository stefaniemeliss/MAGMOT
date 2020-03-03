## analysis of MAGMOT DATA

#### setups ####
#empty work space, load libraries and functions
rm(list=ls())

# define necessary directories
mainDir <- "~/Dropbox/Reading/PhD/Magictricks/fmri_study"
subDirData <- "Data"
version <- "MAGMOT"
dataDir <- file.path(mainDir, subDirData) #"~/Dropbox/Reading/PhD/Magic tricks/fmri_study/Data"
preprocessedDir <- file.path(dataDir, "preprocessed") #"~/Dropbox/Reading/PhD/Magic tricks/behavioural_study/data_fin/MagicBehavioural_preprocessed"
codedDir <- file.path(dataDir, "magicmemory_fmri", "coded") #"~/Dropbox/Reading/PhD/Magic tricks/behavioural_study/data_fin/coded/preprocessing"

analysisDir <- file.path(mainDir, "Analysis")
ratingsDir <- file.path(analysisDir, "Ratings")
memoryDir <- file.path(ratingsDir, "Memory")
brainDir <- file.path(ratingsDir, "BrainBehaviour")
questDir <- file.path(ratingsDir, "Questionnaires")
tricksDir <- file.path(analysisDir, "Tricks")

pooledDir <-  "~/Dropbox/Reading/PhD/Magictricks/pooled_analyses"

pooled <- 1

# check whether these directories exist, if not create them
ifelse(!dir.exists(ratingsDir), dir.create(ratingsDir), FALSE) 
ifelse(!dir.exists(memoryDir), dir.create(memoryDir), FALSE) 
ifelse(!dir.exists(questDir), dir.create(questDir), FALSE) 
ifelse(!dir.exists(brainDir), dir.create(brainDir), FALSE) 
ifelse(!dir.exists(analysisDir), dir.create(analysisDir), FALSE) 
ifelse(!dir.exists(tricksDir), dir.create(tricksDir), FALSE) 


memoryLevels <- c("cuedRecallStrict", "cuedRecallLenient", 
                  "recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf", 
                  "rememberedStrictAboveAvg", "rememberedLenientAboveAvg", "rememberedStrictHigh", "rememberedLenientHigh")
memoryLabels <- c("cuedRecallStrict", "cuedRecallLenient", 
                  "allConf", "highConf", "aboveAvgConf", 
                  "rememberedStrictAboveAvg", "rememberedLenientAboveAvg", "rememberedStrictHigh", "rememberedLenientHigh")

#helper functions and packages #
source("~/Dropbox/Reading/Codes and functions/R/errorbars.R")
source("~/Dropbox/Reading/Codes and functions/R/rbindcolumns.R")
library(ggplot2)
library(lmerTest)
library(reshape2)

### read in data sets ###
setwd(preprocessedDir)
dfWide <- xlsx::read.xlsx("wide_MAGMOT.xlsx", sheetName = "Sheet1")
dfLong <- xlsx::read.xlsx("long_MAGMOT.xlsx", sheetName = "Sheet1")


#####################################################################################################################################
############################################### ANALYSIS BASED ON DATA IN LONG FORMAT ############################################### 
#####################################################################################################################################

########## 1. get descriptives for curiosity and memory for whole sample as well as for each group individually ########## 
if (exists("dfLong") == F) {
  setwd(preprocessedDir)
  dfLong <- xlsx::read.xlsx("long_MAGMOT.xlsx", sheetName = "Sheet1")#,  showWarnings = FALSE)
}

# define dependent variables
workspace <- list.files(path = file.path(codedDir), pattern = "_CP.csv") # check whether the data is coded yet or not
if(length(workspace) == 0) { # if data is not coded yet, only look at recognition performance
  dependentVariables <- c("responseCuriosity", "rtCuriosity", "rtAnswer",
                          "recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
                          "confidence", "confidenceCorrectTrials")
} else {
  dependentVariables <- c("responseCuriosity", "rtCuriosity", "rtAnswer",
                          "cuedRecallStrict", "cuedRecallLenient", 
                          "recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
                          "rememberedStrictAboveAvg", "rememberedLenientAboveAvg", "rememberedStrictHigh", "rememberedLenientHigh",
                          "confidence", "confidenceCorrectTrials")
}

# for each of the variables defined above, compute the descriptive statistics
for(DV in 1:length(dependentVariables)) {
  
  print(dependentVariables[DV])
  
  # all
  descriptive <- psych::describe(dfLong[,dependentVariables[DV] ])
  row.names(descriptive) <- c(paste0(dependentVariables[DV], "_all"))
  # within each group
  descriptive_groupwise <-   by(cbind(dfLong[,dependentVariables[DV]]), dfLong$group, psych::describe)
  descriptive_groupwise_named <- rbind(descriptive_groupwise[[1]], descriptive_groupwise[[2]])
  row.names(descriptive_groupwise_named) <- c(paste0(dependentVariables[DV], "_", names(descriptive_groupwise)[1]),paste0(dependentVariables[DV], "_", names(descriptive_groupwise)[2]) )
  
  # combine output for whole sample and each group
  if (DV == 1) {
    descriptives <- rbind(descriptive, descriptive_groupwise_named)
  } else {
    temp_descriptives <- rbind(descriptive, descriptive_groupwise_named)
    descriptives <- rbind.all.columns(descriptives, temp_descriptives) #rbind all columns will induce NA if there was initially no data saved in the loop per participant
    rm(temp_descriptives)
  }
  # delete unncessary variables
  rm(descriptive, descriptive_groupwise, descriptive_groupwise_named)
}
# round result and save
descriptives <- round(descriptives, digit = 3)
setwd(ratingsDir)
write.csv(descriptives, paste0("Descriptives_dependentVariables_", version, ".csv"))
if (pooled == 1){
  setwd(pooledDir)
  xlsx::write.xlsx(descriptives, file="Descriptives_dependentVariables.xlsx", sheetName = paste(version), append = T) # note: row.names contain variables
}
rm(descriptives)

########## 2. compute mean memory performance for each magic trick ########## 
# define variables for which the mean per trick should be calculated
indicesPerTrick <- c("answerRT", "responseCuriosity", "rtCuriosity",
                     memoryLevels,
                     "confidence", "confidenceCorrectTrials")
indicesPerTrickMean <- paste0("mean_", indicesPerTrick)

# calculate mean values for each magic trick for each index
for (iMean in seq_along(indicesPerTrickMean)){
  
  # put data from long into wide format
  assign(paste(indicesPerTrick[iMean]), reshape::cast(dfLong, ID~stimID,value=paste0(indicesPerTrick[iMean])))
  
  # calculate mean value for each magic trick
  meansPerTrick <-   colMeans(get(indicesPerTrick[iMean]), na.rm = T) # calculate mean for each magic trick
  
  rm(list = indicesPerTrick[iMean])
  
  # combine all means in one data frame
  if (iMean == 1){
    dfMeans <- data.frame(meansPerTrick)
    names(dfMeans) <- indicesPerTrickMean[iMean]
    dfMeans$stimID <- row.names(dfMeans)
  } else {
    dfMeans_temp <- data.frame(meansPerTrick)
    names(dfMeans_temp) <- indicesPerTrickMean[iMean]
    dfMeans_temp$stimID <- row.names(dfMeans_temp)
    dfMeans <- merge(dfMeans, dfMeans_temp, by = "stimID")
    rm(dfMeans_temp)
  }
}

setwd(file.path(analysisDir))
write.csv(dfMeans, paste0("stimuli_MagicBehavioural_memoryPerformance_", version, ".csv"), row.names = F)
if (pooled == 1){
  setwd(pooledDir)
  xlsx::write.xlsx(dfMeans, file="stimuli_MagicBehavioural_memoryPerformance.xlsx", sheetName = paste(version), row.names = F, append = T)
}
rm(dfMeans, meansPerTrick, indicesPerTrick, indicesPerTrickMean)



########## 3. Compute lmer model predicting memory performance using curiosity as a continous variable and group effect coded ########## 

# define dependent variables 
if(length(workspace) == 0) { # if data is not coded yet, only look at recognition performance
  DV_LME <- c("recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
              "confidence", "confidenceCorrectTrials")
} else {
  DV_LME <- c("cuedRecallStrict", "cuedRecallLenient", 
              "recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
              "rememberedStrictAboveAvg", "rememberedLenientAboveAvg", "rememberedStrictHigh", "rememberedLenientHigh",
              "confidence", "confidenceCorrectTrials")
}

##### 3.1 basic LME predicting memory performance with mean-centered curiosity, effect-coded reward and their interaction #####

# loop over dependent variables to compute LME
for (DV in 1:length(DV_LME)){
  # define model based on scaling of variable
  if (max(dfLong[, DV_LME[DV]], na.rm = T) > 1){ # if DV continuous
    LMEmodel_curiosityContinuous <- lmerTest::lmer(dfLong[, DV_LME[DV]] ~ groupEffectCoded*curiosityGroupMeanCentered + (1+curiosityGroupMeanCentered|ID) + (1|stimID), data = dfLong)
  } else { # if DV binary
    LMEmodel_curiosityContinuous <- glmer(dfLong[, DV_LME[DV]] ~ groupEffectCoded*curiosityGroupMeanCentered + (1+curiosityGroupMeanCentered|ID) + (1|stimID), family = "binomial"(link = 'logit'), data = dfLong)
  }
  # put coefficients into data frame
  summaryCuriosityContinuous <- summary(LMEmodel_curiosityContinuous)
  summaryCuriosityContinuousCoefficients <- as.data.frame(summaryCuriosityContinuous$coefficients)
  # change row names
  row.names(summaryCuriosityContinuousCoefficients) <- c(paste0("LME_",DV_LME[DV], "_Intercept"), paste0("LME_",DV_LME[DV], "_groupEffectCoded"), 
                                                         paste0("LME_",DV_LME[DV], "_curiosityGroupMeanCentered"), 
                                                         paste0("LME_",DV_LME[DV], "_interaction"))
  # put coeffiencts from all DVs into one data frame
  if (DV == 1) {
    LMEresults <- summaryCuriosityContinuousCoefficients
  }else {
    temp_LMEresults <- summaryCuriosityContinuousCoefficients
    LMEresults <- rbind.all.columns(LMEresults, temp_LMEresults) #rbind all columns will induce NA if there was initially no data saved in the loop per participant
    rm(temp_LMEresults)
  }
  rm(summaryCuriosityContinuousCoefficients, summaryCuriosityContinuous, LMEmodel_curiosityContinuous)
}
# round results and save
LMEresults <- round(LMEresults, digits = 5)

setwd(memoryDir)
write.csv(LMEresults, paste0("LME_Results_curiosityByReward_", version, ".csv"))
if (pooled == 1){
  setwd(pooledDir)
  xlsx::write.xlsx(LMEresults, file="LME_Results_curiosityByReward.xlsx", sheetName = paste(version), append = T)
}


##### 3.2 LME with covariates predicting memory performance with mean-centered curiosity, effect-coded reward and their interaction #####

# define control vars
controlVar <- c("nback_accurary", "BAS_rewardresponsiveness", "ableToSee", "AvoidanceTemperament")
subjects  <- as.character(dfWide$fMRI)

setwd(memoryDir)

# loop controlling for different variables
for (cov in 1:length(controlVar)){
  # add cov to dfLong
  for (s in seq_along(subjects)){
    # calculate centered variable
    dfWide[[paste0(controlVar[cov], "_c")]] <- dfWide[[paste0(controlVar[cov])]] - mean(dfWide[[paste0(controlVar[cov])]])
    # take CENTERED cov for subject from dfWide and paste it to dfLong
    dfLong[dfLong$fMRI == subjects[s], paste0(controlVar[cov])] <- dfWide[dfWide$fMRI == subjects[s], paste0(controlVar[cov], "_c")]
  }
  
  # loop over dependent variables to compute LME
  for (DV in 1:length(DV_LME)){
    
    # define model based on scaling of variable
    if (max(dfLong[, DV_LME[DV]], na.rm = T) > 1){
      LMEmodel_curiosityContinuousControlled <- lmerTest::lmer(dfLong[, DV_LME[DV]] ~ get(controlVar[cov]) + groupEffectCoded*curiosityGroupMeanCentered + (1+curiosityGroupMeanCentered|ID) + (1|stimID), data = dfLong)
    }else{
      LMEmodel_curiosityContinuousControlled <- glmer(dfLong[, DV_LME[DV]] ~ get(controlVar[cov]) + groupEffectCoded*curiosityGroupMeanCentered + (1+curiosityGroupMeanCentered|ID) + (1|stimID), family = "binomial"(link = 'logit'), data = dfLong)
    }
    # put coefficients into data frame
    summaryCuriosityContinuousControlled <- summary(LMEmodel_curiosityContinuousControlled)
    summaryCuriosityContinuousCoefficientsControlled <- as.data.frame(summaryCuriosityContinuousControlled$coefficients)
    # change row names of data frame
    row.names(summaryCuriosityContinuousCoefficientsControlled) <- c(paste0("LME_", DV_LME[DV], "_Intercept"), paste0("LME_",DV_LME[DV], "_", controlVar[cov]),
                                                                     paste0("LME_", DV_LME[DV], "_groupEffectCoded"),                                                          paste0("LME_", DV_LME[DV], "_curiosityGroupMeanCentered"),
                                                                     paste0("LME_", DV_LME[DV], "_interaction"))
    # round results and save
    if (DV == 1) {
      LMEresultsControlled <- summaryCuriosityContinuousCoefficientsControlled
    }else {
      temp_LMEresults <- summaryCuriosityContinuousCoefficientsControlled
      LMEresultsControlled <- rbind.all.columns(LMEresultsControlled, temp_LMEresults) #rbind all columns will induce NA if there was initially no data saved in the loop per participant
      rm(temp_LMEresults)
    }
    rm(summaryCuriosityContinuousCoefficients, summaryCuriosityContinuous, LMEmodel_curiosityContinuousControlled)
  }
  # round results and save
  LMEresultsControlled <- round(LMEresultsControlled, digits = 5)
  write.csv(LMEresultsControlled, paste0("LME_Results_curiosityByReward_", version, "_controlledFor_", controlVar[cov], ".csv"))
}

# When there seem to be an interaction effect like nback, you can test a model like memory~nback*group to see if the interaction is indeed significant. If it is significant, there might be a way to show that memory performance is different for s particular subset of participants (e.g. those who are low in the nback task). Given that reward group showed numerically lower memory performance but it is still informative for us.

for (DV in 1:length(DV_LME)){
  # define model based on scaling of variable
  if (max(dfLong[, DV_LME[DV]], na.rm = T) > 1){ # if DV continuous
    LMEmodel_curiosityContinuous <- lmerTest::lmer(dfLong[, DV_LME[DV]] ~ groupEffectCoded*nback_accurary + (1|ID) + (1|stimID), data = dfLong)
  } else { # if DV binary
    LMEmodel_curiosityContinuous <- glmer(dfLong[, DV_LME[DV]] ~ groupEffectCoded*nback_accurary + (1|ID) + (1|stimID), family = "binomial"(link = 'logit'), data = dfLong)
  }
  # put coefficients into data frame
  summaryCuriosityContinuous <- summary(LMEmodel_curiosityContinuous)
  summaryCuriosityContinuousCoefficients <- as.data.frame(summaryCuriosityContinuous$coefficients)
  # change row names
  row.names(summaryCuriosityContinuousCoefficients) <- c(paste0("LME_",DV_LME[DV], "_Intercept"), paste0("LME_",DV_LME[DV], "_groupEffectCoded"), 
                                                         paste0("LME_",DV_LME[DV], "_nback_accurary"), 
                                                         paste0("LME_",DV_LME[DV], "_interaction"))
  # put coeffiencts from all DVs into one data frame
  if (DV == 1) {
    LMEresults <- summaryCuriosityContinuousCoefficients
  }else {
    temp_LMEresults <- summaryCuriosityContinuousCoefficients
    LMEresults <- rbind.all.columns(LMEresults, temp_LMEresults) #rbind all columns will induce NA if there was initially no data saved in the loop per participant
    rm(temp_LMEresults)
  }
  rm(summaryCuriosityContinuousCoefficients, summaryCuriosityContinuous, LMEmodel_curiosityContinuous)
}


##### 3.3 basic LME predicting memory performance using subset of reward group (i.e. those who believed in the manipulation) #####

# create subset of data
dfWide$reward_subset <- ifelse(dfWide$group == "cont", 1, #include everyone from control group
                               ifelse(dfWide$rewardBelief_score > 3 & dfWide$group == "exp", 1, 0)) # and those from reward group believing in manipulation

for (s in seq_along(subjects)){
  # take subset variable for subject from dfWide and paste it to dfLong
  dfLong[dfLong$fMRI == subjects[s], "reward_subset"] <- dfWide[dfWide$fMRI == subjects[s], "reward_subset"]
}

dfLong_subset <- subset(dfLong, dfLong$reward_subset == 1)
# loop over dependent variables to compute LME
for (DV in 1:length(DV_LME)){
  # define model based on scaling of variable
  if (max(dfLong_subset[, DV_LME[DV]], na.rm = T) > 1){ # if DV continuous
    LMEmodel_curiosityContinuous <- lmerTest::lmer(dfLong_subset[, DV_LME[DV]] ~ groupEffectCoded*curiosityGroupMeanCentered + (1+curiosityGroupMeanCentered|ID) + (1|stimID), data = dfLong_subset)
  } else { # if DV binary
    LMEmodel_curiosityContinuous <- glmer(dfLong_subset[, DV_LME[DV]] ~ groupEffectCoded*curiosityGroupMeanCentered + (1+curiosityGroupMeanCentered|ID) + (1|stimID), family = "binomial"(link = 'logit'), data = dfLong_subset)
  }
  # put coefficients into data frame
  summaryCuriosityContinuous <- summary(LMEmodel_curiosityContinuous)
  summaryCuriosityContinuousCoefficients <- as.data.frame(summaryCuriosityContinuous$coefficients)
  # change row names
  row.names(summaryCuriosityContinuousCoefficients) <- c(paste0("LME_subset_",DV_LME[DV], "_Intercept"), paste0("LME_subset_",DV_LME[DV], "_groupEffectCoded"), 
                                                         paste0("LME_subset_",DV_LME[DV], "_curiosityGroupMeanCentered"), 
                                                         paste0("LME_subset_",DV_LME[DV], "_interaction"))
  # put coeffiencts from all DVs into one data frame
  if (DV == 1) {
    LMEresults <- summaryCuriosityContinuousCoefficients
  }else {
    temp_LMEresults <- summaryCuriosityContinuousCoefficients
    LMEresults <- rbind.all.columns(LMEresults, temp_LMEresults) #rbind all columns will induce NA if there was initially no data saved in the loop per participant
    rm(temp_LMEresults)
  }
  rm(summaryCuriosityContinuousCoefficients, summaryCuriosityContinuous)
}
# round results and save
LMEresults <- round(LMEresults, digits = 5)
setwd(memoryDir)
write.csv(LMEresults, paste0("LME_Results_curiosityByReward_", version, "_subset.csv"))



##### 3.4  LME predicting memory performance with mean-centered curiosity only #####

# loop over dependent variables to compute LME
for (DV in 1:length(DV_LME)){
  # define model based on scaling of variable
  if (max(dfLong[, DV_LME[DV]], na.rm = T) > 1){ # if DV continuous
    LMEmodel_curiosityContinuous <- lmerTest::lmer(dfLong[, DV_LME[DV]] ~ curiosityGroupMeanCentered + (1+curiosityGroupMeanCentered|ID) + (1|stimID), data = dfLong)
  } else { # if DV binary
    LMEmodel_curiosityContinuous <- glmer(dfLong[, DV_LME[DV]] ~ curiosityGroupMeanCentered + (1+curiosityGroupMeanCentered|ID) + (1|stimID), family = "binomial"(link = 'logit'), data = dfLong)
  }
  # put coefficients into data frame
  summaryCuriosityContinuous <- summary(LMEmodel_curiosityContinuous)
  summaryCuriosityContinuousCoefficients <- as.data.frame(summaryCuriosityContinuous$coefficients)
  # change row names
  row.names(summaryCuriosityContinuousCoefficients) <- c(paste0("LME_",DV_LME[DV], "_Intercept"),
                                                         paste0("LME_",DV_LME[DV], "_curiosityGroupMeanCentered"))
  # put coeffiencts from all DVs into one data frame
  if (DV == 1) {
    LMEresults <- summaryCuriosityContinuousCoefficients
  }else {
    temp_LMEresults <- summaryCuriosityContinuousCoefficients
    LMEresults <- rbind.all.columns(LMEresults, temp_LMEresults) #rbind all columns will induce NA if there was initially no data saved in the loop per participant
    rm(temp_LMEresults)
  }
  rm(summaryCuriosityContinuousCoefficients, summaryCuriosityContinuous)
}
# round results and save
LMEresults <- round(LMEresults, digits = 5)
setwd(memoryDir)
write.csv(LMEresults, paste0("LME_Results_curiosityOnly_", version, "_subset.csv"))


########## 4. Create barplots to visualise the effects of curiosity and reward on memory performance ########## 
if (exists("dfLong") == F) {
  setwd(preprocessedDir)
  dfLong <- xlsx::read.xlsx("long_MAGMOT.xlsx", sheetName = "Sheet1")#,  showWarnings = FALSE)
} 

# create a dichomotised curiosity variable using mean-cenetred curiosity
dfLong$curiosity_dich <- ifelse(dfLong$curiosity_dich == -1, "below", 
                                ifelse(dfLong$curiosity_dich == 1, "above", NA)) # create curiosity_dichotom as factor
# define grouping variables
groupingVariables <- c("mediansplitCuriosity_MAGMOT", "mediansplitCuriosityWithinSubject", "curiosity_dich")
groupingVariables <- c("mediansplitCuriosityWithinSubject", "curiosity_dich")

# determine dependent variables to loop over (same as above)
DV_barplot <- c("cuedRecallStrict", "cuedRecallLenient", 
                "recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
                "rememberedStrictAboveAvg", "rememberedLenientAboveAvg", "rememberedStrictHigh", "rememberedLenientHigh",
                "confidence", "confidenceCorrectTrials", 
                "confidenceGroupMeanCentered", "confidenceGroupMeanCenteredCorrectTrials")

# determine colours to be used in histogram
cols <- c("above" = "#F8766D", "below" = "#00BFC4", "all tricks" = "grey", "median" = "#C77CFF")

# set directory where to save plots
setwd(memoryDir)

# loop over dependent variables to create bar plots
for (DV in DV_barplot){
  
  # loop over grouping variables for continuous curiosity variable
  for (g in 1:length(groupingVariables)){
    # create a data frame containing the data for each group divided regarding curiosity ratings
    outputGroup <- summarySEwithin(dfLong, measurevar=DV, betweenvars="group", withinvars=groupingVariables[g], idvar="ID", na.rm = T)
    levels(outputGroup$group) <- c("No reward", "Reward")
    names(outputGroup)[names(outputGroup) == groupingVariables[g]] <- "cutoff"
    
    # remove NAs (necessary for curiosity_dichotomous)
    outputGroup <- na.omit(outputGroup) 
    
    # add grouping variable as column
    outputGroup$groupingvar <- groupingVariables[g]
    
    # create a data frame containing the data for each group regardless of curiosity ratings
    outputAll <- summarySE(dfLong, measurevar=DV, groupvars="group", na.rm=T,
                           conf.interval=.95, .drop=TRUE)
    levels(outputAll$group) <-  c("No reward", "Reward")
    outputAll$cutoff <- rep("all tricks", 2)    
    outputAll$groupingvar <- groupingVariables[g]
    
    # combine these two data frames and use it for plotting
    if (g == 1){
      output <- rbind.all.columns(outputGroup, outputAll)
    } else {
      outputGroup <- rbind.all.columns(outputGroup, outputAll)
      output <- rbind.all.columns(output, outputGroup)
    }
  }
  output <- output[!output$cutoff =="all tricks", ]  # this line can be commented out to also illustrate effect regardless of group
  
  
  # create bar graph
  graph <- ggplot(output, aes(x=group, y=get(DV), fill=cutoff)) +
    geom_bar(stat="identity", position="dodge") + geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=get(DV)-se, ymax=get(DV)+se)) +
    scale_x_discrete(limits=c("No reward", "Reward")) + 
    labs(x="Between Group manipulation", y="Performance index", fill = "Curiosity category", title = paste(version, ":", DV ))  +
    theme(axis.text=element_text(size=20), axis.title=element_text(size=20, face="bold"), title=element_text(size =20, face="bold"), legend.title = element_text(size=20), legend.text = element_text(size = 20)) +
    theme_classic() + scale_fill_manual(values = cols) +
    facet_grid(. ~ groupingvar)     
  
  # modify bar graph depending on dependent variable
  if (DV %in% c("cuedRecallStrict", "cuedRecallLenient", "rememberedStrictAboveAvg", "rememberedLenientAboveAvg", "rememberedStrictHigh", "rememberedLenientHigh")){
    graph <- graph + coord_cartesian(ylim = c(0, 1))
  }
  if (DV %in% c("recognition", "recognitionAboveMeanConf", "recognitionConfLevel_4_5_6")){
    graph <- graph + coord_cartesian(ylim = c(0, 1)) + geom_hline(yintercept = 0.25, linetype="dashed", color = "black")
  }
  if (DV %in% c("confidence", "confidenceCorrectTrials")){
    graph <- graph + coord_cartesian(ylim = c(0, 6))
  }
  if (DV %in% c("confidenceGroupMeanCentered", "confidenceGroupMeanCenteredCorrectTrials")){
    graph <- graph + coord_cartesian(ylim = c(-1, 1))
  }
  print(graph)
  print(paste0("Bargraph_", DV, ".jpeg"))
  ggsave(paste0("Bargraph_", DV, ".jpeg"))
}



########## 5. Create histograms to further investigate the relation between reward, curiosity and memory ########## 

# define variables
DV_hist <- c("cuedRecallStrict", "cuedRecallLenient", 
             "recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
             "rememberedStrictAboveAvg", "rememberedLenientAboveAvg", "rememberedStrictHigh", "rememberedLenientHigh")
# create data frame and recode memory as factors
output <- dfLong[,c("ID", "group", "curiosityGroupMeanCentered", DV_hist)]
output$group <- ifelse(output$group == "cont", "No reward", "Reward")


# loop over dependent variables to create histograms 
for (DV in DV_hist){
  
  # recode memory performance
  output[[paste0(DV)]] <- ifelse(output[[paste0(DV)]] == 0, "forgotten", "remembered")
  
  # histogram plot: facet: group, col: memory
  graph <- ggplot(output, aes(x=curiosityGroupMeanCentered, col = get(DV), fill = get(DV))) + 
    geom_histogram(aes(y=..density..), binwidth = 0.1, alpha=0.2, position="dodge") +
    geom_density(alpha=.1) +
    theme_classic() + 
    scale_color_brewer(palette="Dark2", limits = c("remembered", "forgotten")) + scale_fill_brewer(palette="Dark2", limits = c("remembered", "forgotten")) +
    facet_grid(group ~ .) +
    coord_cartesian(xlim = c(-5, 5), ylim = c(0, 1)) +
    labs(x="curiosity group mean centered", y="Density", title = paste(version, DV)) +
    theme(legend.position="bottom") + theme(legend.title = element_blank()) + 
    theme(axis.text=element_text(size=14), axis.title=element_text(size=16, face="bold"), title=element_text(size =20, face="bold"), strip.text = element_text(size = 16)) 
  
  print(graph)
  print(paste0("Histogram_", DV, ".jpeg"))
  ggsave(paste0("Histogram_", DV, ".jpeg"))
  
  # histogram plot: facet: memory, col: group
  graph <- ggplot(output, aes(x=curiosityGroupMeanCentered, col = group, fill = group)) + 
    geom_histogram(aes(y=..density..), binwidth = 0.1, alpha=0.2, position="dodge") +
    geom_density(alpha=.1) +
    theme_classic() + 
    facet_grid(get(DV) ~ .) +
    coord_cartesian(xlim = c(-5, 5), ylim = c(0, 1)) +
    labs(x="curiosity group mean centered", y="Density", title = paste(version, DV)) +
    theme(legend.position="bottom") + theme(legend.title = element_blank()) + 
    theme(axis.text=element_text(size=14), axis.title=element_text(size=16, face="bold"), title=element_text(size =20, face="bold"), strip.text = element_text(size = 16)) 
  
  print(graph)
  print(paste0("Histogram_", DV, "_2.jpeg"))
  ggsave(paste0("Histogram_", DV, "_2.jpeg"))
}


#####################################################################################################################################
############################################### ANALYSIS BASED ON DATA IN WIDE FORMAT ############################################### 
#####################################################################################################################################

########## 1. Look at demographics of the sample ########## 

# age
output <- psych::describeBy(dfWide[,"age"], group=dfWide$group)
output <- as.data.frame(rbind(output$cont, output$exp))
output$mot <- rep(c("cont","exp"), each = 1)

outg <- ggplot(output, aes(mot, mean, fill = mot))
outg + geom_bar(stat="identity", position="dodge") + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9))  + scale_x_discrete(limits=c("exp","cont")) + labs(x="Experimental condition", y="Age", fill = "Experimental Condition", title = paste("demogs I", version)) + theme_classic() + scale_fill_discrete(guide=FALSE)

t.test(dfWide$age ~ dfWide$group)

# gender
rm(output)
output <- plyr::count(dfWide, vars = c("gender","group"))

# output <- count(df, c('gender','cond'))
outg <- ggplot(output, aes(group, freq, fill = gender))
outg + geom_bar(stat="identity", position="fill")+ scale_x_discrete(limits=c("cont","exp")) + labs(x="Experimental condition", y="Frequency", fill = "Gender", title = paste("demogs II", version)) + theme_classic() 


########## 2. Look at the effectiveness of reward manipulation ########## 

# check whether responses ppt gave are plausible
dfWide$fMRI[!is.na(dfWide$rewardBelief_score)] # missing value for MAGMOT_06, value for MAGMOT_35
dfWide$rewardBelief_score[dfWide$fMRI == "MAGMOT_35"] <- NA
# manipulation check
t.test(dfWide$rewardBelief_score, alternative = "greater", mu = 3)

dfWide$rewardEffort
cor_rewardEffort <- as.data.frame(t(cor(dfWide$rewardEffort, dfWide[,c("BIS", "BAS_rewardresponsiveness", "BAS_drive", "BAS_funseeking", 
                                                                       "NeedForCognition", "FearOfFailure", "ApproachTemperament", "AvoidanceTemperament", 
                                                                       "TraitCuriosity", "StateCuriosity", 
                                                                       "intrinsicMotivation", "taskEngagement", "interest", "boredom", "effort", "pressure", 
                                                                       "compliance", "ableToSee", 
                                                                       "rewardBelief_score",
                                                                       "corsiSpan", "nback_hitrate", "nback_falsealarmrate", "nback_accurary",
                                                                       "cuedRecallStrict_abs", "cuedRecallLenient_abs",
                                                                       "allConf_abs", "highConf_abs", "aboveAvgConf_abs",
                                                                       "rememberedStrictAboveAvg_abs", "rememberedLenientAboveAvg_abs", "rememberedStrictHigh_abs", "rememberedLenientHigh_abs",
                                                                       "meanConfidence", "meanConfidenceCorrectTrials")], 
                                        use = "pairwise.complete.obs")))
names(cor_rewardEffort) <- "cor_rewardEffort"



cor_rewardBelief <- as.data.frame(t(cor(dfWide$rewardBelief_score, dfWide[,c("BIS", "BAS_rewardresponsiveness", "BAS_drive", "BAS_funseeking", 
                                                                             "NeedForCognition", "FearOfFailure", "ApproachTemperament", "AvoidanceTemperament", 
                                                                             "TraitCuriosity", "StateCuriosity", 
                                                                             "intrinsicMotivation", "taskEngagement", "interest", "boredom", "effort", "pressure", 
                                                                             "compliance", "ableToSee", 
                                                                             "rewardBelief_score",
                                                                             "corsiSpan", "nback_hitrate", "nback_falsealarmrate", "nback_accurary",
                                                                             "cuedRecallStrict_abs", "cuedRecallLenient_abs",
                                                                             "allConf_abs", "highConf_abs", "aboveAvgConf_abs",
                                                                             "rememberedStrictAboveAvg_abs", "rememberedLenientAboveAvg_abs", "rememberedStrictHigh_abs", "rememberedLenientHigh_abs",
                                                                             "meanConfidence", "meanConfidenceCorrectTrials")], 
                                        use = "pairwise.complete.obs")))
names(cor_rewardBelief) <- "cor_rewardBelief"

merge(cor_rewardBelief, cor_rewardEffort, by = "row.names")

hist(dfWide$rewardBelief_score, xlim = c(1,6), main="Did you believe that you would receive a bonus payment based on your performance?",
     xlab = "Level of agreement (1 = diagree, 6 = agree)")

hist(dfWide$rewardEffort, xlim = c(1,7))

########## 3. Create bar plots for  mean values (with SE) of the questionnaire scores for each group ########## 

setwd(ratingsDir)
rm(output)

# define scales
scales <- names(dfWide[,c("BIS", "BAS_rewardresponsiveness", "BAS_drive", "BAS_funseeking", 
                          "NeedForCognition", "FearOfFailure", "ApproachTemperament", "AvoidanceTemperament", 
                          "TraitCuriosity", "StateCuriosity", 
                          "intrinsicMotivation", "taskEngagement", "interest", "boredom", "effort", "pressure", 
                          "compliance", "ableToSee", 
                          "rewardBelief_score",
                          "corsiSpan", "nback_hitrate", "nback_falsealarmrate", "nback_accurary")])
psych::describe(dfWide[,scales])
by(cbind(dfWide[,scales]), dfWide$group, psych::describe)

# combine different measurements / scales
IMI <- names(dfWide[,c("intrinsicMotivation", "taskEngagement", "interest", "boredom", "effort", "pressure")])
BISBAS <- names(dfWide[,c("BIS", "BAS_rewardresponsiveness", "BAS_drive", "BAS_funseeking")])
MCI <- names(dfWide[,c("TraitCuriosity", "StateCuriosity")])
others <- names(dfWide[,c("NeedForCognition", "FearOfFailure", "ApproachTemperament", "AvoidanceTemperament")])
postQest <- names(dfWide[,c("compliance", "ableToSee", "memoryTestKnown", "memoryIntention", "rewardBelief_score", "rewardEffort")])
neuro <- names(dfWide[,c("corsiSpan", "nback_hitrate", "nback_falsealarmrate", "nback_accurary")])

allQ <- c("IMI", "BISBAS", "MCI", "others", "neuro", "postQest")

setwd(questDir)
# plot average responses for each scale seperately for each group
for (q in seq_along(allQ)){
  # create data frame used by ggplot
  output <- by(cbind(dfWide[,get(allQ[q])]), dfWide$group, psych::describe)
  rating <- as.data.frame(rbind(output$cont, output$exp))
  rating$mot <- rep(c("intrinsic","extrinsic"), each = length(get(allQ[q])))
  
  # bar plot
  outg <- ggplot(rating, aes(vars, mean, fill = mot))
  outg <- outg + geom_bar(stat="identity", position="dodge") + 
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9))  + 
    labs(x="Scales", y="Rating", fill = "Experimental Condition", title = paste(allQ[q], "questionnaire",version)) + 
    theme_classic() + #coord_cartesian(ylim = c(0, 7)) +
    scale_x_discrete(limits=get(allQ[q])) +
    theme(legend.position="bottom") 
  print(outg)
  print(paste0("QuestionnaireDataByGroup_", allQ[q],".jpeg")) # save plot
  ggsave(paste0("QuestionnaireDataByGroup_", allQ[q],".jpeg")) # save plot
}

rm(output, rating, outg)

########## 4. Compute t-tests and effect sizes for between-group differences in questionnaire scores ########## 
scales <- names(dfWide[,c("intrinsicMotivation", "taskEngagement", "interest", "boredom", "effort", "pressure", 
                          "BIS", "BAS_rewardresponsiveness", "BAS_drive", "BAS_funseeking", 
                          "NeedForCognition", "FearOfFailure", "ApproachTemperament", "AvoidanceTemperament", 
                          "TraitCuriosity", "StateCuriosity", 
                          "compliance", "ableToSee", 
                          "corsiSpan", "nback_hitrate", "nback_falsealarmrate", "nback_accurary")])

for(scale in 1:length(scales)) {
  print(scales[scale])
  # compute t-test for group difference
  ttest <- t.test(dfWide[,scales[scale]]~dfWide$group)
  t.stats <- as.data.frame(t(round(c(ttest$statistic, ttest$p.value), digits = 3)))
  names(t.stats) <- c("tValue", "pValue(t)")
  
  # compute wilcox
  wilcox <- wilcox.test(dfWide[,scales[scale]]~dfWide$group) 
  attributes(wilcox)
  w.stats <- as.data.frame(t(round(c(wilcox$statistic, wilcox$p.value), digits = 3)))
  names(w.stats) <- c("W", "pValue(W)")
  
  # merge t-test and wilconxon's test
  t.stats <- merge(t.stats, w.stats)
  
  # compute mean for each group
  means <- tapply(dfWide[,scales[scale]], dfWide$group, mean, na.rm = T)
  means <- as.data.frame(t(means))
  means <- merge(t.stats, means)
  
  # compute effect size
  if (means$exp != means$cont) {
    data <- dfWide[,c("group", scales[scale])]
    psych::cohen.d(data, "group")
    d <- psych::cohen.d(data, "group")
    cohen <- as.data.frame(d$cohen.d)
    means <- merge(means, cohen)
  }
  
  # put all in a data frame
  row.names(means) <- paste(scales[scale])
  if (scale == 1) {
    effectsizesScales <- means
  } else {
    temp_effectsizesScales <- means
    effectsizesScales <- rbind.all.columns(effectsizesScales, temp_effectsizesScales) #rbind all columns will induce NA if there was initially no data saved in the loop per participant
    rm(temp_effectsizesScales)
  }
}
rm(data, cohen, means, t.stats, ttest, d, w.stats, wilcox)
setwd(ratingsDir)
write.csv(effectsizesScales, paste0("effectsizesScales_", version, ".csv"))
if (pooled == 1){
  setwd(pooledDir)
  xlsx::write.xlsx(effectsizesScales, file="effectsizesScales.xlsx", sheetName = paste(version), append = T)
}


########## 5. Check whether any of the differences in the questionnaires relates to memory performance and whether there is a group effect  ########## 

# create intercorrelation table

cor_vars <- c("BIS", "BAS_rewardresponsiveness", "BAS_drive", "BAS_funseeking", 
              "NeedForCognition", "FearOfFailure", "ApproachTemperament", "AvoidanceTemperament", 
              "TraitCuriosity", "StateCuriosity", 
              "intrinsicMotivation", "taskEngagement", "interest", "boredom", "effort", "pressure", 
              "compliance", "ableToSee", "rewardBelief_score",
              "corsiSpan", "nback_hitrate", "nback_falsealarmrate", "nback_accurary",
              "cuedRecallStrict_abs", "cuedRecallLenient_abs",
              "allConf_abs", "highConf_abs", "aboveAvgConf_abs",
              "rememberedStrictAboveAvg_abs", "rememberedLenientAboveAvg_abs", "rememberedStrictHigh_abs", "rememberedLenientHigh_abs",
              "meanConfidence", "meanConfidenceCorrectTrials")

cor_all <- as.data.frame(t(cor(dfWide[,cor_vars], use = "pairwise.complete.obs")))

cor_cont <- as.data.frame(t(cor(dfWide[dfWide$group == "cont", cor_vars],  use = "pairwise.complete.obs")))

cor_exp <- as.data.frame(t(cor(dfWide[dfWide$group == "exp", cor_vars],  use = "pairwise.complete.obs")))

cor_all <- round(cor_all, digits = 3)
cor_cont <- round(cor_cont, digits = 3)
cor_exp <- round(cor_exp, digits = 3)

setwd(ratingsDir)
xlsx::write.xlsx(cor_all, file=paste0("Intercorrelation_", version, ".xlsx"), sheetName = "cor_all", row.names = F) 
xlsx::write.xlsx(cor_cont, file=paste0("Intercorrelation_", version, ".xlsx"), sheetName = "cor_cont", row.names = F, append = T) 
xlsx::write.xlsx(cor_exp, file=paste0("Intercorrelation_", version, ".xlsx"), sheetName = "cor_exp", row.names = F, append = T) 


########## 6. Compute summary statistics for measurements of memory (subject sum scores) ########## 

recognitionPerformanceVars <- c("recognitionConfLevel_1", "recognitionConfLevel_above_1", "recognitionConfLevel_1_2", "recognitionConfLevel_1_2_3", "recognitionConfLevel_2",
                                "recognitionConfLevel_above_2", "recognitionConfLevel_3", "recognitionConfLevel_above_3", "recognitionConfLevel_3_4", "recognitionConfLevel_4", "recognitionConfLevel_above_4",
                                "recognitionConfLevel_5", "recognitionConfLevel_above_5", "recognitionConfLevel_5_6", "recognitionConfLevel_6",
                                "meanConfidence", "meanConfidenceCorrectTrials")

recognitionPerformanceVars <- c(recognitionPerformanceVars, paste0(memoryLabels, "_abs"))
recognitionPerformanceVars <- c(recognitionPerformanceVars, paste0(memoryLabels, "_rel"))
recognitionPerformanceVars <- c(recognitionPerformanceVars, paste0("curBen_cont_",memoryLabels))
recognitionPerformanceVars <- c(recognitionPerformanceVars, paste0("curBen_dich_",memoryLabels))
recognitionPerformanceVars <- c(recognitionPerformanceVars, paste0("curBen_rel_",memoryLabels))
recognitionPerformanceVars <- c(recognitionPerformanceVars, paste0("curCor_",memoryLabels))
recognitionPerformanceVars <- c(recognitionPerformanceVars, paste0("curBeta_",memoryLabels))

# create a data frame that shows the summary statistics for each score, for the whole population
dataWideRecognitionPerformance <- dfWide[,recognitionPerformanceVars]
df_mean <- data.frame(apply(dataWideRecognitionPerformance, 2, mean, na.rm = T), recognitionPerformanceVars)
df_sd <- data.frame(apply(dataWideRecognitionPerformance, 2, sd, na.rm = T), recognitionPerformanceVars)
df_min <- data.frame(apply(dataWideRecognitionPerformance, 2, min, na.rm = T), recognitionPerformanceVars)
df_max <- data.frame(apply(dataWideRecognitionPerformance, 2, max, na.rm = T), recognitionPerformanceVars)

descriptivesRecognitionPerformance <- merge(df_mean, df_sd, by = "recognitionPerformanceVars")
descriptivesRecognitionPerformance <- merge(descriptivesRecognitionPerformance, df_min, by = "recognitionPerformanceVars")
descriptivesRecognitionPerformance <- merge(descriptivesRecognitionPerformance, df_max, by = "recognitionPerformanceVars")
rm(df_mean, df_sd, df_min, df_max)

names(descriptivesRecognitionPerformance) <- c("recognitionPerformanceVar", "mean", "sd", "min", "max")

setwd(memoryDir)
write.csv(descriptivesRecognitionPerformance, paste0("subjectSumscore_performance_memory_", version, ".csv"), row.names = F)
if (pooled == 1){
  setwd(pooledDir)
  xlsx::write.xlsx(descriptivesRecognitionPerformance, file="subjectSumscore_performance_memory.xlsx", sheetName = paste(version), row.names = F, append = T)
}
rm(dataWideRecognitionPerformance, descriptivesRecognitionPerformance)


########## 7. Compute t-tests and effect sizes for between-group differences in memory scores ########## 


# define dependent variables (i.e. sum scores)

DV_wide <- c(paste0(memoryLabels, "_abs")) # absolute sum scores
DV_wide <- c(DV_wide, paste0("curBeta_",memoryLabels)) # betas
DV_wide <- c(DV_wide, paste0("curCor_",memoryLabels)) # correlation
for (mem in 1:length(memoryLevels)) { # benefits
  DV_wide <- c(DV_wide, paste0("curBen_cont_", memoryLabels[mem], collapse = ", "))
  DV_wide <- c(DV_wide, paste0("curBen_dich_", memoryLabels[mem], collapse = ", "))
  DV_wide <- c(DV_wide, paste0("curBen_rel_", memoryLabels[mem], collapse = ", "))
}
DV_wide <- c(DV_wide, paste("RSFC_VTAHPC_diff"),  paste("RSFC_VTAHPC_diff_spearman")) # brain data


# for all dependent variables compute two-sample t-test and calculate effect size
for(DV in 1:length(DV_wide)) {
  
  
  print(DV_wide[DV])
  # compute t-test for group difference
  ttest <- t.test(dfWide[,DV_wide[DV]]~dfWide$group)
  t.stats <- as.data.frame(t(round(c(ttest$statistic, ttest$p.value), digits = 3)))
  names(t.stats) <- c("tValue", "pValue(t)")
  
  # compute wilcox
  wilcox <- wilcox.test(dfWide[,DV_wide[DV]]~dfWide$group) 
  attributes(wilcox)
  w.stats <- as.data.frame(t(round(c(wilcox$statistic, wilcox$p.value), digits = 3)))
  names(w.stats) <- c("W", "pValue(W)")
  
  # merge t-test and wilconxon's test
  t.stats <- merge(t.stats, w.stats)
  
  # compute mean for each group
  means <- tapply(dfWide[,DV_wide[DV]], dfWide$group, mean, na.rm = T)
  means <- as.data.frame(t(means))
  means <- merge(t.stats, means)
  
  # compute effect size
  if (means$cont != means$exp) {
    data <- dfWide[,c("group", DV_wide[DV])]
    d <- psych::cohen.d(data, "group")
    cohen <- as.data.frame(d$cohen.d)
    means <- merge(means, cohen)
  }
  
  # put all in a data frame
  row.names(means) <- paste(DV_wide[DV])
  if (DV == 1) {
    effectsizesMemory <- means
  } else {
    temp_effectsizesMemory <- means
    effectsizesMemory <- rbind.all.columns(effectsizesMemory, temp_effectsizesMemory) #rbind all columns will induce NA if there was initially no data saved in the loop per participant
    rm(temp_effectsizesMemory)
  }
}
# delete unnecssary variables, round values and save dataframe
rm(data, cohen, means, t.stats, ttest, d, w.stats, wilcox)
effectsizesMemory <- round(effectsizesMemory, digits = 3)
setwd(memoryDir)
write.csv(effectsizesMemory, paste0("effectsizesMemory_", version, ".csv"))
if (pooled == 1){
  setwd(pooledDir)
  xlsx::write.xlsx(effectsizesMemory, file="effectsizesMemory.xlsx", sheetName = paste(version), append = T)
}
effectsizes <- rbind.all.columns(effectsizesScales, effectsizesMemory)
setwd(ratingsDir)
write.csv(effectsizes, paste0("effectsizesAll_", version, ".csv"))

########## 8. Create violin plots for sum scores of memory measures in each group ########## 
# subset dataWide
output <- dfWide[,c("ID", "group", DV_wide)]
output$group <- ifelse(output$group == "cont", "No reward", "Reward")
setwd(memoryDir)

threshold_remembered1 <- c("rememberedStrictAboveAvg_abs", "curBeta_rememberedStrictAboveAvg", "curCor_rememberedStrictAboveAvg",
                           "curBen_cont_rememberedStrictAboveAvg", "curBen_dich_rememberedStrictAboveAvg", "curBen_rel_rememberedStrictAboveAvg")
threshold_remembered2 <- c("rememberedLenientAboveAvg_abs", "curBeta_rememberedLenientAboveAvg", "curCor_rememberedLenientAboveAvg", 
                           "curBen_cont_rememberedLenientAboveAvg", "curBen_dich_rememberedLenientAboveAvg", "curBen_rel_rememberedLenientAboveAvg")
threshold_remembered3 <- c("rememberedStrictHigh_abs", "curBeta_rememberedStrictHigh", "curCor_rememberedStrictHigh",
                           "curBen_cont_rememberedStrictHigh", "curBen_dich_rememberedStrictHigh", "curBen_rel_rememberedStrictHigh")
threshold_remembered4 <- c("rememberedLenientHigh_abs", "curBeta_rememberedLenientHigh", "curCor_rememberedLenientHigh",
                           "curBen_cont_rememberedLenientHigh", "curBen_dich_rememberedLenientHigh", "curBen_rel_rememberedLenientHigh")

# define variables to plot data for
plotVars <- c("recall", "recollection", # absolute values 
              "betaRecall", "betaRecollection", # beta values
              "absoluteBenefitRecall_dich", "absoluteBenefitRecollection_dich", # dichotomised curiosity benefit
              "absoluteBenefitRecall_cont", "absoluteBenefitRecollection_cont", # continuous curiosity benefit
              "relativeBenefitRecall", "relativeBenefitRecollection", # relative curiosity benefit
              "correlationRecall", "correlationRecollection", # correlation
              "recognition", "betaRecognition", "correlationRecognition",
              "absoluteBenefitRecognition_dich", "absoluteBenefitRecognition_cont", "relativeBenefitRecognition",
              "remembered", "betaRemembered", "correlationRemembered",
              "absoluteBenefitRemembered_dich", "absoluteBenefitRemembered_cont", "relativeBenefitRemembered")

# define thresholds for the facet.grid "criteria"
threshold_recoll1 <- c("cuedRecallStrict_abs", "aboveAvgConf_abs", 
                       "curBeta_cuedRecallStrict", "curBeta_aboveAvgConf", 
                       "curBen_dich_cuedRecallStrict", "curBen_dich_aboveAvgConf",
                       "curBen_cont_cuedRecallStrict", "curBen_cont_highConf", 
                       "curBen_rel_cuedRecallStrict", "curBen_rel_aboveAvgConf", 
                       "curCor_cuedRecallStrict", "curCor_aboveAvgConf")

threshold_recoll2 <- c("cuedRecallLenient_abs", "highConf_abs", 
                       "curBeta_cuedRecallLenient", "curBeta_highConf", 
                       "curBen_dich_cuedRecallLenient", "curBen_dich_highConf",
                       "curBen_cont_cuedRecallLenient", "curBen_cont_aboveAvgConf",
                       "curBen_rel_cuedRecallLenient", "curBen_rel_highConf", 
                       "curCor_cuedRecallLenient", "curCor_highConf")

# define thresholds for unthresholded recognition
threshold_recog <- c("allConf_abs", "curBeta_allConf", "curCor_allConf", "curBen_dich_allConf", "curBen_cont_allConf", "curBen_rel_allConf")

# define thresholds and labels for facet.grid "criteria" * "method"
threshold_remembered1 <- c("rememberedStrictAboveAvg_abs", "curBeta_rememberedStrictAboveAvg", "curCor_rememberedStrictAboveAvg",
                           "curBen_cont_rememberedStrictAboveAvg", "curBen_dich_rememberedStrictAboveAvg", "curBen_rel_rememberedStrictAboveAvg")
threshold_remembered2 <- c("rememberedLenientAboveAvg_abs", "curBeta_rememberedLenientAboveAvg", "curCor_rememberedLenientAboveAvg", 
                           "curBen_cont_rememberedLenientAboveAvg", "curBen_dich_rememberedLenientAboveAvg", "curBen_rel_rememberedLenientAboveAvg")
threshold_remembered3 <- c("rememberedStrictHigh_abs", "curBeta_rememberedStrictHigh", "curCor_rememberedStrictHigh",
                           "curBen_cont_rememberedStrictHigh", "curBen_dich_rememberedStrictHigh", "curBen_rel_rememberedStrictHigh")
threshold_remembered4 <- c("rememberedLenientHigh_abs", "curBeta_rememberedLenientHigh", "curCor_rememberedLenientHigh",
                           "curBen_cont_rememberedLenientHigh", "curBen_dich_rememberedLenientHigh", "curBen_rel_rememberedLenientHigh")

# set variables to loop through
pp <- 0
ppp <- 0

# loop over all variables to plot
for(p in 1:length(plotVars)) {
  
  plot <- plotVars[p]
  print(plot)
  
  if (plot %in% c("recall", "recollection", # absolute values 
                  "betaRecall", "betaRecollection", # beta values
                  "absoluteBenefitRecall_dich", "absoluteBenefitRecollection_dich", # dichotomised curiosity benefit
                  "absoluteBenefitRecall_cont", "absoluteBenefitRecollection_cont", # continuous curiosity benefit
                  "relativeBenefitRecall", "relativeBenefitRecollection", # relative curiosity benefit
                  "correlationRecall", "correlationRecollection")){ # correlation
    # create data frame in long format
    columns <- c("ID", "group", threshold_recoll1[p], threshold_recoll2[p])
    output_plot <- output[,columns]
    output_plot <- reshape2::melt(output_plot, id=c("ID","group"))
    names(output_plot) <- c("ID", "group", "criteria", "performance")
    
  } else if (plot %in% c("recognition", "betaRecognition", "relativeBenefitRecognition", "correlationRecognition",
                         "absoluteBenefitRecognition_dich", "absoluteBenefitRecognition_cont", "relativeBenefitRecognition")){
    # create data frame in long format
    pp <- pp+1
    columns <- c("ID", "group", threshold_recog[pp])
    output_plot <- output[,columns]
    output_plot <- reshape2::melt(output_plot, id=c("ID","group"))
    names(output_plot) <- c("ID", "group", "criteria", "performance")
    
  } else if (plot %in% c("remembered", "betaRemembered", "correlationRemembered",
                         "absoluteBenefitRemembered_dich", "absoluteBenefitRemembered_cont", "relativeBenefitRemembered")){
    # create data frame in long format
    ppp <- ppp+1
    columns <- c("ID", "group", threshold_remembered1[ppp], threshold_remembered2[ppp], threshold_remembered3[ppp], threshold_remembered4[ppp])
    output_plot <- output[,columns]
    output_plot <- reshape2::melt(output_plot, id=c("ID","group"))
    names(output_plot) <- c("ID", "group", "criteria", "performance")
    output_plot$method <- ifelse(output_plot$criteria ==  threshold_remembered1[ppp] | output_plot$criteria ==  threshold_remembered3[ppp], "strict recall", "lenient recall")
    output_plot$criteria <- ifelse(output_plot$criteria == threshold_remembered1[ppp] | output_plot$criteria == threshold_remembered2[ppp], "above average conf", "high conf")
  } 
  
  # Basic violin plot
  # graph <- ggplot(get(paste0("output_",plot)), aes(x=group, y=performance, fill = group)) + 
  graph <- ggplot(get(paste0("output_plot")), aes(x=group, y=performance, fill = group)) + 
    geom_violin(trim=FALSE) +
    geom_boxplot(width=0.1) +
    geom_jitter(size = 1, shape=1, position=position_jitter(0.2)) +
    theme_classic() +
    labs(x="Experimental Condition", y="Sum score", title = paste(version, plot)) +
    theme(legend.position="none") +
    theme(axis.text=element_text(size=20), axis.title=element_text(size=20, face="bold"), title=element_text(size =20, face="bold")) 
  # change facet grid depending on dependent variable
  if (plot %in% c("recall", "recollection", # absolute values 
                  "betaRecall", "betaRecollection", # beta values
                  "absoluteBenefitRecall_dich", "absoluteBenefitRecollection_dich", # dichotomised curiosity benefit
                  "absoluteBenefitRecall_cont", "absoluteBenefitRecollection_cont", # continuous curiosity benefit
                  "relativeBenefitRecall", "relativeBenefitRecollection", # relative curiosity benefit
                  "correlationRecall", "correlationRecollection")){
    graph <- graph + facet_grid(. ~ criteria) 
  } else if (plot %in% c( "remembered", "betaRemembered", "correlationRemembered",
                          "absoluteBenefitRemembered_dich", "absoluteBenefitRemembered_cont", "relativeBenefitRemembered")){
    graph <- graph + facet_grid(criteria ~ method)
  }
  # change y axis depending on variable
  if (plot %in% c("recall", "recognition", "recollection", "remembered")){
    graph <- graph + coord_cartesian(ylim = c(-5, 41))
  } else if (plot %in% c("absoluteBenefitRecall_dich", "absoluteBenefitRecollection_dich", # dichotomised curiosity benefit
                         "absoluteBenefitRecall_cont", "absoluteBenefitRecollection_cont", # continuous curiosity benefit
                         "absoluteBenefitRecognition_dich", "absoluteBenefitRecognition_cont",
                         "absoluteBenefitRemembered_dich", "absoluteBenefitRemembered_cont")){
    graph <- graph + coord_cartesian(ylim = c(-23, 23))
  } else if (plot %in% c("relativeBenefitRecall", "relativeBenefitRecollection",
                         "relativeBenefitRecognition", "relativeBenefitRemembered",
                         "correlationRecall", "correlationRecollection",
                         "correlationRecognition", "correlationRemembered")){
    graph <- graph + coord_cartesian(ylim = c(-1, 1))
  } else if (plot %in% c("betaRecall", "betaRecognition", "betaRemembered", "betaRecollection")){
    graph <- graph + coord_cartesian(ylim = c(-0.4, 0.4))
  }
  
  
  print(graph)
  print(paste0("Violonplot_", plot, ".jpeg"))
  ggsave(paste0("Violonplot_", plot, ".jpeg"))
  
}


########## 9. Look at the change in memory performance between blocks over time ########## 

blocks <- c(1,2,3,4)
blockstring <- c("_firstBlock", "_secondBlock", "_thirdBlock", "")


DV_wide <- c(paste0(memoryLabels, "_rel")) # absolute sum scores
DV_wide <- c(DV_wide, paste0("curCor_",memoryLabels)) # correlation
for (mem in 1:length(memoryLevels)) { # benefits
  DV_wide <- c(DV_wide, paste0("curBen_rel_", memoryLabels[mem], collapse = ", "))
}

# create a data frame that includes the mean + CI for each group and the effect size of the group difference + CI
# aim: data set with col: dependentVar; measure (effect size/group mean);  group (int / ext / effect); value; upper; lower
# rows: dependent variables


for (block in blocks){
  for(DV in 1:length(DV_wide)) {
    paste0(DV_wide[DV])
    # compute mean for each block and group
    means <- as.data.frame(tapply(dfWide[[paste0(DV_wide[DV], blockstring[block])]], dfWide$group, mean, na.rm = T)) # maybe not transpose
    names(means) <- "effect" # need to be the same as 
    means$grouping <- row.names(means)
    
    # compute SD for each block and grouping
    sds <- as.data.frame(tapply(dfWide[[paste0(DV_wide[DV], blockstring[block])]], dfWide$group, sd, na.rm = T)) # maybe not transpose
    names(sds) <- "sd" # need to be the same as 
    sds$grouping <- row.names(sds)
    
    # compute SE
    sds$se[sds$grouping == "cont"] <- (sds$sd[sds$grouping == "cont"]) / (nrow(dfWide[dfWide$group=="cont",]))
    sds$se[sds$grouping == "exp"] <- (sds$sd[sds$grouping == "exp"]) / (nrow(dfWide[dfWide$group=="exp",]))
    
    # merge means and sds
    df <- merge(means, sds, by = "grouping")
    
    # compute confidence interval
    df$upper <- df$effect + 1.96 * df$se
    df$lower <- df$effect - 1.96 * df$se
    
    # calculate effect size
    if (df$effect[df$grouping=="cont"] != df$effect[df$grouping=="exp"]) {
      data <- dfWide[,c("group", paste0(DV_wide[DV], blockstring[block]))]
      d <- psych::cohen.d(data, "group")
      cohen <- as.data.frame(d$cohen.d)
      cohen$grouping <- "effect"
      df <- rbind.all.columns(df, cohen)
      rm(d)
    } else {
      cohen <- data.frame("lower" = 0, "effect" = 0, "upper" = 0, "grouping" = "effect")
      df <- rbind.all.columns(df, cohen)
    }
    rm(cohen)
    
    # combine data from all DVs
    df$dependentVar <- paste(DV_wide[DV])
    
    if (DV == 1) {
      effectsizesMemoryBlock <- df
    } else {
      temp_effectsizesMemory <- df
      effectsizesMemoryBlock <- rbind.all.columns(effectsizesMemoryBlock, temp_effectsizesMemory) #rbind all columns will induce NA if there was initially no data saved in the loop per participant
      rm(temp_effectsizesMemory)
    }
  }
  
  # for those occasions where mean in both groups was identical, add zeros:
  effectsizesMemoryBlock$lower[is.na(effectsizesMemoryBlock$lower)] <- 0
  effectsizesMemoryBlock$upper[is.na(effectsizesMemoryBlock$upper)] <- 0
  effectsizesMemoryBlock$effect[is.na(effectsizesMemoryBlock$effect)] <- 0
  
  # add block
  effectsizesMemoryBlock$blockNumber <- block
  
  # delete row.names
  row.names(effectsizesMemoryBlock) <- NULL
  
  # rbind variables
  if (block == 1){
    effectsizesMemoryBlock_long <- effectsizesMemoryBlock
  }  else {
    effectsizesMemoryBlock_temp <- effectsizesMemoryBlock
    effectsizesMemoryBlock_long <- rbind.all.columns(effectsizesMemoryBlock_long, effectsizesMemoryBlock_temp)
    rm(effectsizesMemoryBlock_temp)
  }
  #assign(paste0("effectsizesMemoryBlock",block), effectsizesMemoryBlock) # if effect sizes for a single block were of interest
  rm(effectsizesMemoryBlock)
}


# add measure to use as facet.grid
effectsizesMemoryBlock_long$measure <- c("group mean", "group mean", "effect size")

# add group discription
effectsizesMemoryBlock_long$grouping <- ifelse(effectsizesMemoryBlock_long$grouping == "cont", "Mean (no reward)",
                                               ifelse(effectsizesMemoryBlock_long$grouping == "exp", "Mean (reward)",
                                                      ifelse(effectsizesMemoryBlock_long$grouping == "effect", "Cohen's d", NA)))
# define colours
cols <- c("Mean (no reward)" = "#F8766D", "Mean (reward)" = "#00BFC4", "Cohen's d" = "black")

# define intercept for horizontal line
hline <- data.frame(measure = c("group mean", "effect size"), intercept = c(NA, 0))

# for each variable in list, create a plot showing the average performance in each group per block as well as the effect size
for(DV in 1:length(DV_wide)) {
  
  # subset the data frame containing the effect sizes per block and select relevant DV only
  output <- subset(effectsizesMemoryBlock_long, effectsizesMemoryBlock_long$dependentVar == DV_wide[DV], col = grouping)
  
  # create graph
  graph <- ggplot(data=output, aes(x=blockNumber, y=effect)) +
    geom_point(aes(col = grouping)) + 
    geom_errorbar(aes(ymin=lower, ymax=upper, col = grouping), width=.025) +
    facet_grid(measure ~., scales = "free_y") +
    theme_classic() + scale_color_manual(values = cols) +
    scale_x_discrete(limits=c("1", "2", "3", "all")) +
    labs(x="Task block", y="Value", group = "", title = paste(version, "Proportial effect over time",DV_wide[DV])) +
    theme(legend.title = element_blank()) +
    theme(axis.text=element_text(size=20), axis.title=element_text(size=20, face="bold"), title=element_text(size =20, face="bold")) 
  
  # print(graph)
  # add dashed line for zero effect size
  graph <-graph + geom_hline(aes(yintercept = intercept), hline, linetype="dashed", color = "grey")
  
  print(graph)  
  print(paste0("Graph_changeOfEffectsizeOverBlocks_", DV_wide[DV], ".jpeg"))
  # save file
  ggsave(file = paste0("Graph_changeOfEffectsizeOverBlocks_", DV_wide[DV], ".jpeg"), graph)
  
}


########## 10. Test changes if RSFC between HPC and VTA for significance ########## 
# t-test for groupwise differences in RSFC change
t.test(dfWide$RSFC_VTAHPC_diff, alternative = c("greater"))
t.test(dfWide$RSFC_VTAHPC_diff[dfWide$group == "cont"], alternative = c("greater"))
t.test(dfWide$RSFC_VTAHPC_diff[dfWide$group == "exp"], alternative = c("greater"))
t.test(dfWide$RSFC_VTAHPC_diff ~ dfWide$group)

t.test(dfWide$RSFC_VTAHPC_diff_spearman, alternative = c("greater"))
t.test(dfWide$RSFC_VTAHPC_diff_spearman[dfWide$group == "cont"], alternative = c("greater"))
t.test(dfWide$RSFC_VTAHPC_diff_spearman[dfWide$group == "exp"], alternative = c("greater"))
t.test(dfWide$RSFC_VTAHPC_diff_spearman ~ dfWide$group) 

# create data in long format
df_long_pearson <- melt(dfWide, id=c("ID", "group"), measure = c("RSFC_VTAHPC_run1_z", "RSFC_VTAHPC_run2_z"))
names(df_long_pearson) <- c("ID", "group", "run", "RSFC_pearson")
df_long_pearson$run <- ifelse(df_long_pearson$run == "RSFC_VTAHPC_run1_z", "1", "2")

df_long_spearman <- melt(dfWide, id=c("ID", "group"), measure = c("RSFC_VTAHPC_run1_z_spearman", "RSFC_VTAHPC_run2_z_spearman"))
names(df_long_spearman) <- c("ID", "group", "run", "RSFC_spearman")
df_long_spearman$run <- ifelse(df_long_spearman$run == "RSFC_VTAHPC_run1_z_spearman", "1", "2")

# mixed model anova
df_long_pearson$ID <- factor(df_long_pearson$ID)
summary(aov(RSFC_pearson~ group*run + Error(ID/run), data = df_long_pearson)) 

df_long_spearman$ID <- factor(df_long_spearman$ID)
summary(aov(RSFC_spearman~ group*run + Error(ID/run), data = df_long_spearman)) 

# ancova predicting RSFC in post learning rest by group and pre learning rest
summary(lm(RSFC_VTAHPC_run2_z~group + RSFC_VTAHPC_run1_z, dfWide))

summary(lm(RSFC_VTAHPC_run2_z_spearman~group + RSFC_VTAHPC_run1_z_spearman, dfWide))

# create data frame for spaghetti plot
output <- merge(df_long_pearson, df_long_spearman, by = c("ID", "group", "run"))
output$group <- ifelse(output$group == "cont", "No reward", "Reward")

output <- melt(output, id=c("ID", "group", "run"), measure = c("RSFC_pearson", "RSFC_spearman"))
names(output) <- c("ID", "group", "run", "coefficient", "RSFC")

# spaghetti plot
graph <- ggplot(data = output, aes(x = run, y = RSFC, group = ID, col = group)) + 
  geom_line() +
  facet_grid(. ~ coefficient) +
  stat_summary(aes(group = group), geom = "point", fun.y = mean, size = 6) +
  theme_classic() +
  labs(x="resting state run", y="RSFC", title = "HPC-VTA RSFC") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20, face="bold"), title=element_text(size =20, face="bold"))
print(graph)
setwd(brainDir)
print(paste0("Spaghettiplot_HPCVTA_RSFC.jpeg")) # save plot
ggsave(paste0("Spaghettiplot_HPCVTA_RSFC.jpeg")) # save plot

# create data frame for violin plot
output <- melt(dfWide, id=c("ID", "group"), measure = c("RSFC_VTAHPC_diff", "RSFC_VTAHPC_diff_spearman"))
names(output) <- c("ID", "group", "coefficient", "RSFC_change")
output$group <- ifelse(output$group == "cont", "No reward", "Reward")
output$coefficient <- ifelse(output$coefficient == "RSFC_VTAHPC_diff", "Pearson", "Spearman")


# violin plot of difference scores
graph <- ggplot(output, aes(x=group, y=RSFC_change, fill = group)) + 
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.1) +
  geom_jitter(size = 1, shape=1, position=position_jitter(0.2)) +
  theme_classic() +
  labs(x="Experimental Condition", y="Change post > pre", title = "Difference in HPC-VTA RSFC") +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20, face="bold"), title=element_text(size =20, face="bold")) +
  coord_cartesian(ylim = c(-1, 1)) + 
  facet_grid(. ~ coefficient)
print(graph)
print(paste0("Violinplot_HPCVTA_RSFC_changes.jpeg")) # save plot
ggsave(paste0("Violinplot_HPCVTA_RSFC_changes.jpeg")) # save plot




########## 11. Create table showing the correlation between RSFC changes and memory measurements ########## 

DV_wide <- c(paste0(memoryLabels, "_abs")) # absolute sum scores
DV_wide <- c(DV_wide, paste0("curBeta_",memoryLabels)) # betas
DV_wide <- c(DV_wide, paste0("curCor_",memoryLabels)) # correlation
for (mem in 1:length(memoryLevels)) { # benefits
  DV_wide <- c(DV_wide, paste0("curBen_cont_", memoryLabels[mem], collapse = ", "))
  DV_wide <- c(DV_wide, paste0("curBen_dich_", memoryLabels[mem], collapse = ", "))
  DV_wide <- c(DV_wide, paste0("curBen_rel_", memoryLabels[mem], collapse = ", "))
}

# for all dependent variables compute correlation between changes in RSFC and memory for whole sample as well as in each group
df_cor <- data.frame(var   = character(length(DV_wide)),
                     pearson = numeric(length(DV_wide)),
                     p_value_pearson = numeric(length(DV_wide)),
                     pearson_noR = numeric(length(DV_wide)),
                     p_value_pearson_noR = numeric(length(DV_wide)),
                     pearson_R = numeric(length(DV_wide)),
                     p_value_pearson_R = numeric(length(DV_wide)),
                     corrDiff_pearson = numeric(length(DV_wide)),
                     p_value_corrDiff_pearson = numeric(length(DV_wide)),
                     spearman = numeric(length(DV_wide)),
                     p_value_spearman = numeric(length(DV_wide)),
                     spearman_noR = numeric(length(DV_wide)),
                     p_value_spearman_noR = numeric(length(DV_wide)),
                     spearman_R = numeric(length(DV_wide)),
                     p_value_spearman_R = numeric(length(DV_wide)),
                     corrDiff_spearman = numeric(length(DV_wide)),
                     p_value_corrDiff_spearman = numeric(length(DV_wide)),
                     stringsAsFactors=FALSE)
for(DV in 1:length(DV_wide)) {
  # paste variable
  df_cor$var[DV] <- paste(DV_wide[DV])
  # compute correlarion with FC computed with pearson
  cor <- cor.test(dfWide[[paste(DV_wide[DV])]], dfWide$RSFC_VTAHPC_diff)
  df_cor$pearson[DV] <- cor$estimate # correlation
  df_cor$p_value_pearson[DV] <- cor$p.value # p value
  cor_noR <- cor.test(dfWide[[paste(DV_wide[DV])]][dfWide$group=="cont"], dfWide$RSFC_VTAHPC_diff[dfWide$group=="cont"])
  df_cor$pearson_noR[DV] <- cor_noR$estimate # correlation
  df_cor$p_value_pearson_noR[DV] <-cor_noR$p.value # p value
  cor_R <- cor.test(dfWide[[paste(DV_wide[DV])]][dfWide$group=="exp"], dfWide$RSFC_VTAHPC_diff[dfWide$group=="exp"])
  df_cor$pearson_R[DV] <- cor_R$estimate # correlation
  df_cor$p_value_pearson_R[DV] <-cor_R$p.value # p value
  # check for group diff in correlation
  corDiff <- cocor::cocor.indep.groups(cor_noR$estimate, cor_R$estimate, 25, 25, alternative = "two.sided",
                            test = "all", alpha = 0.05, conf.level = 0.95, null.value = 0,
                            data.name = NULL, var.labels = NULL, return.htest = FALSE)
  attributes(corDiff)
  
  df_cor$corrDiff_pearson[DV] <- corDiff@diff
  df_cor$p_value_corrDiff_pearson[DV] <- corDiff@fisher1925$p.value
  
  # compute correlarion with FC computed with spearman (using spearman)
  cor <- cor.test(dfWide[[paste(DV_wide[DV])]], dfWide$RSFC_VTAHPC_diff_spearman, method = "spearman")
  df_cor$spearman[DV] <- cor$estimate # correlation
  df_cor$p_value_spearman[DV] <- cor$p.value # p value
  cor_noR <- cor.test(dfWide[[paste(DV_wide[DV])]][dfWide$group=="cont"], dfWide$RSFC_VTAHPC_diff_spearman[dfWide$group=="cont"], method = "spearman")
  df_cor$spearman_noR[DV] <- cor_noR$estimate # correlation
  df_cor$p_value_spearman_noR[DV] <-cor_noR$p.value # p value
  cor_R <- cor.test(dfWide[[paste(DV_wide[DV])]][dfWide$group=="exp"], dfWide$RSFC_VTAHPC_diff_spearman[dfWide$group=="exp"], method = "spearman")
  df_cor$spearman_R[DV] <- cor_R$estimate # correlation
  df_cor$p_value_spearman_R[DV] <-cor_R$p.value # p value
  # check for group diff in correlation
  corDiff <- cocor::cocor.indep.groups(cor_noR$estimate, cor_R$estimate, 25, 25, alternative = "two.sided",
                                       test = "all", alpha = 0.05, conf.level = 0.95, null.value = 0,
                                       data.name = NULL, var.labels = NULL, return.htest = FALSE)
  attributes(corDiff)
  
  df_cor$corrDiff_spearman[DV] <- corDiff@diff
  df_cor$p_value_corrDiff_spearman[DV] <- corDiff@fisher1925$p.value
}
# 
df_cor[, 2:17] <- round(df_cor[, 2:17], digits = 5)
setwd(brainDir)
write.csv(df_cor, file = "Correlation_RSFC_changes_memory.csv", row.names = F)

########## 12. Create scatter plot to show relationship between HPC-VTA RSFC and memory measurements ########## 

# select variables to be plotted
output <- dfWide[,c("ID", "group", DV_wide, "RSFC_VTAHPC_diff", "RSFC_VTAHPC_diff_spearman")]
output$group <- ifelse(output$group == "cont", "No reward", "Reward")

# set variables to loop through
pp <- 0
ppp <- 0

# loop over all variables to plot
for(p in 1:length(plotVars)) {
  
  plot <- plotVars[p]
  print(plot)
  
  if (plot %in% c("recall", "recollection", # absolute values 
                  "betaRecall", "betaRecollection", # beta values
                  "absoluteBenefitRecall_dich", "absoluteBenefitRecollection_dich", # dichotomised curiosity benefit
                  "absoluteBenefitRecall_cont", "absoluteBenefitRecollection_cont", # continuous curiosity benefit
                  "relativeBenefitRecall", "relativeBenefitRecollection", # relative curiosity benefit
                  "correlationRecall", "correlationRecollection")){ # correlation
    # create data frame in long format
    # first for person
    columns <- c("ID", "group", threshold_recoll1[p], threshold_recoll2[p], "RSFC_VTAHPC_diff")
    output_pearson <- output[,columns]
    output_pearson <- reshape2::melt(output_pearson, id=c("ID","group", "RSFC_VTAHPC_diff"))
    names(output_pearson) <- c("ID", "group", "Pearson",  "criteria", "performance")
    
    # then for spearman
    columns <- c("ID", "group", threshold_recoll1[p], threshold_recoll2[p], "RSFC_VTAHPC_diff_spearman")
    output_spearman <- output[,columns]
    output_spearman <- reshape2::melt(output_spearman, id=c("ID","group", "RSFC_VTAHPC_diff_spearman"))
    names(output_spearman) <- c("ID", "group", "Spearman",  "criteria", "performance")
    
    # then merge them
    output_plot <-  merge(output_pearson, output_spearman, by = c("ID", "group", "criteria", "performance"))
    rm(output_pearson, output_spearman)
    output_plot <- reshape2::melt(output_plot, id=c("ID","group", "criteria", "performance"))
    names(output_plot) <- c("ID", "group", "criteria", "performance", "coefficient", "RSFC_change")
    
  } else if (plot %in% c("recognition", "betaRecognition", "relativeBenefitRecognition", "correlationRecognition",
                         "absoluteBenefitRecognition_dich", "absoluteBenefitRecognition_cont", "relativeBenefitRecognition")){
    # create data frame in long format
    pp <- pp+1
    columns <- c("ID", "group", threshold_recog[pp], "RSFC_VTAHPC_diff", "RSFC_VTAHPC_diff_spearman")
    output_plot <- output[,columns]
    output_plot <- reshape2::melt(output_plot, id=c("ID","group", threshold_recog[pp]))
    names(output_plot) <-  c("ID", "group", "performance",  "coefficient", "RSFC_change")
    output_plot$coefficient <- ifelse(output_plot$coefficient == "RSFC_VTAHPC_diff", "Pearson", "Spearman")
    
  } else if (plot %in% c("remembered", "betaRemembered", "correlationRemembered",
                         "absoluteBenefitRemembered_dich", "absoluteBenefitRemembered_cont", "relativeBenefitRemembered")){
    # create data frame in long format
    # first for person
    ppp <- ppp+1
    columns <- c("ID", "group", threshold_remembered1[ppp], threshold_remembered2[ppp], threshold_remembered3[ppp], threshold_remembered4[ppp], "RSFC_VTAHPC_diff")
    output_pearson <- output[,columns]
    output_pearson <- reshape2::melt(output_pearson, id=c("ID","group", "RSFC_VTAHPC_diff"))
    names(output_pearson) <- c("ID", "group", "Pearson",  "criteria", "performance")
    # then for spearman
    columns <- c("ID", "group", threshold_remembered1[ppp], threshold_remembered2[ppp], threshold_remembered3[ppp], threshold_remembered4[ppp], "RSFC_VTAHPC_diff_spearman")
    output_spearman <- output[,columns]
    output_spearman <- reshape2::melt(output_spearman, id=c("ID","group", "RSFC_VTAHPC_diff_spearman"))
    names(output_spearman) <- c("ID", "group", "Spearman",  "criteria", "performance")
    # then merge them
    output_plot <-  merge(output_pearson, output_spearman, by = c("ID", "group", "criteria", "performance"))
    rm(output_pearson, output_spearman)
    output_plot <- reshape2::melt(output_plot, id=c("ID","group", "criteria", "performance"))
    names(output_plot) <- c("ID", "group", "criteria", "performance", "coefficient", "RSFC_change")
    output_plot$method <- ifelse(output_plot$criteria ==  threshold_remembered1[ppp] | output_plot$criteria ==  threshold_remembered3[ppp], "strict recall", "lenient recall")
    output_plot$criteria <- ifelse(output_plot$criteria == threshold_remembered1[ppp] | output_plot$criteria == threshold_remembered2[ppp], "above average conf", "high conf")
    
  } 
  
  # create scatter plot
  graph <- ggplot(get(paste0("output_plot")), aes(x=RSFC_change, y=performance, col = group)) + 
    geom_point() +
    geom_smooth(method=lm, aes(fill=group), alpha = 0.3, fullrange = T) +
    theme_classic() +
    labs(x="Difference in HPC-VTA RSFC (post > pre)", y="Performance index", title = paste("RSFC changes and", plot), fill = "Group", col  = "Group") +
    theme(axis.text=element_text(size=20), axis.title=element_text(size=20, face="bold"), title=element_text(size =20, face="bold")) #+
  # change facet grid depending on dependent variable
  if (plot %in% c("recall", "recollection", # absolute values 
                  "betaRecall", "betaRecollection", # beta values
                  "absoluteBenefitRecall_dich", "absoluteBenefitRecollection_dich", # dichotomised curiosity benefit
                  "absoluteBenefitRecall_cont", "absoluteBenefitRecollection_cont", # continuous curiosity benefit
                  "relativeBenefitRecall", "relativeBenefitRecollection", # relative curiosity benefit
                  "correlationRecall", "correlationRecollection")){
    graph <- graph + facet_grid(coefficient ~ criteria) 
  } else if (plot %in% c( "remembered", "betaRemembered", "correlationRemembered",
                          "absoluteBenefitRemembered_dich", "absoluteBenefitRemembered_cont", "relativeBenefitRemembered")){
    graph <- graph + facet_grid(coefficient ~ criteria+method)
  }  else if (plot %in% c("recognition", "betaRecognition", "relativeBenefitRecognition", "correlationRecognition",
                          "absoluteBenefitRecognition_dich", "absoluteBenefitRecognition_cont", "relativeBenefitRecognition")){
    graph <- graph + facet_grid(coefficient ~ .) 
  }
  # change y axis depending on variable
  if (plot %in% c("recall", "recognition", "recollection", "remembered")){
    graph <- graph + coord_cartesian(ylim = c(-5, 41))
  } else if (plot %in% c("absoluteBenefitRecall_dich", "absoluteBenefitRecollection_dich", # dichotomised curiosity benefit
                         "absoluteBenefitRecall_cont", "absoluteBenefitRecollection_cont", # continuous curiosity benefit
                         "absoluteBenefitRecognition_dich", "absoluteBenefitRecognition_cont",
                         "absoluteBenefitRemembered_dich", "absoluteBenefitRemembered_cont")){
    graph <- graph + coord_cartesian(ylim = c(-23, 23))
  } else if (plot %in% c("relativeBenefitRecall", "relativeBenefitRecollection",
                         "relativeBenefitRecognition", "relativeBenefitRemembered",
                         "correlationRecall", "correlationRecollection",
                         "correlationRecognition", "correlationRemembered")){
    graph <- graph + coord_cartesian(ylim = c(-1, 1))
  } else if (plot %in% c("betaRecall", "betaRecognition", "betaRemembered", "betaRecollection")){
    graph <- graph + coord_cartesian(ylim = c(-0.4, 0.4))
  }
  
  # print graph and save
  print(graph)
  print(paste0("Scatterplot_", plot, ".jpeg"))
  ggsave(paste0("Scatterplot_", plot, ".jpeg"))
  
}
