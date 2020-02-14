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

# check whether these directories exist, if not create them
ifelse(!dir.exists(ratingsDir), dir.create(ratingsDir), FALSE) 
ifelse(!dir.exists(memoryDir), dir.create(memoryDir), FALSE) 
ifelse(!dir.exists(questDir), dir.create(questDir), FALSE) 
ifelse(!dir.exists(brainDir), dir.create(brainDir), FALSE) 

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
                          "cuedRecallLenient", "cuedRecallStrict", 
                          "recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
                          "rememberedLenient", "rememberedStrict",
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
xlsx::write.xlsx(descriptives, file="Descriptives_dependentVariables.xlsx", sheetName = "Sheet1")
rm(descriptives)


########## 2. Compute lmer model predicting memory performance using curiosity as a continous variable and group effect coded ########## 

# define dependent variables 
workspace <- list.files(path = file.path(codedDir), pattern = "_CP.csv") # check whether the data is coded yet or not
if(length(workspace) == 0) { # if data is not coded yet, only look at recognition performance
  dependentVariables <- c("recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
                          "confidence", "confidenceCorrectTrials")
} else {
  dependentVariables <- c("cuedRecallLenient", "cuedRecallStrict",
                          "recognition","recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
                          "rememberedLenient", "rememberedStrict",
                          "confidence", "confidenceCorrectTrials")
}

##### 2.1 basic LME predicting memory performance with mean-centered curiosity, effect-coded reward and their interaction #####

# loop over dependent variables to compute LME
for (DV in 1:length(dependentVariables)){
  # define model based on scaling of variable
  if (max(dfLong[, dependentVariables[DV]], na.rm = T) > 1){ # if DV continuous
    LMEmodel_curiosityContinuous <- lmerTest::lmer(dfLong[, dependentVariables[DV]] ~ groupEffectCoded*curiosityGroupMeanCentered + (1+curiosityGroupMeanCentered|ID) + (1|stimID), data = dfLong)
  } else { # if DV binary
    LMEmodel_curiosityContinuous <- glmer(dfLong[, dependentVariables[DV]] ~ groupEffectCoded*curiosityGroupMeanCentered + (1+curiosityGroupMeanCentered|ID) + (1|stimID), family = "binomial"(link = 'logit'), data = dfLong)
  }
  # put coefficients into data frame
  summaryCuriosityContinuous <- summary(LMEmodel_curiosityContinuous)
  summaryCuriosityContinuousCoefficients <- as.data.frame(summaryCuriosityContinuous$coefficients)
  # change row names
  row.names(summaryCuriosityContinuousCoefficients) <- c(paste0("LME_",dependentVariables[DV], "_Intercept"), paste0("LME_",dependentVariables[DV], "_groupEffectCoded"), 
                                                         paste0("LME_",dependentVariables[DV], "_curiosityGroupMeanCentered"), 
                                                         paste0("LME_",dependentVariables[DV], "_interaction"))
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
xlsx::write.xlsx(LMEresults, file=paste0("LME_Results_", version, "_CuriosityRatingsAndMemoryScores_byCuriosityAsContinuousVariable.xlsx"), sheetName = "Sheet1")


##### 2.2 LME with covariates predicting memory performance with mean-centered curiosity, effect-coded reward and their interaction #####

# define control vars
controlVar <- c("nback_accurary", "BAS_rewardresponsiveness", "ableToSee")
subjects  <- as.character(dfWide$fMRI)

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
  for (DV in 1:length(dependentVariables)){
    
    # define model based on scaling of variable
    if (max(dfLong[, dependentVariables[DV]], na.rm = T) > 1){
      LMEmodel_curiosityContinuousControlled <- lmerTest::lmer(dfLong[, dependentVariables[DV]] ~ get(controlVar[cov]) + groupEffectCoded*curiosityGroupMeanCentered + (1+curiosityGroupMeanCentered|ID) + (1|stimID), data = dfLong)
    }else{
      LMEmodel_curiosityContinuousControlled <- glmer(dfLong[, dependentVariables[DV]] ~ get(controlVar[cov]) + groupEffectCoded*curiosityGroupMeanCentered + (1+curiosityGroupMeanCentered|ID) + (1|stimID), family = "binomial"(link = 'logit'), data = dfLong)
    }
    # put coefficients into data frame
    summaryCuriosityContinuousControlled <- summary(LMEmodel_curiosityContinuousControlled)
    summaryCuriosityContinuousCoefficientsControlled <- as.data.frame(summaryCuriosityContinuousControlled$coefficients)
    # change row names of data frame
    row.names(summaryCuriosityContinuousCoefficientsControlled) <- c(paste0("LME_", dependentVariables[DV], "_Intercept"), paste0("LME_",dependentVariables[DV], "_", controlVar[cov]),
                                                                     paste0("LME_", dependentVariables[DV], "_groupEffectCoded"),
                                                                     paste0("LME_", dependentVariables[DV], "_curiosityGroupMeanCentered"),
                                                                     paste0("LME_", dependentVariables[DV], "_interaction"))
    # round results and save
    if (DV == 1) {
      LMEresultsControlled <- summaryCuriosityContinuousCoefficientsControlled
    }else {
      temp_LMEresults <- summaryCuriosityContinuousCoefficientsControlled
      LMEresultsControlled <- rbind.all.columns(LMEresultsControlled, temp_LMEresults) #rbind all columns will induce NA if there was initially no data saved in the loop per participant
      rm(temp_LMEresults)
    }
    rm(summaryCuriosityContinuousCoefficientsControlled, summaryCuriosityContinuousControlled)
  }
  # round results and save
  LMEresultsControlled <- round(LMEresultsControlled, digits = 5)
  xlsx::write.xlsx(LMEresultsControlled, file=paste0("LME_Results_", version, "_CuriosityRatingsAndMemoryScores_byCuriosityAsContinuousVariable_controlledFor_", controlVar[cov], ".xlsx"), sheetName = "Sheet1")
  
}

# When there seem to be an interaction effect like nback, you can test a model like memory~nback*group to see if the interaction is indeed significant. If it is significant, there might be a way to show that memory performance is different for s particular subset of participants (e.g. those who are low in the nback task). Given that reward group showed numerically lower memory performance but it is still informative for us.

for (DV in 1:length(dependentVariables)){
  # define model based on scaling of variable
  if (max(dfLong[, dependentVariables[DV]], na.rm = T) > 1){ # if DV continuous
    LMEmodel_curiosityContinuous <- lmerTest::lmer(dfLong[, dependentVariables[DV]] ~ groupEffectCoded*nback_accurary + (1|ID) + (1|stimID), data = dfLong)
  } else { # if DV binary
    LMEmodel_curiosityContinuous <- glmer(dfLong[, dependentVariables[DV]] ~ groupEffectCoded*nback_accurary + (1|ID) + (1|stimID), family = "binomial"(link = 'logit'), data = dfLong)
  }
  # put coefficients into data frame
  summaryCuriosityContinuous <- summary(LMEmodel_curiosityContinuous)
  summaryCuriosityContinuousCoefficients <- as.data.frame(summaryCuriosityContinuous$coefficients)
  # change row names
  row.names(summaryCuriosityContinuousCoefficients) <- c(paste0("LME_",dependentVariables[DV], "_Intercept"), paste0("LME_",dependentVariables[DV], "_groupEffectCoded"), 
                                                         paste0("LME_",dependentVariables[DV], "_nback_accurary"), 
                                                         paste0("LME_",dependentVariables[DV], "_interaction"))
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


##### 2.3 basic LME predicting memory performance using subset of reward group (i.e. those who believed in the manipulation) #####

# create subset of data
dfWide$reward_subset <- ifelse(dfWide$group == "int", 1, #include everyone from control group
                               ifelse(dfWide$rewardBelief_score > 3 & dfWide$group == "ext", 1, 0)) # and those from reward group believing in manipulation

for (s in seq_along(subjects)){
  # take subset variable for subject from dfWide and paste it to dfLong
  dfLong[dfLong$fMRI == subjects[s], "reward_subset"] <- dfWide[dfWide$fMRI == subjects[s], "reward_subset"]
}

dfLong_subset <- subset(dfLong, dfLong$reward_subset == 1)
# loop over dependent variables to compute LME
for (DV in 1:length(dependentVariables)){
  # define model based on scaling of variable
  if (max(dfLong_subset[, dependentVariables[DV]], na.rm = T) > 1){ # if DV continuous
    LMEmodel_curiosityContinuous <- lmerTest::lmer(dfLong_subset[, dependentVariables[DV]] ~ groupEffectCoded*curiosityGroupMeanCentered + (1+curiosityGroupMeanCentered|ID) + (1|stimID), data = dfLong_subset)
  } else { # if DV binary
    LMEmodel_curiosityContinuous <- glmer(dfLong_subset[, dependentVariables[DV]] ~ groupEffectCoded*curiosityGroupMeanCentered + (1+curiosityGroupMeanCentered|ID) + (1|stimID), family = "binomial"(link = 'logit'), data = dfLong_subset)
  }
  # put coefficients into data frame
  summaryCuriosityContinuous <- summary(LMEmodel_curiosityContinuous)
  summaryCuriosityContinuousCoefficients <- as.data.frame(summaryCuriosityContinuous$coefficients)
  # change row names
  row.names(summaryCuriosityContinuousCoefficients) <- c(paste0("LME_subset_",dependentVariables[DV], "_Intercept"), paste0("LME_subset_",dependentVariables[DV], "_groupEffectCoded"), 
                                                         paste0("LME_subset_",dependentVariables[DV], "_curiosityGroupMeanCentered"), 
                                                         paste0("LME_subset_",dependentVariables[DV], "_interaction"))
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
xlsx::write.xlsx(LMEresults, file=paste0("LME_Results_", version, "_subset_CuriosityRatingsAndMemoryScores_byCuriosityAsContinuousVariable.xlsx"), sheetName = "Sheet1")



########## 3. Create barplots to visualise the effects of curiosity and reward on memory performance ########## 
if (exists("dfLong") == F) {
  setwd(preprocessedDir)
  dfLong <- xlsx::read.xlsx("long_MAGMOT.xlsx", sheetName = "Sheet1")#,  showWarnings = FALSE)
} 

# create a dichomotised curiosity variable using mean-cenetred curiosity
dfLong$curiosity_dichotom <- ifelse(dfLong$curiosity_dichotom == -1, "below", 
                                    ifelse(dfLong$curiosity_dichotom == 1, "above", NA)) # create curiosity_dichotom as factor
# define grouping variables
groupingVariables <- c("mediansplitCuriosity_MAGMOT", "mediansplitCuriosityWithinSubject", "curiosity_dichotom")
groupingVariables <- c("mediansplitCuriosityWithinSubject", "curiosity_dichotom")

# determine dependent variables to loop over (same as above)
dependentVariables <- c("cuedRecallLenient", "cuedRecallStrict",
                        "recognition","recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
                        "rememberedLenient", "rememberedStrict",
                        "confidence", "confidenceCorrectTrials",
                        "confidenceGroupMeanCentered", "confidenceGroupMeanCenteredCorrectTrials")

# determine colours to be used in histogram
cols <- c("above" = "#F8766D", "below" = "#00BFC4", "all tricks" = "grey", "median" = "#C77CFF")

# set directory where to save plots
setwd(memoryDir)

# loop over dependent variables to create bar plots
for (DV in dependentVariables){
  
  # loop over grouping variables for continuous curiosity variable
  for (g in 1:length(groupingVariables)){
    # create a data frame containing the data for each group divided regarding curiosity ratings
    outputGroup <- summarySEwithin(dfLong, measurevar=DV, betweenvars="group", withinvars=groupingVariables[g], idvar="ID", na.rm = T)
    levels(outputGroup$group) <- c("Reward","No reward")
    names(outputGroup)[names(outputGroup) == groupingVariables[g]] <- "cutoff"
    
    # remove NAs (necessary for curiosity_dichotomous)
    outputGroup <- na.omit(outputGroup) 
    
    # add grouping variable as column
    outputGroup$groupingvar <- groupingVariables[g]
    
    # create a data frame containing the data for each group regardless of curiosity ratings
    outputAll <- summarySE(dfLong, measurevar=DV, groupvars="group", na.rm=T,
                           conf.interval=.95, .drop=TRUE)
    levels(outputAll$group) <- c("Reward","No reward")
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
  if (DV %in% c("cuedRecallStrict", "cuedRecallLenient", "rememberedStrict", "rememberedLenient")){
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



########## 4. Create histograms to further investigate the relation between reward, curiosity and memory ########## 

# define variables
dependentVariables <- c("cuedRecallLenient", "cuedRecallStrict",
                        "recognition","recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
                        "rememberedLenient", "rememberedStrict")
# create data frame and recode memory as factors
output <- dfLong[,c("ID", "group", "curiosityGroupMeanCentered", dependentVariables)]
output$group <- ifelse(output$group == "int", "No reward", "Reward")
output$cuedRecallLenient <- ifelse(output$cuedRecallLenient == 0, "forgotten", "remembered")
output$cuedRecallStrict <- ifelse(output$cuedRecallStrict == 1, "remembered", "forgotten")
output$recognition <- ifelse(output$recognition == 1, "remembered", "forgotten")
output$recognitionConfLevel_4_5_6 <- ifelse(output$recognitionConfLevel_4_5_6 == 1, "remembered", "forgotten")
output$recognitionAboveMeanConf <- ifelse(output$recognitionAboveMeanConf == 1, "remembered", "forgotten")
output$rememberedLenient <- ifelse(output$rememberedLenient == 1, "remembered", "forgotten")
output$rememberedStrict <- ifelse(output$rememberedStrict == 1, "remembered", "forgotten")

# loop over dependent variables to create histograms 
for (DV in dependentVariables){
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
output <- as.data.frame(rbind(output$int, output$ext))
output$mot <- rep(c("int","ext"), each = 1)

outg <- ggplot(output, aes(mot, mean, fill = mot))
outg + geom_bar(stat="identity", position="dodge") + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9))  + scale_x_discrete(limits=c("ext","int")) + labs(x="Experimental condition", y="Age", fill = "Experimental Condition", title = paste("demogs I", version)) + theme_classic() + scale_fill_discrete(guide=FALSE)

t.test(dfWide$age ~ dfWide$group)

# gender
rm(output)
output <- plyr::count(dfWide, vars = c("gender","group"))

# output <- count(df, c('gender','cond'))
outg <- ggplot(output, aes(group, freq, fill = gender))
outg + geom_bar(stat="identity", position="fill")+ scale_x_discrete(limits=c("int","ext")) + labs(x="Experimental condition", y="Frequency", fill = "Gender", title = paste("demogs II", version)) + theme_classic() 


########## 2. Look at the effectiveness of reward manipulation ########## 

# check whether responses ppt gave are plausible
dfWide$fMRI[!is.na(dfWide$rewardBelief_score)] # missing value for MAGMOT_06, value for MAGMOT_35
dfWide$rewardBelief_score[dfWide$fMRI == "MAGMOT_35"] <- NA
# manipulation check
t.test(dfWide$rewardBelief_score, alternative = "greater", mu = 3)
cor.test(dfWide$rewardBelief_score, dfWide$recognitionAboveMeanConf, use = "pairwise.complete.obs")

dfWide$rewardEffort
cor_rewardEffort <- as.data.frame(t(cor(dfWide$rewardEffort, dfWide[,c("BIS", "BAS_rewardresponsiveness", "BAS_drive", "BAS_funseeking", 
                                                                       "NeedForCognition", "FearOfFailure", "ApproachTemperament", "AvoidanceTemperament", 
                                                                       "TraitCuriosity", "StateCuriosity", 
                                                                       "intrinsicMotivation", "taskEngagement", "interest", "boredom", "effort", "pressure", 
                                                                       "compliance", "ableToSee", 
                                                                       "rewardBelief_score",
                                                                       "corsiSpan", "nback_hitrate", "nback_falsealarmrate", "nback_accurary",
                                                                       "cuedRecallStrict", "cuedRecallLenient",
                                                                       "recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
                                                                       "rememberedStrict", "rememberedLenient", 
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
                                                                             "cuedRecallStrict", "cuedRecallLenient",
                                                                             "recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
                                                                             "rememberedStrict", "rememberedLenient", 
                                                                             "meanConfidence", "meanConfidenceCorrectTrials")], 
                                        use = "pairwise.complete.obs")))
names(cor_rewardBelief) <- "cor_rewardBelief"


cor(dfWide$rewardBelief_score, dfWide[,scales], use = "pairwise.complete.obs")

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
  rating <- as.data.frame(rbind(output$int, output$ext))
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
scales <- names(dfWide[,c("BIS", "BAS_rewardresponsiveness", "BAS_drive", "BAS_funseeking", 
                          "NeedForCognition", "FearOfFailure", "ApproachTemperament", "AvoidanceTemperament", 
                          "TraitCuriosity", "StateCuriosity", 
                          "intrinsicMotivation", "taskEngagement", "interest", "boredom", "effort", "pressure", 
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
  if (means$ext != means$int) {
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
write.csv(effectsizesScales, paste0("effectsizesScalesNeuro_", version, ".csv"))


########## 5. Check whether any of the differences in the questionnaires relates to memory performance and whetehr there is a group effect  ########## 

# able to see
cor_ableToSee <- as.data.frame(t(cor(dfWide$ableToSee, dfWide[,c("BIS", "BAS_rewardresponsiveness", "BAS_drive", "BAS_funseeking", 
                                                                 "NeedForCognition", "FearOfFailure", "ApproachTemperament", "AvoidanceTemperament", 
                                                                 "TraitCuriosity", "StateCuriosity", 
                                                                 "intrinsicMotivation", "taskEngagement", "interest", "boredom", "effort", "pressure", 
                                                                 "compliance", 
                                                                 "rewardBelief_score",
                                                                 "corsiSpan", "nback_hitrate", "nback_falsealarmrate", "nback_accurary",
                                                                 "cuedRecallStrict", "cuedRecallLenient",
                                                                 "recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
                                                                 "rememberedStrict", "rememberedLenient", 
                                                                 "meanConfidence", "meanConfidenceCorrectTrials")], 
                                     use = "pairwise.complete.obs")))
names(cor_ableToSee) <- "cor_ableToSee"
cor.test(dfWide$recognitionAboveMeanConf, dfWide$ableToSee)
cor.test(dfWide$recognitionAboveMeanConf[dfWide$group =="ext"], dfWide$ableToSee[dfWide$group =="ext"])
cor.test(dfWide$recognitionAboveMeanConf[dfWide$group =="int"], dfWide$ableToSee[dfWide$group =="int"])
cor.test(dfWide$recognitionAboveMeanConf, dfWide$compliance)

# reward responsiveness
cor.test(dfWide$recognitionAboveMeanConf, dfWide$BAS_rewardresponsiveness)
cor.test(dfWide$recognitionAboveMeanConf[dfWide$group =="ext"], dfWide$BAS_rewardresponsiveness[dfWide$group =="ext"])
cor.test(dfWide$recognitionAboveMeanConf[dfWide$group =="int"], dfWide$BAS_rewardresponsiveness[dfWide$group =="int"])

# nback
cor.test(dfWide$recognitionAboveMeanConf, dfWide$nback_accurary)
cor.test(dfWide$recognitionAboveMeanConf[dfWide$group =="ext"], dfWide$nback_accurary[dfWide$group =="ext"])
cor.test(dfWide$recognitionAboveMeanConf[dfWide$group =="int"], dfWide$nback_accurary[dfWide$group =="int"])

# create intercorrelation table

cor_all <- as.data.frame(t(cor(dfWide[,c("BIS", "BAS_rewardresponsiveness", "BAS_drive", "BAS_funseeking", 
                                         "NeedForCognition", "FearOfFailure", "ApproachTemperament", "AvoidanceTemperament", 
                                         "TraitCuriosity", "StateCuriosity", 
                                         "intrinsicMotivation", "taskEngagement", "interest", "boredom", "effort", "pressure", 
                                         "compliance", "ableToSee", "rewardBelief_score",
                                         "corsiSpan", "nback_hitrate", "nback_falsealarmrate", "nback_accurary",
                                         "cuedRecallStrict", "cuedRecallLenient",
                                         "recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
                                         "rememberedStrict", "rememberedLenient", 
                                         "meanConfidence", "meanConfidenceCorrectTrials")], 
                               use = "pairwise.complete.obs")))
cor_int <- as.data.frame(t(cor(dfWide[dfWide$group == "int", c("BIS", "BAS_rewardresponsiveness", "BAS_drive", "BAS_funseeking", 
                                                               "NeedForCognition", "FearOfFailure", "ApproachTemperament", "AvoidanceTemperament", 
                                                               "TraitCuriosity", "StateCuriosity", 
                                                               "intrinsicMotivation", "taskEngagement", "interest", "boredom", "effort", "pressure", 
                                                               "compliance", "ableToSee", "rewardBelief_score",
                                                               "corsiSpan", "nback_hitrate", "nback_falsealarmrate", "nback_accurary",
                                                               "cuedRecallStrict", "cuedRecallLenient",
                                                               "recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
                                                               "rememberedStrict", "rememberedLenient", 
                                                               "meanConfidence", "meanConfidenceCorrectTrials")], 
                               use = "pairwise.complete.obs")))
cor_ext <- as.data.frame(t(cor(dfWide[dfWide$group == "ext", c("BIS", "BAS_rewardresponsiveness", "BAS_drive", "BAS_funseeking", 
                                                               "NeedForCognition", "FearOfFailure", "ApproachTemperament", "AvoidanceTemperament", 
                                                               "TraitCuriosity", "StateCuriosity", 
                                                               "intrinsicMotivation", "taskEngagement", "interest", "boredom", "effort", "pressure", 
                                                               "compliance", "ableToSee", "rewardBelief_score",
                                                               "corsiSpan", "nback_hitrate", "nback_falsealarmrate", "nback_accurary",
                                                               "cuedRecallStrict", "cuedRecallLenient",
                                                               "recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
                                                               "rememberedStrict", "rememberedLenient", 
                                                               "meanConfidence", "meanConfidenceCorrectTrials")], 
                               use = "pairwise.complete.obs")))
cor_all <- round(cor_all, digits = 3)
cor_ext <- round(cor_ext, digits = 3)

names(cor_ableToSee) <- "cor_ableToSee"

########## 6. Compute t-tests and effect sizes for between-group differences in memory scores ########## 

# define dependent variables (i.e. sum scores)
dependentVariablesWide <- c("cuedRecallStrict", "cuedRecallLenient",
                            "recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
                            "rememberedStrict", "rememberedLenient", 
                            "meanConfidence", "meanConfidenceCorrectTrials",
                            #betas
                            "curiosityBeta_cuedRecallStrict", "curiosityBeta_cuedRecallLenient",
                            "curiosityBeta_allConf", "curiosityBeta_highConf", "curiosityBeta_aboveAvgConf",
                            "curiosityBeta_rememberedStrict", "curiosityBeta_rememberedLenient",
                            # benefit
                            "curiosityBenefit_cuedRecallStrict", "curiosityBenefit_cuedRecallStrict_dichotom", "curiosityBenefit_cuedRecallStrict_perc",
                            "curiosityBenefit_cuedRecallLenient", "curiosityBenefit_cuedRecallLenient_dichotom", "curiosityBenefit_cuedRecallLenient_perc",
                            "curiosityBenefit_allConf", "curiosityBenefit_allConf_dichotom", "curiosityBenefit_allConf_perc",
                            "curiosityBenefit_highConf", "curiosityBenefit_highConf_dichotom", "curiosityBenefit_highConf_perc",
                            "curiosityBenefit_aboveAvgConf", "curiosityBenefit_aboveAvgConf_dichotom", "curiosityBenefit_aboveAvgConf_perc",
                            "curiosityBenefit_rememberedStrict", "curiosityBenefit_rememberedStrict_dichotom", "curiosityBenefit_rememberedStrict_perc",
                            "curiosityBenefit_rememberedLenient", "curiosityBenefit_rememberedLenient_dichotom", "curiosityBenefit_rememberedLenient_perc",
                            # correlations
                            "curiosityCorrelation_cuedRecallStrict", "curiosityCorrelation_cuedRecallLenient", 
                            "curiosityCorrelation_allConf", "curiosityCorrelation_highConf", "curiosityCorrelation_aboveAvgConf", 
                            "curiosityCorrelation_rememberedStrict", "curiosityCorrelation_rememberedLenient",                 
                            # brain data
                            "RSFC_VTAHPC_diff", "RSFC_VTAHPC_diff_spearman")


# for all dependent variables compute two-sample t-test and calculate effect size
for(DV in 1:length(dependentVariablesWide)) {
  
  # compute t-test for group difference
  ttest <- t.test(dfWide[,dependentVariablesWide[DV]]~dfWide$group)
  t.stats <- as.data.frame(t(round(c(ttest$statistic, ttest$p.value), digits = 3)))
  names(t.stats) <- c("tValue", "pValue(t)")
  
  # compute wilcox
  wilcox <- wilcox.test(dfWide[,dependentVariablesWide[DV]]~dfWide$group) 
  attributes(wilcox)
  w.stats <- as.data.frame(t(round(c(wilcox$statistic, wilcox$p.value), digits = 3)))
  names(w.stats) <- c("W", "pValue(W)")
  
  # merge t-test and wilconxon's test
  t.stats <- merge(t.stats, w.stats)
  
  # compute mean for each group
  means <- tapply(dfWide[,dependentVariablesWide[DV]], dfWide$group, mean, na.rm = T)
  means <- as.data.frame(t(means))
  means <- merge(t.stats, means)
  
  # compute effect size
  if (means$ext != means$int) {
    data <- dfWide[,c("group", dependentVariablesWide[DV])]
    d <- psych::cohen.d(data, "group")
    cohen <- as.data.frame(d$cohen.d)
    means <- merge(means, cohen)
  }
  
  # put all in a data frame
  row.names(means) <- paste(dependentVariablesWide[DV])
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
effectsizes <- rbind.all.columns(effectsizesScales, effectsizesMemory)
setwd(ratingsDir)
write.csv(effectsizes, paste0("effectsizesAll_", version, ".csv"))

########## 7. Create violin plots for sum scores of memory measures in each group ########## 

# define dependent variables (i.e. sum scores)
dependentVariablesWide <- c("cuedRecallStrict", "cuedRecallLenient",
                            "recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
                            "rememberedStrict", "rememberedLenient", 
                            "meanConfidence", "meanConfidenceCorrectTrials",
                            #betas
                            "curiosityBeta_cuedRecallStrict", "curiosityBeta_cuedRecallLenient",
                            "curiosityBeta_allConf", "curiosityBeta_highConf", "curiosityBeta_aboveAvgConf",
                            "curiosityBeta_rememberedStrict", "curiosityBeta_rememberedLenient",
                            # benefit
                            "curiosityBenefit_cuedRecallStrict", "curiosityBenefit_cuedRecallStrict_dichotom", "curiosityBenefit_cuedRecallStrict_perc",
                            "curiosityBenefit_cuedRecallLenient", "curiosityBenefit_cuedRecallLenient_dichotom", "curiosityBenefit_cuedRecallLenient_perc",
                            "curiosityBenefit_allConf", "curiosityBenefit_allConf_dichotom", "curiosityBenefit_allConf_perc",
                            "curiosityBenefit_highConf", "curiosityBenefit_highConf_dichotom", "curiosityBenefit_highConf_perc",
                            "curiosityBenefit_aboveAvgConf", "curiosityBenefit_aboveAvgConf_dichotom", "curiosityBenefit_aboveAvgConf_perc",
                            "curiosityBenefit_rememberedStrict", "curiosityBenefit_rememberedStrict_dichotom", "curiosityBenefit_rememberedStrict_perc",
                            "curiosityBenefit_rememberedLenient", "curiosityBenefit_rememberedLenient_dichotom", "curiosityBenefit_rememberedLenient_perc",
                            # correlations
                            "curiosityCorrelation_cuedRecallStrict", "curiosityCorrelation_cuedRecallLenient", 
                            "curiosityCorrelation_allConf", "curiosityCorrelation_highConf", "curiosityCorrelation_aboveAvgConf", 
                            "curiosityCorrelation_rememberedStrict", "curiosityCorrelation_rememberedLenient",                 
                            # brain data
                            "RSFC_VTAHPC_diff", "RSFC_VTAHPC_diff_spearman")


output <- dfWide[,c("ID", "group", dependentVariablesWide)]
output$group <- ifelse(output$group == "int", "No reward", "Reward")
setwd(memoryDir)



# create long data frame for: 
#cuedRecall: lenient vs strict
output_recall <- output[,c("ID", "group", "cuedRecallStrict", "cuedRecallLenient")]
output_recall <- reshape2::melt(output_recall, id=c("ID","group"))
names(output_recall) <- c("ID", "group", "criteria", "performance")
#recognition
output_recognition <- output[,c("ID", "group", "recognition")]
output_recognition <- reshape2::melt(output_recognition, id=c("ID","group"))
names(output_recognition) <- c("ID", "group", "criteria", "performance")
#recollection: high vs aboveAvg
output_recollection <- output[,c("ID", "group", "recognitionAboveMeanConf", "recognitionConfLevel_4_5_6")]
output_recollection <- reshape2::melt(output_recollection, id=c("ID","group"))
names(output_recollection) <- c("ID", "group", "criteria", "performance")
#remembered: lenient vs strict
output_remembered <- output[,c("ID", "group", "rememberedStrict", "rememberedLenient")]
output_remembered <- reshape2::melt(output_remembered, id=c("ID","group"))
names(output_remembered) <- c("ID", "group", "criteria", "performance")

#beta recall: lenient vs strict
output_betaRecall <- output[,c("ID", "group", "curiosityBeta_cuedRecallStrict", "curiosityBeta_cuedRecallLenient")]
output_betaRecall <- reshape2::melt(output_betaRecall, id=c("ID","group"))
names(output_betaRecall) <- c("ID", "group", "criteria", "performance")
#beta familarity
output_betaRecognition <- output[,c("ID", "group", "curiosityBeta_allConf")]
output_betaRecognition <- reshape2::melt(output_betaRecognition, id=c("ID","group"))
names(output_betaRecognition) <- c("ID", "group", "criteria", "performance")
#beta recollection: high vs aboveAvg
output_betaRecollection <- output[,c("ID", "group", "curiosityBeta_aboveAvgConf", "curiosityBeta_highConf")]
output_betaRecollection <- reshape2::melt(output_betaRecollection, id=c("ID","group"))
names(output_betaRecollection) <- c("ID", "group", "criteria", "performance")
#beta recall: lenient vs strict
output_betaRemembered <- output[,c("ID", "group", "curiosityBeta_rememberedStrict", "curiosityBeta_rememberedLenient")]
output_betaRemembered <- reshape2::melt(output_betaRemembered, id=c("ID","group"))
names(output_betaRemembered) <- c("ID", "group", "criteria", "performance")

#benefit recall: 2 (lenient vs strict) * 2 (cont vs dichotom)
output_benefitRecall <- output[,c("ID", "group","curiosityBenefit_cuedRecallStrict", "curiosityBenefit_cuedRecallStrict_dichotom",
                                  "curiosityBenefit_cuedRecallLenient", "curiosityBenefit_cuedRecallLenient_dichotom")]
output_benefitRecall <- reshape2::melt(output_benefitRecall, id=c("ID","group"))
names(output_benefitRecall) <- c("ID", "group", "criteria", "performance")
output_benefitRecall$method <- ifelse(output_benefitRecall$criteria == "curiosityBenefit_cuedRecallStrict" | output_benefitRecall$criteria == "curiosityBenefit_cuedRecallLenient", "continuous", "dichotom")
output_benefitRecall$criteria <- ifelse(output_benefitRecall$criteria == "curiosityBenefit_cuedRecallStrict" | output_benefitRecall$criteria == "curiosityBenefit_cuedRecallStrict_dichotom", "cuedRecallStrict", "cuedRecallLenient")
#benefit recognition: cont vs dichotom
output_benefitRecognition <- output[,c("ID", "group","curiosityBenefit_allConf", "curiosityBenefit_allConf_dichotom")]
output_benefitRecognition <- reshape2::melt(output_benefitRecognition, id=c("ID","group"))
names(output_benefitRecognition) <- c("ID", "group", "criteria", "performance")
output_benefitRecognition$method <- ifelse(output_benefitRecognition$criteria == "curiosityBenefit_allConf", "continuous", "dichotom")
output_benefitRecognition$criteria <- "curiosityBenefit_allConf"
#benefit recognition: 2 (high vs aboveAvg) * 2 (cont vs dichotom)
output_benefitRecollection <- output[,c("ID", "group","curiosityBenefit_aboveAvgConf", "curiosityBenefit_aboveAvgConf_dichotom",
                                        "curiosityBenefit_highConf", "curiosityBenefit_highConf_dichotom")]
output_benefitRecollection <- reshape2::melt(output_benefitRecollection, id=c("ID","group"))
names(output_benefitRecollection) <- c("ID", "group", "criteria", "performance")
output_benefitRecollection$method <- ifelse(output_benefitRecollection$criteria == "curiosityBenefit_aboveAvgConf" | output_benefitRecollection$criteria == "curiosityBenefit_highConf", "continuous", "dichotom")
output_benefitRecollection$criteria <- ifelse(output_benefitRecollection$criteria == "curiosityBenefit_aboveAvgConf" | output_benefitRecollection$criteria == "curiosityBenefit_aboveAvgConf_dichotom", "aboveAvgConf", "highConf")
#benefit Remembered: 2 (lenient vs strict) * 2 (cont vs dichotom)
output_benefitRemembered <- output[,c("ID", "group","curiosityBenefit_rememberedStrict", "curiosityBenefit_rememberedStrict_dichotom",
                                      "curiosityBenefit_rememberedLenient", "curiosityBenefit_rememberedLenient_dichotom")]
output_benefitRemembered <- reshape2::melt(output_benefitRemembered, id=c("ID","group"))
names(output_benefitRemembered) <- c("ID", "group", "criteria", "performance")
output_benefitRemembered$method <- ifelse(output_benefitRemembered$criteria == "curiosityBenefit_rememberedStrict" | output_benefitRemembered$criteria == "curiosityBenefit_rememberedLenient", "continuous", "dichotom")
output_benefitRemembered$criteria <- ifelse(output_benefitRemembered$criteria == "curiosityBenefit_rememberedStrict" | output_benefitRemembered$criteria == "curiosityBenefit_rememberedStrict_dichotom", "rememberedStrict", "rememberedLenient")

#percBenefit recall: lenient vs strict
output_percBenefitRecall <- output[,c("ID", "group", "curiosityBenefit_cuedRecallStrict_perc", "curiosityBenefit_cuedRecallLenient_perc")]
output_percBenefitRecall <- reshape2::melt(output_percBenefitRecall, id=c("ID","group"))
names(output_percBenefitRecall) <- c("ID", "group", "criteria", "performance")
#percBenefit familarity
output_percBenefitRecognition <- output[,c("ID", "group", "curiosityBenefit_allConf_perc")]
output_percBenefitRecognition <- reshape2::melt(output_percBenefitRecognition, id=c("ID","group"))
names(output_percBenefitRecognition) <- c("ID", "group", "criteria", "performance")
#percBenefit recollection: high vs aboveAvg
output_percBenefitRecollection <- output[,c("ID", "group", "curiosityBenefit_aboveAvgConf_perc", "curiosityBenefit_highConf_perc")]
output_percBenefitRecollection <- reshape2::melt(output_percBenefitRecollection, id=c("ID","group"))
names(output_percBenefitRecollection) <- c("ID", "group", "criteria", "performance")
#percBenefit recall: lenient vs strict
output_percBenefitRemembered <- output[,c("ID", "group", "curiosityBenefit_rememberedStrict_perc", "curiosityBenefit_rememberedLenient_perc")]
output_percBenefitRemembered <- reshape2::melt(output_percBenefitRemembered, id=c("ID","group"))
names(output_percBenefitRemembered) <- c("ID", "group", "criteria", "performance")



# define which variables to plot
plotVars <- c("recall", "recognition","recollection", "remembered", 
              "betaRecall", "betaRecognition", "betaRecollection", "betaRemembered", 
              "benefitRecall", "benefitRecognition", "benefitRecollection", "benefitRemembered",
              "percBenefitRecall", "percBenefitRecognition", "percBenefitRecollection", "percBenefitRemembered")

for (plot in plotVars){
  
  # Basic violin plot
  graph <- ggplot(get(paste0("output_",plot)), aes(x=group, y=performance, fill = group)) + 
    geom_violin(trim=FALSE) +
    geom_boxplot(width=0.1) +
    geom_jitter(size = 1, shape=1, position=position_jitter(0.2)) +
    theme_classic() +
    labs(x="Experimental Condition", y="Sum score", title = paste(version, plot)) +
    theme(legend.position="none") +
    theme(axis.text=element_text(size=20), axis.title=element_text(size=20, face="bold"), title=element_text(size =20, face="bold")) 
  # change facet grid depending on dependent variable
  if (plot %in% c("recall", "recollection", "betaRecall", "betaRecollection", "remembered", "betaRemembered")){
    graph <- graph + facet_grid(. ~ criteria) 
    if (plot %in% c("recall", "recollection")){
      graph + coord_cartesian(ylim = c(-5, 41))
    } else {
      graph + coord_cartesian(ylim = c(-.25, 0.25))
    }
  }
  # change x y-axis depending on dependent variable
  if (plot %in% c("benefitRecall", "benefitRecognition", "benefitRemembered", "benefitRecollection")){
    graph <- graph + facet_grid(criteria ~ method) +
      coord_cartesian(ylim = c(-23, 23))
  }
  
  print(graph)
  print(paste0("Violonplot_", plot, ".jpeg"))
  ggsave(paste0("Violonplot_", plot, ".jpeg"))
}  


################### 8. PROBLEM, for future stef to fix ###################

# check whether there is a difference in the effect size depending on task block
# NOTE: SOMETHING IS NOT WORKING HERE!!!!! HOWEVER, the effect sizes do not seem to change over blocks, hence tiredness can be excluded
blocks <- c(1,2,3,4)
blockstring <- c("_firstBlock", "_secondBlock", "_thirdBlock", "")
dependentVariablesWide <- c("cuedRecallStrict_perc", "cuedRecallLenient_perc",
                            "recognition_perc", "recognitionConfLevel_4_5_6_perc", "recognitionAboveMeanConf_perc",
                            "rememberedStrict_perc", "rememberedLenient_perc", 
                            "meanConfidence", "meanConfidenceCorrectTrials",
                            #betas
                            "curiosityBeta_cuedRecallStrict", "curiosityBeta_cuedRecallLenient",
                            "curiosityBeta_allConf", "curiosityBeta_highConf", "curiosityBeta_aboveAvgConf",
                            "curiosityBeta_rememberedStrict", "curiosityBeta_rememberedLenient",
                            # benefit
                            "curiosityBenefit_cuedRecallStrict_perc",
                            "curiosityBenefit_cuedRecallLenient_perc",
                            "curiosityBenefit_allConf_perc",
                            "curiosityBenefit_highConf_perc",
                            "curiosityBenefit_aboveAvgConf_perc",
                            "curiosityBenefit_rememberedStrict_perc",
                            "curiosityBenefit_rememberedLenient_perc",
                            # correlations
                            "curiosityCorrelation_cuedRecallStrict", "curiosityCorrelation_cuedRecallLenient", 
                            "curiosityCorrelation_allConf", "curiosityCorrelation_highConf", "curiosityCorrelation_aboveAvgConf", 
                            "curiosityCorrelation_rememberedStrict", "curiosityCorrelation_rememberedLenient")

dependentVariablesWide <- c("cuedRecallStrict_perc", "cuedRecallLenient_perc",
                            "recognition_perc", "recognitionConfLevel_4_5_6_perc", "recognitionAboveMeanConf_perc",
                            "rememberedStrict_perc", "rememberedLenient_perc") # note: I maybe want to look at them in proportions, i.e. out of 12 or 36


# create a data frame that includes the mean + CI for each group and the effect size of the group difference + CI
# aim: data set with col: dependentVar; measure (effect size/group mean);  group (int / ext / effect); value; upper; lower
# rows: dependent variables


for (block in blocks){
  for(DV in 1:length(dependentVariablesWide)) {
    # compute mean for each block and group
    means <- as.data.frame(tapply(dfWide[[paste0(dependentVariablesWide[DV], blockstring[block])]], dfWide$group, mean, na.rm = T)) # maybe not transpose
    names(means) <- "effect" # need to be the same as 
    means$grouping <- row.names(means)
    
    # compute SD for each block and grouping
    sds <- as.data.frame(tapply(dfWide[[paste0(dependentVariablesWide[DV], blockstring[block])]], dfWide$group, sd, na.rm = T)) # maybe not transpose
    names(sds) <- "sd" # need to be the same as 
    sds$grouping <- row.names(sds)
    
    # compute SE
    sds$se[sds$grouping == "int"] <- (sds$sd[sds$grouping == "int"]) / (nrow(dfWide[dfWide$group=="int",]))
    sds$se[sds$grouping == "ext"] <- (sds$sd[sds$grouping == "ext"]) / (nrow(dfWide[dfWide$group=="ext",]))
    
    # merge means and sds
    df <- merge(means, sds, by = "grouping")
    
    # compute confidence interval
    df$upper <- df$effect + 1.96 * df$se
    df$lower <- df$effect - 1.96 * df$se
    
    # calculate effect size
    if (df$effect[df$grouping=="int"] != df$effect[df$grouping=="ext"]) {
      data <- dfWide[,c("group", paste0(dependentVariablesWide[DV], blockstring[block]))]
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
    df$dependentVar <- paste(dependentVariablesWide[DV])
    
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
effectsizesMemoryBlock_long$grouping <- ifelse(effectsizesMemoryBlock_long$grouping == "int", "Mean (no reward)",
                                               ifelse(effectsizesMemoryBlock_long$grouping == "ext", "Mean (reward)",
                                                      ifelse(effectsizesMemoryBlock_long$grouping == "effect", "Cohen's d", NA)))
# define colours
cols <- c("Mean (no reward)" = "#F8766D", "Mean (reward)" = "#00BFC4", "Cohen's d" = "black")

# for each variable in list, create a plot showing the average performance in each group per block as well as the effect size
for(DV in 1:length(dependentVariablesWide)) {
  
  # subset the data frame containing the effect sizes per block and select relevant DV only
  output <- subset(effectsizesMemoryBlock_long, effectsizesMemoryBlock_long$dependentVar == dependentVariablesWide[DV], col = grouping)
  
  # create graph
  graph <- ggplot(data=output, aes(x=blockNumber, y=effect)) +
    geom_point(aes(col = grouping)) + 
    geom_errorbar(aes(ymin=lower, ymax=upper, col = grouping), width=.025) +
    facet_grid(measure ~., scales = "free_y") +
    theme_classic() + scale_color_manual(values = cols) +
    scale_x_discrete(limits=c("1", "2", "3", "all")) +
    labs(x="Task block", y="Value", group = "", title = paste(version, "Proportial effect over time",dependentVariablesWide[DV])) +
    theme(legend.title = element_blank()) +
    theme(axis.text=element_text(size=20), axis.title=element_text(size=20, face="bold"), title=element_text(size =20, face="bold")) 
  
  print(graph)
  print(paste0("Graph_changeOfEffectsizeOverBlocks_", dependentVariablesWide[DV], ".jpeg"))
  # save file
  ggsave(file = paste0("Graph_changeOfEffectsizeOverBlocks_", dependentVariablesWide[DV], ".jpeg"), graph)
  
}


########## 9. Test changes if RSFC between HPC and VTA for significance ########## 
# t-test for groupwise differences in RSFC change
t.test(dfWide$RSFC_VTAHPC_diff)
t.test(dfWide$RSFC_VTAHPC_diff ~ dfWide$group)
t.test(dfWide$RSFC_VTAHPC_diff_spearman)
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
output$group <- ifelse(output$group == "int", "No reward", "Reward")

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
output$group <- ifelse(output$group == "int", "No reward", "Reward")
output$coefficient <- ifelse(output$coefficient == "RSFC_VTAHPC_diff", "Pearson", "Spearman")


# violin plot of difference scores
graph <- ggplot(output, aes(x=group, y=RSFC_change, fill = group)) + 
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.1) +
  geom_jitter(size = 1, shape=1, position=position_jitter(0.2)) +
  theme_classic() +
  labs(x="Experimental Condition", y="Sum score", title = "Difference in HPC-VTA RSFC") +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20, face="bold"), title=element_text(size =20, face="bold")) +
  coord_cartesian(ylim = c(-1, 1)) + 
  facet_grid(. ~ coefficient)
print(graph)
print(paste0("Violinplot_HPCVTA_RSFC_changes.jpeg")) # save plot
ggsave(paste0("Violinplot_HPCVTA_RSFC_changes.jpeg")) # save plot




########## 10. Create table showing the correlation between RSFC changes and memory measurements ########## 

dependentVariablesWide <- c("cuedRecallStrict", "cuedRecallLenient",
                            "recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
                            "rememberedStrict", "rememberedLenient", 
                            #betas
                            "curiosityBeta_cuedRecallStrict", "curiosityBeta_cuedRecallLenient",
                            "curiosityBeta_allConf", "curiosityBeta_highConf", "curiosityBeta_aboveAvgConf",
                            "curiosityBeta_rememberedStrict", "curiosityBeta_rememberedLenient",
                            # benefit
                            "curiosityBenefit_cuedRecallStrict", "curiosityBenefit_cuedRecallStrict_dichotom",
                            "curiosityBenefit_cuedRecallLenient", "curiosityBenefit_cuedRecallLenient_dichotom",
                            "curiosityBenefit_allConf", "curiosityBenefit_allConf_dichotom", 
                            "curiosityBenefit_highConf", "curiosityBenefit_highConf_dichotom",
                            "curiosityBenefit_aboveAvgConf", "curiosityBenefit_aboveAvgConf_dichotom", 
                            "curiosityBenefit_rememberedStrict", "curiosityBenefit_rememberedStrict_dichotom",
                            "curiosityBenefit_rememberedLenient", "curiosityBenefit_rememberedLenient_dichotom",
                            # correlations
                            "curiosityCorrelation_cuedRecallStrict", "curiosityCorrelation_cuedRecallLenient", 
                            "curiosityCorrelation_allConf", "curiosityCorrelation_highConf", "curiosityCorrelation_aboveAvgConf", 
                            "curiosityCorrelation_rememberedStrict", "curiosityCorrelation_rememberedLenient")


# for all dependent variables compute correlation between changes in RSFC and memory for whole sample as well as in each group
df_cor <- data.frame(var   = character(length(dependentVariablesWide)),
                     pearson = numeric(length(dependentVariablesWide)),
                     pearson_noR = numeric(length(dependentVariablesWide)),
                     pearson_R = numeric(length(dependentVariablesWide)),
                     spearman = numeric(length(dependentVariablesWide)),
                     spearman_noR = numeric(length(dependentVariablesWide)),
                     spearman_R = numeric(length(dependentVariablesWide)),
                     stringsAsFactors=FALSE)
for(DV in 1:length(dependentVariablesWide)) {
  # paste variable
  df_cor$var[DV] <- paste(dependentVariablesWide[DV])
  # compute correlarion with FC computed with pearson
  df_cor$pearson[DV] <- cor(dfWide[[paste(dependentVariablesWide[DV])]], dfWide$RSFC_VTAHPC_diff)
  df_cor$pearson_noR[DV] <- cor(dfWide[[paste(dependentVariablesWide[DV])]][dfWide$group=="int"], dfWide$RSFC_VTAHPC_diff[dfWide$group=="int"])
  df_cor$pearson_R[DV] <- cor(dfWide[[paste(dependentVariablesWide[DV])]][dfWide$group=="ext"], dfWide$RSFC_VTAHPC_diff[dfWide$group=="ext"])
  # compute correlarion with FC computed with spearman
  df_cor$spearman[DV] <- cor(dfWide[[paste(dependentVariablesWide[DV])]], dfWide$RSFC_VTAHPC_diff_spearman)
  df_cor$spearman_noR[DV] <- cor(dfWide[[paste(dependentVariablesWide[DV])]][dfWide$group=="int"], dfWide$RSFC_VTAHPC_diff[dfWide$group=="int"])
  df_cor$spearman_R[DV] <- cor(dfWide[[paste(dependentVariablesWide[DV])]][dfWide$group=="ext"], dfWide$RSFC_VTAHPC_diff[dfWide$group=="ext"])
}
# 
setwd(brainDir)
write.csv(df_cor, file = "Correlation_RSFC_changes_memory.csv", row.names = F)

########## 11. Create scatter plot to show relationship between HPC-VTA RSFC and memory measurements ########## 

# select variables to be plotted
output <- dfWide[,c("ID", "group", dependentVariablesWide, "RSFC_VTAHPC_diff", "RSFC_VTAHPC_diff_spearman")]
output$group <- ifelse(output$group == "int", "No reward", "Reward")

## create long data frame for: 
# cuedRecall: lenient vs strict
output_recall_p <- output[,c("ID", "group", "cuedRecallStrict", "cuedRecallLenient", "RSFC_VTAHPC_diff")]
output_recall_p <- reshape2::melt(output_recall_p, id=c("ID","group", "RSFC_VTAHPC_diff"))
names(output_recall_p) <- c("ID", "group", "RSFC_VTAHPC_diff_p",  "criteria", "performance")

output_recall_s <- output[,c("ID", "group", "cuedRecallStrict", "cuedRecallLenient", "RSFC_VTAHPC_diff_spearman")]
output_recall_s <- reshape2::melt(output_recall_s, id=c("ID","group", "RSFC_VTAHPC_diff_spearman"))
names(output_recall_s) <- c("ID", "group", "RSFC_VTAHPC_diff_s",  "criteria", "performance")

output_recall <- merge(output_recall_p, output_recall_s, by = c("ID", "group", "criteria", "performance"))
rm(output_recall_p, output_recall_s)
output_recall <- reshape2::melt(output_recall, id=c("ID","group", "criteria", "performance"))
names(output_recall) <- c("ID", "group", "criteria", "performance", "coefficient", "RSFC_change")
output_recall$coefficient <- ifelse(output_recall$coefficient == "RSFC_VTAHPC_diff_p", "Pearson", "Spearman")

# recognition
output_recognition <- output[,c("ID", "group", "recognition", "RSFC_VTAHPC_diff", "RSFC_VTAHPC_diff_spearman")]
output_recognition <- reshape2::melt(output_recognition, id=c("ID","group", "recognition"))
names(output_recognition) <-  c("ID", "group", "performance",  "coefficient", "RSFC_change")
output_recognition$coefficient <- ifelse(output_recognition$coefficient == "RSFC_VTAHPC_diff", "Pearson", "Spearman")

#recollection: high vs aboveAvg
output_recollection_p <- output[,c("ID", "group", "recognitionAboveMeanConf", "recognitionConfLevel_4_5_6", "RSFC_VTAHPC_diff")]
output_recollection_p <- reshape2::melt(output_recollection_p, id=c("ID","group", "RSFC_VTAHPC_diff"))
names(output_recollection_p) <- c("ID", "group", "RSFC_VTAHPC_diff_p",  "criteria", "performance")

output_recollection_s <- output[,c("ID", "group", "recognitionAboveMeanConf", "recognitionConfLevel_4_5_6", "RSFC_VTAHPC_diff_spearman")]
output_recollection_s <- reshape2::melt(output_recollection_s, id=c("ID","group", "RSFC_VTAHPC_diff_spearman"))
names(output_recollection_s) <- c("ID", "group", "RSFC_VTAHPC_diff_s",  "criteria", "performance")

output_recollection <- merge(output_recollection_p, output_recollection_s, by = c("ID", "group", "criteria", "performance"))
rm(output_recollection_p, output_recollection_s)
output_recollection <- reshape2::melt(output_recollection, id=c("ID","group", "criteria", "performance"))
names(output_recollection) <- c("ID", "group", "criteria", "performance", "coefficient", "RSFC_change")
output_recollection$coefficient <- ifelse(output_recollection$coefficient == "RSFC_VTAHPC_diff_p", "Pearson", "Spearman")

# remembered: lenient vs strict
output_remembered_p <- output[,c("ID", "group", "rememberedStrict", "rememberedLenient", "RSFC_VTAHPC_diff")]
output_remembered_p <- reshape2::melt(output_remembered_p, id=c("ID","group", "RSFC_VTAHPC_diff"))
names(output_remembered_p) <- c("ID", "group", "RSFC_VTAHPC_diff_p",  "criteria", "performance")

output_remembered_s <- output[,c("ID", "group", "rememberedStrict", "rememberedLenient", "RSFC_VTAHPC_diff_spearman")]
output_remembered_s <- reshape2::melt(output_remembered_s, id=c("ID","group", "RSFC_VTAHPC_diff_spearman"))
names(output_remembered_s) <- c("ID", "group", "RSFC_VTAHPC_diff_s",  "criteria", "performance")

output_remembered <- merge(output_remembered_p, output_remembered_s, by = c("ID", "group", "criteria", "performance"))
rm(output_remembered_p, output_remembered_s)
output_remembered <- reshape2::melt(output_remembered, id=c("ID","group", "criteria", "performance"))
names(output_remembered) <- c("ID", "group", "criteria", "performance", "coefficient", "RSFC_change")
output_remembered$coefficient <- ifelse(output_remembered$coefficient == "RSFC_VTAHPC_diff_p", "Pearson", "Spearman")

# beta recall: lenient vs strict
output_betaRecall_p <- output[,c("ID", "group", "curiosityBeta_cuedRecallStrict", "curiosityBeta_cuedRecallLenient", "RSFC_VTAHPC_diff")]
output_betaRecall_p <- reshape2::melt(output_betaRecall_p, id=c("ID","group", "RSFC_VTAHPC_diff"))
names(output_betaRecall_p) <- c("ID", "group", "RSFC_VTAHPC_diff_p",  "criteria", "performance")

output_betaRecall_s <- output[,c("ID", "group", "curiosityBeta_cuedRecallStrict", "curiosityBeta_cuedRecallLenient", "RSFC_VTAHPC_diff_spearman")]
output_betaRecall_s <- reshape2::melt(output_betaRecall_s, id=c("ID","group", "RSFC_VTAHPC_diff_spearman"))
names(output_betaRecall_s) <- c("ID", "group", "RSFC_VTAHPC_diff_s",  "criteria", "performance")

output_betaRecall <- merge(output_betaRecall_p, output_betaRecall_s, by = c("ID", "group", "criteria", "performance"))
rm(output_betaRecall_p, output_betaRecall_s)
output_betaRecall <- reshape2::melt(output_betaRecall, id=c("ID","group", "criteria", "performance"))
names(output_betaRecall) <- c("ID", "group", "criteria", "performance", "coefficient", "RSFC_change")
output_betaRecall$coefficient <- ifelse(output_betaRecall$coefficient == "RSFC_VTAHPC_diff_p", "Pearson", "Spearman")

# beta recognition
output_betaRecognition <- output[,c("ID", "group", "curiosityBeta_allConf", "RSFC_VTAHPC_diff", "RSFC_VTAHPC_diff_spearman")]
output_betaRecognition <- reshape2::melt(output_betaRecognition, id=c("ID","group", "curiosityBeta_allConf"))
names(output_betaRecognition) <-  c("ID", "group", "performance",  "coefficient", "RSFC_change")
output_betaRecognition$coefficient <- ifelse(output_betaRecognition$coefficient == "RSFC_VTAHPC_diff", "Pearson", "Spearman")

# beta recollection: high vs aboveAvg
output_betaRecollection_p <- output[,c("ID", "group", "curiosityBeta_aboveAvgConf", "curiosityBeta_highConf", "RSFC_VTAHPC_diff")]
output_betaRecollection_p <- reshape2::melt(output_betaRecollection_p, id=c("ID","group", "RSFC_VTAHPC_diff"))
names(output_betaRecollection_p) <- c("ID", "group", "RSFC_VTAHPC_diff_p",  "criteria", "performance")

output_betaRecollection_s <- output[,c("ID", "group", "curiosityBeta_aboveAvgConf", "curiosityBeta_highConf", "RSFC_VTAHPC_diff_spearman")]
output_betaRecollection_s <- reshape2::melt(output_betaRecollection_s, id=c("ID","group", "RSFC_VTAHPC_diff_spearman"))
names(output_betaRecollection_s) <- c("ID", "group", "RSFC_VTAHPC_diff_s",  "criteria", "performance")

output_betaRecollection <- merge(output_betaRecollection_p, output_betaRecollection_s, by = c("ID", "group", "criteria", "performance"))
rm(output_betaRecollection_p, output_betaRecollection_s)
output_betaRecollection <- reshape2::melt(output_betaRecollection, id=c("ID","group", "criteria", "performance"))
names(output_betaRecollection) <- c("ID", "group", "criteria", "performance", "coefficient", "RSFC_change")
output_betaRecollection$coefficient <- ifelse(output_betaRecollection$coefficient == "RSFC_VTAHPC_diff_p", "Pearson", "Spearman")

# beta remembered: lenient vs strict
output_betaRemembered_p <- output[,c("ID", "group", "curiosityBeta_rememberedStrict", "curiosityBeta_rememberedLenient", "RSFC_VTAHPC_diff")]
output_betaRemembered_p <- reshape2::melt(output_betaRemembered_p, id=c("ID","group", "RSFC_VTAHPC_diff"))
names(output_betaRemembered_p) <- c("ID", "group", "RSFC_VTAHPC_diff_p",  "criteria", "performance")

output_betaRemembered_s <- output[,c("ID", "group", "curiosityBeta_rememberedStrict", "curiosityBeta_rememberedLenient", "RSFC_VTAHPC_diff_spearman")]
output_betaRemembered_s <- reshape2::melt(output_betaRemembered_s, id=c("ID","group", "RSFC_VTAHPC_diff_spearman"))
names(output_betaRemembered_s) <- c("ID", "group", "RSFC_VTAHPC_diff_s",  "criteria", "performance")

output_betaRemembered <- merge(output_betaRemembered_p, output_betaRemembered_s, by = c("ID", "group", "criteria", "performance"))
rm(output_betaRemembered_p, output_betaRemembered_s)
output_betaRemembered <- reshape2::melt(output_betaRemembered, id=c("ID","group", "criteria", "performance"))
names(output_betaRemembered) <- c("ID", "group", "criteria", "performance", "coefficient", "RSFC_change")
output_betaRemembered$coefficient <- ifelse(output_betaRemembered$coefficient == "RSFC_VTAHPC_diff_p", "Pearson", "Spearman")

# correlation recall: lenient vs strict
output_correlationRecall_p <- output[,c("ID", "group", "curiosityCorrelation_cuedRecallStrict", "curiosityCorrelation_cuedRecallLenient", "RSFC_VTAHPC_diff")]
output_correlationRecall_p <- reshape2::melt(output_correlationRecall_p, id=c("ID","group", "RSFC_VTAHPC_diff"))
names(output_correlationRecall_p) <- c("ID", "group", "RSFC_VTAHPC_diff_p",  "criteria", "performance")

output_correlationRecall_s <- output[,c("ID", "group", "curiosityCorrelation_cuedRecallStrict", "curiosityCorrelation_cuedRecallLenient", "RSFC_VTAHPC_diff_spearman")]
output_correlationRecall_s <- reshape2::melt(output_correlationRecall_s, id=c("ID","group", "RSFC_VTAHPC_diff_spearman"))
names(output_correlationRecall_s) <- c("ID", "group", "RSFC_VTAHPC_diff_s",  "criteria", "performance")

output_correlationRecall <- merge(output_correlationRecall_p, output_correlationRecall_s, by = c("ID", "group", "criteria", "performance"))
rm(output_correlationRecall_p, output_correlationRecall_s)
output_correlationRecall <- reshape2::melt(output_correlationRecall, id=c("ID","group", "criteria", "performance"))
names(output_correlationRecall) <- c("ID", "group", "criteria", "performance", "coefficient", "RSFC_change")
output_correlationRecall$coefficient <- ifelse(output_correlationRecall$coefficient == "RSFC_VTAHPC_diff_p", "Pearson", "Spearman")

# correlation recognition
output_correlationRecognition <- output[,c("ID", "group", "curiosityCorrelation_allConf", "RSFC_VTAHPC_diff", "RSFC_VTAHPC_diff_spearman")]
output_correlationRecognition <- reshape2::melt(output_correlationRecognition, id=c("ID","group", "curiosityCorrelation_allConf"))
names(output_correlationRecognition) <-  c("ID", "group", "performance",  "coefficient", "RSFC_change")
output_correlationRecognition$coefficient <- ifelse(output_correlationRecognition$coefficient == "RSFC_VTAHPC_diff", "Pearson", "Spearman")

# correlation recollection: high vs aboveAvg
output_correlationRecollection_p <- output[,c("ID", "group", "curiosityCorrelation_aboveAvgConf", "curiosityCorrelation_highConf", "RSFC_VTAHPC_diff")]
output_correlationRecollection_p <- reshape2::melt(output_correlationRecollection_p, id=c("ID","group", "RSFC_VTAHPC_diff"))
names(output_correlationRecollection_p) <- c("ID", "group", "RSFC_VTAHPC_diff_p",  "criteria", "performance")

output_correlationRecollection_s <- output[,c("ID", "group", "curiosityCorrelation_aboveAvgConf", "curiosityCorrelation_highConf", "RSFC_VTAHPC_diff_spearman")]
output_correlationRecollection_s <- reshape2::melt(output_correlationRecollection_s, id=c("ID","group", "RSFC_VTAHPC_diff_spearman"))
names(output_correlationRecollection_s) <- c("ID", "group", "RSFC_VTAHPC_diff_s",  "criteria", "performance")

output_correlationRecollection <- merge(output_correlationRecollection_p, output_correlationRecollection_s, by = c("ID", "group", "criteria", "performance"))
rm(output_correlationRecollection_p, output_correlationRecollection_s)
output_correlationRecollection <- reshape2::melt(output_correlationRecollection, id=c("ID","group", "criteria", "performance"))
names(output_correlationRecollection) <- c("ID", "group", "criteria", "performance", "coefficient", "RSFC_change")
output_correlationRecollection$coefficient <- ifelse(output_correlationRecollection$coefficient == "RSFC_VTAHPC_diff_p", "Pearson", "Spearman")

# correlation remembered: lenient vs strict
output_correlationRemembered_p <- output[,c("ID", "group", "curiosityCorrelation_rememberedStrict", "curiosityCorrelation_rememberedLenient", "RSFC_VTAHPC_diff")]
output_correlationRemembered_p <- reshape2::melt(output_correlationRemembered_p, id=c("ID","group", "RSFC_VTAHPC_diff"))
names(output_correlationRemembered_p) <- c("ID", "group", "RSFC_VTAHPC_diff_p",  "criteria", "performance")

output_correlationRemembered_s <- output[,c("ID", "group", "curiosityCorrelation_rememberedStrict", "curiosityCorrelation_rememberedLenient", "RSFC_VTAHPC_diff_spearman")]
output_correlationRemembered_s <- reshape2::melt(output_correlationRemembered_s, id=c("ID","group", "RSFC_VTAHPC_diff_spearman"))
names(output_correlationRemembered_s) <- c("ID", "group", "RSFC_VTAHPC_diff_s",  "criteria", "performance")

output_correlationRemembered <- merge(output_correlationRemembered_p, output_correlationRemembered_s, by = c("ID", "group", "criteria", "performance"))
rm(output_correlationRemembered_p, output_correlationRemembered_s)
output_correlationRemembered <- reshape2::melt(output_correlationRemembered, id=c("ID","group", "criteria", "performance"))
names(output_correlationRemembered) <- c("ID", "group", "criteria", "performance", "coefficient", "RSFC_change")
output_correlationRemembered$coefficient <- ifelse(output_correlationRemembered$coefficient == "RSFC_VTAHPC_diff_p", "Pearson", "Spearman")

# benefit recall: 2 (lenient vs strict) * 2 (cont vs dichotom) * 2 (corr: pearson vs spearman)
output_benefitRecall_p <- output[,c("ID", "group","curiosityBenefit_cuedRecallStrict", "curiosityBenefit_cuedRecallStrict_dichotom",
                                    "curiosityBenefit_cuedRecallLenient", "curiosityBenefit_cuedRecallLenient_dichotom", "RSFC_VTAHPC_diff")]
output_benefitRecall_p <- reshape2::melt(output_benefitRecall_p, id=c("ID","group", "RSFC_VTAHPC_diff"))
names(output_benefitRecall_p) <- c("ID", "group", "Pearson", "criteria", "performance")
output_benefitRecall_p$method <- ifelse(output_benefitRecall_p$criteria == "curiosityBenefit_cuedRecallStrict" | output_benefitRecall_p$criteria == "curiosityBenefit_cuedRecallLenient", "continuous", "dichotom")
output_benefitRecall_p$criteria <- ifelse(output_benefitRecall_p$criteria == "curiosityBenefit_cuedRecallStrict" | output_benefitRecall_p$criteria == "curiosityBenefit_cuedRecallStrict_dichotom", "cuedRecallStrict", "cuedRecallLenient")

output_benefitRecall_s <- output[,c("ID", "group","curiosityBenefit_cuedRecallStrict", "curiosityBenefit_cuedRecallStrict_dichotom",
                                    "curiosityBenefit_cuedRecallLenient", "curiosityBenefit_cuedRecallLenient_dichotom", "RSFC_VTAHPC_diff_spearman")]
output_benefitRecall_s <- reshape2::melt(output_benefitRecall_s, id=c("ID","group", "RSFC_VTAHPC_diff_spearman"))
names(output_benefitRecall_s) <- c("ID", "group", "Spearman", "criteria", "performance")
output_benefitRecall_s$method <- ifelse(output_benefitRecall_s$criteria == "curiosityBenefit_cuedRecallStrict" | output_benefitRecall_s$criteria == "curiosityBenefit_cuedRecallLenient", "continuous", "dichotom")
output_benefitRecall_s$criteria <- ifelse(output_benefitRecall_s$criteria == "curiosityBenefit_cuedRecallStrict" | output_benefitRecall_s$criteria == "curiosityBenefit_cuedRecallStrict_dichotom", "cuedRecallStrict", "cuedRecallLenient")

output_benefitRecall <- merge(output_benefitRecall_p, output_benefitRecall_s, by = c("ID", "group", "criteria", "performance", "method"))
rm(output_benefitRecall_p, output_benefitRecall_s)
output_benefitRecall <- melt(output_benefitRecall, id=c("ID", "group", "criteria", "performance", "method"))
names(output_benefitRecall) <- c("ID", "group", "criteria", "performance", "method", "coefficient", "RSFC_change")

# benefit recall: 2 (cont vs dichotom) * 2 (corr: pearson vs spearman)
output_benefitRecall_p <- output[,c("ID", "group","curiosityBenefit_cuedRecallStrict", "curiosityBenefit_cuedRecallStrict_dichotom",
                                    "curiosityBenefit_cuedRecallLenient", "curiosityBenefit_cuedRecallLenient_dichotom", "RSFC_VTAHPC_diff")]
output_benefitRecall_p <- reshape2::melt(output_benefitRecall_p, id=c("ID","group", "RSFC_VTAHPC_diff"))
names(output_benefitRecall_p) <- c("ID", "group", "Pearson", "criteria", "performance")
output_benefitRecall_p$method <- ifelse(output_benefitRecall_p$criteria == "curiosityBenefit_cuedRecallStrict" | output_benefitRecall_p$criteria == "curiosityBenefit_cuedRecallLenient", "continuous", "dichotom")
output_benefitRecall_p$criteria <- ifelse(output_benefitRecall_p$criteria == "curiosityBenefit_cuedRecallStrict" | output_benefitRecall_p$criteria == "curiosityBenefit_cuedRecallStrict_dichotom", "cuedRecallStrict", "cuedRecallLenient")

output_benefitRecall_s <- output[,c("ID", "group","curiosityBenefit_cuedRecallStrict", "curiosityBenefit_cuedRecallStrict_dichotom",
                                    "curiosityBenefit_cuedRecallLenient", "curiosityBenefit_cuedRecallLenient_dichotom", "RSFC_VTAHPC_diff_spearman")]
output_benefitRecall_s <- reshape2::melt(output_benefitRecall_s, id=c("ID","group", "RSFC_VTAHPC_diff_spearman"))
names(output_benefitRecall_s) <- c("ID", "group", "Spearman", "criteria", "performance")
output_benefitRecall_s$method <- ifelse(output_benefitRecall_s$criteria == "curiosityBenefit_cuedRecallStrict" | output_benefitRecall_s$criteria == "curiosityBenefit_cuedRecallLenient", "continuous", "dichotom")
output_benefitRecall_s$criteria <- ifelse(output_benefitRecall_s$criteria == "curiosityBenefit_cuedRecallStrict" | output_benefitRecall_s$criteria == "curiosityBenefit_cuedRecallStrict_dichotom", "cuedRecallStrict", "cuedRecallLenient")

output_benefitRecall <- merge(output_benefitRecall_p, output_benefitRecall_s, by = c("ID", "group", "criteria", "performance", "method"))
rm(output_benefitRecall_p, output_benefitRecall_s)
output_benefitRecall <- melt(output_benefitRecall, id=c("ID", "group", "criteria", "performance", "method"))
names(output_benefitRecall) <- c("ID", "group", "criteria", "performance", "method", "coefficient", "RSFC_change")

# benefit recognition: 2 (cont vs dichotom) * 2 (corr: pearson vs spearman)
output_benefitRecognition_p <- output[,c("ID", "group","curiosityBenefit_allConf", "curiosityBenefit_allConf_dichotom", "RSFC_VTAHPC_diff")]
output_benefitRecognition_p <- reshape2::melt(output_benefitRecognition_p, id=c("ID","group", "RSFC_VTAHPC_diff"))
names(output_benefitRecognition_p) <- c("ID", "group", "Pearson", "criteria", "performance")
output_benefitRecognition_p$method <- ifelse(output_benefitRecognition_p$criteria == "curiosityBenefit_allConf", "continuous", "dichotom")
output_benefitRecognition_p$criteria <- "curiosityBenefit_allConf"

output_benefitRecognition_s <- output[,c("ID", "group","curiosityBenefit_allConf", "curiosityBenefit_allConf_dichotom", "RSFC_VTAHPC_diff_spearman")]
output_benefitRecognition_s <- reshape2::melt(output_benefitRecognition_s, id=c("ID","group", "RSFC_VTAHPC_diff_spearman"))
names(output_benefitRecognition_s) <- c("ID", "group", "Spearman", "criteria", "performance")
output_benefitRecognition_s$method <- ifelse(output_benefitRecognition_s$criteria == "curiosityBenefit_allConf", "continuous", "dichotom")
output_benefitRecognition_s$criteria <- "curiosityBenefit_allConf"

output_benefitRecognition <- merge(output_benefitRecognition_p, output_benefitRecognition_s, by = c("ID", "group", "criteria", "performance", "method"))
rm(output_benefitRecognition_p, output_benefitRecognition_s)
output_benefitRecognition <- melt(output_benefitRecognition, id=c("ID", "group", "criteria", "performance", "method"))
names(output_benefitRecognition) <- c("ID", "group", "criteria", "performance", "method", "coefficient", "RSFC_change")

# #benefit recollection: 2 (high vs aboveAvg) * 2 (cont vs dichotom) * 2 (corr: pearson vs spearman)
output_benefitRecollection_p <- output[,c("ID", "group","curiosityBenefit_aboveAvgConf", "curiosityBenefit_aboveAvgConf_dichotom",
                                          "curiosityBenefit_highConf", "curiosityBenefit_highConf_dichotom", "RSFC_VTAHPC_diff")]
output_benefitRecollection_p <- reshape2::melt(output_benefitRecollection_p, id=c("ID","group", "RSFC_VTAHPC_diff"))
names(output_benefitRecollection_p) <- c("ID", "group", "Pearson", "criteria", "performance")
output_benefitRecollection_p$method <- ifelse(output_benefitRecollection_p$criteria == "curiosityBenefit_aboveAvgConf" | output_benefitRecollection_p$criteria == "curiosityBenefit_highConf", "continuous", "dichotom")
output_benefitRecollection_p$criteria <- ifelse(output_benefitRecollection_p$criteria == "curiosityBenefit_aboveAvgConf" | output_benefitRecollection_p$criteria == "curiosityBenefit_aboveAvgConf_dichotom", "aboveAvgConf", "highConf")

output_benefitRecollection_s <- output[,c("ID", "group","curiosityBenefit_aboveAvgConf", "curiosityBenefit_aboveAvgConf_dichotom",
                                          "curiosityBenefit_highConf", "curiosityBenefit_highConf_dichotom", "RSFC_VTAHPC_diff_spearman")]
output_benefitRecollection_s <- reshape2::melt(output_benefitRecollection_s, id=c("ID","group", "RSFC_VTAHPC_diff_spearman"))
names(output_benefitRecollection_s) <- c("ID", "group", "Spearman", "criteria", "performance")
output_benefitRecollection_s$method <- ifelse(output_benefitRecollection_s$criteria == "curiosityBenefit_aboveAvgConf" | output_benefitRecollection_s$criteria == "curiosityBenefit_highConf", "continuous", "dichotom")
output_benefitRecollection_s$criteria <- ifelse(output_benefitRecollection_s$criteria == "curiosityBenefit_aboveAvgConf" | output_benefitRecollection_s$criteria == "curiosityBenefit_aboveAvgConf_dichotom", "aboveAvgConf", "highConf")

output_benefitRecollection <- merge(output_benefitRecollection_p, output_benefitRecollection_s, by = c("ID", "group", "criteria", "performance", "method"))
rm(output_benefitRecollection_p, output_benefitRecollection_s)
output_benefitRecollection <- melt(output_benefitRecollection, id=c("ID", "group", "criteria", "performance", "method"))
names(output_benefitRecollection) <- c("ID", "group", "criteria", "performance", "method", "coefficient", "RSFC_change")


# benefit remembered: 2 (cont vs dichotom) * 2 (corr: pearson vs spearman)
output_benefitRemembered_p <- output[,c("ID", "group","curiosityBenefit_rememberedStrict", "curiosityBenefit_rememberedStrict_dichotom",
                                        "curiosityBenefit_rememberedLenient", "curiosityBenefit_rememberedLenient_dichotom", "RSFC_VTAHPC_diff")]
output_benefitRemembered_p <- reshape2::melt(output_benefitRemembered_p, id=c("ID","group", "RSFC_VTAHPC_diff"))
names(output_benefitRemembered_p) <- c("ID", "group", "Pearson", "criteria", "performance")
output_benefitRemembered_p$method <- ifelse(output_benefitRemembered_p$criteria == "curiosityBenefit_rememberedStrict" | output_benefitRemembered_p$criteria == "curiosityBenefit_rememberedLenient", "continuous", "dichotom")
output_benefitRemembered_p$criteria <- ifelse(output_benefitRemembered_p$criteria == "curiosityBenefit_rememberedStrict" | output_benefitRemembered_p$criteria == "curiosityBenefit_rememberedStrict_dichotom", "rememberedStrict", "rememberedLenient")

output_benefitRemembered_s <- output[,c("ID", "group","curiosityBenefit_rememberedStrict", "curiosityBenefit_rememberedStrict_dichotom",
                                        "curiosityBenefit_rememberedLenient", "curiosityBenefit_rememberedLenient_dichotom", "RSFC_VTAHPC_diff_spearman")]
output_benefitRemembered_s <- reshape2::melt(output_benefitRemembered_s, id=c("ID","group", "RSFC_VTAHPC_diff_spearman"))
names(output_benefitRemembered_s) <- c("ID", "group", "Spearman", "criteria", "performance")
output_benefitRemembered_s$method <- ifelse(output_benefitRemembered_s$criteria == "curiosityBenefit_rememberedStrict" | output_benefitRemembered_s$criteria == "curiosityBenefit_rememberedLenient", "continuous", "dichotom")
output_benefitRemembered_s$criteria <- ifelse(output_benefitRemembered_s$criteria == "curiosityBenefit_rememberedStrict" | output_benefitRemembered_s$criteria == "curiosityBenefit_rememberedStrict_dichotom", "rememberedStrict", "rememberedLenient")

output_benefitRemembered <- merge(output_benefitRemembered_p, output_benefitRemembered_s, by = c("ID", "group", "criteria", "performance", "method"))
rm(output_benefitRemembered_p, output_benefitRemembered_s)
output_benefitRemembered <- melt(output_benefitRemembered, id=c("ID", "group", "criteria", "performance", "method"))
names(output_benefitRemembered) <- c("ID", "group", "criteria", "performance", "method", "coefficient", "RSFC_change")

# define which variables to plot
plotVars <- c("recall", "recognition","recollection", "betaRecall", "betaRecognition", "betaRecollection", "benefitRecall", "benefitRecollection", "benefitRecognition")
plotVars <- c("recall", "recollection", "remembered", "betaRecall", "betaRecollection", "betaRemembered",
              "correlationRecall", "correlationRecollection", "correlationRemembered",
              "recognition", "betaRecognition", "correlationRecognition",
              "benefitRecall", "benefitRecollection", "benefitRecognition", "benefitRemembered")

for (plot in plotVars){
  #
  graph <- ggplot(get(paste0("output_",plot)), aes(x=RSFC_change, y=performance, col = group)) + 
    geom_point() +
    geom_smooth(method=lm, aes(fill=group), alpha = 0.3, fullrange = T) +
    theme_classic() +
    labs(x="Difference in HPC-VTA RSFC (post > pre)", y="Performance index", title = paste("RSFC changes and", plot), fill = "Group", col  = "Group") +
    theme(axis.text=element_text(size=20), axis.title=element_text(size=20, face="bold"), title=element_text(size =20, face="bold")) #+
  # change facet grid depending on dependent variable
  if (plot %in% c("recall", "recollection", "remembered", "betaRecall", "betaRecollection", "betaRemembered",
                  "correlationRecall", "correlationRecollection", "correlationRemembered")){
    graph <- graph + facet_grid(coefficient ~ criteria) 
  } else if (plot %in% c("recognition", "betaRecognition", "correlationRecognition")) {
    graph <- graph + facet_grid(coefficient ~ .) 
  }
  #change facet grid and coordinate system
  if (plot %in% c("benefitRecall", "benefitRecognition", "benefitRecollection", "benefitRemembered")){
    graph <- graph + facet_grid(coefficient ~ criteria + method)
    graph <- graph + coord_cartesian(ylim = c(-18, 18))#, xlim = c(-0.5, 0.5)) 
  }
  # change coordinate system depending on dependent var
  if (plot %in% c("recall", "recollection", "recognition")){
    graph <- graph + coord_cartesian(ylim = c(0, 36))#, xlim = c(-0.5, 0.5)) 
  } 
  print(graph) 
  print(paste0("Scatterplot_", plot, ".jpeg"))
  ggsave(paste0("Scatterplot_", plot, ".jpeg"))
}  


