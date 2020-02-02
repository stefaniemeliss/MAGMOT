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
questDir <- file.path(ratingsDir, "Questionnaires")
tricksDir <- file.path(analysisDir, "Tricks")

# check whether these directories exist, if not create them
ifelse(!dir.exists(ratingsDir), dir.create(ratingsDir), FALSE) 
ifelse(!dir.exists(memoryDir), dir.create(memoryDir), FALSE) 
ifelse(!dir.exists(questDir), dir.create(questDir), FALSE) 

#helper functions and packages #
source("~/Dropbox/Reading/Codes and functions/R/errorbars.R")
source("~/Dropbox/Reading/Codes and functions/R/rbindcolumns.R")
library(ggplot2)

### read in data sets ###
setwd(preprocessedDir)
dfWide <- xlsx::read.xlsx("wide_MAGMOT.xlsx", sheetName = "Sheet1")
dfLong <- xlsx::read.xlsx("long_MAGMOT.xlsx", sheetName = "Sheet1")

#### data set handling wide format ####
scales <- names(dfWide[,c("BIS", "BAS_rewardresponsiveness", "BAS_drive", "BAS_funseeking", "NeedForCognition", "FearOfFailure", "ApproachTemperament", "AvoidanceTemperament", "TraitCuriosity", "StateCuriosity", "intrinsicMotivation", "taskEngagement", "interest", "boredom", "effort", "pressure", "corsiSpan", "nback_hitrate", "nback_falsealarmrate", "nback_accurary")])
psych::describe(dfWide[,scales])
by(cbind(dfWide[,scales]), dfWide$group, psych::describe)


### demogs ####
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

#### plots post questionnaire ####
setwd(ratingsDir)
rm(output)

names(dfWide) [names(dfWide) == "post23"] <- "compliance"
names(dfWide) [names(dfWide) == "post24"] <- "ableToSee"

# combine different measurements / scales
IMI <- names(dfWide[,c("intrinsicMotivation", "taskEngagement", "interest", "boredom", "effort", "pressure")])
BISBAS <- names(dfWide[,c("BIS", "BAS_rewardresponsiveness", "BAS_drive", "BAS_funseeking")])
MCI <- names(dfWide[,c("TraitCuriosity", "StateCuriosity")])
others <- names(dfWide[,c("NeedForCognition", "FearOfFailure", "ApproachTemperament", "AvoidanceTemperament")])
postQest <- names(dfWide[,c("compliance", "ableToSee", "memoryTestKnown", "memoryIntention")])
neuro <- names(dfWide[,c("corsiSpan", "nback_hitrate", "nback_falsealarmrate", "nback_accurary")])


allQ <- c("IMI", "BISBAS", "MCI", "others", "postQest", "neuro")
setwd(questDir)
# plot average responses for each scale seperately for each group
for (q in seq_along(allQ)){
  output <- by(cbind(dfWide[,get(allQ[q])]), dfWide$group, psych::describe)
  rating <- as.data.frame(rbind(output$int, output$ext))
  rating$mot <- rep(c("intrinsic","extrinsic"), each = length(get(allQ[q])))
  
  outg <- ggplot(rating, aes(vars, mean, fill = mot))
  outg <- outg + geom_bar(stat="identity", position="dodge") + 
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9))  + 
    labs(x="Scales", y="Rating", fill = "Experimental Condition", title = paste(allQ[q], "questionnaire",version)) + 
    theme_classic() + #coord_cartesian(ylim = c(0, 7)) +
    scale_x_discrete(limits=get(allQ[q])) +
    theme(legend.position="bottom")
  print(outg)
  ggsave(paste0("QuestionnaireDataByGroup_", allQ[q],".jpeg")) # save plot
}

rm(output, rating, outg)

# effect sizes and t-tests
for(scale in 1:length(scales)) {
  # compute t-test for group difference
  ttest <- t.test(dfWide[,scales[scale]]~dfWide$group)
  t.stats <- as.data.frame(t(round(c(ttest$statistic, ttest$p.value), digits = 3)))
  names(t.stats) <- c("tValue", "pValue")
  
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
rm(data, cohen, means, t.stats, ttest, d)
write.csv(effectsizesScales, paste0("effectsizesScalesNeuro_", version, ".csv"))


###### lme #####
if (exists("dfLong") == F) {
  setwd(preprocessedDir)
  dfLong <- xlsx::read.xlsx("long_MAGMOT.xlsx", sheetName = "Sheet1")#,  showWarnings = FALSE)
}

# define grouping variables
groupingVariables <- c("mediansplitCuriosity_MAGMOT", "mediansplitCuriosityWithinSubject", "curiosity_dichotom")

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
                            "confidence", "confidenceCorrectTrials")
}


#### get descriptives for the whole sample as well as for each group individually
for(DV in 1:length(dependentVariables)) {
  print(dependentVariables[DV])
  # all
  descriptive <- psych::describe(dfLong[,dependentVariables[DV] ])
  row.names(descriptive) <- c(paste0(dependentVariables[DV], "_all"))
  #
  descriptive_groupwise <-   by(cbind(dfLong[,dependentVariables[DV]]), dfLong$group, psych::describe)
  descriptive_groupwise_named <- rbind(descriptive_groupwise[[1]], descriptive_groupwise[[2]])
  row.names(descriptive_groupwise_named) <- c(paste0(dependentVariables[DV], "_", names(descriptive_groupwise)[1]),paste0(dependentVariables[DV], "_", names(descriptive_groupwise)[2]) )
  
  
  if (DV == 1) {
    descriptives <- rbind(descriptive, descriptive_groupwise_named)
  } else {
    temp_descriptives <- rbind(descriptive, descriptive_groupwise_named)
    descriptives <- rbind.all.columns(descriptives, temp_descriptives) #rbind all columns will induce NA if there was initially no data saved in the loop per participant
    rm(temp_descriptives)
  }

  rm(descriptive, descriptive_groupwise, descriptive_groupwise_named)
}
descriptives <- round(descriptives, digit = 3)
setwd(ratingsDir)
xlsx::write.xlsx(descriptives, file="Descriptives_dependentVariables.xlsx", sheetName = "Sheet1")
rm(descriptives)

# lmer model using curiosity as a continous variable
workspace <- list.files(path = file.path(codedDir), pattern = "_CP.csv") # check whether the data is coded yet or not
if(length(workspace) == 0) { # if data is not coded yet, only look at recognition performance
  dependentVariables <- c("recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
                          "confidence", "confidenceCorrectTrials")
} else {
  dependentVariables <- c("cuedRecallLenient", "cuedRecallStrict",
                          "recognition","recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
                          "confidence", "confidenceCorrectTrials",
                          "confidenceGroupMeanCentered", "confidenceGroupMeanCenteredCorrectTrials")
}

library(lmerTest)

# loop over dependent variables
for (DV in 1:length(dependentVariables)){
  # if DV continuous
  if (max(dfLong[, dependentVariables[DV]], na.rm = T) > 1){
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


# look at the effects controlling for nback [same as above, just different model]
for (DV in 1:length(dependentVariables)){
  
  if (max(dfLong[, dependentVariables[DV]], na.rm = T) > 1){
    LMEmodel_curiosityContinuousControlled <- lmerTest::lmer(dfLong[, dependentVariables[DV]] ~ nback_accurary + groupEffectCoded*curiosityGroupMeanCentered + (1+curiosityGroupMeanCentered|ID) + (1|stimID), data = dfLong)
  }else{
    LMEmodel_curiosityContinuousControlled <- glmer(dfLong[, dependentVariables[DV]] ~ nback_accurary + groupEffectCoded*curiosityGroupMeanCentered + (1+curiosityGroupMeanCentered|ID) + (1|stimID), family = "binomial"(link = 'logit'), data = dfLong)
  }
  summaryCuriosityContinuousControlled <- summary(LMEmodel_curiosityContinuousControlled)
  summaryCuriosityContinuousCoefficientsControlled <- as.data.frame(summaryCuriosityContinuousControlled$coefficients)
  
  row.names(summaryCuriosityContinuousCoefficientsControlled) <- c(paste0("LME_",dependentVariables[DV], "_Intercept"), paste0("LME_",dependentVariables[DV], "_nBackAccuracy"), 
                                                                   paste0("LME_",dependentVariables[DV], "_groupEffectCoded"), 
                                                         paste0("LME_",dependentVariables[DV], "_curiosityGroupMeanCentered"), 
                                                         paste0("LME_",dependentVariables[DV], "_interaction"))
  
  if (DV == 1) {
    LMEresultsControlled <- summaryCuriosityContinuousCoefficientsControlled
  }else {
    temp_LMEresults <- summaryCuriosityContinuousCoefficientsControlled
    LMEresultsControlled <- rbind.all.columns(LMEresultsControlled, temp_LMEresults) #rbind all columns will induce NA if there was initially no data saved in the loop per participant
    rm(temp_LMEresults)
  }
  rm(summaryCuriosityContinuousCoefficientsControlled, summaryCuriosityContinuousControlled)
}
LMEresultsControlled <- round(LMEresultsControlled, digits = 5)
xlsx::write.xlsx(LMEresultsControlled, file=paste0("LME_Results_", version, "_CuriosityRatingsAndMemoryScores_byCuriosityAsContinuousVariable_controlledForNback.xlsx"), sheetName = "Sheet1")


#### bar plots ####

dfLong$curiosity_dichotom <- ifelse(dfLong$curiosity_dichotom == -1, "low", 
                                    ifelse(dfLong$curiosity_dichotom == 1, "high", NA)) # create curiosity_dichotom as factor
groupingVariables <- c("mediansplitCuriosity_MAGMOT", "mediansplitCuriosityWithinSubject", "curiosity_dichotom")
groupingVariables <- c("mediansplitCuriosityWithinSubject", "curiosity_dichotom")

if (exists("dfLong") == F) {
  setwd(preprocessedDir)
  dfLong <- xlsx::read.xlsx("long_MagicBehavioural.xlsx", sheetName = "Sheet1")#,  showWarnings = FALSE)
} 
# determine dependent variables to loop over
dependentVariables <- c("cuedRecallLenient", "cuedRecallStrict",
                          "recognition","recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
                          "confidence", "confidenceCorrectTrials",
                          "confidenceGroupMeanCentered", "confidenceGroupMeanCenteredCorrectTrials")

cols <- c("above" = "#F8766D", "below" = "#00BFC4", "all tricks" = "grey", "median" = "#C77CFF",
          "high" = "#F8766D", "low" = "#00BFC4")

setwd(memoryDir)
for (g in 1:length(groupingVariables)){

  for (DV in dependentVariables){
    # create a data frame containing the data for each group divided regarding curiosity ratings
    outputGroup <- summarySEwithin(dfLong, measurevar=DV, betweenvars="group", withinvars=groupingVariables[g], idvar="ID", na.rm = T)
    levels(outputGroup$group) <- c("Reward","No reward")
    
    # remove NAs (necessary for curiosity_dichotomous)
    outputGroup <- na.omit(outputGroup) 
    
    # create a data frame containing the data for each group regardless of curiosity ratings
    outputAll <- summarySE(dfLong, measurevar=DV, groupvars="group", na.rm=T,
                           conf.interval=.95, .drop=TRUE)
    levels(outputAll$group) <- c("Reward","No reward")
    outputAll[[paste0(groupingVariables[g])]] <- rep("all tricks", 2)
    
    # combine these two data frames and use it for plotting
    output <- rbind.all.columns(outputGroup, outputAll)
    output <- outputGroup # this line can be commented out to also illustrate effect regardless of group
    
    # create bar graph
    graph <- ggplot(output, aes(x=group, y=get(DV), fill=get(groupingVariables[g]))) +
      geom_bar(stat="identity", position="dodge") + geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=get(DV)-se, ymax=get(DV)+se)) +
      scale_x_discrete(limits=c("No reward", "Reward")) + 
      labs(x="Between Group manipulation", y="Performance index", fill = "Curiosity split", title = paste(version, ":", DV, "by",groupingVariables[g] ))  +
      theme(axis.text=element_text(size=20), axis.title=element_text(size=20, face="bold"), title=element_text(size =20, face="bold"), legend.title = element_text(size=20), legend.text = element_text(size = 20)) +
      theme_classic() + scale_fill_manual(values = cols)
    # change x y-axis depending on dependent variable
    if (DV %in% c("cuedRecallStrict", "cuedRecallLenient")){
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
    print(paste0("Bargraph_2x2_", DV, "_by", groupingVariables[g], ".jpeg"))
    ggsave(paste0("Bargraph_2x2_", DV, "_by", groupingVariables[g], ".jpeg"))
  }
}

###### create histograms ######
# define variables
dependentVariables <- c("cuedRecallLenient", "cuedRecallStrict",
                        "recognition","recognitionConfLevel_4_5_6", "recognitionAboveMeanConf")
# create data frame and recode memory as factors
output <- dfLong[,c("ID", "group", "curiosityGroupMeanCentered", dependentVariables)]
output$group <- ifelse(output$group == "int", "No reward", "Reward")
output$cuedRecallLenient <- ifelse(output$cuedRecallLenient == 0, "forgotten", "remembered")
output$cuedRecallStrict <- ifelse(output$cuedRecallStrict == 1, "remembered", "forgotten")
output$recognition <- ifelse(output$recognition == 1, "remembered", "forgotten")
output$recognitionConfLevel_4_5_6 <- ifelse(output$recognitionConfLevel_4_5_6 == 1, "remembered", "forgotten")
output$recognitionAboveMeanConf <- ifelse(output$recognitionAboveMeanConf == 1, "remembered", "forgotten")

for (DV in dependentVariables){
  # histogram plot
  graph <- ggplot(output, aes(x=curiosityGroupMeanCentered, col = get(DV), fill = get(DV))) + 
    geom_histogram(aes(y=..density..), binwidth = 0.1, alpha=0.2, position="dodge") +
    geom_density(alpha=.1) +
    theme_classic() + 
    scale_color_brewer(palette="Dark2", limits = c("remembered", "forgotten")) + scale_fill_brewer(palette="Dark2", limits = c("remembered", "forgotten")) +
    facet_grid(group ~ .) +
    coord_cartesian(xlim = c(-5, 5), ylim = c(0, 1)) +
    labs(x="curiosity group mean centered", y="Density", title = paste(version, DV)) +
    theme(legend.position="bottom") +
    theme(axis.text=element_text(size=14), axis.title=element_text(size=16, face="bold"), title=element_text(size =20, face="bold"), strip.text = element_text(size = 16)) 
  
  print(graph)
  print(paste0("Histogram_", DV, ".jpeg"))
  ggsave(paste0("Histogram_", DV, ".jpeg"))
}


#### determine effect sizes ####
source("~/Dropbox/Reading/Codes and functions/R/rbindcolumns.R")

dependentVariablesWide <- c("cuedRecallStrict", "cuedRecallLenient",
                            "recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
                            "meanConfidence", "meanConfidenceCorrectTrials",
                            "curiosityBeta_cuedRecallStrict", "curiosityBeta_cuedRecallLenient",
                            "curiosityBeta_allConf", "curiosityBeta_highConf", "curiosityBeta_aboveAvgConf",
                            "curiosityBenefit_cuedRecallStrict", "curiosityBenefit_cuedRecallStrict_dichotom",
                            "curiosityBenefit_allConf", "curiosityBenefit_allConf_dichotom", 
                            "curiosityBenefit_highConf", "curiosityBenefit_highConf_dichotom",
                            "curiosityBenefit_aboveAvgConf", "curiosityBenefit_aboveAvgConf_dichotom")

for(DV in 1:length(dependentVariablesWide)) {
  
  # compute t-test for group difference
  ttest <- t.test(dfWide[,dependentVariablesWide[DV]]~dfWide$group)
  t.stats <- as.data.frame(t(round(c(ttest$statistic, ttest$p.value), digits = 3)))
  names(t.stats) <- c("tValue", "pValue")
  
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
rm(data, cohen, means, t.stats, ttest, d)
effectsizesMemory <- round(effectsizesMemory, digits = 3)
write.csv(effectsizesMemory, paste0("effectsizesMemory_", version, ".csv"))
effectsizes <- rbind.all.columns(effectsizesScales, effectsizesMemory)
setwd(ratingsDir)
write.csv(effectsizes, paste0("effectsizesAll_", version, ".csv"))

##### create violin plots ####
output <- dfWide[,c("ID", "group", dependentVariablesWide)]
output$group <- ifelse(output$group == "int", "No reward", "Reward")

setwd(memoryDir)
for (DV in dependentVariablesWide){
  
  # Basic violin plot
  graph <- ggplot(output, aes(x=group, y=get(DV), fill = group)) + 
    geom_violin(trim=FALSE) +
    geom_boxplot(width=0.1) +
    geom_jitter(size = 1, shape=1, position=position_jitter(0.2)) +
    theme_classic() +
    labs(x="Experimental Condition", y="Sum score", title = paste(version, DV)) +
    theme(legend.position="none") +
    theme(axis.text=element_text(size=20), axis.title=element_text(size=20, face="bold"), title=element_text(size =20, face="bold")) 
  # change x y-axis depending on dependent variable
  if (DV %in% c("cuedRecallStrict", "cuedRecallLenient", "recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf")){
    graph <- graph + coord_cartesian(ylim = c(-5, 41))
  }
  if (DV %in% c("confidence", "confidenceCorrectTrials")){
    graph <- graph + coord_cartesian(ylim = c(0, 6))
  }
  if (DV %in% c("curiosityBeta_cuedRecallStrict", "curiosityBeta_cuedRecallLenient", "curiosityBeta_allConf", "curiosityBeta_highConf", "curiosityBeta_aboveAvgConf")){
    graph <- graph + coord_cartesian(ylim = c(-.25, 0.25))
  }
  if (DV %in% c("curiosityBenefit_cuedRecallStrict", "curiosityBenefit_cuedRecallStrict_dichotom", 
                "curiosityBenefit_allConf", "curiosityBenefit_allConf_dichotom",  
                "curiosityBenefit_highConf", "curiosityBenefit_highConf_dichotom",        
                "curiosityBenefit_aboveAvgConf", "curiosityBenefit_aboveAvgConf_dichotom")){
    graph <- graph + coord_cartesian(ylim = c(-23, 23))
  }
  
  print(graph)
  print(paste0("Violonplot_", DV, ".jpeg"))
  ggsave(paste0("Violonplot_", DV, ".jpeg"))
}


# ################### PROBLEM, for future stef to fix ################### 
# 
# # check whether there is a difference in the effect size depending on task block
# # NOTE: SOMETHING IS NOT WORKING HERE!!!!! HOWEVER, the effect sizes do not seem to change over blocks, hence tiredness can be excluded
# blocks <- c(1,2,3)
# blockstring <- c("_firstBlock", "_secondBlock", "_thirdBlock")
# dependentVariablesWide <- c("cuedRecallStrict", "cuedRecallLenient",
#                             "recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
#                             "meanConfidence", "meanConfidenceCorrectTrials")
# 
# for (block in blocks){
#   for(DV in 1:length(dependentVariablesWide)) {
#     # compute mean for each block and group
#     means <- tapply(dfWide[[paste0(dependentVariablesWide[DV], blockstring[block])]], dfWide$group, mean, na.rm = T)
#     means <- as.data.frame(t(means))
#     
#     # calculate effect size
#     if (means$ext != means$int) {
#       data <- dfWide[,c("group", paste0(dependentVariablesWide[DV], blockstring[block]))]
#       d <- psych::cohen.d(data, "group")
#       cohen <- as.data.frame(d$cohen.d)
#       means <- merge(cohen, means)
#       rm(cohen, d)
#     }
#     
#     # combine data from all DVs
#     row.names(means) <- paste(dependentVariablesWide[DV])
#     
#     if (DV == 1) {
#       effectsizesMemoryBlock <- means
#     } else {
#       temp_effectsizesMemory <- means
#       effectsizesMemoryBlock <- rbind.all.columns(effectsizesMemoryBlock, temp_effectsizesMemory) #rbind all columns will induce NA if there was initially no data saved in the loop per participant
#       rm(temp_effectsizesMemory)
#     }
#   }
#   # change column names
#   names(effectsizesMemoryBlock)[names(effectsizesMemoryBlock)=="lower"] <- c(paste0("lower", blockstring[block]))
#   names(effectsizesMemoryBlock)[names(effectsizesMemoryBlock)=="effect"] <- c(paste0("effect", blockstring[block]))
#   names(effectsizesMemoryBlock)[names(effectsizesMemoryBlock)=="upper"] <- c(paste0("upper", blockstring[block]))
#   names(effectsizesMemoryBlock)[names(effectsizesMemoryBlock)=="int"] <- c(paste0("int", blockstring[block]))
#   names(effectsizesMemoryBlock)[names(effectsizesMemoryBlock)=="ext"] <- c(paste0("ext", blockstring[block]))
#   
#   assign(paste0("effectsizesMemoryBlock",block), effectsizesMemoryBlock)
#   rm(effectsizesMemoryBlock)
#   
# }
# 
# effectsizesMemoryBlock <- merge(effectsizesMemoryBlock1, effectsizesMemoryBlock2, by = "row.names")
# row.names(effectsizesMemoryBlock) <- effectsizesMemoryBlock$Row.names
# effectsizesMemoryBlock$Row.names <- NULL
# effectsizesMemoryBlock <- merge(effectsizesMemoryBlock, effectsizesMemoryBlock3, by = "row.names")
# row.names(effectsizesMemoryBlock) <- effectsizesMemoryBlock$Row.names
# effectsizesMemoryBlock$Row.names <- NULL
# 
# effectsizesMemoryBlock[,c("effect_firstBlock", "effect_secondBlock", "effect_thirdBlock", "ext_firstBlock", "ext_secondBlock", "ext_thirdBlock", "int_firstBlock", "int_secondBlock", "int_thirdBlock")]
# write.csv(effectsizesMemoryBlock, paste0("effectsizesMemoryPerBlock_", version, ".csv"))
# # effectsizesMemoryBlock<- effectsizesMemoryBlock[,c("effect_firstBlock", "effect_secondBlock", "effect_thirdBlock")]
# # write.csv(effectsizesMemoryBlock, paste0("effectsizesMemoryPerBlock_", version, ".csv"))
# 
# # spaghetti plot
# library(gridExtra)
# setwd(ratingsDir)
# output <- effectsizesMemoryBlock[,c("effect_firstBlock", "effect_secondBlock", "effect_thirdBlock", "ext_firstBlock", "ext_secondBlock", "ext_thirdBlock", "int_firstBlock", "int_secondBlock", "int_thirdBlock")]
# 
# output <- as.data.frame(t(output))
# output$rowname <-   row.names(output)
# output$block <- rep(c("first", "second", "third"), 3)
# output$whatIsMeasured <- rep(c("effectsize", "meanExpGroup", "meanContGroup"), each = 3)
# output$whatIsMeasured2 <- rep(c("effectsize", "group mean", "group mean"), each = 3)
# 
# 
# for (DV in 1:length(dependentVariablesWide)){
#     graph <- ggplot(data=output, aes(x=block, y=get(dependentVariablesWide[DV] ), colour = whatIsMeasured, group = whatIsMeasured)) +
#       geom_point() + geom_line() + theme_bw() +
#       labs(x="Task block", y="", colour = "", title = paste(dependentVariablesWide[DV])) +
#       facet_grid(rows = vars(whatIsMeasured2), scales = "free_y")  + coord_cartesian(ylim = c(-1, 8))
#     print(graph)
#     ggsave(paste0("Graph_changeOfEffectsizeOverBlocks_", dependentVariablesWide[DV], ".jpeg"))
# }
# 
# outputEffect <- subset(output, output$whatIsMeasured2 == "effectsize")
# outputGroupMean <- subset(output, output$whatIsMeasured2 == "group mean")
# 
# for (DV in 1:length(dependentVariablesWide)){
#   # create one graph for the effect size (ranging from -1 to 1)
#   graphEffect <- ggplot(data=outputEffect, aes(x=block, y=get(dependentVariablesWide[DV] ), group = whatIsMeasured2)) +
#     geom_point(color='black') + geom_line(color='black') + theme_bw() +
#     theme(legend.position="bottom") +
#     labs(x="Task block", y="Cohen's d", group = "", title = paste("Effect size",dependentVariablesWide[DV], version)) + 
#     coord_cartesian(ylim = c(-1, 1)) +
#     scale_y_continuous(labels = scales::scientific)
#   
#   # then create a graph for the group average
#   graphGroupMean <- ggplot(data=outputGroupMean, aes(x=block, y=get(dependentVariablesWide[DV] ), colour = whatIsMeasured, group = whatIsMeasured)) +
#     geom_point() + geom_line() + theme_bw() +
#     labs(x="Task block", y="Mean value", colour = "", title = paste("Group means",dependentVariablesWide[DV], version)) + 
#     coord_cartesian(ylim = c(floor(min(outputGroupMean[,dependentVariablesWide[DV]]))-1, floor(min(outputGroupMean[,dependentVariablesWide[DV]]))+3)) +
#     theme(legend.position="bottom") +
#     scale_y_continuous(labels = scales::scientific)
# 
#   # then combine those two in a new graphic using gridExtra
#   gridExtra::grid.arrange(graphEffect, graphGroupMean, nrow = 2) 
#   graph <- gridExtra::arrangeGrob(graphEffect, graphGroupMean, nrow=2) #generates graph
#   ggsave(file = paste0("Graph_changeOfEffectsizeOverBlocks_", dependentVariablesWide[DV], "_v2.jpeg"), graph)
#   
#   rm(graphEffect, graphGroupMean, graph)
# }