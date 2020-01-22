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
codedDir <- file.path(dataDir, "coded", "preprocessed") #"~/Dropbox/Reading/PhD/Magic tricks/behavioural_study/data_fin/coded/preprocessing"
analysisDir <- file.path(mainDir, "Analysis")
ratingsDir <- file.path(analysisDir, "Ratings")
tricksDir <- file.path(analysisDir, "Tricks")

# check whether these directories exist, if not create them
ifelse(!dir.exists(ratingsDir), dir.create(ratingsDir), FALSE) 

#helper functions and packages #

source("~/Dropbox/Reading/Codes and functions/R/errorbars.R")
source("~/Dropbox/Reading/Codes and functions/R/rbindcolumns.R")
# library(xlsx)
# library(dplyr)
# library(psych)
library(ggplot2)

### read in data sets ###
setwd(preprocessedDir)
dfWide <- xlsx::read.xlsx("wide_MAGMOT.xlsx", sheetName = "Sheet1")
#dfLong <- read.xlsx("long_MagicBehavioural.xlsx", sheetName = "Sheet1")

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

IMI <- names(dfWide[,c("intrinsicMotivation", "taskEngagement", "interest", "boredom", "effort", "pressure")])
BISBAS <- names(dfWide[,c("BIS", "BAS_rewardresponsiveness", "BAS_drive", "BAS_funseeking")])
MCI <- names(dfWide[,c("TraitCuriosity", "StateCuriosity")])
others <- names(dfWide[,c("NeedForCognition", "FearOfFailure", "ApproachTemperament", "AvoidanceTemperament")])
postQest <- names(dfWide[,c("compliance", "ableToSee", "memoryTestKnown", "memoryIntention")])
neuro <- names(dfWide[,c("corsiSpan", "nback_hitrate", "nback_falsealarmrate", "nback_accurary")])


allQ <- c("IMI", "BISBAS", "MCI", "others", "postQest", "neuro")

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
  ggsave(paste0("QuestionnaireDataByGroup_", allQ[q],".jpeg"))
}

rm(output, rating, outg)

# effect sizes and t-tests
for(scale in 1:length(scales)) {
  #print(scales[scale])
  ttest <- t.test(dfWide[,scales[scale]]~dfWide$group)
  #print(ttest)
  t.stats <- as.data.frame(t(round(c(ttest$statistic, ttest$p.value), digits = 3)))
  names(t.stats) <- c("tValue", "pValue")
  
  means <- tapply(dfWide[,scales[scale]], dfWide$group, mean, na.rm = T)
  means <- as.data.frame(t(means))
  means <- merge(t.stats, means)
  #print(means)
  
  if (means$ext != means$int) {
    data <- dfWide[,c("group", scales[scale])]
    
    psych::cohen.d(data, "group")
    
    d <- psych::cohen.d(data, "group")
    cohen <- as.data.frame(d$cohen.d)
    #print(cohen)
    
    means <- merge(means, cohen)
  }
  
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
groupingVariables <- c("mediansplitCuriosity_MAGMOT", "mediansplitCuriosityWithinSubject")
#groupingVariables <- c("mediansplitCuriosity_MAGMOT")


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


#### get descriptives
for(DV in 1:length(dependentVariables)) {
  print(dependentVariables[DV])
  descriptive <- psych::describe(dfLong[,dependentVariables[DV] ])
  row.names(descriptive) <- c(paste0(dependentVariables[DV], "_all"))
  
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

if(length(workspace) == 0) { # if data is not coded yet, only look at recognition performance
  dependentVariables <- c("recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
                          "confidence", "confidenceCorrectTrials")
} else {
  dependentVariables <- c("cuedRecallLenient", "cuedRecallStrict",
                          "recognition","recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
                          "confidence", "confidenceCorrectTrials")
}

library(lmerTest)


DV <- 5 # look at confidence for correct trials
DV <- 2 # look at high confidence recognition
DV <- 1 # look at recognition
LMEmodel_curiosityContinuous <- glmer(dfLong[, dependentVariables[DV]] ~ 1 + groupEffectCoded + curiosityGroupMeanCentered + rewardByCuriosity + (1 | ID), data = dfLong, family = binomial(link = 'logit'))


LMEmodel_curiosityContinuous_recogAboveAvgConf <- glmer(recognitionAboveMeanConf ~ 1 + groupEffectCoded + curiosityGroupMeanCentered + rewardByCuriosity + (1 | ID), data = dfLong, family = binomial(link = 'logit'))
summary(LMEmodel_curiosityContinuous_recogAboveAvgConf)
exp(2*LMEmodel_curiosityContinuous_recogAboveAvgConf@beta[2]) # transform beta estimate for log-Odds from effect coding interpretation to actual change
exp(LMEmodel_curiosityContinuous_recogAboveAvgConf@beta[2]) # transform beta estimate for log-Odds from effect coding interpretation to actual change
sjPlot::tab_model(LMEmodel_curiosityContinuous_recogAboveAvgConf)
sjPlot::tab_model(LMEmodel_curiosityContinuous_recogAboveAvgConf, transform = NULL)


LMEmodel_curiosityContinuous_recogHighConf <- glmer(recognitionConfLevel_4_5_6 ~ 1 + groupEffectCoded + curiosityGroupMeanCentered + rewardByCuriosity + (1 | ID), data = dfLong, family = binomial(link = 'logit'))
summary(LMEmodel_curiosityContinuous_recogHighConf)
exp(2*LMEmodel_curiosityContinuous_recogHighConf@beta[2]) # transform beta estimate for log-Odds from effect coding interpretation to actual change
exp(LMEmodel_curiosityContinuous_recogHighConf@beta[2]) # transform beta estimate for log-Odds from effect coding interpretation to actual change
sjPlot::tab_model(LMEmodel_curiosityContinuous_recogHighConf)
sjPlot::tab_model(LMEmodel_curiosityContinuous_recogHighConf, transform = NULL)

LMEmodel_curiosityContinuous_recog <- glmer(recognition ~ 1 + groupEffectCoded + curiosityGroupMeanCentered + rewardByCuriosity + (1 | ID), data = dfLong, family = binomial(link = 'logit'))
summary(LMEmodel_curiosityContinuous_recog)
sjPlot::tab_model(LMEmodel_curiosityContinuous_recog)
sjPlot::tab_model(LMEmodel_curiosityContinuous_recog, transform = NULL)

LMEmodel_curiosityContinuous_confidenceCorrect <- lmerTest::lmer(confidenceCorrectTrials ~ 1 + groupEffectCoded + curiosityGroupMeanCentered + rewardByCuriosity + (1 | ID), data = dfLong)
summary(LMEmodel_curiosityContinuous_confidenceCorrect)
sjPlot::tab_model(LMEmodel_curiosityContinuous_confidenceCorrect)
sjPlot::tab_model(LMEmodel_curiosityContinuous_confidenceCorrect, transform = NULL)

LMEmodel_curiosityContinuous_confidence <- lmerTest::lmer(confidence ~ 1 + groupEffectCoded + curiosityGroupMeanCentered + rewardByCuriosity + (1 | ID), data = dfLong)
summary(LMEmodel_curiosityContinuous_confidence)
sjPlot::tab_model(LMEmodel_curiosityContinuous_confidence)
sjPlot::tab_model(LMEmodel_curiosityContinuous_confidence, transform = NULL)

LMEmodel_curiosityContinuous_confidence_recogAsCov <- lmerTest::lmer(confidence ~ 1 + recognition + groupEffectCoded + curiosityGroupMeanCentered + rewardByCuriosity + (1 | ID), data = dfLong)
summary(LMEmodel_curiosityContinuous_confidence_recogAsCov)
sjPlot::tab_model(LMEmodel_curiosityContinuous_confidence_recogAsCov)
sjPlot::tab_model(LMEmodel_curiosityContinuous_confidence_recogAsCov, transform = NULL, show.stat = T, show.std = T)
sjstats::std_beta(LMEmodel_curiosityContinuous_confidence_recogAsCov)



LMEmodel_recognition <- glmer(recognition ~ 1 + confidenceGroupMeanCentered +  groupEffectCoded*curiosityGroupMeanCentered + (1 | ID), data = dfLong, family = binomial(link = 'logit'))
summary(LMEmodel_recognition)
sjPlot::tab_model(LMEmodel_recognition)

LMEmodel_highConfRec <- glmer(recognitionConfLevel_4_5_6 ~ 1 + confidenceGroupMeanCentered + groupEffectCoded*curiosityGroupMeanCentered + (1 | ID), data = dfLong, family = binomial(link = 'logit'))
summary(LMEmodel_highConfRec)
sjPlot::tab_model(LMEmodel_highConfRec)


cor.test(dfLong$curiosityGroupMeanCentered, dfLong$confidenceGroupMeanCentered, use = "pairwise.complete.obs")
plot(dfLong$curiosityGroupMeanCentered, dfLong$confidenceGroupMeanCentered, main = paste(version))
abline(lm(dfLong$curiosityGroupMeanCentered ~ dfLong$confidenceGroupMeanCentered), col="red")


for (DV in 1:length(dependentVariables)){
  
  if (max(dfLong[, dependentVariables[DV]], na.rm = T) > 1){
    LMEmodel_curiosityContinuous <- lmerTest::lmer(dfLong[, dependentVariables[DV]] ~ groupEffectCoded*curiosityGroupMeanCentered + (1+curiosityGroupMeanCentered|ID) + (1|stimID), data = dfLong)
  }else{
    LMEmodel_curiosityContinuous <- glmer(dfLong[, dependentVariables[DV]] ~ groupEffectCoded*curiosityGroupMeanCentered + (1+curiosityGroupMeanCentered|ID) + (1|stimID), family = "binomial"(link = 'logit'), data = dfLong)
  }
  summaryCuriosityContinuous <- summary(LMEmodel_curiosityContinuous)
  summaryCuriosityContinuousCoefficients <- as.data.frame(summaryCuriosityContinuous$coefficients)
  
  row.names(summaryCuriosityContinuousCoefficients) <- c(paste0("LME_",dependentVariables[DV], "_Intercept"), paste0("LME_",dependentVariables[DV], "_groupEffectCoded"), 
                                                         paste0("LME_",dependentVariables[DV], "_curiosityGroupMeanCentered"), 
                                                         paste0("LME_",dependentVariables[DV], "_interaction"))
  
  if (DV == 1) {
    LMEresults <- summaryCuriosityContinuousCoefficients
  }else {
    temp_LMEresults <- summaryCuriosityContinuousCoefficients
    LMEresults <- rbind.all.columns(LMEresults, temp_LMEresults) #rbind all columns will induce NA if there was initially no data saved in the loop per participant
    rm(temp_LMEresults)
  }
  rm(summaryCuriosityContinuousCoefficients, summaryCuriosityContinuous)
}
LMEresults <- round(LMEresults, digits = 5)

setwd(ratingsDir)
xlsx::write.xlsx(LMEresults, file=paste0("LME_Results_", version, "_CuriosityRatingsAndMemoryScores_byCuriosityAsContinuousVariable.xlsx"), sheetName = "Sheet1")

# look at the effects controlling for nback 
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
setwd(preprocessedDir)
if (exists("groupingVariables") == F) {
  groupingVariables <- c("mediansplitCuriosity_MAGMOT", "mediansplitCuriosityWithinSubject")
}
if (exists("dfLong") == F) {
  setwd(preprocessedDir)
  dfLong <- xlsx::read.xlsx("long_MagicBehavioural.xlsx", sheetName = "Sheet1")#,  showWarnings = FALSE)
}
if (exists("dependentVariables") == F) {
  dependentVariables <- c("recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
                          "confidence", "confidenceCorrectTrials")
}

cols <- c("above" = "#F8766D", "below" = "#00BFC4", "all tricks" = "grey", "median" = "#C77CFF")

setwd(ratingsDir)
for (g in 1:length(groupingVariables)){

  for (DV in dependentVariables){
    # create a data frame containing the data for each group divided regarding curiosity ratings
    outputGroup <- summarySEwithin(dfLong, measurevar=DV, betweenvars="group", withinvars=groupingVariables[g], idvar="ID", na.rm = T)
    levels(outputGroup$group) <- c("Reward","No reward")
    
    # create a data frame containing the data for each group regardless of curiosity ratings
    outputAll <- summarySE(dfLong, measurevar=DV, groupvars="group", na.rm=T,
                           conf.interval=.95, .drop=TRUE)
    levels(outputAll$group) <- c("Reward","No reward")
    outputAll[[paste0(groupingVariables[g])]] <- rep("all tricks", 2)
    
    # combine these two data frames and use it for plotting
    output <- rbind.all.columns(outputGroup, outputAll)
    
    output <- outputGroup
    
    # create bar graph
    graph <- ggplot(output, aes(x=group, y=get(DV), fill=get(groupingVariables[g]))) +
      geom_bar(stat="identity", position="dodge") + geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=get(DV)-se, ymax=get(DV)+se)) +
      scale_x_discrete(limits=c("No reward", "Reward")) + 
      labs(x="Between Group manipulation", y="Performance index", fill = "Curiosity median split", title = paste("Memory performance in",version, ":", DV, "by",groupingVariables[g] ))  +
      theme(axis.text=element_text(size=20), axis.title=element_text(size=20, face="bold"), title=element_text(size =20, face="bold"), legend.title = element_text(size=20), legend.text = element_text(size = 20)) +
      theme_classic() + scale_fill_manual(values = cols)
    if (DV %in% c("cuedRecallStrict", "cuedRecallLenient")){
      graph <- graph + coord_cartesian(ylim = c(0, 1))
    }
    if (DV %in% c("recognition", "recognitionConfLevel_1", "recognitionConfLevel_2", "recognitionConfLevel_3", "recognitionConfLevel_4", "recognitionConfLevel_5", "recognitionConfLevel_6",
                  "recognitionConfLevel_1_2", "recognitionConfLevel_3_4", "recognitionConfLevel_5_6", "recognitionConfLevel_1_2_3", "recognitionConfLevel_4_5_6")){
      graph <- graph + coord_cartesian(ylim = c(0, 1)) + geom_hline(yintercept = 0.25, linetype="dashed", color = "black")
    }
    if (DV %in% c("confidence", "confidenceCorrectTrials")){
      graph <- graph + coord_cartesian(ylim = c(0, 6))
      
    }
    print(graph)
    print(paste0("Graph_2x2_", DV, "_by", groupingVariables[g], ".jpeg"))
    #ggsave(paste0("Graph_2x2_", DV, "_by", groupingVariables[g], "_v2.jpeg"))
    ggsave(paste0("Graph_2x2_", DV, "_by", groupingVariables[g], ".jpeg"))
  }
}

tapply(dfLong$recognition, list(dfLong$group, dfLong$mediansplitCuriosity_MAGMOT), mean)
tapply(dfLong$recognitionConfLevel_4_5_6, list(dfLong$group, dfLong$mediansplitCuriosity_MAGMOT), mean)
tapply(dfLong$recognitionConfLevel_4_5_6, list(dfLong$group, as.factor(dfLong$responseCuriosity)), mean)


## determine effect sizes
source("~/Dropbox/Reading/Codes and functions/R/rbindcolumns.R")

dependentVariablesWide <- c("cuedRecallStrict", "cuedRecallLenient",
                        "recognition", 
                        "recognitionConfLevel_6", "recognitionConfLevel_5_6", "recognitionConfLevel_4_5_6")
dependentVariablesWide <- c("recognition", "meanConfidence", "meanConfidenceCorrectTrials",
                            "recognitionConfLevel_6", "recognitionConfLevel_5_6", "recognitionConfLevel_4_5_6")
dependentVariablesWide <- c("recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
                            "meanConfidence", "meanConfidenceCorrectTrials")

for(DV in 1:length(dependentVariablesWide)) {
  
  
  ttest <- t.test(dfWide[,dependentVariablesWide[DV]]~dfWide$group)
  #print(ttest)
  t.stats <- as.data.frame(t(round(c(ttest$statistic, ttest$p.value), digits = 3)))
  names(t.stats) <- c("tValue", "pValue")
  
  
  means <- tapply(dfWide[,dependentVariablesWide[DV]], dfWide$group, mean, na.rm = T)
  # deviations <- tapply(dfWide[,dependentVariablesWide[DV]], dfWide$group, sd, na.rm = T)
  means <- as.data.frame(t(means))
  # deviations <- as.data.frame(t(deviations))
  means <- merge(t.stats, means)
  
  if (means$ext != means$int) {
    
    data <- dfWide[,c("group", dependentVariablesWide[DV])]
    
    d <- psych::cohen.d(data, "group")
    cohen <- as.data.frame(d$cohen.d)
    
    means <- merge(means, cohen)
    
  }
  row.names(means) <- paste(dependentVariablesWide[DV])
  if (DV == 1) {
    effectsizesMemory <- means
  } else {
    temp_effectsizesMemory <- means
    effectsizesMemory <- rbind.all.columns(effectsizesMemory, temp_effectsizesMemory) #rbind all columns will induce NA if there was initially no data saved in the loop per participant
    rm(temp_effectsizesMemory)
  }
}
rm(data, cohen, means, t.stats, ttest, d)

effectsizesMemory <- round(effectsizesMemory, digits = 3)

write.csv(effectsizesMemory, paste0("effectsizesMemory_", version, ".csv"))

effectsizes <- rbind.all.columns(effectsizesScales, effectsizesMemory)
write.csv(effectsizes, paste0("effectsizesAll_", version, ".csv"))


# check whether there is a difference in the effect size depending on task block
blocks <- c(1,2,3)
blockstring <- c("_firstBlock", "_secondBlock", "_thirdBlock")

for (block in blocks){
  for(DV in 1:length(dependentVariablesWide)) {
    means <- tapply(dfWide[[paste0(dependentVariablesWide[DV], blockstring[block])]], dfWide$group, mean, na.rm = T)
    print(paste0(dependentVariablesWide[DV], blockstring[block]))
    print(means)
    means <- as.data.frame(t(means))
    
    if (means$ext != means$int) {
      
      data <- dfWide[,c("group", paste0(dependentVariablesWide[DV], blockstring[block]))]
      d <- psych::cohen.d(data, "group")
      cohen <- as.data.frame(d$cohen.d)
      means <- merge(cohen, means)
      rm(cohen, d)
    }
    
    row.names(means) <- paste(dependentVariablesWide[DV])
    
    if (DV == 1) {
      effectsizesMemoryBlock <- means
    } else {
      temp_effectsizesMemory <- means
      effectsizesMemoryBlock <- rbind.all.columns(effectsizesMemoryBlock, temp_effectsizesMemory) #rbind all columns will induce NA if there was initially no data saved in the loop per participant
      rm(temp_effectsizesMemory)
    }
  }
  names(effectsizesMemoryBlock)[names(effectsizesMemoryBlock)=="lower"] <- c(paste0("lower", blockstring[block]))
  names(effectsizesMemoryBlock)[names(effectsizesMemoryBlock)=="effect"] <- c(paste0("effect", blockstring[block]))
  names(effectsizesMemoryBlock)[names(effectsizesMemoryBlock)=="upper"] <- c(paste0("upper", blockstring[block]))
  names(effectsizesMemoryBlock)[names(effectsizesMemoryBlock)=="int"] <- c(paste0("int", blockstring[block]))
  names(effectsizesMemoryBlock)[names(effectsizesMemoryBlock)=="ext"] <- c(paste0("ext", blockstring[block]))
  
  assign(paste0("effectsizesMemoryBlock",block), effectsizesMemoryBlock)
  rm(effectsizesMemoryBlock)
  
}

effectsizesMemoryBlock <- merge(effectsizesMemoryBlock1, effectsizesMemoryBlock2, by = "row.names")
row.names(effectsizesMemoryBlock) <- effectsizesMemoryBlock$Row.names
effectsizesMemoryBlock$Row.names <- NULL
effectsizesMemoryBlock <- merge(effectsizesMemoryBlock, effectsizesMemoryBlock3, by = "row.names")
row.names(effectsizesMemoryBlock) <- effectsizesMemoryBlock$Row.names
effectsizesMemoryBlock$Row.names <- NULL

effectsizesMemoryBlock[,c("effect_firstBlock", "effect_secondBlock", "effect_thirdBlock", "ext_firstBlock", "ext_secondBlock", "ext_thirdBlock", "int_firstBlock", "int_secondBlock", "int_thirdBlock")]
write.csv(effectsizesMemoryBlock, paste0("effectsizesMemoryPerBlock_", version, ".csv"))
# effectsizesMemoryBlock<- effectsizesMemoryBlock[,c("effect_firstBlock", "effect_secondBlock", "effect_thirdBlock")]
# write.csv(effectsizesMemoryBlock, paste0("effectsizesMemoryPerBlock_", version, ".csv"))

# spaghetti plot
library(gridExtra)
setwd(ratingsDir)
output <- effectsizesMemoryBlock[,c("effect_firstBlock", "effect_secondBlock", "effect_thirdBlock", "ext_firstBlock", "ext_secondBlock", "ext_thirdBlock", "int_firstBlock", "int_secondBlock", "int_thirdBlock")]

output <- as.data.frame(t(output))
output$rowname <-   row.names(output)
output$block <- rep(c("first", "second", "third"), 3)
output$whatIsMeasured <- rep(c("effectsize", "meanExpGroup", "meanContGroup"), each = 3)
output$whatIsMeasured2 <- rep(c("effectsize", "group mean", "group mean"), each = 3)


for (DV in 1:length(dependentVariablesWide)){
    graph <- ggplot(data=output, aes(x=block, y=get(dependentVariablesWide[DV] ), colour = whatIsMeasured, group = whatIsMeasured)) +
      geom_point() + geom_line() + theme_bw() +
      labs(x="Task block", y="", colour = "", title = paste(dependentVariablesWide[DV])) +
      facet_grid(rows = vars(whatIsMeasured2), scales = "free_y")  + coord_cartesian(ylim = c(-1, 8))
    print(graph)
    ggsave(paste0("Graph_changeOfEffectsizeOverBlocks_", dependentVariablesWide[DV], ".jpeg"))
}

outputEffect <- subset(output, output$whatIsMeasured2 == "effectsize")
outputGroupMean <- subset(output, output$whatIsMeasured2 == "group mean")

for (DV in 1:length(dependentVariablesWide)){
  # create one graph for the effect size (ranging from -1 to 1)
  graphEffect <- ggplot(data=outputEffect, aes(x=block, y=get(dependentVariablesWide[DV] ), group = whatIsMeasured2)) +
    geom_point(color='black') + geom_line(color='black') + theme_bw() +
    theme(legend.position="bottom") +
    labs(x="Task block", y="Cohen's d", group = "", title = paste("Effect size",dependentVariablesWide[DV], version)) + 
    coord_cartesian(ylim = c(-1, 1)) +
    scale_y_continuous(labels = scales::scientific)
  
  # then create a graph for the group average
  graphGroupMean <- ggplot(data=outputGroupMean, aes(x=block, y=get(dependentVariablesWide[DV] ), colour = whatIsMeasured, group = whatIsMeasured)) +
    geom_point() + geom_line() + theme_bw() +
    labs(x="Task block", y="Mean value", colour = "", title = paste("Group means",dependentVariablesWide[DV], version)) + 
    coord_cartesian(ylim = c(floor(min(outputGroupMean[,dependentVariablesWide[DV]]))-1, floor(min(outputGroupMean[,dependentVariablesWide[DV]]))+3)) +
    theme(legend.position="bottom") +
    scale_y_continuous(labels = scales::scientific)

  # then combine those two in a new graphic using gridExtra
  gridExtra::grid.arrange(graphEffect, graphGroupMean, nrow = 2) 
  graph <- gridExtra::arrangeGrob(graphEffect, graphGroupMean, nrow=2) #generates graph
  ggsave(file = paste0("Graph_changeOfEffectsizeOverBlocks_", dependentVariablesWide[DV], "_v2.jpeg"), graph)
  
  rm(graphEffect, graphGroupMean, graph)
}




dfLong$recognisedAndRecalledStrict <- ifelse(dfLong$recognition == dfLong$cuedRecallStrict, 1, 0)
dfLong$recognisedHighConfAndRecalledStrict <- ifelse(dfLong$recognitionConfLevel_4_5_6 == dfLong$cuedRecallStrict, 1, 0)
dfLong$recognisedAndRecalledLenient <- ifelse(dfLong$recognition == dfLong$cuedRecallLenient, 1, 0)
dfLong$recognisedHighConfAndRecalledLenient <- ifelse(dfLong$recognitionConfLevel_4_5_6 == dfLong$cuedRecallLenient, 1, 0)
mean(dfLong$recognisedAndRecalledStrict, na.rm=T)
mean(dfLong$recognisedAndRecalledLenient, na.rm=T)
mean(dfLong$recognisedHighConfAndRecalledStrict, na.rm=T)
mean(dfLong$recognisedHighConfAndRecalledLenient, na.rm=T)


dfLong$memoryStrict <- ifelse(dfLong$recognitionConfLevel_4_5_6 == 1 & dfLong$cuedRecallStrict == 1, "recognisedHighConfAndRecalledStrict", 
                              #ifelse(dfLong$recognition == 1 & dfLong$cuedRecallStrict == 1, "recognisedAndRecalledStrict",
                                                    #ifelse(dfLong$recognition == 0  & dfLong$cuedRecallStrict == 1, "notRecognisedButRecalledStrict",
                                                                  ifelse(dfLong$recognitionConfLevel_4_5_6 == 0 & dfLong$cuedRecallStrict == 1, "notRecognisedHighButAndRecalledStrict",
                                                                                #ifelse(dfLong$recognition == 1 & dfLong$cuedRecallStrict == 0, "recognisedButNotRecalledStrict",
                                                                                              ifelse(dfLong$recognitionConfLevel_4_5_6 == 1 & dfLong$cuedRecallStrict == 0, "recognisedHighConfButNotRecalledStrict", 
                                                                                                            #ifelse(dfLong$recognition == 0 & dfLong$cuedRecallStrict == 0, "notRecognisedAndNotRecalledStrict",
                                                                                                                          ifelse(dfLong$recognitionConfLevel_4_5_6 == 0 & dfLong$cuedRecallStrict == 0, "notRecognisedHighConfAndNotRecalledStrict", 
                                                                                                                                 NA))))#))))

dfLong$memoryLenient <- ifelse(dfLong$recognitionConfLevel_4_5_6 == 1 & dfLong$cuedRecallLenient == 1, "recognisedHighConfAndRecalledLenient",
                               #ifelse(dfLong$recognition == 1 & dfLong$cuedRecallLenient == 1, "recognisedAndRecalledLenient",
                                                   #ifelse(dfLong$recognition == 0 & dfLong$cuedRecallLenient == 1, "notRecognisedButRecalledLenient",
                                                                 ifelse(dfLong$recognitionConfLevel_4_5_6 == 0 & dfLong$cuedRecallLenient == 1, "notRecognisedHighConfButRecalledLenient",
                                                                               #ifelse(dfLong$recognition == 1 & dfLong$cuedRecallLenient == 0, "recognisedButNotRecalledLenient",
                                                                                             ifelse(dfLong$recognitionConfLevel_4_5_6 == 1 & dfLong$cuedRecallLenient == 0, "recognisedHighConfButNotRecalledLenient",
                                                                                                           #ifelse(dfLong$recognition == 0 & dfLong$cuedRecallLenient == 0, "notRecognisedAndNotRecalledLenient",
                                                                                                                         ifelse(dfLong$recognitionConfLevel_4_5_6 == 0 & dfLong$cuedRecallLenient == 0, "notRecognisedHighConfAndNotRecalledLenient",
                                                                                                                                       NA))))#))))
gmodels::CrossTable(dfLong$memoryLenient)
gmodels::CrossTable(dfLong$memoryStrict)

                               
problematic <- subset(dfLong, dfLong$recognisedHighConfAndRecalledLenient == 0)
problematic <- problematic[,c("username", "stimID", "ID", "group", "groupEffectCoded", "cuedRecallLenient", "recognition", "recognitionConfLevel_4_5_6", "recognisedHighConfAndRecalledLenient")]
