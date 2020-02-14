## analysis of the pilot data collected on MTurk may 2018

#### setups ####
#empty work space, load libraries and functions
rm(list=ls())

# define necessary directories
mainDir <- "~/Dropbox/Reading/PhD/Magictricks/behavioural_study"
subDirData <- "data_kittenv2"
version <- "kittenv2"
dataDir <- file.path(mainDir, subDirData) #"~/Dropbox/Reading/PhD/Magic tricks/behavioural_study/data_fin"
groupDir <- file.path(dataDir, "MagicBehavioural_") #"~/Dropbox/Reading/PhD/Magic tricks/behavioural_study/data_fin/MagicBehavioural_"
contDir <- file.path(dataDir, "MagicBehavioural_cont") #"~/Dropbox/Reading/PhD/Magic tricks/behavioural_study/data_fin/MagicBehavioural_cont"
expDir <- file.path(dataDir, "MagicBehavioural_exp") #"~/Dropbox/Reading/PhD/Magic tricks/behavioural_study/data_fin/MagicBehavioural_exp"
memoryDir <- file.path(dataDir, "MagicBehavioural_memory") #"~/Dropbox/Reading/PhD/Magic tricks/behavioural_study/data_fin/MagicBehavioural_memory"
preprocessedDir <- file.path(dataDir, "MagicBehavioural_preprocessed") #"~/Dropbox/Reading/PhD/Magic tricks/behavioural_study/data_fin/MagicBehavioural_preprocessed"
codedDir <- file.path(dataDir, "coded", "preprocessed") #"~/Dropbox/Reading/PhD/Magic tricks/behavioural_study/data_fin/coded/preprocessing"
analysisDir <- file.path(dataDir, "Analysis")
ratingsDir <- file.path(analysisDir, "Ratings")
tricksDir <- file.path(analysisDir, "Tricks")

# check whether these directories exist, if not create them
ifelse(!dir.exists(ratingsDir), dir.create(ratingsDir), FALSE) 

#helper functions and packages #

source("~/Dropbox/Reading/Codes and functions/R/errorbars.R")
library(xlsx)
library(dplyr)
library(psych)
library(ggplot2)

### read in data sets ###
setwd(preprocessedDir)
dfWide <- read.xlsx(paste0("wide_MagicBehavioural_", version, ".xlsx"), sheetName = "Sheet1")
dfLong <- read.xlsx(paste0("long_MagicBehavioural_", version, ".xlsx"), sheetName = "Sheet1")


#### data set handling wide format ####

# compute scale scores

###intrinsic motivation items
# Post1	It was fun to do the experiment.
# Post2	It was boring to do the experiment. ### (R) note: post2_score is already recoded!
# Post3	It was enjoyable to do the experiment.

###task engagement items
# Post4	I was totally absorbed in the experiment.
# Post5	I lost track of time.
# Post6	I concentrated on the experiment.

###interest items
# Post7	The task was interesting.
# Post8	I liked the experiment.
# Post9	I found working on the task interesting.

###boredom items
# Post10	The experiment bored me.
# Post11	I found the experiment fairly dull.
# Post12	I got bored.

####effort/importance
# Post13	I put a lot of effort into this.
# Post14	I didn't try very hard to do well at this activity. ### (R) note: post14_score is already recoded!
# Post15	I tried very hard on this activity.
# Post16	It was important to me to do well at this task.
# Post17	I didn't put much energy into this. ### (R) note: post17_score is already recoded!

###pressure/tension
# Post18	I did not feel nervous at all while doing this. ### (R) note: post18_score is already recoded!
# Post19 I felt very tense while doing this activity.
# Post20	I was very relaxed in doing this experiment. ### (R) note: post20_score is already recoded!
# Post21	I was anxious while working on this task.
# Post22	I felt pressured while doing this task.

###others
# Post23	I tried to find out how many people will be able to find the solution.
# Post24	The amount of magic tricks presented was
# Post25	There were no problems with the internet connection while I participated in the experiment.
# Post26	I was able to see the magic tricks properly.

# Post_comment1	Did you like the experiment? Why? Why not?
# Post_comment2	What did you do while watching the videos?
# Post_comment3	What do you think is the hypothesis behind the experiment?
# Post_comment4	Is there anything else you would like us to know?

scales <- names(dfWide[,c("intrinsicMotivation", "taskEngagement", "interest", "boredom", "effort", "pressure")])
describe(dfWide[,scales])
by(cbind(dfWide[,scales]), dfWide$group, describe)

### demogs ####
# age
output <- describeBy(dfWide[,"age"], group=dfWide$group)
output <- as.data.frame(rbind(output$cont, output$exp))
output$mot <- rep(c("int","ext"), each = 1)

outg <- ggplot(output, aes(mot, mean, fill = mot))
outg + geom_bar(stat="identity", position="dodge") + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9))  + scale_x_discrete(limits=c("ext","int")) + labs(x="Experimental condition", y="Age", fill = "Experimental Condition", title = paste("demogs I", version)) + theme_classic() + scale_fill_discrete(guide=FALSE)

# gender
rm(output)
output <- plyr::count(dfWide, vars = c("gender","group"))

# output <- count(df, c('gender','cond'))
outg <- ggplot(output, aes(group, freq, fill = gender))
outg + geom_bar(stat="identity", position="fill")+ scale_x_discrete(limits=c("cont","exp")) + labs(x="Experimental condition", y="Frequency", fill = "Gender", title = paste("demogs II", version)) + theme_classic() 

#### plots post questionnaire ####
setwd(ratingsDir)
library(ggplot2)
output <- by(cbind(dfWide[,scales]), dfWide$group, describe)
rating <- as.data.frame(rbind(output$cont, output$exp))
rating$mot <- rep(c("intrinsic","extrinsic"), each = 6)

outg <- ggplot(rating, aes(vars, mean, fill = mot))
outg + geom_bar(stat="identity", position="dodge") + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9))  + 
  labs(x="Post task questionnaire", y="Rating", fill = "Experimental Condition", title = paste("all post questions",version)) + 
  theme_classic() + coord_cartesian(ylim = c(0, 7)) +
  scale_x_discrete(limits=paste(scales)) +
  theme(legend.position="bottom")
ggsave("postMainByGroup.jpeg")


output <- by(cbind(dfWide[, c("compliance", "tooManyVids", "problemsInternet", "ableToSee")]),dfWide$group, describe)
output <- as.data.frame(rbind(output$cont, output$exp))
output$mot <- rep(c("intrinsic","extrinsic"), each = 4)
output$vars <- as.factor(output$vars)
levels(output$vars) <- c("task compliance","too many magic tricks","problems internet", "video display")

outg <- ggplot(output, aes(vars, mean, fill = mot))
outg + geom_bar(stat="identity", position="dodge") + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9))  + 
  labs(x="Dependent variable", y="Rating", fill = "Group", title = paste("ratings about experiment", version)) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20, face="bold"), title=element_text(size =20, face="bold"), legend.title = element_text(size=20), legend.text = element_text(size = 20)) + 
  theme_classic() + coord_cartesian(ylim = c(0, 7)) +
  theme(legend.position="bottom") +
  scale_x_discrete(limits=c("task compliance","too many magic tricks","problems internet", "video display")) +
ggsave("postMainByGroup2.jpeg")

rm(outg, output, rating)

#### anovas post questionnaire #####
source("~/Dropbox/Reading/Codes and functions/R/rbindcolumns.R")

# nothing significant

model <- lapply(scales, function(x) {
  lm(substitute(i~groupEffectCoded, list(i = as.name(x))), data = dfWide)})
lapply(model, summary)
test <-t.test(dfWide$pressure~dfWide$group)
test$estimate[1] - test$estimate[2]  
cohen.d(dfWide[,c("group", "pressure")], "group")
t.test(dfWide$taskEngagement~dfWide$group)

# effect sizes
for(scale in 1:length(scales)) {
  
  means <- tapply(dfWide[,scales[scale]], dfWide$group, mean, na.rm = T)
  means <- as.data.frame(t(means))
  # print(means)
  
  data <- dfWide[,c("group", scales[scale])]
  
  d <- cohen.d(data, "group")
  cohen <- as.data.frame(d$cohen.d)
  # print(cohen)
  
  cohen <- merge(cohen, means)
  row.names(cohen) <- paste(scales[scale])
  
  if (scale == 1) {
    effectsizesScales <- cohen
  } else {
    temp_effectsizesScales <- cohen
    effectsizesScales <- rbind.all.columns(effectsizesScales, temp_effectsizesScales) #rbind all columns will induce NA if there was initially no data saved in the loop per participant
    rm(temp_effectsizesScales)
  }
  
}
write.csv(effectsizesScales, paste0("effectsizesScales_", version, ".csv"))


###### lme #####
source("~/Dropbox/Reading/Codes and functions/R/rbindcolumns.R")

if (exists("dfLong") == F) {
  setwd(preprocessedDir)
  dfLong <- xlsx::read.xlsx("long_MagicBehavioural.xlsx", sheetName = "Sheet1")#,  showWarnings = FALSE)
}
groupingVariables <- c("mediansplitCuriosityAya", "mediansplitCuriositySample", "mediansplitCuriosityWithinSubject")
# mediansplitCuriosityAya = median split based an Aya's data, used for initial selection/definition
# mediansplitCuriositySample = median split based on sample data ACROSS subjects
# mediansplitCuriosityWithinSubject =  median split based on sample data within each subject

#dfLong$mediansplitCuriosityAya <- ifelse(dfLong$mediansplitCuriosityAya == "above", 1, -1)


workspace <- list.files(path = file.path(codedDir), pattern = "_CP.csv") # check whether the data is coded yet or not
if(length(workspace) == 0) { # if data is not coded yet, only look at recognition performance
  dependentVariables <- c("curiosity", "curiosityRT", "decision", "decisionRT",
                          "recognition", "recognitionConfLevel_1", "recognitionConfLevel_2", "recognitionConfLevel_3", "recognitionConfLevel_4", "recognitionConfLevel_5", "recognitionConfLevel_6",
                          "recognitionConfLevel_1_2", "recognitionConfLevel_3_4", "recognitionConfLevel_5_6", "recognitionConfLevel_1_2_3", "recognitionConfLevel_4_5_6",
                          "confidence", "confidenceGroupMeanCentered", "confidenceCorrectTrials")
  } else {
    dependentVariables <- c("curiosity", "curiosityRT", "decision", "decisionRT",
                            "cuedRecallLenient", "cuedRecallStrict", 
                            "recognition", "recognitionConfLevel_1", "recognitionConfLevel_2", "recognitionConfLevel_3", "recognitionConfLevel_4", "recognitionConfLevel_5", "recognitionConfLevel_6",
                            "recognitionConfLevel_1_2", "recognitionConfLevel_3_4", "recognitionConfLevel_5_6", "recognitionConfLevel_1_2_3", "recognitionConfLevel_4_5_6",
                            "confidence", "confidenceGroupMeanCentered", "confidenceCorrectTrials")
}


#### get descriptives
for(DV in 1:length(dependentVariables)) {
  # print(dependentVariables[DV])
  descriptive <- describe(dfLong[,dependentVariables[DV] ])
  row.names(descriptive) <- c(paste0(dependentVariables[DV], "_all"))
  
  descriptive_groupwise <-   by(cbind(dfLong[,dependentVariables[DV]]), dfLong$group, describe)
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

# lmer model using curiosity as a categorial variable 

if(length(workspace) == 0) { # if data is not coded yet, only look at recognition performance
  dependentVariables <- c("recognition","recognitionConfLevel_4", "recognitionConfLevel_5", "recognitionConfLevel_6",
                          "recognitionConfLevel_5_6", "recognitionConfLevel_4_5_6",
                          "confidence", "confidenceGroupMeanCentered", "confidenceCorrectTrials")
} else {
  dependentVariables <- c("cuedRecallLenient", "cuedRecallStrict",
                          "recognition","recognitionConfLevel_4", "recognitionConfLevel_5", "recognitionConfLevel_6",
                          "recognitionConfLevel_5_6", "recognitionConfLevel_4_5_6",
                          "confidence", "confidenceGroupMeanCentered", "confidenceCorrectTrials")
}

library(lmerTest)

for (DV in 1:length(dependentVariables)){
  
  # only looking at group, ignoring curiosity
  if (max(dfLong[, dependentVariables[DV]], na.rm = T) > 1){
    LMEmodel_byGroupOnly <- lmer(dfLong[, dependentVariables[DV]] ~ 1 + group + (1 | ID), data = dfLong)
    LMEmodel_byGroupOnly <- lmer(dfLong[, dependentVariables[DV]] ~ 1 + group + (1+curiosityGroupMeanCentered|ID) + (1|stimID), data = dfLong)
  }else{
    LMEmodel_byGroupOnly <- glmer(dfLong[, dependentVariables[DV]] ~ 1 + group + (1 | ID), data = dfLong, family = binomial(link = 'logit'))
    LMEmodel_byGroupOnly <- glmer(dfLong[, dependentVariables[DV]] ~ 1 + group + (1+curiosityGroupMeanCentered|ID) + (1|stimID), data = dfLong, family = binomial(link = 'logit'))
  }
  testGroupOnly <- anova(LMEmodel_byGroupOnly)
  summaryGroupOnly <- summary(LMEmodel_byGroupOnly)
  summaryGroupOnlyCoefficients <- as.data.frame(summaryGroupOnly$coefficients)
  
  row.names(testGroupOnly) <- c(paste0("LME_",dependentVariables[DV], "_byGroupOnly"))
  row.names(summaryGroupOnlyCoefficients) <- c( paste0("LME_",dependentVariables[DV], "_byGroupOnly_Intercept"), paste0("LME_",dependentVariables[DV], "_byGroupOnly_RewardExp"))
   
  # looking at group and curiosity using mediansplitCuriosityAya
  if (max(dfLong[, dependentVariables[DV]], na.rm = T) > 1){
    LMEmodel_byGroupAndMediansplitCuriosityAya <- lmer(dfLong[, dependentVariables[DV]] ~ 1 + group + mediansplitCuriosityAya + group:mediansplitCuriosityAya + (1 | ID), data = dfLong)
    LMEmodel_byGroupAndMediansplitCuriosityAya <- lmer(dfLong[, dependentVariables[DV]] ~ 1 + group + mediansplitCuriosityAya + group:mediansplitCuriosityAya + (1+curiosityGroupMeanCentered|ID) + (1|stimID), data = dfLong)
  }else{
    LMEmodel_byGroupAndMediansplitCuriosityAya <- glmer(dfLong[, dependentVariables[DV]] ~ 1 + group + mediansplitCuriosityAya + group:mediansplitCuriosityAya + (1 | ID), data = dfLong, family = binomial(link = 'logit'))
    LMEmodel_byGroupAndMediansplitCuriosityAya <- glmer(dfLong[, dependentVariables[DV]] ~ 1 + group + mediansplitCuriosityAya + group:mediansplitCuriosityAya + (1+curiosityGroupMeanCentered|ID) + (1|stimID), data = dfLong, family = binomial(link = 'logit'))
  }
  testMediansplitCuriosityAya <- anova(LMEmodel_byGroupAndMediansplitCuriosityAya)
  summaryMediansplitCuriosityAya <- summary(LMEmodel_byGroupAndMediansplitCuriosityAya)
  summaryMediansplitCuriosityAyaCoefficients <- as.data.frame(summaryMediansplitCuriosityAya$coefficients)
  
  row.names(testMediansplitCuriosityAya) <- c(paste0("LME_",dependentVariables[DV], "_byGroupAndMediansplitCuriosityAya_MainEffectReward"), 
                                              paste0("LME_",dependentVariables[DV], "_byGroupAndMediansplitCuriosityAya_MainEffectCuriosity"), 
                                              paste0("LME_",dependentVariables[DV], "_byGroupAndMediansplitCuriosityAya_InteractionEffect")) 
  row.names(summaryMediansplitCuriosityAyaCoefficients) <- c(paste0("LME_",dependentVariables[DV],"_byGroupAndMediansplitCuriosityAya_Intercept"), 
                                                             paste0("LME_",dependentVariables[DV], "_byGroupAndMediansplitCuriosityAya_RewardExp"), 
                                                             paste0("LME_",dependentVariables[DV], "_byGroupAndMediansplitCuriosityAya_BelowMedian"), 
                                                             paste0("LME_",dependentVariables[DV], "_byGroupAndMediansplitCuriosityAya_InteractionEffect"))
  
  # looking at group and curiosity using mediansplitCuriositySample
  if (max(dfLong[, dependentVariables[DV]], na.rm = T) > 1){
    LMEmodel_byGroupAndMediansplitCuriositySample <- lmer(dfLong[, dependentVariables[DV]] ~ 1 + group + mediansplitCuriositySample + group:mediansplitCuriositySample + (1 | ID), data = dfLong)
    LMEmodel_byGroupAndMediansplitCuriositySample <- lmer(dfLong[, dependentVariables[DV]] ~ 1 + group + mediansplitCuriositySample + group:mediansplitCuriositySample + (1+curiosityGroupMeanCentered|ID) + (1|stimID), data = dfLong)
  }else{
    LMEmodel_byGroupAndMediansplitCuriositySample <- glmer(dfLong[, dependentVariables[DV]] ~ 1 + group + mediansplitCuriositySample + group:mediansplitCuriositySample + (1 | ID), data = dfLong, family = binomial(link = 'logit'))
    LMEmodel_byGroupAndMediansplitCuriositySample <- glmer(dfLong[, dependentVariables[DV]] ~ 1 + group + mediansplitCuriositySample + group:mediansplitCuriositySample + (1+curiosityGroupMeanCentered|ID) + (1|stimID), data = dfLong, family = binomial(link = 'logit'))
  }
  
  testMediansplitCuriositySample <- anova(LMEmodel_byGroupAndMediansplitCuriositySample)
  summaryMediansplitCuriositySample <- summary(LMEmodel_byGroupAndMediansplitCuriositySample)
  summaryMediansplitCuriositySampleCoefficients <- as.data.frame(summaryMediansplitCuriositySample$coefficients)
  
  row.names(testMediansplitCuriositySample) <- c(paste0("LME_",dependentVariables[DV], "_byGroupAndMediansplitCuriositySample_MainEffectReward"), 
                                                 paste0("LME_",dependentVariables[DV], "_byGroupAndMediansplitCuriositySample_MainEffectCuriosity"), 
                                                 paste0("LME_",dependentVariables[DV], "_byGroupAndMediansplitCuriositySample_InteractionEffect"))
  row.names(summaryMediansplitCuriositySampleCoefficients) <- c(paste0("LME_",dependentVariables[DV],"_byGroupAndMediansplitCuriositySample_Intercept"), 
                                                                paste0("LME_",dependentVariables[DV], "_byGroupAndMediansplitCuriositySample_RewardExp"), 
                                                                paste0("LME_",dependentVariables[DV], "_byGroupAndMediansplitCuriositySample_BelowMedian"), 
                                                                paste0("LME_",dependentVariables[DV], "_byGroupAndMediansplitCuriositySample_InteractionEffect"))
  
  # looking at group and curiosity using mediansplitCuriosityWithinSubject
  if (max(dfLong[, dependentVariables[DV]], na.rm = T) > 1){
    LMEmodel_byGroupAndMediansplitCuriosityWithinSubject <- lmer(dfLong[, dependentVariables[DV]] ~ 1 + group + mediansplitCuriosityWithinSubject + group:mediansplitCuriosityWithinSubject + (1 | ID), data = dfLong)
    LMEmodel_byGroupAndMediansplitCuriosityWithinSubject <- lmer(dfLong[, dependentVariables[DV]] ~ 1 + group + mediansplitCuriosityWithinSubject + group:mediansplitCuriosityWithinSubject + (1+curiosityGroupMeanCentered|ID) + (1|stimID), data = dfLong)
  }else{
    LMEmodel_byGroupAndMediansplitCuriosityWithinSubject <- glmer(dfLong[, dependentVariables[DV]] ~ 1 + group + mediansplitCuriosityWithinSubject + group:mediansplitCuriosityWithinSubject + (1 | ID), data = dfLong, family = binomial(link = 'logit'))
    LMEmodel_byGroupAndMediansplitCuriosityWithinSubject <- glmer(dfLong[, dependentVariables[DV]] ~ 1 + group + mediansplitCuriosityWithinSubject + group:mediansplitCuriosityWithinSubject + (1+curiosityGroupMeanCentered|ID) + (1|stimID), data = dfLong, family = binomial(link = 'logit'))
  }
  
  testMediansplitCuriosityWithinSubject <- anova(LMEmodel_byGroupAndMediansplitCuriosityWithinSubject)
  summaryMediansplitCuriosityWithinSubject <- summary(LMEmodel_byGroupAndMediansplitCuriosityWithinSubject)
  summaryMediansplitCuriosityWithinSubjectCoefficients <- as.data.frame(summaryMediansplitCuriosityWithinSubject$coefficients)
  
  row.names(testMediansplitCuriosityWithinSubject) <- c(paste0("LME_",dependentVariables[DV], "_byGroupAndMediansplitCuriosityWithinSubject_MainEffectReward"), 
                                                        paste0("LME_",dependentVariables[DV], "_byGroupAndMediansplitCuriosityWithinSubject_MainEffectCuriosity"), 
                                                        paste0("LME_",dependentVariables[DV], "_byGroupAndMediansplitCuriosityWithinSubject_InteractionEffect"))
  row.names(summaryMediansplitCuriosityWithinSubjectCoefficients) <- c(paste0("LME_",dependentVariables[DV],"_byGroupAndMediansplitCuriosityWithinSubject_Intercept"), 
                                                                       paste0("LME_",dependentVariables[DV], "_byGroupAndMediansplitCuriosityWithinSubject_RewardExp"), 
                                                                       paste0("LME_",dependentVariables[DV], "_byGroupAndMediansplitCuriosityWithinSubject_BelowMedian"), 
                                                                       paste0("LME_",dependentVariables[DV], "_byGroupAndMediansplitCuriosityWithinSubject_Median"), 
                                                                       paste0("LME_",dependentVariables[DV], "_byGroupAndMediansplitCuriosityWithinSubject_ExpBelow"),
                                                                       paste0("LME_",dependentVariables[DV], "_byGroupAndMediansplitCuriosityWithinSubject_ExpMedian"))
  
  test <- rbind.all.columns(testGroupOnly, testMediansplitCuriosityAya)
  test <- rbind.all.columns(test, testMediansplitCuriositySample)
  test <- rbind.all.columns(test, testMediansplitCuriosityWithinSubject)
  
  summary <- rbind.all.columns(summaryGroupOnlyCoefficients, summaryMediansplitCuriosityAyaCoefficients)
  summary <- rbind.all.columns(summary, summaryMediansplitCuriositySampleCoefficients)
  summary <- rbind.all.columns(summary, summaryMediansplitCuriosityWithinSubjectCoefficients)
  
  summaryShort <- summaryMediansplitCuriositySampleCoefficients
  
  if (DV == 1) {
    LMEresults <- test
    LMEresults <- summary
    LMEresultsShort <- summaryShort
  }
  else {
    temp_LMEresults <- test
    temp_LMEresults <- summary
    temp_LMEresultsShort <- summaryShort
    LMEresults <- rbind.all.columns(LMEresults, temp_LMEresults) #rbind all columns will induce NA if there was initially no data saved in the loop per participant
    LMEresultsShort <- rbind.all.columns(LMEresultsShort, temp_LMEresultsShort) #rbind all columns will induce NA if there was initially no data saved in the loop per participant
    rm(temp_LMEresults, temp_LMEresultsShort)
  }
  rm(LMEmodel_byGroupOnly, LMEmodel_byGroupAndMediansplitCuriosityAya, LMEmodel_byGroupAndMediansplitCuriositySample, LMEmodel_byGroupAndMediansplitCuriosityWithinSubject)
  rm(testGroupOnly, testMediansplitCuriosityAya, testMediansplitCuriositySample, testMediansplitCuriosityWithinSubject)
  rm(summaryGroupOnly,summaryGroupOnlyCoefficients, summaryMediansplitCuriosityAya, summaryMediansplitCuriosityAyaCoefficients, summaryMediansplitCuriositySample, summaryMediansplitCuriositySampleCoefficients, summaryMediansplitCuriosityWithinSubject, summaryMediansplitCuriosityWithinSubjectCoefficients)
}
rm(test)
rm(summary)
rm(summaryShort)
LMEresults <- round(LMEresults, digits = 5)
LMEresultsShort <- round(LMEresultsShort, digits = 5)

setwd(ratingsDir)
xlsx::write.xlsx(LMEresults, file="LME_Results_CuriosityRatingsAndMemoryScores_byRewardAndCuriosity_v2.xlsx", sheetName = "Sheet1")
xlsx::write.xlsx(LMEresultsShort, file="LME_Results_CuriosityRatingsAndMemoryScores_byRewardAndCuriosity_MedianSplitSample.xlsx", sheetName = "Sheet1")


# lmer model using curiosity as a continous variable
DV <- 11 # look at confidence for correct trials
DV <- 8 # look at high confidence recognition
DV <- 3 # look at recognition
LMEmodel_curiosityContinuous <- glmer(dfLong[, dependentVariables[DV]] ~ 1 + groupEffectCoded + curiosityGroupMeanCentered + rewardByCuriosity + (1 | ID), data = dfLong, family = binomial(link = 'logit'))
LMEmodel_curiosityContinuous <- glmer(dfLong[, dependentVariables[DV]] ~ 1 + groupEffectCoded + curiosityGroupMeanCentered + rewardByCuriosity + (1+curiosityGroupMeanCentered|ID) + (1|stimID), data = dfLong, family = binomial(link = 'logit'))



LMEmodel_curiosityContinuous_recogHighConf <- glmer(recognitionConfLevel_4_5_6 ~ 1 + groupEffectCoded + curiosityGroupMeanCentered + rewardByCuriosity + (1 | ID), data = dfLong, family = binomial(link = 'logit'))
LMEmodel_curiosityContinuous_recogHighConf <- glmer(recognitionConfLevel_4_5_6 ~ 1 + groupEffectCoded + curiosityGroupMeanCentered + rewardByCuriosity + (1+curiosityGroupMeanCentered|ID) + (1|stimID), data = dfLong, family = binomial(link = 'logit'))
summary(LMEmodel_curiosityContinuous_recogHighConf)
exp(2*LMEmodel_curiosityContinuous_recogHighConf@beta[2]) # transform beta estimate for log-Odds from effect coding interpretation to actual change
exp(LMEmodel_curiosityContinuous_recogHighConf@beta[2]) # transform beta estimate for log-Odds from effect coding interpretation to actual change
sjPlot::tab_model(LMEmodel_curiosityContinuous_recogHighConf)
sjPlot::tab_model(LMEmodel_curiosityContinuous_recogHighConf, transform = NULL)

LMEmodel_curiosityContinuous_recog <- glmer(recognition ~ 1 + groupEffectCoded + curiosityGroupMeanCentered + rewardByCuriosity + (1 | ID), data = dfLong, family = binomial(link = 'logit'))
LMEmodel_curiosityContinuous_recog <- glmer(recognition ~ 1 + groupEffectCoded + curiosityGroupMeanCentered + rewardByCuriosity + (1+curiosityGroupMeanCentered|ID) + (1|stimID), data = dfLong, family = binomial(link = 'logit'))
summary(LMEmodel_curiosityContinuous_recog)
sjPlot::tab_model(LMEmodel_curiosityContinuous_recog)
sjPlot::tab_model(LMEmodel_curiosityContinuous_recog, transform = NULL)

LMEmodel_curiosityContinuous_confidenceCorrect <- lmerTest::lmer(confidenceCorrectTrials ~ 1 + groupEffectCoded + curiosityGroupMeanCentered + rewardByCuriosity + (1 | ID), data = dfLong)
LMEmodel_curiosityContinuous_confidenceCorrect <- lmerTest::lmer(confidenceCorrectTrials ~ 1 + groupEffectCoded + curiosityGroupMeanCentered + rewardByCuriosity + (1+curiosityGroupMeanCentered|ID) + (1|stimID), data = dfLong)
summary(LMEmodel_curiosityContinuous_confidenceCorrect)
sjPlot::tab_model(LMEmodel_curiosityContinuous_confidenceCorrect)
sjPlot::tab_model(LMEmodel_curiosityContinuous_confidenceCorrect, transform = NULL)

LMEmodel_curiosityContinuous_confidence <- lmerTest::lmer(confidence ~ 1 + groupEffectCoded + curiosityGroupMeanCentered + rewardByCuriosity + (1 | ID), data = dfLong)
LMEmodel_curiosityContinuous_confidence <- lmerTest::lmer(confidence ~ 1 + groupEffectCoded + curiosityGroupMeanCentered + rewardByCuriosity + (1+curiosityGroupMeanCentered|ID) + (1|stimID), data = dfLong)
summary(LMEmodel_curiosityContinuous_confidence)
sjPlot::tab_model(LMEmodel_curiosityContinuous_confidence)
sjPlot::tab_model(LMEmodel_curiosityContinuous_confidence, transform = NULL)

LMEmodel_curiosityContinuous_confidence_recogAsCov <- lmerTest::lmer(confidence ~ 1 + recognition + groupEffectCoded + curiosityGroupMeanCentered + rewardByCuriosity + (1 | ID), data = dfLong)
LMEmodel_curiosityContinuous_confidence_recogAsCov <- lmerTest::lmer(confidence ~ 1 + recognition + groupEffectCoded + curiosityGroupMeanCentered + rewardByCuriosity + (1+curiosityGroupMeanCentered|ID) + (1|stimID), data = dfLong)
summary(LMEmodel_curiosityContinuous_confidence_recogAsCov)
sjPlot::tab_model(LMEmodel_curiosityContinuous_confidence_recogAsCov)
sjPlot::tab_model(LMEmodel_curiosityContinuous_confidence_recogAsCov, transform = NULL, show.stat = T, show.std = T)
sjstats::std_beta(LMEmodel_curiosityContinuous_confidence_recogAsCov)



LMEmodel_curiosityContinuous_recall <- glmer(cuedRecallStrict ~ 1 + groupEffectCoded + curiosityGroupMeanCentered + rewardByCuriosity + (1 | ID), data = dfLong, family = binomial(link = 'logit'))
LMEmodel_curiosityContinuous_recall <- glmer(cuedRecallStrict ~ 1 + groupEffectCoded + curiosityGroupMeanCentered + rewardByCuriosity + (1+curiosityGroupMeanCentered|ID) + (1|stimID), data = dfLong, family = binomial(link = 'logit'))
summary(LMEmodel_curiosityContinuous_recall)
sjPlot::tab_model(LMEmodel_curiosityContinuous_recall)


LMEmodel <- glmer(recognition ~ 1 + confidenceGroupMeanCentered + curiosityGroupMeanCentered + confidenceGroupMeanCentered:curiosityGroupMeanCentered + (1 | ID), data = dfLong, family = binomial(link = 'logit'))
LMEmodel <- glmer(recognition ~ 1 + confidenceGroupMeanCentered + curiosityGroupMeanCentered + confidenceGroupMeanCentered:curiosityGroupMeanCentered + (1+curiosityGroupMeanCentered|ID) + (1|stimID), data = dfLong, family = binomial(link = 'logit'))
summary(LMEmodel)
sjPlot::tab_model(LMEmodel)

LMEmodel_highConfRec <- glmer(recognitionConfLevel_4_5_6 ~ 1 + confidenceGroupMeanCentered + curiosityGroupMeanCentered + confidenceGroupMeanCentered:curiosityGroupMeanCentered + (1 | ID), data = dfLong, family = binomial(link = 'logit'))
LMEmodel_highConfRec <- glmer(recognitionConfLevel_4_5_6 ~ 1 + confidenceGroupMeanCentered + curiosityGroupMeanCentered + confidenceGroupMeanCentered:curiosityGroupMeanCentered + (1+curiosityGroupMeanCentered|ID) + (1|stimID), data = dfLong, family = binomial(link = 'logit'))
summary(LMEmodel_highConfRec)
sjPlot::tab_model(LMEmodel_highConfRec)


cor.test(dfLong$curiosityGroupMeanCentered, dfLong$confidenceGroupMeanCentered, use = "pairwise.complete.obs")
plot(dfLong$curiosityGroupMeanCentered, dfLong$confidenceGroupMeanCentered, main = paste(version))
abline(lm(dfLong$curiosityGroupMeanCentered ~ dfLong$confidenceGroupMeanCentered), col="red")


for (DV in 1:length(dependentVariables)){
  
  if (max(dfLong[, dependentVariables[DV]], na.rm = T) > 1){
    LMEmodel_curiosityContinuous <- lmerTest::lmer(dfLong[, dependentVariables[DV]] ~ 1 + groupEffectCoded + curiosityGroupMeanCentered + rewardByCuriosity + (1 | ID), data = dfLong)
    LMEmodel_curiosityContinuous <- lmerTest::lmer(dfLong[, dependentVariables[DV]] ~ 1 + groupEffectCoded + curiosityGroupMeanCentered + rewardByCuriosity + (1+curiosityGroupMeanCentered|ID) + (1|stimID), data = dfLong)
  }else{
    LMEmodel_curiosityContinuous <- glmer(dfLong[, dependentVariables[DV]] ~ 1 + groupEffectCoded + curiosityGroupMeanCentered + rewardByCuriosity + (1 | ID), data = dfLong, family = binomial(link = 'logit'))
    LMEmodel_curiosityContinuous <- glmer(dfLong[, dependentVariables[DV]] ~ 1 + groupEffectCoded + curiosityGroupMeanCentered + rewardByCuriosity + (1+curiosityGroupMeanCentered|ID) + (1|stimID), data = dfLong, family = binomial(link = 'logit'))
  }
  summaryCuriosityContinuous <- summary(LMEmodel_curiosityContinuous)
  summaryCuriosityContinuousCoefficients <- as.data.frame(summaryCuriosityContinuous$coefficients)
  
  row.names(summaryCuriosityContinuousCoefficients) <- c(paste0("LME_",dependentVariables[DV], "_Intercept"), paste0("LME_",dependentVariables[DV], "_groupEffectCoded"), 
                                                         paste0("LME_",dependentVariables[DV], "_curiosityGroupMeanCentered"), 
                                                         paste0("LME_",dependentVariables[DV], "_rewardByCuriosity"))
  
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


#### bar plots ####
library(ggplot2)
source('~/Dropbox/Reading/Codes and functions/R/errorbars.R')
setwd(preprocessedDir)
if (exists("groupingVariables") == F) {
  groupingVariables <- c("mediansplitCuriosityAya", "mediansplitCuriositySample", "mediansplitCuriosityWithinSubject")
}
if (exists("dfLong") == F) {
  setwd(preprocessedDir)
  dfLong <- xlsx::read.xlsx("long_MagicBehavioural.xlsx", sheetName = "Sheet1")#,  showWarnings = FALSE)
}
if (exists("dependentVariables") == F) {
  dependentVariables <- c("cuedRecallLenient", "cuedRecallStrict",
                          "recognition", "recognitionConfLevel_1", "recognitionConfLevel_2", "recognitionConfLevel_3", "recognitionConfLevel_4", "recognitionConfLevel_5", "recognitionConfLevel_6",
                          "recognitionConfLevel_1_2", "recognitionConfLevel_3_4", "recognitionConfLevel_5_6", "recognitionConfLevel_1_2_3", "recognitionConfLevel_4_5_6",
                          "confidence", "confidenceCorrectTrials")
}

cols <- c("above" = "#F8766D", "below" = "#00BFC4", "all tricks" = "grey", "median" = "#C77CFF")

setwd(ratingsDir)
for (g in 1:length(groupingVariables)){
  
  for (DV in dependentVariables){
    # create a data frame containing the data for each group divided regarding curiosity ratings
    outputGroup <- summarySEwithin(dfLong, measurevar=DV, betweenvars="group", withinvars=groupingVariables[g], idvar="ID", na.rm = T)
    levels(outputGroup$group) <- c("No reward","Reward")
    
    # # create a data frame containing the data for each group regardless of curiosity ratings
    # outputAll <- summarySE(dfLong, measurevar=DV, groupvars="group", na.rm=T,
    #                        conf.interval=.95, .drop=TRUE)
    # levels(outputAll$group) <- c("No reward","Reward")
    # outputAll[[paste0(groupingVariables[g])]] <- rep("all tricks", 2)
    
    # combine these two data frames and use it for plotting
    output <- outputGroup
    
    # output <- outputGroup
    
    # create bar graph
    graph <- ggplot(output, aes(x=group, y=get(DV), fill=get(groupingVariables[g]))) +
      geom_bar(stat="identity", position="dodge") + geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=get(DV)-1.96*se, ymax=get(DV)+1.96*se)) +
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
    ggsave(paste0("Graph_2x2_", DV, "_by", groupingVariables[g], ".jpeg"))
  }
}

for (g in 1:length(groupingVariables)){

  for (DV in dependentVariables){
    # create a data frame containing the data for each group divided regarding curiosity ratings
    outputGroup <- summarySEwithin(dfLong, measurevar=DV, betweenvars="group", withinvars=groupingVariables[g], idvar="ID", na.rm = T)
    levels(outputGroup$group) <- c("No reward","Reward")
    
    # create a data frame containing the data for each group regardless of curiosity ratings
    outputAll <- summarySE(dfLong, measurevar=DV, groupvars="group", na.rm=T,
                           conf.interval=.95, .drop=TRUE)
    levels(outputAll$group) <- c("No reward","Reward")
    outputAll[[paste0(groupingVariables[g])]] <- rep("all tricks", 2)
    
    # combine these two data frames and use it for plotting
    output <- rbind.all.columns(outputGroup, outputAll)
    
    # output <- outputGroup
    
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
    print(paste0("Graph_2x2_", DV, "_by", groupingVariables[g], "_v2.jpeg"))
    ggsave(paste0("Graph_2x2_", DV, "_by", groupingVariables[g], "_v2.jpeg"))
  }
}

## determine effect sizes
source("~/Dropbox/Reading/Codes and functions/R/rbindcolumns.R")

dependentVariablesWide <- c("cuedRecallStrict", "cuedRecallLenient",
                        "recognition", 
                        "recognitionConfLevel_6", "recognitionConfLevel_5_6", "recognitionConfLevel_4_5_6")
dependentVariablesWide <- c("recognition", "meanConfidence", "meanConfidenceCorrectTrials",
                            "recognitionConfLevel_6", "recognitionConfLevel_5_6", "recognitionConfLevel_4_5_6")
dependentVariablesWide <- c("cuedRecallStrict", "cuedRecallLenient", "recognition", "recognitionConfLevel_1", "recognitionConfLevel_2", "recognitionConfLevel_3", "recognitionConfLevel_4", "recognitionConfLevel_5", "recognitionConfLevel_6",
                            "recognitionConfLevel_1_2", "recognitionConfLevel_3_4", "recognitionConfLevel_5_6", "recognitionConfLevel_1_2_3", "recognitionConfLevel_4_5_6",
                            "meanConfidence", "meanConfidenceCorrectTrials")

for(DV in 1:length(dependentVariablesWide)) {
  means <- tapply(dfWide[,dependentVariablesWide[DV]], dfWide$group, mean, na.rm = T)
  # deviations <- tapply(dfWide[,dependentVariablesWide[DV]], dfWide$group, sd, na.rm = T)
  means <- as.data.frame(t(means))
  # deviations <- as.data.frame(t(deviations))

  data <- dfWide[,c("group", dependentVariablesWide[DV])]
  
  d <- cohen.d(data, "group")
  cohen <- as.data.frame(d$cohen.d)
  rm(d)
  
  cohen <- merge(cohen, means)
  rm(means)
  row.names(cohen) <- paste(dependentVariablesWide[DV])
  if (DV == 1) {
    effectsizesMemory <- cohen
  } else {
    temp_effectsizesMemory <- cohen
    effectsizesMemory <- rbind.all.columns(effectsizesMemory, temp_effectsizesMemory) #rbind all columns will induce NA if there was initially no data saved in the loop per participant
    rm(temp_effectsizesMemory)
  }
  rm(cohen)
}
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
    means <- as.data.frame(t(means))
    
    data <- dfWide[,c("group", paste0(dependentVariablesWide[DV], blockstring[block]))]
    d <- cohen.d(data, "group")
    cohen <- as.data.frame(d$cohen.d)
    rm(d)
    cohen <- merge(cohen, means)
    rm(means)
    row.names(cohen) <- paste(dependentVariablesWide[DV])
    if (DV == 1) {
      effectsizesMemoryBlock <- cohen
    } else {
      temp_effectsizesMemory <- cohen
      effectsizesMemoryBlock <- rbind.all.columns(effectsizesMemoryBlock, temp_effectsizesMemory) #rbind all columns will induce NA if there was initially no data saved in the loop per participant
      rm(temp_effectsizesMemory)
    }
    rm(cohen)
  }
  names(effectsizesMemoryBlock)[names(effectsizesMemoryBlock)=="lower"] <- c(paste0("lower", blockstring[block]))
  names(effectsizesMemoryBlock)[names(effectsizesMemoryBlock)=="effect"] <- c(paste0("effect", blockstring[block]))
  names(effectsizesMemoryBlock)[names(effectsizesMemoryBlock)=="upper"] <- c(paste0("upper", blockstring[block]))
  names(effectsizesMemoryBlock)[names(effectsizesMemoryBlock)=="cont"] <- c(paste0("cont", blockstring[block]))
  names(effectsizesMemoryBlock)[names(effectsizesMemoryBlock)=="exp"] <- c(paste0("exp", blockstring[block]))
  
  assign(paste0("effectsizesMemoryBlock",block), effectsizesMemoryBlock)
  rm(effectsizesMemoryBlock)
  
}

effectsizesMemoryBlock <- merge(effectsizesMemoryBlock1, effectsizesMemoryBlock2, by = "row.names")
row.names(effectsizesMemoryBlock) <- effectsizesMemoryBlock$Row.names
effectsizesMemoryBlock$Row.names <- NULL
effectsizesMemoryBlock <- merge(effectsizesMemoryBlock, effectsizesMemoryBlock3, by = "row.names")
row.names(effectsizesMemoryBlock) <- effectsizesMemoryBlock$Row.names
effectsizesMemoryBlock$Row.names <- NULL

effectsizesMemoryBlock[,c("effect_firstBlock", "effect_secondBlock", "effect_thirdBlock", "exp_firstBlock", "exp_secondBlock", "exp_thirdBlock", "cont_firstBlock", "cont_secondBlock", "cont_thirdBlock")]
write.csv(effectsizesMemoryBlock, paste0("effectsizesMemoryPerBlock_", version, ".csv"))
# effectsizesMemoryBlock<- effectsizesMemoryBlock[,c("effect_firstBlock", "effect_secondBlock", "effect_thirdBlock")]
# write.csv(effectsizesMemoryBlock, paste0("effectsizesMemoryPerBlock_", version, ".csv"))

# spaghetti plot
library(gridExtra)
setwd(ratingsDir)
output <- effectsizesMemoryBlock[,c("effect_firstBlock", "effect_secondBlock", "effect_thirdBlock", "exp_firstBlock", "exp_secondBlock", "exp_thirdBlock", "cont_firstBlock", "cont_secondBlock", "cont_thirdBlock")]

output <- as.data.frame(t(output))
output$rowname <-   row.names(output)
output$block <- rep(c("first", "second", "third"), 3)
output$whatIsMeasured <- rep(c("effectsize", "meanExpGroup", "meanContGroup"), each = 3)
output$whatIsMeasured2 <- rep(c("effectsize", "group mean", "group mean"), each = 3)

dependentVariablesWideMemory <- c("recognition", "recognitionConfLevel_4", "recognitionConfLevel_5", "recognitionConfLevel_6",
                                  "recognitionConfLevel_5_6", "recognitionConfLevel_4_5_6", 
                                  "meanConfidence", "meanConfidenceCorrectTrials")
for (DV in 1:length(dependentVariablesWideMemory)){
    graph <- ggplot(data=output, aes(x=block, y=get(dependentVariablesWideMemory[DV] ), colour = whatIsMeasured, group = whatIsMeasured)) +
      geom_point() + geom_line() + theme_bw() +
      labs(x="Task block", y="", colour = "", title = paste(dependentVariablesWideMemory[DV])) +
      facet_grid(rows = vars(whatIsMeasured2), scales = "free_y")  + coord_cartesian(ylim = c(0, 5))
    print(graph)
    ggsave(paste0("Graph_changeOfEffectsizeOverBlocks_", dependentVariablesWideMemory[DV], ".jpeg"))
}

outputEffect <- subset(output, output$whatIsMeasured2 == "effectsize")
outputGroupMean <- subset(output, output$whatIsMeasured2 == "group mean")

for (DV in 1:length(dependentVariablesWideMemory)){
  # create one graph for the effect size (ranging from -1 to 1)
  graphEffect <- ggplot(data=outputEffect, aes(x=block, y=get(dependentVariablesWideMemory[DV] ), group = whatIsMeasured2)) +
    geom_point(color='black') + geom_line(color='black') + theme_bw() +
    theme(legend.position="bottom") +
    labs(x="Task block", y="Cohen's d", group = "", title = paste("Effect size",dependentVariablesWideMemory[DV], version)) + 
    coord_cartesian(ylim = c(-1, 1)) +
    scale_y_continuous(labels = scales::scientific)
  
  # then create a graph for the group average
  graphGroupMean <- ggplot(data=outputGroupMean, aes(x=block, y=get(dependentVariablesWideMemory[DV] ), colour = whatIsMeasured, group = whatIsMeasured)) +
    geom_point() + geom_line() + theme_bw() +
    labs(x="Task block", y="Mean value", colour = "", title = paste("Group means",dependentVariablesWideMemory[DV], version)) + 
    coord_cartesian(ylim = c(floor(min(outputGroupMean[,dependentVariablesWideMemory[DV]]))-1, floor(min(outputGroupMean[,dependentVariablesWideMemory[DV]]))+3)) +
    theme(legend.position="bottom") +
    scale_y_continuous(labels = scales::scientific)

  # then combine those two in a new graphic using gridExtra
  grid.arrange(graphEffect, graphGroupMean, nrow = 2) 
  graph <- arrangeGrob(graphEffect, graphGroupMean, nrow=2) #generates graph
  ggsave(file = paste0("Graph_changeOfEffectsizeOverBlocks_", dependentVariablesWideMemory[DV], "_v2.jpeg"), graph)
  
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
