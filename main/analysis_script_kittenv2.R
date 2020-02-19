## analysis of the pilot data collected on MTurk may 2018

#### setups ####
#empty work space, load libraries and functions
rm(list=ls())

# define necessary directories
mainDir <- "~/Dropbox/Reading/PhD/Magictricks/behavioural_study"
subDirData <- "data_kittenv2"
dataDir <- file.path(mainDir, subDirData) #"~/Dropbox/Reading/PhD/Magic tricks/behavioural_study/data_fin"
preprocessedDir <- file.path(dataDir, "MagicBehavioural_preprocessed") #"~/Dropbox/Reading/PhD/Magic tricks/behavioural_study/data_fin/MagicBehavioural_preprocessed"
codedDir <- file.path(dataDir, "MagicBehavioural_memory", "preprocessed")

analysisDir <- file.path(dataDir, "Analysis")
ratingsDir <- file.path(analysisDir, "Ratings")
memoryDir <- file.path(ratingsDir, "Memory")


# check whether these directories exist, if not create them
ifelse(!dir.exists(ratingsDir), dir.create(ratingsDir), FALSE) 
ifelse(!dir.exists(memoryDir), dir.create(memoryDir), FALSE) 


# define version 
version <- "kittenv2"
# define block names
blockstring <- c("_firstBlock", "_secondBlock", "_thirdBlock", "")

memoryLevels <- c("cuedRecallStrict", "cuedRecallLenient", 
                  "recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf", 
                  "rememberedStrictAboveAvg", "rememberedLenientAboveAvg", "rememberedStrictHigh", "rememberedLenientHigh")
memoryLabels <- c("cuedRecallStrict", "cuedRecallLenient", 
                  "allConf", "highConf", "aboveAvgConf", 
                  "rememberedStrictAboveAvg", "rememberedLenientAboveAvg", "rememberedStrictHigh", "rememberedLenientHigh")

#helper functions and packages #

source("~/Dropbox/Reading/Codes and functions/R/errorbars.R")
source("~/Dropbox/Reading/Codes and functions/R/rbindcolumns.R")

library(xlsx)
library(psych)
library(ggplot2)

###### read in data sets ###### 
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
output$mot <- rep(c("cont","exp"), each = 1)

outg <- ggplot(output, aes(mot, mean, fill = mot))
outg + geom_bar(stat="identity", position="dodge") + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9))  + scale_x_discrete(limits=c("exp","cont")) + labs(x="Experimental condition", y="Age", fill = "Experimental Condition", title = paste("demogs I", version)) + theme_classic() + scale_fill_discrete(guide=FALSE)

# gender
rm(output)
output <- plyr::count(dfWide, vars = c("gender","group"))

# output <- count(df, c('gender','cond'))
outg <- ggplot(output, aes(group, freq, fill = gender))
outg + geom_bar(stat="identity", position="fill")+ scale_x_discrete(limits=c("cont","exp")) + labs(x="Experimental condition", y="Frequency", fill = "Gender", title = paste("demogs II", version)) + theme_classic() 

#### plots post questionnaire ####
setwd(ratingsDir)

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

########## Compute t-tests and effect sizes for between-group differences in questionnaire scores ########## 

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
write.csv(effectsizesScales, paste0("effectsizesScalesNeuro_", version, ".csv"))

###### lme #####

workspace <- list.files(path = file.path(codedDir), pattern = "_CP.csv") # check whether the data is coded yet or not

if(length(workspace) == 0) { # if data is not coded yet, only look at recognition performance
  dependentVariables <- c("curiosity", "curiosityRT", "decision", "decisionRT",
                          "recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
                          "confidence", "confidenceCorrectTrials")
} else {
  dependentVariables <- c("curiosity", "curiosityRT", "decision", "decisionRT",
                          "cuedRecallStrict", "cuedRecallLenient", 
                          "recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
                          "rememberedStrictAboveAvg", "rememberedLenientAboveAvg", "rememberedStrictHigh", "rememberedLenientHigh",
                          "confidence", "confidenceCorrectTrials")
}


#### get descriptives #### 
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

####### lmer model using curiosity as a continous variable ######

if(length(workspace) == 0) { # if data is not coded yet, only look at recognition performance
  dependentVariables <- c("recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
                          "confidence", "confidenceCorrectTrials")
} else {
  dependentVariables <- c("cuedRecallStrict", "cuedRecallLenient", 
                          "recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
                          "rememberedStrictAboveAvg", "rememberedLenientAboveAvg", "rememberedStrictHigh", "rememberedLenientHigh",
                          "confidence", "confidenceCorrectTrials")
}

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


########## Create barplots to visualise the effects of curiosity and reward on memory performance ########## 

# create a dichomotised curiosity variable using mean-cenetred curiosity
dfLong$curiosity_dich <- ifelse(dfLong$curiosity_dich == -1, "below", 
                                ifelse(dfLong$curiosity_dich == 1, "above", NA)) # create curiosity_dichotom as factor
# define grouping variables
groupingVariables <- c("mediansplitCuriositySample", "mediansplitCuriosityWithinSubject", "curiosity_dich")
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

#loop over dependent variables to create bar plots
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


########## Create histograms to further investigate the relation between reward, curiosity and memory ########## 

# define variables
DV_hist <- c("cuedRecallStrict", "cuedRecallLenient", 
             "recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
             "rememberedStrictAboveAvg", "rememberedLenientAboveAvg", "rememberedStrictHigh", "rememberedLenientHigh")
# create data frame and recode memory as factors
output <- dfLong[,c("ID", "group", "curiosityGroupMeanCentered", DV_hist)]
output$group <- ifelse(output$group == "cont", "No reward", "Reward")

# remove IDs without memory data
output <- subset(output, output$ID != "exp36")
output <- subset(output, output$ID != "exp24")

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



########## 6. Compute t-tests and effect sizes for between-group differences in memory scores ########## 


# define dependent variables (i.e. sum scores)

DV_wide <- c(paste0(memoryLabels, "_abs")) # absolute sum scores
DV_wide <- c(DV_wide, paste0("curBeta_",memoryLabels)) # betas
DV_wide <- c(DV_wide, paste0("curCor_",memoryLabels)) # correlation
for (mem in 1:length(memoryLevels)) { # benefits
  DV_wide <- c(DV_wide, paste0("curBen_cont_", memoryLabels[mem], collapse = ", "))
  DV_wide <- c(DV_wide, paste0("curBen_dich_", memoryLabels[mem], collapse = ", "))
  DV_wide <- c(DV_wide, paste0("curBen_rel_", memoryLabels[mem], collapse = ", "))
}


# for all dependent variables compute two-sample t-test and calculate effect size
for(DV in 1:length(DV_wide)) {
  
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
effectsizes <- rbind.all.columns(effectsizesScales, effectsizesMemory)
setwd(ratingsDir)
write.csv(effectsizes, paste0("effectsizesAll_", version, ".csv"))


########## 7. Create violin plots for sum scores of memory measures in each group ########## 
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


########## 8. Look at the change in memory performance between blocks over time ########## 

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


