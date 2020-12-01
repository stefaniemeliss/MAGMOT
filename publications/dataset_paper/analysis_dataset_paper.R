### this script contains the code used for the dataset paper 

#### setups ####
#empty work space, load libraries and functions
rm(list=ls())

# define necessary directories
analysisDir <- getwd()

filename_tables <- "tables_dataset_paper.xlsx"

# delete output files
ifelse(file.exists(filename_tables), file.remove(filename_tables), FALSE)

# helper functions and packages 
devtools::source_url("https://github.com/stefaniemeliss/MAGMOT/blob/master/functions/errorbars.R?raw=TRUE")
devtools::source_url("https://github.com/stefaniemeliss/MAGMOT/blob/master/functions/rbindcolumns.R?raw=TRUE")

library(dplyr)

# define version --> double-check whether this is necessary or not
version <- "MAGMOT"
version_official <- "fmri"

memoryLevels <- c("cuedRecallStrict", "cuedRecallLenient", 
                  "recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf", 
                  "rememberedStrictAboveAvg", "rememberedLenientAboveAvg", "rememberedStrictHigh", "rememberedLenientHigh")

memoryLabels <- c("cuedRecallStrict", "cuedRecallLenient", 
                  "allConf", "highConf", "aboveAvgConf", 
                  "rememberedStrictAboveAvg", "rememberedLenientAboveAvg", "rememberedStrictHigh", "rememberedLenientHigh")


### downlaod data sets from OSF ###
osfr::osf_auth() # log into OSF
project <- osfr::osf_retrieve_node("fhqb7")
target_dir <- osfr::osf_ls_files(project, pattern = "data") # looks at all files and directories in the project and defines the match with "data"
sub_dir <- osfr::osf_mkdir(target_dir, path = paste0(version_official)) # add folder in OSF data dir

# download all files
osfr::osf_ls_files(sub_dir, pattern = "MAGMOT") %>%
  osfr::osf_download(conflicts = "overwrite")

target_dir <- osfr::osf_ls_files(project, pattern = "stim") # looks at all files and directories in the project and defines the match with "data"

osfr::osf_ls_files(target_dir, pattern = ".csv") %>%
  osfr::osf_download(conflicts = "overwrite")


### read in data sets ###

# data in wide format
scores <- read.csv(paste0(version, "_scores.csv"), stringsAsFactors = F)
other_information <- read.csv(paste0(version, "_other_information.csv"), stringsAsFactors = F)
demographics <- read.csv(paste0(version, "_demographics", ".csv"))
dfWide <- merge(demographics, scores, by = c("ID", "BIDS"))
dfWide <- merge(dfWide, other_information, by = c("ID", "BIDS"))
# data in long format
dfLong <- read.csv(paste0(version, "_experimental_data.csv"))

# stimuli related files
ozono <- read.csv("~/Dropbox/Reading/PhD/Magictricks/stimuli/Ozono_et_al_2020_Detailed_information_about_MagicCATs.csv", stringsAsFactors = F)

duration <- read.table("~/Dropbox/Reading/PhD/Magictricks/stimuli/duration_magictrickfiles.txt", stringsAsFactors = F, header = T)

marker <- read.csv("~/Dropbox/Reading/PhD/Magictricks/stimuli/marker_magictrickfiles.csv", stringsAsFactors = F)

memory_test <- read.csv("recognition_memory_test_incl_deviations_in_pilot.csv", stringsAsFactors = F)
memory_test <- memory_test[,c("stimID", "Recognition.option.1", "Recognition.option.2", "Recognition.option.3", "Recognition.option.4", "Different.wording.in.pilot")]

########## 1. Look at demographics of the sample ########## 

# aim:  descriptive statistics (Elaborate each by group if you have more than one group)
# - Age (mean, SD, range)
# - Sex assigned at birth (absolute count or relative frequencies)
# - Race & ethnicity
#   - White - British
#   - Other White
#   - Asian or Asian British - Bangladeshi
#   - Asian or Asian British - Indian
#   - Asian or Asian British - Pakistani
#   - Asian or Asian British - Chinese
#   - Other Asian background
#   - Black or Black British - Carribbean
#   - Other Black background
#   - Mixed - White and Asian
#   - Mixed - White and Black African
#   - Mixed - White and Black Carribbean
#   - Other Mixed background
#   - Other Ethnic background
#   - Not known
#   - Information refused
# - Education (highest completed), SES
#   - Primary school
#   - GCSEs or equivalent
#   - A-Levels or equivalent
#   - University undergraduate program
#   - University post-graduate program
#   - Doctoral degree
#   - How many years of education have you received, including primary school? (12 = HS diploma, 15 = Bachelor's degree)


# age
psych::describe(dfWide$age)
age <- psych::describeBy(dfWide[,"age"], group=dfWide$group)
age <- as.data.frame(rbind(age$cont, age$exp))
age$mot <- rep(c("cont","exp"), each = 1)

# sex
plyr::count(dfWide, vars =  "sex")
N_per_group <- plyr::count(dfWide, vars = "group")
sex <- plyr::count(dfWide, vars = c("sex","group"))

# ethnicity
ethnicity <- plyr::count(dfWide, vars = c("ethnicity","group"))
dfWide$ethnicity_BAME <- ifelse(dfWide$ethnicity == "White British" | dfWide$ethnicity == "Other White", "White", "BAME")
BAME <- plyr::count(dfWide, vars = c("ethnicity_BAME","group"))

# education
dfWide$yearsOfEducation_num <- ifelse(dfWide$yearsOfEducation == "17 years", 17,
                                  ifelse(dfWide$yearsOfEducation == "20 (18.5 excluding time as phd student)", 20, as.numeric(as.character(dfWide$yearsOfEducation))))

psych::describe(dfWide$yearsOfEducation_num)
edu <- psych::describeBy(dfWide[,"yearsOfEducation_num"], group=dfWide$group)
edu <- as.data.frame(rbind(edu$cont, edu$exp))
edu$mot <- rep(c("cont","exp"), each = 1)

degree <- plyr::count(dfWide, vars = c("education","group"))

# working memory
psych::describe(dfWide$corsiSpan)
corsi <- psych::describeBy(dfWide[,"corsiSpan"], group=dfWide$group)
corsi <- as.data.frame(rbind(corsi$cont, corsi$exp))
corsi$mot <- rep(c("cont","exp"), each = 1)

dfWide$nback_accurary <- dfWide$nback_accurary*100 # convert to percent
psych::describe(dfWide$nback_accurary)
nback <- psych::describeBy(dfWide[,"nback_accurary"], group=dfWide$group)
nback <- as.data.frame(rbind(nback$cont, nback$exp))
nback$mot <- rep(c("cont","exp"), each = 1)


# CREATE DEMOGRAPHICS TABLE FOR PAPER (i.e. Table 1)
demogs <- data.frame() # 1. col: variable, 2. col = control, 3. col = experimental
i = 0
varCol <- 1
contCol <- 2
expCol <- 3

i = i+1
demogs[i,varCol] <- "Age"
demogs[i,contCol] <- paste0(age$mean[1], " (", round(age$sd[1], 2), ")") # control group
demogs[i,expCol] <- paste0(age$mean[2], " (", round(age$sd[2], 2), ")") # experimental group

# this currently produces error, but this is due to missing data in sex
i = i+1
demogs[i,varCol] <- "Sex assigned at birth (% female)"
demogs[i,contCol] <- paste0(sex$freq[sex$sex == "female" & sex$group == "cont"]/N_per_group$freq[N_per_group$group == "cont"]*100) # control group
demogs[i,expCol] <- paste0(sex$freq[sex$sex == "female" & sex$group == "exp"]/N_per_group$freq[N_per_group$group == "exp"]*100) # experimental group

i = i+1
demogs[i,varCol] <- "Ethnicity (% BAME)"
demogs[i,contCol] <- paste0(BAME$freq[BAME$ethnicity_BAME == "BAME" & BAME$group == "cont"]/N_per_group$freq[N_per_group$group == "cont"]*100) # control group
demogs[i,expCol] <- paste0(BAME$freq[BAME$ethnicity_BAME == "BAME" & BAME$group == "exp"]/N_per_group$freq[N_per_group$group == "exp"]*100) # experimental group

i = i+1
demogs[i,varCol] <- "Years of Education"
demogs[i,contCol] <- paste0(edu$mean[1], " (", round(edu$sd[1], 2), ")") # control group
demogs[i,expCol] <- paste0(edu$mean[2], " (", round(edu$sd[2], 2), ")") # experimental group

i = i+1
demogs[i,varCol] <- "Corsi span"
demogs[i,contCol] <- paste0(corsi$mean[1], " (", round(corsi$sd[1], 2), ")") # control group
demogs[i,expCol] <- paste0(corsi$mean[2], " (", round(corsi$sd[2], 2), ")") # experimental group

i = i+1
demogs[i,varCol] <- "n-back (% accuracy)"
demogs[i,contCol] <- paste0(round(nback$mean[1],2), " (", round(nback$sd[1], 2), ")") # control group
demogs[i,expCol] <- paste0(round(nback$mean[2],2), " (", round(nback$sd[2], 2), ")") # experimental group


names(demogs) <- c("", "Control Group", "Experimental Group")

# save demogs
xlsx::write.xlsx(demogs, file=filename_tables, sheetName = "Table_1", append = T, row.names = F) # note: row.names contain variables

########## 2. Determine average testing lengths etc. ########## 

# duration of task blocks 
psych::describe(dfWide$durInMins_firstBlock)
psych::describe(dfWide$durInMins_secondBlock)
psych::describe(dfWide$durInMins_thirdBlock)

# duration of whole experiment
psych::describe(dfWide$durInMins)

block_durInMins <- c(dfWide$durInMins_firstBlock, dfWide$durInMins_secondBlock, dfWide$durInMins_thirdBlock)
psych::describe(block_durInMins)

block_durInSecs <- c(dfWide$durInSecs_firstBlock, dfWide$durInSecs_secondBlock, dfWide$durInSecs_thirdBlock)
psych::describe(block_durInSecs)

# time between experiment and memory tests
psych::describe(dfWide$daysBetweenExpAndMemory)

# average trial length
psych::describe(dfLong$durationTrial)

# duration of memory assessment
psych::describe(dfWide$durMemory)

### information used in trial structure figure - actual presentation times
psych::describe(dfLong$displayVidDuration) # stimulus presentation (this is the actual magic trick + 6 seconds mock video)
psych::describe(dfLong$displayBlankDuration) # blank between end of magic trick video and start of fixation
psych::describe(dfLong$fixationPostVidDuration) # fixation between magic trick and estimate rating
psych::describe(dfLong$rtAnswer) # reaction time estimate rating
psych::describe(dfLong$rtCuriosity) # reaction time estimate rating
psych::describe(dfLong$fixationPostCuriosityDuration) # fixation between curiosity rating and next magic trick

########## 3. Descriptives of magic tricks ########## 

# determine the IDs of all stimuli used
stim_exp <- as.character(unique(dfLong$stimID))
stim_prac <- c("H4_long", "K23")
stim_all <- c(stim_exp, stim_prac)
stim_all <- as.data.frame(stim_all)
names(stim_all) <- "stimID"
stim_all$stimID <- as.character(stim_all$stimID)

# prepare in information from Ozono et al 2020
ozono$X <- NULL #remove col
names(ozono) <- c("stimID", "Name", "Credit", "Phenomena Category", "Materials", "Length", "Subtitle", "Description")
ozono$stimID <- gsub("Short", "short", ozono$stimID)
ozono$stimID <- gsub("Long", "long", ozono$stimID)

# subset ozono to only include relevant tricks
stim_all <- merge(stim_all, ozono, by.x = "stimID")

# categorise stimuli whether they were used for experiment or practice
stim_all$experiment <- ifelse(stim_all$stimID == stim_prac[1] | stim_all$stimID == stim_prac[2], "practice", "experiment")

# convert length
stim_all$Length_secs <- as.numeric(gsub("00:00:", "", stim_all$Length))

# combine this with the actual length of the video files that were used
duration$stimID <- gsub("^\\d\\d_", "", duration$vidFileName)
duration$stimID <- gsub("_combined_small.mp4", "", duration$stimID)
stim_all <- merge(stim_all, duration, by.x = c("stimID", "experiment"))

# combine with memory test information
stim_all <- merge(stim_all, memory_test, by.x = c("stimID"), all.x = T)


# look at average video length
psych::describe(stim_all$vidFileDuration_withoutMock)
psych::describe(stim_all$vidFileDuration)

# look at average video length of tricks used in experiment only
psych::describe(stim_all$vidFileDuration_withoutMock[stim_all$experiment == "experiment"])
psych::describe(stim_all$vidFileDuration[stim_all$experiment == "experiment"])

vidFileDuration <- stim_all$vidFileDuration[stim_all$experiment == "experiment"]


### information used in trial structure figure - intended presentation times
psych::describe(stim_all$vidFileDuration[stim_all$experiment == "experiment"]) # stimulus presentation (this is the actual magic trick + 6 seconds mock video)
psych::describe(dfLong$displayBlankDuration) # blank between end of magic trick video and start of fixation
psych::describe(dfLong$jitterVideo_trial[dfLong$ID == 1]) # jitter between magic trick and estimate rating (constant across subjects)
psych::describe(dfLong$rtAnswer) # reaction time estimate rating
psych::describe(dfLong$rtCuriosity) # reaction time estimate rating
psych::describe(dfLong$jitterRating_trial[dfLong$ID == 1]) # jitter between curiosity rating and next magic trick (constant across subjects)

# CREATE STIMULI TABLE FOR PAPER (i.e. Online-only Table 1)

# add cue image
stim_all$cue <- paste0(stim_all$stimID, "_cue.png")

# marker information
stim_all <- merge(stim_all, marker, by = "stimID", all.x = T)

# select relevant columns
stim_all$duration <- paste0(stim_all$vidFileDuration, " (", stim_all$vidFileDuration_withoutMock, ")")
stim_info <- stim_all[,c("experiment", "stimID", "Name", "Credit", #"vidFileName", "duration", "cue",
                        "Phenomena Category", "Materials", "Description", 
                        "Recognition.option.1", "Recognition.option.2", "Recognition.option.3", "Recognition.option.4", "Different.wording.in.pilot" )]
# order by file name
stim_info <- stim_info[order(stim_info$experiment, stim_info$stimID),] 

# rename columns
names(stim_info) <- c("Occurance", "Stimulus ID", "Name", "Credit", #"Video file name", "Video duration (without mock video)", "Cue image file name",
                     "Phenomena category", "Materials", "Description", 
                     "Recognition option 1", "Recognition option 2", "Recognition option 3", "Recognition option 4", "Different wording in pilot study")

# write file
xlsx::write.xlsx(stim_info, file=filename_tables, sheetName = "Online-only Table_1", append = T, row.names = F, showNA = F) # note: row.names contain variables



# CREATE MARKER TABLE FOR PAPER (i.e. Online-only Table 2)

# create object containing the timings of each magic trick marker and transform this into long format
timing <- stim_all[,c("stimID",  "vidFileName", "duration", "cueImage", 
                      "momentOfSurprise_1",
                      "momentOfSurprise_2",
                      "momentOfSurprise_3", 
                      "momentOfSurprise_4",
                      "momentOfSurprise_5",
                      "momentOfSurprise_6", 
                      "momentOfSurprise_7",
                      "additionalMarker_momentOfSurprise_1", 
                      "additionalMarker_momentOfSurprise_2", 
                      "mainMomentOfSurprise")]
timing_long <- reshape2::melt(timing, id.vars=c("stimID",  "vidFileName", "duration",  "mainMomentOfSurprise"), value.name = "timing")

# create object containing the description of each magic trick marker and transform this into long format
description <- stim_all[,c("stimID",  "vidFileName", "duration", "cue",
                           "momentOfSurprise_1_DESCRIPTION",
                           "momentOfSurprise_2_DESCRIPTION", 
                           "momentOfSurprise_3_DESCRIPTION",
                           "momentOfSurprise_4_DESCRIPTION",
                           "momentOfSurprise_5_DESCRIPTION",
                           "momentOfSurprise_6_DESCRIPTION",
                           "momentOfSurprise_7_DESCRIPTION",
                           "ADDITIONAL.marker.for.moment.1.Description",
                           "ADDITIONAL.marker.for.moment.2.Description",
                           "mainMomentOfSurprise")]

description_long <- reshape2::melt(description, id.vars=c("stimID",  "vidFileName", "duration", "mainMomentOfSurprise"), value.name = "description")

# recode the names of description_long$variable and merge both data sets
description_long$variable <- ifelse(description_long$variable == "momentOfSurprise_1_DESCRIPTION", "momentOfSurprise_1",
                                    ifelse(description_long$variable == "momentOfSurprise_2_DESCRIPTION", "momentOfSurprise_2",
                                           ifelse(description_long$variable == "momentOfSurprise_3_DESCRIPTION", "momentOfSurprise_3",
                                                  ifelse(description_long$variable == "momentOfSurprise_4_DESCRIPTION", "momentOfSurprise_4",
                                                         ifelse(description_long$variable == "momentOfSurprise_5_DESCRIPTION", "momentOfSurprise_5",
                                                                ifelse(description_long$variable == "momentOfSurprise_6_DESCRIPTION", "momentOfSurprise_6",
                                                                       ifelse(description_long$variable == "momentOfSurprise_7_DESCRIPTION", "momentOfSurprise_7", 
                                                                              ifelse(description_long$variable == "ADDITIONAL.marker.for.moment.1.Description", "additionalMarker_momentOfSurprise_1",
                                                                                     ifelse(description_long$variable == "ADDITIONAL.marker.for.moment.2.Description", "additionalMarker_momentOfSurprise_2",
                                                                                            ifelse(description_long$variable == "cue", "cueImage", 
                                                                                                   NA))))))))))

marker_long <- merge(timing_long, description_long, by =c("stimID",  "vidFileName", "duration", "mainMomentOfSurprise", "variable") )

# recode factor levels of variable and sort data frame
marker_long$var <- factor(marker_long$variable, levels = c("cueImage", "momentOfSurprise_1", "additionalMarker_momentOfSurprise_1", "momentOfSurprise_2", "additionalMarker_momentOfSurprise_2", 
                                                           "momentOfSurprise_3", "momentOfSurprise_4", "momentOfSurprise_5", "momentOfSurprise_6", "momentOfSurprise_7"))
marker_long <- marker_long[order(marker_long$stimID, marker_long$var),] 

# rephrase variables
marker_long$Marker <- ifelse(marker_long$variable == "momentOfSurprise_1", "1. moment of suprise",
                             ifelse(marker_long$variable ==  "momentOfSurprise_2", "2. moment of suprise",
                                    ifelse(marker_long$variable == "momentOfSurprise_3", "3. moment of suprise",
                                           ifelse(marker_long$variable == "momentOfSurprise_4", "4. moment of suprise",
                                                  ifelse(marker_long$variable == "momentOfSurprise_5", "5. moment of suprise",
                                                         ifelse(marker_long$variable == "momentOfSurprise_6", "6. moment of suprise",
                                                                ifelse(marker_long$variable == "momentOfSurprise_7",  "7. moment of suprise",
                                                                       ifelse(marker_long$variable == "additionalMarker_momentOfSurprise_1", "additional marker for 1. moment of suprise",
                                                                              ifelse(marker_long$variable == "additionalMarker_momentOfSurprise_2", "additional marker for 2. moment of suprise",
                                                                                     ifelse(marker_long$variable == "cueImage", "Time stamp of cue image",
                                                                                            NA))))))))))
marker_long$Notes <- ifelse(marker_long$mainMomentOfSurprise == "moment of surprise 1", "Main moment of surprise is the first moment of surprise",
                            ifelse(marker_long$mainMomentOfSurprise == "moment of surprise 2", "Main moment of surprise is the second moment of surprise",
                                   ifelse(marker_long$mainMomentOfSurprise == "moment of surprise 3", "Main moment of surprise is the third moment of surprise",
                                          NA)))

# remove unnecessary rows
marker_long$subset <- ifelse(is.na(marker_long$timing) == F, 1, 0)
marker_long <- subset(marker_long, marker_long$subset == 1)

# remove unnecessary columns
marker_long <- marker_long[,c("stimID", "vidFileName", "duration", "Marker", "var", "timing", "description", "Notes")]

# overwrite unncessary cells
marker_long$stimID <- ifelse(marker_long$var == "cueImage", as.character(marker_long$stimID), NA)
marker_long$vidFileName <- ifelse(marker_long$var == "cueImage", as.character(marker_long$vidFileName), NA)
marker_long$duration <- ifelse(marker_long$var == "cueImage", as.character(marker_long$duration), NA)
marker_long$Notes <- ifelse(marker_long$var == "cueImage", as.character(marker_long$Notes), NA)

#rename columns
names(marker_long) <- c("Stimulus ID", "Video file name", "Video duration (without mock video)", 
                        "Marker", "Variable name in dataset", "Timing", "Description", "Notes")

# write file
xlsx::write.xlsx(marker_long, file=filename_tables, sheetName = "Online-only Table_2", append = T, row.names = F, showNA = F) # note: row.names contain variables


### DURATION OF MAGIC TRICKs
# transform long data into wide data
displayVidDuration <- reshape::cast(dfLong, ID~stimID,value="displayVidDuration")
vidFileDuration  <- reshape::cast(stim_all, ~stimID,value="vidFileDuration")
rownames(vidFileDuration) <- "vidFileDuration"

# compute mean, sd, min, max and range for the display duration of each magic trick
mean_displayVidDuration <- mapply(mean, displayVidDuration[-1])
sd_displayVidDuration <- mapply(sd, displayVidDuration[-1])
min_displayVidDuration <- mapply(min, displayVidDuration[-1])
max_displayVidDuration <- mapply(max, displayVidDuration[-1])
range_displayVidDuration <- max_displayVidDuration - min_displayVidDuration

# put all information into one data frame
descript_displayVidDuration <- rbind(mean_displayVidDuration, sd_displayVidDuration)
descript_displayVidDuration <- rbind(descript_displayVidDuration, min_displayVidDuration)
descript_displayVidDuration <- rbind(descript_displayVidDuration, max_displayVidDuration)
descript_displayVidDuration <- rbind(descript_displayVidDuration, range_displayVidDuration)

# merge this with actual file length
descript_displayVidDuration <- as.data.frame(descript_displayVidDuration)
descript_displayVidDuration <- rbind.match.columns(descript_displayVidDuration, vidFileDuration)
descript_displayVidDuration <- as.data.frame(t(descript_displayVidDuration))

# compute difference between vidFileDurattion and displayVidDuration
descript_displayVidDuration$diff_file_mean <- descript_displayVidDuration$mean_displayVidDuration - descript_displayVidDuration$vidFileDuration
descript_displayVidDuration$diff_file_min <- descript_displayVidDuration$min_displayVidDuration - descript_displayVidDuration$vidFileDuration
descript_displayVidDuration$diff_file_max <- descript_displayVidDuration$max_displayVidDuration - descript_displayVidDuration$vidFileDuration

descript_displayVidDuration$diff_mean_min <- descript_displayVidDuration$mean_displayVidDuration - descript_displayVidDuration$min_displayVidDuration
descript_displayVidDuration$diff_mean_max <- descript_displayVidDuration$max_displayVidDuration - descript_displayVidDuration$mean_displayVidDuration

psych::describe(descript_displayVidDuration$diff_file_mean)
psych::describe(descript_displayVidDuration$diff_file_min)
psych::describe(descript_displayVidDuration$diff_file_max)

psych::describe(descript_displayVidDuration$diff_mean_min)
psych::describe(descript_displayVidDuration$diff_mean_max)

# latencies
psych::describe(dfLong$fixationPostVidDuration - dfLong$jitterVideo_trial) # fixation after video
psych::describe(dfLong$fixationPostCuriosityDuration - dfLong$jitterRating_trial) # fixation after curiosity rating



