#empty work space, load libraries and functions
rm(list=ls())

compareDataCollections <- 0
doISCprep <- 0
debug <- 1

####################################################################################################################################
##################################################  COMPARISON DATA COLLECTIONS   ##################################################
####################################################################################################################################

if (compareDataCollections == 1){
  
  # this is the third data collection investigating the effects of reward and curiosity on memory
  # the following code is used to look at whether the magic tricks got similar curiosity ratings and memory performances across data collections 
  
  
  # define all three data collections
  datasets <- c("MAGMOT", "kittenv2", "fin")
  
  # read in the data in long format for all three data collections
  long_MAGMOT <- xlsx::read.xlsx("~/Dropbox/Reading/PhD/Magictricks/fmri_study/Data/preprocessed/long_MAGMOT.xlsx", sheetName = "Sheet1")
  #names(long_MAGMOT)[names(long_MAGMOT)=="responseCuriosity"] <- "curiosity" # to make things easier, rename responseCuriosity
  long_fin <- xlsx::read.xlsx("~/Dropbox/Reading/PhD/Magictricks/behavioural_study/data_fin/MagicBehavioural_preprocessed/long_MagicBehavioural_fin.xlsx", sheetName = "Sheet1")
  long_kittenv2 <- xlsx::read.xlsx("~/Dropbox/Reading/PhD/Magictricks/behavioural_study/data_kittenv2/MagicBehavioural_preprocessed/long_MagicBehavioural_kittenv2.xlsx", sheetName = "Sheet1")
  
  # dedine the indices of interest
  indicesPerTrick <- c("curiosity", "curiosityGroupMeanCentered",
                       "recognition", "recognitionConfLevel_4_5_6",
                       "confidence", "confidenceCorrectTrials")
  indicesPerTrickMean <- c( "meanCuriosity", "meanCuriosityGroupMeanCentered",
                            "meanRecognition", "meanRecognitionConfLevel_4_5_6", 
                            "meanConfidence", "meanConfidenceCorrectTrials" )
  
  # create a df for each of the indices per trick for each data collection
  for (i in seq_along(indicesPerTrick)){
    assign(paste0(indicesPerTrick[i],"_fin"), reshape::cast(long_fin, ID~stimID,value=paste0(indicesPerTrick[i])))
    assign(paste0(indicesPerTrick[i],"_kittenv2"), reshape::cast(long_kittenv2, ID~stimID,value=paste0(indicesPerTrick[i])))
    assign(paste0(indicesPerTrick[i],"_MAGMOT"), reshape::cast(long_MAGMOT, ID~stimID,value=paste0(indicesPerTrick[i])))
  }
  
  # create list of tricks
  tricks <- as.character(levels(long_MAGMOT$stimID)) 
  
  # calculate mean values for each magic trick for each index
  for (iMean in seq_along(indicesPerTrickMean)){
    for (data in seq_along(datasets)){
      assign(paste0(indicesPerTrickMean[iMean],"_",datasets[data]), numeric(length(tricks))) # creates object meanCuriosity num[1:36]
      currentMean <- numeric(length(tricks))
      
      currentIndex <- indicesPerTrick[iMean] # get current index
      currentIndex <- paste0(indicesPerTrick[iMean], "_", datasets[data]) # get current index
      meansPerTrick <-   colMeans(get(currentIndex), na.rm = T) # calculate mean for each magic trick
      
      # combine all means in one data frame
      if (iMean == 1 && data == 1){
        dfMeans <- data.frame(meansPerTrick)
        names(dfMeans) <- paste0(indicesPerTrickMean[iMean], "_", datasets[data])
        dfMeans$stimID <- row.names(dfMeans)
      } else {
        dfMeans_temp <- data.frame(meansPerTrick)
        names(dfMeans_temp) <- paste0(indicesPerTrickMean[iMean], "_", datasets[data])
        dfMeans_temp$stimID <- row.names(dfMeans_temp)
        dfMeans <- merge(dfMeans, dfMeans_temp, by = "stimID")
        rm(dfMeans_temp)
      }
    }
  }
  
  long_MAGMOT$vidDurCalc <- long_MAGMOT$displayVidOffset - long_MAGMOT$mockOffset
  
  dur_MAGMOT <- reshape::cast(long_MAGMOT, ID~stimID,value="vidDurCalc")
  avgDur_MAGMOT <- colMeans(dur_MAGMOT, na.rm = T) 
  
  dfMeans_temp <- data.frame(avgDur_MAGMOT)
  names(dfMeans_temp) <- "avgVidDur_MAGMOT"
  dfMeans_temp$stimID <- row.names(dfMeans_temp)
  dfMeans <- merge(dfMeans, dfMeans_temp, by = "stimID")
  rm(dfMeans_temp, dur_MAGMOT, avgDur_MAGMOT)
  
  # remove all variables no longer needd
  rm(list=ls(pattern = "_MAGMOT"))
  rm(list=ls(pattern = "_fin"))
  rm(list=ls(pattern = "_kittenv2"))
  
  # compute standardised curiosity ratings, get the rank for each magictrick for each data collection and calculate median splits
  dfMeans$meanCuriosityStandardised_MAGMOT <- (dfMeans$meanCuriosity_MAGMOT - mean(dfMeans$meanCuriosity_MAGMOT)) / sd(dfMeans$meanCuriosity_MAGMOT)
  dfMeans$mediansplitCuriosity_MAGMOT <- ifelse(dfMeans$meanCuriosity_MAGMOT > median(dfMeans$meanCuriosity_MAGMOT), "above", 
                                                ifelse(dfMeans$meanCuriosity_MAGMOT < median(dfMeans$meanCuriosity_MAGMOT), "below", "median")) 
  #dfMeans$mediansplitCuriosityGroupMeanCentered_MAGMOT <- ifelse(dfMeans$meanCuriosityGroupMeanCentered_MAGMOT > median(dfMeans$meanCuriosityGroupMeanCentered_MAGMOT), "above", 
  #                                              ifelse(dfMeans$meanCuriosityGroupMeanCentered_MAGMOT < median(dfMeans$meanCuriosityGroupMeanCentered_MAGMOT), "below", "median")) 
  dfMeans$rankedCuriosity_MAGMOT <- rank(dfMeans$meanCuriosityGroupMeanCentered_MAGMOT)
  
  dfMeans$meanCuriosityStandardised_kittenv2 <- (dfMeans$meanCuriosity_kittenv2 - mean(dfMeans$meanCuriosity_kittenv2)) / sd(dfMeans$meanCuriosity_kittenv2)
  dfMeans$mediansplitCuriosity_kittenv2 <- ifelse(dfMeans$meanCuriosity_kittenv2 > median(dfMeans$meanCuriosity_kittenv2), "above", 
                                                  ifelse(dfMeans$meanCuriosity_kittenv2 < median(dfMeans$meanCuriosity_kittenv2), "below", "median")) 
  dfMeans$rankedCuriosity_kittenv2 <- rank(dfMeans$meanCuriosity_kittenv2)
  
  dfMeans$meanCuriosityStandardised_fin <- (dfMeans$meanCuriosity_fin - mean(dfMeans$meanCuriosity_fin)) / sd(dfMeans$meanCuriosity_fin)
  dfMeans$mediansplitCuriosity_fin <- ifelse(dfMeans$meanCuriosity_fin > median(dfMeans$meanCuriosity_fin), "above", 
                                             ifelse(dfMeans$meanCuriosity_fin < median(dfMeans$meanCuriosity_fin), "below", "median")) 
  dfMeans$rankedCuriosity_fin <- rank(dfMeans$meanCuriosity_fin)
  
  # compare mediansplits across different data collections
  dfMeans$differentSplits_MAGMOT_kittenv2 <- ifelse(dfMeans$mediansplitCuriosity_MAGMOT != dfMeans$mediansplitCuriosity_kittenv2, "different", "same") 
  dfMeans$stimID[dfMeans$differentSplits_MAGMOT_kittenv2 == "different"]
  dfMeans$differentSplits_MAGMOT_fin <- ifelse(dfMeans$mediansplitCuriosity_MAGMOT != dfMeans$mediansplitCuriosity_fin, "different", "same") 
  dfMeans$stimID[dfMeans$differentSplits_MAGMOT_fin == "different"]
  dfMeans$differentSplits_kittenv2_fin <- ifelse(dfMeans$mediansplitCuriosity_fin != dfMeans$mediansplitCuriosity_kittenv2, "different", "same") 
  dfMeans$stimID[dfMeans$differentSplits_kittenv2_fin == "different"]
  
  # save this overview
  setwd("~/Dropbox/Reading/PhD/Magictricks/fmri_study/Analysis/Tricks/")
  xlsx::write.xlsx(dfMeans, file="MAGMOT_recognitionAndCuriosity_perTrick.xlsx", sheetName = "Sheet1", row.names = F)
  
}



####################################################################################################################################
############################################################  SET UPS   ############################################################
####################################################################################################################################


source("~/Dropbox/Reading/Codes and functions/R/rbindcolumns.R")

# define core variables
version = "fmri"
blockstring <- c("_firstBlock", "_secondBlock", "_thirdBlock", "")
SME_outcome <- c("bothRemembered", "bothForgotten", "differentResponses")
confidence_levels <- c("allConf", "aboveAvgConf", "highConf")
TR = 2

# define whether data on VM should be overwritten
overwrite = "no"

# define whether you want feedback printed to the console
feedback = "no"


# define necessary directories
mainDir <- "~/Dropbox/Reading/PhD/Magictricks/fmri_study"
preDir <- file.path(mainDir, "Data", "MAGMOT_pre")
postDir <- file.path(mainDir, "Data", "MAGMOT_post")
brainDir <- file.path(mainDir, "Data", "brain_activation")
dataMemoryDir <- file.path(mainDir, "Data", "magicmemory_fmri", "decrypted")
codedDir <- file.path(mainDir, "Data", "magicmemory_fmri", "coded")
preprocessedDir <- file.path(mainDir, "Data", "preprocessed")
preprocessedQuestDir <- file.path(preprocessedDir, "quest")
preprocessedMemoryDir <- file.path(preprocessedDir, "memory")
preprocessedLongDir <- file.path(preprocessedDir, "long")

preprocessedconcatDir <- file.path(preprocessedDir, "BIDS_eventfiles")
preprocessedEventsRootDir <- file.path(mainDir, "derivatives", "magictrickwatching", "concat")

dirVM <- '/Users/stefaniemeliss/cinn/2018/MAGMOT'
BIDSdirVM <- file.path(dirVM, "MAGMOT_BIDS")
preprocessedEventsRootDirVM <-  file.path(dirVM, "derivatives", "magictrickwatching")


# check whether directories for preprocessed data exist, if so empty them and recreate
ifelse(dir.exists(preprocessedDir), unlink(preprocessedDir, recursive = TRUE), FALSE)
ifelse(dir.exists(preprocessedEventsRootDir), unlink(preprocessedEventsRootDir, recursive = TRUE), FALSE)
ifelse(!dir.exists(preprocessedDir), dir.create(preprocessedDir), FALSE)
ifelse(!dir.exists(preprocessedQuestDir), dir.create(preprocessedQuestDir), FALSE)
ifelse(!dir.exists(preprocessedMemoryDir), dir.create(preprocessedMemoryDir), FALSE)
ifelse(!dir.exists(preprocessedLongDir), dir.create(preprocessedLongDir), FALSE)
ifelse(!dir.exists(preprocessedEventsRootDir), dir.create(preprocessedEventsRootDir), FALSE)
ifelse(!dir.exists(preprocessedconcatDir), dir.create(preprocessedconcatDir), FALSE)

# read in data on curuiosity mean split
if (exists("dfMeans") == F) {
  setwd(preprocessedDir)
  dfMeans <- xlsx::read.xlsx("~/Dropbox/Reading/PhD/Magictricks/fmri_study/Analysis/Tricks/MAGMOT_recognitionAndCuriosity_perTrick.xlsx", sheetName = "Sheet1")#,  showWarnings = FALSE)
}
dfMeans <- dfMeans[,c("stimID", "meanCuriosity_MAGMOT", "meanCuriosityStandardised_MAGMOT", "mediansplitCuriosity_MAGMOT", "avgVidDur_MAGMOT")]

####################################################################################################################################
####################################################  PROCESS PSYTOOLKIT DATA   ####################################################
####################################################################################################################################


#################################################### read in MAGMOT_pre ####################################################
setwd(preDir)
#MAGMOT_pre <- read.table("data.csv", sep = ",", na.strings = "", fill = T, row.names = NULL)
#names(MAGMOT_pre) <- names(MAGMOT_pre)[-1]
#MAGMOT_pre[dim(MAGMOT_pre)[2]] <- NULL
data_pilot <- xlsx::read.xlsx("data_pilot.xlsx", sheetName="data") # data of participant 1 to 6 
data_rest <- xlsx::read.xlsx("data.xlsx", sheetName="data") # data of participant 7 onwards
MAGMOT_pre <- rbind.all.columns(data_pilot, data_rest)
rm(data_pilot, data_rest)

# remove incomplete data sets
MAGMOT_pre <- subset(MAGMOT_pre, !is.na(MAGMOT_pre$TIME_end)) # subset data to only include complete data sets
MAGMOT_pre <- subset(MAGMOT_pre, !is.na(MAGMOT_pre$corsi.1)) # subset data to only include complete data sets
MAGMOT_pre <- subset(MAGMOT_pre, nchar(as.character(MAGMOT_pre$id.1)) <= 2) # subset data to not include test data sets

# change participant IDs
MAGMOT_pre$id.1 <- gsub("O", "", MAGMOT_pre$id.1) #replace any O in the file names
MAGMOT_pre$id.1 <- ifelse(nchar(MAGMOT_pre$id.1)==1, paste0("0", MAGMOT_pre$id.1), MAGMOT_pre$id.1) # add a zeo if not present in front of one digit numbers

# delete doubled/reassigned IDs
participantToDelete1 <- as.character(MAGMOT_pre$participant[(MAGMOT_pre$id.1 == 16 & MAGMOT_pre$DOB.1 != "25 08 1995")]) # ppt born on 25/08/95 should be kept
participantToDelete2 <- as.character(MAGMOT_pre$participant[(MAGMOT_pre$id.1 == 37 & MAGMOT_pre$DOB.1 != "03/01/1996")]) # ppt born on 03/01/96 should be kept
participantToDelete3 <- as.character(MAGMOT_pre$participant[(MAGMOT_pre$id.1 == "09" & MAGMOT_pre$DOB.1 != "19/03/1998")]) # ppt born on 19/03/98 should be kept
participantToDelete4 <- as.character(MAGMOT_pre$participant[(MAGMOT_pre$id.1 == 23 & MAGMOT_pre$DOB.1 == "17/12/1991")]) # ppt born on 17/12/91 should be deleted
participantToDelete5 <- as.character(MAGMOT_pre$participant[(MAGMOT_pre$id.1 == 32 & MAGMOT_pre$DOB.1 == "10/03/1997")]) # ppt born on 10/03/97 should be deleted
participantToDelete6 <- as.character(MAGMOT_pre$participant[(MAGMOT_pre$id.1 == 46 & MAGMOT_pre$DOB.1 == "15/04/1996")]) # ppt born on 15/04/96 should be deleted
participantToDelete7 <- as.character(MAGMOT_pre$participant[(MAGMOT_pre$id.1 == 48 & MAGMOT_pre$DOB.1 != "15/11/1999")]) # ppt born on 15/11/99 should be deleted
participantToDelete8 <- as.character(MAGMOT_pre$participant[(MAGMOT_pre$id.1 == 50 & MAGMOT_pre$DOB.1 != "05/09/1993")]) # ppt born on 05/09/93 should be kept
participantToDelete9 <- as.character(MAGMOT_pre$participant[(MAGMOT_pre$id.1 == 12 & MAGMOT_pre$DOB.1 == "28/02/1998")]) # ppt withdraw
MAGMOT_pre <- subset(MAGMOT_pre, MAGMOT_pre$participant != participantToDelete1 & MAGMOT_pre$participant != participantToDelete2 & MAGMOT_pre$participant != participantToDelete3 & MAGMOT_pre$participant != participantToDelete4 &
                       MAGMOT_pre$participant != participantToDelete5 & MAGMOT_pre$participant != participantToDelete6 & MAGMOT_pre$participant != participantToDelete7 & MAGMOT_pre$participant != participantToDelete8 & MAGMOT_pre$participant != participantToDelete9)
rm(list=ls(pattern = "participantToDelete"))
sort(MAGMOT_pre$id.1)

### recode demographic information and then delete the old variable
MAGMOT_pre$gender <- ifelse(MAGMOT_pre$gender.1 == 1, "male", # gender
                            ifelse(MAGMOT_pre$gender.1 == 2, "female",
                                   ifelse(MAGMOT_pre$gender.1 == 3, "different", NA)))
MAGMOT_pre$gender.1 <- NULL
MAGMOT_pre$ethnicity <- ifelse(MAGMOT_pre$ethnicity.1 == 1, "White British", # ethnicity
                               ifelse(MAGMOT_pre$ethnicity.1 == 2, "Other White",
                                      ifelse(MAGMOT_pre$ethnicity.1 == 3, "Asian or Asian British - Bangladeshi",
                                             ifelse(MAGMOT_pre$ethnicity.1 == 4, "Asian or Asian British - Indian",
                                                    ifelse(MAGMOT_pre$ethnicity.1 == 5, "Asian or Asian British - Pakistani",
                                                           ifelse(MAGMOT_pre$ethnicity.1 == 6, "Asian or Asian British - Chinese",
                                                                  ifelse(MAGMOT_pre$ethnicity.1 == 7, "Other Asian background",
                                                                         ifelse(MAGMOT_pre$ethnicity.1 == 8, "Black or Black British - Carribbean",
                                                                                ifelse(MAGMOT_pre$ethnicity.1 == 9, "Other Black background",
                                                                                       ifelse(MAGMOT_pre$ethnicity.1 == 10, "Mixed - White and Asian",
                                                                                              ifelse(MAGMOT_pre$ethnicity.1 == 11, "Mixed - White and Black African",
                                                                                                     ifelse(MAGMOT_pre$ethnicity.1 == 12, "Mixed - White and Black Carribbean",
                                                                                                            ifelse(MAGMOT_pre$ethnicity.1 == 13, "Other Mixed background",
                                                                                                                   ifelse(MAGMOT_pre$ethnicity.1 == 14, "Other Ethnic background",
                                                                                                                          ifelse(MAGMOT_pre$ethnicity.1 == 15, "Not known",
                                                                                                                                 ifelse(MAGMOT_pre$ethnicity.1 == 16, "Information refused",
                                                                                                                                        NA))))))))))))))))
MAGMOT_pre$ethnicity.1 <- NULL
MAGMOT_pre$education <- ifelse(MAGMOT_pre$education.1 == 1, "Primary school", # education
                               ifelse(MAGMOT_pre$education.1 == 2, "GCSEs or equivalent",
                                      ifelse(MAGMOT_pre$education.1 == 3, "A-Levels or equivalent",
                                             ifelse(MAGMOT_pre$education.1 == 4, "University undergraduate program",
                                                    ifelse(MAGMOT_pre$education.1 == 5, "University post-graduate program",
                                                           ifelse(MAGMOT_pre$education.1 == 6, "Doctoral degree",
                                                                  NA))))))
MAGMOT_pre$education.1 <- NULL
MAGMOT_pre$employment <-  ifelse(MAGMOT_pre$employment.1 == 1, "Unemployed", # employment
                                 ifelse(MAGMOT_pre$employment.1 == 2, "Self-employed part-time",
                                        ifelse(MAGMOT_pre$employment.1 == 3, "Self-employed full-time",
                                               ifelse(MAGMOT_pre$employment.1 == 4, "Part-time employment within organisation/company",
                                                      ifelse(MAGMOT_pre$employment.1 == 5, "Full-time employment within organisation/company",
                                                             ifelse(MAGMOT_pre$employment.1 == 6, "Full-time student",
                                                                    ifelse(MAGMOT_pre$employment.1 == 7, "Part-time student",
                                                                           NA)))))))
MAGMOT_pre$employment.1 <- NULL
MAGMOT_pre$handedness <- ifelse(MAGMOT_pre$handedness.1 == 1, "Left", # handedness
                                ifelse(MAGMOT_pre$handedness.1 == 2, "Right",
                                       ifelse(MAGMOT_pre$handedness.1 == 3, "Both",
                                              NA)))
MAGMOT_pre$handedness.1 <- NULL
MAGMOT_pre$vision <-  ifelse(MAGMOT_pre$vision.1 == 1, "normal", # vision
                             ifelse(MAGMOT_pre$vision.1 == 2 & MAGMOT_pre$vision_contactlenses.1 == 1, "corrected with contacts",
                                    ifelse(MAGMOT_pre$vision.1 == 2 & MAGMOT_pre$vision_contactlenses.1 == 2, "corrected with glasses",
                                           NA)))
MAGMOT_pre$vision.1 <- NULL
MAGMOT_pre$vision_contactlenses.1 <- NULL
MAGMOT_pre$english <- ifelse(MAGMOT_pre$english_native.1 == 1, "native speaker", # english proficiency
                             ifelse(MAGMOT_pre$english_native.1 == 2, "non-native speaker",
                                    NA))
MAGMOT_pre$english_native.1 <- NULL

### rename variables in MAGMOT_pre
names(MAGMOT_pre)[names(MAGMOT_pre)=="id.1"] <- "ID"
names(MAGMOT_pre)[names(MAGMOT_pre)=="age.1"] <- "age"
names(MAGMOT_pre)[names(MAGMOT_pre)=="education_years.1"] <- "yearsOfEducation"
names(MAGMOT_pre)[names(MAGMOT_pre)=="student_subject.1"] <- "studySubject"
names(MAGMOT_pre)[names(MAGMOT_pre)=="english_age.1"] <- "AgeEnglishAcquisition"
names(MAGMOT_pre)[names(MAGMOT_pre)=="health.1"] <- "health"
names(MAGMOT_pre)[names(MAGMOT_pre)=="neurodisorders.1"] <- "neurodisorders"
names(MAGMOT_pre)[names(MAGMOT_pre)=="participant"] <- "preFile"
names(MAGMOT_pre)[names(MAGMOT_pre)=="TIME_start"] <- "startPre"
names(MAGMOT_pre)[names(MAGMOT_pre)=="TIME_end"] <- "endPre"
names(MAGMOT_pre)[names(MAGMOT_pre)=="DOB.1"] <- "DOB"


#### MRI screening: any yes?
items <- c("screening_MRI.1", "screening_MRI.2", "screening_MRI.3", "screening_MRI.4", "screening_MRI.5", "screening_MRI.6", "screening_MRI.7", "screening_MRI.8", "screening_MRI.9", "screening_MRI.10",
           "screening_MRI.11", "screening_MRI.12", "screening_MRI.13", "screening_MRI.14", "screening_MRI.15", "screening_MRI.16", "screening_MRI.17", "screening_MRI.18", "screening_MRI.19")
# compute summary score for MRI screening answers
for (item in seq_along(items)){
  if (item == 1){
    MAGMOT_pre$screening_MRI <- MAGMOT_pre[[paste0(items[item])]]
  } else {
    MAGMOT_pre$screening_MRI <- MAGMOT_pre$screening_MRI + MAGMOT_pre[[paste0(items[item])]]
  }
  MAGMOT_pre[[paste0(items[item])]] <- ifelse(MAGMOT_pre[[paste0(items[item])]] == 0, "no",
                                              ifelse(MAGMOT_pre[[paste0(items[item])]] == 1, "yes",
                                                     NA))
}
rm(item, items)

#### BIS-BAS
# 1 - very true for me; 2 - somewhat true for me; 3 - somewhat false for me; 4 - very false for me
items <- c("BISBAS.2", "BISBAS.20") # two items need recoding
for (item in items){ # recode items
  item_score <-  paste0(item, "R")
  MAGMOT_pre[,item_score] <- ifelse(MAGMOT_pre[,item] == 1, 4,
                                    ifelse(MAGMOT_pre[,item] == 2, 3,
                                           ifelse(MAGMOT_pre[,item] == 3, 2,
                                                  ifelse(MAGMOT_pre[,item] == 4, 1,
                                                         NA))))
  rm(item_score)
}
rm(item, items)

# compute scales
MAGMOT_pre$BIS <- MAGMOT_pre$BISBAS.15 +  MAGMOT_pre$BISBAS.1 +  MAGMOT_pre$BISBAS.12 +  MAGMOT_pre$BISBAS.2R + MAGMOT_pre$BISBAS.17 + MAGMOT_pre$BISBAS.20R #BIS items: 15, 1, 8, 12, 2 (R), 17, 20 (R)
MAGMOT_pre$BAS_rewardresponsiveness <-  MAGMOT_pre$BISBAS.7 +  MAGMOT_pre$BISBAS.4 +  MAGMOT_pre$BISBAS.16 +  MAGMOT_pre$BISBAS.6 + MAGMOT_pre$BISBAS.13 #BAS reward responsiveness: 7, 4, 16, 6, 13
MAGMOT_pre$BAS_drive <- MAGMOT_pre$BISBAS.9 +  MAGMOT_pre$BISBAS.3 +  MAGMOT_pre$BISBAS.11 +  MAGMOT_pre$BISBAS.19 #BAS drive: 9, 3, 11, 19
MAGMOT_pre$BAS_funseeking <- MAGMOT_pre$BISBAS.10 +  MAGMOT_pre$BISBAS.18 +  MAGMOT_pre$BISBAS.5 +  MAGMOT_pre$BISBAS.14 #BAS fun seeking: 10, 18, 5, 14

### Need for Cognition
# +4 = very strong agreement; +3 = strong agreement; +2 = moderate agreement; +1 = slight agreement; 0 = neither agreement nor disagreement;
# -1 = slight disagreement; -2 = moderate disagreement; -3 = strong disagreement; -4 = very strong disagreement
items <- c("NeedForCognition.3", "NeedForCognition.4", "NeedForCognition.5", "NeedForCognition.8", "NeedForCognition.9", "NeedForCognition.12", "NeedForCognition.16", "NeedForCognition.17") # items that need recoding
for (item in items){ # recode items
  item_score <-  paste0(item, "R")
  MAGMOT_pre[,item_score] <- ifelse(MAGMOT_pre[,item] == 4, -4,
                                    ifelse(MAGMOT_pre[,item] == 3, -3,
                                           ifelse(MAGMOT_pre[,item] == 2, -2,
                                                  ifelse(MAGMOT_pre[,item] == 1, -1,
                                                         ifelse(MAGMOT_pre[,item] == 0, 0,
                                                                ifelse(MAGMOT_pre[,item] == -1, 1,
                                                                       ifelse(MAGMOT_pre[,item] == -2, 2,
                                                                              ifelse(MAGMOT_pre[,item] == -3, 3,
                                                                                     ifelse(MAGMOT_pre[,item] == -4, 4,
                                                                                            NA)))))))))
  rm(item_score)
}
rm(item, items)

# compute scale
MAGMOT_pre$NeedForCognition <- MAGMOT_pre$NeedForCognition.1 + MAGMOT_pre$NeedForCognition.2 + MAGMOT_pre$NeedForCognition.3R + MAGMOT_pre$NeedForCognition.4R + MAGMOT_pre$NeedForCognition.5R + MAGMOT_pre$NeedForCognition.6 +
  MAGMOT_pre$NeedForCognition.7 + MAGMOT_pre$NeedForCognition.8R + MAGMOT_pre$NeedForCognition.9R + MAGMOT_pre$NeedForCognition.10 + MAGMOT_pre$NeedForCognition.11 + MAGMOT_pre$NeedForCognition.12R + MAGMOT_pre$NeedForCognition.13 +
  MAGMOT_pre$NeedForCognition.14 + MAGMOT_pre$NeedForCognition.15 + MAGMOT_pre$NeedForCognition.16R + MAGMOT_pre$NeedForCognition.17R

### Fear of Failure compute scale
# no information available on whether it includes items that need recoding
MAGMOT_pre$FearOfFailure <- MAGMOT_pre$FearOfFailure.1 + MAGMOT_pre$FearOfFailure.2  + MAGMOT_pre$FearOfFailure.3 + MAGMOT_pre$FearOfFailure.4 + MAGMOT_pre$FearOfFailure.5 +
  MAGMOT_pre$FearOfFailure.6 + MAGMOT_pre$FearOfFailure.7 + MAGMOT_pre$FearOfFailure.8 + MAGMOT_pre$FearOfFailure.9

### Approach and Avaoidance Temperament compute scales
MAGMOT_pre$ApproachTemperament <- MAGMOT_pre$ApproachAndAvoidanceTemperament.2 + MAGMOT_pre$ApproachAndAvoidanceTemperament.4 + MAGMOT_pre$ApproachAndAvoidanceTemperament.5 + MAGMOT_pre$ApproachAndAvoidanceTemperament.8 +
  MAGMOT_pre$ApproachAndAvoidanceTemperament.10 + MAGMOT_pre$ApproachAndAvoidanceTemperament.11 #Approach temperament = item2+item4+item5+item8+item10+item11
MAGMOT_pre$AvoidanceTemperament <-  MAGMOT_pre$ApproachAndAvoidanceTemperament.1 + MAGMOT_pre$ApproachAndAvoidanceTemperament.3 + MAGMOT_pre$ApproachAndAvoidanceTemperament.6 + MAGMOT_pre$ApproachAndAvoidanceTemperament.7  +
  MAGMOT_pre$ApproachAndAvoidanceTemperament.9 + MAGMOT_pre$ApproachAndAvoidanceTemperament.12 #Avoidance temperament = item1+item3+item6+item7+item9+item12

### Trait Curiosity: compute scale
MAGMOT_pre$TraitCuriosity <- MAGMOT_pre$TraitCuriosity.1 + MAGMOT_pre$TraitCuriosity.2 + MAGMOT_pre$TraitCuriosity.3 + MAGMOT_pre$TraitCuriosity.4 + MAGMOT_pre$TraitCuriosity.5 + MAGMOT_pre$TraitCuriosity.6 + MAGMOT_pre$TraitCuriosity.7 +
  MAGMOT_pre$TraitCuriosity.8 + MAGMOT_pre$TraitCuriosity.9 + MAGMOT_pre$TraitCuriosity.10 + MAGMOT_pre$TraitCuriosity.11 + MAGMOT_pre$TraitCuriosity.12 + MAGMOT_pre$TraitCuriosity.13 + MAGMOT_pre$TraitCuriosity.14 +
  MAGMOT_pre$TraitCuriosity.15 + MAGMOT_pre$TraitCuriosity.16 + MAGMOT_pre$TraitCuriosity.17 + MAGMOT_pre$TraitCuriosity.18 + MAGMOT_pre$TraitCuriosity.19 + MAGMOT_pre$TraitCuriosity.20

### reduce MAGMOT_pre to relevant variables
MAGMOT_pre <- MAGMOT_pre[, c("ID", "preFile", "startPre", "endPre", "corsi.1", "X2nback.1", "age", "DOB",   "gender", "ethnicity", "education", "yearsOfEducation", 
                             "employment", "studySubject", "english", "AgeEnglishAcquisition", "handedness", "vision", "health", "neurodisorders", "screening_MRI",
                             "BIS", "BAS_rewardresponsiveness", "BAS_drive", "BAS_funseeking", "NeedForCognition", "FearOfFailure", "ApproachTemperament", "AvoidanceTemperament", "TraitCuriosity")]


#################################################### get CORSI data ####################################################
# check whether there are any files missing
MAGMOT_pre$ID[is.na(MAGMOT_pre$corsi.1)]

# convert corsi.1  to character
MAGMOT_pre$corsi.1 <- as.character(MAGMOT_pre$corsi.1)

# create list of corsi files
corsi_filelist <- as.character(MAGMOT_pre$corsi.1)

# create empty data frame corsiData
corsiData <- data.frame(matrix(nrow=length(corsi_filelist), ncol=1))
names(corsiData) <- c("corsi.1")

# replace subject data in created data frame
for (i in seq_along(corsi_filelist)){
  
  if(is.na(corsi_filelist[i]) == T | corsi_filelist[i] == "missing"){ # not all participants have complete data
    print(paste("There is a corsi file missing for ID", MAGMOT_pre$ID[i] ))
  } else {
    corsi <- read.table(paste(corsi_filelist[i]))
    #V1: highest corsi span thus far, V2: number of items to remember, V3: corret/incorrect, V4: random presentation
    
    corsiSpanMax <- corsi[nrow(corsi),1]
    #print(corsi_filelist[i])
    corsiData[i, "corsi.1"] <- corsi_filelist[i]
    corsiData[i, "corsiSpan"] <- corsiSpanMax
    rm(corsi, corsiSpanMax)
  }
}

# delete all rows in corsiData with NAs (these are the participants that do not have data)
corsiData <- na.omit(corsiData)

# merge data from MAGMOT_pre and corsi
MAGMOT <- merge(MAGMOT_pre, corsiData, by = "corsi.1", all.x = T)
rm(corsiData, corsi_filelist, i)

#################################################### get NBACK data ####################################################
# check whether there are any files missing
MAGMOT_pre$ID[is.na(MAGMOT_pre$X2nback.1)]

# convert corsi.1  to character
MAGMOT_pre$X2nback.1 <- as.character(MAGMOT_pre$X2nback.1)

# create list of nback files
nback_filelist <- as.character(MAGMOT_pre$X2nback.1)

# create empty data frame nbackData
nbackData <- data.frame(matrix(nrow=length(nback_filelist), ncol=1))
names(nbackData) <- c("X2nback.1")

# replace subject data in created data frame
for (i in seq_along(nback_filelist)){
  if(is.na(nback_filelist[i]) == T){ # not all participants have complete data
    print(paste("There is a nback file missing for ID", MAGMOT_pre$ID[i] ))
  } else {
    nback <- read.table(paste(nback_filelist[i]))
    # V1: name of block,V2: correct (1=correct, 2=wrong, 3=too slow), V3: which key was pressed, V4: RT in ms
    # V5: random number used for conditions (1=same as 2-back, 2-5 other letter), V6: trial number (per block), V7: current letter, V8: letter in previous trial
    nback <- subset(nback, nback$V1 != "training") # subset nback to only include data from task block
    
    nbackData[i, "X2nback.1"] <- nback_filelist[i]
    nbackData[i, "NBacks"] <- sum(nback$V5 == 1)
    nbackData[i, "nonNBacks"] <- sum(nback$V5 != 1)
    
    # subset data set so that it only contains the data of each of signal detection theory (SDT) cells
    nbackData[i, "nback_hits"] <-  nrow(subset(nback, nback$V2 == 1 & nback$V5 == 1))
    nbackData[i, "nback_misses_inclTooSlow"] <-  nrow( subset(nback, nback$V2 != 1 & nback$V5 == 1))
    nbackData[i, "nback_misses_exclTooSlow"] <-  nrow(subset(nback, nback$V2 == 2 & nback$V5 == 1))
    nbackData[i, "nback_correctrejections"] <-  nrow(subset(nback, nback$V2 == 1 & nback$V5 > 1))
    nbackData[i, "nback_falsealarms_inclTooSlow"] <-  nrow(subset(nback, nback$V2 != 1 & nback$V5 > 1)) # wrong/too slow response when letter is not the same as 2-back
    nbackData[i, "nback_falsealarms_exclTooSlow"] <-  nrow(subset(nback, nback$V2 == 2 & nback$V5 > 1)) # wrong response when letter is not the same as 2-back
    
    # compute some SDT measures
    nbackData[i, "nback_hitrate"] <-  nbackData[i, "nback_hits"] / nbackData[i, "NBacks"]  # percentage hits in actual nBack trials
    nbackData[i, "nback_falsealarmrate"] <-  nbackData[i, "nback_falsealarms_exclTooSlow"] / nbackData[i, "nonNBacks"] # percentage false alarms in actual no nBack trials
    nbackData[i, "nback_accurary"] <-  nbackData[i, "nback_hitrate"] - nbackData[i, "nback_falsealarmrate"] # hit rate corrected for false alarms
    
    rm(nback)
  }
}

# delete all rows in nbackData with NAs (these are the participants that do not have data)
nbackData <- na.omit(nbackData)

# merge data from MAGMOT (= MAGMOT_pre and corsi) and nback
MAGMOT <- merge(MAGMOT, nbackData, by = "X2nback.1", all.x = T)
rm(nbackData, nback_filelist, MAGMOT_pre,i)

# rename corsi and nback
names(MAGMOT)[names(MAGMOT)=="X2nback.1"] <- "nbackFile"
names(MAGMOT)[names(MAGMOT)=="corsi.1"] <- "corsiFile"

# reorder MAGMOT
MAGMOT <- MAGMOT[,c(grep("ID", colnames(MAGMOT)), grep("corsiFile", colnames(MAGMOT)), grep("nbackFile", colnames(MAGMOT)), 4:ncol(MAGMOT))]


#################################################### read in MAGMOT_post ####################################################
setwd(postDir)
MAGMOT_post <- xlsx::read.xlsx("data.xlsx", sheetName="data")
MAGMOT_post <- subset(MAGMOT_post, !is.na(MAGMOT_post$TIME_end)) # subset data to only include complete data sets
MAGMOT_post <-  subset(MAGMOT_post, nchar(as.character(MAGMOT_post$id.1)) <= 2) # remove test data

# change participant IDs
MAGMOT_post$id.1 <- gsub("O", "", MAGMOT_post$id.1) #replace any O in the file names
MAGMOT_post$id.1 <- ifelse(nchar(MAGMOT_post$id.1)==1, paste0("0", MAGMOT_post$id.1), MAGMOT_post$id.1) # add a zeo if not present in front of one digit numbers

# delete reassigned ID's
participantToDelete1 <- as.character(MAGMOT_post$participant[(MAGMOT_post$id.1 == "09" & MAGMOT_post$TIME_start == "2019-05-13-19-20")]) # ppt born on 19/03/98 should be kept
participantToDelete2 <- as.character(MAGMOT_post$participant[(MAGMOT_post$id.1 == 23 & MAGMOT_post$TIME_start == "2019-05-22-11-51")]) # ppt born on 17/12/91 should be deleted
participantToDelete3 <- as.character(MAGMOT_post$participant[(MAGMOT_post$id.1 == 32 & MAGMOT_post$TIME_start == "2019-06-03-15-42")]) # ppt born on 10/03/97 should be deleted
participantToDelete4 <- as.character(MAGMOT_post$participant[(MAGMOT_post$id.1 == 46 & MAGMOT_post$TIME_start == "2019-06-26-19-36")]) # ppt born on 15/04/96 should be deleted
participantToDelete5 <- as.character(MAGMOT_post$participant[(MAGMOT_post$id.1 == 48 & MAGMOT_post$TIME_start == "2019-07-03-12-48")]) # ppt born on 15/11/99 should be deleted

MAGMOT_post <- subset(MAGMOT_post, MAGMOT_post$participant != participantToDelete1 & MAGMOT_post$participant != participantToDelete2 & MAGMOT_post$participant != participantToDelete3 & MAGMOT_post$participant != participantToDelete4 & MAGMOT_post$participant != participantToDelete5)
rm(list=ls(pattern = "participantToDelete"))
### recode demographic information and then delete the old variable
MAGMOT_post$group <- ifelse(MAGMOT_post$group.1 == 1, "int", # group
                            ifelse(MAGMOT_post$group.1 == 2, "ext", NA))
MAGMOT_post$group.1 <- NULL
MAGMOT_post$groupEffectCoded <-  ifelse(MAGMOT_post$group == "ext", 1, -1)

MAGMOT_post$alcohol <- ifelse(MAGMOT_post$alcohol.1 == 1, "alcohol in the last 24 hours", # alcohol
                              ifelse(MAGMOT_post$alcohol.1 == 2, "no alcohol in the last 24 hours", NA))
MAGMOT_post$alcohol.1 <- NULL

### recode and rename checklist items
items <- c("checklist.1", "checklist.2", "checklist.3", "checklist.4", "checklist.5", "checklist.6", "checklist.7", "checklist.8", "checklist.9", "checklist.10", "checklist.11")
itemsRenamed <-c("eyetracking","fieldmap","preLearningRest","taskBlock1","taskBlock2","taskBlock3","postLearningRest","T1w","eyetrackingData","taskData","questionnaireData")
for (item in seq_along(items)){
  item_new <-  paste0(itemsRenamed[item])
  MAGMOT_post[,item_new] <- ifelse(MAGMOT_post[[paste0(items[item])]] == 0, "no",
                                   ifelse(MAGMOT_post[[paste0(items[item])]] == 1, "yes",
                                          ifelse(MAGMOT_post[[paste0(items[item])]] == 999, "noInformation",
                                                 NA)))
  MAGMOT_post[[paste0(items[item])]] <- NULL
  rm(item_new)
}
rm(item, items, itemsRenamed)

### rename variables
names(MAGMOT_post)[names(MAGMOT_post)=="id.1"] <- "ID"
names(MAGMOT_post)[names(MAGMOT_post)=="sleep_lastnight.1"] <- "sleepLastNight"
names(MAGMOT_post)[names(MAGMOT_post)=="sleep_average.1"] <- "sleepAverage"
names(MAGMOT_post)[names(MAGMOT_post)=="alcohol_amount.1"] <- "alcoholAmount"
names(MAGMOT_post)[names(MAGMOT_post)=="rewardEffort.1"] <- "rewardEffort"
names(MAGMOT_post)[names(MAGMOT_post)=="rewardExpectations.1"] <- "rewardExpectations"
names(MAGMOT_post)[names(MAGMOT_post)=="TIME_start"] <- "startPost"
names(MAGMOT_post)[names(MAGMOT_post)=="TIME_end"] <- "endPost"
names(MAGMOT_post)[names(MAGMOT_post)=="participant"] <- "postFile"
names(MAGMOT_post)[names(MAGMOT_post)=="comments_participant.1"] <- "comment_ppt1"
names(MAGMOT_post)[names(MAGMOT_post)=="comments_participant.2"] <- "comment_ppt2"
names(MAGMOT_post)[names(MAGMOT_post)=="comments_participant.3"] <- "comment_ppt3"
names(MAGMOT_post)[names(MAGMOT_post)=="comments_experimenter.1"] <- "comment_exp"

### State Curiosity: compute scale
MAGMOT_post$StateCuriosity <- MAGMOT_post$StateCuriosity.1 + MAGMOT_post$StateCuriosity.2 + MAGMOT_post$StateCuriosity.3 + MAGMOT_post$StateCuriosity.4 + MAGMOT_post$StateCuriosity.5 + MAGMOT_post$StateCuriosity.6 + MAGMOT_post$StateCuriosity.7 +
  MAGMOT_post$StateCuriosity.8 + MAGMOT_post$StateCuriosity.9 + MAGMOT_post$StateCuriosity.10 + MAGMOT_post$StateCuriosity.11 + MAGMOT_post$StateCuriosity.12 + MAGMOT_post$StateCuriosity.13 + MAGMOT_post$StateCuriosity.14 +
  MAGMOT_post$StateCuriosity.15 + MAGMOT_post$StateCuriosity.16 + MAGMOT_post$StateCuriosity.17 + MAGMOT_post$StateCuriosity.18 + MAGMOT_post$StateCuriosity.19 + MAGMOT_post$StateCuriosity.20

# select relevant rows from MAGMOT_post
MAGMOT_post <- MAGMOT_post[,c("ID", "postFile", "startPost", "endPost", 
               "group", "groupEffectCoded", "StateCuriosity", 
               "sleepLastNight", "sleepAverage", 
               "alcohol", "alcoholAmount", "rewardEffort", "rewardExpectations", 
               "comment_ppt1", "comment_ppt2", "comment_ppt3", "comment_exp",
               "eyetracking", "fieldmap", "preLearningRest", "taskBlock1", "taskBlock2", "taskBlock3", "postLearningRest", "T1w",
               "eyetrackingData", "taskData", "questionnaireData" 
               )]

# merge data from MAGMOT (= MAGMOT_pre, corsi, nback) and MAGMOT_post
MAGMOT <- merge(MAGMOT, MAGMOT_post, by = "ID", all = T)
rm(MAGMOT_post)


########################################################################################################################################
####################################################  PROCESS THE DATA SUBJECTWISE  ####################################################
########################################################################################################################################
subjects  <- as.character(MAGMOT$ID)

for (s in seq_along(subjects)){
  
  # create BIDS string
  # if (nchar(subjects[s])==1){ # add a zeo if not present in front of one digit numbers
  #   subjects[s] <- paste0("0", subjects[s])
  # }
  if (MAGMOT$group[MAGMOT$ID == subjects[s]] == "ext"){ # define BIDS name depending on group
    BIDSstring = paste0("sub-experimental0", subjects[s])
  } else if (MAGMOT$group[MAGMOT$ID == subjects[s]] == "int"){
    BIDSstring = paste0("sub-control0", subjects[s])
  }
  print(paste("processing data for",BIDSstring))
  
  # create directory to save preprocessed events.tsv files
  preprocessedEventsSubjDir <- file.path(preprocessedEventsRootDir, BIDSstring)
  ifelse(!dir.exists(preprocessedEventsSubjDir), dir.create(preprocessedEventsSubjDir), FALSE)
  
  
  #################################################### read in questionnaire data  ####################################################
  MRIdataDir <- file.path(mainDir, "PsychToolBox_script", "behavioural_data", paste0(subjects[s]))
  setwd(MRIdataDir)
  quest <- rmatio::read.mat(paste0("MAGMOT_", subjects[s], "_questionnaire.mat"))
  respMatQuest<- quest$respMatQuest
  quest <- data.frame(matrix(unlist(respMatQuest), nrow=length(respMatQuest), byrow=T),stringsAsFactors=FALSE)
  
  # rename variables questionnaire
  names(quest) <- c("ID","fMRI","group","motivation","q_trial", "question",
                    "startValueQuestion","clicksQuestion","endValue","displayQuestionOnset","timestampQuestion","rtQuestion")
  
  # assign names to question
  quest$IMI <- ifelse(quest$question == "It was fun to do the experiment.", "post1",
                      ifelse(quest$question == "It was boring to do the experiment.", "post2",
                             ifelse(quest$question == "It was enjoyable to do the experiment.", "post3",
                                    ifelse(quest$question == "I was totally absorbed in the experiment.", "post4",
                                           ifelse(quest$question == "I lost track of time.", "post5",
                                                  ifelse(quest$question == "I concentrated on the experiment.", "post6",
                                                         ifelse(quest$question == "The task was interesting.", "post7",
                                                                ifelse(quest$question == "I liked the experiment.", "post8",
                                                                       ifelse(quest$question == "I found working on the task interesting.", "post9",
                                                                              ifelse(quest$question == "The experiment bored me.", "post10",
                                                                                     ifelse(quest$question == "I found the experiment fairly dull.", "post11",
                                                                                            ifelse(quest$question == "I got bored." , "post12",
                                                                                                   ifelse(quest$question == "I put a lot of effort into this.", "post13",
                                                                                                          ifelse(quest$question == "I did not try very hard\\nto do well at this activity.", "post14",
                                                                                                                 ifelse(quest$question == "I tried very hard on this activity.", "post15",
                                                                                                                        ifelse(quest$question == "It was important to me to do well at this task.", "post16",
                                                                                                                               ifelse(quest$question == "I did not put much energy into this.", "post17",
                                                                                                                                      ifelse(quest$question == "I did not feel nervous at all while doing this.", "post18",
                                                                                                                                             ifelse(quest$question == "I felt very tense while doing this activity." , "post19",
                                                                                                                                                    ifelse(quest$question == "I was very relaxed in doing this experiment.", "post20",
                                                                                                                                                           ifelse(quest$question == "I was anxious while working on this task." , "post21",
                                                                                                                                                                  ifelse(quest$question == "I felt pressured while doing this task.", "post22",
                                                                                                                                                                         ifelse(quest$question == "I tried to find out how many people\\nwill be able to find the solution.", "post23",
                                                                                                                                                                                ifelse(quest$question == "I was able to see the magic tricks properly.", "post24",
                                                                                                                                                                                       NA))))))))))))))))))))))))
  
  
  
  
  # recode items questionnaire
  items <- c("post2", "post14", "post17", "post18", "post20") #items to recode
  quest$score <- quest$endValue
  for (item in items){
    quest$score[quest$IMI == item] <- ifelse(quest$endValue[quest$IMI == item] == 7, 1,
                                             ifelse(quest$endValue[quest$IMI == item] == 6, 2,
                                                    ifelse(quest$endValue[quest$IMI == item] == 5, 3,
                                                           ifelse(quest$endValue[quest$IMI == item] == 4, 4,
                                                                  ifelse(quest$endValue[quest$IMI == item] == 3, 5,
                                                                         ifelse(quest$endValue[quest$IMI == item] == 2, 6,
                                                                                ifelse(quest$endValue[quest$IMI == item] == 1, 7, NA)))))))
    
  }
  rm(item, items)
  
  # save questionnaire data in long format
  setwd(preprocessedQuestDir)
  filenameQuest <- paste0(BIDSstring,"_questionnaire.tsv")
  write.table(quest, file = filenameQuest, sep = "\t", row.names = F, quote = F)
  
  # reduce questionnaire data set and convert data in wide format
  questLong <- quest[,c("ID", "fMRI", "group", "motivation", "IMI", "score")]
  questWide <- reshape2::dcast(questLong, ID + fMRI + group + motivation ~ IMI, value.var="score")
  
  # rbind the questWide files of each subject to a data frame
  if(s == 1){
    questDataWide <- questWide
  } else {
    temp_questDataWide <-  questWide
    questDataWide <- rbind.all.columns(questDataWide, temp_questDataWide) #rbind all columns will induce NA if there was initially no data saved in the loop per participant
    rm(temp_questDataWide)
  }
  if (s == length(subjects)){
    questDataWide[,grep("post", colnames(questDataWide))] <- as.numeric(unlist(questDataWide[,grep("post", colnames(questDataWide))])) # make sure the scores are numeric
    questDataWide$intrinsicMotivation <- (questDataWide$post1 + questDataWide$post2 + questDataWide$post3)/3 ###intrinsic motivation items
    questDataWide$taskEngagement <- (questDataWide$post4 + questDataWide$post5 + questDataWide$post6)/3 ###task engagement items
    questDataWide$interest <- (questDataWide$post7 + questDataWide$post8 + questDataWide$post9)/3 ###interest items
    questDataWide$boredom <- (questDataWide$post10 + questDataWide$post11 + questDataWide$post12)/3 ###boredom items
    questDataWide$effort <- (questDataWide$post13 + questDataWide$post14 + questDataWide$post15 + questDataWide$post16 + questDataWide$post17)/5 ####effort/importance
    questDataWide$pressure <- (questDataWide$post18 + questDataWide$post19 + questDataWide$post20 + questDataWide$post21 + questDataWide$post22)/5 ###pressure/tension
  }
  
  # remove unneccsary values and data
  if (debug == 0){
    rm(quest, questLong, questWide, respMatQuest)
  }
  
  #################################################### read in task data  ####################################################
  MRIdataDir <- file.path(mainDir, "PsychToolBox_script", "behavioural_data", paste0(subjects[s]))
  setwd(MRIdataDir)
  taskfile <- list.files(pattern = glob2rx("MAGMOT*task_inclTimingCorrection*txt"))
  ptbdata <- read.table(file = taskfile, sep = "\t", header = T)
  
  # change participant IDs
  names(ptbdata)[names(ptbdata)=="subject"] <- "ID"
  ptbdata$ID <- gsub("O", "", ptbdata$ID) #replace any O in the file names
  ptbdata$ID <- ifelse(nchar(ptbdata$ID)==1, paste0("0", ptbdata$ID), ptbdata$ID) # add a zeo if not present in front of one digit numbers
  ptbdata$BIDS <- BIDSstring
  ptbdata$groupEffectCoded <-  ifelse(ptbdata$group == "ext", 1, -1)
  
  
  # change stimID of respMatTask to match collector
  names(ptbdata)[names(ptbdata)=="stimID"] <- "vidFileName"
  ptbdata$stimID <- gsub("^\\d\\d_", "", ptbdata$vidFileName)
  ptbdata$stimID <- gsub("_combined_small.mp4", "", ptbdata$stimID)
  
  # add curiosity median split for MAGMOT sample (computed at beginning of script based on old MAGMOT_long.xlsx data set)
  ptbdata <- merge(ptbdata, dfMeans, by = "stimID")
  
  # recode curiosity: replace all curiosity ratings that have exceeded the timeout 
  # NOTE::: this means we should CHANGE/REPLACE responseCuriosity with  curiosity in the following, NOT DONE YET
  ptbdata$curiosity <- ifelse(ptbdata$timestampCuriosity > ptbdata$displayCuriosityOnset + ptbdata$timeoutCuriosity - 3*ptbdata$timingCorrection, NA, ptbdata$responseCuriosity)
  ptbdata$curiosity_RT <- ifelse(ptbdata$timestampCuriosity > ptbdata$displayCuriosityOnset + ptbdata$timeoutCuriosity - 3*ptbdata$timingCorrection, NA, ptbdata$rtCuriosity)
  curiosityNAs <- sum(is.na(ptbdata$curiosity))
  if (feedback == "yes"){
    print(paste("number of NAs in curiosity", curiosityNAs))
  }
  
  # compute median split (WITHIN SUBJECT) for curiousity rating
  ptbdata$mediansplitCuriosityWithinSubject <- ifelse(ptbdata$responseCuriosity > median(ptbdata$responseCuriosity), "above", 
                                                      ifelse(ptbdata$responseCuriosity < median(ptbdata$responseCuriosity), "below", "median"))
  ptbdata$mediansplitCuriosityWithinSubject_updated <- ifelse(ptbdata$curiosity > median(ptbdata$curiosity), "above", 
                                                              ifelse(ptbdata$curiosity < median(ptbdata$curiosity), "below", "median"))
  
  # compute group-mean centered curiosity rating
  ptbdata$curiosityGroupMeanCentered <- ptbdata$responseCuriosity - mean(ptbdata$responseCuriosity, na.rm = T)
  ptbdata$curiosityGroupMeanCentered_updated <- ptbdata$curiosity - mean(ptbdata$curiosity, na.rm = T)
  
  # compute reward by curiosity interaction
  ptbdata$rewardByCuriosity <- ptbdata$curiosityGroupMeanCentered * ptbdata$groupEffectCoded
  ptbdata$rewardByCuriosity_updated <- ptbdata$curiosityGroupMeanCentered_updated * ptbdata$groupEffectCoded
  
  # define number of blocks
  numBlocks <- max(ptbdata$block)
  
  # order MEMO according to trial number
  ptbdata <- ptbdata[order(ptbdata$trial),] 
  
  # define number of acquisitions
  numAcqTotal <- length(which(is.na(ptbdata$endBlock) == F))
  if (numAcqTotal > numBlocks) { # variable acq is added for each participant rather than just for those
    print(paste (subjects[s], "has", numAcqTotal, "acquisitions"))
  }
  
  # for each of the task blocks, create the variable acq reflecting the scanner runs of that block
  for (b in 1:numBlocks) {
    
    # subset to only look at the data for each run
    run <- subset(ptbdata, ptbdata$block == b)
    
    # how many scanner runs relate to this block?
    numAcqBlock <- length(which(is.na(run$endBlock) == F))
    
    # add the variable acq to the data set
    endOfAcqTrials <- which(is.na(run$endBlock) == F) # determine number of endBlock values (= end of a scanning block)
    for (t in 1:length(endOfAcqTrials)) {
      # define first and last trial of each scanner acquisition
      if (t == 1){
        startOfAcq <- 1
      } else {
        startOfAcq <- endOfAcq + 1
      }
      endOfAcq <- endOfAcqTrials[t]
      for (row in 1:nrow(run)) {
        # if the current row falls within start and end trial, add the acquistion as t
        if (row >= startOfAcq & row <= endOfAcq) {
          run$acq[row] <- t 
        }
      }
    }
    if (b == 1){
      task <- run
    } else {
      task <- rbind(task, run)
    }
  }
  
  # rearrange culumns in data frame: acq is currently the last row, move it next to block
  col_acq <- which(names(task) == "acq")
  col_block <- which(names(task) == "block")
  task <- task[,c(1:col_block, col_acq, (col_block+1):(ncol(task)-1))]
  
  
  
  #################################################### read in memory data  ####################################################
  setwd(dataMemoryDir)
  
  f <- paste0("magicmemory_", version, "_", subjects[s], ".csv")
  if (file.exists(f)){
    memory <- read.csv(f, header = T)
    memory$memoryFile <- f
  } else { # if file does not exist, print into console and check other spelling of file
    print(paste0("magicmemory_", version, "_", subjects[s], ".csv does not exist"))
    
    subjectsAlt <- gsub("0", "", subjects[s])
    f <- paste0("magicmemory_", version, "_", subjectsAlt, ".csv")
    if (file.exists(f)){
      memory <- read.csv(f, header = T)
      memory$memoryFile <- f
    } else { # if file does not exist, print into console and use filler file to create NAs
      print(paste0("magicmemory_", version, "_", subjectsAlt, ".csv does not exist either"))
      rm(subjectsAlt)
      
      # CREATE MEMORY FILLER!! actually, that might not even be necessary
    }
  }
  
  # overwrite username to ensure that it is the same across all data sets
  memory$username <- subjects[s]
  
  
  if("post_0_trial_end_date" %in% colnames(memory)){  # check whether the memory data set has information about when it has started/finished
    memory$startMemory <- memory$post_0_trial_end_date[2]
    memory$endMemory <- memory$post_0_trial_end_date[dim(memory)[2]]
  } else { # if not, add a note
    memory$startMemory <- "check manually"
    memory$endMemory <- "check manually"
  }
  
  
  #################################################### process the RECALL memory data: select relevant rows and columns, change item counter
  cuedRecall <- subset(memory, memory$trial.type == "magic_recall")
  cuedRecall$itemRecall <- cuedRecall$item-1
  cuedRecall$trialRecall <- c(1:dim(cuedRecall)[1])
  cuedRecall <- cuedRecall[,c("username", "trialRecall", "itemRecall", "stimid", "trickAnswer")]
  names(cuedRecall) <- c("ID", "trialRecall", "itemRecall", "stimID", "description")
  
  # pre code answers for recall task: in the experiment, participants are asked to insert "no recall" in the field if they cannot recall the magic trick
  # if that has occured, recall performance on that magic trick is coded with 0
  for (j in 1:nrow(cuedRecall)){
    if (is.na(cuedRecall$description[j]) == T) {
      cuedRecall$cuedRecallStrict[j] = NaN
      cuedRecall$cuedRecallLenient[j] = NaN
    } else if (cuedRecall$description[j] == "no recall") { # check different spellings of NO RECALL
      cuedRecall$cuedRecallStrict[j] = 0
      cuedRecall$cuedRecallLenient[j] = 0
    } else if (cuedRecall$description[j] == "No recall") {
      cuedRecall$cuedRecallStrict[j] = 0
      cuedRecall$cuedRecallLenient[j] = 0
    } else if (cuedRecall$description[j] == "No Recall") {
      cuedRecall$cuedRecallStrict[j] = 0
      cuedRecall$cuedRecallLenient[j] = 0
    } else if (cuedRecall$description[j] == " No recall") {
      cuedRecall$cuedRecallStrict[j] = 0
      cuedRecall$cuedRecallLenient[j] = 0
    } else if (cuedRecall$description[j] == "NO recall") {
      cuedRecall$cuedRecallStrict[j] = 0
      cuedRecall$cuedRecallLenient[j] = 0
    } else { # if none of the spelling matches, insert NA --> manual coding necessary
      cuedRecall$cuedRecallStrict[j] = NA
      cuedRecall$cuedRecallLenient[j] = NA
    }
  }
  
  # save recall memory data if the participants initially participated in the memory part
  if (file.exists(file.path(dataMemoryDir, f))) { # f = magicmemory_fmri.csv
    # save recall file
    setwd(preprocessedMemoryDir)
    filenameRecall <- paste0(BIDSstring,"_cuedRecall.csv")
    write.table(cuedRecall, file = filenameRecall, sep = ",", row.names = F)
    
    # check whether a coded recall file exists; NOTE: THE STRING MIGHT CHANGE!
    setwd(codedDir)
    #f_coded <- list.files(pattern = glob2rx("sub*cuedRecall_CP.csv"))# the suffix _SM indicated that I have coded the remaining answers in the recall task
    f_coded <-  paste0(BIDSstring,"_cuedRecall_CP.csv") # the suffix _SM indicated that I have coded the remaining answers in the recall task
    if (file.exists(f_coded)) { #if the recall performance has already been coded, the recall df gets overwritten
      cuedRecall_coded <- read.csv(f_coded, header = T)
      cuedRecall_coded <- cuedRecall_coded[1:nrow(ptbdata),] # reduce cuedRecall_coded to only those rows that relate to stimuli
      jj <- 0 #define a variable jj as zero to be used in the next statement
      for (j in 1:nrow(cuedRecall_coded)){ # Cristina has inserted some FALSE as coding, this will be replaced with NA to make sure that the script still runs
        if (is.na(cuedRecall_coded$cuedRecallStrict[j] == TRUE)) {
          jj <- jj+1 #update jj
          if(jj == 1){
            print(paste(f_coded, "has NA in cuedRecallStrict"))
          }
        } else if(cuedRecall_coded$cuedRecallStrict[j] == "FALSE") {
          cuedRecall_coded$cuedRecallStrict[j] = NA
          cuedRecall_coded$cuedRecallLenient[j] = NA
        }
      }
      names(cuedRecall_coded)[names(cuedRecall_coded)=="Username"] <- "ID" # change Username to ID
      cuedRecall_coded$ID <- subjects[s]
      cuedRecall <- cuedRecall_coded # overwrite cuedRecall
    }
  } else { # if there is no _cuedRecall_CP file ALTHOUGH there was a memory data file, print into console
    print(paste(f_coded, "does not exist"))
  }
  
  ####################################################  process RECOGNITION memory data: select relevant rows and columns, change item counter
  recognition <- subset(memory, memory$trial.type == "magic_recognition")
  recognition$itemRecognition <- recognition$item-1
  recognition$trialRecognition <- c(1:dim(recognition)[1])
  
  # check whether the answer selected by participants is the correct answer
  recognition$option1 <- as.character(recognition$option1)
  recognition$Recognition_chosen <- as.character(recognition$Recognition_chosen)
  recognition$recognition <- ifelse(is.na(recognition$Recognition_chosen) == TRUE, NA,
                                    ifelse(recognition$option1 == recognition$Recognition_chosen, 1, 0))  # in the collector sheet, the correct answer is always in the option1 column
  
  # compute group-mean centered confidence rating
  recognition$confidenceGroupMeanCentered <- recognition$Confidence - mean(recognition$Confidence, na.rm = T)
  
  # reduce and rename data set
  recognition <- recognition[,c("username", "trialRecognition", "itemRecognition", "stimid", "Recognition_chosen", "recognition_RT", "recognition" , "Confidence", "confidenceGroupMeanCentered", "confidence_RT")]
  names(recognition) <- c("ID", "trialRecognition", "itemRecognition", "stimID", "answer", "answerRT", "recognition", "confidence", "confidenceGroupMeanCentered", "confidenceRT")
  
  # compute scores of recognition performance
  recognition$confidenceCorrectTrials <- ifelse(recognition$recognition == 1, recognition$confidence, NA)
  recognition$confidenceGroupMeanCenteredCorrectTrials <- ifelse(recognition$recognition == 1, recognition$confidenceGroupMeanCentered, NA)
  recognition$recognitionAboveMeanConf <- ifelse(recognition$recognition == 1 & recognition$confidenceGroupMeanCentered > 0, 1, 0)

  for (k in 1:6) { #confidence ranges from 1 to 6, potentially code can be made more flexible by using min(data$confidence) and max(data$confidence)
    recognition[[paste0("recognitionConfLevel_", k)]] <- ifelse(recognition$confidence == k & recognition$recognition == 1, 1, 0)
    if (k < 6) {
      recognition[[paste0("recognitionConfLevel_above_", k)]] <- ifelse(recognition$confidence > k & recognition$recognition == 1, 1, 0)
    }
    
    if (k == 1 || k == 3 || k == 5) {
      recognition[[paste0("recognitionConfLevel_", k, "_", k+1)]] <- ifelse(recognition$confidence == k  & recognition$recognition == 1 | recognition$confidence == k+1  & recognition$recognition == 1, 1, 0)
    }
    if (k == 1 || k == 4) {
      recognition[[paste0("recognitionConfLevel_", k, "_", k+1, "_", k+2)]] <- ifelse(recognition$confidence == k  & recognition$recognition == 1 | recognition$confidence == k+1  & recognition$recognition == 1 | recognition$confidence == k+2  & recognition$recognition == 1, 1, 0)
    }
  }
  
  # save recognition data
  setwd(dataMemoryDir)
  if (file.exists(f)) { # if the participants initially participated in the memory part, their preprocessed recall data is saved
    setwd(preprocessedMemoryDir)
    filenameRecog <- paste0(BIDSstring,"_recognition.csv")
    write.table(recognition, file = filenameRecog, sep = ",", row.names = F)
  }
  
  # MERGE DATA FROM TASK AND MEMORY TEST
  MEMO <- merge(task, cuedRecall, by = c("ID", "stimID"))
  MEMO <- merge(MEMO, recognition, by = c("ID", "stimID"))
  if (debug == 0){
    rm(cuedRecall, recognition, task)
  }
  
  # compute whether a trick has been remembered based on Hasson et al. (2008).  
  # They classified an event to be remembered if there was a correct answer using either recall or high confidence recognition.
  MEMO$rememberedStrict <- ifelse(MEMO$cuedRecallStrict == 1 | MEMO$recognitionAboveMeanConf == 1, 1, 0)
  MEMO$rememberedLenient <- ifelse(MEMO$cuedRecallLenient == 1 | MEMO$recognitionAboveMeanConf == 1, 1, 0)
  
  # calculate curiosity-driven memory memory benefit (continouos)
  MEMO$curiosityBenefit_cuedRecallStrict <- MEMO$curiosityGroupMeanCentered*MEMO$cuedRecallStrict
  MEMO$curiosityBenefit_cuedRecallLenient <- MEMO$curiosityGroupMeanCentered*MEMO$cuedRecallLenient
  MEMO$curiosityBenefit_allConf <- MEMO$curiosityGroupMeanCentered*MEMO$recognition
  MEMO$curiosityBenefit_highConf <- MEMO$curiosityGroupMeanCentered*MEMO$recognitionConfLevel_4_5_6
  MEMO$curiosityBenefit_aboveAvgConf <- MEMO$curiosityGroupMeanCentered*MEMO$recognitionAboveMeanConf
  MEMO$curiosityBenefit_rememberedStrict <- MEMO$curiosityGroupMeanCentered*MEMO$rememberedStrict
  MEMO$curiosityBenefit_rememberedLenient <- MEMO$curiosityGroupMeanCentered*MEMO$rememberedLenient

  # calculate curiosity-driven memory memory benefit (dichotomous)
  MEMO$curiosity_dichotom <- ifelse(MEMO$curiosityGroupMeanCentered > 0, 1,
                                    ifelse(MEMO$curiosityGroupMeanCentered < 0, -1, NA))
  MEMO$curiosityBenefit_cuedRecallStrict_dichotom <- MEMO$curiosity_dichotom*MEMO$cuedRecallStrict
  MEMO$curiosityBenefit_cuedRecallLenient_dichotom <- MEMO$curiosity_dichotom*MEMO$cuedRecallLenient
  MEMO$curiosityBenefit_allConf_dichotom <- MEMO$curiosity_dichotom*MEMO$recognition
  MEMO$curiosityBenefit_highConf_dichotom <- MEMO$curiosity_dichotom*MEMO$recognitionConfLevel_4_5_6
  MEMO$curiosityBenefit_aboveAvgConf_dichotom <- MEMO$curiosity_dichotom*MEMO$recognitionAboveMeanConf
  MEMO$curiosityBenefit_rememberedStrict_dichotom <- MEMO$curiosity_dichotom*MEMO$rememberedStrict
  MEMO$curiosityBenefit_rememberedLenient_dichotom <- MEMO$curiosity_dichotom*MEMO$rememberedLenient
  
  # save data in long format
  setwd(preprocessedLongDir)
  filenameLong <- paste0(BIDSstring,"_long.csv")
  write.table(MEMO, file = filenameLong, sep = ",", row.names = F)
  
  
  ########### CREATE EVENTS.TSV FILES   ########### 
  
  # for each of the task blocks, create the variable acq reflecting the scanner runs of that block
  for (b in 1:numBlocks) {
    
    # subset to only look at the data for each run
    run <- subset(MEMO, MEMO$block == b)
    
    # for each of the acquisitions within the block, save a separate .tsv file FOR BIDS
    for (a in 1:max(run$acq)){
      
      # subset to only look at the data for each acq
      run_acq <- subset(run, run$acq == a)
      
      # replace NAs with n/a
      run_acq$cuedRecallStrict <- ifelse(!is.na(run_acq$cuedRecallStrict), run_acq$cuedRecallStrict, "n/a")
      run_acq$cuedRecallLenient <- ifelse(!is.na(run_acq$cuedRecallLenient), run_acq$cuedRecallLenient, "n/a")
      
      # pick relevant onset variables and create data in long format
      onset <- run_acq[,c("trial","vidFileName", "displayVidOnset", "mockOffset",  "displayVidOffset", "fixationPostVidOnset",
                          "displayAnswerOnset",  "responseAnswer", "timestampAnswer", "fixationPostAnswerOnset",
                          "displayCuriosityOnset", "responseCuriosity", "timestampCuriosity",  "fixationPostCuriosityOnset",
                          "cuedRecallStrict", "cuedRecallLenient", "recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf", "confidence")]
      onset$timestampVidOnset <- onset$displayVidOnset
      onset$timestampMockOffset <- onset$mockOffset
      onset$timestampPostVidFixation <- onset$fixationPostVidOnset
      # onset <- reshape2::melt(onset, id.vars=c("vidFileName", "trial", "mockOffset", "displayVidOffset", "timestampPostVidFixation", "vidDurCalc", "cuedRecallStrict",  "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
      #                                          "responseAnswer", "timestampAnswer", "responseCuriosity", "timestampCuriosity"), value.name = "onset")
      onset <- reshape2::melt(onset, id.vars=c("vidFileName", "trial", "timestampVidOnset", "timestampMockOffset", "displayVidOffset", "timestampPostVidFixation",
                                               "cuedRecallStrict",  "cuedRecallLenient", "recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",  "confidence",
                                               "responseAnswer", "timestampAnswer", "responseCuriosity", "timestampCuriosity"), value.name = "onset")
      
      # recode variable (for later merging)
      # Levels: displayVidOnset mockOffset fixationPostVidOnset displayAnswerOnset fixationPostAnswerOnset displayCuriosityOnset fixationPostCuriosityOnset
      levels(onset$variable) <- c("displayMockVid", "displayVid", "fixationPostVid",
                                  "displayAnswer", "fixationPostAnswer", "displayCuriosity", "fixationPostCuriosity")
      
      # pick relevant duration variables and create data in long format
      duration <- run_acq[,c("trial","vidFileName", "avgVidDur_MAGMOT",
                             "displayVidOnset", "mockOffset", "displayVidOffset",
                             "displayVidDuration", "fixationPostVidDuration",
                             "displayAnswerDuration", "fixationPostAnswerDuration", "displayCuriosityDuration", "fixationPostCuriosityDuration")]
      duration$displayMockVidDuration <- duration$mockOffset - duration$displayVidOnset
      duration$displayVidDuration <- duration$displayVidOffset - duration$mockOffset
      duration <- duration[,c("trial", "vidFileName", "avgVidDur_MAGMOT", "displayMockVidDuration", "displayVidDuration", "fixationPostVidDuration",
                              "displayAnswerDuration", "fixationPostAnswerDuration", "displayCuriosityDuration", "fixationPostCuriosityDuration")]
      
      duration <- reshape2::melt(duration, id.vars=c("trial","vidFileName","avgVidDur_MAGMOT"), value.name = "duration")
      # recode variable (for later merging)
      # Levels: displayMockVidDuration displayVidDuration fixationPostVidDuration displayAnswerDuration fixationPostAnswerDuration displayCuriosityDuration fixationPostCuriosityDuration
      levels(duration$variable) <-  c("displayMockVid", "displayVid", "fixationPostVid",
                                      "displayAnswer", "fixationPostAnswer", "displayCuriosity", "fixationPostCuriosity")
      
      # merge onset and duration
      run_BIDS <- merge(onset, duration, by = c("trial","vidFileName", "variable"))
      
      
      # add ppt estimate on how many people are able to find the solution to the magic trick and the response time
      run_BIDS$response <- ifelse(run_BIDS$variable == "displayAnswer", as.character(run_BIDS$responseAnswer),
                                  ifelse(run_BIDS$variable == "displayCuriosity", run_BIDS$responseCuriosity, "not applicable"))
      run_BIDS$response_timestamp <- ifelse(run_BIDS$variable == "displayAnswer", run_BIDS$timestampAnswer,
                                            ifelse(run_BIDS$variable == "displayCuriosity", run_BIDS$timestampCuriosity, "not applicable"))
      
      # code memory performance and trial type (remembered vs forgotten)
      run_BIDS$allConfRecognition <- ifelse(run_BIDS$variable == "displayVid", run_BIDS$recognition, "not applicable")
      run_BIDS$highConfRecognition <- ifelse(run_BIDS$variable == "displayVid", run_BIDS$recognitionConfLevel_4_5_6, "not applicable")
      run_BIDS$aboveAverageConfRecognition <- ifelse(run_BIDS$variable == "displayVid", run_BIDS$recognitionAboveMeanConf, "not applicable")
      run_BIDS$confidence <- ifelse(run_BIDS$variable == "displayVid", run_BIDS$confidence, "not applicable")
      run_BIDS$cuedRecallStrict <- ifelse(run_BIDS$variable == "displayVid", run_BIDS$cuedRecallStrict, "not applicable")
      run_BIDS$cuedRecallLenient <- ifelse(run_BIDS$variable == "displayVid", run_BIDS$cuedRecallLenient, "not applicable")
      # run_BIDS$trial_type <-  ifelse(run_BIDS$highConfRecognition == 1 | run_BIDS$cuedRecall == 1, "remembered",
      #                                ifelse(run_BIDS$highConfRecognition == 0 & run_BIDS$cuedRecall == 0, "forgotten",
      #                                       ifelse(run_BIDS$cuedRecall ==  "n/a" & run_BIDS$highConfRecognition == 0, "n/a",
      #                                              ifelse(run_BIDS$highConfRecognition == "not applicable" & run_BIDS$cuedRecall == "not applicable", "not applicable",
      #                                                     "undefined"))))
      
      # note: this has to change after machine learning, currently trial type is only based on recognition performance
      run_BIDS$trial_type_cuedRecallStrict <-  ifelse(run_BIDS$cuedRecallStrict == 1, "remembered", 
                                                       ifelse(run_BIDS$cuedRecallStrict == 0, "forgotten",
                                                              ifelse(run_BIDS$cuedRecallStrict == "not applicable", "not applicable",
                                                                     "undefined")))
      
      run_BIDS$trial_type_cuedRecallLenient <-  ifelse(run_BIDS$cuedRecallLenient == 1, "remembered", 
                                                     ifelse(run_BIDS$cuedRecallLenient == 0, "forgotten",
                                                            ifelse(run_BIDS$cuedRecallLenient == "not applicable", "not applicable",
                                                                   "undefined")))

      run_BIDS$trial_type_allConfRecognition <-  ifelse(run_BIDS$allConfRecognition == 1, "remembered", 
                                                       ifelse(run_BIDS$allConfRecognition == 0, "forgotten",
                                                              ifelse(run_BIDS$allConfRecognition == "not applicable", "not applicable",
                                                                     "undefined")))
      
      # note: this has to change after machine learning, currently trial type is only based on recognition performance
      run_BIDS$trial_type_highConfRecognition <-  ifelse(run_BIDS$highConfRecognition == 1, "remembered",
                                                         ifelse(run_BIDS$highConfRecognition == 0, "forgotten",
                                                                ifelse(run_BIDS$highConfRecognition == "not applicable", "not applicable",
                                                                       "undefined")))
      
      # note: this has to change after machine learning, currently trial type is only based on recognition performance
      run_BIDS$trial_type_aboveAverageConfRecognition <-  ifelse(run_BIDS$aboveAverageConfRecognition == 1, "remembered",
                                                                 ifelse(run_BIDS$aboveAverageConfRecognition == 0, "forgotten",
                                                                        ifelse(run_BIDS$aboveAverageConfRecognition == "not applicable", "not applicable",
                                                                               "undefined")))
      
      # add not applicable whenever necessary
      run_BIDS$timestampVidOnset <- ifelse(run_BIDS$variable == "displayVid", run_BIDS$timestampVidOnset, "not applicable")
      run_BIDS$timestampMockOffset <- ifelse(run_BIDS$variable == "displayVid", run_BIDS$timestampMockOffset, "not applicable")
      run_BIDS$displayVidOffset <- ifelse(run_BIDS$variable == "displayVid", run_BIDS$displayVidOffset, "not applicable")
      run_BIDS$timestampPostVidFixation <- ifelse(run_BIDS$variable == "displayVid", run_BIDS$timestampPostVidFixation, "not applicable")
      run_BIDS$avgVidDur_MAGMOT <- ifelse(run_BIDS$variable == "displayVid", run_BIDS$avgVidDur_MAGMOT, "not applicable")
      
      run_BIDS$responseAnswer <- ifelse(run_BIDS$variable == "displayAnswer", run_BIDS$responseAnswer, "not applicable")
      run_BIDS$timestampAnswer <- ifelse(run_BIDS$variable == "displayAnswer", run_BIDS$timestampAnswer, "not applicable")
      
      run_BIDS$responseCuriosity <- ifelse(run_BIDS$variable == "displayCuriosity", run_BIDS$responseCuriosity, "not applicable")
      run_BIDS$timestampCuriosity <- ifelse(run_BIDS$variable == "displayCuriosity", run_BIDS$timestampCuriosity, "not applicable")
      
      # change col names of data frame according to BIDS specification
      names(run_BIDS)[names(run_BIDS)=="vidFileName"] <- "stim_file"
      names(run_BIDS)[names(run_BIDS)=="variable"] <- "event"
      run_BIDS <-  run_BIDS[order(run_BIDS$onset),]
      run_BIDS$onset <- round(run_BIDS$onset, digits = 3)
      run_BIDS$duration <- round(run_BIDS$duration, digits = 3)
      
      events_BIDS <- run_BIDS[, c("onset", "duration", "trial", "stim_file",  "event", "response", "response_timestamp",  
                                  "trial_type_cuedRecallStrict", "trial_type_cuedRecallLenient",
                                  "trial_type_allConfRecognition", "trial_type_highConfRecognition", "trial_type_aboveAverageConfRecognition")]
      
      # save file for BIDS
      if (max(run$acq) > 1) { # include acq in filename if there was more than one acq in a run
        BIDSfilename <- paste0(BIDSstring, "_task-magictrickwatching_acq-", a, "_run-", b,"_events.tsv")
      } else {
        BIDSfilename <- paste0(BIDSstring, "_task-magictrickwatching_run-", b,"_events.tsv")
      }
      
      setwd(preprocessedEventsSubjDir)
      write.table(events_BIDS, file=BIDSfilename, quote=FALSE, sep="\t", row.names = FALSE, na = "n/a")
      write.table(events_BIDS, file=file.path(preprocessedconcatDir, BIDSfilename), quote=FALSE, sep="\t", row.names = FALSE, na = "n/a")
      
      
      # if connected to VM, save it there too
      if (dir.exists(dirVM)==T & overwrite == "yes"){
        setwd(file.path(BIDSdirVM, BIDSstring, 'func'))
        write.table(events_BIDS, file=BIDSfilename, quote=FALSE, sep="\t", row.names = FALSE, na = "n/a")
        setwd(preprocessedEventsSubjDir)
      } else {
        if (feedback == "yes"){
          print(paste("Trying to save",  BIDSfilename, ", but not connected to study drive!!!"))
        }
      }
      
      # crate file for concatenation of BOLD data: mockOffset and onsets have to be translated
      run_BIDS$run <- b
      run_BIDS$acq <- a
      names(run_BIDS)[names(run_BIDS)=="timestampVidOnset"] <- "vid_per_run"
      names(run_BIDS)[names(run_BIDS)=="timestampMockOffset"] <- "mock_per_run"
      names(run_BIDS)[names(run_BIDS)=="onset"] <- "onset_per_run"
      
      if(b == 1 && a == 1){
        run_BIDS$onset <-  run_BIDS$onset_per_run
        run_BIDS$vid <-  as.numeric(run_BIDS$vid_per_run)
        run_BIDS$mock <-  as.numeric(run_BIDS$mock_per_run)
        BIDS <- run_BIDS
      } else {
        run_BIDS$onset <- run_BIDS$onset_per_run + totalDurationScan
        run_BIDS$vid <- as.numeric(run_BIDS$vid_per_run) + totalDurationScan
        run_BIDS$mock <- as.numeric(run_BIDS$mock_per_run) + totalDurationScan
        tempBIDS <-  run_BIDS
        BIDS <- rbind.all.columns(BIDS, tempBIDS) #rbind all columns will induce NA if there was initially no data saved in the loop per participant
        rm(tempBIDS)
      }
      
      # define the end of the last event
      whichLast <- which.max(run_BIDS$onset_per_run)
      lastOnset <- run_BIDS[[ whichLast, "onset"]]
      lastDuration <- run_BIDS[[ whichLast, "duration"]]
      
      durationRun <- run_BIDS[[ whichLast, "onset_per_run"]] + run_BIDS[[ whichLast, "duration"]]
      
      totalDurationScan <- lastOnset + lastDuration
      if (feedback == "yes"){
        print(paste0(BIDSfilename, "-vidOnly.tsv has a total duration of ", totalDurationScan, "s (", totalDurationScan/2, " TR)"))
      }
      scaninfo <- data.frame(BIDSstring, gsub("_events.tsv", "", BIDSfilename), durationRun, totalDurationScan)
      
      if (s == 1 && b == 1){ # for the first run of the first subject
        scaninfoAll <-  scaninfo
      } else {
        temp_scaninfo <-  scaninfo
        scaninfoAll <- rbind.all.columns(scaninfoAll, temp_scaninfo) #rbind all columns will induce NA if there was initially no data saved in the loop per participant
        rm(temp_scaninfo)
      }
      
      if (b == 3 && a == max(run$acq)){ # after the last run is finished
        if (debug == 0){
          rm(scaninfo, totalDurationScan, whichLast, lastOnset, lastDuration)
        }
      }
    }
    
    # save a file containg onsets, durations, stimID, and events with stimID and events as numbers to use it for the concat script in Matlab
    # save a file only containing onsets and durations of the actual video
    BIDS[(BIDS$event == "displayVid"), "responseCuriosity"] <-     BIDS[(BIDS$event == "displayCuriosity"), "responseCuriosity"] # overwrite the "not applicable" with the actual curiosity rating
    BIDS_concat <- subset(BIDS, BIDS$event == "displayVid") # BIDS contains all onsets and durations for all 36 trials as well as the memory performance associated with all of them
    BIDS_concat$duration_mock <- BIDS_concat$mock-BIDS_concat$vid
    BIDS_concat$duration_vid <- as.numeric(BIDS_concat$displayVidOffset)-as.numeric(BIDS_concat$vid_per_run)
    BIDS_concat$duration_vid_withoutMock <-as.numeric(BIDS_concat$displayVidOffset)-as.numeric(BIDS_concat$mock_per_run)
    BIDS_concat$duration_vid_postFixation <-as.numeric(BIDS_concat$timestampPostVidFixation)-as.numeric(BIDS_concat$vid_per_run)
    BIDS_concat$duration_vid_withoutMock_postFixation <-as.numeric(BIDS_concat$timestampPostVidFixation)-as.numeric(BIDS_concat$mock_per_run)
    
    #names(BIDS_concat)
    BIDS_concat <- BIDS_concat[,c("vid", "mock", "duration_vid", "duration_vid_withoutMock", "duration_vid_withoutMock_postFixation", "avgVidDur_MAGMOT", "trial", "stim_file", "responseCuriosity", 
                                  "trial_type_cuedRecallStrict", "trial_type_cuedRecallLenient", "trial_type_allConfRecognition", 
                                  "trial_type_highConfRecognition", "trial_type_aboveAverageConfRecognition",  "confidence",  "run", "acq")]
    
    if (feedback == "yes"){
      print(paste("total Duration is", sum(BIDS_concat$duration)))
    }
    write.table(BIDS_concat, file=paste0(BIDSstring, "_task-magictrickwatching_concat.tsv"), quote=FALSE, sep="\t", row.names = FALSE, na = "n/a")
    write.table(BIDS_concat, file=file.path(preprocessedconcatDir, paste0(BIDSstring, "_task-magictrickwatching_concat.tsv")), quote=FALSE, sep="\t", row.names = FALSE, na = "n/a")
    
    
    
    
    # if connected to VM, save it there too
    if (dir.exists(dirVM)==T & overwrite == "yes"){
      preprocessedEventsSubjDirVM <- file.path(preprocessedEventsRootDirVM, "concat", BIDSstring)
      ifelse(!dir.exists(preprocessedEventsSubjDirVM), dir.create(preprocessedEventsSubjDirVM), FALSE)
      setwd(preprocessedEventsSubjDirVM)
      write.table(BIDS_concat, file=paste0(BIDSstring, "_task-magictrickwatching_concat.tsv"), quote=FALSE, sep="\t", row.names = FALSE, na = "n/a")
      setwd(preprocessedEventsSubjDir)
    } else {
      if (feedback == "yes" & overwrite == "yes"){
        print("Trying to save concat.tsv, but not connected to study drive!!!")
      }
    }
  }
  
  # delete no longer needed variables
  if (debug == 0){
    rm(duration, onset, run, run_BIDS, BIDS, events_BIDS)
  }
  
  
  
  
  #################################################### process questionnaire data collected during memory part
  postMemory <-subset(memory, trial.type == "surveycat")   # select relevant rows and columns of memory task questions data
  postMemory <- postMemory[,c("username","startMemory", "endMemory", "memoryFile",
                              "survey_sleep_response","survey_sleep_hours",
                              "survey_test_known_response", "survey_memory_intention_response", "survey_reward_belief_response",
                              "survey_magictrick_experience_response", "survey_connection_response", "survey_comment_response")]
  names(postMemory) <- c("ID", "startMemory", "endMemory", "memoryFile",
                         "sleepBeforeMemoryTest","sleepHours", "memoryTestKnown", "memoryIntention", "rewardBelief", "magictrickExperience", "connection", "comment")
  
  # add curiosity NAs to the data set in wide format
  postMemory$curiosityNAs <- curiosityNAs
  
  
  ########### add columns to postMemory data
  postMemory$BIDS <- BIDSstring
  if (file.exists(file.path(dataMemoryDir,f))) { # check whether there is  data from the memory test at all; if so compute sum scores for recognition performance

    # subset the data depending on block
    for (BLOCK in 1:(max(MEMO$block)+1)) {
      data_subset <- subset(MEMO, MEMO$block == BLOCK)
      if(BLOCK == 4){ # as well as for the data set in total
        data_subset <- MEMO
      }
      if (feedback == "yes"){
        print(paste("rows for each of the blocks:",nrow(data_subset)))
      }
      
      # sum up the scores for recall task
      if (file.exists(file.path(codedDir,f_coded))) {  # if there is no data, NA will be added when rbinding all information across subjects
        postMemory[[paste0("cuedRecallStrict", blockstring[BLOCK])]] <- sum(MEMO$cuedRecallStrict, na.rm = T) #please note that this needs to be changed as it is not looking at any form of coded data
        postMemory[[paste0("cuedRecallLenient", blockstring[BLOCK])]] <- sum(MEMO$cuedRecallLenient, na.rm = T) #please note that this needs to be changed as it is not looking at any form of coded data
      }
      
      
      
      # sum up the scores for the recognition task, once in total and once seperated for the different levels of confidence
      postMemory[[paste0("responseCuriosity", blockstring[BLOCK])]] <- mean(data_subset$responseCuriosity, na.rm = T)
      postMemory[[paste0("curiosity", blockstring[BLOCK])]] <- mean(data_subset$curiosity, na.rm = T)
      postMemory[[paste0("recognition", blockstring[BLOCK])]] <- sum(data_subset$recognition, na.rm = T)
      postMemory[[paste0("recognitionAboveMeanConf", blockstring[BLOCK])]] <- sum(data_subset$recognitionAboveMeanConf, na.rm = T)
      postMemory[[paste0("recognitionContConf", blockstring[BLOCK])]] <- sum(data_subset$recognitionContConf, na.rm = T)
      
      postMemory[[paste0("meanConfidence", blockstring[BLOCK])]]  <- mean(data_subset$confidence, na.rm = T)
      temp_data <- subset(data_subset, data_subset$recognition == 1)
      postMemory[[paste0("meanConfidenceCorrectTrials", blockstring[BLOCK])]]   <- mean(temp_data$confidence, na.rm = T)
      
      for (k in 1:6) { #confidence ranges from 1 to 6, potentially code can be made more flexible by using min(data$confidence) and max(data$confidence)
        temp_data <- subset(data_subset, data_subset$confidence == k)
        postMemory[[paste0("recognitionConfLevel_", k, blockstring[BLOCK])]] <- sum(temp_data$recognition, na.rm = T)
        if (k < 6) {
          temp_data_above <- subset(data_subset, data_subset$confidence > k)
          postMemory[[paste0("recognitionConfLevel_above_", k, blockstring[BLOCK])]] <-sum(temp_data_above$recognition, na.rm = T) 
          rm(temp_data_above)
        }
        
        # sum up the scores for the recognition task pooled
        if (k == 1 || k == 3 || k == 5) {
          temp_data_plus1 <- subset(data_subset, data_subset$confidence == (k+1))
          postMemory[[paste0("recognitionConfLevel_", k, "_", k+1, blockstring[BLOCK])]] <- sum(temp_data$recognition, na.rm = T) + sum(temp_data_plus1$recognition, na.rm = T)
          rm(temp_data_plus1)
        }
        if (k == 1 || k == 4) {
          temp_data_plus1 <- subset(data_subset, data_subset$confidence == (k+1))
          temp_data_plus2 <- subset(data_subset, data_subset$confidence == (k+2))
          postMemory[[paste0("recognitionConfLevel_", k, "_", k+1, "_", k+2, blockstring[BLOCK])]] <- sum(temp_data$recognition, na.rm = T) + sum(temp_data_plus1$recognition, na.rm = T) + sum(temp_data_plus2$recognition, na.rm = T)
          rm(temp_data_plus1)
          rm(temp_data_plus2)
        }
        rm(temp_data)
      }
    }
    if (debug == 0){
      rm(data_subset)
    }
  }
  
  # sum up curiosity-driven memory memory benefit (continouos) for subjects
  postMemory$curiosityBenefit_cuedRecallStrict <-  sum(MEMO$curiosityBenefit_cuedRecallStrict, na.rm = T)
  postMemory$curiosityBenefit_cuedRecallLenient <-  sum(MEMO$curiosityBenefit_cuedRecallLenient, na.rm = T)
  postMemory$curiosityBenefit_allConf <-  sum(MEMO$curiosityBenefit_allConf, na.rm = T)
  postMemory$curiosityBenefit_highConf <-  sum(MEMO$curiosityBenefit_highConf, na.rm = T)
  postMemory$curiosityBenefit_aboveAvgConf <- sum(MEMO$curiosityBenefit_aboveAvgConf, na.rm = T)
  postMemory$curiosityBenefit_rememberedStrict <- sum(MEMO$curiosityBenefit_rememberedStrict, na.rm = T)
  postMemory$curiosityBenefit_rememberedLenient <- sum(MEMO$curiosityBenefit_rememberedLenient, na.rm = T)
  
  # sum up curiosity-driven memory memory benefit (dichotomous) for subjects
  postMemory$curiosityBenefit_cuedRecallStrict_dichotom <-  sum(MEMO$curiosityBenefit_cuedRecallStrict_dichotom, na.rm = T)
  postMemory$curiosityBenefit_cuedRecallLenient_dichotom <-  sum(MEMO$curiosityBenefit_cuedRecallLenient_dichotom, na.rm = T)
  postMemory$curiosityBenefit_allConf_dichotom <-  sum(MEMO$curiosityBenefit_allConf_dichotom, na.rm = T)
  postMemory$curiosityBenefit_highConf_dichotom <-  sum(MEMO$curiosityBenefit_highConf_dichotom, na.rm = T)
  postMemory$curiosityBenefit_aboveAvgConf_dichotom <-  sum(MEMO$curiosityBenefit_aboveAvgConf_dichotom, na.rm = T)
  postMemory$curiosityBenefit_rememberedStrict_dichotom <- sum(MEMO$curiosityBenefit_rememberedStrict_dichotom, na.rm = T)
  postMemory$curiosityBenefit_rememberedLenient_dichotom <- sum(MEMO$curiosityBenefit_rememberedLenient_dichotom, na.rm = T)
  
  # rbind the postMemory files of each subject to a data frame
  if(s == 1){
    postMemoryWide <- postMemory
  } else {
    temp_postMemoryWide <-  postMemory
    postMemoryWide <- rbind.all.columns(postMemoryWide, temp_postMemoryWide) #rbind all columns will induce NA if there was initially no data saved in the loop per participant
    rm(temp_postMemoryWide)
  }
  if (debug == 0){
    rm(postMemory, memory)
  }
  
  
  #################################################### at the end of the loop, merge data sets ####################################################
  if (s == length(subjects)){
    MAGMOT <- merge(MAGMOT, questDataWide, by = c("ID", "group"), all = T)
    MAGMOT <- merge(MAGMOT, postMemoryWide, by = c("ID"), all = T)
    if (debug == 0){
      rm(postMemoryWide, questDataWide)
    }
    
    # compute summary statistics for measurements of memory
    workspace <- list.files(path = file.path(codedDir), pattern = "_CP.csv") # check whether the data is coded yet or not
    if(length(workspace) == 0) { # if data is not coded yet, only look at recognition performance
      recognitionPerformanceVars <- c("recognition", "recognitionConfLevel_1", "recognitionConfLevel_2", "recognitionConfLevel_3", "recognitionConfLevel_4", "recognitionConfLevel_5", "recognitionConfLevel_6",
                                      "recognitionConfLevel_1_2", "recognitionConfLevel_3_4", "recognitionConfLevel_5_6", "recognitionConfLevel_1_2_3", "recognitionConfLevel_4_5_6")
    } else {
      recognitionPerformanceVars <- c("cuedRecallLenient", "cuedRecallStrict",
                                      "recognition", "recognitionConfLevel_1", "recognitionConfLevel_2", "recognitionConfLevel_3", "recognitionConfLevel_4", "recognitionConfLevel_5", "recognitionConfLevel_6",
                                      "recognitionConfLevel_1_2", "recognitionConfLevel_3_4", "recognitionConfLevel_5_6", "recognitionConfLevel_1_2_3", "recognitionConfLevel_4_5_6")
    }
    
    # create a data frame that shows the summary statistics for each score, for the whole population
    dataWideRecognitionPerformance <- MAGMOT[,recognitionPerformanceVars]
    mean <- data.frame(apply(dataWideRecognitionPerformance, 2, mean,  na.rm = T), recognitionPerformanceVars)
    sd <- data.frame(apply(dataWideRecognitionPerformance, 2, sd,  na.rm = T), recognitionPerformanceVars)
    min <- data.frame(apply(dataWideRecognitionPerformance, 2, min,  na.rm = T), recognitionPerformanceVars)
    max <- data.frame(apply(dataWideRecognitionPerformance, 2, max,  na.rm = T), recognitionPerformanceVars)
    
    descriptivesRecognitionPerformance <- merge(mean, sd, by = "recognitionPerformanceVars")
    descriptivesRecognitionPerformance <- merge(descriptivesRecognitionPerformance, min, by = "recognitionPerformanceVars")
    descriptivesRecognitionPerformance <- merge(descriptivesRecognitionPerformance, max, by = "recognitionPerformanceVars")
    
    names(descriptivesRecognitionPerformance) <- c("recognitionPerformanceVar", "mean", "sd", "min", "max")
    
    setwd(preprocessedDir)
    xlsx::write.xlsx(descriptivesRecognitionPerformance, file= paste0("memoryPerformance_MAGMOT_N", length(subjects), "_", format(Sys.time(), "%Y-%m-%d"), ".xlsx"), sheetName = "Sheet1", row.names = F)
    
    # write information about scan durations
    names(scaninfoAll) <- c("ID", "scan", "duration_run_seconds", "duration_scan_seconds")
    scaninfoAll$duration_run_seconds <- round(scaninfoAll$duration_run_seconds, digits = 0)
    scaninfoAll$duration_scan_seconds <- round(scaninfoAll$duration_scan_seconds, digits = 0)
    scaninfoAll$duration_run_TR <- scaninfoAll$duration_run_seconds/TR
    scaninfoAll$duration_scan_TR <- scaninfoAll$duration_scan_seconds/TR
    setwd(preprocessedEventsRootDir)
    write.table(scaninfoAll, file="MAGMOT_informationAboutScanDuration.tsv", quote=FALSE, sep="\t", row.names = FALSE, na = "n/a")
    
    if (dir.exists(dirVM)==T & overwrite == "yes"){
      setwd(file.path(preprocessedEventsRootDirVM, 'scripts'))
      write.table(scaninfoAll, file="MAGMOT_informationAboutScanDuration.tsv", quote=FALSE, sep="\t", row.names = FALSE, na = "n/a")
      setwd(preprocessedEventsSubjDir)
    } else {
      if (feedback == "yes" & overwrite == "yes"){
        print("Trying to save information about scan duration, but not connected to study drive!!!")
      }
    }
    
    
    #### combine single _long.csv files to one file ####
    setwd(file.path(preprocessedDir, "long"))
    file_list <- list.files(pattern = "long.csv")
    
    for (n in seq_along(file_list)){
      if(n == 1){
        dataLong <- read.csv(file_list[n], header=T)
      }
      else {
        temp_datalang <- read.csv(file_list[n], header=T)
        dataLong <- rbind.all.columns(dataLong, temp_datalang) #rbind all columns will induce NA if there was initially no data saved in the loop per participant
        rm(temp_datalang)
      }
    }
    setwd(preprocessedDir)
    xlsx::write.xlsx(dataLong, file="long_MAGMOT.xlsx", sheetName = "Sheet1", row.names = F) 
    
    # compute the glmer models to extract the slopes
    library(lme4)
    LMEmodel_cuedRecallStrict <- glmer(cuedRecallStrict ~ groupEffectCoded*curiosityGroupMeanCentered + (1+curiosityGroupMeanCentered|ID) + (1|stimID) , family = "binomial"(link = 'logit'), data = dataLong)
    MAGMOT$curiosityBeta_cuedRecallStrict <- coef(LMEmodel_cuedRecallStrict)$ID$curiosityGroupMeanCentered
    MAGMOT$curiosityBeta_cuedRecallStrict_c <-  MAGMOT$curiosityBeta_cuedRecallStrict - mean(MAGMOT$curiosityBeta_cuedRecallStrict)   
    
    LMEmodel_cuedRecallLenient <- glmer(cuedRecallLenient ~ groupEffectCoded*curiosityGroupMeanCentered + (1+curiosityGroupMeanCentered|ID) + (1|stimID) , family = "binomial"(link = 'logit'), data = dataLong)
    MAGMOT$curiosityBeta_cuedRecallLenient <- coef(LMEmodel_cuedRecallLenient)$ID$curiosityGroupMeanCentered
    MAGMOT$curiosityBeta_cuedRecallLenient_c <-  MAGMOT$curiosityBeta_cuedRecallLenient - mean(MAGMOT$curiosityBeta_cuedRecallLenient)   
    
    LMEmodel_recogAboveAvgConf <- glmer(recognitionAboveMeanConf ~ groupEffectCoded*curiosityGroupMeanCentered + (1+curiosityGroupMeanCentered|ID) + (1|stimID) , family = "binomial"(link = 'logit'), data = dataLong)
    MAGMOT$curiosityBeta_aboveAvgConf <- coef(LMEmodel_recogAboveAvgConf)$ID$curiosityGroupMeanCentered
    MAGMOT$curiosityBeta_aboveAvgConf_c <-  MAGMOT$curiosityBeta_aboveAvgConf - mean(MAGMOT$curiosityBeta_aboveAvgConf)
    
    LMEmodel_recogHighConf <- lme4::glmer(recognitionConfLevel_4_5_6 ~ groupEffectCoded*curiosityGroupMeanCentered + (1+curiosityGroupMeanCentered|ID) + (1|stimID) , family = "binomial"(link = 'logit'), data = dataLong)
    MAGMOT$curiosityBeta_highConf <- coef(LMEmodel_recogHighConf)$ID$curiosityGroupMeanCentered
    MAGMOT$curiosityBeta_highConf_c <-  MAGMOT$curiosityBeta_highConf - mean(MAGMOT$curiosityBeta_highConf)
    
    LMEmodel_recogAllConf <- glmer(recognition ~ groupEffectCoded*curiosityGroupMeanCentered + (1+curiosityGroupMeanCentered|ID) + (1|stimID) , family = "binomial"(link = 'logit'), data = dataLong)
    MAGMOT$curiosityBeta_allConf <- coef(LMEmodel_recogAllConf)$ID$curiosityGroupMeanCentered
    MAGMOT$curiosityBeta_allConf_c <-  MAGMOT$curiosityBeta_allConf - mean(MAGMOT$curiosityBeta_allConf)
    
    LMEmodel_rememberedStrict <- glmer(rememberedStrict ~ groupEffectCoded*curiosityGroupMeanCentered + (1+curiosityGroupMeanCentered|ID) + (1|stimID) , family = "binomial"(link = 'logit'), data = dataLong)
    print(summary(LMEmodel_rememberedStrict))
    MAGMOT$curiosityBeta_rememberedStrict <- coef(LMEmodel_rememberedStrict)$ID$curiosityGroupMeanCentered
    MAGMOT$curiosityBeta_rememberedStrict_c <-  MAGMOT$curiosityBeta_rememberedStrict - mean(MAGMOT$curiosityBeta_rememberedStrict)
    
    LMEmodel_rememberedLenient <- glmer(rememberedLenient ~ groupEffectCoded*curiosityGroupMeanCentered + (1+curiosityGroupMeanCentered|ID) + (1|stimID) , family = "binomial"(link = 'logit'), data = dataLong)
    print(summary(LMEmodel_rememberedLenient))
    MAGMOT$curiosityBeta_rememberedLenient <- coef(LMEmodel_rememberedLenient)$ID$curiosityGroupMeanCentered
    MAGMOT$curiosityBeta_rememberedLenient_c <-  MAGMOT$curiosityBeta_rememberedLenient - mean(MAGMOT$curiosityBeta_rememberedLenient)
    
    
    # add RSFC estimates between HPC & VTA (Pearson)
    setwd(brainDir)
    RSFC <- read.delim2("RSFC_VTA-HPC.txt") # read in data
    names(RSFC)[names(RSFC)=="RSFC_run.1"] <- "RSFC_VTAHPC_run1" # change col name
    names(RSFC)[names(RSFC)=="RSFC_run.2"] <- "RSFC_VTAHPC_run2" # change col name
    RSFC$RSFC_VTAHPC_run1 <- as.numeric(as.character(RSFC$RSFC_VTAHPC_run1)) # change from factor to numeric
    RSFC$RSFC_VTAHPC_run2 <- as.numeric(as.character(RSFC$RSFC_VTAHPC_run2)) # change from factor to numeric
    
    # fisher z transform correlations and compute RSFC change
    RSFC$RSFC_VTAHPC_run1_z <- atanh(RSFC$RSFC_VTAHPC_run1) # fisher z transform correlations
    RSFC$RSFC_VTAHPC_run2_z <- atanh(RSFC$RSFC_VTAHPC_run2) # fisher z transform correlations
    RSFC$RSFC_VTAHPC_diff <- RSFC$RSFC_VTAHPC_run2_z - RSFC$RSFC_VTAHPC_run1_z
    
    MAGMOT <- merge(MAGMOT, RSFC, by = "BIDS")
    
    # add RSFC estimates between HPC & VTA (Spearman)
    RSFC <- read.delim2("RSFC_VTA-HPC_spearman.txt") # read in data
    names(RSFC)[names(RSFC)=="RSFC_run.1_spearman"] <- "RSFC_VTAHPC_run1_spearman" # change col name
    names(RSFC)[names(RSFC)=="RSFC_run.2_spearman"] <- "RSFC_VTAHPC_run2_spearman" # change col name
    RSFC$RSFC_VTAHPC_run1_spearman <- as.numeric(as.character(RSFC$RSFC_VTAHPC_run1_spearman)) # change from factor to numeric
    RSFC$RSFC_VTAHPC_run2_spearman <- as.numeric(as.character(RSFC$RSFC_VTAHPC_run2_spearman)) # change from factor to numeric
    
    # fisher z transform correlations and compute RSFC change
    RSFC$RSFC_VTAHPC_run1_z_spearman <- atanh(RSFC$RSFC_VTAHPC_run1_spearman) # fisher z transform correlations
    RSFC$RSFC_VTAHPC_run2_z_spearman <- atanh(RSFC$RSFC_VTAHPC_run2_spearman) # fisher z transform correlations
    RSFC$RSFC_VTAHPC_diff_spearman <- RSFC$RSFC_VTAHPC_run2_z_spearman - RSFC$RSFC_VTAHPC_run1_z_spearman
    
    MAGMOT <- merge(MAGMOT, RSFC, by = "BIDS")
    
    
    # par(mfrow=c(1,3))
    # boxplot(MAGMOT$curiosityBeta_aboveAvgConf, main = "curiosityBeta_aboveAvgConf")
    # boxplot(MAGMOT$curiosityBeta_highConf, main = "curiosityBeta_highConf")
    # boxplot(MAGMOT$curiosityBeta_allConf, main = "curiosityBeta_allConf")
    
    setwd(preprocessedDir)
    xlsx::write.xlsx(MAGMOT, file="wide_MAGMOT.xlsx", sheetName = "Sheet1", row.names = F)
    
    
    #################### as a last step, create the files we need for concatenation ####################
    
    if (doISCprep == 1){
      
      # aim: create:
      # "sub-control001sub-control002_sub-control001_bothRemembered_avgConf_SME_concat.tsv"
      # "sub-control001sub-control002_sub-control002_bothRemembered_avgConf_SME_concat.tsv"
      # "sub-control001sub-control002_sub-control001_bothForgotten_avgConf_SME_concat.tsv"
      # "sub-control001sub-control002_sub-control002_bothForgotten_avgConf_SME_concat.tsv"
      # "sub-control001sub-control002_sub-control001_differentResponses_avgConf_SME_concat.tsv"
      # "sub-control001sub-control002_sub-control002_differentResponses_avgConf_SME_concat.tsv"
      
      # "sub-control001sub-control002_sub-control001_bothRemembered_highConf_SME_concat.tsv"
      # "sub-control001sub-control002_sub-control002_bothRemembered_highConf_SME_concat.tsv"
      # "sub-control001sub-control002_sub-control001_bothForgotten_highConf_SME_concat.tsv"
      # "sub-control001sub-control002_sub-control002_bothForgotten_highConf_SME_concat.tsv"
      # "sub-control001sub-control002_sub-control001_differentResponses_highConf_SME_concat.tsv"
      # "sub-control001sub-control002_sub-control002_differentResponses_highConf_SME_concat.tsv"
      
      # get a list with all subject BIDS strings
      subjectsCorr <- levels(dataLong$BIDS)
      
      #noMemoryYet <- c("sub-experimental012")
      #subjectsCorr <- subjectsCorr[!(subjectsCorr %in% noMemoryYet)]
      subjectsToCorrelate <- subjectsCorr
      
      N <- 0.5*length(subjectsCorr)*length(subjectsToCorrelate)
      
      # create empty df for 3dISC -dataTable
      dataTable_allConf <- data.frame()
      dataTable_highConf <- data.frame()
      dataTable_aboveAvgConf <- data.frame()
      dataTable_ISC <- data.frame()
      dataTable_ISC_dummy <- data.frame()
      dataTable <- data.frame()
      x <- 0
      
      if ( identical(subjectsToCorrelate, character(0)) == F){
        
        for(s in seq_along(subjectsCorr)){
          
          #ss = s+1
          subjectsToCorrelate <-subjectsCorr[-c(1:s)]
          
          # read in concat for subjectsCorr[s]
          preprocessedEventsSubjDir_s <- file.path(preprocessedEventsRootDir, subjectsCorr[s])
          setwd(preprocessedEventsSubjDir_s)
          file_s <- list.files(pattern = "concat")
          
          events_s <- read.delim(file = file_s, header = T, sep="\t", na = "n/a")
          
          # rename vars
          events_s <- events_s[,c("vid", "mock", "duration_vid", "duration_vid_withoutMock", "avgVidDur_MAGMOT", "stim_file", "responseCuriosity", 
                                  "trial_type_cuedRecallStrict", "trial_type_cuedRecallLenient",
                                  "trial_type_allConfRecognition", "trial_type_highConfRecognition", "trial_type_aboveAverageConfRecognition", "confidence")]
          events_s$responseCuriosity <- events_s$responseCuriosity - mean(events_s$responseCuriosity, na.rm = T) #mean center curiosity
          
          names(events_s)[names(events_s)=="vid"] <- paste0("onset_", subjectsCorr[s])
          names(events_s)[names(events_s)=="duration_vid"] <- paste0("duration_vid_", subjectsCorr[s])
          names(events_s)[names(events_s)=="mock"] <- paste0("mock_", subjectsCorr[s])
          names(events_s)[names(events_s)=="duration_vid_withoutMock"] <- paste0("duration_vid_withoutMock_", subjectsCorr[s])
          names(events_s)[names(events_s)=="responseCuriosity"] <- paste0("responseCuriosity_", subjectsCorr[s])
          names(events_s)[names(events_s)=="trial_type_cuedRecallStrict"] <- paste0("trial_type_cuedRecallStrict_", subjectsCorr[s])
          names(events_s)[names(events_s)=="trial_type_cuedRecallLenient"] <- paste0("trial_type_cuedRecallLenient_", subjectsCorr[s])
          names(events_s)[names(events_s)=="trial_type_allConfRecognition"] <- paste0("trial_type_allConfRecognition_", subjectsCorr[s])
          names(events_s)[names(events_s)=="trial_type_highConfRecognition"] <- paste0("trial_type_highConfRecognition_", subjectsCorr[s])
          names(events_s)[names(events_s)=="trial_type_aboveAverageConfRecognition"] <- paste0("trial_type_aboveAverageConfRecognition_", subjectsCorr[s])
          events_s$confidence <- events_s$confidence - mean(events_s$confidence, na.rm = T) #mean center confidence
          names(events_s)[names(events_s)=="confidence"] <- paste0("confidence_", subjectsCorr[s])
          
          for(ss in seq_along(subjectsToCorrelate)){
            
            # read in concat for subjectsToCorrelate[ss]
            preprocessedEventsSubjDir_ss <- file.path(preprocessedEventsRootDir, paste(subjectsToCorrelate[ss]))
            setwd(preprocessedEventsSubjDir_ss)
            file_ss <- list.files(pattern = "concat")
            
            events_ss <- read.delim(file = file_ss, header = T, sep="\t", na = "n/a")
            
            # rename vars
            events_ss <- events_ss[,c("vid", "mock", "duration_vid", "duration_vid_withoutMock", "avgVidDur_MAGMOT", "stim_file", "responseCuriosity", 
                                      "trial_type_cuedRecallStrict", "trial_type_cuedRecallLenient",
                                      "trial_type_allConfRecognition", "trial_type_highConfRecognition", "trial_type_aboveAverageConfRecognition", "confidence")]
            events_ss$responseCuriosity <- events_ss$responseCuriosity - mean(events_ss$responseCuriosity, na.rm = T) #mean center curiosity
            
            names(events_ss)[names(events_ss)=="vid"] <- paste0("onset_", subjectsToCorrelate[ss])
            names(events_ss)[names(events_ss)=="duration_vid"] <- paste0("duration_vid_", subjectsToCorrelate[ss])
            names(events_ss)[names(events_ss)=="mock"] <- paste0("mock_", subjectsToCorrelate[ss])
            names(events_ss)[names(events_ss)=="duration_vid_withoutMock"] <- paste0("duration_vid_withoutMock_", subjectsToCorrelate[ss])
            names(events_ss)[names(events_ss)=="responseCuriosity"] <- paste0("responseCuriosity_", subjectsToCorrelate[ss])
            names(events_ss)[names(events_ss)=="trial_type_cuedRecallStrict"] <- paste0("trial_type_cuedRecallStrict_", subjectsToCorrelate[ss])
            names(events_ss)[names(events_ss)=="trial_type_cuedRecallLenient"] <- paste0("trial_type_cuedRecallLenient_", subjectsToCorrelate[ss])
            names(events_ss)[names(events_ss)=="trial_type_allConfRecognition"] <- paste0("trial_type_allConfRecognition_", subjectsToCorrelate[ss])
            names(events_ss)[names(events_ss)=="trial_type_highConfRecognition"] <- paste0("trial_type_highConfRecognition_", subjectsToCorrelate[ss])
            names(events_ss)[names(events_ss)=="trial_type_aboveAverageConfRecognition"] <- paste0("trial_type_aboveAverageConfRecognition_", subjectsToCorrelate[ss])
            events_ss$confidence <- events_ss$confidence - mean(events_ss$confidence, na.rm = T) #mean center confidence
            names(events_ss)[names(events_ss)=="confidence"] <- paste0("confidence_", subjectsToCorrelate[ss])
            
            # marge both subjects
            events <- merge(events_s, events_ss, by = c("stim_file", "avgVidDur_MAGMOT"))
            
            # define match in memory performance
            events$behavParcel_cuedRecallStrict <- ifelse(events[[paste0("trial_type_cuedRecallStrict_",subjectsCorr[s])]] == "remembered" & events[[paste0("trial_type_cuedRecallStrict_",subjectsToCorrelate[ss])]] == "remembered", "bothRemembered",
                                                 ifelse(events[[paste0("trial_type_cuedRecallStrict_",subjectsCorr[s])]] == "forgotten" & events[[paste0("trial_type_cuedRecallStrict_",subjectsToCorrelate[ss])]] == "forgotten", "bothForgotten",
                                                        "differentResponses"))
            events$behavParcel_cuedRecallLenient <- ifelse(events[[paste0("trial_type_cuedRecallLenient_",subjectsCorr[s])]] == "remembered" & events[[paste0("trial_type_cuedRecallLenient_",subjectsToCorrelate[ss])]] == "remembered", "bothRemembered",
                                                 ifelse(events[[paste0("trial_type_cuedRecallLenient_",subjectsCorr[s])]] == "forgotten" & events[[paste0("trial_type_cuedRecallLenient_",subjectsToCorrelate[ss])]] == "forgotten", "bothForgotten",
                                                        "differentResponses"))
            events$behavParcel_allConf <- ifelse(events[[paste0("trial_type_allConfRecognition_",subjectsCorr[s])]] == "remembered" & events[[paste0("trial_type_allConfRecognition_",subjectsToCorrelate[ss])]] == "remembered", "bothRemembered",
                                                 ifelse(events[[paste0("trial_type_allConfRecognition_",subjectsCorr[s])]] == "forgotten" & events[[paste0("trial_type_allConfRecognition_",subjectsToCorrelate[ss])]] == "forgotten", "bothForgotten",
                                                        "differentResponses"))

            events$behavParcel_highConf <- ifelse(events[[paste0("trial_type_highConfRecognition_",subjectsCorr[s])]] == "remembered" & events[[paste0("trial_type_highConfRecognition_",subjectsToCorrelate[ss])]] == "remembered", "bothRemembered",
                                                  ifelse(events[[paste0("trial_type_highConfRecognition_",subjectsCorr[s])]] == "forgotten" & events[[paste0("trial_type_highConfRecognition_",subjectsToCorrelate[ss])]] == "forgotten", "bothForgotten",
                                                         "differentResponses"))
            events$behavParcel_aboveAvgConf <- ifelse(events[[paste0("trial_type_aboveAverageConfRecognition_",subjectsCorr[s])]] == "remembered" & events[[paste0("trial_type_aboveAverageConfRecognition_",subjectsToCorrelate[ss])]] == "remembered", "bothRemembered",
                                                      ifelse(events[[paste0("trial_type_aboveAverageConfRecognition_",subjectsCorr[s])]] == "forgotten" & events[[paste0("trial_type_aboveAverageConfRecognition_",subjectsToCorrelate[ss])]] == "forgotten", "bothForgotten",
                                                             "differentResponses"))
            
            # create an effect coded memory variable
            events[[paste0("trial_type_cuedRecallStrict_",subjectsCorr[s], "_effectCoded")]] <- ifelse( events[[paste0("trial_type_cuedRecallStrict_",subjectsCorr[s])]] == "remembered", 1, -1)
            events[[paste0("trial_type_cuedRecallStrict_",subjectsToCorrelate[ss], "_effectCoded")]]  <- ifelse( events[[paste0("trial_type_cuedRecallStrict_",subjectsToCorrelate[ss])]] == "remembered", 1, -1)
            events[[paste0("trial_type_cuedRecallLenient_",subjectsCorr[s], "_effectCoded")]] <- ifelse( events[[paste0("trial_type_cuedRecallLenient_",subjectsCorr[s])]] == "remembered", 1, -1)
            events[[paste0("trial_type_cuedRecallLenient_",subjectsToCorrelate[ss], "_effectCoded")]]  <- ifelse( events[[paste0("trial_type_cuedRecallLenient_",subjectsToCorrelate[ss])]] == "remembered", 1, -1)
            
            events[[paste0("trial_type_allConfRecognition_",subjectsCorr[s], "_effectCoded")]] <- ifelse( events[[paste0("trial_type_allConfRecognition_",subjectsCorr[s])]] == "remembered", 1, -1)
            events[[paste0("trial_type_allConfRecognition_",subjectsToCorrelate[ss], "_effectCoded")]]  <- ifelse( events[[paste0("trial_type_allConfRecognition_",subjectsToCorrelate[ss])]] == "remembered", 1, -1)
            events[[paste0("trial_type_highConfRecognition_",subjectsCorr[s], "_effectCoded")]] <- ifelse( events[[paste0("trial_type_highConfRecognition_",subjectsCorr[s])]] == "remembered", 1, -1)
            events[[paste0("trial_type_highConfRecognition_",subjectsToCorrelate[ss], "_effectCoded")]]  <- ifelse( events[[paste0("trial_type_highConfRecognition_",subjectsToCorrelate[ss])]] == "remembered", 1, -1)
            events[[paste0("trial_type_aboveAverageConfRecognition_",subjectsCorr[s], "_effectCoded")]] <- ifelse( events[[paste0("trial_type_aboveAverageConfRecognition_",subjectsCorr[s])]] == "remembered", 1, -1)
            events[[paste0("trial_type_aboveAverageConfRecognition_",subjectsToCorrelate[ss], "_effectCoded")]]  <- ifelse( events[[paste0("trial_type_aboveAverageConfRecognition_",subjectsToCorrelate[ss])]] == "remembered", 1, -1)
            # multiply memory and confidence
            # events[[paste0("trial_type_allConfRecognition_",subjectsCorr[s], "_dummyCoded")]] <- ifelse( events[[paste0("trial_type_allConfRecognition_",subjectsCorr[s])]] == "remembered", 1, 0)
            # events[[paste0("trial_type_allConfRecognition_",subjectsToCorrelate[ss], "_dummyCoded")]]  <- ifelse( events[[paste0("trial_type_allConfRecognition_",subjectsToCorrelate[ss])]] == "remembered", 1, 0)
            # events[[paste0("trial_type_contConfRecognition_",subjectsCorr[s], "_dummyCoded")]] <- events[[paste0("trial_type_allConfRecognition_",subjectsCorr[s], "_dummyCoded")]] * events[[ paste0("confidence_", subjectsCorr[s])]]
            # events[[paste0("trial_type_contConfRecognition_",subjectsToCorrelate[ss], "_dummyCoded")]]  <- events[[paste0("trial_type_allConfRecognition_",subjectsToCorrelate[ss], "_dummyCoded")]] * events[[ paste0("confidence_",subjectsToCorrelate[ss])]]
            
            # create effect coded group variable
            events[[paste0("reward_",subjectsCorr[s], "_effectCoded")]] <- ifelse(grepl("cont", subjectsCorr[s]), -1, 1)
            events[[paste0("reward_",subjectsToCorrelate[ss], "_effectCoded")]] <- ifelse(grepl("cont", subjectsToCorrelate[ss]), -1, 1)
            
            # create curiosity-dependent memory variable --> product approach doesn't work
            # events[[paste0("mecu_allConf_",subjectsCorr[s], "_effectCoded")]] <- events[[paste0("trial_type_allConfRecognition_",subjectsCorr[s], "_effectCoded")]] * events[[ paste0("responseCuriosity_", subjectsCorr[s])]]
            # events[[paste0("mecu_allConf_",subjectsToCorrelate[ss], "_effectCoded")]] <- events[[paste0("trial_type_allConfRecognition_",subjectsToCorrelate[ss], "_effectCoded")]] * events[[ paste0("responseCuriosity_",subjectsToCorrelate[ss])]]
            # events[[paste0("mecu_highConf_",subjectsCorr[s], "_effectCoded")]] <- events[[paste0("trial_type_highConfRecognition_",subjectsCorr[s], "_effectCoded")]] * events[[ paste0("responseCuriosity_", subjectsCorr[s])]]
            # events[[paste0("mecu_highConf_",subjectsToCorrelate[ss], "_effectCoded")]] <- events[[paste0("trial_type_highConfRecognition_",subjectsToCorrelate[ss], "_effectCoded")]] * events[[ paste0("responseCuriosity_", subjectsToCorrelate[ss])]]
            # events[[paste0("mecu_aboveAvgConf_",subjectsCorr[s], "_effectCoded")]] <- events[[paste0("trial_type_aboveAverageConfRecognition_",subjectsCorr[s], "_effectCoded")]] * events[[ paste0("responseCuriosity_", subjectsCorr[s])]]
            # events[[paste0("mecu_aboveAvgConf_",subjectsToCorrelate[ss], "_effectCoded")]] <- events[[paste0("trial_type_aboveAverageConfRecognition_",subjectsToCorrelate[ss], "_effectCoded")]] * events[[ paste0("responseCuriosity_", subjectsToCorrelate[ss])]]
            # events[[paste0("mecu_contConf_",subjectsCorr[s], "_dummyCoded")]] <- events[[paste0("trial_type_contConfRecognition_",subjectsCorr[s], "_dummyCoded")]] * events[[ paste0("responseCuriosity_", subjectsCorr[s])]]
            # events[[paste0("mecu_contConf_",subjectsToCorrelate[ss], "_dummyCoded")]] <- events[[paste0("trial_type_contConfRecognition_",subjectsToCorrelate[ss], "_dummyCoded")]] * events[[ paste0("responseCuriosity_",subjectsToCorrelate[ss])]]
            
            # create reward-dependent memory variable --> product approach doesn't work
            # events[[paste0("memo_allConf_",subjectsCorr[s], "_effectCoded")]] <- events[[paste0("trial_type_allConfRecognition_",subjectsCorr[s], "_effectCoded")]] * events[[ paste0("reward_", subjectsCorr[s], "_effectCoded")]]
            # events[[paste0("memo_allConf_",subjectsToCorrelate[ss], "_effectCoded")]] <- events[[paste0("trial_type_allConfRecognition_",subjectsToCorrelate[ss], "_effectCoded")]] * events[[ paste0("reward_",subjectsToCorrelate[ss], "_effectCoded")]]
            # events[[paste0("memo_highConf_",subjectsCorr[s], "_effectCoded")]] <- events[[paste0("trial_type_highConfRecognition_",subjectsCorr[s], "_effectCoded")]] * events[[ paste0("reward_", subjectsCorr[s], "_effectCoded")]]
            # events[[paste0("memo_highConf_",subjectsToCorrelate[ss], "_effectCoded")]] <- events[[paste0("trial_type_highConfRecognition_",subjectsToCorrelate[ss], "_effectCoded")]] * events[[ paste0("reward_", subjectsToCorrelate[ss], "_effectCoded")]]
            # events[[paste0("memo_aboveAvgConf_",subjectsCorr[s], "_effectCoded")]] <- events[[paste0("trial_type_aboveAverageConfRecognition_",subjectsCorr[s], "_effectCoded")]] * events[[ paste0("reward_", subjectsCorr[s], "_effectCoded")]]
            # events[[paste0("memo_aboveAvgConf_",subjectsToCorrelate[ss], "_effectCoded")]] <- events[[paste0("trial_type_aboveAverageConfRecognition_",subjectsToCorrelate[ss], "_effectCoded")]] * events[[ paste0("reward_", subjectsToCorrelate[ss], "_effectCoded")]]
            # events[[paste0("memo_contConf_",subjectsCorr[s], "_dummyCoded")]] <- events[[paste0("trial_type_contConfRecognition_",subjectsCorr[s], "_dummyCoded")]] * events[[ paste0("reward_", subjectsCorr[s], "_effectCoded")]]
            # events[[paste0("memo_contConf_",subjectsToCorrelate[ss], "_dummyCoded")]] <- events[[paste0("trial_type_contConfRecognition_",subjectsToCorrelate[ss], "_dummyCoded")]] * events[[ paste0("reward_",subjectsToCorrelate[ss], "_effectCoded")]]
            
            ### here we start a loop for each member of the pair
            pair <- c(subjectsCorr[s], subjectsToCorrelate[ss])
            # create an index variable
            x <- x+1
            
            # fill in information to dataTable and dataTable_ISC
            dataTable[x,1] <- subjectsCorr[s]
            dataTable[x,2] <- subjectsToCorrelate[ss]          
            
            dataTable_ISC[x,1] <- subjectsCorr[s]
            dataTable_ISC[x,2] <- subjectsToCorrelate[ss]
            dataTable_ISC[x,3] <- ifelse(grepl("cont", subjectsCorr[s]) & grepl("cont", subjectsToCorrelate[ss]), "G11", 
                                         ifelse(grepl("cont", subjectsCorr[s]) & grepl("exp", subjectsToCorrelate[ss]), "G12", 
                                                ifelse(grepl("exp", subjectsCorr[s]) & grepl("exp", subjectsToCorrelate[ss]), "G22",NA )))
            
            dataTable_ISC[x,4] <- paste0("ISC_",subjectsCorr[s],subjectsToCorrelate[ss],"_magictrickwatching_z.nii.gz")
            if(x < N) {
              dataTable_ISC[x,5] <- '\\' #add back slash at end of the row
            }
            
            # fill in information to dataTable_ISC_dummy
            dataTable_ISC_dummy[x,1] <- subjectsCorr[s]
            dataTable_ISC_dummy[x,2] <- subjectsToCorrelate[ss]
            #we adopt deviation coding for the two groups by replacing two groups G1 (cont) and G2 (exp) with -0.5 and 0.5. 
            #Then add up the two values for each row (each subject pair), resulting in three possible values of 1, -1 and 0.
            dataTable_ISC_dummy[x,3] <- ifelse(grepl("cont", subjectsCorr[s]) & grepl("cont", subjectsToCorrelate[ss]), -1, 
                                               ifelse(grepl("cont", subjectsCorr[s]) & grepl("exp", subjectsToCorrelate[ss]), 0, 
                                                      ifelse(grepl("exp", subjectsCorr[s]) & grepl("exp", subjectsToCorrelate[ss]), 1,NA )))
            
            # here we need to add contConfRecognition !!
            
            #curiosity and curiosity interaction
            dataTable_ISC_dummy[x,4] <- cor(events[, paste0("responseCuriosity_", subjectsCorr[s])], events[, paste0("responseCuriosity_", subjectsToCorrelate[ss])] )
            dataTable_ISC_dummy[x,5] <- dataTable_ISC_dummy[x,4] * dataTable_ISC_dummy[x,3] 
            dataTable_ISC_dummy[x,6] <- cor(events[, paste0("confidence_", subjectsCorr[s])], events[, paste0("confidence_", subjectsToCorrelate[ss])] )
            dataTable_ISC_dummy[x,7] <- dataTable_ISC_dummy[x,6] * dataTable_ISC_dummy[x,3] 
            
            #memory and reward-memory interaction
            dataTable_ISC_dummy[x,8] <- cor(events[, paste0("trial_type_cuedRecallStrict_", subjectsCorr[s], "_effectCoded")], events[, paste0("trial_type_cuedRecallStrict_", subjectsToCorrelate[ss], "_effectCoded")])
            dataTable_ISC_dummy[x,9] <- dataTable_ISC_dummy[x,8] * dataTable_ISC_dummy[x,3] 
            dataTable_ISC_dummy[x,10] <- cor(events[, paste0("trial_type_cuedRecallLenient_", subjectsCorr[s], "_effectCoded")], events[, paste0("trial_type_cuedRecallLenient_", subjectsToCorrelate[ss], "_effectCoded")])
            dataTable_ISC_dummy[x,11] <- dataTable_ISC_dummy[x,10] * dataTable_ISC_dummy[x,3] 
            
            dataTable_ISC_dummy[x,12] <- cor(events[, paste0("trial_type_allConfRecognition_", subjectsCorr[s], "_effectCoded")], events[, paste0("trial_type_allConfRecognition_", subjectsToCorrelate[ss], "_effectCoded")])
            dataTable_ISC_dummy[x,13] <- dataTable_ISC_dummy[x,12] * dataTable_ISC_dummy[x,3] 
            dataTable_ISC_dummy[x,14] <- cor(events[, paste0("trial_type_highConfRecognition_", subjectsCorr[s], "_effectCoded")], events[, paste0("trial_type_highConfRecognition_", subjectsToCorrelate[ss], "_effectCoded")])
            dataTable_ISC_dummy[x,15] <- dataTable_ISC_dummy[x,14] * dataTable_ISC_dummy[x,3] 
            dataTable_ISC_dummy[x,16] <- cor(events[, paste0("trial_type_aboveAverageConfRecognition_", subjectsCorr[s], "_effectCoded")], events[, paste0("trial_type_aboveAverageConfRecognition_", subjectsToCorrelate[ss], "_effectCoded")])
            dataTable_ISC_dummy[x,17] <- dataTable_ISC_dummy[x,16] * dataTable_ISC_dummy[x,3]
            # dataTable_ISC_dummy[x,14] <- cor(events[, paste0("trial_type_contConfRecognition_", subjectsCorr[s], "_dummyCoded")], events[, paste0("trial_type_contConfRecognition_", subjectsToCorrelate[ss], "_dummyCoded")])
            # dataTable_ISC_dummy[x,15] <- dataTable_ISC_dummy[x,14] * dataTable_ISC_dummy[x,3]
            
            # #curiosity-dependent memory and curiosity-dependent memory interaction  --> product approach doesn't work
            # dataTable_ISC_dummy[x,14] <- cor(events[, paste0("mecu_allConf_", subjectsCorr[s], "_effectCoded")], events[, paste0("mecu_allConf_", subjectsToCorrelate[ss], "_effectCoded")])
            # dataTable_ISC_dummy[x,15] <- dataTable_ISC_dummy[x,14] * dataTable_ISC_dummy[x,3] 
            # dataTable_ISC_dummy[x,16] <- cor(events[, paste0("mecu_highConf_", subjectsCorr[s], "_effectCoded")], events[, paste0("mecu_highConf_", subjectsToCorrelate[ss], "_effectCoded")])
            # dataTable_ISC_dummy[x,17] <- dataTable_ISC_dummy[x,16] * dataTable_ISC_dummy[x,3] 
            # dataTable_ISC_dummy[x,18] <- cor(events[, paste0("mecu_aboveAvgConf_", subjectsCorr[s], "_effectCoded")], events[, paste0("mecu_aboveAvgConf_", subjectsToCorrelate[ss], "_effectCoded")])
            # dataTable_ISC_dummy[x,19] <- dataTable_ISC_dummy[x,18] * dataTable_ISC_dummy[x,3] 
            # dataTable_ISC_dummy[x,20] <- cor(events[, paste0("mecu_contConf_", subjectsCorr[s], "_dummyCoded")], events[, paste0("mecu_contConf_", subjectsToCorrelate[ss], "_dummyCoded")])
            # dataTable_ISC_dummy[x,21] <- dataTable_ISC_dummy[x,20] * dataTable_ISC_dummy[x,3] 
            
            # #reward-dependent memory and reward-dependent memory interaction  --> product approach doesn't work
            # dataTable_ISC_dummy[x,22] <- cor(events[, paste0("memo_allConf_", subjectsCorr[s], "_effectCoded")], events[, paste0("memo_allConf_", subjectsToCorrelate[ss], "_effectCoded")])
            # dataTable_ISC_dummy[x,23] <- dataTable_ISC_dummy[x,22] * dataTable_ISC_dummy[x,3] 
            # dataTable_ISC_dummy[x,24] <- cor(events[, paste0("memo_highConf_", subjectsCorr[s], "_effectCoded")], events[, paste0("memo_highConf_", subjectsToCorrelate[ss], "_effectCoded")])
            # dataTable_ISC_dummy[x,25] <- dataTable_ISC_dummy[x,24] * dataTable_ISC_dummy[x,3] 
            # dataTable_ISC_dummy[x,26] <- cor(events[, paste0("memo_aboveAvgConf_", subjectsCorr[s], "_effectCoded")], events[, paste0("memo_aboveAvgConf_", subjectsToCorrelate[ss], "_effectCoded")])
            # dataTable_ISC_dummy[x,27] <- dataTable_ISC_dummy[x,26] * dataTable_ISC_dummy[x,3] 
            # dataTable_ISC_dummy[x,28] <- cor(events[, paste0("memo_contConf_", subjectsCorr[s], "_dummyCoded")], events[, paste0("memo_contConf_", subjectsToCorrelate[ss], "_dummyCoded")])
            # dataTable_ISC_dummy[x,29] <- dataTable_ISC_dummy[x,28] * dataTable_ISC_dummy[x,3] 
            
            #curiosity-driven memory benefit and curiosity-driven memory benefit interaction
            dataTable_ISC_dummy[x,18] <- (MAGMOT$curiosityBeta_cuedRecallStrict_c[MAGMOT$ID == subjects[s]] - MAGMOT$curiosityBeta_cuedRecallStrict_c[MAGMOT$BIDS == subjectsToCorrelate[ss]])/2
            dataTable_ISC_dummy[x,19] <- dataTable_ISC_dummy[x,18] * dataTable_ISC_dummy[x,3] 
            dataTable_ISC_dummy[x,20] <- (MAGMOT$curiosityBeta_cuedRecallLenient_c[MAGMOT$ID == subjects[s]] - MAGMOT$curiosityBeta_cuedRecallLenient_c[MAGMOT$BIDS == subjectsToCorrelate[ss]])/2
            dataTable_ISC_dummy[x,21] <- dataTable_ISC_dummy[x,20] * dataTable_ISC_dummy[x,3] 
            #dataTable_ISC_dummy[x,30] <- 1-abs(MAGMOT$curiosityBeta_allConf[MAGMOT$ID == subjects[s]] - MAGMOT$curiosityBeta_allConf[MAGMOT$BIDS == subjectsToCorrelate[ss]])
            dataTable_ISC_dummy[x,22] <- (MAGMOT$curiosityBeta_allConf_c[MAGMOT$ID == subjects[s]] - MAGMOT$curiosityBeta_allConf_c[MAGMOT$BIDS == subjectsToCorrelate[ss]])/2
            dataTable_ISC_dummy[x,23] <- dataTable_ISC_dummy[x,22] * dataTable_ISC_dummy[x,3] 
            # dataTable_ISC_dummy[x,32] <- 1-abs(MAGMOT$curiosityBeta_highConf[MAGMOT$ID == subjects[s]] - MAGMOT$curiosityBeta_highConf[MAGMOT$BIDS == subjectsToCorrelate[ss]])
            dataTable_ISC_dummy[x,24] <- (MAGMOT$curiosityBeta_highConf_c[MAGMOT$ID == subjects[s]] - MAGMOT$curiosityBeta_highConf_c[MAGMOT$BIDS == subjectsToCorrelate[ss]])/2
            dataTable_ISC_dummy[x,25] <- dataTable_ISC_dummy[x,24] * dataTable_ISC_dummy[x,3] 
            # dataTable_ISC_dummy[x,34] <- 1-abs(MAGMOT$curiosityBeta_aboveAvgConf[MAGMOT$ID == subjects[s]] - MAGMOT$curiosityBeta_aboveAvgConf[MAGMOT$BIDS == subjectsToCorrelate[ss]])
            dataTable_ISC_dummy[x,26] <- (MAGMOT$curiosityBeta_aboveAvgConf_c[MAGMOT$ID == subjects[s]] - MAGMOT$curiosityBeta_aboveAvgConf_c[MAGMOT$BIDS == subjectsToCorrelate[ss]])/2
            dataTable_ISC_dummy[x,27] <- dataTable_ISC_dummy[x,26] * dataTable_ISC_dummy[x,3] 
            
            # # determine the unique contribution of curiosity, memory and confidence
            dataTable_ISC_dummy[,28] <- NA #residuals(uniqueCurAboveAvgConf)
            dataTable_ISC_dummy[,29] <-  NA #dataTable_ISC_dummy[x,24] * dataTable_ISC_dummy[x,3] 
            dataTable_ISC_dummy[,30] <- NA #residuals(uniqueMemAboveAvgConf)
            dataTable_ISC_dummy[,31] <- NA #dataTable_ISC_dummy[x,26] * dataTable_ISC_dummy[x,3] 
            dataTable_ISC_dummy[,32] <- NA #residuals(uniqueConf)
            dataTable_ISC_dummy[,33] <- NA #dataTable_ISC_dummy[x,28] * dataTable_ISC_dummy[x,3] 
            
            
            dataTable_ISC_dummy[x,34] <- paste0("ISC_",subjectsCorr[s],subjectsToCorrelate[ss],"_magictrickwatching_z.nii.gz")
            if(x < N) {
              dataTable_ISC_dummy[x,35] <- '\\' #add back slash at end of the row
            }
            
            # HERE DELETE ALL "_effectCoded" variables in events!!!!!!!!!!!!
            
            for (p in seq_along(pair)) {
              
              # disentangle events for pair[p]
              SME_events <- events[,grep(pair[p], colnames(events))] #picks all columns relating to one of the subjects in the pair
              SME_events$stim_file <- events$stim_file
              SME_events$avgVidDur_MAGMOT <- events$avgVidDur_MAGMOT
              SME_events$behavParcel_allConf <- events$behavParcel_allConf
              SME_events$behavParcel_highConf <- events$behavParcel_highConf
              SME_events$behavParcel_aboveAvgConf <- events$behavParcel_aboveAvgConf
              
              names(SME_events)[names(SME_events)=="behavParcel_allConf"] = paste0("behavParcel_allConf_",subjectsCorr[s],subjectsToCorrelate[ss])
              names(SME_events)[names(SME_events)=="behavParcel_highConf"] = paste0("behavParcel_highConf_",subjectsCorr[s],subjectsToCorrelate[ss])
              names(SME_events)[names(SME_events)=="behavParcel_aboveAvgConf"] = paste0("behavParcel_aboveAvgConf_",subjectsCorr[s],subjectsToCorrelate[ss])
              
              for (o in seq_along(SME_outcome)){
                
                for (c in seq_along(confidence_levels)){
                  
                  # subset data depnding on memory performance
                  SME_events_outcome <-  subset(SME_events, SME_events[[paste0("behavParcel_", confidence_levels[c], "_", subjectsCorr[s],subjectsToCorrelate[ss])]] == paste0(SME_outcome[o]))
                  if (feedback == "yes"){
                    print(paste("behavParcel", confidence_levels[c], subjectsCorr[s],subjectsToCorrelate[ss], "has", dim(SME_events_outcome)[1], "rows"))
                  }
                  
                  # save data for each of the subjects
                  if(dim(SME_events_outcome)[1]>0){
                    preprocessedEventsPairDir <- file.path(preprocessedEventsRootDir, paste0(subjectsCorr[s],subjectsToCorrelate[ss]))
                    ifelse(!dir.exists(preprocessedEventsPairDir), dir.create(preprocessedEventsPairDir), FALSE)
                    setwd(preprocessedEventsPairDir)
                    
                    #write.table(SME_events_outcome, file = paste0(subjectsCorr[s], subjectsToCorrelate[ss], "_", pair[p], "_", SME_outcome[o], "_", confidence_levels[c], "_SME_concat.tsv"), quote=FALSE, sep="\t", row.names = FALSE, na = "n/a")
                    
                    
                    # # if study drive is connected, save it there too
                    # if (dir.exists(dirVM)==T && overwrite == "yes"){
                    #   preprocessedEventsPairDirVM <- file.path(preprocessedEventsRootDirVM,"concat", paste0(subjectsCorr[s],subjectsToCorrelate[ss]))
                    #   ifelse(!dir.exists(preprocessedEventsPairDirVM), dir.create(preprocessedEventsPairDirVM), FALSE)
                    #   setwd(preprocessedEventsPairDirVM)
                    #   if (SME_outcome[o] != "differentResponses"){
                    #     write.table(SME_events_outcome, file = paste0(subjectsCorr[s], subjectsToCorrelate[ss], "_", pair[p], "_", SME_outcome[o], "_", confidence_levels[c], "_SME_concat.tsv"), quote=FALSE, sep="\t", row.names = FALSE, na = "n/a")
                    #     setwd(preprocessedEventsPairDir)                      
                    #   }
                    #   if (o == 1 && p == 1 && c == 1){
                    #     print(paste0("saving SME files for pair ",subjectsCorr[s],subjectsToCorrelate[ss],  " to VM..."))
                    #   }
                    # } else if (dir.exists(dirVM)==F && o == 1){
                    #   print(paste0("trying to save SME files for pair ",subjectsCorr[s],subjectsToCorrelate[ss],  " to VM, but not connected to study drive"))
                    # } else if (overwrite == "no"){
                    #   print(paste0("trying to save SME files for pair ",subjectsCorr[s],subjectsToCorrelate[ss],  " to VM, but cannot overwrite files on VM"))
                    # }
                    
                    
                    if (p == 1){
                      if(SME_outcome[o] == "bothRemembered"){
                        
                        dataTable[x,c+2] <- dim(SME_events_outcome)[1]
                        
                        if (confidence_levels[c] == "highConf"){
                          
                          # fill in information to dataTable
                          dataTable_highConf[x,1] <- subjectsCorr[s]
                          dataTable_highConf[x,2] <- subjectsToCorrelate[ss]
                          dataTable_highConf[x,3] <- ifelse(grepl("cont", subjectsCorr[s]) & grepl("cont", subjectsToCorrelate[ss]), "G11", 
                                                            ifelse(grepl("cont", subjectsCorr[s]) & grepl("exp", subjectsToCorrelate[ss]), "G12", 
                                                                   ifelse(grepl("exp", subjectsCorr[s]) & grepl("exp", subjectsToCorrelate[ss]), "G22",NA )))
                          dataTable_highConf[x,4] <- paste0("ISC_",subjectsCorr[s],subjectsToCorrelate[ss],"_SME_highConf.nii.gz")
                          dataTable_highConf[x,5] <- '\\' #add back slash at end of the row
                          
                        } else if(confidence_levels[c] == "aboveAvgConf"){
                          
                          dataTable_aboveAvgConf[x,1] <- subjectsCorr[s]
                          dataTable_aboveAvgConf[x,2] <- subjectsToCorrelate[ss]
                          dataTable_aboveAvgConf[x,3] <- ifelse(grepl("cont", subjectsCorr[s]) & grepl("cont", subjectsToCorrelate[ss]), "G11", 
                                                                ifelse(grepl("cont", subjectsCorr[s]) & grepl("exp", subjectsToCorrelate[ss]), "G12", 
                                                                       ifelse(grepl("exp", subjectsCorr[s]) & grepl("exp", subjectsToCorrelate[ss]), "G22",NA )))
                          dataTable_aboveAvgConf[x,4] <- paste0("ISC_",subjectsCorr[s],subjectsToCorrelate[ss],"_SME_aboveAvgConf.nii.gz")
                          dataTable_aboveAvgConf[x,5] <- '\\' #add back slash at end of the row
                        }
                      } else if(SME_outcome[o] == "bothForgotten"){ # for allConf we need to look at bothForgotten as some pairs don't have any magictricks that they both have forgotten
                        
                        dataTable[x,c+5] <- dim(SME_events_outcome)[1]
                        
                        if(confidence_levels[c] == "allConf"){
                          
                          # fill in information to dataTable
                          dataTable_allConf[x,1] <- subjectsCorr[s]
                          dataTable_allConf[x,2] <- subjectsToCorrelate[ss]
                          dataTable_allConf[x,3] <- ifelse(grepl("cont", subjectsCorr[s]) & grepl("cont", subjectsToCorrelate[ss]), "G11", 
                                                           ifelse(grepl("cont", subjectsCorr[s]) & grepl("exp", subjectsToCorrelate[ss]), "G12", 
                                                                  ifelse(grepl("exp", subjectsCorr[s]) & grepl("exp", subjectsToCorrelate[ss]), "G22",NA )))
                          dataTable_allConf[x,4] <- paste0("ISC_",subjectsCorr[s],subjectsToCorrelate[ss],"_SME_allConf.nii.gz")
                          dataTable_allConf[x,5] <- '\\' #add back slash at end of the row
                          
                        }
                        
                      } #else {
                      #print(paste0("pair ",subjectsCorr[s],subjectsToCorrelate[ss],  " does not have any ", SME_outcome[o], "_", confidence_levels[c], " magictricks." ))
                      #}
                    }
                  }
                }
              }
            }
          } # end of subjCorrel
        }
      }
      
      
      
      # process 3dISC -dataTable
      dataTable_allConf <- na.omit(dataTable_allConf) # remove na rows
      dataTable_aboveAvgConf <- na.omit(dataTable_aboveAvgConf) # remove na rows
      dataTable_highConf <- na.omit(dataTable_highConf) # remove na rows
      
      dataTable_ISC[nrow(dataTable_ISC),ncol(dataTable_ISC)] <- NA # no // for last column
      dataTable_ISC_dummy[nrow(dataTable_ISC_dummy),ncol(dataTable_ISC_dummy)] <- NA # no // for last column
      dataTable_allConf[nrow(dataTable_allConf),ncol(dataTable_allConf)] <- NA  # no // for last column
      dataTable_highConf[nrow(dataTable_highConf),ncol(dataTable_highConf)] <- NA  # no // for last column
      dataTable_aboveAvgConf[nrow(dataTable_aboveAvgConf),ncol(dataTable_aboveAvgConf)] <- NA  # no // for last column
      
      names(dataTable_ISC) <- c("Subj1", "Subj2", "grp", "InputFile", "\\")
      
      names(dataTable_ISC_dummy) <- c("Subj1", "Subj2", "grp", "corrCuriosity", "grpCorrCuriosity", "corrConfidence", "grpCorrConfidence",
                                      "corrRecallStrict", "grpCorrRecallStrict", "corrRecallLenient", "grpCorrRecallLenient",
                                      "corrAllConf", "grpCorrAllConf", "corrHighConf", "grpCorrHighConf", "corrAboveAvgConf", "grpCorrAboveAvgConf", #"corrContConf", "grpCorrContConf", 
                                      #"mecuAllConf", "grpMecuAllConf", "mecuHighConf", "grpMecuHighConf", "mecuAboveAvgConf", "grpMecuAboveAvgConf", "mecuContConf", "grpMecuContConf", 
                                      #"memoAllConf", "grpMemoAllConf", "memoHighConf", "grpMemoHighConf", "memoAboveAvgConf", "grpMemoAboveAvgConf", "memoContConf", "grpMemoContConf", 
                                      "cuBetaRecallStrict", "grpCuBetaRecallStrict", "cuBetaRecallLenient", "grpCuBetaRecallLenient",
                                      "cuBetaAllConf", "grpCuBetaAllConf", "cuBetaHighConf", "grpCuBetaHighConf", "cuBetaAboveAvgConf", "grpCuBetaAboveAvgConf", 
                                      "uniqueCurAboveAvgConf", "grpUniqueCurAboveAvgConf", "uniqueMemAboveAvgConf", "grpUniqueMemAboveAvgConf", "uniqueConfidence", "grpUniqueConfidence",
                                      "InputFile", "\\")
      
      # determine the unique contribution of curiosity, memory and confidence
      uniqueCurAboveAvgConf <- lm(dataTable_ISC_dummy$corrCuriosity ~ dataTable_ISC_dummy$corrAboveAvgConf)
      plot(predict(uniqueCurAboveAvgConf), rstandard(uniqueCurAboveAvgConf))
      dataTable_ISC_dummy[,28] <- residuals(uniqueCurAboveAvgConf)
      dataTable_ISC_dummy[,29] <- dataTable_ISC_dummy[,24] * dataTable_ISC_dummy[,3] 
      
      uniqueMemAboveAvgConf <- lm(dataTable_ISC_dummy$corrAboveAvgConf ~ dataTable_ISC_dummy$corrCuriosity)
      plot(predict(uniqueMemAboveAvgConf), rstandard(uniqueMemAboveAvgConf))
      dataTable_ISC_dummy[,30] <- residuals(uniqueMemAboveAvgConf)
      dataTable_ISC_dummy[,31] <- dataTable_ISC_dummy[,26] * dataTable_ISC_dummy[,3] 
      
      uniqueConf <- lm(dataTable_ISC_dummy$corrConfidence ~ dataTable_ISC_dummy$corrAllConf)
      plot(predict(uniqueConf), rstandard(uniqueConf))
      dataTable_ISC_dummy[,32] <- residuals(uniqueConf)
      dataTable_ISC_dummy[,33] <- dataTable_ISC_dummy[,28] * dataTable_ISC_dummy[,3] 

      names(dataTable_aboveAvgConf) <- c("Subj1", "Subj2", "grp", "InputFile", "\\")
      names(dataTable_highConf) <- c("Subj1", "Subj2", "grp", "InputFile", "\\")
      names(dataTable_allConf) <- c("Subj1", "Subj2", "grp", "InputFile", "\\")
      names(dataTable) <- c("Subj1", "Subj2", paste0(confidence_levels, "_remembered"),  paste0(confidence_levels, "_forgotten"))
      
      
      # round values
      for (c in (seq_along(colnames(dataTable_ISC_dummy)))){
        currentCol <- colnames(dataTable_ISC_dummy)[c]
        if(is.numeric(dataTable_ISC_dummy[names(dataTable_ISC_dummy)==currentCol]) == T){
          dataTable_ISC_dummy[names(dataTable_ISC_dummy)==currentCol] <- round(dataTable_ISC_dummy[names(dataTable_ISC_dummy)==currentCol], digits = 1)
        }
      }
      
      #save files
      setwd(preprocessedEventsRootDir)
      write.table(dataTable_ISC, file="dataTable_magictrickwatching.txt", quote=FALSE, sep="\t", row.names = FALSE, na = "")
      write.table(dataTable_ISC_dummy, file="dataTable_magictrickwatching_memo.txt", quote=FALSE, sep="\t", row.names = FALSE, na = "")
      write.table(dataTable_aboveAvgConf, file="dataTable_aboveAvgConf.txt", quote=FALSE, sep="\t", row.names = FALSE, na = "")
      write.table(dataTable_highConf, file="dataTable_highConf.txt", quote=FALSE, sep="\t", row.names = FALSE, na = "")
      write.table(dataTable_allConf, file="dataTable_allConf.txt", quote=FALSE, sep="\t", row.names = FALSE, na = "")
      write.table(dataTable, file="dataTable_pairwise_memoryScore.csv", quote=FALSE, sep=",", row.names = FALSE, na = "NA")
      
      
      if (dir.exists(dirVM)==T & overwrite == "yes"){
        setwd(file.path(dirVM, 'derivatives', 'magictrickwatching', 'analyses'))
        write.table(dataTable, file=file.path("ISC_magictrickwatching","dataTable.txt"), quote=FALSE, sep="\t", row.names = FALSE, na = "")
        write.table(dataTable_aboveAvgConf, file=file.path("ISC_SME_aboveAvgConf","dataTable_aboveAvgConf"), quote=FALSE, sep="\t", row.names = FALSE, na = "")
        write.table(dataTable_highConf, file=file.path("ISC_SME_highConf","dataTable_highConf"), quote=FALSE, sep="\t", row.names = FALSE, na = "")
        write.table(dataTable_allConf, file=file.path("ISC_SME_allConf","dataTable_allConf"), quote=FALSE, sep="\t", row.names = FALSE, na = "")
      } else {
        if (feedback == "yes" & overwrite == "yes"){
          print("Trying to save information for 3dISC -dataTable, but not connected to study drive!!!")
        }
      }
    }
  }
}

#corrCoefficients <- dataTable_ISC_dummy[, c("corrCuriosity", "corrConfidence", "corrAllConf", "corrHighConf", "corrAboveAvgConf", "cuBetaAllConf", "cuBetaHighConf", "cuBetaAboveAvgConf")]
#corrCoefficients <- dataTable_ISC_dummy[, c("corrCuriosity", "corrConfidence", "corrAllConf", "corrAboveAvgConf", "cuBetaAllConf", "cuBetaAboveAvgConf")]
#cor(corrCoefficients)