#empty work space, load libraries and functions
rm(list=ls())

doISCprep <- 0
debug <- 1




####################################################################################################################################
############################################################  SET UPS   ############################################################
####################################################################################################################################

devtools::source_url("https://github.com/stefaniemeliss/MAGMOT/blob/master/functions/rbindcolumns.R?raw=TRUE")
library(lme4)
options(scipen=999) # suppresses scientific notation


# define core variables
blockstring <- c("_firstBlock", "_secondBlock", "_thirdBlock", "")
SME_outcome <- c("bothRemembered", "bothForgotten", "differentResponses")
SME_outcome <- c("bothRemembered", "bothForgotten")
TR = 2
memoryLevels <- c("cuedRecallStrict", "cuedRecallLenient", 
                  "recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf", 
                  "rememberedStrictAboveAvg", "rememberedLenientAboveAvg", "rememberedStrictHigh", "rememberedLenientHigh")
memoryLabels <- c("cuedRecallStrict", "cuedRecallLenient", 
                  "allConf", "highConf", "aboveAvgConf", 
                  "rememberedStrictAboveAvg", "rememberedLenientAboveAvg", "rememberedStrictHigh", "rememberedLenientHigh")

# define version 
version <- "MAGMOT"
version_official <- "fmri"

# define whether data on VM should be overwritten
overwrite = "no"

# define whether you want feedback printed to the console
feedback = "no"


# define necessary directories
mainDir <- "~/Dropbox/Reading/PhD/Magictricks/fmri_study"
preDir <- file.path(mainDir, "Data", "MAGMOT_pre")
postDir <- file.path(mainDir, "Data", "MAGMOT_post")
PTBDir <- file.path(mainDir, "Data", "PTB")
brainDir <- file.path(mainDir, "Data", "brain_activation")
dataMemoryDir <- file.path(mainDir, "Data", "magicmemory_fmri", "decrypted")
codedDir <- file.path(mainDir, "Data", "magicmemory_fmri", "coded")
preprocessedDir <- file.path(mainDir, "Data", "preprocessed")
preprocessedDir <- "~/Dropbox/Reading/PhD/Magictricks/preprocessed_data/fmri" #"~/Dropbox/Reading/PhD/Magic tricks/behavioural_study/data_fin/MagicBehavioural_preprocessed"

preprocessedQuestDir <- file.path(preprocessedDir, "quest")
preprocessedMemoryDir <- file.path(preprocessedDir, "memory")
preprocessedLongDir <- file.path(preprocessedDir, "long")
preprocessedShareDir <- file.path(preprocessedDir, "share")

#preprocessedEventsRootDir <- file.path(preprocessedDir, "BIDS_eventfiles") # needed?!
preprocessedEventsRootDir <- file.path(mainDir, "derivatives", "magictrickwatching", "concat")
preprocessedEventsRootDir <- file.path(preprocessedDir, "concat") 
preprocessedEventsRootDir <- file.path(preprocessedDir, "BIDS_eventfiles")
concatRootDir <- file.path(preprocessedDir, "concat") 

dataTableDir <- file.path(preprocessedDir, "dataTable") 


# check whether directories for preprocessed data exist, if so empty them and recreate
ifelse(dir.exists(preprocessedDir), unlink(preprocessedDir, recursive = TRUE), FALSE)
ifelse(!dir.exists(preprocessedDir), dir.create(preprocessedDir), FALSE)
ifelse(!dir.exists(preprocessedQuestDir), dir.create(preprocessedQuestDir), FALSE)
ifelse(!dir.exists(preprocessedMemoryDir), dir.create(preprocessedMemoryDir), FALSE)
ifelse(!dir.exists(preprocessedLongDir), dir.create(preprocessedLongDir), FALSE)
ifelse(!dir.exists(preprocessedShareDir), dir.create(preprocessedShareDir), FALSE)
ifelse(!dir.exists(preprocessedEventsRootDir), dir.create(preprocessedEventsRootDir), FALSE)
ifelse(!dir.exists(concatRootDir), dir.create(concatRootDir), FALSE)
ifelse(!dir.exists(dataTableDir), dir.create(dataTableDir), FALSE)
#ifelse(!dir.exists(preprocessedEventsRootDir), dir.create(preprocessedEventsRootDir), FALSE)

# read in data on curuiosity mean split
if (exists("dfMeans") == F) {
  setwd(preprocessedDir)
  dfMeans <- xlsx::read.xlsx("~/Dropbox/Reading/PhD/Magictricks/fmri_study/Analysis/Tricks/MAGMOT_recognitionAndCuriosity_perTrick.xlsx", sheetName = "Sheet1")#, showWarnings = FALSE)
}
dfMeans <- dfMeans[,c("stimID", "meanCuriosity_MAGMOT", "meanCuriosityStandardised_MAGMOT", "mediansplitCuriosity_MAGMOT", "avgVidDur_MAGMOT")]

####################################################################################################################################
####################################################  PROCESS PSYTOOLKIT DATA   ####################################################
####################################################################################################################################


#################################################### read in MAGMOT_pre ####################################################
setwd(preDir)

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
MAGMOT_pre$sex <- ifelse(MAGMOT_pre$sex.1 == 1, "male", 
                         ifelse(MAGMOT_pre$sex.1 == 2, "female", NA))
MAGMOT_pre$sex.1 <- NULL
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
names(MAGMOT_pre)[names(MAGMOT_pre)=="english_age.1"] <- "ageEnglishAcquisition"
names(MAGMOT_pre)[names(MAGMOT_pre)=="health.1"] <- "health"
names(MAGMOT_pre)[names(MAGMOT_pre)=="neurodisorders.1"] <- "neurodisorders"
names(MAGMOT_pre)[names(MAGMOT_pre)=="participant"] <- "preFile"
names(MAGMOT_pre)[names(MAGMOT_pre)=="DOB.1"] <- "DOB"

### compute duration of pre assessment
MAGMOT_pre$startPre <- as.POSIXct(MAGMOT_pre$TIME_start, format = "%Y-%m-%d-%H-%M") # convert date-time format
MAGMOT_pre$endPre <- as.POSIXct(MAGMOT_pre$TIME_end, format = "%Y-%m-%d-%H-%M") # convert date-time format
MAGMOT_pre$durPre <- difftime(MAGMOT_pre$endPre, MAGMOT_pre$startPre, units='mins')

### inclusion check
items <- c("inclusioncheck.1", "inclusioncheck.2", "inclusioncheck.3", "inclusioncheck.4", "inclusioncheck.5", "inclusioncheck.6", "inclusioncheck.7")

for (item in seq_along(items)){
  # recode
  MAGMOT_pre[[paste0(items[item])]] <- ifelse(MAGMOT_pre[[paste0(items[item])]] == 0, "no",
                                              ifelse(MAGMOT_pre[[paste0(items[item])]] == 1, "yes",
                                                     NA))
}
rm(item, items)

### MRI screening: any yes?
items <- c("screening_MRI.1", "screening_MRI.2", "screening_MRI.3", "screening_MRI.4", "screening_MRI.5", "screening_MRI.6", "screening_MRI.7", "screening_MRI.8", "screening_MRI.9", "screening_MRI.10",
           "screening_MRI.11", "screening_MRI.12", "screening_MRI.13", "screening_MRI.14", "screening_MRI.15", "screening_MRI.16", "screening_MRI.17", "screening_MRI.18", "screening_MRI.19")

# compute summary score for MRI screening answers
for (item in seq_along(items)){
  # if (item == 1){
  #   MAGMOT_pre$screening_MRI <- MAGMOT_pre[[paste0(items[item])]]
  # } else {
  #   MAGMOT_pre$screening_MRI <- MAGMOT_pre$screening_MRI + MAGMOT_pre[[paste0(items[item])]]
  # }
  # recode answers
  MAGMOT_pre[[paste0(items[item])]] <- ifelse(MAGMOT_pre[[paste0(items[item])]] == 0, "no",
                                              ifelse(MAGMOT_pre[[paste0(items[item])]] == 1, "yes",
                                                     NA))
}
rm(item, items)

#### health screening
items <- c("health_current.1", "health_current.2", "health_current.3", "health_current.4", 
           "health_ever.1", "health_ever.2", "health_ever.3", "health_ever.4", "health_ever.5", "health_ever.6",
           "health_ever.7", "health_ever.8", "health_ever.9", "health_ever.10", "health_ever.11", "health_ever.12",
           "health_ever.13", "health_ever.14", "health_ever.15", "health_ever.16", "health_ever.17", "health_ever.18",
           "health_ever.19", "health_ever.20", "health_ever.21", "health_ever.22", "health_ever.23", "health_ever.24")
# compute summary score for health screening answers
for (item in seq_along(items)){
  # if (item == 1){
  #   MAGMOT_pre$health_check <- MAGMOT_pre[[paste0(items[item])]]
  # } else {
  #   MAGMOT_pre$health_check <- MAGMOT_pre$screening_MRI + MAGMOT_pre[[paste0(items[item])]]
  # }
  # recode answers
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
MAGMOT_pre$BISBAS_inhibition <- MAGMOT_pre$BISBAS.15 +  MAGMOT_pre$BISBAS.1 +  MAGMOT_pre$BISBAS.12 +  MAGMOT_pre$BISBAS.2R + MAGMOT_pre$BISBAS.17 + MAGMOT_pre$BISBAS.20R #BIS items: 15, 1, 8, 12, 2 (R), 17, 20 (R)
MAGMOT_pre$BISBAS_rewardresponsiveness <-  MAGMOT_pre$BISBAS.7 +  MAGMOT_pre$BISBAS.4 +  MAGMOT_pre$BISBAS.16 +  MAGMOT_pre$BISBAS.6 + MAGMOT_pre$BISBAS.13 #BAS reward responsiveness: 7, 4, 16, 6, 13
MAGMOT_pre$BISBAS_drive <- MAGMOT_pre$BISBAS.9 +  MAGMOT_pre$BISBAS.3 +  MAGMOT_pre$BISBAS.11 +  MAGMOT_pre$BISBAS.19 #BAS drive: 9, 3, 11, 19
MAGMOT_pre$BISBAS_funseeking <- MAGMOT_pre$BISBAS.10 +  MAGMOT_pre$BISBAS.18 +  MAGMOT_pre$BISBAS.5 +  MAGMOT_pre$BISBAS.14 #BAS fun seeking: 10, 18, 5, 14

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

### create data frame with all the raw questionnaire data
cols <- grepl("ID",names(MAGMOT_pre)) | grepl("inclusioncheck.",names(MAGMOT_pre)) | grepl("screening_MRI.",names(MAGMOT_pre)) | grepl("BISBAS.",names(MAGMOT_pre)) |  grepl("NeedForCognition.",names(MAGMOT_pre)) |  
  grepl("FearOfFailure.",names(MAGMOT_pre))  |  grepl("ApproachAndAvoidanceTemperament.",names(MAGMOT_pre)) |  grepl("health_",names(MAGMOT_pre))  |  grepl("TraitCuriosity.",names(MAGMOT_pre))
questionnaire_raw <- MAGMOT_pre[,  cols ]
rm(cols)


### reduce MAGMOT_pre to relevant variables
MAGMOT_pre <- MAGMOT_pre[, c("ID", "preFile", "startPre", "endPre", "durPre", "corsi.1", "X2nback.1", "age", "DOB", "sex",  "gender", "ethnicity", "education", "yearsOfEducation", 
                             "employment", "studySubject", "english", "ageEnglishAcquisition", "handedness", "vision", "health", "neurodisorders",
                             "BISBAS_inhibition", "BISBAS_rewardresponsiveness", "BISBAS_drive", "BISBAS_funseeking", "NeedForCognition", "FearOfFailure", "ApproachTemperament", "AvoidanceTemperament", "TraitCuriosity")]


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
    nbackData[i, "nback_hits"] <-  nrow(subset(nback, nback$V2 == 1 & nback$V5 == 1)) # letter equal to what was presented two trials ago and ppt indicated this correctly
    nbackData[i, "nback_misses_inclTooSlow"] <-  nrow( subset(nback, nback$V2 != 1 & nback$V5 == 1)) # letter equal to what was presented two trials ago and ppt did not indicated this OR reacted too slow
    nbackData[i, "nback_misses_exclTooSlow"] <-  nrow(subset(nback, nback$V2 == 2 & nback$V5 == 1)) # letter equal to what was presented two trials ago and ppt did not indicated this
    nbackData[i, "nback_correctrejections"] <-  nrow(subset(nback, nback$V2 == 1 & nback$V5 > 1)) # letter not equal to what was presented two trials ago and ppt indicated this correctly
    nbackData[i, "nback_falsealarms_inclTooSlow"] <-  nrow(subset(nback, nback$V2 != 1 & nback$V5 > 1)) # letter not equal to what was presented two trials ago and ppt indicated did not this correctly OR reacted too slow
    nbackData[i, "nback_falsealarms_exclTooSlow"] <-  nrow(subset(nback, nback$V2 == 2 & nback$V5 > 1)) # letter not equal to what was presented two trials ago and ppt indicated did not this correctly
    
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
participantToDelete1 <- as.character(MAGMOT_post$participant[(MAGMOT_post$id.1 == "09" & MAGMOT_post$TIME_start == "2019-05-13-19-20")]) # ppt assessed on 13/05/19 should be removed 
participantToDelete2 <- as.character(MAGMOT_post$participant[(MAGMOT_post$id.1 == 23 & MAGMOT_post$TIME_start == "2019-05-22-11-51")]) # ppt assessed on 13/05/22 should be removed
participantToDelete3 <- as.character(MAGMOT_post$participant[(MAGMOT_post$id.1 == 32 & MAGMOT_post$TIME_start == "2019-06-03-15-42")]) # ppt assessed on 13/06/03 should be removed
participantToDelete4 <- as.character(MAGMOT_post$participant[(MAGMOT_post$id.1 == 46 & MAGMOT_post$TIME_start == "2019-06-26-19-36")]) # ppt assessed on 13/06/26 should be removed
participantToDelete5 <- as.character(MAGMOT_post$participant[(MAGMOT_post$id.1 == 48 & MAGMOT_post$TIME_start == "2019-07-03-12-48")]) # ppt assessed on 13/07/03 should be removed

MAGMOT_post <- subset(MAGMOT_post, MAGMOT_post$participant != participantToDelete1 & MAGMOT_post$participant != participantToDelete2 & MAGMOT_post$participant != participantToDelete3 & MAGMOT_post$participant != participantToDelete4 & MAGMOT_post$participant != participantToDelete5)
rm(list=ls(pattern = "participantToDelete"))
### recode demographic information and then delete the old variable
MAGMOT_post$group <- ifelse(MAGMOT_post$group.1 == 1, "cont", # group
                            ifelse(MAGMOT_post$group.1 == 2, "exp", NA))
MAGMOT_post$group.1 <- NULL
MAGMOT_post$groupEffectCoded <-  ifelse(MAGMOT_post$group == "exp", 1, -1)

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
names(MAGMOT_post)[names(MAGMOT_post)=="participant"] <- "postFile"
names(MAGMOT_post)[names(MAGMOT_post)=="comments_participant.1"] <- "comment_task1"
names(MAGMOT_post)[names(MAGMOT_post)=="comments_participant.2"] <- "comment_task2"
names(MAGMOT_post)[names(MAGMOT_post)=="comments_participant.3"] <- "comment_task3"
names(MAGMOT_post)[names(MAGMOT_post)=="comments_experimenter.1"] <- "comment_exp"

### compute duration of post assessment
MAGMOT_post$startPost <- as.POSIXct(MAGMOT_post$TIME_start, format = "%Y-%m-%d-%H-%M") # convert date-time format
MAGMOT_post$endPost <- as.POSIXct(MAGMOT_post$TIME_end, format = "%Y-%m-%d-%H-%M") # convert date-time format
MAGMOT_post$durPost <- difftime(MAGMOT_post$endPost, MAGMOT_post$startPost, units='mins')

### State Curiosity: compute scale
MAGMOT_post$StateCuriosity <- MAGMOT_post$StateCuriosity.1 + MAGMOT_post$StateCuriosity.2 + MAGMOT_post$StateCuriosity.3 + MAGMOT_post$StateCuriosity.4 + MAGMOT_post$StateCuriosity.5 + MAGMOT_post$StateCuriosity.6 + MAGMOT_post$StateCuriosity.7 +
  MAGMOT_post$StateCuriosity.8 + MAGMOT_post$StateCuriosity.9 + MAGMOT_post$StateCuriosity.10 + MAGMOT_post$StateCuriosity.11 + MAGMOT_post$StateCuriosity.12 + MAGMOT_post$StateCuriosity.13 + MAGMOT_post$StateCuriosity.14 +
  MAGMOT_post$StateCuriosity.15 + MAGMOT_post$StateCuriosity.16 + MAGMOT_post$StateCuriosity.17 + MAGMOT_post$StateCuriosity.18 + MAGMOT_post$StateCuriosity.19 + MAGMOT_post$StateCuriosity.20

### create data frame with all the raw questionnaire data
cols <- grepl("ID",names(MAGMOT_post)) | grepl("StateCuriosity.",names(MAGMOT_post))
questionnaire_raw_2 <- MAGMOT_post[,  cols ]

questionnaire_raw <- merge(questionnaire_raw, questionnaire_raw_2, by = "ID")
rm(cols, questionnaire_raw_2)

# select relevant rows from MAGMOT_post
MAGMOT_post <- MAGMOT_post[,c("ID", "postFile", "startPost", "endPost", "durPost",
                              "group", "groupEffectCoded", "StateCuriosity", 
                              "sleepLastNight", "sleepAverage", 
                              "alcohol", "alcoholAmount", "rewardEffort", "rewardExpectations", 
                              "comment_task1", "comment_task2", "comment_task3", "comment_exp",
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
  if (MAGMOT$group[MAGMOT$ID == subjects[s]] == "exp"){ # define BIDS name depending on group
    BIDSstring = paste0("sub-experimental0", subjects[s])
  } else if (MAGMOT$group[MAGMOT$ID == subjects[s]] == "cont"){
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
  
  # change group from int/ext to cont/exp
  quest$group <- ifelse(quest$group == "int", "cont", "exp")
  
  # assign names to question
  quest$IMI <- ifelse(quest$question == "It was fun to do the experiment.", "post1",
                      ifelse(quest$question == "It was boring to do the experiment.", "post2", #reverse
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
                                                                                                          ifelse(quest$question == "I did not try very hard\\nto do well at this activity.", "post14", #reverse
                                                                                                                 ifelse(quest$question == "I tried very hard on this activity.", "post15",
                                                                                                                        ifelse(quest$question == "It was important to me to do well at this task.", "post16",
                                                                                                                               ifelse(quest$question == "I did not put much energy into this.", "post17", #reverse
                                                                                                                                      ifelse(quest$question == "I did not feel nervous at all while doing this.", "post18", #reverse
                                                                                                                                             ifelse(quest$question == "I felt very tense while doing this activity." , "post19",
                                                                                                                                                    ifelse(quest$question == "I was very relaxed in doing this experiment.", "post20", #reverse
                                                                                                                                                           ifelse(quest$question == "I was anxious while working on this task." , "post21",
                                                                                                                                                                  ifelse(quest$question == "I felt pressured while doing this task.", "post22",
                                                                                                                                                                         ifelse(quest$question == "I tried to find out how many people\\nwill be able to find the solution.", "post23",
                                                                                                                                                                                ifelse(quest$question == "I was able to see the magic tricks properly.", "post24",
                                                                                                                                                                                       NA))))))))))))))))))))))))
  
  # create item name that overlaps with the naming convention from the other questionnaires
  quest$item <- gsub("post", "PostExpAssessment.", quest$IMI)
  
  # recode items questionnaire
  items <- c("post2", "post14", "post17", "post18", "post20") #items to recode
  quest$score_raw <- quest$endValue
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
  questLong <- quest[,c("ID", "fMRI", "group", "motivation", "IMI", "score", "item", "score_raw")]
  questWide <- reshape2::dcast(questLong, ID + fMRI + group + motivation ~ IMI, value.var="score")
  questWide_raw <- reshape2::dcast(questLong, ID + fMRI + group + motivation ~ item, value.var="score_raw")
  
  # rbind the questWide files of each subject to a data frame
  if(s == 1){
    questDataWide <- questWide
    questDataWide_raw <- questWide_raw
  } else {
    temp_questDataWide <-  questWide
    temp_questDataWide_raw <-  questWide_raw
    questDataWide <- rbind.all.columns(questDataWide, temp_questDataWide) #rbind all columns will induce NA if there was initially no data saved in the loop per participant
    questDataWide_raw <- rbind.all.columns(questDataWide_raw, temp_questDataWide_raw) #rbind all columns will induce NA if there was initially no data saved in the loop per participant
    rm(temp_questDataWide, temp_questDataWide_raw)
  }
  
  if (s == length(subjects)){
    # compute IMI scores
    questDataWide[,grep("post", colnames(questDataWide))] <- as.numeric(unlist(questDataWide[,grep("post", colnames(questDataWide))])) # make sure the scores are numeric
    questDataWide$IMI_intrinsicMotivation <- (questDataWide$post1 + questDataWide$post2 + questDataWide$post3)/3 ###intrinsic motivation items
    questDataWide$IMI_taskEngagement <- (questDataWide$post4 + questDataWide$post5 + questDataWide$post6)/3 ###task engagement items
    questDataWide$IMI_interest <- (questDataWide$post7 + questDataWide$post8 + questDataWide$post9)/3 ###interest items
    questDataWide$IMI_boredom <- (questDataWide$post10 + questDataWide$post11 + questDataWide$post12)/3 ###boredom items
    questDataWide$IMI_effort <- (questDataWide$post13 + questDataWide$post14 + questDataWide$post15 + questDataWide$post16 + questDataWide$post17)/5 ####effort/importance
    questDataWide$IMI_pressure <- (questDataWide$post18 + questDataWide$post19 + questDataWide$post20 + questDataWide$post21 + questDataWide$post22)/5 ###pressure/tension
    
    # rename two variables
    names(questDataWide) [names(questDataWide) == "post23"] <- "compliance"
    names(questDataWide) [names(questDataWide) == "post24"] <- "ableToSee"
    
    # combine IMI values with questionnaire_raw
    questionnaire_raw <- merge(questionnaire_raw, questDataWide_raw, by = "ID")
  }
  
  # remove unneccsary values and data
  if (debug == 0){
    rm(quest, questLong, questWide, questWide_raw, respMatQuest)
    
  }
  
  #################################################### read in task data  ####################################################
  setwd(PTBDir)
  taskfile <- paste0("MAGMOT_", subjects[s], "_ptbdata.txt")
  ptbdata <- read.table(file = taskfile, sep = "\t", header = T)
  
  # change participant IDs
  names(ptbdata)[names(ptbdata)=="subject"] <- "ID"
  ptbdata$ID <- gsub("O", "", ptbdata$ID) #replace any O in the file names
  ptbdata$ID <- ifelse(nchar(ptbdata$ID)==1, paste0("0", ptbdata$ID), ptbdata$ID) # add a zeo if not present in front of one digit numbers
  ptbdata$BIDS <- BIDSstring
  ptbdata$groupEffectCoded <-  ifelse(ptbdata$group == "ext", 1, -1)
  ptbdata$group <-  ifelse(ptbdata$group == "ext", "exp", "cont") # note: behavioural studies have cont vs exp as croup names
  
  # change stimID of respMatTask to match collector
  names(ptbdata)[names(ptbdata)=="stimID"] <- "vidFileName"
  ptbdata$stimID <- gsub("^\\d\\d_", "", ptbdata$vidFileName)
  ptbdata$stimID <- gsub("_combined_small.mp4", "", ptbdata$stimID)
  
  # convert timings
  ptbdata$endExperiment_raw <- as.POSIXct(ptbdata$endExperiment, format = "%d-%m-%Y_%H:%M")
  ptbdata$endExperiment_LDN <- as.POSIXct(ptbdata$endExperiment, format = "%d-%m-%Y_%H:%M", tz = "Europe/London")
  ptbdata$endExperiment_UTC <- format(ptbdata$endExperiment_LDN, tz="GMT", usetz=TRUE)
  
  if (is.na(ptbdata$endPractice) == T) {
    ptbdata$endPractice_raw <- NA
    ptbdata$endPractice_LDN <- NA
    ptbdata$endPractice_UTC <- NA
  } else {
    ptbdata$endPractice_raw <- as.POSIXct(ptbdata$endPractice, format = "%d-%m-%Y_%H:%M")
    ptbdata$endPractice_LDN <- as.POSIXct(ptbdata$endPractice, format = "%d-%m-%Y_%H:%M", tz = "Europe/London")
    ptbdata$endPractice_UTC <- format(ptbdata$endPractice_LDN, tz="GMT", usetz=TRUE)
  }

  # add curiosity median split for MAGMOT sample (computed at beginning of script based on old MAGMOT_long.xlsx data set)
  ptbdata <- merge(ptbdata, dfMeans, by = "stimID")
  
  # recode curiosity: replace all curiosity ratings that have exceeded the timeout 
  # NOTE::: this means we should CHANGE/REPLACE responseCuriosity with  curiosity in the following
  ptbdata$curiosity <- ifelse(ptbdata$timestampCuriosity > ptbdata$displayCuriosityOnset + ptbdata$timeoutCuriosity - 3*ptbdata$timingCorrection, NA, ptbdata$responseCuriosity)
  ptbdata$curiosity_RT <- ifelse(ptbdata$timestampCuriosity > ptbdata$displayCuriosityOnset + ptbdata$timeoutCuriosity - 3*ptbdata$timingCorrection, NA, ptbdata$rtCuriosity)
  ptbdata$curiosity_tooSlow <- ifelse(ptbdata$timestampCuriosity > ptbdata$displayCuriosityOnset + ptbdata$timeoutCuriosity - 3*ptbdata$timingCorrection, 1, 0)
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
  
  # dichotimise curiosity
  ptbdata$curiosity_dich <- ifelse(ptbdata$curiosityGroupMeanCentered > 0, 1,
                                   ifelse(ptbdata$curiosityGroupMeanCentered < 0, -1, NA))
  
  # compute reward by curiosity interaction
  ptbdata$rewardByCuriosity <- ptbdata$curiosityGroupMeanCentered * ptbdata$groupEffectCoded
  ptbdata$rewardByCuriosity_updated <- ptbdata$curiosityGroupMeanCentered_updated * ptbdata$groupEffectCoded
  
  # process answer to the question how many people will be able to find a solution to the magic trick --> variable "estimate"
  ptbdata$answer_tooSlow <- ifelse(ptbdata$timestampAnswer > ptbdata$displayAnswerOnset + ptbdata$timeoutAnswer - 3*ptbdata$timingCorrection, 1, 0)
  ptbdata$responseEstimate <- ifelse(ptbdata$responseAnswer == "0 to 10", 0,
                                     ifelse(ptbdata$responseAnswer == "11 to 20", 1,
                                            ifelse(ptbdata$responseAnswer == "21 to 30", 2,
                                                   ifelse(ptbdata$responseAnswer == "31 or more", 3, NaN))))
  ptbdata$rtEstimate <- ptbdata$rtAnswer
  
  # calculate how long a blank screen was displayed for between the end of the magic trick video and the start of the fixation
  ptbdata$displayBlankDuration <- ptbdata$fixationPostVidOnset - ptbdata$displayVidOffset
  
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
  
  f <- paste0("magicmemory_", version_official, "_", subjects[s], ".csv")
  if (file.exists(f)){
    memory <- read.csv(f, header = T)
    memory$memoryFile <- f
  } else { # if file does not exist, print into console and check other spelling of file
    print(paste0("magicmemory_", version_official, "_", subjects[s], ".csv does not exist"))
    
    subjectsAlt <- gsub("0", "", subjects[s])
    f <- paste0("magicmemory_", version_official, "_", subjectsAlt, ".csv")
    if (file.exists(f)){
      memory <- read.csv(f, header = T)
      memory$memoryFile <- f
    } else { # if file does not exist, print into console and use filler file to create NAs
      print(paste0("magicmemory_", version_official, "_", subjectsAlt, ".csv does not exist either"))
      rm(subjectsAlt)
    }
  }
  
  # overwrite username to ensure that it is the same across all data sets
  memory$username <- subjects[s]
  
  # get time stamps
  if("post_0_trial_end_date" %in% colnames(memory)){  # check whether the memory data set has information about when it has started/finished
    
    memory$startMemory <- memory$post_0_trial_end_date[2] # note: no post_0_trial_start_date available
    memory$endMemory <- memory$post_0_trial_end_date[dim(memory)[1]]
    
    # memory$startMemory_UTC <- as.POSIXct(min(memory$post_0_trial_end_ms, na.rm = T)/1000, origin = "1970-01-01", tz = "GMT") # this command produces the same output as memory$startMemory <- memory$post_0_trial_end_date[2]
    memory$startMemory_UTC <- as.POSIXct(min(memory$post_0_trial_start_ms, na.rm = T)/1000, origin = "1970-01-01", tz = "GMT")
    memory$endMemory_UTC <- as.POSIXct(max(memory$post_0_trial_end_ms, na.rm = T)/1000, origin = "1970-01-01", tz = "GMT")
    
    ### compute duration of memory assessment
    memory$durMemory <- difftime(memory$endMemory_UTC, memory$startMemory_UTC, units='mins')
    
  } else { # if not, add a note
    memory$startMemory <- "check manually"
    memory$startMemory_UTC <- "check manually"
    memory$endMemory <- "check manually"
    memory$endMemory_UTC <- "check manually"
    memory$durMemory <- "check manually"
  }
  
  
  #################################################### process the RECALL memory data: select relevant rows and columns, change item counter
  cuedRecall <- subset(memory, memory$trial.type == "magic_recall")
  cuedRecall$itemRecall <- cuedRecall$item-1
  cuedRecall$trialRecall <- c(1:dim(cuedRecall)[1])
  cuedRecall <- cuedRecall[,c("username", "trialRecall", "itemRecall", "stimid", "trickAnswer")]
  names(cuedRecall) <- c("ID", "trialRecall", "itemRecall", "stimID", "responseRecall")
  
  # pre code answers for recall task: in the experiment, participants are asked to insert "no recall" in the field if they cannot recall the magic trick
  # if that has occured, recall performance on that magic trick is coded with 0
  for (j in 1:nrow(cuedRecall)){
    if (is.na(cuedRecall$responseRecall[j]) == T) {
      cuedRecall$cuedRecallStrict[j] = NaN
      cuedRecall$cuedRecallLenient[j] = NaN
    } else if (cuedRecall$responseRecall[j] == "no recall") { # check different spellings of NO RECALL
      cuedRecall$cuedRecallStrict[j] = 0
      cuedRecall$cuedRecallLenient[j] = 0
    } else if (cuedRecall$responseRecall[j] == "No recall") {
      cuedRecall$cuedRecallStrict[j] = 0
      cuedRecall$cuedRecallLenient[j] = 0
    } else if (cuedRecall$responseRecall[j] == "No Recall") {
      cuedRecall$cuedRecallStrict[j] = 0
      cuedRecall$cuedRecallLenient[j] = 0
    } else if (cuedRecall$responseRecall[j] == " No recall") {
      cuedRecall$cuedRecallStrict[j] = 0
      cuedRecall$cuedRecallLenient[j] = 0
    } else if (cuedRecall$responseRecall[j] == "NO recall") {
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
    write.csv(cuedRecall, file = filenameRecall, row.names = F)
    
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
      names(cuedRecall_coded)[names(cuedRecall_coded)=="description"] <- "responseRecall" # change description to responseRecall
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
  names(recognition) <- c("ID", "trialRecognition", "itemRecognition", "stimID", "responseRecognition", "rtRecognition", "recognition", "responseConfidence", "confidenceGroupMeanCentered", "rtConfidence")
  
  # compute scores of recognition performance
  recognition$responseConfidenceCorrectTrials <- ifelse(recognition$recognition == 1, recognition$responseConfidence, NA)
  recognition$confidenceGroupMeanCenteredCorrectTrials <- ifelse(recognition$recognition == 1, recognition$confidenceGroupMeanCentered, NA)
  recognition$recognitionAboveMeanConf <- ifelse(recognition$recognition == 1 & recognition$confidenceGroupMeanCentered > 0, 1, 0)
  
  for (k in 1:6) { #confidence ranges from 1 to 6, potentially code can be made more flexible by using min(data$confidence) and max(data$confidence)
    recognition[[paste0("recognitionConfLevel_", k)]] <- ifelse(recognition$responseConfidence == k & recognition$recognition == 1, 1, 0)
    if (k < 6) {
      recognition[[paste0("recognitionConfLevel_above_", k)]] <- ifelse(recognition$responseConfidence > k & recognition$recognition == 1, 1, 0)
    }
    
    if (k == 1 || k == 3 || k == 5) {
      recognition[[paste0("recognitionConfLevel_", k, "_", k+1)]] <- ifelse(recognition$responseConfidence == k  & recognition$recognition == 1 | recognition$responseConfidence == k+1  & recognition$recognition == 1, 1, 0)
    }
    if (k == 1 || k == 4) {
      recognition[[paste0("recognitionConfLevel_", k, "_", k+1, "_", k+2)]] <- ifelse(recognition$responseConfidence == k  & recognition$recognition == 1 | recognition$responseConfidence == k+1  & recognition$recognition == 1 | recognition$responseConfidence == k+2  & recognition$recognition == 1, 1, 0)
    }
  }
  
  # save recognition data
  setwd(dataMemoryDir)
  if (file.exists(f)) { # if the participants initially participated in the memory part, their preprocessed recall data is saved
    setwd(preprocessedMemoryDir)
    filenameRecog <- paste0(BIDSstring,"_recognition.csv")
    write.table(recognition, file = filenameRecog, sep = ",", row.names = F)
    write.csv(recognition, file = filenameRecog, row.names = F)
  }
  
  # MERGE DATA FROM TASK AND MEMORY TEST
  MEMO <- merge(task, cuedRecall, by = c("ID", "stimID"))
  MEMO <- merge(MEMO, recognition, by = c("ID", "stimID"))
  if (debug == 0){
    rm(cuedRecall, cuedRecall_coded, recognition, task)
  }
  
  # compute whether a trick has been remembered based on Hasson et al. (2008).  
  # They classified an event to be remembered if there was a correct answer using either recall or high confidence recognition.
  MEMO$rememberedStrictAboveAvg<- ifelse(MEMO$cuedRecallStrict == 1 | MEMO$recognitionAboveMeanConf == 1, 1, 0)
  MEMO$rememberedLenientAboveAvg<- ifelse(MEMO$cuedRecallLenient == 1 | MEMO$recognitionAboveMeanConf == 1, 1, 0)
  MEMO$rememberedStrictHigh <- ifelse(MEMO$cuedRecallStrict == 1 | MEMO$recognitionConfLevel_4_5_6 == 1, 1, 0)
  MEMO$rememberedLenientHigh <- ifelse(MEMO$cuedRecallLenient == 1 | MEMO$recognitionConfLevel_4_5_6 == 1, 1, 0)
  
  for (mem in 1:length(memoryLevels)) {
    # calculate curiosity-driven memory memory benefit (continuous, absolute) --> range: [min(curiosityGroupMeanCentered); max(curiosityGroupMeanCentered)]
    MEMO[[paste0("curBen_cont_", memoryLabels[mem])]] <- MEMO$curiosityGroupMeanCentered * MEMO[[paste0(memoryLevels[mem])]]
    # calculate curiosity-driven memory memory benefit (dichotomous, absolute)  --> range: [-1; 1]
    MEMO[[paste0("curBen_dich_", memoryLabels[mem])]] <- MEMO$curiosity_dich * MEMO[[paste0(memoryLevels[mem])]]
  }
  
  # sort MEMO in order of trials of main experiment
  MEMO <- MEMO[order(MEMO$trial),]
  
  # save data in long format
  setwd(preprocessedLongDir)
  filenameLong <- paste0(BIDSstring,"_long.csv")
  write.csv(MEMO, file = filenameLong, row.names = F)
  
  
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
      onset <- run_acq[,c("trial","vidFileName", "displayVidOnset", "mockOffset", "displayVidOffset", "fixationPostVidOnset",
                          "displayAnswerOnset", "responseAnswer", "timestampAnswer", "fixationPostAnswerOnset",
                          "displayCuriosityOnset", "responseCuriosity", "timestampCuriosity", "fixationPostCuriosityOnset", 
                          "responseRecall", "cuedRecallStrict", "cuedRecallLenient",
                          "responseRecognition", "responseConfidence",
                          "recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf", 
                          "rememberedStrictAboveAvg", "rememberedLenientAboveAvg", "rememberedStrictHigh", "rememberedLenientHigh")]
      
      onset$timestampVidOnset <- onset$displayVidOnset
      onset$timestampMockOffset <- onset$mockOffset
      onset$timestampPostVidFixation <- onset$fixationPostVidOnset
      # onset <- reshape2::melt(onset, id.vars=c("vidFileName", "trial", "mockOffset", "displayVidOffset", "timestampPostVidFixation", "vidDurCalc", "cuedRecallStrict", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf",
      #                                          "responseAnswer", "timestampAnswer", "responseCuriosity", "timestampCuriosity"), value.name = "onset")
      onset <- reshape2::melt(onset, id.vars=c("vidFileName", "trial", "timestampVidOnset", "timestampMockOffset", "displayVidOffset", "timestampPostVidFixation",
                                               "responseRecall", "cuedRecallStrict", "cuedRecallLenient", 
                                               "responseRecognition", "responseConfidence",
                                               "recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf", 
                                               "rememberedStrictAboveAvg", "rememberedLenientAboveAvg", "rememberedStrictHigh", "rememberedLenientHigh",
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
      
      for (mem in 1:length(memoryLevels)) {
        # code memory performance and trial type (remembered vs forgotten)
        run_BIDS[[paste0(memoryLabels[mem])]] <- ifelse(run_BIDS$variable == "displayVid", run_BIDS[[paste0(memoryLevels[mem])]], "not applicable")
        # Define trial types for each memory cut off seperately 
        run_BIDS[[paste0("trial_type_", memoryLabels[mem])]] <-  ifelse(run_BIDS[[paste0(memoryLabels[mem])]] == 1, "remembered", 
                                                                        ifelse(run_BIDS[[paste0(memoryLabels[mem])]] == 0, "forgotten",
                                                                               ifelse(run_BIDS[[paste0(memoryLabels[mem])]] == "not applicable", "not applicable",
                                                                                      "undefined")))
      }
      
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
      
      run_BIDS$response_recall <-  ifelse(run_BIDS$variable == "displayVid", as.character(run_BIDS$responseRecall), "not applicable")
      run_BIDS$response_recognition <-  ifelse(run_BIDS$variable == "displayVid", run_BIDS$responseRecognition, "not applicable")
      run_BIDS$response_confidence <-  ifelse(run_BIDS$variable == "displayVid", run_BIDS$responseConfidence, "not applicable")
      
      # change col names of data frame according to BIDS specification
      names(run_BIDS)[names(run_BIDS)=="vidFileName"] <- "stim_file"
      names(run_BIDS)[names(run_BIDS)=="variable"] <- "event"
      run_BIDS <-  run_BIDS[order(run_BIDS$onset),]
      run_BIDS$onset <- round(run_BIDS$onset, digits = 3)
      run_BIDS$duration <- round(run_BIDS$duration, digits = 3)
      
      events_BIDS <- run_BIDS[, c("onset", "duration", "trial", "stim_file", "event", "response", "response_timestamp", 
                                  "response_recall","trial_type_cuedRecallStrict", "trial_type_cuedRecallLenient",
                                  "response_recognition", "response_confidence",
                                  "trial_type_allConf", "trial_type_highConf", "trial_type_aboveAvgConf",
                                  "trial_type_rememberedStrictAboveAvg", "trial_type_rememberedLenientAboveAvg", "trial_type_rememberedStrictHigh", "trial_type_rememberedLenientHigh")]
      
      # save file for BIDS
      if (max(run$acq) > 1) { # include acq in filename if there was more than one acq in a run
        BIDSfilename <- paste0(BIDSstring, "_task-magictrickwatching_acq-", a, "_run-", b,"_events.tsv")
      } else {
        BIDSfilename <- paste0(BIDSstring, "_task-magictrickwatching_run-", b,"_events.tsv")
      }
      
      setwd(preprocessedEventsSubjDir)
      write.table(events_BIDS, file=BIDSfilename, quote = F, sep="\t", row.names = F, na = "n/a")
      #write.table(events_BIDS, file=file.path(preprocessedEventsRootDir, BIDSfilename), quote=FALSE, sep="\t", row.names = FALSE, na = "n/a")
      
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
    BIDS[(BIDS$event == "displayVid"), "responseCuriosity"] <- BIDS[(BIDS$event == "displayCuriosity"), "responseCuriosity"] # overwrite the "not applicable" with the actual curiosity rating
    BIDS_concat <- subset(BIDS, BIDS$event == "displayVid") # BIDS contains all onsets and durations for all 36 trials as well as the memory performance associated with all of them
    BIDS_concat$duration_mock <- BIDS_concat$mock-BIDS_concat$vid
    BIDS_concat$duration_vid <- as.numeric(BIDS_concat$displayVidOffset)-as.numeric(BIDS_concat$vid_per_run)
    BIDS_concat$duration_vid_withoutMock <-as.numeric(BIDS_concat$displayVidOffset)-as.numeric(BIDS_concat$mock_per_run)
    BIDS_concat$duration_vid_postFixation <-as.numeric(BIDS_concat$timestampPostVidFixation)-as.numeric(BIDS_concat$vid_per_run)
    BIDS_concat$duration_vid_withoutMock_postFixation <-as.numeric(BIDS_concat$timestampPostVidFixation)-as.numeric(BIDS_concat$mock_per_run)
    
    #names(BIDS_concat)
    BIDS_concat <- BIDS_concat[,c("vid", "mock", "duration_vid", "duration_vid_withoutMock", "duration_vid_withoutMock_postFixation", "avgVidDur_MAGMOT", "trial", "stim_file", "responseCuriosity", 
                                  "trial_type_cuedRecallStrict", "trial_type_cuedRecallLenient",
                                  "trial_type_allConf", "trial_type_highConf", "trial_type_aboveAvgConf",
                                  "trial_type_rememberedStrictAboveAvg", "trial_type_rememberedLenientAboveAvg", "trial_type_rememberedStrictHigh", "trial_type_rememberedLenientHigh", 
                                  "responseConfidence", "run", "acq")]
    
    if (feedback == "yes"){
      print(paste("total Duration is", sum(BIDS_concat$duration)))
    }
    #write.table(BIDS_concat, file=paste0(BIDSstring, "_task-magictrickwatching_concat.tsv"), quote=FALSE, sep="\t", row.names = FALSE, na = "n/a")
    write.table(BIDS_concat, file=file.path(concatRootDir, paste0(BIDSstring, "_task-magictrickwatching_concat.tsv")), quote=FALSE, sep="\t", row.names = FALSE, na = "n/a")
    
  }
  
  # delete no longer needed variables
  if (debug == 0){
    rm(duration, onset, run, run_BIDS, BIDS, events_BIDS)
  }
  
  
  ########### process questionnaire data collected during memory part   ########### 
  
  postMemory <-subset(memory, trial.type == "surveycat")   # select relevant rows and columns of memory task questions data
  postMemory <- postMemory[,c("username","startMemory", "startMemory_UTC", "endMemory", "endMemory_UTC", "durMemory", "memoryFile",
                              "survey_sleep_response","survey_sleep_hours",
                              "survey_test_known_response", "survey_memory_intention_response", "survey_reward_belief_response",
                              "survey_magictrick_experience_response", "survey_connection_response", "survey_comment_response")]
  names(postMemory) <- c("ID", "startMemory", "startMemory_UTC", "endMemory", "endMemory_UTC", "durMemory", "memoryFile",
                         "sleepBeforeMemoryTest","sleepHours", "memoryTestKnown", "memoryIntention", "rewardBelief", "magictrickExperience", "connection", "comment_memory")
  
  # recode rewardBelief
  postMemory$rewardBelief_score <- ifelse(postMemory$rewardBelief == "Not applicable", NA,
                                          ifelse(postMemory$rewardBelief == "Definitely agree ", 6,
                                                 ifelse(postMemory$rewardBelief == "Somehow agree", 5,
                                                        ifelse(postMemory$rewardBelief == "Slightly agree", 4,
                                                               ifelse(postMemory$rewardBelief == "Slightly disagree", 3,
                                                                      ifelse(postMemory$rewardBelief == "Somehow disagree", 2,
                                                                             ifelse(postMemory$rewardBelief == "Definitely disagree", 1, 0)))))))
  
  # add curiosity NAs to the data set in wide format
  postMemory$curiosityNAs <- curiosityNAs
  
  # add endExperiment and endPractice
  postMemory$endExperiment_raw <- ptbdata$endExperiment_raw[1]
  postMemory$endExperiment_LDN <- ptbdata$endExperiment_LDN[1]
  postMemory$endExperiment_UTC <- ptbdata$endExperiment_UTC[1]
  postMemory$endPractice_raw <- ptbdata$endPractice_raw[1]
  postMemory$endPractice_LDN <- ptbdata$endPractice_LDN[1]
  postMemory$endPractice_UTC <- ptbdata$endPractice_UTC[1]
  
  # calculate time span between experiment and memory test
  postMemory$durScanningSession <- difftime(postMemory$endExperiment_UTC, postMemory$endPractice_UTC, units = "mins")
  postMemory$durScanningSession <- ifelse( postMemory$durScanningSession > 180, NA,  postMemory$durScanningSession) # replace any values bigger than 3 hours with NA
  postMemory$daysBetweenExpAndMemory <- difftime(postMemory$startMemory_UTC, postMemory$endExperiment_UTC)
  
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
      
      # duration of each block
      postMemory[[paste0("durInSecs", blockstring[BLOCK])]] <- sum(data_subset$endBlock, na.rm = T) # using sum instead of max necessary because there were 2 acq for run2 in ppt 16
      postMemory[[paste0("durInMins", blockstring[BLOCK])]] <- sum(data_subset$endBlock, na.rm = T)/60
      
      # average curiosity
      postMemory[[paste0("responseCuriosity", blockstring[BLOCK])]] <- mean(data_subset$responseCuriosity, na.rm = T)
      postMemory[[paste0("curiosity", blockstring[BLOCK])]] <- mean(data_subset$curiosity, na.rm = T)
      
      # # sum up mean centered curiosity (continous) over trials 
      # postMemory[[paste0("highCuriosityTrials", blockstring[BLOCK])]] <- sum(data_subset$curiosityGroupMeanCentered > 0, na.rm = T) # sum up amount of tricks that are subject-wise high curiosity tricks
      # postMemory[[paste0("highCuriosityTrials_perc", blockstring[BLOCK])]] <- sum(data_subset$curiosityGroupMeanCentered > 0, na.rm = T) / dim(data_subset)[1] # percentage og high curiosity trials
      # postMemory[[paste0("lowCuriosityTrials", blockstring[BLOCK])]] <- sum(data_subset$curiosityGroupMeanCentered < 0, na.rm = T) # sum up amount of tricks that are subject-wise high curiosity tricks
      # postMemory[[paste0("lowCuriosityTrials_perc", blockstring[BLOCK])]] <- sum(data_subset$curiosityGroupMeanCentered < 0, na.rm = T) / dim(data_subset)[1] # percentage og high curiosity trials
      # 
      # # sum up mean centered curiosity (dichotomuous) over trials
      # postMemory[[paste0("highCuriosityTrials_dichotom", blockstring[BLOCK])]] <- sum(data_subset$curiosity_dich > 0, na.rm = T) # sum up amount of tricks that are subject-wise high curiosity tricks
      # postMemory[[paste0("highCuriosityTrials_dichotom_perc", blockstring[BLOCK])]] <- sum(data_subset$curiosity_dich > 0, na.rm = T) / dim(data_subset)[1] # percentage og high curiosity trials
      # postMemory[[paste0("lowCuriosityTrials_dichotom", blockstring[BLOCK])]] <- sum(data_subset$curiosity_dich < 0, na.rm = T) # sum up amount of tricks that are subject-wise high curiosity tricks
      # postMemory[[paste0("lowCuriosityTrials_dichotom_perc", blockstring[BLOCK])]] <- sum(data_subset$curiosity_dich < 0, na.rm = T) / dim(data_subset)[1] # percentage og high curiosity trials
      
      for (mem in 1:length(memoryLevels)) {
        # sum up scores for all memory levels
        postMemory[[paste0(memoryLabels[mem], "_abs", blockstring[BLOCK])]] <- sum(data_subset[[paste0(memoryLevels[mem])]], na.rm = T) 
        postMemory[[paste0(memoryLabels[mem], "_rel", blockstring[BLOCK])]] <- sum(data_subset[[paste0(memoryLevels[mem])]], na.rm = T)  / dim(data_subset)[1]
        # sum up curiosity-driven memory memory benefit (continuous, absolute)
        postMemory[[paste0("curBen_cont_", memoryLabels[mem], blockstring[BLOCK])]] <- sum(data_subset[[paste0("curBen_cont_", memoryLabels[mem])]], na.rm = T) 
        # sum up curiosity-driven memory memory benefit (dichotomuous, absolute)
        postMemory[[paste0("curBen_dich_", memoryLabels[mem], blockstring[BLOCK])]] <- sum(data_subset[[paste0("curBen_dich_", memoryLabels[mem])]], na.rm = T) 
        # calculate relative curiosity benefit
        postMemory[[paste0("curBen_rel_", memoryLabels[mem], blockstring[BLOCK])]] <- (sum(data_subset[[paste0("curBen_cont_", memoryLabels[mem])]] > 0, na.rm = T) / sum(data_subset$curiosityGroupMeanCentered > 0, na.rm = T)) - (sum(data_subset[[paste0("curBen_cont_", memoryLabels[mem])]] < 0, na.rm = T) / sum(data_subset$curiosityGroupMeanCentered < 0, na.rm = T))
        # calculate correlation between curiosity and memory
        postMemory[[paste0("curCor_",  memoryLabels[mem], blockstring[BLOCK])]] <- cor(data_subset$curiosityGroupMeanCentered, data_subset[[paste0(memoryLevels[mem])]], use = "pairwise.complete.obs") 
        
        # LM to predict memory using curiosity
        glm <- glm(data_subset[[paste0(memoryLevels[mem])]] ~ data_subset$curiosityGroupMeanCentered, family=binomial(link='logit'))
        postMemory[[paste0("lmBeta_", memoryLabels[mem])]] <- glm$coefficients[2]
      }
      
      # average mean confidence
      postMemory[[paste0("meanConfidence", blockstring[BLOCK])]]  <- mean(data_subset$responseConfidence, na.rm = T)
      temp_data <- subset(data_subset, data_subset$recognition == 1)
      postMemory[[paste0("meanConfidenceCorrectTrials", blockstring[BLOCK])]]   <- mean(temp_data$responseConfidence, na.rm = T)
      
      for (k in 1:6) { #confidence ranges from 1 to 6, potentially code can be made more flexible by using min(data$confidence) and max(data$confidence)
        temp_data <- subset(data_subset, data_subset$responseConfidence == k)
        postMemory[[paste0("recognitionConfLevel_", k, blockstring[BLOCK])]] <- sum(temp_data$recognition, na.rm = T)
        postMemory[[paste0("recognitionConfLevel_", k, "_perc", blockstring[BLOCK])]] <- sum(temp_data$recognition, na.rm = T) / dim(data_subset)[1]
        if (k < 6) {
          temp_data_above <- subset(data_subset, data_subset$responseConfidence > k)
          postMemory[[paste0("recognitionConfLevel_above_", k, blockstring[BLOCK])]] <-sum(temp_data_above$recognition, na.rm = T)
          postMemory[[paste0("recognitionConfLevel_above_", k, "_perc", blockstring[BLOCK])]] <-sum(temp_data_above$recognition, na.rm = T) / dim(data_subset)[1]
          rm(temp_data_above)
        }
        
        # sum up the scores for the recognition task pooled
        if (k == 1 || k == 3 || k == 5) {
          temp_data_plus1 <- subset(data_subset, data_subset$responseConfidence == (k+1))
          postMemory[[paste0("recognitionConfLevel_", k, "_", k+1, blockstring[BLOCK])]] <- sum(temp_data$recognition, na.rm = T) + sum(temp_data_plus1$recognition, na.rm = T)
          postMemory[[paste0("recognitionConfLevel_", k, "_", k+1, "_perc", blockstring[BLOCK])]] <- (sum(temp_data$recognition, na.rm = T) + sum(temp_data_plus1$recognition, na.rm = T))  / dim(data_subset)[1]
          rm(temp_data_plus1)
        }
        if (k == 1 || k == 4) {
          temp_data_plus1 <- subset(data_subset, data_subset$responseConfidence == (k+1))
          temp_data_plus2 <- subset(data_subset, data_subset$responseConfidence == (k+2))
          postMemory[[paste0("recognitionConfLevel_", k, "_", k+1, "_", k+2, blockstring[BLOCK])]] <- sum(temp_data$recognition, na.rm = T) + sum(temp_data_plus1$recognition, na.rm = T) + sum(temp_data_plus2$recognition, na.rm = T)
          postMemory[[paste0("recognitionConfLevel_", k, "_", k+1, "_", k+2, "_perc", blockstring[BLOCK])]] <- (sum(temp_data$recognition, na.rm = T) + sum(temp_data_plus1$recognition, na.rm = T) + sum(temp_data_plus2$recognition, na.rm = T))  / dim(data_subset)[1]
          rm(temp_data_plus1)
          rm(temp_data_plus2)
        }
        rm(temp_data)
      }
    } # end of loop over blockString
    
    if (debug == 0){
      rm(data_subset)
    }
  }
  
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
    
    ### long format data ###
    
    setwd(file.path(preprocessedDir, "long"))
    file_list <- list.files(pattern = "long.csv")
    # combine single _long.csv files to one file
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
    # save the final file
    setwd(preprocessedShareDir)
    xlsx::write.xlsx(dataLong, file=paste0("long_MagicBehavioural_", version_official, ".xlsx"), sheetName = "Sheet1", row.names = F) 
    write.csv(dataLong, file=paste0("long_MagicBehavioural_", version_official, ".csv"), row.names = FALSE, na = "NA")   
    
    ### wide format data ###
    
    # merge data in wide format
    MAGMOT <- merge(MAGMOT, questDataWide, by = c("ID", "group"), all = T)
    MAGMOT <- merge(MAGMOT, postMemoryWide, by = c("ID"), all = T)
    if (debug == 0){
      rm(postMemoryWide, questDataWide)
    }
    
    # compute the glmer models to extract the slopes
    for (mem in 1:length(memoryLevels)) {
      # model with curiosity, reward and their interaction
      LMEmodel <- glmer(dataLong[, memoryLevels[mem]] ~ groupEffectCoded*curiosityGroupMeanCentered + (1+curiosityGroupMeanCentered|ID) + (1|stimID), family = "binomial"(link = 'logit'), data = dataLong)
      MAGMOT[[paste0("curBeta_", memoryLabels[mem])]] <-  coef(LMEmodel)$ID$curiosityGroupMeanCentered
      MAGMOT[[paste0("curBeta_c_", memoryLabels[mem])]] <-  MAGMOT[[paste0("curBeta_", memoryLabels[mem])]] - mean(MAGMOT[[paste0("curBeta_", memoryLabels[mem])]])
      
      # model with curiosity only
      LMEmodel2 <- glmer(dataLong[, memoryLevels[mem]] ~ curiosityGroupMeanCentered + (1+curiosityGroupMeanCentered|ID) + (1|stimID), family = "binomial"(link = 'logit'), data = dataLong)
      MAGMOT[[paste0("curBetaOnly_", memoryLabels[mem])]] <-  coef(LMEmodel2)$ID$curiosityGroupMeanCentered
      MAGMOT[[paste0("curBetaOnly_c_", memoryLabels[mem])]] <-  MAGMOT[[paste0("curBetaOnly_", memoryLabels[mem])]] - mean(MAGMOT[[paste0("curBetaOnly_", memoryLabels[mem])]])
    }
    
    # add RSFC estimates between HPC & VTA (Pearson)
    setwd(brainDir)
    RSFC <- read.delim2("RSFC_VTA-HPC_pearson.txt") # read in data
    names(RSFC)[names(RSFC)=="RSFC_run.1"] <- "RSFC_VTAHPC_run1" # change col name
    names(RSFC)[names(RSFC)=="RSFC_run.2"] <- "RSFC_VTAHPC_run2" # change col name
    RSFC$RSFC_VTAHPC_run1 <- as.numeric(as.character(RSFC$RSFC_VTAHPC_run1)) # change from factor to numeric
    RSFC$RSFC_VTAHPC_run2 <- as.numeric(as.character(RSFC$RSFC_VTAHPC_run2)) # change from factor to numeric
    
    # fisher z transform correlations and compute RSFC change
    RSFC$RSFC_VTAHPC_run1_z <- atanh(RSFC$RSFC_VTAHPC_run1) # fisher z transform correlations
    RSFC$RSFC_VTAHPC_run2_z <- atanh(RSFC$RSFC_VTAHPC_run2) # fisher z transform correlations
    RSFC$RSFC_VTAHPC_diff <- RSFC$RSFC_VTAHPC_run2_z - RSFC$RSFC_VTAHPC_run1_z
    
    # merge wide format data
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
    
    # merge wide format data
    MAGMOT <- merge(MAGMOT, RSFC, by = "BIDS")
    
    # save data
    setwd(preprocessedShareDir)
    xlsx::write.xlsx(MAGMOT, file=paste0("wide_MagicBehavioural_", version_official, ".xlsx"), sheetName = "Sheet1", row.names = F) 
    write.csv(MAGMOT, file=paste0("wide_MagicBehavioural_", version_official, ".csv"), row.names = FALSE, na = "NA")
    
    ### create files that for the dataset paper ###
    # demographics
    demographics <- MAGMOT[,c("ID", "BIDS", "group", "age", "DOB", "sex", "gender", "ethnicity", "english", "ageEnglishAcquisition", 
                              "education", "yearsOfEducation", "employment", "studySubject", "handedness", "vision", "health")]
    write.csv(demographics, file=paste0( version, "_demographics", ".csv"), row.names = FALSE, na = "NA")
    
    # scores
    scores <- MAGMOT[,c("ID", "BIDS", 
                        # working memory
                        "corsiSpan", "nback_hits", "nback_misses_inclTooSlow", "nback_misses_exclTooSlow", 
                        "nback_correctrejections", "nback_falsealarms_inclTooSlow", "nback_falsealarms_exclTooSlow", 
                        "nback_hitrate", "nback_falsealarmrate", "nback_accurary",
                        # questionnaires
                        "BISBAS_inhibition", "BISBAS_rewardresponsiveness", "BISBAS_drive", "BISBAS_funseeking",
                        "NeedForCognition", "FearOfFailure", "ApproachTemperament", "AvoidanceTemperament", "TraitCuriosity", "StateCuriosity",
                        "IMI_intrinsicMotivation", "IMI_interest", "IMI_taskEngagement", "IMI_boredom", "IMI_pressure", "IMI_effort",
                        # memory performance
                        "cuedRecallStrict_abs", "cuedRecallStrict_rel", "cuedRecallLenient_abs", "cuedRecallLenient_rel", 
                        "allConf_abs", "allConf_rel",
                        "highConf_abs", "highConf_rel", "aboveAvgConf_abs", "aboveAvgConf_rel",
                        "rememberedStrictHigh_abs", "rememberedLenientHigh_rel", "rememberedLenientHigh_abs", "rememberedLenientHigh_rel",
                        "rememberedStrictAboveAvg_abs", "rememberedStrictAboveAvg_rel", "rememberedLenientAboveAvg_abs", "rememberedLenientAboveAvg_rel")]
    write.csv(scores, file=paste0( version, "_scores.csv"), row.names = FALSE, na = "NA")
    
    # raw_quest_data
    questionnaires <- c("inclusioncheck", "health_current", "health_ever", "screening_MRI", 
                        "BISBAS", "NeedForCognition", "FearOfFailure", "ApproachAndAvoidanceTemperament", "TraitCuriosity", 
                        "StateCuriosity", "PostExpAssessment")
    numItems <- c(7,4,24, 19, 20, 18, 9, 12, 18, 20, 24)
   
    for (q in seq_along(questionnaires)){ # determine item list
      for (ii in 1:numItems[q]){
        if (ii == 1){
          itemList <- paste0(questionnaires[q], ".", ii)
        } else {
          itemList <-  c(itemList, paste0(questionnaires[q], ".", ii, collapse = ", "))
        }
      }
      if (q == 1){
        itemListALL <- itemList
      } else {
        itemListALL <-  c(itemListALL, itemList)
      }
    }
    raw_quest_data <- questionnaire_raw[, c("ID", itemListALL)]
    write.csv(raw_quest_data, file=paste0( version, "_raw_quest_data.csv"), row.names = FALSE, na = "NA")
    
    # other information
    duration_info <- MAGMOT[,c(grepl("ID",names(MAGMOT)) | grepl("dur",names(MAGMOT)))]
    other_information <- MAGMOT[,c("ID", "BIDS", "ableToSee", "compliance", "sleepLastNight", "sleepAverage", "alcohol", "alcoholAmount", "rewardEffort", "rewardExpectations",
                                   "comment_task1", "comment_task2", "comment_task3", 
                                   "sleepBeforeMemoryTest", "sleepHours", "memoryTestKnown", "memoryIntention", "rewardBelief", "magictrickExperience", "connection", "comment_memory", "daysBetweenExpAndMemory")]
    other_information <- merge(other_information, duration_info, by = c("ID", "BIDS"))
    write.csv(other_information, file=paste0( version, "_other_information.csv"), row.names = FALSE, na = "NA")
    
    # experimental_data
    experimental_data <- dataLong[,c("ID", "BIDS", "group", "orderNumber", "block", "acq", "startBlock", "endBlock", 
                                     
                                     "timingCorrection", "jitterVideo_trial", "jitterRating_trial",
                                     
                                     "vidFileName", "trial", "tTrialStart", "tTrialEnd", "durationTrial", "fixationInitialDuration",
                                     "displayVidOnset", "displayVidOffset", "displayAnswerDuration", "displayBlankDuration", "fixationPostVidOnset", "fixationPostVidDuration",
                                     
                                     "displayAnswerOnset", "displayAnswerDuration",	"timeoutAnswer", "responseAnswer", "timestampAnswer", "timestampAnswerWhite", "rtAnswer", 
                                     "fixationPostAnswerOnset", "fixationPostAnswerDuration","betweenRatingFixation",
                                     "displayCuriosityOnset", "displayCuriosityDuration", "timeoutCuriosity", "responseCuriosity", "timestampCuriosity", "timestampCuriosityWhite",
                                     "rtCuriosity", "startValueCuriosity", "clicksCuriosity", "fixationPostCuriosityOnset", "fixationPostCuriosityDuration",
                                     "curiosity_tooSlow", # ifelse(ptbdata$timestampCuriosity > ptbdata$displayCuriosityOnset + ptbdata$timeoutCuriosity - 3*ptbdata$timingCorrection, 1, 0)
                                     "mockOffset", "cueImage", "momentOfSurprise_1", "momentOfSurprise_2", "momentOfSurprise_3", "momentOfSurprise_4", "momentOfSurprise_5", "momentOfSurprise_6", "momentOfSurprise_7",	
                                     "additionalMarker_momentOfSurprise_1", "additionalMarker_momentOfSurprise_2",
     
                                     "trialRecall", "responseRecall", "cuedRecallStrict", "cuedRecallLenient", "Flagging", "Comments", 
                                     "trialRecognition", "responseRecognition", "rtRecognition", 
                                     "responseConfidence", "rtConfidence",
                                     "recognition", "recognitionAboveMeanConf", "recognitionConfLevel_4_5_6", 
                                     "rememberedStrictAboveAvg", "rememberedLenientAboveAvg", "rememberedStrictHigh", "rememberedLenientHigh"
                                     )]
    write.csv(experimental_data, file=paste0( version, "_experimental_data.csv"), row.names = FALSE, na = "NA")
    
    
    ### upload file to OSF
    osfr::osf_auth() # log into OSF
    project <- osfr::osf_retrieve_node("fhqb7")
    target_dir <- osfr::osf_ls_files(project, pattern = "data") # looks at all files and directories in the project and defines the match with "data"
    sub_dir <- osfr::osf_mkdir(target_dir, path = paste0(version_official)) # add folder in OSF data dir
    # check whether file already exists - this is necessary due to a bug in the package
    file_exists <- osfr::osf_ls_files(sub_dir, pattern = "MagicBehavioural") # check whether file already exists
    while (dim(file_exists)[1] > 0){ #s delete files if they exists. use while loop because only the first row will be used
      osfr::osf_rm(file_exists, recurse = T, verbose = FALSE, check = F)
      file_exists <- osfr::osf_ls_files(sub_dir, pattern = "MagicBehavioural")
    }
    file_exists <- osfr::osf_ls_files(sub_dir, pattern = paste(version)) # check whether file already exists
    while (dim(file_exists)[1] > 0){ #s delete files if they exists. use while loop because only the first row will be used
      osfr::osf_rm(file_exists, recurse = T, verbose = FALSE, check = F)
      file_exists <- osfr::osf_ls_files(sub_dir, pattern = paste(version))
    }
    # upload all files in this directory
    osfr::osf_upload(sub_dir, path = ".", recurse = TRUE, conflicts = "overwrite")

    ### write information about scan durations
    names(scaninfoAll) <- c("ID", "scan", "duration_run_seconds", "duration_scan_seconds")
    scaninfoAll$duration_run_seconds <- round(scaninfoAll$duration_run_seconds, digits = 0)
    scaninfoAll$duration_scan_seconds <- round(scaninfoAll$duration_scan_seconds, digits = 0)
    scaninfoAll$duration_run_TR <- scaninfoAll$duration_run_seconds/TR
    scaninfoAll$duration_scan_TR <- scaninfoAll$duration_scan_seconds/TR
    setwd(preprocessedEventsRootDir)
    write.table(scaninfoAll, file="MAGMOT_informationAboutScanDuration.tsv", quote=FALSE, sep="\t", row.names = FALSE, na = "n/a")
    
    # upload file to OSF
    # check whether file already exists - this is necessary due to a bug in the package
    file_exists <- osfr::osf_ls_files(sub_dir, pattern = "MAGMOT_informationAboutScanDuration.tsv")
    if (dim(file_exists)[1] > 0){ # delete file if it exists
      osfr::osf_rm(file_exists, recurse = T, verbose = FALSE, check = F)
    }
    osfr::osf_upload(sub_dir, path = "MAGMOT_informationAboutScanDuration.tsv", conflicts = "overwrite") 

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
      subjectsToCorrelate <- subjectsCorr
      
      N <- 0.5*length(subjectsCorr)*length(subjectsToCorrelate)
      
      # create empty df for 3dISC -dataTable
      for (mem in seq_along(memoryLabels)){
        assign(paste0("dataTable_",memoryLabels[mem]), data.frame()) # if effect sizes for a single block were of interest  
      }
      
      dataTable_ISC <- data.frame()
      dataTable_ISC_dummy <- data.frame()
      dataTable <- data.frame()
      fillerRow <- data.frame()
      x <- 0
      
      if ( identical(subjectsToCorrelate, character(0)) == F){
        
        for(s in seq_along(subjectsCorr)){
          
          #ss = s+1
          subjectsToCorrelate <-subjectsCorr[-c(1:s)]
          
          # read in concat for subjectsCorr[s]
          # concatSubjDir_s <- file.path(concatRootDir, subjectsCorr[s])
          # setwd(preprocessedEventsSubjDir_s)
          setwd(concatRootDir)
          #file_s <- list.files(pattern = "concat")
          file_s <- list.files(pattern = paste0(subjectsCorr[s], "_task-magictrickwatching_concat.tsv"))
          
          events_s <- read.delim(file = file_s, header = T, sep="\t", na = "n/a")
          
          # pick relevant columns
          events_s <- events_s[,c("vid", "mock", "duration_vid", "duration_vid_withoutMock", "avgVidDur_MAGMOT", "stim_file", "responseCuriosity", 
                                  "trial_type_cuedRecallStrict", "trial_type_cuedRecallLenient",
                                  "trial_type_allConf", "trial_type_highConf", "trial_type_aboveAvgConf", 
                                  "trial_type_rememberedStrictAboveAvg", "trial_type_rememberedLenientAboveAvg", "trial_type_rememberedStrictHigh", "trial_type_rememberedLenientHigh", 
                                  "responseConfidence")]
          # mean-center curiosity and confidence
          events_s$responseCuriosity <- events_s$responseCuriosity - mean(events_s$responseCuriosity, na.rm = T) #mean center curiosity
          events_s$responseConfidence <- events_s$responseConfidence - mean(events_s$responseConfidence, na.rm = T) #mean center confidence
          
          # change column names so that they include the subject ID
          names(events_s)[names(events_s)=="vid"] <- "onset"
          names(events_s) <- paste0(names(events_s), "_", subjectsCorr[s])
          names(events_s)[names(events_s)== paste0("stim_file_", subjectsCorr[s])] <- "stim_file"
          names(events_s)[names(events_s)== paste0("avgVidDur_MAGMOT_", subjectsCorr[s])] <- "avgVidDur_MAGMOT"
          
          for(ss in seq_along(subjectsToCorrelate)){
            
            # read in concat for subjectsToCorrelate[ss]
            # concatSubjDir_s <- file.path(concatRootDir, subjectsCorr[s])
            # setwd(preprocessedEventsSubjDir_s)
            setwd(concatRootDir)
            #file_s <- list.files(pattern = "concat")
            file_ss <- list.files(pattern = paste0(subjectsToCorrelate[ss], "_task-magictrickwatching_concat.tsv"))
            
            # preprocessedEventsSubjDir_ss <- file.path(preprocessedEventsRootDir, paste(subjectsToCorrelate[ss]))
            # setwd(preprocessedEventsSubjDir_ss)
            # file_ss <- list.files(pattern = "concat")
            
            events_ss <- read.delim(file = file_ss, header = T, sep="\t", na = "n/a")
            
            # pick relevant columns
            events_ss <- events_ss[,c("vid", "mock", "duration_vid", "duration_vid_withoutMock", "avgVidDur_MAGMOT", "stim_file", "responseCuriosity", 
                                      "trial_type_cuedRecallStrict", "trial_type_cuedRecallLenient",
                                      "trial_type_allConf", "trial_type_highConf", "trial_type_aboveAvgConf", 
                                      "trial_type_rememberedStrictAboveAvg", "trial_type_rememberedLenientAboveAvg", "trial_type_rememberedStrictHigh", "trial_type_rememberedLenientHigh", 
                                      "responseConfidence")]
            # mean-center curiosity and confidence
            events_ss$responseCuriosity <- events_ss$responseCuriosity - mean(events_ss$responseCuriosity, na.rm = T) #mean center curiosity
            events_ss$responseConfidence <- events_ss$responseConfidence - mean(events_ss$responseConfidence, na.rm = T) #mean center confidence
            
            # change column names so that they include the subject ID
            names(events_ss)[names(events_ss)=="vid"] <- "onset"
            names(events_ss) <- paste0(names(events_ss), "_", subjectsToCorrelate[ss])
            names(events_ss)[names(events_ss)== paste0("stim_file_", subjectsToCorrelate[ss])] <- "stim_file"
            names(events_ss)[names(events_ss)== paste0("avgVidDur_MAGMOT_", subjectsToCorrelate[ss])] <- "avgVidDur_MAGMOT"
            
            # marge both subjects
            events <- merge(events_s, events_ss, by = c("stim_file", "avgVidDur_MAGMOT"))
            
            for (mem in 1:length(memoryLevels)) {
              # behavioural parcellation: define match in memory performance
              events[[paste0("behavParcel_", memoryLabels[mem])]] <- ifelse(events[[paste0("trial_type_",  memoryLabels[mem], "_",subjectsCorr[s])]] == "remembered" & events[[paste0("trial_type_", memoryLabels[mem], "_", subjectsToCorrelate[ss])]] == "remembered", "bothRemembered",
                                                                            ifelse(events[[paste0("trial_type_",  memoryLabels[mem], "_",subjectsCorr[s])]] == "forgotten" & events[[paste0("trial_type_", memoryLabels[mem], "_", subjectsToCorrelate[ss])]] == "forgotten", "bothForgotten",
                                                                                   "differentResponses"))
              # create an effect coded memory variable
              events[[paste0("trial_type_", memoryLabels[mem],"_",subjectsCorr[s], "_effectCoded")]] <- ifelse(events[[paste0("trial_type_",  memoryLabels[mem], "_",subjectsCorr[s])]] == "remembered", 1, -1)
              events[[paste0("trial_type_", memoryLabels[mem],"_",subjectsToCorrelate[ss], "_effectCoded")]] <- ifelse(events[[paste0("trial_type_",  memoryLabels[mem], "_",subjectsToCorrelate[ss])]] == "remembered", 1, -1)
            }
            
            # create effect coded group variable
            events[[paste0("reward_",subjectsCorr[s], "_effectCoded")]] <- ifelse(grepl("cont", subjectsCorr[s]), -1, 1)
            events[[paste0("reward_",subjectsToCorrelate[ss], "_effectCoded")]] <- ifelse(grepl("cont", subjectsToCorrelate[ss]), -1, 1)
            
            ### here we start a loop for each member of the pair
            pair <- c(subjectsCorr[s], subjectsToCorrelate[ss])
            # create an index variable
            x <- x+1
            
            ### fill in information to dataTable (group coding) ###
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
            
            ### fill in information to dataTable_ISC_dummy (group coding) ###
            dataTable_ISC_dummy[x,1] <- subjectsCorr[s]
            dataTable_ISC_dummy[x,2] <- subjectsToCorrelate[ss]
            #we adopt deviation coding for the two groups by replacing two groups G1 (cont) and G2 (exp) with -0.5 and 0.5. 
            #Then add up the two values for each row (each subject pair), resulting in three possible values of 1, -1 and 0.
            dataTable_ISC_dummy[x,3] <- ifelse(grepl("cont", subjectsCorr[s]) & grepl("cont", subjectsToCorrelate[ss]), -1, 
                                               ifelse(grepl("cont", subjectsCorr[s]) & grepl("exp", subjectsToCorrelate[ss]), 0, 
                                                      ifelse(grepl("exp", subjectsCorr[s]) & grepl("exp", subjectsToCorrelate[ss]), 1,NA )))
            
            #curiosity and curiosity interaction
            dataTable_ISC_dummy[x,4] <- cor(events[, paste0("responseCuriosity_", subjectsCorr[s])], events[, paste0("responseCuriosity_", subjectsToCorrelate[ss])] )
            dataTable_ISC_dummy[x,5] <- dataTable_ISC_dummy[x,4] * dataTable_ISC_dummy[x,3] 
            dataTable_ISC_dummy[x,6] <- cor(events[, paste0("responseConfidence_", subjectsCorr[s])], events[, paste0("responseConfidence_", subjectsToCorrelate[ss])] )
            dataTable_ISC_dummy[x,7] <- dataTable_ISC_dummy[x,6] * dataTable_ISC_dummy[x,3] 
            
            addCol <- 0
            for (mem in 1:length(memoryLevels)) {
              # memory and reward-memory interaction
              memoCol <- 8 + (10*addCol) # add more columns
              
              dataTable_ISC_dummy[x,memoCol] <- cor(events[, paste0("trial_type_", memoryLabels[mem],"_",subjectsCorr[s], "_effectCoded")], events[, paste0("trial_type_", memoryLabels[mem],"_",subjectsToCorrelate[ss], "_effectCoded")])
              dataTable_ISC_dummy[x,memoCol+1] <-  dataTable_ISC_dummy[x,memoCol] * dataTable_ISC_dummy[x,3] # third col has group information
              
              # curiosity beta and curiosity beta interaction
              betaCol <- memoCol + 2
              dataTable_ISC_dummy[x,betaCol] <- (MAGMOT[[paste0("curBeta_c_", memoryLabels[mem])]][MAGMOT$BIDS == subjectsCorr[s]] + MAGMOT[[paste0("curBeta_c_", memoryLabels[mem])]][MAGMOT$BIDS == subjectsToCorrelate[ss]])/2
              dataTable_ISC_dummy[x,betaCol+1] <- dataTable_ISC_dummy[x,betaCol] * dataTable_ISC_dummy[x,3] 
              
              # curiosity only beta and curiosity only beta interaction
              betaOnlyCol <- betaCol + 2
              
              dataTable_ISC_dummy[x,betaOnlyCol] <- (MAGMOT[[paste0("curBetaOnly_c_", memoryLabels[mem])]][MAGMOT$BIDS == subjectsCorr[s]] + MAGMOT[[paste0("curBetaOnly_c_", memoryLabels[mem])]][MAGMOT$BIDS == subjectsToCorrelate[ss]])/2
              dataTable_ISC_dummy[x,betaOnlyCol+1] <- dataTable_ISC_dummy[x,betaOnlyCol] * dataTable_ISC_dummy[x,3] 
              
              # curiosity benefit and curiosity benefit interaction (RELATIVE)
              benCol <- betaOnlyCol + 2
              
              dataTable_ISC_dummy[x,benCol] <- (MAGMOT[[paste0("curBen_rel_", memoryLabels[mem])]][MAGMOT$BIDS == subjectsCorr[s]] + MAGMOT[[paste0("curBen_rel_", memoryLabels[mem])]][MAGMOT$BIDS == subjectsToCorrelate[ss]])/2
              dataTable_ISC_dummy[x,benCol+1] <- dataTable_ISC_dummy[x,benCol] * dataTable_ISC_dummy[x,3]              
              
              # curiosity correlation and curiosity correlation interaction (RELATIVE)
              corCol <- benCol + 2
              
              dataTable_ISC_dummy[x,corCol] <- (MAGMOT[[paste0("curCor_", memoryLabels[mem])]][MAGMOT$BIDS == subjectsCorr[s]] + MAGMOT[[paste0("curCor_", memoryLabels[mem])]][MAGMOT$BIDS == subjectsToCorrelate[ss]])/2
              dataTable_ISC_dummy[x,corCol+1] <- dataTable_ISC_dummy[x,corCol] * dataTable_ISC_dummy[x,3]                
              
              addCol <- addCol + 1 # updates multiplier
            }
            
            # add file name and back slashs
            nextCol <- corCol + 2 # looks at current number of columns in object and adds 1
            
            dataTable_ISC_dummy[x,nextCol] <- paste0("ISC_",subjectsCorr[s],subjectsToCorrelate[ss],"_magictrickwatching_z.nii.gz")
            if(x < N) {
              dataTable_ISC_dummy[x,nextCol+1] <- '\\' #add back slash at end of the row
            }
            
            ### create SME dataTables ###
            for (p in seq_along(pair)) {
              
              # disentangle events for pair[p]
              SME_events <- events[,grep(pair[p], colnames(events))] #picks all columns relating to one of the subjects in the pair
              colToDelete <- grep("*effectCoded", colnames(SME_events)) # define columns with effect coded performance
              SME_events <- SME_events[, - colToDelete] # delete columns with effect coded performance
              # add file information
              SME_events$stim_file <- events$stim_file
              SME_events$avgVidDur_MAGMOT <- events$avgVidDur_MAGMOT
              # add behavioural parcellation
              behavParc <- events[,grep("behavParcel", colnames(events))] #picks all columns relating to behavioural parcellation
              names(behavParc) <- paste0(names(behavParc), "_",subjectsCorr[s],subjectsToCorrelate[ss])
              SME_events <- merge(SME_events, behavParc, by = "row.names")
              
              xx <- 2
              # loop over all memory labels to create SME tables
              for (mem in seq_along(memoryLabels)){
                
                # reset variable for if statement
                anyRemembered <- "no"
                anyForgotten <- "no"
                
                # for remembered and forgotten
                for (o in seq_along(SME_outcome)){
                  
                  # subset data depnding on memory performance
                  SME_events_outcome <-  subset(SME_events, SME_events[[paste0("behavParcel_", memoryLabels[mem], "_", subjectsCorr[s],subjectsToCorrelate[ss])]] == paste0(SME_outcome[o]))
                  if (feedback == "yes"){
                    print(paste("behavParcel", memoryLabels[mem], subjectsCorr[s],subjectsToCorrelate[ss], "has", dim(SME_events_outcome)[1], "rows"))
                  }
                  
                  # save data for each of the subjects
                  if(dim(SME_events_outcome)[1]>0){
                    concatPairDir <- file.path(concatRootDir, paste0(subjectsCorr[s],subjectsToCorrelate[ss]))
                    ifelse(!dir.exists(concatPairDir), dir.create(concatPairDir), FALSE)
                    setwd(concatPairDir)
                    write.table(SME_events_outcome, file = paste0(subjectsCorr[s], subjectsToCorrelate[ss], "_", pair[p], "_", SME_outcome[o], "_", memoryLabels[mem], "_SME_concat.tsv"), quote=FALSE, sep="\t", row.names = FALSE, na = "n/a")
                    
                    # for the first subject in each for
                    if (p == 1){
                      if(SME_outcome[o] == "bothRemembered"){
                        
                        anyRemembered <- "yes" # this variable needs to be created to use in if statement
                        dataTable[x,xx+o] <- dim(SME_events_outcome)[1]
                        
                      } else if(SME_outcome[o] == "bothForgotten"){ # for allConf we need to look at bothForgotten as some pairs don't have any magictricks that they both have forgotten
                        
                        anyForgotten <- "yes" # this variable needs to be created to use in if statement
                        dataTable[x,xx+o] <- dim(SME_events_outcome)[1]
                        xx <- xx+2
                      } 
                      if (anyRemembered == "yes" && anyForgotten == "yes") {
                        # create filler table
                        fillerRow[1,1] <- subjectsCorr[s]
                        fillerRow[1,2] <- subjectsToCorrelate[ss] 
                        fillerRow[1,3] <- ifelse(grepl("cont", subjectsCorr[s]) & grepl("cont", subjectsToCorrelate[ss]), -1, 
                                                 ifelse(grepl("cont", subjectsCorr[s]) & grepl("exp", subjectsToCorrelate[ss]), 0, 
                                                        ifelse(grepl("exp", subjectsCorr[s]) & grepl("exp", subjectsToCorrelate[ss]), 1,NA )))
                        fillerRow[1,4] <- paste0("ISC_",subjectsCorr[s],subjectsToCorrelate[ss],"_SME_", memoryLabels[mem], ".nii.gz")
                        fillerRow[1,5] <- '\\' #add back slash at end of the row
                        
                        # rbind fillerRow to dataTable
                        if(s == 1 && ss == 1){
                          assign(paste0("dataTable_", memoryLabels[mem]), fillerRow)
                        } else if(s < 50){
                          fillerTable <- rbind(get(paste0("dataTable_",memoryLabels[mem])), fillerRow) # if effect sizes for a single block were of interest
                          assign(paste0("dataTable_", memoryLabels[mem]), fillerTable) 
                        }
                        
                        # reset variable for if statement
                        anyRemembered <- "no"
                        anyForgotten <- "no"
                        # clear filler row
                        fillerRow <- data.frame()
                        
                      }
                    }
                    
                  } else { # end of "if(dim(SME_events_outcome)[1]>0)"
                    print(paste0("pair ",subjectsCorr[s],subjectsToCorrelate[ss], " does not have any ", SME_outcome[o], "_", memoryLabels[mem], " magictricks." ))
                  }
                }  # end of "for (o in seq_along(SME_outcome))"
              } # end of "for (mem in seq_along(memoryLabels))"
            } # end of "for (p in seq_along(pair))"
          } # end of subjectsToCorrelate
        } # end of subjectsCorr
      } 
      
      ############# once all subject pairs are processed ############# 
      
      ### process dataTable_ISC ###
      dataTable_ISC[nrow(dataTable_ISC),ncol(dataTable_ISC)] <- NA # no // for last column
      dataTable_ISC_dummy[nrow(dataTable_ISC_dummy),ncol(dataTable_ISC_dummy)] <- NA # no // for last column
      rm(fillerTable)
      for (mem in seq_along(memoryLabels)){
        fillerTable <- get(paste0("dataTable_",memoryLabels[mem]))
        fillerTable[nrow(fillerTable),ncol(fillerTable)] <- NA # no // for last column
        names(fillerTable) <- c("Subj1", "Subj2", "grp", "InputFile", "\\")
        assign(paste0("dataTable_",memoryLabels[mem]), fillerTable) 
        
        if (mem == 1){
          dataTablenames <- c(paste0(memoryLabels[mem], "_remembered"))
        } else {
          dataTablenames <-  c(dataTablenames, paste0(memoryLabels[mem], "_remembered", collapse = ", "))
        }
        dataTablenames <-  c(dataTablenames, paste0(memoryLabels[mem], "_forgotten", collapse = ", "))
      }
      
      names(dataTable_ISC) <- c("Subj1", "Subj2", "grp", "InputFile", "\\")
      names(dataTable) <- c("Subj1", "Subj2", dataTablenames)
      
      ### process dataTable_ISC_dummy ###
      # create vector with all column names for ISC table
      for (mem in 1:length(memoryLevels)) {
        if (mem == 1){
          ISC_table_names <- c(paste0("corr_", memoryLabels[mem]))
        } else {
          ISC_table_names <-  c(ISC_table_names, paste0("corr_", memoryLabels[mem], collapse = ", "))
        }
        ISC_table_names <-  c(ISC_table_names, paste0("grCorr_", memoryLabels[mem], collapse = ", "))
        ISC_table_names <-  c(ISC_table_names, paste0("curBeta_", memoryLabels[mem], collapse = ", "))
        ISC_table_names <-  c(ISC_table_names, paste0("grCurBeta_", memoryLabels[mem], collapse = ", "))
        ISC_table_names <-  c(ISC_table_names, paste0("curBetaOnly_", memoryLabels[mem], collapse = ", "))
        ISC_table_names <-  c(ISC_table_names, paste0("grCurBetaOnly_", memoryLabels[mem], collapse = ", "))
        ISC_table_names <-  c(ISC_table_names, paste0("curBen_", memoryLabels[mem], collapse = ", "))
        ISC_table_names <-  c(ISC_table_names, paste0("grCurBen_", memoryLabels[mem], collapse = ", "))
        ISC_table_names <-  c(ISC_table_names, paste0("curCor_", memoryLabels[mem], collapse = ", "))
        ISC_table_names <-  c(ISC_table_names, paste0("grCurCor_", memoryLabels[mem], collapse = ", "))
      }
      # round values
      for (cc in (seq_along(colnames(dataTable_ISC_dummy)))){
        currentCol <- colnames(dataTable_ISC_dummy)[cc]
        if(is.numeric(dataTable_ISC_dummy[,names(dataTable_ISC_dummy)==currentCol]) == T){
          dataTable_ISC_dummy[names(dataTable_ISC_dummy)==currentCol] <- round(dataTable_ISC_dummy[names(dataTable_ISC_dummy)==currentCol], digits = 5)
        }
      }
      # add column names to ISC table dummy
      names(dataTable_ISC_dummy) <- c("Subj1", "Subj2", "grp", "corr_curiosity", "grCorr_curiosity", "corr_confidence", "grCorr_confidence",
                                      ISC_table_names, "InputFile", "\\")
      
      ### create data table with unique effects of curiosity, memory and their interaction ###
      dataTable_ISC_controlled <- dataTable_ISC_dummy[, c("Subj1", "Subj2", "grp")] # copy some columns from dataTable_ISC_dummy
      
      # unique confidence
      uniqueConfidence <- lm(dataTable_ISC_dummy$corr_confidence ~ dataTable_ISC_dummy$corr_allConf)
      dataTable_ISC_controlled[,4] <- round(residuals(uniqueConfidence), digits = 5) # residuals from models
      dataTable_ISC_controlled[,5] <- dataTable_ISC_controlled[,4] * dataTable_ISC_controlled[,3] # interaction
      
      addCol2 <- 0
      for (mem in 1:length(memoryLevels)) {
        # unique curiosity
        curCol <- 6 + (12*addCol2) # add more columns
        uniqueCuriosity <- lm(dataTable_ISC_dummy$corr_curiosity ~ dataTable_ISC_dummy[[paste0("corr_", memoryLabels[mem])]])
        dataTable_ISC_controlled[,curCol] <- round(residuals(uniqueCuriosity), digits = 5) # residuals from models
        dataTable_ISC_controlled[,curCol+1] <- dataTable_ISC_controlled[,curCol] * dataTable_ISC_controlled[,3] # interaction
        
        # unique memory
        memCol <- curCol + 2 # add more columns
        uniqueMemory <- lm(dataTable_ISC_dummy[[paste0("corr_", memoryLabels[mem])]] ~ dataTable_ISC_dummy$corr_curiosity)
        dataTable_ISC_controlled[,memCol] <- round(residuals(uniqueMemory), digits = 5) # residuals from models
        dataTable_ISC_controlled[,memCol+1] <- dataTable_ISC_controlled[,memCol] * dataTable_ISC_controlled[,3] # interaction
        
        # unique curiosity beta
        uBetaCol <- memCol + 2 # add more columns
        uniqueBeta <- lm(dataTable_ISC_dummy[[paste0("curBeta_", memoryLabels[mem])]]  ~ dataTable_ISC_dummy[[paste0("corr_", memoryLabels[mem])]] + dataTable_ISC_dummy$corr_curiosity)
        dataTable_ISC_controlled[,uBetaCol] <- round(residuals(uniqueBeta), digits = 5) # residuals from models
        dataTable_ISC_controlled[,uBetaCol+1] <- dataTable_ISC_controlled[,uBetaCol] * dataTable_ISC_controlled[,3] # interaction
        
        # unique curiosity beta only
        uBetaOnlyCol <- uBetaCol + 2 # add more columns
        uniqueBetaOnly <- lm(dataTable_ISC_dummy[[paste0("curBetaOnly_", memoryLabels[mem])]]  ~ dataTable_ISC_dummy[[paste0("corr_", memoryLabels[mem])]] + dataTable_ISC_dummy$corr_curiosity)
        dataTable_ISC_controlled[,uBetaOnlyCol] <- round(residuals(uniqueBetaOnly), digits = 5) # residuals from models
        dataTable_ISC_controlled[,uBetaOnlyCol+1] <- dataTable_ISC_controlled[,uBetaOnlyCol] * dataTable_ISC_controlled[,3] # interaction
        
        # unique curiosity benefit
        uBenCol <- uBetaOnlyCol + 2 # add more columns
        uniqueBenefit <- lm(dataTable_ISC_dummy[[paste0("curBen_", memoryLabels[mem])]]  ~ dataTable_ISC_dummy[[paste0("corr_", memoryLabels[mem])]] + dataTable_ISC_dummy$corr_curiosity)
        dataTable_ISC_controlled[,uBenCol] <- round(residuals(uniqueBenefit), digits = 5) # residuals from models
        dataTable_ISC_controlled[,uBenCol+1] <- dataTable_ISC_controlled[,uBenCol] * dataTable_ISC_controlled[,3] # interaction
        
        # unique curiosity correlation
        uCorCol <- uBenCol + 2 # add more columns
        uniqueCorrelation <- lm(dataTable_ISC_dummy[[paste0("curCor_", memoryLabels[mem])]]  ~ dataTable_ISC_dummy[[paste0("corr_", memoryLabels[mem])]] + dataTable_ISC_dummy$corr_curiosity)
        dataTable_ISC_controlled[,uCorCol] <- round(residuals(uniqueCorrelation), digits = 5) # residuals from models
        dataTable_ISC_controlled[,uCorCol+1] <- dataTable_ISC_controlled[,uCorCol] * dataTable_ISC_controlled[,3] # interaction
        
        addCol2 <- addCol2 + 1 # updates multiplier
      }
      
      # add file name and back slashs
      nextCol2 <- uCorCol + 2 # looks at current number of columns in object and adds 1
      dataTable_ISC_controlled[,nextCol2] <- dataTable_ISC_dummy[, "InputFile"]
      dataTable_ISC_controlled[,nextCol2+1] <- '\\' #add back slash at end of the row
      dataTable_ISC_controlled[nrow(dataTable_ISC_controlled),ncol(dataTable_ISC_controlled)] <- NA # no // for last column
      
      for (mem in 1:length(memoryLevels)) {
        if (mem == 1){
          ISC_table_names <- c(paste0("uniqueCur_", memoryLabels[mem]))
        } else {
          ISC_table_names <-  c(ISC_table_names, paste0("uniqueCur_", memoryLabels[mem], collapse = ", "))
        }
        ISC_table_names <-  c(ISC_table_names, paste0("grUniqueCur_", memoryLabels[mem], collapse = ", "))
        ISC_table_names <-  c(ISC_table_names, paste0("uniqueMem_", memoryLabels[mem], collapse = ", "))
        ISC_table_names <-  c(ISC_table_names, paste0("grUniqueMem_", memoryLabels[mem], collapse = ", "))
        ISC_table_names <-  c(ISC_table_names, paste0("uniqueBeta_", memoryLabels[mem], collapse = ", "))
        ISC_table_names <-  c(ISC_table_names, paste0("grUniqueBeta_", memoryLabels[mem], collapse = ", "))
        ISC_table_names <-  c(ISC_table_names, paste0("uniqueBetaOnly_", memoryLabels[mem], collapse = ", "))
        ISC_table_names <-  c(ISC_table_names, paste0("grUniqueBetaOnly_", memoryLabels[mem], collapse = ", "))
        ISC_table_names <-  c(ISC_table_names, paste0("uniqueBen_", memoryLabels[mem], collapse = ", "))
        ISC_table_names <-  c(ISC_table_names, paste0("grUniqueBen_", memoryLabels[mem], collapse = ", "))
        ISC_table_names <-  c(ISC_table_names, paste0("uniqueCor_", memoryLabels[mem], collapse = ", "))
        ISC_table_names <-  c(ISC_table_names, paste0("grUniqueCor_", memoryLabels[mem], collapse = ", "))
      }
      # add column names to dataTable_ISC_controlled
      names(dataTable_ISC_controlled) <- c("Subj1", "Subj2", "grp", "uniqueConf", "grpUniqueConf",
                                           ISC_table_names, "InputFile", "\\")
      
      ### save all dataTables ###
      setwd(dataTableDir)
      write.table(dataTable_ISC, file="dataTable_magictrickwatching.txt", quote=FALSE, sep="\t", row.names = FALSE, na = "")
      write.table(dataTable_ISC_dummy, file="dataTable_magictrickwatching_memo.txt", quote=FALSE, sep="\t", row.names = FALSE, na = "")
      write.table(dataTable_ISC_controlled, file="dataTable_magictrickwatching_unique.txt", quote=FALSE, sep="\t", row.names = FALSE, na = "")
      write.table(dataTable, file="dataTable_pairwise_memoryScore.csv", quote=FALSE, sep=",", row.names = FALSE, na = "NA")
      
      for (mem in seq_along(memoryLabels)){
        write.table(get(paste0("dataTable_", memoryLabels[mem])), file=paste0("dataTable_",   memoryLabels[mem], ".txt"), quote=FALSE, sep="\t", row.names = FALSE, na = "")
      }
      
    } # end doISCprep
  }
}

