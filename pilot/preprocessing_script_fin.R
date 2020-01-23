# this script preprocesses all single participant files saved by collector. 
# it first reads in the data from the main task, saves some of it in wide format (i.e. questionnaire data, demographic input, etc), 
# wheres it also creates a long format with the information for each magic trick seperately sorted in a row.
# It also handles the memory test data and adds it to the created wide and long format if such data exists.
# Final outputs are: 
#   for each participant data in wide and long format
#   for the whole sample data in wide and long format

#### setups ####

#empty work space, load libraries and functions
rm(list=ls())
library(xlsx)
source("~/Dropbox/Reading/Codes and functions/R/rbindcolumns.R")

# define necessary directories
mainDir <- "~/Dropbox/Reading/PhD/Magictricks/behavioural_study"
subDirData <- "data_fin"
dataDir <- file.path(mainDir, subDirData) #"~/Dropbox/Reading/PhD/Magic tricks/behavioural_study/data_fin"
groupDir <- file.path(dataDir, "MagicBehavioural_") #"~/Dropbox/Reading/PhD/Magic tricks/behavioural_study/data_fin/MagicBehavioural_"
contDir <- file.path(dataDir, "MagicBehavioural_cont") #"~/Dropbox/Reading/PhD/Magic tricks/behavioural_study/data_fin/MagicBehavioural_cont"
expDir <- file.path(dataDir, "MagicBehavioural_exp") #"~/Dropbox/Reading/PhD/Magic tricks/behavioural_study/data_fin/MagicBehavioural_exp"
memoryDir <- file.path(dataDir, "MagicBehavioural_memory") #"~/Dropbox/Reading/PhD/Magic tricks/behavioural_study/data_fin/MagicBehavioural_memory"
preprocessedDir <- file.path(dataDir, "MagicBehavioural_preprocessed") #"~/Dropbox/Reading/PhD/Magic tricks/behavioural_study/data_fin/MagicBehavioural_preprocessed"
codedDir <- file.path(dataDir, "coded", "preprocessing") #"~/Dropbox/Reading/PhD/Magic tricks/behavioural_study/data_fin/coded/preprocessing"

# check whether these directories exist, if not create them
ifelse(!dir.exists(dataDir), dir.create(dataDir), FALSE)
ifelse(!dir.exists(contDir), dir.create(contDir), FALSE)
ifelse(!dir.exists(expDir), dir.create(expDir), FALSE)
ifelse(!dir.exists(memoryDir), dir.create(memoryDir), FALSE)
ifelse(!dir.exists(preprocessedDir), dir.create(preprocessedDir), FALSE) 

# read in the file that contains information about the stimuli used
setwd(file.path(mainDir, "stimuli"))
info_old <- read.xlsx("magic_selection.xlsx", sheetName = "preselection", showWarnings = FALSE) # recomment this file!
info_old <- subset(info_old, info_old$stimID != "H4_long" ) #practice trial
info_old <- subset(info_old, info_old$stimID != "K23") #practice trial
info_old <- info_old[,c("stimID", "length")]
info <- read.xlsx("stimuli_MagicBehavioural_fin_2019-03-27.xlsx", sheetName = "Sheet1") # mistake in code detected on 2019-03-27, previously 2018-12-19
info <- info[,c("stimID", "meanCuriosityStandardisedAya", "mediansplitCuriosityAya", "meanCuriosity", "meanCuriosityStandardisedSample", "mediansplitCuriositySample")]
names(info) <- c("stimID", "meanCuriosityStandardisedAya", "mediansplitCuriosityAya", "meanCuriositySample", "meanCuriosityStandardisedSample", "mediansplitCuriositySample")

info <- merge(info, info_old, by = "stimID")

# define experimental groups
group <- c("cont", "exp")

# define version 
version <- "fin"

# define block names
blockstring <- c("_firstBlock", "_secondBlock", "_thirdBlock", "")

##### LOOP to preprocess all data files #####

# start to loop over groups
for (l in seq_along(group)) { 
  
  # define subject list from prolific export file
  setwd(dataDir)
  report <- read.csv(paste0("prolific_export_MagicBehavioural_", group[l], "_", version, ".csv"), header = T)
  report <- subset(report, report$status == "AWAITING REVIEW" | report$status == "APPROVED") #NOTE: CHANGE!!
  
  subjects <-  report$participant_id # we are reading in data of all subjects and have fillers for the exp data + output that data does not exist
  subjects <- as.character(subjects)

  # start loop over subjects
  for (s in seq_along(subjects)){
    
    # read in responses.csv files for video watching part
    # if.exists() statement to check whether a file is available or not
    setwd(paste0(groupDir, group[l]))
    dataset <- paste0("MagicBehavioural_", group[l], "_", version, "_", subjects[s], "_responses.csv")
    if (file.exists(dataset)){
      exp1 <- read.csv(paste0("MagicBehavioural_", group[l], "_", version, "_", subjects[s], "_responses.csv"), header = T)

      exp1$Username <- gsub(" ", "", exp1$Username) #replace any spaces in the file names
      exp1$group <- group[l] #add group and ID
      exp1$groupEffectCoded <-ifelse(exp1$group == "exp", 1, -1) 
      exp1$ID <- paste0(group[l],s)
    }
    else { # if file does not exist, print into console and use filler file to create NAs
      print(paste0("MagicBehavioural_", group[l], "_", version, "_", subjects[s], "_responses.csv does not exist"))
      
      exp1 <- read.csv(file.path(dataDir, "filler_mainExp.csv"), header = T)
      exp1$group <- group[l]
      exp1$groupEffectCoded <-ifelse(exp1$group == "exp", 1, -1) 
      exp1$Username <- subjects[s]
      exp1$ID <-  paste0(group[l],s)
    }
    
    if (max(exp1$Trial) < 9){ # check whether participants have actually gone further than just the instructions; if not, print into console and use filler
      if (file.exists(dataset)){
        print(paste0("MagicBehavioural_", group[l], "_", version, "_", subjects[s], "_responses.csv is too short"))
      }
      exp1 <- read.csv(file.path(dataDir, "filler_mainExp.csv"), header = T)
      exp1$group <- group[l]
      exp1$groupEffectCoded <-ifelse(exp1$group == "exp", 1, -1) 
      exp1$Username <- subjects[s]
      exp1$ID <-  paste0(group[l],s)
    }
    
    # read in responses.csv files for memory part: check whether a file for a subject exists
    # uses a filler file in case there is no responses file for the memory test and prints into console
    setwd(memoryDir)
    f <- paste0("MagicBehavioural_memory_", version, "_", subjects[s], "_responses.csv")
    f_coded <- paste0(subjects[s], "_cuedRecall_CP.csv") # the suffix _SM indicated that I have coded the remaining answers in the recall task, suffix _CP means that Cristina has done it
    # note: this will have to be changed to _CP as Cristina is coded them now
    if (file.exists(f)){
      memory <- read.csv(f, header = T)
    }
    else { # if file does not exist, print into console and use filler file to create NAs
      print(paste0("MagicBehavioural_memory_", version, "_", subjects[s], "_responses.csv does not exist"))
      
      memory <- read.csv2(file.path(dataDir, "filler_memory.csv"), header = T)
      memory$group <- group[l]
      memory$groupEffectCoded <-ifelse(memory$group == "exp", 1, -1) 
      memory$Username <- subjects[s]
      memory$ID <-  paste0(group[l],s)
    }
    
    if (suppressWarnings(max(memory$Trial, na.rm = T)) < 5){ # check whether participants have gone further than instructions
      if (file.exists(f)){
        print(paste0("MagicBehavioural_memory_", version, "_", subjects[s], "_responses.csv does not have recall and recognition data"))
      }
      memory <- read.csv(file.path(dataDir, "filler_memory.csv"), header = T)
      memory$group <- group[l]
      memory$groupEffectCoded <-ifelse(memory$group == "exp", 1, -1) 
      memory$Username <- subjects[s]
      memory$ID <-  paste0(group[l],s)
    }
    
    memory$Username <- gsub(" ", "", memory$Username) #replace any spaces in the file names
    
    
    if("Resp_Start_Time_Date" %in% colnames(memory)){  # check whether the main data set has information about when it has started/finished 
      memory$startMemory <- memory$Resp_Start_Time_Date[1]
      memory$startMemoryNumeric <- memory$Resp_Start_Time_Date_Numeric[1]
      memory$endMemory <- memory$Resp_End_Time_Date[dim(exp1)[1]]
      memory$endMemoryNumeric <- memory$Resp_End_Time_Date_Numeric[dim(exp1)[1]]
    } else { # if not, add a note
      memory$startMemory <- "check manually"
      memory$startMemoryNumeric <- "check manually"
      memory$endMemory <- "check manually"
      memory$endMemoryNumeric <- "check manually" }
    
    if("Resp_Start_Time_Date" %in% colnames(exp1)){  # check whether the memory data set has information about when it has started/finished 
      exp1$startMain <- exp1$Resp_Start_Time_Date[1]
      exp1$startMainNumeric <- exp1$Resp_Start_Time_Date_Numeric[1]
      exp1$endMain <- exp1$Resp_End_Time_Date[dim(exp1)[1]]
      exp1$endMainNumeric <- exp1$Resp_End_Time_Date_Numeric[dim(exp1)[1]]
    } else { # if not, add a note
      exp1$startMain <- "check manually"
      exp1$startMainNumeric <- "check manually"
      exp1$endMain <- "check manually"
      exp1$endMainNumeric <- "check manually"}
    
    main <- exp1 # transfer exp1 to main
    main<-subset(main, Proc_comment == "experiment") # chose only the rows corresponding to the experiment itself
    
    # process the data from the first part of the experiment: change variable names
    names(main)[names(main)=="Stim_stimID"] <- "stimID"  
    names(main)[names(main)=="Resp_Curiosity"] <- "curiosity" 
    names(main)[names(main)=="Resp_curiosity_RT"] <- "curiosityRT" 
    names(main)[names(main)=="Resp_Answer"] <- "decision"
    names(main)[names(main)=="Resp_answer_RT"] <- "decisionRT"
    names(main)[names(main)=="Trial"] <- "trialMain"
    names(main)[names(main)=="Proc_item"] <- "itemMain"
    
    # process the data from the first part of the experiment: change numerics for item counter to delete instructions etc.
    main$itemMain <- main$itemMain-1
    main$trialMain<- main$trialMain-7
    
    # merge the response file of each participant with the info for each magic trick
    main <- merge(main, info, by = "stimID", all.y = T) # include all.y = T to produce NA if a participant has not completed the task
    main$Username <- subjects[s]
    main$group <- group[l]
    main$groupEffectCoded <-ifelse(main$group == "exp", 1, -1) 
    main$ID <-  paste0(group[l],s)
    
    # within subject calculations
    setwd(paste0(groupDir, group[l]))
    if (file.exists(dataset)){ # check whether data set exists
      # compute median split (WITHIN SUBJECT) for curiousity rating
      main$mediansplitCuriosityWithinSubject <- ifelse(main$curiosity > median(main$curiosity), "above", 
                                                       ifelse(main$curiosity < median(main$curiosity), "below", "median"))
      # compute group-mean centered curiosity rating
      main$curiosityGroupMeanCentered <- main$curiosity - mean(main$curiosity, na.rm = T)
      # compute reward by curiosity interaction
      main$rewardByCuriosity <- main$curiosityGroupMeanCentered * main$groupEffectCoded
    } else { # create NAs if no data
      main$mediansplitCuriosityWithinSubject <- NA 
      main$curiosityGroupMeanCentered  <- NA }
    
    
    # reduce main to relevant variables only
    main <- main[, c("Username", "ID", "group", "groupEffectCoded", "stimID", "length", "meanCuriosityStandardisedAya", "mediansplitCuriosityAya", 
                     "meanCuriositySample", "meanCuriosityStandardisedSample", "mediansplitCuriositySample", "mediansplitCuriosityWithinSubject",
                     "trialMain", "itemMain", "decision", "decisionRT", "curiosity", "curiosityGroupMeanCentered", "curiosityRT", "rewardByCuriosity")]
    
    # create variable for block
    main$blockMain <- ifelse(main$trialMain <= 12, 1, 
                             ifelse(main$trialMain > 24 , 3, 2))
    
    # save preprocessed data of video watching part in a csv file
    ifelse(!dir.exists(paste0(groupDir, group[l], "/preprocessed/")), dir.create(paste0(groupDir, group[l], "/preprocessed/")), FALSE)
    setwd(paste0(groupDir, group[l], "/preprocessed/"))
    write.csv(main, file = paste0(subjects[s], "_main.csv"), row.names = F)
    
    # process the recall memory data: select relevant rows and columns, change item counter
    cuedRecall <- subset(memory, memory$Proc_trial.type == "magic_recall")
    cuedRecall <- cuedRecall[,c("Username", "Trial", "Proc_item", "Stim_stimID", "Resp_trickAnswer")]
    names(cuedRecall) <- c("Username", "trialRecall", "itemRecall", "stimID", "description")
    cuedRecall$itemRecall <- cuedRecall$itemRecall-1
    cuedRecall$trialRecall <- cuedRecall$trialRecall-3
    
    # if participants have gone further than just reading the instructions for the memory test
    if (suppressWarnings(max(memory$Trial, na.rm = T)) > 5){
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
        } else { # if none of the spelling matches, insert NA --> manual coding necessary
          cuedRecall$cuedRecallStrict[j] = NA
          cuedRecall$cuedRecallLenient[j] = NA
        }
      }
    } else { # if none of the spelling matches, insert NA --> manual coding necessary
      cuedRecall$cuedRecallStrict = NA
      cuedRecall$cuedRecallLenient = NA
    }
    
    
    setwd(memoryDir)
    if (file.exists(f)) { # if the participants initially participated in the memory part, their preprocessed recall data is saved
      ifelse(!dir.exists(file.path(memoryDir, "preprocessed")), dir.create(file.path(memoryDir, "preprocessed")), FALSE)
      setwd(file.path(memoryDir, "preprocessed"))
      write.csv(cuedRecall, file = paste0(subjects[s], "_cuedRecall.csv"), row.names = F)    
      
      setwd(codedDir)
      f_coded <- paste0(subjects[s], "_cuedRecall_CP.csv") # the suffix _SM indicated that I have coded the remaining answers in the recall task
      
      if (file.exists(f_coded)) { #if the recall performance has already been coded, the recall df gets overwritten
        cuedRecall_coded <- read.csv(f_coded)
        cuedRecall_coded <- cuedRecall_coded[1:nrow(info),] # reduce cuedRecall_coded to only those rows that relate to stimuli
        cuedRecall <- cuedRecall_coded # overwrite cuedRecall
        
      } else { # if there is no _cuedRecall_CP file ALTHOUGH there was a memory data file, print into console
        if (suppressWarnings(max(memory$Trial, na.rm = T)) > 5){
          print(paste0(subjects[s], "_cuedRecall_CP.csv does not exist"))
        }
      }
    }
    
    # process recognition memory data: select relevant rows and convert variables to characters
    if (suppressWarnings(max(memory$Trial, na.rm = T)) < 41){  # check whether participants have done more than just the recall
      if (suppressWarnings(max(memory$Trial, na.rm = T)) < 5){
        recognition <- subset(memory, memory$Proc_trial.type == "magic_recognition")
      }
      if (suppressWarnings(max(memory$Trial, na.rm = T)) > 5){
        print(paste0("MagicBehavioural_memory_", version, "_", subjects[s], "_responses.csv does not have recognition data"))
        
        recognition <- read.csv(file.path(dataDir, "filler_recognition.csv"), header = T)
      }
    } else {
      recognition <- subset(memory, memory$Proc_trial.type == "magic_recognition")
    }
    
    recognition$group <- group[l]
    recognition$groupEffectCoded <-ifelse(recognition$group == "exp", 1, -1) 
    recognition$Username <- subjects[s]
    recognition$ID <-  paste0(group[l],s)
    
    recognition$option1 <- as.character(recognition$Stim_option1)
    recognition$answer <- as.character(recognition$Resp_Answer)
    
    # in the collector sheet, the correct answer is always in the option1 column
    # check whether the answer selected by participants is the correct answer
    recognition$recognition <- NA
    recognition$recognition <- ifelse(is.na(recognition$answer) == TRUE, NA, 
                                      ifelse(recognition$option1 == recognition$answer, 1, 0))
    
    # compute group-mean centered confidence rating
    recognition$confidenceGroupMeanCentered <- recognition$Resp_Confidence - mean(recognition$Resp_Confidence, na.rm = T)
    print(paste(subjects[s], mean(recognition$confidenceGroupMeanCentered)))
    
    # process the memory data: select relevant columns, change item counter
    recognition <- recognition[,c("Username", "Trial", "Proc_item", "Stim_stimID", "answer", "Resp_Confidence", "confidenceGroupMeanCentered", "recognition" )]
    names(recognition) <- c("Username", "trialRecognition", "itemRecognition", "stimID", "answer", "confidence", "confidenceGroupMeanCentered", "recognition")
    recognition$itemRecognition <- recognition$itemRecognition-1
    recognition$trialRecognition <- recognition$trialRecognition-40
    
    # compute scores of recognition performance
    recognition$confidenceCorrectTrials <- ifelse(recognition$recognition == 1, recognition$confidence, NA)
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
    
    
    setwd(memoryDir)
    if (file.exists(f)) { # if the participants initially participated in the memory part, their preprocessed recognition data is saved
      ifelse(!dir.exists(file.path(memoryDir, "preprocessed")), dir.create(file.path(memoryDir, "preprocessed")), FALSE)
      setwd(file.path(memoryDir, "preprocessed"))
      write.csv(recognition, file = paste0(subjects[s], "_recognition.csv"), row.names = F)
    }
    
    # merge all data (video watching, recall & recognition) together
    # all.x = T makes sure that NA are produced if a pps has not completed the full task
    recall <- merge(cuedRecall, recognition, by = c("Username", "stimID"), all.x = T)
    data <- merge(main, recall, by = c("Username", "stimID"), all.x = T)
    if (is.factor(data$cuedRecallLenient) == T){
      print(paste0(subjects[s], "_cuedRecall_CP.csv has cuedRecallLenient as factor"))
      data$cuedRecallLenient <- as.numeric(data$cuedRecallLenient)
      data$cuedRecallLenient <- ifelse(data$cuedRecallLenient == 1, 0,
                                       ifelse(data$cuedRecallLenient == 2, 1, NA))
    }
    if (is.factor(data$cuedRecallStrict) == T){
      print(paste0(subjects[s], "_cuedRecall_CP.csv has cuedRecallStrict as factor"))
      
      data$cuedRecallStrict <- as.numeric(data$cuedRecallStrict)
      data$cuedRecallStrict <- ifelse(data$cuedRecallStrict == 1, 0,
                                      ifelse(data$cuedRecallStrict == 2, 1, NA))
    }
    
    # save data per participant in long format, this includes data from movie watching and memory test
    ifelse(!dir.exists(file.path(preprocessedDir, "long")), dir.create(file.path(preprocessedDir, "long")), FALSE)
    setwd(file.path(preprocessedDir, "long"))
    write.csv(data, file = paste0(subjects[s], "_long.csv"), row.names = F)
    
    # process demographic information collected during movie watch part
    demogs <-subset(exp1, Proc_comment == "demogs details") 
    demogs <- demogs[,c("Username", "ID", "group", "groupEffectCoded", "Resp_survey_gender", "Resp_survey_age", "Resp_survey_education_years", "Resp_survey_education_level", "Resp_survey_profession", "Resp_survey_profession_student", 
                        "Resp_survey_ethnicity", "Resp_survey_sleep", "Resp_survey_sleep_average", "Resp_survey_health", "Resp_survey_english", "Resp_survey_handedness", "Resp_survey_neuro_disorders", "Resp_survey_contactdetails",
                        "startMain", "startMainNumeric", "endMain", "endMainNumeric")]
    names(demogs) <- c("Username", "ID", "group", "groupEffectCoded", "gender", "age", "education_years", "education_level", "profession", "profession_student", 
                       "ethnicity", "sleep_main", "sleep_average", "health", "survey_english", "handedness", "neuro_disorders", "contactdetails",
                       "startMain", "startMainNumeric", "endMain", "endMainNumeric")
    
    # process task motivation questions  collected during movie watch part
    postMain <-subset(exp1, Proc_comment == "post")
    
    if (dim(postMain)[1] == 0) { #if a participant has not completed the experiment, read in filler and print into console
      postMain <- read.csv(file = file.path(dataDir, "filler_postMain.csv"), header = T)
      postMain$group <- group[l]
      postMain$groupEffectCoded <- ifelse(postMain$group == "exp", 1, -1)
      postMain$Username <- subjects[s]
      postMain$ID <-  paste0(group[l],s)
      
      print(paste0("MagicBehavioural_", group[l], "_", version, "_", subjects[s], "_responses.csv does not have questionnaire data"))
      
    }
    
    # select relevant columns of task motivation data
    postMain <- postMain[,c("Username", "ID", "group", "groupEffectCoded", "Resp_survey_post1", "Resp_survey_post1_score", "Resp_survey_post2", "Resp_survey_post2_score",
                            "Resp_survey_post3", "Resp_survey_post3_score", "Resp_survey_post4", "Resp_survey_post4_score",
                            "Resp_survey_post5", "Resp_survey_post5_score", "Resp_survey_post6", "Resp_survey_post6_score",
                            "Resp_survey_post7", "Resp_survey_post7_score", "Resp_survey_post8", "Resp_survey_post8_score",
                            "Resp_survey_post9", "Resp_survey_post9_score", "Resp_survey_post10", "Resp_survey_post10_score",
                            "Resp_survey_post11", "Resp_survey_post11_score", "Resp_survey_post12", "Resp_survey_post12_score",
                            "Resp_survey_post13", "Resp_survey_post13_score", "Resp_survey_post14", "Resp_survey_post14_score",
                            "Resp_survey_post15", "Resp_survey_post15_score", "Resp_survey_post16", "Resp_survey_post16_score",
                            "Resp_survey_post17", "Resp_survey_post17_score", "Resp_survey_post18", "Resp_survey_post18_score",
                            "Resp_survey_post19", "Resp_survey_post19_score", "Resp_survey_post20", "Resp_survey_post20_score",
                            "Resp_survey_post21", "Resp_survey_post21_score", "Resp_survey_post22", "Resp_survey_post22_score",
                            "Resp_survey_post23", "Resp_survey_post23_score", "Resp_survey_post24", "Resp_survey_post24_score",
                            "Resp_survey_post25", "Resp_survey_post25_score",
                            "Resp_survey_post_comment1", "Resp_survey_post_comment2", "Resp_survey_post_comment3", "Resp_survey_post_comment4")]
    names(postMain) <- c("Username", "ID", "group", "groupEffectCoded", "post1", "post1_score", "post2", "post2_score", "post3", "post3_score", "post4", "post4_score",
                         "post5", "post5_score", "post6", "post6_score", "post7", "post7_score", "post8", "post8_score", "post9", "post9_score",
                         "post10", "post10_score", "post11", "post11_score", "post12", "post12_score", "post13", "post13_score", "post14", "post14_score",
                         "post15", "post15_score", "post16", "post16_score", "post17", "post17_score", "post18", "post18_score", "post19", "post19_score",
                         "post20", "post20_score", "post21", "post21_score", "post22", "post22_score", "post23", "post23_score", "post24", "post24_score",
                         "post25", "post25_score", "post_comment1", "post_comment2", "post_comment3", "post_comment4")
    
    
    # compute item scores
    items <- c("post1", "post3", "post4", "post5", "post6", "post7", "post8", "post9", "post10", "post11", "post12", "post13", "post15", "post16", "post19", "post21", "post22", "post23", "post25")
    itemsReversed <- c("post2", "post14", "post17", "post18", "post20")
    
    
    for (item in items){
      item_score <-  paste0(item, "_score")
      postMain[,item_score] <- ifelse(postMain[,item] == "Strongly agree", 7, 
                                      ifelse(postMain[,item] == "Somehow agree", 6, 
                                             ifelse(postMain[,item] == "Slightly agree", 5, 
                                                    ifelse(postMain[,item] == "Neither agree nor disagree", 4, 
                                                           ifelse(postMain[,item] == "Slightly disagree", 3, 
                                                                  ifelse(postMain[,item] == "Somehow disagree", 2, 
                                                                         ifelse(postMain[,item] == "Strongly disagree", 1, NA)))))))
      rm(item_score)
    }
    
    for (item in itemsReversed){
      item_score <-  paste0(item, "_score")
      postMain[,item_score] <- ifelse(postMain[,item] == "Strongly agree", 1, 
                                      ifelse(postMain[,item] == "Somehow agree", 2, 
                                             ifelse(postMain[,item] == "Slightly agree", 3, 
                                                    ifelse(postMain[,item] == "Neither agree nor disagree", 4, 
                                                           ifelse(postMain[,item] == "Slightly disagree", 5, 
                                                                  ifelse(postMain[,item] == "Somehow disagree", 6, 
                                                                         ifelse(postMain[,item] == "Strongly disagree", 7, NA)))))))
      rm(item_score)
    }
    
    # compute scales
    postMain$intrinsicMotivation <- (postMain$post1_score + postMain$post2_score + postMain$post3_score)/3
    postMain$taskEngagement <- (postMain$post4_score + postMain$post5_score + postMain$post6_score)/3
    postMain$interest <- (postMain$post7_score + postMain$post8_score + postMain$post9_score)/3
    postMain$boredom <- (postMain$post10_score + postMain$post11_score + postMain$post12_score)/3
    postMain$effort <- (postMain$post13_score + postMain$post14_score + postMain$post15_score + postMain$post16_score + postMain$post17_score)/5
    postMain$pressure <- (postMain$post18_score + postMain$post19_score + postMain$post20_score + postMain$post21_score)/4
    
    
    # process task m questions collected during memory part
    postMemory <-subset(memory, Proc_comment == "questionnaire")
    
    if (dim(postMemory)[1] == 0) { #if a participant has not completed the experiment, read in filler and print into console
      postMemory <- read.csv(file = file.path(dataDir, "filler_postMemory.csv"), header = T)
      postMemory$group <- group[l]
      postMemory$groupEffectCoded <-  ifelse(postMain$group == "exp", 1, -1)
      postMemory$Username <- subjects[s]
      postMemory$ID <-  paste0(group[l],s)
      
      print(paste0("MagicBehavioural_memory_", version, "_", subjects[s], "_responses.csv does not have questionnaire data"))
    }
    
    # select relevant columns of memory task questions data
    postMemory <- postMemory[,c("Username", "ID","startMemory", "startMemoryNumeric", "endMemory", "endMemoryNumeric",
                                "Resp_survey_sleep","Resp_survey_sleep_hours", 
                                "Resp_survey_test_known", "Resp_survey_test_known_score", "Resp_survey_memory_intention", "Resp_survey_memory_intention_score", 
                                "Resp_survey_connection", "Resp_survey_connection_score", "Resp_survey_comment")]
    names(postMemory) <- c("Username", "ID", "startMemory", "startMemoryNumeric", "endMemory", "endMemoryNumeric", 
                           "sleep_memory","sleep_hours", "test_known", "test_known_score", "memory_intention", "memory_intention_score", "connection", "connection_score", "comment")
    
    # merge all the data from demogs, task motivation and post memory test questions
    subjectData <- merge(demogs, postMain, by = c("Username", "ID", "group", "groupEffectCoded"))
    subjectData <- merge(subjectData, postMemory, by = c("Username"))
    
    # ADD column for complete data: check data
    subjectData$completeData <- ifelse(is.na(subjectData$connection) == F, "completeDataset",
                                       ifelse(any(is.na(data$confidence)) == F, "completedRecognitionData",
                                              ifelse(any(is.na(data$description)) == F, "completedRecallData",
                                                     ifelse(is.na(subjectData$post_comment3) == F, "completedMovieWatchingDataInclQuestionnaire",
                                                            ifelse(any(is.na(data$curiosity)) == F, "completedMovieWatchingData",
                                                                   "noData")))))
    
    # check whether there is coded data from the memory test at all; if so compute sum scores for recall performance
    # if there is no data, NA will be added when rbinding all information across subjects
    setwd(codedDir)
    if (file.exists(f_coded)) {
      
      for (block in 1:(max(data$blockMain)+1)) {
        data_subset <- subset(data, data$blockMain == block)
        
        if(block == 4){ # as well as for the data set in total
          data_subset <- data
        }
        
        subjectData[[paste0("cuedRecallLenient", blockstring[block])]] <- sum(data_subset$cuedRecallLenient, na.rm = T) #please note that this needs to be changed as it is not looking at any form of coded data
        subjectData[[paste0("cuedRecallStrict", blockstring[block])]] <- sum(data_subset$cuedRecallStrict, na.rm = T) #please note that this needs to be changed as it is not looking at any form of coded data
        
      }
    }
    
    # check whether there is  data from the memory test at all; if so compute sum scores for recognition performance
    # if there is no data, NA will be added when rbinding all information across subjects
    setwd(memoryDir)
    if (file.exists(f)) {
      
      # for each of the blocks
      for (block in 1:(max(data$blockMain)+1)) {
        data_subset <- subset(data, data$blockMain == block)
        
        if(block == 4){ # as well as for the data set in total
          data_subset <- data
        }
        
        # sum up the scores for the recognition task, once in total and once seperated for the different levels of confidence
        subjectData[[paste0("recognition", blockstring[block])]] <- sum(data_subset$recognition, na.rm = T)
        subjectData[[paste0("recognitionAboveMeanConf", blockstring[block])]] <- sum(data_subset$recognitionAboveMeanConf, na.rm = T)
        subjectData[[paste0("meanConfidence", blockstring[block])]]  <- mean(data_subset$confidence, na.rm = T)
        temp_data <- subset(data_subset, data_subset$recognition == 1)
        subjectData[[paste0("meanConfidenceCorrectTrials", blockstring[block])]]   <- mean(temp_data$confidence, na.rm = T)
        
        for (k in 1:6) { #confidence ranges from 1 to 6, potentially code can be made more flexible by using min(data$confidence) and max(data$confidence)
          temp_data <- subset(data_subset, data_subset$confidence == k)
          subjectData[[paste0("recognitionConfLevel_", k, blockstring[block])]] <- sum(temp_data$recognition, na.rm = T)
          
          if (k < 6) {
            temp_data_above <- subset(data_subset, data_subset$confidence > k)
            subjectData[[paste0("recognitionConfLevel_above_", k, blockstring[block])]] <-sum(temp_data_above$recognition, na.rm = T) 
            rm(temp_data_above)
          }
          
          # sum up the scores for the recognition task pooled
          if (k == 1 || k == 3 || k == 5) {
            temp_data_plus1 <- subset(data_subset, data_subset$confidence == (k+1))
            subjectData[[paste0("recognitionConfLevel_", k, "_", k+1, blockstring[block])]] <- sum(temp_data$recognition, na.rm = T) + sum(temp_data_plus1$recognition, na.rm = T)
            rm(temp_data_plus1)
          }
          if (k == 1 || k == 4) {
            temp_data_plus1 <- subset(data_subset, data_subset$confidence == (k+1))
            temp_data_plus2 <- subset(data_subset, data_subset$confidence == (k+2))
            subjectData[[paste0("recognitionConfLevel_", k, "_", k+1, "_", k+2, blockstring[block])]] <- sum(temp_data$recognition, na.rm = T) + sum(temp_data_plus1$recognition, na.rm = T) + sum(temp_data_plus2$recognition, na.rm = T)
            rm(temp_data_plus1)
            rm(temp_data_plus2)
          }
          if (suppressWarnings(max(memory$Trial, na.rm = T)) < 41){
            subjectData$recognition <- NA
            subjectData[[paste0("recognitionConfLevel_", k)]] <- NA
            if (k == 1 || k == 3 || k == 5) {subjectData[[paste0("recognitionConfLevel_", k, "_", k+1, blockstring[block])]] <- NA}
            if (k == 1 || k == 4) { subjectData[[paste0("recognitionConfLevel_", k, "_", k+1, "_", k+2, blockstring[block])]] <- NA}
          }
          rm(temp_data)
        }
      }
    }
    
    # save the file containing all information about demogs, task motivation / post memory test questions and memory performance
    ifelse(!dir.exists(file.path(preprocessedDir, "wide")), dir.create(file.path(preprocessedDir, "wide")), FALSE)
    setwd(file.path(preprocessedDir, "wide"))
    write.csv(subjectData, file = paste0(subjects[s], "_wide.csv"), row.names = F)
    
    # ´removing all variables created, comment out for debugging!
    rm(cuedRecall)
    if (file.exists(file.path(codedDir, f_coded))) {
      rm(cuedRecall_coded)
    }
    rm(data,demogs,exp1,main,memory,postMain,postMemory,recall,recognition,subjectData)
  }
  
}

#### combine single _wide.csv files to one file ####
setwd(file.path(preprocessedDir, "wide"))
file_list <- list.files(pattern = "wide.csv")

for (m in seq_along(file_list)){
  if (m == 1) {
    dataWide <- read.csv(file_list[m], header=T)
  }
  else {
    temp_datawide <- read.csv(file_list[m], header=T)
    dataWide <- rbind.all.columns(dataWide, temp_datawide) #rbind all columns will induce NA if there was initially no data saved in the loop per participant
    rm(temp_datawide)
  }
}
setwd(preprocessedDir)
write.xlsx(dataWide, file=paste0("wide_MagicBehavioural_", version, ".xlsx"), sheetName = "Sheet1", row.names = F) 


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
dataLong$confidenceCorrectTrials <- ifelse(dataLong$recognition == 1, dataLong$confidence, NA)

write.xlsx(dataLong, file=paste0("long_MagicBehavioural_", version, ".xlsx"), sheetName = "Sheet1", row.names = F) 

#### compute summary statistics for measurements of memory ####
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
dataWideRecognitionPerformance <- dataWide[,recognitionPerformanceVars]
mean <- data.frame(apply(dataWideRecognitionPerformance, 2, mean,  na.rm = T), recognitionPerformanceVars)
sd <- data.frame(apply(dataWideRecognitionPerformance, 2, sd,  na.rm = T), recognitionPerformanceVars)
min <- data.frame(apply(dataWideRecognitionPerformance, 2, min,  na.rm = T), recognitionPerformanceVars)
max <- data.frame(apply(dataWideRecognitionPerformance, 2, max,  na.rm = T), recognitionPerformanceVars)

descriptivesRecognitionPerformance <- merge(mean, sd, by = "recognitionPerformanceVars")
descriptivesRecognitionPerformance <- merge(descriptivesRecognitionPerformance, min, by = "recognitionPerformanceVars")
descriptivesRecognitionPerformance <- merge(descriptivesRecognitionPerformance, max, by = "recognitionPerformanceVars")

names(descriptivesRecognitionPerformance) <- c("randomisations", "mean", "sd", "min", "max")

setwd(preprocessedDir)
write.xlsx(descriptivesRecognitionPerformance, file= paste0("memoryPerformance_MagicBehavioural_", version, "_", format(Sys.time(), "%Y-%m-%d"), ".xlsx"), sheetName = "Sheet1", row.names = F) 


#### create a file for each magic trick seperately ####
# purpose: to check the coded answers for consistency
tricks <- info$stimID
for (t in seq_along(tricks)){
  trickdata <- subset(dataLong, dataLong$stimID == tricks[t])
  write.csv(trickdata, file = paste0("~/Dropbox/Reading/PhD/Magic tricks/behavioural_study/data_fin/Analysis/Tricks/MagicBehavioural_memory_", version, "_", tricks[t], ".csv"), row.names = F)
}

#### create a seperate wide format file with magic tricks as columns and participants as rows ####
# purpose of this is to have an document that shows whether the items are easy/difficult enough
# check whether this can be done better, but so far, it does the trick
library(reshape)

indicesPerTrick <- c("stimID","decision", "decisionRT", "curiosity", "curiosityRT", "cuedRecallLenient", "cuedRecallStrict", 
                     "recognition", "recognitionConfLevel_6", "recognitionConfLevel_5_6", "recognitionConfLevel_4_5_6",
                     "confidence", "confidenceCorrectTrials")

# rm(decision, decisionRT, curiosity, curiosityRT, cuedRecallLenient, cuedRecallStrict, 
#                         recognitionAll, recognitionConfLevel_6, recognitionConfLevel_5_6, recognitionConfLevel_4_5_6,
#                         confidence, confidenceCorrectTrials)
dataLongTrick <- dataLong[,indicesPerTrick]

indicesPerTrick <- indicesPerTrick[-1]


for (i in seq_along(indicesPerTrick)){
  # 
  # if (indicesPerTrick[i] == "recognition") {
  #   assign(paste0(indicesPerTrick[i], "All"), cast(dataLong, ID~stimID,value=paste0(indicesPerTrick[i])))
  #   indicesPerTrick[i] = "recognitionAll"
  # } else {
    assign(paste(indicesPerTrick[i]), cast(dataLong, ID~stimID,value=paste0(indicesPerTrick[i])))
  #}
  
}

tricks <- as.character(info$stimID)

# dfMeans <- data.frame(stimID=tricks, meanDecision = "", meanDecisionRT = "", meanCuriosity = "", meanCuriosityRT = "", 
#                       meanCuedRecallLenient = "",  meanCuedRecallStrict = "", meanRecognition = "", 
#                       meanRecognitionConfLevel_6 = "", meanRecognitionConfLevel_5_6 = "", meanRecognitionConfLevel_4_5_6 = "",
#                       meanConfidence = "", meanConfidenceCorrectTrials = "")
# indicesPerTrickMean <- names(dfMeans)
# indicesPerTrickMean <- indicesPerTrickMean[-1]
indicesPerTrickMean <- c("meanDecision", "meanDecisionRT", "meanCuriosity", "meanCuriosityRT", "meanCuedRecallLenient", "meanCuedRecallStrict", 
                         "meanRecognition", "meanRecognitionConfLevel_6", "meanRecognitionConfLevel_5_6",
                         "meanRecognitionConfLevel_4_5_6", "meanConfidence", "meanConfidenceCorrectTrials" )

# calculate mean values for each magic trick for each index
for (iMean in seq_along(indicesPerTrickMean)){
  assign(paste(indicesPerTrickMean[iMean]), numeric(length(tricks)))
  currentMean <- numeric(length(tricks))
  
  currentIndex <- indicesPerTrick[iMean] # get current index
  meansPerTrick <-   colMeans(get(currentIndex), na.rm = T) # calculate mean for each magic trick
  
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

dfMeans$problematic <- ifelse(dfMeans$meanRecognition == 1, "too easy", 
                              ifelse(dfMeans$meanRecognition == 0, "too difficult", "all good")) 

dfMeans$stimID[dfMeans$problematic == "too easy"]
dfMeans$stimID[dfMeans$problematic == "too difficult"]

dfMeans$stimID[dfMeans$meanCuedRecallLenient < 0.1]
dfMeans$stimID[dfMeans$meanCuedRecallStrict < 0.1]

dfMeans$stimID[dfMeans$meanRecognition < 0.25]
dfMeans$stimID[dfMeans$meanConfidenceCorrectTrials < 4]

recognitionTrials <- read.csv("~/Dropbox/Apps/Open-Collector/experiments/MagicBehavioural_memory_fin/Stimuli/Stimuli.csv")

dfStimuli <- merge(dfMeans, recognitionTrials)
dfStimuli <- merge(dfStimuli, info, by = "stimID")

dfStimuli$meanCuriosityStandardisedSample <- (dfStimuli$meanCuriosity - mean(dfStimuli$meanCuriosity)) / sd(dfStimuli$meanCuriosity)

dfStimuli$mediansplitCuriositySample <- ifelse(dfStimuli$meanCuriositySample > median(dfStimuli$meanCuriositySample), "above", 
                                               ifelse(dfStimuli$meanCuriositySample < median(dfStimuli$meanCuriositySample), "below", "median")) 


dfStimuli$differentSplits <- ifelse(dfStimuli$mediansplitCuriosityAya != dfStimuli$mediansplitCuriositySample, "different", "same") 
dfStimuli$stimID[dfStimuli$differentSplits == "different"]
dfStimuli$order.curious[dfStimuli$differentSplits == "different"]

dfStimuli <- dfStimuli[,c("stimID","meanDecision","meanDecisionRT","meanCuriosity","meanCuriosityRT",
                          "meanCuriosityStandardisedAya","mediansplitCuriosityAya",
                          "meanCuriosityStandardisedSample","mediansplitCuriositySample",
                          "meanCuedRecallLenient","meanCuedRecallStrict","meanRecognition",
                          "meanRecognitionConfLevel_6","meanRecognitionConfLevel_5_6","meanRecognitionConfLevel_4_5_6",
                          "option1","option2","option3","option4")]

setwd(file.path(mainDir, "stimuli"))
write.xlsx(dfMeans, file= paste0("stimuli_MagicBehavioural_memoryPerformance_",version, "_", format(Sys.time(), "%Y-%m-%d"), ".xlsx"), sheetName = "Sheet1", row.names = F) 

setwd(file.path(mainDir, "stimuli"))
write.xlsx(dfStimuli, file= paste0("stimuli_MagicBehavioural_",version, "_", format(Sys.time(), "%Y-%m-%d"), ".xlsx"), sheetName = "Sheet1", row.names = F) 

