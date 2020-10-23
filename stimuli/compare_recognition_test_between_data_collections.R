rm(list = ls())

setwd("~/Dropbox/Reading/PhD/Magictricks/stimuli/")

# download data from OSF
osfr::osf_auth() # log into OSF
project <- osfr::osf_retrieve_node("fhqb7")
target_dir <- osfr::osf_ls_files(project, pattern = "stimuli") # looks at all files and directories in the project and defines the match with "data"
osfr::osf_ls_files(target_dir, pattern = ".xlsx") %>%
  osfr::osf_download(conflicts = "overwrite")

fmri <- xlsx::read.xlsx("recognition_memory_tests.xlsx", sheetName = "fmri", stringsAsFactors=FALSE)
main <- xlsx::read.xlsx("recognition_memory_tests.xlsx", sheetName = "main", stringsAsFactors=FALSE)
pilot <- xlsx::read.xlsx("recognition_memory_tests.xlsx", sheetName = "pilot", stringsAsFactors=FALSE)

# read in information from Ozono et al 2020
ozono <- read.csv("Ozono_et_al_2020_Detailed_information_about_MagicCATs.csv", stringsAsFactors = F)
ozono$X <- NULL #remove col
names(ozono) <- c("stimID", "Name", "Credit", "Phenomena Category", "Materials", "Length", "Subtitle", "Description")
ozono$stimID <- gsub("Short", "short", ozono$stimID)
ozono$stimID <- gsub("Long", "long", ozono$stimID)

# compare stimids --> if 36, then they are all the same!
sum(fmri$stimid == main$stimid)
sum(pilot$stimid == main$stimid)
sum(fmri$stimid == pilot$stimid)

# compare dropbox links
sum(fmri$image == main$image)
sum(pilot$image == main$image)
sum(fmri$image == pilot$image)

# compare option1
sum(fmri$option1 == main$option1) # 36
sum(pilot$option1 == main$option1) # 34
sum(fmri$option1 == pilot$option1) # 34
    
# compare option2
sum(fmri$option2 == main$option2) # 36
sum(pilot$option2 == main$option2) # 34
sum(fmri$option2 == pilot$option2) # 34

# compare option3
sum(fmri$option3 == main$option3) # 36
sum(pilot$option3 == main$option3) # 36
sum(fmri$option3 == pilot$option3) # 36

# compare option4
sum(fmri$option4 == main$option4) # 36
sum(pilot$option4 == main$option4) # 35
sum(fmri$option4 == pilot$option4) # 35

# rename columns
names(fmri) <- c("stimid", "image", "option1_fmri", "option2_fmri", "option3_fmri", "option4_fmri")
names(main) <- c("stimid", "image", "option1_main", "option2_main", "option3_main", "option4_main")
names(pilot) <- c("stimid", "image", "option1_pilot", "option2_pilot", "option3_pilot", "option4_pilot")

# merge all together
all <- merge(fmri, main, by = c("stimid", "image"))
all <- merge(all, pilot, by = c("stimid", "image"))

# compare the wording for each option between the data collections
all$option1_same_fmri_main <- all$option1_fmri == all$option1_main
all$option1_same_fmri_pilot <- all$option1_fmri == all$option1_pilot
all$option1_same_main_pilot <- all$option1_pilot == all$option1_main

all$option2_same_fmri_main <- all$option2_fmri == all$option2_main
all$option2_same_fmri_pilot <- all$option2_fmri == all$option2_pilot
all$option2_same_main_pilot <- all$option2_pilot == all$option2_main

all$option3_same_fmri_main <- all$option3_fmri == all$option3_main
all$option3_same_fmri_pilot <- all$option3_fmri == all$option3_pilot
all$option3_same_main_pilot <- all$option3_pilot == all$option3_main

all$option4_same_fmri_main <- all$option4_fmri == all$option4_main
all$option4_same_fmri_pilot <- all$option4_fmri == all$option4_pilot
all$option4_same_main_pilot <- all$option4_pilot == all$option4_main

# check data
colSums(all[grepl("same", names(all))])
# this shows that there are no differences in the memory test between fmri and main, but that there are some deviations between pilot and the other two 


# create dataframe memory test
memory <- xlsx::read.xlsx("recognition_memory_tests.xlsx", sheetName = "fmri", stringsAsFactors=FALSE)


# copy the wording from pilot if it differs from the wording in fmri or main
memory$comments_option1 <- ifelse(all$option1_same_fmri_pilot == FALSE | all$option1_same_main_pilot == FALSE, all$option1_pilot, "")
memory$comments_option2 <- ifelse(all$option2_same_fmri_pilot == FALSE | all$option2_same_main_pilot == FALSE, all$option2_pilot, "")
memory$comments_option3 <- ifelse(all$option3_same_fmri_pilot == FALSE | all$option3_same_main_pilot == FALSE, all$option3_pilot, "")
memory$comments_option4 <- ifelse(all$option4_same_fmri_pilot == FALSE | all$option4_same_main_pilot == FALSE, all$option4_pilot, "")

# merge comments
memory$deviation_in_pilot <- paste0(memory$comments_option1, memory$comments_option2, memory$comments_option3, memory$comments_option4)
memory$deviation_in_pilot <- gsub(".M", ". M", memory$deviation_in_pilot)
names(memory) <- c("stimID", "Cue image", "Recognition option 1", "Recognition option 2", "Recognition option 3", "Recognition option 4", "comments_option1", "comments_option2", "comments_option3", "comments_option4", "Different wording in pilot")
  
# merge this file with the information from Ozono et al 2020
memory <- merge(memory, ozono, by.x = "stimID", all.x = T)

# select necessary columns only
memory <- memory[,c("stimID", "Name", "Credit", "Phenomena Category", "Materials", "Length", "Description","Recognition option 1", "Recognition option 2", "Recognition option 3", "Recognition option 4", "Different wording in pilot")]

# save and upload file
write.csv(memory, file = "recognition_memory_test_incl_deviations_in_pilot.csv", row.names = F)
file_exists <- osfr::osf_ls_files(target_dir, pattern = "recognition_memory_test_incl_deviations_in_pilot.csv")
if (dim(file_exists)[1] > 0){ # delete file if it exists
  osfr::osf_rm(file_exists, recurse = T, verbose = FALSE, check = F)
}
osfr::osf_upload(target_dir, path = "recognition_memory_test_incl_deviations_in_pilot.csv", conflicts = "overwrite")





option1 <- all[, c("stimid", "image", "option1_fmri", "option1_main", "option1_pilot")]
option2 <- all[, c("stimid", "image", "option2_fmri", "option2_main", "option2_pilot")]
option3 <- all[, c("stimid", "image", "option3_fmri", "option3_main", "option3_pilot")]
option4 <- all[, c("stimid", "image", "option4_fmri", "option4_main", "option4_pilot")]

# create df in long format for each option
option1_long <- reshape2::melt(option1, id.vars=c("stimid", "image"), value.name = "option1")
names(option1_long)[names(option1_long) == "variable"] <- "exp"
option1_long$exp <- gsub("option1_", "", option1_long$exp)
option2_long <- reshape2::melt(option2, id.vars=c("stimid", "image"), value.name = "option2")
names(option2_long)[names(option2_long) == "variable"] <- "exp"
option2_long$exp <- gsub("option2_", "", option2_long$exp)
option3_long <- reshape2::melt(option3, id.vars=c("stimid", "image"), value.name = "option3")
names(option3_long)[names(option3_long) == "variable"] <- "exp"
option3_long$exp <- gsub("option3_", "", option3_long$exp)
option4_long <- reshape2::melt(option4, id.vars=c("stimid", "image"), value.name = "option4")
names(option4_long)[names(option4_long) == "variable"] <- "exp"
option4_long$exp <- gsub("option4_", "", option4_long$exp)

# combine all four options together again
comp <- merge(option1_long, option2_long, by = c("stimid", "image", "exp"))
comp <- merge(comp, option3_long, by = c("stimid", "image", "exp"))
comp <- merge(comp, option4_long, by = c("stimid", "image", "exp"))

# save file
write.csv(comp, file = "comparison_recognotion_test_between_data_collections.csv", row.names = F)
osfr::osf_upload(target_dir, path = "comparison_recognotion_test_between_data_collections.csv", conflicts = "overwrite")


# compare the wording for each option between the data collections
options <- data.frame(numeric(length = dim(fmri)[1]))
options$option1_same_fmri_main<- comp$option1[comp$exp == "fmri"] == comp$option1[comp$exp == "main"]
options$option1_same_fmri_pilot<- comp$option1[comp$exp == "fmri"] == comp$option1[comp$exp == "pilot"]
options$option1_same_main_pilot<- comp$option1[comp$exp == "main"] == comp$option1[comp$exp == "pilot"]

options$option2_same_fmri_main<- comp$option2[comp$exp == "fmri"] == comp$option2[comp$exp == "main"]
options$option2_same_fmri_pilot<- comp$option2[comp$exp == "fmri"] == comp$option2[comp$exp == "pilot"]
options$option2_same_main_pilot<- comp$option2[comp$exp == "main"] == comp$option2[comp$exp == "pilot"]

options$option3_same_fmri_main<- comp$option3[comp$exp == "fmri"] == comp$option3[comp$exp == "main"]
options$option3_same_fmri_pilot<- comp$option3[comp$exp == "fmri"] == comp$option3[comp$exp == "pilot"]
options$option3_same_main_pilot<- comp$option3[comp$exp == "main"] == comp$option3[comp$exp == "pilot"]

options$option4_same_fmri_main<- comp$option4[comp$exp == "fmri"] == comp$option4[comp$exp == "main"]
options$option4_same_fmri_pilot<- comp$option4[comp$exp == "fmri"] == comp$option4[comp$exp == "pilot"]
options$option4_same_main_pilot<- comp$option4[comp$exp == "main"] == comp$option4[comp$exp == "pilot"]
