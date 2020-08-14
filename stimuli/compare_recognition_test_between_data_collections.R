# download data from OSF
osfr::osf_auth() # log into OSF
project <- osfr::osf_retrieve_node("fhqb7")
target_dir <- osfr::osf_ls_files(project, pattern = "stimuli") # looks at all files and directories in the project and defines the match with "data"
osfr::osf_ls_files(target_dir, pattern = ".xlsx") %>%
  osfr::osf_download(conflicts = "overwrite")

fmri <- xlsx::read.xlsx("recognition_memory_tests.xlsx", sheetName = "fmri", stringsAsFactors=FALSE)
main <- xlsx::read.xlsx("recognition_memory_tests.xlsx", sheetName = "main", stringsAsFactors=FALSE)
pilot <- xlsx::read.xlsx("recognition_memory_tests.xlsx", sheetName = "pilot", stringsAsFactors=FALSE)

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
