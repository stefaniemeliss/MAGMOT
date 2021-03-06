##################################################################
### extract information from the txt files created by pyfMRIqc ###
##################################################################

# determine directories
path="/storage/shared/research/cinn/2018/MAGMOT"
BIDS_dir=file.path(path, "MAGMOT_BIDS")
qc_dir=file.path(path, "derivatives", "pyfMRIqc")

# create variable subjects based on folders in fMRIqc directory
setwd(qc_dir)
subjects <- list.dirs(recursive = F)
subjects <- gsub('./', '', subjects)
subjects <- subjects[grepl("sub", subjects)]

# loop over subjects
for (s in seq_along(subjects)){
  
  # go into subjects folder
  setwd(file.path(qc_dir, subjects[s]))
  
  # create variable with file names
  files <- list.files(pattern = ".txt", recursive = T)
  
  # loop over txt files
  for (f in seq_along(files)) {
    
    # read in file
    param <- read.delim(files[f], sep = ":")
    
    # format param
    param <- as.data.frame(t(param)) # transpose
    rownames(param) <- NULL # remove row names
    for (col in 1:dim(param)[2]){
      param[,col] <- as.character(param[,col])  # change all columns from factor to character
    }
    param[1,] <- gsub(' ', '_', param[1,]) # replace spaces
    param[1,] <- gsub('__', '', param[1,]) # remove double underscore
    colnames(param) <- param[1,] # use first row to determine column names
    param <- param[-1,] # delete first row
    
    # determine group
    if (grepl("control", subjects[s])){
      group <- "control"
    } else {
      group <- "experimental"
    }
    
    # determine file name etc
    if (grepl("magictrickwatching_run-1", files[f])){
      task <- "magictrickwatching"
      run <- 1
      BOLD <- "magictrickwatching_run-1"
    } else if (grepl("magictrickwatching_run-2", files[f])){
      task <- "magictrickwatching"
      run <- 2
      BOLD <- "magictrickwatching_run-2"
    } else if (grepl("magictrickwatching_acq-1_run-2", files[f])){
      task <- "magictrickwatching"
      run <- 2
      BOLD <- "magictrickwatching_run-2_acq-1"
    } else if (grepl("magictrickwatching_acq-2_run-2", files[f])){
      task <- "magictrickwatching"
      run <- 2
      BOLD <- "magictrickwatching_run-2_acq-2"
    } else if (grepl("magictrickwatching_run-3", files[f])){
      task <- "magictrickwatching"
      run <- 3
      BOLD <- "magictrickwatching_run-3"
    } else if (grepl("rest_run-1", files[f])){
      task <- "rest"
      run <- 1
      BOLD <- "rest_run-1"
    } else if (grepl("rest_run-2", files[f])){
      task <- "rest"
      run <- 2
      BOLD <- "rest_run-2"
    }  
    
    # determine scan
    scan <- paste0(subjects[s], "_task-", BOLD)
    
    # determine values to extract
    extract <- c("SNR_voxel_MEAN", "SNR_voxel_STD", "SNR_voxel_value_range", "Mean", "Mean_(mask)", "SD", "SD_(mask)", 
                 "Min_Slice_SNR", "Max_Slice_SNR", "Mean_voxel_SNR", "Mean_absolute_Movement", "Max_absolute_Movement", 
                 "Max_relative_Movement", "Relative_movements_(>0.1mm)", "Relative_movements_(>0.5mm)", "Relative_movements_(>voxelsize)") 
    
    
    # create data frame
    df <- data.frame(subject = subjects[s],
                     group = group,
                     task = task,
                     run = run,
                     BOLD = BOLD,
                     scan = scan,
                     param = character(length(extract)),
                     value = numeric(length(extract)),
                     stringsAsFactors=FALSE)
    
    
    
    # extract parameters
    for(v in 1:length(extract)) {
      # paste variable
      df$param[v] <- paste(extract[v])
      df$value[v] <- param[1,paste(extract[v])]
    }
    
    # combine information from different subjects
    if (s == 1 && f == 1){ # for the first run of the first subject
      scanparam <-  df
    } else {
      temp_scanparam <-  df
      scanparam <- rbind(scanparam, temp_scanparam) 
      rm(temp_scanparam)
    }
  }
  
  # when finished with processing of last subject
  if(s == 50){
    
    # make value numeric
    scanparam$range <-ifelse(scanparam$param == "SNR_voxel_value_range", scanparam$value, NA)
    scanparam$value <-ifelse(scanparam$param != "SNR_voxel_value_range", as.numeric(scanparam$value), NA)
    
    # add censor information
    setwd(qc_dir)
    censor <- read.delim("afni_volreg_censor_count.tsv")
    censor$scan <- gsub("_bold_concat.nii", "", censor$scan)
    censor$scan <- gsub("_bold.nii", "", censor$scan)
    
    # transform from wide into long format
    long <- reshape2::melt(censor, id.vars = c("subject", "scan"))
    names(long) <- c("subject", "scan", "param", "value")
    
    # add columns
    long$group <- ifelse(grepl("control", long$subject), "control", "experimental")
    
    long$task <- ifelse(grepl("magictrickwatching", long$scan), "magictrickwatching", "rest")
    long$run <- ifelse(grepl("run-1", long$scan), 1, 
                       ifelse(grepl("run-2", long$scan), 2, 3))
    long$BOLD <- ifelse(grepl("magictrickwatching_run-1", long$scan), "magictrickwatching_run-1", 
                        ifelse(grepl("magictrickwatching_run-2", long$scan), "magictrickwatching_run-2", 
                               ifelse(grepl("magictrickwatching_run-3", long$scan), "magictrickwatching_run-3", 
                                      ifelse(grepl("rest_run-1", long$scan), "rest_run-1", 
                                             ifelse(grepl("rest_run-2", long$scan), "rest_run-2",
                                                    ifelse(grepl("magictrickwatching_acq-1_run-2", long$scan), "magictrickwatching_acq-1_run-2",
                                                           ifelse(grepl("magictrickwatching_acq-2_run-2", long$scan), "magictrickwatching_acq-2_run-2",
                                                                  NA)))))))
    long$range <- NA
    
    # merge
    scanparam <- rbind(scanparam, long)
    
    # order rows
    scanparam <- scanparam[order(scanparam$subject, scanparam$BOLD),]
    
    # save data in csv
    write.csv(scanparam, file = "pyfMRIqc_output.csv", row.names = F)
    
  }
}
