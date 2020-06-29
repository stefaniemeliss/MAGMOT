##################################################################
### analyse information from the txt files created by pyfMRIqc ###
##################################################################


#### setups ####
#empty work space, load libraries and functions
rm(list=ls())

# define necessary directories
qualDir <- getwd()

# load libraries
library(dplyr)
library(ggplot2)
library(ggrepel)
library(grid)

# function to determine outliers
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

### downlaod data from OSF ###
version_official <- "fmri"
osfr::osf_auth() # log into OSF
project <- osfr::osf_retrieve_node("fhqb7")
target_dir <- osfr::osf_ls_files(project, pattern = "data") # looks at all files and directories in the project and defines the match with "data"
sub_dir <- osfr::osf_mkdir(target_dir, path = paste0(version_official)) # add folder in OSF data dir

osfr::osf_ls_files(sub_dir, pattern = "pyfMRIqc_output.csv") %>%
  osfr::osf_download(conflicts = "overwrite")

# load in dataset
scanparam <- read.csv("pyfMRIqc_output.csv")

# create ID column
scanparam$ID <- gsub("sub-control0", "", scanparam$subject)
scanparam$ID <- gsub("sub-experimental0", "", scanparam$ID)
#scanparam$ID <- as.numeric(scanparam$ID)

# remove scientific notation
scanparam$value <-format(scanparam$value, scientific = FALSE) 
scanparam$value <-as.numeric(scanparam$value) 

# make subject a character variable
scanparam$subject <- as.character(scanparam$subject)

# eliminate acq in BOLD
scanparam$acq <- gsub("_acq-1", "", scanparam$BOLD)
scanparam$acq <- gsub("_acq-2", "", scanparam$acq)

# determine DV to analyse
DV_plot <- c("SNR_voxel_MEAN", "SNR_voxel_STD", "Mean", "Mean_(mask)", "SD", "SD_(mask)", 
             "Min_Slice_SNR", "Max_Slice_SNR", "Mean_voxel_SNR", "Mean_absolute_Movement", "Max_absolute_Movement", 
             "Max_relative_Movement", "Relative_movements_(>0.1mm)", "Relative_movements_(>0.5mm)", "Relative_movements_(>voxelsize)",
             "volumes", "threshold_02", "fraction_02", "threshold_03", "fraction_03", "threshold_05", "fraction_05", "threshold_1", "fraction_1")

# determine font sizes
text_size <- 8
axis_text_size <- 12
axis_title_size <- 14
title_size <- 16
strip_text_size <- 14

#### create graphs for these DVs #### 
# a) for all scans --> histogram
# b) per task: magictrickwatching vs rest & per group: control & experimental --> histogram
# c) per acq: magictrickwatching_run-1 vs. magictrickwatching_run-2 vs. magictrickwatching_run-3 vs. rest_run-1 vs. rest_run-2
# a) per subject --> spaghetti plot


# loop over dependent variables to create plots 
for (DV in DV_plot){
  
  # create data frame and recode memory as factors
  print(DV)
  output <- subset(scanparam, scanparam$param == paste0(DV))
  
  # mean = red line; median = green line
  # Create a text
  Mean <- grobTree(textGrob(paste0("Mean = ", round(mean(output$value),2), " (SD = ", round(sd(output$value),2),")"), x=0.6,  y=0.9, hjust=0,
                            gp=gpar(col="red", fontsize=13, fontface="italic")))
  Median <- grobTree(textGrob(paste("Median =", round(median(output$value),2)), x=0.6,  y=0.85, hjust=0,
                              gp=gpar(col="green", fontsize=13, fontface="italic")))
  # (a) histogram for all 250 scans
  all <- ggplot(output, aes(x=value, col = group, fill = group)) + 
    geom_histogram(position = "stack", alpha=0.2, col = "black") +
    geom_vline(aes(xintercept = mean(value)),col='red', size=1) + 
    geom_vline(aes(xintercept = median(value)),col='green', size=1) + 
    theme_classic() +  scale_colour_grey() + scale_fill_grey() +
    labs(x=paste(DV), y="Frequency", title = "Histogram all scans") +
    theme(axis.text=element_text(size=axis_text_size), axis.title=element_text(size=axis_title_size, face="bold"), title=element_text(size =title_size, face="bold"), strip.text = element_text(size = strip_text_size))
  #print(all)
  all + annotation_custom(Mean) + annotation_custom(Median)
  ggsave(paste0("Hist_all_", DV, ".jpeg"))
  
  # (b) histogram splitted by task and group
  # mean = red line; median = green line
  task <- ggplot(output, aes(x=value)) + 
    theme_classic() +  #scale_colour_grey() +
    geom_histogram(position="dodge", alpha=0.2, col = "black") +
    geom_vline(data = plyr::ddply(output, c("task", "group"), summarize, mean = mean(value, na.rm = T)), aes(xintercept=mean), col='red') +
    geom_vline(data = plyr::ddply(output, "task", summarize, median = median(value, na.rm = T)), aes(xintercept=median), col='green') +
    labs(x=paste(DV), y="Frequency", title = "Histogram scans split by task and group") +
    theme(axis.text=element_text(size=axis_text_size), axis.title=element_text(size=axis_title_size, face="bold"), title=element_text(size =title_size, face="bold"), strip.text = element_text(size = strip_text_size))+ 
    facet_grid(group ~ task)
  # add mean and median as values
  xpos <- 0.5*(min(output$value) + max(output$value))
  ypos1 <- max(ggplot_build(task)$data[[1]]$count) - sd(ggplot_build(task)$data[[1]]$count)
  ypos2 <- ypos1 - 0.5*sd(ggplot_build(task)$data[[1]]$count)
  ypos3 <- ypos2 - 0.5*sd(ggplot_build(task)$data[[1]]$count)
  task + 
    geom_text(data = plyr::ddply(output, c("task", "group"), summarize, mean = round(mean(value, na.rm = T), 2)), 
              aes(label=paste("Mean =", mean), x = xpos, y = ypos1), vjust = 0, hjust = 0, col='red') +
    geom_text(data = plyr::ddply(output, c("task", "group"), summarize, sd = round(sd(value, na.rm = T), 2)), 
              aes(label=paste0("(SD = ", sd, ")"), x = xpos, y = ypos2), vjust = 0, hjust = 0, col='red') +
    geom_text(data = plyr::ddply(output, c("task", "group"), summarize, median = round(median(value, na.rm = T), 2)), 
              aes(label=paste("Median =", median), x = xpos, y = ypos3), vjust = 0, hjust = 0, col='green')
  ggsave(paste0("Hist_group_", DV, ".jpeg"))
  
  # create label for outliers in boxplot
  # outlier = x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x)
  output_outliers <- output %>% 
    group_by(acq) %>%
    mutate(outlier = ifelse(is_outlier(value), subject, NA))
  # (c) create boxplot for each acquisition
  boxplot <- ggplot(output_outliers, aes(y = value, x = acq, label = outlier)) +
    geom_boxplot()  + #geom_jitter(color = "red") +
    geom_label_repel(aes(label = outlier, color = group), nudge_x = 0.15, direction = "y", hjust = 0, segment.size = 0.2) +
    theme_classic() + scale_colour_grey() + theme(legend.position="none") +
    labs(y=paste(DV), x="BOLD acquisition", title = "Boxplot for each acquisition") +
    theme(axis.text=element_text(size=axis_text_size), axis.title=element_text(size=axis_title_size, face="bold"), title=element_text(size =title_size, face="bold"), strip.text = element_text(size = strip_text_size))
  print(boxplot)
  ggsave(paste0("Boxplot_", DV, ".jpeg"))
  
  # creeate label for spaghetti plot
  # subject gets labeled if  parameter > mean + 1SD
  output$above <- ifelse((mean(output$value, na.rm = T) + sd(output$value, na.rm = T)) < output$value, output$subject, NA) # add subject in column if value > mean + 1SD
  output$unique <- output$subject %in% unique(output$above) # determine unique subject names
  output$test <- ifelse(output$unique & output$acq == "rest_run-2", output$subject, NA) # add label to last acq only
  # (d) create spaghetti plot for each subject
  subj <- ggplot(output, aes(x=acq, y = value, group=subject, col = ID)) + geom_line() +
    facet_grid(group ~ .) +
    theme_classic() + theme(legend.position="none") +
    scale_x_discrete(limits=c("rest_run-1", "magictrickwatching_run-1", "magictrickwatching_run-2", "magictrickwatching_run-3", "rest_run-2")) +
    labs(y=paste(DV), x="BOLD acquisition", title = "Parameters over time per subject") +
    theme(axis.text=element_text(size=axis_text_size), axis.title=element_text(size=axis_title_size, face="bold"), title=element_text(size =title_size, face="bold"), strip.text = element_text(size = strip_text_size))+
    geom_label_repel(aes(label = test), nudge_x = 1, na.rm = TRUE) 
  print(subj)
  ggsave(paste0("Spaghettiplot_", DV, ".jpeg"))
  
}



