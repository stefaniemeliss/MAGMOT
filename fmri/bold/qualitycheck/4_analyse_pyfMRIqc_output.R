##################################################################
### analyse information from the txt files created by pyfMRIqc ###
##################################################################


#### setups ####
#empty work space, load libraries and functions
rm(list=ls())

# define necessary directories
qualDir <- getwd()

library(dplyr)

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

# remove scientific notation
scanparam$value <-format(scanparam$value, scientific = FALSE) 
scanparam$value <-as.numeric(scanparam$value) 

# determine DV to analyse
DV_plot <- c("SNR_voxel_MEAN", "SNR_voxel_STD", "Mean", "Mean_(mask)", "SD_(mask)", 
             "Min_Slice_SNR", "Max_Slice_SNR", "Mean_voxel_SNR", "Mean_absolute_Movement", "Max_absolute_Movement", 
             "Max_relative_Movement", "Relative_movements_(>0.1mm)", "Relative_movements_(>0.5mm)", "Relative_movements_(>voxelsize)") 

# create graphs for these DVs
# a) for all scans --> histogram
# b) per task: magictrickwatching vs rest --> violin plot
# c) per BOLD: magictrickwatching_run-1 vs. magictrickwatching_run-2 vs. magictrickwatching_run-3 vs. rest_run-1 vs. rest_run-2
# a) per subject --> spaghetti plot


# loop over dependent variables to create histograms 
for (DV in DV_plot){
  # create data frame and recode memory as factors
  output <- subset(scanparam, scanparam$param == paste0(DV))
  
}
  
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
    labs(x="curiosity group mean centered", y="Density", title = paste(version_official, DV)) +
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
    labs(x="curiosity group mean centered", y="Density", title = paste(version_official, DV)) +
    theme(legend.position="bottom") + theme(legend.title = element_blank()) + 
    theme(axis.text=element_text(size=14), axis.title=element_text(size=16, face="bold"), title=element_text(size =20, face="bold"), strip.text = element_text(size = 16)) 
  
  print(graph)
  print(paste0("Histogram_", DV, "_2.jpeg"))
  ggsave(paste0("Histogram_", DV, "_2.jpeg"))
}

