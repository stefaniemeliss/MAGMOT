### this script contains the code used for the dataset paper 

#### setups ####
#empty work space, load libraries and functions
rm(list=ls())

# define necessary directories
analysisDir <- getwd()
ratingsDir <- file.path(analysisDir, "ratings")
memoryDir <- file.path(analysisDir, "memory")
tricksDir <- file.path(analysisDir, "tricks")
brainDir <- file.path(analysisDir, "brain")

setwd('..') # moves up in relative path
rootDir <- getwd()
pooledDir <-  file.path(rootDir, "pooled")
setwd(analysisDir)

# check whether these directories exist, if not create them
ifelse(!dir.exists(ratingsDir), dir.create(ratingsDir), FALSE) 
ifelse(!dir.exists(memoryDir), dir.create(memoryDir), FALSE) 
ifelse(!dir.exists(tricksDir), dir.create(tricksDir), FALSE) 
ifelse(!dir.exists(brainDir), dir.create(brainDir), FALSE) 
ifelse(!dir.exists(pooledDir), dir.create(pooledDir), FALSE) 

# helper functions and packages #
devtools::source_url("https://github.com/stefaniemeliss/MAGMOT/blob/master/functions/errorbars.R?raw=TRUE")
devtools::source_url("https://github.com/stefaniemeliss/MAGMOT/blob/master/functions/rbindcolumns.R?raw=TRUE")

library(lmerTest)
library(psych)
library(ggplot2)
library(dplyr)
library(reshape2)

# define version 
version <- "MAGMOT"
version_official <- "fmri"

memoryLevels <- c("cuedRecallStrict", "cuedRecallLenient", 
                  "recognition", "recognitionConfLevel_4_5_6", "recognitionAboveMeanConf", 
                  "rememberedStrictAboveAvg", "rememberedLenientAboveAvg", "rememberedStrictHigh", "rememberedLenientHigh")

memoryLabels <- c("cuedRecallStrict", "cuedRecallLenient", 
                  "allConf", "highConf", "aboveAvgConf", 
                  "rememberedStrictAboveAvg", "rememberedLenientAboveAvg", "rememberedStrictHigh", "rememberedLenientHigh")

pooled <- 1

dataCoded <- 1




### downlaod data sets from OSF ###
osfr::osf_auth() # log into OSF
project <- osfr::osf_retrieve_node("fhqb7")
target_dir <- osfr::osf_ls_files(project, pattern = "data") # looks at all files and directories in the project and defines the match with "data"
sub_dir <- osfr::osf_mkdir(target_dir, path = paste0(version_official)) # add folder in OSF data dir

osfr::osf_ls_files(sub_dir, pattern = ".xlsx") %>%
  osfr::osf_download(conflicts = "overwrite")

osfr::osf_ls_files(sub_dir, pattern = ".txt") %>%
  osfr::osf_download(conflicts = "overwrite")

### read in data sets ###
dfWide <- xlsx::read.xlsx(paste0("wide_MagicBehavioural_", version_official, ".xlsx"), sheetName = "Sheet1")
dfLong <- xlsx::read.xlsx(paste0("long_MagicBehavioural_", version_official, ".xlsx"), sheetName = "Sheet1")
dfROI <- read.delim("ISC_ROI.txt")