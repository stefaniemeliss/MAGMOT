# this script downloads files from OSF storage.
# the files are needed as input for the Python code that cuts the fMRI data.

library(dplyr)

version_official <- "fmri"
git_dir <- "~/Dropbox/Reading/PhD/Magictricks/git" # this needs to be changed to reflect the directory where the git repository has been downloaded to
setwd(file.path(git_dir, "fmri", "bold", "concat")) # please use the same directory here where the jupyter notebook is saved in!!

osfr::osf_auth() # log into OSF
project <- osfr::osf_retrieve_node("fhqb7")
target_dir <- osfr::osf_ls_files(project, pattern = "data") # looks at all files and directories in the project and defines the match with "data"
sub_dir <- osfr::osf_mkdir(target_dir, path = paste0(version_official)) # add folder in OSF data dir

# download file from OSF
osfr::osf_ls_files(sub_dir, pattern = "informationAboutScanDuration") %>%
  osfr::osf_download(conflicts = "overwrite")

