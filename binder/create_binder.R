# this code creates the binder for the MAGMOT R analysis

# 1. install packages
remotes::install_github("karthik/holepunch")

# 2. Setting up your project as a compendium
library(holepunch)
write_compendium_description(package = "Investigating the effects of monetary reward and curiosity on incidental memory encoding", 
                             description = "In this project, we examine the effects of monetary reward on incidental learning and how this is influenced by experiencing curiosity. This includes behavioural experiments conducted online as well as an fMRI study to investigate the neural underpinnings of the effects. We use dynamic stimuli (e.g. short movie clips) to increase the ecological validity of our research.")
# to write a description, with dependencies. Be sure to fill in placeholder text

write_dockerfile(maintainer = "Stef Meliss") 
# To write a Dockerfile. It will automatically pick the date of the last 
# modified file, match it to that version of R and add it here. You can 
# override this by passing r_date to some arbitrary date
# (but one for which a R version exists).

generate_badge() # This generates a badge for your readme.

# ----------------------------------------------
# At this time 🙌 push the code to GitHub 🙌
# ----------------------------------------------

# And click on the badge or use the function below to get the build 
# ready ahead of time.
build_binder()