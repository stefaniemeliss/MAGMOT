# note: computed jitter has the contraint of having 2 seconds wait at the begin of each run, but no waiting time at the end of each run

rm(list=ls()) # empty workspace

# define working directory
setwd('~/Dropbox/Reading/PhD/Magictricks/fmri_study/jitter/stim_results/')
summary <- read.delim("NSD_sums", header = F, sep = " ")
summary[which.min(summary$V1),]

iter <- 0222

# basic inputs for duration of video clips, ratings and runs
durVideos <- 38
durRatings <- 12
durRun <- 720

# read in the file for the onset timings of the videos
videosOnsets <- read.delim(paste0('stimes.0', iter, '_01_video.1D'), header = F, sep = " ")
videosOnsets$V13 <- NULL # delete V13
videosOnsets <- t(videosOnsets) #transpose
videosOnsets <- data.frame(videosOnsets) # create data frame
names(videosOnsets) <- c("run1", "run2", "run3")
videosOffsets <- videosOnsets + durVideos # determine when videos end

# create identifier
videosOnsets$identifier <- NA 
for (i in 1:nrow(videosOnsets)) {
  videosOnsets[i,"identifier"] <- c(paste0("videosOnset_stim", i))
}
videosOffsets$identifier <- NA
for (i in 1:nrow(videosOffsets)) {
  videosOffsets[i,"identifier"] <- c(paste0("videosOffset_stim", i))
}

# create a data frame containing both, onsets and offsets of the videos
videos <- rbind(videosOnsets, videosOffsets)
videos$stimType  <- "video"
videos$trial  <- c(rep(c(1:12),2))

# read in the file for the onset timings of the videos
ratingsOnsets <- read.delim(paste0('stimes.0', iter, '_02_rating.1D'), header = F, sep = " ")
ratingsOnsets$V13 <- NULL
ratingsOnsets <- t(ratingsOnsets)
ratingsOnsets <- data.frame(ratingsOnsets)
names(ratingsOnsets) <- c("run1", "run2", "run3")
ratingsOffsets <- ratingsOnsets + durRatings # determine when ratings end

# create identifier
ratingsOnsets$identifier <- NA
for (i in 1:nrow(ratingsOnsets)) {
  ratingsOnsets[i,"identifier"] <- c(paste0("ratingsOnset_stim", i))
}
ratingsOffsets$identifier <- NA
for (i in 1:nrow(ratingsOffsets)) {
  ratingsOffsets[i,"identifier"] <- c(paste0("ratingsOffset_stim", i))
}

# create a data frame containing both, onsets and offsets of the videos
ratings <- rbind(ratingsOnsets, ratingsOffsets)
ratings$stimType  <- "rating"
ratings$trial  <- c(rep(c(1:12),2))

# combine information for videos and rating
allTimings <- rbind(videos, ratings)
allTimings <- allTimings[order(allTimings$run1),] # sort it by ascending order of onsets

# compute ISI for all three runs: how much time is in between offset of event i and onset of event i+1?
allTimings$run1_ISI <- NA
for (i in 1:nrow(allTimings)){
  allTimings[i,"run1_ISI"] <- allTimings[i+1,"run1"] - allTimings[i,"run1"]
  if (i == nrow(allTimings)){
    allTimings[i,"run1_ISI"] <- durRun - allTimings[i,"run1"]
  } 
}
allTimings$run2_ISI <- NA
for (i in 1:nrow(allTimings)){
  allTimings[i,"run2_ISI"] <- allTimings[i+1,"run2"] - allTimings[i,"run2"]
  if (i == nrow(allTimings)){
    allTimings[i,"run2_ISI"] <- durRun - allTimings[i,"run2"]
  } 
}
allTimings$run3_ISI <- NA
for (i in 1:nrow(allTimings)){
  allTimings[i,"run3_ISI"] <- allTimings[i+1,"run3"] - allTimings[i,"run3"]
  if (i == nrow(allTimings)){
    allTimings[i,"run3_ISI"] <- durRun - allTimings[i,"run3"]
  } }


# identify jitter: jitter refers to the timing between offset video and onset rating as well as between offset rating and onset video
allTimings$run1_jitter <- ifelse(allTimings$run1_ISI == durVideos & allTimings$stimType == "video", NA, #set to NA if it reflects the duration of the video
                                 ifelse(allTimings$run1_ISI == durRatings & allTimings$stimType == "rating", NA, allTimings$run1_ISI)) #or the duration of the rating
allTimings$run2_jitter <- ifelse(allTimings$run2_ISI == durVideos & allTimings$stimType == "video", NA,
                                 ifelse(allTimings$run2_ISI == durRatings & allTimings$stimType == "rating", NA, allTimings$run2_ISI))
allTimings$run3_jitter <- ifelse(allTimings$run3_ISI == durVideos & allTimings$stimType == "video", NA,
                                 ifelse(allTimings$run3_ISI == durRatings & allTimings$stimType == "rating", NA, allTimings$run3_ISI))

setwd('~/Dropbox/Reading/PhD/Magictricks/fmri_study/PsychToolBox_script/')

videosJitter <- subset(allTimings, allTimings$stimType == "video") #extract timings for video only
videosJitter <- subset(videosJitter, is.na(videosJitter$run1_jitter)==F) # and only those referring to jitter rather than ISIs (= durVid or durRatings)
videoJitter <- c(videosJitter$run1_jitter, videosJitter$run2_jitter,videosJitter$run3_jitter) # combine all three runs into one long vector
videoJitter
hist(videoJitter, ylim = c(0,36), xlim = c(0, 12))
write(videoJitter, file = "videoJitter.tsv", ncolumns = 1, sep = "\t")

#and now do the same for the ratings
ratingsJitter <- subset(allTimings, allTimings$stimType == "rating")
ratingsJitter <- subset(ratingsJitter, is.na(ratingsJitter$run1_jitter)==F)
ratingJitter <- c(ratingsJitter$run1_jitter, ratingsJitter$run2_jitter,ratingsJitter$run3_jitter)
ratingJitter
hist(ratingJitter, ylim = c(0,36), xlim = c(0, 12))
write(ratingJitter, file = "ratingJitter.tsv", ncolumns = 1, sep = "\t")

