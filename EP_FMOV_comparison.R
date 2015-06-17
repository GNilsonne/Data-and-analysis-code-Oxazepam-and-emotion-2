# Script to analyse correlations between outcomes from different experiments in the 2011 Oxazepam study
# Gustav Nilsonne 2015-06-16

# Require packages
require(RCurl) # To read data from GitHub

# Define functions
# panel.cor puts correlation in upper panels of pairs plots, size proportional to correlation
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, use = "pairwise.complete.obs")
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

# Read files
# First demographic data
demDataURL <- getURL("https://raw.githubusercontent.com/GNilsonne/Data-and-analysis-code-Oxazepam-and-emotion/master/demographics.csv", ssl.verifypeer = FALSE)
demData <- read.csv(text = demDataURL)

# Then results files
ep_ratings_url <- getURL("https://raw.githubusercontent.com/GNilsonne/Data-and-analysis-code-Oxazepam-and-emotion/master/IndividualResponsesUnpleasantness.csv", ssl.verifypeer = FALSE)
ep_ratings <- read.csv(text = ep_ratings_url)

ep_scr_url <- getURL("https://raw.githubusercontent.com/GNilsonne/Data-and-analysis-code-Oxazepam-and-emotion/master/IndividualResponsesSCR.csv", ssl.verifypeer = FALSE)
ep_scr <- read.csv(text = ep_scr_url)

ep_emg_url <- getURL("https://raw.githubusercontent.com/GNilsonne/Data-and-analysis-code-Oxazepam-and-emotion/master/IndividualResponsesEMG.csv", ssl.verifypeer = FALSE)
ep_emg <- read.csv(text = ep_emg_url)

ep_hr_url <- getURL("https://raw.githubusercontent.com/GNilsonne/Data-and-analysis-code-Oxazepam-and-emotion/master/IndividualResponsesHeartRate.csv", ssl.verifypeer = FALSE)
ep_hr <- read.csv(text = ep_hr_url)

fmov_emg_url <- getURL("https://raw.githubusercontent.com/GNilsonne/Data-and-analysis-code-Oxazepam-and-emotion-2/master/IndividualResponsesFMOV.csv", ssl.verifypeer = FALSE)
fmov_emg <- read.csv(text = fmov_emg_url)

# Analyse covariance within EP experiment
# Merge data files
ep_ind <- merge(ep_ratings, ep_scr, by = "Subject", all = T)
ep_ind <- merge(ep_ind, ep_emg, by = "Subject", all = T)
ep_ind <- merge(ep_ind, ep_hr, by = "Subject", all = T)

# Inspect covariance
pairs(~ Unpleasantness + sqrtSCR + EMG_corr_mean + hr_mean, data = ep_ind,
      lower.panel = panel.smooth, upper.panel = panel.cor, 
      pch = 20, main = "EP Other High")

# Analyse covariance within FMOV experiment
# Inspect covariance
pairs(~ CorrAngry + CorrHappy + ZygAngry + ZygHappy, data = fmov_emg,
      lower.panel = panel.smooth, upper.panel = panel.cor, 
      pch = 20, main = "FMOV")

# Analyse covariance between FMOV and EP experiments
# I choose to make one compound variable for FMOV responses based on corrugator and zygomatic responses to happy faces, since that was where we found main effects and since they correlated inversely, as expected.
fmov_emg$fmov_compound_z <- scale(fmov_emg$ZygHappy - fmov_emg$CorrHappy)
ind <- merge(ep_ind, fmov_emg[, c("Subject", "CorrHappy", "ZygHappy", "fmov_compound_z")], by = "Subject", all = T)

pairs(~ Unpleasantness + sqrtSCR + EMG_corr_mean + hr_mean + fmov_compound_z, data = ind,
      lower.panel = panel.smooth, upper.panel = panel.cor, 
      pch = 20, main = "EP Other High and FMOV")
