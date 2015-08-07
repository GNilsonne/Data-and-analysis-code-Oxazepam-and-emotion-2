# Script to analyse EMG in the FMOV experiment of the Oxazepam project 2011
# Gustav Nilsonne 2015-06-08

# Require packages
require(RCurl) # To read data from GitHub
require(reshape2)
require(RColorBrewer)
require(nlme)
require(effects)
require(latticeExtra) # To make dotcharts
require(parallel) # To increase processing speed

# Define colors for later
col1 = brewer.pal(8, "Dark2")[1]
col2 = brewer.pal(8, "Dark2")[2]
col3 = brewer.pal(8, "Dark2")[3]
col4 = brewer.pal(8, "Dark2")[4]
col5 = brewer.pal(8, "Dark2")[5]
col6 = brewer.pal(8, "Dark2")[6]
col7 = brewer.pal(8, "Dark2")[7]
col8 = brewer.pal(8, "Dark2")[8]

# Define function to read data
fun_readEMGdata <- function(X){
  
  require(RCurl) # To read data from GitHub
  
  # Read data file
  if(X == 9){
    X <- "09"
  }
  EMGDataURL <- getURL(paste("https://raw.githubusercontent.com/GNilsonne/Data-and-analysis-code-Oxazepam-and-emotion-2/master/FMOV/", X, "_FMOV_preprocessed.txt", sep = ""), ssl.verifypeer = FALSE)
  data <- read.table(text = EMGDataURL, header = FALSE)
  if(length(data) == 8){
    names(data) <- c("SCR", "Raw_EMG_Zyg", "Raw_EMG_corr", "EMG_zyg", "EMG_corr", "Angry", "Happy", "Neutral")
  } else if(length(data) == 10){
    names(data) <- c("Time", "SCR", "Raw_EMG_Zyg", "Raw_EMG_corr", "EMG_zyg", "EKG", "EMG_corr", "Angry", "Happy", "Neutral")
  }
  
  # Verify that data are reasonably orthogonal
  plot(data$EMG_zyg ~ data$EMG_corr, frame.plot = F, main = X)
  
  # Find event onsets
  data$AngryDiff <- c(0, diff(data$Angry, lag = 1))
  data$AngryOnsets <- data$AngryDiff == -5
  data$HappyDiff <- c(0, diff(data$Happy, lag = 1))
  data$HappyOnsets <- data$HappyDiff == -5
  data$NeutralDiff <- c(0, diff(data$Neutral, lag = 1))
  data$NeutralOnsets <- data$NeutralDiff == -5
  
  # Add time vector. Sampling was at 10 millisecond intervals.
  data$millisecond <- 10*(1:length(data$SCR))
  
  # Get EMG data
  # Define data for each event  
  AngryIndex <- which(data$AngryOnsets == TRUE)
  data$AngryWindow <- FALSE
  for(i in AngryIndex){
    if(length(data[, 1]) > (i + 600)){
      data$AngryWindow[(i - 200):(i + 600)] <- TRUE
    } else {
      this_length <- length(data[, 1]) - i
      data$AngryWindow[(i - 200):(i + this_length)] <- TRUE
    }
  }
  HappyIndex <- which(data$HappyOnsets == TRUE)
  data$HappyWindow <- FALSE
  for(i in HappyIndex){
    if(length(data[, 1]) > (i + 600)){
      data$HappyWindow[(i - 200):(i + 600)] <- TRUE
    } else {
      this_length <- length(data[, 1]) - i
      data$HappyWindow[(i - 200):(i + this_length)] <- TRUE
    }
  }
  NeutralIndex <- which(data$NeutralOnsets == TRUE)
  data$NeutralWindow <- FALSE
  for(i in NeutralIndex){
    if(length(data[, 1]) > (i + 600)){
      data$NeutralWindow[(i - 200):(i + 600)] <- TRUE
    } else {
      this_length <- length(data[, 1]) - i
      data$NeutralWindow[(i - 200):(i + this_length)] <- TRUE
    }
  }
  
  # Extract data, downsample
  AngryEMG_corr <- matrix(data$EMG_corr[data$AngryWindow == TRUE], ncol = length(AngryIndex))
  AngryEMG_corrlowess <- apply(AngryEMG_corr, 2, lowess, f = 0.01)
  AngryEMG_corrlowess <- matrix(unlist(AngryEMG_corrlowess), byrow = F, ncol = length(AngryEMG_corrlowess))[802:1602, ]
  
  plot(AngryEMG_corr[, 1], type = "l", col = "gray", frame.plot = F, main = X)
  lines(AngryEMG_corrlowess[, 1])
  lines(y = AngryEMG_corrlowess[, 1][seq(1, NROW(AngryEMG_corrlowess[, 1]), by=10)], x = c(1:801)[seq(1, NROW(AngryEMG_corrlowess[, 1]), by=10)], col = "red", type = "o")
  
  AngryEMG_corrdownsampled <- AngryEMG_corr[seq(1, NROW(AngryEMG_corrlowess), by=10), ]
  
  AngryEMG_zyg <- matrix(data$EMG_zyg[data$AngryWindow == TRUE], ncol = length(AngryIndex))
  AngryEMG_zyglowess <- apply(AngryEMG_zyg, 2, lowess, f = 0.01)
  AngryEMG_zyglowess <- matrix(unlist(AngryEMG_zyglowess), byrow = F, ncol = length(AngryEMG_zyglowess))[802:1602, ]
  AngryEMG_zygdownsampled <- AngryEMG_zyg[seq(1, NROW(AngryEMG_zyglowess), by=10), ]
  
  HappyEMG_corr <- matrix(data$EMG_corr[data$HappyWindow == TRUE], ncol = length(HappyIndex))
  HappyEMG_corrlowess <- apply(HappyEMG_corr, 2, lowess, f = 0.01)
  HappyEMG_corrlowess <- matrix(unlist(HappyEMG_corrlowess), byrow = F, ncol = length(HappyEMG_corrlowess))[802:1602, ]
  HappyEMG_corrdownsampled <- HappyEMG_corr[seq(1, NROW(HappyEMG_corrlowess), by=10), ]
  
  HappyEMG_zyg <- matrix(data$EMG_zyg[data$HappyWindow == TRUE], ncol = length(HappyIndex))
  HappyEMG_zyglowess <- apply(HappyEMG_zyg, 2, lowess, f = 0.01)
  HappyEMG_zyglowess <- matrix(unlist(HappyEMG_zyglowess), byrow = F, ncol = length(HappyEMG_zyglowess))[802:1602, ]
  HappyEMG_zygdownsampled <- HappyEMG_zyg[seq(1, NROW(HappyEMG_zyglowess), by=10), ]
  
  NeutralEMG_corr <- matrix(data$EMG_corr[data$NeutralWindow == TRUE], ncol = length(NeutralIndex))
  NeutralEMG_corrlowess <- apply(NeutralEMG_corr, 2, lowess, f = 0.01)
  NeutralEMG_corrlowess <- matrix(unlist(NeutralEMG_corrlowess), byrow = F, ncol = length(NeutralEMG_corrlowess))[802:1602, ]
  NeutralEMG_corrdownsampled <- NeutralEMG_corr[seq(1, NROW(NeutralEMG_corrlowess), by=10), ]
  
  NeutralEMG_zyg <- matrix(data$EMG_zyg[data$NeutralWindow == TRUE], ncol = length(NeutralIndex))
  NeutralEMG_zyglowess <- apply(NeutralEMG_zyg, 2, lowess, f = 0.01)
  NeutralEMG_zyglowess <- matrix(unlist(NeutralEMG_zyglowess), byrow = F, ncol = length(NeutralEMG_zyglowess))[802:1602, ]
  NeutralEMG_zygdownsampled <- NeutralEMG_zyg[seq(1, NROW(NeutralEMG_zyglowess), by=10), ]
  
  # Make a list for each participant and return it as output
  EMGData <- list(as.integer(X), AngryEMG_zygdownsampled, AngryEMG_corrdownsampled, HappyEMG_zygdownsampled, HappyEMG_corrdownsampled, NeutralEMG_zygdownsampled, NeutralEMG_corrdownsampled)
  return(EMGData)
}

# Read files
# First demographic data etc
demDataURL <- getURL("https://raw.githubusercontent.com/GNilsonne/Data-and-analysis-code-Oxazepam-and-emotion/master/demographics.csv", ssl.verifypeer = FALSE)
demData <- read.csv(text = demDataURL)
# Then the Acqknowledge logfiles containing EMG data
IncludedSubjects <- demData$Subject[demData$Included_FMOV == T]
SubjectsWave1 <- demData$Subject[demData$Included_FMOV == T & demData$Wave == 1]
SubjectsWave2 <- demData$Subject[demData$Included_FMOV == T & demData$Wave == 2]
#EMGDataList <- lapply(IncludedSubjects, FUN = fun_readEMGdata)
# Read Acqknowledge datafiles with parallelization to save time
no_cores <- detectCores() - 1 # Determine number of cores
cl <- makeCluster(no_cores) # Initiate cluster
EMGDataList <- parLapply(cl, IncludedSubjects, fun = fun_readEMGdata)
stopCluster(cl) # Stop cluster

# Analyse data, starting with corrugator activity
# Move all data to one big frame
EMG_corr <- data.frame()

for(i in 1:length(EMGDataList)){
  if(!is.null(EMGDataList[i][[1]])){
    Angry_corr <- melt(EMGDataList[i][[1]][[3]], na.rm = F)
    names(Angry_corr) <- c("time_0.1s", "event_no", "EMG_corr")
    Angry_corr$Stimulus <- "Angry"
    
    Happy_corr <- melt(EMGDataList[i][[1]][[5]], na.rm = F)
    names(Happy_corr) <- c("time_0.1s", "event_no", "EMG_corr")
    Happy_corr$Stimulus <- "Happy"
    
    Neutral_corr <- melt(EMGDataList[i][[1]][[7]], na.rm = F)
    names(Neutral_corr) <- c("time_0.1s", "event_no", "EMG_corr")
    Neutral_corr$Stimulus <- "Neutral"
    
    EMG_corr_temp <- rbind(Angry_corr, Happy_corr, Neutral_corr)
    EMG_corr_temp$subject <- IncludedSubjects[i]
  }
  EMG_corr <- rbind(EMG_corr, EMG_corr_temp)
}

# Convert time scale and set 0 to event onset
EMG_corr$time_s <- EMG_corr$time_0.1s/10
EMG_corr$time_s <- EMG_corr$time_s -2

# Make spaghetti plots
plot(EMG_corr ~ time_s, data = subset(EMG_corr, Stimulus == "Angry"), frame.plot = F, type = 'n', main = "Angry", ylim = c(0, 0.004))
for(i in unique(EMG_corr$subject)){
  for( j in unique(EMG_corr$event_no[EMG_corr$subject == i])){
    lines(EMG_corr ~ time_s, data = EMG_corr[EMG_corr$subject == i & EMG_corr$event_no == j & EMG_corr$Stimulus == "Angry", ], col = col4)
  }
}

plot(EMG_corr ~ time_s, data = subset(EMG_corr, Stimulus == "Happy"), frame.plot = F, type = 'n', main = "Happy", ylim = c(0, 0.004))
for(i in unique(EMG_corr$subject)){
  for( j in unique(EMG_corr$event_no[EMG_corr$subject == i])){
    lines(EMG_corr ~ time_s, data = EMG_corr[EMG_corr$subject == i & EMG_corr$event_no == j & EMG_corr$Stimulus == "Happy", ], col = col7)
  }
}

# Since data have varying baselines between conditions, we index each response to its baseline from -2 to 0 seconds
EMG_corr$EMG_corr_index <- NA
for(i in unique(EMG_corr$subject)){
  for(j in unique(EMG_corr$event_no[EMG_corr$subject == i])){
    for(k in c("Angry", "Happy", "Neutral")){
      EMG_corr$EMG_corr_index[EMG_corr$subject == i & EMG_corr$event_no == j & EMG_corr$Stimulus == k] <- EMG_corr$EMG_corr[EMG_corr$subject == i & EMG_corr$event_no == j & EMG_corr$Stimulus == k]/mean(EMG_corr$EMG_corr[EMG_corr$subject == i & EMG_corr$event_no == j & EMG_corr$Stimulus == k][1:20])
    }
  }
}

# Make spaghetti plots of normalised data
plot(EMG_corr_index ~ time_s, data = subset(EMG_corr, Stimulus == "Angry"), frame.plot = F, type = 'n', main = "Angry", ylim = c(0, 15))
for(i in unique(EMG_corr$subject)){
  for(j in unique(EMG_corr$event_no[EMG_corr$subject == i])){
    lines(EMG_corr_index ~ time_s, data = EMG_corr[EMG_corr$subject == i & EMG_corr$event_no == j & EMG_corr$Stimulus == "Angry", ], col = col4)
  }
}

plot(EMG_corr_index ~ time_s, data = subset(EMG_corr, Stimulus == "Happy"), frame.plot = F, type = 'n', main = "Happy", ylim = c(0, 15))
for(i in unique(EMG_corr$subject)){
  for(j in unique(EMG_corr$event_no[EMG_corr$subject == i])){
    lines(EMG_corr_index ~ time_s, data = EMG_corr[EMG_corr$subject == i & EMG_corr$event_no == j & EMG_corr$Stimulus == "Happy", ], col = col7)
  }
}

plot(EMG_corr_index ~ time_s, data = subset(EMG_corr, Stimulus == "Neutral"), frame.plot = F, type = 'n', main = "Neutral", ylim = c(0, 15))
for(i in unique(EMG_corr$subject)){
  for(j in unique(EMG_corr$event_no[EMG_corr$subject == i])){
    lines(EMG_corr_index ~ time_s, data = EMG_corr[EMG_corr$subject == i & EMG_corr$event_no == j & EMG_corr$Stimulus == "Neutral", ], col = col8)
  }
}

MeanEMGCorrAngry <- aggregate(EMG_corr_index ~ time_s, data = subset(EMG_corr, Stimulus == "Angry"), mean)
MeanEMGCorrHappy <- aggregate(EMG_corr_index ~ time_s, data = subset(EMG_corr, Stimulus == "Happy"), mean)
MeanEMGCorrNeutral <- aggregate(EMG_corr_index ~ time_s, data = subset(EMG_corr, Stimulus == "Neutral"), mean)

MeanEMGCorrAngryWave1 <- aggregate(EMG_corr_index ~ time_s, data = subset(EMG_corr, Stimulus == "Angry" & subject %in% SubjectsWave1), mean)
MeanEMGCorrHappyWave1 <- aggregate(EMG_corr_index ~ time_s, data = subset(EMG_corr, Stimulus == "Happy" & subject %in% SubjectsWave1), mean)
MeanEMGCorrNeutralWave1 <- aggregate(EMG_corr_index ~ time_s, data = subset(EMG_corr, Stimulus == "Neutral" & subject %in% SubjectsWave1), mean)

MeanEMGCorrAngryWave2 <- aggregate(EMG_corr_index ~ time_s, data = subset(EMG_corr, Stimulus == "Angry" & subject %in% SubjectsWave2), mean)
MeanEMGCorrHappyWave2 <- aggregate(EMG_corr_index ~ time_s, data = subset(EMG_corr, Stimulus == "Happy" & subject %in% SubjectsWave2), mean)
MeanEMGCorrNeutralWave2 <- aggregate(EMG_corr_index ~ time_s, data = subset(EMG_corr, Stimulus == "Neutral" & subject %in% SubjectsWave2), mean)


plot(MeanEMGCorrAngry, type = "l", col = col4, frame.plot = F, main = "Corrugator EMG", ylim = c(0.8, 1.2))
lines(MeanEMGCorrHappy, col = col7)
lines(MeanEMGCorrNeutral, col = col8)
legend("topleft", lty = 1, col = c(col4, col7, col8), legend = c("Angry", "Happy", "Neutral"), bty = "n")
abline(v = c(0), lty = 2)

plot(MeanEMGCorrAngry, type = "n", col = col4, frame.plot = F, main = "Corrugator EMG", ylim = c(0.8, 1.2))
lines(lowess(MeanEMGCorrAngry, f = 0.1), col = col4)
lines(lowess(MeanEMGCorrHappy, f = 0.1), col = col7)
lines(lowess(MeanEMGCorrNeutral, f = 0.1), col = col8)
legend("topleft", lty = 1, col = c(col4, col7, col8), legend = c("Angry", "Happy", "Neutral"), bty = "n")
abline(v = c(0), lty = 2)


# Make figures
pdf("Fig_EMG1.pdf", width = 7, height = 5)
plot(MeanEMGCorrAngryWave1, type = "n", col = col4, frame.plot = F, main = "Corrugator EMG, Wave 1", xlab = "Time (s)", ylab = "EMG (ratio)", xaxt = "n", yaxt = "n", ylim = c(0.8, 1.3))
rect(2, 0.000001, 4, 10, col = "gray88", border = NA)
abline(v = c(0, 2, 5), lty = 3)
axis(1)
axis(2, at = c(0.8, 1, 1.2))
lines(lowess(MeanEMGCorrAngryWave1, f = 0.1), col = col4, lwd = 2)
lines(lowess(MeanEMGCorrHappyWave1, f = 0.1), col = col7, lwd = 2)
lines(lowess(MeanEMGCorrNeutralWave1, f= 0.1), col = col8, lwd = 2)
legend("topleft", lty = 1, col = c(col4, col7, col8), legend = c("Angry", "Happy", "Neutral"), lwd = 2, bty = "n")
dev.off()

pdf("Fig_EMG2.pdf", width = 7, height = 5)
plot(MeanEMGCorrAngryWave2, type = "n", col = col4, frame.plot = F, main = "Corrugator EMG, Wave 2", xlab = "Time (s)", ylab = "EMG (ratio)", xaxt = "n", yaxt = "n", ylim = c(0.8, 1.3))
rect(2, 0.000001, 4, 10, col = "gray88", border = NA)
abline(v = c(0, 2, 4), lty = 3)
axis(1)
axis(2, at = c(0.8, 1, 1.2))
lines(lowess(MeanEMGCorrAngryWave2, f = 0.1), col = col4, lwd = 2)
lines(lowess(MeanEMGCorrHappyWave2, f = 0.1), col = col7, lwd = 2)
lines(lowess(MeanEMGCorrNeutralWave2, f = 0.1), col = col8, lwd = 2)
dev.off()

# Average data over interval for each event for statistical modelling
EMGEventData <- data.frame()
for(i in unique(EMG_corr$subject)){
  for(j in unique(EMG_corr$event_no[EMG_corr$subject == i])){
    for(k in c("Angry", "Happy", "Neutral")){
      if(i %in% demData$Subject[demData$Wave == 1]){
        EMG_corr_mean <- mean(EMG_corr$EMG_corr_index[EMG_corr$subject == i & EMG_corr$event_no == j & EMG_corr$Stimulus == k][40:60])
        EMGEventData <- rbind(EMGEventData, data.frame(i, j, k, EMG_corr_mean))
      } else if(i %in% demData$Subject[demData$Wave == 2]){
        EMG_corr_mean <- mean(EMG_corr$EMG_corr_index[EMG_corr$subject == i & EMG_corr$event_no == j & EMG_corr$Stimulus == k][40:60])
        EMGEventData <- rbind(EMGEventData, data.frame(i, j, k, EMG_corr_mean))
      }
    }
  }
}
names(EMGEventData) <- c("Subject", "event_no", "Stimulus", "EMG_corr_mean")

hist(EMGEventData$EMG_corr_mean)
hist(log(EMGEventData$EMG_corr_mean))

EMGEventData$EMG_corr_mean <- log(EMGEventData$EMG_corr_mean)

# Analyse data
EMGEventData <- merge(EMGEventData, demData, by = "Subject")
EMGEventData$Subject <- as.factor(EMGEventData$Subject)
EMGEventData$Stimulus <- relevel(EMGEventData$Stimulus, ref = "Neutral")
EMGEventData$Treatment <- relevel(EMGEventData$Treatment, ref = "Placebo")
EMGEventData$IRI_EC_z <- scale(EMGEventData$IRI_EC)

# Build model
lme1 <- lme(EMG_corr_mean ~ Stimulus*(Treatment + IRI_EC_z) + Wave, data = EMGEventData, random = ~1|Subject, na.action = na.omit)
lme1b <- lme(EMG_corr_mean ~ Stimulus*Treatment + Wave, data = EMGEventData, random = ~1|Subject, na.action = na.omit)
lme1c <- lme(EMG_corr_mean ~ Stimulus*(Treatment + IRI_EC_z) + Wave, data = subset(EMGEventData, Stimulus %in% c("Angry", "Happy")), random = ~1|Subject, na.action = na.omit)

plot(lme1)
summary(lme1)
intervals(lme1)

eff1 <- effect("Stimulus*Treatment", lme1)

# Compare plots to less custom-generated output for verification
plot(effect("Treatment*Stimulus", lme1))

pdf("Fig_EMG3.pdf", width = 7, height = 5)
plot(c(eff1$fit[1], eff1$fit[2], eff1$fit[3]),
     type = "b",
     frame.plot = F,
     ylab = "EMG (log ratio)",
     xlab = "Stimulus type",
     xaxt = "n",
     yaxt = "n",
     xlim = c(1, 3.1),
     ylim = c(-0.2, 0.1),
     col = col1,
     main = "Corrugator EMG",
)
lines(c(1, 1), c(eff1$upper[1], eff1$lower[1]), col = col1)
lines(c(2, 2), c(eff1$upper[2], eff1$lower[2]), col = col1)
lines(c(3, 3), c(eff1$upper[3], eff1$lower[3]), col = col1)
lines(c(1.1, 2.1, 3.1), c(eff1$fit[4], eff1$fit[5], eff1$fit[6]), type = "b", col = col2, pch = 16)
lines(c(1.1, 1.1), c(eff1$upper[4], eff1$lower[4]), col = col2)
lines(c(2.1, 2.1), c(eff1$upper[5], eff1$lower[5]), col = col2)
lines(c(3.1, 3.1), c(eff1$upper[6], eff1$lower[6]), col = col2)
axis(1, at = c(1.05, 2.05, 3.05), labels = c("Neutral", "Angry", "Happy"))
axis(2, at = c(-0.2, -0.1, 0, 0.1))
legend("topright", col = c(col1, col2), pch = c(1, 16), legend = c("Placebo", "Oxazepam"), bty = "n", lty = 1)
dev.off()

# Analyse other rating scales as predictors
# Since rating scales may be collinear, we include them one by one

# IRI-PT
EMGEventData$IRI_PT_z <- scale(EMGEventData$IRI_PT)
lme2 <- lme(EMG_corr_mean ~ Stimulus*(Treatment + IRI_PT_z) + Wave, data = EMGEventData, random = ~1|Subject, na.action = na.omit)
plot(lme2)
summary(lme2)
intervals(lme2)
lme2c <- lme(EMG_corr_mean ~ Stimulus*(Treatment + IRI_PT_z) + Wave, data = subset(EMGEventData, Stimulus %in% c("Angry", "Happy")), random = ~1|Subject, na.action = na.omit)

# IRI-PD
EMGEventData$IRI_PD_z <- scale(EMGEventData$IRI_PD)
lme3 <- lme(EMG_corr_mean ~ Stimulus*(Treatment + IRI_PD_z) + Wave, data = EMGEventData, random = ~1|Subject, na.action = na.omit)
plot(lme3)
summary(lme3)
intervals(lme3)
lme3c <- lme(EMG_corr_mean ~ Stimulus*(Treatment + IRI_PD_z) + Wave, data = subset(EMGEventData, Stimulus %in% c("Angry", "Happy")), random = ~1|Subject, na.action = na.omit)

# IRI-F
EMGEventData$IRI_F_z <- scale(EMGEventData$IRI_F)
lme4 <- lme(EMG_corr_mean ~ Stimulus*(Treatment + IRI_F_z) + Wave, data = EMGEventData, random = ~1|Subject, na.action = na.omit)
plot(lme4)
summary(lme4)
intervals(lme4)
lme4c <- lme(EMG_corr_mean ~ Stimulus*(Treatment + IRI_F_z) + Wave, data = subset(EMGEventData, Stimulus %in% c("Angry", "Happy")), random = ~1|Subject, na.action = na.omit)

# STAI-T
EMGEventData$STAI.T_z <- scale(EMGEventData$STAI.T)
lme5 <- lme(EMG_corr_mean ~ Stimulus*(Treatment + STAI.T_z) + Wave, data = EMGEventData, random = ~1|Subject, na.action = na.omit)
plot(lme5)
summary(lme5)
intervals(lme5)
lme5c <- lme(EMG_corr_mean ~ Stimulus*(Treatment + STAI.T_z) + Wave, data = subset(EMGEventData, Stimulus %in% c("Angry", "Happy")), random = ~1|Subject, na.action = na.omit)

# TAS-20
EMGEventData$TAS.20_z <- scale(EMGEventData$TAS.20)
lme6 <- lme(EMG_corr_mean ~ Stimulus*(Treatment + TAS.20_z) + Wave, data = EMGEventData, random = ~1|Subject, na.action = na.omit)
plot(lme6)
summary(lme6)
intervals(lme6)
lme6c <- lme(EMG_corr_mean ~ Stimulus*(Treatment + TAS.20_z) + Wave, data = subset(EMGEventData, Stimulus %in% c("Angry", "Happy")), random = ~1|Subject, na.action = na.omit)

# PPI-R-SCI
EMGEventData$PPI_SCI_z <- EMGEventData$PPI_1_SCI_R
EMGEventData$PPI_SCI_z[EMGEventData$PPI_1_IR >= 45] <- NA # Exclude participants with too high scores on Inconsistent Responding measure
EMGEventData$PPI_SCI_z <- scale(EMGEventData$PPI_SCI_z)
lme7 <- lme(EMG_corr_mean ~ Stimulus*(Treatment + PPI_SCI_z) + Wave, data = EMGEventData, random = ~1|Subject, na.action = na.omit)
plot(lme7)
summary(lme7)
intervals(lme7)
lme7c <- lme(EMG_corr_mean ~ Stimulus*(Treatment + PPI_SCI_z) + Wave, data = subset(EMGEventData, Stimulus %in% c("Angry", "Happy")), random = ~1|Subject, na.action = na.omit)

# PPI-R-FD
EMGEventData$PPI_FD_z <- EMGEventData$PPI_1_FD_R
EMGEventData$PPI_FD_z[EMGEventData$PPI_1_IR >= 45] <- NA # Exclude participants with too high scores on Inconsistent Responding measure
EMGEventData$PPI_FD_z <- scale(EMGEventData$PPI_FD_z)
lme8 <- lme(EMG_corr_mean ~ Stimulus*(Treatment + PPI_FD_z) + Wave, data = EMGEventData, random = ~1|Subject, na.action = na.omit)
plot(lme8)
summary(lme8)
intervals(lme8)
lme8c <- lme(EMG_corr_mean ~ Stimulus*(Treatment + PPI_FD_z) + Wave, data = subset(EMGEventData, Stimulus %in% c("Angry", "Happy")), random = ~1|Subject, na.action = na.omit)

# PPI-R-C
EMGEventData$PPI_C_z <- EMGEventData$PPI_1_C_R
EMGEventData$PPI_C_z[EMGEventData$PPI_1_IR >= 45] <- NA # Exclude participants with too high scores on Inconsistent Responding measure
EMGEventData$PPI_C_z <- scale(EMGEventData$PPI_C_z)
lme9 <- lme(EMG_corr_mean ~ Stimulus*(Treatment + PPI_C_z) + Wave, data = EMGEventData, random = ~1|Subject, na.action = na.omit)
plot(lme9)
summary(lme9)
intervals(lme9)
lme9c <- lme(EMG_corr_mean ~ Stimulus*(Treatment + PPI_C_z) + Wave, data = subset(EMGEventData, Stimulus %in% c("Angry", "Happy")), random = ~1|Subject, na.action = na.omit)


# Plot effects of rating scales
# Corrugator responses to angry faces
data_corr_angry <- data.frame(scale = "PPI-R-C", beta = intervals(lme9)$fixed[9, 2], lower = intervals(lme9)$fixed[9, 1], upper = intervals(lme9)$fixed[9, 3], group = "PPI", p = round(summary(lme9)$tTable[9, 5], 3))
data_corr_angry <- rbind(data_corr_angry, data.frame(scale = "PPI-R-FD", beta = intervals(lme8)$fixed[9, 2], lower = intervals(lme8)$fixed[9, 1], upper = intervals(lme8)$fixed[9, 3], group = "PPI", p = round(summary(lme8)$tTable[9, 5], 3)))
data_corr_angry <- rbind(data_corr_angry, data.frame(scale = "PPI-R-SCI", beta = intervals(lme7)$fixed[9, 2], lower = intervals(lme7)$fixed[9, 1], upper = intervals(lme7)$fixed[9, 3], group = "PPI", p = round(summary(lme7)$tTable[9, 5], 3)))
data_corr_angry <- rbind(data_corr_angry, data.frame(scale = "TAS-20", beta = intervals(lme6)$fixed[9, 2], lower = intervals(lme6)$fixed[9, 1], upper = intervals(lme6)$fixed[9, 3], group = "TAS", p = round(summary(lme6)$tTable[9, 5], 3)))
data_corr_angry <- rbind(data_corr_angry, data.frame(scale = "STAI-T", beta = intervals(lme5)$fixed[9, 2], lower = intervals(lme5)$fixed[9, 1], upper = intervals(lme5)$fixed[9, 3], group = "STAI", p = round(summary(lme5)$tTable[9, 5], 3)))
data_corr_angry <- rbind(data_corr_angry, data.frame(scale = "IRI-F", beta = intervals(lme4)$fixed[9, 2], lower = intervals(lme4)$fixed[9, 1], upper = intervals(lme4)$fixed[9, 3], group = "IRI", p = round(summary(lme4)$tTable[9, 5], 3)))
data_corr_angry <- rbind(data_corr_angry, data.frame(scale = "IRI-PD", beta = intervals(lme3)$fixed[9, 2], lower = intervals(lme3)$fixed[9, 1], upper = intervals(lme3)$fixed[9, 3], group = "IRI", p = round(summary(lme3)$tTable[9, 5], 3)))
data_corr_angry <- rbind(data_corr_angry, data.frame(scale = "IRI-PT", beta = intervals(lme2)$fixed[9, 2], lower = intervals(lme2)$fixed[9, 1], upper = intervals(lme2)$fixed[9, 3], group = "IRI", p = round(summary(lme2)$tTable[9, 5], 3)))
data_corr_angry <- rbind(data_corr_angry, data.frame(scale = "IRI-EC", beta = intervals(lme1)$fixed[9, 2], lower = intervals(lme1)$fixed[9, 1], upper = intervals(lme1)$fixed[9, 3], group = "IRI", p = round(summary(lme1)$tTable[9, 5], 3)))
data_corr_angry <- data_corr_angry[c(9:1), ] # Reverse order

pdf("Fig_EMG4.pdf", width = 7, height = 5)
par(mar=c(5.1, 5.4, 4, 14))
plot(x = data_corr_angry$beta, y = c(12:9, 7, 5, 3:1), xlab = expression(beta), ylab = "", frame.plot = F, xlim = c(min(data_corr_angry$lower), max(data_corr_angry$upper)), xaxt = "n", yaxt = "n")
title("Corrugator, Angry", line = 2)
abline(v = 0, col = "gray")
axis(1, at = c(-0.04, 0, 0.04), labels = c(-0.04, 0, 0.04))
points(x = data_corr_angry$beta, y = c(12:9, 7, 5, 3:1), pch = 16)
arrows(data_corr_angry$lower, c(12:9, 7, 5, 3:1), data_corr_angry$upper, c(12:9, 7, 5, 3:1), length = 0, lwd = 1.5)
par(las=1)
mtext(side = 2, at = c(12:9, 7, 5, 3:1), text = c("IRI-EC", " IRI-PT", "IRI-PD", "IRI-F", "STAI-T", "TAS-20", "PPI-R-SCI", "PPI-R_FD", "PPI-R-C"), line = 1)
mtext(side = 4, at = 13, text = expression(beta), line = 1)
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = round(data_corr_angry$beta, 2), line = 1)
mtext(side = 4, at = 13, text = "95 % CI", line = 4)
CI <- paste("[", round(data_corr_angry$lower, 2), ", ", round(data_corr_angry$upper, 2), "]", sep = "")
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = CI, line = 4)
mtext(side = 4, at = 13, text = expression(italic(p)), line = 10)
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = data_corr_angry$p, line = 10)
dev.off()

# Corrugator responses to happy faces
data_corr_happy <- data.frame(scale = "PPI-R-C", beta = intervals(lme9)$fixed[10, 2], lower = intervals(lme9)$fixed[10, 1], upper = intervals(lme9)$fixed[10, 3], group = "PPI", p = round(summary(lme9)$tTable[10, 5], 3))
data_corr_happy <- rbind(data_corr_happy, data.frame(scale = "PPI-R-FD", beta = intervals(lme8)$fixed[10, 2], lower = intervals(lme8)$fixed[10, 1], upper = intervals(lme8)$fixed[10, 3], group = "PPI", p = round(summary(lme8)$tTable[10, 5], 3)))
data_corr_happy <- rbind(data_corr_happy, data.frame(scale = "PPI-R-SCI", beta = intervals(lme7)$fixed[10, 2], lower = intervals(lme7)$fixed[10, 1], upper = intervals(lme7)$fixed[10, 3], group = "PPI", p = round(summary(lme7)$tTable[10, 5], 3)))
data_corr_happy <- rbind(data_corr_happy, data.frame(scale = "TAS-20", beta = intervals(lme6)$fixed[10, 2], lower = intervals(lme6)$fixed[10, 1], upper = intervals(lme6)$fixed[10, 3], group = "TAS", p = round(summary(lme6)$tTable[10, 5], 3)))
data_corr_happy <- rbind(data_corr_happy, data.frame(scale = "STAI-T", beta = intervals(lme5)$fixed[10, 2], lower = intervals(lme5)$fixed[10, 1], upper = intervals(lme5)$fixed[10, 3], group = "STAI", p = round(summary(lme5)$tTable[10, 5], 3)))
data_corr_happy <- rbind(data_corr_happy, data.frame(scale = "IRI-F", beta = intervals(lme4)$fixed[10, 2], lower = intervals(lme4)$fixed[10, 1], upper = intervals(lme4)$fixed[10, 3], group = "IRI", p = round(summary(lme4)$tTable[10, 5], 3)))
data_corr_happy <- rbind(data_corr_happy, data.frame(scale = "IRI-PD", beta = intervals(lme3)$fixed[10, 2], lower = intervals(lme3)$fixed[10, 1], upper = intervals(lme3)$fixed[10, 3], group = "IRI", p = round(summary(lme3)$tTable[10, 5], 3)))
data_corr_happy <- rbind(data_corr_happy, data.frame(scale = "IRI-PT", beta = intervals(lme2)$fixed[10, 2], lower = intervals(lme2)$fixed[10, 1], upper = intervals(lme2)$fixed[10, 3], group = "IRI", p = round(summary(lme2)$tTable[10, 5], 3)))
data_corr_happy <- rbind(data_corr_happy, data.frame(scale = "IRI-EC", beta = intervals(lme1)$fixed[10, 2], lower = intervals(lme1)$fixed[10, 1], upper = intervals(lme1)$fixed[10, 3], group = "IRI", p = round(summary(lme1)$tTable[10, 5], 3)))
data_corr_happy <- data_corr_happy[c(9:1), ] # Reverse order

pdf("Fig_EMG5.pdf", width = 7, height = 5)
par(mar=c(5.1, 5.4, 4, 14))
plot(x = data_corr_happy$beta, y = c(12:9, 7, 5, 3:1), xlab = expression(beta), ylab = "", frame.plot = F, xlim = c(min(data_corr_happy$lower), max(data_corr_happy$upper)), xaxt = "n", yaxt = "n")
title("Corrugator, Happy", line = 2)
abline(v = 0, col = "gray")
axis(1, at = c(-0.05, 0, 0.05), labels = c(-0.05, 0, 0.05))
points(x = data_corr_happy$beta, y = c(12:9, 7, 5, 3:1), pch = 16)
arrows(data_corr_happy$lower, c(12:9, 7, 5, 3:1), data_corr_happy$upper, c(12:9, 7, 5, 3:1), length = 0, lwd = 1.5)
par(las=1)
mtext(side = 2, at = c(12:9, 7, 5, 3:1), text = c("IRI-EC", " IRI-PT", "IRI-PD", "IRI-F", "STAI-T", "TAS-20", "PPI-R-SCI", "PPI-R-FD", "PPI-R-C"), line = 1)
mtext(side = 4, at = 13, text = expression(beta), line = 1)
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = round(data_corr_happy$beta, 2), line = 1)
mtext(side = 4, at = 13, text = "95 % CI", line = 4)
CI <- paste("[", round(data_corr_happy$lower, 2), ", ", round(data_corr_happy$upper, 2), "]", sep = "")
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = CI, line = 4)
mtext(side = 4, at = 13, text = expression(italic(p)), line = 10)
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = data_corr_happy$p, line = 10)
dev.off()

# Corrugator responses to angry faces vs happy faces
data_corr_ah <- data.frame(scale = "PPI-R-C", beta = intervals(lme9c)$fixed[7, 2], lower = intervals(lme9c)$fixed[7, 1], upper = intervals(lme9c)$fixed[7, 3], group = "PPI", p = round(summary(lme9c)$tTable[7, 5], 3))
data_corr_ah <- rbind(data_corr_ah, data.frame(scale = "PPI-R-FD", beta = intervals(lme8c)$fixed[7, 2], lower = intervals(lme8c)$fixed[7, 1], upper = intervals(lme8c)$fixed[7, 3], group = "PPI", p = round(summary(lme8c)$tTable[7, 5], 3)))
data_corr_ah <- rbind(data_corr_ah, data.frame(scale = "PPI-R-SCI", beta = intervals(lme7c)$fixed[7, 2], lower = intervals(lme7c)$fixed[7, 1], upper = intervals(lme7c)$fixed[7, 3], group = "PPI", p = round(summary(lme7c)$tTable[7, 5], 3)))
data_corr_ah <- rbind(data_corr_ah, data.frame(scale = "TAS-20", beta = intervals(lme6c)$fixed[7, 2], lower = intervals(lme6c)$fixed[7, 1], upper = intervals(lme6c)$fixed[7, 3], group = "TAS", p = round(summary(lme6c)$tTable[7, 5], 3)))
data_corr_ah <- rbind(data_corr_ah, data.frame(scale = "STAI-T", beta = intervals(lme5c)$fixed[7, 2], lower = intervals(lme5c)$fixed[7, 1], upper = intervals(lme5c)$fixed[7, 3], group = "STAI", p = round(summary(lme5c)$tTable[7, 5], 3)))
data_corr_ah <- rbind(data_corr_ah, data.frame(scale = "IRI-F", beta = intervals(lme4c)$fixed[7, 2], lower = intervals(lme4c)$fixed[7, 1], upper = intervals(lme4c)$fixed[7, 3], group = "IRI", p = round(summary(lme4c)$tTable[7, 5], 3)))
data_corr_ah <- rbind(data_corr_ah, data.frame(scale = "IRI-PD", beta = intervals(lme3c)$fixed[7, 2], lower = intervals(lme3c)$fixed[7, 1], upper = intervals(lme3c)$fixed[7, 3], group = "IRI", p = round(summary(lme3c)$tTable[7, 5], 3)))
data_corr_ah <- rbind(data_corr_ah, data.frame(scale = "IRI-PT", beta = intervals(lme2c)$fixed[7, 2], lower = intervals(lme2c)$fixed[7, 1], upper = intervals(lme2c)$fixed[7, 3], group = "IRI", p = round(summary(lme2c)$tTable[7, 5], 3)))
data_corr_ah <- rbind(data_corr_ah, data.frame(scale = "IRI-EC", beta = intervals(lme1c)$fixed[7, 2], lower = intervals(lme1c)$fixed[7, 1], upper = intervals(lme1c)$fixed[7, 3], group = "IRI", p = round(summary(lme1c)$tTable[7, 5], 3)))
data_corr_ah <- data_corr_ah[c(9:1), ] # Reverse order

pdf("Fig_EMG4c.pdf", width = 7, height = 5)
par(mar=c(5.1, 5.4, 4, 14))
plot(x = data_corr_ah$beta, y = c(12:9, 7, 5, 3:1), xlab = expression(beta), ylab = "", frame.plot = F, xlim = c(min(data_corr_ah$lower), max(data_corr_ah$upper)), xaxt = "n", yaxt = "n")
title("Corrugator, happy vs angry", line = 2)
abline(v = 0, col = "gray")
axis(1, at = c(-0.05, 0, 0.05), labels = c(-0.05, 0, 0.05))
points(x = data_corr_ah$beta, y = c(12:9, 7, 5, 3:1), pch = 16)
arrows(data_corr_ah$lower, c(12:9, 7, 5, 3:1), data_corr_ah$upper, c(12:9, 7, 5, 3:1), length = 0, lwd = 1.5)
par(las=1)
mtext(side = 2, at = c(12:9, 7, 5, 3:1), text = c("IRI-EC", " IRI-PT", "IRI-PD", "IRI-F", "STAI-T", "TAS-20", "PPI-R-SCI", "PPI-R-FD", "PPI-R-C"), line = 1)
mtext(side = 4, at = 13, text = expression(beta), line = 1)
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = round(data_corr_ah$beta, 2), line = 1)
mtext(side = 4, at = 13, text = "95 % CI", line = 4)
CI <- paste("[", round(data_corr_ah$lower, 2), ", ", round(data_corr_ah$upper, 2), "]", sep = "")
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = CI, line = 4)
mtext(side = 4, at = 13, text = expression(italic(p)), line = 10)
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = data_corr_ah$p, line = 10)
dev.off()


# Analyze zygomatic responses
# Move all data to one big frame
EMG_zyg <- data.frame()

for(i in 1:length(EMGDataList)){
  if(!is.null(EMGDataList[i][[1]])){
    Angry_zyg <- melt(EMGDataList[i][[1]][[2]], na.rm = F)
    names(Angry_zyg) <- c("time_0.1s", "event_no", "EMG_zyg")
    Angry_zyg$Stimulus <- "Angry"
    
    Happy_zyg <- melt(EMGDataList[i][[1]][[4]], na.rm = F)
    names(Happy_zyg) <- c("time_0.1s", "event_no", "EMG_zyg")
    Happy_zyg$Stimulus <- "Happy"
    
    Neutral_zyg <- melt(EMGDataList[i][[1]][[6]], na.rm = F)
    names(Neutral_zyg) <- c("time_0.1s", "event_no", "EMG_zyg")
    Neutral_zyg$Stimulus <- "Neutral"
    
    EMG_zyg_temp <- rbind(Angry_zyg, Happy_zyg, Neutral_zyg)
    EMG_zyg_temp$subject <- IncludedSubjects[i]
  }
  EMG_zyg <- rbind(EMG_zyg, EMG_zyg_temp)
}

# Convert time scale and set 0 to event onset
EMG_zyg$time_s <- EMG_zyg$time_0.1s/10
EMG_zyg$time_s <- EMG_zyg$time_s -2

# Make spaghetti plots
plot(EMG_zyg ~ time_s, data = subset(EMG_zyg, Stimulus == "Angry"), frame.plot = F, type = 'n', main = "Angry", ylim = c(0, 0.004))
for(i in unique(EMG_zyg$subject)){
  for( j in unique(EMG_zyg$event_no[EMG_zyg$subject == i])){
    lines(EMG_zyg ~ time_s, data = EMG_zyg[EMG_zyg$subject == i & EMG_zyg$event_no == j & EMG_zyg$Stimulus == "Angry", ], col = col4)
  }
}

plot(EMG_zyg ~ time_s, data = subset(EMG_zyg, Stimulus == "Happy"), frame.plot = F, type = 'n', main = "Happy", ylim = c(0, 0.004))
for(i in unique(EMG_zyg$subject)){
  for( j in unique(EMG_zyg$event_no[EMG_zyg$subject == i])){
    lines(EMG_zyg ~ time_s, data = EMG_zyg[EMG_zyg$subject == i & EMG_zyg$event_no == j & EMG_zyg$Stimulus == "Happy", ], col = col7)
  }
}

# Since data have varying baselines between conditions, we index each response to its baseline from -2 to 0 seconds
EMG_zyg$EMG_zyg_index <- NA
for(i in unique(EMG_zyg$subject)){
  for(j in unique(EMG_zyg$event_no[EMG_zyg$subject == i])){
    for(k in c("Angry", "Happy", "Neutral")){
      EMG_zyg$EMG_zyg_index[EMG_zyg$subject == i & EMG_zyg$event_no == j & EMG_zyg$Stimulus == k] <- EMG_zyg$EMG_zyg[EMG_zyg$subject == i & EMG_zyg$event_no == j & EMG_zyg$Stimulus == k]/mean(EMG_zyg$EMG_zyg[EMG_zyg$subject == i & EMG_zyg$event_no == j & EMG_zyg$Stimulus == k][1:20])
    }
  }
}

# Make spaghetti plots of normalised data
plot(EMG_zyg_index ~ time_s, data = subset(EMG_zyg, Stimulus == "Angry"), frame.plot = F, type = 'n', main = "Angry", ylim = c(0, 15))
for(i in unique(EMG_zyg$subject)){
  for(j in unique(EMG_zyg$event_no[EMG_zyg$subject == i])){
    lines(EMG_zyg_index ~ time_s, data = EMG_zyg[EMG_zyg$subject == i & EMG_zyg$event_no == j & EMG_zyg$Stimulus == "Angry", ], col = col4)
  }
}

plot(EMG_zyg_index ~ time_s, data = subset(EMG_zyg, Stimulus == "Happy"), frame.plot = F, type = 'n', main = "Happy", ylim = c(0, 15))
for(i in unique(EMG_zyg$subject)){
  for(j in unique(EMG_zyg$event_no[EMG_zyg$subject == i])){
    lines(EMG_zyg_index ~ time_s, data = EMG_zyg[EMG_zyg$subject == i & EMG_zyg$event_no == j & EMG_zyg$Stimulus == "Happy", ], col = col7)
  }
}

plot(EMG_zyg_index ~ time_s, data = subset(EMG_zyg, Stimulus == "Neutral"), frame.plot = F, type = 'n', main = "Neutral", ylim = c(0, 15))
for(i in unique(EMG_zyg$subject)){
  for(j in unique(EMG_zyg$event_no[EMG_zyg$subject == i])){
    lines(EMG_zyg_index ~ time_s, data = EMG_zyg[EMG_zyg$subject == i & EMG_zyg$event_no == j & EMG_zyg$Stimulus == "Neutral", ], col = col8)
  }
}



MeanEMGzygAngry <- aggregate(EMG_zyg_index ~ time_s, data = subset(EMG_zyg, Stimulus == "Angry"), mean)
MeanEMGzygHappy <- aggregate(EMG_zyg_index ~ time_s, data = subset(EMG_zyg, Stimulus == "Happy"), mean)
MeanEMGzygNeutral <- aggregate(EMG_zyg_index ~ time_s, data = subset(EMG_zyg, Stimulus == "Neutral"), mean)

MeanEMGzygAngryWave1 <- aggregate(EMG_zyg_index ~ time_s, data = subset(EMG_zyg, Stimulus == "Angry" & subject %in% SubjectsWave1), mean)
MeanEMGzygHappyWave1 <- aggregate(EMG_zyg_index ~ time_s, data = subset(EMG_zyg, Stimulus == "Happy" & subject %in% SubjectsWave1), mean)
MeanEMGzygNeutralWave1 <- aggregate(EMG_zyg_index ~ time_s, data = subset(EMG_zyg, Stimulus == "Neutral" & subject %in% SubjectsWave1), mean)

MeanEMGzygAngryWave2 <- aggregate(EMG_zyg_index ~ time_s, data = subset(EMG_zyg, Stimulus == "Angry" & subject %in% SubjectsWave2), mean)
MeanEMGzygHappyWave2 <- aggregate(EMG_zyg_index ~ time_s, data = subset(EMG_zyg, Stimulus == "Happy" & subject %in% SubjectsWave2), mean)
MeanEMGzygNeutralWave2 <- aggregate(EMG_zyg_index ~ time_s, data = subset(EMG_zyg, Stimulus == "Neutral" & subject %in% SubjectsWave2), mean)


plot(MeanEMGzygAngry, type = "l", col = col4, frame.plot = F, main = "zygomatic EMG", ylim = c(0.8, 1.5))
lines(MeanEMGzygHappy, col = col7)
lines(MeanEMGzygNeutral, col = col8)
legend("topleft", lty = 1, col = c(col4, col7, col8), legend = c("Angry", "Happy", "Neutral"), bty = "n")
abline(v = c(0, 2), lty = 2)

plot(MeanEMGzygAngry, type = "n", col = col4, frame.plot = F, main = "zygomatic EMG", ylim = c(0.8, 1.5))
lines(lowess(MeanEMGzygAngry, f = 0.1), col = col4)
lines(lowess(MeanEMGzygHappy, f = 0.1), col = col7)
lines(lowess(MeanEMGzygNeutral, f = 0.1), col = col8)
legend("topleft", lty = 1, col = c(col4, col7, col8), legend = c("Angry", "Happy", "Neutral"), bty = "n")
abline(v = c(0, 2), lty = 2)


# Make figures
pdf("Fig_EMG6.pdf", width = 7, height = 5)
plot(MeanEMGzygAngryWave1, type = "n", col = col4, frame.plot = F, main = "Zygomatic EMG, Wave 1", xlab = "Time (s)", ylab = "EMG (ratio)", xaxt = "n", yaxt = "n", ylim = c(0.8, 1.5))
rect(2, 0.000001, 4, 10, col = "gray88", border = NA)
abline(v = c(0, 2, 5), lty = 3)
axis(1)
axis(2, at = c(0.8, 1, 1.2, 1.4))
lines(lowess(MeanEMGzygAngryWave1, f = 0.1), col = col4, lwd = 2)
lines(lowess(MeanEMGzygHappyWave1, f = 0.1), col = col7, lwd = 2)
lines(lowess(MeanEMGzygNeutralWave1, f = 0.1), col = col8, lwd = 2)
legend("topleft", lty = 1, col = c(col4, col7, col8), legend = c("Angry", "Happy", "Neutral"), bty = "n", lwd = 2)
dev.off()

pdf("Fig_EMG7.pdf", width = 7, height = 5)
plot(MeanEMGzygAngryWave2, type = "n", col = col4, frame.plot = F, main = "Zygomatic EMG, Wave 2", xlab = "Time (s)", ylab = "EMG (ratio)", xaxt = "n", yaxt = "n", ylim = c(0.8, 1.5))
rect(2, 0.000001, 4, 10, col = "gray88", border = NA)
abline(v = c(0, 2, 4), lty = 3)
axis(1)
axis(2, at = c(0.8, 1, 1.2, 1.4))
lines(lowess(MeanEMGzygAngryWave2, f = 0.1), col = col4, lwd = 2)
lines(lowess(MeanEMGzygHappyWave2, f = 0.1), col = col7, lwd = 2)
lines(lowess(MeanEMGzygNeutralWave2, f = 0.1), col = col8, lwd = 2)
dev.off()


# Average data over interval for each event for statistical modelling
EMGEventData2 <- data.frame()
for(i in unique(EMG_zyg$subject)){
  for(j in unique(EMG_zyg$event_no[EMG_zyg$subject == i])){
    for(k in c("Angry", "Happy", "Neutral")){
      if(i %in% demData$Subject[demData$Wave == 1]){
        EMG_zyg_mean <- mean(EMG_zyg$EMG_zyg_index[EMG_zyg$subject == i & EMG_zyg$event_no == j & EMG_zyg$Stimulus == k][40:60])
        EMGEventData2 <- rbind(EMGEventData2, data.frame(i, j, k, EMG_zyg_mean))
      } else if(i %in% demData$Subject[demData$Wave == 2]){
        EMG_zyg_mean <- mean(EMG_zyg$EMG_zyg_index[EMG_zyg$subject == i & EMG_zyg$event_no == j & EMG_zyg$Stimulus == k][40:60])
        EMGEventData2 <- rbind(EMGEventData2, data.frame(i, j, k, EMG_zyg_mean))
      }
    }
  }
}
names(EMGEventData2) <- c("Subject", "event_no", "Stimulus", "EMG_zyg_mean")

hist(EMGEventData2$EMG_zyg_mean)
hist(log(EMGEventData2$EMG_zyg_mean))

EMGEventData2$EMG_zyg_mean <- log(EMGEventData2$EMG_zyg_mean)

# Analyse data
EMGEventData2 <- merge(EMGEventData2, demData, by = "Subject")
EMGEventData2$Subject <- as.factor(EMGEventData2$Subject)
EMGEventData2$Stimulus <- relevel(EMGEventData2$Stimulus, ref = "Neutral")
EMGEventData2$Treatment <- relevel(EMGEventData2$Treatment, ref = "Placebo")
EMGEventData2$IRI_EC_z <- scale(EMGEventData2$IRI_EC)


# Build model
lmezyg1 <- lme(EMG_zyg_mean ~ Stimulus*(Treatment + IRI_EC_z) + Wave, data = EMGEventData2, random = ~1|Subject, na.action = na.omit)
lmezyg1c <- lme(EMG_zyg_mean ~ Stimulus*(Treatment + IRI_EC_z) + Wave, data = subset(EMGEventData2, Stimulus %in% c("Angry", "Happy")), random = ~1|Subject, na.action = na.omit)

plot(lmezyg1)
summary(lmezyg1)
intervals(lmezyg1)

effzyg1 <- effect("Treatment*Stimulus", lmezyg1)

# Compare plots to less custom-generated output for verification
plot(effect("Treatment*Stimulus", lmezyg1))
plot(effect("Stimulus*Treatment", lmezyg1b))

# Plot results
pdf("Fig_EMG8.pdf", height = 5, width = 7)
plot(c(effzyg1$fit[1], effzyg1$fit[3], effzyg1$fit[5]),
     type = "b",
     frame.plot = F,
     ylab = "EMG (ratio)",
     xlab = "Stimulus type",
     xaxt = "n",
     yaxt = "n",
     xlim = c(1, 3.1),
     ylim = c(-0.2, 0.15),
     col = col1,
     main = "Zygomatic EMG"
)
lines(c(1, 1), c(effzyg1$upper[1], effzyg1$lower[1]), col = col1)
lines(c(2, 2), c(effzyg1$upper[3], effzyg1$lower[3]), col = col1)
lines(c(3, 3), c(effzyg1$upper[5], effzyg1$lower[5]), col = col1)
lines(c(1.1, 2.1, 3.1), c(effzyg1$fit[2], effzyg1$fit[4], effzyg1$fit[6]), type = "b", col = col2, pch = 16)
lines(c(1.1, 1.1), c(effzyg1$upper[2], effzyg1$lower[2]), col = col2)
lines(c(2.1, 2.1), c(effzyg1$upper[4], effzyg1$lower[4]), col = col2)
lines(c(3.1, 3.1), c(effzyg1$upper[6], effzyg1$lower[6]), col = col2)
axis(1, at = c(1.05, 2.05, 3.05), labels = c("Neutral", "Angry", "Happy"))
axis(2, at = c(-0.2, -0.1, 0, 0.1))
legend("topleft", col = c(col1, col2), pch = c(1, 16), legend = c("Placebo", "Oxazepam"), bty = "n", lty = 1)
dev.off()

# Analyse other rating scales as predictors
# Since rating scales may be collinear, we include them one by one

# IRI-PT
EMGEventData2$IRI_PT_z <- scale(EMGEventData2$IRI_PT)
lmezyg2 <- lme(EMG_zyg_mean ~ Stimulus*(Treatment + IRI_PT_z) + Wave, data = EMGEventData2, random = ~1|Subject, na.action = na.omit)
plot(lmezyg2)
summary(lmezyg2)
intervals(lmezyg2)
lmezyg2c <- lme(EMG_zyg_mean ~ Stimulus*(Treatment + IRI_PT_z) + Wave, data = subset(EMGEventData2, Stimulus %in% c("Angry", "Happy")), random = ~1|Subject, na.action = na.omit)

# IRI-PD
EMGEventData2$IRI_PD_z <- scale(EMGEventData2$IRI_PD)
lmezyg3 <- lme(EMG_zyg_mean ~ Stimulus*(Treatment + IRI_PD_z) + Wave, data = EMGEventData2, random = ~1|Subject, na.action = na.omit)
plot(lmezyg3)
summary(lmezyg3)
intervals(lmezyg3)
lmezyg3c <- lme(EMG_zyg_mean ~ Stimulus*(Treatment + IRI_PD_z) + Wave, data = subset(EMGEventData2, Stimulus %in% c("Angry", "Happy")), random = ~1|Subject, na.action = na.omit)

# IRI-F
EMGEventData2$IRI_F_z <- scale(EMGEventData2$IRI_F)
lmezyg4 <- lme(EMG_zyg_mean ~ Stimulus*(Treatment + IRI_F_z) + Wave, data = EMGEventData2, random = ~1|Subject, na.action = na.omit)
plot(lmezyg4)
summary(lmezyg4)
intervals(lmezyg4)
lmezyg4c <- lme(EMG_zyg_mean ~ Stimulus*(Treatment + IRI_F_z) + Wave, data = subset(EMGEventData2, Stimulus %in% c("Angry", "Happy")), random = ~1|Subject, na.action = na.omit)

# STAI-T
EMGEventData2$STAI.T_z <- scale(EMGEventData2$STAI.T)
lmezyg5 <- lme(EMG_zyg_mean ~ Stimulus*(Treatment + STAI.T_z) + Wave, data = EMGEventData2, random = ~1|Subject, na.action = na.omit)
plot(lmezyg5)
summary(lmezyg5)
intervals(lmezyg5)
lmezyg5c <- lme(EMG_zyg_mean ~ Stimulus*(Treatment + STAI.T_z) + Wave, data = subset(EMGEventData2, Stimulus %in% c("Angry", "Happy")), random = ~1|Subject, na.action = na.omit)

# TAS-20
EMGEventData2$TAS.20_z <- scale(EMGEventData2$TAS.20)
lmezyg6 <- lme(EMG_zyg_mean ~ Stimulus*(Treatment + TAS.20_z) + Wave, data = EMGEventData2, random = ~1|Subject, na.action = na.omit)
plot(lmezyg6)
summary(lmezyg6)
intervals(lmezyg6)
lmezyg6c <- lme(EMG_zyg_mean ~ Stimulus*(Treatment + TAS.20_z) + Wave, data = subset(EMGEventData2, Stimulus %in% c("Angry", "Happy")), random = ~1|Subject, na.action = na.omit)

# PPI-R-SCI
EMGEventData2$PPI_SCI_z <- EMGEventData2$PPI_1_SCI_R
EMGEventData2$PPI_SCI_z[EMGEventData2$PPI_1_IR >= 45] <- NA # Exclude participants with too high scores on Inconsistent Responding measure
EMGEventData2$PPI_SCI_z <- scale(EMGEventData2$PPI_SCI_z)
lmezyg7 <- lme(EMG_zyg_mean ~ Stimulus*(Treatment + PPI_SCI_z) + Wave, data = EMGEventData2, random = ~1|Subject, na.action = na.omit)
plot(lmezyg7)
summary(lmezyg7)
intervals(lmezyg7)
lmezyg7c <- lme(EMG_zyg_mean ~ Stimulus*(Treatment + PPI_SCI_z) + Wave, data = subset(EMGEventData2, Stimulus %in% c("Angry", "Happy")), random = ~1|Subject, na.action = na.omit)

# PPI-R-FD
EMGEventData2$PPI_FD_z <- EMGEventData2$PPI_1_FD_R
EMGEventData2$PPI_FD_z[EMGEventData2$PPI_1_IR >= 45] <- NA # Exclude participants with too high scores on Inconsistent Responding measure
EMGEventData2$PPI_FD_z <- scale(EMGEventData2$PPI_FD_z)
lmezyg8 <- lme(EMG_zyg_mean ~ Stimulus*(Treatment + PPI_FD_z) + Wave, data = EMGEventData2, random = ~1|Subject, na.action = na.omit)
plot(lmezyg8)
summary(lmezyg8)
intervals(lmezyg8)
lmezyg8c <- lme(EMG_zyg_mean ~ Stimulus*(Treatment + PPI_FD_z) + Wave, data = subset(EMGEventData2, Stimulus %in% c("Angry", "Happy")), random = ~1|Subject, na.action = na.omit)

# PPI-R-C
EMGEventData2$PPI_C_z <- EMGEventData2$PPI_1_C_R
EMGEventData2$PPI_C_z[EMGEventData2$PPI_1_IR >= 45] <- NA # Exclude participants with too high scores on Inconsistent Responding measure
EMGEventData2$PPI_C_z <- scale(EMGEventData2$PPI_C_z)
lmezyg9 <- lme(EMG_zyg_mean ~ Stimulus*(Treatment + PPI_C_z) + Wave, data = EMGEventData2, random = ~1|Subject, na.action = na.omit)
plot(lmezyg9)
summary(lmezyg9)
intervals(lmezyg9)
lmezyg9c <- lme(EMG_zyg_mean ~ Stimulus*(Treatment + PPI_C_z) + Wave, data = subset(EMGEventData2, Stimulus %in% c("Angry", "Happy")), random = ~1|Subject, na.action = na.omit)

# Plot effects of rating scales
# Zygomatic responses to angry faces
data_zyg_angry <- data.frame(scale = "PPI-R-C", beta = intervals(lmezyg9)$fixed[9, 2], lower = intervals(lmezyg9)$fixed[9, 1], upper = intervals(lmezyg9)$fixed[9, 3], group = "PPI", p = round(summary(lmezyg9)$tTable[9, 5], 3))
data_zyg_angry <- rbind(data_zyg_angry, data.frame(scale = "PPI-R-FD", beta = intervals(lmezyg8)$fixed[9, 2], lower = intervals(lmezyg8)$fixed[9, 1], upper = intervals(lmezyg8)$fixed[9, 3], group = "PPI", p = round(summary(lmezyg8)$tTable[9, 5], 3)))
data_zyg_angry <- rbind(data_zyg_angry, data.frame(scale = "PPI-R-SCI", beta = intervals(lmezyg7)$fixed[9, 2], lower = intervals(lmezyg7)$fixed[9, 1], upper = intervals(lmezyg7)$fixed[9, 3], group = "PPI", p = round(summary(lmezyg7)$tTable[9, 5], 3)))
data_zyg_angry <- rbind(data_zyg_angry, data.frame(scale = "TAS-20", beta = intervals(lmezyg6)$fixed[9, 2], lower = intervals(lmezyg6)$fixed[9, 1], upper = intervals(lmezyg6)$fixed[9, 3], group = "TAS", p = round(summary(lmezyg6)$tTable[9, 5], 3)))
data_zyg_angry <- rbind(data_zyg_angry, data.frame(scale = "STAI-T", beta = intervals(lmezyg5)$fixed[9, 2], lower = intervals(lmezyg5)$fixed[9, 1], upper = intervals(lmezyg5)$fixed[9, 3], group = "STAI", p = round(summary(lmezyg5)$tTable[9, 5], 3)))
data_zyg_angry <- rbind(data_zyg_angry, data.frame(scale = "IRI-F", beta = intervals(lmezyg4)$fixed[9, 2], lower = intervals(lmezyg4)$fixed[9, 1], upper = intervals(lmezyg4)$fixed[9, 3], group = "IRI", p = round(summary(lmezyg4)$tTable[9, 5], 3)))
data_zyg_angry <- rbind(data_zyg_angry, data.frame(scale = "IRI-PD", beta = intervals(lmezyg3)$fixed[9, 2], lower = intervals(lmezyg3)$fixed[9, 1], upper = intervals(lmezyg3)$fixed[9, 3], group = "IRI", p = round(summary(lmezyg3)$tTable[9, 5], 3)))
data_zyg_angry <- rbind(data_zyg_angry, data.frame(scale = "IRI-PT", beta = intervals(lmezyg2)$fixed[9, 2], lower = intervals(lmezyg2)$fixed[9, 1], upper = intervals(lmezyg2)$fixed[9, 3], group = "IRI", p = round(summary(lmezyg2)$tTable[9, 5], 3)))
data_zyg_angry <- rbind(data_zyg_angry, data.frame(scale = "IRI-EC", beta = intervals(lmezyg1)$fixed[9, 2], lower = intervals(lmezyg1)$fixed[9, 1], upper = intervals(lmezyg1)$fixed[9, 3], group = "IRI", p = round(summary(lmezyg1)$tTable[9, 5], 3)))
data_zyg_angry <- data_zyg_angry[c(9:1), ] # Reverse order

pdf("Fig_EMG9.pdf", width = 7, height = 5)
par(mar=c(5.1, 5.4, 4, 14))
plot(x = data_zyg_angry$beta, y = c(12:9, 7, 5, 3:1), xlab = expression(beta), ylab = "", frame.plot = F, xlim = c(min(data_zyg_angry$lower), max(data_zyg_angry$upper)), xaxt = "n", yaxt = "n")
title("Zygomatic, Angry", line = 2)
abline(v = 0, col = "gray")
axis(1, at = c(-0.05, 0, 0.05), labels = c(-0.05, 0, 0.05))
points(x = data_zyg_angry$beta, y = c(12:9, 7, 5, 3:1), pch = 16)
arrows(data_zyg_angry$lower, c(12:9, 7, 5, 3:1), data_zyg_angry$upper, c(12:9, 7, 5, 3:1), length = 0, lwd = 1.5)
par(las=1)
mtext(side = 2, at = c(12:9, 7, 5, 3:1), text = c("IRI-EC", " IRI-PT", "IRI-PD", "IRI-F", "STAI-T", "TAS-20", "PPI-R-SCI", "PPI-R_FD", "PPI-R-C"), line = 1)
mtext(side = 4, at = 13, text = expression(beta), line = 1)
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = round(data_zyg_angry$beta, 2), line = 1)
mtext(side = 4, at = 13, text = "95 % CI", line = 4)
CI <- paste("[", round(data_zyg_angry$lower, 2), ", ", round(data_zyg_angry$upper, 2), "]", sep = "")
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = CI, line = 4)
mtext(side = 4, at = 13, text = expression(italic(p)), line = 10)
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = data_zyg_angry$p, line = 10)
dev.off()

# Zygomatic responses to happy faces
data_zyg_happy <- data.frame(scale = "PPI-R-C", beta = intervals(lmezyg9)$fixed[10, 2], lower = intervals(lmezyg9)$fixed[10, 1], upper = intervals(lmezyg9)$fixed[10, 3], group = "PPI", p = round(summary(lmezyg9)$tTable[10, 5], 3))
data_zyg_happy <- rbind(data_zyg_happy, data.frame(scale = "PPI-R-FD", beta = intervals(lmezyg8)$fixed[10, 2], lower = intervals(lmezyg8)$fixed[10, 1], upper = intervals(lmezyg8)$fixed[10, 3], group = "PPI", p = round(summary(lmezyg8)$tTable[10, 5], 3)))
data_zyg_happy <- rbind(data_zyg_happy, data.frame(scale = "PPI-R-SCI", beta = intervals(lmezyg7)$fixed[10, 2], lower = intervals(lmezyg7)$fixed[10, 1], upper = intervals(lmezyg7)$fixed[10, 3], group = "PPI", p = round(summary(lmezyg7)$tTable[10, 5], 3)))
data_zyg_happy <- rbind(data_zyg_happy, data.frame(scale = "TAS-20", beta = intervals(lmezyg6)$fixed[10, 2], lower = intervals(lmezyg6)$fixed[10, 1], upper = intervals(lmezyg6)$fixed[10, 3], group = "TAS", p = round(summary(lmezyg6)$tTable[10, 5], 3)))
data_zyg_happy <- rbind(data_zyg_happy, data.frame(scale = "STAI-T", beta = intervals(lmezyg5)$fixed[10, 2], lower = intervals(lmezyg5)$fixed[10, 1], upper = intervals(lmezyg5)$fixed[10, 3], group = "STAI", p = round(summary(lmezyg5)$tTable[10, 5], 3)))
data_zyg_happy <- rbind(data_zyg_happy, data.frame(scale = "IRI-F", beta = intervals(lmezyg4)$fixed[10, 2], lower = intervals(lmezyg4)$fixed[10, 1], upper = intervals(lmezyg4)$fixed[10, 3], group = "IRI", p = round(summary(lmezyg4)$tTable[10, 5], 3)))
data_zyg_happy <- rbind(data_zyg_happy, data.frame(scale = "IRI-PD", beta = intervals(lmezyg3)$fixed[10, 2], lower = intervals(lmezyg3)$fixed[10, 1], upper = intervals(lmezyg3)$fixed[10, 3], group = "IRI", p = round(summary(lmezyg3)$tTable[10, 5], 3)))
data_zyg_happy <- rbind(data_zyg_happy, data.frame(scale = "IRI-PT", beta = intervals(lmezyg2)$fixed[10, 2], lower = intervals(lmezyg2)$fixed[10, 1], upper = intervals(lmezyg2)$fixed[10, 3], group = "IRI", p = round(summary(lmezyg2)$tTable[10, 5], 3)))
data_zyg_happy <- rbind(data_zyg_happy, data.frame(scale = "IRI-EC", beta = intervals(lmezyg1)$fixed[10, 2], lower = intervals(lmezyg1)$fixed[10, 1], upper = intervals(lmezyg1)$fixed[10, 3], group = "IRI", p = round(summary(lmezyg1)$tTable[10, 5], 3)))
data_zyg_happy <- data_zyg_happy[c(9:1), ] # Reverse order

pdf("Fig_EMG10.pdf", width = 7, height = 5)
par(mar=c(5.1, 5.4, 4, 14))
plot(x = data_zyg_happy$beta, y = c(12:9, 7, 5, 3:1), xlab = expression(beta), ylab = "", frame.plot = F, xlim = c(min(data_zyg_happy$lower), max(data_zyg_happy$upper)), xaxt = "n", yaxt = "n")
title("Zygomatic, happy", line = 2)
abline(v = 0, col = "gray")
axis(1, at = c(-0.05, 0, 0.05), labels = c(-0.05, 0, 0.05))
points(x = data_zyg_happy$beta, y = c(12:9, 7, 5, 3:1), pch = 16)
arrows(data_zyg_happy$lower, c(12:9, 7, 5, 3:1), data_zyg_happy$upper, c(12:9, 7, 5, 3:1), length = 0, lwd = 1.5)
par(las=1)
mtext(side = 2, at = c(12:9, 7, 5, 3:1), text = c("IRI-EC", " IRI-PT", "IRI-PD", "IRI-F", "STAI-T", "TAS-20", "PPI-R-SCI", "PPI-R_FD", "PPI-R-C"), line = 1)
mtext(side = 4, at = 13, text = expression(beta), line = 1)
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = round(data_zyg_happy$beta, 2), line = 1)
mtext(side = 4, at = 13, text = "95 % CI", line = 4)
CI <- paste("[", round(data_zyg_happy$lower, 2), ", ", round(data_zyg_happy$upper, 2), "]", sep = "")
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = CI, line = 4)
mtext(side = 4, at = 13, text = expression(italic(p)), line = 10)
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = data_zyg_happy$p, line = 10)
dev.off()

# Zygomatic responses to angry faces
data_zyg_ah <- data.frame(scale = "PPI-R-C", beta = intervals(lmezyg9c)$fixed[7, 2], lower = intervals(lmezyg9c)$fixed[7, 1], upper = intervals(lmezyg9c)$fixed[7, 3], group = "PPI", p = round(summary(lmezyg9c)$tTable[7, 5], 3))
data_zyg_ah <- rbind(data_zyg_ah, data.frame(scale = "PPI-R-FD", beta = intervals(lmezyg8c)$fixed[7, 2], lower = intervals(lmezyg8c)$fixed[7, 1], upper = intervals(lmezyg8c)$fixed[7, 3], group = "PPI", p = round(summary(lmezyg8c)$tTable[7, 5], 3)))
data_zyg_ah <- rbind(data_zyg_ah, data.frame(scale = "PPI-R-SCI", beta = intervals(lmezyg7c)$fixed[7, 2], lower = intervals(lmezyg7c)$fixed[7, 1], upper = intervals(lmezyg7c)$fixed[7, 3], group = "PPI", p = round(summary(lmezyg7c)$tTable[7, 5], 3)))
data_zyg_ah <- rbind(data_zyg_ah, data.frame(scale = "TAS-20", beta = intervals(lmezyg6c)$fixed[7, 2], lower = intervals(lmezyg6c)$fixed[7, 1], upper = intervals(lmezyg6c)$fixed[7, 3], group = "TAS", p = round(summary(lmezyg6c)$tTable[7, 5], 3)))
data_zyg_ah <- rbind(data_zyg_ah, data.frame(scale = "STAI-T", beta = intervals(lmezyg5c)$fixed[7, 2], lower = intervals(lmezyg5c)$fixed[7, 1], upper = intervals(lmezyg5c)$fixed[7, 3], group = "STAI", p = round(summary(lmezyg5c)$tTable[7, 5], 3)))
data_zyg_ah <- rbind(data_zyg_ah, data.frame(scale = "IRI-F", beta = intervals(lmezyg4c)$fixed[7, 2], lower = intervals(lmezyg4c)$fixed[7, 1], upper = intervals(lmezyg4c)$fixed[7, 3], group = "IRI", p = round(summary(lmezyg4c)$tTable[7, 5], 3)))
data_zyg_ah <- rbind(data_zyg_ah, data.frame(scale = "IRI-PD", beta = intervals(lmezyg3c)$fixed[7, 2], lower = intervals(lmezyg3c)$fixed[7, 1], upper = intervals(lmezyg3c)$fixed[7, 3], group = "IRI", p = round(summary(lmezyg3c)$tTable[7, 5], 3)))
data_zyg_ah <- rbind(data_zyg_ah, data.frame(scale = "IRI-PT", beta = intervals(lmezyg2c)$fixed[7, 2], lower = intervals(lmezyg2c)$fixed[7, 1], upper = intervals(lmezyg2c)$fixed[7, 3], group = "IRI", p = round(summary(lmezyg2c)$tTable[7, 5], 3)))
data_zyg_ah <- rbind(data_zyg_ah, data.frame(scale = "IRI-EC", beta = intervals(lmezyg1c)$fixed[7, 2], lower = intervals(lmezyg1c)$fixed[7, 1], upper = intervals(lmezyg1c)$fixed[7, 3], group = "IRI", p = round(summary(lmezyg1c)$tTable[7, 5], 3)))
data_zyg_ah <- data_zyg_ah[c(9:1), ] # Reverse order

pdf("Fig_EMG9c.pdf", width = 7, height = 5)
par(mar=c(5.1, 5.4, 4, 14))
plot(x = data_zyg_ah$beta, y = c(12:9, 7, 5, 3:1), xlab = expression(beta), ylab = "", frame.plot = F, xlim = c(min(data_zyg_ah$lower), max(data_zyg_ah$upper)), xaxt = "n", yaxt = "n")
title("Zygomatic, happy vs angry", line = 2)
abline(v = 0, col = "gray")
axis(1, at = c(-0.05, 0, 0.05), labels = c(-0.05, 0, 0.05))
points(x = data_zyg_ah$beta, y = c(12:9, 7, 5, 3:1), pch = 16)
arrows(data_zyg_ah$lower, c(12:9, 7, 5, 3:1), data_zyg_ah$upper, c(12:9, 7, 5, 3:1), length = 0, lwd = 1.5)
par(las=1)
mtext(side = 2, at = c(12:9, 7, 5, 3:1), text = c("IRI-EC", " IRI-PT", "IRI-PD", "IRI-F", "STAI-T", "TAS-20", "PPI-R-SCI", "PPI-R-FD", "PPI-R-C"), line = 1)
mtext(side = 4, at = 13, text = expression(beta), line = 1)
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = round(data_zyg_ah$beta, 2), line = 1)
mtext(side = 4, at = 13, text = "95 % CI", line = 4)
CI <- paste("[", round(data_zyg_ah$lower, 2), ", ", round(data_zyg_ah$upper, 2), "]", sep = "")
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = CI, line = 4)
mtext(side = 4, at = 13, text = expression(italic(p)), line = 10)
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = data_zyg_ah$p, line = 10)
dev.off()

# Write regression output tables
write.csv(summary(lme1b)$tTable, file = "Corr_unadjusted.csv")
write.csv(summary(lme1)$tTable, file = "Corr_IRI_EC.csv")

# Make new figures for final publication
pdf("Fig_EMG1a.pdf", height = 8, width = 4*1.61)
par(mfrow=c(2,1))
par(mar=c(4.5, 4.5, 0, 0))
plot(MeanEMGCorrAngryWave1, type = "n", col = col4, frame.plot = F, main = "", xlab = "", ylab = "EMG (ratio)", xaxt = "n", yaxt = "n", ylim = c(0.85, 1.2))
rect(2, 0.000001, 4, 10, col = "gray88", border = NA)
abline(v = c(0, 2, 6), lty = 3, lwd = 2)
axis(1)
axis(2, at = c(0.9, 1, 1.1, 1.2))
lines(lowess(MeanEMGCorrAngryWave1, f = 0.1), col = col4, lwd = 2, lty = 6)
lines(lowess(MeanEMGCorrHappyWave1, f = 0.1), col = col7, lwd = 2, lty = 1)
lines(lowess(MeanEMGCorrNeutralWave1, f= 0.1), col = col8, lwd = 2, lty = 2)
legend("topleft", col = c(col4, col7, col8), legend = c("Angry", "Happy", "Neutral"), lty = c(6, 1, 2), lwd = 2, bty = "n")

plot(MeanEMGCorrAngryWave2, type = "n", col = col4, frame.plot = F, main = "", xlab = "Time (s)", ylab = "EMG (ratio)", xaxt = "n", yaxt = "n", ylim = c(0.85, 1.2))
rect(2, 0.000001, 4, 10, col = "gray88", border = NA)
abline(v = c(0, 2, 4), lty = 3, lwd = 2)
axis(1)
axis(2, at = c(0.9, 1, 1.1, 1.2))
lines(lowess(MeanEMGCorrAngryWave2, f = 0.1), col = col4, lwd = 2, lty = 6)
lines(lowess(MeanEMGCorrHappyWave2, f = 0.1), col = col7, lwd = 2, lty = 1)
lines(lowess(MeanEMGCorrNeutralWave2, f = 0.1), col = col8, lwd = 2, lty = 2)
dev.off()

pdf("Fig_EMG2a.pdf", height = 4, width = 4*1.61)
par(mar=c(4.5, 4.5, 0, 0))
plot(c(eff1$fit[1], eff1$fit[2], eff1$fit[3]),
     type = "b",
     frame.plot = F,
     ylab = "EMG (log ratio)",
     xlab = "Stimulus type",
     xaxt = "n",
     yaxt = "n",
     xlim = c(1, 3.1),
     ylim = c(-0.2, 0.1),
     col = col1,
     main = "",
)
lines(c(1, 1), c(eff1$upper[1], eff1$lower[1]), col = col1)
lines(c(2, 2), c(eff1$upper[2], eff1$lower[2]), col = col1)
lines(c(3, 3), c(eff1$upper[3], eff1$lower[3]), col = col1)
lines(c(1.1, 2.1, 3.1), c(eff1$fit[4], eff1$fit[5], eff1$fit[6]), type = "b", col = col2, pch = 16)
lines(c(1.1, 1.1), c(eff1$upper[4], eff1$lower[4]), col = col2)
lines(c(2.1, 2.1), c(eff1$upper[5], eff1$lower[5]), col = col2)
lines(c(3.1, 3.1), c(eff1$upper[6], eff1$lower[6]), col = col2)
axis(1, at = c(1.05, 2.05, 3.05), labels = c("Neutral", "Angry", "Happy"))
axis(2, at = c(-0.2, -0.1, 0, 0.1))
legend("topright", col = c(col1, col2), pch = c(1, 16), legend = c("Placebo", "Oxazepam"), bty = "n", lty = 1)
dev.off()


pdf("Fig_EMG3a.pdf", height = 8, width = 4*1.61)
par(mfrow=c(2,1))
par(mar=c(4.5, 4.5, 0, 0))
plot(MeanEMGzygAngryWave1, type = "n", col = col4, frame.plot = F, xlab = "", ylab = "EMG (ratio)", xaxt = "n", yaxt = "n", ylim = c(0.9, 1.45))
rect(2, 0.000001, 4, 10, col = "gray88", border = NA)
abline(v = c(0, 2, 6), lty = 3, lwd = 2)
axis(1)
axis(2, at = c(0.9, 1, 1.1, 1.2, 1.3, 1.4))
lines(lowess(MeanEMGzygAngryWave1, f = 0.1), col = col4, lwd = 2, lty = 6)
lines(lowess(MeanEMGzygHappyWave1, f = 0.1), col = col7, lwd = 2, lty = 1)
lines(lowess(MeanEMGzygNeutralWave1, f = 0.1), col = col8, lwd = 2, lty = 2)
legend("topleft", lty = c(6, 1, 2), col = c(col4, col7, col8), legend = c("Angry", "Happy", "Neutral"), bty = "n", lwd = 2)

plot(MeanEMGzygAngryWave2, type = "n", col = col4, frame.plot = F, xlab = "Time (s)", ylab = "EMG (ratio)", xaxt = "n", yaxt = "n", ylim = c(0.9, 1.45))
rect(2, 0.000001, 4, 10, col = "gray88", border = NA)
abline(v = c(0, 2, 4), lty = 3, lwd = 2)
axis(1)
axis(2, at = c(0.9, 1, 1.1, 1.2, 1.3, 1.4))
lines(lowess(MeanEMGzygAngryWave2, f = 0.1), col = col4, lwd = 2, lty = 6)
lines(lowess(MeanEMGzygHappyWave2, f = 0.1), col = col7, lwd = 2, lty = 1)
lines(lowess(MeanEMGzygNeutralWave2, f = 0.1), col = col8, lwd = 2, lty = 2)
dev.off()

pdf("Fig_EMG4a.pdf", height = 5, width = 7)
par(mar=c(4.5, 4.5, 0, 0))
plot(c(effzyg1$fit[1], effzyg1$fit[3], effzyg1$fit[5]),
     type = "b",
     frame.plot = F,
     ylab = "EMG (ratio)",
     xlab = "Stimulus type",
     xaxt = "n",
     yaxt = "n",
     xlim = c(1, 3.1),
     ylim = c(-0.2, 0.15),
     col = col1,
     main = ""
)
lines(c(1, 1), c(effzyg1$upper[1], effzyg1$lower[1]), col = col1)
lines(c(2, 2), c(effzyg1$upper[3], effzyg1$lower[3]), col = col1)
lines(c(3, 3), c(effzyg1$upper[5], effzyg1$lower[5]), col = col1)
lines(c(1.1, 2.1, 3.1), c(effzyg1$fit[2], effzyg1$fit[4], effzyg1$fit[6]), type = "b", col = col2, pch = 16)
lines(c(1.1, 1.1), c(effzyg1$upper[2], effzyg1$lower[2]), col = col2)
lines(c(2.1, 2.1), c(effzyg1$upper[4], effzyg1$lower[4]), col = col2)
lines(c(3.1, 3.1), c(effzyg1$upper[6], effzyg1$lower[6]), col = col2)
axis(1, at = c(1.05, 2.05, 3.05), labels = c("Neutral", "Angry", "Happy"))
axis(2, at = c(-0.2, -0.1, 0, 0.1))
legend("topleft", col = c(col1, col2), pch = c(1, 16), legend = c("Placebo", "Oxazepam"), bty = "n", lty = 1)
dev.off()


# Yet a new version of figures for publication
pdf("Fig_EMG1b.pdf", height = 7.0866*1.4/3, width = 7.0866*1.4)
par(mfrow=c(1,3))
par(mar=c(4.5, 4.5, 1, 1))
plot(MeanEMGCorrAngryWave1, type = "n", col = col4, frame.plot = F, main = "", xlab = "Time (s)", ylab = "EMG (ratio)", xaxt = "n", yaxt = "n", ylim = c(0.85, 1.2))
rect(2, 0.000001, 4, 10, col = "gray88", border = NA)
abline(v = c(0, 2, 5), lty = 3)
legend("topleft", lty = 1, col = c(col4, col7, col8), legend = c("Angry", "Happy", "Neutral"), lwd = 2, box.lty = 0, bg = "white")
axis(1)
axis(2, at = c(0.9, 1, 1.1, 1.2))
lines(lowess(MeanEMGCorrAngryWave1, f = 0.1), col = col4, lwd = 2)
lines(lowess(MeanEMGCorrHappyWave1, f = 0.1), col = col7, lwd = 2)
lines(lowess(MeanEMGCorrNeutralWave1, f= 0.1), col = col8, lwd = 2)

plot(MeanEMGCorrAngryWave2, type = "n", col = col4, frame.plot = F, main = "", xlab = "Time (s)", ylab = "EMG (ratio)", xaxt = "n", yaxt = "n", ylim = c(0.85, 1.2))
rect(2, 0.000001, 4, 10, col = "gray88", border = NA)
abline(v = c(0, 2, 4), lty = 3)
legend("topleft", lty = 1, col = c(col4, col7, col8), legend = c("Angry", "Happy", "Neutral"), lwd = 2, box.lty = 0, bg = "white")
axis(1)
axis(2, at = c(0.9, 1, 1.1, 1.2))
lines(lowess(MeanEMGCorrAngryWave2, f = 0.1), col = col4, lwd = 2)
lines(lowess(MeanEMGCorrHappyWave2, f = 0.1), col = col7, lwd = 2)
lines(lowess(MeanEMGCorrNeutralWave2, f = 0.1), col = col8, lwd = 2)

plot(c(eff1$fit[1], eff1$fit[2], eff1$fit[3]),
     type = "b",
     frame.plot = F,
     ylab = "EMG (log ratio)",
     xlab = "Stimulus type",
     xaxt = "n",
     yaxt = "n",
     xlim = c(1, 3.1),
     ylim = c(-0.2, 0.1),
     col = col1,
     main = "",
)
lines(c(1, 1), c(eff1$upper[1], eff1$lower[1]), col = col1)
lines(c(2, 2), c(eff1$upper[2], eff1$lower[2]), col = col1)
lines(c(3, 3), c(eff1$upper[3], eff1$lower[3]), col = col1)
lines(c(1.1, 2.1, 3.1), c(eff1$fit[4], eff1$fit[5], eff1$fit[6]), type = "b", col = col2, pch = 16)
lines(c(1.1, 1.1), c(eff1$upper[4], eff1$lower[4]), col = col2)
lines(c(2.1, 2.1), c(eff1$upper[5], eff1$lower[5]), col = col2)
lines(c(3.1, 3.1), c(eff1$upper[6], eff1$lower[6]), col = col2)
axis(1, at = c(1.05, 2.05, 3.05), labels = c("Neutral", "Angry", "Happy"))
axis(2, at = c(-0.2, -0.1, 0, 0.1))
legend("topright", col = c(col1, col2), pch = c(1, 16), legend = c("Placebo", "Oxazepam"), bty = "n", lty = 1)
dev.off()