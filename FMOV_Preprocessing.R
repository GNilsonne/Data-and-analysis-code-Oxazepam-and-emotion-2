# Script to preprocess EMG data from the FMOV experiment of the Oxazepam project 2011
# Gustav Nilsonne 2015-06-09

# Require packages
require(RCurl) # To read data from GitHub
require(reshape2) # To reshape data
require(psych) # To convert r to z

# Define function to read data
fun_readEMGdata <- function(X){
  
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
  EMGData <- list(as.integer(X), AngryEMG_zygdownsampled, AngryEMG_corrdownsampled, HappyEMG_zygdownsampled, HappyEMG_corrdownsampled, NeutralEMG_zygdownsampled, NeutralEMG_corrdownsampled, cor.test(data$EMG_zyg, data$EMG_corr)$estimate)
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
EMGDataList <- lapply(IncludedSubjects, FUN = fun_readEMGdata)

# Extract data on correlation between corrugator and zygomatic activity
r_vec <- c() # Initialise vector
for(i in 1:length(EMGDataList)){
  if(!is.null(EMGDataList[i][[1]])){
    r <- EMGDataList[[i]][[8]]
  }
  r_vec <- c(r_vec, r)
}
DataCorrelationZygomaticCorrugator <- data.frame(Subject = IncludedSubjects, r = r_vec)
DataCorrelationZygomaticCorrugator$z <- fisherz(DataCorrelationZygomaticCorrugator$r)
write.csv(DataCorrelationZygomaticCorrugator, file = "CorrelationZygomaticCorrugator.csv", row.names = F)

# Move corrugator data to one big data frame
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

# Since data have varying baselines between conditions, we index each response to its baseline from -2 to 0 seconds
EMG_corr$EMG_corr_index <- NA
for(i in unique(EMG_corr$subject)){
  for(j in unique(EMG_corr$event_no[EMG_corr$subject == i])){
    for(k in c("Angry", "Happy", "Neutral")){
      EMG_corr$EMG_corr_index[EMG_corr$subject == i & EMG_corr$event_no == j & EMG_corr$Stimulus == k] <- EMG_corr$EMG_corr[EMG_corr$subject == i & EMG_corr$event_no == j & EMG_corr$Stimulus == k]/mean(EMG_corr$EMG_corr[EMG_corr$subject == i & EMG_corr$event_no == j & EMG_corr$Stimulus == k][1:20])
    }
  }
}

# Get mean time courses
MeanEMGCorrAngry <- aggregate(EMG_corr_index ~ time_s, data = subset(EMG_corr, Stimulus == "Angry"), mean)
MeanEMGCorrHappy <- aggregate(EMG_corr_index ~ time_s, data = subset(EMG_corr, Stimulus == "Happy"), mean)
MeanEMGCorrNeutral <- aggregate(EMG_corr_index ~ time_s, data = subset(EMG_corr, Stimulus == "Neutral"), mean)

MeanEMGCorrAngryWave1 <- aggregate(EMG_corr_index ~ time_s, data = subset(EMG_corr, Stimulus == "Angry" & subject %in% SubjectsWave1), mean)
MeanEMGCorrHappyWave1 <- aggregate(EMG_corr_index ~ time_s, data = subset(EMG_corr, Stimulus == "Happy" & subject %in% SubjectsWave1), mean)
MeanEMGCorrNeutralWave1 <- aggregate(EMG_corr_index ~ time_s, data = subset(EMG_corr, Stimulus == "Neutral" & subject %in% SubjectsWave1), mean)

MeanEMGCorrAngryWave2 <- aggregate(EMG_corr_index ~ time_s, data = subset(EMG_corr, Stimulus == "Angry" & subject %in% SubjectsWave2), mean)
MeanEMGCorrHappyWave2 <- aggregate(EMG_corr_index ~ time_s, data = subset(EMG_corr, Stimulus == "Happy" & subject %in% SubjectsWave2), mean)
MeanEMGCorrNeutralWave2 <- aggregate(EMG_corr_index ~ time_s, data = subset(EMG_corr, Stimulus == "Neutral" & subject %in% SubjectsWave2), mean)

MeanTimecoursesCorrugator <- cbind(MeanEMGCorrAngryWave1, MeanEMGCorrHappyWave1, MeanEMGCorrNeutralWave1, MeanEMGCorrAngryWave2, MeanEMGCorrHappyWave2, MeanEMGCorrNeutralWave2)
MeanTimecoursesCorrugator <- MeanTimecoursesCorrugator[, -c(3, 5, 7, 9, 11)]
names(MeanTimecoursesCorrugator) <- c("time_s", "AngryWave1", "HappyWave1", "NeutralWave1", "AngryWave2", "HappyWave2", "NeutralWave2")
write.csv(MeanTimecoursesCorrugator, file = "MeanTimecoursesCorrugator.csv", row.names = F)

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

write.csv(EMGEventData, file = "CorrugatorEventData.csv", row.names = F)


# Analyze zygomatic responses
# Move data to one big frame
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

# Since data have varying baselines between conditions, we index each response to its baseline from -2 to 0 seconds
EMG_zyg$EMG_zyg_index <- NA
for(i in unique(EMG_zyg$subject)){
  for(j in unique(EMG_zyg$event_no[EMG_zyg$subject == i])){
    for(k in c("Angry", "Happy", "Neutral")){
      EMG_zyg$EMG_zyg_index[EMG_zyg$subject == i & EMG_zyg$event_no == j & EMG_zyg$Stimulus == k] <- EMG_zyg$EMG_zyg[EMG_zyg$subject == i & EMG_zyg$event_no == j & EMG_zyg$Stimulus == k]/mean(EMG_zyg$EMG_zyg[EMG_zyg$subject == i & EMG_zyg$event_no == j & EMG_zyg$Stimulus == k][1:20])
    }
  }
}

# Get mean time courses
MeanEMGzygAngry <- aggregate(EMG_zyg_index ~ time_s, data = subset(EMG_zyg, Stimulus == "Angry"), mean)
MeanEMGzygHappy <- aggregate(EMG_zyg_index ~ time_s, data = subset(EMG_zyg, Stimulus == "Happy"), mean)
MeanEMGzygNeutral <- aggregate(EMG_zyg_index ~ time_s, data = subset(EMG_zyg, Stimulus == "Neutral"), mean)

MeanEMGzygAngryWave1 <- aggregate(EMG_zyg_index ~ time_s, data = subset(EMG_zyg, Stimulus == "Angry" & subject %in% SubjectsWave1), mean)
MeanEMGzygHappyWave1 <- aggregate(EMG_zyg_index ~ time_s, data = subset(EMG_zyg, Stimulus == "Happy" & subject %in% SubjectsWave1), mean)
MeanEMGzygNeutralWave1 <- aggregate(EMG_zyg_index ~ time_s, data = subset(EMG_zyg, Stimulus == "Neutral" & subject %in% SubjectsWave1), mean)

MeanEMGzygAngryWave2 <- aggregate(EMG_zyg_index ~ time_s, data = subset(EMG_zyg, Stimulus == "Angry" & subject %in% SubjectsWave2), mean)
MeanEMGzygHappyWave2 <- aggregate(EMG_zyg_index ~ time_s, data = subset(EMG_zyg, Stimulus == "Happy" & subject %in% SubjectsWave2), mean)
MeanEMGzygNeutralWave2 <- aggregate(EMG_zyg_index ~ time_s, data = subset(EMG_zyg, Stimulus == "Neutral" & subject %in% SubjectsWave2), mean)

MeanTimecoursesZygomatic <- cbind(MeanEMGzygAngryWave1, MeanEMGzygHappyWave1, MeanEMGzygNeutralWave1, MeanEMGzygAngryWave2, MeanEMGzygHappyWave2, MeanEMGzygNeutralWave2)
MeanTimecoursesZygomatic <- MeanTimecoursesZygomatic[, -c(3, 5, 7, 9, 11)]
names(MeanTimecoursesZygomatic) <- c("time_s", "AngryWave1", "HappyWave1", "NeutralWave1", "AngryWave2", "HappyWave2", "NeutralWave2")
write.csv(MeanTimecoursesZygomatic, file = "MeanTimecoursesZygomatic.csv", row.names = F)

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

write.csv(EMGEventData2, file = "ZygomaticEventData.csv", row.names = F)