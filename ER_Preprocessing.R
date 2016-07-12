# Script to preprocess EMG and heart rate data from the ER experiment of the Oxazepam project 2011
# Gustav Nilsonne 2015-06-18

# Require packages
require(RCurl) # To read data from GitHub
require(reshape2) # To reshape data
require(parallel) # To increase processing speed

# Define function to read data
fun_read_data <- function(X){
  
  # Read data file
  require(RCurl)
  url <- getURL(paste("https://raw.githubusercontent.com/GNilsonne/Data-and-analysis-code-Oxazepam-and-emotion-2/master/ER/Fp", X, "_ER_processed.txt", sep = ""), ssl.verifypeer = FALSE)
  data <- read.csv(text = url, header = FALSE)
  names(data) <- c("SCR", "Raw_EMG_Zyg", "Raw_EMG_corr", "EMG_zyg", "EKG", "EMG_corr", "down_neg", "up_neg", "down_neu", "up_neu", "")
  
  # Find event onsets
  data$down_neg_diff <- c(0, diff(data$down_neg, lag = 1))
  data$down_neg_onsets <- data$down_neg_diff == -5
  data$down_neu_diff <- c(0, diff(data$down_neu, lag = 1))
  data$down_neu_onsets <- data$down_neu_diff == -5
  data$up_neg_diff <- c(0, diff(data$up_neg, lag = 1))
  data$up_neg_onsets <- data$up_neg_diff == -5
  data$up_neu_diff <- c(0, diff(data$up_neu, lag = 1))
  data$up_neu_onsets <- data$up_neu_diff == -5
  
  # Add time vector. Sampling was at 10 millisecond intervals.
  data$millisecond <- 10*(1:length(data$SCR))
  
  # Define time windows for each event  
  down_neg_index <- which(data$down_neg_onsets == TRUE)
  data$down_neg_window <- FALSE
  for(i in down_neg_index){
    if(length(data[, 1]) > (i + 400)){
      data$down_neg_window[(i - 400):(i + 400)] <- TRUE
    } else {
      this_length <- length(data[, 1]) - i
      data$down_neg_window[(i - 400):(i + down_neg_window)] <- TRUE
    }
  }
  down_neu_index <- which(data$down_neu_onsets == TRUE)
  data$down_neu_window <- FALSE
  for(i in down_neu_index){
    if(length(data[, 1]) > (i + 400)){
      data$down_neu_window[(i - 400):(i + 400)] <- TRUE
    } else {
      this_length <- length(data[, 1]) - i
      data$down_neu_window[(i - 400):(i + down_neu_window)] <- TRUE
    }
  }
  up_neg_index <- which(data$up_neg_onsets == TRUE)
  data$up_neg_window <- FALSE
  for(i in up_neg_index){
    if(length(data[, 1]) > (i + 400)){
      data$up_neg_window[(i - 400):(i + 400)] <- TRUE
    } else {
      this_length <- length(data[, 1]) - i
      data$up_neg_window[(i - 400):(i + up_neg_window)] <- TRUE
    }
  }
  up_neu_index <- which(data$up_neu_onsets == TRUE)
  data$up_neu_window <- FALSE
  for(i in up_neu_index){
    if(length(data[, 1]) > (i + 400)){
      data$up_neu_window[(i - 400):(i + 400)] <- TRUE
    } else {
      this_length <- length(data[, 1]) - i
      data$up_neu_window[(i - 400):(i + up_neu_window)] <- TRUE
    }
  }
  
  # Get SCR data
  # Extract data, downsample
  down_neg_SCR <- matrix(data$SCR[data$down_neg_window == TRUE], ncol = length(down_neg_index))
  down_neg_SCRlowess <- apply(down_neg_SCR, 2, lowess, f = 0.01)
  down_neg_SCRlowess <- matrix(unlist(down_neg_SCRlowess), byrow = F, ncol = length(down_neg_SCRlowess))[802:1602, ]
  plot(down_neg_SCR[, 1], type = "l", col = "gray", frame.plot = F, main = X)
  lines(down_neg_SCRlowess[, 1])
  lines(y = down_neg_SCRlowess[, 1][seq(1, NROW(down_neg_SCRlowess[, 1]), by=10)], x = c(1:801)[seq(1, NROW(down_neg_SCRlowess[, 1]), by=10)], col = "red", type = "o")
  down_neg_SCR_downsampled <- down_neg_SCR[seq(1, NROW(down_neg_SCRlowess), by=10), ]
  
  down_neu_SCR <- matrix(data$SCR[data$down_neu_window == TRUE], ncol = length(down_neu_index))
  down_neu_SCRlowess <- apply(down_neu_SCR, 2, lowess, f = 0.01)
  down_neu_SCRlowess <- matrix(unlist(down_neu_SCRlowess), byrow = F, ncol = length(down_neu_SCRlowess))[802:1602, ]
  down_neu_SCR_downsampled <- down_neu_SCR[seq(1, NROW(down_neu_SCRlowess), by=10), ]
  
  up_neg_SCR <- matrix(data$SCR[data$up_neg_window == TRUE], ncol = length(up_neg_index))
  up_neg_SCRlowess <- apply(up_neg_SCR, 2, lowess, f = 0.01)
  up_neg_SCRlowess <- matrix(unlist(up_neg_SCRlowess), byrow = F, ncol = length(up_neg_SCRlowess))[802:1602, ]
  up_neg_SCR_downsampled <- up_neg_SCR[seq(1, NROW(up_neg_SCRlowess), by=10), ]
  
  up_neu_SCR <- matrix(data$SCR[data$up_neu_window == TRUE], ncol = length(up_neu_index))
  up_neu_SCRlowess <- apply(up_neu_SCR, 2, lowess, f = 0.01)
  up_neu_SCRlowess <- matrix(unlist(up_neu_SCRlowess), byrow = F, ncol = length(up_neu_SCRlowess))[802:1602, ]
  up_neu_SCR_downsampled <- up_neu_SCR[seq(1, NROW(up_neu_SCRlowess), by=10), ]
  
  # Get EMG data
  # Extract data, downsample
  down_neg_EMG_corr <- matrix(data$EMG_corr[data$down_neg_window == TRUE], ncol = length(down_neg_index))
  down_neg_EMG_corrlowess <- apply(down_neg_EMG_corr, 2, lowess, f = 0.01)
  down_neg_EMG_corrlowess <- matrix(unlist(down_neg_EMG_corrlowess), byrow = F, ncol = length(down_neg_EMG_corrlowess))[802:1602, ]
  #plot(down_neg_EMG_corr[, 1], type = "l", col = "gray", frame.plot = F, main = X)
  #lines(down_neg_EMG_corrlowess[, 1])
  #lines(y = down_neg_EMG_corrlowess[, 1][seq(1, NROW(down_neg_EMG_corrlowess[, 1]), by=10)], x = c(1:801)[seq(1, NROW(down_neg_EMG_corrlowess[, 1]), by=10)], col = "red", type = "o")
  down_neg_EMG_corr_downsampled <- down_neg_EMG_corr[seq(1, NROW(down_neg_EMG_corrlowess), by=10), ]
  
  down_neu_EMG_corr <- matrix(data$EMG_corr[data$down_neu_window == TRUE], ncol = length(down_neu_index))
  down_neu_EMG_corrlowess <- apply(down_neu_EMG_corr, 2, lowess, f = 0.01)
  down_neu_EMG_corrlowess <- matrix(unlist(down_neu_EMG_corrlowess), byrow = F, ncol = length(down_neu_EMG_corrlowess))[802:1602, ]
  down_neu_EMG_corr_downsampled <- down_neu_EMG_corr[seq(1, NROW(down_neu_EMG_corrlowess), by=10), ]
  
  up_neg_EMG_corr <- matrix(data$EMG_corr[data$up_neg_window == TRUE], ncol = length(up_neg_index))
  up_neg_EMG_corrlowess <- apply(up_neg_EMG_corr, 2, lowess, f = 0.01)
  up_neg_EMG_corrlowess <- matrix(unlist(up_neg_EMG_corrlowess), byrow = F, ncol = length(up_neg_EMG_corrlowess))[802:1602, ]
  up_neg_EMG_corr_downsampled <- up_neg_EMG_corr[seq(1, NROW(up_neg_EMG_corrlowess), by=10), ]
  
  up_neu_EMG_corr <- matrix(data$EMG_corr[data$up_neu_window == TRUE], ncol = length(up_neu_index))
  up_neu_EMG_corrlowess <- apply(up_neu_EMG_corr, 2, lowess, f = 0.01)
  up_neu_EMG_corrlowess <- matrix(unlist(up_neu_EMG_corrlowess), byrow = F, ncol = length(up_neu_EMG_corrlowess))[802:1602, ]
  up_neu_EMG_corr_downsampled <- up_neu_EMG_corr[seq(1, NROW(up_neu_EMG_corrlowess), by=10), ]
  
  # Get HR data
  require(quantmod) # To find peaks for heartbeats
  if(length(data$EKG) > 0){
    HeartBeats <- findPeaks(data$EKG[data$EKG > 0.22], thresh = 0.0001) # This threshold seems to work all right
    HeartBeatTimes <- data$millisecond[data$EKG > 0.22][HeartBeats]
    HeartBeatIntervals <- (c(HeartBeatTimes, 0) - c(0, HeartBeatTimes))/1000 # Convert to seconds
    
    HeartBeatTimes <- HeartBeatTimes[HeartBeatIntervals > 0.25] # Remove heartbeats that were unrealisticaly close to the last one, they are artefacts of the signal or of the peak finding algorithm
    HeartBeatIntervals <- (c(HeartBeatTimes, 0) - c(0, HeartBeatTimes))/1000 # Convert to seconds
    
    HeartRate <- 60/HeartBeatIntervals
    
    plot(EKG ~ millisecond, data = data[2000:6000, c("EKG", "millisecond")], type = "l", main = X, frame.plot = F, col = "gray")
    points(HeartBeatTimes, rep(0.5, length(HeartBeatTimes)), col = "red", pch = "|")
    
    data$HeartRate <- NA
    index <- data$millisecond %in% HeartBeatTimes
    HeartBeatLocations <- which(index == T)
    for(i in 1:(length(HeartBeatLocations)-1)){
      data$HeartRate[HeartBeatLocations[i]:HeartBeatLocations[i+1]] <- HeartRate[i]
    }
    
    data$HeartRate[data$HeartRate < 40] <- NA # Heart rate below 40 is thought to arise from missing beats, set to NA 
    data$HeartRate[data$HeartRate > 200] <- NA
  }
  down_neg_hr <- matrix(data$HeartRate[data$down_neg_window == TRUE], ncol = length(down_neg_index))
  down_neg_hr <- down_neg_hr[seq(1, NROW(down_neg_hr), by=10), ]
  down_neu_hr <- matrix(data$HeartRate[data$down_neu_window == TRUE], ncol = length(down_neu_index))
  down_neu_hr <- down_neu_hr[seq(1, NROW(down_neu_hr), by=10), ]
  up_neg_hr <- matrix(data$HeartRate[data$up_neg_window == TRUE], ncol = length(up_neg_index))
  up_neg_hr <- up_neg_hr[seq(1, NROW(up_neg_hr), by=10), ]
  up_neu_hr <- matrix(data$HeartRate[data$up_neu_window == TRUE], ncol = length(up_neu_index))
  up_neu_hr <- up_neu_hr[seq(1, NROW(up_neu_hr), by=10), ]
  
  # Make a list for each participant and return it as output
  data <- list(as.integer(X), down_neg_EMG_corr_downsampled, down_neu_EMG_corr_downsampled, up_neg_EMG_corr_downsampled, up_neu_EMG_corr_downsampled, down_neg_hr, down_neu_hr, up_neg_hr, up_neu_hr, down_neg_SCR_downsampled, down_neu_SCR_downsampled, up_neg_SCR_downsampled, up_neu_SCR_downsampled)
  return(data)
}

# Read files
# First demographic data etc
demDataURL <- getURL("https://raw.githubusercontent.com/GNilsonne/Data-and-analysis-code-Oxazepam-and-emotion/master/demographics.csv", ssl.verifypeer = FALSE)
demData <- read.csv(text = demDataURL)
# Then the Acqknowledge logfiles containing EMG data
IncludedSubjects <- demData$Subject[demData$Included_ER == T & demData$Wave == 2] 

# Read Acqknowledge datafiles with parallelization to save time
no_cores <- detectCores() - 1 # Determine number of cores
cl <- makeCluster(no_cores) # Initiate cluster
data_list <- parLapply(cl, IncludedSubjects, fun = fun_read_data)
stopCluster(cl) # Stop cluster


# Move SCR data to one big data frame
SCR <- data.frame()

for(i in 1:length(data_list)){
  if(!is.null(data_list[i][[1]])){
    down_neg_SCR <- melt(data_list[i][[1]][[10]], na.rm = F)
    names(down_neg_SCR) <- c("time_0.1s", "event_no", "SCR")
    down_neg_SCR$instruction <- "Downregulate"
    down_neg_SCR$valence <- "Negative"
    
    down_neu_SCR <- melt(data_list[i][[1]][[11]], na.rm = F)
    names(down_neu_SCR) <- c("time_0.1s", "event_no", "SCR")
    down_neu_SCR$instruction <- "Downregulate"
    down_neu_SCR$valence <- "Neutral"
    
    up_neg_SCR <- melt(data_list[i][[1]][[12]], na.rm = F)
    names(up_neg_SCR) <- c("time_0.1s", "event_no", "SCR")
    up_neg_SCR$instruction <- "Upregulate"
    up_neg_SCR$valence <- "Negative"
    
    up_neu_SCR <- melt(data_list[i][[1]][[13]], na.rm = F)
    names(up_neu_SCR) <- c("time_0.1s", "event_no", "SCR")
    up_neu_SCR$instruction <- "Upregulate"
    up_neu_SCR$valence <- "Neutral"
    
    SCR_temp <- rbind(down_neg_SCR, down_neu_SCR, up_neg_SCR, up_neu_SCR)
    SCR_temp$subject <- IncludedSubjects[i]
  }
  SCR <- rbind(SCR, SCR_temp)
}

# Convert time scale and set 0 to event onset
SCR$time_s <- SCR$time_0.1s/10
SCR$time_s <- SCR$time_s -2

# Since data have varying baselines between conditions, we index each response to its baseline from -2 to 0 seconds
SCR$SCR_index <- NA
for(i in unique(SCR$subject)){
  for(j in unique(SCR$event_no[SCR$subject == i])){
    for(k in c("Downregulate", "Upregulate")){
      for(l in c("Neutral", "Negative")){
        SCR$SCR_index[SCR$subject == i & SCR$event_no == j & SCR$instruction == k & SCR$valence == l] <- SCR$SCR[SCR$subject == i & SCR$event_no == j & SCR$instruction == k & SCR$valence == l]/mean(SCR$SCR[SCR$subject == i & SCR$event_no == j & SCR$instruction == k & SCR$valence == l][1:20])
      }
    }
  }
}

# Get mean time courses
mean_down_neg_SCR <- aggregate(SCR_index ~ time_s, data = SCR[SCR$instruction == "Downregulate" & SCR$valence == "Negative", ], mean)
mean_down_neu_SCR <- aggregate(SCR_index ~ time_s, data = SCR[SCR$instruction == "Downregulate" & SCR$valence == "Neutral", ], mean)
mean_up_neg_SCR <- aggregate(SCR_index ~ time_s, data = SCR[SCR$instruction == "Upregulate" & SCR$valence == "Negative", ], mean)
mean_up_neu_SCR <- aggregate(SCR_index ~ time_s, data = SCR[SCR$instruction == "Upregulate" & SCR$valence == "Neutral", ], mean)

mean_timecourses_SCR <- cbind(mean_down_neg_SCR, mean_down_neu_SCR[, "SCR_index"], mean_up_neg_SCR[, "SCR_index"], mean_up_neu_SCR[, "SCR_index"])
names(mean_timecourses_SCR) <- c("time_s", "down_neg", "down_neu", "up_neg", "up_neu")
write.csv(mean_timecourses_SCR, file = "MeanTimecoursesSCR_ER.csv", row.names = F)

# Plot time courses
plot(mean_down_neg_SCR, type = "n", frame.plot = F, ylim = c(0.998, 1.02))
abline(v = c(0, 2, 3, 5), col = "gray")
lines(mean_down_neg_SCR, col = "red")
lines(mean_down_neu_SCR, col = "pink")
lines(mean_up_neg_SCR, col = "blue")
lines(mean_up_neu_SCR, col = "lightblue")

# Average data over interval for each event for statistical modelling
SCR_event_data <- data.frame()
for(i in unique(SCR$subject)){
  for(j in unique(SCR$event_no[SCR$subject == i])){
    for(k in c("Downregulate", "Upregulate")){
      for(l in c("Neutral", "Negative")){
        SCR_mean_instruction <- mean(SCR$SCR_index[SCR$subject == i & SCR$event_no == j & SCR$instruction == k & SCR$valence == l][21:40])
        SCR_mean_stimulus <- mean(SCR$SCR_index[SCR$subject == i & SCR$event_no == j & SCR$instruction == k & SCR$valence == l][51:70])
        SCR_event_data <- rbind(SCR_event_data, cbind(data.frame(i, j, k, l), SCR_mean_instruction, SCR_mean_stimulus))
      }
    }
  }
}
names(SCR_event_data) <- c("Subject", "event_no", "instruction", "valence", "SCR_mean_instruction", "SCR_mean_stimulus")

hist(SCR_event_data$SCR_mean_stimulus)
hist(log(SCR_event_data$SCR_mean_stimulus))

SCR_event_data$log_SCR_mean_instruction <- log(SCR_event_data$SCR_mean_instruction)
SCR_event_data$log_SCR_mean_stimulus <- log(SCR_event_data$SCR_mean_stimulus)

write.csv(SCR_event_data, file = "SCREventData_ER.csv", row.names = F)


# Move corrugator data to one big data frame
EMG_corr <- data.frame()

for(i in 1:length(data_list)){
  if(!is.null(data_list[i][[1]])){
    down_neg_corr <- melt(data_list[i][[1]][[2]], na.rm = F)
    names(down_neg_corr) <- c("time_0.1s", "event_no", "EMG_corr")
    down_neg_corr$instruction <- "Downregulate"
    down_neg_corr$valence <- "Negative"
    
    down_neu_corr <- melt(data_list[i][[1]][[3]], na.rm = F)
    names(down_neu_corr) <- c("time_0.1s", "event_no", "EMG_corr")
    down_neu_corr$instruction <- "Downregulate"
    down_neu_corr$valence <- "Neutral"
    
    up_neg_corr <- melt(data_list[i][[1]][[4]], na.rm = F)
    names(up_neg_corr) <- c("time_0.1s", "event_no", "EMG_corr")
    up_neg_corr$instruction <- "Upregulate"
    up_neg_corr$valence <- "Negative"
    
    up_neu_corr <- melt(data_list[i][[1]][[5]], na.rm = F)
    names(up_neu_corr) <- c("time_0.1s", "event_no", "EMG_corr")
    up_neu_corr$instruction <- "Upregulate"
    up_neu_corr$valence <- "Neutral"
    
    EMG_corr_temp <- rbind(down_neg_corr, down_neu_corr, up_neg_corr, up_neu_corr)
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
    for(k in c("Downregulate", "Upregulate")){
      for(l in c("Neutral", "Negative")){
        EMG_corr$EMG_corr_index[EMG_corr$subject == i & EMG_corr$event_no == j & EMG_corr$instruction == k & EMG_corr$valence == l] <- EMG_corr$EMG_corr[EMG_corr$subject == i & EMG_corr$event_no == j & EMG_corr$instruction == k & EMG_corr$valence == l]/mean(EMG_corr$EMG_corr[EMG_corr$subject == i & EMG_corr$event_no == j & EMG_corr$instruction == k & EMG_corr$valence == l][1:20])
      }
    }
  }
}

# Get mean time courses
mean_down_neg_corr <- aggregate(EMG_corr_index ~ time_s, data = EMG_corr[EMG_corr$instruction == "Downregulate" & EMG_corr$valence == "Negative", ], mean)
mean_down_neu_corr <- aggregate(EMG_corr_index ~ time_s, data = EMG_corr[EMG_corr$instruction == "Downregulate" & EMG_corr$valence == "Neutral", ], mean)
mean_up_neg_corr <- aggregate(EMG_corr_index ~ time_s, data = EMG_corr[EMG_corr$instruction == "Upregulate" & EMG_corr$valence == "Negative", ], mean)
mean_up_neu_corr <- aggregate(EMG_corr_index ~ time_s, data = EMG_corr[EMG_corr$instruction == "Upregulate" & EMG_corr$valence == "Neutral", ], mean)

mean_timecourses_corrugator <- cbind(mean_down_neg_corr, mean_down_neu_corr[, "EMG_corr_index"], mean_up_neg_corr[, "EMG_corr_index"], mean_up_neu_corr[, "EMG_corr_index"])
names(mean_timecourses_corrugator) <- c("time_s", "down_neg", "down_neu", "up_neg", "up_neu")
write.csv(mean_timecourses_corrugator, file = "MeanTimecoursesCorrugator_ER.csv", row.names = F)

# Plot time courses
plot(mean_down_neg_corr, type = "n", frame.plot = F, ylim = c(0.95, 1.50))
abline(v = c(0, 2, 3, 5), col = "gray")
lines(mean_down_neg_corr, col = "red")
lines(mean_down_neu_corr, col = "pink")
lines(mean_up_neg_corr, col = "blue")
lines(mean_up_neu_corr, col = "lightblue")

# Average data over interval for each event for statistical modelling
emg_event_data <- data.frame()
for(i in unique(EMG_corr$subject)){
  for(j in unique(EMG_corr$event_no[EMG_corr$subject == i])){
    for(k in c("Downregulate", "Upregulate")){
      for(l in c("Neutral", "Negative")){
        emg_corr_mean_instruction <- mean(EMG_corr$EMG_corr_index[EMG_corr$subject == i & EMG_corr$event_no == j & EMG_corr$instruction == k & EMG_corr$valence == l][21:40])
        emg_corr_mean_stimulus <- mean(EMG_corr$EMG_corr_index[EMG_corr$subject == i & EMG_corr$event_no == j & EMG_corr$instruction == k & EMG_corr$valence == l][41:60])
        emg_event_data <- rbind(emg_event_data, cbind(data.frame(i, j, k, l), emg_corr_mean_instruction, emg_corr_mean_stimulus))
      }
    }
  }
}
names(emg_event_data) <- c("Subject", "event_no", "instruction", "valence", "emg_corr_mean_instruction", "emg_corr_mean_stimulus")

hist(emg_event_data$emg_corr_mean_stimulus)
hist(log(emg_event_data$emg_corr_mean_stimulus))

emg_event_data$log_emg_corr_mean_instruction <- log(emg_event_data$emg_corr_mean_instruction)
emg_event_data$log_emg_corr_mean_stimulus <- log(emg_event_data$emg_corr_mean_stimulus)

write.csv(emg_event_data, file = "CorrugatorEventData_ER.csv", row.names = F)

# Heart rate data

# Move heart rate data to one big data frame
hr <- data.frame()
for(i in 1:length(data_list)){
  if(!is.null(data_list[i][[1]])){
    down_neg_hr <- melt(data_list[i][[1]][[6]], na.rm = F)
    names(down_neg_hr) <- c("time_0.1s", "event_no", "hr")
    down_neg_hr$instruction <- "Downregulate"
    down_neg_hr$valence <- "Negative"
    
    down_neu_hr <- melt(data_list[i][[1]][[7]], na.rm = F)
    names(down_neu_hr) <- c("time_0.1s", "event_no", "hr")
    down_neu_hr$instruction <- "Downregulate"
    down_neu_hr$valence <- "Neutral"
    
    up_neg_hr <- melt(data_list[i][[1]][[8]], na.rm = F)
    names(up_neg_hr) <- c("time_0.1s", "event_no", "hr")
    up_neg_hr$instruction <- "Upregulate"
    up_neg_hr$valence <- "Negative"
    
    up_neu_hr <- melt(data_list[i][[1]][[9]], na.rm = F)
    names(up_neu_hr) <- c("time_0.1s", "event_no", "hr")
    up_neu_hr$instruction <- "Upregulate"
    up_neu_hr$valence <- "Neutral"
    
    hr_temp <- rbind(down_neg_hr, down_neu_hr, up_neg_hr, up_neu_hr)
    hr_temp$subject <- IncludedSubjects[i]
  }
  hr <- rbind(hr, hr_temp)
}

# Exclude subject 50 due to poor recording
hr <- hr[hr$subject != 50, ]

# Convert time scale and set 0 to event onset
hr$time_s <- hr$time_0.1s/10
hr$time_s <- hr$time_s -2

# Since data have varying baselines between conditions, we index each response to its baseline from -2 to 0 seconds
# Urry 2009 have indexed additively, I add it for comparability even though indexing might be better 
hr$hr_index <- NA
hr$hr_rel <- NA
for(i in unique(hr$subject)){
  for(j in unique(hr$event_no[hr$subject == i])){
    for(k in c("Downregulate", "Upregulate")){
      for(l in c("Neutral", "Negative")){
        hr$hr_index[hr$subject == i & hr$event_no == j & hr$instruction == k & hr$valence == l] <- hr$hr[hr$subject == i & hr$event_no == j & hr$instruction == k & hr$valence == l]/mean(hr$hr[hr$subject == i & hr$event_no == j & hr$instruction == k & hr$valence == l][1:20])
        hr$hr_rel[hr$subject == i & hr$event_no == j & hr$instruction == k & hr$valence == l] <- hr$hr[hr$subject == i & hr$event_no == j & hr$instruction == k & hr$valence == l] - mean(hr$hr[hr$subject == i & hr$event_no == j & hr$instruction == k & hr$valence == l][1:20])
      }
    }
  }
}

# Get mean time courses
mean_down_neg_hr <- aggregate(hr_index ~ time_s, data = hr[hr$instruction == "Downregulate" & hr$valence == "Negative", ], mean)
mean_down_neu_hr <- aggregate(hr_index ~ time_s, data = hr[hr$instruction == "Downregulate" & hr$valence == "Neutral", ], mean)
mean_up_neg_hr <- aggregate(hr_index ~ time_s, data = hr[hr$instruction == "Upregulate" & hr$valence == "Negative", ], mean)
mean_up_neu_hr <- aggregate(hr_index ~ time_s, data = hr[hr$instruction == "Upregulate" & hr$valence == "Neutral", ], mean)
mean_down_neg_hr_rel <- aggregate(hr_rel ~ time_s, data = hr[hr$instruction == "Downregulate" & hr$valence == "Negative", ], mean)
mean_down_neu_hr_rel <- aggregate(hr_rel ~ time_s, data = hr[hr$instruction == "Downregulate" & hr$valence == "Neutral", ], mean)
mean_up_neg_hr_rel <- aggregate(hr_rel ~ time_s, data = hr[hr$instruction == "Upregulate" & hr$valence == "Negative", ], mean)
mean_up_neu_hr_rel <- aggregate(hr_rel ~ time_s, data = hr[hr$instruction == "Upregulate" & hr$valence == "Neutral", ], mean)

mean_timecourses_hr <- cbind(mean_down_neg_hr, mean_down_neu_hr[, "hr_index"], mean_up_neg_hr[, "hr_index"], mean_up_neu_hr[, "hr_index"], mean_down_neg_hr_rel[, "hr_rel"], mean_down_neu_hr_rel[, "hr_rel"], mean_up_neg_hr_rel[, "hr_rel"], mean_up_neu_hr_rel[, "hr_rel"])
names(mean_timecourses_hr) <- c("time_s", "down_neg", "down_neu", "up_neg", "up_neu", "down_neg_rel", "down_neu_rel", "up_neg_rel", "up_neu_rel")
write.csv(mean_timecourses_hr, file = "MeanTimecoursesHR_ER.csv", row.names = F)

# Plot time courses
plot(mean_down_neg_hr, type = "n", frame.plot = F, ylim = c(0.94, 1.02))
abline(v = c(0, 2, 3, 5), col = "gray")
lines(mean_down_neg_hr, col = "red")
lines(mean_down_neu_hr, col = "pink")
lines(mean_up_neg_hr, col = "blue")
lines(mean_up_neu_hr, col = "lightblue")

plot(mean_down_neg_hr_rel, type = "n", frame.plot = F, ylim = c(-5, 1))
abline(v = c(0, 2, 3, 5), col = "gray")
lines(mean_down_neg_hr_rel, col = "red")
lines(mean_down_neu_hr_rel, col = "pink")
lines(mean_up_neg_hr_rel, col = "blue")
lines(mean_up_neu_hr_rel, col = "lightblue")

# Average data over interval for each event for statistical modelling
hr_event_data <- data.frame()
for(i in unique(hr$subject)){
  for(j in unique(hr$event_no[hr$subject == i])){
    for(k in c("Downregulate", "Upregulate")){
      for(l in c("Neutral", "Negative")){
        hr_mean_stimulus <- mean(hr$hr_index[hr$subject == i & hr$event_no == j & hr$instruction == k & hr$valence == l][41:80])
        hr_mean_stimulus2 <- mean(hr$hr_index[hr$subject == i & hr$event_no == j & hr$instruction == k & hr$valence == l][61:80])
        hr_mean_stimulus3 <- mean(hr$hr_index[hr$subject == i & hr$event_no == j & hr$instruction == k & hr$valence == l][61:70])
        hr_mean_stimulus4 <- mean(hr$hr_index[hr$subject == i & hr$event_no == j & hr$instruction == k & hr$valence == l][51:70])
        hr_mean_stimulus_rel <- mean(hr$hr_rel[hr$subject == i & hr$event_no == j & hr$instruction == k & hr$valence == l][41:80])
        hr_mean_stimulus2_rel <- mean(hr$hr_rel[hr$subject == i & hr$event_no == j & hr$instruction == k & hr$valence == l][61:80])
        hr_mean_stimulus3_rel <- mean(hr$hr_rel[hr$subject == i & hr$event_no == j & hr$instruction == k & hr$valence == l][61:70])
        hr_mean_stimulus4_rel <- mean(hr$hr_rel[hr$subject == i & hr$event_no == j & hr$instruction == k & hr$valence == l][51:70])
        hr_event_data <- rbind(hr_event_data, cbind(data.frame(i, j, k, l), hr_mean_stimulus, hr_mean_stimulus2, hr_mean_stimulus3, hr_mean_stimulus4, hr_mean_stimulus_rel, hr_mean_stimulus2_rel, hr_mean_stimulus3_rel, hr_mean_stimulus4_rel))
      }
    }
  }
}
names(hr_event_data) <- c("Subject", "event_no", "instruction", "valence", "hr_mean_stimulus", "hr_mean_stimulus2", "hr_mean_stimulus3", "hr_mean_stimulus4", "hr_mean_stimulus_rel", "hr_mean_stimulus2_rel", "hr_mean_stimulus3_rel", "hr_mean_stimulus4_rel")

hist(hr_event_data$hr_mean_stimulus)
hist(sqrt(hr_event_data$hr_mean_stimulus))

hist(hr_event_data$hr_mean_stimulus_rel)
hist(sqrt(hr_event_data$hr_mean_stimulus_rel))

hist(hr_event_data$hr_mean_stimulus4)
hist(sqrt(hr_event_data$hr_mean_stimulus4))

hr_event_data$sqrt_hr_mean_stimulus <- sqrt(hr_event_data$hr_mean_stimulus)
hr_event_data$sqrt_hr_mean_stimulus2 <- sqrt(hr_event_data$hr_mean_stimulus2)
hr_event_data$sqrt_hr_mean_stimulus3 <- sqrt(hr_event_data$hr_mean_stimulus3)
hr_event_data$sqrt_hr_mean_stimulus4 <- sqrt(hr_event_data$hr_mean_stimulus4)

write.csv(hr_event_data, file = "HeartRateEventData_ER.csv", row.names = F)
