# Script to analyse EMG in the FMOV experiment of the Oxazepam project 2011
# Gustav Nilsonne 2015-03-04

# Works up to and including participant 85

# Require packages
require(RCurl) # To read data from GitHub
require(reshape2)
require(RColorBrewer)
require(nlme)
require(effects)
require(latticeExtra) # To make dotcharts

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
EMGDataList <- lapply(IncludedSubjects, FUN = fun_readEMGdata)


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
pdf("Fig_EMG1.pdf", width = 4, height = 4)
plot(MeanEMGCorrAngryWave1, type = "n", col = col4, frame.plot = F, main = "Corrugator EMG, Wave 1", xlab = "Time (s)", ylab = "EMG (ratio)", xaxt = "n", yaxt = "n", ylim = c(0.8, 1.3))
rect(2, 0.000001, 4, 10, col = "gray88", border = NA)
abline(v = c(0, 2, 5), lty = 3)
axis(1)
axis(2, at = c(0.8, 1, 1.2))
lines(MeanEMGCorrAngryWave1, col = col4, lwd = 2)
lines(MeanEMGCorrHappyWave1, col = col7, lwd = 2)
lines(MeanEMGCorrNeutralWave1, col = col8, lwd = 2)
legend("topleft", lty = 1, col = c(col4, col7, col8), legend = c("Angry", "Happy", "Neutral"), bty = "n", lwd = 2)
dev.off()

pdf("Fig_EMG3.pdf", width = 4, height = 4)
plot(MeanEMGCorrAngryWave2, type = "n", col = col4, frame.plot = F, main = "Corrugator EMG, Wave 2", xlab = "Time (s)", ylab = "EMG (ratio)", xaxt = "n", yaxt = "n", ylim = c(0.8, 1.3))
rect(2, 0.000001, 4, 10, col = "gray88", border = NA)
abline(v = c(0, 2, 4), lty = 3)
axis(1)
axis(2, at = c(0.8, 1, 1.2))
lines(MeanEMGCorrAngryWave2, col = col4, lwd = 2)
lines(MeanEMGCorrHappyWave2, col = col7, lwd = 2)
lines(MeanEMGCorrNeutralWave2, col = col8, lwd = 2)
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
EMGEventData$Stimulus <- ordered(EMGEventData$Stimulus, levels = c("Angry", "Neutral", "Happy"))
EMGEventData$IRI_EC_z <- scale(EMGEventData$IRI_EC)

# Build model
lme1 <- lme(EMG_corr_mean ~ Treatment*Stimulus + Wave, data = EMGEventData, random = ~1|Subject, na.action = na.omit)
lme2 <- lme(EMG_corr_mean ~ Treatment*Stimulus + Wave + IRI_EC_z, data = subset(EMGEventData, Stimulus %in% c("Angry", "Neutral")), random = ~1|Subject, na.action = na.omit)

plot(lme1)
summary(lme1)
intervals(lme1)

eff1 <- effect("Treatment*Stimulus", lme1)
eff2 <- effect("Treatment*Stimulus", lme2)

# Compare plots to less custom-generated output for verification
plot(effect("Treatment*Stimulus", lme1))
plot(effect("Treatment*Stimulus", lme2))

pdf("Fig_EMG5.pdf", width = 4, height = 4)
plot(c(eff1$fit[2], eff1$fit[4], eff1$fit[6]),
     type = "b",
     frame.plot = F,
     ylab = "EMG (ratio)",
     xlab = "Stimulus type",
     xaxt = "n",
     yaxt = "n",
     xlim = c(1, 3.1),
     ylim = c(-0.2, 0.1),
     col = col1,
     main = "Corrugator EMG"
)
lines(c(1, 1), c(eff1$upper[2], eff1$lower[2]), col = col1)
lines(c(2, 2), c(eff1$upper[4], eff1$lower[4]), col = col1)
lines(c(3, 3), c(eff1$upper[6], eff1$lower[6]), col = col1)
lines(c(1.1, 2.1, 3.1), c(eff1$fit[1], eff1$fit[3], eff1$fit[5]), type = "b", col = col2, pch = 16)
lines(c(1.1, 1.1), c(eff1$upper[1], eff1$lower[1]), col = col2)
lines(c(2.1, 2.1), c(eff1$upper[3], eff1$lower[3]), col = col2)
lines(c(3.1, 3.1), c(eff1$upper[5], eff1$lower[5]), col = col2)
axis(1, at = c(1.05, 2.05, 3.05), labels = c("Angry", "Neutral", "Happy"))
axis(2, at = c(-0.2, -0.1, 0, 0.1))
legend("topright", col = c(col1, col2), pch = c(1, 16), legend = c("Placebo", "Oxazepam"), bty = "n", lty = 1)
dev.off()

# Analyse other rating scales as predictors
# Since rating scales may be collinear, we include them one by one

# IRI-EC
EMGEventData$IRI_EC_z <- scale(EMGEventData$IRI_EC)

# IRI-PT
EMGEventData$IRI_PT_z <- scale(EMGEventData$IRI_PT)
EMGEventData$IRI_PT_z_OtherHigh <- EMGEventData$IRI_PT_z * EMGEventData$OtherHigh
EMGEventData$IRI_PT_z_OtherHigh[EMGEventData$IRI_PT_z_OtherHigh > 0 & !is.na(EMGEventData$IRI_PT_z_OtherHigh)] <- EMGEventData$IRI_PT_z_OtherHigh[EMGEventData$IRI_PT_z_OtherHigh > 0 & !is.na(EMGEventData$IRI_PT_z_OtherHigh)] - mean(EMGEventData$IRI_PT_z_OtherHigh[EMGEventData$IRI_PT_z_OtherHigh > 0 & !is.na(EMGEventData$IRI_PT_z_OtherHigh)], na.rm = TRUE)
lme2 <- lme(EMG_corr_mean ~ Treatment*Stimulus*Condition + IRI_PT_z + IRI_PT_z_OtherHigh + Wave, data = EMGEventData, random = ~1|Subject, na.action = na.omit)
plot(lme2)
summary(lme2)
intervals(lme2)

# IRI-PD
EMGEventData$IRI_PD_z <- scale(EMGEventData$IRI_PD)
EMGEventData$IRI_PD_z_OtherHigh <- EMGEventData$IRI_PD_z * EMGEventData$OtherHigh
EMGEventData$IRI_PD_z_OtherHigh[EMGEventData$IRI_PD_z_OtherHigh > 0 & !is.na(EMGEventData$IRI_PD_z_OtherHigh)] <- EMGEventData$IRI_PD_z_OtherHigh[EMGEventData$IRI_PD_z_OtherHigh > 0 & !is.na(EMGEventData$IRI_PD_z_OtherHigh)] - mean(EMGEventData$IRI_PD_z_OtherHigh[EMGEventData$IRI_PD_z_OtherHigh > 0 & !is.na(EMGEventData$IRI_PD_z_OtherHigh)], na.rm = TRUE)
lme3 <- lme(EMG_corr_mean ~ Treatment*Stimulus*Condition + IRI_PD_z + IRI_PD_z_OtherHigh + Wave, data = EMGEventData, random = ~1|Subject, na.action = na.omit)
plot(lme3)
summary(lme3)
intervals(lme3)

# IRI-F
EMGEventData$IRI_F_z <- scale(EMGEventData$IRI_F)
EMGEventData$IRI_F_z_OtherHigh <- EMGEventData$IRI_F_z * EMGEventData$OtherHigh
EMGEventData$IRI_F_z_OtherHigh[EMGEventData$IRI_F_z_OtherHigh > 0 & !is.na(EMGEventData$IRI_F_z_OtherHigh)] <- EMGEventData$IRI_F_z_OtherHigh[EMGEventData$IRI_F_z_OtherHigh > 0 & !is.na(EMGEventData$IRI_F_z_OtherHigh)] - mean(EMGEventData$IRI_F_z_OtherHigh[EMGEventData$IRI_F_z_OtherHigh > 0 & !is.na(EMGEventData$IRI_F_z_OtherHigh)], na.rm = TRUE)
lme4 <- lme(EMG_corr_mean ~ Treatment*Stimulus*Condition + IRI_F_z + IRI_F_z_OtherHigh + Wave, data = EMGEventData, random = ~1|Subject, na.action = na.omit)
plot(lme4)
summary(lme4)
intervals(lme4)

# STAI-T
EMGEventData$STAI.T_z <- scale(EMGEventData$STAI.T)
EMGEventData$STAI.T_z_OtherHigh <- EMGEventData$STAI.T_z * EMGEventData$OtherHigh
EMGEventData$STAI.T_z_OtherHigh[EMGEventData$STAI.T_z_OtherHigh > 0 & !is.na(EMGEventData$STAI.T_z_OtherHigh)] <- EMGEventData$STAI.T_z_OtherHigh[EMGEventData$STAI.T_z_OtherHigh > 0 & !is.na(EMGEventData$STAI.T_z_OtherHigh)] - mean(EMGEventData$STAI.T_z_OtherHigh[EMGEventData$STAI.T_z_OtherHigh > 0 & !is.na(EMGEventData$STAI.T_z_OtherHigh)], na.rm = TRUE)
lme5 <- lme(EMG_corr_mean ~ Treatment*Stimulus*Condition + STAI.T_z + STAI.T_z_OtherHigh + Wave, data = EMGEventData, random = ~1|Subject, na.action = na.omit)
plot(lme5)
summary(lme5)
intervals(lme5)

# TAS-20
EMGEventData$TAS.20_z <- scale(EMGEventData$TAS.20)
EMGEventData$TAS.20_z_OtherHigh <- EMGEventData$TAS.20_z * EMGEventData$OtherHigh
EMGEventData$TAS.20_z_OtherHigh[EMGEventData$TAS.20_z_OtherHigh > 0 & !is.na(EMGEventData$TAS.20_z_OtherHigh)] <- EMGEventData$TAS.20_z_OtherHigh[EMGEventData$TAS.20_z_OtherHigh > 0 & !is.na(EMGEventData$TAS.20_z_OtherHigh)] - mean(EMGEventData$TAS.20_z_OtherHigh[EMGEventData$TAS.20_z_OtherHigh > 0 & !is.na(EMGEventData$TAS.20_z_OtherHigh)], na.rm = TRUE)
lme6 <- lme(EMG_corr_mean ~ Treatment*Stimulus*Condition + TAS.20_z + TAS.20_z_OtherHigh + Wave, data = EMGEventData, random = ~1|Subject, na.action = na.omit)
plot(lme6)
summary(lme6)
intervals(lme6)

# PPI-R-SCI
EMGEventData$PPI_SCI_z <- EMGEventData$PPI_1_SCI_R
EMGEventData$PPI_SCI_z[EMGEventData$PPI_1_IR >= 45] <- NA # Exclude participants with too high scores on Inconsistent Responding measure
EMGEventData$PPI_SCI_z <- scale(EMGEventData$PPI_SCI_z)
EMGEventData$PPI_SCI_z_OtherHigh <- EMGEventData$PPI_SCI_z * EMGEventData$OtherHigh
EMGEventData$PPI_SCI_z_OtherHigh[EMGEventData$PPI_SCI_z_OtherHigh > 0 & !is.na(EMGEventData$PPI_SCI_z_OtherHigh)] <- EMGEventData$PPI_SCI_z_OtherHigh[EMGEventData$PPI_SCI_z_OtherHigh > 0 & !is.na(EMGEventData$PPI_SCI_z_OtherHigh)] - mean(EMGEventData$PPI_SCI_z_OtherHigh[EMGEventData$PPI_SCI_z_OtherHigh > 0 & !is.na(EMGEventData$PPI_SCI_z_OtherHigh)], na.rm = TRUE)
lme7 <- lme(EMG_corr_mean ~ Treatment*Stimulus*Condition + PPI_SCI_z + PPI_SCI_z_OtherHigh + Wave, data = EMGEventData, random = ~1|Subject, na.action = na.omit)
plot(lme7)
summary(lme7)
intervals(lme7)

# PPI-R-FD
EMGEventData$PPI_FD_z <- EMGEventData$PPI_1_FD_R
EMGEventData$PPI_FD_z[EMGEventData$PPI_1_IR >= 45] <- NA # Exclude participants with too high scores on Inconsistent Responding measure
EMGEventData$PPI_FD_z <- scale(EMGEventData$PPI_FD_z)
EMGEventData$PPI_FD_z_OtherHigh <- EMGEventData$PPI_FD_z * EMGEventData$OtherHigh
EMGEventData$PPI_FD_z_OtherHigh[EMGEventData$PPI_FD_z_OtherHigh > 0 & !is.na(EMGEventData$PPI_FD_z_OtherHigh)] <- EMGEventData$PPI_FD_z_OtherHigh[EMGEventData$PPI_FD_z_OtherHigh > 0 & !is.na(EMGEventData$PPI_FD_z_OtherHigh)] - mean(EMGEventData$PPI_FD_z_OtherHigh[EMGEventData$PPI_FD_z_OtherHigh > 0 & !is.na(EMGEventData$PPI_FD_z_OtherHigh)], na.rm = TRUE)
lme8 <- lme(EMG_corr_mean ~ Treatment*Stimulus*Condition + PPI_FD_z + PPI_FD_z_OtherHigh + Wave, data = EMGEventData, random = ~1|Subject, na.action = na.omit)
plot(lme8)
summary(lme8)
intervals(lme8)

# PPI-R-C
EMGEventData$PPI_C_z <- EMGEventData$PPI_1_C_R
EMGEventData$PPI_C_z[EMGEventData$PPI_1_IR >= 45] <- NA # Exclude participants with too high scores on Inconsistent Responding measure
EMGEventData$PPI_C_z <- scale(EMGEventData$PPI_C_z)
EMGEventData$PPI_C_z_OtherHigh <- EMGEventData$PPI_C_z * EMGEventData$OtherHigh
EMGEventData$PPI_C_z_OtherHigh[EMGEventData$PPI_C_z_OtherHigh > 0 & !is.na(EMGEventData$PPI_C_z_OtherHigh)] <- EMGEventData$PPI_C_z_OtherHigh[EMGEventData$PPI_C_z_OtherHigh > 0 & !is.na(EMGEventData$PPI_C_z_OtherHigh)] - mean(EMGEventData$PPI_C_z_OtherHigh[EMGEventData$PPI_C_z_OtherHigh > 0 & !is.na(EMGEventData$PPI_C_z_OtherHigh)], na.rm = TRUE)
lme9 <- lme(EMG_corr_mean ~ Treatment*Stimulus*Condition + PPI_C_z + PPI_C_z_OtherHigh + Wave, data = EMGEventData, random = ~1|Subject, na.action = na.omit)
plot(lme9)
summary(lme9)
intervals(lme9)

# Make plots to compare effects for different scales
# Put data in new frame
data_main <- data.frame(scale = "PPI-R-C", beta = intervals(lme9)$fixed[5, 2], lower = intervals(lme9)$fixed[5, 1], upper = intervals(lme9)$fixed[5, 3], group = "PPI")
data_main <- rbind(data_main, data.frame(scale = "PPI-R-FD", beta = intervals(lme8)$fixed[5, 2], lower = intervals(lme8)$fixed[5, 1], upper = intervals(lme8)$fixed[5, 3], group = "PPI"))
data_main <- rbind(data_main, data.frame(scale = "PPI-R-SCI", beta = intervals(lme7)$fixed[5, 2], lower = intervals(lme7)$fixed[5, 1], upper = intervals(lme7)$fixed[5, 3], group = "PPI"))
data_main <- rbind(data_main, data.frame(scale = "TAS-20", beta = intervals(lme6)$fixed[5, 2], lower = intervals(lme6)$fixed[5, 1], upper = intervals(lme6)$fixed[5, 3], group = "TAS"))
data_main <- rbind(data_main, data.frame(scale = "STAI-T", beta = intervals(lme5)$fixed[5, 2], lower = intervals(lme5)$fixed[5, 1], upper = intervals(lme5)$fixed[5, 3], group = "STAI"))
data_main <- rbind(data_main, data.frame(scale = "IRI-F", beta = intervals(lme4)$fixed[5, 2], lower = intervals(lme4)$fixed[5, 1], upper = intervals(lme4)$fixed[5, 3], group = "IRI"))
data_main <- rbind(data_main, data.frame(scale = "IRI-PD", beta = intervals(lme3)$fixed[5, 2], lower = intervals(lme3)$fixed[5, 1], upper = intervals(lme3)$fixed[5, 3], group = "IRI"))
data_main <- rbind(data_main, data.frame(scale = "IRI-PT", beta = intervals(lme2)$fixed[5, 2], lower = intervals(lme2)$fixed[5, 1], upper = intervals(lme2)$fixed[5, 3], group = "IRI"))
data_main <- rbind(data_main, data.frame(scale = "IRI-EC", beta = intervals(lme1)$fixed[5, 2], lower = intervals(lme1)$fixed[5, 1], upper = intervals(lme1)$fixed[5, 3], group = "IRI"))

# Make plot
pdf("Fig_EMG7.pdf", width = 4, height = 4)
axis.L <- function(side, ..., line.col){
  if (side %in% c("bottom", "left")) {
    col <- trellis.par.get("axis.text")$col
    axis.default(side, ..., line.col = col)
    if (side == "bottom")
      grid::grid.lines(y = 0)
  }
}

sty <- list()
sty$axis.line$col <- NA
sty$strip.border$col <- NA
sty$strip.background$col <- NA
segplot(scale ~ lower + upper, data = data_main, 
        centers = beta, 
        lwd = 2,
        draw.bands = FALSE,
        col = c(col8, col8, col8, col3, col6, col5, col5, col5, col5),
        par.settings = sty,
        axis=axis.L,
        xlab = "Beta, 95% CI",
        main = "D. Corrugator EMG",
        panel = function (x,y,z,...){
          panel.segplot(x,y,z,...)
          panel.abline(v=0,lty=2)
        })
dev.off()

# Put data in new frame
data_emp <- data.frame(scale = "PPI-R-C", beta = intervals(lme9)$fixed[6, 2], lower = intervals(lme9)$fixed[6, 1], upper = intervals(lme9)$fixed[6, 3], group = "PPI")
data_emp <- rbind(data_emp, data.frame(scale = "PPI-R-FD", beta = intervals(lme8)$fixed[6, 2], lower = intervals(lme8)$fixed[6, 1], upper = intervals(lme8)$fixed[6, 3], group = "PPI"))
data_emp <- rbind(data_emp, data.frame(scale = "PPI-R-SCI", beta = intervals(lme7)$fixed[6, 2], lower = intervals(lme7)$fixed[6, 1], upper = intervals(lme7)$fixed[6, 3], group = "PPI"))
data_emp <- rbind(data_emp, data.frame(scale = "TAS-20", beta = intervals(lme6)$fixed[6, 2], lower = intervals(lme6)$fixed[6, 1], upper = intervals(lme6)$fixed[6, 3], group = "TAS"))
data_emp <- rbind(data_emp, data.frame(scale = "STAI-T", beta = intervals(lme5)$fixed[6, 2], lower = intervals(lme5)$fixed[6, 1], upper = intervals(lme5)$fixed[6, 3], group = "STAI"))
data_emp <- rbind(data_emp, data.frame(scale = "IRI-F", beta = intervals(lme4)$fixed[6, 2], lower = intervals(lme4)$fixed[6, 1], upper = intervals(lme4)$fixed[6, 3], group = "IRI"))
data_emp <- rbind(data_emp, data.frame(scale = "IRI-PD", beta = intervals(lme3)$fixed[6, 2], lower = intervals(lme3)$fixed[6, 1], upper = intervals(lme3)$fixed[6, 3], group = "IRI"))
data_emp <- rbind(data_emp, data.frame(scale = "IRI-PT", beta = intervals(lme2)$fixed[6, 2], lower = intervals(lme2)$fixed[6, 1], upper = intervals(lme2)$fixed[6, 3], group = "IRI"))
data_emp <- rbind(data_emp, data.frame(scale = "IRI-EC", beta = intervals(lme1)$fixed[6, 2], lower = intervals(lme1)$fixed[6, 1], upper = intervals(lme1)$fixed[6, 3], group = "IRI"))

# Make plot
pdf("Fig_EMG8.pdf", width = 4, height = 4)
axis.L <- function(side, ..., line.col){
  if (side %in% c("bottom", "left")) {
    col <- trellis.par.get("axis.text")$col
    axis.default(side, ..., line.col = col)
    if (side == "bottom")
      grid::grid.lines(y = 0)
  }
}

sty <- list()
sty$axis.line$col <- NA
sty$strip.border$col <- NA
sty$strip.background$col <- NA
segplot(scale ~ lower + upper, data = data_emp, 
        centers = beta, 
        lwd = 2,
        draw.bands = FALSE,
        col = c(col8, col8, col8, col3, col6, col5, col5, col5, col5),
        par.settings = sty,
        axis=axis.L,
        xlab = "Beta, 95% CI",
        main = "D. Corrugator EMG",
        panel = function (x,y,z,...){
          panel.segplot(x,y,z,...)
          panel.abline(v=0,lty=2)
        })
dev.off()

# Analyse zygomatic responses
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
pdf("Fig_EMG1.pdf", width = 4, height = 4)
plot(MeanEMGzygAngryWave1, type = "n", col = col4, frame.plot = F, main = "Zygomatic EMG, Wave 1", xlab = "Time (s)", ylab = "EMG (ratio)", xaxt = "n", yaxt = "n", ylim = c(0.8, 1.5))
rect(2, 0.000001, 4, 10, col = "gray88", border = NA)
abline(v = c(0, 2, 5), lty = 3)
axis(1)
axis(2, at = c(0.8, 1, 1.2, 1.4))
lines(MeanEMGzygAngryWave1, col = col4, lwd = 2)
lines(MeanEMGzygHappyWave1, col = col7, lwd = 2)
lines(MeanEMGzygNeutralWave1, col = col8, lwd = 2)
legend("topleft", lty = 1, col = c(col4, col7, col8), legend = c("Angry", "Happy", "Neutral"), bty = "n", lwd = 2)
dev.off()

pdf("Fig_EMG3.pdf", width = 4, height = 4)
plot(MeanEMGzygAngryWave2, type = "n", col = col4, frame.plot = F, main = "Zygomatic EMG, Wave 2", xlab = "Time (s)", ylab = "EMG (ratio)", xaxt = "n", yaxt = "n", ylim = c(0.8, 1.5))
rect(2, 0.000001, 4, 10, col = "gray88", border = NA)
abline(v = c(0, 2, 4), lty = 3)
axis(1)
axis(2, at = c(0.8, 1, 1.2, 1.4))
lines(MeanEMGzygAngryWave2, col = col4, lwd = 2)
lines(MeanEMGzygHappyWave2, col = col7, lwd = 2)
lines(MeanEMGzygNeutralWave2, col = col8, lwd = 2)
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
EMGEventData2$Stimulus <- ordered(EMGEventData2$Stimulus, levels = c("Angry", "Neutral", "Happy"))

# Build model
EMGEventData2$IRI_EC_z <- scale(EMGEventData2$IRI_EC)

lmezyg1 <- lme(EMG_zyg_mean ~ Treatment*Stimulus + Wave, data = EMGEventData2, random = ~1|Subject, na.action = na.omit)
lmezyg2 <- lme(EMG_zyg_mean ~ Treatment + Wave, data = subset(EMGEventData2, Stimulus %in% c("Angry")), random = ~1|Subject, na.action = na.omit)

plot(lmezyg1)
summary(lmezyg1)
intervals(lmezyg1)

summary(lmezyg2)

effzyg1 <- effect("Treatment*Stimulus", lmezyg1)
#effzyg2 <- effect("Treatment*Stimulus", lmezyg2)

# Compare plots to less custom-generated output for verification
plot(effect("Treatment*Stimulus", lmezyg1))
plot(effect("Treatment*Stimulus", lmezyg2))

# Plot results
plot(c(effzyg1$fit[2], effzyg1$fit[4], effzyg1$fit[6]),
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
lines(c(1, 1), c(effzyg1$upper[2], effzyg1$lower[2]), col = col1)
lines(c(2, 2), c(effzyg1$upper[4], effzyg1$lower[4]), col = col1)
lines(c(3, 3), c(effzyg1$upper[6], effzyg1$lower[6]), col = col1)
lines(c(1.1, 2.1, 3.1), c(effzyg1$fit[1], effzyg1$fit[3], effzyg1$fit[5]), type = "b", col = col2, pch = 16)
lines(c(1.1, 1.1), c(effzyg1$upper[1], effzyg1$lower[1]), col = col2)
lines(c(2.1, 2.1), c(effzyg1$upper[3], effzyg1$lower[3]), col = col2)
lines(c(3.1, 3.1), c(effzyg1$upper[5], effzyg1$lower[5]), col = col2)
axis(1, at = c(1.05, 2.05, 3.05), labels = c("Angry", "Neutral", "Happy"))
axis(2, at = c(-0.2, -0.1, 0, 0.1))
legend("topleft", col = c(col1, col2), pch = c(1, 16), legend = c("Placebo", "Oxazepam"), bty = "n", lty = 1)

