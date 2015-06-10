# Script to analyse EMG in the FMOV experiment of the Oxazepam project 2011
# Gustav Nilsonne 2015-06-09

# Require packages
require(RCurl) # To read data from GitHub
require(RColorBrewer)
require(nlme)
require(effects)
require(psych)

# Define colors
col1 = brewer.pal(8, "Dark2")[1]
col2 = brewer.pal(8, "Dark2")[2]
col3 = brewer.pal(8, "Dark2")[3]
col4 = brewer.pal(8, "Dark2")[4]
col5 = brewer.pal(8, "Dark2")[5]
col6 = brewer.pal(8, "Dark2")[6]
col7 = brewer.pal(8, "Dark2")[7]
col8 = brewer.pal(8, "Dark2")[8]

# Read files

# First demographic data etc
demDataURL <- getURL("https://raw.githubusercontent.com/GNilsonne/Data-and-analysis-code-Oxazepam-and-emotion/master/demographics.csv", ssl.verifypeer = FALSE)
demData <- read.csv(text = demDataURL)

# Then the preprocessed EMG data
CorrelationZygomaticCorrugatorURL <- getURL("https://raw.githubusercontent.com/GNilsonne/Data-and-analysis-code-Oxazepam-and-emotion-2/master/CorrelationZygomaticCorrugator.csv", ssl.verifypeer = FALSE)
CorrelationZygomaticCorrugator <- read.csv(text = CorrelationZygomaticCorrugatorURL)

MeanTimecoursesCorrugatorURL <- getURL("https://raw.githubusercontent.com/GNilsonne/Data-and-analysis-code-Oxazepam-and-emotion-2/master/MeanTimecoursesCorrugator.csv", ssl.verifypeer = FALSE)
MeanTimecoursesCorrugator <- read.csv(text = MeanTimecoursesCorrugatorURL)

CorrugatorEventDataURL <- getURL("https://raw.githubusercontent.com/GNilsonne/Data-and-analysis-code-Oxazepam-and-emotion-2/master/CorrugatorEventData.csv", ssl.verifypeer = FALSE)
CorrugatorEventData <- read.csv(text = CorrugatorEventDataURL)

MeanTimecoursesZygomaticURL <- getURL("https://raw.githubusercontent.com/GNilsonne/Data-and-analysis-code-Oxazepam-and-emotion-2/master/MeanTimecoursesZygomatic.csv", ssl.verifypeer = FALSE)
MeanTimecoursesZygomatic <- read.csv(text = MeanTimecoursesZygomaticURL)

ZygomaticEventDataURL <- getURL("https://raw.githubusercontent.com/GNilsonne/Data-and-analysis-code-Oxazepam-and-emotion-2/master/ZygomaticEventData.csv", ssl.verifypeer = FALSE)
ZygomaticEventData <- read.csv(text = ZygomaticEventDataURL)


# Analyse data

# Correlation between zygomatic and corrugator activity
fisherz2r(mean(CorrelationZygomaticCorrugator$z))
fisherz2r(sd(CorrelationZygomaticCorrugator$z))

# Corrugator activity
# Make figures
pdf("Fig_EMG1.pdf", width = 7, height = 5)
plot(MeanTimecoursesCorrugator$AngryWave1 ~ MeanTimecoursesCorrugator$time_s, type = "n", col = col4, frame.plot = F, main = "Corrugator EMG, Wave 1", xlab = "Time (s)", ylab = "EMG (ratio)", xaxt = "n", yaxt = "n", ylim = c(0.85, 1.2))
rect(2, 0.000001, 4, 10, col = "gray88", border = NA)
abline(v = c(0, 2, 5), lty = 3)
axis(1)
axis(2, at = c(0.9, 1, 1.1, 1.2))
lines(lowess(x = MeanTimecoursesCorrugator$time_s, y = MeanTimecoursesCorrugator$AngryWave1, f = 0.1), col = col4, lwd = 2)
lines(lowess(x = MeanTimecoursesCorrugator$time_s, y = MeanTimecoursesCorrugator$HappyWave1, f = 0.1), col = col7, lwd = 2)
lines(lowess(x = MeanTimecoursesCorrugator$time_s, y = MeanTimecoursesCorrugator$NeutralWave1, f= 0.1), col = col8, lwd = 2)
legend("topleft", lty = 1, col = c(col4, col7, col8), legend = c("Angry", "Happy", "Neutral"), lwd = 2, bty = "n")
dev.off()

pdf("Fig_EMG2.pdf", width = 7, height = 5)
plot(MeanEMGCorrAngryWave2, type = "n", col = col4, frame.plot = F, main = "Corrugator EMG, Wave 2", xlab = "Time (s)", ylab = "EMG (ratio)", xaxt = "n", yaxt = "n", ylim = c(0.85, 1.2))
rect(2, 0.000001, 4, 10, col = "gray88", border = NA)
abline(v = c(0, 2, 4), lty = 3)
axis(1)
axis(2, at = c(0.9, 1, 1.1, 1.2))
lines(lowess(x = MeanTimecoursesCorrugator$time_s, y = MeanTimecoursesCorrugator$AngryWave2, f = 0.1), col = col4, lwd = 2)
lines(lowess(x = MeanTimecoursesCorrugator$time_s, y = MeanTimecoursesCorrugator$HappyWave2, f = 0.1), col = col7, lwd = 2)
lines(lowess(x = MeanTimecoursesCorrugator$time_s, y = MeanTimecoursesCorrugator$NeutralWave2, f= 0.1), col = col8, lwd = 2)
dev.off()

# Prepare data for statistical modelling
CorrugatorEventData <- merge(CorrugatorEventData, demData, by = "Subject")
CorrugatorEventData$Subject <- as.factor(CorrugatorEventData$Subject)
CorrugatorEventData$Stimulus <- relevel(CorrugatorEventData$Stimulus, ref = "Neutral")
CorrugatorEventData$Treatment <- relevel(CorrugatorEventData$Treatment, ref = "Placebo")
CorrugatorEventData$IRI_EC_z <- scale(CorrugatorEventData$IRI_EC)

# Build model
lme1 <- lme(EMG_corr_mean ~ Stimulus*(Treatment + IRI_EC_z) + Wave, data = CorrugatorEventData, random = ~1|Subject, na.action = na.omit)
lme1b <- lme(EMG_corr_mean ~ Stimulus*Treatment + Wave, data = CorrugatorEventData, random = ~1|Subject, na.action = na.omit)

plot(lme1)
summary(lme1)
intervals(lme1)

eff1 <- effect("Stimulus*Treatment", lme1b)

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

# Compare plot to less custom-generated output for verification
plot(effect("Stimulus*Treatment", lme1b))

# Analyse other rating scales as predictors
# Since rating scales may be collinear, we include them one by one

# IRI-PT
CorrugatorEventData$IRI_PT_z <- scale(CorrugatorEventData$IRI_PT)
lme2 <- lme(EMG_corr_mean ~ Stimulus*(Treatment + IRI_PT_z) + Wave, data = CorrugatorEventData, random = ~1|Subject, na.action = na.omit)
plot(lme2)
summary(lme2)
intervals(lme2)

# IRI-PD
CorrugatorEventData$IRI_PD_z <- scale(CorrugatorEventData$IRI_PD)
lme3 <- lme(EMG_corr_mean ~ Stimulus*(Treatment + IRI_PD_z) + Wave, data = CorrugatorEventData, random = ~1|Subject, na.action = na.omit)
plot(lme3)
summary(lme3)
intervals(lme3)

# IRI-F
CorrugatorEventData$IRI_F_z <- scale(CorrugatorEventData$IRI_F)
lme4 <- lme(EMG_corr_mean ~ Stimulus*(Treatment + IRI_F_z) + Wave, data = CorrugatorEventData, random = ~1|Subject, na.action = na.omit)
plot(lme4)
summary(lme4)
intervals(lme4)

# STAI-T
CorrugatorEventData$STAI.T_z <- scale(CorrugatorEventData$STAI.T)
lme5 <- lme(EMG_corr_mean ~ Stimulus*(Treatment + STAI.T_z) + Wave, data = CorrugatorEventData, random = ~1|Subject, na.action = na.omit)
plot(lme5)
summary(lme5)
intervals(lme5)

# TAS-20
CorrugatorEventData$TAS.20_z <- scale(CorrugatorEventData$TAS.20)
lme6 <- lme(EMG_corr_mean ~ Stimulus*(Treatment + TAS.20_z) + Wave, data = CorrugatorEventData, random = ~1|Subject, na.action = na.omit)
plot(lme6)
summary(lme6)
intervals(lme6)

# PPI-R-SCI
CorrugatorEventData$PPI_SCI_z <- CorrugatorEventData$PPI_1_SCI_R
CorrugatorEventData$PPI_SCI_z[CorrugatorEventData$PPI_1_IR >= 45] <- NA # Exclude participants with too high scores on Inconsistent Responding measure
CorrugatorEventData$PPI_SCI_z <- scale(CorrugatorEventData$PPI_SCI_z)
lme7 <- lme(EMG_corr_mean ~ Stimulus*(Treatment + PPI_SCI_z) + Wave, data = CorrugatorEventData, random = ~1|Subject, na.action = na.omit)
plot(lme7)
summary(lme7)
intervals(lme7)

# PPI-R-FD
CorrugatorEventData$PPI_FD_z <- CorrugatorEventData$PPI_1_FD_R
CorrugatorEventData$PPI_FD_z[CorrugatorEventData$PPI_1_IR >= 45] <- NA # Exclude participants with too high scores on Inconsistent Responding measure
CorrugatorEventData$PPI_FD_z <- scale(CorrugatorEventData$PPI_FD_z)
lme8 <- lme(EMG_corr_mean ~ Stimulus*(Treatment + PPI_FD_z) + Wave, data = CorrugatorEventData, random = ~1|Subject, na.action = na.omit)
plot(lme8)
summary(lme8)
intervals(lme8)

# PPI-R-C
CorrugatorEventData$PPI_C_z <- CorrugatorEventData$PPI_1_C_R
CorrugatorEventData$PPI_C_z[CorrugatorEventData$PPI_1_IR >= 45] <- NA # Exclude participants with too high scores on Inconsistent Responding measure
CorrugatorEventData$PPI_C_z <- scale(CorrugatorEventData$PPI_C_z)
lme9 <- lme(EMG_corr_mean ~ Stimulus*(Treatment + PPI_C_z) + Wave, data = CorrugatorEventData, random = ~1|Subject, na.action = na.omit)
plot(lme9)
summary(lme9)
intervals(lme9)


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

# Analyse effects with angry as reference level
#CorrugatorEventData$Stimulus <- relevel(CorrugatorEventData$Stimulus, ref = "Angry")
#lme1c <- lme(EMG_corr_mean ~ Stimulus*(Treatment + IRI_EC_z) + Wave, data = CorrugatorEventData, random = ~1|Subject, na.action = na.omit)
#summary(lme1c)


# Analyze zygomatic responses

# Make figures
pdf("Fig_EMG6.pdf", width = 7, height = 5)
plot(MeanTimecoursesZygomatic$AngryWave1 ~ MeanTimecoursesZygomatic$time_s, type = "n", col = col4, frame.plot = F, main = "Zygomatic EMG, Wave 1", xlab = "Time (s)", ylab = "EMG (ratio)", xaxt = "n", yaxt = "n", ylim = c(0.9, 1.45))
rect(2, 0.000001, 4, 10, col = "gray88", border = NA)
abline(v = c(0, 2, 5), lty = 3)
axis(1)
axis(2, at = c(1, 1.2, 1.4))
lines(lowess(x = MeanTimecoursesZygomatic$time_s, y = MeanTimecoursesZygomatic$AngryWave1, f = 0.1), col = col4, lwd = 2)
lines(lowess(x = MeanTimecoursesZygomatic$time_s, y = MeanTimecoursesZygomatic$HappyWave1, f = 0.1), col = col7, lwd = 2)
lines(lowess(x = MeanTimecoursesZygomatic$time_s, y = MeanTimecoursesZygomatic$NeutralWave1, f= 0.1), col = col8, lwd = 2)
legend("topleft", lty = 1, col = c(col4, col7, col8), legend = c("Angry", "Happy", "Neutral"), bty = "n", lwd = 2)
dev.off()

pdf("Fig_EMG7.pdf", width = 7, height = 5)
plot(MeanTimecoursesZygomatic$AngryWave2 ~ MeanTimecoursesZygomatic$time_s, type = "n", col = col4, frame.plot = F, main = "Zygomatic EMG, Wave 2", xlab = "Time (s)", ylab = "EMG (ratio)", xaxt = "n", yaxt = "n", ylim = c(0.9, 1.45))
rect(2, 0.000001, 4, 10, col = "gray88", border = NA)
abline(v = c(0, 2, 4), lty = 3)
axis(1)
axis(2, at = c(1, 1.2, 1.4))
lines(lowess(x = MeanTimecoursesZygomatic$time_s, y = MeanTimecoursesZygomatic$AngryWave2, f = 0.1), col = col4, lwd = 2)
lines(lowess(x = MeanTimecoursesZygomatic$time_s, y = MeanTimecoursesZygomatic$HappyWave2, f = 0.1), col = col7, lwd = 2)
lines(lowess(x = MeanTimecoursesZygomatic$time_s, y = MeanTimecoursesZygomatic$NeutralWave2, f= 0.1), col = col8, lwd = 2)
dev.off()

# Prepare data for statistical modelling
ZygomaticEventData <- merge(ZygomaticEventData, demData, by = "Subject")
ZygomaticEventData$Subject <- as.factor(ZygomaticEventData$Subject)
ZygomaticEventData$Stimulus <- relevel(ZygomaticEventData$Stimulus, ref = "Neutral")
ZygomaticEventData$Treatment <- relevel(ZygomaticEventData$Treatment, ref = "Placebo")
ZygomaticEventData$IRI_EC_z <- scale(ZygomaticEventData$IRI_EC)

# Build model
lmezyg1 <- lme(EMG_zyg_mean ~ Stimulus*(Treatment + IRI_EC_z) + Wave, data = ZygomaticEventData, random = ~1|Subject, na.action = na.omit)
lmezyg1b <- lme(EMG_zyg_mean ~ Stimulus*Treatment + Wave, data = ZygomaticEventData, random = ~1|Subject, na.action = na.omit)

plot(lmezyg1)
summary(lmezyg1)
intervals(lmezyg1)

summary(lmezyg1b)

effzyg1 <- effect("Stimulus*Treatment", lmezyg1b)

# Compare plots to less custom-generated output for verification
plot(effect("Stimulus*Treatment", lmezyg1))
plot(effect("Stimulus*Treatment", lmezyg1b))

# Plot results
pdf("Fig_EMG8.pdf", height = 5, width = 7)
plot(c(effzyg1$fit[1], effzyg1$fit[2], effzyg1$fit[3]),
     type = "b",
     frame.plot = F,
     ylab = "EMG (log ratio)",
     xlab = "Stimulus type",
     xaxt = "n",
     yaxt = "n",
     xlim = c(1, 3.1),
     ylim = c(-0.2, 0.15),
     col = col1,
     main = "Zygomatic EMG"
)
lines(c(1, 1), c(effzyg1$upper[1], effzyg1$lower[1]), col = col1)
lines(c(2, 2), c(effzyg1$upper[2], effzyg1$lower[2]), col = col1)
lines(c(3, 3), c(effzyg1$upper[3], effzyg1$lower[3]), col = col1)
lines(c(1.1, 2.1, 3.1), c(effzyg1$fit[4], effzyg1$fit[5], effzyg1$fit[6]), type = "b", col = col2, pch = 16)
lines(c(1.1, 1.1), c(effzyg1$upper[4], effzyg1$lower[4]), col = col2)
lines(c(2.1, 2.1), c(effzyg1$upper[5], effzyg1$lower[5]), col = col2)
lines(c(3.1, 3.1), c(effzyg1$upper[6], effzyg1$lower[6]), col = col2)
axis(1, at = c(1.05, 2.05, 3.05), labels = c("Neutral", "Angry", "Happy"))
axis(2, at = c(-0.2, -0.1, 0, 0.1))
legend("topleft", col = c(col1, col2), pch = c(1, 16), legend = c("Placebo", "Oxazepam"), bty = "n", lty = 1)
dev.off()

# Analyse other rating scales as predictors
# Since rating scales may be collinear, we include them one by one

# IRI-PT
ZygomaticEventData$IRI_PT_z <- scale(ZygomaticEventData$IRI_PT)
lmezyg2 <- lme(EMG_zyg_mean ~ Stimulus*(Treatment + IRI_PT_z) + Wave, data = ZygomaticEventData, random = ~1|Subject, na.action = na.omit)
plot(lmezyg2)
summary(lmezyg2)
intervals(lmezyg2)

# IRI-PD
ZygomaticEventData$IRI_PD_z <- scale(ZygomaticEventData$IRI_PD)
lmezyg3 <- lme(EMG_zyg_mean ~ Stimulus*(Treatment + IRI_PD_z) + Wave, data = ZygomaticEventData, random = ~1|Subject, na.action = na.omit)
plot(lmezyg3)
summary(lmezyg3)
intervals(lmezyg3)

# IRI-F
ZygomaticEventData$IRI_F_z <- scale(ZygomaticEventData$IRI_F)
lmezyg4 <- lme(EMG_zyg_mean ~ Stimulus*(Treatment + IRI_F_z) + Wave, data = ZygomaticEventData, random = ~1|Subject, na.action = na.omit)
plot(lmezyg4)
summary(lmezyg4)
intervals(lmezyg4)

# STAI-T
ZygomaticEventData$STAI.T_z <- scale(ZygomaticEventData$STAI.T)
lmezyg5 <- lme(EMG_zyg_mean ~ Stimulus*(Treatment + STAI.T_z) + Wave, data = ZygomaticEventData, random = ~1|Subject, na.action = na.omit)
plot(lmezyg5)
summary(lmezyg5)
intervals(lmezyg5)

# TAS-20
ZygomaticEventData$TAS.20_z <- scale(ZygomaticEventData$TAS.20)
lmezyg6 <- lme(EMG_zyg_mean ~ Stimulus*(Treatment + TAS.20_z) + Wave, data = ZygomaticEventData, random = ~1|Subject, na.action = na.omit)
plot(lmezyg6)
summary(lmezyg6)
intervals(lmezyg6)

# PPI-R-SCI
ZygomaticEventData$PPI_SCI_z <- ZygomaticEventData$PPI_1_SCI_R
ZygomaticEventData$PPI_SCI_z[ZygomaticEventData$PPI_1_IR >= 45] <- NA # Exclude participants with too high scores on Inconsistent Responding measure
ZygomaticEventData$PPI_SCI_z <- scale(ZygomaticEventData$PPI_SCI_z)
lmezyg7 <- lme(EMG_zyg_mean ~ Stimulus*(Treatment + PPI_SCI_z) + Wave, data = ZygomaticEventData, random = ~1|Subject, na.action = na.omit)
plot(lmezyg7)
summary(lmezyg7)
intervals(lmezyg7)

# PPI-R-FD
ZygomaticEventData$PPI_FD_z <- ZygomaticEventData$PPI_1_FD_R
ZygomaticEventData$PPI_FD_z[ZygomaticEventData$PPI_1_IR >= 45] <- NA # Exclude participants with too high scores on Inconsistent Responding measure
ZygomaticEventData$PPI_FD_z <- scale(ZygomaticEventData$PPI_FD_z)
lmezyg8 <- lme(EMG_zyg_mean ~ Stimulus*(Treatment + PPI_FD_z) + Wave, data = ZygomaticEventData, random = ~1|Subject, na.action = na.omit)
plot(lmezyg8)
summary(lmezyg8)
intervals(lmezyg8)

# PPI-R-C
ZygomaticEventData$PPI_C_z <- ZygomaticEventData$PPI_1_C_R
ZygomaticEventData$PPI_C_z[ZygomaticEventData$PPI_1_IR >= 45] <- NA # Exclude participants with too high scores on Inconsistent Responding measure
ZygomaticEventData$PPI_C_z <- scale(ZygomaticEventData$PPI_C_z)
lmezyg9 <- lme(EMG_zyg_mean ~ Stimulus*(Treatment + PPI_C_z) + Wave, data = ZygomaticEventData, random = ~1|Subject, na.action = na.omit)
plot(lmezyg9)
summary(lmezyg9)
intervals(lmezyg9)

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

# Write regression output tables
write.csv(summary(lme1b)$tTable, file = "Corr_unadjusted.csv")
write.csv(summary(lme1)$tTable, file = "Corr_IRI_EC.csv")
write.csv(summary(lme2)$tTable, file = "Corr_IRI_PT.csv")
write.csv(summary(lme3)$tTable, file = "Corr_IRI_PD.csv")
write.csv(summary(lme4)$tTable, file = "Corr_IRI_F.csv")
write.csv(summary(lme5)$tTable, file = "Corr_STAI_T.csv")
write.csv(summary(lme6)$tTable, file = "Corr_TAS_20.csv")
write.csv(summary(lme7)$tTable, file = "Corr_PPI_R_SCI.csv")
write.csv(summary(lme8)$tTable, file = "Corr_PPI_R_FD.csv")
write.csv(summary(lme9)$tTable, file = "Corr_PPI_R_C.csv")

write.csv(summary(lmezyg1b)$tTable, file = "Zyg_unadjusted.csv")
write.csv(summary(lmezyg1)$tTable, file = "Zyg_IRI_EC.csv")
write.csv(summary(lmezyg2)$tTable, file = "Zyg_IRI_PT.csv")
write.csv(summary(lmezyg3)$tTable, file = "Zyg_IRI_PD.csv")
write.csv(summary(lmezyg4)$tTable, file = "Zyg_IRI_F.csv")
write.csv(summary(lmezyg5)$tTable, file = "Zyg_STAI_T.csv")
write.csv(summary(lmezyg6)$tTable, file = "Zyg_TAS_20.csv")
write.csv(summary(lmezyg7)$tTable, file = "Zyg_PPI_R_SCI.csv")
write.csv(summary(lmezyg8)$tTable, file = "Zyg_PPI_R_FD.csv")
write.csv(summary(lmezyg9)$tTable, file = "Zyg_PPI_R_C.csv")

# Analyse happy vs angry conditions separately
lme1br <- lme(EMG_corr_mean ~ Stimulus*Treatment + Wave, data = CorrugatorEventData[CorrugatorEventData$Stimulus != "Neutral", ], random = ~1|Subject, na.action = na.omit)
summary(lme1br)
write.csv(summary(lme1br)$tTable, file = "Corr_AngryHappy.csv")

lmezyg1br <- lme(EMG_zyg_mean ~ Stimulus*Treatment + Wave, data = ZygomaticEventData[ZygomaticEventData$Stimulus != "Neutral", ], random = ~1|Subject, na.action = na.omit)
summary(lmezyg1br)
write.csv(summary(lmezyg1br)$tTable, file = "Zyg_AngryHappy.csv")

# Analyse waves 1 and 2 separately
lme1bw1 <- lme(EMG_corr_mean ~ Stimulus*Treatment, data = CorrugatorEventData[CorrugatorEventData$Wave == 1, ], random = ~1|Subject, na.action = na.omit)
lme1bw2 <- lme(EMG_corr_mean ~ Stimulus*Treatment, data = CorrugatorEventData[CorrugatorEventData$Wave == 2, ], random = ~1|Subject, na.action = na.omit)
summary(lme1bw1)
summary(lme1bw2)
write.csv(summary(lme1bw1)$tTable, file = "Corr_wave1.csv")
write.csv(summary(lme1bw2)$tTable, file = "Corr_wave2.csv")

lmezyg1bw1 <- lme(EMG_zyg_mean ~ Stimulus*Treatment, data = ZygomaticEventData[ZygomaticEventData$Wave == 1, ], random = ~1|Subject, na.action = na.omit)
lmezyg1bw2 <- lme(EMG_zyg_mean ~ Stimulus*Treatment, data = ZygomaticEventData[ZygomaticEventData$Wave == 2, ], random = ~1|Subject, na.action = na.omit)
summary(lmezyg1bw1)
summary(lmezyg1bw2)
write.csv(summary(lmezyg1bw1)$tTable, file = "Zyg_wave1.csv")
write.csv(summary(lmezyg1bw2)$tTable, file = "Zyg_wave2.csv")

# Calculate a response index for each participant and write
Corr_agg <- aggregate(CorrugatorEventData[, c("EMG_corr_mean", "Subject", "Stimulus")], list(Subject = CorrugatorEventData$Subject, Stimulus = CorrugatorEventData$Stimulus), FUN = "mean")
Zyg_agg <- aggregate(ZygomaticEventData[, c("EMG_zyg_mean", "Subject", "Stimulus")], list(Subject = ZygomaticEventData$Subject, Stimulus = ZygomaticEventData$Stimulus), FUN = "mean")
IndividualResponses <- data.frame(Subject = Corr_agg$Subject[Corr_agg$Stimulus == "Angry"], 
                                  CorrAngry = Corr_agg$EMG_corr_mean[Corr_agg$Stimulus == "Angry"] - Corr_agg$EMG_corr_mean[Corr_agg$Stimulus == "Neutral"],
                                  CorrHappy = Corr_agg$EMG_corr_mean[Corr_agg$Stimulus == "Happy"] - Corr_agg$EMG_corr_mean[Corr_agg$Stimulus == "Neutral"],
                                  ZygAngry = Zyg_agg$EMG_zyg_mean[Zyg_agg$Stimulus == "Angry"] - Zyg_agg$EMG_zyg_mean[Zyg_agg$Stimulus == "Neutral"],
                                  ZygHappy = Zyg_agg$EMG_zyg_mean[Zyg_agg$Stimulus == "Happy"] - Zyg_agg$EMG_zyg_mean[Zyg_agg$Stimulus == "Neutral"])
write.csv(IndividualResponses, file = "IndividualResponsesFMOV.csv", row.names = FALSE)
