# Script to analyse skin conductance responses in the FMOV experiment of the Oxazepam project 2011
# Gustav Nilsonne 2015-07-17

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

# Then the preprocessed SCR data
MeanTimecourses <- read.csv("MeanTimecoursesSCR_ER.csv")
SCR_event_data <- read.csv("SCREventData_ER.csv")


# Analyse data
# Make figures
pdf("Fig_SCR_ER1.pdf", width = 7, height = 5)
plot(MeanTimecourses$down_neg ~ MeanTimecourses$time_s, type = "n", frame.plot = F, main = "Skin conductance responses, Neutral stimuli", xlab = "Time (s)", ylab = "EMG (ratio)", xaxt = "n", yaxt = "n", ylim = c(0.998, 1.02))
rect(3, 0.000001, 5, 10, col = "gray88", border = NA)
abline(v = c(0, 2, 3), lty = 3)
axis(1)
axis(2, at = c(1, 1.01, 1.02))
lines(lowess(x = MeanTimecourses$time_s, y = MeanTimecourses$down_neu, f = 0.05), col = col3, lwd = 2)
lines(lowess(x = MeanTimecourses$time_s, y = MeanTimecourses$up_neu, f= 0.05), col = col4, lwd = 2, lty = 5)
legend("topleft", lty = c(1, 5), col = c(col3, col4), legend = c("Downregulate", "Upregulate"), lwd = 2, bty = "n")
dev.off()

pdf("Fig_SCR_ER2.pdf", width = 7, height = 5)
plot(MeanTimecourses$down_neg ~ MeanTimecourses$time_s, type = "n", frame.plot = F, main = "Skin conductance responses, Negative stimuli", xlab = "Time (s)", ylab = "EMG (ratio)", xaxt = "n", yaxt = "n", ylim = c(0.998, 1.02))
rect(3, 0.000001, 5, 10, col = "gray88", border = NA)
abline(v = c(0, 2, 3), lty = 3)
axis(1)
axis(2, at = c(1, 1.01, 1.02))
lines(lowess(x = MeanTimecourses$time_s, y = MeanTimecourses$down_neg, f = 0.05), col = col3, lwd = 2)
lines(lowess(x = MeanTimecourses$time_s, y = MeanTimecourses$up_neg, f= 0.05), col = col4, lwd = 2, lty = 5)
legend("topleft", lty = c(1, 5), col = c(col3, col4), legend = c("Downregulate", "Upregulate"), lwd = 2, bty = "n")
dev.off()

# Compare to others' data
# TODO
#Urry2009 <- read.csv("~/BeteendeforsokEmpati2011/Data-and-analysis-code-Oxazepam-and-emotion-2/Urry2009.csv", header=FALSE, dec=",")
#Urry2009$time <- seq(0, 7.75, 0.25)
#Urry2009$corr_emg_ratio <- 1 + Urry2009$V1

#pdf("Fig_EMG_ER3.pdf", width = 7, height = 5)
#plot(Urry2009$corr_emg_ratio ~ Urry2009$time, type = "l", frame.plot= F, lwd = 2, main = "Upregulate negative", xlab = "time (s)", ylab = "EMG (ratio)", col = col5)
#abline(v = c(0, 2, 3), lty = 3)
#lines(lowess(x = MeanTimecourses$time_s, y = MeanTimecourses$up_neg, f= 0.05), col = col4, lwd = 2)
#legend("bottomright", lty = 1, lwd = 2, col = c(col5, col4), legend = c("Urry2009", "Nilsonne2015"), bty = "n")
#dev.off()


# Prepare data for statistical modelling
SCR_event_data <- merge(SCR_event_data, demData, by = "Subject")
SCR_event_data$Subject <- as.factor(SCR_event_data$Subject)
SCR_event_data$Treatment <- relevel(SCR_event_data$Treatment, ref = "Placebo")
SCR_event_data$instruction <- relevel(SCR_event_data$instruction, ref = "Downregulate")
SCR_event_data$valence <- relevel(SCR_event_data$valence, ref = "Neutral")
contrasts(SCR_event_data$instruction) <- c(-0.5, 0.5) 
contrasts(SCR_event_data$valence) <- c(-0.5, 0.5) 
contrasts(SCR_event_data$Treatment) <- c(-0.5, 0.5) 


# Build model
lme1 <- lme(log_SCR_mean_stimulus ~ instruction*valence*Treatment, data = SCR_event_data, random = ~1|Subject, na.action = na.omit)
lme1b <- lme(log_SCR_mean_stimulus ~ instruction*valence, data = SCR_event_data, random = ~1|Subject, na.action = na.omit)
lme1c <- lme(log_SCR_mean_stimulus ~ instruction*Treatment, data = SCR_event_data[SCR_event_data$valence == "Negative", ], random = ~1|Subject, na.action = na.omit)

plot(lme1)
summary(lme1)
intervals(lme1)

eff1 <- effect("instruction*valence*Treatment", lme1)

pdf("Fig_SCR_ER4.pdf", width = 7, height = 5)
plot(c(eff1$fit[1], eff1$fit[2]),
     type = "b",
     frame.plot = F,
     ylab = "SCR (log ratio)",
     xlab = "Instruction",
     xaxt = "n",
     yaxt = "n",
     xlim = c(1, 2.1),
     ylim = c(-0.01, 0.02),
     col = col1,
     main = "Skin conductance responses, neutral stimuli"
)
lines(c(1, 1), c(eff1$upper[1], eff1$lower[1]), col = col1)
lines(c(2, 2), c(eff1$upper[2], eff1$lower[2]), col = col1)
lines(c(1.1, 2.1), c(eff1$fit[5], eff1$fit[6]), type = "b", col = col2, pch = 16)
lines(c(1.1, 1.1), c(eff1$upper[5], eff1$lower[5]), col = col2)
lines(c(2.1, 2.1), c(eff1$upper[6], eff1$lower[6]), col = col2)
axis(1, at = c(1.05, 2.05), labels = c("Downregulate", "Upregulate"))
axis(2, at = c(-0.01, 0, 0.01, 0.02))
legend("topleft", col = c(col1, col2), pch = c(1, 16), legend = c("Placebo", "Oxazepam"), bty = "n", lty = 1)
dev.off()

pdf("Fig_SCR_ER5.pdf", width = 7, height = 5)
plot(c(eff1$fit[3], eff1$fit[4]),
     type = "b",
     frame.plot = F,
     ylab = "SCR (log ratio)",
     xlab = "Instruction",
     xaxt = "n",
     yaxt = "n",
     xlim = c(1, 2.1),
     ylim = c(-0.01, 0.02),
     col = col1,
     main = "Skin conductance responses, negative stimuli"
)
lines(c(1, 1), c(eff1$upper[3], eff1$lower[3]), col = col1)
lines(c(2, 2), c(eff1$upper[4], eff1$lower[4]), col = col1)
lines(c(1.1, 2.1), c(eff1$fit[7], eff1$fit[8]), type = "b", col = col2, pch = 16)
lines(c(1.1, 1.1), c(eff1$upper[7], eff1$lower[7]), col = col2)
lines(c(2.1, 2.1), c(eff1$upper[8], eff1$lower[8]), col = col2)
axis(1, at = c(1.05, 2.05), labels = c("Downregulate", "Upregulate"))
axis(2, at = c(-0.01, 0, 0.01, 0.02))
legend("topleft", col = c(col1, col2), pch = c(1, 16), legend = c("Placebo", "Oxazepam"), bty = "n", lty = 1)
dev.off()

# Compare plot to less custom-generated output for verification
plot(effect("instruction*valence*Treatment", lme1))

# Analyse other rating scales as predictors
# Since rating scales may be collinear, we include them one by one
# Since oxazepam had no significant effect, we remove it from the model

# IRI-EC
emg_event_data$IRI_EC_z <- scale(emg_event_data$IRI_EC)
lme1d <- lme(log_SCR_mean_stimulus ~ instruction*valence + IRI_EC_z*instruction + IRI_EC_z*valence, data = emg_event_data, random = ~1|Subject, na.action = na.omit)
plot(lme1d)
summary(lme1d)

# IRI-PT
emg_event_data$IRI_PT_z <- scale(emg_event_data$IRI_PT)
lme2 <- lme(log_SCR_mean_stimulus ~ instruction*valence + IRI_PT_z*instruction + IRI_PT_z*valence, data = emg_event_data, random = ~1|Subject, na.action = na.omit)
plot(lme2)
summary(lme2)

# IRI-PD
emg_event_data$IRI_PD_z <- scale(emg_event_data$IRI_PD)
lme3 <- lme(log_SCR_mean_stimulus ~ instruction*valence + IRI_PD_z*instruction + IRI_PD_z*valence, data = emg_event_data, random = ~1|Subject, na.action = na.omit)
plot(lme3)
summary(lme3)

# IRI-F
emg_event_data$IRI_F_z <- scale(emg_event_data$IRI_F)
lme4 <- lme(log_SCR_mean_stimulus ~ instruction*valence + IRI_F_z*instruction + IRI_F_z*valence, data = emg_event_data, random = ~1|Subject, na.action = na.omit)
plot(lme4)
summary(lme4)

# STAI-T
emg_event_data$STAI.T_z <- scale(emg_event_data$STAI.T)
lme5 <- lme(log_SCR_mean_stimulus ~ instruction*valence + STAI.T_z*instruction + STAI.T_z*valence, data = emg_event_data, random = ~1|Subject, na.action = na.omit)
plot(lme5)
summary(lme5)

# TAS-20
emg_event_data$TAS.20_z <- scale(emg_event_data$TAS.20)
lme6 <- lme(log_SCR_mean_stimulus ~ instruction*valence + TAS.20_z*instruction + TAS.20_z*valence, data = emg_event_data, random = ~1|Subject, na.action = na.omit)
plot(lme6)
summary(lme6)

# PPI-R-SCI
emg_event_data$PPI_SCI_z <- emg_event_data$PPI_1_SCI_R
emg_event_data$PPI_SCI_z[emg_event_data$PPI_1_IR >= 45] <- NA # Exclude participants with too high scores on Inconsistent Responding measure
emg_event_data$PPI_SCI_z <- scale(emg_event_data$PPI_SCI_z)
lme7 <- lme(log_SCR_mean_stimulus ~ instruction*valence + PPI_SCI_z*instruction + PPI_SCI_z*valence, data = emg_event_data, random = ~1|Subject, na.action = na.omit)
plot(lme7)
summary(lme7)

# PPI-R-FD
emg_event_data$PPI_FD_z <- emg_event_data$PPI_1_FD_R
emg_event_data$PPI_FD_z[emg_event_data$PPI_1_IR >= 45] <- NA # Exclude participants with too high scores on Inconsistent Responding measure
emg_event_data$PPI_FD_z <- scale(emg_event_data$PPI_FD_z)
lme8 <- lme(log_SCR_mean_stimulus ~ instruction*valence + PPI_FD_z*instruction + PPI_FD_z*valence, data = emg_event_data, random = ~1|Subject, na.action = na.omit)
plot(lme8)
summary(lme8)

# PPI-R-C
emg_event_data$PPI_C_z <- emg_event_data$PPI_1_C_R
emg_event_data$PPI_C_z[emg_event_data$PPI_1_IR >= 45] <- NA # Exclude participants with too high scores on Inconsistent Responding measure
emg_event_data$PPI_C_z <- scale(emg_event_data$PPI_C_z)
lme9 <- lme(log_SCR_mean_stimulus ~ instruction*valence + PPI_C_z*instruction + PPI_C_z*valence, data = emg_event_data, random = ~1|Subject, na.action = na.omit)
plot(lme9)
summary(lme9)

# Plot effects of rating scales
# Main effect on corrugator responses
data_corr_main <- data.frame(scale = "PPI-R-C", beta = intervals(lme9)$fixed[4, 2], lower = intervals(lme9)$fixed[4, 1], upper = intervals(lme9)$fixed[4, 3], group = "PPI", p = round(summary(lme9)$tTable[4, 5], 3))
data_corr_main <- rbind(data_corr_main, data.frame(scale = "PPI-R-FD", beta = intervals(lme8)$fixed[4, 2], lower = intervals(lme8)$fixed[4, 1], upper = intervals(lme8)$fixed[4, 3], group = "PPI", p = round(summary(lme8)$tTable[4, 5], 3)))
data_corr_main <- rbind(data_corr_main, data.frame(scale = "PPI-R-SCI", beta = intervals(lme7)$fixed[4, 2], lower = intervals(lme7)$fixed[4, 1], upper = intervals(lme7)$fixed[4, 3], group = "PPI", p = round(summary(lme7)$tTable[4, 5], 3)))
data_corr_main <- rbind(data_corr_main, data.frame(scale = "TAS-20", beta = intervals(lme6)$fixed[4, 2], lower = intervals(lme6)$fixed[4, 1], upper = intervals(lme6)$fixed[4, 3], group = "TAS", p = round(summary(lme6)$tTable[4, 5], 3)))
data_corr_main <- rbind(data_corr_main, data.frame(scale = "STAI-T", beta = intervals(lme5)$fixed[4, 2], lower = intervals(lme5)$fixed[4, 1], upper = intervals(lme5)$fixed[4, 3], group = "STAI", p = round(summary(lme5)$tTable[4, 5], 3)))
data_corr_main <- rbind(data_corr_main, data.frame(scale = "IRI-F", beta = intervals(lme4)$fixed[4, 2], lower = intervals(lme4)$fixed[4, 1], upper = intervals(lme4)$fixed[4, 3], group = "IRI", p = round(summary(lme4)$tTable[4, 5], 3)))
data_corr_main <- rbind(data_corr_main, data.frame(scale = "IRI-PD", beta = intervals(lme3)$fixed[4, 2], lower = intervals(lme3)$fixed[4, 1], upper = intervals(lme3)$fixed[4, 3], group = "IRI", p = round(summary(lme3)$tTable[4, 5], 3)))
data_corr_main <- rbind(data_corr_main, data.frame(scale = "IRI-PT", beta = intervals(lme2)$fixed[4, 2], lower = intervals(lme2)$fixed[4, 1], upper = intervals(lme2)$fixed[4, 3], group = "IRI", p = round(summary(lme2)$tTable[4, 5], 3)))
data_corr_main <- rbind(data_corr_main, data.frame(scale = "IRI-EC", beta = intervals(lme1d)$fixed[4, 2], lower = intervals(lme1d)$fixed[4, 1], upper = intervals(lme1d)$fixed[4, 3], group = "IRI", p = round(summary(lme1d)$tTable[4, 5], 3)))
data_corr_main <- data_corr_main[c(9:1), ] # Reverse order

pdf("Fig_EMG_ER6.pdf", width = 7, height = 5)
par(mar=c(5.1, 5.4, 4, 14))
plot(x = data_corr_main$beta, y = c(12:9, 7, 5, 3:1), xlab = expression(beta), ylab = "", frame.plot = F, xlim = c(min(data_corr_main$lower), max(data_corr_main$upper)), xaxt = "n", yaxt = "n")
title("Corrugator EMG, Main effects of personality measures", line = 2)
abline(v = 0, col = "gray")
axis(1, at = c(-0.1, 0, 0.1), labels = c(-0.1, 0, 0.1))
points(x = data_corr_main$beta, y = c(12:9, 7, 5, 3:1), pch = 16)
arrows(data_corr_main$lower, c(12:9, 7, 5, 3:1), data_corr_main$upper, c(12:9, 7, 5, 3:1), length = 0, lwd = 1.5)
par(las=1)
mtext(side = 2, at = c(12:9, 7, 5, 3:1), text = c("IRI-EC", " IRI-PT", "IRI-PD", "IRI-F", "STAI-T", "TAS-20", "PPI-R-SCI", "PPI-R-FD", "PPI-R-C"), line = 1)
mtext(side = 4, at = 13, text = expression(beta), line = 1)
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = round(data_corr_main$beta, 2), line = 1)
mtext(side = 4, at = 13, text = "95 % CI", line = 4)
CI <- paste("[", round(data_corr_main$lower, 2), ", ", round(data_corr_main$upper, 2), "]", sep = "")
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = CI, line = 4)
mtext(side = 4, at = 13, text = expression(italic(p)), line = 10)
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = data_corr_main$p, line = 10)
dev.off()

# Interaction with valence
data_corr_valence <- data.frame(scale = "PPI-R-C", beta = intervals(lme9)$fixed[7, 2], lower = intervals(lme9)$fixed[7, 1], upper = intervals(lme9)$fixed[7, 3], group = "PPI", p = round(summary(lme9)$tTable[7, 5], 3))
data_corr_valence <- rbind(data_corr_valence, data.frame(scale = "PPI-R-FD", beta = intervals(lme8)$fixed[7, 2], lower = intervals(lme8)$fixed[7, 1], upper = intervals(lme8)$fixed[7, 3], group = "PPI", p = round(summary(lme8)$tTable[7, 5], 3)))
data_corr_valence <- rbind(data_corr_valence, data.frame(scale = "PPI-R-SCI", beta = intervals(lme7)$fixed[7, 2], lower = intervals(lme7)$fixed[7, 1], upper = intervals(lme7)$fixed[7, 3], group = "PPI", p = round(summary(lme7)$tTable[7, 5], 3)))
data_corr_valence <- rbind(data_corr_valence, data.frame(scale = "TAS-20", beta = intervals(lme6)$fixed[7, 2], lower = intervals(lme6)$fixed[7, 1], upper = intervals(lme6)$fixed[7, 3], group = "TAS", p = round(summary(lme6)$tTable[7, 5], 3)))
data_corr_valence <- rbind(data_corr_valence, data.frame(scale = "STAI-T", beta = intervals(lme5)$fixed[7, 2], lower = intervals(lme5)$fixed[7, 1], upper = intervals(lme5)$fixed[7, 3], group = "STAI", p = round(summary(lme5)$tTable[7, 5], 3)))
data_corr_valence <- rbind(data_corr_valence, data.frame(scale = "IRI-F", beta = intervals(lme4)$fixed[7, 2], lower = intervals(lme4)$fixed[7, 1], upper = intervals(lme4)$fixed[7, 3], group = "IRI", p = round(summary(lme4)$tTable[7, 5], 3)))
data_corr_valence <- rbind(data_corr_valence, data.frame(scale = "IRI-PD", beta = intervals(lme3)$fixed[7, 2], lower = intervals(lme3)$fixed[7, 1], upper = intervals(lme3)$fixed[7, 3], group = "IRI", p = round(summary(lme3)$tTable[7, 5], 3)))
data_corr_valence <- rbind(data_corr_valence, data.frame(scale = "IRI-PT", beta = intervals(lme2)$fixed[7, 2], lower = intervals(lme2)$fixed[7, 1], upper = intervals(lme2)$fixed[7, 3], group = "IRI", p = round(summary(lme2)$tTable[7, 5], 3)))
data_corr_valence <- rbind(data_corr_valence, data.frame(scale = "IRI-EC", beta = intervals(lme1d)$fixed[7, 2], lower = intervals(lme1d)$fixed[7, 1], upper = intervals(lme1d)$fixed[7, 3], group = "IRI", p = round(summary(lme1d)$tTable[7, 5], 3)))
data_corr_valence <- data_corr_valence[c(9:1), ] # Reverse order

pdf("Fig_EMG5.pdf", width = 7, height = 5)
par(mar=c(5.1, 5.4, 4, 14))
plot(x = data_corr_valence$beta, y = c(12:9, 7, 5, 3:1), xlab = expression(beta), ylab = "", frame.plot = F, xlim = c(min(data_corr_valence$lower), max(data_corr_valence$upper)), xaxt = "n", yaxt = "n")
title("Corrugator EMG, interaction with stimulus valence", line = 2)
abline(v = 0, col = "gray")
axis(1, at = c(-0.08, 0, 0.08), labels = c(-0.08, 0, 0.08))
points(x = data_corr_valence$beta, y = c(12:9, 7, 5, 3:1), pch = 16)
arrows(data_corr_valence$lower, c(12:9, 7, 5, 3:1), data_corr_valence$upper, c(12:9, 7, 5, 3:1), length = 0, lwd = 1.5)
par(las=1)
mtext(side = 2, at = c(12:9, 7, 5, 3:1), text = c("IRI-EC", " IRI-PT", "IRI-PD", "IRI-F", "STAI-T", "TAS-20", "PPI-R-SCI", "PPI-R-FD", "PPI-R-C"), line = 1)
mtext(side = 4, at = 13, text = expression(beta), line = 1)
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = round(data_corr_valence$beta, 2), line = 1)
mtext(side = 4, at = 13, text = "95 % CI", line = 4)
CI <- paste("[", round(data_corr_valence$lower, 2), ", ", round(data_corr_valence$upper, 2), "]", sep = "")
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = CI, line = 4)
mtext(side = 4, at = 13, text = expression(italic(p)), line = 10)
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = data_corr_valence$p, line = 10)
dev.off()

# Interaction with instruction
data_corr_instruction <- data.frame(scale = "PPI-R-C", beta = intervals(lme9)$fixed[6, 2], lower = intervals(lme9)$fixed[6, 1], upper = intervals(lme9)$fixed[6, 3], group = "PPI", p = round(summary(lme9)$tTable[6, 5], 3))
data_corr_instruction <- rbind(data_corr_instruction, data.frame(scale = "PPI-R-FD", beta = intervals(lme8)$fixed[6, 2], lower = intervals(lme8)$fixed[6, 1], upper = intervals(lme8)$fixed[6, 3], group = "PPI", p = round(summary(lme8)$tTable[6, 5], 3)))
data_corr_instruction <- rbind(data_corr_instruction, data.frame(scale = "PPI-R-SCI", beta = intervals(lme7)$fixed[6, 2], lower = intervals(lme7)$fixed[6, 1], upper = intervals(lme7)$fixed[6, 3], group = "PPI", p = round(summary(lme7)$tTable[6, 5], 3)))
data_corr_instruction <- rbind(data_corr_instruction, data.frame(scale = "TAS-20", beta = intervals(lme6)$fixed[6, 2], lower = intervals(lme6)$fixed[6, 1], upper = intervals(lme6)$fixed[6, 3], group = "TAS", p = round(summary(lme6)$tTable[6, 5], 3)))
data_corr_instruction <- rbind(data_corr_instruction, data.frame(scale = "STAI-T", beta = intervals(lme5)$fixed[6, 2], lower = intervals(lme5)$fixed[6, 1], upper = intervals(lme5)$fixed[6, 3], group = "STAI", p = round(summary(lme5)$tTable[6, 5], 3)))
data_corr_instruction <- rbind(data_corr_instruction, data.frame(scale = "IRI-F", beta = intervals(lme4)$fixed[6, 2], lower = intervals(lme4)$fixed[6, 1], upper = intervals(lme4)$fixed[6, 3], group = "IRI", p = round(summary(lme4)$tTable[6, 5], 3)))
data_corr_instruction <- rbind(data_corr_instruction, data.frame(scale = "IRI-PD", beta = intervals(lme3)$fixed[6, 2], lower = intervals(lme3)$fixed[6, 1], upper = intervals(lme3)$fixed[6, 3], group = "IRI", p = round(summary(lme3)$tTable[6, 5], 3)))
data_corr_instruction <- rbind(data_corr_instruction, data.frame(scale = "IRI-PT", beta = intervals(lme2)$fixed[6, 2], lower = intervals(lme2)$fixed[6, 1], upper = intervals(lme2)$fixed[6, 3], group = "IRI", p = round(summary(lme2)$tTable[6, 5], 3)))
data_corr_instruction <- rbind(data_corr_instruction, data.frame(scale = "IRI-EC", beta = intervals(lme1d)$fixed[6, 2], lower = intervals(lme1d)$fixed[6, 1], upper = intervals(lme1d)$fixed[6, 3], group = "IRI", p = round(summary(lme1d)$tTable[6, 5], 3)))
data_corr_instruction <- data_corr_instruction[c(9:1), ] # Reverse order

pdf("Fig_EMG5.pdf", width = 7, height = 5)
par(mar=c(5.1, 5.4, 4, 14))
plot(x = data_corr_instruction$beta, y = c(12:9, 7, 5, 3:1), xlab = expression(beta), ylab = "", frame.plot = F, xlim = c(min(data_corr_instruction$lower), max(data_corr_instruction$upper)), xaxt = "n", yaxt = "n")
title("Corrugator EMG, interaction with instruction", line = 2)
abline(v = 0, col = "gray")
axis(1, at = c(-0.1, 0, 0.1), labels = c(-0.1, 0, 0.1))
points(x = data_corr_instruction$beta, y = c(12:9, 7, 5, 3:1), pch = 16)
arrows(data_corr_instruction$lower, c(12:9, 7, 5, 3:1), data_corr_instruction$upper, c(12:9, 7, 5, 3:1), length = 0, lwd = 1.5)
par(las=1)
mtext(side = 2, at = c(12:9, 7, 5, 3:1), text = c("IRI-EC", " IRI-PT", "IRI-PD", "IRI-F", "STAI-T", "TAS-20", "PPI-R-SCI", "PPI-R-FD", "PPI-R-C"), line = 1)
mtext(side = 4, at = 13, text = expression(beta), line = 1)
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = round(data_corr_instruction$beta, 2), line = 1)
mtext(side = 4, at = 13, text = "95 % CI", line = 4)
CI <- paste("[", round(data_corr_instruction$lower, 2), ", ", round(data_corr_instruction$upper, 2), "]", sep = "")
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = CI, line = 4)
mtext(side = 4, at = 13, text = expression(italic(p)), line = 10)
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = data_corr_instruction$p, line = 10)
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