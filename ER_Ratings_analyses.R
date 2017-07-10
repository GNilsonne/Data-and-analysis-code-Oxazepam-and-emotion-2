# Script to analyse ratings in the ER experiment of the Oxazepam project 2011
# Gustav Nilsonne 2015-07-16

# Require packages
require(RCurl) # To read data from GitHub
require(RColorBrewer)
require(nlme)
require(effects)

# Define colors for later
col1 = brewer.pal(8, "Dark2")[1]
col2 = brewer.pal(8, "Dark2")[2]
col3 = brewer.pal(8, "Dark2")[3]
col4 = brewer.pal(8, "Dark2")[4]
col5 = brewer.pal(8, "Dark2")[5]
col6 = brewer.pal(8, "Dark2")[6]
col7 = brewer.pal(8, "Dark2")[7]
col8 = brewer.pal(8, "Dark2")[8]

# Read data
RatingsData <- read.csv("RatingsData_ER.csv")

demDataURL <- getURL("https://raw.githubusercontent.com/GNilsonne/Data-and-analysis-code-Oxazepam-and-emotion/master/demographics.csv", ssl.verifypeer = FALSE)
demData <- read.csv(text = demDataURL)

data <- merge(RatingsData, demData, by = "Subject")

# Set up contrasts
data$Subject <- as.factor(data$Subject)
data$Treatment <- as.factor(data$Treatment)
data$Treatment <- relevel(data$Treatment, ref = "Placebo")
data$instruction <- relevel(data$instruction, ref = "Downregulate")
data$valence <- relevel(data$valence, ref = "Neutral")
contrasts(data$instruction) <- c(-0.5, 0.5) 
contrasts(data$valence) <- c(-0.5, 0.5) 
contrasts(data$Treatment) <- c(-0.5, 0.5) 

# Plot data
boxplot(RatedUnpleasantness ~ valence*instruction*Treatment, data = data, border = c(col3, col3, col4, col4, col3, col3, col4, col4), names = c("Neu Pl", "Neg Pl", "Neu Pl", "Neg Pl", "Neu Ox", "Neg Ox", "Neu Ox", "Neg Ox"))
legend("topleft", lty = 1, col = c(col3, col4), legend = c("Downregulate", "Upregulate"))

# Analyse data
lme1 <- lme(RatedUnpleasantness ~ instruction*valence*Treatment, data = data, random = ~1|Subject, na.action = na.omit)
plot(lme1)
summary(lme1)
intervals(lme1)

eff1 <- effect("instruction*valence*Treatment", lme1)

pdf("Fig_Ratings_ER1.pdf", width = 7, height = 5)
plot(c(eff1$fit[1], eff1$fit[2]),
     type = "b",
     frame.plot = F,
     ylab = "VAS rating",
     xlab = "Instruction",
     xaxt = "n",
     yaxt = "n",
     xlim = c(1, 2.1),
     ylim = c(0, 60),
     col = col1,
     main = "Rated unpleasantness, neutral stimuli"
)
lines(c(1, 1), c(eff1$upper[1], eff1$lower[1]), col = col1)
lines(c(2, 2), c(eff1$upper[2], eff1$lower[2]), col = col1)
lines(c(1.1, 2.1), c(eff1$fit[5], eff1$fit[6]), type = "b", col = col2, pch = 16)
lines(c(1.1, 1.1), c(eff1$upper[5], eff1$lower[5]), col = col2)
lines(c(2.1, 2.1), c(eff1$upper[6], eff1$lower[6]), col = col2)
axis(1, at = c(1.05, 2.05), labels = c("Downregulate", "Upregulate"))
axis(2, at = c(0, 20, 40, 60))
legend("topleft", col = c(col1, col2), pch = c(1, 16), legend = c("Placebo", "Oxazepam"), bty = "n", lty = 1)
dev.off()

pdf("Fig_Ratings_ER2.pdf", width = 7, height = 5)
plot(c(eff1$fit[3], eff1$fit[4]),
     type = "b",
     frame.plot = F,
     ylab = "VAS rating",
     xlab = "Instruction",
     xaxt = "n",
     yaxt = "n",
     xlim = c(1, 2.1),
     ylim = c(0, 60),
     col = col1,
     main = "Rated unpleasantness, negative stimuli"
)
lines(c(1, 1), c(eff1$upper[3], eff1$lower[3]), col = col1)
lines(c(2, 2), c(eff1$upper[4], eff1$lower[4]), col = col1)
lines(c(1.1, 2.1), c(eff1$fit[7], eff1$fit[8]), type = "b", col = col2, pch = 16)
lines(c(1.1, 1.1), c(eff1$upper[7], eff1$lower[7]), col = col2)
lines(c(2.1, 2.1), c(eff1$upper[8], eff1$lower[8]), col = col2)
axis(1, at = c(1.05, 2.05), labels = c("Downregulate", "Upregulate"))
axis(2, at = c(0, 20, 40, 60))
legend("topleft", col = c(col1, col2), pch = c(1, 16), legend = c("Placebo", "Oxazepam"), bty = "n", lty = 1)
dev.off()

# Compare plots to less custom-generated output for verification
plot(effect("instruction*valence*Treatment", lme1))

# Analyse other rating scales as predictors
# Since rating scales may be collinear, we include them one by one
# Since oxazepam had no significant effect, we remove it from the model

# IRI-EC
data$IRI_EC_z <- scale(data$IRI_EC)
lme1d <- lme(RatedUnpleasantness ~ instruction*valence*Treatment + IRI_EC_z*instruction + IRI_EC_z*valence, data = data, random = ~1|Subject, na.action = na.omit)
plot(lme1d)
summary(lme1d)

# IRI-PT
data$IRI_PT_z <- scale(data$IRI_PT)
lme2 <- lme(RatedUnpleasantness ~ instruction*valence*Treatment + IRI_PT_z*instruction + IRI_PT_z*valence, data = data, random = ~1|Subject, na.action = na.omit)
plot(lme2)
summary(lme2)

# IRI-PD
data$IRI_PD_z <- scale(data$IRI_PD)
lme3 <- lme(RatedUnpleasantness ~ instruction*valence*Treatment + IRI_PD_z*instruction + IRI_PD_z*valence, data = data, random = ~1|Subject, na.action = na.omit)
plot(lme3)
summary(lme3)

# IRI-F
data$IRI_F_z <- scale(data$IRI_F)
lme4 <- lme(RatedUnpleasantness ~ instruction*valence*Treatment + IRI_F_z*instruction + IRI_F_z*valence, data = data, random = ~1|Subject, na.action = na.omit)
plot(lme4)
summary(lme4)

# STAI-T
data$STAI.T_z <- scale(data$STAI.T)
lme5 <- lme(RatedUnpleasantness ~ instruction*valence*Treatment + STAI.T_z*instruction + STAI.T_z*valence, data = data, random = ~1|Subject, na.action = na.omit)
plot(lme5)
summary(lme5)

# TAS-20
data$TAS.20_z <- scale(data$TAS.20)
lme6 <- lme(RatedUnpleasantness ~ instruction*valence*Treatment + TAS.20_z*instruction + TAS.20_z*valence, data = data, random = ~1|Subject, na.action = na.omit)
plot(lme6)
summary(lme6)

# PPI-R-SCI
data$PPI_SCI_z <- data$PPI_1_SCI_R
data$PPI_SCI_z[data$PPI_1_IR >= 45] <- NA # Exclude participants with too high scores on Inconsistent Responding measure
data$PPI_SCI_z <- scale(data$PPI_SCI_z)
lme7 <- lme(RatedUnpleasantness ~ instruction*valence*Treatment + PPI_SCI_z*instruction + PPI_SCI_z*valence, data = data, random = ~1|Subject, na.action = na.omit)
plot(lme7)
summary(lme7)

# PPI-R-FD
data$PPI_FD_z <- data$PPI_1_FD_R
data$PPI_FD_z[data$PPI_1_IR >= 45] <- NA # Exclude participants with too high scores on Inconsistent Responding measure
data$PPI_FD_z <- scale(data$PPI_FD_z)
lme8 <- lme(RatedUnpleasantness ~ instruction*valence*Treatment + PPI_FD_z*instruction + PPI_FD_z*valence, data = data, random = ~1|Subject, na.action = na.omit)
plot(lme8)
summary(lme8)

# PPI-R-C
data$PPI_C_z <- data$PPI_1_C_R
data$PPI_C_z[data$PPI_1_IR >= 45] <- NA # Exclude participants with too high scores on Inconsistent Responding measure
data$PPI_C_z <- scale(data$PPI_C_z)
lme9 <- lme(RatedUnpleasantness ~ instruction*valence*Treatment + PPI_C_z*instruction + PPI_C_z*valence, data = data, random = ~1|Subject, na.action = na.omit)
plot(lme9)
summary(lme9)

# Plot effects of rating scales
# Main effect on ratings
data_main <- data.frame(scale = "PPI-R-C", beta = intervals(lme9)$fixed[4, 2], lower = intervals(lme9)$fixed[4, 1], upper = intervals(lme9)$fixed[4, 3], group = "PPI", p = round(summary(lme9)$tTable[4, 5], 3))
data_main <- rbind(data_main, data.frame(scale = "PPI-R-FD", beta = intervals(lme8)$fixed[4, 2], lower = intervals(lme8)$fixed[4, 1], upper = intervals(lme8)$fixed[4, 3], group = "PPI", p = round(summary(lme8)$tTable[4, 5], 3)))
data_main <- rbind(data_main, data.frame(scale = "PPI-R-SCI", beta = intervals(lme7)$fixed[4, 2], lower = intervals(lme7)$fixed[4, 1], upper = intervals(lme7)$fixed[4, 3], group = "PPI", p = round(summary(lme7)$tTable[4, 5], 3)))
data_main <- rbind(data_main, data.frame(scale = "TAS-20", beta = intervals(lme6)$fixed[4, 2], lower = intervals(lme6)$fixed[4, 1], upper = intervals(lme6)$fixed[4, 3], group = "TAS", p = round(summary(lme6)$tTable[4, 5], 3)))
data_main <- rbind(data_main, data.frame(scale = "STAI-T", beta = intervals(lme5)$fixed[4, 2], lower = intervals(lme5)$fixed[4, 1], upper = intervals(lme5)$fixed[4, 3], group = "STAI", p = round(summary(lme5)$tTable[4, 5], 3)))
data_main <- rbind(data_main, data.frame(scale = "IRI-F", beta = intervals(lme4)$fixed[4, 2], lower = intervals(lme4)$fixed[4, 1], upper = intervals(lme4)$fixed[4, 3], group = "IRI", p = round(summary(lme4)$tTable[4, 5], 3)))
data_main <- rbind(data_main, data.frame(scale = "IRI-PD", beta = intervals(lme3)$fixed[4, 2], lower = intervals(lme3)$fixed[4, 1], upper = intervals(lme3)$fixed[4, 3], group = "IRI", p = round(summary(lme3)$tTable[4, 5], 3)))
data_main <- rbind(data_main, data.frame(scale = "IRI-PT", beta = intervals(lme2)$fixed[4, 2], lower = intervals(lme2)$fixed[4, 1], upper = intervals(lme2)$fixed[4, 3], group = "IRI", p = round(summary(lme2)$tTable[4, 5], 3)))
data_main <- rbind(data_main, data.frame(scale = "IRI-EC", beta = intervals(lme1d)$fixed[4, 2], lower = intervals(lme1d)$fixed[4, 1], upper = intervals(lme1d)$fixed[4, 3], group = "IRI", p = round(summary(lme1d)$tTable[4, 5], 3)))
data_main <- data_main[c(9:1), ] # Reverse order

pdf("Fig_Ratings_ER3.pdf", width = 7, height = 5)
par(mar=c(5.1, 5.4, 4, 14))
plot(x = data_main$beta, y = c(12:9, 7, 5, 3:1), xlab = expression(beta), ylab = "", frame.plot = F, xlim = c(min(data_main$lower), max(data_main$upper)), xaxt = "n", yaxt = "n")
title("Rated Unpleasantness, Main effects of personality measures", line = 2)
abline(v = 0, col = "gray")
axis(1, at = c(-10, 0, 10), labels = c(-10, 0, 10))
points(x = data_main$beta, y = c(12:9, 7, 5, 3:1), pch = 16)
arrows(data_main$lower, c(12:9, 7, 5, 3:1), data_main$upper, c(12:9, 7, 5, 3:1), length = 0, lwd = 1.5)
par(las=1)
mtext(side = 2, at = c(12:9, 7, 5, 3:1), text = c("IRI-EC", " IRI-PT", "IRI-PD", "IRI-F", "STAI-T", "TAS-20", "PPI-R-SCI", "PPI-R-FD", "PPI-R-C"), line = 1)
mtext(side = 4, at = 13, text = expression(beta), line = 1)
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = round(data_main$beta, 2), line = 1)
mtext(side = 4, at = 13, text = "95 % CI", line = 4)
CI <- paste("[", round(data_main$lower, 2), ", ", round(data_main$upper, 2), "]", sep = "")
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = CI, line = 4)
mtext(side = 4, at = 13, text = expression(italic(p)), line = 10)
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = data_main$p, line = 10)
dev.off()

# Interaction with valence
data_valence <- data.frame(scale = "PPI-R-C", beta = intervals(lme9)$fixed[10, 2], lower = intervals(lme9)$fixed[10, 1], upper = intervals(lme9)$fixed[10, 3], group = "PPI", p = round(summary(lme9)$tTable[10, 5], 3))
data_valence <- rbind(data_valence, data.frame(scale = "PPI-R-FD", beta = intervals(lme8)$fixed[10, 2], lower = intervals(lme8)$fixed[10, 1], upper = intervals(lme8)$fixed[10, 3], group = "PPI", p = round(summary(lme8)$tTable[10, 5], 3)))
data_valence <- rbind(data_valence, data.frame(scale = "PPI-R-SCI", beta = intervals(lme7)$fixed[10, 2], lower = intervals(lme7)$fixed[10, 1], upper = intervals(lme7)$fixed[10, 3], group = "PPI", p = round(summary(lme7)$tTable[10, 5], 3)))
data_valence <- rbind(data_valence, data.frame(scale = "TAS-20", beta = intervals(lme6)$fixed[10, 2], lower = intervals(lme6)$fixed[10, 1], upper = intervals(lme6)$fixed[10, 3], group = "TAS", p = round(summary(lme6)$tTable[10, 5], 3)))
data_valence <- rbind(data_valence, data.frame(scale = "STAI-T", beta = intervals(lme5)$fixed[10, 2], lower = intervals(lme5)$fixed[10, 1], upper = intervals(lme5)$fixed[10, 3], group = "STAI", p = round(summary(lme5)$tTable[10, 5], 3)))
data_valence <- rbind(data_valence, data.frame(scale = "IRI-F", beta = intervals(lme4)$fixed[10, 2], lower = intervals(lme4)$fixed[10, 1], upper = intervals(lme4)$fixed[10, 3], group = "IRI", p = round(summary(lme4)$tTable[10, 5], 3)))
data_valence <- rbind(data_valence, data.frame(scale = "IRI-PD", beta = intervals(lme3)$fixed[10, 2], lower = intervals(lme3)$fixed[10, 1], upper = intervals(lme3)$fixed[10, 3], group = "IRI", p = round(summary(lme3)$tTable[10, 5], 3)))
data_valence <- rbind(data_valence, data.frame(scale = "IRI-PT", beta = intervals(lme2)$fixed[10, 2], lower = intervals(lme2)$fixed[10, 1], upper = intervals(lme2)$fixed[10, 3], group = "IRI", p = round(summary(lme2)$tTable[10, 5], 3)))
data_valence <- rbind(data_valence, data.frame(scale = "IRI-EC", beta = intervals(lme1d)$fixed[10, 2], lower = intervals(lme1d)$fixed[10, 1], upper = intervals(lme1d)$fixed[10, 3], group = "IRI", p = round(summary(lme1d)$tTable[10, 5], 3)))
data_valence <- data_valence[c(9:1), ] # Reverse order

pdf("Fig_Ratings_ER4.pdf", width = 7, height = 5)
par(mar=c(5.1, 5.4, 4, 14))
plot(x = data_valence$beta, y = c(12:9, 7, 5, 3:1), xlab = expression(beta), ylab = "", frame.plot = F, xlim = c(min(data_valence$lower), max(data_valence$upper)), xaxt = "n", yaxt = "n")
title("Rated unpleasantness, interaction with stimulus valence", line = 2)
abline(v = 0, col = "gray")
axis(1, at = c(-6, 0, 6), labels = c(-6, 0, 6))
points(x = data_valence$beta, y = c(12:9, 7, 5, 3:1), pch = 16)
arrows(data_valence$lower, c(12:9, 7, 5, 3:1), data_valence$upper, c(12:9, 7, 5, 3:1), length = 0, lwd = 1.5)
par(las=1)
mtext(side = 2, at = c(12:9, 7, 5, 3:1), text = c("IRI-EC", " IRI-PT", "IRI-PD", "IRI-F", "STAI-T", "TAS-20", "PPI-R-SCI", "PPI-R-FD", "PPI-R-C"), line = 1)
mtext(side = 4, at = 13, text = expression(beta), line = 1)
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = round(data_valence$beta, 2), line = 1)
mtext(side = 4, at = 13, text = "95 % CI", line = 4)
CI <- paste("[", round(data_valence$lower, 2), ", ", round(data_valence$upper, 2), "]", sep = "")
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = CI, line = 4)
mtext(side = 4, at = 13, text = expression(italic(p)), line = 10)
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = data_valence$p, line = 10)
dev.off()

# Interaction with instruction
data_instruction <- data.frame(scale = "PPI-R-C", beta = intervals(lme9)$fixed[9, 2], lower = intervals(lme9)$fixed[9, 1], upper = intervals(lme9)$fixed[9, 3], group = "PPI", p = round(summary(lme9)$tTable[9, 5], 3))
data_instruction <- rbind(data_instruction, data.frame(scale = "PPI-R-FD", beta = intervals(lme8)$fixed[9, 2], lower = intervals(lme8)$fixed[9, 1], upper = intervals(lme8)$fixed[9, 3], group = "PPI", p = round(summary(lme8)$tTable[9, 5], 3)))
data_instruction <- rbind(data_instruction, data.frame(scale = "PPI-R-SCI", beta = intervals(lme7)$fixed[9, 2], lower = intervals(lme7)$fixed[9, 1], upper = intervals(lme7)$fixed[9, 3], group = "PPI", p = round(summary(lme7)$tTable[9, 5], 3)))
data_instruction <- rbind(data_instruction, data.frame(scale = "TAS-20", beta = intervals(lme6)$fixed[9, 2], lower = intervals(lme6)$fixed[9, 1], upper = intervals(lme6)$fixed[9, 3], group = "TAS", p = round(summary(lme6)$tTable[9, 5], 3)))
data_instruction <- rbind(data_instruction, data.frame(scale = "STAI-T", beta = intervals(lme5)$fixed[9, 2], lower = intervals(lme5)$fixed[9, 1], upper = intervals(lme5)$fixed[9, 3], group = "STAI", p = round(summary(lme5)$tTable[9, 5], 3)))
data_instruction <- rbind(data_instruction, data.frame(scale = "IRI-F", beta = intervals(lme4)$fixed[9, 2], lower = intervals(lme4)$fixed[9, 1], upper = intervals(lme4)$fixed[9, 3], group = "IRI", p = round(summary(lme4)$tTable[9, 5], 3)))
data_instruction <- rbind(data_instruction, data.frame(scale = "IRI-PD", beta = intervals(lme3)$fixed[9, 2], lower = intervals(lme3)$fixed[9, 1], upper = intervals(lme3)$fixed[9, 3], group = "IRI", p = round(summary(lme3)$tTable[9, 5], 3)))
data_instruction <- rbind(data_instruction, data.frame(scale = "IRI-PT", beta = intervals(lme2)$fixed[9, 2], lower = intervals(lme2)$fixed[9, 1], upper = intervals(lme2)$fixed[9, 3], group = "IRI", p = round(summary(lme2)$tTable[9, 5], 3)))
data_instruction <- rbind(data_instruction, data.frame(scale = "IRI-EC", beta = intervals(lme1d)$fixed[9, 2], lower = intervals(lme1d)$fixed[9, 1], upper = intervals(lme1d)$fixed[9, 3], group = "IRI", p = round(summary(lme1d)$tTable[9, 5], 3)))
data_instruction <- data_instruction[c(9:1), ] # Reverse order

pdf("Fig_Ratings_ER5.pdf", width = 7, height = 5)
par(mar=c(5.1, 5.4, 4, 14))
plot(x = data_instruction$beta, y = c(12:9, 7, 5, 3:1), xlab = expression(beta), ylab = "", frame.plot = F, xlim = c(min(data_instruction$lower), max(data_instruction$upper)), xaxt = "n", yaxt = "n")
title("Rated unpleasantness, interaction with instruction", line = 2)
abline(v = 0, col = "gray")
axis(1, at = c(-2, 0, 2), labels = c(-2, 0, 2))
points(x = data_instruction$beta, y = c(12:9, 7, 5, 3:1), pch = 16)
arrows(data_instruction$lower, c(12:9, 7, 5, 3:1), data_instruction$upper, c(12:9, 7, 5, 3:1), length = 0, lwd = 1.5)
par(las=1)
mtext(side = 2, at = c(12:9, 7, 5, 3:1), text = c("IRI-EC", " IRI-PT", "IRI-PD", "IRI-F", "STAI-T", "TAS-20", "PPI-R-SCI", "PPI-R-FD", "PPI-R-C"), line = 1)
mtext(side = 4, at = 13, text = expression(beta), line = 1)
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = round(data_instruction$beta, 2), line = 1)
mtext(side = 4, at = 13, text = "95 % CI", line = 4)
CI <- paste("[", round(data_instruction$lower, 2), ", ", round(data_instruction$upper, 2), "]", sep = "")
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = CI, line = 4)
mtext(side = 4, at = 13, text = expression(italic(p)), line = 10)
mtext(side = 4, at = c(12:9, 7, 5, 3:1), text = data_instruction$p, line = 10)
dev.off()

# Write regression output tables
write.csv(summary(lme1)$tTable, file = "Ratings_main.csv")
write.csv(summary(lme1d)$tTable, file = "Ratings_IRI_EC.csv")
write.csv(summary(lme2)$tTable, file = "Ratings_IRI_PT.csv")
write.csv(summary(lme3)$tTable, file = "Ratings_IRI_PD.csv")
write.csv(summary(lme4)$tTable, file = "Ratings_IRI_F.csv")
write.csv(summary(lme5)$tTable, file = "Ratings_STAI_T.csv")
write.csv(summary(lme6)$tTable, file = "Ratings_TAS_20.csv")
write.csv(summary(lme7)$tTable, file = "Ratings_PPI_R_SCI.csv")
write.csv(summary(lme8)$tTable, file = "Ratings_PPI_R_FD.csv")
write.csv(summary(lme9)$tTable, file = "Ratings_PPI_R_C.csv")
