# Script to analyse response times from the oxazepam and emotion project
# Gustav Nilsonne 2015-08-07

# Require packages
library(RCurl) # To read data from GitHub
library(nlme) # To build mixed-effects models
library(effects) # To get confidence intervals on estimates
library(RColorBrewer) # To get good diverging colors for graphs

# Define colors for later
col1 = brewer.pal(3, "Dark2")[1]
col2 = brewer.pal(3, "Dark2")[2]

# Read data
RTDataURL <- getURL("https://raw.githubusercontent.com/GNilsonne/Data-and-analysis-code-Oxazepam-and-emotion/master/Reaction_time_data.csv", ssl.verifypeer = FALSE)
RTData <- read.csv(text = RTDataURL)

demDataURL <- getURL("https://raw.githubusercontent.com/GNilsonne/Data-and-analysis-code-Oxazepam-and-emotion/master/demographics.csv", ssl.verifypeer = FALSE)
demData <- read.csv(text = demDataURL)

# Inverse transform
RTData$Inv_RT <- 1/RTData$ResponseTime_ms

# Keep only data included in experiment
RTData <- merge(RTData, demData[, c("Subject", "Included_FMOV")], by = "Subject")
RTData <- RTData[RTData$Included_FMOV == TRUE, ]

# Verify normality of distribution
plot(density(RTData$ResponseTime_ms))
plot(density(RTData$Inv_RT))

# Test effects
RTData$FirstOrSecondTest <- as.factor(RTData$FirstOrSecondTest)
lme1 <- lme(Inv_RT ~ Treatment * FirstOrSecondTest + Wave, data = RTData, random = ~1|Subject)
summary(lme1)

# Inspect residuals
plot(lme1)

# Get estimates
intervals(lme1)
1/intervals(lme1)$fixed
1/intervals(lme1)$fixed[1, 2] - 1/(abs(intervals(lme1)$fixed) + abs(intervals(lme1)$fixed[1, 2]))

# Plot effects with transformation back to original scale
eff1 <- effect("Treatment * FirstOrSecondTest", lme1)

pdf("Fig_RT.pdf", height = 4, width = 4*1.61)
par(mar = c(4,4,1,1))
plot(1/eff1$fit[c(2, 4)],
     frame.plot = F,
     xaxt = "n",
     yaxt = "n",
     type = "b",
     xlab = "",
     ylab = "ms",
     xlim = c(1, 2.1),
     ylim = c(297, 347),
     pch = 1,
     #main = "A. Reaction times",
     col = col1)
lines(c(1.1,2.1), 1/eff1$fit[c(1, 3)], type = "b", col = col2, pch = 16)
lines(c(1, 1), c((1/eff1$upper[2]), (1/eff1$lower[2])), col = col1)
lines(c(2, 2), c((1/eff1$upper[4]), (1/eff1$lower[4])), col = col1)
lines(c(1.1, 1.1), c((1/eff1$upper[1]), (1/eff1$lower[1])), col = col2)
lines(c(2.1, 2.1), c((1/eff1$upper[3]), (1/eff1$lower[3])), col = col2)
axis(1, labels = c("Before", "After"), at = c(1.05, 2.05))
axis(2, at = c(300, 320, 340))
legend("topleft", col = c(col1, col2), pch = c(1, 16), legend = c("Placebo", "Oxazepam"), bty = "n", lwd = 1)
dev.off()
