# Script to analyse demographic data and rating scales from the oxazepam and emotion project
# Gustav Nilsonne 2015-08-07

# Require packages
library(RCurl) # To read data from GitHub
library(nlme) # To build mixed-effects models
library(effects) # To get confidence intervals on estimates
library(RColorBrewer) # To get good diverging colors for graphs

# Define colors for later
col1 = brewer.pal(3, "Dark2")[1]
col2 = brewer.pal(3, "Dark2")[2]

add.alpha <- function(col, alpha=1){ ## Function to add an alpha value to a colour, from: http://www.magesblog.com/2013/04/how-to-change-alpha-value-of-colours-in.html
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

col3 <- add.alpha(col1, alpha = 0.2)

# Read data
demDataURL <- getURL("https://raw.githubusercontent.com/GNilsonne/Data-and-analysis-code-Oxazepam-and-emotion/master/demographics.csv", ssl.verifypeer = FALSE)
demData <- read.csv(text = demDataURL)

# Descriptive analyses, n per group
length(demData$Subject[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"])
length(demData$Subject[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"])
length(demData$Subject[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"])
length(demData$Subject[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"])

# Descriptive analyses, IRI
mean(demData$IRI_EC[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
sd(demData$IRI_EC[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
mean(demData$IRI_EC[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
sd(demData$IRI_EC[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
mean(demData$IRI_EC[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
sd(demData$IRI_EC[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
mean(demData$IRI_EC[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
sd(demData$IRI_EC[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)

mean(demData$IRI_PT[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
sd(demData$IRI_PT[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
mean(demData$IRI_PT[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
sd(demData$IRI_PT[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
mean(demData$IRI_PT[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
sd(demData$IRI_PT[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
mean(demData$IRI_PT[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
sd(demData$IRI_PT[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)

mean(demData$IRI_PD[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
sd(demData$IRI_PD[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
mean(demData$IRI_PD[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
sd(demData$IRI_PD[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
mean(demData$IRI_PD[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
sd(demData$IRI_PD[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
mean(demData$IRI_PD[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
sd(demData$IRI_PD[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)

mean(demData$IRI_F[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
sd(demData$IRI_F[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
mean(demData$IRI_F[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
sd(demData$IRI_F[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
mean(demData$IRI_F[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
sd(demData$IRI_F[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
mean(demData$IRI_F[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
sd(demData$IRI_F[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)

# Descriptive analyses, TAS-20
mean(demData$TAS.20[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
sd(demData$TAS.20[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
mean(demData$TAS.20[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
sd(demData$TAS.20[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
mean(demData$TAS.20[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
sd(demData$TAS.20[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
mean(demData$TAS.20[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
sd(demData$TAS.20[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)

mean(demData$Difficulty.identifying.feelings[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
sd(demData$Difficulty.identifying.feelings[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
mean(demData$Difficulty.identifying.feelings[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
sd(demData$Difficulty.identifying.feelings[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
mean(demData$Difficulty.identifying.feelings[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
sd(demData$Difficulty.identifying.feelings[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
mean(demData$Difficulty.identifying.feelings[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
sd(demData$Difficulty.identifying.feelings[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)

mean(demData$Difficulty.describing.feelings[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
sd(demData$Difficulty.describing.feelings[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
mean(demData$Difficulty.describing.feelings[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
sd(demData$Difficulty.describing.feelings[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
mean(demData$Difficulty.describing.feelings[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
sd(demData$Difficulty.describing.feelings[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
mean(demData$Difficulty.describing.feelings[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
sd(demData$Difficulty.describing.feelings[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)

mean(demData$Externally.oriented.thinking[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
sd(demData$Externally.oriented.thinking[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
mean(demData$Externally.oriented.thinking[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
sd(demData$Externally.oriented.thinking[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
mean(demData$Externally.oriented.thinking[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
sd(demData$Externally.oriented.thinking[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
mean(demData$Externally.oriented.thinking[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
sd(demData$Externally.oriented.thinking[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)

# Descriptive analyses, STAI-T
mean(demData$STAI.T[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
sd(demData$STAI.T[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
mean(demData$STAI.T[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
sd(demData$STAI.T[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
mean(demData$STAI.T[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
sd(demData$STAI.T[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
mean(demData$STAI.T[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
sd(demData$STAI.T[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)

# Descriptive analyses, PPI-R
mean(demData$PPI_1_SCI_R[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
sd(demData$PPI_1_SCI_R[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
mean(demData$PPI_1_SCI_R[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
sd(demData$PPI_1_SCI_R[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
mean(demData$PPI_1_SCI_R[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
sd(demData$PPI_1_SCI_R[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
mean(demData$PPI_1_SCI_R[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
sd(demData$PPI_1_SCI_R[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)

mean(demData$PPI_1_FD_R[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
sd(demData$PPI_1_FD_R[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
mean(demData$PPI_1_FD_R[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
sd(demData$PPI_1_FD_R[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
mean(demData$PPI_1_FD_R[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
sd(demData$PPI_1_FD_R[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
mean(demData$PPI_1_FD_R[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
sd(demData$PPI_1_FD_R[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)

mean(demData$PPI_1_C_R[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
sd(demData$PPI_1_C_R[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
mean(demData$PPI_1_C_R[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
sd(demData$PPI_1_C_R[demData$Wave == 1 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
mean(demData$PPI_1_C_R[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
sd(demData$PPI_1_C_R[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Placebo"], na.rm = T)
mean(demData$PPI_1_C_R[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)
sd(demData$PPI_1_C_R[demData$Wave == 2 & demData$Included_FMOV == 1 & demData$Treatment == "Oxazepam"], na.rm = T)

# IRI, test-retest
demData$IRIdiff <- demData$IRI_retest_EC - demData$IRI_EC
mean(demData$IRIdiff[demData$Included_FMOV == TRUE], na.rm = T)
sd(demData$IRIdiff[demData$Included_FMOV == TRUE], na.rm = T)

demData$IRIdiff2 <- demData$IRI_scrambled_EC - demData$IRI_EC
t.test(IRIdiff2 ~ Treatment, data = demData[demData$Included_FMOV == TRUE, ])

# Analyse effect of oxazepam on rated state anxiety
# Make dataframe for mixed-effects model
STAISData <- rbind(demData[, c("Subject", "Treatment", "Wave", "Included_FMOV", "STAI.S", "STAI.S.Scrambled")], demData[, c("Subject", "Treatment", "Wave", "Included_FMOV", "STAI.S", "STAI.S.Scrambled")]) 
STAISData <- STAISData[STAISData$Included_FMOV == T, ] # Remove participants not included in this experiment
STAISData$FirstOrSecond <- c(rep.int(1, 0.5*length(STAISData$Subject)), rep.int(2, 0.5*length(STAISData$Subject)))
STAISData$STAIS <- NA # Make new column for STAI-S rating, then fill it with values for the first and second ratings, respectively
STAISData$STAIS[STAISData$FirstOrSecond == 1] <- STAISData$STAI.S[STAISData$FirstOrSecond == 1]
STAISData$STAIS[STAISData$FirstOrSecond == 2] <- STAISData$STAI.S.Scrambled[STAISData$FirstOrSecond == 2]
STAISData$Treatment <- relevel(STAISData$Treatment, ref = "Placebo")

lme1 <- lme(STAIS ~ Treatment * FirstOrSecond + Wave, data = STAISData, random = ~1|Subject, na.action = na.omit)
summary(lme1)

# Inspect residuals
plot(lme1)

# Get estimates
intervals(lme1)

# Plot effects
eff1 <- effect("Treatment * FirstOrSecond", lme1)

pdf("Fig_STAIS.pdf", height = 4, width = 4*1.61)
par(mar = c(4,4,1,1))
plot(eff1$fit[c(2, 4)],
     frame.plot = F,
     xaxt = "n",
     yaxt = "n",
     type = "b",
     xlab = "",
     ylab = "STAI-S score",
     xlim = c(1, 2.1),
     ylim = c(32, 38),
     #main = "B. State anxiety",
     col = col1)
lines(c(1.1,2.1), eff1$fit[c(1, 3)], type = "b", col = col2, pch = 16)
lines(c(1, 1), c((eff1$upper[2]), (eff1$lower[2])), col = col1)
lines(c(2, 2), c((eff1$upper[4]), (eff1$lower[4])), col = col1)
lines(c(1.1, 1.1), c((eff1$upper[1]), (eff1$lower[1])), col = col2)
lines(c(2.1, 2.1), c((eff1$upper[3]), (eff1$lower[3])), col = col2)
axis(1, labels = c("Before", "After"), at = c(1.05, 2.05))
axis(2, at = c(32, 34, 36, 38))
#legend("top", col = c("blue", "red"), pch = c(1, 16), legend = c("Placebo", "Oxazepam"), bty = "n")
dev.off()


# Analyse participants' guesses of treatment group membership
demData$Guessed.group <- factor(demData$Guessed.group, levels = c("Placebo", "Likely_placebo", "Equivocal", "Likely_oxa", "Oxazepam"), ordered = TRUE)
demData$Guessed.group[demData$Included_FMOV == 0] <- NA

pdf("Fig_Blinding.pdf", height = 4, width = 4*1.61)
par(mar = c(4,4,1,1))
barplot(t(matrix(c(table(demData$Guessed.group[demData$Treatment == "Placebo"]), table(demData$Guessed.group[demData$Treatment == "Oxazepam"])), nr = 5)), 
        beside = TRUE, 
        names.arg = c("Placebo", "Likely placebo", "Equivocal", "Likely oxa", "Oxazepam"), 
        xlab = "Guessed group",
        ylab = "n",
        yaxt = "n",
        col = c(col3, col2),
        border = c(col1, col2),
        #main = "D. Efficacy of blinding",
        lwd = 5)
axis(2, at = c(0, 2, 4, 6, 8))
legend(c(0, 9), legend = c("Placebo", "Oxazepam"), fill = c(col3, col2), border = c(col1, col2), bty = "n")
dev.off()

test1 <- wilcox.test(as.numeric(Guessed.group) ~ Treatment, data = demData, alternative = "greater", paired = F, conf.int = T)
test1
