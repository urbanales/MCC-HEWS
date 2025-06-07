################################################################################
# SCRIPT NAME: 02.stage2_meta_regression.R
#
# DESCRIPTION:
#   Performs the second-stage meta-analysis of city-specific heat effects estimated
#   in the first stage. Models the modification of heat effects over time and across
#   regions, and analyzes the impact of Heat Prevention Plans (HPP). Produces
#   region-specific and policy-specific risk curves.
#
# INPUTS:
#   - data/tmeanperpar.csv: City-period coefficients from stage 1.
#   - data/avgtmeansum.csv: City-period temperature summaries.
#
# OUTPUTS:
#   - Main meta-regression model object (mainmod)
#   - Plots of exposure-response curves
#
# USAGE:
#   Run after completing stage 1 and saving required input files. Ensure all input files
#   are present in the "data/" directory. Run this script in R.
#
# AUTHORs: Aleš Urban, Veronika Huber, Pierre Masselot, Antonio Gasparrini
# DATE: June 2025
################################################################################

#===============================================================================
# 1. LOAD REQUIRED PACKAGES
#===============================================================================
library(mixmeta)     # For meta-analysis and mixed models
library(dlnm)        # For distributed lag non-linear models
library(splines)     # For spline-based modeling
library(FluMoDL)     # MCC-related modeling tools
library(nlme)        # Linear mixed-effects models
library(scales)      # For plotting
library(dplyr)       # Data manipulation
library(ggplot2)     # Graphics
library(lubridate)   # Date manipulation
library(readr)       # CSV reading
library(tidyverse)   # Tidy tools
library(lmtest)      # Diagnostic testing

#===============================================================================
# 2. LOAD AND PREPARE DATA
#===============================================================================

# Load city-specific coefficients and variances from first stage
tmeanperpar <- read.csv("data/tmeanperpar.csv", sep = ",")

# Filter out rows with missing data 
matrix <- subset(tmeanperpar, !is.na(hws))

# Extract coefficient and variance matrices
coef <- as.matrix(matrix[, grep("coef", names(matrix))])
vcov <- as.matrix(matrix[, grep("vcov", names(matrix))])
cityinfo <- matrix[, 1:21]

# Load mean summer temperature distribution to define exposure basis
avgtmeansum <- read.csv("data/avgtmeansum.csv")
tmean <- avgtmeansum$tmean
knots <- tmean[avgtmeansum$perc %in% paste0(c(50,90), ".0%")]
bvar <- onebasis(tmean, fun="bs", degree=2, knots=knots)

# Define percentiles and labels for plotting
xperc <- c(0,1,5,25,50,75,95,99,100)
xval <- tmean[avgtmeansum$perc %in% paste0(xperc, ".0%")]

#===============================================================================
# 3. FIT META-REGRESSION MODELS 
#===============================================================================

# Model 0: HPP (hws) only
# mod0 <- mixmeta(coef ~ hws, vcov, data=cityinfo, method="ml",
#                 random=~I(year-2005)|cityname, bscov="diag")

# Model 3: Year and Region interaction (selected model)
mod3 <- update(mod0, coef ~ hws + I(year - 2005) * Region)

# Model 4: Region-specific HPP effects
# mod4 <- update(mod0, coef ~ hws:Region + I(year - 2005) * Region)

# Model 5: HPP class-specific effects
# mod5 <- update(mod0, coef ~ hws:as.factor(hwclass) + I(year - 2005) * Region)

#===============================================================================
# 3.1 SENSITIVTY ANALYSIS
#===============================================================================

# AIC(mod0,mod1,mod2,mod3,mod4,mod5) 
# 
# # TEST OF A REGION EFFECT
# lrtest(mod0, mod1)
# 
# # TEST OF A TREND
# lrtest(mod1, mod2)
# 
# # TEST OF REGION-SPECIFIC TRENDS
# lrtest(mod2, mod3)
# 
# # TEST NOW HWS EFFECT 
# # waldtest(mod3, "hws")
# # drop1(mod3, test="Chisq")
# # lrtest(mod3, update(mod3, .~.-hws))
# 
# # TEST OF REGION-SPECIFIC HWS EFFECTS
# lrtest(mod3, mod4)
# 
# # TEST OF CLASS-SPECIFIC HWS EFFECTS
# lrtest(mod3, mod5)
# lrtest(mod5, mod6)
# 
# lrtest(mod0,mod1,mod2,mod3,mod4,mod5)

# drop1(mod0,test="Chisq")
# drop1(mod1,test="Chisq")
# drop1(mod2,test="Chisq")
# drop1(mod3,test="Chisq")
# drop1(mod4,test="Chisq")
# drop1(mod5,test="Chisq")
# 
# 
# logLik(mod0)
# logLik(mod1)
# logLik(mod2)
# logLik(mod3)
# logLik(mod4)
# logLik(mod5)
# 
# summary(mod0)$i2stat[1]
# summary(mod1)$i2stat[1]
# summary(mod2)$i2stat[1]
# summary(mod3)$i2stat[1]
# summary(mod4)$i2stat[1]
# summary(mod5)$i2stat[1]
# 
# summary(mod0)$qstat$Q[1]
# summary(mod1)$qstat$Q[1]
# summary(mod2)$qstat$Q[1]
# summary(mod3)$qstat$Q[1]
# summary(mod4)$qstat$Q[1]
# summary(mod5)$qstat$Q[1]


# Select final model
mainmod <- mod3

#===============================================================================
# 4. PREDICT AND PLOT REGION-SPECIFIC ERFs
#===============================================================================

# Define regions
regions <- c("Northern Europe", "Eastern Europe", "Western Europe", "Southern Europe")

# Prepare summary RR output matrix
RRsumlist <- array(NA, c(17,5))
RRsumlist[1,] <- c("Region", "HPPclass", "RR", "low", "high")
radek <- 0  # Row index tracker

# Set up plot panel
par(mfrow=c(2,2), mex=0.8, mgp=c(1,1,0), cex=0.8)

#===============================================================================
# FIGURE 3 - ERFs FOR FACTUAL AND COUNTER-FACTUAL SCENARIOS
#===============================================================================

# Loop over regions
for (r in seq_along(regions)) {
  
  # Subset data for region
  matrix0 <- subset(matrix, Region == regions[r])
  cityinfo0 <- matrix0[, 1:23]
  coef0 <- as.matrix(matrix0[, grep("coef", names(matrix0))])
  vcov0 <- as.matrix(matrix0[, grep("vcov", names(matrix0))])
  
  # Estimate overall curve for region (used to define centering)
  model0 <- mixmeta(coef0, vcov0, data=cityinfo0, method="ml", bscov="diag")
  cen <- tmean[which.min(bvar %*% coef(model0))]  # Minimum risk temperature
  
  # Predict fitted values with and without HPP
  datapred1 <- data.frame(Region = regions[r], hws = 0, year = 2015)
  datapred2 <- data.frame(Region = regions[r], hws = 1, year = 2015)
  
  pred1 <- predict(mainmod, datapred1, vcov=TRUE)
  pred2 <- predict(mainmod, datapred2, vcov=TRUE)
  
  cp1 <- crosspred(bvar, coef=pred1$fit, vcov=pred1$vcov, cen=cen,
                   model.link="log", at=tmean)
  cp2 <- crosspred(bvar, coef=pred2$fit, vcov=pred2$vcov, cen=cen,
                   model.link="log", at=tmean)
  
  # Plot temperature-response curves
  plot(cp1, ylim=c(0.8, 3), xlab="Temperature percentile", ylab="RR",
       lab=c(6,5,7), las=1, lwd=2, xaxt="n", ci="area", mgp=c(2.5,1,0),
       col=6, ci.arg=list(col=alpha(6, 0.3)), main=regions[r])
  
  lines(cp2, lwd=2, col=3, ci="area", ci.arg=list(col=alpha(3, 0.3)))
  axis(1, at=xval, labels=paste0(xperc, "%"), cex.axis=1)
  abline(v=cp1$cen, lty=2, col=grey(0.7))
  abline(v=xval[8], lty=2, col="red")
  legend("topleft", c("counterfactual (w/o HPP)", "factual (w/ HPP)"),
         lwd=2, col=c(6,3), bty="n", inset=0.02)
}

#===============================================================================
# 5. BLUP CALCULATION FROM META-REGRESSION MODEL
#===============================================================================

# Predict fitted coefficients from the main model
predcoeff <- predict(mainmod, vcov=TRUE)
names(predcoeff) <- paste(cityinfo$cityname, cityinfo$year)

# Best Linear Unbiased Predictions (BLUPs) - full (fixed + random effects)
blup_fix <- blup(mainmod, vcov=TRUE, type="outcome")
names(blup_fix) <- paste(cityinfo$cityname, cityinfo$year)

# BLUPs - random effects only (residuals)
blup_rnd <- blup(mainmod, vcov=TRUE, type="residual")
names(blup_rnd) <- paste(cityinfo$cityname, cityinfo$year)

