################################################################################
# SCRIPT NAME: 01.stage1_period_specific_erfs.R
#
# DESCRIPTION:
#   Performs the first-stage analysis for the MCC-HEWS project. This script fits
#   city-specific distributed lag non-linear models (DLNMs) on summer temperature
#   and mortality data to estimate excess risk functions (ERFs) for each city and
#   time period. It saves coefficients and summaries for use in later meta-analysis.
#
#   NOTE: The input file "MCCdata_20230830.RData" cannot be publicly provided due to
#   data sharing restrictions. This script demonstrates the methods and steps used 
#   to obtain the ERF coefficients in stage 1; to reproduce the exact results, 
#   researchers can contact the corresponding author (Aleš Urban, urban@ufa.cas.cz)
#   for information on accessing the data used for this study and for the R code. 
#   The ERF  coefficients based on this dataset are provided in tmeanperpar.csv.
#
# INPUTS:
#   - data/MCCdata_20230830.RData: Main data object with city-wise data lists (not provided).
#   - data/eu_cities.csv: Metadata about EU cities.
#   - data/hws_ind.csv: national HPP class by year.
#
# OUTPUTS:
#   - R objects: tmeanpar, tmeanperpar, tmeansum, avgtmeansum (can be saved as CSV)
#
# USAGE:
#   Ensure all input files are present in the "data/" directory. Run this script in R.
#
# AUTHOR: Urban et al. The effectiveness of heat prevention plans to reduce heat-related mortality in Europe
# DATE: June 2025
################################################################################

# --- Setup Environment ---
# Place your input files in the data/ folder.

library(dlnm)
library(mixmeta)
library(tsModel)
library(splines)
library(lubridate)
library(readr)

# --- Load Data ---
load("data/MCCdata_20230830.RData")
eu_cities <- read.csv2("data/eu_cities.csv", sep = ";")
hws_ind <- read_delim("data/hws_ind.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)

eucities <- dlist[eu_cities$city]

# --- Initialize Result Containers ---
n_cities <- nrow(eu_cities)
tmeanparlist <- vector("list", n_cities)
tmeanperparlist <- vector("list", n_cities)
tmeansumlist <- vector("list", n_cities)

i<-1

# --- Main Loop ---
for (i in seq_len(n_cities)) {
  
  cat("Processing:", eu_cities$cityname[i], "\n")
  
  data <- eucities[[i]]
  data$date <- as.Date(data$date, format="%d.%m.%Y")
  data$year <- year(data$date)
  
  # Match HWS index
  df <- data.frame(year = hws_ind$year, hws = hws_ind[[eu_cities$countryname[i]]])
  data$hwclass <- NA
  data$hws <- NA
  
  for (y in df$year) {
    idx <- data$year == y
    val <- df$hws[df$year == y]
    data$hwclass[idx] <- val
    if (!is.na(val)) {
      data$hws[idx] <- as.integer(val > 0)
    }
  }
  
  # Subset to summer months
  datasum <- subset(data, month(date) %in% 5:9)
  all <- datasum[,7]  # assuming death is column 7
  datasum$tempdev <- datasum$tmean - quantile(datasum$tmean, 0.95, na.rm=TRUE)
  
  # --- DLNM Model: Full Summer ---
  cbtmean <- crossbasis(datasum$tmean, lag=10,
                        argvar=list(fun="bs", degree=2,
                                    knots=quantile(datasum$tmean, c(0.5, 0.9), na.rm=TRUE)),
                        arglag=list(knots=logknots(10,2)),
                        group=year(datasum$date))
  
  model <- glm(all ~ cbtmean + dow + ns(yday(date), df=4)*factor(year(date)),
               data=datasum, family=quasipoisson)
  redpred <- crossreduce(cbtmean, model, cen=mean(datasum$tmean, na.rm=TRUE))
  
  # Store Coefficients
  ncoef <- length(coef(redpred))
  par <- c(coef(redpred), vechMat(vcov(redpred)))
  names(par) <- c(paste0("coef", seq_len(ncoef)),
                  paste0("vcov", seq_len(ncoef*(ncoef+1)/2)))
  tmeanparlist[[i]] <- data.frame(eu_cities[i, c("city", "cityname", "countryname", "Region", "pop",
                                                 "gdp", "lat", "long", "kgclzone", "MeanST", "IQR", "Range",
                                                 "hwscore","hwclass")], t(par), row.names=i)
  
  # --- By Period Model ---
  yearlist <- list(
    a=1990:1992, b=1993:1995, c=1996:1998, d=1999:2001,
    e=2002:2004, f=2005:2007, g=2008:2010, h=2011:2013,
    i=2014:2016, j=2017:2019)
  
  # Remove unused periods
  yearlist <- yearlist[sapply(yearlist, function(x) all(x %in% year(datasum$date)))]
  
  parlist <- lapply(yearlist, function(x) {
    model <- glm(all ~ cbtmean + dow + ns(yday(datasum$date), df=4),
                 data=datasum, family=quasipoisson, subset=year(datasum$date) %in% x)
    redpred <- crossreduce(cbtmean, model, cen=mean(datasum$tmean, na.rm=TRUE))
    t(c(coef(redpred), vechMat(vcov(redpred))))
  })
  
  # Aggregate period-based metrics
  heatlist <- sapply(yearlist, function(x) {
    sum(subset(datasum, year(date) %in% x & tempdev > 0, na.rm=TRUE)$tempdev)
  })
  
  tmeanlist <- sapply(yearlist, function(x) {
    mean(subset(datasum, year(date) %in% x)$tmean, na.rm=TRUE)
  })
  
  hwclist <- sapply(yearlist, function(x) {
    as.integer(median(subset(datasum, year(date) %in% x)$hwclass, na.rm=TRUE))
  })
  
  hwslist <- sapply(yearlist, function(x) {
    as.integer(mean(subset(datasum, year(date) %in% x)$hws, na.rm=TRUE))
  })
  
  deathlist <- sapply(yearlist, function(x) {
    sum(subset(datasum, year(date) %in% x)[,7], na.rm=TRUE)
  })
  
  par <- do.call(rbind, parlist)
  colnames(par) <- names(tmeanparlist[[i]])[-(1:14)]
  
  tmeanperparlist[[i]] <- data.frame(
    eu_cities[i, c("city", "cityname", "countryname", "Region", "pop", "gdp", "lat", "long",
                   "kgclzone", "MeanST", "IQR", "Range", "hwscore","hwclass")],
    period = sapply(yearlist, function(x) paste(range(x), collapse="-")),
    year = sapply(yearlist, mean),
    hwi = heatlist,
    hws = hwslist,
    hwc = hwclist,
    die = deathlist,
    tmean = tmeanlist,
    par,
    row.names=paste0(i, ".", seq_along(yearlist))
  )
  
  # --- Temperature Distribution (Summer Only) ---
  per <- c(0:9/10, 1:99, 991:1000/10)/100
  tmeansumlist[[i]] <- quantile(datasum$tmean, per, na.rm=TRUE)
}

# --- Final Aggregation ---
tmeanpar <- do.call(rbind, tmeanparlist)
tmeanperpar <- do.call(rbind, tmeanperparlist)
tmeansum <- do.call(rbind, tmeansumlist)
avgtmeansum <- data.frame(perc=names(tmeansumlist[[1]]),
                          tmean=apply(do.call(cbind, tmeansumlist), 1, mean))

# --- Write Output ---
# write.csv(tmeanpar, "data/tmeanpar.csv", row.names=FALSE)
# write.csv(tmeanperpar, "data/tmeanperpar.csv", row.names=FALSE)
# write.csv(avgtmeansum, "data/avgtmeansum.csv", row.names=FALSE)
# write.csv(tmeansum, "data/tmeansum.csv", row.names=tmeanpar$cityname)
