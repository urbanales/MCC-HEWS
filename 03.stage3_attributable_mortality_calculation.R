################################################################################
# SCRIPT NAME: 03.stage3_attributable_mortality_calculation.R
#
# DESCRIPTION:
#   This script performs the third-stage analysis: using BLUPs from the second-stage
#   meta-regression to calculate Attributable Numbers (AN) and Attributable Fractions (AF)
#   across cities and countries, with uncertainty quantified via simulation.
#
# INPUTS:
#   - Output from stage 2: mainmod, blup_rnd
#   - City/country data objects from stage 1 (e.g., eu_cities, dlist, hws_ind)
#
#   NOTE: The input file "dlist" cannot be publicly provided due to
#   data sharing restrictions. This script demonstrates the methods and steps used 
#   to obtain the AF and AN estimates in stage 3; to reproduce the exact results, 
#   researchers can contact the corresponding author (Aleš Urban, urban@ufa.cas.cz)
#   for information on accessing the data used for this study and for the R code. 
#   The outputs based on this stage are provided in "data".
#
# OUTPUTS:
#   - Attributable mortality estimates -
#       * af_table_countries_hws0, af_table_countries_hws1
#   - Simulation results for uncertainty quantification -  
#       * ansimlist_0_hws0, ansimlist_1_hws0, ansimlist_0_hws1, ansimlist_1_hws1
#       * afsimlist_0_hws0, afsimlist_1_hws0, afsimlist_0_hws1, afsimlist_1_hws1
#
# USAGE:
#   This script must be run twice: once for the counterfactual scenario (hws = 0), and
#   once for the factual scenario (hws = 1). Adjust the relevant variable accordingly.
#
# AUTHORs: Aleš Urban, Veronika Huber, Pierre Masselot, Antonio Gasparrini
# DATE: June 2025
################################################################################

#===============================================================================
# 1. LOAD REQUIRED PACKAGES
#===============================================================================

library(MASS)        # For multivariate normal simulation
library(dlnm)        # For distributed lag non-linear models
library(dplyr)       # Data manipulation
library(lubridate)   # Date handling
library(pracma)      # For ndims

#===============================================================================
# 2. LOAD REQUIRED DATA
#===============================================================================
# Assumes mainmod and blup_rnd already available 
# Assumes eu_cities, dlist, and hws_ind are already loaded from stage 1

eu_cities <- read.csv2("data/eu_cities.csv", sep = ";")
hws_ind <- read_delim("data/hws_ind.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)


#===============================================================================
# 3. SIMULATE COEFFICIENTS FROM MAIN META-REGRESSION MODEL
#===============================================================================
nsim <- 1000  # Consider using at least 1000 for final analysis
set.seed(12345)
metacoefsim <- mvrnorm(nsim, coef(mainmod), vcov(mainmod))

#===============================================================================
# 4. PREPARE STRUCTURES FOR STORAGE
#===============================================================================
country <- levels(as.factor(eu_cities$countryname))
citynames <- eu_cities$cityname

country <- levels(as.factor(eu_cities$countryname))

###############!!!
# Rename all objects with "hws1" to "hws0" to get outputs for both scenarios
###############!!!

ansimlist_0_hws0 <- ansimlist_1_hws0 <- 
  afsimlist_0_hws0 <- afsimlist_1_hws0 <- 
  matrix(nrow = length(country), ncol = nsim,
         dimnames = list(country, c(1:nsim)))

ancitysimlist_0_hws0 <- ancitysimlist_1_hws0 <- 
  afcitysimlist_0_hws0 <- afcitysimlist_1_hws0 <- 
  matrix(nrow = length(eu_cities$city), ncol = nsim,
         dimnames = list(eu_cities$cityname, c(1:nsim)))

# Tables for storing results
af_table_countries_hws0 <- array(NA, c(length(country)*100, 11))
colnames(af_table_countries_hws0) <- c("region","hwsclass","country","year","hws","AN","ANlow","ANhigh","AF","AFlow","AFhigh")
row_index <- 0

af_table_cities_hws0 <- array(NA,c(length(eu_cities$city)*10+1,12))
colnames(af_table_cities_hws0) <- c("region","hwsclass","country","city","year","hws","AN","ANlow","ANhigh","AF","AFlow","AFhigh")
row <- 0
r <- 0


#===============================================================================
# 6. EXPANDED LOOP: COUNTRY -> CITY -> PERIOD
#===============================================================================

af_table_countries_hws0 <- array(NA,c(length(country)*100,11))
colnames(af_table_countries_hws0) <- c("region","hwsclass","country","year","hws","AN","ANlow","ANhigh","AF","AFlow","AFhigh")
row_index <- 0

af_table_cities_hws0 <- array(NA,c(length(eu_cities$city)*10+1,12))
colnames(af_table_cities_hws0) <- c("region","hwsclass","country","city","year","hws","AN","ANlow","ANhigh","AF","AFlow","AFhigh")
row <- 0
r <- 0

for (r in 1:length(country)){
  
  # I select the cities for the region as an example
  cities <- eu_cities[(which(eu_cities$countryname == country[r])),]
  eucities_select <- print(eucities[cities$city])
  
  
  #----- Now we can loop across cities and years
  
  # Objects to store
  yearlistII <- list(a=1990:1992,b=1993:1995,c=1996:1998,
                     d=1999:2001,e=2002:2004,f=2005:2007,
                     g=2008:2010,h=2011:2013,i=2014:2016,
                     j=2017:2019)  
  
  anlist_hws0 <- aflist_hws0 <- hwslist_hws0 <- matrix(nrow = nrow(cities), ncol = length(yearlistII), 
                                                       dimnames = list(cities$city, names(yearlistII)))
  ansimlist_hws0 <- afsimlist_hws0 <- array(dim = c(nrow(cities), length(yearlistII), nsim),
                                            dimnames = list(cities$city, names(yearlistII), NULL))
  
  # Loop across cities
  city <- 1
  for(city in seq(nrow(cities))) {
    
    cat(cities$cityname[city],"")
    
    # SELECT VARIABLES
    
    data <- eucities_select[[city]]
    vars <- c("date","dow","death","temp")
    data$date <- as.Date(data$date,format="%d.%m.%Y")
    rok <- as.factor(hws_ind$year)
    df <- data.frame(hws_ind$year,hws_ind[,which(colnames(hws_ind)==cities$countryname[city])])
    names(df)<-c("year","hws")
    
    data$hwclass <- NA
    data$hws <- NA
    
    y <- 0
    
    for (y in 1:length(levels(rok))){
      
      data$hwclass[which(data$year==rok[y])] <- df$hws[which(df$year==rok[y])]
      if (is.na(df$hws[which(df$year==rok[y])])){
        y=y+1}
      else if (df$hws[which(df$year==rok[y])] > 0){
        data$hws[which(data$year==rok[y])] <- 1
      } else {
        data$hws[which(data$year==rok[y])] <- 0
      }
      
    }
    
    # SUBSET FOR SUMMER-ONLY
    datasum <- subset(data, month(date) %in% 5:9)
    all <- datasum[,7]
    datasum$tempdev<- (datasum$tmean - quantile(datasum$tmean,0.95,na.rm=T))
    
    ################################################################################
    # ANALYSIS OF TEMPERATURE - ALL-CAUSE MORTALITY (SUMMER-ONLY, BY PERIOD)
    
    # DEFINE THE PERIODS
    
    yearlist <- list(a=1990:1992,b=1993:1995,c=1996:1998,
                     d=1999:2001,e=2002:2004,f=2005:2007,
                     g=2008:2010,h=2011:2013,i=2014:2016,
                     j=2017:2019)  
    
    j <- 0
    
    for (j in c("a","b","c","d","e","f","g","h","i","j")){
      
      if (FALSE %in% (yearlist[[j]] %in% year(datasum$date))) {
        yearlist<-yearlist[names(yearlist) %in% j == FALSE]
      }
    }
    
    yearlist
    
    y<-1
    
    # PERFORM MODEL BY PERIOD
    for (y in seq(length(yearlist))){
      
      
      predper = c(seq(0,100))
      indlab = predper%in%c(0,50,75,90,95,100)
      temp = datasum$tmean
      summary(temp)
      
      dat0 <- subset(datasum,datasum$year==yearlist[[y]][1]|datasum$year==yearlist[[y]][2]|datasum$year==yearlist[[y]][3])
      all0 <- dat0[,7]
      temp0 <- as.numeric(dat0$tmean)
      summary(temp0)
      
      tempper0 = c(seq(min(temp0,na.rm=T),max(temp0,na.rm=T),0.1))
      Tmeancountry0 = quantile(tempper0, predper/100,na.rm=T)
      
      
      # Here we can just create the exposure dimension since we have coefficients for the overall cumulative ERF
      btmean0 <- onebasis(temp0, fun="bs", degree=2,
                          knots=quantile(temp, c(50,90)/100, na.rm=T))
      
      #----- We have to prepare simulated coefficients
      
      ###############!!!
      # Set hws = 0 or hws = 1 here before running the script for each scenario!
      ###############!!!
      
      # Data for prediction 
      datapred1 <- data.frame(Region=cities$Region[1],hws=0, year=yearlist[[y]][2]) # Expanded to all years
      
      # We need to build the model matrix to predict (adapted from predict.mixmeta function)
      tt <- delete.response(mainmod$terms)
      contr <- mixmeta:::getContrXlev(tt, mainmod$contrasts)
      xlev <- mixmeta:::getContrXlev(tt, mainmod$xlev)
      mf <- model.frame(tt, datapred1, xlev = xlev)
      modmat1 <- model.matrix(tt, mf, contr)
      
      # Predict the point estimate for fixed effect
      pred1 <- predict(mainmod, datapred1)
      
      # Predict the point estimate for fixed effect (using the model matrix)
      pred1_alt <- coef(mainmod) %*% t(modmat1 %x% diag(mainmod$dim$k)) # Kronecker product because multivariate meta-analysis
      pred1; pred1_alt # They are identical, good
      
      # Now we can compute simulated spline coefficients for the region
      pred1_sim <- metacoefsim %*% t(modmat1 %x% diag(mainmod$dim$k))
      
      
      #----- We can compute the point estimate for AN
      
      # We need to add the BLUP residual 
      pred1blup <- pred1 + blup_rnd[[paste(cities$cityname[city],yearlist[[y]][2])]]$blup
      
      # Find the MMT and center basis
      mmt <- median(temp0, na.rm=T) #temp0[which.min(btmean0 %*% pred1blup)]
      cenbasis <- onebasis(mmt, fun="bs", degree=2,
                           knots=quantile(temp, c(50,90)/100, na.rm=T), Boundary.knots = range(temp0, na.rm=T))
      btmean0cen <- scale(btmean0, center = cenbasis, scale = F)
      
      # index for days above 95th percentile
      heatind <- temp0 >= quantile(temp, 0.95, na.rm=T)
      
      # We can now compute the AN
      anday <- (1 - exp(-btmean0cen %*% pred1blup)) * all0
      anlist_hws0[city, names(yearlist)[y]] <- ((sum(anday[heatind], na.rm = T)/cities$pop[city])*100000)/3
      
      # And the AF (add [heatind] if computed for hot days only)
      aflist_hws0[city, names(yearlist)[y]] <- sum(anday[heatind], na.rm = T) / sum(all0,na.rm=T)*100
      
      hwslist_hws0[city, names(yearlist)[y]] <- round(mean(dat0$hws))
      
      #----- Now we have to do the same for all simulated coefficients
      sim <- 1
      for (sim in seq_len(nsim)){
        
        # Coefficients (we have to add the BLUP residuals one more time)
        sim1blup <- pred1_sim[sim,] + blup_rnd[[paste(cities$cityname[city],yearlist[[y]][2])]]$blup
        
        # We can now compute AN as before
        andaysim <- (1 - exp(-btmean0cen %*% sim1blup)) * all0
        ansimlist_hws0[city, names(yearlist)[y], sim] <- ((sum(andaysim[heatind], na.rm = T)/cities$pop[city])*100000)/3
        
        # And the AF
        afsimlist_hws0[city, names(yearlist)[y], sim] <- sum(andaysim[heatind], na.rm = T) / 
          sum(all0,na.rm=T)*100
      }
    } 
  }
  
  #----- Confidence intervals at city/period level
  
  # Compute CI as quantiles of the simulated ANs
  AN <- anlist_hws0
  ANlow <- apply(ansimlist_hws0, 1:2, quantile, .025, na.rm = T)
  ANhigh <- apply(ansimlist_hws0, 1:2, quantile, .975, na.rm = T)
  
  # We can look at an example (Gothenburg)
  anlist_hws0[1,]; ANlow[1,]; ANhigh[1,]
  
  # Same for AF
  AF <- aflist_hws0
  AFlow <- apply(afsimlist_hws0, 1:2, quantile, .025, na.rm = T)
  AFhigh <- apply(afsimlist_hws0, 1:2, quantile, .975, na.rm = T)
  
  aflist_hws0[1,]; AFlow[1,]; AFhigh[1,]
  
  HWS <- hwslist_hws0
  
  j<-1
  for (j in 1:length(cities$cityname)){
    df1 <- data.frame(cities$Region[j],cities$hwc[j],cities$countryname[j],cities$cityname[j],as.numeric(sapply(yearlistII, mean)),HWS[j,],as.numeric(AN[j,]),as.numeric(ANlow[j,]),as.numeric(ANhigh[j,]),
                      as.numeric(AF[j,]),as.numeric(AFlow[j,]),as.numeric(AFhigh[j,]))
    k<-0
    for (k in 1:length(df1)){
      af_table_cities_hws0[row+1:length(yearlistII),k] <- df1[,k]
    }
    row <- row+length(yearlistII)
  }
  # ----- Now for regions we have to sum city-level ANs
  
  HWShws0 <- round(colMeans(hwslist_hws0, na.rm = T))
  
  # Point estimate
  ANhws0 <- colMeans(anlist_hws0, na.rm = T)
  
  ANhws00 <- if(length(anlist_hws0[,HWShws0=="0"])==1){
    mean(anlist_hws0[,HWShws0=="0"])} else if (pracma::ndims(anlist_hws0[,HWShws0=="0"])==1){
      mean(anlist_hws0[,HWShws0=="0"])} else {
        mean(colMeans(anlist_hws0[,HWShws0=="0"], na.rm = T))
      }
  
  ANhws01 <- if(length(anlist_hws0[,HWShws0=="1"])==1){
    mean(anlist_hws0[,HWShws0=="1"])} else if (pracma::ndims(anlist_hws0[,HWShws0=="1"])==1){
      mean(anlist_hws0[,HWShws0=="1"])} else {
        mean(colMeans(anlist_hws0[,HWShws0=="1"], na.rm = T))
      }
  
  # Confidence intervals, we first have to sum ANs for all cities for each simulation
  anhws0sim <- apply(ansimlist_hws0, 2:3, mean, na.rm = T) # sum all cities
  ANhws0high <- apply(anhws0sim, 1, quantile, .975, na.rm = T) # Quantiles on sum
  ANhws0low <- apply(anhws0sim, 1, quantile, .025, na.rm = T) # Quantiles on sum
  
  ANhws0; ANhws0low; ANhws0high
  
  anhws0sim0 <- apply(ansimlist_hws0, 2:3, mean, na.rm = T) # sum all cities
  # anhws0sim0 <- anhws0sim0[HWShws0=="0",]         # select periods
  anhws0sim0 <- if(pracma::ndims(anhws0sim0[HWShws0=="0",])==1){
    anhws0sim0[HWShws0=="0",]} else {
      apply(anhws0sim0[HWShws0=="0",], 2, mean, na.rm = T)}  # sum of periods
  ANhws0high0 <- quantile(anhws0sim0,.975, na.rm = T) # Quantiles of sums
  ANhws0low0 <- quantile(anhws0sim0,.025, na.rm = T) # Quantiles of sums
  
  ansimlist_0_hws0[country[r],] <- anhws0sim0 # safe simulations
  
  ANhws00; ANhws0low0; ANhws0high0
  
  anhws0sim1 <- apply(ansimlist_hws0, 2:3, mean, na.rm = T) # sum all cities
  anhws0sim1 <- if (pracma::ndims(anhws0sim1[HWShws0=="1",])==1){
    anhws0sim1[HWShws0=="1",]} else {
      apply(anhws0sim1[HWShws0=="1",], 2, mean, na.rm = T)}  # sum of periods
  ANhws0high1 <- quantile(anhws0sim1,.975, na.rm = T) # Quantiles on sum
  ANhws0low1 <- quantile(anhws0sim1,.025, na.rm = T) # Quantiles on sum
  
  ansimlist_1_hws0[country[r],] <- anhws0sim1 # safe simulations
  
  ANhws01; ANhws0low1; ANhws0high1
  
  
  #----- Now for regions we have to average city-level AFs
  
  # Point estimate
  AFhws0 <- colMeans(aflist_hws0, na.rm = T)
  AFhws00 <- if(pracma::ndims(aflist_hws0[,HWShws0=="0"])==1){mean(aflist_hws0[,HWShws0=="0"])} else {mean(colMeans(aflist_hws0[,HWShws0=="0"], na.rm = T))}
  AFhws01 <- if(pracma::ndims(aflist_hws0[,HWShws0=="1"])==1){mean(aflist_hws0[,HWShws0=="1"])} else {mean(colMeans(aflist_hws0[,HWShws0=="1"], na.rm = T))}
  
  # Confidence intervals, we first have to mean AFs for all cities for each simulation
  afhws0sim <- apply(afsimlist_hws0, 2:3, mean, na.rm = T) # mean all cities
  AFhws0high <- apply(afhws0sim, 1, quantile, .975, na.rm = T) # Quantiles on mean
  AFhws0low <- apply(afhws0sim, 1, quantile, .025, na.rm = T) # Quantiles on mean
  
  AFhws0; AFhws0low; AFhws0high
  
  afhws0sim0 <- apply(afsimlist_hws0, 2:3, mean, na.rm = T) # mean all cities
  afhws0sim0 <- if (pracma::ndims(afhws0sim0[HWShws0=="0",])==1){
    afhws0sim0[HWShws0=="0",]} else {
      apply(afhws0sim0[HWShws0=="0",], 2, mean, na.rm = T)}  # mean of periods
  AFhws0high0 <- quantile(afhws0sim0, .975, na.rm = T) # Quantiles on mean
  AFhws0low0 <- quantile(afhws0sim0, .025, na.rm = T) # Quantiles on mean
  
  afsimlist_0_hws0[country[r],] <- afhws0sim0 # safe simulations
  
  AFhws00; AFhws0low0; AFhws0high0
  
  afhws0sim1 <- apply(afsimlist_hws0, 2:3, mean, na.rm = T) # mean all cities
  afhws0sim1 <- if (pracma::ndims(afhws0sim1[HWShws0=="1",])==1){
    afhws0sim1[HWShws0=="1",]} else {
      apply(afhws0sim1[HWShws0=="1",], 2, mean, na.rm = T)}  # mean of periods
  AFhws0high1 <- quantile(afhws0sim1,.975, na.rm = T) # Quantiles on mean
  AFhws0low1 <- quantile(afhws0sim1,.025, na.rm = T) # Quantiles on mean
  
  afsimlist_1_hws0[country[r],] <- afhws0sim1 # safe simulations
  
  AFhws01; AFhws0low1; AFhws0high1
  
  df2 <- data.frame(cities$Region[1],cities$hwc[1],country[r],as.numeric(sapply(yearlistII, mean)),HWShws0,as.numeric(ANhws0),as.numeric(ANhws0low),as.numeric(ANhws0high),
                    as.numeric(AFhws0),as.numeric(AFhws0low),as.numeric(AFhws0high))
  
  i<-0
  for (i in 1:length(df2)){
    af_table_countries_hws0[row_index+1:length(yearlistII),i] <- df2[,i]
  }
  row_index <- row_index+length(yearlistII)
  
  row_index <- row_index+1
  af_table_countries_hws0[row_index,] <- c(cities$Region[1],cities$hwc[1],country[r],"period1","0",as.numeric(ANhws00),as.numeric(ANhws0low0),as.numeric(ANhws0high0),
                                       as.numeric(AFhws00),as.numeric(AFhws0low0),as.numeric(AFhws0high0))
  row_index <- row_index+1
  af_table_countries_hws0[row_index,] <- c(cities$Region[1],cities$hwc[1],country[r],"period2","1",as.numeric(ANhws01),as.numeric(ANhws0low1),as.numeric(ANhws0high1),
                                       as.numeric(AFhws01),as.numeric(AFhws0low1),as.numeric(AFhws0high1))
}

#===============================================================================
# 7. EXPORT FINAL OUTPUTS (UNCOMMENT TO SAVE)
#===============================================================================
# write.csv(af_table_countries_hws0, "data/af_table_countries_hws0.csv", row.names=FALSE)
# write.csv(ansimlist_0_hws0, "data/ansimlist_0_hws0.csv", row.names=FALSE)
# write.csv(ansimlist_1_hws0, "data/ansimlist_1_hws0.csv", row.names=FALSE)


