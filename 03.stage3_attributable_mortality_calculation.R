
################################################################################
# DESCRIPTION:
# This script performs the third-stage analysis: using BLUPs from the second-stage
# meta-regression to calculate Attributable Numbers (AN) and Attributable Fractions (AF)
# across cities and countries, with uncertainty quantified via simulation.

# IMPORTANT:
# This script needs to be run twice: once for the counterfactual scenario (hws = 0),
# and once for the factual scenario (hws = 1), adjusting the relevant variable in accordingly.
# I.e. this requires changing the hws parameter in line 200
# All objects including "hws0" from line 60 further need to be renamed.
################################################################################

#===============================================================================
# 1. LOAD REQUIRED PACKAGES
#===============================================================================
library(MASS)        # For multivariate normal simulation
library(dplyr)       # Data manipulation
library(lubridate)   # Date handling
library(pracma)      # For ndims

#===============================================================================
# 2. LOAD REQUIRED DATA
#===============================================================================
# Assumes mainmod and blup_rnd already available 
# Assumes eu_cities, dlist, and hws_ind are already loaded from stage 1

check_required_objects <- function(required) {
  missing <- required[!required %in% ls(envir = .GlobalEnv)]
  if (length(missing) == 0) {
    message("All required objects are present.")
  } else {
    warning("Missing objects: ", paste(missing, collapse = ", "))
  }
}

required_objects <- c("mainmod","blup_rnd",
                      "eu_cities", "dlist",
                      "hws_ind")

check_required_objects(required_objects)


#===============================================================================
# 3. SIMULATE COEFFICIENTS FROM MAIN META-REGRESSION MODEL
#===============================================================================
nsim <- 100  # Consider using at least 1000 for final analysis
set.seed(12345)
metacoefsim <- mvrnorm(nsim, coef(mainmod), vcov(mainmod))

#===============================================================================
# 4. PREPARE STRUCTURES FOR STORAGE
#===============================================================================
country <- levels(as.factor(eu_cities$countryname))
citynames <- eu_cities$cityname

country <- levels(as.factor(eu_cities$countryname))

ansimlist_0_hws1 <- ansimlist_1_hws1 <- 
  afsimlist_0_hws1 <- afsimlist_1_hws1 <- 
  matrix(nrow = length(country), ncol = nsim,
         dimnames = list(country, c(1:nsim)))

ancitysimlist_0_hws1 <- ancitysimlist_1_hws1 <- 
  afcitysimlist_0_hws1 <- afcitysimlist_1_hws1 <- 
  matrix(nrow = length(eu_cities$city), ncol = nsim,
         dimnames = list(eu_cities$cityname, c(1:nsim)))

# Tables for storing results
af_table_countries_hws1 <- array(NA, c(length(country)*100, 11))
colnames(af_table_countries_hws1) <- c("region","hwsclass","country","year","hws","AN","ANlow","ANhigh","AF","AFlow","AFhigh")
radek <- 0

af_table_cities_hws1 <- array(NA,c(length(eu_cities$city)*10+1,12))
colnames(af_table_cities_hws1) <- c("region","hwsclass","country","city","year","hws","AN","ANlow","ANhigh","AF","AFlow","AFhigh")
row <- 0
r <- 0


#===============================================================================
# 6. EXPANDED LOOP: COUNTRY -> CITY -> PERIOD
#===============================================================================

af_table_countries_hws1 <- array(NA,c(length(country)*100,11))
colnames(af_table_countries_hws1) <- c("region","hwsclass","country","year","hws","AN","ANlow","ANhigh","AF","AFlow","AFhigh")
radek <- 0

af_table_cities_hws1 <- array(NA,c(length(eu_cities$city)*10+1,12))
colnames(af_table_cities_hws1) <- c("region","hwsclass","country","city","year","hws","AN","ANlow","ANhigh","AF","AFlow","AFhigh")
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
  
  anlist_hws1 <- aflist_hws1 <- hwslist_hws1 <- matrix(nrow = nrow(cities), ncol = length(yearlistII), 
                                                       dimnames = list(cities$city, names(yearlistII)))
  ansimlist_hws1 <- afsimlist_hws1 <- array(dim = c(nrow(cities), length(yearlistII), nsim),
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
      
      # Data for prediction 
      datapred1 <- data.frame(Region=cities$Region[1],hws=1, year=yearlist[[y]][2]) # Expanded to all years
      
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
      anlist_hws1[city, names(yearlist)[y]] <- ((sum(anday[heatind], na.rm = T)/cities$pop[city])*100000)/3
      
      # And the AF (add [heatind] if computed for hot days only)
      aflist_hws1[city, names(yearlist)[y]] <- sum(anday[heatind], na.rm = T) / sum(all0,na.rm=T)*100
      
      hwslist_hws1[city, names(yearlist)[y]] <- round(mean(dat0$hws))
      
      #----- Now we have to do the same for all simulated coefficients
      sim <- 1
      for (sim in seq_len(nsim)){
        
        # Coefficients (we have to add the BLUP residuals one more time)
        sim1blup <- pred1_sim[sim,] + blup_rnd[[paste(cities$cityname[city],yearlist[[y]][2])]]$blup
        
        # We can now compute AN as before
        andaysim <- (1 - exp(-btmean0cen %*% sim1blup)) * all0
        ansimlist_hws1[city, names(yearlist)[y], sim] <- ((sum(andaysim[heatind], na.rm = T)/cities$pop[city])*100000)/3
        
        # And the AF
        afsimlist_hws1[city, names(yearlist)[y], sim] <- sum(andaysim[heatind], na.rm = T) / 
          sum(all0,na.rm=T)*100
      }
    } 
  }
  
  #----- Confidence intervals at city/period level
  
  # Compute CI as quantiles of the simulated ANs
  AN <- anlist_hws1
  ANlow <- apply(ansimlist_hws1, 1:2, quantile, .025, na.rm = T)
  ANhigh <- apply(ansimlist_hws1, 1:2, quantile, .975, na.rm = T)
  
  # We can look at an example (Gothenburg)
  anlist_hws1[1,]; ANlow[1,]; ANhigh[1,]
  
  # Same for AF
  AF <- aflist_hws1
  AFlow <- apply(afsimlist_hws1, 1:2, quantile, .025, na.rm = T)
  AFhigh <- apply(afsimlist_hws1, 1:2, quantile, .975, na.rm = T)
  
  aflist_hws1[1,]; AFlow[1,]; AFhigh[1,]
  
  HWS <- hwslist_hws1
  
  j<-1
  for (j in 1:length(cities$cityname)){
    df1 <- data.frame(cities$Region[j],cities$hwc[j],cities$countryname[j],cities$cityname[j],as.numeric(sapply(yearlistII, mean)),HWS[j,],as.numeric(AN[j,]),as.numeric(ANlow[j,]),as.numeric(ANhigh[j,]),
                      as.numeric(AF[j,]),as.numeric(AFlow[j,]),as.numeric(AFhigh[j,]))
    k<-0
    for (k in 1:length(df1)){
      af_table_cities_hws1[row+1:length(yearlistII),k] <- df1[,k]
    }
    row <- row+length(yearlistII)
  }
  # ----- Now for regions we have to sum city-level ANs
  
  HWShws1 <- round(colMeans(hwslist_hws1, na.rm = T))
  
  # Point estimate
  ANhws1 <- colMeans(anlist_hws1, na.rm = T)
  
  ANhws10 <- if(length(anlist_hws1[,HWShws1=="0"])==1){
    mean(anlist_hws1[,HWShws1=="0"])} else if (pracma::ndims(anlist_hws1[,HWShws1=="0"])==1){
      mean(anlist_hws1[,HWShws1=="0"])} else {
        mean(colMeans(anlist_hws1[,HWShws1=="0"], na.rm = T))
      }
  
  ANhws11 <- if(length(anlist_hws1[,HWShws1=="1"])==1){
    mean(anlist_hws1[,HWShws1=="1"])} else if (pracma::ndims(anlist_hws1[,HWShws1=="1"])==1){
      mean(anlist_hws1[,HWShws1=="1"])} else {
        mean(colMeans(anlist_hws1[,HWShws1=="1"], na.rm = T))
      }
  
  # Confidence intervals, we first have to sum ANs for all cities for each simulation
  anhws1sim <- apply(ansimlist_hws1, 2:3, mean, na.rm = T) # sum all cities
  ANhws1high <- apply(anhws1sim, 1, quantile, .975, na.rm = T) # Quantiles on sum
  ANhws1low <- apply(anhws1sim, 1, quantile, .025, na.rm = T) # Quantiles on sum
  
  ANhws1; ANhws1low; ANhws1high
  
  anhws1sim0 <- apply(ansimlist_hws1, 2:3, mean, na.rm = T) # sum all cities
  # anhws1sim0 <- anhws1sim0[HWShws1=="0",]         # select periods
  anhws1sim0 <- if(pracma::ndims(anhws1sim0[HWShws1=="0",])==1){
    anhws1sim0[HWShws1=="0",]} else {
      apply(anhws1sim0[HWShws1=="0",], 2, mean, na.rm = T)}  # sum of periods
  ANhws1high0 <- quantile(anhws1sim0,.975, na.rm = T) # Quantiles of sums
  ANhws1low0 <- quantile(anhws1sim0,.025, na.rm = T) # Quantiles of sums
  
  ansimlist_0_hws1[country[r],] <- anhws1sim0 # safe simulations
  
  ANhws10; ANhws1low0; ANhws1high0
  
  anhws1sim1 <- apply(ansimlist_hws1, 2:3, mean, na.rm = T) # sum all cities
  anhws1sim1 <- if (pracma::ndims(anhws1sim1[HWShws1=="1",])==1){
    anhws1sim1[HWShws1=="1",]} else {
      apply(anhws1sim1[HWShws1=="1",], 2, mean, na.rm = T)}  # sum of periods
  ANhws1high1 <- quantile(anhws1sim1,.975, na.rm = T) # Quantiles on sum
  ANhws1low1 <- quantile(anhws1sim1,.025, na.rm = T) # Quantiles on sum
  
  ansimlist_1_hws1[country[r],] <- anhws1sim1 # safe simulations
  
  ANhws11; ANhws1low1; ANhws1high1
  
  
  #----- Now for regions we have to average city-level AFs
  
  # Point estimate
  AFhws1 <- colMeans(aflist_hws1, na.rm = T)
  AFhws10 <- if(pracma::ndims(aflist_hws1[,HWShws1=="0"])==1){mean(aflist_hws1[,HWShws1=="0"])} else {mean(colMeans(aflist_hws1[,HWShws1=="0"], na.rm = T))}
  AFhws11 <- if(pracma::ndims(aflist_hws1[,HWShws1=="1"])==1){mean(aflist_hws1[,HWShws1=="1"])} else {mean(colMeans(aflist_hws1[,HWShws1=="1"], na.rm = T))}
  
  # Confidence intervals, we first have to mean AFs for all cities for each simulation
  afhws1sim <- apply(afsimlist_hws1, 2:3, mean, na.rm = T) # mean all cities
  AFhws1high <- apply(afhws1sim, 1, quantile, .975, na.rm = T) # Quantiles on mean
  AFhws1low <- apply(afhws1sim, 1, quantile, .025, na.rm = T) # Quantiles on mean
  
  AFhws1; AFhws1low; AFhws1high
  
  afhws1sim0 <- apply(afsimlist_hws1, 2:3, mean, na.rm = T) # mean all cities
  afhws1sim0 <- if (pracma::ndims(afhws1sim0[HWShws1=="0",])==1){
    afhws1sim0[HWShws1=="0",]} else {
      apply(afhws1sim0[HWShws1=="0",], 2, mean, na.rm = T)}  # mean of periods
  AFhws1high0 <- quantile(afhws1sim0, .975, na.rm = T) # Quantiles on mean
  AFhws1low0 <- quantile(afhws1sim0, .025, na.rm = T) # Quantiles on mean
  
  afsimlist_0_hws1[country[r],] <- afhws1sim0 # safe simulations
  
  AFhws10; AFhws1low0; AFhws1high0
  
  afhws1sim1 <- apply(afsimlist_hws1, 2:3, mean, na.rm = T) # mean all cities
  afhws1sim1 <- if (pracma::ndims(afhws1sim1[HWShws1=="1",])==1){
    afhws1sim1[HWShws1=="1",]} else {
      apply(afhws1sim1[HWShws1=="1",], 2, mean, na.rm = T)}  # mean of periods
  AFhws1high1 <- quantile(afhws1sim1,.975, na.rm = T) # Quantiles on mean
  AFhws1low1 <- quantile(afhws1sim1,.025, na.rm = T) # Quantiles on mean
  
  afsimlist_1_hws1[country[r],] <- afhws1sim1 # safe simulations
  
  AFhws11; AFhws1low1; AFhws1high1
  
  df2 <- data.frame(cities$Region[1],cities$hwc[1],country[r],as.numeric(sapply(yearlistII, mean)),HWShws1,as.numeric(ANhws1),as.numeric(ANhws1low),as.numeric(ANhws1high),
                    as.numeric(AFhws1),as.numeric(AFhws1low),as.numeric(AFhws1high))
  
  i<-0
  for (i in 1:length(df2)){
    af_table_countries_hws1[radek+1:length(yearlistII),i] <- df2[,i]
  }
  radek <- radek+length(yearlistII)
  
  radek <- radek+1
  af_table_countries_hws1[radek,] <- c(cities$Region[1],cities$hwc[1],country[r],"period1","0",as.numeric(ANhws10),as.numeric(ANhws1low0),as.numeric(ANhws1high0),
                                       as.numeric(AFhws10),as.numeric(AFhws1low0),as.numeric(AFhws1high0))
  radek <- radek+1
  af_table_countries_hws1[radek,] <- c(cities$Region[1],cities$hwc[1],country[r],"period2","1",as.numeric(ANhws11),as.numeric(ANhws1low1),as.numeric(ANhws1high1),
                                       as.numeric(AFhws11),as.numeric(AFhws1low1),as.numeric(AFhws1high1))
}

#===============================================================================
# 7. EXPORT FINAL OUTPUTS (UNCOMMENT TO SAVE)
#===============================================================================
# write.csv(af_table_countries_hws0, "data/af_table_countries_hws0.csv", row.names=FALSE)

#===============================================================================
# END OF SCRIPT
#===============================================================================
