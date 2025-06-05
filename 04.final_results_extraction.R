
################################################################################
# SCRIPT NAME: 04.stage4_final_analysis_and_reporting.R
#
# DESCRIPTION:
#   This script aggregates results from previous stages, produces summary tables,
#   and generates final figures for reporting. It compares counterfactual and
#   factual scenarios of heat warning systems (HWS) on heat-attributable mortality.
#
# INPUTS:
#   - Outputs from 03.mcc_hews_attributable_mortality_calculation.R:
#       * af_table_countries_hws0, af_table_countries_hws1
#       * ansimlist_0_hws0, ansimlist_1_hws0, ansimlist_0_hws1, ansimlist_1_hws1
#       * afsimlist_0_hws0, afsimlist_1_hws0, afsimlist_0_hws1, afsimlist_1_hws1
#   - City/country metadata: eu_cities, data/eu_countries.csv
#
# OUTPUTS:
#   - Summary tables (country, region, HWS class, EU level)
#   - Plots/Figures for main text and supplementary materials
#
# USAGE:
#   Ensure all required input objects are loaded into the R environment.
#   Update 'path_out' variable for output location as needed.
#
# AUTHORs: Aleš Urban, Veronika Huber, Pierre Masselot, Antonio Gasparrini
# DATE: June 2025
################################################################################

#===============================================================================
# 1. REQUIRED PACKAGES
#===============================================================================

library(scales)
library(dplyr)
library(ggplot2)
library(lubridate)
library(readr)
library(tidyverse)
library(ggh4x)

#===============================================================================
# 2. LOAD REQUIRED DATA
#===============================================================================
# Requires objects from "03.mcc_hews_attributable_mortality_calculation.R"

check_required_objects <- function(required) {
  missing <- required[!required %in% ls(envir = .GlobalEnv)]
  if (length(missing) == 0) {
    message("All required objects are present.")
  } else {
    warning("Missing objects: ", paste(missing, collapse = ", "))
  }
}

required_objects <- c("af_table_countries_hws0","af_table_countries_hws1",
                      "ansimlist_0_hws0", "ansimlist_1_hws0",
                      "ansimlist_0_hws1", "ansimlist_1_hws1",
                      "afsimlist_0_hws0", "afsimlist_1_hws0",
                      "afsimlist_0_hws1", "afsimlist_1_hws1")

check_required_objects(required_objects)

#===============================================================================
# 3. PREPARE AND MERGE TABLES
#===============================================================================

aftable <- data.frame(af_table_countries_hws0, af_table_countries_hws1) %>%
  mutate(across(starts_with("AN"), as.numeric),
         across(starts_with("AF"), as.numeric)) %>%
  dplyr::select(-c(12:16)) %>%
  na.omit()

eucountries <- read.csv2("data/eu_countries.csv", sep = ";")

# Helper: update simulation lists with region/class info
sim_lists <- c("ansimlist_0_hws0", "ansimlist_1_hws0", "ansimlist_0_hws1", "ansimlist_1_hws1",
               "afsimlist_0_hws0", "afsimlist_1_hws0", "afsimlist_0_hws1", "afsimlist_1_hws1")

for (obj in sim_lists) {
  tmp <- as.data.frame(get(obj))
  tmp$class <- eucountries$HWSclass
  tmp$region <- eucountries$Region
  assign(obj, tmp)
}

regions <- c("Northern Europe","Eastern Europe","Western Europe","Southern Europe")
hwclass <- c(1,2,3)

#===============================================================================
# 4. OUTPUT DATA
#===============================================================================

#===============================================================================
# 4.1 COUNTRY-LEVEL SUMMARY
#===============================================================================

af_aftable_countries <- array(NA, c(length(levels(as.factor(eu_cities$countryname))) * 6, 12))
colnames(af_aftable_countries) <- c("region","hwsclass","country","year","hws","AN","ANlow","ANhigh","AF","AFlow","AFhigh","scenario")
row_index <- 1

for (j in hwclass) {
  aftable_sub <- aftable[aftable$hwsclass == j, ]
  country_list <- levels(as.factor(aftable_sub$country))
  
  for (country in country_list) {
    get_af_row <- function(df, year, scenario_label, cols = 1:11) {
      row <- df[df$year == year & df$country == country, cols]
      if (nrow(row) > 0) return(c(unlist(row), scenario_label)) else return(rep(NA, 12))
    }
    af_aftable_countries[row_index, ] <- get_af_row(aftable_sub, "period1", "period1")
    row_index <- row_index + 1
    af_aftable_countries[row_index, ] <- get_af_row(aftable_sub, "period2", "period2_cft")
    row_index <- row_index + 1
    af_aftable_countries[row_index, ] <- get_af_row(aftable_sub, "period2", "period2_fct", c(1:5, 12:17))
    row_index <- row_index + 1
    
    # Index values
    get_row_index <- function(scenario) which(af_aftable_countries[, 3] == country & af_aftable_countries[, 12] == scenario)
    get_af_ratio <- function(val1, val2) (val1 - val2) / val2 * 100
    quant_diff <- function(v1, v2, prob) quantile(v1 - v2, prob, na.rm = TRUE)
    quant_ratio <- function(v1, v2, prob) quantile((v1 - v2) / v2, prob, na.rm = TRUE) * 100
    
    idx_p1 <- get_row_index("period1")
    idx_p2_cft <- get_row_index("period2_cft")
    idx_p2_fct <- get_row_index("period2_fct")
    region_val <- eu_cities$Region[which(eu_cities$countryname == country)][1]
    
    af_aftable_countries[row_index, ] <- c(
      region_val, j, country, "HAF ratio_cft", "-",
      as.numeric(af_aftable_countries[idx_p2_cft, 6]) - as.numeric(af_aftable_countries[idx_p1, 6]),
      quant_diff(ansimlist_1_hws0[country, 1:100], ansimlist_0_hws0[country, 1:100], 0.025),
      quant_diff(ansimlist_1_hws0[country, 1:100], ansimlist_0_hws0[country, 1:100], 0.975),
      get_af_ratio(as.numeric(af_aftable_countries[idx_p2_cft, 9]), as.numeric(af_aftable_countries[idx_p1, 9])),
      quant_ratio(afsimlist_1_hws0[country, 1:100], afsimlist_0_hws0[country, 1:100], 0.025),
      quant_ratio(afsimlist_1_hws0[country, 1:100], afsimlist_0_hws0[country, 1:100], 0.975),
      "-"
    )
    row_index <- row_index + 1
    
    af_aftable_countries[row_index, ] <- c(
      region_val, j, country, "HAF ratio_fct", "-",
      as.numeric(af_aftable_countries[idx_p2_fct, 6]) - as.numeric(af_aftable_countries[idx_p1, 6]),
      quant_diff(ansimlist_1_hws1[country, 1:100], ansimlist_0_hws0[country, 1:100], 0.025),
      quant_diff(ansimlist_1_hws1[country, 1:100], ansimlist_0_hws0[country, 1:100], 0.975),
      get_af_ratio(as.numeric(af_aftable_countries[idx_p2_fct, 9]), as.numeric(af_aftable_countries[idx_p1, 9])),
      quant_ratio(afsimlist_1_hws1[country, 1:100], afsimlist_0_hws0[country, 1:100], 0.025),
      quant_ratio(afsimlist_1_hws1[country, 1:100], afsimlist_0_hws0[country, 1:100], 0.975),
      "-"
    )
    row_index <- row_index + 1
    
    af_aftable_countries[row_index, ] <- c(
      region_val, j, country, "HAF ratio_cft-fct", "-",
      as.numeric(af_aftable_countries[idx_p2_fct, 6]) - as.numeric(af_aftable_countries[idx_p2_cft, 6]),
      quant_diff(ansimlist_1_hws1[country, 1:100], ansimlist_1_hws0[country, 1:100], 0.025),
      quant_diff(ansimlist_1_hws1[country, 1:100], ansimlist_1_hws0[country, 1:100], 0.975),
      get_af_ratio(as.numeric(af_aftable_countries[idx_p2_fct, 9]), as.numeric(af_aftable_countries[idx_p2_cft, 9])),
      quant_ratio(afsimlist_1_hws1[country, 1:100], afsimlist_1_hws0[country, 1:100], 0.025),
      quant_ratio(afsimlist_1_hws1[country, 1:100], afsimlist_1_hws0[country, 1:100], 0.975),
      "-"
    )
    row_index <- row_index + 1
  }
}

# Final output for inspection
View(af_aftable_countries)

#===============================================================================
# 4.2 REGION-LEVEL SUMMARY
#===============================================================================

regions <- c("Northern Europe","Eastern Europe","Western Europe","Southern Europe")
af_aftable_regions <- array(NA,c(length(regions)*6,8))
colnames(af_aftable_regions) <- c("region","AN","ANlow","ANhigh","AF","AFlow","AFhigh","scenario")
row_index <- 1

i <- 1

for (i in 1:length(regions)){
  
  af_aftable_regions[row_index,]<-c(regions[i],colMeans(aftable[which(aftable$year=="period1" & aftable$region==regions[i]),c(6:11)]),
                                "period1")
  
  row_index <- row_index+1
  
  af_aftable_regions[row_index,]<-(c(regions[i],colMeans(aftable[which(aftable$year=="period2" & aftable$region==regions[i]),c(6:11)]),
                                 "period2_cft"))
  
  row_index <- row_index+1
  
  af_aftable_regions[row_index,]<-(c(regions[i],colMeans(aftable[which(aftable$year=="period2" & aftable$region==regions[i]),c(12:17)]),
                                 "period2_fct"))
  row_index <- row_index+1
  
  af_aftable_regions[row_index,]=c(regions[i],
                               as.numeric(af_aftable_regions[,2][which(af_aftable_regions[,1]==regions[i] & af_aftable_regions[,8]=="period2_cft")])-as.numeric(af_aftable_regions[,2][which(af_aftable_regions[,1]==regions[i] & af_aftable_regions[,8]=="period1")]),
                               "-","-",
                               ((as.numeric(af_aftable_regions[,5][which(af_aftable_regions[,1]==regions[i] & af_aftable_regions[,8]=="period2_cft")])-as.numeric(af_aftable_regions[,5][which(af_aftable_regions[,1]==regions[i] & af_aftable_regions[,8]=="period1")]))/
                                  as.numeric(af_aftable_regions[,5][which(af_aftable_regions[,1]==regions[i] & af_aftable_regions[,8]=="period1")])),
                               "-","-","HAF ratio_cft")
  row_index <- row_index+1
  
  af_aftable_regions[row_index,]=c(regions[i],
                               as.numeric(af_aftable_regions[,2][which(af_aftable_regions[,1]==regions[i] & af_aftable_regions[,8]=="period2_fct")])-as.numeric(af_aftable_regions[,2][which(af_aftable_regions[,1]==regions[i] & af_aftable_regions[,8]=="period1")]),
                               "-","-",
                               ((as.numeric(af_aftable_regions[,5][which(af_aftable_regions[,1]==regions[i] & af_aftable_regions[,8]=="period2_fct")])-as.numeric(af_aftable_regions[,5][which(af_aftable_regions[,1]==regions[i] & af_aftable_regions[,8]=="period1")]))/
                                  as.numeric(af_aftable_regions[,5][which(af_aftable_regions[,1]==regions[i] & af_aftable_regions[,8]=="period1")])),
                               "-","-","HAF ratio_fct")
  row_index <- row_index+1
  
  af_aftable_regions[row_index,]=c(regions[i],
                               as.numeric(af_aftable_regions[,2][which(af_aftable_regions[,1]==regions[i] & af_aftable_regions[,8]=="period2_fct")])-as.numeric(af_aftable_regions[,2][which(af_aftable_regions[,1]==regions[i] & af_aftable_regions[,8]=="period2_cft")]),
                               quantile(colMeans(ansimlist_1_hws1[which(ansimlist_1_hws1$region==regions[i]),c(1:100)])-colMeans(ansimlist_1_hws0[which(ansimlist_1_hws0$region==regions[i]),c(1:100)]),.025, na.rm = T),
                               quantile(colMeans(ansimlist_1_hws1[which(ansimlist_1_hws1$region==regions[i]),c(1:100)])-colMeans(ansimlist_1_hws0[which(ansimlist_1_hws0$region==regions[i]),c(1:100)]),.975, na.rm = T),
                               
                               (as.numeric(af_aftable_regions[,5][which(af_aftable_regions[,1]==regions[i] & af_aftable_regions[,8]=="period2_fct")])-as.numeric(af_aftable_regions[,5][which(af_aftable_regions[,1]==regions[i] & af_aftable_regions[,8]=="period2_cft")]))/
                                 as.numeric(af_aftable_regions[,5][which(af_aftable_regions[,1]==regions[i] & af_aftable_regions[,8]=="period2_cft")])*100,
                               quantile((colMeans(afsimlist_1_hws1[which(afsimlist_1_hws1$region==regions[i]),c(1:100)])-colMeans(afsimlist_1_hws0[which(afsimlist_1_hws0$region==regions[i]),c(1:100)]))/
                                          colMeans(afsimlist_1_hws0[which(afsimlist_1_hws0$region==regions[i]),c(1:100)]),.025, na.rm = T)*100,
                               quantile((colMeans(afsimlist_1_hws1[which(afsimlist_1_hws1$region==regions[i]),c(1:100)])-colMeans(afsimlist_1_hws0[which(afsimlist_1_hws0$region==regions[i]),c(1:100)]))/
                                          colMeans(afsimlist_1_hws0[which(afsimlist_1_hws0$region==regions[i]),c(1:100)]),.975, na.rm = T)*100,
                               "HAF ratio_cft-fct")
  row_index <- row_index+1
  
  
  
}

View(af_aftable_regions)

#===============================================================================
# 4.3 HPPclass SUMMARY
#===============================================================================

hwclass <- c(1,2,3)
af_aftable_hwsclass <- array(NA,c(length(hwclass)*6,8))
colnames(af_aftable_hwsclass) <- c("class","AN","ANlow","ANhigh","AF","AFlow","AFhigh","scenario")
row_index <- 1

i <- 1

for (i in 1:length(hwclass)){
  
  af_aftable_hwsclass[row_index,]<-c(hwclass[i],colMeans(aftable[which(aftable$year=="period1" & aftable$hwsclass==hwclass[i]),c(6:11)]),
                                 "period1")
  
  row_index <- row_index+1
  
  af_aftable_hwsclass[row_index,]<-(c(hwclass[i],colMeans(aftable[which(aftable$year=="period2" & aftable$hwsclass==hwclass[i]),c(6:11)]),
                                  "period2_cft"))
  
  row_index <- row_index+1
  
  af_aftable_hwsclass[row_index,]<-(c(hwclass[i],colMeans(aftable[which(aftable$year=="period2" & aftable$hwsclass==hwclass[i]),c(12:17)]),
                                  "period2_fct"))
  row_index <- row_index+1
  
  af_aftable_hwsclass[row_index,]=c(hwclass[i],
                                as.numeric(af_aftable_hwsclass[,2][which(af_aftable_hwsclass[,1]==hwclass[i] & af_aftable_hwsclass[,8]=="period2_cft")])-as.numeric(af_aftable_hwsclass[,2][which(af_aftable_hwsclass[,1]==hwclass[i] & af_aftable_hwsclass[,8]=="period1")]),
                                "-","-",
                                ((as.numeric(af_aftable_hwsclass[,5][which(af_aftable_hwsclass[,1]==hwclass[i] & af_aftable_hwsclass[,8]=="period2_cft")])-as.numeric(af_aftable_hwsclass[,5][which(af_aftable_hwsclass[,1]==hwclass[i] & af_aftable_hwsclass[,8]=="period1")]))/
                                   as.numeric(af_aftable_hwsclass[,5][which(af_aftable_hwsclass[,1]==hwclass[i] & af_aftable_hwsclass[,8]=="period1")])),
                                "-","-","HAF ratio_cft")
  row_index <- row_index+1
  
  af_aftable_hwsclass[row_index,]=c(hwclass[i],
                                as.numeric(af_aftable_hwsclass[,2][which(af_aftable_hwsclass[,1]==hwclass[i] & af_aftable_hwsclass[,8]=="period2_fct")])-as.numeric(af_aftable_hwsclass[,2][which(af_aftable_hwsclass[,1]==hwclass[i] & af_aftable_hwsclass[,8]=="period1")]),
                                "-","-",
                                ((as.numeric(af_aftable_hwsclass[,5][which(af_aftable_hwsclass[,1]==hwclass[i] & af_aftable_hwsclass[,8]=="period2_fct")])-as.numeric(af_aftable_hwsclass[,5][which(af_aftable_hwsclass[,1]==hwclass[i] & af_aftable_hwsclass[,8]=="period1")]))/
                                   as.numeric(af_aftable_hwsclass[,5][which(af_aftable_hwsclass[,1]==hwclass[i] & af_aftable_hwsclass[,8]=="period1")])),
                                "-","-","HAF ratio_fct")
  row_index <- row_index+1
  
  af_aftable_hwsclass[row_index,]=c(hwclass[i],
                                as.numeric(af_aftable_hwsclass[,2][which(af_aftable_hwsclass[,1]==hwclass[i] & af_aftable_hwsclass[,8]=="period2_fct")])-as.numeric(af_aftable_hwsclass[,2][which(af_aftable_hwsclass[,1]==hwclass[i] & af_aftable_hwsclass[,8]=="period2_cft")]),
                                quantile(colMeans(ansimlist_1_hws1[which(ansimlist_1_hws1$class==hwclass[i]),c(1:100)])-colMeans(ansimlist_1_hws0[which(ansimlist_1_hws0$class==hwclass[i]),c(1:100)]),.025, na.rm = T),
                                quantile(colMeans(ansimlist_1_hws1[which(ansimlist_1_hws1$class==hwclass[i]),c(1:100)])-colMeans(ansimlist_1_hws0[which(ansimlist_1_hws0$class==hwclass[i]),c(1:100)]),.975, na.rm = T),
                                
                                (as.numeric(af_aftable_hwsclass[,5][which(af_aftable_hwsclass[,1]==hwclass[i] & af_aftable_hwsclass[,8]=="period2_fct")])-as.numeric(af_aftable_hwsclass[,5][which(af_aftable_hwsclass[,1]==hwclass[i] & af_aftable_hwsclass[,8]=="period2_cft")]))/
                                  as.numeric(af_aftable_hwsclass[,5][which(af_aftable_hwsclass[,1]==hwclass[i] & af_aftable_hwsclass[,8]=="period2_cft")])*100,
                                quantile((colMeans(afsimlist_1_hws1[which(afsimlist_1_hws1$class==hwclass[i]),c(1:100)])-colMeans(afsimlist_1_hws0[which(afsimlist_1_hws0$class==hwclass[i]),c(1:100)]))/
                                           colMeans(afsimlist_1_hws0[which(afsimlist_1_hws0$class==hwclass[i]),c(1:100)]),.025, na.rm = T)*100,
                                quantile((colMeans(afsimlist_1_hws1[which(afsimlist_1_hws1$class==hwclass[i]),c(1:100)])-colMeans(afsimlist_1_hws0[which(afsimlist_1_hws0$class==hwclass[i]),c(1:100)]))/
                                           colMeans(afsimlist_1_hws0[which(afsimlist_1_hws0$class==hwclass[i]),c(1:100)]),.975, na.rm = T)*100,
                                "HAF ratio_cft-fct")
  row_index <- row_index+1
  
  
  
}

View(af_aftable_hwsclass)

#===============================================================================
# 4.4 EU-level SUMMARY
#===============================================================================

# Preallocate output array to store EU-level AF-related values
af_aftable_EU <- array(NA,c(6,8))
colnames(af_aftable_EU) <- c("region","AN","ANlow","ANhigh","AF","AFlow","AFhigh","scenario")
row_index <- 1

# i <- 1
# for (i in 1:length(hwclass)){

af_aftable_EU[row_index,]<-c("EU",colMeans(aftable[which(aftable$year=="period1"),c(6:11)]),
                         "period1")

row_index <- row_index+1

af_aftable_EU[row_index,]<-(c("EU",colMeans(aftable[which(aftable$year=="period2"),c(6:11)]),
                          "period2_cft"))

row_index <- row_index+1

af_aftable_EU[row_index,]<-(c("EU",colMeans(aftable[which(aftable$year=="period2"),c(12:17)]),
                          "period2_fct"))
row_index <- row_index+1

af_aftable_EU[row_index,]=c("EU",
                        as.numeric(af_aftable_EU[,2][which(af_aftable_EU[,8]=="period2_cft")])-as.numeric(af_aftable_EU[,2][which(af_aftable_EU[,8]=="period1")]),
                        "-","-",
                        ((as.numeric(af_aftable_EU[,5][which(af_aftable_EU[,8]=="period2_cft")])-as.numeric(af_aftable_EU[,5][which(af_aftable_EU[,8]=="period1")]))/
                           as.numeric(af_aftable_EU[,5][which(af_aftable_EU[,8]=="period1")])),
                        "-","-","HAF ratio_cft")
row_index <- row_index+1

af_aftable_EU[row_index,]=c("EU",
                        as.numeric(af_aftable_EU[,2][which(af_aftable_EU[,8]=="period2_fct")])-as.numeric(af_aftable_EU[,2][which(af_aftable_EU[,8]=="period1")]),
                        "-","-",
                        ((as.numeric(af_aftable_EU[,5][which(af_aftable_EU[,8]=="period2_fct")])-as.numeric(af_aftable_EU[,5][which(af_aftable_EU[,8]=="period1")]))/
                           as.numeric(af_aftable_EU[,5][which(af_aftable_EU[,8]=="period1")])),
                        "-","-","HAF ratio_fct")
row_index <- row_index+1

af_aftable_EU[row_index,]=c("EU",
                        as.numeric(af_aftable_EU[,2][which(af_aftable_EU[,8]=="period2_fct")])-as.numeric(af_aftable_EU[,2][which(af_aftable_EU[,8]=="period2_cft")]),
                        quantile(colMeans(ansimlist_1_hws1[,c(1:100)])-colMeans(ansimlist_1_hws0[,c(1:100)]),.025, na.rm = T),
                        quantile(colMeans(ansimlist_1_hws1[,c(1:100)])-colMeans(ansimlist_1_hws0[,c(1:100)]),.975, na.rm = T),
                        
                        (as.numeric(af_aftable_EU[,5][which(af_aftable_EU[,8]=="period2_fct")])-as.numeric(af_aftable_EU[,5][which(af_aftable_EU[,8]=="period2_cft")]))/
                          as.numeric(af_aftable_EU[,5][which(af_aftable_EU[,8]=="period2_cft")])*100,
                        quantile((colMeans(afsimlist_1_hws1[,c(1:100)])-colMeans(afsimlist_1_hws0[,c(1:100)]))/
                                   colMeans(afsimlist_1_hws0[,c(1:100)]),.025, na.rm = T)*100,
                        quantile((colMeans(afsimlist_1_hws1[,c(1:100)])-colMeans(afsimlist_1_hws0[,c(1:100)]))/
                                   colMeans(afsimlist_1_hws0[,c(1:100)]),.975, na.rm = T)*100,
                        "HAF ratio_cft-fct")
row_index <- row_index+1



# }

View(af_aftable_EU)

#===============================================================================
# FIGURE 4
#===============================================================================

# plot AF

af_aftable_countries <- af_aftable_countries %>%
  as.data.frame() %>%
  mutate(across(c(AN, ANlow, ANhigh, AF, AFlow, AFhigh), as.numeric))


# af_aftable_sub <- subset(af_aftable_regions,af_aftable_regions$country=="-")
af_aftable_sub <- subset(af_aftable_countries,af_aftable_countries$scenario=="period1" | af_aftable_countries$scenario=="period2_cft" | af_aftable_countries$scenario=="period2_fct")
af_aftable_sub$region <- factor(af_aftable_sub$region,levels = c("Northern Europe","Eastern Europe",  
                                                                 "Western Europe", "Southern Europe"))
af_aftable_sub$hwsclass <- factor(af_aftable_sub$hwsclass,levels = c("1","2","3"))
af_aftable_sub$country <- as.factor(af_aftable_sub$country)

my_colors <- RColorBrewer::brewer.pal(4, "Blues")[4:2]

AF <- ggplot(af_aftable_sub,aes(x=country, y=AF,ymin=AFlow, ymax=AFhigh, fill=scenario))+
  geom_col(position = "dodge")+
  geom_errorbar(aes(ymin = AFlow, ymax = AFhigh), 
                position = position_dodge(0.9), width = .3)+
  scale_fill_manual(name = "scenario", labels = c("pre-implementation","counterfactual","factual"),values = my_colors)+
  labs(y="HAF [%]",x="Region")+
  ylim(-0.5, 4)+
  facet_wrap(.~region,ncol=2,
             scales = "free")+
  theme_bw()+
  theme(
    panel.grid.minor.y = element_blank(),
    strip.text = element_text(
      size = rel(1.5), face = "bold"),
    plot.title = element_text(size = rel(1.2),face = "bold"),
    axis.title.y = element_text(vjust= 1.0,size = rel(1.5)),
    axis.text.y = element_text(size = rel(1.2)),
    axis.title.x = element_text(vjust= 1.0,size = rel(1.5)),
    axis.text.x = element_text(size = rel(1.2)),
    legend.text = element_text(size = rel(1.5)),
    legend.title = element_text(size = rel(1.5),face="bold"),
    legend.position = "top"
    
  )
print(AF)

# pdf(paste(path_out,"HAF_regions_scenarios.pdf",sep=""),width=10,height=8)
# print(AF)
# dev.off()

#===============================================================================
# FIGURE 5
#===============================================================================

# PLOTS AF ratio period 2 fct vs period 2 cft

af_aftable_sub <- subset(af_aftable_countries,af_aftable_countries$year=="HAF ratio_cft-fct")
af_aftable_sub$region <- factor(af_aftable_sub$region,levels = c("Northern Europe","Eastern Europe",  
                                                                 "Western Europe", "Southern Europe","All"))
af_aftable_sub$hwsclass <- factor(af_aftable_sub$hwsclass,levels = c("1","2","3","-"))
af_aftable_sub$country <- as.factor(af_aftable_sub$country)

af_aftable_sub <- rbind(af_aftable_sub, c("All","-","Overall", af_aftable_EU[6,8],"-",af_aftable_EU[6,2:7],"-"))

af_aftable_sub <- af_aftable_sub %>%
  mutate(across(c(AN, ANlow, ANhigh, AF, AFlow, AFhigh), as.numeric))


library(ggh4x)

strip <- strip_themed(background_y = elem_list_rect(fill = c("cornflowerblue","orange","chartreuse3","coral3","grey85")))

AFratio <- ggplot(af_aftable_sub,aes(x=country, y=AF, ymin=AFlow,ymax=AFhigh, colour = region))+
  geom_pointrange(size=1,linewidth=0.9,position=position_dodge(0.9),aes(shape=region))+
  # geom_point(size=4,position=position_dodge(0.9),aes(shape=region),colour = "black")+
  scale_shape_manual(values = c(15,19,17,18,8),name="Region") +
  # scale_fill_manual(values = c("cornflowerblue","orange","chartreuse3","coral3"),name="Region")+
  scale_colour_manual(values = c("cornflowerblue","orange","chartreuse3","coral3","black"),name="Region")+
  labs(y="HAF change [%]",x="")+ 
  ylim(-75, 35)+
  scale_x_discrete(limits = rev)+
  geom_hline(yintercept=0)+
  coord_flip()+
  facet_grid2(region~., space="free",      #labeller = to_string,
              scales = "free",strip=strip)+
  theme_bw()+
  theme(
    panel.grid.minor.y = element_blank(),
    strip.text = element_text(size = rel(1.2), face = "bold"),
    # strip.background = element_rect(fill=c("cornflowerblue","orange","chartreuse3","coral3")),
    plot.title = element_text(size = rel(1.6),face = "bold"),
    axis.title.y = element_text(vjust= 1.0,size = rel(1.6)),
    axis.text.y = element_text(size = rel(1.6)),
    axis.title.x = element_text(vjust= 1.0,size = rel(1.6)),
    axis.text.x = element_text(size = rel(1.6)),
    legend.position="none",
    legend.text = element_text(size=rel(1.6)),
    legend.title = element_text(size=rel(1.6),face = "bold")
  )

print(AFratio)

#===============================================================================
# FIGURES S4 
#===============================================================================

par(mfrow=c(3,3))

aftable$country <- as.factor(aftable$country)
country <- levels(aftable$country)

r <- 1
for (r in 1:length(country)){
  
  aftable_trends <- aftable[which(aftable$country==country[r]),]
  aftable_trends <- aftable_trends[-which(aftable_trends$year=="period1" | aftable_trends$year=="period2"),]
  
  
  maxyear <- max(aftable_trends$year[which(aftable_trends$country==country[r] & aftable_trends$hws=="0")])
  minyear <- min(aftable_trends$year[which(aftable_trends$country==country[r] & aftable_trends$hws=="1")])
  
  
  plot(1991:2018, seq(1991:2018), type="n", ylim=c(-0.5,6),axes=F,
       ylab="HAF [%]", xlab="subperiod", bty="l", las=1, mgp=c(2.5,1,0), cex.axis=0.8,
       main=country[r])
  
  axis(1,at=seq(from=1991,to=2018,by=3),cex.axis=1.0)
  axis(2,at=seq(from=0,to=8,by=2),cex.axis=1.0)
  
  # counterfactual
  
  points(as.numeric(aftable_trends$year[c(which(aftable_trends$year==maxyear),which(aftable_trends$year>maxyear))]), c(aftable_trends$AF[which(aftable_trends$year==maxyear)],aftable_trends$AF[which(aftable_trends$year>maxyear)]),
         type="o", col=3, pch=19,lty="dashed")
  arrows(as.numeric(aftable_trends$year[c(which(aftable_trends$year==maxyear),which(aftable_trends$year>maxyear))]), c(aftable_trends$AFlow[which(aftable_trends$year==maxyear)],aftable_trends$AFlow[which(aftable_trends$year>maxyear)]),
         as.numeric(aftable_trends$year[c(which(aftable_trends$year==maxyear),which(aftable_trends$year>maxyear))]), c(aftable_trends$AFhigh[which(aftable_trends$year==maxyear)],aftable_trends$AFhigh[which(aftable_trends$year>maxyear)]),
         col=alpha(3,0.5),
         code=3, angle=90, length=0.05, lwd=1)
  
  # factual
  
  points(as.numeric(aftable_trends$year[c(which(aftable_trends$year<maxyear),which(aftable_trends$year==maxyear))]), c(aftable_trends$AF[which(aftable_trends$year<maxyear)],aftable_trends$AF[which(aftable_trends$year==maxyear)]),
         type="o", col=6, pch=19,lty="solid")
  arrows(as.numeric(aftable_trends$year[c(which(aftable_trends$year<maxyear),which(aftable_trends$year==maxyear))]), c(aftable_trends$AFlow[which(aftable_trends$year<maxyear)],aftable_trends$AFlow[which(aftable_trends$year==maxyear)]),
         as.numeric(aftable_trends$year[c(which(aftable_trends$year<maxyear),which(aftable_trends$year==maxyear))]), c(aftable_trends$AFhigh[which(aftable_trends$year<maxyear)],aftable_trends$AFhigh[which(aftable_trends$year==maxyear)]),
         col=alpha(6,0.5),
         code=3, angle=90, length=0.05, lwd=1)
  
  points(as.numeric(aftable_trends$year[c(which(aftable_trends$year==maxyear),which(aftable_trends$year>maxyear))]), c(aftable_trends$AF[which(aftable_trends$year==maxyear)],aftable_trends$AF.1[which(aftable_trends$year>maxyear)]),
         type="o", col=6, pch=19,lty="solid")
  arrows(as.numeric(aftable_trends$year[c(which(aftable_trends$year==maxyear),which(aftable_trends$year>maxyear))]), c(aftable_trends$AFlow[which(aftable_trends$year==maxyear)],aftable_trends$AFlow.1[which(aftable_trends$year>maxyear)]),
         as.numeric(aftable_trends$year[c(which(aftable_trends$year==maxyear),which(aftable_trends$year>maxyear))]), c(aftable_trends$AFhigh[which(aftable_trends$year==maxyear)],aftable_trends$AFhigh.1[which(aftable_trends$year>maxyear)]),
         col=alpha(6,0.5),
         code=3, angle=90, length=0.05, lwd=1)
  
  
  
  abline(h=0,col="grey75")
  legend("topright", c("factual","counterfactual"), pch=19, col=c(6,3), bty="n",cex=0.9)
  
}
