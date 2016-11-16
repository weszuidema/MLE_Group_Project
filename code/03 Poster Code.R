#' Code for generating analysis for poster project
#' Wesley Zuidema and Ellen Ahlness
#' MLE Poster Project

# Load library packages
library(ggplot2)
library(dplyr)


# Load data
load("data/wvs_data.rdata")
wvs_data <- WVS_Longitudinal_1981_2014_R_v2015_04_18
rm(WVS_Longitudinal_1981_2014_R_v2015_04_18)

# Rename variables we are using to make them easier to read
wvs_data$immigrationpref <- wvs_data$E143
wvs_data$changeacceptance <- wvs_data$E047
wvs_data$age <- wvs_data$X003
wvs_data$incomelevel <- wvs_data$X047
wvs_data$female <- wvs_data$X001
wvs_data$maritalstatus <- wvs_data$X007

wvs_data$wave <- as.factor(wvs_data$S002)
wvs_data$country <- as.factor(wvs_data$S003)

# Create binary variable to indicate support of far right party
# Have data for 20 - Andorra, 36 - Australia, 348 - Hungary, 360 - Indonesia, 
# 484 - Mexico, 840 - USA, and 858 - Uraguay
wvs_data$partypref <- wvs_data$E256
wvs_filter <- filter(wvs_data, country == 348, partypref > 7)
unique(wvs_filter$partypref)
ggplot(wvs_filter) + geom_jitter(aes(x = country, y = partypref, col = country))

# Recode data to remove non-responses


# Build model and estimate with least squares first for reference
model_ls <-   ~ immigrationpref + changeacceptance + age + incomelevel + female + maritalstatus
