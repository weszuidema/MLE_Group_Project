#' Title: WVS Data Cleanup
#' Author: Wesley Zuidema and Ellen Ahlness
#' Project: MLE Group Project
################################################################################

# Load Library packages. Use install.packages() for any package not already installed
library(dplyr)
library(tidyr)



################################################################################
# Read and extract data into R
################################################################################

# Checks to see if an rdata file already exists. If not, extracts data from zipped file.

# These are just paths that could be changed if necessary
file_dir <- "rawdata/"
exdir_dir <- "rawdata"
file_name <- "wvslongdata.rdata"
zip_name <- "rawdata.zip"
search_pattern <- ".rdata"
rawdata <- paste(file_dir, file_name, sep = "")

output_dir <- "data/"
output_name <- "wvs_data.rdata"
output_full <- paste(output_dir, output_name, sep = "")

if (!file.exists(rawdata)) {
  unzip(paste(file_dir, zip_name, sep = ""), exdir = exdir_dir)
  search_file <- list.files(file_dir, pattern = search_pattern)
  file.rename(paste(file_dir, search_file, sep = ""), rawdata)
}

# Read data into R
load(rawdata)


# Save whole dataset, as well as countries of interest
if (!file.exists(output_full)) {
  
save(WVS_Longitudinal_1981_2014_R_v2015_04_18, file = output_full)
  
}

USA <- filter(WVS_Longitudinal_1981_2014_R_v2015_04_18, grepl("840", S003))
UK <- filter(WVS_Longitudinal_1981_2014_R_v2015_04_18, grepl("826", S003))

save(USA, UK, file = "data/wvs_states.rdata")
