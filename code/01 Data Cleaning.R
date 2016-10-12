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
file_name <- "wvslongdata.rdata"
zip_name <- "rawdata.zip"
search_pattern <- ".rdata"
rawdata <- paste(file_dir, file_name, sep = "")

if (!file.exists(rawdata)) {
  unzip(paste(file_dir, zip_name, sep = ""), exdir = "rawdata")
  search_file <- list.files(file_dir, pattern = search_pattern)
  file.rename(paste(file_dir, search_file, sep = ""), rawdata)
}

# Read data into R
data_orig <- load(rawdata