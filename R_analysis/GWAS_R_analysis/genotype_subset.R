#!/usr/bin/env Rscript

library(data.table)
library(tidyverse)

#working directory
wd <- "/home/namulil/lstm_projects/funestus_llineup/R_analysis"
setwd("/home/namulil/lstm_projects/funestus_llineup/R_analysis")


# Load the genotyping data
gt <- fread('.funestus_llineup/R_analysis/gt_glm.csv', na.strings = '')
#subset gt into 10k rows and save as rds
nrows_subset <- 100000
subset_indices <- seq(1,nrow(gt),nrows_subset)

# iterate through each subset
for (i in 1:length(subset_indices)){
  # calculate the start and end indices for the current subset
  start_index <- subset_indices[i]
  end_index <- min(start_index+nrows_subset - 1, nrow(gt))
  
  #subset df
  subset_gt <- gt[start_index:end_index,]
  filename <- paste0('gt_',i,'.rds')
  #save
  saveRDS(subset_gt, paste('./glm/genotype_subset/', filename, sep=''))
  # print progress
  cat('Subset', i, 'saved as ',filename, '\n')
}
