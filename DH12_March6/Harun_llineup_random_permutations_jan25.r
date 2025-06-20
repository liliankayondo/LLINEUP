library(data.table)
library(stringr)
library(RColorBrewer)
library(magrittr)

setwd('~/lstm_projects/funestus_llineup/notebooks/DH12_permutations')

# complete llineup_seq data
## loading csv (df_samples) created from combing af1000 metadata to the funestus llineup extra meta data
meta <- fread('permutation_df_samples.csv',na.strings = '') ## fread used to read in data

# before screening sample_id column, I want to have a look at the columns in my data
names(meta)

meta <- unique(meta, by='sample_id') ## uniques used to filter data

length(meta$LLIN.actual)  # checking length of my metadata

#Remove first degree related samples fro KING relatedness analysis
## analysis from Sanjay at start showed no relatedness so commenting next codes out
#samples.remove <- c('VBS50531-6645STDY11194268','VBS50533-6645STDY11194270','VBS50528-6645STDY11194265')
#meta<- meta[!sample_id %in% samples.remove]


# since i was getting error when trying to relabel llin.actual column, I am going to first try and drop row with empty entries
# Remove rows where LLIN.actual is empty or NA
meta <- meta[!(meta$LLIN.actual == "" | is.na(meta$LLIN.actual)), ]
# checking size of my meta after the drop
length(meta$LLIN.actual)
## looks like this was successful have 1045 from 1150


# tring to change llin.actual entries from numerics to categories 
meta$LLIN.actual <- factor(meta$LLIN.actual, 
                         levels = c(1, 2), 
                         labels = c("PBO", "NonPBO"))


## need to recode my rounds from baseline to pre and 25 months to post and rest as intermediates in my data
meta$Control_phase <- factor(meta$RND, 
                             levels = c(1,2,3,4,5), 
                             labels = c("Pre", "Intermediate", "Intermediate", "Intermediate", "Post"))



#meta$LLIN.actual <- gsub("Non-PBO", "NonPBO", meta$LLIN.actual) commeting out Haruns code cuz my LLIN.actual colum has numerics not factors 

# Have a column to indicate population (location/insecticide)
meta$population <- with(meta, paste(LLIN.actual, sep = '_')) #randomise by round at each intervention seperately
#meta$population <- with(meta, paste('uganda')) #randomise by round the whole population

# Reorder and remove columns
column.order <- c('sample_id', 'location', 'Control_phase','population', 'admin1_iso')
meta <- meta[, ..column.order]

#select round 
meta <- meta[Control_phase %in% c('Pre','Post'), ] # ensure in your metadata you have a column where you have grouped your intervention to pre, intermediate and post so that you can select pre and post

# A function that will return a list of length k, where each k is a random shuffling of the input vector
shuffle <- function(x, k){
  data.table(replicate(k, sample(x, length(x), replace = F)))
}
# Create 1000 or more randomisations of the phenotype labels, stratified by population, and add them to the control
# table
set.seed(42)
num.randomisations <- 1000
replicate.names <- paste('r', str_pad(1:num.randomisations, nchar(as.character(num.randomisations)), pad = 0), sep = '')

meta[, (replicate.names) := shuffle(Control_phase, num.randomisations), by = population]

# Write the table to file
fwrite(meta, '~/lstm_projects/funestus_llineup/notebooks/DH12_permutations/net_round_randomisations_pop.csv', sep = '\t', quote = F)

