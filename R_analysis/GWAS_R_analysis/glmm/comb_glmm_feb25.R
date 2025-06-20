#!/usr/bin/env Rscript
## bash code for runnng script(Rscript ./comb_glmm_feb25.R > comb_glmm.log 2>&1)
library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(glmmTMB)
library(grid)
library(gridExtra)
library(gtable)
#library(tidyverse)
library(future.apply)

#working directory
wd <- '/home/namulil/lstm_projects/R_analysis'

setwd(wd)

# Load the genotyping data
# fetch genotype subset  files 
dir.genotype_subset <-  './glm/genotype_subset'
genotype.files <- list.files(dir.genotype_subset)
genotype.files <- genotype.files [grepl('.*rds', genotype.files)]
dir.glm.result <- './glm /glmm_comb_feb25' 


#read subset genotype files  and run gwas
for (c in genotype.files){
  file.name <- sub('.rds', '', c)
  setwd(dir.genotype_subset)
  gt  <- readRDS(c)
  setwd(wd)
  # loading metadata by self 
  metadata <- fread ('/home/namulil/lstm_projects/R_analysis/metadata.csv', na.strings = '')# load your metadata
  
  #Prep data-Select only pbo samples in round 1 and 5 or any net and rounds you want
  sample.names <- colnames(gt)[3:ncol(gt)]
  #metadata_pbo<- metadata[sample_id	%in% sample.names & sex_call=='F' & RND %in% c(1,5)&LLIN.actual==2]
  metadata_comb<- metadata[sample_id %in% sample.names & sex_call=='F' & RND %in% c(1,5)]
  metadata <- metadata_comb %>% distinct(sample_id, .keep_all = TRUE)#Remove duplicated row in meta
  gt.id <- gt$`Unnamed: 0`
  #gt <- gt[,1:2 := NULL]
  gt <- gt %>% select(metadata$sample_id)
  setcolorder(gt, metadata$sample_id)
  
  
  # We will turn the gt to haplotypes
  gen2hap <- function(n){
    outvec <- rep(1,2)
    outvec[seq_len(2-n)] <- 0
    outvec
  }
  
  gen2hap.vector <- function(N){
    unlist(lapply(N, gen2hap))
  }
  
  #run in parallel
  # multi cores processing
  n.cores <- 30
  # Set the number of cores that future_apply will use for all its operations
  plan(tweak(multisession, workers = n.cores))
  
  #allow variables larger than default 500mb
  options(future.globals.maxSize=+Inf)
  #genotype to haplotype
  start.time <- Sys.time()
  ug.hap <- future_apply(gt, 1, gen2hap.vector) %>% data.frame()
  
  end.time <- Sys.time()
  hap.conversion.time <- end.time-start.time
  cat('Haplotype conversion', c, 'done after',hap.conversion.time, '\n')
  
  colnames(ug.hap) <- gt.id
  
  #double metadata because Turning genotype to haps duplicated each snp position
  
  model.meta<- metadata[,.(sample_id= rep(sample_id, each =2),
                           RND = rep(RND, each = 2),
                           HSD = rep(HSD, each = 2),
                           EA = rep(EA, each = 2),
                           HHID.rnd = rep(HHID.rnd, each = 2),
                           Wave = rep(Wave, each = 2),
                           LLIN.actual = rep(LLIN.actual, each = 2),
                           Location = rep(admin1_iso, each = 2),
                           wgs.sample.id = rep(sample_id, each = 2)
  )]
  
  #Running the model
  rnd<- as.factor(model.meta$RND)
  hsd<- as.factor(model.meta$HSD)
  ea<- as.factor(model.meta$EA)
  hhid <- as.factor(model.meta$HHID.rnd)
  wave <- as.factor(model.meta$Wave)
  location<- as.factor(model.meta$Location)
  llin <- as.factor(model.meta$LLIN.actual)
  
  #Extract p values in coefficients
  get.term.details <- function(model, test, model.term, test.term){
    p <- test[['Pr(>Chi)']][which(attributes(test)$row.names == test.term)]
    coefficient <- coef(model)[model.term]
    c(P = signif(p, 2), coef = signif(coefficient, 2))
  }
  
  #glm model function
run.model <- function(markers, ug.hap){
  
  model <- glmmTMB(get(markers) ~ rnd+ llin +(1|hsd), 
                   data = ug.hap, family = 'binomial')
  test1 <- drop1(model, test = 'Chisq')
  #print(test1)
  
  test.summary = get.term.details(model,test1,'rnd','rnd')
  test.summary
}
  markers <- gt.id
  
  # multi cores processing
  n.cores <- 30
  ##plan(tweak(multisession, workers = n.cores))
  #allow variables larger than default 500mb
  options(future.globals.maxSize=+Inf)
  
  #Run model
  start.time <- Sys.time()
  model.test <- do.call(rbind,future_lapply(markers, run.model,ug.hap))
  model.test <- as.data.frame(model.test) 
  model.test$snp.id <- gt.id
  saveRDS(model.test,paste(dir.glm.result,'/',file.name,'_glm.rds',sep = ''))
  cat('Subset comb_glmm', c, 'saved as ',file.name, '_comb_glmmfeb.rds','\n')
  end.time <- Sys.time()
  model_100k.time <- end.time-start.time
  cat('Glmm of', c, 'started at', start.time, 'and finished after',end.time, 'taking',model_100k.time, 'to complete', '\n')
}


