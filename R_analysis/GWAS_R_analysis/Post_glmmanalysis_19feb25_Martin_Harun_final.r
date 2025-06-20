library(tidyverse)
library(qqman)
library(fdrtool)
library(dplyr)
#merge the pvalues
wd <- '/home/lilian/Desktop/D/LSTM/Data analysis'

setwd(wd)
# fetch npbo corr results files 
dir.glm.results <- './Glmm_Feb 25/glmm_comb_feb25'
glm.files <- list.files(dir.glm.results)
glm.files <- glm.files [grepl('.rds', glm.files)]
glm.list <- list()

setwd(dir.glm.results)

#merge all ug haplotype npbo files
for (g in glm.files){
  glm.list [[g]] <- readRDS(g)
}


df <- do.call(rbind, glm.list)
df$coef.NA = NULL
df <- na.omit(df)
p <- df$P
snp.id <- df$snp.id

#fdr

fdr <- data.frame(fdrtool(p, statistic = 'pvalue', plot = F))
#df <- data.frame(snp_id,p)
#df.adjusted <- data.frame(snp.id,p=fdr$qval)
df_ad <- data.frame(snp.id,p=fdr$pval)

#significant filter
#sig <-df %>% filter(p<0.0001)
sig <-df_ad %>% filter(p<0.05)
sig <- arrange(sig, p)
#sig2 <- arrange(sig, p)
snpsOfInterest <- sig$snp.id

# MANHATAN
df[c('CHR', 'BP')] <- str_split_fixed(df$snp.id, ':', 2)

df<- df%>% mutate(CHR=fct_recode(CHR, "1"="2RL", "2"="3RL", "3"="X"))

 #checking if chromosome is a factor
#is.factor(df$CHR)
 #converting chromosome to a numeric
df$CHR <- as.numeric(df$CHR)
df$BP <- as.numeric(df$BP)

#df <- df %>% mutate(CHR=fct_recode(CHR, "2RL" = "1", "3RL" = "2", "X" = "3"))
#df_n <- as.data.frame(apply(df[,-1], 2, as.numeric))
## trun column bp and p to numerics if not alread

#df <- df %>%
  #mutate(across(c(p, BP))  # Convert p and BP to numeric

## check if mutate to numeric has been succesful: 
str(df)


#df_n$snps.id <- df$snp.id
#df_n<- df_n[,c(4,1,2,3)]
## omitting nas
df <- na.omit(df)
## visualising in my enviroment 
manhattan(df, snp = 'snp.id', chr = 'CHR', bp = 'BP', p = 'P', 
          chrlabs = c("2RL", "3RL", "X"),
          suggestiveline = F, genomewideline = F, highlight = snpsOfInterest)


## saving svg file
svg('glmm gwas.svg-p.05', height =7, width=9)
manhattan(df, snp = 'snp.id',chr = 'CHR', bp='BP', p='P', chrlabs = c("2RL", "3RL", "X"),
          suggestiveline = F,genomewideline = F, highlight = snpsOfInterest)
dev.off()
