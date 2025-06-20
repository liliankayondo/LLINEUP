setwd("/home/lilian/Desktop/D/LSTM/Data analysis")
library(readxl)
data_raw <- read_excel("./filtered_data_usedforanaysis_2.11.24.xlsx")


#Libraries
library(tidyverse)
library(jpeg)
library(ggplot2)
library (dplyr)
library(glmmTMB)


# select a few columns to use in the data

data <- data_raw %>% dplyr::select(RND, EA, HH, LLIN.actual, MolSpecies)

# filtering out funestus 
data2 <- data %>% dplyr::filter(MolSpecies== "FUN")

# filterng na's in the LLIN.actual column
data_PBO <- data %>% dplyr::filter(!is.na(LLIN.actual))

# grouping my coluns of intrest to count fun densties
data_PBO_summary <- data_PBO %>%
  group_by(EA, HH, LLIN.actual, RND) %>%
  summarise(MolSpecies_count = n())

# changing rnd column to factor
data_PBO_summary$RND <- factor(data_PBO_summary$RND, levels = c("Baseline", "6 months", "12 months", "18 months", "25 months"))
ggplot(data_PBO_summary, aes(y = MolSpecies_count, 
                             x = factor(RND), 
                             fill = factor(LLIN.actual))) +
  geom_boxplot() +
  theme_bw() +
  scale_y_continuous(name = "Mean number of An. funestus per Household per cluster") +
  scale_x_discrete(name = "") +
  ggtitle("An. funestus densities across the trial period") +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),   
    panel.grid.minor = element_blank()      
  )



## adding plot function to see data distribution
ggplot(data_PBO_summary, aes(x = MolSpecies_count)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of An. funestus Counts", 
       x = "Mosquito Count per Household", 
       y = "Frequency")




## log transforming the data before doing normalisation test again 



## checking if my mean and variance differ to dtermine which model to use
mean_val <- mean(data_PBO_summary$MolSpecies_count)
var_val <- var(data_PBO_summary$MolSpecies_count)
mean_val; var_val

## using negative binomial since my data is over dispersed, ie mean is much less than variance
library(MASS)
#nb_model <- glm.nb(MolSpecies_count ~ RND + LLIN.actual, data = data_PBO_summary)
#summary(nb_model)


### trying gwas model


nb_model <- glmmTMB(MolSpecies_count ~ RND+LLIN.actual+(1|HH), 
                 data = data_PBO_summary, family= nbinom2)
summary(nb_model)

     ## check if model fits
library(performance)
check_overdispersion(nb_model)

