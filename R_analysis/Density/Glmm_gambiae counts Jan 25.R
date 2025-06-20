setwd("/home/lilian/Desktop/D/LSTM/Data analysis")
library(readxl)
data_raw <- read_excel("./PBO_Ento_ALL_molecular_data_V4_10-03-21_merged-GRAPHS.xlsx")


#Libraries
library(tidyverse)
library(jpeg)
library(ggplot2)
library (dplyr)


# select a few columns to use in the data

data <- data_raw %>% dplyr::select(RND, EA, HH, LLIN.actual, MolSpecies)

# filtering out funestus 
data2 <- data %>% dplyr::filter(MolSpecies %in% c("S", "M")) %>%
  dplyr::filter(!is.na(LLIN.actual))

# filterng na's in the LLIN.actual column
data_PBO <- data2 %>% dplyr::filter(!is.na(LLIN.actual))

# grouping my coluns of intrest to creat a column with count  densties
data_PBO_summary <- data_PBO %>%
  group_by(EA, HH, LLIN.actual, RND, MolSpecies) %>%
  summarise(MolSpecies_count = n(), .groups = "drop")

# Replace numeric values in RND with descriptive labels- code from ai
data_PBO_summary2 <- data_PBO_summary %>%
  mutate(RND = case_when(
    RND == 1 ~ "Baseline",
    RND == 2 ~ "6 months",
    RND == 3 ~ "12 months",
    RND == 4 ~ "18 months",
    RND == 5 ~ "25 months",
  ))

### Doing the same as above for llins.actual
data_PBO_summary3 <- data_PBO_summary2 %>%
  mutate(LLIN.actual = case_when(
    LLIN.actual == 1 ~ "Non-PBO",
    LLIN.actual == 2 ~ "PBO"
  ))
## leaving them into an oder as surveys
data_PBO_summary3$RND <- factor(data_PBO_summary3$RND, 
                                levels = c("Baseline", "6 months", "12 months", "18 months", "25 months"))


## adding plot function to see data distribution
ggplot(data_PBO_summary3, aes(x = MolSpecies_count)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of An. gambiae Counts", 
       x = "Mosquito Count per Household", 
       y = "Frequency")




## log transforming the data before doing normalisation test again 

data_PBO_summary3$log_MolSpecies_count <- log(data_PBO_summary3$MolSpecies_count + 1)
ggplot(data_PBO_summary3, aes(x = log_MolSpecies_count)) +
  geom_histogram(binwidth = 0.2, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Log-Transformed Distribution of An. gambiae Counts", 
       x = "Log(Mosquito Count + 1)", 
       y = "Frequency")



## checking if my mean and variance differ to dtermine which model to use
mean_val <- mean(data_PBO_summary3$MolSpecies_count)
var_val <- var(data_PBO_summary3$MolSpecies_count)
mean_val; var_val

## using negative bionmial since my data is over dispersed, ie mean is much less than variance
library(MASS)
nb_model <-glmmTMB(MolSpecies_count ~ RND+LLIN.actual+(1|HH), 
        data = data_PBO_summary3,family= nbinom2)
summary(nb_model)









