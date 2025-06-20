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

data_PBO_summary$RND <- factor(data_PBO_summary$RND, 
                               levels = c("Baseline", "6 months", "12 months", "18 months", "25 months"))

# Fix the factor levels
# Fix the factor levels
ggplot(data_PBO_summary, aes(y = MolSpecies_count, 
                              x = factor(RND), 
                              fill = factor(LLIN.actual))) +  # remove colour from here
  geom_boxplot(outlier.shape = 16, outlier.size = 2, whisker.width = 0,
               aes(colour = factor(LLIN.actual)),  # assign colour here for plot only
               alpha = 0.8) +
  
  scale_y_continuous(name = "Mean number of An. funestus per Household per cluster") +
  scale_x_discrete(name = "") +
  coord_cartesian(ylim = c(0, 75)) +
  
  scale_fill_manual(
    values = c("Non-PBO" = "#c28374", "PBO" = "#52607b"),
    name = "LLIN Type"
  ) +
  scale_colour_manual(
    values = c("Non-PBO" = "#c28374", "PBO" = "#52607b"),
    guide = "none"  # hide colour legend
  ) +
  
  ggtitle("An. funestus densities across the trial period") +
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank() 
  )


