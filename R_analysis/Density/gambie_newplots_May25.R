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


# Fix the factor levels
ggplot(data_PBO_summary3, aes(y = MolSpecies_count, 
                             x = factor(RND), 
                             fill = factor(LLIN.actual))) +  # remove colour from here
  geom_boxplot(outlier.shape = 16, outlier.size = 2, whisker.width = 0,
               aes(colour = factor(LLIN.actual)),  # assign colour here for plot only
               alpha = 0.8) +
  
  scale_y_continuous(name = "Mean number of An. gambiae per Household per cluster") +
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
  
  ggtitle("An. gambiae densities across the trial period") +
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank() 
  )
