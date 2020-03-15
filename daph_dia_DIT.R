#=============================================================================
#Information theory and population dynamics 
#=============================================================================
# This is R code to take population time series of two competing populations
# from a controlled (lab) experiment. The data are time series of Daphnia and
# Diaphanosoma species at 5 different Temperature treatment levels. The 
# different mesocosms also correspond with either an invasion or single species
# resident treatment. 
#
# This analysis measures several specific dynamic information theoretic 
# features of the time series: 
#   1. Excess entropy -- interpreted as the memory required/stored to produce
#      a particular output, and also equated to the complexity of the process.
#   2. Transient information -- information that is transferred between 
#      different sub-processes (i.e. population time series) of the system. 
#=============================================================================

#=============================================================================
#Load libraries
#=============================================================================
library(tidyverse)
source( )
#=============================================================================
#Load data
#=============================================================================
#Coreen's experiment data
m1=read.csv(file= "mesocosm_experiment_coreen_forbes.csv") 
#m1_inv = m1[1:72,] #Just take the invasion experiments

m1_data_long = m1 %>% 
  gather(key = day, value = N, day1:day83) %>% 
  dplyr::select(-X) %>% 
  mutate(day_n = as.numeric(gsub(x = day, pattern = "[^0-9.-]", replacement = "")))

m1_data_long$replicate.number= factor(m1_data_long$replicate.number)
#=============================================================================
#Plot the data
#=============================================================================
m1_data_long %>% 
  ungroup() %>% 
  filter(!is.na(N)) %>% 
  ggplot(aes(x = day_n, y = N, color = species, group = interaction(species,replicate.number)))+
  geom_line()+
  facet_grid(temperature~invade.mono)+
  scale_y_log10()+
  theme_bw()+
  scale_color_manual(values = c("dodgerblue1", "red1"), name = NA)+
  xlab("day")+
  ylab("population size")+
  theme(strip.background = element_rect(colour=NA, fill=NA))
ggsave("./Downloads/Coreen_plot.pdf", width = 8, height = 10)
Collapse
