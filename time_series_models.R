
#=============================================================================
#Load libraries
#=============================================================================
library(tidyverse)
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

#=============================================================================
# Use nlme to fit an AR model to each treatment time series. 
# Temporal autocorrelation is treated using corARMA correlation structure. 
# Mesocosm is a random effect. 
#=============================================================================
ntreatments =  length(unique(m1_data_long$temperature))
treats = unique(m1_data_long$mesocosm.ID)
invader = unique(m1_data_long$invade.mono)

all_models = list(matrix(0,ntreatments*2,1)) #To store all of the fitted models

#Loop over invader
for i(in 1:2) {  
	#Loop over treatments 
	for (n in 1: (ntreatments)){
		#Since the treatments are regular, set the start/end position of
		#the subset qith:  
		pos1 = (n-1)*2 + n
		pos2 = pos1+2
		u_treats = treats[pos1:pos2]
		u_invader = invader[i]
		mydata = subset(m1_data_long, mesocosm.ID == u_treats &  invade.mono == u_invader  )
		all_models[n] = list (lme(height ~ type * time, random= ~ 1|box/plant, 
			correlation=corARMA(0.2, form=~time|box/plant, p=1, q=0), data=mydata))
	}
}






