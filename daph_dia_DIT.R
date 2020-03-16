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
library(fields)
source("./functions/info_theory_functions.R")
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

#=============================================================================
#Collect basic data on number of treatments, species, etc.
#=============================================================================
ntreatments =  length(unique(m1_data_long$temperature))
treats = unique(m1_data_long$mesocosm.ID)
invader = unique(m1_data_long$species)
invasions_per = length(treats)/4 #Treatment entries correponding to each spp 
inv_day = 28 #The first day of attempted invasion
no_reps = 18 #The number of replicated mesocosms total per resident/invader

#=============================================================================
# Variables for data collection
#=============================================================================
out1 = list(matrix(0,ntreatments*2,1))

#=============================================================================
# First, loop over both diaphanosoma and daphnia when diaphanasoma is the 
# invader
#=============================================================================
#The resident and invader: 
i = 1
  #Loop over treatments 
  u_invader = invader[-i]
  u_res = invader[i]

  for (n in 1: (ntreatments)){
    #Since the treatments are regular, set the start/end position of
    #the subset with:  
    pos1 = (n-1)*2 + n
    pos2 = pos1+2
    u_treats = treats[pos1:pos2]
    index1 = n+(i-1)*ntreatments

    #=============================================================================
    #Make new resident and invader data sets.
    #=============================================================================
    #Make a new resident data set to fit the growth rate function with nlme/nls 
    #The new data set includes a column for n(t+1)/n(t) and a columnf for the time 
    #interval across subsequent measurements.  
    res_tmp = subset(m1_data_long, mesocosm.ID == u_treats &  species == u_res   )
        res_tmp = na.exclude(res_tmp) #NLME won't work with NAs 
        #Arrange the data by replicate number, then add a new column for the delta N
        res_tmp = res_tmp %>% 
            arrange(replicate.number) 
    res_tmp = res_tmp %>% mutate(N_res = N)

    #Make a new invader data set to fit the growth rate function with nlme/nls 
    #The new data set includes a column for n(t+1)/n(t) and a columnf for the time 
    #interval across subsequent measurements.  
    inv_tmp = subset(m1_data_long, mesocosm.ID == u_treats &  species == u_invader   )
        inv_tmp = na.exclude(inv_tmp) #NLME won't work with NAs 
        #inv_tmp$N[is.na(inv_tmp$N)] = 0 #Replace NAs with 0  
        #Arrange the data by replicate number, then add a new column for the delta N
        inv_tmp = inv_tmp %>% 
            arrange(replicate.number) 
    inv_tmp = inv_tmp %>% mutate(N_inv = N)
    #This line adds the resident densities on the matching days from the matching 
    #replicates to the data table inv_tmp
    inv_tmp = inv_tmp %>% left_join( select( res_tmp, replicate.number,day_n,N_res),
     by= (c( "day_n" = "day_n", "replicate.number"="replicate.number" )) )
   
    #Create a dummy data set in the few cases that experiment was essentially 
    #unsuccesful: 
    if(dim(inv_tmp)[1] <=0 ){
      res_tmp = res_tmp %>% left_join( select( inv_tmp, replicate.number,day_n,N_inv),
       by= (c( "day_n" = "day_n", "replicate.number"="replicate.number" )) )
          res_tmp$N_inv[is.na(res_tmp$N_inv)] = 0 #Replace NAs with 0  
      mydata_res = subset(res_tmp, day_n >= inv_day  )
      mydata = cbind(mydata_res$N_inv,mydata_res$N)
      out1[n] = list(as.matrix(na.exclude(mydata)))

    } else {
      mydata_inv = subset(inv_tmp, day_n >= inv_day  )
      mydata = cbind(mydata_inv$N,mydata_inv$N_res)
      out1[n] = list(as.matrix(na.exclude(mydata)))
    }
    
}


#=============================================================================
  # Information processing networks
  #=============================================================================
  ## This code takes the population time-series counts output by the ODEs and 
  ## calculates Excess Entropy, Active Information Storage, and Transfer Entropy.
  ## Each quantity is calculated at both the average and local level.  
  #=============================================================================
  # This function gives:
  # EE_mean   Average mutual information per species
  # AI_mean   Average active information per species
  # TE_mean   Average transfer entropy per species
  # 
  # EE_local    Local mutual information per species
  # AI_local    Local active information per species
  # TE_local    Local transfer entropy per species
  #=============================================================================
  nt1 = 1
  nt2 = tl
  f1 = 100 #scaling factor
  di_web[w] = list(get_info_dynamics(pop_ts = floor(f1*out1[[w]]$out[nt1:tl,2:(nspp+1)]), 
    k=k,with_blocks=FALSE))

  ## This code takes the population time-series counts output by the ODEs and 
  ## calculates the average Transfer Entropy from each species to every other 
  ## species. The goal is to get an overview of the major information pathways 
  ## in the web.   
  #=============================================================================
  # This function gives:
  # te_web    Average transfer entropy per species as a pairwise matrix
  #=============================================================================
  te_web[w] = list( get_te_web( pop_ts = floor(f1*out1[[w]]$out[nt1:tl,2:(nspp+1)]), 
    k=k) )