#=============================================================================
# R code to create to explore the Dynamic Information Theoretic properties of 
# a two-species consumer-resource system. 
# Data are taken from an experiment where a Daphnia and Diaphanosoma species were 
# grown in mesocosms at temperatures ranging from 20 - 32C and fed algae daily. 
#
# The theoretical models include different versions of the MacArthur RC model 
# where parameters are fit from the experiments. 
#
# The DIT metrics are chosen to quantify the overall complexity of the competitive
# dynamics across the different temperature treatments. 
#
#=============================================================================
#Load libraries
#=============================================================================
library(deSolve)
library(tidyverse)
library(lubridate)
library(mgcv)
library(fields)
library(zoo)
library(plyr)
library(gridExtra)
#Custom functions:
source("./functions/food_web_functions.R")
source("./functions/info_theory_functions.R")
source("./functions/resource_fitting.R")
#=============================================================================
#STAGE 1: Data loading and preprocessing: 
#=============================================================================
#=============================================================================
#Load data
#=============================================================================
#Coreen's experiment data
#Preprocessing m1_...
m1=read.csv(file= "tequila_master_data.csv") 
algae1 = read.csv(file= "algae_export.csv")[,1:3]

#Fix dates in both data sets:
m1$date = ymd(m1$date)
algae1$date =  parse_date_time(algae1$date, "md")
year(algae1$date) [year(algae1$date) == 0000] = 2019
year(algae1$date) [year(algae1$date) == 2020 & 
    month(algae1$date) > 9  & month(algae1$date) <=12  ] = 2019
algae1$date = ymd (algae1$date)

#There are dates where multiple individuals were measured. Replace individual
#measurements with a group average so that there is only one sample per day:

cols2group=colnames(m1)[c(2:10,12,15,17:18)] #Which columns to group
cols2summarize = colnames(m1)[c(14,16)] #Which columns to summarize 
m1_data_long = m1 %>% 
       group_by_at(cols2group) %>%
       summarise_at(cols2summarize, .funs= mean)
m1_data_long= ungroup(m1_data_long)
m1_data_long$algae_abundance = as.numeric(m1_data_long$algae_abundance)
colnames(m1_data_long)[colnames(m1_data_long) == "zooplankton_abundance"] = "N" 

#Join algae to m1_...
m1_data_long= m1_data_long %>% 
  left_join( algae1) 

#Fix NAs in algal abundance
m1_data_long$algae_cells_mL[is.na(m1_data_long$algae_cells_mL)] = mean(algae1$algae_cells_mL,na.rm=T)
m1_data_long$algae_mL_media[is.na(m1_data_long$algae_mL_media)] = mean(algae1$algae_mL_media,na.rm=T)

#Fix NA dates. These come up on the half days, so can be replaced with lag of dates
m1_data_long$N[ is.na(m1_data_long$date)] = lag(m1_data_long$N)[ is.na(m1_data_long$date)]
m1_data_long$algae_cells_mL[ is.na(m1_data_long$date)] = lag(m1_data_long$algae_cells_mL)[ is.na(m1_data_long$date)]
m1_data_long$date[ is.na(m1_data_long$date)] = lag(m1_data_long$date)[ is.na(m1_data_long$date)]
#=============================================================================
#Plot the data
#=============================================================================
m1_data_long %>% 
  ungroup() %>% 
  filter(!is.na(N)) %>% 
  ggplot(aes(x = day_n, y = N, color = species, group = interaction(species,replicate_number)))+
  geom_line()+
  facet_grid(temperature~invade_monoculture)+
  #scale_y_log10()+
  theme_bw()+
  scale_color_manual(values = c("dodgerblue1", "red1"), name = NA)+
  xlab("day")+
  ylab("population size")+
  theme(strip.background = element_rect(colour=NA, fill=NA))



#=============================================================================
#Replace missing measurements of algal abundance. This uses s GAMM to 
#predict on the basis of density per-species, as a funciton of temperature, 
#with mesocosm_id as a Random Effect.
#Note that the data are fit with the log link function (log transformed) because 
#values cannot be <0. 
#=============================================================================
m1_data_long$algae_abundance[is.na(m1_data_long$algae_abundance)] = 
	alg_replace(m1_data_long)$algae_abundance 
#=============================================================================
#Plot the data
#=============================================================================
m1_data_long %>% 
ggplot(aes(y = algae_abundance, x = N, color = species, group = interaction(species,replicate_number)))+
  geom_point()+
  facet_grid(temperature~species)

#=============================================================================
#Collect/set basic data on number of treatments, species, etc.
#=============================================================================
#These are set manually: 
inv_day = 28 #The first day of attempted invasion
inv_end = 50#The last day of invasion conditions
no_reps = 18 #The number of replicated mesocosms total per resident/invader
nspp=2

#These are from data: 
temps = (unique(m1_data_long$temperature))
ntemps =  length(unique(m1_data_long$temperature))
mesos = unique(m1_data_long$mesocosm_id)
nmesos =  length(mesos)
rspecies = unique(m1_data_long$species)
invader = unique(m1_data_long$invade_monoculture)
invasions_per = length(mesos)/4 #Treatment entries correponding to each spp 

#=============================================================================
#Make the single species data set by isolating mesocosms and portions of time
#series where only a single species is present. Use both the monoculture data 
#and the data for the resident pre-invasion for  resource consumption fits, 
#intraspecific competition, intrinsic growth
#=============================================================================
mesos1 = subset(m1_data_long,species == "daphnia" & invade_monoculture == "monoculture" )
mesos2 = subset(m1_data_long,species == "daphnia" & invade_monoculture == "dia invade" )
mesos2 = mesos2[ mesos2$day_n<=inv_day, ]
mesos_daph = rbind(mesos1,mesos2)

mesos1 = subset(m1_data_long,species == "diaphanosoma" & invade_monoculture == "monoculture" )
mesos2 = subset(m1_data_long,species == "diaphanosoma" & invade_monoculture == "daph invade" )
mesos2 = mesos2[ mesos2$day_n<=inv_day, ]
mesos_dia = rbind(mesos1,mesos2)
rm(mesos1,mesos2)

#Make the invasion data sets by isolating portions of time series where a species
#invades. 
mesos_inv2 = subset(m1_data_long, invade_monoculture == "daph invade" )
mesos_inv2 = mesos_inv2[ mesos_inv2$day_n>=inv_day & mesos_inv2$day_n <= inv_end , ]
mesos_inv_daph = mesos_inv2 

mesos_inv2 = subset(m1_data_long,invade_monoculture == "dia invade" )
mesos_inv2 = mesos_inv2[ mesos_inv2$day_n>=inv_day & mesos_inv2$day_n <= inv_end , ]
mesos_inv_dia= mesos_inv2 
rm(mesos_inv2)

#=============================================================================
#STAGE 2: Fit resource consumption rates, intrinsic growth rates, and 
#phenomenological competition coefficients from time series of the data.
#=============================================================================

#Fit the resource consumption model per temperature. 
cl_daph = vector("list", 6); cl_daph_plot = NULL; cl_daph_pred = NULL
cl_dia = vector("list", 6); cl_dia_plot = NULL; cl_dia_pred = NULL

#Fit the consumer growth rate per temperature 
cR_daph = vector("list", 6); cR_daph_pred = NULL;
cR_dia = vector("list", 6);cR_dia_pred = NULL;

#Fit competition models per temperature 
#Intraspecific
igr_daph = vector("list", 6); lvii_daph = vector("list", 6);igr_daph_pred = NULL
lvii_daph_pred = NULL
igr_dia = vector("list", 6);lvii_dia = vector("list", 6);igr_dia_pred = NULL;
lvii_dia_pred = NULL;
#Interspecific
igrij_daph = vector("list", 6);lvij_daph = vector("list", 6);igrij_daph_pred = NULL;
lvij_daph_pred = NULL;
igrij_dia = vector("list", 6); lvij_dia = vector("list", 6); igrij_dia_pred = NULL;
lvij_dia_pred = NULL;



#par(mfrow=c(6,1),oma = c(5,4,0,0) + 0.1,mar = c(0,0,1,1) + 0.1)

algae_start = 1E7 #This was the target value
#For comparison with target value: on average, it's pretty damn close (9.8E6 vs 10E6)
#algae_start = mean(m1_data_long$algae_cells_mL*m1_data_long$algae_mL_media,na.rm=T)

maxN = 500 #max(m1_data_long$N,na.rm=T)
maxT = max(m1_data_long$day_n,na.rm=T)

for(t in 1:ntemps) { 

}
