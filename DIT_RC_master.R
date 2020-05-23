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
#STAGE 2: Fit resource consumption rates, intrinsic growth rates, and 
#phenomenological competition coefficients from time series of the data.
# Choose which model to fit for each of the different fits by uncommenting from
# the list below.
# Current favorite is starred. 
#=============================================================================
# Resource consumption models
#=============================================================================
algae_start = 1E7 #This was the target value of algal cells
#1. Type 2 (saturating) response
#alg2 = formula (Adiff ~  a1/(1+b1*exp(c1*N) ) ); a1 = algae_start

#2. Decaying logistic/exponential response
#alg2 = formula (algae_abundance/algae_start ~  1/(1+b1*exp(c1*N) ) )

#3. *** Decaying exponential response ***
alg2 = formula (algae_abundance/algae_start~  exp(c1*N+b1) )

#3b. Alternative decaying exponential 
#alg2 = formula (Adiff ~  exp(c1*N+b1) )

#=============================================================================
# Intrinsic growth rate rate models -- as a function of algal consumption 
#=============================================================================
#1. Exponential as a function of Adiff
#cR = formula (N~exp(c1*Adiff+b1) )

#2. Exponential, after accountinf for species-specific consumption rate:
#	For this one, p1 is a constant that must be set as 
#	 p1 = abs(coef(cl_daph[[w]])[1]) or  p1 = abs(coef(cl_dia[[w]])[1])
#cR = formula (N~exp(c1*p1*Adiff+b1) )

#3. *** Exponential as a function of algal abundance. ***
cR = formula ( N ~  exp(c1*algae_abundance+b1)  )

#=============================================================================
# Intrinsic growth rate rate models -- as a function of time
# These are used for a cleaner fit of competition coefficients. 
#=============================================================================
#1. Exponential over initial growth phase. This requires a short time window, 
#   i.e. use inv_end to remove equilibrium densities.
# igr = formula (N ~  a1*exp(c1*day_n ) )

#2. *** Saturating logistic function. ***    
#   Note: a1 (the carrying capacity) could be set as a constant instead of 
#   being fit if there are convergence issues. 
igr = formula (N ~  a1/(1+b1*exp(c1*day_n ) ) )

#3. and 4., The interspecific invasion complements to 1. and 2.: 
#3. 
#igrij = formula (N_inv ~  (b1*exp(c1*day_n ) ) )

#4. ***    ***
igrij = formula (N_inv ~  a1/(1+b1*exp(c1*day_n ) ) )

#=============================================================================
# Competition models
# These fall into two approaches. The first is used in combo with the fitted
# IGR models above, in which case they are just linear models. 
# The second ignores the IGR models and fits directly from the data. 
# While the second is more conceptually ideal, it seems to fail mostly with 
# these data. 
#=============================================================================
#=============================================================================
# First approach, in combo with IGR models
#=============================================================================
# Intraspecific competition models
f_lvii = Ndiff ~ N
# Interspecific competition models
f_lvij = Ndiff ~ N

#=============================================================================
#Second approach
#With all of these models, ri may be considered as fixed or as a free
#parameter. 
#For interspecific, the intraspecific competition coefficient can be included
#in the model. In this case, it works best to treat it as a fixed parameter
#and set it equal to the fit obtained from the intraspecific step. 
#=============================================================================
####################################
#Intraspecific competition models

#1. Log transformed, Lotka Volterra 
#f_lvii = formula ( I(log(Ndiff)) ~  ri*(ki-aii*N) )

#2. Untransformed Lotka Volterra 
#f_lvii = formula ( Ndiff ~  ri*(1-aii*N) )

#3. *** Leslie-Gower ***
#f_lvii = formula ( Ndiff ~  ri/(1+aii*N) )

####################################
#Interspecific competition models

#1. Lotka Volterra, aij only 
#f_lvij = formula ( Ndiff ~  ri*(1-aij*N_res) )

#2. Leslie Gower, aij only. 
#f_lvij = formula ( Ndiff ~  ri/(1+aij*N_res) )

#3. LV, aii and aij
#f_lvij = formula ( Ndiff ~  ri/(1+aii*N_inv+aij*N_res) )

#4. LG, aii and aij
#f_lvij = formula ( Ndiff ~  ri*(1-aii*N_inv-aij*N_res) )

#=============================================================================
#Variables:

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

maxN = 500 #max(m1_data_long$N,na.rm=T)
maxT = max(m1_data_long$day_n,na.rm=T)

for(t in 1:ntemps) { 
  #=============================================================================
  #Resource consumption and intraspecific competition
  #=============================================================================
  #Make the single species data set by isolating mesocosms and portions of time
  #series where only a single species is present. Use both the monoculture data 
  #and the data for the resident pre-invasion for resource consumption fits, 
  #intraspecific competition, intrinsic growth
  #
  #Get the data for each temperature. 
  #=============================================================================
  daph_tmp = get_new_mono(m1_data_long, rspecies[1], temps[t], inv_day )
  dia_tmp = get_new_mono(m1_data_long, rspecies[2], temps[t], inv_day )

  



  #Make the invasion data sets by isolating portions of time series where a species
  #invades. 
  daph_inv_tmp = get_new_inv(m1_data_long, rspecies[1], temps[t], inv_day, inv_end )
  dia_inv_tmp = get_new_inv(m1_data_long, rspecies[2], temps[t], inv_day, inv_end )

  #Output
  print(mean(as.data.frame(subset(daph_inv_tmp,day_n <=34))$Ndiff,na.rm=T) )
  print(mean(as.data.frame(subset(dia_inv_tmp,day_n <=34))$Ndiff,na.rm=T) )
  
  #Keep building these for plotting: 
  cl_daph_plot = rbind( cl_daph_plot, daph_tmp )
  cl_dia_plot = rbind( cl_dia_plot, dia_tmp )



}
