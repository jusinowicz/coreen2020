#=============================================================================
#File: Fit resource utilization functions for Coreen's competition experiment
#=============================================================================
# This is R script to fit resource utilization cuves for Daphnia and 
# Diaphanisoma, where the resource is algae. The fitting is based on an 
# underlying assumption of MacArthur consumer dynamics. 
#=============================================================================
#=============================================================================
#Load libraries
#=============================================================================
library(tidyverse)
library(lubridate)
library(mgcv)
library(gamm4)
library(fields)
source("./functions/info_theory_functions.R")
#=============================================================================
#Load data
#=============================================================================
#Coreen's experiment data

#Preprocessing m1_...
m1=read.csv(file= "tequila_master_data.csv") 
m1$date = ymd(m1$date) #Fix date

#There are certain days with data across individuals from a 
#mesocosm. This is to take all of those measurements and combine
#them into one measurement per day per mesocosm. This will
#use the average: 
cols2group=colnames(m1)[c(2:10,12,15,17:18)] 
cols2summarize = colnames(m1)[c(14,16)]
treats = unique(m1$mesocosm_id)
m1_data_long = m1 %>% 
     group_by_at(cols2group) %>%
     summarise_at(cols2summarize, .funs= mean)
m1_data_long= ungroup(m1_data_long)
m1_data_long$algae_abundance = as.numeric(m1_data_long$algae_abundance)


algae1 = read.csv(file= "algae_export.csv")[,1:3]
#Fix dates in algae:
algae1$date =  parse_date_time(algae1$date, "md")
year(algae1$date) [year(algae1$date) == 0000] = 2019
year(algae1$date) [year(algae1$date) == 2020 & 
    month(algae1$date) > 9  & month(algae1$date) <=12  ] = 2019

algae1$date = ymd (algae1$date)

#Join algae to m1_...
m1_data_long= m1_data_long %>% 
  left_join( algae1) 

colnames(m1_data_long)[colnames(m1_data_long) == "zooplankton_abundance"] = "N"

#=============================================================================
#Fit a GAMM to the data to predict algal abundance on the basis of density per-species, 
#as a funciton of temperature, with mesocosm_id as a Random Effect
#=============================================================================
#Note that the data are fir with the log link function (log transformed) because 
#values cannot be <0. 
alg_gam = gam( algae_abundance ~ s(N)+s(temperature,k=5)+s(species, bs="re")+s(mesocosm_id,bs="re"),
  family=Gamma(link='log'), data=m1_data_long )

#Use the fitted GAMM to interpolate NAs for algal abundance in the data set
alg_newdata = subset(m1_data_long, is.na(algae_abundance) )
alg_replace = as.vector(predict.gam(alg_gam, newdata=alg_newdata,type= "response") )
alg_newdata$algae_abundance = alg_replace
m1_data_long$algae_abundance[is.na(m1_data_long$algae_abundance)] = alg_newdata$algae_abundance

m1_data_long %>% 
ggplot(aes(y = algae_abundance, x = N, color = species, group = interaction(species,replicate_number)))+
  geom_point()
ggsave("./algae_projection1.pdf", width = 8, height = 10)

#=============================================================================
#Plot the data
#=============================================================================
m1_data_long %>% 
  ungroup() %>% 
  filter(!is.na(N)) %>% 
  ggplot(aes(x = day_n, y = N, color = species, group = interaction(species,replicate_number)))+
  geom_line()+
  facet_grid(temperature~invade_monoculture)+
  scale_y_log10()+
  theme_bw()+
  scale_color_manual(values = c("dodgerblue1", "red1"), name = NA)+
  xlab("day")+
  ylab("population size")+
  theme(strip.background = element_rect(colour=NA, fill=NA))

#=============================================================================
#Collect basic data on number of treatments, species, etc.
#=============================================================================
nspp=2
temps = (unique(m1_data_long$temperature))
ntreatments =  length(unique(m1_data_long$temperature))
mesos = unique(m1_data_long$mesocosm_id)
nmesos =  length(mesos)
rspecies = unique(m1_data_long$species)
invader = unique(m1_data_long$invade_monoculture)
invasions_per = length(treats)/4 #Treatment entries correponding to each spp 
inv_day = 28 #The first day of attempted invasion
no_reps = 18 #The number of replicated mesocosms total per resident/invader

#=============================================================================
#Fit a model at each treatment level (temperature)
#=============================================================================
#Use both the monoculture data and the data for the resident pre-invasion
mesos1 = subset(m1_data_long,species == "daphnia" & invade_monoculture == "monoculture" )
mesos2 = subset(m1_data_long,species == "daphnia" & invade_monoculture == "dia invade" )
mesos2 = mesos2[ mesos2$day_n<=inv_day, ]
mesos_daph = rbind(mesos1,mesos2)

mesos1 = subset(m1_data_long,species == "diaphanosoma" & invade_monoculture == "monoculture" )
mesos2 = subset(m1_data_long,species == "diaphanosoma" & invade_monoculture == "daph invade" )
mesos2 = mesos2[ mesos2$day_n<=inv_day, ]
mesos_dia = rbind(mesos1,mesos2)

#Fit the resource consumption model as a function of temperature:  



#Fit the resource consumption model per temperature. 
cl_daph = vector("list", 6)
cl_dia = vector("list", 6)

for(t in 1:6) { 

  #Get the data for each temperature. 
  daph_tmp = subset(mesos_daph, temperature == temps[t])
  dia_tmp = subset(mesos_dia, temperature == temps[t])

  #Add the rate of population growth, Ndiff
  daph_tmp = daph_tmp %>%    
            arrange(replicate_number) %>%
            mutate(Ndiff = lead(N-lag(N),)/lead(day_n-lag(day_n),) ) #"lead" lines up the result 
  daph_tmp$Ndiff[is.infinite(daph_tmp$Ndiff)] = NA

  dia_tmp = dia_tmp %>%    
            arrange(replicate_number) %>%
            mutate(Ndiff = lead(N-lag(N),)/lead(day_n-lag(day_n),) ) #"lead" lines up the result 
  dia_tmp$Ndiff[is.infinite(dia_tmp$Ndiff)] = NA

  #The basic MacArthur model is a linear consumption rate so just fit with a GLMM
  # cl_daph[[t]] = gam( Ndiff ~ +s(temperature,k=5)+s(mesocosm_id,bs="re"),family=Gamma(link='log'), data=daph_tmp)
  # cl_dia[[t]] = gam( Ndiff ~ +s(temperature,k=5)+s(mesocosm_id,bs="re"),family=Gamma(link='log'), data=dia_tmp)

  cl_daph[[t]] = gam( Ndiff ~ algae_abundance +s(mesocosm_id,bs="re"), data=daph_tmp)
  cl_dia[[t]] = gam( Ndiff ~ algae_abundance +s(mesocosm_id,bs="re"), data=dia_tmp)


}