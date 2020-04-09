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
algae1$algae_cells_mL[is.na(algae1$algae_cells_mL)] = mean(algae1$algae_cells_mL,na.rm=T)
algae1$algae_mL_media[is.na(algae1$algae_mL_media)] = mean(algae1$algae_mL_media,na.rm=T)

#Join algae to m1_...
m1_data_long= m1_data_long %>% 
  left_join( algae1) 

colnames(m1_data_long)[colnames(m1_data_long) == "zooplankton_abundance"] = "N"

#Fix NAs in algal abundance
m1_data_long$algae_cells_mL[is.na(m1_data_long$algae_cells_mL)] = mean(algae1$algae_cells_mL,na.rm=T)
m1_data_long$algae_mL_media[is.na(m1_data_long$algae_mL_media)] = mean(algae1$algae_mL_media,na.rm=T)

#Fix NA dates. These come up on the half days, so can be replaced with lag of dates
m1_data_long$N[ is.na(m1_data_long$date)] = lag(m1_data_long$N)[ is.na(m1_data_long$date)]
m1_data_long$algae_cells_mL[ is.na(m1_data_long$date)] = lag(m1_data_long$algae_cells_mL)[ is.na(m1_data_long$date)]
#m1_data_long$algae_abundance[ is.na(m1_data_long$date)] = lag(m1_data_long$algae_abundance)[ is.na(m1_data_long$date)]
m1_data_long$date[ is.na(m1_data_long$date)] = lag(m1_data_long$date)[ is.na(m1_data_long$date)]


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
  geom_point()+
  facet_grid(temperature~species)

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
#Are these effectively consumption rates? 
m1_data_long %>% 
ggplot(aes(y =N, x =mean(algae_cells_mL,na.rm=T)- algae_abundance, color = species, group = interaction(species,replicate_number)))+
  geom_point()+
  facet_grid(temperature~species)

#Use both the monoculture data and the data for the resident pre-invasion
# inv_day = 1000
mesos1 = subset(m1_data_long,species == "daphnia" & invade_monoculture == "monoculture" )
mesos2 = subset(m1_data_long,species == "daphnia" & invade_monoculture == "dia invade" )
mesos2 = mesos2[ mesos2$day_n<=inv_day, ]
mesos_daph = mesos1 #rbind(mesos1,mesos2)

mesos1 = subset(m1_data_long,species == "diaphanosoma" & invade_monoculture == "monoculture" )
mesos2 = subset(m1_data_long,species == "diaphanosoma" & invade_monoculture == "daph invade" )
mesos2 = mesos2[ mesos2$day_n<=inv_day, ]
mesos_dia = mesos1 #rbind(mesos1,mesos2)

#Fit the resource consumption model as a function of temperature:  




#Fit the resource consumption model per temperature. 
cl_daph = vector("list", 6)
cl_dia = vector("list", 6)

cl_daph_plot = NULL
cl_dia_plot = NULL

cl_daph_pred = NULL
cl_dia_pred = NULL

par(mfrow=c(6,1),oma = c(5,4,0,0) + 0.1,mar = c(0,0,1,1) + 0.1)

algae_start = mean(m1_data_long$algae_cells_mL,na.rm=T)

for(t in 1:6) { 

  #Get the data for each temperature. 
  daph_tmp = subset(mesos_daph, temperature == temps[t])
  dia_tmp = subset(mesos_dia, temperature == temps[t])

  #Add the rate of population growth, Ndiff
  daph_tmp = daph_tmp %>%    
            arrange(replicate_number) %>%
            mutate(Ndiff = (N-lag(N))/(day_n-lag(day_n)) ) #"lead" lines up the result 
  daph_tmp$Ndiff[is.infinite(daph_tmp$Ndiff)] = NA

  dia_tmp = dia_tmp %>%    
            arrange(replicate_number) %>%
            mutate(Ndiff = (N-lag(N))/(day_n-lag(day_n)) ) #"lead" lines up the result 
  dia_tmp$Ndiff[is.infinite(dia_tmp$Ndiff)] = NA


  #Add the rate of algal consumption, Adiff
  daph_tmp = daph_tmp %>%    
            arrange(replicate_number) %>%
            mutate(Adiff = algae_start-algae_abundance)  #"lead" lines up the result 
            #mutate(Adiff = exp(N*0.01)+algae_abundance )  #"lead" lines up the result 
  daph_tmp$Adiff[is.infinite(daph_tmp$Adiff)] = NA

  dia_tmp = dia_tmp %>%    
            arrange(replicate_number) %>%
            mutate(Adiff = algae_start-algae_abundance)#"lead" lines up the result 
  dia_tmp$Adiff[is.infinite(dia_tmp$Adiff)] = NA

  #For NLS: Remove NA entries! 
  daph_tmp = daph_tmp[!is.na(daph_tmp$Adiff) & !is.na(daph_tmp$N), ]
  dia_tmp = dia_tmp[!is.na(dia_tmp$Adiff) & !is.na(dia_tmp$N), ]

  cl_daph_plot = rbind( cl_daph_plot, daph_tmp )
  cl_dia_plot = rbind( cl_dia_plot, dia_tmp )

  #The basic MacArthur model is a linear consumption rate so just fit with a GLMM
  # cl_daph[[t]] = gam( Ndiff ~ +s(temperature,k=5)+s(mesocosm_id,bs="re"),family=Gamma(link='log'), data=daph_tmp)
  # cl_dia[[t]] = gam( Ndiff ~ +s(temperature,k=5)+s(mesocosm_id,bs="re"),family=Gamma(link='log'), data=dia_tmp)
  # plot(daph_tmp$algae_abundance, daph_tmp$N)
  # plot(dia_tmp$algae_abundance, dia_tmp$N)

  ########################
  #This works better!!! 
  #Fit the algal consumption rate instead: 
  # cl_daph[[t]] = gam( Adiff ~  N+s(mesocosm_id,bs="re"), family=binomial, data=daph_tmp)
  # cl_dia[[t]] = gam( Adiff ~  N +s(mesocosm_id,bs="re"), famly=binomial, data=dia_tmp)

  #Use NLS to fit a Type 2 (saturating) response
  alg2 = formula (Adiff ~  algae_start/(1+c1*exp(b1*N) ) )

  m1 = lm(I(log(algae_start/Adiff-1))~I(N), data = daph_tmp ) 
  tryCatch({ 
    cl_daph[[t]] = nls( formula= alg2, data = daph_tmp, 
    start=list(b1=(as.numeric(coef(m1)[2])), c1=exp(as.numeric(coef(m1)[1] ) ) ),
    control=nls.control(maxiter = 1000), trace=F ) 

    #Predicted values: 
    s=seq(min(daph_tmp$N),max(daph_tmp$N),1 )
    d_tmp = predict(cl_daph[[t]], list( N = s ) )
    daph_pred_tmp = data.frame( species = rspecies[1], temperature = temps[t], s=s, N_pred = d_tmp )
    cl_daph_pred = rbind(cl_daph_pred, daph_pred_tmp)
  }, error = function(e) {} ) 

  m1 = lm(I(log(algae_start/Adiff-1))~I(N), data = dia_tmp ) 
  tryCatch( {
    cl_dia[[t]] = nls( formula= alg2, data = dia_tmp, 
    start=list(b1=(as.numeric(coef(m1)[2])), c1=exp(as.numeric(coef(m1)[1] ) ) ),
    control=nls.control(maxiter = 1000), trace=F ) 

    #Predicted values: 
    s=seq(min(dia_tmp$N),max(dia_tmp$N),1 )
    d_tmp = predict(cl_dia[[t]], list( N = s ) )
    dia_pred_tmp = data.frame( species = rspecies[2], temperature = temps[t], s=s, N_pred = d_tmp )
    cl_dia_pred = rbind(cl_dia_pred, dia_pred_tmp)
  }, error = function(e) {} ) 


  #To plot in the loop: 
  plot((daph_tmp$N),(daph_tmp$Adiff),col="blue",
    ylim=c(min( c(min(daph_tmp$Adiff),min(dia_tmp$Adiff)) ), max(c(max(daph_tmp$Adiff),max(dia_tmp$Adiff) ) ) ),
    xlim=c(min( c(min(daph_tmp$N),min(dia_tmp$N)) ), max(c(max(daph_tmp$N),max(dia_tmp$N) ) ) )
     )
  s=seq(min(daph_tmp$N),max(daph_tmp$N),1 )
  lines(s, predict(cl_daph[[t]], list( N = s ) ), col = "blue")

  points((dia_tmp$N),(dia_tmp$Adiff),col="red")
  s=seq(min(dia_tmp$N),max(dia_tmp$N),1 )
  lines(s, predict(cl_dia[[t]], list( N = s ) ), col = "red")



  ############################  
  #The inverse relationship: 
  #alg2 = formula (Adiff ~ b/(1+b*N) )
  # alg2 = formula (N ~  exp(b1*Adiff+c1)  )
  # m1 = lm(I(log(N+1))~I((Adiff)), data = daph_tmp ) 
  # cl_daph[[t]] = nls( formula= alg2, data = daph_tmp, 
  #   start=list( b1=(as.numeric(coef(m1)[2])), c1=(as.numeric(coef(m1)[1] ) ) ),
  #   control=nls.control(maxiter = 1000), trace=T )

  # m1 = lm(I(log(N+1))~I((Adiff)), data = dia_tmp ) 
  # cl_dia[[t]] = nls( formula= alg2, data = dia_tmp, 
  #   start=list( b1=(as.numeric(coef(m1)[2])), c1=(as.numeric(coef(m1)[1] ) ) ),
  #   control=nls.control(maxiter = 1000), trace=T )

  # plot((daph_tmp$Adiff),(daph_tmp$N))
  # s=seq(min(daph_tmp$Adiff),max(daph_tmp$Adiff),10 )
  # lines(s, predict(cl_daph[[t]], list(Adiff= s ) ), col = "green")

  # plot((dia_tmp$Adiff),(dia_tmp$N ) )
  # s=seq(min(dia_tmp$Adiff),max(dia_tmp$Adiff),10 )
  # lines(s, predict(cl_dia[[t]], list(Adiff= s ) ), col = "green")


}

cl_plot = rbind(cl_daph_plot,cl_dia_plot)
cl_pred = rbind(cl_daph_pred,cl_dia_pred)

ggplot(cl_plot, aes(x = N, y =Adiff, color = species) ) + 
  geom_point( )+ facet_grid(temperature~species)+ 
  geom_line(data= cl_pred, mapping= aes(x = s, y =N_pred, color=species) )+
  facet_grid(temperature~species)+
  xlab("Zooplankton abundance ")+
  ylab("Algal consumption rate")+
  theme(strip.background = element_rect(colour=NA, fill=NA))
ggsave("./algal_consump_diaDaph.pdf", width = 8, height = 10)