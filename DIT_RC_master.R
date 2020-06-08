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
inv_day = 27 #The first day of attempted invasion
inv_end = 50#The last day of invasion conditions
no_reps = 18 #The number of replicated mesocosms total per resident/invader
nEspp=2

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
alg2 = formula (algae_abundance/algae_start ~ exp(c1*N+b1) )
alg2_prms= c("c1","b1" )
alg_lm=formula(I(log(algae_abundance/algae_start ))~N )

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
cR_prms= c("c1","b1" )
cR_lm=formula(I(log(N+1))~algae_abundance)

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
igr_prms = c("a1","b1","c1")

#3. and 4., The interspecific invasion complements to 1. and 2.: 
#3. 
#igrij = formula (N_inv ~  (b1*exp(c1*day_n ) ) )

#4. ***    ***
igrij = formula (N_inv ~  a1/(1+b1*exp(c1*day_n ) ) )
igrij_prms = c("a1","b1","c1")
#igrij_prms = c("b1","c1") #Treat a1 as fixed 


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
#f_lvii = Ndiff ~ N 
# Interspecific competition models
#f_lvij = Ndiff ~ N 

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

#2b. LV model without fitting ri, constrained intercept (using offset=1: see below)
#f_lvii = Ndiff ~ N +0
f_lvii = Ndiff ~ N 


#3. *** Leslie-Gower ***
#f_lvii = formula ( Ndiff ~  ri/(1+aii*N) )

####################################
#Interspecific competition models

#1. Lotka Volterra, aij only 
#f_lvij = formula ( Ndiff ~  ri*(1-aij*N_res) )

#1b.. LV model without fitting ri, constrained intercept (using offset=1: see below)
f_lvij = Ndiff ~ N_res +0

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

for(t in 1:ntemps) { 
	#=============================================================================
	#Resource consumption and intraspecific competition
	#=============================================================================
	#Make the single species data set by isolating mesocosms and portions of time
	#series where only a single species is present. Use both the monoculture data 
	#and the data for the resident pre-invasion for resource consumption fits, 
	#intraspecific competition, intrinsic growth
	#=============================================================================
	#Get the data for each temperature. 
	daph_tmp = get_new_mono(m1_data_long, rspecies[1], temps[t], inv_day )
	dia_tmp = get_new_mono(m1_data_long, rspecies[2], temps[t], inv_day )


	#Resource consumption rates:
	m1 = lm(alg_lm, daph_tmp) #Use for starting values for NLS in prm_start
	cl_daph[[t]] = get_mod_fit( mod_data =daph_tmp, mod_fit = alg2, mod_prms = alg2_prms,
					prm_start = c((as.numeric(coef(m1)[2])), (as.numeric(coef(m1)[1] ) ) ), 
					mod_x = "N", lm_mod = alg_lm  ) 
	if(!is.null(cl_daph[[t]])) { 
	daph_pred_tmp = data.frame( species = rspecies[1], temperature = temps[t], 
					N=cl_daph[[t]]$new_fit$N, N_pred = cl_daph[[t]]$new_fit$N_pred )
    cl_daph_pred = rbind(cl_daph_pred, daph_pred_tmp)
	}

	m1 = lm(alg_lm, dia_tmp)
	cl_dia[[t]] = get_mod_fit( mod_data =dia_tmp, mod_fit = alg2, mod_prms = alg2_prms,
					prm_start = c((as.numeric(coef(m1)[2])), (as.numeric(coef(m1)[1] ) ) ), 
					mod_x = "N", lm_mod = alg_lm  ) 
	if(!is.null(cl_dia[[t]])) { 
	dia_pred_tmp = data.frame( species = rspecies[2], temperature = temps[t], 
					N=cl_dia[[t]]$new_fit$N, N_pred = cl_dia[[t]]$new_fit$N_pred )	 
    cl_dia_pred = rbind(cl_dia_pred, dia_pred_tmp)
	}

	#Intrinsic growth f(algae)
	m1 = lm(cR_lm, daph_tmp)
	cR_daph[[t]] =get_mod_fit( mod_data =daph_tmp, mod_fit = cR, mod_prms = cR_prms,
					prm_start = c((as.numeric(coef(m1)[2])), (as.numeric(coef(m1)[1] ) ) ), 
					mod_x = "algae_abundance", lm_mod = cR_lm  ) 
	if(!is.null(cR_daph[[t]])) { 
	daph_pred_tmp = data.frame( species = rspecies[1], temperature = temps[t], 
					algae_abundance=cR_daph[[t]]$new_fit$algae_abundance, 
					N_pred = cR_daph[[t]]$new_fit$N_pred )
	cR_daph_pred = rbind(cR_daph_pred, daph_pred_tmp)
	}
	
	m1 = lm(cR_lm, dia_tmp)
	cR_dia[[t]] = get_mod_fit( mod_data =dia_tmp, mod_fit = cR, mod_prms = cR_prms,
					prm_start = c((as.numeric(coef(m1)[2])), (as.numeric(coef(m1)[1] ) ) ), 
					mod_x = "algae_abundance", lm_mod = cR_lm  )
	if(!is.null(cR_dia[[t]])) { 
	dia_pred_tmp = data.frame( species = rspecies[2], temperature = temps[t], 
					algae_abundance=cR_dia[[t]]$new_fit$algae_abundance, 
					N_pred = cR_dia[[t]]$new_fit$N_pred )
	cR_dia_pred = rbind(cR_dia_pred, dia_pred_tmp)
	}

	#Intrinsic growth rate f(time)
	#Adjust the days to make both data sets line up
  	daph_tmpb=daph_tmp
  	daph_tmpb$day_n[daph_tmpb$day_n>inv_day] = daph_tmpb$day_n[daph_tmpb$day_n>inv_day]-inv_day 
  	
  	dia_tmpb=dia_tmp
  	dia_tmpb$day_n[dia_tmpb$day_n>inv_day] = dia_tmpb$day_n[dia_tmpb$day_n>inv_day]-inv_day 

	test1=NULL
	test1 = get_mod_fit( mod_data =daph_tmpb, mod_fit = igr, mod_prms = igr_prms,
					prm_start = c(max(daph_tmpb$N,na.rm=T), 1, 0.05 ), mod_x = "day_n")
 	if(!is.null(test1)) { 
 	igr_daph[[t]] =test1
	igr_daph_tmp = data.frame( species = rspecies[1], temperature = temps[t],
		 			day_n=igr_daph[[t]]$new_fit$day_n,  
		 			N = igr_daph[[t]]$new_fit$N_pred )
    igr_daph_pred = rbind(igr_daph_pred, igr_daph_tmp)
	}
	
	test1 = NULL
	test1= get_mod_fit( mod_data =dia_tmpb, mod_fit = igr, mod_prms = igr_prms,
					prm_start = c(max(dia_tmpb$N,na.rm=T), 1, 0.05 ), mod_x = "day_n")
	if(!is.null(test1)) { 
	igr_dia[[t]] =test1
	igr_dia_tmp = data.frame( species = rspecies[2], temperature = temps[t],
		 			day_n=igr_dia[[t]]$new_fit$day_n,  
		 			N = igr_dia[[t]]$new_fit$N_pred )
    igr_dia_pred = rbind(igr_dia_pred, igr_dia_tmp)
	}
	
	#Intraspecific competition, when fitting to the above model fits igr_specis 
	#instead of directly to the data. 
	#Add Ndiff
	# igr_daph_tmp$Ndiff = (igr_daph_tmp$N-lag(igr_daph_tmp$N))/
 #                        (igr_daph_tmp$day_n-lag(igr_daph_tmp$day_n))*1/igr_daph_tmp$N
	# igr_dia_tmp$Ndiff = (igr_dia_tmp$N-lag(igr_dia_tmp$N))/
 #                        (igr_dia_tmp$day_n-lag(igr_dia_tmp$day_n))*1/igr_dia_tmp$N
	# lvii_daph[[t]] = lm(f_lvii, igr_daph_tmp )
	# lvii_dia[[t]] = lm(f_lvii,igr_dia_tmp )

	#Intraspecific competition, when fitting only to the "invasion" portion of data
 	# daph_tmpb_inv = subset(daph_tmpb, day_n <= (inv_end-inv_day) )
 	# dia_tmpb_inv = subset(dia_tmpb, day_n <= (inv_end-inv_day) )

 	#Assume that only positive increments reflect growth? 
	daph_tmpb_inv  = subset(daph_tmpb , Ndiff>0)
	dia_tmpb_inv= subset(dia_tmpb, Ndiff>0)

	lvii_daph[[t]] = lm(f_lvii, daph_tmpb_inv)#, offset = rep(1,length(daph_tmpb_inv$N))  )
	lvii_daph[[t]]$new_fit$N = seq(1, max(daph_tmpb_inv$N,na.rm=T) ) 
	lvii_daph[[t]]$new_fit$N_pred = lvii_daph[[t]]$new_fit$N* coef(lvii_daph[[t]])[2] + coef(lvii_daph[[t]])[1]
	#lvii_daph[[t]]$new_fit$N_pred = lvii_daph[[t]]$new_fit$N* coef(lvii_daph[[t]])+1
	lvii_daph_tmp = data.frame( species = rspecies[1], temperature = temps[t],
		 			N=lvii_daph[[t]]$new_fit$N,  
		 			Ndiff= lvii_daph[[t]]$new_fit$N_pred )
    lvii_daph_pred = rbind(lvii_daph_pred, lvii_daph_tmp)

	lvii_dia[[t]] = lm(f_lvii,dia_tmpb_inv)#, offset = rep(1,length(dia_tmpb_inv$N))  )
	lvii_dia[[t]]$new_fit$N = seq(1, max(dia_tmpb_inv$N,na.rm=T) ) 
	lvii_dia[[t]]$new_fit$N_pred = lvii_dia[[t]]$new_fit$N* coef(lvii_dia[[t]])[2] + coef(lvii_dia[[t]])[1]
	#lvii_dia[[t]]$new_fit$N_pred = lvii_dia[[t]]$new_fit$N* coef(lvii_dia[[t]]) + 1
	lvii_dia_tmp = data.frame( species = rspecies[2], temperature = temps[t],
		 			N=lvii_dia[[t]]$new_fit$N,  
		 			Ndiff= lvii_dia[[t]]$new_fit$N_pred )
    lvii_dia_pred = rbind(lvii_dia_pred, lvii_dia_tmp)
	#=============================================================================
	#Interspecific competition
	#=============================================================================
	#First, try fitting a logistic model to the growth rate: 
	daph_inv_tmp = get_new_inv(m1_data_long, rspecies[1], temps[t], inv_day, inv_end)
	dia_inv_tmp = get_new_inv(m1_data_long, rspecies[2], temps[t], inv_day, inv_end )

	test1=NULL
	test1 = get_mod_fit( mod_data =daph_inv_tmp, mod_fit = igrij, mod_prms = igrij_prms,
					prm_start = c(max(daph_inv_tmp$N_inv,na.rm=T), 2, -0.05 ), mod_x = "day_n"
					)
	if(!is.null(test1)) { 
	igrij_daph[[t]]=test1
	igrij_daph_tmp = data.frame( species = rspecies[1], temperature = temps[t],
		 			day_n=igrij_daph[[t]]$new_fit$day_n,  
		 			N= igrij_daph[[t]]$new_fit$N_pred )
    igrij_daph_pred = rbind(igrij_daph_pred, igrij_daph_tmp)
	}
	
	test1=NULL
	test1 = get_mod_fit( mod_data =dia_inv_tmp, mod_fit = igrij, mod_prms = igrij_prms,
					prm_start = c( max(dia_inv_tmp$N_inv,na.rm=T), 2, -0.05 ), mod_x = "day_n"
					)
	if(!is.null(igrij_dia[[t]])) { 
	igrij_dia[[t]] = test1
	igrij_dia_tmp = data.frame( species = rspecies[2], temperature = temps[t],
		 			day_n=igrij_dia[[t]]$new_fit$day_n,  
		 			N= igrij_dia[[t]]$new_fit$N_pred )
    igrij_dia_pred = rbind(igrij_dia_pred, igrij_dia_tmp)
	}

	#Interspecific competition, when fitting to the above model fits igrij_specis 
	#instead of directly to the data. 
	#Add Ndiff
	# igrij_daph_tmp$Ndiff = (igrij_daph_tmp$N-lag(igrij_daph_tmp$N))/
 #                        (igrij_daph_tmp$day_n-lag(igrij_daph_tmp$day_n))*1/igrij_daph_tmp$N
	# igrij_dia_tmp$Ndiff = (igrij_dia_tmp$N-lag(igrij_dia_tmp$N))/
 #                        (igrij_dia_tmp$day_n-lag(igrij_dia_tmp$day_n))*1/igrij_dia_tmp$N
	# lvij_daph[[t]] = lm(f_lvij, igrij_daph_tmp )
	# lvij_dia[[t]] = lm(f_lvij,igrij_dia_tmp )

	####Fit a linear model with a fixed intercept to the invasion growth rate:
	daph_inv_tmp = get_new_inv(m1_data_long, rspecies[1], temps[t], inv_day, inv_end)
	dia_inv_tmp = get_new_inv(m1_data_long, rspecies[2], temps[t], inv_day, inv_end )

	#Assume that only positive increments reflect growth? 
	daph_inv_tmp = subset(daph_inv_tmp, Ndiff>0)
	dia_inv_tmp = subset(dia_inv_tmp, Ndiff>0)

	lvij_daph[[t]] = lm(f_lvij, daph_inv_tmp, offset = rep(1,length(daph_inv_tmp$N_res))  )
	lvij_daph[[t]]$new_fit$N_res = seq(1, max(daph_inv_tmp$N_res,na.rm=T) ) 
	lvij_daph[[t]]$new_fit$N_pred = lvij_daph[[t]]$new_fit$N_res* coef(lvij_daph[[t]]) + 1
	lvij_daph_tmp = data.frame( species = rspecies[1], temperature = temps[t],
		 			N_res=lvij_daph[[t]]$new_fit$N_res,  
		 			Ndiff= lvij_daph[[t]]$new_fit$N_pred )
    lvij_daph_pred = rbind(lvij_daph_pred, lvij_daph_tmp)
	
	lvij_dia[[t]] = lm(f_lvij,dia_inv_tmp, offset = rep(1,length(dia_inv_tmp$N_res))  )
	lvij_dia[[t]]$new_fit$N_res = seq(1, max(dia_inv_tmp$N_res,na.rm=T) ) 
	lvij_dia[[t]]$new_fit$N_pred = lvij_dia[[t]]$new_fit$N_res* coef(lvij_dia[[t]]) + 1
	lvij_dia_tmp = data.frame( species = rspecies[2], temperature = temps[t],
		 			N_res=lvij_dia[[t]]$new_fit$N_res,  
		 			Ndiff= lvij_dia[[t]]$new_fit$N_pred )
    lvij_dia_pred = rbind(lvij_dia_pred, lvij_dia_tmp)

	#Output
	print(mean(as.data.frame(subset(daph_inv_tmp,day_n <=34))$Ndiff,na.rm=T) )
	print(mean(as.data.frame(subset(dia_inv_tmp,day_n <=34))$Ndiff,na.rm=T) )

	#Keep building these for plotting: 
	cl_daph_plot = rbind( cl_daph_plot, daph_tmp )
	cl_dia_plot = rbind( cl_dia_plot, dia_tmp )


}

#Since invasion treatment 6 was nonexistant: 
lvij_dia[[6]]=lvij_daph[[6]]
lvij_dia[[6]]$coefficients[1] = 1E-5 #Make it effectively nothing 

#Consumption functions
cl_plot = rbind(cl_daph_plot,cl_dia_plot)
cl_pred = rbind(cl_daph_pred,cl_dia_pred)

#Growth rate functions
cR_pred = rbind(cR_daph_pred,cR_dia_pred)

#competition
lvii_pred = rbind(lvii_daph_pred, lvii_dia_pred)
lvij_pred = rbind(lvij_daph_pred, lvij_dia_pred)


#=============================================================================
#Plotting for Stage 2: 
#=============================================================================
#Consumption functions: 
#Separate panels
#Pick the appropriate first line: 
#ggplot(cl_plot, aes(x = N, y =algae_abundance, color = species) )+  #2
ggplot(cl_plot, aes(x = N, y =algae_abundance/algae_start, color = species) )+ #3
#ggplot(cl_plot, aes(x = N, y =Adiff, color = species) )+ #3
  geom_point( )+ facet_grid(temperature~species)+ 
  geom_line(data= cl_pred, mapping= aes(x = N, y =N_pred, color=species) ) + 
  facet_grid(temperature~species)+
  xlab("Zooplankton abundance ")+
  ylab("Algal consumption rate")+
  theme(strip.background = element_rect(colour=NA, fill=NA))
#ggsave("./algal_consump3_diaDaph.pdf", width = 8, height = 10)

#=============================================================================
#Intrinsic growth rate functions: 
#Separate panels
#Pick the appropriate first line: 
#ggplot(cl_plot, aes(x = day_n, y =N, color = species) ) + #1. 
#ggplot(cl_plot, aes(x = Adiff, y =N, color = species) ) + #2. 
ggplot(cl_plot, aes(x = algae_abundance, y =N, color = species) ) + #2. 
  geom_point( )+ facet_grid(temperature~species)+ 
  geom_line(data= cR_pred, mapping= aes(x = algae_abundance, y =N_pred, color=species) )+
  facet_grid(temperature~species)+ #xlim( min(cl_plot$Adiff), max(cl_plot$Adiff))+
  xlab("Zooplankton abundance ")+
  ylab("Time")+
  theme(strip.background = element_rect(colour=NA, fill=NA))
#ggsave("./intrinsicR_diaDaph2.pdf", width = 8, height = 10)

#=============================================================================
#Competition models 
#Separate panels
#Pick the appropriate first line: 
#ggplot(cl_plot, aes(x = day_n, y =N, color = species) ) + #1. 
#ggplot(cl_plot, aes(x = Adiff, y =N, color = species) ) + #2. 
ggplot(cl_plot, aes(x = N, y =Ndiff, color = species) ) + #2. 
  geom_point( )+ 
  geom_line(data= lvii_pred, mapping= aes(x = N, y =Ndiff, color=species) )+
  facet_grid(temperature~species) + #xlim( min(cl_plot$Adiff), max(cl_plot$Adiff))+
  xlab("Zooplankton abundance ")+
  ylab("Growth rate")+  xlim(0,60)+
  theme(strip.background = element_rect(colour=NA, fill=NA))
#ggsave("./intrinsicR_diaDaph2.pdf", width = 8, height = 10)

ggplot(cl_plot, aes(x = N, y =Ndiff, color = species) ) + #2. 
  geom_point( )+
  geom_line(data= lvij_pred, mapping= aes(x = N_res, y =Ndiff, color=species) )+
  facet_grid(temperature~species) + #xlim( min(cl_plot$Adiff), max(cl_plot$Adiff))+
  xlab("Zooplankton abundance ") +
  ylab("Growth rate")+  xlim(0,60)+
  theme(strip.background = element_rect(colour=NA, fill=NA))

#=============================================================================
#STAGE 3: Theoretical simulations of population dynamics. 
# Generate population data using a theoretical Resource-Consumer (RC) model. 
# Check the file [[]] for different modeling approaches. 
#
# The current favorite approach is to use the basic MacArthur RC model with 
# parameters fitted from STAGE 2 at each temperature level. Since the
# experiment uses only one resource -- algae -- but still allows coexistence, 
# a fictitious 2nd resource is modeled which is more heavily consumed by 
# Diaphanosoma. The goal is to make the consumption rates of both resources 
# to match up with the phenomenological competition coefficients, especially 
# where coexistence is possible at temps[[5]].  
#=============================================================================
#=============================================================================
#Set parameters for simulation
#=============================================================================
#Length and time steps of each model run
tend = 100
delta1 = 0.01
tl=tend/delta1+1

###The maximum block depth for dynamic info metrics (larger is more accurate, but
#slower and could cause crashing if too large)
k= 2

###Build a series of scenarios going from simple to more complex dynamics
#Number of food webs to generate
nwebs = length(temps)

###Output of each web
out1 = vector("list",nwebs)
#Invasion scenario
out_inv1 = vector("list",nwebs)

#Competition coefficients.
aii_all = matrix(0,nwebs,2)
aij_all = matrix(0,nwebs,2)
aii_rc = aii_all
aij_rc = aii_all


###Species numbers
#Assume 2 trophic levels unless otherwise specified.
nRsp = 2 #Algae
nCsp = 2 #Spp 1 is Daphnia, Spp 2 is Diaphanosoma
nPsp = 1 #This is actually 0 --> Just a dummy predator
nspp = nRsp+nCsp+nPsp

#Phenomenological coexistence conditions
rho1 = matrix(0,6,2)
the1 = matrix(0,6,2)
kappa1 = matrix(0,6,2)

#=============================================================================
# Outer loop. First run to equilibrate the population dynamics
#=============================================================================
par(mfrow=c(6,1),oma = c(5,4,0,0) + 0.1,mar = c(0,0,1,1) + 0.1)

for (w in 1:nwebs){ 
	#Main list of species parameters for the model
	spp_prms = NULL

	#############################
	#Resources: 
	#This isn't the actual carrying capacity, but the amount that's added every 
	#time step and we assume it's the maximum possible. 
	spp_prms$Kr = matrix(10E6, nRsp, 1) #carrying capacity -- Constant algal addition

	#Random algal/consumer fluctuations. 
	#Algae
	c1 = 1 #10E6
	amp1 = 1E-5#0.25 # #100000 #1/exp(1) 
	spp_prms

	############################
	#Consumers: 
	#
	###Random fluctuations in the consumers
	daph_tmp = subset(m1_data_long, species == "daphnia" & temperature == temps[[w]]) #Daphnia at temp
	dia_tmp = subset(m1_data_long, species == "diaphanosoma" & temperature == temps[[w]] ) #Daphnia at temp

	#Approximate a carrying capacity to get the right magnitude of fluctuations: 
	daph_tmp_a=daph_tmp[daph_tmp$day_n<50, ]
	dia_tmp_a=dia_tmp[dia_tmp$day_n<50, ]
	spp_prms$Kc = matrix(c(mean(daph_tmp_a$N,na.rm=T),2*mean(dia_tmp_a$N,na.rm=T)), nCsp, 1) #carrying capacities, approximately matching data
	Kc_sd = matrix(c(sqrt(var(daph_tmp_a$N,na.rm=T)),sqrt(var(dia_tmp_a$N,na.rm=T))), nCsp, 1)
	
	c2 = 1
	amp2 = 1E-5#mean(spp_prms$Kc/Kc_sd) #  #Take an average coefficient of variation. 
	res_R = c(amp1,c1,amp2,c2)

	
	###Consumption rates: 
	#The coefficient that we want could be in one of two places, depending on whether the 
	#growth curve (in resource_fits1.R) was fit by an LM or NLS. 
	if(class(cl_daph[[w]]$fit_mod) == "nls" ){ daph_C = coef(cl_daph[[w]]$fit_mod)[1]} else {daph_C = coef(cl_daph[[w]]$fit_mod)[2] }
	if(class(cl_dia[[w]]$fit_mod) == "nls" ){ dia_C = coef(cl_dia[[w]]$fit_mod)[1]} else {dia_C = coef(cl_dia[[w]]$fit_mod)[2] }
	spp_prms$cC = t(matrix(abs(c(daph_C,dia_C)),nRsp,nCsp ))

	#Intrinsic growth rates: 
	#The coefficient that we want could be in one of two places, depending on whether the 
	#growth curve (in resource_fits1.R) was fit by an LM or NLS. 
	if(class(cR_daph[[w]]$fit_mod) == "nls" ){ daph_C = coef(cR_daph[[w]]$fit_mod)[1]} else {daph_C = coef(cR_daph[[w]]$fit_mod)[2] }
	if(class(cR_dia[[w]]$fit_mod) == "nls" ){ dia_C = coef(cR_dia[[w]]$fit_mod)[1]} else {dia_C = coef(cR_dia[[w]]$fit_mod)[2] }
	spp_prms$rC = t(matrix(abs(c(daph_C,dia_C)),nRsp,nCsp ) )#/spp_prms$cC 

	#The competition coefficients in the MacArthur RC model are given by cil*cjl*ri*Kl/Rl: 
	#cil/cjl are the consumption rates of resource l by species i/j, ri is species i's growth
	#rate based on consumption and Kl/Rl is the carrying capacity/growth rate of resource l. 
	#We use this relationship and set these alphas = phenomenological alphas then solve for Rl so
	#coefficients match.  
	aii_rc[w,] = (10E6*spp_prms$cC[1,]^2*spp_prms$rC[1,])
	aij_rc[w,] = (10E6*spp_prms$cC[1,]*spp_prms$cC[1,2:1]*spp_prms$rC[1,])

	#Phenomenological alphas.
  	aii_all[w,1] = abs(coef(lvii_daph[[w]])[2])
  	aii_all[w,2] = abs(coef(lvii_dia[[w]])[2])
  	aij_all[w,1] = abs(coef(lvij_daph[[w]])[1])
  	aij_all[w,2] = abs(coef(lvij_dia[[w]])[1])

  	#if(w<5){aij_all[w,2]=2*aij_all[w,2]}

	# The trick now is to figure out how to adjust the RC alphas to match the phenomenological
	# alphas by determining an appropriate R1 and the properties of a fictitious resource 2.  
	# The approach here starts by assuming aii (for Daphnia) is correct and building from there.

	#1.) Determine the R1 that makes aii correct and assume C1 only benefits from resource 1: 
	spp_prms$rR =aii_rc[w,]/aii_all[w,]
	spp_prms$rC[2,1] =0 

	#2.) Correct cC[1,2] (Diaphanosoma) based on aij
	spp_prms$cC[1,2] = spp_prms$rR[1]/spp_prms$Kr[1] * aij_all[w,1]/(spp_prms$cC[1,1]*spp_prms$rC[1,1])

	#3.) This leads to a new ajj for Diaphanosoma:
	aii_rc[w,2]=(10E6*spp_prms$cC[1,2]^2*spp_prms$rC[1,2])/spp_prms$rR[1]

	#4.) Which differs from the phenomenological ajj by:
	ajj_diff = abs(aii_rc[w,2] - aii_all[w,2] )

	#5.) Make this the basis for a fictitious second resource that makes ajj correct now. 
	# Adjust C2's consumption so that species fitnesses are very close in value. 
	# *Assume C2's consumption rate of R2 is the same as R1
	# *Assume the relative fitness ratio is equivalent to the fitted intrinsic growth rates
	# *Make rR2 = 1 and assume that Kr2 is the free parameter. 
	#spp_prms$rR[2] = spp_prms$rR[1]/10 
	spp_prms$rC[2,2] = spp_prms$rC[1,2]
	if (!is.null(igr_daph[[w]]$fit_mod) & !is.null(igr_dia[[w]]$fit_mod) ){ 
		fit_rat= coef(igr_dia[[w]]$fit_mod)[["c1"]]/coef(igr_daph[[w]]$fit_mod)[["c1"]] 
	} else { fit_rat = 0.75}
	fit_tmp1 = spp_prms$cC[1,1]*spp_prms$rC[1,1]*(spp_prms$Kr[1])-1
	fit_tmp2 = fit_tmp1*fit_rat
	spp_prms$cC[2,2] = (fit_tmp2 + 1) / (spp_prms$rC[2,2]*spp_prms$Kr[2] )
	#spp_prms$cC[2,2] = spp_prms$cC[1,1] -spp_prms$cC[1,2]
	spp_prms$Kr[2] = (spp_prms$rR[2]*ajj_diff)/( spp_prms$cC[2,2]^2*spp_prms$rC[2,2] )

	# *Make K2 = K1 and assume that rR2 is the free parameter. This keeps C1 and C2 
	# about equal in fitness. 
	# *Assume that the consumption of R2 makes C2's total consumption of resources
	# about equivalent to C1. 
	# NOTE: This version is potentially problematic when rR2 is very small because
	# it then violates the timescale separation that makes all of this work.  
	# spp_prms$Kr[2] = spp_prms$Kr[1]
	# spp_prms$cC[2,2] = spp_prms$cC[1,1] -spp_prms$cC[1,2]
	# spp_prms$rR[2] = ( spp_prms$cC[2,2]^2*spp_prms$rC[2,2] )/ajj_diff



	#6.) This creates a new aji for Dia that must be corrected for by giving Daphnia a positive
	# consumption rate of R2 (which then must be removed from Daphnia's growth equation by making
	# its rC2=0 )
	spp_prms$cC[2,1]= (aij_all[w,2]-(spp_prms$cC[1,2]*spp_prms$cC[1,1]*spp_prms$rC[1,2]*spp_prms$Kr[1]/spp_prms$rR[1]) )*
	spp_prms$rR[2]/(spp_prms$Kr[2]*spp_prms$cC[2,2]*spp_prms$rC[2,2])  
	
	#7.) Check: 
	aii_check = t(spp_prms$cC^2*spp_prms$rC)%*%(spp_prms$Kr/spp_prms$rR) 
	aij_check = t(spp_prms$cC*spp_prms$cC[,2:1]*spp_prms$rC)%*%(spp_prms$Kr/spp_prms$rR) 

	#Misc.
	spp_prms$eFc = matrix(1,nCsp,nRsp) # just make the efficiency for everything 1 for now
	spp_prms$monoculture = matrix(1, nCsp, 1) #matrix(rnorm(nCsp,0.6,0.1), nCsp, 1) #mortality rates	
	#To give species their actual growth rates: 
	# spp_prms$eFc[1] = abs(coef(igr_daph[[w]]$fit_mod)["c1"])
	# spp_prms$eFc[2] = abs(coef(igr_dia[[w]]$fit_mod)["c1"])

	#This parameterizes the CR model to match to the LV in the standard form of 
	# (1 - competition), where 1 = k in Chesson 1990. 
	 spp_prms$muC=t(spp_prms$cC*spp_prms$rC)%*%(spp_prms$Kr)-1
	# aii_intercept = c( abs(coef(lvii_daph[[w]])[1]), abs(coef(lvii_dia[[w]])[1]))
	# spp_prms$muC=t(spp_prms$cC*spp_prms$rC)%*%(spp_prms$Kr)-aii_intercept

	spp_prms$aii = aii_all[w,]
	spp_prms$aij = aij_all[w,]


	########Predators: These are just dummy variables for now
	spp_prms$rP =  matrix(0.0, 1, 1) #matrix(rnorm(nPsp,0.5,0), nPsp, 1) #intrisic growth
	spp_prms$eFp = matrix(1,1,nCsp) # just make the efficiency for everything 1 for now
	spp_prms$muP = matrix(0.0, 1, 1)#matrix(rnorm(nPsp,0.6,0), nPsp, 1)  #mortality rates
	#Consumption rates: 
	spp_prms$cP = matrix(c(0.0,0.0),nCsp,1)
	# spp_prms$aii = matrix(c(0.10,0.08),1,2)
	# spp_prms$aij = matrix(c(0.10,0.08),1,2)

	#Phenomenological coexistence conditions
	rho1[,1] = aij_all[,1]/(sqrt(aii_all[,1]*aii_all[,2]) )
	rho1[,2] = aij_all[,2]/(sqrt(aii_all[,1]*aii_all[,2]))
	the1[,1] = sqrt(aii_all[,1]/aii_all[,2])
	the1[,2] = sqrt(aii_all[,2]/aii_all[,1])
	kappa1[w,] = spp_prms$muC + 1
	#kappa1[w,] = spp_prms$muC + aii_intercept

	#=============================================================================
	# This function gives: 
	# out 		The time series for of population growth for each species in the web
	#			This can be set to just give the final 2 time steps of the web with
	#			"final = TRUE"
	# spp_prms	The parameters of all species in the food web
	#=============================================================================
	#
	winit = matrix(c(spp_prms$Kr[1], spp_prms$Kr[2], 4,4,0))
	tryCatch( {out1[w] = list(food_web_dynamics (spp_list = c(nRsp,nCsp,nPsp), spp_prms = spp_prms, 
		tend, delta1, winit = winit, res_R = res_R,final = FALSE ))}, error = function(e){}) 

	# print(out1[[w]]$out[tl,])

	# #Look at basic plots on the fly: 
	# plot(out1[[w]]$out[,5],t="l",ylim=c(0,max(out1[[w]]$out[tl,4:5])))
	# points(dia_tmp$day_n, dia_tmp$N)
	# lines( out1[[w]]$out[,4],col="red",t="l")
	# points(daph_tmp$day_n,daph_tmp$N, col="red")

	#=============================================================================
	# Inner loop: Mutual invasion:remove each species and track the dynamics
	#=============================================================================
	a_temp = NULL
	for (s in 3:4){

		out_temp =NULL
		out_temp2 =NULL
		inv_spp = s
		#winit =  out1[[w]]$out[tl,2:(nspp+1)]
		winit = matrix(c(spp_prms$Kr[1], spp_prms$Kr[2], 1/spp_prms$aii[1],1/spp_prms$aii[2],0))
		winit[inv_spp] = 0

		#Equilibrate new community
		tryCatch( {out_temp = (food_web_dynamics (spp_list = c(nRsp,nCsp,nPsp), spp_prms = spp_prms, 
		tend, delta1, winit = winit, res_R = res_R,final = FALSE ))}, error = function(e){}) 

		#Invade
		#ti = which(out_temp$out[ ,nRsp+3] == max(out_temp$out[,nRsp+3] ) )
		ti = tl
		winit =  out_temp$out[ti,2:(nspp+1)]
		winit[inv_spp] = 4

		#Invade new community
		tryCatch( {out_temp2 = (food_web_dynamics (spp_list = c(nRsp,nCsp,nPsp), spp_prms = spp_prms, 
		tend, delta1, winit = winit, res_R = res_R,final = FALSE ))}, error = function(e){}) 
		
		a_temp =rbind(a_temp, out_temp$out, out_temp2$out)

	}

	out_inv1 [[w]] = a_temp
	#Look at basic plots on the fly: 
	# plot(out_inv1 [[w]][30000:31000,1], out_inv1 [[w]][30000:31000,4], t="l", col="red")
	# points(daph_tmp$day_n,daph_tmp$N, col="red")
	# lines(out_inv1 [[w]][10000:11000,1], out_inv1 [[w]][10000:11000,5])	
	# points(dia_tmp$day_n, dia_tmp$N)

	#Look at basic plots on the fly: 
	plot(out_inv1 [[w]][,5],t="l",ylim=c(0,max(out_inv1 [[w]][,4:5])))
	lines( out_inv1 [[w]][,4],col="red")

}



#=============================================================================
# STAGE 4: Dynamic Information Theoretic analysis of real and simulated data. 
#
# Take both population time-series counts output by the ODEs and the actual
# time series data and calculate Excess Entropy, Active Information Storage, 
# and Transfer Entropy.
# Each quantity is calculated at both the average and local level.  
#=============================================================================
#Variables: 
#Set k, the block size: 
k=2
f1R=1 #scaling term for converting fractional populations in simulations
f1S=1 #scaling term for converting fractional populations in simulations


#=============================================================================
# Variables for data collection
#=============================================================================
#Real data: 
#=============================================================================
#Main data frame per mesocosm
out1R = vector("list",nmesos)

#Dynamic information metrics calculated from the (discretized) time series 
di_webR =  vector("list",nmesos)

#Track the average transfer entropy and separable information between each pair of 
#species as a way to build a network of information flow through the network. 
te_webR =  vector("list",nmesos)
si_webR = vector("list",nmesos) 

#The ensemble version of the AI
aiE_webR = vector("list",nmesos)

#=============================================================================
#Sim data: 
#=============================================================================
#Dynamic information metrics calculated from the (discretized) time series 
di_webS = vector("list",nwebs)

#Track the average transfer entropy and separable information between each pair of 
#species as a way to build a network of information flow through the network. 
te_webS = vector("list",nwebs)
si_webS =vector("list",nwebs) 

#The ensemble version of the AI
aiE_webS = vector("list",nwebs)


#=============================================================================
#Real data: Loop over all mesocosms
#=============================================================================
for (i in 1:nmesos) {

	mesosi = mesos[i]
	out1R[[i]] = get_info_ts(m1_data_long, mesosi )
	ti = which(temps == unique(out1R[[i]]$temperature)) #Get the temperature index
	#Take the two-species time series
	pop_ts = floor(f1R*out1R[[i]][c("N_res","N_inv")])

	nt1 = 1
	nt2 = dim(pop_ts)[1]
	if(nt2 <=k){ k = 1}

	if(nt1 != nt2) { 

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
	di_webR[i] = list(get_info_dynamics(pop_ts = pop_ts , k=k,with_blocks=FALSE))

	## This code takes the population time-series counts output by the ODEs and 
	## calculates the average Transfer Entropy from each species to every other 
	## species. The goal is to get an overview of the major information pathways 
	## in the web.   
	#=============================================================================
	# This function gives:
	# te_web    Average transfer entropy per species as a pairwise matrix
	#=============================================================================
	te_webR[i] = list( get_te_web( pop_ts = pop_ts, k=k) )

	#=============================================================================
	# This function gives:
	# aiE_web    The AI of the entire ensemble, treated as a single time series. 
	#=============================================================================
	aiE_webR[i] = list( get_ais (  series1 = pop_ts, k=k, ensemble = FALSE)    )

	#=============================================================================  
	#Build these back out into a data frame that includes all of the mesocosm,
	#treatment, and species information. 

	#Add the DIT to the data frames. There will be 11 new columns. 
	ncnames = c("N_res","N_inv","aiE","te1","te2","ee1","ee2","ai1","ai2","si1","si2",
	  "res_spp", "inv_spp", "res_K", "inv_K","res_aii", "res_aij", "inv_aii","inv_aij")
	DIT_tmp = data.frame(matrix(0,nt2,19) )
	colnames(DIT_tmp) = ncnames
	DIT_tmp[,1:2] = as.matrix(pop_ts)/f1R
	DIT_tmp[(k+1):nt2,3] = aiE_webR[[i]]$local
	DIT_tmp[(k+1):nt2,4:5] = di_webR[[i]]$te_local
	DIT_tmp[(k*2):nt2,6:7] = di_webR[[i]]$ee_local
	DIT_tmp[(k+1):nt2,8:9] = di_webR[[i]]$ai_local
	DIT_tmp[(k+1):nt2,10:11] = di_webR[[i]]$si_local
	DIT_tmp[,12] = factor( levels = levels(rspecies))
	DIT_tmp[,13] = factor( levels = levels(rspecies))

	#Add some of the population parameters to the data.frame: 
	t1 = as.character(unique(out1R[[i]]$species)[1])
	tm1 = which(temps == unique(out1R[[i]]$temperature))
	if( out1R[[i]]$invade_monoculture[1] == "monoculture") {
		DIT_tmp[,12] = unique(out1R[[i]]$species)[1]
		DIT_tmp[,13] = unique(out1R[[i]]$species)[1]
		
		
		if(length(grep("daph",t1)) >0 ) { 
			#the carrying capacity from the fitted IGR
			DIT_tmp[,14] = coef(igr_daph[[ti]]$fit_mod)["a1"]
			DIT_tmp[,15] = coef(igr_daph[[ti]]$fit_mod)["a1"]
			DIT_tmp[,16] = out1[[tm1]]$spp_prms$aii[1]
			DIT_tmp[,17] = out1[[tm1]]$spp_prms$aij[1]
			DIT_tmp[,18] = out1[[tm1]]$spp_prms$aii[1]
			DIT_tmp[,19] = out1[[tm1]]$spp_prms$aij[1]
	

		} else{ 
			DIT_tmp[,14] = coef(igr_dia[[ti]]$fit_mod)["a1"]
			DIT_tmp[,15] = coef(igr_dia[[ti]]$fit_mod)["a1"]
			DIT_tmp[,16] = out1[[tm1]]$spp_prms$aii[2]
			DIT_tmp[,17] = out1[[tm1]]$spp_prms$aij[2]
			DIT_tmp[,18] = out1[[tm1]]$spp_prms$aii[2]
			DIT_tmp[,19] = out1[[tm1]]$spp_prms$aij[2]
		}
	} else {
	  if (out1R[[i]]$invade_monoculture[1] == "daph invade"){
		DIT_tmp[,12] = (rspecies)[2]
		DIT_tmp[,13] = (rspecies)[1]
		DIT_tmp[,14] = coef(igr_dia[[ti]]$fit_mod)["a1"]
		DIT_tmp[,15] = coef(igr_daph[[ti]]$fit_mod)["a1"]
		DIT_tmp[,16] = out1[[tm1]]$spp_prms$aii[2]
		DIT_tmp[,17] = out1[[tm1]]$spp_prms$aij[2]
		DIT_tmp[,18] = out1[[tm1]]$spp_prms$aii[1]
		DIT_tmp[,19] = out1[[tm1]]$spp_prms$aij[1]
	  } else {
	    DIT_tmp[,12] = (rspecies)[1] 
	    DIT_tmp[,13] = (rspecies)[2]
	   	DIT_tmp[,14] = coef(igr_daph[[ti]]$fit_mod)["a1"]
	  	DIT_tmp[,15] = coef(igr_dia[[ti]]$fit_mod)["a1"]
		DIT_tmp[,16] = out1[[tm1]]$spp_prms$aii[1]
		DIT_tmp[,17] = out1[[tm1]]$spp_prms$aij[1]
		DIT_tmp[,18] = out1[[tm1]]$spp_prms$aii[2]
		DIT_tmp[,19] = out1[[tm1]]$spp_prms$aij[2]
	  }
	}

	#DIT_tmp = as.data.frame(DIT_tmp) 
	out1R[[i]] = cbind(out1R[[i]],DIT_tmp) #%>% left_join(DIT_tmp)

	}
}

#Take all of the new data and combine it into one data frame: 
m1_DIT = bind_rows(out1R, .id = "column_label")
m1_DIT=m1_DIT %>%
 mutate(alg_per_Nres = (10E6-algae_abundance)/N_res^2)  #"lead" lines up the result
m1_DIT=m1_DIT %>%
 mutate(alg_per_Ninv = (10E6-algae_abundance)/N_inv^2)  #"lead" lines up the result
m1_DIT=m1_DIT %>%
 mutate(alg_per_Ndiff = (10E6-algae_abundance)/lead( N-lag(N),)) 
m1_DIT=m1_DIT %>%
 mutate(alg_perzoo = (10E6-algae_abundance)/(N_inv+N_res)^2)  #"lead" lines up the result

m1_DIT$alg_per_Nres[is.infinite(m1_DIT$alg_per_Nres)] = NA
m1_DIT$alg_per_Ninv[is.infinite(m1_DIT$alg_per_Ninv)] = NA
m1_DIT$alg_per_Ndiff[is.infinite(m1_DIT$alg_per_Ndiff)] = NA
m1_DIT$alg_perzoo[is.infinite(m1_DIT$alg_perzoo)] = NA


#=============================================================================
#Sim data: Loop over all temperature treatments
#=============================================================================
for (w in 1:nwebs){ 

	pop_ts = ceiling(f1S*out_inv1[[w]][,2:nspp])
	nt1 = 1
	nt2 = dim(pop_ts)[1] 
	if(nt2 <=k){ k = 1}

	if(nt1 != nt2) { 

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
	di_webS[w] = list(get_info_dynamics(pop_ts = pop_ts , k=k,with_blocks=TRUE))

	## This code takes the population time-series counts output by the ODEs and 
	## calculates the average Transfer Entropy from each species to every other 
	## species. The goal is to get an overview of the major information pathways 
	## in the web.   
	#=============================================================================
	# This function gives:
	# te_web    Average transfer entropy per species as a pairwise matrix
	#=============================================================================
	te_webS[w] = list( get_te_web( pop_ts = pop_ts, k=k) )

	#=============================================================================
	# This function gives:
	# aiE_web    The AI of the entire ensemble, treated as a single time series. 
	#=============================================================================
	aiE_webS[w] = list( get_ais (  series1 = pop_ts[,3:4], k=k, ensemble = TRUE)    )

	}
}

#=============================================================================
#Post-loop processing to create a single data frame for plotting and stats
#=============================================================================
##Recursively flatten this into a data frame. 
# In order to match the real data, ai1 is always the resident and ai2 is the invader 
# (and follow that pattern for all of the DIT metrics). "3" corresponds to the algae. 
# 
mDIT_tmp = data.frame(matrix( nrow=0, ncol =38 ) ) 
ncnames = c("Algae", "N_res","N_inv", "res_spp","inv_spp", "aiE","te1","te2", "te3","ee1","ee2",
	"ee3","ai1","ai2","ai3","si1","si2","si3","temperature","invade_monoculture","nspp","alg_per_Nres",
	"alg_per_Ninv","zoo","alg_perzoo","res_aii_rt","res_aij_rt","inv_aii_rt",
	"inv_aij_rt", "m_alg_perInv", "m_alg_perzoo", "res_K", "inv_K","res_aii", "res_aij", 
	"inv_aii","inv_aij", "replicate_number" )

colnames(mDIT_tmp) = ncnames
#mesocosms = factor(c("A","B")) #A is Daphnia invader, B is Daphnia resident
nspps2 = factor(c(1,2)) #1 is the pre-invasion phase, 2 is post-invasion phase
nt21 = length(seq(0,tend,delta1))
nt22=nt2-k
for (f in 1:nwebs){
	spp_prms = out1[[f]]$spp_prms
	DIT_tmp = data.frame(matrix(0,nt2,38 ))
	colnames(DIT_tmp) = ncnames

    DIT_tmp[,1] = as.matrix(out_inv1[[f]][,2])
    DIT_tmp[1:(nt2/2),2] = as.matrix(out_inv1[[f]][1:(nt2/2),5]) #First Daph as invader
    DIT_tmp[1:(nt2/2),3] = as.matrix(out_inv1[[f]][1:(nt2/2),4])
    DIT_tmp[(nt2/2+1):nt2,2] = as.matrix(out_inv1[[f]][(nt2/2+1):nt2,4]) #Second Dia as invader
    DIT_tmp[(nt2/2+1):nt2,3] = as.matrix(out_inv1[[f]][(nt2/2+1):nt2,5])
    
    DIT_tmp[,4] = factor( levels = levels(rspecies)) #Invader/resident IDs
    DIT_tmp[,5] = factor( levels = levels(rspecies))
	DIT_tmp[1:(nt2/2),4] = (rspecies)[2]
    DIT_tmp[1:(nt2/2),5] = (rspecies)[1]
   	DIT_tmp[(nt2/2+1):nt2,4] = (rspecies)[1]
    DIT_tmp[(nt2/2+1):nt2,5] = (rspecies)[2]

    DIT_tmp[(k+1):nt2,6] = aiE_webS[[f]]$local

    #Because resident corresponds to te1 (ai1, etc.), and which species is 
    #resident switches:
    DIT_tmp[ (k+1):(nt21*2),7:8] = di_webS[[f]]$te_local[1:(nt21*2-k),4:3]
    DIT_tmp[ (nt21*2+k+1):(nt21*4),7:8] = di_webS[[f]]$te_local[(nt21*2+1):(nt21*4-k),3:4]
    DIT_tmp[(k+1):(nt2), 9 ] = di_webS[[f]]$te_local[,1] #Algae

	DIT_tmp[ (k+1):(nt21*2-k),10:11] = di_webS[[f]]$ee_local[1:(nt21*2-2*k),4:3]
    DIT_tmp[ (nt21*2+1):(nt21*4-2*k),10:11] = di_webS[[f]]$ee_local[(nt21*2+1):(nt21*4-2*k),3:4]
    DIT_tmp[ (k*2):(nt2), 12 ] = di_webS[[f]]$ee_local[,1] #Algae

  	DIT_tmp[ (k+1):(nt21*2),13:14] = di_webS[[f]]$ai_local[1:(nt21*2-k),4:3]
    DIT_tmp[ (nt21*2+k+1):(nt21*4),13:14] = di_webS[[f]]$ai_local[(nt21*2+1):(nt21*4-k),3:4]
    DIT_tmp[(k+1):(nt2), 15 ] = di_webS[[f]]$ai_local[,1] #Algae

  	DIT_tmp[ (k+1):(nt21*2),16:17] = di_webS[[f]]$si_local[1:(nt21*2-k),4:3]
    DIT_tmp[ (nt21*2+k+1):(nt21*4),16:17] = di_webS[[f]]$si_local[(nt21*2+1):(nt21*4-k),3:4]
    DIT_tmp[(k+1):(nt2), 18 ] = di_webS[[f]]$si_local[,1] #Algae

    #Temperatures
    DIT_tmp[,19] = temps[f]

    #Treatment type (invasion vs. monoculture) and number of species
	spp1 = c(1:(nt2/4),(nt2/2+1):(nt2/2+nt2/4)) #Sinlge species
    spp2 = c((nt2/4+1):(nt2/2),(nt2/2+nt2/4+1):nt2) #Invasion, 2 spp
    DIT_tmp[,20] = factor(levels = levels(invader))
    DIT_tmp[1:nt2/2,20] = invader[2]
    DIT_tmp[(nt2/2+1):nt2,20] = invader[3]
    DIT_tmp[spp1,20] =  invader[1]

    DIT_tmp[,21] = factor(levels = levels(nspps2))
    DIT_tmp[spp1,21] = nspps2[1]
    DIT_tmp[spp2,21] = nspps2[2]


	#Algal consumption rates. 
	#Resident
	DIT_tmp[1:(nt2/2),22] = (spp_prms$cC[1,2]*DIT_tmp$Algae[1:(nt2/2)])/DIT_tmp[1:(nt2/2),2]#Dia
	DIT_tmp[(nt2/2+1):nt2,22]  = (spp_prms$cC[1,1]*DIT_tmp$Algae[(nt2/2+1):nt2])/DIT_tmp[(nt2/2+1):nt2,2]#Daph
	DIT_tmp[,22] [is.infinite(DIT_tmp$alg_per_Nres)] = NA
	
	#Invader
    DIT_tmp[1:(nt2/2),23] = (spp_prms$cC[1,1]*DIT_tmp$Algae[1:(nt2/2)])/DIT_tmp[1:(nt2/2),3]#Daph
	DIT_tmp[(nt2/2+1):nt2,23] = (spp_prms$cC[1,2]*DIT_tmp$Algae[(nt2/2+1):nt2])/DIT_tmp[(nt2/2+1):nt2,3]#Dia
	DIT_tmp[,23][is.infinite(DIT_tmp$alg_per_Ninv)] = NA


	#Total Zoo
	DIT_tmp[,24] = DIT_tmp$N_res+DIT_tmp$N_inv
	DIT_tmp[,25] = ((spp_prms$cC[1,1]+spp_prms$cC[1,2])*DIT_tmp$Algae)/DIT_tmp$zoo
	DIT_tmp[,25] [is.infinite(DIT_tmp$alg_perzoo)] = NA


	#Real-time competition for resources: 
	#Resident:
	DIT_tmp[1:(nt2/2),26] = rowSums( out_inv1[[f]][1:(nt2/2),2:3] * matrix( 
		t(spp_prms$cC[,2]^2*spp_prms$rC[,2]) * (1/spp_prms$rR), length(1:(nt2/2)), 2, byrow=T ) )
	DIT_tmp[1:(nt2/2),27] = rowSums( out_inv1[[f]][1:(nt2/2),2:3] * matrix( 
		t(spp_prms$cC[,2]*spp_prms$cC[,1]*spp_prms$rC[,2]) * (1/spp_prms$rR), length(1:(nt2/2)), 2, byrow=T ) )
	
	DIT_tmp[(nt2/2+1):nt2,26] = rowSums( out_inv1[[f]][(nt2/2+1):nt2,2:3] * matrix( 
		t(spp_prms$cC[,1]^2*spp_prms$rC[,1]) * (1/spp_prms$rR), length(1:(nt2/2)), 2, byrow=T ) )
	DIT_tmp[(nt2/2+1):nt2,27] = rowSums( out_inv1[[f]][(nt2/2+1):nt2,2:3] * matrix( 
		t(spp_prms$cC[,2]*spp_prms$cC[,1]*spp_prms$rC[,1]) * (1/spp_prms$rR), length(1:(nt2/2)), 2, byrow=T ) )
	
	#Invader: 
	DIT_tmp[1:(nt2/2),28] = rowSums( out_inv1[[f]][1:(nt2/2),2:3] * matrix( 
		t(spp_prms$cC[,1]^2*spp_prms$rC[,1]) * (1/spp_prms$rR), length(1:(nt2/2)), 2, byrow=T ) )
	DIT_tmp[1:(nt2/2),29] = rowSums( out_inv1[[f]][1:(nt2/2),2:3] * matrix( 
		t(spp_prms$cC[,2]*spp_prms$cC[,1]*spp_prms$rC[,1]) * (1/spp_prms$rR), length(1:(nt2/2)), 2, byrow=T ) )
	
	DIT_tmp[(nt2/2+1):nt2,28] = rowSums( out_inv1[[f]][(nt2/2+1):nt2,2:3] * matrix( 
		t(spp_prms$cC[,2]^2*spp_prms$rC[,2]) * (1/spp_prms$rR), length(1:(nt2/2)), 2, byrow=T ) )
	DIT_tmp[(nt2/2+1):nt2,29] = rowSums( out_inv1[[f]][(nt2/2+1):nt2,2:3] * matrix( 
		t(spp_prms$cC[,2]*spp_prms$cC[,1]*spp_prms$rC[,2]) * (1/spp_prms$rR), length(1:(nt2/2)), 2, byrow=T ) )


	DIT_tmp[1:(nt2/2),30] = mean(DIT_tmp$m_alg_perInv[1:(nt2/2)])
	DIT_tmp[(nt2/2+1):nt2,30] = mean(DIT_tmp$m_alg_perInv[(nt2/2+1):nt2])

	DIT_tmp[,31] = mean(DIT_tmp$alg_perzoo,na.rm=T )
	
	#Get the carrying capacity for each scenario: 
	DIT_tmp[1:(nt2/4),32] = mean(DIT_tmp$N_res[1:(nt2/4)],na.rm=T)
	DIT_tmp[(nt2/4+1):nt2/2,32] = mean(DIT_tmp$N_res[(nt2/4+1):nt2/2],na.rm=T)
	DIT_tmp[(nt2/2+1):(nt2*3/4),32] = mean(DIT_tmp$N_res[(nt2/2+1):(nt2*3/4)],na.rm=T)
	DIT_tmp[(nt2*3/4+1):(nt2),32] = mean(DIT_tmp$N_res[(nt2*3/4+1):(nt2)],na.rm=T)

	DIT_tmp[1:(nt2/4),33] = mean(DIT_tmp$N_inv[1:(nt2/4)],na.rm=T)
	DIT_tmp[(nt2/4+1):nt2/2,33] = mean(DIT_tmp$N_inv[(nt2/4+1):nt2/2],na.rm=T)
	DIT_tmp[(nt2/2+1):(nt2*3/4),33] = mean(DIT_tmp$N_inv[(nt2/2+1):(nt2*3/4)],na.rm=T)
	DIT_tmp[(nt2*3/4+1):(nt2),33] = mean(DIT_tmp$N_inv[(nt2*3/4+1):(nt2)],na.rm=T)

	#Resident competition
    DIT_tmp[1:(nt2/2),34] = spp_prms$aii[2]#Dia
	DIT_tmp[1:(nt2/2),35] = spp_prms$aij[2]#Dia
	DIT_tmp[(nt2/2+1):nt2,34] = spp_prms$aii[1]#Daph
	DIT_tmp[(nt2/2+1):nt2,35] = spp_prms$aij[1]#Daph

	#Invader competition
   	DIT_tmp[1:(nt2/2),36] = spp_prms$aii[1]#Daph
	DIT_tmp[1:(nt2/2),37] = spp_prms$aij[1]#Daph
	DIT_tmp[(nt2/2+1):nt2,36] = spp_prms$aii[2]#Dia
	DIT_tmp[(nt2/2+1):nt2,37] = spp_prms$aij[2]#Dia

	DIT_tmp[,38] = paste(f,DIT_tmp$res_spp,sep="")

    DIT_tmp = as.data.frame(DIT_tmp) 
    mDIT_tmp = rbind(mDIT_tmp,DIT_tmp)
}

#Add two time columns: time1 is absolute time from beginning of simulation, 
#day_n is adjusted relative to the invasion event to match up with the 
#experimental time scale
inv_day = 28
inv_day_n = ((tend*1/delta1+1) - inv_day )
mDIT = cbind(time1 = matrix(seq(0,tend*2+delta1,delta1),dim(mDIT_tmp)[1],1), 
	day_n =matrix(seq(0,(tend*2)*1/delta1+1,1)-inv_day_n,dim(mDIT_tmp)[1],1) , mDIT_tmp)

#=============================================================================
#Saving
#=============================================================================
save(file="daphDia_DIT_13100_nVar_k2f1.var", "m1_DIT","out1R","di_webR",
	"te_webR","si_webR", "aiE_webR" , "mDIT", "out1","out_inv1","di_webS",
	"te_webS","si_webS", "aiE_webS")

save(file="daphDia_DIT_13100_NoVar_k2f1.var", "m1_DIT","out1R","di_webR",
	"te_webR","si_webR", "aiE_webR" , "mDIT", "out1","out_inv1","di_webS",
	"te_webS","si_webS", "aiE_webS")

#=============================================================================
#Plots in combination with the real data using all data
# m1_DIT is the real data
# mDIT is the simulation data
#=============================================================================

m1_DIT_sub = subset(m1_DIT, invade_monoculture != "monoculture") #Daphnia at 28 C
mDIT_sub = subset(mDIT, invade_monoculture != "monoculture")

m1_DIT_sub = m1_DIT  #Daphnia at 28 C
mDIT_sub = mDIT

# m1_DIT_sub = subset(m1_DIT, invade_monoculture != "monoculture" & temperature == 28 & inv_spp == "daphnia") #Daphnia at 28 C
# mDIT_sub = subset(mDIT, nspp == 2 & temperature == 28 & inv_spp == "daphnia")

# m1_DIT_sub = subset(m1_DIT, invade_monoculture != "monoculture" & temperature == 28 & inv_spp == "diaphanosoma") #Daphnia at 28 C
# mDIT_sub = subset(mDIT, nspp == 2 & temperature == 28 & inv_spp == "diaphanosoma")



p1= ggplot()+ 
		geom_line(data=mDIT_sub, mapping= aes(x = time1-tend+1, y =N_res,  color = res_spp, group = interaction(res_spp,replicate_number) ) )+  
		geom_line( data=mDIT_sub, mapping=aes(x = time1-tend+1, y =N_inv,  color = inv_spp, group = interaction(inv_spp,replicate_number) ) ) +  
		geom_point(data=m1_DIT_sub, mapping= aes(x = day_n-inv_day+1, y =N_res,  color = res_spp, group = interaction(res_spp,replicate_number) ) )+  
		geom_point(data=m1_DIT_sub, mapping= aes(x = day_n-inv_day+1, y =N_inv,  color = inv_spp, group = interaction(inv_spp,replicate_number) ) )+  
		facet_grid(temperature~invade_monoculture) +ylim(0,2E2) +xlim(0,90) +
		#scale_color_discrete(name ="", labels = c("Experiment", "Simulation" ) )+
		ylab("Population")+
		xlab("Day")+
		theme(axis.title.x=element_blank(),axis.text.x = element_blank(), axis.ticks = element_blank())


p2 = ggplot()+ 
		geom_line(data=mDIT_sub, mapping= aes(x = time1-tend+1, y =ai1,  color = res_spp, group = interaction(res_spp,replicate_number) ) )+  
		geom_line( data=mDIT_sub, mapping=aes(x = time1-tend+1, y =ai2,  color = inv_spp, group = interaction(inv_spp,replicate_number) ) ) +  
		geom_point(data=m1_DIT_sub, mapping= aes(x = day_n-inv_day+1, y =ai1,  color = res_spp, group = interaction(res_spp,replicate_number) ) )+  
		geom_point(data=m1_DIT_sub, mapping= aes(x = day_n-inv_day+1, y =ai2,  color = inv_spp, group = interaction(inv_spp,replicate_number) ) )+  
		facet_grid(temperature~invade_monoculture) + xlim(0,90)    +
		#scale_color_discrete(name ="", labels = c("Experiment", "Simulation" ) )+
		ylab("Bits (AI) ")+
		xlab("Day")+
		theme(axis.title.x=element_blank(),axis.text.x = element_blank(), axis.ticks = element_blank())

p2 = ggplot()+ 
		geom_line(data=mDIT_sub, mapping= aes(x = time1-tend+1, y =aiE,  color = res_spp, group = interaction(res_spp,replicate_number) ) )+  
		geom_line( data=mDIT_sub, mapping=aes(x = time1-tend+1, y =aiE,  color = inv_spp, group = interaction(inv_spp,replicate_number) ) ) +  
		geom_point(data=m1_DIT_sub, mapping= aes(x = day_n-inv_day+1, y =aiE,  color = res_spp, group = interaction(res_spp,replicate_number) ) )+  
		geom_point(data=m1_DIT_sub, mapping= aes(x = day_n-inv_day+1, y =aiE,  color = inv_spp, group = interaction(inv_spp,replicate_number) ) )+  
		facet_grid(temperature~invade_monoculture) +xlim(0,90)    +
		#scale_color_discrete(name ="", labels = c("Experiment", "Simulation" ) )+
		ylab("Bits (AI) ")+
		xlab("Day")+
		theme(axis.title.x=element_blank(),axis.text.x = element_blank(), axis.ticks = element_blank())


p3 = ggplot()+ 
		geom_line(data=mDIT_sub, mapping= aes(x = time1-tend+1, y =alg_per_Nres,  color = res_spp, group = interaction(res_spp,replicate_number) ) )+  
		geom_line( data=mDIT_sub, mapping=aes(x = time1-tend+1, y =alg_per_Ninv,  color = inv_spp, group = interaction(inv_spp,replicate_number) ) ) +  
		geom_point(data=m1_DIT_sub, mapping= aes(x = day_n-inv_day+1, y =alg_per_Nres,  color = res_spp, group = interaction(res_spp,replicate_number) ) )+  
		geom_point(data=m1_DIT_sub, mapping= aes(x = day_n-inv_day+1, y =alg_per_Ninv,  color = inv_spp, group = interaction(inv_spp,replicate_number) ) )+  
		facet_grid(temperature~invade_monoculture) +ylim(0,2E6) +xlim(0,90) +
		#scale_color_discrete(name ="", labels = c("Experiment", "Simulation" ) )+
		xlab("Day")+
		ylab("Algae consumed per N ")

gp1 = ggplotGrob(p1)
gp2 = ggplotGrob(p2)
gp3 = ggplotGrob(p3)
grid.draw(rbind(gp1, gp2,gp3))


#To save as a PDF
pdf("./exp_v_simAll1.pdf",width = 8, height = 10)
grid.draw(rbind(gp1, gp2,gp3))
dev.off()

#=============================================================================
m1_DIT_sub = subset(m1_DIT, invade_monoculture != "monoculture") #Daphnia at 28 C
mDIT_sub = subset(mDIT, invade_monoculture != "monoculture")

m1_DIT_sub$res_spp = revalue(m1_DIT_sub$res_spp, c("daphnia" = "daphnia_E", "diaphanosoma" = "dia_E" ) )
m1_DIT_sub$inv_spp = revalue(m1_DIT_sub$inv_spp, c("daphnia" = "daphnia_E", "diaphanosoma" = "dia_E" ) )


#Consumption vs. complexity 1
ggplot()+ 
	geom_point( data=mDIT_sub, shape = 0, mapping= aes(x = ai1, y =alg_per_Nres ,  color = res_spp , group = interaction(res_spp,replicate_number) ) ) +
	geom_point( data=mDIT_sub, shape = 0, mapping=aes(x = ai2, y =alg_per_Ninv,  color = inv_spp, group = interaction(inv_spp,replicate_number) ) ) +  
	geom_point( data=m1_DIT_sub,shape = 18, mapping= aes(x = ai1, y =alg_per_Nres,  color = res_spp, group = interaction(res_spp,replicate_number) ) ) +
	geom_point( data=m1_DIT_sub, shape = 18, mapping=aes(x = ai2, y =alg_per_Ninv,  color = inv_spp, group = interaction(inv_spp,replicate_number) ) ) +  
	facet_grid(temperature~invade_monoculture)  +ylim(0,1E6)
	#+ylim(0,0.001)

#Consumption vs. complexity 2
ggplot()+ 
	geom_point( data=mDIT_sub, shape = 0, mapping= aes(x = aiE, y =alg_per_Nres,  color = res_spp, group = interaction(res_spp,replicate_number) ) ) +
	geom_point( data=mDIT_sub, shape = 0, mapping=aes(x = aiE, y =alg_per_Ninv,  color = inv_spp, group = interaction(inv_spp,replicate_number) ) ) +  
	geom_point( data=m1_DIT_sub,shape = 18, mapping= aes(x = aiE, y =alg_per_Nres,  color = res_spp, group = interaction(res_spp,replicate_number) ) ) +
	geom_point( data=m1_DIT_sub, shape = 18, mapping=aes(x = aiE, y =alg_per_Ninv,  color = inv_spp, group = interaction(inv_spp,replicate_number) ) ) +  
	facet_grid(temperature~invade_monoculture)  +ylim(0,5E5)


#Distance from Equilbrium vs. complexity 1
ggplot()+ 
	geom_point( data=mDIT_sub, shape = 0, mapping= aes(y = (N_res-res_K), x =ai1,  color = res_spp, group = interaction(res_spp,replicate_number) ) ) +
	geom_point( data=mDIT_sub, shape = 0, mapping=aes(y = (N_inv-inv_K), x =ai2,  color = inv_spp, group = interaction(inv_spp,replicate_number) ) ) +  
	geom_point( data=m1_DIT_sub,shape = 18, mapping= aes(y =  (N_res-res_K), x =ai1,  color = res_spp, group = interaction(res_spp,replicate_number) ) ) +
	geom_point( data=m1_DIT_sub, shape = 18, mapping=aes(y = (N_inv-inv_K), x =ai2,  color = inv_spp, group = interaction(inv_spp,replicate_number) ) ) +  
	facet_grid(temperature~invade_monoculture)  +ylim(-100,100) +xlim(-5,15)

#Distance from Equilbrium vs. complexity 2
ggplot()+ 
	geom_point( data=mDIT_sub, shape = 0, mapping= aes(y = (N_res-res_K), x =aiE,  color = res_spp, group = interaction(res_spp,replicate_number) ) ) +
	geom_point( data=mDIT_sub, shape = 0, mapping=aes(y = (N_inv-inv_K), x =aiE,  color = inv_spp, group = interaction(inv_spp,replicate_number) ) ) +  
	geom_point( data=m1_DIT_sub,shape = 18, mapping= aes(y =  (N_res-res_K), x =aiE,  color = res_spp, group = interaction(res_spp,replicate_number) ) ) +
	geom_point( data=m1_DIT_sub, shape = 18, mapping=aes(y = (N_inv-inv_K), x =aiE,  color = inv_spp, group = interaction(inv_spp,replicate_number) ) ) +  
	facet_grid(temperature~invade_monoculture)  +ylim(-300,300) +xlim(-5,15)

#Consumption vs. distance from Equilibrium -- this gives such a nice smooth fit!
ggplot()+ 
	geom_point( data=mDIT_sub, shape = 0, mapping= aes(x = (N_res-res_K), y =alg_per_Nres,  color = res_spp, group = interaction(res_spp,replicate_number) ) ) +
	geom_point( data=mDIT_sub, shape = 0, mapping=aes(x = (N_inv-inv_K), y =alg_per_Ninv,  color = inv_spp, group = interaction(inv_spp,replicate_number) ) ) +  
	geom_point( data=m1_DIT_sub,shape = 18, mapping= aes(x =  (N_res-res_K), y =alg_per_Nres,  color = res_spp, group = interaction(res_spp,replicate_number) ) ) +
	geom_point( data=m1_DIT_sub, shape = 18, mapping=aes(x = (N_inv-inv_K), y =alg_per_Ninv,  color = inv_spp, group = interaction(inv_spp,replicate_number) ) ) +  
	facet_grid(temperature~invade_monoculture)   +ylim(0,8E4)

#Competition vs. complexity 2
ggplot()+ 
	geom_point( data=mDIT_sub,  mapping= aes(x = (inv_aij_rt*res_aij_rt)/(res_aii_rt*inv_aii_rt), y =aiE,  color = res_spp, group = interaction(res_spp,replicate_number) ) ) +
	geom_point( data=mDIT_sub,  mapping=aes(x = (inv_aij_rt*res_aij_rt)/(res_aii_rt*inv_aii_rt), y =aiE,  color = inv_spp, group = interaction(inv_spp,replicate_number) ) ) +  
	# geom_point( data=m1_DIT_sub,shape = 18, mapping= aes(y =  (N_res-res_K), x =aiE,  color = res_spp, group = interaction(res_spp,replicate_number) ) ) +
	# geom_point( data=m1_DIT_sub, shape = 18, mapping=aes(y = (N_inv-inv_K), x =aiE,  color = inv_spp, group = interaction(inv_spp,replicate_number) ) ) +  
	facet_grid(temperature~invade_monoculture)



#=============================================================================
#=============================================================================
m1_DIT_sub = subset(m1_DIT, invade_monoculture != "monoculture") #Daphnia at 28 C
mDIT_sub = subset(mDIT, invade_monoculture != "monoculture")

m1_DIT_sub$res_spp = revalue(m1_DIT_sub$res_spp, c("daphnia" = "daphnia_E", "diaphanosoma" = "dia_E" ) )
m1_DIT_sub$inv_spp = revalue(m1_DIT_sub$inv_spp, c("daphnia" = "daphnia_E", "diaphanosoma" = "dia_E" ) )

#Some summary numbers
#For deterministic sim: 
aie_sum=mDIT_sub %>% 
dplyr::group_by(temperature,invade_monoculture,replicate_number, res_spp,inv_spp,nspp,res_aii,res_aij,inv_aii,inv_aij) %>%
dplyr::summarize(aiE_mean = aiE[time1==tend*2], #mean(aiE,na.rm=T),
	ai1_mean = ai1[time1==tend*2],# mean(ai1,na.rm=T), 
	ai2_mean = ai2[time1==tend*2], #mean(ai2,na.rm=T),
	alg_Nres_mean = mean(alg_per_Nres,na.rm=T),
	alg_Ninv_mean = mean(alg_per_Ninv,na.rm=T)) %>% as.data.frame

#For stochastic sim: 
mDIT_sub2 = subset(mDIT_sub, mDIT_sub$time1 >= 3/4*(tend*2) & mDIT_sub$time1 <= (tend*2)) 
#mDIT_sub2 = subset(mDIT_sub, mDIT_sub$time1 >= tend+tend/4 & mDIT_sub$time1 <= (tend+tend/2)) 

aie_sum=mDIT_sub2 %>% 
dplyr::group_by(temperature,invade_monoculture,replicate_number, res_spp,inv_spp,nspp,res_aii,res_aij,inv_aii,inv_aij) %>%
dplyr::summarize(aiE_mean = mean(aiE,na.rm=T),
	ai1_mean =mean(ai1,na.rm=T), 
	ai2_mean = mean(ai2,na.rm=T),
	alg_Nres_mean = mean(alg_per_Nres,na.rm=T),
	alg_Ninv_mean = mean(alg_per_Ninv,na.rm=T)) %>% as.data.frame

#For real data: 
aie_sumR=m1_DIT_sub %>% 
dplyr::group_by(temperature,invade_monoculture,replicate_number, res_spp,inv_spp,res_aii,res_aij,inv_aii,inv_aij) %>%
dplyr::summarize(aiE_mean = mean(aiE,na.rm=T),
	ai1_mean = mean(ai1,na.rm=T), 
	ai2_mean = mean(ai2,na.rm=T),
	alg_Nres_mean = mean(alg_per_Nres,na.rm=T),
	alg_Ninv_mean = mean(alg_per_Ninv,na.rm=T)) %>% as.data.frame


#Competition vs. complexity 2
ggplot()+ 
	geom_point( data=mDIT_sub,  mapping= aes(x = (inv_aij_rt*res_aij_rt)/(res_aii_rt*inv_aii_rt), y =aiE,  color = res_spp, group = interaction(res_spp,replicate_number) ) ) +
	geom_point( data=mDIT_sub,  mapping=aes(x = (inv_aij_rt*res_aij_rt)/(res_aii_rt*inv_aii_rt), y =aiE,  color = inv_spp, group = interaction(inv_spp,replicate_number) ) ) +  
	geom_point(data=aie_sum,  mapping= aes(x =(inv_aij*res_aij)/(inv_aii*res_aii) , y =aiE_mean,  color = "1", group=interaction(res_spp,replicate_number) ),size=3) +
	geom_point(data=aie_sum,  mapping= aes(x =(inv_aij*res_aij)/(inv_aii*res_aii) , y =aiE_mean,  color = "2", group=interaction(inv_spp,replicate_number) ),size=3 ) +
	geom_point( data=m1_DIT_sub,  mapping= aes(x = (inv_aij*res_aij)/(res_aii*inv_aii), y =aiE,  color = res_spp, group = interaction(res_spp,replicate_number) ) ) +
	geom_point( data=m1_DIT_sub,  mapping=aes(x = (inv_aij*res_aij)/(res_aii*inv_aii), y =aiE,  color = inv_spp, group = interaction(inv_spp,replicate_number) ) ) +  
	geom_point(data=aie_sumR,  mapping= aes(x =(inv_aij*res_aij)/(inv_aii*res_aii) , y =aiE_mean,  color = "3", group=interaction(res_spp,replicate_number) ),size=3) +
	geom_point(data=aie_sumR,  mapping= aes(x =(inv_aij*res_aij)/(inv_aii*res_aii) , y =aiE_mean,  color = "4", group=interaction(inv_spp,replicate_number) ),size=3 ) +
	facet_grid(temperature~invade_monoculture)

#Complexity vs. temp.
ggplot()+ 
	geom_point( data=aie_sum,  mapping= aes(y =aiE_mean, x =temperature,  color = interaction( res_spp, nspp ) ) )  +
	geom_point( data=aie_sum,  mapping=aes(y =aiE_mean, x =temperature,  color = interaction( inv_spp, nspp) ) ) +  
	geom_point( data=aie_sumR,  mapping= aes(y =aiE_mean, x =temperature,  color = interaction( res_spp, nspp ) ) )  +
	geom_point( data=aie_sumR,  mapping=aes(y =aiE_mean, x =temperature,  color = interaction( inv_spp, nspp) ) ) +  
	facet_grid(invade_monoculture~.) 

ggplot()+ 
	geom_point( data=aie_sum,  mapping= aes(y =ai1_mean, x =temperature,  color = interaction( res_spp, nspp ) ) )  +
	geom_point( data=aie_sum,  mapping=aes(y =ai2_mean, x =temperature,  color = interaction( inv_spp, nspp) ) ) +  
	geom_point( data=aie_sumR,  mapping= aes(y =ai1_mean, x =temperature,  color = interaction( res_spp, nspp ) ) )  +
	geom_point( data=aie_sumR,  mapping=aes(y =ai2_mean, x =temperature,  color = interaction( inv_spp, nspp) ) ) +  
	facet_grid(invade_monoculture~.) 

#complexity vs. competition
ggplot()+ 
	geom_point( data=aie_sum,  mapping= aes(x =(inv_aij*res_aij)/(inv_aii*res_aii) , y =aiE_mean,  color = interaction( res_spp, nspp ) ) )  +
	geom_point( data=aie_sum,  mapping=aes(x = (inv_aij*res_aij)/(inv_aii*res_aii), y =aiE_mean,  color = interaction( inv_spp, nspp) ) ) +  
	geom_point( data=aie_sumR,  mapping= aes(x =(inv_aij*res_aij)/(inv_aii*res_aii) , y =aiE_mean,  color = interaction( res_spp, nspp ) ) )  +
	geom_point( data=aie_sumR,  mapping=aes(x = (inv_aij*res_aij)/(inv_aii*res_aii), y =aiE_mean,  color = interaction( inv_spp, nspp) ) ) +  
	facet_grid(invade_monoculture~.) 

#Competition vs. temp
ggplot()+ 
	geom_point( data=aie_sum,  mapping= aes(y =inv_aij*res_aij/(inv_aii*res_aii) , x =temperature,  color = interaction( res_spp, nspp ) ) )  +
	geom_point( data=aie_sum,  mapping=aes(y = inv_aij*res_aij/(inv_aii*res_aii), x =temperature,  color = interaction( inv_spp, nspp) ) ) +  
	geom_point( data=aie_sumR,  mapping= aes(y =inv_aij*res_aij/(inv_aii*res_aii) , x =temperature,  color = interaction( res_spp, nspp ) ) )  +
	geom_point( data=aie_sumR,  mapping=aes(y = inv_aij*res_aij/(inv_aii*res_aii), x =temperature,  color = interaction( inv_spp, nspp) ) ) +  
	facet_grid(invade_monoculture~.)



ggplot()+ 
	geom_point( data=aie_sum,  mapping= aes(x =inv_aij*res_aij/(inv_aii*res_aii) , y=alg_Nres_mean,  color = interaction( res_spp, nspp ) ) )  +
	geom_point( data=aie_sum,  mapping=aes(x = inv_aij*res_aij/(inv_aii*res_aii), y =alg_Ninv_mean,  color = interaction( inv_spp, nspp) ) ) +  
	#geom_point( data=m1_DIT_sub,shape = 18, mapping= aes(x = aiE, y =alg_per_Nres,  color = res_spp, group = interaction(res_spp,replicate_number) ) ) +
	#geom_point( data=m1_DIT_sub, shape = 18, mapping=aes(x = aiE, y =alg_per_Ninv,  color = inv_spp, group = interaction(inv_spp,replicate_number) ) ) +  
	facet_grid(.~invade_monoculture) + ylim(0,8E4)




#=============================================================================
m1_DIT_sub = subset(m1_DIT, inv_spp == "daphnia" & temperature == 28 & invade_monoculture == "daph invade") #Daphnia at 28 C
mDIT_sub = subset(mDIT, inv_spp == "daphnia" & temperature == 28 & invade_monoculture == "daph invade") #Daphnia invasion 

ggplot()+ 
	geom_point( data=mDIT_sub, shape = 0, mapping= aes(x = ai2, y =alg_per_Nres,  color = time1))
	geom_point( data=mDIT_sub, shape = 0, mapping=aes(x = ai2, y =alg_per_Ninv,  color = time1)) +  , group = interaction(inv_spp,replicate_number) ) )

m1_DIT_sub = subset(m1_DIT, inv_spp == "daphnia"& invade_monoculture != "monoculture"& temperature == temps[[3]]) 
mDIT_sub = subset(mDIT, inv_spp == "daphnia"& nspp == 2 & temperature == temps[[3]])

mDIT_sub = subset(mDIT, nspp == 2 & temperature == temps[[2]])


aa=mDIT_sub[mDIT_sub$time1<40,]
aa2=subset(aa, res_spp == "diaphanosoma")
plot(subset(aa, res_spp == "diaphanosoma")$ai1)



m1_DIT_sub = subset(m1_DIT, invade_monoculture == "monoculture"& temperature == 28) 
mDIT_sub = subset(mDIT, nspp == 1 & temperature == 28)

ggplot()+ 
		geom_line(data=mDIT_sub, mapping= aes(x = time1, y =N_res,  color = res_spp, group = interaction(res_spp,replicate_number) ) )+  
		geom_point(data=m1_DIT_sub, mapping= aes(x = day_n, y =N_res,  color = res_spp, group = interaction(res_spp,replicate_number) ) )+  
		scale_color_discrete(name ="", labels = c("Experiment", "Simulation" ) )+
		ylab("Population")+
		xlab("Day")+
		theme(axis.title.x=element_blank(),axis.text.x = element_blank(), axis.ticks = element_blank())


