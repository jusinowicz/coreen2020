
#=============================================================================
#File 2.A: Leslie-Gower population model, NLS
#=============================================================================
# This R script differs from the previous in that we define our own model 
# and then use NLS to fit it. Specifically, we fit the Leslie-Gower model. 
# The model is fit to residents and invader populations seperately. 
# The intrinsic rate of reproduction, %lambda, is fit in two different ways: 
# 	1. by allowing nls to fit it as an intercept in the model. 
# 
# Fitting with NLS is just meant as a first pass to explore the data. I 
# probably would not trust any of the parameter values that it spits out 
# at this point! I.e., don't believe lambda, alpha_ii or alpha_ij for a lot
# of the model fits. Trust the NLME fits in the complimentary model. 
#
# As this code is currently written, NLS can barely handle fitting most of 
# the replicates. It often seems to fail when the data should be "easiest"
# to interpret, and maybe this makes sense. For example, it is bad at fitting
# the model parameters when Daphnia is the resident and quickly reaches its 
# equilibrium abundance: most of the variation is likely due to (AR) noise. 
#
# First, loop over both diaphanosoma and daphnia when diaphanasoma is the 
# invader
# Second, loop over both diaphanosoma and daphnia when daphnia is the invader
#=============================================================================

#=============================================================================
#Load libraries
#=============================================================================
library(tidyverse)
library(nlme)
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
#Define the growth rate functions for the LG population model(s)
#=============================================================================
#The model for the resident with only the resident species
LG_res = formula (Ndiff_res ~ lambda/(1+alpha_ii*N_res) + 0.5)
#The model for the invader with only the invader species
LG_inv = formula (Ndiff_inv ~ lambda/(1+alpha_ij*N_res) + 0.5)
#The model with both species (this tends to have convergence problems)
LG_2sp = formula (Ndiff_inv ~ lambda/(1+alpha_ii*N_inv+alpha_ij*N_res)+0.5)

#=============================================================================
# First, loop over both diaphanosoma and daphnia when diaphanasoma is the 
# invader
#=============================================================================
all_models_dia_invade = vector( "list", ntreatments*2) #To store all of the fitted models
all_models_daph_res = vector( "list", ntreatments*2) #To store all of the fitted models

model_fit_dia_invade = vector( "list", ntreatments*2) #The predicted IGR based on the fitted LME 
model_fit_daph_res = vector( "list", ntreatments*2)

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
        		arrange(replicate.number) %>%
       			mutate(Ndiff_res = lead( N/lag(N),)) #"lead" lines up the result 
		res_tmp$Ndiff_res[res_tmp$day_n == max(res_tmp$day_n )] =NA #Remove last day
		res_tmp = res_tmp %>% mutate(tdiff =day_n-lag(day_n)) #Size of time step
		res_tmp = res_tmp %>% mutate(N_res = N)
        res_tmp$tdiff[res_tmp$tdiff<0] = NA #Remove negative time steps
		#Make this a grouped data object:
		res_data = groupedData(Ndiff_res~N_res+day_n|replicate.number, data = res_tmp)
	
		#Make a new invader data set to fit the growth rate function with nlme/nls 
		#The new data set includes a column for n(t+1)/n(t) and a columnf for the time 
		#interval across subsequent measurements.  
		inv_tmp = subset(m1_data_long, mesocosm.ID == u_treats &  species == u_invader   )
        inv_tmp = na.exclude(inv_tmp) #NLME won't work with NAs 
       	#Arrange the data by replicate number, then add a new column for the delta N
        inv_tmp = inv_tmp %>% 
        		arrange(replicate.number) %>%
       			mutate(Ndiff_inv = lead( N/lag(N),)) #"lead" lines up the result 
		inv_tmp$Ndiff_inv[inv_tmp$day_n == max(inv_tmp$day_n )] =NA #Remove last day
		inv_tmp = inv_tmp %>% mutate(tdiff =day_n-lag(day_n)) #Size of time step
		inv_tmp = inv_tmp %>% mutate(N_inv = N)
		#This line adds the resident densities on the matching days from the matching 
		#replicates to the data table inv_tmp
		inv_tmp = inv_tmp %>% left_join( select( res_tmp, replicate.number,day_n,N_res),
		 by= (c( "day_n" = "day_n", "replicate.number"="replicate.number" )) )
        inv_tmp$tdiff[inv_tmp$tdiff<0] = NA #Remove negative time steps
		#Make this a grouped data object:
		inv_data = groupedData(Ndiff_inv~N_res+N_inv+day_n|replicate.number, data = inv_tmp)
	
		#================================================================================
		# Use NLS to test basic model assumptions and convergence, fine-tune fits, etc.
		# This is essentially what is nested within NLME, it just does not allow for 
		# any fixed or random effects to be properly accounted for, in either fitting or
		# errors. 
		#================================================================================
		rep_use = 6 #Change this in order to pull data from a specific replicate
		###Resident
		#Pull out one replicate AND only use the N that are pre-invasion: 
		#Note, if you use data that are post invasion to fit the single-species resident 
		#model, LG_1sp, then alpha_ii will take on weird values. 
		#Choose 1: 
		mydata1 = subset(res_data,replicate.number==rep_use & day_n < (inv_day+12) )
		#mydata1 = subset(res_data, day_n < (inv_day+12) )

		tryCatch( { m1=nls(LG_res, data = mydata1, start=list(lambda = 1.1, alpha_ii=0.5) )

		#Model parameters: 
		lam1 =  m1$m$getPars()[1]
		aii =  m1$m$getPars()[2]
		p1 = ( lam1/(1+aii*mydata1$N_res)+0.5 )
		#p1=predict(m1, data = mydata1) #Same thing
		
		#To shift things back to the population scale ( we're fitting the growth
		#rate function and not population counts directly)
		#Note, this should be exactly the same as N_res in res_data: 
		actN_res = mydata1$Ndiff_res*(mydata1$N_res)
		#To get the predicted N_res, multiply the NLS fit of our G by N: 
		predN_res = p1*(mydata1$N_res )
		predN_res2 = p1*(subset(res_data,replicate.number==rep_use)$N_res)

		#Basic plots of data and fits: 
		par(mfrow = c(2,1))
		#First plot
		plot(actN_res, ylim = c(0,max( c( actN_res,predN_res ),na.rm=T ) ),xlab = "Time",
			ylab = "Resident density" )
		points(predN_res,col="blue")
		#Second, over full time series
		plot(subset(res_data,replicate.number==rep_use)$N_res , ylim = c(0,max( subset(res_data,replicate.number==rep_use)$N_res,
			na.rm=T ) ), xlab = "Time",	ylab = "Resident density" )
		points(predN_res2, col ="blue")

		#Save model fits and predicted values
		all_models_daph_res[[index1]] = m1
		model_fit1 = data.frame(day_n =mydata1$day_n, N_res = predN_res )
	    model_fit_daph_res[[index1]] = model_fit1
		}, error = function(e) {} ) 

		###Invader
		#Pull out one replicate AND only use the N that are post-invasion but pre invader
		#equilibrium (this is very subjective, but also tends to be fairly robust) : 
		#Choose 1: 
		#mydata1 = subset(inv_data,replicate.number==6 & day_n > inv_day & day_n < (day_n+15) )
		#mydata1 = subset(inv_data,replicate.number==rep_use & day_n >= inv_day  )
		mydata1 = subset(inv_data, day_n >= inv_day  )

		tryCatch( { m1=nls(LG_inv, data = mydata1, start=list(lambda = 1.1, alpha_ij=0.1) )
		#Model parameters: 
		lam1 =  m1$m$getPars()[1]
		aij =  m1$m$getPars()[2]
		p1 = ( lam1/(1+aij*mydata1$N_res)+0.5 )
		#p1=predict(m1, data = mydata1) #Same thing
		
		#To shift things back to the population scale ( we're fitting the growth
		#rate function and not population counts directly)
		#Note, this should be exactly the same as N_res in res_data: 
		actN_inv = mydata1$Ndiff_inv*(mydata1$N_inv)
		#To get the predicted N_res, multiply the NLS fit of our G by N: 
		predN_inv = p1*(mydata1$N_inv )
		predN_inv2 = p1*(subset(inv_data,replicate.number==rep_use)$N_inv)

		#Basic plots of data and fits: 
		par(mfrow = c(2,1))
		#First plot
		plot(actN_inv, ylim = c(0,max( c( actN_inv,predN_inv ),na.rm=T ) ),xlab = "Time",
			ylab = "Invader density" )
		points(predN_inv,col="blue")
		#Second, over full time series
		plot(subset(inv_data,replicate.number==rep_use)$N_inv, ylim = c(0,max( subset(inv_data,replicate.number==rep_use)$N_inv ,
			na.rm=T ) ), xlab = "Time",	ylab = "Invader density" )
		points(predN_inv2, col ="blue")

		#Save model fits and predicted values
		all_models_dia_invade[[index1]] = m1
		model_fit1 = data.frame(day_n =mydata1$day_n, N_inv = predN_inv )
	    model_fit_dia_invade[[index1]] = model_fit1
		}, error = function(e) {} ) 


}


#=============================================================================
# Second, loop over both diaphanosoma and daphnia when daphnia is the invader
#=============================================================================
all_models_daph_invade = vector( "list", ntreatments*2) #To store all of the fitted models
all_models_dia_res = vector( "list", ntreatments*2) #To store all of the fitted models

model_fit_daph_invade = vector( "list", ntreatments*2) #The predicted IGR based on the fitted LME 
model_fit_dia_res = vector( "list", ntreatments*2)

#The resident and invader: 
i=2  
	#Loop over treatments 
	u_invader = invader[-i]
	u_res = invader[i]

	for (n in 1: (ntreatments)){
		#Since the treatments are regular, set the start/end position of
		#the subset qith:  
		pos1 = (n-1)*2 + n+invasions_per
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
        		arrange(replicate.number) %>%
       			mutate(Ndiff_res = lead( N/lag(N),)) #"lead" lines up the result 
		res_tmp$Ndiff_res[res_tmp$day_n == max(res_tmp$day_n )] =NA #Remove last day
		res_tmp = res_tmp %>% mutate(tdiff =day_n-lag(day_n)) #Size of time step
		res_tmp = res_tmp %>% mutate(N_res = N)
        res_tmp$tdiff[res_tmp$tdiff<0] = NA #Remove negative time steps
		#Make this a grouped data object:
		res_data = groupedData(Ndiff_res~N_res+day_n|replicate.number, data = res_tmp)
	
		#Make a new invader data set to fit the growth rate function with nlme/nls 
		#The new data set includes a column for n(t+1)/n(t) and a columnf for the time 
		#interval across subsequent measurements.  
		inv_tmp = subset(m1_data_long, mesocosm.ID == u_treats &  species == u_invader   )
        inv_tmp = na.exclude(inv_tmp) #NLME won't work with NAs 
       	#Arrange the data by replicate number, then add a new column for the delta N
        inv_tmp = inv_tmp %>% 
        		arrange(replicate.number) %>%
       			mutate(Ndiff_inv = lead( N/lag(N),)) #"lead" lines up the result 
		inv_tmp$Ndiff_inv[inv_tmp$day_n == max(inv_tmp$day_n )] =NA #Remove last day
		inv_tmp = inv_tmp %>% mutate(tdiff =day_n-lag(day_n)) #Size of time step
		inv_tmp = inv_tmp %>% mutate(N_inv = N)
		#This line adds the resident densities on the matching days from the matching 
		#replicates to the data table inv_tmp
		inv_tmp = inv_tmp %>% left_join( select( res_tmp, replicate.number,day_n,N_res),
		 by= (c( "day_n" = "day_n", "replicate.number"="replicate.number" )) )
        inv_tmp$tdiff[inv_tmp$tdiff<0] = NA #Remove negative time steps
		#Make this a grouped data object:
		inv_data = groupedData(Ndiff_inv~N_res+N_inv+day_n|replicate.number, data = inv_tmp)
	
		#================================================================================
		# Use NLS to test basic model assumptions and convergence, fine-tune fits, etc.
		# This is essentially what is nested within NLME, it just does not allow for 
		# any fixed or random effects to be properly accounted for, in either fitting or
		# errors. 
		#================================================================================
		rep_use = 6 #Change this in order to pull data from a specific replicate
		###Resident
		#Pull out one replicate AND only use the N that are pre-invasion: 
		#Note, if you use data that are post invasion to fit the single-species resident 
		#model, LG_1sp, then alpha_ii will take on weird values. 
		#Choose 1: 
		#mydata1 = subset(res_data,replicate.number==rep_use & day_n < (inv_day+12) )
		mydata1 = subset(res_data, day_n < (inv_day+12) )

		tryCatch( { m1=nls(LG_res, data = mydata1, start=list(lambda = 1.1, alpha_ii=0.5) )
		#Model parameters: 
		lam1 =  m1$m$getPars()[1]
		aii =  m1$m$getPars()[2]
		p1 = ( lam1/(1+aii*mydata1$N_res)+0.5 )
		#p1=predict(m1, data = mydata1) #Same thing
		
		#To shift things back to the population scale ( we're fitting the growth
		#rate function and not population counts directly)
		#Note, this should be exactly the same as N_res in res_data: 
		actN_res = mydata1$Ndiff_res*(mydata1$N_res)
		#To get the predicted N_res, multiply the NLS fit of our G by N: 
		predN_res = p1*(mydata1$N_res )
		predN_res2 = p1*(subset(res_data,replicate.number==rep_use)$N_res)

		#Basic plots of data and fits: 
		par(mfrow = c(2,1))
		#First plot
		plot(actN_res, ylim = c(0,max( c( actN_res,predN_res ),na.rm=T ) ),xlab = "Time",
			ylab = "Resident density" )
		points(predN_res,col="blue")
		#Second, over full time series
		plot(subset(res_data,replicate.number==rep_use)$N_res , ylim = c(0,max( subset(res_data,replicate.number==rep_use)$N_res,
			na.rm=T ) ), xlab = "Time",	ylab = "Resident density" )
		points(predN_res2, col ="blue")

		#Save model fits and predicted values
		all_models_dia_res[[index1]] = m1
		model_fit1 = data.frame(day_n =mydata1$day_n, N_inv = predN_res )
	    model_fit_dia_res[[index1]] = model_fit1
		}, error = function(e) {} ) 


		###Invader
		#Pull out one replicate AND only use the N that are post-invasion but pre invader
		#equilibrium (this is very subjective, but also tends to be fairly robust) : 
		#Choose 1: 
		#mydata1 = subset(inv_data,replicate.number==6 & day_n > inv_day & day_n < (day_n+15) )
		#mydata1 = subset(inv_data,replicate.number==rep_use & day_n >= inv_day  )
		mydata1 = subset(inv_data, day_n >= inv_day  )

		tryCatch( { m1=nls(LG_inv, data = mydata1, start=list(lambda = 1.1, alpha_ij=0.1) )
	
		#Model parameters: 
		lam1 =  m1$m$getPars()[1]
		aij =  m1$m$getPars()[2]
		p1 = ( lam1/(1+aij*mydata1$N_res)+0.5 )
		#p1=predict(m1, data = mydata1) #Same thing
		
		#To shift things back to the population scale ( we're fitting the growth
		#rate function and not population counts directly)
		#Note, this should be exactly the same as N_res in res_data: 
		actN_inv = mydata1$Ndiff_inv*(mydata1$N_inv)
		#To get the predicted N_res, multiply the NLS fit of our G by N: 
		predN_inv = p1*(mydata1$N_inv )
		predN_inv2 = p1*(subset(inv_data,replicate.number==rep_use)$N_inv)

		#Basic plots of data and fits: 
		par(mfrow = c(2,1))
		#First plot
		plot(actN_inv, ylim = c(0,max( c( actN_inv,predN_inv ),na.rm=T ) ),xlab = "Time",
			ylab = "Invader density" )
		points(predN_inv,col="blue")
		#Second, over full time series
		plot(subset(inv_data,replicate.number==rep_use)$N_res, ylim = c(0,max( subset(inv_data,replicate.number==rep_use)$N_inv ,
			na.rm=T ) ), xlab = "Time",	ylab = "Invader density" )
		points(predN_inv2, col ="blue")

		#Save model fits and predicted values
		all_models_daph_invade[[index1]] = m1
		model_fit1 = data.frame(day_n =mydata1$day_n, N_inv = predN_inv )
	    model_fit_daph_invade[[index1]] = model_fit1
		}, error = function(e) {} ) 


}
