
#=============================================================================
#File 2: Leslie-Gower population model 
#=============================================================================
# This R script differs from the previous in that we define our own model 
# and then use nlme to fit it. Specifically, we fit the Leslie-Gower model. 
# The model is fit to residents and invader populations seperately. 
# The intrinsic rate of reproduction, %lambda, is fit in two different ways: 
# 	1. by assuming that it has been measured in the earliest resident times. 
#	2. by allowing nlme to fit it as an intercept in the model. 
# 
# As previously: Use nlme to fit an AR model to each treatment time series. 
# Temporal autocorrelation is treated using corARMA correlation structure. 
# Mesocosm is a random effect. 
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
all_models_dia_invade = vector( "list", ntreatments) #To store all of the fitted models
all_models_daph_res = vector( "list", ntreatments) #To store all of the fitted models

model_fit_dia_invade = vector( "list", ntreatments) #The predicted IGR based on the fitted LME 
model_fit_daph_res = vector( "list", ntreatments)

#The resident and invader: 
i = 1
	#Loop over treatments 
	u_invader = invader[-i]
	u_res = invader[i]

for (n in 1: (ntreatments)){
		#Since the treatments are regular, set the start/end position of
		#the subset qith:  
		pos1 = (n-1)*2 + n
		pos2 = pos1+2
		u_treats = treats[pos1:pos2]
		index1 = n#+(i-1)*ntreatments

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
        res_tmp = na.exclude(res_tmp) #NLME won't work with NAs 
		
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
        inv_tmp = na.exclude(inv_tmp) #NLME won't work with NAs 
		#Make this a grouped data object:
		inv_data = groupedData(Ndiff_inv~N_res+N_inv+day_n|replicate.number, data = inv_tmp)
	

		#================================================================================
		#Use NLME to fit the model. This should work better than NLS since we can include 
		#AR correlated error and allow for fixed and random effects.
		#================================================================================
		###Resident
		#Pull out one replicate AND only use the N that are pre-invasion: 
		#Note, if you use data that are post invasion to fit the single-species resident 
		#model, LG_1sp, then alpha_ii will take on weird values. 
		#Choose 1: 
		mydata1 = subset(res_data, day_n < (inv_day) )
    	tryCatch( { 

    		all_models_daph_res[index1] = list( 
    			nlme(LG_res, fixed = list(lambda~1,alpha_ii~1), 
				random = list( replicate.number=lambda~1, replicate.number=alpha_ii~1),
        		start = c(lambda = 1.8, alpha_ii=0.5), 
        		correlation=corARMA(0.2, form=~day_n, p=1, q=0),data=mydata1)
    			)

	        #Use the fitted NLME to predict the average across all 3 mesocosms. 
	        #The "level=0" part of the argument will create a column with the across-
	        #replicate average. I have chosen to subset replicate.number == 4, but they
	        #are all identical.
	        #Each entry is a list object with two components: the days on which observations
	        #were made, and the predicted values on those days.  
	        model_fit1 = NULL
	        dn = unique(mydata1$day_n)
	        Nr= subset(predict(all_models_daph_res[[index1]],mydata1,level=0:1),replicate.number==4)$predict.fixed
	        model_fit_daph_res[[index1]] = data.frame(day_n = dn, N_res = Nr*mydata1$N_res )
			###A simple plot
			plot((mydata1$day_n),mydata1$N_res)
			points( unlist(model_fit_daph_res[[index1]]$day_n),unlist(model_fit_daph_res[[index1]]$N_res),col="red")
		

		}, error = function(e) {} ) 

		#Dump model output into a textfile! 
    	daph_res_out = capture.output(summary(all_models_daph_res[[index1]]))
		cat(paste("Daphnia as resident, T =", unique(mydata1$temperature)),  file="nlme_daphnia_res.txt", sep="n", append=TRUE)
		write.table(daph_res_out,file="nlme_daphnia_res.txt",sep = ",", quote = FALSE, row.names = F,append =T)

		#Pull out one replicate AND only use the N that are post-invasion but pre invader
		#equilibrium (this is very subjective, but also tends to be fairly robust) : 
		#Choose 1: 
		#mydata1 = subset(inv_data,replicate.number==6 & day_n > inv_day & day_n < (day_n+15) )
		#mydata1 = subset(inv_data,replicate.number==rep_use & day_n >= inv_day  )
		mydata1 = subset(inv_data, day_n >= inv_day  )

		tryCatch( { 

    		all_models_dia_invade[index1] = list( 
    			nlme(LG_inv, fixed = list(lambda~1,alpha_ij~1), 
				random = list( replicate.number=lambda~1, replicate.number=alpha_ij~1),
        		start = c(lambda = 4, alpha_ij=-0.1), 
        		correlation=corARMA(0.2, form=~day_n, p=1, q=0),data=mydata1)
    			)

	        #Use the fitted NLME to predict the average across all 3 mesocosms. 
	        #The "level=0" part of the argument will create a column with the across-
	        #replicate average. I have chosen to subset replicate.number == 4, but they
	        #are all identical.
	        #Each entry is a list object with two components: the days on which observations
	        #were made, and the predicted values on those days.  
	        model_fit1 = NULL
	        dn = unique(mydata1$day_n)
	        Ni= subset(predict(all_models_dia_invade[[index1]],mydata1,level=0:1),replicate.number==4)$predict.fixed
	        model_fit_dia_invade[[index1]] = data.frame(day_n = dn, N_inv = Ni*mydata1$N_inv )
			###A simple plot
			plot((mydata1$day_n),mydata1$N_inv)
			points( unlist(model_fit_dia_invade[[index1]]$day_n),unlist(model_fit_dia_invade[[index1]]$N_inv),col="red")
		}, error = function(e) {} ) 

		#Dump model output into a textfile! 
		dia_inv_out = capture.output(summary(all_models_dia_invade[[index1]]))
		cat(paste("Dia as invader, T =", unique(mydata1$temperature)), file="nlme_dia_inv.txt", sep="n", append=TRUE)
		write.table(dia_inv_out,file="nlme_dia_inv.txt",sep = ",", quote = FALSE, row.names = F,append =T)



}


#=============================================================================
# Second, loop over both diaphanosoma and daphnia when daphnia is the invader
#=============================================================================
#=============================================================================
# Second, loop over both diaphanosoma and daphnia when daphnia is the invader
#=============================================================================
all_models_daph_invade = vector( "list", ntreatments) #To store all of the fitted models
all_models_dia_res = vector( "list", ntreatments) #To store all of the fitted models

model_fit_daph_invade = vector( "list", ntreatments) #The predicted IGR based on the fitted LME 
model_fit_dia_res = vector( "list", ntreatments)

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
		index1 = n #+(i-1)*ntreatments

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
        res_tmp = na.exclude(res_tmp) #NLME won't work with NAs 
		
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
        inv_tmp = na.exclude(inv_tmp) #NLME won't work with NAs 
		#Make this a grouped data object:
		inv_data = groupedData(Ndiff_inv~N_res+N_inv+day_n|replicate.number, data = inv_tmp)
	

		#================================================================================
		#Use NLME to fit the model. This should work better than NLS since we can include 
		#AR correlated error and allow for fixed and random effects.
		#================================================================================
		###Resident
		#Pull out one replicate AND only use the N that are pre-invasion: 
		#Note: if you use data that are post invasion to fit the single-species resident 
		#model, LG_1sp, then alpha_ii will take on weird values.
		#Note: starting with alpha_ii = -0.1 will work for n=3 
		#Choose 1: 
		mydata1 = subset(res_data, day_n < (inv_day) )
    	tryCatch( { 

    		all_models_dia_res[index1] = list( 
    			nlme(LG_res, fixed = list(lambda~1,alpha_ii~1), 
				random = list( replicate.number=lambda~1, replicate.number=alpha_ii~1),
        		start = c(lambda = 1.8, alpha_ii=0.1), 
        		correlation=corARMA(0.2, form=~day_n, p=1, q=0),data=mydata1)
    			)

	        #Use the fitted NLME to predict the average across all 3 mesocosms. 
	        #The "level=0" part of the argument will create a column with the across-
	        #replicate average. I have chosen to subset replicate.number == 4, but they
	        #are all identical.
	        #Each entry is a list object with two components: the days on which observations
	        #were made, and the predicted values on those days.  
	        model_fit1 = NULL
	        dn = unique(mydata1$day_n)
	        Nr= subset(predict(all_models_dia_res[[index1]],mydata1,level=0:1),replicate.number==4)$predict.fixed
	        model_fit_dia_res[[index1]] = data.frame(day_n = dn, N_res = Nr*mydata1$N_res )
			###A simple plot
			plot((mydata1$day_n),mydata1$N_res)
			points( unlist(model_fit_dia_res[[index1]]$day_n),unlist(model_fit_dia_res[[index1]]$N_res),col="red")
		}, error = function(e) {} ) 

		#Dump model output into a textfile! 
    	dia_res_out = capture.output(summary(all_models_dia_res[[index1]]))
		cat(paste("Dia as resident, T =", unique(mydata1$temperature)),  file="nlme_dia_res.txt", sep="n", append=TRUE)
		write.table(dia_res_out,file="nlme_dia_res.txt",sep = ",", quote = FALSE, row.names = F,append =T)

		###Invader
		#Pull out one replicate AND only use the N that are post-invasion but pre invader
		#equilibrium (this is very subjective, but also tends to be fairly robust) : 
		#Choose 1: 
		mydata1 = subset(inv_data, day_n > inv_day & day_n < (day_n+15) )
		#mydata1 = subset(inv_data, day_n >= inv_day  )

		tryCatch( { 

    		all_models_daph_invade[index1] = list( 
    			nlme(LG_inv, fixed = list(lambda~1,alpha_ij~1), 
				random = list( replicate.number=lambda~1, replicate.number=alpha_ij~1),
        		start = c(lambda = 2, alpha_ij=0.1), 
        		correlation=corARMA(0.2, form=~day_n, p=1, q=0),data=mydata1)
    			)

	        #Use the fitted NLME to predict the average across all 3 mesocosms. 
	        #The "level=0" part of the argument will create a column with the across-
	        #replicate average. I have chosen to subset replicate.number == 4, but they
	        #are all identical.
	        #Each entry is a list object with two components: the days on which observations
	        #were made, and the predicted values on those days.  
	        model_fit1 = NULL
	        dn = unique(mydata1$day_n)
	        Ni= subset(predict(all_models_daph_invade[[index1]],mydata1,level=0:1),replicate.number==4)$predict.fixed
	        model_fit_daph_invade[[index1]] = data.frame(day_n = dn, N_inv = Ni*mydata1$N_inv )
			###A simple plot
			plot((mydata1$day_n),mydata1$N_inv)
			points( unlist(model_fit_daph_invade[[index1]]$day_n),unlist(model_fit_daph_invade[[index1]]$N_inv),col="red")
		}, error = function(e) {} ) 

		#Dump model output into a textfile! 
    	daph_inv_out = capture.output(summary(all_models_daph_invade[[index1]]))
		cat(paste("Daphnia as invader, T =", unique(mydata1$temperature)),  file="nlme_daphnia_inv.txt", sep="n", append=TRUE)
		write.table(daph_inv_out,file="nlme_daphnia_inv.txt",sep = ",", quote = FALSE, row.names = F,append =T)


}
