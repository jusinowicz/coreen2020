
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

#=============================================================================
# First, loop over both diaphanosoma and daphnia when diaphanasoma is the 
# invader
#=============================================================================

all_models_dia_invade = vector( "list", ntreatments*2) #To store all of the fitted models
model_fit_dia_invade = vector( "list", ntreatments*2) #The predicted IGR based on the fitted LME 

#Loop over both diaphanosoma and daphnia when diaphanasoma is the invader
for (i in 1:2) {  
	#Loop over treatments 
	u_invader = invader[i]

	for (n in 1: (ntreatments)){
		#Since the treatments are regular, set the start/end position of
		#the subset qith:  
		pos1 = (n-1)*2 + n
		pos2 = pos1+2
		u_treats = treats[pos1:pos2]
		index1 = n+(i-1)*ntreatments
		#================================================================================
		#The LME approach: simple but wrong
		#================================================================================
		# #Parse out the subset of 3 mesocosms based on the invader ID and Temp combo
		# mydata = subset(m1_data_long, mesocosm.ID == u_treats &  species == u_invader  )
		# #Use LME to fit an AR(1) model with replicate number as a random effect
		# all_models[n] = list (nlme(N ~ day_n, random= ~ 1|replicate.number, 
		# 	correlation=corARMA(0.2, form=~day_n, p=1, q=0), na.action=na.omit, data=mydata))
		# #Now use the fitted LME to predict the average across all 3 mesocosms
		# model_fit[,n] = predict(all_models[[n]],mydata)
		#================================================================================
		#The NLME approach: better, but complicated
		#================================================================================
		#Parse out the subset of 3 mesocosms based on the invader ID and Temp combo
		mydata = subset(m1_data_long, mesocosm.ID == u_treats &  species == u_invader  )
		#Make this a grouped data set for NLME. This groups by replicate number
        mydata2 = groupedData(N~day_n|replicate.number, data = mydata)
        mydata2 = na.exclude(mydata2) #NLME won't work with NAs 
        #Use NLME to fit an AR(1) model with replicate number as a random effect
        # all_models[n] = list(nlme(N ~ SSasymp(day_n, Asym, R0, lrc),fixed = Asym+R0+lrc~1, 
        # 	random =  Asym ~1, correlation=corARMA(0.2, form=~day_n, p=1, q=0),
        # 	start = c(unlist(getInitial(N ~ SSasymp(day_n, Asym, R0, lrc),data=mydata2) )), 
        # 	data=mydata2))

    	tryCatch( { 

    		all_models_dia_invade[index1] = list(nlme(N ~ SSlogis(day_n, Asym, R0, lrc),fixed = Asym+R0+lrc~1, 
        	random =  Asym ~1, correlation=corARMA(0.2, form=~day_n, p=1, q=0),
        	start = c(unlist(getInitial(N ~ SSlogis(day_n, Asym, R0, lrc),data=mydata2) )), 
        	data=mydata2))

	        #Now use the fitted NLME to predict the average across all 3 mesocosms. 
	        #The "level=0" part of the argument will create a column with the across-
	        #replicate average. I have chosen to subset replicate.number == 4, but they
	        #are all identical.
	        #Each entry is a list object with two components: the days on which observations
	        #were made, and the predicted values on those days.  
	        model_fit1 = NULL
	        model_fit1$days = list(unique(mydata2$day_n))
	        model_fit1$N = list(subset(predict(all_models_dia_invade[[index1]],mydata2,level=0:1),replicate.number==4)$predict.fixed)
	        model_fit_dia_invade[[index1]] = model_fit1
			###A simple plot
			plot((mydata2$day_n),mydata2$N)
			lines( unlist(model_fit_dia_invade[[index1]]$days),unlist(model_fit_dia_invade[[index1]]$N))
		}, error = function(e) {} ) 
	}
}


#=============================================================================
# Second, loop over both diaphanosoma and daphnia when daphnia is the invader
#=============================================================================
all_models_dap_invade = vector( "list", ntreatments*2) #To store all of the fitted models
model_fit_dap_invade = vector( "list", ntreatments*2) #The predicted IGR based on the fitted LME 

#Loop over both diaphanosoma and daphnia when daphnia is the invader
for (i in 1:2) {  
	#Loop over treatments 
	u_invader = invader[i]

	for (n in 1: (ntreatments)){
		#Since the treatments are regular, set the start/end position of
		#the subset qith:  
		pos1 = (n-1)*2 + n+invasions_per
		pos2 = pos1+2
		u_treats = treats[pos1:pos2]
		index1 = n+(i-1)*ntreatments
		#================================================================================
		#The LME approach: simple but wrong
		#================================================================================
		# #Parse out the subset of 3 mesocosms based on the invader ID and Temp combo
		# mydata = subset(m1_data_long, mesocosm.ID == u_treats &  species == u_invader  )
		# #Use LME to fit an AR(1) model with replicate number as a random effect
		# all_models[n] = list (nlme(N ~ day_n, random= ~ 1|replicate.number, 
		# 	correlation=corARMA(0.2, form=~day_n, p=1, q=0), na.action=na.omit, data=mydata))
		# #Now use the fitted LME to predict the average across all 3 mesocosms
		# model_fit[,n] = predict(all_models[[n]],mydata)
		#================================================================================
		#The NLME approach: better, but complicated
		#================================================================================
		#Parse out the subset of 3 mesocosms based on the invader ID and Temp combo
		mydata = subset(m1_data_long, mesocosm.ID == u_treats &  species == u_invader  )
		#Make this a grouped data set for NLME. This groups by replicate number
        mydata2 = groupedData(N~day_n|replicate.number, data = mydata)
        mydata2 = na.exclude(mydata2) #NLME won't work with NAs 
        #Use NLME to fit an AR(1) model with replicate number as a random effect
        # all_models[n] = list(nlme(N ~ SSasymp(day_n, Asym, R0, lrc),fixed = Asym+R0+lrc~1, 
        # 	random =  Asym ~1, correlation=corARMA(0.2, form=~day_n, p=1, q=0),
        # 	start = c(unlist(getInitial(N ~ SSasymp(day_n, Asym, R0, lrc),data=mydata2) )), 
        # 	data=mydata2))

    	tryCatch( { 

    		all_models_dap_invade[index1] = list(nlme(N ~ SSlogis(day_n, Asym, R0, lrc),fixed = Asym+R0+lrc~1, 
        	random =  Asym ~1, correlation=corARMA(0.2, form=~day_n, p=1, q=0),
        	start = c(unlist(getInitial(N ~ SSlogis(day_n, Asym, R0, lrc),data=mydata2) )), 
        	data=mydata2))

	        #Now use the fitted NLME to predict the average across all 3 mesocosms. 
	        #The "level=0" part of the argument will create a column with the across-
	        #replicate average. I have chosen to subset replicate.number == 4, but they
	        #are all identical.
	        #Each entry is a list object with two components: the days on which observations
	        #were made, and the predicted values on those days.  
	        model_fit1 = NULL
	        model_fit1$days = list(unique(mydata2$day_n))
	        model_fit1$N = list(subset(predict(all_models_dap_invade[[index1]],mydata2,level=0:1),replicate.number==4)$predict.fixed)
	        model_fit_dap_invade[[index1]] = model_fit1
			###A simple plot
			plot((mydata2$day_n),mydata2$N)
			lines( unlist(model_fit_dap_invade[[index1]]$days),unlist(model_fit_dap_invade[[index1]]$N))
		}, error = function(e) {} ) 
	}
}
