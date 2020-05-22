#=============================================================================
# R code to create to explore the Information Theoretic properties of 
# simple food webs. This creates a food web with an underlying dynamic
# model. 
#
# This version is meant to emulate the competitive dynamics of Daphnia and 
# Diaphanisoma measured in Coreen Forbes' 2019-2020 experiment.   	  
# 1. Food-web includes resource and two competitor species.
#	 A. Resource is a fluctuating supply source. 
#		1. Competition between consumers and resources emerges from consumption
#		2. Parameters at each level can be made a function of temperature. 
#	 B. Resources can be stochastic due to environmental fluctuations. 
# 2. Generate a series of foodwebs building from simple to complex structure
# 3. Perturb the food web in one of two ways: 
#	 A. Remove a species
#	 B. Add a pulse or press perturbation to a species growth rate
# 3. Use information theory to track the resulting food-web structures. 
# 4. This file has a lot of code for visualizing output of both the foodweb 
#	 its information theoretic properties after the main loop. 
#=============================================================================
#=============================================================================
# load libraries
#=============================================================================
library(deSolve)
library(fields)
library(gridExtra)
library(tidyverse)
library(plyr)
source("./functions/food_web_functions.R")
source("./functions/info_theory_functions.R")

#=============================================================================
# Outer loop. Set the number of trials and determine how to generate 
# combinations of species and parameters. 
#=============================================================================

#Length and time steps of each model run
tend = 20
delta1 = 0.01
tl=tend/delta1

#The maximum block depth for dynamic info metrics (larger is more accurate, but
#slower and could cause crashing if too large)
k= 2

###Build a series of scenarios going from simple to more complex dynamics
#Number of food webs to generate
nwebs = 1
# scenarios = list(matrix(0,nwebs,1))

###
#Output of each web
out1 = vector("list",nwebs)
#Invasion scenario
out_inv1 = vector("list",nwebs)

#Dynamic information metrics calculated from the (discretized) time series 
di_web = vector("list",nwebs)

#Track the average transfer entropy and separable information between each pair of 
#species as a way to build a network of information flow through the network. 
te_web = vector("list",nwebs)
si_web =vector("list",nwebs) 

#The ensemble version of the AI
aiE_web = vector("list",nwebs)

#Random resources:
 c1 = 10E6
 amp1 = 00000 #1/exp(1)
#Random consumers
 c2 = 1
 amp2 = 0.0 #1/exp(1)
 res_R = c(amp1,c1,amp2,c2)

 # c = 0
 # amp = 0 #1/exp(1)
 # res_R = c(amp,c)

#or 
# res_R = NULL

# scenarios[[1]] = list(nRsp = 1, nCsp =0, nPsp = 0)
# scenarios[[2]] = list(nRsp = 3, nCsp =0, nPsp = 0)
# scenarios[[3]] = list(nRsp = 3, nCsp =0, nPsp = 0)
# scenarios[[4]] = list(nRsp = 2, nCsp =1, nPsp = 0)
# scenarios[[5]] = list(nRsp = 1, nCsp =1, nPsp = 1)

#The structure of this code is based on taking an initial food web and 
#going through a series of perturbations. The "w" index now corresponds
#to each of the perturbations. 
w = 1 
#for (w in 1:nwebs){ 
print(w)

#Assume 2 trophic levels unless otherwise specified.
nRsp = 1 #Algae
nCsp = 2 #Spp 1 is Daphnia, Spp 2 is Diaphanosoma
nPsp = 1 #This is actually 0 --> Just a dummy predator
nspp = nRsp+nCsp+nPsp

#Randomly generate the species parameters for the model as well: 
spp_prms = NULL
#Resource: Nearly identical resource dynamics: 
spp_prms$rR = matrix(10E6, nRsp, 1) #intrinsic growth -- Not used here
spp_prms$Kr = matrix(10E6, nRsp, 1) #carrying capacity -- Constant algal addition

############################
#Consumers: 
# In the single species exponential growth phase the growth rate of each consumer rC
# is found by fitting an exponential model (e.g. cR_daph [[temp]] in resource_fits1.R).
# However, this exponent must ultimately reflect the consumption rate of algae, which
# means there is a conversion factor missing such that consumption_alga * X = exponent. 
# Solving for X should give the intrinsic growth rate rC. From
# resource_fits1.R, this means taking coef(cR_daph[[2]])[1]/coef(cl_daph[[2]])[1]. 
#
#spp_prms$rC = matrix( c(5,11) , nCsp, 1)  #intrisic growth - Solved fits in resource_fits1.R
spp_prms$rC = matrix( c(6.86e-5, 2.46e-4) , nCsp, 1)  #intrisic growth -  Direct fits in resource_fits1.R
#coef(cl_daph[[5]])[1]/(coef(cR_daph[[5]])[1]*10e6)
spp_prms$eFc = matrix(1,nCsp,nRsp) # just make the efficiency for everything 1 for now
spp_prms$muC = matrix(0.0, nCsp, 1) #matrix(rnorm(nCsp,0.6,0.1), nCsp, 1) #mortality rates
#Consumption rates: 
spp_prms$cC = matrix(c(0.035,0.015),nRsp,nCsp)
spp_prms$Kc = matrix(c(45,350), nCsp, 1) #carrying capacities, approximately matching data

# #Predators: These are just dummy variables for now
spp_prms$rP =  matrix(0.0, 1, 1) #matrix(rnorm(nPsp,0.5,0), nPsp, 1) #intrisic growth
spp_prms$eFp = matrix(1,1,nCsp) # just make the efficiency for everything 1 for now
spp_prms$muP = matrix(0.0, 1, 1)#matrix(rnorm(nPsp,0.6,0), nPsp, 1)  #mortality rates
#Consumption rates: 
spp_prms$cP = matrix(c(0.0,0.0),nCsp,1)

for (w in 1:nwebs){ 
	#=============================================================================
	# Outer loop. First run to equilibrate the population dynamics
	#=============================================================================
	#=============================================================================
	# This function gives: 
	# out 		The time series for of population growth for each species in the web
	#			This can be set to just give the final 2 time steps of the web with
	#			"final = TRUE"
	# spp_prms	The parameters of all species in the food web
	#=============================================================================
	#
	winit = matrix(c(10E6,1,1,0))
	tryCatch( {out1[w] = list(food_web_dynamics (spp_list = c(nRsp,nCsp,nPsp), spp_prms = spp_prms, 
		tend, delta1, winit = winit, res_R = res_R,final = FALSE ))}, error = function(e){}) 

	out1[[w]]$out[tl,]

	#=============================================================================
	# Inner loop: Mutual invasion:remove each species and track the dynamics
	#=============================================================================
	a_temp = NULL
	for (s in 1:2){

		out_temp =NULL
		out_temp2 =NULL
		inv_spp = s+1
		winit =  out1[[w]]$out[tl,2:(nspp+1)]
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


		#Competition from competitor
		#plot(out_temp2$out[1500:2000,(s+1)])
		#ts1=log(out_temp2$out[1500:2000,(s+1)])
		
		#"Competition" from predator
		# plot(out_temp2$out[50:300,(s+1)])
		# ts1=log(out_temp2$out[50:300,(s+1)])
		# ts2=1:length(ts1)
		# summary(lm(ts1~ts2))
	}

	out_inv1 [[w]] = a_temp
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
	  #Set k, the block size: 
	  k=2
	  #Get the populations of both species
	  f1=1 #scaling term
	  pop_ts = ceiling(f1*out_inv1[[w]])
	  
	  nt1 = 1
	  nt2 = dim(pop_ts)[1]
	  if(nt2 <=k){ k = 1}

	  if(nt1 != nt2) { 

	    f1 = 1 #scaling factor
	    di_web[w] = list(get_info_dynamics(pop_ts = pop_ts , k=k,with_blocks=TRUE))

	    ## This code takes the population time-series counts output by the ODEs and 
	    ## calculates the average Transfer Entropy from each species to every other 
	    ## species. The goal is to get an overview of the major information pathways 
	    ## in the web.   
	    #=============================================================================
	    # This function gives:
	    # te_web    Average transfer entropy per species as a pairwise matrix
	    #=============================================================================
	    te_web[w] = list( get_te_web( pop_ts = pop_ts, k=k) )

	    #=============================================================================
	    # This function gives:
	    # aiE_web    The AI of the entire ensemble, treated as a single time series. 
	    #=============================================================================
	    aiE_web[w] = list( get_ais (  series1 = pop_ts, k=k, ensemble = TRUE)    )

	    #=============================================================================  
	    #Build these back out into a data frame that includes all of the mesocosm,
	    #treatment, and species information. 

	    #Add the DIT to the data frames. There will be 11 new columns. 
	    # DIT_tmp = matrix(0,nt2,11)
	    # ncnames = c("N_res","N_inv","aiE","te1","te2","ee1","ee2","ai1","ai2","si1","si2")
	    # DIT_tmp[,1:2] = as.matrix(pop_ts)
	    # DIT_tmp[(k+1):nt2,3] = aiE_web[[index1]]$local
	    # DIT_tmp[(k+1):nt2,4:5] = di_web[[index1]]$te_local
	    # DIT_tmp[(k*2):nt2,6:7] = di_web[[index1]]$ee_local
	    # DIT_tmp[(k+1):nt2,8:9] = di_web[[index1]]$ai_local
	    # DIT_tmp[(k+1):nt2,10:11] = di_web[[index1]]$si_local
	    # colnames(DIT_tmp) = ncnames
	    # DIT_tmp = as.data.frame(DIT_tmp) 
	    # out1[[index1]] = out1[[index1]] %>% left_join(DIT_tmp)

	}
}

#=============================================================================
#Post-loop processing to create a single data frame for plotting and stats
#=============================================================================
##Recursively flatten this into a data frame. 
mDIT_tmp = data.frame(matrix( nrow=0, ncol =( (nRsp+nCsp)*5+4) )) 
ncnames = c("Algae", "Daphnia","Diaphanosoma","aiE","te1","te2", "te3","ee1","ee2",
	"ee3","ai1","ai2","ai3","si1","si2","si3","run","mesocosm","nspp")
colnames(mDIT_tmp) = ncnames
mesocosms = factor(c("A","B")) #A is Daphnia invader, B is Daphnia resident
nspp = factor(c(1,2)) #1 is the pre-invasion phase, 2 is post-invasion phase
for (f in 1:nwebs){
	 DIT_tmp = data.frame(matrix(0,nt2,( (nRsp+nCsp)*5+4) ))
     DIT_tmp[,1:3] = as.matrix(out_inv1[[f]][,2:4])
     DIT_tmp[(k+1):nt2,4] = aiE_web[[f]]$local
     DIT_tmp[(k+1):nt2,5:7] = di_web[[f]]$te_local[,2:4]
     DIT_tmp[(k*2):nt2,8:10] = di_web[[f]]$ee_local[,2:4]
     DIT_tmp[(k+1):nt2,11:13] = di_web[[f]]$ai_local[,2:4]
     DIT_tmp[(k+1):nt2,14:16] = di_web[[f]]$si_local[,2:4]
     DIT_tmp[,17] = f

     DIT_tmp[,18] = factor(levels = levels(mesocosms))
     DIT_tmp[1:nt2/2,18] = mesocosms[1]
     DIT_tmp[(nt2/2+1):nt2,18] = mesocosms[2]

     spp1 = c(1:(nt2/4),(nt2/2+1):(nt2/2+nt2/4))
     spp2 = c((nt2/4+1):(nt2/2),(nt2/2+nt2/4+1):nt2)
     DIT_tmp[,19] = factor(levels = levels(nspp))
     DIT_tmp[spp1,19] = nspp[1]
     DIT_tmp[spp2,19] = nspp[2]

     colnames(DIT_tmp) = ncnames
     DIT_tmp = as.data.frame(DIT_tmp) 
     mDIT_tmp = rbind(mDIT_tmp,DIT_tmp)
}

#Add the time column
#mDIT = cbind(time = matrix(seq(0,tend*2+delta1,delta1),dim(mDIT_tmp)[1],1), mDIT_tmp)
mDIT = cbind(time = matrix(seq(0,tend,delta1),dim(mDIT_tmp)[1],1), mDIT_tmp)
a1 = exp((mDIT$Algae-lag(mDIT$Algae))/(
	(mDIT$Daphnia+mDIT$Diaphanosoma)-( lag(mDIT$Daphnia)+lag(mDIT$Diaphanosoma) )))
a1[is.na(a1)] = 0
a1[is.infinite(a1)] = 0

#=============================================================================
#Saving
#=============================================================================
#save(file="daphDia_DIT_100.var", "mDIT", "out1","out_inv1","di_web",
#	"te_web","si_web", "aiE_web")
# save(file="daphDia_DIT_20_v01.var", "mDIT", "out1","out_inv1","di_web",
# 	"te_web","si_web", "aiE_web")
#save(file="daphDia_DIT_20_v00.var", "mDIT", "out1","out_inv1","di_web",
#	"te_web","si_web", "aiE_web")
#=============================================================================
#Loading
#=============================================================================
#Load saved files
#load("daphDia_DIT_100.var")
#load("daphDia_DIT_20_v00B.var")
#load("daphDia_DIT_20_v01.var")
#=============================================================================
#Add some new columns to the data: 
#=============================================================================
mDIT$alg_perDaph = (c(spp_prms$Kr)-mDIT$Algae)/mDIT$Daphnia
mDIT$alg_perDaph[is.infinite(mDIT$alg_perDaph)] = 0

mDIT$alg_perDia = (c(spp_prms$Kr)-mDIT$Algae)/mDIT$Diaphanosoma
mDIT$alg_perDia[is.infinite(mDIT$alg_perDia)] = 0

mDIT$zoo = mDIT$Daphnia+mDIT$Diaphanosoma

mDIT$alg_perzoo = (c(spp_prms$Kr)-mDIT$Algae)/mDIT$zoo
mDIT$alg_perzoo[is.infinite(mDIT$alg_perzoo)] = 0

#Delta columns: 
mDIT$alg_per_DDaph = (c(spp_prms$Kr)-mDIT$Algae)/c(NA,diff(mDIT$Daphnia))
	#(mDIT$Daphnia *spp_prms$cC[1,1]*spp_prms$rC[1]*mDIT$Algae)
mDIT$alg_per_DDaph[is.infinite(mDIT$alg_per_DDaph)] = 0

mDIT$alg_per_DDia = (c(spp_prms$Kr)-mDIT$Algae)/c(NA,diff(mDIT$Diaphanosoma))
	#(mDIT$Dia *spp_prms$cC[1,2]*spp_prms$rC[2]*mDIT$Algae)
mDIT$alg_per_DDia[is.infinite(mDIT$alg_per_DDia)] = 0

mDIT$alg_per_Dzoo = (c(spp_prms$Kr)-mDIT$Algae)/c(NA,diff(mDIT$zoo))
mDIT$alg_per_Dzoo[is.infinite(mDIT$alg_per_Dzoo)] = 0


mDIT=ddply(mDIT, mesocosm~nspp, transform,
 m_alg_perDaph  = mean(alg_perDaph ), m_alg_perzoo  = mean(alg_perzoo ),
 m_Daph  = mean(Daphnia ), m_zoo  = mean(zoo )) 

# mDIT1 = mDIT
# mDIT2 = mDIT


#=============================================================================
#Plots
#=============================================================================
# ggplot()+ geom_point(data= mDIT, mapping= aes(x = aiE, y =(Algae-lag(Algae))/( 
# 	(Daphnia+Diaphanosoma)-( lag(Daphnia)+lag(Diaphanosoma) )  ) ) )+
#Just pop plots
ggplot()+ geom_point(data= mDIT, mapping= aes(x = time, y =Daphnia,  color = time ) )+  
	geom_point(data= mDIT, mapping= aes(x = time, y =Diaphanosoma,  color = time ) )+  
	facet_grid(mesocosm~nspp)  +  ylim(50000,250000)
ggplot()+ geom_point(data= mDIT[mDIT$time<5,], mapping= aes(x = time, y =Daphnia,  color = time ) )+  
	facet_grid(mesocosm~nspp)   + ylim(50000,250000)

#Pop and aiE
ggplot()+ geom_point(data= mDIT[mDIT$time<2,], mapping= aes(x = time, y =Daphnia,  color = time ) )+  
	geom_line(data= mDIT[mDIT$time<2,], mapping= aes(x = time, y =ai1,  color = time ) )+  
	facet_grid(mesocosm~nspp) 


ggplot()+ geom_line(data= mDIT[mDIT$time<1.5,], mapping= aes(x = time, y =Daphnia,color="1") )+  
		  geom_point(data= mDIT[mDIT$time<1.5,], mapping= aes(x = time, y =ai1,color="2" ) )+ 
		  geom_line(data= mDIT[mDIT$time<1.5,], mapping= aes(x = time, y =(alg_perDaph)/2E5,color="3"))+
		  scale_color_discrete(name ="", labels = c("Daphnia", "AI", "algal cost per Daph"))+
facet_grid(mesocosm~nspp)
ggsave("./Daph_ai1_algalperDaph_think1.pdf", width = 8, height = 10)


ggplot()+geom_point(data= mDIT[mDIT$time<2,], mapping= aes(x = ai1, y =(alg_perDaph)/2E5) )+
facet_grid(mesocosm~nspp) 
ggsave("./ai1_algalperDaph_think1.pdf", width = 8, height = 10)


#When the var in Zoo is 0, aiE is driven by Algae: 
ggplot()+ geom_point(data= mDIT, mapping= aes(x = time, y =Algae, 
 	color = time ) )+  facet_grid(mesocosm~nspp)

#Relative to perturbation of algae consumed per Daphnia
#Raw with mean
#Normal
ggplot()+ geom_point(data= mDIT, mapping= aes(x = aiE, y =(alg_perDaph), 
 	color = time ) )+ 
	geom_line(data= mDIT, mapping= aes(x = aiE, y =m_alg_perDaph, 
 	color = time ))  + facet_grid(mesocosm~nspp) 
#Difference
ggplot()+ geom_point(data= mDIT, mapping= aes(x = aiE, y =(alg_per_DDaph), 
 	color = time ) )+ 
	geom_line(data= mDIT, mapping= aes(x = aiE, y =m_alg_perDaph, 
 	color = time ))  + facet_grid(mesocosm~nspp) 

#Pop vs. consumption
#Normal
ggplot()+ geom_point(data= mDIT, mapping= aes(x = Daphnia, y =(alg_perDaph), 
 	color = time ) )+ facet_grid(mesocosm~nspp) 
#Difference
ggplot()+ geom_point(data= mDIT, mapping= aes(x = Daphnia, y =(alg_per_DDaph), 
 	color = time ) )+ facet_grid(mesocosm~nspp) 

#********************
#aiE vs. consumption
#*********************
ggplot()+ geom_point(data= mDIT2, mapping= aes(x = aiE, y =alg_perDaph, 
 	color = time ) )+
	geom_line(data= mDIT, mapping= aes(x = aiE, y =m_alg_perDaph, 
 	color = time ))+
 	facet_grid(mesocosm~nspp)+  ylim(50000,250000)
#total
ggplot()+ geom_point(data= mDIT[mDIT$time>2,], mapping= aes(x = aiE, y =alg_perzoo, 
 	color = time ) )+
	geom_line(data= mDIT[mDIT$time>2,], mapping= aes(x = aiE, y =m_alg_perzoo, 
 	color = time ))+
 	facet_grid(mesocosm~nspp)  +  ylim(50000,250000)

#multiple
ggplot()+ geom_point(data= mDIT1, mapping= aes(x = aiE, y =alg_perDaph, 
 	color = time ) )+
	geom_line(data= mDIT1, mapping= aes(x = aiE, y =m_alg_perDaph, 
 	color = time ))  + geom_point(data= mDIT2, mapping= aes(x = aiE, y =alg_perDaph, 
 	color = time ) )+  
	geom_line(data= mDIT2, mapping= aes(x = aiE, y =m_alg_perDaph, 
 	color = time ))  +
 	facet_grid(mesocosm~nspp)
 	# +  ylim(50000,250000)
# ggplot()+ geom_point(data= mDIT[mDIT$time>2,], mapping= aes(x = aiE, y =alg_perDaph, 
#  	color = time ) )+  facet_grid(mesocosm~nspp) + 
# 	geom_line(data= mDIT[mDIT$time>2,], mapping= aes(x = aiE, y =m_alg_perDaph, 
#  	color = time ))+ylim(50000,250000)
#Absolute, mean subtracted
ggplot()+ geom_point(data= mDIT[mDIT$time>2,], mapping= aes(x = aiE, y =(alg_perDaph-m_alg_perDaph), 
 	color = time ) )+  facet_grid(mesocosm~nspp)

#************************************************
#More complex the further it is perturbed from K? 
#************************************************
#Absolute, mean subtracted
# ggplot()+ geom_point(data= mDIT[mDIT$time>2,], mapping= aes(x = aiE, y =abs(m_Daph-Daphnia), 
#  	color = time ) )+  facet_grid(mesocosm~nspp)
ggplot()+ geom_point(data= mDIT[mDIT$time>2,], mapping= aes(x = aiE, y =abs(m_Daph-Daphnia), 
 	color = time ) )+  facet_grid(mesocosm~nspp)
ggplot()+ geom_point(data= mDIT[mDIT$time>2,], mapping= aes(x = aiE, y =(spp_prms$Kc[1,1]-Daphnia), 
 	color = time ) )+  facet_grid(mesocosm~nspp)
#Both spp
ggplot()+ geom_point(data= mDIT[mDIT$time>2,], mapping= aes(x = aiE, y =((spp_prms$Kc[1,1]+spp_prms$Kc[2,1])-zoo), 
 	color = time ) )+  facet_grid(mesocosm~nspp)

#Look at N vs. AIE
ggplot()+ geom_point(data= mDIT, mapping= aes(x =Daphnia, y = ai1, 
 	color = time ) )+  facet_grid(mesocosm~nspp)
ggplot()+ geom_point(data= mDIT[mDIT$time>2,], mapping= aes(x =Daphnia, y = aiE, 
 	color = time ) )+  facet_grid(mesocosm~nspp) + 
	# geom_line(data= mDIT[mDIT$time>2,], mapping= aes(x = aiE, y =m_Daph, 
 # 	color = time ))+ylim(50000,250000)

ggplot()+ geom_point(data= mDIT[mDIT$time>2,], mapping= aes(x = aiE, y =Daphnia, 
 	color = time ) )+  facet_grid(mesocosm~nspp) + 
	geom_line(data= mDIT[mDIT$time>2,], mapping= aes(x = aiE, y =m_Daph, 
 	color = time ))+ylim(50000,250000)

#Fot total zooplankton: 
#Consumption
ggplot()+ geom_point(data= mDIT[mDIT$time>2,], mapping= aes(x = aiE, y =(alg_perzoo-m_alg_perzoo), 
 	color = time ) )+  facet_grid(mesocosm~nspp)
#Zoo biomass:
ggplot()+ geom_point(data= mDIT, mapping= aes(x = aiE, y =abs(m_zoo-zoo), 
 	color = time ) )+  facet_grid(mesocosm~nspp)

ggplot()+ geom_point(data= mDIT1, mapping= aes(x = aiE, y =abs(m_zoo-zoo), 
 	color = time ) )+ geom_point(data= mDIT2, mapping= aes(x = aiE, y =abs(m_zoo-zoo), 
 	color = time ) ) +facet_grid(mesocosm~nspp)


ggplot()+ geom_point(data= mDIT, mapping= aes(x = aiE, y =(Algae-lag(Algae))/( 
	(Daphnia- lag(Daphnia) )  ) ) )+
ggplot()+ geom_point(data= mDIT, mapping= aes(x = aiE, y =(c(spp_prms$Kr)-Algae)/( 
	(Daphnia )  ),  color = interaction(run,mesocosm,nspp ) ) )+
ggplot()+ geom_point(data= mDIT, mapping= aes(x = time, y =(c(spp_prms$Kr)-Algae)/( 
	(Daphnia )  ),  color = run ) )+  facet_grid(mesocosm~nspp)   + ylim(70000,90000)



ggplot()+ geom_point(data= mDIT[mDIT$time>2,], mapping= aes(x = aiE, y =abs( ( (c(spp_prms$Kr)-Algae)/Daphnia ) -
	mean( (c(spp_prms$Kr)-Algae)/Daphnia )  ),  color = time ) )+  facet_grid(mesocosm~nspp)   + ylim(50000,250000)



ggplot()+ geom_line(data= mDIT, mapping= aes(x = time, y =Daphnia,  
	color = run ) )+  facet_grid(mesocosm~nspp)   +
geom_line(data= mDIT, mapping= aes(x = time, y =aiE,  
	color = run ) )+  facet_grid(mesocosm~nspp) +xlim(0,10)
ggplot()+ geom_line(data= mDIT, mapping= aes(x = time, y =Daphnia,
	color = interaction(run,mesocosm,nspp ) ) )   +
	geom_point(data= mDIT, mapping= aes(x = time, y =aiE ) )+
  xlab("Active information (bits) ")+
  ylab("Algal consumption per individual")+
  theme(strip.background = element_rect(colour=NA, fill=NA))
ggsave("./aiE_algalperN_theory.pdf", width = 8, height = 10)


#=============================================================================
#Plots in combination with the real data
# m1_DIT is the real data
# mDIT is the simulation data
#=============================================================================
m1_DIT_sub = subset(m1_DIT, species == "daphnia" & temperature == 28 & invade_monoculture == "daph invade") #Daphnia at 28 C
mDIT_sub = subset(mDIT, mesocosm == "A" & nspp == 2) #Daphnia invasion 

p1 = ggplot()+ geom_line(data=m1_DIT_sub, mapping= aes(x= (day_n-min(day_n)+1), y=N, group = mesocosm_id,color="1" )) +
	geom_line(data=mDIT_sub[mDIT_sub$time<3.5,], mapping= aes(x = time*100, y =Daphnia,  color = "2") )+
	scale_color_discrete(name ="", labels = c("Experiment", "Simulation" ) )+
  	ylab("Population (Daphnia")+
  	theme(axis.title.x=element_blank(),axis.text.x = element_blank(), axis.ticks = element_blank())
p2 = ggplot() + geom_point(data=m1_DIT_sub, mapping= aes(x= (day_n-min(day_n)+1), y=ai1, color="1" )) +
	geom_point(data= mDIT_sub[mDIT_sub$time<3.5,], mapping= aes(x = time*100, y =ai2,  color = "2" ) ) +
	scale_color_discrete(name ="", labels = c("Experiment", "Simulation" ) )+
  	ylab("Bits (AI) ")+
  	theme(axis.title.x=element_blank(),axis.text.x = element_blank(), axis.ticks = element_blank())
p3 = ggplot()+ geom_line(data=m1_DIT_sub, mapping= aes(x= (day_n-min(day_n)+1), y=alg_per_N, group = mesocosm_id,,color="1" )) +
	geom_line(data= mDIT_sub[mDIT_sub$time<3.5,], mapping= aes(x = time*100, y =(alg_perDaph), color="2") )+
	scale_color_discrete(name ="", labels = c("Experiment", "Simulation" ) )+
	xlab("Day")+
  	ylab("Algae consumed per N ")

gp1 = ggplotGrob(p1)
gp2 = ggplotGrob(p2)
gp3 = ggplotGrob(p3)
grid.draw(rbind(gp1, gp2,gp3))

#To save as a PDF
pdf("./Daph_exp_v_simv1_a1_all_long.pdf",width = 8, height = 10)
grid.draw(rbind(gp1, gp2,gp3))
dev.off()


#################
#efficiency, 3 panels

m1_DIT_sub = subset(m1_DIT, species == "daphnia" & temperature == 28 ) #Daphnia at 28 C
mDIT_sub = subset(mDIT, mesocosm == "A" & nspp == 2) #Daphnia invasion 

p1 = ggplot()+ geom_point(data=m1_DIT_sub, mapping= aes(x= aiE, y=alg_per_N, group = mesocosm_id,color="1" )) +
	geom_point(data= mDIT_sub[mDIT_sub$time<3.5,], mapping= aes(x = aiE, y =(alg_perDaph), color="2") ) +
	scale_color_discrete(name ="", labels = c("Experiment", "Simulation" ) )+
	xlab("Bits (AI)")+
  	ylab("Algae consumed per N ")

p2 = ggplot()+ geom_point(data=m1_DIT_sub, mapping= aes(x= aiE, y=abs(mean(N) - N), group = mesocosm_id,color="1" )) +
	geom_point(data= mDIT_sub[mDIT_sub$time<3.5,], mapping= aes(x = aiE, y =abs(mean(Daphnia)-Daphnia), color="2") ) +
	scale_color_discrete(name ="", labels = c("Experiment", "Simulation" ) )+
	xlab("Bits (AI)")+
  	ylab("abs(Distance from K) ")

p3 = ggplot()+ geom_point(data=m1_DIT_sub, mapping= aes(x= N, y=alg_per_N, group = mesocosm_id,color="1" )) +
	geom_point(data= mDIT_sub[mDIT_sub$time<3.5,], mapping= aes(x = Daphnia, y =(alg_perDaph), color="2") ) +
	scale_color_discrete(name ="", labels = c("Experiment", "Simulation" ) )+
	xlab("Population N")+
  	ylab("Algae consumed per N ")

gp1 = ggplotGrob(p1)
gp2 = ggplotGrob(p2)
gp3 = ggplotGrob(p3)
grid.draw(rbind(gp1, gp2,gp3))
	

p3 = ggplot()+ geom_point(data=m1_DIT_sub, mapping= aes(x= N, y=ai1, group = mesocosm_id,color="1" )) +
	geom_point(data= mDIT_sub[mDIT_sub$time<3.5,], mapping= aes(x = Daphnia, y =ai2, color="2") ) +
	scale_color_discrete(name ="", labels = c("Experiment", "Simulation" ) )+
	xlab("Day")+
  	ylab("Algae consumed per N ")

ggplot()+ geom_point(data= mDIT[mDIT$time>2,], mapping= aes(x = aiE, y =abs(m_Daph-Daphnia), 
 	color = time ) )+  facet_grid(mesocosm~nspp)
ggplot()+ geom_point(data= mDIT[mDIT$time>2,], mapping= aes(x = aiE, y =(spp_prms$Kc[1,1]-Daphnia), 
 	color = time ) )

####
#For Diaphanosoma:
m1_DIT_sub = subset(m1_DIT, species == "diaphanosoma" & temperature == 28& invade_monoculture == "dia invade") #Daphnia at 28 C
mDIT_sub = subset(mDIT, mesocosm == "B" & nspp == 2) #Daphnia invasion 

p1 = ggplot()+ geom_line(data=m1_DIT_sub, mapping= aes(x= (day_n-min(day_n)+1), y=N, group = mesocosm_id,color="1" )) +
	geom_line(data=mDIT_sub[mDIT_sub$time<3.5,], mapping= aes(x = time*100, y =Diaphanosoma,  color = "2") )+
	scale_color_discrete(name ="", labels = c("Experiment", "Simulation" ) )+
  	ylab("Population (Daphnia")+
  	theme(axis.title.x=element_blank(),axis.text.x = element_blank(), axis.ticks = element_blank())
p2 = ggplot() + geom_point(data=m1_DIT_sub, mapping= aes(x= (day_n-min(day_n)+1), y=ai2, color="1" )) +
	geom_point(data= mDIT_sub[mDIT_sub$time<3.5,], mapping= aes(x = time*100, y =ai3,  color = "2" ) ) +
	scale_color_discrete(name ="", labels = c("Experiment", "Simulation" ) )+
  	ylab("Bits (AI) ")+
  	theme(axis.title.x=element_blank(),axis.text.x = element_blank(), axis.ticks = element_blank())
p3 = ggplot()+ geom_line(data=m1_DIT_sub, mapping= aes(x= (day_n-min(day_n)+1), y=alg_per_N, group = mesocosm_id,,color="1" )) +
	geom_line(data= mDIT_sub[mDIT_sub$time<3.5,], mapping= aes(x = time*100, y =(alg_perDaph), color="2") )+
	scale_color_discrete(name ="", labels = c("Experiment", "Simulation" ) )+
	xlab("Day")+
  	ylab("Algae consumed per N ")

gp1 = ggplotGrob(p1)
gp2 = ggplotGrob(p2)
gp3 = ggplotGrob(p3)
grid.draw(rbind(gp1, gp2,gp3))

#To save as a PDF
pdf("./Dia_exp_v_simv1_a1_all_long.pdf",width = 8, height = 10)
grid.draw(rbind(gp1, gp2,gp3))
dev.off()



