#=============================================================================
#Implementation of a lotka-volterra competition model 
#=============================================================================
# load libraries
#=============================================================================
library(tidyverse)
library(deSolve)
library(gridExtra)

#=============================================================================
# Define the population dynamics through the following function
#=============================================================================
# This is model v.1, where the intrinsic growth rate is specified as ri, so it
# is obvious. 
lvs = function(times,sp,parms){
	with( as.list(c(parms, sp )),
		{		
		
			nspp = parms$nspp 
			Ni = matrix(sp[1:nspp], nspp, 1)

			dNi = Ni 
			for( i in 1:nspp){
				#Lotka-volterra consumption with a saturating mortality term 
				dNi[i] = Ni[i] * ( ri[i] *(1 - t(alphas[i,])%*%Ni )  )
			}

	  	list( c(dNi) )
		})	

}  

# This is model v.2, which corresponds to eq. 3 in Chesson 1990, with bi and ki. 
lvs2 = function(times,sp,parms2 ){
	with( as.list(c(parms, sp )),
		{		
		
			nspp = parms$nspp 
			Ni = matrix(sp[1:nspp], nspp, 1)

			dNi = Ni 
			for( i in 1:nspp){
				#Lotka-volterra consumption with a saturating mortality term 
				dNi[i] = Ni[i] * (bi[i] * (ki[i] - t(alphas2[i,])%*%Ni ) )
			}

	  	list( c(dNi) )
		})	

}  

#This is model v.3, which still corresponds to eq. 3 in Chesson 1990, but 
#where the model has been written out fully nstead of in terms of the 
#aggregate parameters ki and alphas. 
# This is model v.2, which corresponds to eq. 3 in Chesson 1990, with bi and ki. 
#NOT IMPLEMENTED YET
# lvs3 = function(times,sp,parms2 ){
# 	with( as.list(c(parms, sp )),
# 		{		
		
# 			nspp = parms$nspp 
# 			Ni = matrix(sp[1:nspp], nspp, 1)

# 			dNi = Ni 
# 			for( i in 1:nspp){
# 				#Lotka-volterra consumption with a saturating mortality term 
# 				dNi[i] = Ni[i] * (bi[i] * (cil[i]*wl[i]*Ki -mi[i] - t(alphas2[i,])%*%Ni ) )
# 			}

# 	  	list( c(dNi) )
# 		})	

# }  

#=============================================================================
# Set values of parameters
#=============================================================================
#How many species to start? 
nspp = 2

spp_prms = NULL
#spp_prms$ri = matrix(rpois(nspp,10), nspp, 1) #intrinsic growth
spp_prms$ri = matrix(c(5,2), nspp, 1) #intrinsic growth

#model 2: There are a few ways to set this, just to compare with 
#model 1. 

#This scales bi and ki to match ri and alphas in model 1
spp_prms$ki =  matrix(c(3,3), nspp, 1)
spp_prms$bi =  spp_prms$ri/spp_prms$ki
 
#Make ki equal to ri, and set bi to 1
# spp_prms$bi =  matrix( c(1,1),2, 1)
# spp_prms$ki =  spp_prms$ri

#Make bi equal to ri, and set ki to 1
# spp_prms$ki =  matrix(c(1,1), nspp, 1)
# spp_prms$bi =  spp_prms$ri

#Competition coefficients: 
#spp_prms$alphas = matrix( runif(nspp^2,  0.1, 0.9),nspp,nspp )
#diag(spp_prms$alphas) = matrix(rnorm(nspp, 0.99,0.01),nspp,1) #Set the intraspecific alphas = 1
spp_prms$alphas = matrix( c(1,0.8,0.8,1),nspp,nspp )
spp_prms$alphas2 = matrix(spp_prms$ki,nspp,nspp)*spp_prms$alphas

#Pass all of these parameters as a list
parms = list(
	nspp=nspp, ri = spp_prms$ri, alphas2 = spp_prms$alphas2, 
	alphas = spp_prms$alphas, bi = spp_prms$bi, ki = spp_prms$ki
 )

#=============================================================================
# Run the model with initial conditions 
#=============================================================================
tend = 100
delta1 = 0.1
times  = seq(from = 0, to = tend, by = delta1)
tl = length(times)
minit = c( matrix(0.001,nspp,1) )

#This is model1
lv_out = dede(y=minit, times=times, func=lvs, parms=parms, atol = 1e-9)
lv_out = as.data.frame(lv_out)

#This is model2
lv_out2 = dede(y=minit, times=times, func=lvs2, parms=parms, atol = 1e-9)
lv_out2 = as.data.frame(lv_out2)


#=============================================================================
# Plots!
#=============================================================================
lv_long =lv_out %>% gather( species, N, 2:(nspp+1) )
lv_long2 =lv_out2 %>% gather( species, N, 2:(nspp+1) )

#Total number of species still around at each time: 
lv_long = lv_long %>% group_by(time) %>% mutate(spp_tot = sum(N>0.001) )
lv_long2 = lv_long2 %>% group_by(time) %>% mutate(spp_tot = sum(N>0.001) )

#Their collective biomass: 
lv_long = lv_long %>% group_by(time) %>% mutate(biomass = sum(N[N>0.001] ) )
lv_long2 = lv_long2%>% group_by(time) %>% mutate(biomass = sum(N[N>0.001] ) )


#Each species' time trajectory
p1=ggplot()+ geom_line( data = lv_long, aes ( x = time, y = N, color = species)  )+ 
ylab("Population")+
theme(axis.text.x=element_blank(), axis.title.x=element_blank(), legend.position = "none") 

p2=ggplot()+ geom_line( data = lv_long2, aes ( x = time, y = N, color = species)  )+ 
ylab("Population")+
theme(axis.text.x=element_blank(), axis.title.x=element_blank(), legend.position = "none") 

# #Total number of species
# p2=ggplot()+ geom_line( data = lv_long, aes ( x = time, y = spp_tot)  )+
# ylab("Number of species")+
# theme(axis.text.x=element_blank(), axis.title.x=element_blank(), legend.position = "none") 

# #Total biomass of community
# p3=ggplot()+ geom_line( data = lv_long, aes ( x = time, y = biomass)  )+
# ylab("Biomass")

p4 = grid.arrange(p1,p2, nrow = 2 )
p4