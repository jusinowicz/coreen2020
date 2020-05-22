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
source("./info_theory_functions.R")
#=============================================================================
#Data sorting and preprocessing: 
#=============================================================================

#=============================================================================
#alg_replace()
#Fit a GAMM to the data to predict algal abundance on the basis of density 
#per-species, as a funciton of temperature, with mesocosm_id as a Random Effect.
#Note that the data are fit with the log link function (log transformed) because 
#values cannot be <0. 
#=============================================================================
alg_replace = function(m1_data_long){   

  alg_gam = gam( algae_abundance ~ s(N)+s(temperature,k=5)+s(species, bs="re")+s(mesocosm_id,bs="re"),
    family=Gamma(link='log'), data=m1_data_long )

  #Use the fitted GAMM to interpolate NAs for algal abundance in the data set
  alg_newdata = subset(m1_data_long, is.na(algae_abundance) )
  alg_replace = as.vector(predict.gam(alg_gam, newdata=alg_newdata,type= "response") )
  alg_newdata$algae_abundance = alg_replace
  return(alg_newdata)

}

#################
#There are two different approaches to fitting this. For now, just 
#comment/uncomment based on the approach: 
#
#1. Increasing logistic model of total consumption vs. N
#algae_start = mean(m1_data_long$algae_cells_mL,na.rm=T)

#2 Decreasing logsitic model of algae abundance vs. N
#3 Decreasing exponential model of algae abundance / algae_start vs. N
#algae_start = max(m1_data_long$algae_abundance,na.rm=T)

for(t in 1:ntemps) { 

  #=============================================================================
  #Intraspecific and resource consumption:

  #Get the data for each temperature. 
  daph_tmp = subset(mesos_daph, temperature == temps[t])
  dia_tmp = subset(mesos_dia, temperature == temps[t])


  #Remove NA entries! 
  # daph_tmp = daph_tmp[!is.na(daph_tmp$Adiff) & !is.na(daph_tmp$N), ]
  # dia_tmp = dia_tmp[!is.na(dia_tmp$Adiff) & !is.na(dia_tmp$N), ]

  daph_tmp = daph_tmp[!is.na(daph_tmp$N), ]
  dia_tmp = dia_tmp[ !is.na(dia_tmp$N), ]

  #Add the rate of population growth, Ndiff
  dp1 = subset(daph_tmp,!is.na(daph_tmp$N) ) %>%    
            arrange(replicate_number) %>%
            mutate(Ndiff = (N-lag(N))/(day_n-lag(day_n))*1/N )  #"lead" lines up the result
            #mutate(Ndiff = (N/lag(N))/(day_n/lag(day_n)) ) #"lead" lines up the result 
  dp1$Ndiff[is.infinite(dp1$Ndiff)] = NA

  daph_tmp = daph_tmp %>% left_join(dp1)

  da1 = subset(dia_tmp,!is.na(dia_tmp$N) ) %>%    
            arrange(replicate_number) %>%
            mutate(Ndiff = (N-lag(N))/(day_n-lag(day_n))*1/N ) #"lead" lines up the result
            #mutate(Ndiff = (N/lag(N))/(day_n/lag(day_n)) ) #"lead" lines up the result 
 
  da1$Ndiff[is.infinite(da1$Ndiff)] = NA

  dia_tmp = dia_tmp %>% left_join(da1)


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

  #=============================================================================
  #Interspecific competition
  daph_inv_tmp = subset(mesos_inv_daph, temperature == temps[t])
  dia_inv_tmp = subset(mesos_inv_dia, temperature == temps[t])

  #Creates specific columns designating invader and resident: 
  dp2i = subset(daph_inv_tmp, species == "daphnia")
  dp2r = subset(daph_inv_tmp, species == "diaphanosoma")
  colnames(dp2i)[colnames(dp2i)=="N"] = "N_inv"
  colnames(dp2r)[colnames(dp2r)=="N"] = "N_res"
  daph_inv_tmp = dp2i %>% 
    left_join(dp2r[,c("day_n","mesocosm_id","N_res")], by=c("day_n","mesocosm_id"))
  rm(dp2i,dp2r)

  #Add the rate of population growth, Ndiff
  dp1 = subset(daph_inv_tmp,!is.na(daph_inv_tmp$N_inv) ) %>%    
            arrange(replicate_number) %>%
            mutate(Ndiff = (N_inv-lag(N_inv))/(day_n-lag(day_n))*1/N_inv) #"lead" lines up the result 
            #mutate(Ndiff = (N/lag(N))/(day_n/lag(day_n)) ) #"lead" lines up the result 

  dp1$Ndiff[is.infinite(dp1$Ndiff)] = NA
  daph_inv_tmp = daph_inv_tmp %>% left_join(dp1)
  rm(dp1)

  #Creates specific columns designating invader and resident: 
  dp2i = subset(dia_inv_tmp, species == "daphnia")
  dp2r = subset(dia_inv_tmp, species == "diaphanosoma")
  colnames(dp2i)[colnames(dp2i)=="N"] = "N_res"
  colnames(dp2r)[colnames(dp2r)=="N"] = "N_inv"
  dia_inv_tmp = dp2r %>% 
    left_join(dp2i[,c("day_n","mesocosm_id","N_res")], by=c("day_n","mesocosm_id"))
  rm(dp2i,dp2r)

  #Add the rate of population growth, Ndiff
  da1 = subset(dia_inv_tmp,!is.na(dia_inv_tmp$N_inv) ) %>%    
            arrange(replicate_number) %>%
            mutate(Ndiff = (N_inv-lag(N_inv))/(day_n-lag(day_n))*1/N_inv ) #"lead" lines up the result 
            #mutate(Ndiff = (N/lag(N))/(day_n/lag(day_n)) ) #"lead" lines up the result 
  da1$Ndiff[is.infinite(da1$Ndiff)] = NA
  dia_inv_tmp = dia_inv_tmp %>% left_join(da1)
  rm(da1)

  print(mean(as.data.frame(subset(daph_inv_tmp,day_n <=34))$Ndiff,na.rm=T) )
  print(mean(as.data.frame(subset(dia_inv_tmp,day_n <=34))$Ndiff,na.rm=T) )
  #For plotting: 
  cl_daph_plot = rbind( cl_daph_plot, daph_tmp )
  cl_dia_plot = rbind( cl_dia_plot, dia_tmp )

  #The basic MacArthur model is a linear consumption rate so just fit with a GLMM
  # cl_daph[[t]] = gam( Ndiff ~ +s(temperature,k=5)+s(mesocosm_id,bs="re"),family=Gamma(link='log'), data=daph_tmp)
  # cl_dia[[t]] = gam( Ndiff ~ +s(temperature,k=5)+s(mesocosm_id,bs="re"),family=Gamma(link='log'), data=dia_tmp)
  # plot(daph_tmp$algae_abundance, daph_tmp$N)
  # plot(dia_tmp$algae_abundance, dia_tmp$N)

  #=============================================================================
  #This works better!!! 
  #Fit the algal consumption rate instead: 
  # cl_daph[[t]] = gam( Adiff ~  N+s(mesocosm_id,bs="re"), family=binomial, data=daph_tmp)
  # cl_dia[[t]] = gam( Adiff ~  N +s(mesocosm_id,bs="re"), famly=binomial, data=dia_tmp)

 

  # #1. Use NLS to fit a Type 2 (saturating) response
  # alg2 = formula (Adiff ~  algae_start/(1+c1*exp(b1*N) ) )

  # m1 = lm(I(log(algae_start/Adiff-1))~I(N), data = daph_tmp ) 
  # tryCatch({ 
  #   cl_daph[[t]] = nls( formula= alg2, data = daph_tmp, 
  #   start=list(b1=(as.numeric(coef(m1)[2])), c1=exp(as.numeric(coef(m1)[1] ) ) ),
  #   control=nls.control(maxiter = 1000), trace=F ) 

  #   #Predicted values: 
  #   s=seq(min(daph_tmp$N),max(daph_tmp$N),1 )
  #   d_tmp = predict(cl_daph[[t]], list( N = s ) )
  #   daph_pred_tmp = data.frame( species = rspecies[1], temperature = temps[t], s=s, N_pred = d_tmp )
  #   cl_daph_pred = rbind(cl_daph_pred, daph_pred_tmp)
  # }, error = function(e) {} ) 

  # m1 = lm(I(log(algae_start/Adiff-1))~I(N), data = dia_tmp ) 
  # tryCatch( {
  #   cl_dia[[t]] = nls( formula= alg2, data = dia_tmp, 
  #   start=list(b1=(as.numeric(coef(m1)[2])), c1=exp(as.numeric(coef(m1)[1] ) ) ),
  #   control=nls.control(maxiter = 1000), trace=F ) 

  #   #Predicted values: 
  #   s=seq(min(dia_tmp$N),max(dia_tmp$N),1 )
  #   d_tmp = predict(cl_dia[[t]], list( N = s ) )
  #   dia_pred_tmp = data.frame( species = rspecies[2], temperature = temps[t], s=s, N_pred = d_tmp )
  #   cl_dia_pred = rbind(cl_dia_pred, dia_pred_tmp)
  # }, error = function(e) {} ) 

  ######Note: Approaches 2 and 3 give essentially the same fit of consumption rates c1
 # #2. Use NLS to fit the decaying logistic/exponential response
 # alg2 = formula (algae_abundance/algae_start ~  1/(1+a1*exp(c1*N) ) )

 #  m1 = lm(I(log(algae_start/algae_abundance))~I(N), data = daph_tmp ) 
 #  tryCatch({ 
 #    cl_daph[[t]] = nls( formula= alg2, data = daph_tmp, 
 #    start=list(c1=(as.numeric(coef(m1)[2])), a1=exp(as.numeric(coef(m1)[1] ) ) ),
 #    control=nls.control(maxiter = 1000), trace=F ) 

 #    #Predicted values: 
 #    #s=seq(min(daph_tmp$N),max(daph_tmp$N),1 )
 #    s=seq(0,maxN,1 )
 #    d_tmp = predict(cl_daph[[t]], list( N = s ) )
 #    daph_pred_tmp = data.frame( species = rspecies[1], temperature = temps[t], s=s, N_pred = d_tmp )
 #    cl_daph_pred = rbind(cl_daph_pred, daph_pred_tmp)
 #  }, error = function(e) {} ) 

 #  m1 = lm(I(log(algae_start/algae_abundance))~I(N), data = dia_tmp ) 
 #  tryCatch( {
 #    cl_dia[[t]] = nls( formula= alg2, data = dia_tmp, 
 #    start=list(c1=(as.numeric(coef(m1)[2])), a1=exp(as.numeric(coef(m1)[1] ) ) ),
 #    control=nls.control(maxiter = 1000), trace=F ) 

 #    #Predicted values: 
 #    #s=seq(min(dia_tmp$N),max(dia_tmp$N),1 )
 #    s=seq(0,maxN,1 )
 #    d_tmp = predict(cl_dia[[t]], list( N = s ) )
 #    dia_pred_tmp = data.frame( species = rspecies[2], temperature = temps[t], s=s, N_pred = d_tmp )
 #    cl_dia_pred = rbind(cl_dia_pred, dia_pred_tmp)
 #  }, error = function(e) {} ) 


# #3. Use NLS to fit a decaying exponential response
  alg2 = formula (algae_abundance/algae_start~  exp(c1*N+a1) )

  m1 = lm(I(log(algae_abundance/algae_start ))~I(N), data = daph_tmp ) 
  tryCatch({ 
    cl_daph[[t]] = nls( formula= alg2, data = daph_tmp, 
    start=list(c1=(as.numeric(coef(m1)[2])),a1=(as.numeric(coef(m1)[1] ) )   ),
    control=nls.control(maxiter = 1000), trace=F ) 

    #Predicted values: 
    #s=seq(min(daph_tmp$N),max(daph_tmp$N),1 )
    s=seq(0,maxN,1 )
    d_tmp = predict(cl_daph[[t]], list( N = s ) )
    daph_pred_tmp = data.frame( species = rspecies[1], temperature = temps[t], s=s, N_pred = d_tmp )
    cl_daph_pred = rbind(cl_daph_pred, daph_pred_tmp)
  }, error = function(e) {} ) 

  m1 = lm(I(log(algae_abundance/algae_start))~I(N), data = dia_tmp ) 
  tryCatch( {
    cl_dia[[t]] = nls( formula= alg2, data = dia_tmp, 
    start=list(c1=(as.numeric(coef(m1)[2])),a1=(as.numeric(coef(m1)[1] ) )   ),
    control=nls.control(maxiter = 1000), trace=F ) 

    #Predicted values: 
    #s=seq(min(dia_tmp$N),max(dia_tmp$N),1 )
    s=seq(0,maxN,1 )
    d_tmp = predict(cl_dia[[t]], list( N = s ) )
    dia_pred_tmp = data.frame( species = rspecies[2], temperature = temps[t], s=s, N_pred = d_tmp )
    cl_dia_pred = rbind(cl_dia_pred, dia_pred_tmp)
  }, error = function(e) {} ) 

#3b. Use NLS to fit a decaying exponential response
  # alg2 = formula (Adiff ~  exp(c1*N+a1) )

  # m1 = lm(I(log(Adiff ))~I(N), data = daph_tmp ) 
  # tryCatch({ 
  #   cl_daph[[t]] = nls( formula= alg2, data = daph_tmp, 
  #   start=list(c1=(as.numeric(coef(m1)[2])),a1=(as.numeric(coef(m1)[1] ) )   ),
  #   control=nls.control(maxiter = 1000), trace=F ) 

  #   #Predicted values: 
  #   #s=seq(min(daph_tmp$N),max(daph_tmp$N),1 )
  #   s=seq(0,maxN,1 )
  #   d_tmp = predict(cl_daph[[t]], list( N = s ) )
  #   daph_pred_tmp = data.frame( species = rspecies[1], temperature = temps[t], s=s, N_pred = d_tmp )
  #   cl_daph_pred = rbind(cl_daph_pred, daph_pred_tmp)
  # }, error = function(e) {} ) 

  # m1 = lm(I(log(Adiff ))~I(N), data = dia_tmp ) 
  # tryCatch( {
  #   cl_dia[[t]] = nls( formula= alg2, data = dia_tmp, 
  #   start=list(c1=(as.numeric(coef(m1)[2])),a1=(as.numeric(coef(m1)[1] ) )   ),
  #   control=nls.control(maxiter = 1000), trace=F ) 

  #   #Predicted values: 
  #   #s=seq(min(dia_tmp$N),max(dia_tmp$N),1 )
  #   s=seq(0,maxN,1 )
  #   d_tmp = predict(cl_dia[[t]], list( N = s ) )
  #   dia_pred_tmp = data.frame( species = rspecies[2], temperature = temps[t], s=s, N_pred = d_tmp )
  #   cl_dia_pred = rbind(cl_dia_pred, dia_pred_tmp)
  # }, error = function(e) {} ) 

  ############################  
  #The inverse relationship between N and Adiff, related to approach 1: 
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

  #=============================================================================
  #Get the intrinsic growth rate of each species
  # cR = formula (N ~  I(a1*exp(c1*day_n ) ) )

  # c1 = lm(I(log(N+1))~I(day_n ), data = daph_tmp ) 
  # tryCatch({ 
  #   cR_daph[[t]] = nls( formula= cR, data = daph_tmp, 
  #   start=list(c1=(as.numeric(coef(c1)[2])), a1=exp(as.numeric(coef(c1)[1] ) ) ),
  #   control=nls.control(maxiter = 1000), algorithm="port", trace=F ) 

  #   #Predicted values: 
  #   s=seq(min(daph_tmp$day_n),max(daph_tmp$day_n),1 )
  #   #s=seq(0,maxT,1 )
  #   d_tmp = predict(cR_daph[[t]], list( day_n = s ) )
  #   daph_pred_tmp = data.frame( species = rspecies[1], temperature = temps[t], s=s, N_pred = d_tmp )
  #   cR_daph_pred = rbind(cR_daph_pred, daph_pred_tmp)
  # }, error = function(e) {} ) 

  # c1 = lm(I(log(N+1))~I(day_n ), data = dia_tmp ) 
  # tryCatch({ 
  #   cR_dia[[t]] = nls( formula= cR, data = dia_tmp, 
  #   start=list(c1=(as.numeric(coef(c1)[2])), a1=exp(as.numeric(coef(c1)[1] ) ) ),
  #   control=nls.control(maxiter = 1000), algorithm="port", trace=F ) 

  #   #Predicted values: 
  #   s=seq(min(dia_tmp$day_n),max(dia_tmp$day_n),1 )
  #   #s=seq(0,maxT,1 )
  #   d_tmp = predict(cR_dia[[t]], list( day_n = s ) )
  #   dia_pred_tmp = data.frame( species = rspecies[2], temperature = temps[t], s=s, N_pred = d_tmp )
  #   cR_dia_pred = rbind(cR_dia_pred, dia_pred_tmp)
  # }, error = function(e) {} ) 

  #=============================================================================

  #1. Get the intrinsic growth rate of each species as a function of algal density
  # cR = formula (N~  I(exp(c1*Adiff+a1) ) )

  # c1 = lm(I(log(N+1))~I(Adiff), data = daph_tmp ) 
  # tryCatch({ 
  #   cR_daph[[t]] = nls( formula= cR, data = daph_tmp, 
  #   start=list(c1=(as.numeric(coef(c1)[2])), a1=(as.numeric(coef(c1)[1] ) ) ),
  #   control=nls.control(maxiter = 1000), algorithm="port", trace=F ) 
  # }, error = function(e) {} ) 
  #   #Predicted values: 
  #   s=seq(min(daph_tmp$Adiff),max(daph_tmp$Adiff),1 )
  #   #s=seq(0,maxT,1 )

  #   if( is.null(cR_daph[[t]]) ) {
  #       cR_daph[[t]] = c1
  #       d_tmp = exp(predict(cR_daph[[t]], list( Adiff= s ) ) )
  #   }else {
  #       d_tmp = predict(cR_daph[[t]], list( Adiff= s ) )
  #   }

  #   daph_pred_tmp = data.frame( species = rspecies[1], temperature = temps[t], s=s, N_pred = d_tmp )
  #   cR_daph_pred = rbind(cR_daph_pred, daph_pred_tmp)
 

  # c1 = lm(I(log(N+1))~I(Adiff), data = dia_tmp ) 
  # tryCatch({ 
  #   cR_dia[[t]] = nls( formula= cR, data = dia_tmp, 
  #   start=list(c1=(as.numeric(coef(c1)[2])), a1=(as.numeric(coef(c1)[1] ) ) ),
  #   control=nls.control(maxiter = 1000), algorithm="port", trace=F ) 
  # }, error = function(e) {} ) 
  #   #Predicted values: 
  #   s=seq(min(dia_tmp$Adiff),max(dia_tmp$Adiff),1 )
  #   #s=seq(0,maxT,1 )
  #   if( is.null(cR_dia[[t]]) ) {
  #       cR_dia[[t]] = c1
  #       d_tmp = exp(predict(cR_dia[[t]], list( Adiff= s ) ) )
  #   }else {
  #       d_tmp = predict(cR_dia[[t]], list( Adiff= s ) )
  #   }
  #   dia_pred_tmp = data.frame( species = rspecies[2], temperature = temps[t], s=s, N_pred = d_tmp )
  #   cR_dia_pred = rbind(cR_dia_pred, dia_pred_tmp)


# #2. Get the intrinsic growth rate of each species as a function of algal density
#   cR = formula (N~  I(exp(c1*p1*Adiff+a1) ) )

#   p1 = abs(coef(cl_daph[[w]])[1])
#   c1 = lm(I(log(N+1))~I(p1*Adiff), data = daph_tmp ) 
#   tryCatch({ 
#     cR_daph[[t]] = nls( formula= cR, data = daph_tmp, 
#     start=list(c1=(as.numeric(coef(c1)[2])), a1=(as.numeric(coef(c1)[1] ) ) ),
#     control=nls.control(maxiter = 1000), algorithm="port", trace=F ) 
#   }, error = function(e) {} ) 
#     #Predicted values: 
#     s=seq(min(daph_tmp$Adiff),max(daph_tmp$Adiff),1 )
#     #s=seq(0,maxT,1 )

#     if( is.null(cR_daph[[t]]) ) {
#         cR_daph[[t]] = c1
#         d_tmp = exp(predict(cR_daph[[t]], list( Adiff= s ) ) )
#     }else {
#         d_tmp = predict(cR_daph[[t]], list( Adiff= s ) )
#     }

#     daph_pred_tmp = data.frame( species = rspecies[1], temperature = temps[t], s=s, N_pred = d_tmp )
#     cR_daph_pred = rbind(cR_daph_pred, daph_pred_tmp)
 


#   p1 = abs(coef(cl_dia[[w]])[1])
#   c1 = lm(I(log(N+1))~I(p1*Adiff), data = dia_tmp ) 
#   tryCatch({ 
#     cR_dia[[t]] = nls( formula= cR, data = dia_tmp, 
#     start=list(c1=(as.numeric(coef(c1)[2])), a1=(as.numeric(coef(c1)[1] ) ) ),
#     control=nls.control(maxiter = 1000), algorithm="port", trace=F ) 
#   }, error = function(e) {} ) 
#     #Predicted values: 
#     s=seq(min(dia_tmp$Adiff),max(dia_tmp$Adiff),1 )
#     #s=seq(0,maxT,1 )
#     if( is.null(cR_dia[[t]]) ) {
#         cR_dia[[t]] = c1
#         d_tmp = exp(predict(cR_dia[[t]], list( Adiff= s ) ) )
#     }else {
#         d_tmp = predict(cR_dia[[t]], list( Adiff= s ) )
#     }
#     dia_pred_tmp = data.frame( species = rspecies[2], temperature = temps[t], s=s, N_pred = d_tmp )
#     cR_dia_pred = rbind(cR_dia_pred, dia_pred_tmp)

#3. Get the intrinsic growth rate of each species as a function of algal density
  cR = formula ( N ~  exp(c1*algae_abundance+a1)  )

  c1 = lm(I(log(N+1))~I(algae_abundance), data = daph_tmp ) 
  tryCatch({ 
    cR_daph[[t]] = nls( formula= cR, data = daph_tmp, 
    start=list(c1=(as.numeric(coef(c1)[2])), a1=(as.numeric(coef(c1)[1] ) ) ),
    control=nls.control(maxiter = 1000), algorithm="port", trace=F ) 
  }, error = function(e) {} ) 
    #Predicted values: 
    s=seq(min(daph_tmp$algae_abundance),max(daph_tmp$algae_abundance),1 )
    #s=seq(0,maxT,1 )

    if( is.null(cR_daph[[t]]) ) {
        cR_daph[[t]] = c1
        d_tmp = exp(predict(cR_daph[[t]], list( algae_abundance= s ) ) )
    }else {
        d_tmp = predict(cR_daph[[t]], list( algae_abundance= s ) )
    }

    daph_pred_tmp = data.frame( species = rspecies[1], temperature = temps[t], s=s, N_pred = d_tmp )
    cR_daph_pred = rbind(cR_daph_pred, daph_pred_tmp)
 

  c1 = lm(I(log(N+1))~I(algae_abundance), data = dia_tmp ) 
  tryCatch({ 
    cR_dia[[t]] = nls( formula= cR, data = dia_tmp, 
    start=list(c1=(as.numeric(coef(c1)[2])), a1=(as.numeric(coef(c1)[1] ) ) ),
    control=nls.control(maxiter = 1000), algorithm="port", trace=F ) 
  }, error = function(e) {} ) 
    #Predicted values: 
    s=seq(min(dia_tmp$algae_abundance),max(dia_tmp$algae_abundance),1 )
    #s=seq(0,maxT,1 )
    if( is.null(cR_dia[[t]]) ) {
        cR_dia[[t]] = c1
        d_tmp = exp(predict(cR_dia[[t]], list( algae_abundance= s ) ) )
    }else {
        d_tmp = predict(cR_dia[[t]], list( algae_abundance= s ) )
    }
    dia_pred_tmp = data.frame( species = rspecies[2], temperature = temps[t], s=s, N_pred = d_tmp )
    cR_dia_pred = rbind(cR_dia_pred, dia_pred_tmp)

  #=============================================================================
  #Fit intraspecific competition coefficients directly using the single species and invasion
  #scenarios.
  #Get the intrinsic growth rate of each species
  # igr = formula (N ~  I(a1*exp(c1*day_n ) ) )
  igr = formula (N ~  a1/(1+b1*exp(c1*day_n ) ) )
 
  #Adjust the days to make both data sets line up
  daph_tmpb=daph_tmp
  daph_tmpb$day_n[daph_tmpb$day_n>inv_day] = daph_tmpb$day_n[daph_tmpb$day_n>inv_day]-inv_day 

  #Fit the model
  tryCatch({ 
    igr_daph[[t]] = nls( formula= igr, data = daph_tmpb, 
    #start=list(c1=(as.numeric(coef(igr1)[2])), a1=exp(as.numeric(coef(igr1)[1] ) ) ),
    start=list(c1=.05, a1=max(daph_tmpb$N),b1=1 ),
    control=nls.control(maxiter = 1000), algorithm="port", trace=F ) 

    #Predicted values: 
    s=seq(min(daph_tmp$day_n),max(daph_tmp$day_n),1 )
    #s=seq(0,maxT,1 )
    d_tmp = predict(igr_daph[[t]], list( day_n = s ) )
    igr_daph_tmp = data.frame( species = rspecies[1], temperature = temps[t], s=s, N_pred = d_tmp )
    igr_daph_pred = rbind(igr_daph_pred, igr_daph_tmp)
  }, error = function(e) {} ) 

  #Make new columns out of the fitted data
  igr_daph_tmp$N = igr_daph_tmp$N_pred
  igr_daph_tmp$Ndiff = (igr_daph_tmp$N_pred-lag(igr_daph_tmp$N_pred))/
                        (igr_daph_tmp$s-lag(igr_daph_tmp$s))*1/igr_daph_tmp$N_pred

  #The fitted data should be straightforward: 
  lvii_daph[[t]] = lm(data=igr_daph_tmp, Ndiff ~ N)
  s=seq(1,max(daph_tmp$N),1 )
  d_tmp = predict(lvii_daph[[t]], list( N = s ) )
  lvii_daph_tmp = data.frame( species = rspecies[1], temperature = temps[t], s=s, N_pred = d_tmp )
  lvii_daph_pred = rbind(lvii_daph_pred, lvii_daph_tmp)


  #Dia 
  #Adjust the days to make both data sets line up
  dia_tmpb=dia_tmp
  dia_tmpb$day_n[dia_tmpb$day_n>inv_day] = dia_tmpb$day_n[dia_tmpb$day_n>inv_day]-inv_day 

  #Fit the model
  tryCatch({ 
    igr_dia[[t]] = nls( formula= igr, data = dia_tmpb, 
    #start=list(c1=(as.numeric(coef(igr1)[2])), a1=exp(as.numeric(coef(igr1)[1] ) ) ),
    start=list(c1=.05, a1=max(dia_tmpb$N),b1=1 ),
    control=nls.control(maxiter = 1000), algorithm="port", trace=F ) 

    #Predicted values: 
    s=seq(min(dia_tmp$day_n),max(dia_tmp$day_n),1 )
    #s=seq(0,maxT,1 )
    d_tmp = predict(igr_dia[[t]], list( day_n = s ) )
    igr_dia_tmp = data.frame( species = rspecies[2], temperature = temps[t], s=s, N_pred = d_tmp )
    igr_dia_pred = rbind(igr_dia_pred, igr_dia_tmp)
  }, error = function(e) {} ) 

  #Make new columns out of the fitted data
  igr_dia_tmp$N = igr_dia_tmp$N_pred
  igr_dia_tmp$Ndiff = (igr_dia_tmp$N_pred-lag(igr_dia_tmp$N_pred))/
                        (igr_dia_tmp$s-lag(igr_dia_tmp$s))*1/igr_dia_tmp$N_pred

  #The fitted data should be straightforward: 
  lvii_dia[[t]] = lm(data=igr_dia_tmp, Ndiff ~ N)
  s=seq(1,max(dia_tmp$N),1 )
  d_tmp = predict(lvii_dia[[t]], list( N = s ) )
  lvii_dia_tmp = data.frame( species = rspecies[2], temperature = temps[t], s=s, N_pred = d_tmp )
  lvii_dia_pred = rbind(lvii_dia_pred, lvii_dia_tmp)
  #=============================================================================
  #Fit interspecific competition coefficients directly using the invasion
  #scenarios.
  #igrij = formula (N_inv ~  (b1*exp(c1*day_n ) ) )
  igrij = formula (N_inv ~  a1/(1+b1*exp(c1*day_n ) ) )
 
  a1 = max(daph_inv_tmp$N_inv,na.rm=T )
  #Fit the model
  tryCatch({ 
    igrij_daph[[t]] = nls( formula= igrij, data = daph_inv_tmp, 
    #start=list(c1=(as.numeric(coef(igr1)[2])), a1=exp(as.numeric(coef(igr1)[1] ) ) ),
    start=list(c1=-.05, b1=2 ),
    control=nls.control(maxiter = 1000), algorithm="port", trace=F ) 

    #Predicted values: 
    s=seq(min(daph_inv_tmp$day_n),max(daph_tmp$day_n),1 )
    #s=seq(0,maxT,1 )
    d_tmp = predict(igrij_daph[[t]], list( day_n = s ) )
    igrij_daph_tmp = data.frame( species = rspecies[1], temperature = temps[t], s=s, N_pred = d_tmp )
    igrij_daph_pred = rbind(igrij_daph_pred, igrij_daph_tmp)
  }, error = function(e) {} ) 

  #Make new columns out of the fitted data
  igrij_daph_tmp$N = igrij_daph_tmp$N_pred
  igrij_daph_tmp$Ndiff = (igrij_daph_tmp$N_pred-lag(igrij_daph_tmp$N_pred))/
                        (igrij_daph_tmp$s-lag(igrij_daph_tmp$s))*1/igrij_daph_tmp$N_pred

  #The fitted data should be straightforward: 
  lvij_daph[[t]] = lm(data=igrij_daph_tmp, Ndiff ~ N)
  s=seq(1,max(daph_tmp$N),1 )
  d_tmp = predict(lvij_daph[[t]], list( N = s ) )
  lvij_daph_tmp = data.frame( species = rspecies[1], temperature = temps[t], s=s, N_pred = d_tmp )
  lvij_daph_pred = rbind(lvij_daph_pred, lvij_daph_tmp)

  #Dia

  a1 = max(dia_inv_tmp$N_inv,na.rm=T )
  #Fit the model
  tryCatch({ 
    igrij_dia[[t]] = nls( formula= igrij, data = dia_inv_tmp, 
    #start=list(c1=(as.numeric(coef(igr1)[2])), a1=exp(as.numeric(coef(igr1)[1] ) ) ),
    start=list(c1=-.05, b1=2 ),
    control=nls.control(maxiter = 1000), algorithm="port", trace=F ) 

    #Predicted values: 
    s=seq(min(dia_inv_tmp$day_n),max(dia_tmp$day_n),1 )
    #s=seq(0,maxT,1 )
    d_tmp = predict(igrij_dia[[t]], list( day_n = s ) )
    igrij_dia_tmp = data.frame( species = rspecies[2], temperature = temps[t], s=s, N_pred = d_tmp )
    igrij_dia_pred = rbind(igrij_dia_pred, igrij_dia_tmp)
  }, error = function(e) {} ) 

  #Make new columns out of the fitted data
  igrij_dia_tmp$N = igrij_dia_tmp$N_pred
  igrij_dia_tmp$Ndiff = (igrij_dia_tmp$N_pred-lag(igrij_dia_tmp$N_pred))/
                        (igrij_dia_tmp$s-lag(igrij_dia_tmp$s))*1/igrij_dia_tmp$N_pred

  #The fitted data should be straightforward: 
  lvij_dia[[t]] = lm(data=igrij_dia_tmp, Ndiff ~ N)
  s=seq(1,max(dia_tmp$N),1 )
  d_tmp = predict(lvij_dia[[t]], list( N = s ) )
  lvij_dia_tmp = data.frame( species = rspecies[2], temperature = temps[t], s=s, N_pred = d_tmp )
  lvij_dia_pred = rbind(lvij_dia_pred, lvij_dia_tmp)
  
  #=============================================================================
#   #Fit intraspecific competition coefficients directly using the single species and invasion
#   #scenarios.
#   #1. Lotka Volterra

#   #f_lvii = formula ( I(log(Ndiff)) ~  ri*(ki-aii*N) )
#   #f_lvii = formula ( Ndiff ~  ri*(1-aii*N) )
#   f_lvii = formula ( Ndiff ~  ri/(1+aii*N) )


#   #Daphnia
#   lm_ii = lm(data=igr_daph_tmp, Ndiff ~ N)
#   #ri=coef(lm_ii)[1]
#   #ri=1
#   tryCatch({ 
#     lvii_daph[[t]] = nls( formula= f_lvii, data = igr_daph_tmp, 
#     start=list(ri=2, aii=0.5 ),
#     control=nls.control(maxiter = 1000), algorithm="port", trace=F ) 
#   }, error = function(e) {} ) 
#     #lvii_daph[[t]]$ri = ri
#     #Predicted values: 
#     s=seq(1,max(daph_tmp$N),1 )
#     #s=seq(0,maxT,1 )
#     if( is.null(lvii_daph[[t]]) ) {
#         f_lvii = formula ( Ndiff ~  ri*(1-aii*N) )
#         lvii_daph[[t]] =  nls( formula= f_lvii, data = daph_tmp, 
#         start=list(ri=2, aii=0.5 ),
#         control=nls.control(maxiter = 1000), algorithm="port", trace=F ) 
#     }

#     d_tmp = predict(lvii_daph[[t]], list( N = s ) )
#     lvii_daph_tmp = data.frame( species = rspecies[1], temperature = temps[t], s=s, N_pred = d_tmp )
#     lvii_daph_pred = rbind(lvii_daph_pred, lvii_daph_tmp)

#   #Dia
#   #lm_ii = lm(data=dia_tmp, Ndiff ~ N)
#   #ri=coef(lm_ii)[1]
#   #ri=1
#   tryCatch({ 
#     lvii_dia[[t]] = nls( formula= f_lvii, data = dia_tmp, 
#     start=list(ri=2, aii=0.5 ),
#     control=nls.control(maxiter = 1000), algorithm="port", trace=F ) 
#   }, error = function(e) {} ) 
#     #lvii_dia[[t]]$ri = ri
#     #Predicted values: 
#     s=seq(1,max(dia_tmp$N),1 )
#     if( is.null(lvii_dia[[t]]) ) {
#         f_lvii = formula ( Ndiff ~  ri*(1-aii*N) )
#         lvii_dia[[t]] =  nls( formula= f_lvii, data = dia_tmp, 
#         start=list(ri=2, aii=0.5 ),
#         control=nls.control(maxiter = 1000), algorithm="port", trace=F ) 
#     }

#     d_tmp = predict(lvii_dia[[t]], list( N = s ) )
#     lvii_dia_tmp = data.frame( species = rspecies[2], temperature = temps[t], s=s, N_pred = d_tmp )
#     lvii_dia_pred = rbind(lvii_dia_pred, lvii_dia_tmp)

# #=============================================================================
#   #Fit interspecific competition coefficients directly using the invasion
#   #scenarios.



#   f_lvij = formula ( Ndiff ~  ri*(1-aij*N_res) )
#   #f_lvij = formula ( Ndiff ~  ri/(1+aij*N_res) )
#   #f_lvij = formula ( Ndiff ~  ri/(1+aii*N_inv+aij*N_res) )
#   #f_lvij = formula ( Ndiff ~  ri*(1-aii*N_inv-aij*N_res) )


#   #Daphnia
#   #lm_ij = lm(data=daph_inv_tmp, Ndiff ~ N_res)
#   #ri=lvii_daph[[t]]$ri
#   #ri = subset(lvii_daph_pred,s == 1 & temperature == temps[t] )$N_pred[1]
#   #ri=1
#   aii = coef(lvii_daph[[t]])[2]
#   # ri = coef(lvii_daph[[t]])[1]
#   tryCatch({ 
#     lvij_daph[[t]] = nls( formula= f_lvij, data = daph_inv_tmp, 
#     start=list(ri=coef(lvii_daph[[t]])[1], aij=0.5 ),
#     control=nls.control(maxiter = 1000), algorithm="port", trace=F ) 
#   }, error = function(e) {} ) 
#     lvij_daph[[t]]$ri = ri
#     #Predicted values: 
#     s=seq(1,max(daph_inv_tmp$N_res,na.rm=T),1 )
#     #s=seq(0,maxT,1 )
#     # if( is.null(cR_dia[[t]]) ) {
#     #     cR_dia[[t]] = c1
#     #     d_tmp = exp(predict(cR_dia[[t]], list( algae_abundance= s ) ) )
#     # }else {
#     #     d_tmp = predict(cR_dia[[t]], list( algae_abundance= s ) )
#     # }
#     d_tmp = predict(lvij_daph[[t]], list( N_res = s ) )
#     lvij_daph_tmp = data.frame( species = rspecies[1], temperature = temps[t], s=s, N_pred = d_tmp )
#     lvij_daph_pred = rbind(lvij_daph_pred, lvij_daph_tmp)

#   #Dia

#   # c1 = lm(I(log(N+1))~I(day_n ), data = dia_tmp ) 
#   # tryCatch({ 
#   #   cR_dia[[t]] = nls( formula= cR, data = dia_tmp, 
#   #   start=list(c1=(as.numeric(coef(c1)[2])), a1=exp(as.numeric(coef(c1)[1] ) ) ),
#   #   control=nls.control(maxiter = 1000), algorithm="port", trace=F ) 

#   #   #Predicted values: 
#   #   s=seq(min(dia_tmp$day_n),max(dia_tmp$day_n),1 )
#   #   #s=seq(0,maxT,1 )
#   #   d_tmp = predict(cR_dia[[t]], list( day_n = s ) )
#   #   dia_pred_tmp = data.frame( species = rspecies[2], temperature = temps[t], s=s, N_pred = d_tmp )
#   #   cR_dia_pred = rbind(cR_dia_pred, dia_pred_tmp)
#   # }, error = function(e) {} ) 

#   #1. Lotka Volterra
#   #Use the growth rate from the intraspecific fits
#   #lm_ij = lm(data=dia_inv_tmp, Ndiff ~ N_res)
#   #ri=lvii_dia[[t]]$ri
#   #ri = subset(lvii_dia_pred,s == 1 & temperature == temps[t] )$N_pred[1]
#   #ri=1
#   aii = coef(lvii_dia[[t]])[2]
#   tryCatch({ 
#     lvij_dia[[t]] = nls( formula= f_lvij, data = dia_inv_tmp, 
#     start=list(ri =coef(lvii_dia[[t]])[1], aij=0.5 ),
#     control=nls.control(maxiter = 1000), algorithm="port", trace=F ) 
#   }, error = function(e) {} ) 
#     lvij_dia[[t]]$ri = ri
#     #Predicted values: 
#     s=seq(1,max(dia_inv_tmp$N_res,na.rm=T),1 )
#     #s=seq(0,maxT,1 )
#     # if( is.null(cR_dia[[t]]) ) {
#     #     cR_dia[[t]] = c1
#     #     d_tmp = exp(predict(cR_dia[[t]], list( algae_abundance= s ) ) )
#     # }else {
#     #     d_tmp = predict(cR_dia[[t]], list( algae_abundance= s ) )
#     # }
#     d_tmp = predict(lvij_dia[[t]], list( N_res = s ) )
#     lvij_dia_tmp = data.frame( species = rspecies[2], temperature = temps[t], s=s, N_pred = d_tmp )
#     lvij_dia_pred = rbind(lvij_dia_pred, lvij_dia_tmp)



}

#Consumption functions
cl_plot = rbind(cl_daph_plot,cl_dia_plot)
cl_pred = rbind(cl_daph_pred,cl_dia_pred)

#Growth rate functions
cR_pred = rbind(cR_daph_pred,cR_dia_pred)

#competition
lvii_pred = rbind(lvii_daph_pred, lvii_dia_pred)
lvij_pred = rbind(lvij_daph_pred, lvij_dia_pred)

aii_all = matrix(0,6,2)
aij_all = matrix(0,6,2)
ci_all = matrix(0,6,2)
wi_all = matrix(0,6,2)

for(n in 1:6){ 
  aii_all[n,1] = abs(coef(lvii_daph[[n]])[2])
  aii_all[n,2] = abs(coef(lvii_dia[[n]])[2])
  aij_all[n,1] = abs(coef(lvij_daph[[n]])[2])
  aij_all[n,2] = abs(coef(lvij_dia[[n]])[2])
  ci_all[n,1] = abs(coef(cl_daph[[n]])[1])
  ci_all[n,2] = abs(coef(cl_dia[[n]])[1])
  wi_all[n,1] = abs(coef(cR_daph[[n]])[1])
  wi_all[n,2] = abs(coef(cR_dia[[n]])[1])
}

aiib_all = ci_all^2*wi_all
aijb_all = ci_all*ci_all[,2:1]*wi_all

rho1 = matrix(0,6,2)
rho1b = matrix(0,6,2)
the1 = matrix(0,6,2)
the1b = matrix(0,6,2)

rho1[,1] = aij_all[,1]/(sqrt(aii_all[,1]*aii_all[,2]) )
rho1[,2] = aij_all[,2]/(sqrt(aii_all[,1]*aii_all[,2]))
rho1b[,1] = aijb_all[,1]/(sqrt(aiib_all[,1]*aiib_all[,2]))
rho1b[,2] = aijb_all[,2]/(sqrt(aiib_all[,1]*aiib_all[,2]))

the1[,1] = sqrt(aii_all[,1]/aii_all[,2])
the1[,2] = sqrt(aii_all[,2]/aii_all[,1])
the1b[,1] = sqrt(aiib_all[,1]/aiib_all[,2])
the1b[,2] = sqrt(aiib_all[,2]/aiib_all[,1])


#=============================================================================
#Consumption functions: 
#Separate panels
#Pick the appropriate first line: 
#ggplot(cl_plot, aes(x = N, y =algae_abundance, color = species) )+  #2
ggplot(cl_plot, aes(x = N, y =algae_abundance/algae_start, color = species) )+ #3
#ggplot(cl_plot, aes(x = N, y =Adiff, color = species) )+ #3
  geom_point( )+ facet_grid(temperature~species)+ 
  geom_line(data= cl_pred, mapping= aes(x = s, y =N_pred, color=species) ) + 
  facet_grid(temperature~species)+
  xlab("Zooplankton abundance ")+
  ylab("Algal consumption rate")+
  theme(strip.background = element_rect(colour=NA, fill=NA))
ggsave("./algal_consump3_diaDaph.pdf", width = 8, height = 10)

#By temp, species combined.
ggplot(cl_plot, aes(x = N, y =algae_abundance, color = species) ) + 
  geom_point( )+ facet_grid(temperature~.)+ 
  geom_line(data= cl_pred, mapping= aes(x = s, y =N_pred, color=species) )+
  facet_grid(temperature~.)+
  xlab("Zooplankton abundance ")+
  ylab("Algal consumption rate")+
  theme(strip.background = element_rect(colour=NA, fill=NA))
ggsave("./algal_consump2_samediaDaph.pdf", width = 8, height = 10)

#By species, temp combined.
ggplot(data= cl_pred, mapping= aes(x = s, y =N_pred, color=(temperature),group=(temperature) ) )+
  scale_fill_gradient(low = "yellow", high = "red", na.value = NA)+
  geom_line()+
  facet_grid(.~species)+
  xlab("Zooplankton abundance ")+
  ylab("Algal consumption rate")+
  theme(strip.background = element_rect(colour=NA, fill=NA))
ggsave("./algal_consump2_tempdiaDaph.pdf", width = 8, height = 10)

#=============================================================================
#Intrinsic growth rate functions: 
#Separate panels
#Pick the appropriate first line: 
#ggplot(cl_plot, aes(x = day_n, y =N, color = species) ) + #1. 
#ggplot(cl_plot, aes(x = Adiff, y =N, color = species) ) + #2. 
ggplot(cl_plot, aes(x = algae_abundance, y =N, color = species) ) + #2. 
  geom_point( )+ facet_grid(temperature~species)+ 
  geom_line(data= cR_pred, mapping= aes(x = s, y =N_pred, color=species) )+
  facet_grid(temperature~species)+ #xlim( min(cl_plot$Adiff), max(cl_plot$Adiff))+
  xlab("Zooplankton abundance ")+
  ylab("Time")+
  theme(strip.background = element_rect(colour=NA, fill=NA))
ggsave("./intrinsicR_diaDaph2.pdf", width = 8, height = 10)

#By temp, species combined.
ggplot(cl_plot, aes(x = day_n, y =N, color = species) ) + 
  geom_point( )+ facet_grid(temperature~.)+ 
  geom_line(data= cR_pred, mapping= aes(x = s, y =N_pred, color=species) )+
  facet_grid(temperature~.)+
  xlab("Zooplankton abundance ")+
  ylab("Time")+
  theme(strip.background = element_rect(colour=NA, fill=NA))
ggsave("./intrinsicR_samediaDaph.pdf", width = 8, height = 10)

#By species, temp combined.
ggplot(data= cR_pred, mapping= aes(x = s, y =N_pred, color=(temperature),group=(temperature) ) )+
  scale_fill_gradient(low = "yellow", high = "red", na.value = NA)+
  geom_line()+
  facet_grid(.~species)+
  xlab("Zooplankton abundance ")+
  ylab("Time")+
  theme(strip.background = element_rect(colour=NA, fill=NA))
ggsave("./intrinsicR_tempdiaDaph.pdf", width = 8, height = 10)

#=============================================================================
#Competition models 
#Separate panels
#Pick the appropriate first line: 
#ggplot(cl_plot, aes(x = day_n, y =N, color = species) ) + #1. 
#ggplot(cl_plot, aes(x = Adiff, y =N, color = species) ) + #2. 
ggplot(cl_plot, aes(x = N, y =Ndiff, color = species) ) + #2. 
  geom_point( )+ facet_grid(temperature~species) + 
  geom_line(data= lvii_pred, mapping= aes(x = s, y =N_pred, color=species) )+
  facet_grid(temperature~species)+ #xlim( min(cl_plot$Adiff), max(cl_plot$Adiff))+
  xlab("Zooplankton abundance ")+
  ylab("Growth rate")+  xlim(0,60)+
  theme(strip.background = element_rect(colour=NA, fill=NA))
#ggsave("./intrinsicR_diaDaph2.pdf", width = 8, height = 10)

ggplot(cl_plot, aes(x = N, y =Ndiff, color = species) ) + #2. 
  geom_point( )+ facet_grid(temperature~species) + 
  geom_line(data= lvij_pred, mapping= aes(x = s, y =N_pred, color=species) )+
  facet_grid(temperature~species)+ #xlim( min(cl_plot$Adiff), max(cl_plot$Adiff))+
  xlab("Zooplankton abundance ")+
  ylab("Growth rate")+  xlim(0,60)+
  theme(strip.background = element_rect(colour=NA, fill=NA))