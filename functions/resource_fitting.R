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
library(fields)
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

#=============================================================================
#get_new_mono()
# Build the needed sub data sets, which are broken down by temperature and 
# species and include a new column for the growth rate. This function is only
# meaningful for this particular data set. 
# Adds Ndiff, which is dN/dt * 1/N
# Adds Adiff, which is starting algae - ending algae. 
#=============================================================================
get_new_mono = function (m1_data_long, species1, temperature,inv_day ) {
    all_species = unique(m1_data_long$species)
    all_invasions = unique(m1_data_long$invade_monoculture)

    if(species1 == "daphnia") {invasion = all_invasions[3]} else {invasion = all_invasions[2]}

    mesos1 = subset(m1_data_long,species == species1 & invade_monoculture == "monoculture" )
    mesos2 = subset(m1_data_long,species == species1 & invade_monoculture == invasion )
    mesos2 = mesos2[ mesos2$day_n<=inv_day, ]
    mesos_species = rbind(mesos1,mesos2)
    species_tmp = subset(mesos_species, temperature == temps[t])

    #Add the rate of population growth, Ndiff
    dp1 = subset(species_tmp,!is.na(species_tmp$N) ) %>%    
            arrange(replicate_number) %>%
            mutate(Ndiff = (N-lag(N))/(day_n-lag(day_n))*1/N )  #"lead" lines up the result
    dp1$Ndiff[is.infinite(dp1$Ndiff)] = NA
    species_tmp = species_tmp %>% left_join(dp1)

    #Add the rate of algal consumption, Adiff
    species_tmp = species_tmp %>%    
            arrange(replicate_number) %>%
            mutate(Adiff = algae_start-algae_abundance)  #"lead" lines up the result 
    species_tmp$Adiff[is.infinite(species_tmp$Adiff)] = NA

    return(species_tmp)

}

#=============================================================================
#get_new_inv()
# Build the needed sub data sets, which are broken down by temperature and 
# species and include a new column for the growth rate. This function is only
# meaningful for this particular data set. 
# Adds Ndiff, which is dN/dt * 1/N
# Adds Adiff, which is starting algae - ending algae. 
#=============================================================================
get_new_inv = function (m1_data_long, species1, temperature,inv_day,inv_end ) {
  
  all_species = unique(m1_data_long$species)
  all_invasions = unique(m1_data_long$invade_monoculture)
  species2 = all_species[all_species != species1 ] #The resident 

  if(species1 == "daphnia") {invasion = all_invasions[2]} else {invasion = all_invasions[3]}
  mesos2 = subset(m1_data_long, invade_monoculture == invasion )
  mesos2 =  mesos2[ mesos2$day_n>=inv_day & mesos2$day_n <= inv_end , ]
  species_tmp = subset(mesos2, temperature == temps[t])

  #Creates specific columns designating invader and resident: 
  dp2i = subset(species_tmp, species == species1)
  dp2r = subset(species_tmp, species == species2)
  colnames(dp2i)[colnames(dp2i)=="N"] = "N_inv"
  colnames(dp2r)[colnames(dp2r)=="N"] = "N_res"
  species_inv_tmp = dp2i %>% 
    left_join(dp2r[,c("day_n","mesocosm_id","N_res")], by=c("day_n","mesocosm_id"))

  #Add the rate of population growth, Ndiff
  dp1 = subset(species_inv_tmp,!is.na(species_inv_tmp$N_inv) ) %>%    
            arrange(replicate_number) %>%
            mutate(Ndiff = (N_inv-lag(N_inv))/(day_n-lag(day_n))*1/N_inv) #"lead" lines up the result 
  dp1$Ndiff[is.infinite(dp1$Ndiff)] = NA
  species_inv_tmp = species_inv_tmp %>% left_join(dp1)
  

  return(species_inv_tmp)

}

#=============================================================================
#Function fitting
#=============================================================================


get_mod_fit = function ( mod_data, mod_fit, mod_prms, prm_start, mod_x, 
  lm_mod = NULL, fixed = NULL) {

  mpl = length ( mod_prms)
  start1 = vector("list",mpl)

  #Make the list of initial model values. 
  for ( n in 1:mpl ){ 
    start1[[n]] = prm_start[n]
    names(start1)[n] = mod_prms[n]
  }

  #Check to see if any fixed parameters have been passed. If so, assign them. 
  if( !is.null(fixed)) { 
    nfixed = length(fixed)
    for ( n in 1:nfixed) {
      assign(paste(names(fixed)[n]), unlist(c(fixed[n],recursive =T ) ) ) 
    }
  }

  #Fit the model with NLS
  fit_mod = NULL
  tryCatch({ 
    fit_mod = nls( formula= mod_fit, data = mod_data, 
    start=start1, control=nls.control(maxiter = 1000), trace=F ) 
  }, error = function(e) {} ) 

  #If NLS failed, try it once more with a different algorithm
  if(is.null(fit_mod) ){ 
    tryCatch({ 
    fit_mod = nls( formula= mod_fit, data = mod_data, 
      start=start1, control=nls.control(maxiter = 1000),algorithm="port", trace=F ) 
    }, error = function(e) {} ) 
  }

  #Fitting data:
  #Sequence of data to predict over
  y_tmp =seq(min(mod_data[[paste(mod_x)]],na.rm=T),
    max(mod_data[[paste(mod_x)]],na.rm=T),1 )

  fit_dat = list(y_tmp)
  names(fit_dat)[1] = mod_x

  #See whether both NLS attempts failed. If lm_mod has been provided then attempt 
  #a linear model fit. Otherwise, do not try to fit the data and return with
  #nothing. 

  if(is.null(fit_mod) ){ 
    if( !is.null(lm_mod) ) {
      #Fit the linear model
      fit_mod = lm(lm_mod, data= mod_data )
      #Predicted values: 
      new_fit = predict(fit_mod, fit_dat )
      if( grep("log", lm_mod)>0) {   
        new_fit = exp(new_fit)
      }
      new_fit = data.frame( cbind(fit_dat[[1]],new_fit) )
      names(new_fit)[1] = paste(mod_x)
      names(new_fit)[2] = "N_pred"
      fit_mod = list(fit_mod = fit_mod, new_fit=new_fit)
    } else { 
      #Fail
      return(fit_mod)
    }

  }else{  
    #Predicted values: 
    new_fit = predict(fit_mod, fit_dat )
    new_fit = data.frame( cbind(fit_dat[[1]],new_fit) )
    names(new_fit)[1] = paste(mod_x)
    names(new_fit)[2] = "N_pred"
    fit_mod = list(fit_mod = fit_mod, new_fit=new_fit)
    return(fit_mod)
  }
 
}