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
# get_mod_fit ()
# Function fitting
# This is a wrapper function for fitting NLS models, predicting values based
# on the fitted model, and error controlling when the models fail. If an NLS 
# fit fails with default conditions a second algorithm is attempted, and if this 
# also fails then a linear model is attempted. 
# mod_data            The input data set 
# mod_fit             The formula of the model to fit 
# mod_prms            A list of model parameters 
# prm_start           The starting values of parameters to be fit
# mod_x               The X value used in the model fitting, for prediction
# lm_mod              If supplied, use this linear model structure if NLS fails
# fixed               Use this to treat a parameter in the supplied mod_fit as 
#                     a fixed parameter. 
#=============================================================================

get_mod_fit = function ( mod_data, mod_fit, mod_prms, prm_start, mod_x, 
  lm_mod = NULL, fixed = NULL) {


  #Check to see if any fixed parameters have been passed. If so, assign them. 
  if( !is.null(fixed)) { 
    nfixed = length(fixed)
    for ( n in 1:nfixed) {
      assign(paste(names(fixed)[n]), unlist(c(fixed[n],recursive =T ) ) ) 
    }
    #Change the length of the parameter vector
    prm_start = prm_start[-(mod_prms == names(fixed))]
    mod_prms = mod_prms[-(mod_prms == names(fixed))] 

  }

  #Initial parameters
  mpl = length ( mod_prms)
  start1 = vector("list",mpl)

  #Make the list of initial model values. 
  for ( n in 1:mpl ){ 
    start1[[n]] = prm_start[n]
    names(start1)[n] = mod_prms[n]
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
  y_tmp =seq(1, max(mod_data[[paste(mod_x)]],na.rm=T),1 )

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


#=============================================================================
# get_info_ts() 
# Data processing to prepare data from the real empirical data to be passed
# to the info theory functions (in info_theory_functions.R). 
#
#
#=============================================================================

get_info_ts = function (m1_data_long, mesosi ){ 
  #=============================================================================
  #Make new resident and invader data sets.
  #=============================================================================
  #Make a new resident data set to fit the growth rate function with nlme/nls 
  #The new data set includes a column for n(t+1)/n(t) and a columnf for the time 
  #interval across subsequent measurements.  
  m1_tmp = subset(m1_data_long, mesocosm_id == mesosi)
  sp_tmp = substr(as.character(unique(mesosi)),1,3) #Resident
  r_sp = rspecies[grepl(sp_tmp,rspecies,fixed=T)]

  res_tmp = subset(m1_tmp, species == r_sp  ) %>% arrange(day_n)
  res_tmp = res_tmp[!is.na(res_tmp$N),] #NLME won't work with NAs 
  #Arrange the data by replicate number, then add a new column for the delta N
  res_tmp = res_tmp %>% 
      arrange(replicate_number)%>%
      mutate(Ndiff_res = lead( N-lag(N),)/lead(day_n-lag(day_n))*1/N) #"lead" lines up the result 
     
  #Add a new column for the change in algal biomass
  res_tmp = res_tmp %>%
      mutate(Ndiff_alg = (lag(algae_cells_mL)- algae_abundance ) )  #"lead" lines up the result

  res_tmp$Ndiff_res[res_tmp$day_n == max(res_tmp$day_n,na.rm=T )] =NA #Remove last day
  res_tmp = res_tmp %>% mutate(tdiff =day_n-lag(day_n)) #Size of time step
  res_tmp$tdiff[res_tmp$tdiff<0] = 1 #Remove negative time steplag(N)s
  res_tmp = res_tmp %>% mutate(N_res = N)
  colnames(res_tmp)[colnames(res_tmp)=="zooplankton_length_cm"] = "zooplankton_length_cm_res"

  #Make a new invader data set to fit the growth rate function with nlme/nls 
  #The new data set includes a column for n(t+1)/n(t) and a columnf for the time 
  #interval across subsequent measurements.  
  i_sp = rspecies[!grepl(sp_tmp,rspecies,fixed=T)]
  inv_tmp = subset(m1_tmp, species == i_sp )%>% arrange(day_n)
  inv_tmp = inv_tmp[!is.na(inv_tmp$N),]
  #inv_tmp$N[is.na(inv_tmp$N)] = 0 #Replace NAs with 0  
  #Arrange the data by replicate number, then add a new column for the delta N
  inv_tmp = inv_tmp %>% 
      arrange(replicate_number)%>%
      mutate(Ndiff_inv = lead( N-lag(N),)/lead(day_n-lag(day_n))*1/N) #"lead" lines up the result 
  
  #Add a new column for the change in algal biomass
  inv_tmp = inv_tmp %>%
      mutate(Ndiff_alg = (lag(algae_cells_mL)- algae_abundance) )  #"lead" lines up the result

  inv_tmp$Ndiff_inv[inv_tmp$day_n == max(inv_tmp$day_n,na.rm=T )] =NA #Remove last day
  inv_tmp = inv_tmp %>% mutate(tdiff =day_n-lag(day_n)) #Size of time step 
  inv_tmp = inv_tmp %>% mutate(N_inv = N)
  colnames(inv_tmp)[colnames(inv_tmp)=="zooplankton_length_cm"] = "zooplankton_length_cm_inv"

  #This line adds the resident densities on the matching days from the matching 
  #replicates to the data table inv_tmp
  inv_tmp = inv_tmp %>% left_join( 
    select( res_tmp, replicate_number, day_n, zooplankton_length_cm_res, N_res, Ndiff_res),
    by= (c( "day_n" = "day_n", "replicate_number"="replicate_number" )) )
  res_tmp = res_tmp %>% left_join( 
     select( inv_tmp, replicate_number,day_n, zooplankton_length_cm_inv, N_inv, Ndiff_inv),
     by= (c( "day_n" = "day_n", "replicate_number"="replicate_number" )) )
  inv_tmp = inv_tmp[!is.na(inv_tmp$N_inv), ] #Omit NAs
  res_tmp = res_tmp[!is.na(res_tmp$N_res) , ]

  #Create a dummy data set in the few cases that experiment was essentially 
  #unsuccesful: 
  if(dim(inv_tmp)[1] <=0  ){
    
    #New population, diff,algal, and body size data
    #mydata_res = subset(res_tmp, day_n >= inv_day  ) #Only days after invasion
    mydata_res = res_tmp
    mydata_res$N_inv = 0
    out1[[i]] = mydata_res[!is.na(mydata_res$N_res), ]

    #Take the time series and interpolate missing days:
    new_days = data.frame(day_n = seq(min(out1[[i]]$day_n),max(out1[[i]]$day_n,1) ),
             species = out1[[i]]$species[1], replicate_number = out1[[i]]$replicate_number[1],  
             invade_monoculture = out1[[i]]$invade_monoculture[1], 
             temperature= out1[[i]]$temperature[1], mesocosm_id =out1[[i]]$mesocosm_id[1] )

    out_tmp = new_days %>% left_join(out1[[i]][,7:24])
    #Try an NA interpolation on every column (this doesn't always work, hence the 
    #obnoxious tryCatch for loop)  
    for (j in 7:23){  tryCatch({ out_tmp[,j] = na.approx(out_tmp[,j])}, error = function(e) {} )}
    return(out_tmp)

  } else if (dim(res_tmp)[1] <=0  ){
    #New population, diff,algal, and body size data
    #mydata_inv = subset(inv_tmp, day_n >= inv_day  )#Only days after invasion
    mydata_inv = inv_tmp
    mydata_inv$N_res = 0
    out1[[i]] = mydata_inv[!is.na(mydata_res$N_inv), ]

    #Take the time series and interpolate missing days:
    new_days = data.frame(day_n = seq(min(out1[[i]]$day_n),max(out1[[i]]$day_n,1) ),
             species = out1[[i]]$species[1], replicate_number = out1[[i]]$replicate_number[1],  
             invade_monoculture = out1[[i]]$invade_monoculture[1], 
             temperature= out1[[i]]$temperature[1], mesocosm_id =out1[[i]]$mesocosm_id[1] )

    out_tmp = new_days %>% left_join(out1[[i]][,7:24])
    #Try an NA interpolation on every column (this doesn't always work, hence the 
    #obnoxious tryCatch for loop)  
    for (j in 7:23){  tryCatch({ out_tmp[,j] = na.approx(out_tmp[,j])}, error = function(e) {} )}
    return(out_tmp)

  } else {
    #New population, diff,algal, and body size data
    #mydata_inv = subset(inv_tmp, day_n >= inv_day  )#Only days after invasion
    mydata_inv= inv_tmp
    out1[[i]] = mydata_inv[!is.na(mydata_inv$N_inv) & 
      !is.na(mydata_inv$N_res) , ]
    
    #Take the time series and interpolate missing days:
    new_days = data.frame(day_n = seq(min(out1[[i]]$day_n),max(out1[[i]]$day_n,1) ),
             species = out1[[i]]$species[1], replicate_number = out1[[i]]$replicate_number[1],  
             invade_monoculture = out1[[i]]$invade_monoculture[1], 
             temperature= out1[[i]]$temperature[1], mesocosm_id =out1[[i]]$mesocosm_id[1] )

    out_tmp = new_days %>% left_join(out1[[i]][,7:24])
    #Try an NA interpolation on every column (this doesn't always work, hence the 
    #obnoxious tryCatch for loop)  
    for (j in 7:23){  tryCatch({ out_tmp[,j] = na.approx(out_tmp[,j])}, error = function(e) {} )}
    return(out_tmp)

  } 

}
