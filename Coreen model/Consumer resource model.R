###########################MacArthur consumer resource schoolfield fits
library(tidyverse)



##Functions for temp dependencies from best fit models

sharpeschoolfull_modified_inverse <- function(temp, r_tref, e, el, tl, eh, th, tref){
  tref <- 273.15 + tref
  k <- 8.62e-05
  boltzmann.term <- r_tref*exp(e/k * (1/tref - 1/(temp + 273.15)))
  inactivation.term <- (1 + exp(-el/k * (1/(tl + 273.15) - 1/(temp + 273.15))) + exp(eh/k * (1/(th + 273.15) - 1/(temp + 273.15))))
  return(boltzmann.term * inactivation.term)
}


sharpeschoolhigh_1981 <- function(temp, r_tref, e, eh, th, tref){
  tref <- 273.15 + tref
  k <- 8.62e-05
  boltzmann.term <- r_tref*exp(e/k * (1/tref - 1/(temp + 273.15)))
  inactivation.term <- 1/(1 + exp(eh/k * (1/(th + 273.15) - 1/(temp + 273.15))))
  return(boltzmann.term * inactivation.term)
}




##### Read in the parameter estimates 

mortality.parameters <- read.csv("Coreen model/estimates.for.mortality.parameters.csv")


mortality_daphnia <- mortality.parameters%>%
  filter(species == "daphnia")%>%
  filter(parameter.type == "mortality")%>%
  filter(best.model == "yes")%>%
  column_to_rownames(var = 'term')

tref <- 20

mortality_diaphanosoma <- mortality.parameters%>%
  filter(species == "diaphanosoma")%>%
  filter(parameter.type == "mortality")%>%
  filter(best.model == "yes")%>%
  column_to_rownames(var = 'term')



grazing.parameters <- read.csv("Coreen model/estimates.for.grazing.parameters.csv")


###Daphnia consumption

grazing_daphnia <- grazing.parameters%>%
  filter(species == "daphnia")%>%
  filter(parameter.type == "grazing")%>%
  filter(best.model == "yes")%>%
  column_to_rownames(var = 'term')

###Diaphanosoma consumption

grazing_diaphanosoma <- grazing.parameters%>%
  filter(species == "diaphanosoma")%>%
  filter(parameter.type == "grazing")%>%
  filter(best.model == "yes")%>%
  column_to_rownames(var = 'term')


###Conversion factor: Assimilation efficiency of herbivores = 0.8 from Anderson 2004
###Here we will use carbon in one algae: carbon in one zooplankton
##See grazing rate for chapter 2 script for values
##dry weights
scenedesmus.wt = 48.51664^-6  #ug C per cell

#daphnia.wt =52.42142 #ug per average individual from grazing rate experiment

conversion.factor = 0.8
daphnia.wt= 8.39518
diaphanosoma.wt = 1.250867


daphna.conversion.efficiency = conversion.factor * scenedesmus.wt/daphnia.wt
diaphanosoma.conversion.efficiency = conversion.factor * scenedesmus.wt/ diaphanosoma.wt



temp_dependences_MacArthur <- function(temp = temp, ref_temp2 = 1, tref = 20, 
                                       #r_EaN = 0.5, 
                                       
                                       ###resource consumption parameters                        
                                       ##Diaphanosoma                        
                        r_tref_dia_consumption = grazing_diaphanosoma["r_tref", "estimate"],
                        e_dia_consumption = grazing_diaphanosoma["e", "estimate"],
                        #el_dia_consumption<- grazing_diaphanosoma["el", "estimate"]
                        #tl_dia_consumption<- grazing_diaphanosoma["tl", "estimate"]
                        eh_dia_consumption = grazing_diaphanosoma["eh", "estimate"],
                        th_dia_consumption = grazing_diaphanosoma["th", "estimate"],                        
                                       ##Daphnia
                        r_tref_daph_consumption = grazing_daphnia["r_tref", "estimate"],
                        e_daph_consumption = grazing_daphnia["e", "estimate"],
                        #el_daph_consumption<- grazing_daphnia["el", "estimate"]
                        #tl_daph_consumption<- grazing_daphnia["tl", "estimate"]
                        eh_daph_consumption = grazing_daphnia["eh", "estimate"],
                        th_daph_consumption = grazing_daphnia["th", "estimate"],
                                       
                                       ##Mortality rates (scale parameter)
                                       ##Daphnia
                        r_tref_daph_mort = mortality_daphnia["r_tref", "estimate"],
                        e_daph_mort = mortality_daphnia["e", "estimate"],
                        el_daph_mort = mortality_daphnia["el", "estimate"],
                        tl_daph_mort = mortality_daphnia["tl", "estimate"],
                        eh_daph_mort = mortality_daphnia["eh", "estimate"],
                        th_daph_mort = mortality_daphnia["th", "estimate"],
                                       
                        #Diaphanosoma
                        r_tref_dia_mort = mortality_diaphanosoma["r_tref", "estimate"],
                        e_dia_mort = mortality_diaphanosoma["e", "estimate"],
                        el_dia_mort = mortality_diaphanosoma["el", "estimate"],
                        tl_dia_mort = mortality_diaphanosoma["tl", "estimate"],
                        eh_dia_mort = mortality_diaphanosoma["eh", "estimate"],
                        th_dia_mort = mortality_diaphanosoma["th", "estimate"]
                                       
                                       
                          
                        #KR = 384.0671,            ## resource carrying capacity number of cells per day is 10000000 per mesocosm, ug carbon is 384.0671
                                       
                                      
){
  
  # resource supply rate- Should this have an r and K still?
  #KR = KR
  # resource growth rates
  #rR = arrhenius_function(temp = temp, E = 0, b1 = 1) ##Try making this large so resources are replenished?
  rR = 1
  # resource carrying capacity
  #rK = arrhenius_function(temp = temp, E = 0, b1 = 2) 
  ###Let's make this in terms of the number of cells per day since rsource consumption is in these units?
  rK = 10000000
  
  c_daph_R = sharpeschoolhigh_1981(temp = temp, r_tref = r_tref_daph_consumption, e = e_daph_consumption, eh = eh_daph_consumption, th = th_daph_consumption, tref = 20)
  
  
  c_dia_R = sharpeschoolhigh_1981(temp = temp, r_tref = r_tref_dia_consumption, e = e_dia_consumption, eh = eh_dia_consumption, th = th_dia_consumption, tref = 20)
  
  
  # vij = conversion factor that converts resource into biomass of consumer
  #v_daph_R = arrhenius_function(temp = temp, E = 0, b1 = 0.6) ####Dell assumes no effect oftemperature
  
  v_daph_R = daphna.conversion.efficiency
  v_dia_R = diaphanosoma.conversion.efficiency
  #v_dia_R = arrhenius_function(temp = temp, E = 0, b1 = 0.4)
  
  # mortality rates
  ##For units to work this must be proportion of population that dies per day I think
  
  m_daph = sharpeschoolfull_modified_inverse(temp = temp, r_tref = r_tref_daph_mort, e = e_daph_mort, el = el_daph_mort, tl = tl_daph_mort, eh = eh_daph_mort, th = th_daph_mort, tref = 20)
  
  m_dia = sharpeschoolfull_modified_inverse(temp = temp, r_tref = r_tref_dia_mort, e = e_dia_mort, el = el_dia_mort, tl = tl_dia_mort, eh = eh_dia_mort, th = th_dia_mort, tref = 20)
  
  
  g1 = v_daph_R * c_daph_R * rR * rK  - (m_daph) ### growth rate of daphnia
  g2 = v_dia_R * c_dia_R * rR * rK  - (m_daph) ### growth rate of diaphanosoma
  
  
  # Absolute competition coefficients
  beta11 = v_daph_R * c_daph_R * (rR*rK) * c_daph_R  ### intra
  beta12 = v_daph_R * c_daph_R * (rR*rK) * c_dia_R  ### inter
  beta22 = v_dia_R * c_dia_R * (rR*rK) * c_dia_R  ### intra
  beta21 = v_dia_R * c_dia_R * (rR*rK) * c_daph_R  ### inter
  
  ###Fixed betas
  ###These seem wrong
  
  # Relative competition coefficients 
  a11 = beta11 / g1
  a21 = beta21 / g2
  a22 = beta22 / g2
  a12 = beta12 / g1
  
  
  ###Invasion growth rates? See Song 2019
  
  r1 = g1*(1-(a12/a22))
  r2 = g2*(1-(a21/a11))
  
  # MCT parameters
  rho <- sqrt((a12*a21)/(a11*a22)) #niche overlap
  stabil_potential <- 1 - rho #stabilizing potential
  fit_ratio <- sqrt((a11*a12)/(a22*a21))  #fitness ratio
  coexist <- rho < fit_ratio &  fit_ratio < 1/rho
  
  # report results
  results <-  data.frame(temp = temp, 
                         a11 = a11, a12 = a12, a22 = a22, a21 = a21, g1 = g1, g2 = g2, r1 = r1, r2 = r2, beta12 = beta12, beta22= beta22, beta21 = beta21, beta11 = beta11,
                         stabil_potential = stabil_potential, fit_ratio = fit_ratio, rho = rho, coexist = coexist, 
                         m_daph = m_daph, m_dia = m_dia, rR = rR, rK = rK, 
                         c_daph_R = c_daph_R, c_dia_R = c_dia_R)
  
  return(results)
}




temp <- seq(4, 36, by = 0.1)


results <- as.data.frame(temp_dependences_MacArthur(temp))






##Plot some results
#Yikes

ggplot(results)+
  geom_line(aes(x = temp, y = g1*.1), colour = "blue", data = results)+
  geom_line(aes(x = temp, y = g2*.1), colour = "red", data = results)


ggplot(results)+
  geom_line(aes(x = temp, y = r1), colour = "blue", data = results)+
  geom_line(aes(x = temp, y = r2), colour = "red", data = results)+
  ylim(-1,3.5)



###Competition coefficients-muy yikes
ggplot(results)+
  geom_line(aes(x = temp, y = a11), colour = "blue", data = results)+
  geom_line(aes(x = temp, y = a22), colour = "red", data = results)+
  geom_line(aes(x = temp, y = a21), colour = "orange", data = results)+
  geom_line(aes(x = temp, y = a12), colour = "purple", data = results)+
  ylim(-2, 30000)

ggplot(results)+
  geom_line(aes(x = temp, y = beta11), colour = "blue", data = results)+
  geom_line(aes(x = temp, y = beta22), colour = "red", data = results)+
  geom_line(aes(x = temp, y = beta21), colour = "orange", data = results)+
  geom_line(aes(x = temp, y = beta12), colour = "purple", data = results)+
  ylim(-2, 30000)


ggplot(aes(x = temp, y = stabil_potential), data = results) + 
  geom_point()



ggplot(aes(x = temp, y = fit_ratio), data = results) + 
  geom_point()




