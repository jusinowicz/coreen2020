#=============================================================================
#Information theory and population dynamics 
#=============================================================================
# This is R code to take population time series of two competing populations
# from a controlled (lab) experiment. The data are time series of Daphnia and
# Diaphanosoma species at 5 different Temperature treatment levels. The 
# different mesocosms also correspond with either an invasion or single species
# resident treatment. 
#
# This analysis measures several specific dynamic information theoretic 
# features of the time series: 
#   1. Excess entropy -- interpreted as the memory required/stored to produce
#      a particular output, and also equated to the complexity of the process.
#   2. Transient information -- information that is transferred between 
#      different sub-processes (i.e. population time series) of the system. 
#=============================================================================

#=============================================================================
#Load libraries
#=============================================================================
library(tidyverse)
library(lubridate)
library(mgcv)
library(gamm4)
library(fields)
source("./functions/info_theory_functions.R")
#=============================================================================
#Load data
#=============================================================================
#Coreen's experiment data
#m1=read.csv(file= "mesocosm_experiment_coreen_forbes.csv") 
#m1_inv = m1[1:72,] #Just take the invasion experiments

#This bit is now depricated with the new data sheet 
# m1_data_long = m1 %>% 
#   gather(key = day, value = N, day1:day83) %>% 
#   dplyr::select(-X) %>% 
#   mutate(day_n = as.numeric(gsub(x = day, pattern = "[^0-9.-]", replacement = "")))

# m1_data_long$replicate_number= factor(m1_data_long$replicate_number)

#Preprocessing m1_...
m1=read.csv(file= "tequila_master_data.csv") 
m1$date = ymd(m1$date) #Fix date

#There are certain days with data across individuals from a 
#mesocosm. This is to take all of those measurements and combine
#them into one measurement per day per mesocosm. This will
#use the average: 
cols2group=colnames(m1)[c(2:10,12,15,17:18)] 
cols2summarize = colnames(m1)[c(14,16)]
treats = unique(m1$mesocosm_id)
m1_data_long = m1 %>% 
     group_by_at(cols2group) %>%
     summarise_at(cols2summarize, .funs= mean)
m1_data_long= ungroup(m1_data_long)
m1_data_long$algae_abundance = as.numeric(m1_data_long$algae_abundance)


algae1 = read.csv(file= "algae_export.csv")[,1:3]
#Fix dates in algae:
algae1$date =  parse_date_time(algae1$date, "md")
year(algae1$date) [year(algae1$date) == 0000] = 2019
year(algae1$date) [year(algae1$date) == 2020 & 
    month(algae1$date) > 9  & month(algae1$date) <=12  ] = 2019

algae1$date = ymd (algae1$date)
algae1$algae_cells_mL[is.na(algae1$algae_cells_mL)] = mean(algae1$algae_cells_mL)
algae1$algae_mL_media[is.na(algae1$algae_mL_media)] = mean(algae1$algae_mL_media)

#Join algae to m1_...
m1_data_long= m1_data_long %>% 
  left_join( algae1) 

colnames(m1_data_long)[colnames(m1_data_long) == "zooplankton_abundance"] = "N"

#write.csv (file="tequila_master_data_algae.csv", m1_data_long)

#=============================================================================
#Fit a GAMM to the data to predict algal abundance on the basis of density per-species, 
#as a funciton of temperature, with mesocosm_id as a Random Effect
#=============================================================================
#Note that the data are fir with the log link function (log transformed) because 
#values cannot be <0. 
alg_gam = gam( algae_abundance ~ s(N)+s(temperature,k=5)+s(species, bs="re")+s(mesocosm_id,bs="re"),
  family=Gamma(link='log'), data=m1_data_long )

#Use the fitted GAMM to interpolate NAs for algal abundance in the data set
alg_newdata = subset(m1_data_long, is.na(algae_abundance) )
alg_replace = as.vector(predict.gam(alg_gam, newdata=alg_newdata,type= "response") )
alg_newdata$algae_abundance = alg_replace
m1_data_long$algae_abundance[is.na(m1_data_long$algae_abundance)] = alg_newdata$algae_abundance

m1_data_long %>% 
ggplot(aes(y = algae_abundance, x = N, color = species, group = interaction(species,replicate_number)))+
  geom_point()
ggsave("./algae_projection1.pdf", width = 8, height = 10)

#=============================================================================
#Plot the data
#=============================================================================
m1_data_long %>% 
  ungroup() %>% 
  filter(!is.na(N)) %>% 
  ggplot(aes(x = day_n, y = N, color = species, group = interaction(species,replicate_number)))+
  geom_line()+
  facet_grid(temperature~invade_monoculture)+
  scale_y_log10()+
  theme_bw()+
  scale_color_manual(values = c("dodgerblue1", "red1"), name = NA)+
  xlab("day")+
  ylab("population size")+
  theme(strip.background = element_rect(colour=NA, fill=NA))
ggsave("./Coreen_plot.pdf", width = 8, height = 10)
Collapse

#=============================================================================
#Collect basic data on number of treatments, species, etc.
#=============================================================================
nspp=2
temps = (unique(m1_data_long$temperature))
ntreatments =  length(unique(m1_data_long$temperature))
mesos = unique(m1_data_long$mesocosm_id)
nmesos =  length(mesos)
rspecies = unique(m1_data_long$species)
invader = unique(m1_data_long$invade_monoculture)
invasions_per = length(treats)/4 #Treatment entries correponding to each spp 
inv_day = 28 #The first day of attempted invasion
no_reps = 18 #The number of replicated mesocosms total per resident/invader

#=============================================================================
# Variables for data collection
#=============================================================================
out1 = vector("list",nmesos)
out1_diff = vector("list",nmesos)
out1_alg = vector("list",nmesos)
out1_bs = vector("list",nmesos)

#Dynamic information metrics calculated from the (discretized) time series 
di_web =  vector("list",nmesos)
#Track the average transfer entropy and separable information between each pair of 
#species as a way to build a network of information flow through the network. 
te_web =  vector("list",nmesos)
si_web = vector("list",nmesos) 

#The ensemble version of the AI
aiE_web = vector("list",nmesos)

#=============================================================================
# First, loop over both diaphanosoma and daphnia when diaphanasoma is the 
# invader
#=============================================================================

#Loop over all mesocosms
for (i in 1:nmesos) {

  index1=i
  #=============================================================================
  #Make new resident and invader data sets.
  #=============================================================================
  #Make a new resident data set to fit the growth rate function with nlme/nls 
  #The new data set includes a column for n(t+1)/n(t) and a columnf for the time 
  #interval across subsequent measurements.  
  m1_tmp = subset(m1_data_long, mesocosm_id == mesos[i])
  sp_tmp = substr(as.character(unique(mesos[i])),1,3) #Resident
  r_sp = rspecies[grepl(sp_tmp,rspecies,fixed=T)]

  res_tmp = subset(m1_tmp, species == r_sp  ) %>% arrange(day_n)
  res_tmp = res_tmp[!is.na(res_tmp$N),] #NLME won't work with NAs 
  #Arrange the data by replicate number, then add a new column for the delta N
  res_tmp = res_tmp %>% 
      arrange(replicate_number)%>%
      mutate(Ndiff_res = lead( N-lag(N),)/lead(day_n-lag(day_n))) #"lead" lines up the result 
     
  #Add a new column for the change in algal biomass
  res_tmp = res_tmp %>%
      mutate(Ndiff_alg = (lag(algae_cells_mL)- algae_abundance ) )  #"lead" lines up the result

  res_tmp$Ndiff_res[res_tmp$day_n == max(res_tmp$day_n )] =NA #Remove last day
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
      mutate(Ndiff_inv = lead( N-lag(N),)/lead(day_n-lag(day_n))) #"lead" lines up the result 
  
  #Add a new column for the change in algal biomass
  inv_tmp = inv_tmp %>%
      mutate(Ndiff_alg = (lag(algae_cells_mL)- algae_abundance) )  #"lead" lines up the result

  inv_tmp$Ndiff_inv[inv_tmp$day_n == max(inv_tmp$day_n )] =NA #Remove last day
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
    out1[[index1]] = mydata_res[!is.na(mydata_res$N_res), ]
  

  } else if (dim(res_tmp)[1] <=0  ){
    #New population, diff,algal, and body size data
    #mydata_inv = subset(inv_tmp, day_n >= inv_day  )#Only days after invasion
    mydata_inv = inv_tmp
    mydata_inv$N_res = 0
    out1[[index1]] = mydata_inv[!is.na(mydata_res$N_inv), ]

  } else {
    #New population, diff,algal, and body size data
    #mydata_inv = subset(inv_tmp, day_n >= inv_day  )#Only days after invasion
    mydata_inv= inv_tmp
    out1[[index1]] = mydata_inv[!is.na(mydata_inv$N_inv) & 
      !is.na(mydata_inv$N_res) , ]
  } 
  
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
  pop_ts = floor(f1*out1[[index1]][c("N_res","N_inv")])
  
  nt1 = 1
  nt2 = dim(pop_ts)[1]
  if(nt2 <=k){ k = 1}

  if(nt1 != nt2) { 

    f1 = 1 #scaling factor
    di_web[index1] = list(get_info_dynamics(pop_ts = pop_ts , k=k,with_blocks=TRUE))

    ## This code takes the population time-series counts output by the ODEs and 
    ## calculates the average Transfer Entropy from each species to every other 
    ## species. The goal is to get an overview of the major information pathways 
    ## in the web.   
    #=============================================================================
    # This function gives:
    # te_web    Average transfer entropy per species as a pairwise matrix
    #=============================================================================
    te_web[index1] = list( get_te_web( pop_ts = pop_ts, k=k) )

    #=============================================================================
    # This function gives:
    # aiE_web    The AI of the entire ensemble, treated as a single time series. 
    #=============================================================================
    aiE_web[index1] = list( get_ais (  series1 = pop_ts, k=k, ensemble = TRUE)    )

    #=============================================================================  
    #Build these back out into a data frame that includes all of the mesocosm,
    #treatment, and species information. 

    #Add the DIT to the data frames. There will be 11 new columns. 
    DIT_tmp = matrix(0,nt2,11)
    ncnames = c("N_res","N_inv","aiE","te1","te2","ee1","ee2","ai1","ai2","si1","si2")
    DIT_tmp[,1:2] = as.matrix(pop_ts)
    DIT_tmp[(k+1):nt2,3] = aiE_web[[index1]]$local
    DIT_tmp[(k+1):nt2,4:5] = di_web[[index1]]$te_local
    DIT_tmp[(k*2):nt2,6:7] = di_web[[index1]]$ee_local
    DIT_tmp[(k+1):nt2,8:9] = di_web[[index1]]$ai_local
    DIT_tmp[(k+1):nt2,10:11] = di_web[[index1]]$si_local
    colnames(DIT_tmp) = ncnames
    DIT_tmp = as.data.frame(DIT_tmp) 
    out1[[index1]] = out1[[index1]] %>% left_join(DIT_tmp)

  }

 


}

#Take all of the new data and combine it into one data frame: 
m1_DIT = bind_rows(out1, .id = "column_label")
m1_DIT=m1_DIT %>%
 mutate(alg_per_N = algae_abundance/N)  #"lead" lines up the result
m1_DIT=m1_DIT %>%
 mutate(alg_per_Ndiff = algae_abundance/lead( N-lag(N),)) 
m1_DIT$alg_per_N[is.infinite(m1_DIT$alg_per_N)] = NA
m1_DIT$alg_per_Ndiff[is.infinite(m1_DIT$alg_per_Ndiff)] = NA


#=============================================================================
#Fit GAMMs to the DIT results, grouping everything into one model with 
#species and mesocosm as random effects.Fit fit the hump-shaped 
#relationships between the DIT and the rate of algal consumption per individual. 
#The next plot is the better one 
#=============================================================================
N_tmp = 0:max(m1_DIT$N)
aiE_tmp = seq(0,max(m1_DIT$aiE),by= 0.1)
DIT_newdata = crossing(species=rspecies, temperature = temps, mesocosm_id= 1,aiE=aiE_tmp)

#Fit the GAMM
aiE_gam = gam( alg_per_N ~ s(aiE,k=3)+s(temperature,k=5)+s(species, bs="re")+s(mesocosm_id,bs="re"),
  data=m1_DIT )

#Use the fitted GAMM to plot the resulting relationships
N_plot = as.vector(predict.gam(aiE_gam, newdata=DIT_newdata, exclude = "s(mescosom_id)", type= "response") )
aiEp =  cbind(DIT_newdata, N_plot)
aiEp %>% 
  #filter(temperature == temps[1])%>%
  ggplot(aes(x = aiE, y =N_plot, color = species))+
  geom_line()+
  geom_point(data= m1_DIT, mapping= aes(x = aiE, y =alg_per_N, color = species) ) +
  facet_grid(temperature~species)

#=============================================================================
#Fit GAMMs to each temperature treatment X species separately to account for 
#differing curvature at each treatment level. Fit fit the hump-shaped 
#relationships between the DIT and the rate of algal consumption per individual. 
#============================================================================= 
#Loop through temperature treatments
m1_DIT_daph = vector("list", 6)
m1_DIT_dia = vector("list", 6)

m1_daph_plot = NULL
m1_dia_plot = NULL

for(t in 1:6) { 
  

  #Pull out each species for each temperature
  #########
  #Daphnia
  m1_DIT_daph_tmp = subset (m1_DIT, temperature == temps[t] & species == rspecies[1] )
  
  #Fit the GAMM
  aiE_daph_gam = gam( alg_per_N~ s(aiE,k=3)+s(mesocosm_id,bs="re"), data=m1_DIT_daph_tmp  )
  
  m1_DIT_daph [[t]] = aiE_daph_gam

  #Create the dummy data set for plotting
  N_tmp = 0:max(m1_DIT_daph_tmp$alg_per_N,na.rm=T)
  aiE_tmp = seq(0,max(m1_DIT_daph_tmp$aiE,na.rm=T),by= 0.1)
  DIT_daph_newdata = crossing(species=rspecies[1], temperature = temps[t], 
      mesocosm_id= 1,aiE=aiE_tmp)
  DIT_daph_newdata = subset(DIT_daph_newdata, species ==rspecies[1] )

  #Use the fitted GAMM to plot the resulting relationships and save the data: 
  N_plot = as.vector(predict.gam(aiE_daph_gam, newdata=DIT_daph_newdata, exclude = "s(mescosom_id)", 
    type= "response") )
  aiEp_daph =  cbind(DIT_daph_newdata, N_plot)
  m1_daph_plot = rbind( m1_daph_plot , aiEp_daph)

  ggplot(aiEp_daph, aes(x = aiE, y =N_plot, color = species) ) + 
  geom_line( )+  
  geom_point(data= m1_DIT_daph_tmp, mapping= aes(x = aiE, y =alg_per_N, color = species) )


  #########
  #Diaphanasoma
  m1_DIT_dia_tmp = subset (m1_DIT, temperature == temps[t] & species == rspecies[2] )

  #Fit the GAMM
  tryCatch( {aiE_dia_gam = gam( alg_per_N~ s(aiE,k=3)+s(mesocosm_id,bs="re"), data=m1_DIT_dia_tmp  )
  m1_DIT_dia [[t]] = aiE_dia_gam}, error = function(e){})

  #Create the dummy data set for plotting
  N_tmp = 0:max(m1_DIT_dia_tmp$alg_per_N,na.rm=T)
  aiE_tmp = seq(0,max(m1_DIT_dia_tmp$aiE,na.rm=T),by= 0.1)
  DIT_dia_newdata = crossing(species=rspecies[2], temperature = temps[t], 
      mesocosm_id= 1,aiE=aiE_tmp)
  DIT_dia_newdata = subset(DIT_dia_newdata, species ==rspecies[2] )

  #Use the fitted GAMM to plot the resulting relationships and save the data: 
  N_plot = as.vector(predict.gam(aiE_dia_gam, newdata=DIT_dia_newdata, exclude = "s(mescosom_id)", 
    type= "response") )
  aiEp_dia =  cbind(DIT_dia_newdata, N_plot)
  m1_dia_plot = rbind( m1_dia_plot , aiEp_dia)

  ggplot(aiEp_dia, aes(x = aiE, y =N_plot, color = species) ) + 
  geom_line( )+  
  geom_point(data= m1_DIT_dia_tmp, mapping= aes(x = aiE, y =alg_per_N, color = species) )


}

m1_DIT_plot = rbind(m1_daph_plot, m1_dia_plot)

ggplot(m1_DIT_plot, aes(x = aiE, y =N_plot, color = species) ) + 
  geom_line( )+ facet_grid(temperature~species)+ 
  geom_point(data= m1_DIT, mapping= aes(x = aiE, y =alg_per_N, color = species) )+
  facet_grid(temperature~species)+
  xlab("Active information (bits) ")+
  ylab("Algal consumption per individual")+
  theme(strip.background = element_rect(colour=NA, fill=NA))
ggsave("./aiE_algalperN.pdf", width = 8, height = 10)


# ggplot(m1_DIT_plot, aes(x = aiE, y =N_plot, color=interaction(temperature,species)), 
#   group=interaction(temperature,species) )  + 
#   geom_line( )+ facet_grid(.~species)+ 
#   geom_point(data= m1_DIT, mapping= aes(x = aiE, y =alg_per_N, color=interaction(temperature,species), 
#   group=interaction(temperature,species)) )+  facet_grid(.~species)+ 
#   xlab("Active information (bits) ")+
#   ylab("Algal consumption per individual")+
#   theme(strip.background = element_rect(colour=NA, fill=NA))
# ggsave("./aiE_algalperN_same.pdf", width = 8, height = 10)


#=============================================================================
# Plot of complexity (Excess Entropy) against ?????Cost????? per temperature
# treatment. Meant to look for something similar to Moore's curves of 
# integrated circuit complexity vs. manufacturing cost from 1965 paper. 
#=============================================================================
m1_DIT %>% 
  ungroup() %>% 
  filter(!is.na(N)) %>% 
  ggplot(aes(x = aiE, y =N, color = species, group = interaction(species,replicate_number)))+
  geom_point()+
  facet_grid(temperature~species)+
  #scale_y_log10()+
  theme_bw()+
  scale_color_manual(values = c("dodgerblue1", "red1"), name = NA)+
  xlab("Active information (bits) ")+
  ylab("Algal consumption per individual")+
  theme(strip.background = element_rect(colour=NA, fill=NA))
ggsave("./aiE_algalperN1.pdf", width = 8, height = 10)



#=============================================================================
# Plot each of the average information theoretic metrics as a bar graph
#=============================================================================

for(w in 1:(ntreatments*2)) {

fig.name = paste("average_dynamics_daphDia1k3",w,".pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

layout.matrix=matrix(c(1:4), nrow = 2, ncol = 2)
layout(mat = layout.matrix,
       heights = c(5,5), # Heights of the rows
       widths = c(5,5)) # Widths of columns

#layout.show(4)

barplot(di_web[[w]]$ee_means,cex.lab =1.3, beside = TRUE,ylab="Bits of information", xlab = "")
barplot(di_web[[w]]$ai_means,cex.lab =1.3, beside = TRUE,ylab="Bits of information", xlab = "Species #")
barplot(di_web[[w]]$te_means,cex.lab =1.3, beside = TRUE,ylab="", xlab = "")
barplot(di_web[[w]]$si_means,cex.lab =1.3, beside = TRUE,ylab="", xlab = "Species #")

dev.off()

}

#=============================================================================
# Plot a local DIT quantity alongside the population dynamics
#=============================================================================
#Active information and population dynamics
ai_out_all = data.frame( matrix(ncol = 6,nrow=0) ) 
colnames(ai_out_all) = c("temp", "group", "time", "species", "ai", "out1")
for(w in 2:(ntreatments)) {
  

  #First half of the treatments: 
  k=2
  nt1 = k
  if(nrow(out1[[w]])<=k){ k = 1; nt1=k}
  nt2 = dim(out1[[w]])[1]-1

  ai_comb = di_web[[w]]$ai_local
  out1_comb = out1[[w]][nt1:nt2,]/(0.5*max(out1[[w]][nt1:nt2,] ) ) #Scale by max? 
  #Turn ai into a data frame with headings matching ai_all
  colnames(ai_comb) = rspecies
  ai_comb=data.frame(ai_comb)
  ai_comb$temp = temps[w] #Add temperature 
  ai_comb$group = c("dia")
  ai_comb=ai_comb%>%mutate(time=row_number())
  ai_temp = ai_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=ai, rspecies)

  #Turn out1 into a data frame with headings matching ai_all
  colnames(out1_comb) = rspecies
  out1_comb=data.frame(out1_comb)
  out1_comb$temp = temps[w] #Add temperature 
  out1_comb$group = c("dia")
  out1_comb=out1_comb%>%mutate(time=row_number())
  out1_temp = out1_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=out1, rspecies)
  
  all_temp1 = ai_temp%>% left_join(out1_temp,all.x=T)

  #Second half of the treatments: 
  w2 = w+ntreatments
  k=2
  nt1 = k
  if(nrow(out1[[w2]])<=k){ k = 1; nt1=k}
  nt2 = dim(out1[[w2]])[1]-1

  ai_comb = di_web[[w2]]$ai_local
  out1_comb = out1[[w2]][nt1:nt2,]/(0.5*max(out1[[w2]][nt1:nt2,] ) ) #Scale by max? 
 
  #Turn ai into a data frame with headings matching ai_all
  colnames(ai_comb) = rspecies
  ai_comb=data.frame(ai_comb)
  ai_comb$temp = temps[w] #Add temperature 
  ai_comb$group = c("daph")
  ai_comb=ai_comb%>%mutate(time=row_number())
  ai_temp = ai_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=ai, rspecies)

  #Turn out1 into a data frame with headings matching ai_all
  colnames(out1_comb) = rspecies
  out1_comb=data.frame(out1_comb)
  out1_comb$temp = temps[w] #Add temperature 
  out1_comb$group = c("daph")
  out1_comb=out1_comb%>%mutate(time=row_number())
  out1_temp = out1_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=out1, rspecies)
  
  all_temp2 = ai_temp%>% left_join(out1_temp,all.x=T)

  all_temp = all_temp1%>% full_join(all_temp2)

  ai_out_all=rbind(ai_out_all,all_temp)
      
}

#Make the plot
ggplot()+ 
  geom_point(data = ai_out_all, aes(x = time, y = ai, color = interaction(group,"AI"), group = interaction(species,group)))+
  geom_line(data = ai_out_all, aes(x = time, y = ai, color = interaction(group,"AI"), group = interaction(species,group)))+
  geom_line(data = ai_out_all, aes(x = time, y = out1, color =interaction(group,"Pop"),  group = interaction(species,group)))+
  #scale_y_continuous(sec.axis = sec_axis(~.*5, name = "AI"))+
  facet_grid(temp~species)+
  #scale_y_log10()+
  theme_bw()+
  #scale_color_manual(values = c("dodgerblue1", "red1"), name = "")+
  xlab("Time")+
  ylab("Population/active information")+
  theme(strip.background = element_rect(colour=NA, fill=NA))
ggsave("./time_localAI_plot1.pdf", width = 8, height = 10)


#Transfer entropy and population dynamics
te_out_all = data.frame( matrix(ncol = 6,nrow=0) ) 
colnames(te_out_all) = c("temp", "group", "time", "species", "te", "out1")
for(w in 2:(ntreatments)) {
  

  #First half of the treatments: 
  k=2
  nt1 = k
  if(nrow(out1[[w]])<=k){ k = 1; nt1=k}
  nt2 = dim(out1[[w]])[1]-1

  te_comb = di_web[[w]]$te_local
  out1_comb = out1[[w]][nt1:nt2,]/(0.5*max(out1[[w]][nt1:nt2,] ) ) #Scale by max? 
  #Turn te into a data frame with headings matching te_all
  colnames(te_comb) = rspecies
  te_comb=data.frame(te_comb)
  te_comb$temp = temps[w] #Add temperature 
  te_comb$group = c("dia")
  te_comb=te_comb%>%mutate(time=row_number())
  te_temp = te_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=te, rspecies)

  #Turn out1 into a data frame with headings matching te_all
  colnames(out1_comb) = rspecies
  out1_comb=data.frame(out1_comb)
  out1_comb$temp = temps[w] #Add temperature 
  out1_comb$group = c("dia")
  out1_comb=out1_comb%>%mutate(time=row_number())
  out1_temp = out1_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=out1, rspecies)
  
  all_temp1 = te_temp%>% left_join(out1_temp,all.x=T)

  #Second half of the treatments: 
  w2 = w+ntreatments
  k=2
  nt1 = k
  if(nrow(out1[[w2]])<=k){ k = 1; nt1=k}
  nt2 = dim(out1[[w2]])[1]-1

  te_comb = di_web[[w2]]$te_local
  out1_comb = out1[[w2]][nt1:nt2,]/(0.5*max(out1[[w2]][nt1:nt2,] ) ) #Scale by max? 
 
  #Turn te into a data frame with headings matching te_all
  colnames(te_comb) = rspecies
  te_comb=data.frame(te_comb)
  te_comb$temp = temps[w] #Add temperature 
  te_comb$group = c("daph")
  te_comb=te_comb%>%mutate(time=row_number())
  te_temp = te_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=te, rspecies)

  #Turn out1 into a data frame with headings matching te_all
  colnames(out1_comb) = rspecies
  out1_comb=data.frame(out1_comb)
  out1_comb$temp = temps[w] #Add temperature 
  out1_comb$group = c("daph")
  out1_comb=out1_comb%>%mutate(time=row_number())
  out1_temp = out1_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=out1, rspecies)
  
  all_temp2 = te_temp%>% left_join(out1_temp,all.x=T)

  all_temp = all_temp1%>% full_join(all_temp2)

  te_out_all=rbind(te_out_all,all_temp)
      
}

#Make the plot
ggplot()+ 
  geom_line(data = te_out_all, aes(x = time, y = te, color = interaction(group,"te"), group = interaction(species,group)))+
  geom_point(data = te_out_all, aes(x = time, y = te, color = interaction(group,"te"), group = interaction(species,group)))+
  geom_line(data = te_out_all, aes(x = time, y = out1, color =interaction(group,"Pop"),  group = interaction(species,group)))+
  #scale_y_continuous(sec.axis = sec_axis(~.*5, name = "te"))+
  facet_grid(temp~species)+
  #scale_y_log10()+
  theme_bw()+
  #scale_color_manual(values = c("dodgerblue1", "red1"), name = "")+
  xlab("Time")+
  ylab("Population/transfer entropy")+
  theme(strip.background = element_rect(colour=NA, fill=NA))
ggsave("./time_localTE_plot1.pdf", width = 8, height = 10)



#=============================================================================
# Plot of complexity (Excess Entropy) against ?????Cost????? per temperature
# treatment. Meant to look for something similar to Moore's curves of 
# integrated circuit complexity vs. manufacturing cost from 1965 paper. 
#=============================================================================
#Active information vs. number of individuals 
ggplot()+ 
  geom_point(data = ai_out_all, aes(x = ai, y = out1, color = interaction(group,"AI"), group = interaction(species,group)))+
  #scale_y_continuous(sec.axis = sec_axis(~.*5, name = "AI"))+
  facet_grid(temp~species)+
  #scale_y_log10()+
  theme_bw()+
  #scale_color_manual(values = c("dodgerblue1", "red1"), name = "")+
  xlab("active information (bits)")+
  ylab("Population")+
  theme(strip.background = element_rect(colour=NA, fill=NA))
ggsave("./pop_localAI_plot1.pdf", width = 8, height = 10)

#=============================================================================
#Transfer entropy vs. number of individuals 
ggplot()+ 
  geom_point(data = te_out_all, aes(x = te, y = out1, color = interaction(group,"TE"), group = interaction(species,group)))+
  #scale_y_continuous(sec.axis = sec_axis(~.*5, name = "AI"))+
  facet_grid(temp~species)+
  #scale_y_log10()+
  theme_bw()+
  #scale_color_manual(values = c("dodgerblue1", "red1"), name = "")+
  xlab("transfer entropy (bits)")+
  ylab("Population")+
  theme(strip.background = element_rect(colour=NA, fill=NA))
ggsave("./pop_localTE_plot1.pdf", width = 8, height = 10)


#=============================================================================
#Ensemble active information and algal consumption: 
ai_alg_all = data.frame( matrix(ncol = 6,nrow=0) ) 
colnames(ai_alg_all) = c("temp", "group", "time", "species", "aiE", "out1_alg")
for(w in 1:nmesos) {
  

  #First half of the treatments: 
  k=2
  nt1 = k-1
  if(nrow(out1_alg[[w]])<=k){ k = 1; nt1=k}
  nt2 = dim(out1_alg[[w]])[1]-1
  din = nt2-nt1+1

  ai_comb = aiE_web[[w]]$local[1:din,]
  out1_alg_comb = out1_alg[[w]][nt1:nt2,]/(0.5*max(out1_alg[[w]][nt1:nt2,] ) ) #Scale by max? 
  #Turn ai into a data frame with headings matching ai_all
  colnames(ai_comb) = rspecies
  ai_comb=data.frame(ai_comb)
  ai_comb$temp = temps[w] #Add temperature 
  ai_comb$group = c("dia")
  ai_comb=ai_comb%>%mutate(time=row_number())
  ai_temp = ai_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=aiE, rspecies)

  #Turn out1_alg into a data frame with headings matching ai_all
  colnames(out1_alg_comb) = rspecies
  out1_alg_comb=data.frame(out1_alg_comb)
  out1_alg_comb$temp = temps[w] #Add temperature 
  out1_alg_comb$group = c("dia")
  out1_alg_comb=out1_alg_comb%>%mutate(time=row_number())
  out1_alg_temp = out1_alg_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=out1_alg, rspecies)
  
  all_temp1 = ai_temp%>% left_join(out1_alg_temp,all.x=T)

  #Second half of the treatments: 
  w2 = w+ntreatments
  k=2
  nt1 = k-1
  if(nrow(out1_alg[[w2]])<=k){ k = 1; nt1=k}
  nt2 = dim(out1_alg[[w2]])[1]-1
  din = nt2-nt1+1

  ai_comb = di_web[[w2]]$ai_local[1:din,]
  out1_alg_comb = out1_alg[[w2]][nt1:nt2,]/(0.5*max(out1_alg[[w2]][nt1:nt2,] ) ) #Scale by max? 
 
  #Turn ai into a data frame with headings matching ai_all
  colnames(ai_comb) = rspecies
  ai_comb=data.frame(ai_comb)
  ai_comb$temp = temps[w] #Add temperature 
  ai_comb$group = c("daph")
  ai_comb=ai_comb%>%mutate(time=row_number())
  ai_temp = ai_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=aiE, rspecies)

  #Turn out1_alg into a data frame with headings matching ai_all
  colnames(out1_alg_comb) = rspecies
  out1_alg_comb=data.frame(out1_alg_comb)
  out1_alg_comb$temp = temps[w] #Add temperature 
  out1_alg_comb$group = c("daph")
  out1_alg_comb=out1_alg_comb%>%mutate(time=row_number())
  out1_alg_temp = out1_alg_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=out1_alg, rspecies)
  
  all_temp2 = ai_temp%>% left_join(out1_alg_temp,all.x=T)

  all_temp = all_temp1%>% full_join(all_temp2)

  ai_alg_all=rbind(ai_alg_all,all_temp)
      
}

#Make the plot
ggplot()+ 
  geom_point(data = ai_alg_all, aes(x = ai, y = out1_alg, color = interaction(group,"AI"), group = interaction(species,group)))+
  #scale_y_continuous(sec.axis = sec_axis(~.*5, name = "AI"))+
  facet_grid(temp~species)+
  #scale_y_log10()+
  theme_bw()+
  #scale_color_manual(values = c("dodgerblue1", "red1"), name = "")+
  xlab("active information (bits) ")+
  ylab("1st derivative (individuals)")+
  theme(strip.background = element_rect(colour=NA, fill=NA))
ggsave("./alg_localAI_plot1.pdf", width = 8, height = 10)


#=============================================================================
#Active information and slope of population growth curve
ai_outdiff_all = data.frame( matrix(ncol = 6,nrow=0) ) 
colnames(ai_outdiff_all) = c("temp", "group", "time", "species", "ai", "out1_diff")
for(w in 2:(ntreatments)) {
  

  #First half of the treatments: 
  k=2
  nt1 = k-1
  if(nrow(out1_diff[[w]])<=k){ k = 1; nt1=k}
  nt2 = dim(out1_diff[[w]])[1]-1
  din = nt2-nt1+1

  ai_comb = di_web[[w]]$ai_local[1:din,]
  out1_diff_comb = out1_diff[[w]][nt1:nt2,]/(0.5*max(out1_diff[[w]][nt1:nt2,] ) ) #Scale by max? 
  #Turn ai into a data frame with headings matching ai_all
  colnames(ai_comb) = rspecies
  ai_comb=data.frame(ai_comb)
  ai_comb$temp = temps[w] #Add temperature 
  ai_comb$group = c("dia")
  ai_comb=ai_comb%>%mutate(time=row_number())
  ai_temp = ai_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=ai, rspecies)

  #Turn out1_diff into a data frame with headings matching ai_all
  colnames(out1_diff_comb) = rspecies
  out1_diff_comb=data.frame(out1_diff_comb)
  out1_diff_comb$temp = temps[w] #Add temperature 
  out1_diff_comb$group = c("dia")
  out1_diff_comb=out1_diff_comb%>%mutate(time=row_number())
  out1_diff_temp = out1_diff_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=out1_diff, rspecies)
  
  all_temp1 = ai_temp%>% left_join(out1_diff_temp,all.x=T)

  #Second half of the treatments: 
  w2 = w+ntreatments
  k=2
  nt1 = k-1
  if(nrow(out1_diff[[w2]])<=k){ k = 1; nt1=k}
  nt2 = dim(out1_diff[[w2]])[1]-1
  din = nt2-nt1+1

  ai_comb = di_web[[w2]]$ai_local[1:din,]
  out1_diff_comb = out1_diff[[w2]][nt1:nt2,]/(0.5*max(out1_diff[[w2]][nt1:nt2,] ) ) #Scale by max? 
 
  #Turn ai into a data frame with headings matching ai_all
  colnames(ai_comb) = rspecies
  ai_comb=data.frame(ai_comb)
  ai_comb$temp = temps[w] #Add temperature 
  ai_comb$group = c("daph")
  ai_comb=ai_comb%>%mutate(time=row_number())
  ai_temp = ai_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=ai, rspecies)

  #Turn out1_diff into a data frame with headings matching ai_all
  colnames(out1_diff_comb) = rspecies
  out1_diff_comb=data.frame(out1_diff_comb)
  out1_diff_comb$temp = temps[w] #Add temperature 
  out1_diff_comb$group = c("daph")
  out1_diff_comb=out1_diff_comb%>%mutate(time=row_number())
  out1_diff_temp = out1_diff_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=out1_diff, rspecies)
  
  all_temp2 = ai_temp%>% left_join(out1_diff_temp,all.x=T)

  all_temp = all_temp1%>% full_join(all_temp2)

  ai_outdiff_all=rbind(ai_outdiff_all,all_temp)
      
}

#Make the plot
ggplot()+ 
  geom_point(data = ai_outdiff_all, aes(x = ai, y = out1_diff, color = interaction(group,"AI"), group = interaction(species,group)))+
  #scale_y_continuous(sec.axis = sec_axis(~.*5, name = "AI"))+
  facet_grid(temp~species)+
  #scale_y_log10()+
  theme_bw()+
  #scale_color_manual(values = c("dodgerblue1", "red1"), name = "")+
  xlab("active information (bits) ")+
  ylab("1st derivative (individuals)")+
  theme(strip.background = element_rect(colour=NA, fill=NA))
ggsave("./diff_localAI_plot1.pdf", width = 8, height = 10)


#=============================================================================
#Transfer entropy and population dynamics
te_outdiff_all = data.frame( matrix(ncol = 6,nrow=0) ) 
colnames(te_outdiff_all) = c("temp", "group", "time", "species", "te", "out1_diff")
for(w in 2:(ntreatments)) {
  

 #First half of the treatments: 
  k=2
  nt1 = k-1
  if(nrow(out1_diff[[w]])<=k){ k = 1; nt1=k}
  nt2 = dim(out1_diff[[w]])[1]-1
  din = nt2-nt1+1


  te_comb = di_web[[w]]$te_local[1:din,]
  out1_diff_comb = out1_diff[[w]][nt1:nt2,]/(0.5*max(out1_diff[[w]][nt1:nt2,] ) ) #Scale by max? 
  #Turn te into a data frame with headings matching te_all
  colnames(te_comb) = rspecies
  te_comb=data.frame(te_comb)
  te_comb$temp = temps[w] #Add temperature 
  te_comb$group = c("dia")
  te_comb=te_comb%>%mutate(time=row_number())
  te_temp = te_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=te, rspecies)

  #Turn out1_diff into a data frame with headings matching te_all
  colnames(out1_diff_comb) = rspecies
  out1_diff_comb=data.frame(out1_diff_comb)
  out1_diff_comb$temp = temps[w] #Add temperature 
  out1_diff_comb$group = c("dia")
  out1_diff_comb=out1_diff_comb%>%mutate(time=row_number())
  out1_diff_temp = out1_diff_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=out1_diff, rspecies)
  
  all_temp1 = te_temp%>% left_join(out1_diff_temp,all.x=T)

  #Second half of the treatments: 
  w2 = w+ntreatments
  k=2
  nt1 = k-1
  if(nrow(out1_diff[[w2]])<=k){ k = 1; nt1=k}
  nt2 = dim(out1_diff[[w2]])[1]-1
  din = nt2-nt1+1

  te_comb = di_web[[w2]]$te_local[1:din,]
  out1_diff_comb = out1_diff[[w2]][nt1:nt2,]/(0.5*max(out1_diff[[w2]][nt1:nt2,] ) ) #Scale by max? 
 
  #Turn te into a data frame with headings matching te_all
  colnames(te_comb) = rspecies
  te_comb=data.frame(te_comb)
  te_comb$temp = temps[w] #Add temperature 
  te_comb$group = c("daph")
  te_comb=te_comb%>%mutate(time=row_number())
  te_temp = te_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=te, rspecies)

  #Turn out1_diff into a data frame with headings matching te_all
  colnames(out1_diff_comb) = rspecies
  out1_diff_comb=data.frame(out1_diff_comb)
  out1_diff_comb$temp = temps[w] #Add temperature 
  out1_diff_comb$group = c("daph")
  out1_diff_comb=out1_diff_comb%>%mutate(time=row_number())
  out1_diff_temp = out1_diff_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=out1_diff, rspecies)
  
  all_temp2 = te_temp%>% left_join(out1_diff_temp,all.x=T)

  all_temp = all_temp1%>% full_join(all_temp2)

  te_outdiff_all=rbind(te_outdiff_all,all_temp)
      
}

#Make the plot
ggplot()+ 
  geom_point(data = te_outdiff_all, aes(x = te, y = out1_diff, color = interaction(group,"TE"), group = interaction(species,group)))+

   #scale_y_continuous(sec.axis = sec_axis(~.*5, name = "te"))+
  facet_grid(temp~species)+
  #scale_y_log10()+
  theme_bw()+
  #scale_color_manual(values = c("dodgerblue1", "red1"), name = "")+
  xlab("active information (bits) ")+
  ylab("1st derivative (individuals)")+
  theme(strip.background = element_rect(colour=NA, fill=NA))
ggsave("./diff_localTE_plot1.pdf", width = 8, height = 10)



#============================================================================================
#============================================================================================




# fig.name = paste("moore_curve_DaphDia1.pdf",sep="")
# pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)
ee_all = NULL
pops_all = NULL
ee_all2 = NULL
pops_all2 = NULL


layout.matrix=matrix(c(1:(ntreatments*4)), nrow = ntreatments, ncol = 4)
layout(mat = layout.matrix,
       heights = c(matrix(5,ntreatments*2,1)), # Heights of the rows
       widths = c(matrix(5,4,1)) ) # Widths of columns

#Active information vs population
for(w in 2:(ntreatments*2)) {
  k=2
  nt1 = k
  if(nrow(out1[[w]])<=k){ k = 1; nt1=k}
  nt2 = dim(out1[[w]])[1]-1
  plot(di_web[[w]]$ai_local[,1],out1[[w]][nt1:nt2,1],xlim = c(0,8),ylim=c(0,500) )
  plot(di_web[[w]]$ai_local[,2],out1[[w]][nt1:nt2,2],xlim = c(0,8),ylim=c(0,500)   )
}


#Active information vs slope
ai_all = data.frame( matrix(ncol = 5,nrow=0) ) 
colnames(ai_all) = c("temp", "time", "species", "ai", "out1")
for(w in 2:(ntreatments)) {
  

  #First half of the treatments: 
  k=2
  nt1 = k-1
  if(nrow(out1_diff[[w]])<=k){ k = 1; nt1=k}
  nt2 = dim(out1_diff[[w]])[1]-1
  din = nt2-nt1+1

  ai_comb = di_web[[w]]$ai_local[1:din,]
  out1_comb = out1_diff[[w]][nt1:nt2,]

  #Second half of the treatments: 
  w2 = w+ntreatments
  k=2
  nt1 = k-1
  if(nrow(out1_diff[[w2]])<=k){ k = 1; nt1=k}
  nt2 = dim(out1_diff[[w2]])[1]-1
  din = nt2-nt1+1

  #Turn ai into a data frame with headings matching ai_all
  ai_comb = rbind( ai_comb, di_web[[w2]]$ai_local[1:din,])
  colnames(ai_comb) = rspecies
  ai_comb=data.frame(ai_comb)
  ai_comb$temp = temps[w] #Add temperature 
  ai_comb=ai_comb%>%mutate(time=row_number())
  ai_temp = ai_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=ai, rspecies)

  #Turn out1 into a data frame with headings matching ai_all
  out1_comb = rbind( out1_comb, out1_diff[[w2]][nt1:nt2,])
  colnames(out1_comb) = rspecies
  out1_comb=data.frame(out1_comb)
  out1_comb$temp = temps[w] #Add temperature 
  out1_comb=out1_comb%>%mutate(time=row_number())
  out1_temp = out1_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=out1, rspecies)
  
  all_temp = ai_temp%>% left_join(out1_temp)

  ai_all=rbind(ai_all,all_temp)
      
}

#Make the plot
ai_all %>% 
  ungroup() %>% 
  #filter(!is.na(N)) %>% 
  ggplot(aes(x = ai, y = out1, color = species))+
  geom_point()+
  facet_grid(temp~species)+
  #scale_y_log10()+
  theme_bw()+
  scale_color_manual(values = c("dodgerblue1", "red1"), name = "species")+
  xlab("active information")+
  ylab("1st derivative of growth")+
  theme(strip.background = element_rect(colour=NA, fill=NA))
ggsave("./AI_plot1.pdf", width = 8, height = 10)
Collapse


#Active information vs slope
for(w in 2:(ntreatments*2)) {
  k=2
  nt1 = k-1
  if(nrow(out1_diff[[w]])<=k){ k = 1; nt1=k}
  nt2 = dim(out1_diff[[w]])[1]-1
  din = nt2-nt1+1
  plot(di_web[[w]]$ai_local[1:din,1],out1_diff[[w]][nt1:nt2,1],xlim = c(0,8),ylim=c(0,500) )
  plot(di_web[[w]]$ai_local[1:din,2],out1_diff[[w]][nt1:nt2,2],xlim = c(0,8),ylim=c(0,500)   )
}


#Excess Entropy
for(w in 2:(ntreatments*2)) {
  k=2
  nt1 = k
  if(nrow(out1[[w]])<=k){ k = 1; nt1=k}
  nt2 = dim(out1[[w]])[1]-k
  plot(di_web[[w]]$ee_local[,1],out1[[w]][nt1:nt2,1],xlim = c(0,8),ylim=c(0,500) )
  plot(di_web[[w]]$ee_local[,2],out1[[w]][nt1:nt2,2],xlim = c(0,8),ylim=c(0,500)   )
}


for(w in 1:(ntreatments*2)) {
  k=3
  nt1 = k
  if(nrow(out1[[w]])<=k){ k = 1; nt1=k}
  nt2 = dim(out1[[w]])[1]-k
  ee_all = rbind(ee_all,di_web[[w]]$ee_local )
  pops_all = rbind(pops_all,out1[[w]][nt1:nt2,] )
}


layout.matrix=matrix(c(1:4), nrow = 2, ncol = 2)
layout(mat = layout.matrix,
       heights = c(5,5), # Heights of the rows
       widths = c(5,5)) # Widths of columns

#layout.show(4)

barplot(di_web[[w]]$ee_means,cex.lab =1.3, beside = TRUE,ylab="Bits of information", xlab = "")
barplot(di_web[[w]]$ai_means,cex.lab =1.3, beside = TRUE,ylab="Bits of information", xlab = "Species #")
barplot(di_web[[w]]$te_means,cex.lab =1.3, beside = TRUE,ylab="", xlab = "")
barplot(di_web[[w]]$si_means,cex.lab =1.3, beside = TRUE,ylab="", xlab = "Species #")

dev.off()
