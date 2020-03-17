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
library(fields)
source("./functions/info_theory_functions.R")
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
nspp=2
temps = (unique(m1_data_long$temperature))
ntreatments =  length(unique(m1_data_long$temperature))
treats = unique(m1_data_long$mesocosm.ID)
invader = unique(m1_data_long$species)
invasions_per = length(treats)/4 #Treatment entries correponding to each spp 
inv_day = 28 #The first day of attempted invasion
no_reps = 18 #The number of replicated mesocosms total per resident/invader

#=============================================================================
# Variables for data collection
#=============================================================================
out1 = vector("list",ntreatments*2)
out1_diff = vector("list",ntreatments*2)
#Dynamic information metrics calculated from the (discretized) time series 
di_web =  vector("list",ntreatments*2)
#Track the average transfer entropy and separable information between each pair of 
#species as a way to build a network of information flow through the network. 
te_web =  vector("list",ntreatments*2)
si_web = vector("list",ntreatments*2) 

#=============================================================================
# First, loop over both diaphanosoma and daphnia when diaphanasoma is the 
# invader
#=============================================================================
#The resident and invader: 
for (i in 1:2) {  
  #Loop over treatments 
  u_invader = invader[-i]
  u_res = invader[i]

  for (n in 1: (ntreatments)){
    #Since the treatments are regular, set the start/end position of
    #the subset with:  
    if(i==1){ pos1 = (n-1)*2 + n } else {pos1 = (n-1)*2 + n+invasions_per}
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
            arrange(replicate.number)%>%
            mutate(Ndiff_res = lead( N-lag(N),)/lead(day_n-lag(day_n))) #"lead" lines up the result 
    res_tmp$Ndiff_res[res_tmp$day_n == max(res_tmp$day_n )] =NA #Remove last day
    res_tmp = res_tmp %>% mutate(tdiff =day_n-lag(day_n)) #Size of time step
    res_tmp$tdiff[res_tmp$tdiff<0] = 1 #Remove negative time steplag(N)s
    res_tmp = res_tmp %>% mutate(N_res = N)


    #Make a new invader data set to fit the growth rate function with nlme/nls 
    #The new data set includes a column for n(t+1)/n(t) and a columnf for the time 
    #interval across subsequent measurements.  
    inv_tmp = subset(m1_data_long, mesocosm.ID == u_treats &  species == u_invader   )
        inv_tmp = na.exclude(inv_tmp) #NLME won't work with NAs 
        #inv_tmp$N[is.na(inv_tmp$N)] = 0 #Replace NAs with 0  
        #Arrange the data by replicate number, then add a new column for the delta N
        inv_tmp = inv_tmp %>% 
            arrange(replicate.number)%>%
            mutate(Ndiff_inv = lead( N-lag(N),)/lead(day_n-lag(day_n))) #"lead" lines up the result 
    inv_tmp$Ndiff_inv[inv_tmp$day_n == max(inv_tmp$day_n )] =NA #Remove last day
    inv_tmp = inv_tmp %>% mutate(tdiff =day_n-lag(day_n)) #Size of time step 
    inv_tmp = inv_tmp %>% mutate(N_inv = N)
    #This line adds the resident densities on the matching days from the matching 
    #replicates to the data table inv_tmp
    inv_tmp = inv_tmp %>% left_join( select( res_tmp, replicate.number,day_n,N_res,Ndiff_res),
     by= (c( "day_n" = "day_n", "replicate.number"="replicate.number" )) )
    res_tmp = res_tmp %>% left_join( select( inv_tmp, replicate.number,day_n,N_inv,Ndiff_inv),
       by= (c( "day_n" = "day_n", "replicate.number"="replicate.number" )) )
    res_tmp$N_inv[is.na(res_tmp$N_inv)] = 0 #Replace NAs with 0  

    #Create a dummy data set in the few cases that experiment was essentially 
    #unsuccesful: 
    if(dim(inv_tmp)[1] <=0  ){
      
      mydata_res = subset(res_tmp, day_n >= inv_day  )
      mydata_res$Ndiff_res[is.na(mydata_res$Ndiff_res)] = 0 #Replace NAs with 0 

      #Population data
      mydata = matrix(0,dim(mydata_res)[1],nspp)
      colnames(mydata) = invader
      mydata[,u_invader] = mydata_res$N_inv
      mydata[,u_res] = mydata_res$N_res
      out1[index1] = list(as.matrix(na.exclude(mydata)))
      
      #Population diff data
      mydata_diff = matrix(0,dim(mydata_res)[1],nspp)
      colnames(mydata_diff) = invader
      mydata_diff[,u_invader] = mydata_res$Ndiff_res
      mydata_diff[,u_res] = mydata_res$Ndiff_res
      out1_diff[index1] = list(as.matrix(na.exclude(mydata_diff)))
      

    } else if (dim(res_tmp)[1] <=0  ){
      mydata_inv = subset(inv_tmp, day_n >= inv_day  )

      #Population data
      mydata_inv$N_res[is.na(mydata_inv$N_res)] = 0 #Replace NAs with 0 
      mydata = matrix(0,dim(mydata_inv)[1],nspp)
      colnames(mydata) = invader
      mydata[,u_invader] = mydata_inv$N_inv
      mydata[,u_res] = mydata_inv$N_res
      out1[index1] = list(as.matrix(na.exclude(mydata)))

      #Population diff data
      mydata_inv$Ndiff_res[is.na(mydata_inv$Ndiff_res)] = 0 #Replace NAs with 0 
      mydata_diff = matrix(0,dim(mydata_inv)[1],nspp)
      colnames(mydata_diff) = invader
      mydata_diff[,u_invader] = mydata_inv$Ndiff_inv
      mydata_diff[,u_res] = mydata_inv$Ndiff_res
      out1_diff[index1] = list(as.matrix(na.exclude(mydata_diff)))


    } else {

      mydata_inv = subset(inv_tmp, day_n >= inv_day  )

      #Population data
      mydata = matrix(0,dim(mydata_inv)[1],nspp)
      colnames(mydata) = invader
      mydata[,u_invader] = mydata_inv$N_inv
      mydata[,u_res] = mydata_inv$N_res
      out1[index1] = list(as.matrix(na.exclude(mydata)))

      #Population diff data
      mydata_diff = matrix(0,dim(mydata_inv)[1],nspp)
      colnames(mydata_diff) = invader
      mydata_diff[,u_invader] = mydata_inv$Ndiff_inv
      mydata_diff[,u_res] = mydata_inv$Ndiff_res
      out1_diff[index1] = list(as.matrix(na.exclude(mydata_diff)))
      
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
    if(nrow(out1[[index1]])<=k){ k = 1}
    #

    nt1 = 1
    nt2 = dim(out1[[index1]])[1]
    f1 = 1 #scaling factor
    di_web[index1] = list(get_info_dynamics(pop_ts = floor(f1*out1[[index1]][nt1:nt2,]), 
      k=k,with_blocks=TRUE))

    ## This code takes the population time-series counts output by the ODEs and 
    ## calculates the average Transfer Entropy from each species to every other 
    ## species. The goal is to get an overview of the major information pathways 
    ## in the web.   
    #=============================================================================
    # This function gives:
    # te_web    Average transfer entropy per species as a pairwise matrix
    #=============================================================================
    te_web[index1] = list( get_te_web( pop_ts = floor(f1*out1[[index1]][nt1:nt2,]), 
      k=k) )

  }
}

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
  colnames(ai_comb) = invader
  ai_comb=data.frame(ai_comb)
  ai_comb$temp = temps[w] #Add temperature 
  ai_comb$group = c("dia")
  ai_comb=ai_comb%>%mutate(time=row_number())
  ai_temp = ai_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=ai, invader)

  #Turn out1 into a data frame with headings matching ai_all
  colnames(out1_comb) = invader
  out1_comb=data.frame(out1_comb)
  out1_comb$temp = temps[w] #Add temperature 
  out1_comb$group = c("dia")
  out1_comb=out1_comb%>%mutate(time=row_number())
  out1_temp = out1_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=out1, invader)
  
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
  colnames(ai_comb) = invader
  ai_comb=data.frame(ai_comb)
  ai_comb$temp = temps[w] #Add temperature 
  ai_comb$group = c("daph")
  ai_comb=ai_comb%>%mutate(time=row_number())
  ai_temp = ai_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=ai, invader)

  #Turn out1 into a data frame with headings matching ai_all
  colnames(out1_comb) = invader
  out1_comb=data.frame(out1_comb)
  out1_comb$temp = temps[w] #Add temperature 
  out1_comb$group = c("daph")
  out1_comb=out1_comb%>%mutate(time=row_number())
  out1_temp = out1_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=out1, invader)
  
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
  ylab("Population")+
  theme(strip.background = element_rect(colour=NA, fill=NA))
ggsave("./pop_localAI_plot1.pdf", width = 8, height = 10)


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
  colnames(te_comb) = invader
  te_comb=data.frame(te_comb)
  te_comb$temp = temps[w] #Add temperature 
  te_comb$group = c("dia")
  te_comb=te_comb%>%mutate(time=row_number())
  te_temp = te_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=te, invader)

  #Turn out1 into a data frame with headings matching te_all
  colnames(out1_comb) = invader
  out1_comb=data.frame(out1_comb)
  out1_comb$temp = temps[w] #Add temperature 
  out1_comb$group = c("dia")
  out1_comb=out1_comb%>%mutate(time=row_number())
  out1_temp = out1_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=out1, invader)
  
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
  colnames(te_comb) = invader
  te_comb=data.frame(te_comb)
  te_comb$temp = temps[w] #Add temperature 
  te_comb$group = c("daph")
  te_comb=te_comb%>%mutate(time=row_number())
  te_temp = te_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=te, invader)

  #Turn out1 into a data frame with headings matching te_all
  colnames(out1_comb) = invader
  out1_comb=data.frame(out1_comb)
  out1_comb$temp = temps[w] #Add temperature 
  out1_comb$group = c("daph")
  out1_comb=out1_comb%>%mutate(time=row_number())
  out1_temp = out1_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=out1, invader)
  
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
  ylab("Population")+
  theme(strip.background = element_rect(colour=NA, fill=NA))
ggsave("./pop_localTE_plot1.pdf", width = 8, height = 10)

#=============================================================================
# Plot of complexity (Excess Entropy) against ?????Cost????? per temperature
# treatment. Meant to look for something similar to Moore's curves of 
# integrated circuit complexity vs. manufacturing cost from 1965 paper. 
#=============================================================================
#Active information and population dynamics
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
  colnames(ai_comb) = invader
  ai_comb=data.frame(ai_comb)
  ai_comb$temp = temps[w] #Add temperature 
  ai_comb$group = c("dia")
  ai_comb=ai_comb%>%mutate(time=row_number())
  ai_temp = ai_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=ai, invader)

  #Turn out1_diff into a data frame with headings matching ai_all
  colnames(out1_diff_comb) = invader
  out1_diff_comb=data.frame(out1_diff_comb)
  out1_diff_comb$temp = temps[w] #Add temperature 
  out1_diff_comb$group = c("dia")
  out1_diff_comb=out1_diff_comb%>%mutate(time=row_number())
  out1_diff_temp = out1_diff_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=out1_diff, invader)
  
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
  colnames(ai_comb) = invader
  ai_comb=data.frame(ai_comb)
  ai_comb$temp = temps[w] #Add temperature 
  ai_comb$group = c("daph")
  ai_comb=ai_comb%>%mutate(time=row_number())
  ai_temp = ai_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=ai, invader)

  #Turn out1_diff into a data frame with headings matching ai_all
  colnames(out1_diff_comb) = invader
  out1_diff_comb=data.frame(out1_diff_comb)
  out1_diff_comb$temp = temps[w] #Add temperature 
  out1_diff_comb$group = c("daph")
  out1_diff_comb=out1_diff_comb%>%mutate(time=row_number())
  out1_diff_temp = out1_diff_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=out1_diff, invader)
  
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
  xlab("Time")+
  ylab("Population")+
  theme(strip.background = element_rect(colour=NA, fill=NA))
ggsave("./diff_localAI_plot1.pdf", width = 8, height = 10)


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
  colnames(te_comb) = invader
  te_comb=data.frame(te_comb)
  te_comb$temp = temps[w] #Add temperature 
  te_comb$group = c("dia")
  te_comb=te_comb%>%mutate(time=row_number())
  te_temp = te_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=te, invader)

  #Turn out1_diff into a data frame with headings matching te_all
  colnames(out1_diff_comb) = invader
  out1_diff_comb=data.frame(out1_diff_comb)
  out1_diff_comb$temp = temps[w] #Add temperature 
  out1_diff_comb$group = c("dia")
  out1_diff_comb=out1_diff_comb%>%mutate(time=row_number())
  out1_diff_temp = out1_diff_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=out1_diff, invader)
  
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
  colnames(te_comb) = invader
  te_comb=data.frame(te_comb)
  te_comb$temp = temps[w] #Add temperature 
  te_comb$group = c("daph")
  te_comb=te_comb%>%mutate(time=row_number())
  te_temp = te_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=te, invader)

  #Turn out1_diff into a data frame with headings matching te_all
  colnames(out1_diff_comb) = invader
  out1_diff_comb=data.frame(out1_diff_comb)
  out1_diff_comb$temp = temps[w] #Add temperature 
  out1_diff_comb$group = c("daph")
  out1_diff_comb=out1_diff_comb%>%mutate(time=row_number())
  out1_diff_temp = out1_diff_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=out1_diff, invader)
  
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
  xlab("Time")+
  ylab("Population")+
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
  colnames(ai_comb) = invader
  ai_comb=data.frame(ai_comb)
  ai_comb$temp = temps[w] #Add temperature 
  ai_comb=ai_comb%>%mutate(time=row_number())
  ai_temp = ai_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=ai, invader)

  #Turn out1 into a data frame with headings matching ai_all
  out1_comb = rbind( out1_comb, out1_diff[[w2]][nt1:nt2,])
  colnames(out1_comb) = invader
  out1_comb=data.frame(out1_comb)
  out1_comb$temp = temps[w] #Add temperature 
  out1_comb=out1_comb%>%mutate(time=row_number())
  out1_temp = out1_comb%>%
    mutate(time=row_number())%>%
      gather(key=species, value=out1, invader)
  
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
