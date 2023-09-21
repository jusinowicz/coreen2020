#=============================================================================
# Download and prepare data
#Right now this is set up with source "av" = Alpha Vantage: Our standard API 
#call frequency is 5 calls per minute and 500 calls per day. Please visit 
#https://www.alphavantage.co/premium/ if you would like to target a higher 
#API call frequency.
#My API key: M6EXVU1REBBKKLJK
#=============================================================================

yourKey = "M6EXVU1REBBKKLJK"
#Download multiple data sets for analysis:
symbols = c("SPY","IWM","QQQ")
ns = length(symbols)
getSymbols(symbols, src="av", api.key="yourKey", output.size="full",
  periodicity="intraday") 

#Merge them into one xts object with a shared time index.
c1 = paste(symbols, collapse =",")
c2 = paste("merge(",c1,")") #Create a function call with all symbols
sl1 = eval(parse(text=c2)) #Use eval(parse(text= )) to interpret the text

#Interpolate NAs so that analyses are meaningful
sl1.int = na.approx(sl1)

#Turn an xts object into a data.frame with zoo
sl2 = fortify.zoo(sl1.int)


#=============================================================================
# Dynamic Information Theoretic analysis of data. 
#
# Calculate Excess Entropy, Active Information Storage, 
# and Transfer Entropy.
# Each quantity is calculated at both the average and local level.  
#=============================================================================
#Variables: 
#Set k, the block size: 
k=2



#=============================================================================
# Variables for data collection
#=============================================================================
#Real data: 
#=============================================================================
#Main data frame per ticker
out1R = vector("list",ns)

#Dynamic information metrics calculated from the time series 
di_all =  vector("list",ns)

#Track the average transfer entropy and separable information between each pair of 
#stocks as a way to build a network of information flow through the network. 
te_all =  vector("list",ns)
si_all= vector("list",ns) 

#The ensemble version of the AI
aiE_all = vector("list",ns)
