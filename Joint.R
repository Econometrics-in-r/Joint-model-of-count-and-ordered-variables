rm(list = ls())
require(stats4)
require(maxLik)
require(randtoolbox)

# reading and storing data in a dataframe
dataset <- read.csv(file.choose(),header=T)
#Sample Size
N <- nrow(dataset) 
#Dependent variable (crash counts in this example); change the variable as required
DVar <- dataset$Crash 

# subcounts of the overall count variable
subc1 <- dataset[,'Fatal']#(subcount of first category: fatal crashes in this example); change the variable as required
subc2 <- dataset[,'Major']#(subcount of second category: major crashes in this example); change the variable as required
subc3 <- dataset[,'Minor']#(subcount of third category: minor crashes in this example); change the variable as required

#Matrix of fractions of dependent variable (for three categories; adjust as needed)
w=matrix(NA,nrow = N,ncol = 3)

for (i in 1:N) {
  if ((subc1[i] + subc2[i]	+ subc3[i])==0){
    w[i,1] <- 0
    w[i,2] <- 0
    w[i,3] <- 0  
  }
  else {
    w[i,1] <- subc1[i]/(subc1[i] + subc2[i]	+ subc3[i])
    w[i,2] <- subc2[i]/(subc1[i] + subc2[i]	+ subc3[i])
    w[i,3] <- subc3[i]/(subc1[i] + subc2[i]	+ subc3[i])
  }}

# #standardization, comment out if not used
# namesVar = c("AADT") # name of variables to be standardized
# meanVar = NULL
# sdVar = NULL
# for (namei in namesVar){
#   meanVar = c(meanVar,mean(dataset[,namei]))
#   sdVar = c(sdVar,sd(dataset[,namei]))            
#   names(meanVar)[length(meanVar)] = namei
#   names(sdVar)[length(sdVar)] = namei
#   dataset[,namei] =  (dataset[,namei] - meanVar[namei])/sdVar[namei]
#}

# Halton Draws 
preparedraws=function()
{
  d=1
  while(d<(length(normaldraws)+1))
  {
    draws1[,normaldraws[d]]<<- qnorm(draws1[,normaldraws[d]])
    d=d+1
  }
}

Ndraws=500      # set number of draws 
dim1=2    # define number of random parameters in the count model
dim2=2    # define number of random parameters in the ordered model
dimensions=dim1+dim2    # define number of random parameters in the model
# generate draws (using Halton)
draws1=as.matrix(halton(Ndraws*N,dimensions))

# assign names to individual sets of draws - need one entry per dimension
colnames(draws1)=c("HRbeta1","HRbeta2","HRbeta3","HRbeta4")
# define whether any draws should be transformed to Normals, which is also needed for e.g. lognormals (leave empty if not)
normaldraws=c("HRbeta1","HRbeta2","HRbeta3","HRbeta4")

# preparing draws for estimation - this may take a while
preparedraws()

# fixing parameters across grouped observations i.e. grouped random parameters
# comment out if there is no panel
block = length(unique(dataset[,'ID']))
ngroup = length(unique(dataset[,'Group']))
for (i in 1:Ndraws){
  tempInd = ((i-1)*block*ngroup) + (1:block)
  for (ii in 2:ngroup){
    draws1[tempInd+(ii-1)*block,] = draws1[tempInd,]
  }
}

## data preparation
# separating the variables with fixed parameters in the count model
dataF1 =  as.matrix(data.frame(1,log(dataset$Length)))
# separating the variables with fixed parameters in the ordered model
dataF2 =  as.matrix(data.frame(dataset$SEAL,dataset$radius,dataset$FuncDum))

# separating the variables with random parameters 
dataR1 = as.matrix(data.frame(log(dataset$AADT),dataset$LWIDTH))
# separating the variables with random parameters 
dataR2 = as.matrix(data.frame(dataset$curv,dataset$NumLane))

dataR11=NULL
dataR22=NULL
Dvar2 = NULL
for(i in 1:Ndraws){
  dataR11=rbind(dataR11,dataR1)
  dataR22=rbind(dataR22,dataR2)
  Dvar2 = c(Dvar2,DVar)
}

draws1 = draws1[,1:dimensions]

# placeholder for probabilities
prob = matrix(NA,nrow = N*Ndraws, ncol = S)

# Likelihood function
LL <- function(params){  
  disp <- params[1] # dispersion parameter of the NB distribution
  Fbeta1 <- params[2:3] # Fixed parameters in the mean Function of count model
  MRbeta1 <- params[4:5]  # Mean of Random parameters in the mean function of count model
  SDRbeta1 <- params[6:7]  # Std of Random parameters in the mean function of count model
  
  Fbeta2 <- params[8:10] # Fixed parameters in the mean Function of ordered model (excluding the constant)
  MRbeta2 <- params[11:12]  # Mean of Random parameters in the mean function of ordered model
  SDRbeta2 <- params[13:14]  # Std of Random parameters in the mean function of ordered model
  cutpoint1 = params[15] #first threshold of ordered model
  cutpoint2 = params[16] #second threshold of ordered model

  #############################################################################################################################
  # count model #
  # vector of indipendent variables with fixed parameters
  offset1 = rep.int(dataF1%*%as.matrix(Fbeta1,ncol=1),Ndraws)
  # simulating random parameters from their means and standard deviation
  beta1 = t( t(draws1[,1:dim1])*SDRbeta1 + MRbeta1 )
  # constructing the mean function
  mu <- exp(offset1+rowSums(dataR11*beta1))
  
  #############################################################################################################################
  # ordered model #
  
  # vector of indipendent variables with fixed parameters
  offset2 = rep.int(dataF2%*%as.matrix(Fbeta2,ncol=1),Ndraws)
  # simulating random parameters from their means and standard deviation
  beta2 = t( t(draws1[,dim1+1:dim2])*SDRbeta2 + MRbeta2)
  # constructing the mean function
  MEQ = offset2 + rowSums(dataR22*beta2) + offset1#observed correlation term
  
  prob[,1] <- plogis(cutpoint1 - MEQ)
  prob[,2] <- plogis(cutpoint2 - MEQ) - prob[,1]
  prob[,3] <- 1 - prob[,1] - prob[,2] 
  
  # cumulative probability functions for normal distribution (ordered probit)
  # prob[,1] <- pnorm(cutpoint1 - MEQ)
  # prob[,2] <- pnorm(cutpoint2 - MEQ) - prob[,1]
  # prob[,3] <- 1 - prob[,1] - prob[,2]
  
  # Lgarithim of probability for each category
  PR1 <- log(rowMeans(matrix(prob[,1], ncol = Ndraws)))
  PR2 <- log(rowMeans(matrix(prob[,2], ncol = Ndraws)))
  PR3 <- log(rowMeans(matrix(prob[,3], ncol = Ndraws)))
  
  # Logarithm of the probability for all fractions
  LNPR <- w[,1]*PR1+w[,2]*PR2+w[,3]*PR3 
  
  # simulated loglikelihood for the joint fractional split model
  CLIK = sum(log(rowMeans(matrix(dnbinom(CrashM,size=disp,mu=mu,log = F), ncol = Ndraws))))
  OLIK=  sum(LNPR)
  loglik <-  CLIK + OLIK

  return(loglik)
}

# initial values for optimization
init <- c(1,#dispersion parameter
          0.01, 0.01,#fixed parameters in count model 
          0.01, 0.01,#rMean of andom parameters in count model
          0.01, 0.01,#SD of random parameters in count model
          0.01, 0.01, 0.01,#fixed parameters in ordered model
          0.01, 0.01,#mean of random parameters in ordered model
          0.01, 0.01,#SD of random parameters in ordered model
          0.01, 10)#cut-off points

# optimization (maximization of likelihood function)
fit1 <- maxLik(LL,start=init,method="BFGS")

summary(fit1)






