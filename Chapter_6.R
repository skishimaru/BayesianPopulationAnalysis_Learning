# Bayesian Population Analysis using WinBUGS
# Chapter 6: Estimation of the Size of a Closed Population from Capture-Recapture Data

library(R2jags) #to run JAGS
library(shinystan) #to run shiny stan
library(tidyverse) #to utilize pipe operators

#NOTE: All example models in this chapter assume closure, which means each species is part of the sampled community at a given point and must be available for detection during all replicate surveys

# 6.2 Generation and Analysis of Simulated Data with Data Augmentation
# 6.2.1 Introduction to Data Augmentation for the Simplest Case: Model M0
#-------------------------------------------------------------------------------
#Detection probability is constant over 2 dimensions of the detection parameter (individuals and time periods)
#Arguments: Population size (N), Detection probability (p), and Number of sampling occasions (T)
#For this model T is assumed constant for all individuals

#Define function to simulate data under M0
data.fn <- function(N= 100, p= 0.5, T= 3){
  yfull <- yobs <- array(NA, dim= c(N, T))
  for (j in 1:T){
    yfull[,j] <- rbinom(n= N, size= 1, prob= p)
  }
  ever.detected <- apply(yfull, 1, max)
  C <- sum(ever.detected)
  yobs <- yfull[ever.detected== 1,]
  cat(C, "out of", N, "animals present were detected")
  return(list(N=N, p= p, C= C, T= T, yfull= yfull, yobs= yobs))
}

data <- data.fn()

str(data) #shows us the list of R objects
#yfull= full capture history matrix of all N animals, including those which were never captured and have an all zero capture history
#yobs= observed data 

#Challenge= the dimension of the parameter for the abundance vector (N) may change at every iteration of the MCMC algorithm 
#Solution, Data Augmentation= Converts the closed population model into an occupancy model and turns the problem of estimating abundance (N) into estimating occupancy
#Achieves this by 1. adding an arbitrary number of 0s to the data set (a large number of potential unobserved individuals with no encounter history) and 2. analyzing a reparameterized version of the original model

#Data Augmentation Initial Example
#Augment data set by 150 potential individuals
nz <- 150
yaug <- rbind(data$yobs, array(0, dim= c(nz, data$T)))

#Specify model in JAGS
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors
  omega ~ dunif(0,1)
  p ~ dunif(0,1)
  
  #Likelihood
  for (i in 1:M){
    z[i] ~ dbern(omega) #inclusion indicators
    for (j in 1:T){
      yaug[i,j] ~ dbern(p.eff[i,j])
      p.eff[i,j] <- z[i]*p # can only be detected if z= 1
    }
  }
  
  #Derived Quantities 
  N <- sum(z[])
}

jags.data <- list(yaug= yaug, M= row(yaug), T= ncol(yaug)) #Bundle data
inits <- function() list(z= rep(1, nrow(yaug)), p= runif(1, 0, 1)) #Initial Values
params <- c("N", "p", "omega") #Parameters monitored

#MCMC Settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

#Call JAGS from R
out <- jags(data  = jags.data,
            inits = inits,
            parameters.to.save = params,
            model.file = jags.model.txt,
            n.chains = nc,
            n.thin= nt,
            n.iter = ni,
            n.burnin = nb)

print(out, dig = 3) #summarize Posteriors

hist(out$BUGSoutput$sims.list$deviance, nclass= 50, col= "gray", main= "Figure 6.1", xlab= "Population size", las= 1, xlim= c(80,150))
abline(v= data$C, lwd= 3)

k<-mcmcplots::as.mcmc.rjags(out)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#Data Augmentation Example: Proof Adding 0s Does Not Effect Estimates of Detection Probability (p) and Population Size (N)
#This example will create figure 6.2 which plots the posterior distributions of population size N under different degrees of data augmentation (black line= observed number of animals, blue line= posterior mean)
#1. Augment data set to 5 potential individuals
nz <- 5
yaug <- rbind(data$yobs, array(0, dim= c(nz, data$T)))
jags.data <- list(yaug= yaug, M= dim(yaug)[1], T= dim(yaug)[2])
out5 <- jags(data  = jags.data,
            inits = inits,
            parameters.to.save = params,
            model.file = jags.model.txt,
            n.chains = nc,
            n.thin= nt,
            n.iter = ni,
            n.burnin = nb) 
print(out5, dig= 3)
par(mfrow= c(3,1))
hist(out5$BUGSoutput$sims.list$deviance, nclass= 30, col= "gray", main= "Augmentation by 5", xlab= "Population size", las= 1, xlim= c(80,140))
abline(v= data$C, col= "black", lwd= 5)
abline(v= mean(out5$BUGSoutput$sims.list$deviance), col= "blue", lwd= 3)

#2. Augment data set to 150 potential individuals
nz <- 150
yaug <- rbind(data$yobs, array(0, dim= c(nz, data$T)))
jags.data <- list(yaug= yaug, M= dim(yaug)[1], T= dim(yaug)[2])
out150 <- jags(data  = jags.data,
             inits = inits,
             parameters.to.save = params,
             model.file = jags.model.txt,
             n.chains = nc,
             n.thin= nt,
             n.iter = ni,
             n.burnin = nb) 
print(out150, dig= 3)
hist(out150$BUGSoutput$sims.list$deviance, nclass= 30, col= "gray", main= "Augmentation by 150", xlab= "Population size", las= 1, xlim= c(80,140))
abline(v= data$C, col= "black", lwd= 5)
abline(v= mean(out150$BUGSoutput$sims.list$deviance), col= "blue", lwd= 3)

#3. Augment data set to 150 potential individuals
nz <- 1500
yaug <- rbind(data$yobs, array(0, dim= c(nz, data$T)))
jags.data <- list(yaug= yaug, M= dim(yaug)[1], T= dim(yaug)[2])
out1500 <- jags(data  = jags.data,
               inits = inits,
               parameters.to.save = params,
               model.file = jags.model.txt,
               n.chains = nc,
               n.thin= nt,
               n.iter = ni,
               n.burnin = nb) 
print(out1500, dig= 3)
hist(out1500$BUGSoutput$sims.list$deviance, nclass= 30, col= "gray", main= "Augmentation by 1500", xlab= "Population size", las= 1, xlim= c(80,140))
abline(v= data$C, col= "black", lwd= 5)
abline(v= mean(out1500$BUGSoutput$sims.list$deviance), col= "blue", lwd= 3)
mtext("Figure 6.2", side= 3, line= -2, outer= T) #adding main title to multiplot

#This exercise shows that data augmentation does not have an effect on estimates but effects efficiency. Data augmentation is more costly in terms of computation time
#NOTE, Diagnosing not enough of data augmentation: Looks at the posterior distribution of N. If it is truncated on the right by your choice of M you need to repeat the analysis with more 0s added
#-------------------------------------------------------------------------------

# 6.2.2 Time Effects: Model Mt
#-------------------------------------------------------------------------------
#Assumes that detection probability (p) varies by occasion because of weather conditions, different traps, different detection devices used, etc.
#The simulations in this example use different detection methods (ex. traps vs observers) during a single occasion and can be treated exactly as time effects in capture-recapture modeling

#Define function to simulate data under Mt
data.fn <- function(N= 100, mean.p= 0.5, T= 3, time.eff= runif(T, -2, 2)){
  yfull <- yobs <- array(NA, dim= c(N,T))
  p.vec <- array(NA, dim= T)
  for (j in 1:T){
    p <- plogis(log(mean.p/(1-mean.p)) + time.eff[j])
    yfull[,j] <- rbinom(n= N, size= 1, prob= p)
    p.vec[j] <- p
  }
  ever.detected <- apply(yfull, 1, max)
  C <- sum(ever.detected)
  yobs <- yfull[ever.detected== 1,]
  cat(C, "out of", N, "animals present were detected.")
  return(list(N= N, p.vec= p.vec, C= C, T= T, yfull= yfull, yobs= yobs))
}

data <- data.fn() #Created a data set

#Analysis in JAGS
#Augment data set
nz <- 150
yaug <- rbind(data$yobs, array(0, dim= c(nz, data$T)))

#Specify model in JAGS
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors
  omega ~ dunif(0,1)
  for (i in 1:T){
    p[i] ~ dunif(0,1)
  }
  
  #Likelihood
  for (i in 1:M){
    z[i] ~ dbern(omega) #inclusion indicators
    for (j in 1:T){
      yaug[i,j] ~ dbern(p.eff[i,j])
      p.eff[i,j] <- z[i]*p # can only be detected if z= 1
    }
  }
  
  #Derived Quantities 
  N <- sum(z[])
}

jags.data <- list(yaug= yaug, M= row(yaug), T= ncol(yaug)) #Bundle data
inits <- function() list(z= rep(1, nrow(yaug)), p= runif(data$T, 0, 1)) #Initial Values
params <- c("N", "p", "omega") #Parameters monitored

#MCMC Settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

#Call JAGS from R
out <- jags(data  = jags.data,
            inits = inits,
            parameters.to.save = params,
            model.file = jags.model.txt,
            n.chains = nc,
            n.thin= nt,
            n.iter = ni,
            n.burnin = nb)

print(out, dig = 3) #summarize Posteriors

hist(out$BUGSoutput$sims.list$deviance, nclass= 50, col= "gray", main= "", xlab= "Population size", las= 1, xlim= c(80,150))
abline(v= data$C, lwd= 3)

k<-mcmcplots::as.mcmc.rjags(out)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan
#-------------------------------------------------------------------------------

# 6.2.3 Behavioral or Memory Effects: Model Mb
#-------------------------------------------------------------------------------
#Modeling immediate trap response: when an individual is captured it has a different detection probability only on the nest following occasion but not after unless if its captured again

#Define function to simulate data under Mb
#in this function we simulate trap response on the probability scale. We denote the capture probability as c or p depending on whether an individual has or has not been captured during the preeceding occasion
data.fn <- function(N= 200, T= 5, p= 0.3, c= 0.4){
  yfull <- yobs <- array(NA, dim= c(N,T))
  p.eff <- array(NA, dim= N)
  #First Capture Occasion
  yfull[,1] <- rbinom(n= N, size= 1, prob= p)
  #Later Capture Occasions
  for (j in 2:T){
    p.eff <- (1 - yfull[,(j-1)])*p + yfull[,(j-1)]*c
    yfull[,j] <- rbinom(n=N, size= 1, prob= p.eff)
  }
  ever.detected <- apply(yfull, 1, max)
  C <- sum(ever.detected)
  yobs <- yfull[ever.detected== 1,]
  cat(C, "out of", N, "animals present were detected.")
  return(list(N= N, p= p, C= C, T= T, yfull= yfull, yobs= yobs))
}

data <- data.fn(N= 200) #Created a data set with trap happiness (p<c)

#Analysis in JAGS
#Augment data set
nz <- 150
yaug <- rbind(data$yobs, array(0, dim= c(nz, data$T)))

#Specify model in JAGS
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors
  omega ~ dunif(0,1)
  p ~ dunif(0,1) #can prob when not caught t-1
  c ~ dunif(0,1) #cap prob when caught at t-1
  
  #Likelihood
  for (i in 1:M){
    z[i] ~ dbern(omega) #inclusion indicators
    #First Occasion
    yaug[i,1] ~ dbern(p.eff[i,1])
    p.eff[i,1] <- z[i]*p
    #All subsequent occasions
    for (j in 2:T){
      yaug[i,j] ~dbern(p.eff[i,j])
      p.eff[i,j] <- z[i]*((1 - yaug[i, (j-1)])*p + yaug[i, (j-1)]*c)
    }
  }
  #Derived Quantities 
  N <- sum(z[])
  trap.response <- c - p
}

jags.data <- list(yaug= yaug, M= row(yaug), T= ncol(yaug)) #Bundle data
inits <- function() list(z= rep(1, nrow(yaug)), p= runif(1, 0, 1)) #Initial Values
params <- c("N", "p", "C", "trap.response", "omega") #Parameters monitored

#MCMC Settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

#Call JAGS from R
out <- jags(data  = jags.data,
            inits = inits,
            parameters.to.save = params,
            model.file = jags.model.txt,
            n.chains = nc,
            n.thin= nt,
            n.iter = ni,
            n.burnin = nb)

print(out, dig = 3) #summarize Posteriors

hist(out$BUGSoutput$sims.list$deviance, nclass= 50, col= "gray", main= "", xlab= "Population size", las= 1, xlim= c(80,150))
abline(v= data$C, lwd= 3)

k<-mcmcplots::as.mcmc.rjags(out)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan
#-------------------------------------------------------------------------------

# 6.2.4 Individual (Random) Effects: The Heterogeneity Model Mh
#-------------------------------------------------------------------------------
#Assumes each individual has its own detection probability and there are individuals latent (random) effects in detection probability
#In the model individual heterogeneity is modeled as random noise around some mean on a logit-transformed scale. The model for noise is the normal distribution, so it is called a logistic normal
#In this model we aggregate capture histories to capture frequencies
#mean.p= average detection probability of individuals in the population, sd= standard deviation of the normal distribution which describes the heterogeneity in the individual detection probability on the logit scale

#Define function to simulate data under Mh
data.fn <- function(N= 100, mean.p= 0.4, T= 5, sd= 1){
  yfull <- yobs <- array(NA, dim= c(N,T))
  mean.lp <- log(mean.p/(1-mean.p))
  p.vec <- plogis(mean.lp + rnorm(N, 0, sd))
  for (i in 1:N){
    yfull[i,] <- rbinom(n= T, size= 1, prob= p.vec[i])
  }
  ever.detected <- apply(yfull, 1, max)
  C <- sum(ever.detected)
  yobs <- yfull[ever.detected== 1,]
  cat(C, "out of", N, "animals present were detected.")
  hist(p.vec, xlim= c(0,1), nclass= 20, col= "gray", main= "Figure 6.3", xlab= "Detection Probability", las= 1)
  return(list(N= N, p.vec= p.vec, mean.lp= mean.lp, C= C, T= T, yfull= yfull, yobs= yobs))
}

data <- data.fn() #Created a data set with trap happiness (p<c)

#Analysis in JAGS
#Aggregate capture histories and augment data set
y <- sort(apply(data$yobs, 1, sum), decreasing= TRUE)
nz <- 300
yaug <- c(y, rep(0, nz))

#Specify model in JAGS
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors
  omega ~ dunif(0,1)
  mean.lp <- logit(mean.p)
  mean.p <- dunif(0,1)
  tau <- 1/(sd*sd)
  sd ~ dunif(0,5)
  
  #Likelihood
  for (i in 1:M){
    z[i] ~ dbern(omega) #inclusion indicators
    logit(p[i]) <- eps[i]
    eps[i] ~ dnorm(mean.lp, tau) I(-16,16) #I() changes the class of an object to be treated "as is", inhibit interpretation
    p.eff[i] <- z[i]*p[i]
    y[i] ~dbin(p.eff[i], T)
  }
  #Derived Quantities 
  N <- sum(z[])
}

jags.data <- list(yaug= yaug, M= length(yaug), T= ncol(data$yobs)) #Bundle data
inits <- function() list(z= rep(1, length(yaug)), sd= runif(1, 0.1, 0.9)) #Initial Values
params <- c("N", "mean.p", "sd", "omega") #Parameters monitored

#MCMC Settings
ni <- 25000
nt <- 2
nb <- 5000
nc <- 3

#Call JAGS from R
out <- jags(data  = jags.data,
            inits = inits,
            parameters.to.save = params,
            model.file = jags.model.txt,
            n.chains = nc,
            n.thin= nt,
            n.iter = ni,
            n.burnin = nb)

print(out, dig = 3) #summarize Posteriors

hist(out$BUGSoutput$sims.list$deviance, nclass= 50, col= "gray", main= "Figure 6.4", xlab= "Population size", las= 1, xlim= c(80,250))
abline(v= data$C, lwd= 3)

k<-mcmcplots::as.mcmc.rjags(out)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan
#We see that precision for N is lower compared to other models. Estimation in individual effects models tend to be more challenging than models without random effects and we need to run them for a longer time to get an adequate posterior sample.
#-------------------------------------------------------------------------------

# 6.2.5 Combined Effects (time + random): Model Mth
#-------------------------------------------------------------------------------
#Example that shows the additive fixed effect of time and random (logistic normal) effects of individuals

#Define function to simulate data under Mh
data.fn <- function(N= 100, T= 5, mean.p= 0.4, time.effects= runif(T, -1, 1), sd= 1){ 
  yfull <- yobs <- p <- array(NA, dim= c(N,T))
  mean.lp <- log(mean.p/(1-mean.p)) #mean p on logit scale
  eps <- rnorm(N, 0, sd) #Individual effects
  for (j in 1:T){
    pp <- p[,j] <- plogis(mean.lp + time.effects[j] + eps)
    yfull[,j] <- rbinom(n= N, size= 1, prob= pp)
  }
  ever.detected <- apply(yfull, 1, max)
  C <- sum(ever.detected)
  yobs <- yfull[ever.detected== 1,]
  cat(C, "out of", N, "animals present were detected.")
  cat("Mean p per occasion:", round(apply(p, 2, mean), 2), "\n")
  par(mfrow= c(2,1))
  plot(plogis(mean.lp + time.effects), xlab= "Occasion", type= "b", main= "Approximate m=Mean p at Each Occasion", ylim= c(0,1))
  hist(plogis(mean.lp + eps), xlim= c(0,1), col= "gray", main= "Approximate Distribution of p at Average Occasion")
  return(list(N= N, mean.lp= mean.lp, time.effects= time.effects, C= C, T= T, yfull= yfull, yobs= yobs))
}

data <- data.fn() #Created a data set with trap happiness (p<c)

#Analysis in JAGS
#Aggregate capture histories and augment data set
nz <- 300
yaug <- rbind(data$yobs, array(0, dim= c(nz, data$T)))

#Specify model in JAGS
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors
  omega ~ dunif(0,1)
  for (j in 1:T){
    mean.lp[j] <- log(mean.p[j]/(1 - mean.p[j])) #Define logit
    mean.p[j] ~ dunif(0,1)
  }
  tau <- 1/(sd*sd)
  sd ~ dunif(0,5)
  
  #Likelihood
  for (i in 1:M){
    z[i] ~ dbern(omega) #inclusion indicators
    eps[i] ~ dnorm(0, tau) I(-16,16) #I() changes the class of an object to be treated "as is", inhibit interpretation
    for (j in 1:T){
      lp[i,j] <- mean.lp[j] + eps[i]
      p[i,j] <- 1/(1 + exp(-lp[i,j])) #Define logit
      p.eff[i,j] <- z[i]*p[i,j]
      y[i,j] ~ dbern(p.eff[i,j])
    }
  }
  #Derived Quantities 
  N <- sum(z[])
}

jags.data <- list(yaug= yaug, M= nrow(yaug), T= ncol(yaug)) #Bundle data
inits <- function() list(z= rep(1, length(yaug)), sd= runif(1, 0.1, 0.9)) #Initial Values
params <- c("N", "mean.p", "mean.lp", "sd", "omega") #Parameters monitored

#MCMC Settings
ni <- 25000
nt <- 2
nb <- 5000
nc <- 3

#Call JAGS from R
out <- jags(data  = jags.data,
            inits = inits,
            parameters.to.save = params,
            model.file = jags.model.txt,
            n.chains = nc,
            n.thin= nt,
            n.iter = ni,
            n.burnin = nb)

print(out, dig = 3) #summarize Posteriors

hist(out$BUGSoutput$sims.list$deviance, nclass= 50, col= "gray", main= "Figure 6.5", xlab= "Population size", las= 1, xlim= c(80,200))
abline(v= data$C, lwd= 3)

k<-mcmcplots::as.mcmc.rjags(out)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan
#-------------------------------------------------------------------------------

# 6.3 Analysis of a Real Data Set: Model Mtbh for Species Richness Estimation
#-------------------------------------------------------------------------------
#Uses point count bird count data from Czech republic. An East-West transect spanning the whole country and conducted a 5 minute point count every 500m from 2004-2005. This process was repeated 5 times within a breeding season
#146 species were detected but this analysis will focus on the Wryneck (point count number= 610)
#Analysis is reduced to simple detection/nondetection data

#Read in data
p610 <- read.table("/Users/shelbie.ishimaru/Documents/GitHub/BayesianPopulationAnalysis_Learning/p610.csv", header= TRUE)
y <- p610[,5:9] #grab counts
y[y > 1] <- 1 #counts to det-nondetections
C <- sum(apply(y, 1, max)) ; print(C) #number of observed species
table(apple(y, 1, sum)) #capture frequencies
#NOTE: This data contains detection history of the 31 species detected at that particular point, but also all species detected anywhere along the entire transect. So the data is naturally augmented. We will not add additional zeros because 115 zeros corresponding to species not detected at point 610 but somewhere else in the transect is enough.

#Example that shows the additive fixed effect of time, behavior, and random (logistic normal) effects of individuals
#Specify model in JAGS
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors
  omega ~ dunif(0,1)
  for (j in 1:T){
    alpha[j] <- log(mean.p[j]/(1 - mean.p[j])) #Define logit
    mean.p[j] ~ dunif(0,1) #Detect intercepts
  }
  gamma ~ dnorm(0, 0.01)
  tau <- 1/(sd*sd)
  sd ~ dunif(0,3)
  
  #Likelihood
  for (i in 1:M){
    z[i] ~ dbern(omega) #inclusion indicators
    eps[i] ~ dnorm(0, tau) I(-16,16) #I() changes the class of an object to be treated "as is", inhibit interpretation
    #First occasion: no term for recapture (gamma)
    y[i,1] ~ dbern(p.eff[i,1])
    p.eff[i,1] <- z[i]*p[i,1]
    p[i,1] <- 1/(1 + exp(-lp[i,1])) #Define logit
    lp[i,1] <- alpha[1] + eps[i]
    #All subsequent occasions: includes recapture term (gamma)
    for (j in 2:T){
      y[i,j] ~ dbern(p.eff[i,j])
      p.eff[i,j] <- z[i]*p[i,j]
      p[i,j] <- 1/(1 + exp(-lp[i,j])) #Define logit
      lp[i,j] <- alpha[j] + eps[i] + gamma * y[i,(j - 1)]
    }
  }
  #Derived Quantities 
  N <- sum(z[])
}

jags.data <- list(y= as.matrix(y), M= nrow(y), T= ncol(y)) #Bundle data
inits <- function() list(z= rep(1, nrow(y)), sd= runif(1, 0.1, 0.9)) #Initial Values
params <- c("N", "mean.p", "gamma", "sd", "omega") #Parameters monitored

#MCMC Settings
ni <- 50000
nt <- 4
nb <- 10000
nc <- 3

#Call JAGS from R
out <- jags(data  = jags.data,
            inits = inits,
            parameters.to.save = params,
            model.file = jags.model.txt,
            n.chains = nc,
            n.thin= nt,
            n.iter = ni,
            n.burnin = nb)

print(out, dig = 3) #summarize Posteriors
#in this model we detected 31 species and estimate there was 42 species, so only 74% were detected
#For a short survey this outcome is reasonable because average detection of a species community ranges from 0.22 - 0.32

par(mfrow= c(1,2))
hist(out$BUGSoutput$sims.list$deviance, nclass= 35, col= "gray", main= "Figure 6.7", xlab= "Community Size", las= 1, xlim= c(30,100), freq= FALSE)
abline(v= C, col= "black", lwd= 3)

k<-mcmcplots::as.mcmc.rjags(out)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#Example: Shows how inference under model M0 often leads to a serious negative bias in the estimate of N where there is heterogeneity among species p
#Specify model in JAGS
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors
  omega ~ dunif(0,1)
  p ~ dunif(0,1)
  
  #Likelihood
  for (i in 1:M){
    z[i] ~ dbern(omega) #inclusion indicators
    for (j in 1:T){
      y[i,j] ~ dbern(p.eff[i,j])
      p.eff[i,j] <- z[i]*p
    }
  }
  #Derived Quantities 
  N <- sum(z[])
}

inits <- function() list(z= rep(1, nrow(y))) #Initial Values
params <- c("N", "p", "omega") #Parameters monitored

#MCMC Settings
ni <- 50000
nt <- 4
nb <- 10000
nc <- 3

#Call JAGS from R
out0 <- jags(data  = jags.data,
            inits = inits,
            parameters.to.save = params,
            model.file = jags.model.txt,
            n.chains = nc,
            n.thin= nt,
            n.iter = ni,
            n.burnin = nb)

print(out0, dig = 3) #summarize Posteriors
#Under this model we estimated 38 instead of 42 species. This illustrates that ignoring individual heterogeneity in detection probability produces underestimates of and too short standard errors for population size (N)

k<-mcmcplots::as.mcmc.rjags(out)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan
#-------------------------------------------------------------------------------

# 6.4 Capture-Recapture Models with Individual Covariables: Model Mt+x
# 6.4.1 Individual Covariate Model for Species Richness Estimation
#-------------------------------------------------------------------------------
#Read in data
p610 <- read.table("/Users/shelbie.ishimaru/Documents/GitHub/BayesianPopulationAnalysis_Learning/p610.csv", header= TRUE)
y <- p610[,5:9] #grab counts
y[y > 1] <- 1 #counts to det-nondetections
ever.observed <- apply(y, 1, max)
wt <- p610$bm[ever.obsessed== 1] #body mass
yy <- as.matrix(y[ever.observed== 1,]) #detection histories
dimnames(yy) <- NULL

mlog <- mean(log(p610$bm^(1/3)))
sdlog <- sd(log(p610$bm^(1/3)))
hist(p610$bm^(1/3), xlim= c(0,30), nclass= 25, freq= FALSE, col= "gray")
lines(density(rlnorm(n= 10^6, meanlog= mlog, sdlog= sdlog)), col= "blue", lwd= 3)

#Augment both data sets
#since we disgarded the data on undetected species we now need to activey augment the data
nz= 150
yaug <- rbind(yy, array(0, dim= c(nz, ncol(yy))))
logwt3 <- c(log(wt^(1/3)), rep(NA, nz))

#Specify model in JAGS
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors
  omega ~ dunif(0,1)
  for (j in 1:T){
    alpha[j] <- log(mean.p[j]/(1 - mean.p[j])) #Define logit
    mean.p[j] ~ dunif(0,1) #Detect intercepts
  }
  beta ~ dnorm(0, 0.01)
  mu.size ~ dnorm(0, 0.01)
  gamma ~ dnorm(0, 0.01)
  tau.size <- 1/(sd.size*sd.size)
  sd.size ~ dunif(0,prior.sd.upper) #provide upper bound as data
  
  #Likelihood
  for (i in 1:M){
    z[i] ~ dbern(omega) #inclusion indicators
    size[i] ~ dnorm(mu.size, tau.size) I(-6,6)
    for (j in 1:T){
      y[i,j] ~ dbern(p.eff[i,j])
      p.eff[i,j] <- z[i]*p[i,j]
      p[i,j] <- 1/(1 + exp(-lp[i,j])) #Define logit
      lp[i,j] <- alpha[j] + beta * size[i]
    }
  }
  #Derived Quantities 
  N <- sum(z[])
}

jags.data <- list(y= yaug, size= logwt3 - mean(logwt3, na.rm= TRUE), M= nrow(yaug), prior.sd.upper= 3) #Bundle data
inits <- function() list(z= rep(1, nrow(yaug)), beta= runif(1, 0, 1), mu.size= rnorm(1, 0, 1)) #Initial Values
params <- c("N", "mean.p", "beta", "omega", "mu.size", "sd.size") #Parameters monitored

#MCMC Settings
ni <- 50000
nt <- 4
nb <- 10000
nc <- 3

#Call JAGS from R
outX <- jags(data  = jags.data,
            inits = inits,
            parameters.to.save = params,
            model.file = jags.model.txt,
            n.chains = nc,
            n.thin= nt,
            n.iter = ni,
            n.burnin = nb)

print(outX, dig = 3) #summarize Posteriors
#in this model we detected 31 species and estimate there was 42 species, so only 74% were detected
#For a short survey this outcome is reasonable because average detection of a species community ranges from 0.22 - 0.32

hist(outX$BUGSoutput$sims.list$deviance, breaks= 100, col= "gray", main= "", xlab= "Community Size", las= 1, xlim= c(30,100), freq= FALSE)
abline(v= 31, col= "black", lwd= 3)

k<-mcmcplots::as.mcmc.rjags(out)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

pred.wt <- seq(5, 2000, length.out= 100) #cov. vals for prediction
pred.wt.st <- log(pred.wt^(1/3))-mean(logwt3, na.rm= TRUE)
#Transform them in the same way as in the analysis
pred.p <- plogis(log(mean(outX$BUGSoutput$mean$mean.p)/(1-mean(outX$BUGSoutput$mean$mean.p)))+ outX$BUGSoutput$mean$beta*pred.wt.st)
#Compute predicted response 
plot(pred.wt, pred.p, type= "l", lwd= 3, col= "blue", las= 1, frame.plot= FALSE, ylim= c(0, 0.5))
#-------------------------------------------------------------------------------

# 6.4.2 Individual Covariate Model for Population Size Estimation
#-------------------------------------------------------------------------------
pinna <- read.table("/Users/shelbie.ishimaru/Documents/GitHub/BayesianPopulationAnalysis_Learning/p610.csv", header= TRUE)
y <- cbind(pinna$d1, pinna$d2)
size <- pinna$width
hist(size, col= "gray", nclass= 50, xlim= c(0,30), freq= FALSE)
lines(density(rnorm(10^6, mean= mean(size), sd= sd(size))), col= "blue", lwd= 3)

#Augment both data sets
nz= 150
yaug <- rbind(y, array(0, dim= c(nz, ncol(y))))
size <- c(size, rep(NA, nz))

#MCMC Settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

#Call JAGS from R
outXX <- jags(data  = jags.data,
             inits = inits,
             parameters.to.save = params,
             model.file = jags.model.txt,
             n.chains = nc,
             n.thin= nt,
             n.iter = ni,
             n.burnin = nb)

print(outXX, dig = 2) #summarize Posteriors
#in this model we estimate that 172 instead of 143 pen shells were available for detection along the survey transects, so 29 were missed by both teams of divers
#detection probability was higher for the first team as expected, and there was a positive relationship with shell width

k<-mcmcplots::as.mcmc.rjags(out)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#Plot posterior for N and prediction of p
par(mfrow= c(1,2), mar= c(4.5, 4, 2, 1))
hist(outXX$BUGSoutput$sims.list$deviance, breaks= 30, col= "gray", main= "Figure 6.9", xlab= "Population size", las= 1, xlim= c(143, 220), freq= FALSE)
abline(v= 143, col= "black", lwd= 3)
pred.size <- seq(0, 30, length.out= 1000) #cov. vals for prediction
pred.size.st <- pred.size - mean(size,na.rm= TRUE) #Transform them
pred.p <- plogis(log(mean(outXX$BUGSoutput$mean$mean.p)/(1-mean(outXX$BUGSoutput$mean$mean.p)))+ outXX$BUGSoutput$mean$beta*pred.size.st)
#Compute predicted response 
plot(pred.size, pred.p, type= "l", lwd= 3, col= "blue", las= 1, frame.plot= FALSE, ylim= c(0, 1), xlab= "Shell width (cm)", ylab= "Predicted detection probability")
#------------------------------------------------------------------------------- 