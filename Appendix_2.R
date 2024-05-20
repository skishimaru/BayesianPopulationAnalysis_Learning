# Bayesian Population Analysis using WinBUGS
# Appendix 2: Two Further Useful Multistate Capture-Recapture Models

library(R2jags) #to run JAGS
library(shinystan) #to run shiny stan
library(tidyverse) #to utilize pipe operators

# 2.1 Estimation of Age-Specific Survival Probabilities
# 2.1.2 Generation of Simulated Data
#-------------------------------------------------------------------------------
#Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
phi.juv <- 0.3
phi.ad <- 0.65
p <- 0.5
n.occasions <- 6
n.states <- 3
n.obs <- 3
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(200, n.occasions) #Juveniles
marked[,2] <- rep(30, n.occasions) #Adults
marked[,3] <- rep(0, n.occasions) #Dead individuals

#Define matrices with survival, transition and recapture probabilities
#These are 4-dimensional matrices, with
#Dimension 1: state of departure
#Dimension 2: state of arrival
#Dimension 3: individual
#Dimension 4: time
# 1. State process matrix
totrel <- sum(marked)*(n.occasions-1)
PSI.STATE <- array(NA, dim=c(n.states, n.states, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.STATE[,,i,t] <- matrix(c(
      0, phi.juv, 1-phi.juv,
      0, phi.ad, 1-phi.ad,
      0, 0, 1 ), nrow = n.states, byrow = TRUE)
  } 
} 

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.OBS[,,i,t] <- matrix(c(
      0, 0, 1,
      0, p, 1-p,
      0, 0, 1 ), nrow = n.states, byrow = TRUE)
  } 
} 

#Execute simulation function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked)
CH <- sim$CH

#Compute vector with occasion of first capture
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

#Recode CH matrix: note, a 0 is not allowed in WinBUGS!
#1 = seen alive as juvenile, 2 = seen alive as adult, 3 = not seen
rCH <- CH #Recoded CH
rCH[rCH==0] <- 3
#-------------------------------------------------------------------------------

# 2.1.3 Analysis of the Model
#-------------------------------------------------------------------------------
#Specify model in JAGS language
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  #---------------------------------------
  #Parameters:
  # phi.juv: juvenile survival probability
  # phi.ad: adult survival probability
  # p: recapture probability
  #---------------------------------------
  #States (S):
  # 1 alive as juvenile
  # 2 alive as adult
  # 3 dead
  #Observations (O):
  # 1 seen as juvenile
  # 2 seen as adult
  # 3 not seen
  #---------------------------------------
  
  #Priors and constraints
  for (t in 1:(n.occasions-1)){
    phi.juv[t] <- mean.phijuv
    phi.ad[t] <- mean.phiad
    p[t] <- mean.p
  }
  mean.phijuv ~ dunif(0, 1)
  mean.phiad ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  
  #Define state-transition and observation matrices
  for (i in 1:nind){
    #Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
      ps[1,i,t,1] <- 0
      ps[1,i,t,2] <- phi.juv[t]
      ps[1,i,t,3] <- 1-phi.juv[t]
      ps[2,i,t,1] <- 0
      ps[2,i,t,2] <- phi.ad[t]
      ps[2,i,t,3] <- 1-phi.ad[t]
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- 1
      #Define probabilities of O(t) given S(t)
      po[1,i,t,1] <- 0
      po[1,i,t,2] <- 0
      po[1,i,t,3] <- 1
      po[2,i,t,1] <- 0
      po[2,i,t,2] <- p[t]
      po[2,i,t,3] <- 1-p[t]
      po[3,i,t,1] <- 0
      po[3,i,t,2] <- 0
      po[3,i,t,3] <- 1
    } 
  } 
  
  #State-space model likelihood
  for (i in 1:nind){
    z[i,f[i]] <- Y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
      #State equation: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
      #Observation equation: draw O(t) given S(t)
      Y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } 
  } 
}

#Bundle data
jags.data <- list(Y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1])

#Initial values
inits <- function(){list(mean.phijuv = runif(1, 0, 1), mean.phiad = runif(1, 0, 1), mean.p = runif(1, 0, 1), z = ch.init(rCH, f))}

#Parameters monitored
params <- c("mean.phijuv", "mean.phiad", "mean.p")

#MCMC settings
ni <- 2000
nt <- 3
nb <- 1000
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

print(out, digits = 3)

k<-mcmcplots::as.mcmc.rjags(out)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan
#-------------------------------------------------------------------------------

# 2.2 Accounting for Immediate Trap Response
# 2.2.2 Generation of Simulated Date
#-------------------------------------------------------------------------------
#Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
phi <- 0.55
pss <- 0.75
pns <- 0.3
n.occasions <- 10
n.states <- 3
n.obs <- 2
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(100, n.occasions) #Alive, seen
marked[,2] <- rep(0, n.occasions) #Alive, not seen
marked[,3] <- rep(0, n.occasions) #Dead

#Define matrices with survival, transition and recapture probabilities
#These are 4-dimensional matrices, with
#Dimension 1: state of departure
#Dimension 2: state of arrival
#Dimension 3: individual
#Dimension 4: time
# 1. State process matrix
totrel <- sum(marked)*(n.occasions-1)
PSI.STATE <- array(NA, dim=c(n.states, n.states, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.STATE[,,i,t] <- matrix(c(
      phi*pss, phi*(1-pss), 1-phi,
      phi*pns, phi*(1-pns), 1-phi,
      0, 0, 1 ), nrow = n.states, byrow = TRUE)
  } 
} 

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.OBS[,,i,t] <- matrix(c(
      1, 0,
      0, 1,
      0, 1 ), nrow = n.states, byrow = TRUE)
  } 
} 

#Execute simulation function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked, unobservable = 2)
CH <- sim$CH

#Compute vector with occasion of first capture
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

#Recode CH matrix: note, a 0 is not allowed in WinBUGS!
#1 = seen, 2 = not seen
rCH <- CH #Recoded CH
rCH[rCH==0] <- 2
#-------------------------------------------------------------------------------

# 2.2.3 Analysis of the Model 
#-------------------------------------------------------------------------------
#Specify model in JAGS language
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  #----------------------------------------------------------
  #Parameters:
  # phi: survival probability
  # pss: recapture probability at t, given captured at t−1
  # pns: recapture probability at t, given not captured at t−1
  #----------------------------------------------------------
  #States (S):
  # 1 alive, seen at t−1
  # 2 alive, not seen at t−1
  # 3 dead
  #Observations (O):
  # 1 seen
  # 2 not seen
  #----------------------------------------------------------
  
  #Priors and constraints
  for (t in 1:(n.occasions-1)){
    phi[t] <- mean.phi
    pss[t] <- mean.pss
    pns[t] <- mean.pns
  }
  mean.phi ~ dunif(0, 1)
  mean.pss ~ dunif(0, 1)
  mean.pns ~ dunif(0, 1)
  
  #Define state-transition and observation matrices
  for (i in 1:nind){
    #Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
      ps[1,i,t,1] <- phi[t]*pss[t]
      ps[1,i,t,2] <- phi[t]*(1-pss[t])
      ps[1,i,t,3] <- 1-phi[t]
      ps[2,i,t,1] <- phi[t]*pns[t]
      ps[2,i,t,2] <- phi[t]*(1-pns[t])
      ps[2,i,t,3] <- 1-phi[t]
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- 1
      #Define probabilities of O(t) given S(t)
      po[1,i,t,1] <- 1
      po[1,i,t,2] <- 0
      po[2,i,t,1] <- 0
      po[2,i,t,2] <- 1
      po[3,i,t,1] <- 0
      po[3,i,t,2] <- 1
    } 
  } 
  
  #State-space model likelihood
  for (i in 1:nind){
    z[i,f[i]] <- Y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
      #State equation: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t−1], i, t−1,])
      #Observation equation: draw O(t) given S(t)
      Y[i,t] ~ dcat(po[z[i,t], i, t−1,])
    } 
  } 
}

#Bundle data
jags.data <- list(Y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1], z = known.state.ms(rCH, 2))

#Initial values
inits <- function(){list(mean.phi = runif(1, 0, 1), mean.pss =runif(1, 0, 1), mean.pns = runif(1, 0, 1), z = ms.init.z(rCH, f))}

#Parameters monitored
params <- c("mean.phi", "mean.pss", "mean.pns")

#MCMC settings
ni <- 20000
nt <- 6
nb <- 10000
nc <- 3
  
#Call JAGS from R
immtrap <- jags(data  = jags.data,
            inits = inits,
            parameters.to.save = params,
            model.file = jags.model.txt,
            n.chains = nc,
            n.thin= nt,
            n.iter = ni,
            n.burnin = nb)

print(immtrap, digits = 3)

k<-mcmcplots::as.mcmc.rjags(immtrap)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan
#-------------------------------------------------------------------------------