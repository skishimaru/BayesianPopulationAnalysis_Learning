# Bayesian Population Analysis using WinBUGS
# Chapter 9: Estimation of Survival and Movement from Capture-Recapture Data Using Multistate Models

library(R2jags) #to run JAGS
library(shinystan) #to run shiny stan
library(tidyverse) #to utilize pipe operators

# 9.2 Estimation of Movement between Two Sites
# 9.2.2 Generation of Simulated Data
#-------------------------------------------------------------------------------
#Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
phiA <- 0.8
phiB <- 0.7
psiAB <- 0.3
psiBA <- 0.5
pA <- 0.7
pB <- 0.4
n.occasions <- 6
n.states <- 3
n.obs <- 3
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(100, n.occasions)
marked[,2] <- rep(60, n.occasions)
marked[,3] <- rep(0, n.occasions)

#Define matrices with survival, transition and recapture probabilities. These are 4-dimensional matrices, with
# Dimension 1: state of departure
# Dimension 2: state of arrival
# Dimension 3: individual
# Dimension 4: time

#1. State process matrix
totrel <- sum(marked)*(n.occasions-1)
PSI.STATE <- array(NA, dim=c(n.states, n.states, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.STATE[,,i,t] <- matrix(c(
      phiA*(1-psiAB), phiA*psiAB, 1-phiA,
      phiB*psiBA, phiB*(1-psiBA), 1-phiB,
      0, 0, 1 ), nrow = n.states,
      byrow = TRUE)
  } 
} 

#2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.OBS[,,i,t] <- matrix(c(
      pA, 0, 1-pA,
      0, pB, 1-pB,
      0, 0, 1 ), nrow = n.states, byrow = TRUE)
  } 
} 

#Define function to simulate multistate capture-recapture data
simul.ms <- function(PSI.STATE, PSI.OBS, marked, unobservable = NA){
  #Unobservable: number of state that is unobservable
  n.occasions <- dim(PSI.STATE)[4] + 1
  CH <- CH.TRUE <- matrix(NA, ncol = n.occasions, nrow = sum(marked))
  #Define a vector with the occasion of marking
  mark.occ <- matrix(0, ncol = dim(PSI.STATE)[1], nrow = sum(marked))
  g <- colSums(marked)
  for (s in 1:dim(PSI.STATE)[1]){
    if (g[s]==0) next # To avoid error message if nothing to replace
    mark.occ[(cumsum(g[1:s])-g[s]+1)[s]:cumsum(g[1:s])[s],s] <-
      rep(1:n.occasions, marked[1:n.occasions,s])
  }
  for (i in 1:sum(marked)){
    for (s in 1:dim(PSI.STATE)[1]){
      if (mark.occ[i,s]==0) next
      first <- mark.occ[i,s]
      CH[i,first] <- s
      CH.TRUE[i,first] <- s
    } #s
    for (t in (first+1):n.occasions){
      #Multinomial trials for state transitions
      if (first==n.occasions) next
      state <- which(rmultinom(1, 1, PSI.STATE[CH.TRUE[i,t-1],,i,t-1])==1)
      CH.TRUE[i,t] <- state
      #Multinomial trials for observation process
      event <- which(rmultinom(1, 1, PSI.OBS[CH.TRUE[i,t],,i,t-1])==1)
      CH[i,t] <- event
    } 
  } 
  #Replace the NA and the highest state number (dead) in the file by 0
  CH[is.na(CH)] <- 0
  CH[CH==dim(PSI.STATE)[1]] <- 0
  CH[CH==unobservable] <- 0
  id <- numeric(0)
  for (i in 1:dim(CH)[1]){
    z <- min(which(CH[i,]!=0))
    ifelse(z==dim(CH)[2], id <- c(id,i), id <- c(id))
  }
  return(list(CH=CH[-id,], CH.TRUE=CH.TRUE[-id,]))
  #CH: capture-histories to be used
  #CH.TRUE: capture-histories with perfect observation
}
#Execute function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked)
CH <- sim$CH

#Compute vector with occasion of first capture
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

#Recode CH matrix: note, a 0 is not allowed in WinBUGS!
#1 = seen alive in A, 2 = seen alive in B, 3 = not seen
rCH <- CH # Recoded CH
rCH[rCH==0] <- 3
#-------------------------------------------------------------------------------

# 9.2.3 Analysis of the Model
#-------------------------------------------------------------------------------
#Specify model in JAGS
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #------------------------------------------------
  #Parameters:
  # phiA: survival probability at site A
  # phiB: survival probability at site B
  # psiAB: movement probability from site A to site B
  # psiBA: movement probability from site B to site A
  # pA: recapture probability at site A
  # pB: recapture probability at site B
  #------------------------------------------------
  #States (S):
  # 1 alive at A
  # 2 alive at B
  # 3 dead
  #Observations (O):
  # 1 seen at A
  # 2 seen at B
  # 3 not seen
  #------------------------------------------------
  
  #Priors and constraints
  for (t in 1:(n.occasions-1)){
    phiA[t] <- mean.phi[1]
    phiB[t] <- mean.phi[2]
    psiAB[t] <- mean.psi[1]
    psiBA[t] <- mean.psi[2]
    pA[t] <- mean.p[1]
    pB[t] <- mean.p[2]
  }
  for (u in 1:2){
    mean.phi[u] ~ dunif(0, 1) #Priors for mean state-spec. survival
    mean.psi[u] ~ dunif(0, 1) #Priors for mean transitions
    mean.p[u] ~ dunif(0, 1) #Priors for mean state-spec. recapture
  }
  
  #Define state-transition and observation matrices
  for (i in 1:nind){
    #Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
      ps[1,i,t,1] <- phiA[t] * (1-psiAB[t])
      ps[1,i,t,2] <- phiA[t] * psiAB[t]
      ps[1,i,t,3] <- 1-phiA[t]
      ps[2,i,t,1] <- phiB[t] * psiBA[t]
      ps[2,i,t,2] <- phiB[t] * (1-psiBA[t])
      ps[2,i,t,3] <- 1-phiB[t]
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- 1
      #Define probabilities of O(t) given S(t)
      po[1,i,t,1] <- pA[t]
      po[1,i,t,2] <- 0
      po[1,i,t,3] <- 1-pA[t]
      po[2,i,t,1] <- 0
      po[2,i,t,2] <- pB[t]
      po[2,i,t,3] <- 1-pB[t]
      po[3,i,t,1] <- 0
      po[3,i,t,2] <- 0
      po[3,i,t,3] <- 1
    } 
  } 
  
  #Likelihood
  for (i in 1:nind){
    #Define latent state at first capture
    z[i,f[i]] <- Y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
      #State process: draw S(t) given S(t−1)
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
      #Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } 
  } 
}

#Function to create known latent states z
known.state.ms <- function(ms, notseen){
  #notseen: label for ‘not seen’
  state <- ms
  state[state==notseen] <- NA
  for (i in 1:dim(ms)[1]){
    m <- min(which(!is.na(state[i,])))
    state[i,m] <- NA
  }
  return(state)
}

#Function to create initial values for unknown z
ms.init.z <- function(ch, f){
  for (i in 1:dim(ch)[1]){ch[i,1:f[i]] <- NA}
  states <- max(ch, na.rm = TRUE)
  known.states <- 1:(states-1)
  v <- which(ch==states)
  ch[-v] <- NA
  ch[v] <- sample(known.states, length(v), replace = TRUE)
  return(ch)
}

#Bundle data
jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1], z = known.state.ms(rCH, 3))

#Initial values
inits <- function(){list(mean.phi = runif(2, 0, 1), mean.psi = runif(2,0, 1), mean.p = runif(2, 0, 1), z = ms.init.z(rCH, f))}

#Parameters monitored
params <- c("mean.phi", "mean.psi", "mean.p")

#MCMC settings
ni <- 10000
nt <- 6
nb <- 2000
nc <- 3

#Call JAGS from R
ms <- jags(data  = jags.data,
              inits = inits,
              parameters.to.save = params,
              model.file = jags.model.txt,
              n.chains = nc,
              n.thin= nt,
              n.iter = ni,
              n.burnin = nb)

print(ms, digits = 3)

k<-mcmcplots::as.mcmc.rjags(ms)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#Create Plot
par(mfrow = c(3, 2), las = 1)
hist(ms$BUGSoutput$sims.list$mean.phi[,1], col = "gray", main = "", xlab = expression(phi[A]), ylim=c(0,1300))
abline(v = phiA, col = "red")
hist(ms$BUGSoutput$sims.list$mean.phi[,2], col = "gray", main = "", xlab = expression(phi[B]), ylim=c(0,1300), ylab="")
abline(v = phiB, col="red")
hist(ms$BUGSoutput$sims.list$mean.psi[,1], col = "gray", main = "", xlab = expression(psi[AB]), ylim=c(0,1300))
abline(v = psiAB, col="red")
hist(ms$BUGSoutput$sims.list$mean.psi[,2], col = "gray", main = "", xlab = expression(psi[BA]), ylab="", ylim=c(0,1300))
abline(v = psiBA, col="red")
hist(ms$BUGSoutput$sims.list$mean.p[,1], col = "gray", main = "", xlab = expression(p[A]), ylim=c(0,1300))
abline(v = pA, col = "red")
hist(ms$BUGSoutput$sims.list$mean.p[,2], col = "gray", main = "", xlab = expression(p[B]), ylab="", ylim=c(0,1300))
abline(v = pB, col = "red")

#Specify alternative model 1 in JAGS
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  #Priors and constraints
  for (t in 1:(n.occasions-1)){
    phiA[t] <- mean.phi[1]
    phiB[t] <- mean.phi[2]
    psiAB[t] <- mean.psi[1]
    psiBA[t] <- mean.psi[2]
    pA[t] <- mean.p[1]
    pB[t] <- mean.p[2]
  }
  for (u in 1:2){
    mean.phi[u] ~ dunif(0, 1) #Priors for mean state-spec. survival
    mean.psi[u] ~ dunif(0, 1) #Priors for mean transitions
    mean.p[u] ~ dunif(0, 1) #Priors for mean state-spec. recapture
  }
  #Define state-transition and observation matrices
  #Define probabilities of state S(t+1) given S(t)
  for (t in 1:(n.occasions-1)){
    ps[1,t,1] <- phiA[t] * (1-psiAB[t])
    ps[1,t,2] <- phiA[t] * psiAB[t]
    ps[1,t,3] <- 1-phiA[t]
    ps[2,t,1] <- phiB[t] * psiBA[t]
    ps[2,t,2] <- phiB[t] * (1-psiBA[t])
    ps[2,t,3] <- 1-phiB[t]
    ps[3,t,1] <- 0
    ps[3,t,2] <- 0
    ps[3,t,3] <- 1
    
    #Define probabilities of O(t) given S(t)
    po[1,t,1] <- pA[t]
    po[1,t,2] <- 0
    po[1,t,3] <- 1-pA[t]
    po[2,t,1] <- 0
    po[2,t,2] <- pB[t]
    po[2,t,3] <- 1-pB[t]
    po[3,t,1] <- 0
    po[3,t,2] <- 0
    po[3,t,3] <- 1
  } 
  
  #Likelihood
  for (i in 1:nind){
    #Define latent state at first capture
    z[i,f[i]] <- Y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
      #State process: draw S(t) given S(t−1)
      z[i,t] ~ dcat(ps[z[i,t-1], t-1,])
      #Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t], t-1,])
    } 
  } 
}

#Specify alternative model 2 in JAGS
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  #Priors and constraints
  phiA ~ dunif(0, 1) #Prior for mean survival in A
  phiB ~ dunif(0, 1) #Prior for mean survival in B
  psiAB ~ dunif(0, 1) #Prior for mean movement from A to B
  psiBA ~ dunif(0, 1) #Prior for mean movement from B to A
  pA ~ dunif(0, 1) #Prior for mean recapture in A
  pB ~ dunif(0, 1) #Prior for mean recapture in B
  
  #Define state-transition and observation matrices
  #Define probabilities of state S(t+1) given S(t)
  ps[1,1] <- phiA * (1-psiAB)
  ps[1,2] <- phiA * psiAB
  ps[1,3] <- 1-phiA
  ps[2,1] <- phiB * psiBA
  ps[2,2] <- phiB * (1-psiBA)
  ps[2,3] <- 1-phiB
  ps[3,1] <- 0
  ps[3,2] <- 0
  ps[3,3] <- 1
  
  #Define probabilities of O(t) given S(t)
  po[1,1] <- pA
  po[1,2] <- 0
  po[1,3] <- 1-pA
  po[2,1] <- 0
  po[2,2] <- pB
  po[2,3] <- 1-pB
  po[3,1] <- 0
  po[3,2] <- 0
  po[3,3] <- 1
  
  #Likelihood
  for (i in 1:nind){
    #Define latent state at first capture
    z[i,f[i]] <- Y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
      #State process: draw S(t) given S(t−1)
      z[i,t] ~ dcat(ps[z[i,t−1],])
      #Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t],])
    } 
  } 
}
#-------------------------------------------------------------------------------

# 9.3 Accounting for Temporary Emigration
# 9.3.2 Generation of Simulated Data
#-------------------------------------------------------------------------------
#Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
phi <- 0.85
psiIO <- 0.2
psiOI <- 0.3
p <- 0.7
n.occasions <- 8
n.states <- 3
n.obs <- 2
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(70, n.occasions) #Present
marked[,2] <- rep(0, n.occasions) #Absent
marked[,3] <- rep(0, n.occasions) #Dead

#Define matrices with survival, transition and recapture probabilities. These are 4-dimensional matrices, with
#Dimension 1: state of departure
# Dimension 2: state of arrival
# Dimension 3: individual
# Dimension 4: time

#1. State process matrix
totrel <- sum(marked)*(n.occasions-1)
PSI.STATE <- array(NA, dim=c(n.states, n.states, totrel,
                             n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.STATE[,,i,t] <- matrix(c(
      phi*(1-psiIO), phi*psiIO, 1-phi,
      phi*psiOI, phi*(1-psiOI), 1-phi,
      0, 0, 1 ), nrow = n.states,
      byrow = TRUE)
  } 
} 

#2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.OBS[,,i,t] <- matrix(c(
      p, 1-p,
      0, 1,
      0, 1 ), nrow = n.states, byrow = TRUE)
  } 
} 

#Execute simulation function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked)
CH <- sim$CH

#Compute vector with occasion of first capture
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)
# Recode CH matrix: note, a 0 is not allowed!
# 1 = seen alive, 2 = not seen
rCH <- CH # Recoded CH
rCH[rCH==0] <- 2
#-------------------------------------------------------------------------------

# 9.3.3 Analysis of the Model
#-------------------------------------------------------------------------------
#Specify model in JAGS
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #--------------------------------
  #Parameters:
  # phi: survival probability
  # psiIO: probability to emigrate
  # psiOI: probability to immigrate
  # p: recapture probability
  #--------------------------------
  #States (S):
  # 1 alive and present
  # 2 alive and absent
  # 3 dead
  #Observations (O):
  # 1 seen
  # 2 not seen
  #--------------------------------
  
  #Priors and constraints
  for (t in 1:(n.occasions-1)){
    phi[t] <- mean.phi
    psiIO[t] <- mean.psiIO
    psiOI[t] <- mean.psiOI
    p[t] <- mean.p
  }
  mean.phi ~ dunif(0, 1) #Prior for mean survival
  mean.psiIO ~ dunif(0, 1) #Prior for mean temp. emigration
  mean.psiOI ~ dunif(0, 1) #Prior for mean temp. immigration
  mean.p ~ dunif(0, 1) #Prior for mean recapture
  
  #Define state-transition and observation matrices
  for (i in 1:nind){
    #Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
      ps[1,i,t,1] <- phi[t] * (1-psiIO[t])
      ps[1,i,t,2] <- phi[t] * psiIO[t]
      ps[1,i,t,3] <- 1-phi[t]
      ps[2,i,t,1] <- phi[t] * psiOI[t]
      ps[2,i,t,2] <- phi[t] * (1-psiOI[t])
      ps[2,i,t,3] <- 1-phi[t]
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- 1
      
      #Define probabilities of O(t) given S(t)
      po[1,i,t,1] <- p[t]
      po[1,i,t,2] <- 1-p[t]
      po[2,i,t,1] <- 0
      po[2,i,t,2] <- 1
      po[3,i,t,1] <- 0
      po[3,i,t,2] <- 1
    } 
  } 
  
  #Likelihood
  for (i in 1:nind){
    #Define latent state at first capture
    z[i,f[i]] <- Y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
      #State process: draw S(t) given S(t−1)
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
      #Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } 
  } 
}

#Bundle data
jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1], z = known.state.ms(rCH, 2))

#Initial values
inits <- function(){list(mean.phi = runif(1, 0, 1), mean.psiIO = runif(1, 0, 1), mean.psiOI = runif(1, 0, 1), mean.p = runif(1, 0, 1), z = ms.init.z(rCH, f))}

# Parameters monitored
params <- c("mean.phi", "mean.psiIO", "mean.psiOI", "mean.p")

# MCMC settings
ni <- 50000
nt <- 10
nb <- 10000
nc <- 3

#Call JAGS from R
tempemi <- jags(data  = jags.data,
           inits = inits,
           parameters.to.save = params,
           model.file = jags.model.txt,
           n.chains = nc,
           n.thin= nt,
           n.iter = ni,
           n.burnin = nb)

print(tempemi, digits = 3)

k<-mcmcplots::as.mcmc.rjags(tempemi)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#Create plot
par(mfrow = c(2, 2), las = 1)
hist(tempemi$BUGSoutput$sims.list$mean.phi, col = "gray", main = "", xlab = expression(phi))
abline(v = phi, col = "red", lwd = 2)
hist(tempemi$BUGSoutput$sims.list$mean.psiIO, col = "gray", main = "", xlab = expression(psi[IO]), ylab = "")
abline(v = psiIO, col = "red", lwd = 2)
hist(tempemi$BUGSoutput$sims.list$mean.psiOI, col = "gray", main = "", xlab = expression(psi[OI]))
abline(v = psiOI, col = "red", lwd = 2)
hist(tempemi$BUGSoutput$sims.list$mean.p, col = "gray", main = "", xlab = expression(p), ylab = "")
abline(v = p, col = "red", lwd = 2)
#-------------------------------------------------------------------------------

# 9.4 Estimation of Age-Specific Probability of First Breeding
# 9.4.2 Generation of Simulated Data
#-------------------------------------------------------------------------------
#Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
phi.1 <- 0.4
phi.2 <- 0.7
phi.ad <- 0.8
alpha.1 <- 0.2
alpha.2 <- 0.6
p.NB <- 0.5
p.B <- 0.7
n.occasions <- 7
n.states <- 5
n.obs <- 4
marked <- matrix(0, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(100, n.occasions) #Releases only as juveniles

#Define matrices with survival, transition and recapture probabilities
#These are 4-dimensional matrices, with
#Dimension 1: state of departure
#Dimension 2: state of arrival
#Dimension 3: individual
#Dimension 4: time

#1. State process matrix
totrel <- sum(marked)*(n.occasions-1)
PSI.STATE <- array(NA, dim=c(n.states, n.states, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.STATE[,,i,t] <- matrix(c(
      0, phi.1*(1-alpha.1), 0, phi.1*alpha.1, 1-phi.1,
      0, 0, phi.2*(1-alpha.2), phi.2*alpha.2, 1-phi.2,
      0, 0, 0, phi.ad, 1-phi.ad,
      0, 0, 0, phi.ad, 1-phi.ad,
      0, 0, 0, 0, 1), nrow = n.states,byrow = TRUE)
  } 
} 

#2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.OBS[,,i,t] <- matrix(c(
      0, 0, 0, 1,
      0, p.NB, 0, 1-p.NB,
      0, p.NB, 0, 1-p.NB,
      0, 0, p.B, 1-p.B,
      0, 0, 0, 1), nrow = n.states, byrow = TRUE)
  } 
} 

#Execute simulation function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked)
CH <- sim$CH

#Compute vector with occasion of first capture
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

#Recode CH matrix: note, a 0 is not allowed!
#1 = seen as juv, 2 = seen no rep, 3 = seen rep, 4 = not seen
rCH <- CH # Recoded CH
rCH[rCH==0] <- 4
#-------------------------------------------------------------------------------

# 9.4.3 Analysis of the Model
#-------------------------------------------------------------------------------
#Specify model in JAGS
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #------------------------------------------------
  #Parameters:
  # phi.1: first year survival probability
  # phi.2: second year survival probability
  # phi.ad: adult survival probability
  # alpha.1: probability to start breeding when 1 year old
  # alpha.2: probability to start breeding when 2 years old
  # p.NB: recapture probability of non-breeders
  # p.B: recapture probability of breeders
  #------------------------------------------------
  #States (S):
  # 1 juvenile
  # 2 not yet breeding at age 1 year
  # 3 not yet breeding at age 2 years
  # 4 breeder
  # 5 dead
  #Observations (O):
  # 1 seen as juvenile
  # 2 seen as not yet breeding
  # 3 seen breeding
  # 4 not seen
  #------------------------------------------------
  
  #Priors and constraints
  for (t in 1:(n.occasions-1)){
    phi.1[t] <- mean.phi1
    phi.2[t] <- mean.phi2
    phi.ad[t] <- mean.phiad
    alpha.1[t] <- mean.alpha1
    alpha.2[t] <- mean.alpha2
    p.NB[t] <- mean.pNB
    p.B[t] <- mean.pB
  }
  mean.phi1 ~ dunif(0, 1) #Prior for mean 1y survival
  mean.phi2 ~ dunif(0, 1) #Prior for mean 2y survival
  mean.phiad ~ dunif(0, 1) #Prior for mean ad survival
  mean.alpha1 ~ dunif(0, 1) #Prior for mean 1y breeding prob.
  mean.alpha2 ~ dunif(0, 1) #Prior for mean 2y breeding prob.
  mean.pNB ~ dunif(0, 1) #Prior for mean recapture non-breeders
  mean.pB ~ dunif(0, 1) #Prior for mean recapture breeders
  
  #Define state-transition and observation matrices
  for (i in 1:nind){
    #Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
      ps[1,i,t,1] <- 0
      ps[1,i,t,2] <- phi.1[t] * (1-alpha.1[t])
      ps[1,i,t,3] <- 0
      ps[1,i,t,4] <- phi.1[t] * alpha.1[t]
      ps[1,i,t,5] <- 1-phi.1[t]
      ps[2,i,t,1] <- 0
      ps[2,i,t,2] <- 0
      ps[2,i,t,3] <- phi.2[t] * (1-alpha.2[t])
      ps[2,i,t,4] <- phi.2[t] * alpha.2[t]
      ps[2,i,t,5] <- 1-phi.2[t]
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- 0
      ps[3,i,t,4] <- phi.ad[t]
      ps[3,i,t,5] <- 1-phi.ad[t]
      ps[4,i,t,1] <- 0
      ps[4,i,t,2] <- 0
      ps[4,i,t,3] <- 0
      ps[4,i,t,4] <- phi.ad[t]
      ps[4,i,t,5] <- 1-phi.ad[t]
      ps[5,i,t,1] <- 0
      ps[5,i,t,2] <- 0
      ps[5,i,t,3] <- 0
      ps[5,i,t,4] <- 0
      ps[5,i,t,5] <- 1
      
      #Define probabilities of O(t) given S(t)
      po[1,i,t,1] <- 0
      po[1,i,t,2] <- 0
      po[1,i,t,3] <- 0
      po[1,i,t,4] <- 1
      po[2,i,t,1] <- 0
      po[2,i,t,2] <- p.NB[t]
      po[2,i,t,3] <- 0
      po[2,i,t,4] <- 1-p.NB[t]
      po[3,i,t,1] <- 0
      po[3,i,t,2] <- p.NB[t]
      po[3,i,t,3] <- 0
      po[3,i,t,4] <- 1-p.NB[t]
      po[4,i,t,1] <- 0
      po[4,i,t,2] <- 0
      po[4,i,t,3] <- p.B[t]
      po[4,i,t,4] <- 1-p.B[t]
      po[5,i,t,1] <- 0
      po[5,i,t,2] <- 0
      po[5,i,t,3] <- 0
      po[5,i,t,4] <- 1
    } 
  } 
  
  #Likelihood
  for (i in 1:nind){
    #Define latent state at first capture
    z[i,f[i]] <- Y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
      #State process: draw S(t) given S(t−1)
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
      #Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } 
  } 
}

#Bundle data
jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1])

#Initial values (note: function ch.init is defined in section 7.3)
inits <- function(){list(mean.phi1 = runif(1, 0, 1), mean.phi2 = runif(1, 0, 1), mean.phiad = runif(1, 0, 1), mean.alpha1 = runif(1, 0, 1), mean.alpha2 = runif(1, 0, 1), mean.pNB = runif(1, 0, 1), mean.pB = runif(1, 0, 1), z = ch.init(rCH, f))}

#Parameters monitored
params <- c("mean.phi1", "mean.phi2", "mean.phiad", "mean.alpha1", "mean.alpha2", "mean.pNB", "mean.pB")

#MCMC settings
ni <- 2000
nt <- 3
nb <- 1000
nc <- 3

#Call JAGS from R
agefirst <- jags(data  = jags.data,
                inits = inits,
                parameters.to.save = params,
                model.file = jags.model.txt,
                n.chains = nc,
                n.thin= nt,
                n.iter = ni,
                n.burnin = nb)

print(agefirst, digits = 3)

k<-mcmcplots::as.mcmc.rjags(agefirst)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#Create plot
par(mfrow = c(3, 3), las = 1)
hist(agefirst$BUGSoutput$sims.list$mean.phi1, col = "gray", main = "", xlab = expression(phi[1]))
abline(v = phi.1, col = "red", lwd = 2)
hist(agefirst$BUGSoutput$sims.list$mean.phi2, col = "gray", main = "", xlab = expression(phi[2]), ylab = "")
abline(v = phi.2, col = "red", lwd = 2)
hist(agefirst$BUGSoutput$sims.list$mean.phiad, col = "gray", main = "", xlab = expression(phi[ad]) , ylab = "")
abline(v = phi.ad, col = "red", lwd = 2)
hist(agefirst$BUGSoutput$sims.list$mean.alpha1, col = "gray", main = "", xlab = expression(alpha[1]))
abline(v = alpha.1, col = "red", lwd = 2)
hist(agefirst$BUGSoutput$sims.list$mean.alpha2, col = "gray", main = "", xlab = expression(alpha[2]) , ylab = "")
abline(v = alpha.2, col = "red", lwd = 2)
plot(0, type = "n", axes = F, ylab = "", xlab = "")
hist(agefirst$BUGSoutput$sims.list$mean.pNB, col = "gray", main = "", xlab = expression(p[NB]))
abline(v = p.NB, col = "red", lwd = 2)
hist(agefirst$BUGSoutput$sims.list$mean.pB, col = "gray", main = "", xlab = expression(p[B]) , ylab = "")
abline(v = p.B, col = "red", lwd = 2)
#-------------------------------------------------------------------------------

# 9.5 Joint Analysis of Capture-Recapture and Mark-Recovery Data
# 9.5.2 Generation of Simulated Data
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

# 9.5.3 Analysis of the Model
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

# 9.6 Estimation of Movement Among Three Sites
# 9.6.2 Generation of Simulated Data
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

# 9.6.3 Analysis of the Model
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

# 9.7 Real-Data Example: the Showy Lady's Slipper
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------