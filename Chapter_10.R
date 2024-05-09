# Bayesian Population Analysis using WinBUGS
# Chapter 10: Estimation of Survival, Recruitment, and Population Size from Capture-Recapture Data Using the Jolly-Seber Model 

library(R2jags) #to run JAGS
library(shinystan) #to run shiny stan
library(tidyverse) #to utilize pipe operators

# 10.3 Fitting the JS Model with Data Augmentation
# 10.3.1 The JS Model as a Restricted Dynamic Occupancy Model
#-------------------------------------------------------------------------------
#Specify model in JAGS
jags.js.rest.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  #Priors and constraints
  for (i in 1:M){
    for (t in 1:(n.occasions-1)){
      phi[i,t] <- mean.phi
    } 
    for (t in 1:n.occasions){
      p[i,t] <- mean.p
    }
  } 
  mean.phi ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  for (t in 1:n.occasions){
    gamma[t] ~ dunif(0, 1)
  } 
  
  #Likelihood
  for (i in 1:M){
    #First occasion
    #State process
    z[i,1] ~ dbern(gamma[1])
    mu1[i] <- z[i,1] * p[i,1]
    #Observation process
    y[i,1] ~ dbern(mu1[i])
    #Subsequent occasions
    for (t in 2:n.occasions){
      #State process
      q[i,t-1] <- 1-z[i,t-1] # Availability for recruitment
      mu2[i,t] <- phi[i,t-1] * z[i,t-1] + gamma[t] * prod(q[i,1:(t-1)])
      z[i,t] ~ dbern(mu2[i,t])
      #Observation process
      mu3[i,t] <- z[i,t] * p[i,t]
      y[i,t] ~ dbern(mu3[i,t])
    } 
  } 
  
  #Calculate derived population parameters
  for (t in 1:n.occasions){
    qgamma[t] <- 1-gamma[t]
  }
  cprob[1] <- gamma[1]
  for (t in 2:n.occasions){
    cprob[t] <- gamma[t] * prod(qgamma[1:(t-1)])
  } 
  psi <- sum(cprob[]) #Inclusion probability
  for (t in 1:n.occasions){
    b[t] <- cprob[t] / psi #Entry probability
  } 
  for (i in 1:M){
    recruit[i,1] <- z[i,1]
    for (t in 2:n.occasions){
      recruit[i,t] <- (1-z[i,t-1]) * z[i,t]
    } 
  } 
  for (t in 1:n.occasions){
    N[t] <- sum(z[1:M,t]) #Actual population size
    B[t] <- sum(recruit[1:M,t]) #Number of entries
  } 
  for (i in 1:M){
    Nind[i] <- sum(z[i,1:n.occasions])
    Nalive[i] <- 1-equals(Nind[i], 0)
  } 
  Nsuper <- sum(Nalive[]) #Superpopulation size
}
#-------------------------------------------------------------------------------

# 10.3.2 The JS Model as a Multistate Model
#-------------------------------------------------------------------------------
#Specify model in JAGS
jags.js.ms.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #-----------------------------------
  #Parameters:
  # phi: survival probability
  # gamma: removal entry probability
  # p: capture probability
  #-----------------------------------
  #States (S):
  # 1 not yet entered
  # 2 alive
  # 3 dead
  #Observations (O):
  # 1 seen
  # 2 not seen
  #-----------------------------------
  
  # Priors and constraints
  for (t in 1:(n.occasions-1)){
    phi[t] <- mean.phi
    gamma[t] ~ dunif(0, 1) #Prior for entry probabilities
    p[t] <- mean.p
  }
  mean.phi ~ dunif(0, 1) #Prior for mean survival
  mean.p ~ dunif(0, 1) #Prior for mean capture
  
  #Define state-transition and observation matrices
  for (i in 1:M){
    #Define probabilities of state S(t+1) given S(t)
    for (t in 1:(n.occasions-1)){
      ps[1,i,t,1] <- 1-gamma[t]
      ps[1,i,t,2] <- gamma[t]
      ps[1,i,t,3] <- 0
      ps[2,i,t,1] <- 0
      ps[2,i,t,2] <- phi[t]
      ps[2,i,t,3] <- 1-phi[t]
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- 1
      #Define probabilities of O(t) given S(t)
      po[1,i,t,1] <- 0
      po[1,i,t,2] <- 1
      po[2,i,t,1] <- p[t]
      po[2,i,t,2] <- 1-p[t]
      po[3,i,t,1] <- 0
      po[3,i,t,2] <- 1
    } 
  } 
  
  #Likelihood
  for (i in 1:M){
    #Define latent state at first occasion
    z[i,1] <- 1 #Make sure that all M individuals are in state 1 at t=1
    for (t in 2:n.occasions){
      #State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
      #Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } 
  } 
  
  #Calculate derived population parameters
  for (t in 1:(n.occasions-1)){
    qgamma[t] <- 1-gamma[t]
  }
  cprob[1] <- gamma[1]
  for (t in 2:(n.occasions-1)){
    cprob[t] <- gamma[t] * prod(qgamma[1:(t-1)])
  } 
  psi <- sum(cprob[]) #Inclusion probability
  for (t in 1:(n.occasions-1)){
    b[t] <- cprob[t] / psi #Entry probability
  } 
  for (i in 1:M){
    for (t in 2:n.occasions){
      al[i,t-1] <- equals(z[i,t], 2)
    } 
    for (t in 1:(n.occasions-1)){
      d[i,t] <- equals(z[i,t]-al[i,t],0)
    } 
    alive[i] <- sum(al[i,])
  } 
  for (t in 1:(n.occasions-1)){
    N[t] <- sum(al[,t]) #Actual population size
    B[t] <- sum(d[,t]) #Number of entries
  } 
  for (i in 1:M){
    w[i] <- 1-equals(alive[i],0)
  } 
  Nsuper <- sum(w[]) #Superpopulation size
}
#-------------------------------------------------------------------------------

# 10.3.3 The Superpopulation Parameterization
#-------------------------------------------------------------------------------
#Specify model in JAGS
jags.js.super.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors and constraints
  for (i in 1:M){
    for (t in 1:(n.occasions-1)){
      phi[i,t] <- mean.phi
    } 
    for (t in 1:n.occasions){
      p[i,t] <- mean.p
    } 
  } 
  mean.phi ~ dunif(0, 1) #Prior for mean survival
  mean.p ~ dunif(0, 1) #Prior for mean capture
  psi ~ dunif(0, 1) #Prior for inclusion probability
  
  #Dirichlet prior for entry probabilities
  for (t in 1:n.occasions){
    beta[t] ~ dgamma(1, 1)
    b[t] <- beta[t] / sum(beta[1:n.occasions])
  }
  
  #Convert entry probs to conditional entry probs
  nu[1] <- b[1]
  for (t in 2:n.occasions){
    nu[t] <- b[t] / (1-sum(b[1:(t-1)]))
  } 
  
  #Likelihood
  for (i in 1:M){
    #First occasion
    #State process
    w[i] ~ dbern(psi) #Draw latent inclusion
    z[i,1] ~ dbern(nu[1])
    #Observation process
    mu1[i] <- z[i,1] * p[i,1] * w[i]
    y[i,1] ~ dbern(mu1[i])
    #Subsequent occasions
    for (t in 2:n.occasions){
      #State process
      q[i,t-1] <- 1-z[i,t-1]
      mu2[i,t] <- phi[i,t-1] * z[i,t-1] + nu[t] * prod(q[i,1:(t-1)])
      z[i,t] ~ dbern(mu2[i,t])
      #Observation process
      mu3[i,t] <- z[i,t] * p[i,t] * w[i]
      y[i,t] ~ dbern(mu3[i,t])
    } 
  } 
  
  #Calculate derived population parameters
  for (i in 1:M){
    for (t in 1:n.occasions){
      u[i,t] <- z[i,t]*w[i] #Deflated latent state (u)
    }
  }
  for (i in 1:M){
    recruit[i,1] <- u[i,1]
    for (t in 2:n.occasions){
      recruit[i,t] <- (1-u[i,t-1]) * u[i,t]
    } 
  } 
  for (t in 1:n.occasions){
    N[t] <- sum(u[1:M,t]) #Actual population size
    B[t] <- sum(recruit[1:M,t]) #Number of entries
  }
  for (i in 1:M){
    Nind[i] <- sum(u[i,1:n.occasions])
    Nalive[i] <- 1-equals(Nind[i], 0)
  } #i
  Nsuper <- sum(Nalive[]) #Superpopulation size
}
#-------------------------------------------------------------------------------

# 10.4 Models with Constant Survival and Time-Dependent Entry
#-------------------------------------------------------------------------------
# Define parameter values
n.occasions <- 7 #Number of capture occasions
N <- 400 #Superpopulation size
phi <- rep(0.7, n.occasions-1) #Survival probabilities
b <- c(0.34, rep(0.11, n.occasions-1)) #Entry probabilities
p <- rep(0.5, n.occasions) #Capture probabilities
PHI <- matrix(rep(phi, (n.occasions-1)*N), ncol = n.occasions-1, nrow = N, byrow = T)
P <- matrix(rep(p, n.occasions*N), ncol = n.occasions, nrow = N, byrow = T)

#Function to simulate capture-recapture data under the JS model
simul.js <- function(PHI, P, b, N){
  B <- rmultinom(1, N, b) #Generate no. of entering ind. per occasion
  n.occasions <- dim(PHI)[2] + 1
  CH.sur <- CH.p <- matrix(0, ncol = n.occasions, nrow = N)
  #Define a vector with the occasion of entering the population
  ent.occ <- numeric()
  for (t in 1:n.occasions){
    ent.occ <- c(ent.occ, rep(t, B[t]))
  }
  
  #Simulate survival
  for (i in 1:N){
    CH.sur[i, ent.occ[i]] <- 1 #Write 1 when ind. enters the pop.
    if (ent.occ[i] == n.occasions) next
    for (t in (ent.occ[i]+1):n.occasions){
      #Bernoulli trial: has individual survived occasion?
      sur <- rbinom(1, 1, PHI[i,t-1])
      ifelse (sur==1, CH.sur[i,t] <- 1, break)
    } 
  } 
  
  #Simulate capture
  for (i in 1:N){
    CH.p[i,] <- rbinom(n.occasions, 1, P[i,])
  } 
  
  #Full capture-recapture matrix
  CH <- CH.sur * CH.p
  
  #Remove individuals never captured
  cap.sum <- rowSums(CH)
  never <- which(cap.sum == 0)
  CH <- CH[-never,]
  Nt <- colSums(CH.sur) #Actual population size
  return(list(CH=CH, B=B, N=Nt))
}

#Execute simulation function
sim <- simul.js(PHI, P, b, N)
CH <- sim$CH
#-------------------------------------------------------------------------------

# 10.4.1 Analysis of the JS Model as a Restricted Occupancy Model
#-------------------------------------------------------------------------------
#Augment the capture-histories by nz pseudo-individuals
nz <- 500
CH.aug <- rbind(CH, matrix(0, ncol = dim(CH)[2], nrow = nz))

#Bundle data
jags.data <- list(y = CH.aug, n.occasions = dim(CH.aug)[2], M = dim(CH.aug)[1])

#Initial values from vogelwarte.ch
#Good initial values for the latent state z are needed to run the model in JAGS. The simplest option that works is to give just a matrix with a 1 at all places.
z.init <- CH.aug
z.init[z.init==0] <- 1
inits <- function(){list(mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1), z = z.init)}  

#Parameters monitored
params <- c("psi", "mean.p", "mean.phi", "b", "Nsuper", "N", "B", "gamma")

# MCMC settings
ni <- 5000
nt <- 3
nb <- 2000
nc <- 3

#Call JAGS from R
js.occ <- jags(data  = jags.data,
           inits = inits,
           parameters.to.save = params,
           model.file = jags.js.rest.txt,
           n.chains = nc,
           n.thin= nt,
           n.iter = ni,
           n.burnin = nb)

print(js.occ, digits = 3)

k<-mcmcplots::as.mcmc.rjags(js.occ)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan
#-------------------------------------------------------------------------------

# 10.4.2 Analysis of the JS Model as a Multistate Model
#-------------------------------------------------------------------------------
#Add dummy occasion
CH.du <- cbind(rep(0, dim(CH)[1]), CH)

# Code to fix error from: https://groups.google.com/g/hmecology/c/S4HO-tnzep8?pli=1 START
my.z.init <- CH.du

first.one <- apply(my.z.init[,1:ncol(CH.du)], 1, function(x) min(which(x == 1)))
last.one  <- apply(my.z.init[,1:ncol(CH.du)], 1, function(x) max(which(x == 1)))

for(i in 1:nrow(my.z.init)) {
  my.z.init[i,     first.one[i]  : last.one[i]        ] = 2
  if(first.one[i] > 1)               my.z.init[i,                1  : (first.one[i] - 1) ] = 1
  if(last.one[i]  < ncol(my.z.init)) my.z.init[i, (last.one[i] + 1) : ncol(my.z.init)    ] = 3
}
# Code to fix error from: https://groups.google.com/g/hmecology/c/S4HO-tnzep8?pli=1 END

#Augment data
nz <- 500
CH.ms <- rbind(CH.du, matrix(0, ncol = dim(CH.du)[2], nrow = nz))

#Recode CH matrix: a 0 is not allowed in WinBUGS!
CH.ms[CH.ms==0] <- 2 #Not seen = 2, seen = 1

#Code to fix error from: https://groups.google.com/g/hmecology/c/S4HO-tnzep8?pli=1
my.z.init.ms <- rbind(my.z.init, matrix(0, ncol = dim(my.z.init)[2], nrow = nz))
my.z.init.ms[my.z.init.ms==0] <- 1

#Bundle data
jags.data <- list(y = CH.ms, n.occasions = dim(CH.ms)[2], M = dim(CH.ms)[1])

inits <- function(){list(mean.phi = runif(1, 0, 1), #Code to fix error from: https://groups.google.com/g/hmecology/c/S4HO-tnzep8?pli=1
                         mean.p = runif(1, 0, 1),
                         z = cbind(rep(NA, dim(my.z.init.ms)[1]), my.z.init.ms[,-1]))}

#Parameters monitored
params <- c("mean.p", "mean.phi", "b", "Nsuper", "N", "B")

# MCMC settings
ni <- 20000
nt <- 3
nb <- 5000
nc <- 3

#Call JAGS from R
js.ms <- jags(data  = jags.data,
               inits = inits,
               parameters.to.save = params,
               model.file = jags.js.ms.txt,
               n.chains = nc,
               n.thin= nt,
               n.iter = ni,
               n.burnin = nb)

print(js.ms, digits = 3)

k<-mcmcplots::as.mcmc.rjags(js.ms)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan
#-------------------------------------------------------------------------------

# 10.4.3 Analysis of the JS Model Under the Superpopulation Parameterization
#-------------------------------------------------------------------------------
#Augment capture-histories by nz pseudo-individuals
nz <- 500
CH.aug <- rbind(CH, matrix(0, ncol = dim(CH)[2], nrow = nz))

#Bundle data
jags.data <- list(y = CH.aug, n.occasions = dim(CH.aug)[2], M = dim(CH.aug)[1])

#Initial values
#Good initial values for the latent state z are needed to run the model in JAGS. The simplest option that works is to give just a matrix with a 1 at all places. Moreover, initial values for the inclusion parameter w should also be given, the simplest option is to provide a vector with 1's.
z.init <- CH.aug
z.init[z.init==0] <- 1
w.init <- rep(1, nrow(CH.aug))

inits <- function(){list(mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1), psi = runif(1, 0, 1), w = w.init, z = z.init)}  
# Parameters monitored
params <- c("psi", "mean.p", "mean.phi", "b", "Nsuper", "N", "B", "nu")

#MCMC settings
ni <- 5000
nt <- 3
nb <- 2000
nc <- 3

#Call JAGS from R
js.super <- jags(data  = jags.data,
              inits = inits,
              parameters.to.save = params,
              model.file = jags.js.super.txt,
              n.chains = nc,
              n.thin= nt,
              n.iter = ni,
              n.burnin = nb)

print(js.super, digits = 3)

k<-mcmcplots::as.mcmc.rjags(js.super)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan
#-------------------------------------------------------------------------------

# 10.4.4 Comparison of Estimates
#-------------------------------------------------------------------------------
#Code to produce Fig. 10.4
par(mfrow = c(1,2), mar = c(5, 6, 2, 1), mgp=c(3.4, 1, 0), las = 1)
plot(density(js.occ$BUGSoutput$sims.list$Nsuper), main = "", xlab = "",
     ylab = "Density", frame = FALSE, lwd = 2, ylim=c(0, 0.023),
     col = "blue")
points(density(js.ms$BUGSoutput$sims.list$Nsuper), type = "l", lty = 2,
       col = "blue", lwd = 2)
points(density(js.super$BUGSoutput$sims.list$Nsuper), type = "l", lty = 3,
       col = "blue", lwd = 2)
abline(v = N, col = "red", lwd = 2)
mtext("Size of superpopulation", 1, line = 3)

b1.lower <- b2.lower <- b3.lower <- b1.upper <- b2.upper <- b3.upper <- numeric()

for (t in 1:n.occasions){
  b1.lower[t] <- quantile(js.occ$BUGSoutput$sims.list$b[,t], 0.025)
  b2.lower[t] <- quantile(js.ms$BUGSoutput$sims.list$b[,t], 0.025)
  b3.lower[t] <- quantile(js.super$BUGSoutput$sims.list$b[,t], 0.025)
  b1.upper[t] <- quantile(js.occ$BUGSoutput$sims.list$b[,t], 0.975)
  b2.upper[t] <- quantile(js.ms$BUGSoutput$sims.list$b[,t], 0.975)
  b3.upper[t] <- quantile(js.super$BUGSoutput$sims.list$b[,t], 0.975)
}

time <- 1:n.occasions
plot(x = time-0.25, y = js.occ$BUGSoutput$mean$b, xlab = "", ylab = "Entry probability", frame = FALSE, las = 1, xlim = c(0.5, 7.5), pch = 16, ylim = c(0, max(c(b1.upper, b2.upper))))
segments(time-0.25, b1.lower, time-0.25, b1.upper)
points(x = time, y = js.ms$BUGSoutput$mean$b[1:7], pch = 1)
segments(time, b2.lower, time, b2.upper)
points(x = time+0.25, y = js.super$BUGSoutput$mean$b, pch = 17)
segments(time+0.25, b3.lower, time+0.25, b3.upper)
points(x = time, y = b, pch = 18, col = "red")
mtext("Year", 1, line = 3)
mtext("Figure 10.4", side= 3, line= -1.5, outer= T) #adding main title to multiplot
#-------------------------------------------------------------------------------

# 10.5 Models with Individual Capture Heterogeneity 
#-------------------------------------------------------------------------------
#Define parameter values
n.occasions <- 8 #Number of capture occasions
N <- 300 #Size of the superpopulation
phi <- rep(0.75, n.occasions-1) #Survival probabilities
b <- c(0.37, rep(0.09, n.occasions-1)) #Entry probabilities
mean.p <- 0.6 #Mean capture probability
var.p <- 1 #Indv. Variance of capture prob.
p <- plogis(rnorm(N, qlogis(mean.p), var.p^0.5))
PHI <- matrix(rep(phi, (n.occasions-1)*N), ncol = n.occasions-1, nrow = N, byrow = T)
P <- matrix(rep(p, n.occasions), ncol = n.occasions, nrow = N, byrow = F)

#Execute simulation function
sim <- simul.js(PHI, P, b, N)
CH <- sim$CH

#Specify model in JAGS
jags.js.super.indran.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors and constraints
  for (i in 1:M){
    for (t in 1:(n.occasions-1)){
      phi[i,t] <- mean.phi
    } 
    for (t in 1:n.occasions){
      logit(p[i,t]) <- mean.lp + epsilon[i]
    } 
  } 
  mean.phi ~ dunif(0, 1) #Prior for mean survival
  mean.lp <- log(mean.p / (1-mean.p))
  mean.p ~ dunif(0, 1) #Prior for mean capture
  for (i in 1:M){
    epsilon[i] ~ dnorm(0, tau);T(-15,15)
  }
  tau <- 1/(sigma*sigma)
  sigma ~ dunif(0, 5) #Prior for sd of indv. variation of p
  sigma2 <- (sigma*sigma)
  psi ~ dunif(0, 1) #Prior for inclusion probability
  
  #Dirichlet prior for entry probabilities
  for (t in 1:n.occasions){
    beta[t] ~ dgamma(1, 1)
    b[t] <- beta[t] / sum(beta[1:n.occasions])
  }
  
  #Convert entry probs to conditional entry probs
  nu[1] <- b[1]
  for (t in 2:n.occasions){
    nu[t] <- b[t] / (1-sum(b[1:(t-1)]))
  }
  
  #Likelihood
  for (i in 1:M){
    #First occasion
    #State process
    w[i] ~ dbern(psi) #Draw latent inclusion
    z[i,1] ~ dbern(nu[1])
    #Observation process
    mu1[i] <- z[i,1] * p[i,1] * w[i]
    y[i,1] ~ dbern(mu1[i])
    #Subsequent occasions
    for (t in 2:n.occasions){
      #State process
      q[i,t-1] <- 1-z[i,t-1]
      mu2[i,t] <- phi[i,t-1] * z[i,t-1] + nu[t] * prod(q[i,1:(t-1)])
      z[i,t] ~ dbern(mu2[i,t])
      #Observation process
      mu3[i,t] <- z[i,t] * p[i,t] * w[i]
      y[i,t] ~ dbern(mu3[i,t])
    } 
  } 
  
  #Calculate derived population parameters
  for (i in 1:M){
    for (t in 1:n.occasions){
      u[i,t] <- z[i,t]*w[i] #Deflated latent state (u)
    }
  }
  for (i in 1:M){
    recruit[i,1] <- u[i,1]
    for (t in 2:n.occasions){
      recruit[i,t] <- (1-u[i,t-1]) * u[i,t]
    } 
  } 
  for (t in 1:n.occasions){
    N[t] <- sum(u[1:M,t]) #Actual population size
    B[t] <- sum(recruit[1:M,t]) #Number of entries
  } 
  for (i in 1:M){
    Nind[i] <- sum(u[i,1:n.occasions])
    Nalive[i] <- 1-equals(Nind[i], 0)
  } 
  Nsuper <- sum(Nalive[]) #Superpopulation size
}

#Augment the capture-histories by nz pseudo-individuals
nz <- 300
CH.aug <- rbind(CH, matrix(0, ncol = dim(CH)[2], nrow = nz))

#Bundle data
jags.data <- list(y = CH.aug, n.occasions = dim(CH.aug)[2], M = dim(CH.aug)[1])

#Initial values
#Good initial values for the latent state z are needed to run the model in JAGS. The simplest option that works is to give just a matrix with a 1 at all places. Moreover, initial values for the inclusion parameter w should also be given, the simplest option is to provide a vector with 1's.
z.init <- CH.aug
z.init[z.init==0] <- 1
w.init <- rep(1, nrow(CH.aug))
inits <- function(){list(mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1), sigma = runif(1, 0, 1), w = w.init, z = z.init)}  

#Parameters monitored
params <- c("sigma2","psi", "mean.p", "mean.phi", "N", "Nsuper", "b", "B")

# MCMC settings
ni <- 20000
nt <- 6
nb <- 5000
nc <- 3

#Call JAGS from R
js.ran <- jags(data  = jags.data,
              inits = inits,
              parameters.to.save = params,
              model.file = jags.js.super.indran.txt,
              n.chains = nc,
              n.thin= nt,
              n.iter = ni,
              n.burnin = nb)

print(js.ran, digits = 3)

k<-mcmcplots::as.mcmc.rjags(js.ran)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan
#-------------------------------------------------------------------------------

# 10.7 Analysis of a Real Data Set: Survival, Recruitment and Population Size of Leisler's Bats
#-------------------------------------------------------------------------------
#Specify model in JAGS
jags.js.tempran.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors and constraints
  for (i in 1:M){
    for (t in 1:(n.occasions-1)){
      logit(phi[i,t]) <- mean.lphi + epsilon[t]
    } 
    for (t in 1:n.occasions){
      p[i,t] <- mean.p
    } 
  } 
  mean.p ~ dunif(0, 1) #Prior for mean capture
  mean.phi ~ dunif(0, 1) #Prior for mean survival
  mean.lphi <- log(mean.phi / (1-mean.phi))
  for (t in 1:(n.occasions-1)){
    epsilon[t] ~ dnorm(0, tau)
  }
  tau <- 1/(sigma*sigma)
  sigma ~ dunif(0, 5) #Prior for sd of indv. variation of phi
  sigma2 <- (sigma*sigma)
  for (t in 1:n.occasions){
    gamma[t] ~ dunif(0, 1)
  } 
  
  #Likelihood
  for (i in 1:M){
    #First occasion
    #State process
    z[i,1] ~ dbern(gamma[1])
    mu1[i] <- z[i,1] * p[i,1]
    #Observation process
    y[i,1] ~ dbern(mu1[i])
    #Subsequent occasions
    for (t in 2:n.occasions){
      #State process
      q[i,t-1] <- 1-z[i,t-1]
      mu2[i,t] <- phi[i,t-1] * z[i,t-1] + gamma[t] * prod(q[i,1:(t-1)])
      z[i,t] ~ dbern(mu2[i,t])
      #Observation process
      mu3[i,t] <- z[i,t] * p[i,t]
      y[i,t] ~ dbern(mu3[i,t])
    } 
  } 
  
  #Calculate derived population parameters
  for (t in 1:n.occasions){
    qgamma[t] <- 1-gamma[t]
  }
  cprob[1] <- gamma[1]
  for (t in 2:n.occasions){
    cprob[t] <- gamma[t] * prod(qgamma[1:(t-1)])
  } 
  psi <- sum(cprob[]) #Inclusion probability
  for (t in 1:n.occasions){
    b[t] <- cprob[t] / psi #Entry probability
  } 
  for (i in 1:M){
    recruit[i,1] <- z[i,1]
    for (t in 2:n.occasions){
      recruit[i,t] <- (1-z[i,t-1]) * z[i,t]
    } 
  } 
  for (t in 1:n.occasions){
    N[t] <- sum(z[1:M,t]) #Actual population size
    B[t] <- sum(recruit[1:M,t]) #Number of entries
  } 
  for (i in 1:M){
    Nind[i] <- sum(z[i,1:n.occasions])
    Nalive[i] <- 1-equals(Nind[i], 0)
  } 
  Nsuper <- sum(Nalive[]) #Size of superpopulation
}

leis <- as.matrix(read.table("/Users/shelbie.ishimaru/Documents/GitHub/BayesianPopulationAnalysis_Learning/leisleri.txt", sep = " ", header = FALSE))
nz <- 300
CH.aug <- rbind(leis, matrix(0, ncol = dim(leis)[2], nrow = nz))

#Bundle data
jags.data <- list(y = CH.aug, n.occasions = dim(CH.aug)[2], M = dim(CH.aug)[1])

#Initial values
#Good initial values for the latent state z are needed to run the model in JAGS. The simplest option that works is to give just a matrix with a 1 at all places.
z.init <- CH.aug
z.init[z.init==0] <- 1
inits <- function(){list(mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1), sigma = runif(1, 0, 1), z = z.init)} 

#Parameters monitored
params <- c("psi", "mean.p", "sigma2", "mean.phi", "N", "Nsuper", "b", "B")

#MCMC settings
ni <- 10000
nt <- 6
nb <- 5000
nc <- 3

#Call JAGS from R
n1 <- jags(data  = jags.data,
               inits = inits,
               parameters.to.save = params,
               model.file = jags.js.tempran.txt,
               n.chains = nc,
               n.thin= nt,
               n.iter = ni,
               n.burnin = nb)

print(n1, digits = 3)

k<-mcmcplots::as.mcmc.rjags(n1)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#Code to produce Fig. 10.6
#Calculate per-capita recruitment
T <- dim(leis)[2]
f <- matrix(NA, ncol = T, nrow = length(n1$BUGSoutput$sims.list$B[,1]))
for (t in 1:(T-1)){
  f[,t] <- n1$BUGSoutput$sims.list$B[,t+1] / n1$BUGSoutput$sims.list$N[,t+1]
}

n.lower <- n.upper <- f.lower <- f.upper <- f.mean <- numeric()

for (t in 1:T){
  n.lower[t] <- quantile(n1$BUGSoutput$sims.list$N[,t], 0.025)
  n.upper[t] <- quantile(n1$BUGSoutput$sims.list$N[,t], 0.975)
}
for (t in 1:(T-1)){
  f.lower[t] <- quantile(f[,t], 0.025)
  f.upper[t] <- quantile(f[,t], 0.975)
  f.mean[t] <- mean(f[,t])
}
par(mfrow = c(1, 2))
plot(n1$BUGSoutput$mean$N, type = "b", pch = 19, ylab = "Population size", xlab = "", cex = 1.5, ylim = c(10, max(n.upper)), xaxt= "none")
axis(1, at = seq(1, T, 2), labels = seq(1990, 2008, 2))
axis(1, at = 1:T, labels = rep("", T), tcl = -0.25)
axis(2, las = 1)
segments(1:T, n.lower, 1:T, n.upper)
plot(f.mean, type = "b", pch = 19, ylab = "Local per capita recruitment", xlab = "", axes = F, cex = 1.5, ylim = c(0, 0.8), xaxt= "none")
axis(1, at = seq(1, (T-1), 2), labels = seq(1991, 2008, 2))
axis(1, at = 1:(T-1), labels = rep("", T-1), tcl = -0.25)
axis(2, las = 1)
segments(1:(T-1), f.lower, 1:(T-1), f.upper)
mtext("Figure 10.6", side= 3, line= -1.5, outer= T) #adding main title to multiplot
#-------------------------------------------------------------------------------