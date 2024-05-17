# Bayesian Population Analysis using WinBUGS
# Chapter 13: Estimation of Occupancy and Species Distributions from Detection/Nondetection Data in Metapopulation Designs Using Site-Occupancy Models

library(R2jags) #to run JAGS
library(shinystan) #to run shiny stan
library(tidyverse) #to utilize pipe operators

# 13.2 What Happens When p<1 and Constant and p is Not Accounted for in a Species Distribution Model?
#-------------------------------------------------------------------------------
nreps <- 10^5 #No. replicates, leave at 10^5 to produce figure 13.1
estimates <- array(NA, dim = c(nreps, 2)) #Array to contain the estimates
R <- 250 #No. sites

for (i in 1:nreps) {
  cat(i, "\n"); flush.console()
  x <- runif(R, 0, 10) #choose covariate values
  state<-rbinom(n = R, size = 1, prob = plogis(-3 + 1 * x)) #Occ. state
  obs <- rbinom(n = R, size = 1, prob = 0.6) * state #Observations
  fm <- glm(obs~x, family = binomial)
  estimates[i,] <- fm$coef
}

par(mfrow = c(3, 1))
hist(estimates[,1], col = "gray", nclass = 50, main = "", xlab = "Intercept estimates", las = 1, ylab = "Density", freq = FALSE)
abline(v = -3, col = "red", lwd = 3) #Truth
hist(estimates[,2], col = "gray", nclass = 50, main = "", xlab = "Slope estimates", xlim = c(0,1), las = 1, ylab = "Density", freq = FALSE)
abline(v = 1, col = "red", lwd = 3) #Truth

plot(1:10, plogis(estimates[1,1] + estimates[1,2] * (1:10)), col ="gray", lwd = 1, ylab = "Occupancy probability", xlab = "Covariate value", type = "l", ylim = c(0, 1), frame.plot = FALSE, las = 1)
samp <- sample(1:nreps, 1000)
for (i in samp){
  lines(1:10, plogis(estimates[i,1] + estimates[i,2] * (1:10)), col = "gray", lwd = 1, type = "l")
} #gray line= random sample of 1000 estimated regression lines
lines(1:10, plogis(-3 + 1 * (1:10)), col = "red", lwd = 3, type = "l") #truth
mtext("Figure 13.1", side= 3, line= -1.5, outer= T) #adding main title
#-------------------------------------------------------------------------------

# 13.3 Generation and Analysis of Simulated Data for Single-Season Occupancy
# 13.3.1 The Simplest Possible Site-Occupancy Model
#-------------------------------------------------------------------------------
#Select sample sizes (spatial and temporal replication)
R <- 200
T <- 3

#Determine process parameters
psi <- 0.8 #Occupancy probability
p <- 0.5 #Detection probability

#Create structure to contain counts
y <- matrix(NA, nrow = R, ncol = T)

#Ecological process: Sample true occurrence (z, yes/no) from a Bernoulli (occurrence probability = psi)
z <- rbinom(n = R, size = 1, prob = psi) #Latent occurrence state

#Observation process: Sample detection/nondetection observations from a Bernoulli(with p) if z=1
for (j in 1:T){
  y[,j] <- rbinom(n = R, size = 1, prob = z * p)
}

#Look at truth and at our imperfect observations
sum(z) #Realized occupancy among 200 surveyed sites
sum(apply(y, 1, max)) #Observed occupancy

#Specify model in BUGS language
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  #Priors
  psi ~ dunif(0, 1)
  p ~ dunif(0, 1)
  
  #Likelihood
  #Ecological model for true occurrence
  for (i in 1:R) {
    z[i] ~ dbern(psi)
    p.eff[i] <- z[i] * p
    #Observation model for replicated detection/nondetection observations
    for (j in 1:T) {
      y[i,j] ~ dbern(p.eff[i])
    }
  }
  
  #Derived quantities
  occ.fs <- sum(z[]) #Number of occupied sites among the 200
}

#Bundle data
jags.data <- list(y = y, R = nrow(y), T = ncol(y))

#Initial values
zst <- apply(y, 1, max) #Observed occurrence as starting values for z
inits <- function() list(z = zst)

#Parameters monitored
params <- c("psi", "p", "occ.fs")

#MCMC settings
ni <- 1200
nt <- 2
nb <- 200
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

# 13.3.2 Site-Occupancy Models with Covariates
#-------------------------------------------------------------------------------
#Define function for generating species distribution data
data.fn <- function(R = 200, T = 3, xmin = -1, xmax = 1, alpha.psi = -1, beta.psi = 3, alpha.p = 1, beta.p = -3) {
  y <- array(dim = c(R, T)) #Array for counts
  
  #Ecological process
  #Covariate values
  X <- sort(runif(n = R, min = xmin, max = xmax))
  
  #Relationship expected occurrence – covariate
  psi <- plogis(alpha.psi + beta.psi * X) #Apply inverse logit
  
  #Add Bernoulli noise: draw occurrence indicator z from Bernoulli(psi)
  z <- rbinom(n = R, size = 1, prob = psi)
  occ.fs <- sum(z) #Finite-sample occupancy (see Royle and Kéry, 2007)
  
  #Observation process
  #Relationship detection prob – covariate
  p <- plogis(alpha.p + beta.p * X)
  
  #Make a 'census'
  p.eff <- z * p
  for (i in 1:T){
    y[,i] <- rbinom(n = R, size = 1, prob = p.eff)
  }
  
  #Naïve regression
  naive.pred <- plogis(predict(glm(apply(y, 1, max) ~ X + I(X^2), family = binomial)))
  
  #Plot features of the simulated system
  par(mfrow = c(2, 2))
  plot(X, psi, main = "Expected occurrence", xlab = "Covariate", ylab = "Occupancy probability", las = 1, type = "l", col = "red", lwd = 3, frame.plot = FALSE)
  plot(X, z, main = "Realised (true) occurrence", xlab = "Covariate", ylab = "Occurrence", las = 1, frame.plot = FALSE, col = "red",)
  plot(X, p, ylim = c(0,1), main = "Detection probability", xlab = "Covariate", ylab = "p", type = "l", lwd = 3, col = "red", las = 1, frame.plot = FALSE)
  plot(X, naive.pred, main = "Detection/nondetection observations \n and conventional SDM", xlab = "Covariate", ylab = "Apparent occupancy", ylim = c(min(y), max(y)), type = "l", lwd = 3, lty = 2, col = "blue", las = 1, frame.plot = FALSE)
  points(rep(X, T), y)
  mtext("Figure 13.2", side= 3, line= -1.5, outer= T) #adding main title
  
  #Return stuff
  return(list(R = R, T = T, X = X, alpha.psi = alpha.psi, beta.psi = beta.psi, alpha.p = alpha.p , beta.p = beta.p, psi = psi, z = z, occ.fs = occ.fs, p = p, y = y))
}

sodata <- data.fn()
str(sodata) #Look at data
summary(glm(apply(y, 1, max) ~ X + I(X^2), family = binomial, data = sodata))

#Specify model in BUGS language
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  #Priors
  alpha.occ ~ dunif(-10, 10)
  beta.occ ~ dunif(-10, 10)
  alpha.p ~ dunif(-10, 10)
  beta.p ~ dunif(-10, 10)
  
  #Likelihood
  for (i in 1:R) {
    
    #True state model for the partially observed true state
    z[i] ~ dbern(psi[i]) #True occupancy z at site i
    logit(psi[i]) <- alpha.occ + beta.occ * X[i]
    for (j in 1:T) {
      #Observation model for the actual observations
      y[i,j] ~ dbern(p.eff[i,j]) #Detection-nondetection at i and j
      p.eff[i,j] <- z[i] * p[i,j]
      logit(p[i,j]) <- alpha.p + beta.p * X[i]
    } 
  } 
  
  #Derived quantities
  occ.fs <- sum(z[]) #Number of occupied sites among those studied
}

#Bundle data
jags.data <- list(y = sodata$y, X = sodata$X, R = nrow(sodata$y), T = ncol(sodata$y))

#Initial values
zst <- apply(sodata$y, 1, max) #Good inits for latent states essential
inits <- function(){list(z = zst, alpha.occ = runif(1, -3, 3), beta.occ = runif(1, -3, 3), alpha.p = runif(1, -3, 3), beta.p = runif(1, -3, 3))}

#Parameters monitored
params <- c("alpha.occ", "beta.occ", "alpha.p", "beta.p", "occ.fs")

#MCMC settings
ni <- 10000
nt <- 8
nb <- 2000
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

TRUTH <- c(sodata$alpha.psi, sodata$beta.psi, sodata$alpha.p, sodata$beta.p, sum(sodata$z))
print(cbind(TRUTH, out$summary[1:5, c(1,2,3,7)]), dig = 3)

sum(apply(sodata$y, 1, sum) > 0) #Apparent number of occupied sites

naive.pred <- plogis(predict(glm(apply(sodata$y, 1, max) ~ X + I(X^2), family = binomial, data = sodata)))
lin.pred2 <- out$BUGSoutput$mean$alpha.occ + out$BUGSoutput$mean$beta.occ * sodata$X
plot(sodata$X, sodata$psi, ylim = c(0, 1), main = "", ylab = "Occupancy probability", xlab = "Covariate", type = "l", lwd = 3, col = "red", las = 1, frame.plot = FALSE)
lines(sodata$X, naive.pred, ylim = c(0 ,1), type = "l", lty = 2, lwd = 3, col = "blue")
lines(sodata$X, plogis(lin.pred2), ylim = c(0, 1), type = "l", lty = 1, lwd = 2, col = "blue")
mtext("Figure 13.3", side= 3, line= -1.5, outer= T) #adding main title
#-------------------------------------------------------------------------------

# 13.4 Analysis of A Real Data Set: Single-Season Occupancy Model
#-------------------------------------------------------------------------------
#Read in the data
data <- read.table("/Users/shelbie.ishimaru/Documents/GitHub/BayesianPopulationAnalysis_Learning/bluebug.txt", header = TRUE)

#Collect the data into suitable structures
y <- as.matrix(data[,4:9]) #as.matrix essential for WinBUGS
y[y>1] <- 1 #Reduce counts to 0/1
edge <- data$forest_edge
dates <- as.matrix(data[,10:15])
hours <- as.matrix(data[,16:21])

#Standardize covariates
mean.date <- mean(dates, na.rm = TRUE)
sd.date <- sd(dates[!is.na(dates)])
DATES <- (dates-mean.date)/sd.date #Standardize date
DATES[is.na(DATES)] <- 0 #Impute zeroes (means)
mean.hour <- mean(hours, na.rm = TRUE)
sd.hour <- sd(hours[!is.na(hours)])
HOURS <- (hours-mean.hour)/sd.hour #Standardize hour
HOURS[is.na(HOURS)] <- 0 #Impute zeroes (means)

#Specify model in BUGS language
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  #Priors
  alpha.psi ~ dnorm(0, 0.01)
  beta.psi ~ dnorm(0, 0.01)
  alpha.p ~ dnorm(0, 0.01)
  beta1.p ~ dnorm(0, 0.01)
  beta2.p ~ dnorm(0, 0.01)
  beta3.p ~ dnorm(0, 0.01)
  beta4.p ~ dnorm(0, 0.01)
  
  #Likelihood
  #Ecological model for the partially observed true state
  for (i in 1:R) {
    z[i] ~ dbern(psi[i]) #True occurrence z at site i
    psi[i] <- 1 / (1 + exp(-lpsi.lim[i]))
    lpsi.lim[i] <- min(999, max(-999, lpsi[i]))
    lpsi[i] <- alpha.psi + beta.psi * edge[i]
    
    #Observation model for the observations
    for (j in 1:T) {
      y[i,j] ~ dbern(mu.p[i,j]) #Detection-nondetection at i and j
      mu.p[i,j] <- z[i] * p[i,j]
      p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))
      lp.lim[i,j] <- min(999, max(-999, lp[i,j]))
      lp[i,j] <- alpha.p + beta1.p * DATES[i,j] + beta2.p *
        pow(DATES[i,j], 2) + beta3.p * HOURS[i,j] + beta4.p *
        pow(HOURS[i,j], 2)
    } 
  } 
  
  #Derived quantities
  occ.fs <- sum(z[]) #Number of occupied sites
  mean.p <- exp(alpha.p) / (1 + exp(alpha.p)) #Average detection
}

#Bundle data
jags.data <- list(y = y, R = nrow(y), T = ncol(y), edge = edge, DATES = DATES, HOURS = HOURS)

#Initial values
zst <- apply(y, 1, max, na.rm = TRUE) # Good starting values crucial
inits <- function(){list(z = zst, alpha.psi=runif(1, -3, 3), alpha.p = runif(1, -3, 3))}

#Parameters monitored
params <- c("alpha.psi", "beta.psi", "mean.p", "occ.fs", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p")

#MCMC settings
ni <- 30000
nt <- 10
nb <- 20000
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

#Posterior distribution of the number of occupied woodpiles in actual sample
hist(out$BUGSoutput$sims.list$occ.fs, nclass = 30, col = "gray", main = "", xlab = "Number of occupied woodpiles (occ.fs)", xlim = c(9, 27))
abline(v = 10, lwd = 2) #The observed number
mtext("Figure 13.5", side= 3, line= -1.5, outer= T) #adding main title

Pstar <- array(NA, dim = c(out$BUGSoutput$n.sims, 10))
x <- cbind(rep(1, 3000), rep(2, 3000), rep(3, 3000), rep(4, 3000), rep(5, 3000), rep(6, 3000), rep(7, 3000), rep(8, 3000), rep(9, 3000), rep(10, 3000))
for (i in 1:out$BUGSoutput$n.sims) {
  for (j in 1:10){
    Pstar[i,j] <- 1 - (1 - out$BUGSoutput$sims.list$mean.p[i])^j
  } 
} 

boxplot(Pstar ~ x, col = "gray", las = 1, ylab = "Pstar", xlab = "Number of surveys", outline = FALSE)
abline(h = 0.95, lty = 2, lwd = 2)
mtext("Figure 13.6", side= 3, line= -1.5, outer= T) #adding main title

par(mfrow = c(2, 1))
hist(plogis(out$BUGSoutput$sims.list$alpha.psi), nclass = 40, col = "gray", main = "Forest interior", xlab = "Occupancy probability", xlim = c(0, 1))
hist(plogis(out$BUGSoutput$sims.list$alpha.psi+ out$BUGSoutput$sims.list$beta.psi), nclass = 40, col = "gray", main = "Forest edge", xlab = "Occupancy probability", xlim = c(0, 1))
mtext("Figure 13.7", side= 3, line= -1.5, outer= T) #adding main title

#Predict effect of time of day with uncertainty
mcmc.sample <- out$BUGSoutput$n.sims
original.date.pred <- seq(0, 60, length.out = 30)
original.hour.pred <- seq(180, 540, length.out = 30)
date.pred <- (original.date.pred - mean.date)/sd.date
hour.pred <- (original.hour.pred - mean.hour)/sd.hour
p.pred.date <- plogis(out$BUGSoutput$mean$alpha.p + out$BUGSoutput$mean$beta1.p * date.pred + out$BUGSoutput$mean$beta2.p * date.pred^2)
p.pred.hour <- plogis(out$BUGSoutput$mean$alpha.p + out$BUGSoutput$mean$beta3.p * hour.pred + out$BUGSoutput$mean$beta4.p * hour.pred^2 )
array.p.pred.hour <- array.p.pred.date <- array(NA, dim = c(length(hour.pred), mcmc.sample))

for (i in 1:mcmc.sample){
  array.p.pred.date[,i] <- plogis(out$BUGSoutput$sims.list$alpha.p[i] + out$BUGSoutput$sims.list$beta1.p[i] * date.pred + out$BUGSoutput$sims.list$beta2.p[i] * date.pred^2)
  array.p.pred.hour[,i] <- plogis(out$BUGSoutput$sims.list$alpha.p[i] + out$BUGSoutput$sims.list$beta3.p[i] * hour.pred + out$BUGSoutput$sims.list$beta4.p[i] * hour.pred^2)
}

#Plot for a subsample of MCMC draws
sub.set <- sort(sample(1:mcmc.sample, size = 200))

par(mfrow = c(2, 1))
plot(original.date.pred, p.pred.date, main = "", ylab = "Detection probability", xlab = "Date (1 = 1 July)", ylim = c(0, 1), type = "l", lwd = 3, frame.plot = FALSE)

for (i in sub.set){
  lines(original.date.pred, array.p.pred.date[,i], type = "l", lwd = 1, col = "gray")
}

lines(original.date.pred, p.pred.date, type = "l", lwd = 3, col = "blue")
plot(original.hour.pred, p.pred.hour, main = "", ylab = "Detection probability", xlab = "Time of day (mins after noon)", ylim = c(0, 1), type = "l", lwd = 3, frame.plot = FALSE)

for (i in sub.set){
  lines(original.hour.pred, array.p.pred.hour[,i], type = "l", lwd = 1, col = "gray")
}

lines(original.hour.pred, p.pred.hour, type = "l", lwd = 3, col = "blue")
mtext("Figure 13.8", side= 3, line= -1.5, outer= T) #adding main title
#-------------------------------------------------------------------------------

# 13.5 Dynamic (Multiseason) Site-Occupancy Models
# 13.5.1 Generation and Analysis of Simulated Data
#-------------------------------------------------------------------------------
data.fn <- function(R = 250, J = 3, K = 10, psi1 = 0.4, range.p = c(0.2, 0.4), range.phi = c(0.6, 0.8), range.gamma = c(0, 0.1)) {  #Function to simulate detection/nondetection data for dynamic site-occ model
  #R = 250, J = 3, K = 10, psi1 = 0.4, range.p = c(0.2, 0.4), range.phi = c(0.6, 0.8), range.gamma = c(0, 0.1)
  #Annual variation in probabilities of patch survival, colonization and
  #detection is specified by the bounds of a uniform distribution.
  #Function arguments:
  #R – Number of sites
  #J – Number of replicate surveys
  #K – Number of years
  #psi1 – occupancy probability in first year
  #range.p – bounds of uniform distribution from which annual p drawn
  #range.psi and range.gamma – same for survival and colonization probability
  
  # Set up some required arrays
  site <- 1:R					# Sites
  year <- 1:K					# Years
  psi <- rep(NA, K)				# Occupancy probability
  muZ <- z <- array(dim = c(R, K))	# Expected and realized occurrence
  y <- array(NA, dim = c(R, J, K))	# Detection histories
  
  # Determine initial occupancy and demographic parameters
  psi[1] <- psi1				# Initial occupancy probability
  p <- runif(n = K, min = range.p[1], max = range.p[2])
  phi <- runif(n = K-1, min = range.phi[1], max = range.phi[2])
  gamma <- runif(n = K-1, min = range.gamma[1], max = range.gamma[2])
  
  # Generate latent states of occurrence
  # First year
  z[,1] <- rbinom(R, 1, psi[1])		# Initial occupancy state
  # Later years
  for(i in 1:R){				# Loop over sites
    for(k in 2:K){				# Loop over years
      muZ[k] <- z[i, k-1]*phi[k-1] + (1-z[i, k-1])*gamma[k-1] # Prob for occ.
      z[i,k] <- rbinom(1, 1, muZ[k])
    }
  }
  
  # Plot realized occupancy
  plot(year, apply(z, 2, mean), type = "l", xlab = "Year", ylab = "Occupancy or Detection prob.", col = "red", xlim = c(0,K+1), ylim = c(0,1), lwd = 2, lty = 1, frame.plot = FALSE, las = 1)
  lines(year, p , type = "l", col = "red", lwd = 2, lty = 2)
  
  # Generate detection/nondetection data
  for(i in 1:R){
    for(k in 1:K){
      prob <- z[i,k] * p[k]
      for(j in 1:J){
        y[i,j,k] <- rbinom(1, 1, prob)
      }
    }
  }
  
  # Compute annual population occupancy
  for (k in 2:K){
    psi[k] <- psi[k-1]*phi[k-1] + (1-psi[k-1])*gamma[k-1]
  }
  
  # Plot apparent occupancy
  psi.app <- apply(apply(y, c(1,3), max), 2, mean)
  lines(year, psi.app, type = "l", col = "black", lwd = 2)
  text(0.85*K, 0.06, labels = "red solid - true occupancy\n red dashed - detection\n black - observed occupancy")
  
  # Return data
  return(list(R = R, J = J, K = K, psi = psi, psi.app = psi.app, z = z, phi = phi, gamma = gamma, p = p, y = y))
}

data <- data.fn(R = 250, J = 3, K = 10, psi1 = 0.6, range.p = c(0.1, 0.9), range.phi = c(0.7, 0.9), range.gamma = c(0.1, 0.5))

#attach(data)
#str(data)

#Specify model in BUGS language
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  # Specify priors
  psi1 ~ dunif(0, 1)
  for (k in 1:(nyear-1)){
    phi[k] ~ dunif(0, 1)
    gamma[k] ~ dunif(0, 1)
    p[k] ~ dunif(0, 1) 
  }
  p[nyear] ~ dunif(0, 1)
  
  # Ecological submodel: Define state conditional on parameters
  for (i in 1:nsite){
    z[i,1] ~ dbern(psi1)
    for (k in 2:nyear){
      muZ[i,k]<- z[i,k-1]*phi[k-1] + (1-z[i,k-1])*gamma[k-1]
      z[i,k] ~ dbern(muZ[i,k])
    } #k
  } #i
  
  # Observation model
  for (i in 1:nsite){
    for (j in 1:nrep){
      for (k in 1:nyear){
        muy[i,j,k] <- z[i,k]*p[k]
        y[i,j,k] ~ dbern(muy[i,j,k])
      } #k
    } #j
  } #i
  
  # Derived parameters: Sample and population occupancy, growth rate and turnover
  psi[1] <- psi1
  n.occ[1]<-sum(z[1:nsite,1])
  for (k in 2:nyear){
    psi[k] <- psi[k-1]*phi[k-1] + (1-psi[k-1])*gamma[k-1]
    n.occ[k] <- sum(z[1:nsite,k])
    growthr[k-1] <- psi[k]/psi[k-1]                         # originally we had growthr[k]. JAGS seem to dislike vectoring going from 2..K.
    turnover[k-1] <- (1 - psi[k-1]) * gamma[k-1]/psi[k]
  }
}

#Bundle data
jags.data <- list(y = y, nsite = dim(y)[1], nrep = dim(y)[2], nyear = dim(y)[3])

# Initial values
zst <- apply(y, c(1,3), max)	# Observed occurrence as inits for z
inits <- function(){ list(z = zst)}

# Parameters monitored
params <- c("psi", "phi", "gamma", "p", "n.occ", "growthr", "turnover")

# MCMC settings
ni <- 2500
nt <- 4
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

print(out, digits = 2)

k<-mcmcplots::as.mcmc.rjags(out)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

psiall <- paste("psi[", 1:K, "]", sep="")
print(cbind(data$psi, out$BUGSoutput$summary[psiall, c(1, 2, 3, 7)]), dig = 3)
phiall <- paste("phi[", 1:(K-1), "]", sep="")
print(cbind(data$phi, out$BUGSoutput$summary[phiall, c(1, 2, 3, 7)]), dig = 3)
gammaall <- paste("gamma[", 1:(K-1), "]", sep="")
print(cbind(data$gamma, out$BUGSoutput$summary[gammaall, c(1, 2, 3, 7)]), dig = 3)
pall <- paste("p[", 1:K, "]", sep="")
print(cbind(data$p, out$BUGSoutput$summary[pall, c(1, 2, 3, 7)]), dig = 3)

plot(1:K, data$psi, type = "l", xlab = "Year", ylab = "Occupancy probability", col = "red", xlim = c(0,K+1), ylim = c(0,1), lwd = 2, lty = 1, frame.plot = FALSE, las = 1)
lines(1:K, data$psi.app, type = "l", col = "black", lwd = 2)
points(1:K, out$BUGSoutput$mean$psi, type = "l", col = "blue", lwd = 2)
segments(1:K, out$BUGSoutput$summary[psiall,3], 1:K, out$BUGSoutput$summary[psiall,7], col = "blue", lwd = 1) 
#-------------------------------------------------------------------------------

# 13.5.2 Dynamic Occupancy Modeling in a Real Data Set
#-------------------------------------------------------------------------------
#Read in the data and put it into 3D array
bdat <- read.table(file = "/Users/shelbie.ishimaru/Documents/GitHub/BayesianPopulationAnalysis_Learning/burnet.txt", header = T)
str(bdat)
y <- array(NA, dim = c(95, 2, 7)) #95 sites, 2 reps, 7 days
for (i in 1:7){
  sel.rows <- bdat$day == i
  y[,,i] <- as.matrix(bdat)[sel.rows, 3:4]
}
str(y)

#Convert counts to detection/nondetection data
y[y>0] <- 1

#Look at the number of sites with detections for each day
tmp <- apply(y, c(1,3), max, na.rm = TRUE)
tmp[tmp == "-Inf"] <- NA
apply(tmp, 2, sum, na.rm = TRUE)

#Bundle data
jags.data <- list(y = y, nsite = dim(y)[1], nrep = dim(y)[2], nyear = dim(y)[3])

#Initial values
inits <- function(){ list(z = apply(y, c(1, 3), max))}

#Parameters monitored
params <- c("psi", "phi", "gamma", "p", "n.occ", "growthr", "turnover")

#MCMC settings
ni <- 5000
nt <- 4
nb <- 1000
nc <- 3

#Call JAGS from R
out1 <- jags(data  = jags.data,
            inits = inits,
            parameters.to.save = params,
            model.file = jags.model.txt,
            n.chains = nc,
            n.thin= nt,
            n.iter = ni,
            n.burnin = nb)

print(out1, digits = 2)

k<-mcmcplots::as.mcmc.rjags(out1)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#Specify model in BUGS language
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  #Specify priors
  psi1 ~ dunif(0, 1)
  for (k in 1:(nyear-1)){
    phi[k] ~ dunif(0, 1)
    gamma[k] ~ dunif(0, 1)
  }
  p ~ dunif(0, 1)
  
  #Both models at once
  for (i in 1:nsite){
    z[i,1] ~ dbern(psi1) # State model 1: Initial state
    for (k in 2:nyear){ # State model 2: State dynamics
      muZ[i,k] <- z[i,k-1]*phi[k-1] + (1-z[i,k-1])*gamma[k-1]
      z[i,k] ~ dbern(muZ[i,k])
      
      #Observation model
      muy[i,k] <- z[i,k]*p
      y[i,k] ~ dbin(muy[i,k], 2)
    } 
  } 
  
  #Derived parameters: Sample and population occupancy, growth rate and turnover
  psi[1] <- psi1
  n.occ[1] <- sum(z[1:nsite,1])
  for (k in 2:nyear){
    psi[k] <- psi[k-1]*phi[k-1] + (1-psi[k-1])*gamma[k-1]
    n.occ[k] <- sum(z[1:nsite,k])
    growthr[k] <- psi[k]/psi[k-1]
    turnover[k-1] <- (1 - psi[k-1]) * gamma[k-1]/psi[k]
  }
}

#Aggregate detections over reps within a day and bundle data
yy <- apply(y, c(1, 3), sum, na.rm = TRUE)
jags.data <- list(y = yy, nsite = dim(yy)[1], nyear = dim(yy)[2])

#Initial values
inits <- function(){list(z = apply(y, c(1, 3), max))}

#Parameters monitored
params <- c("psi", "phi", "gamma", "p", "n.occ", "growthr", "turnover")

#MCMC settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

#Call JAGS from R
out2 <- jags(data  = jags.data,
             inits = inits,
             parameters.to.save = params,
             model.file = jags.model.txt,
             n.chains = nc,
             n.thin= nt,
             n.iter = ni,
             n.burnin = nb)

print(out2, digits = 2)

k<-mcmcplots::as.mcmc.rjags(out2)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

DAY <- cbind(rep(1, out2$BUGSoutput$n.sims), rep(2, out2$BUGSoutput$n.sims), rep(3,out2$BUGSoutput$n.sims), rep(4, out2$BUGSoutput$n.sims), rep(5, out2$BUGSoutput$n.sims), rep(6,out2$BUGSoutput$n.sims), rep(7, out2$BUGSoutput$n.sims))
boxplot(out2$BUGSoutput$sims.list$psi ~ DAY, col = "gray", ylab = "Occupancy probability", xlab = "Day of survey", las = 1, frame.plot = FALSE)
mtext("Figure 13.12", side= 3, line= -1.5, outer= T) #adding main title

apply(apply(y, c(1,3), max), 2, function(x){sum(!is.na(x))}) #comparing sample sizes for eahc day
#-------------------------------------------------------------------------------

# 13.6 Multistate Occupancy Models
#-------------------------------------------------------------------------------
owls <- read.table("/Users/shelbie.ishimaru/Documents/GitHub/BayesianPopulationAnalysis_Learning/owls.txt", header = TRUE)
str(owls)

#Specify model in BUGS language
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  #Priors
  p2 ~ dunif(0, 1)
  psi ~ dunif(0, 1)
  r ~ dunif(0, 1)
  for (i in 1:3) {
    beta[i] ~ dgamma(1, 1) #Induce Dirichlet prior
    p3[i] <- beta[i]/sum(beta[])
  }
  
  #Define state vector
  for (s in 1:R){
    phi[s,1] <- 1 - psi #Prob. of nonoccupation
    phi[s,2] <- psi * (1 - r) #Prob. of occupancy without repro
    phi[s,3] <- psi * r #Prob. of occupancy and repro
  }
  
  #Define observation matrix
  #Order of indices: true state, time, observed state
  for (t in 1:T){
    p[1,t,1] <- 1
    p[1,t,2] <- 0
    p[1,t,3] <- 0
    p[2,t,1] <- 1-p2
    p[2,t,2] <- p2
    p[2,t,3] <- 0
    p[3,t,1] <- p3[1]
    p[3,t,2] <- p3[2]
    p[3,t,3] <- p3[3]
  }
  
  #State-space likelihood
  #State equation: model of true states (z)
  for (s in 1:R){
    z[s] ~ dcat(phi[s,])
  }
  
  #Observation equation
  for (s in 1:R){
    for (t in 1:T){
      y[s,t] ~ dcat(p[z[s],t,])
    } 
  } 
  
  #Derived quantities
  for (s in 1:R){
    occ1[s] <- equals(z[s], 1)
    occ2[s] <- equals(z[s], 2)
    occ3[s] <- equals(z[s], 3)
  }
  n.occ[1] <- sum(occ1[]) #Sites in state 1
  n.occ[2] <- sum(occ2[]) #Sites in state 2
  n.occ[3] <- sum(occ3[]) #Sites in state 3
}

#Bundle data
y <- as.matrix(owls[, 2:6])
y <- y + 1
jags.data <- list(y = y, R = dim(y)[1], T = dim(y)[2])

#Initial values
zst <- apply(y, 1, max, na.rm = TRUE)
zst[zst == "-Inf"] <- 1
inits <- function(){list(z = zst)}

#Parameters monitored
params <- c("p2", "p3", "r", "psi", "n.occ") #Might want to add "z"

#MCMC settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

#Call JAGS from R
out1 <- jags(data  = jags.data,
             inits = inits,
             parameters.to.save = params,
             model.file = jags.model.txt,
             n.chains = nc,
             n.thin= nt,
             n.iter = ni,
             n.burnin = nb)

print(out1, digits = 2)

k<-mcmcplots::as.mcmc.rjags(out1)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#Specify model in BUGS language
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  #Priors
  psi ~ dunif(0, 1)
  r ~ dunif(0,1 )
  for (t in 1:T){
    p2[t] ~ dunif(0, 1)
    for (i in 1:3) {
      beta[i,t] ~ dgamma(1, 1) #Induce Dirichlet prior
      p3[i,t] <- beta[i,t]/sum(beta[,t])
    } 
  } 
  
  #Define state vector
  for (s in 1:R){
    phi[s,1] <- 1 - psi #Prob. of nonoccupation
    phi[s,2] <- psi * (1 - r) #Prob. of occupancy without repro.
    phi[s,3] <- psi * r #Prob. of occupancy and repro.
  }
  
  #Define observation matrix
  #Order of indices: true state, time, observed state
  for (t in 1:T){
    p[1,t,1] <- 1
    p[1,t,2] <- 0
    p[1,t,3] <- 0
    p[2,t,1] <- 1-p2[t]
    p[2,t,2] <- p2[t]
    p[2,t,3] <- 0
    p[3,t,1] <- p3[1,t]
    p[3,t,2] <- p3[2,t]
    p[3,t,3] <- p3[3,t]
  }
  
  #State-space likelihood
  #State equation: model of true states (z)
  for (s in 1:R){
    z[s] ~ dcat(phi[s,])
  }
  
  #Observation equation
  for (s in 1:R){
    for (t in 1:T){
      y[s,t] ~ dcat(p[z[s],t,])
    } 
  } 
  
  #Derived quantities
  for (s in 1:R){
    occ1[s] <- equals(z[s], 1)
    occ2[s] <- equals(z[s], 2)
    occ3[s] <- equals(z[s], 3)
  }
  n.occ[1] <- sum(occ1[]) #Sites in state 1
  n.occ[2] <- sum(occ2[]) #Sites in state 2
  n.occ[3] <- sum(occ3[]) #Sites in state 3
}

#Call JAGS from R
out2 <- jags(data  = jags.data,
             inits = inits,
             parameters.to.save = params,
             model.file = jags.model.txt,
             n.chains = nc,
             n.thin= nt,
             n.iter = ni,
             n.burnin = nb)

print(out2, digits = 2)

k<-mcmcplots::as.mcmc.rjags(out2)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan
#-------------------------------------------------------------------------------