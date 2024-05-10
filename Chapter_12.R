# Bayesian Population Analysis using WinBUGS
# Chapter 12: Estimation of Abundance from Counts in Metapopulation Design Using the Binomial Mixture Model

library(R2jags) #to run JAGS
library(shinystan) #to run shiny stan
library(tidyverse) #to utilize pipe operators

# 12.2 Generation and Analysis of Simulated Data
# 12.2.1 The Simplest Case with Constant Parameters
#-------------------------------------------------------------------------------
#Determine sample sizes (spatial and temporal replication)
R <- 200
T <- 3

#Create structure to contain counts
y <- array(dim = c(R, T))

#Sample abundance from a Poisson(lambda = 2)
N <- rpois(n = R, lambda = 2)

#Sample counts from a Binomial (N, p = 0.5)
for (j in 1:T){
  y[,j] <- rbinom(n = R, size = N, prob = 0.5)
}

#Look at realization of biological and observation processes
cbind(N, y)

#Specify model in BUGS language
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  #Priors
  lambda ~ dgamma(0.005, 0.005) #Standard vague prior for lambda
  #lambda ~ dunif(0, 10) # Other possibility
  p ~ dunif(0, 1)
  
  #Likelihood
  #Biological model for true abundance
  for (i in 1:R) {
    N[i] ~ dpois(lambda)
    #Observation model for replicated counts
    for (j in 1:T) {
      y[i,j] ~ dbin(p, N[i])
    } 
  } 
}

#Bundle data
jags.data <- list(y = y, R = nrow(y), T = ncol(y))

#Initial values
Nst <- apply(y, 1, max) + 1 #This line is important
inits <- function() list(N = Nst)

#Parameters monitored
params <- c("lambda", "p")

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

print(out, digits = 2)

k<-mcmcplots::as.mcmc.rjags(out)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan
#-------------------------------------------------------------------------------

# 12.2.2 Introducing Covariates
#-------------------------------------------------------------------------------
#Define function for generating binom-mix model data
data.fn <- function(R = 200, T = 3, xmin = -1, xmax = 1, alpha0 = 1, alpha1 = 3, beta0 = 0, beta1 = -5){
  #R: number of sites at which counts were made (= number of spatial reps)
  #T: number of times that counts were made at each site
  #(= number of temporal reps)
  #xmin, xmax: define range of the covariate X
  #alpha0 and alpha1: intercept and slope of log-linear regression
  #relating abundance to the site covariate A
  #beta0 and beta1: intercept and slope of logistic-linear regression
  #of detection probability on A
  y <- array(dim = c(R, T)) # Array for counts
  
  #Ecological process
  #Covariate values: sort for ease of presentation
  X <- sort(runif(n = R, min = xmin, max = xmax))
  #Relationship expected abundance – covariate
  lam <- exp(alpha0 + alpha1 * X)
  
  #Add Poisson noise: draw N from Poisson(lambda)
  N <- rpois(n = R, lambda = lam)
  table(N) #Distribution of abundances across sites
  sum(N > 0) / R #Empirical occupancy
  totalN <- sum(N) ; totalN
  
  #Observation process
  #Relationship detection prob – covariate
  p <- plogis(beta0 + beta1 * X)
  #Make a 'census' (i.e., go out and count things)
  for (i in 1:T){
    y[,i] <- rbinom(n = R, size = N, prob = p)
  }
  
  #Naïve regression
  naive.pred <- exp(predict(glm(apply(y, 1, max) ~ X + I(X^2), family = poisson)))
  
  #Plot features of the simulated system
  par(mfrow = c(2, 2))
  plot(X, lam, main = "Expected abundance", xlab = "Covariate", ylab = "lambda", las = 1, type = "l", col = "red", lwd = 3, frame.plot = FALSE)
  plot(X, N, main = "Realised abundance", xlab = "Covariate", ylab = "N", las = 1, frame.plot = FALSE, col = "red", cex = 1.2)
  plot(X, p, ylim = c(0, 1), main = "Detection probability", xlab = "Covariate", ylab = "p", type = "l", col = "red", lwd = 3, las = 1, frame.plot = FALSE)
  plot(X, naive.pred, main = "Actual counts \n and naïve regression", xlab = "Covariate", ylab = "Relative abundance", ylim = c(min(y), max(y)), type = "l", lty = 2, lwd = 4, col = "blue", las = 1, frame.plot = FALSE)
  points(rep(X, T), y, col = "black", cex = 1.2)
  
  #Return stuff
  return(list(R = R, T = T, X = X, alpha0 = alpha0, alpha1 = alpha1, beta0 = beta0, beta1 = beta1, lam = lam, N = N, totalN = totalN, p = p, y = y))
}

data <- data.fn()
str(data)

#Specify model in BUGS language
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  #Priors
  alpha0 ~ dunif(-10, 10)
  alpha1 ~ dunif(-10, 10)
  beta0 ~ dunif(-10, 10)
  beta1 ~ dunif(-10, 10)
  
  #Likelihood
  #Ecological model for true abundance
  for (i in 1:R){
    N[i] ~ dpois(lambda[i])
    log(lambda[i]) <- alpha0 + alpha1 * X[i]
    
    #Observation model for replicated counts
    for (j in 1:T){
      y[i,j] ~ dbin(p[i,j], N[i])
      p[i,j] <- exp(lp[i,j])/(1+exp(lp[i,j]))
      lp[i,j] <- beta0 + beta1 * X[i]
    } 
  } 
  #Derived quantities
  totalN <- sum(N[])
}
  
#Bundle data
y <- data$y
jags.data <- list(y = y, R = nrow(y), T = ncol(y), X = data$X)

#Initial values
Nst <- apply(y, 1, max) + 1 #Important to give good inits for latent N
inits <- function() list(N = Nst, alpha0 = runif(1, -1, 1), alpha1 = runif(1, -1, 1), beta0 = runif(1, -1, 1), beta1 = runif(1, -1, 1))

#Parameters monitored
params <- c("totalN", "alpha0", "alpha1", "beta0", "beta1")

#MCMC settings
ni <- 22000
nt <- 20
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

print(out, digits = 3)

k<-mcmcplots::as.mcmc.rjags(out)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#Plot posteriors
par(mfrow = c(2, 2))
hist(out$sims.list$alpha0, col = "gray", main = "", xlab = "alpha0", las = 1)
abline(v = data$alpha0, lwd = 3, col = "red")
hist(out$sims.list$alpha1, col = "gray", main = "", xlab = "alpha1", las = 1)
abline(v = data$alpha1, lwd = 3, col = "red")
hist(out$sims.list$beta0, col = "gray", main = "", xlab = "beta0", las=1)
abline(v = data$beta0, lwd = 3, col = "red")
hist(out$sims.list$beta1, col = "gray", main = "", xlab = "beta1", las = 1)
abline(v = data$beta1, lwd = 3, col = "red")

#Plot predicted covariate relationship with abundance
plot(data$X, data$N, main = "", xlab = "Covariate", ylab = "Abundance", las = 1, ylim = c(0, max(data$N)), frame.plot = FALSE)
lines(data$X, data$lam, type = "l", col = "red", lwd = 3)
GLM.pred <- exp(predict(glm(apply(data$y, 1, max) ~ X + I(X^2), family = poisson, data = data)))
lines(data$X, GLM.pred, type = "l", lty = 2, col = "blue", lwd = 3)
Nmix.pred <- exp(out$mean$alpha0 + out$mean$alpha1 * data$X)
points(data$X, Nmix.pred, type = "l", col = "blue", lwd = 3)
#-------------------------------------------------------------------------------

# 12.3 Analysis of Real Data: Open-Population Binomial Mixture Models
#-------------------------------------------------------------------------------
#Get the data and put them into 3D array
bdat <- read.table("fritillary.txt", header = TRUE)
y <- array(NA, dim = c(95, 2, 7)) #95 sites, 2 reps, 7 days
for(k in 1:7){
  sel.rows <- bdat$day == k
  y[,,k] <- as.matrix(bdat)[sel.rows, 3:4]
}
y #Look at data set in 3D layout
str(y)

#Have a look at raw data
day.max <- apply(y, c(1, 3), max, na.rm = TRUE) # Max count each site and day
day.max

site.max <- apply(day.max, 1, max, na.rm = TRUE) # Max count each site
site.max

table(site.max) #Frequency distribution of max counts
plot(table(site.max))
table(site.max>0) #Observed occupancy is only 56%

#Sum of observed max as conventional estimator of total abundance
max1 <- apply(y, c(1, 3), max)
obs.max.sum <- apply(max1, 2, sum, na.rm = TRUE)
#-------------------------------------------------------------------------------

# 12.3.1 Simple Poisson Model
#-------------------------------------------------------------------------------
#Specify model in BUGS language
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  #Priors
  for (k in 1:7){
    alpha.lam[k] ~ dnorm(0, 0.01)
    p[k] ~ dunif(0, 1)
  }
  
  #Likelihood
  #Ecological model for true abundance
  for (k in 1:7){ #Loop over days (7)
    lambda[k] <- exp(alpha.lam[k])
    for (i in 1:R){ #Loop over R sites (95)
      N[i,k] ~ dpois(lambda[k]) #Abundance
      
      #Observation model for replicated counts
      for (j in 1:T){ #Loop over temporal reps (2)
        y[i,j,k] ~ dbin(p[k], N[i,k]) #Detection
        
        #Assess model fit using Chi-squared discrepancy
        #Compute fit statistic E for observed data
        eval[i,j,k] <- p[k] * N[i,k] #Expected values
        E[i,j,k] <- pow((y[i,j,k] - eval[i,j,k]),2) / (eval[i,j,k] + 0.5)
        
        #Generate replicate data and compute fit stats for them
        y.new[i,j,k] ~ dbin(p[k], N[i,k])
        E.new[i,j,k] <- pow((y.new[i,j,k] - eval[i,j,k]),2) / (eval[i,j,k] + 0.5)
      } 
    } 
  } 
  #Derived and other quantities
  for (k in 1:7){
    totalN[k] <- sum(N[,k]) #Total pop. size across all sites
    mean.abundance[k] <- exp(alpha.lam[k])
  }
  fit <- sum(E[,,])
  fit.new <- sum(E.new[,,])
}

#Bundle data
R = nrow(y)
T = ncol(y)
jags.data <- list(y = y, R = R, T = T)

#Initial values
Nst <- apply(y, c(1, 3), max) + 1
Nst[is.na(Nst)] <- 1
inits <- function(){list(N = Nst, alpha.lam = runif(7, -1, 1))}

#Parameters monitored
params <- c("totalN", "mean.abundance", "alpha.lam", "p", "fit", "fit.new")

#MCMC settings
ni <- 10000
nt <- 8
nb <- 2000
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

print(out0, digits = 3)

k<-mcmcplots::as.mcmc.rjags(out0)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#Evaluation of fit
plot(out0$sims.list$fit, out0$sims.list$fit.new, main = "", xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", frame.plot = FALSE)
abline(0, 1, lwd = 2, col = "black")

mean(out0$sims.list$fit.new > out0$sims.list$fit)

mean(out0$mean$fit) / mean(out0$mean$fit.new)
#-------------------------------------------------------------------------------

# 12.3.2 Zero-Inflated Poisson Binomial Mixture Model
#-------------------------------------------------------------------------------
#Specify model in BUGS language
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  #Priors
  omega ~ dunif(0, 1)
  for (k in 1:7){
    alpha.lam[k] ~ dnorm(0, 0.01)
    p[k] ~ dunif(0, 1)
  }
  
  #Likelihood
  #Ecological model for true abundance
  for (i in 1:R){ #Loop over R sites (95)
    z[i] ~ dbern(omega) #Latent suitability state
    for (k in 1:7){ #Loop over survey periods (seasons)
      N[i,k] ~ dpois(lam.eff[i,k]) #Latent abundance state
      lam.eff[i,k] <- z[i] * lambda[i,k]
      log(lambda[i,k]) <- alpha.lam[k]
      
      #Observation model for replicated counts
      for (j in 1:T){ #Loop over temporal reps (2)
        y[i,j,k] ~ dbin(p[k], N[i,k]) #Detection
        
        #Assess model fit using Chi-squared discrepancy
        #Compute fit statistic for observed data
        eval[i,j,k] <- p[k] * N[i,k]
        E[i,j,k] <- ((y[i,j,k] - eval[i,j,k])*(y[i,j,k] - eval[i,j,k])) / (eval[i,j,k] + 0.5)
        
        #Generate replicate data and compute fit stats for them
        y.new[i,j,k] ~ dbin(p[k], N[i,k])
        E.new[i,j,k] <- ((y.new[i,j,k] - eval[i,j,k]) * (y.new[i,j,k] - eval[i,j,k])) / (eval[i,j,k]+0.5)
      } 
    } 
  } 
  
  #Derived and other quantities
  for (k in 1:7){
    
    #Estimate total pop. size across all sites
    totalN[k] <- sum(N[,k])
    mean.abundance[k] <- exp(alpha.lam[k])
  }
  fit <- sum(E[,,])
  fit.new <- sum(E.new[,,])
}
  
#Bundle data
R = nrow(y)
T = ncol(y)
jags.data <- list(y = y, R = R, T = T)

#Initial values
Nst <- apply(y, c(1, 3), max) + 1
Nst[is.na(Nst)] <- 1
inits <- function(){list(N = Nst, alpha.lam = runif(7, -1, 1))}

#Parameters monitored
params <- c("omega", "totalN", "alpha.lam", "p", "mean.abundance",
            "fit", "fit.new")

#MCMC settings
ni <- 30000
nt <- 15
nb <- 15000
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

print(out1, digits = 3)

k<-mcmcplots::as.mcmc.rjags(out1)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#Evaluation of fit
plot(out1$sims.list$fit, out1$sims.list$fit.new, main = "", xlab =
       "Discrepancy actual data", ylab = "Discrepancy replicate data",
     frame.plot = FALSE)
abline(0, 1, lwd = 2, col = "black")
mean(out1$sims.list$fit.new > out1$sims.list$fit)

mean(out1$mean$fit) / mean(out1$mean$fit.new)
#-------------------------------------------------------------------------------

# 12.3.3 Binomial Mixture Model with Overdispersion in Both Abundace and Detection
#-------------------------------------------------------------------------------
#Specify model in BUGS language
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  #Priors
  for (k in 1:7){
    alpha.lam[k] ~ dnorm(0, 0.1)
    beta[k] ~ dnorm(0, 0.1)
  }
  
  #Abundance site and detection site-by-day random effects
  for (i in 1:R){
    eps[i] ~ dnorm(0, tau.lam) # Abundance noise
  }
  tau.lam <- 1 / (sd.lam * sd.lam)
  sd.lam ~ dunif(0, 3)
  tau.p <- 1 / (sd.p * sd.p)
  sd.p ~ dunif(0, 3)
  
  #Likelihood
  #Ecological model for true abundance
  for (i in 1:R){ #Loop over R sites (95)
    for (k in 1:7){ #Loop over days (7)
      N[i,k] ~ dpois(lambda[i,k]) #Abundance
      log(lambda[i,k]) <- alpha.lam[k] + eps[i]
      
      #Observation model for replicated counts
      for (j in 1:T){ #Loop over temporal reps (2)
        y[i,j,k] ~ dbin(p[i,j,k], N[i,k]) #Detection
        p[i,j,k] <- 1 / (1 + exp(-lp[i,j,k]))
        lp[i,j,k] ~ dnorm(beta[k], tau.p) #Random delta defined implicitly
        
        #Assess model fit using Chi-squared discrepancy
        #Compute fit statistic for observed data
        eval[i,j,k] <- p[i,j,k] * N[i,k]
        E[i,j,k] <- pow((y[i,j,k] - eval[i,j,k]),2) / (eval[i,j,k]+0.5)
        
        #Generate replicate data and compute fit stats for them
        y.new[i,j,k] ~ dbin(p[i,j,k], N[i,k])
        E.new[i,j,k] <- pow((y.new[i,j,k] - eval[i,j,k]),2) / (eval[i,j,k]+0.5)
      } 
      ik.p[i,k] <- mean(p[i,,k])
    } 
  } 
  
  #Derived and other quantities
  for (k in 1:7){
    totalN[k] <- sum(N[,k]) #Estimate total pop. size across all
    sites
    mean.abundance[k] <- mean(lambda[,k])
    mean.N[k] <- mean(N[,k])
    mean.detection[k] <- mean(ik.p[,k])
  }
  fit <- sum(E[,,])
  fit.new <- sum(E.new[,,])
}

#Bundle data
R = nrow(y)
T = ncol(y)
win.data <- list(y = y, R = R, T = T)

#Initial values
Nst <- apply(y, c(1, 3), max) + 1
Nst[is.na(Nst)] <- 1
inits <- function(){list(N = Nst, alpha.lam = runif(7, -3, 3), beta = runif(7, -3, 3), sd.lam = runif(1, 0, 1), sd.p = runif(1, 0, 1))}

#Parameters monitored
params <- c("totalN", "alpha.lam", "beta", "sd.lam", "sd.p","mean.abundance", "mean.N", "mean.detection", "fit", "fit.new")

#MCMC settings
ni <- 350000
nt <- 300
nb <- 50000
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

print(out2, digits = 3)

k<-mcmcplots::as.mcmc.rjags(out2)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#Evaluation of fit
plot(out2$sims.list$fit, out2$sims.list$fit.new, main = "", xlab =
     "Discrepancy actual data", ylab = "Discrepancy replicate data",
     frame.plot = FALSE, xlim = c(50, 200), ylim = c(50, 200))
abline(0, 1, lwd = 2, col = "black")
mean(out2$sims.list$fit.new > out2$sims.list$fit)

mean(out2$mean$fit) / mean(out2$mean$fit.new)

#Plot
max.day.count <- apply(y, c(1, 3), max, na.rm = TRUE)
max.day.count[max.day.count == "-Inf"] <- NA
mean.max.count <- apply(max.day.count, 2, mean, na.rm = TRUE)
mean.max.count
par(mfrow = c(2, 1))
plot(1:7, mean.max.count, xlab = "Day", ylab = "Mean daily abundance",
     las = 1, ylim = c(0, 16), type = "b", main = "", frame.plot = FALSE,
     pch = 16, lwd = 2)
lines(1:7, out2$summary[24:30,5], type = "b", pch = 16, col = "blue", lwd = 2)
segments(1:7, out2$summary[24:30,3], 1:7, out2$summary[24:30,7], col = "blue")
plot(1:7, out2$summary[38:44,1], xlab = "Day", ylab = "Detection
     probability ", las = 1, ylim = c(0, 1), type = "b", 
     col = "blue",pch = 16, frame.plot = FALSE, lwd = 2)
segments(1:7, out2$summary[38:44,3], 1:7, out2$summary[38:44,7], col = "blue")
#-------------------------------------------------------------------------------