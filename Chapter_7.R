# Bayesian Population Analysis using WinBUGS
# Chapter 7: Estimation of Survival from Capture-Recapture Data Using the Cormack-Jolly-Seber Model

library(R2jags) #to run JAGS
library(shinystan) #to run shiny stan
library(tidyverse) #to utilize pipe operators

# 7.3 Models with Constant Parameters
#-------------------------------------------------------------------------------
#Simple model where survival and recapture are identical for all individuals at all occasions
#Define parameter values
n.occasions <- 6 #Number of capture occasions
marked <- rep(50, n.occasions - 1) #Annual number of newly marked individuals 
phi <- rep(0.65, n.occasions - 1) #Survival of adults
p <- rep(0.4, n.occasions - 1) #Recapture probability

#Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol= n.occasions - 1, nrow= sum(marked))
P <- matrix(p, ncol= n.occasions - 1, nrow= sum(marked))

#Define function to simulate a capture history (CH) matrix
simul.cjs <- function(PHI, P, marked) {
  n.occasions <- dim(PHI)[2] + 1
  CH <- matrix(0, ncol= n.occasions, nrow= sum(marked))
  #Define a vector with the occasion of marking
  mark.occ <- rep(1:length(marked), marked[1:length(marked)])
  #Fill the CH matrix
  for (i in 1:sum(marked)) {
    CH[i, mark.occ[i]] <- 1 #Write 1 at the release occasion
    if (mark.occ[i] == n.occasions) next
    for (t in (mark.occ[i] + 1):n.occasions) {
      #Bernoulli trial: does individual survive an occasion?
      sur <- rbinom(1, 1, PHI[i, t-1])
      if (sur==0) break #if dead move to next individual
      rp <- rbinom(1, 1, P[i, t - 1])
      if (rp==1) CH[i, t] <- 1
    }
  }
return(CH)
}

#Execute function
CH <- simul.cjs(PHI, P, marked)

#Create vector with occasion of marking (f)
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

#Specify model in JAGS
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors and Constraints
  for (i in 1:nind) {
    for (t in f[i]: (n.occasions - 1)) {
      phi[i, t] <- mean.phi
      p[i, t] <-mean.p
    }
  }
  mean.phi ~ dunif(0,1) #Prior for mean survival
  mean.p ~ dunif(0,1) #prior for mean recapture
  
  #Likelihood
  for (i in 1:nind) {
    #Define latent state at first capture
    z[i,f[i]] <- 1
    for (t in (f[i] + 1):n.occasions) {
      #State process
      z[i, t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i, t-1] * z[i, t-1]
      #Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i, t-1] * z[i,t]
    }
  }
}

#Bundle data
jags.data <- list(y= CH, f= f, nind= dim(CH)[1], n.occasions= dim(CH)[2])

#Function to create a matrix of initial values for latent state z
ch.init <- function (ch, f) {
  for(i in 1:dim(ch)[1]) {ch[i, 1:f[i]] <- NA}
  return(ch)
}

#Initial values
known.state.cjs <- function(ch) { #NOT IN BOOK: In JAGS we must provide good initial values for the latent state z
  state <- ch
  for (i in 1:dim(ch)[1]) {
    n1 <- min(which(ch[i,]==1))
    n2 <- max(which(ch[i,]==1))
    state[i, n1:n2] <- 1
    state[i, n1] <- NA
  }
  state[state== 0] <- NA
  return(state)
} #When an individual was observed OR an individual was not observed but was alive for sure: z=1

inits <- function() {list(mean.phi= runif(1, 0, 1), mean.p= runif(1, 0, 1), z= known.state.cjs(CH))}

#Parameters monitored
params <- c("mean.phi", "mean.p")

#MCMC settings
ni <- 10000
nt <- 6
nb <- 5000
nc <- 3

#Call JAGS from R
cjs.c.c <- jags(data  = jags.data,
            inits = inits,
            parameters.to.save = params,
            model.file = jags.model.txt,
            n.chains = nc,
            n.thin= nt,
            n.iter = ni,
            n.burnin = nb)

print(cjs.c.c, dig = 3) #summarize Posteriors

k<-mcmcplots::as.mcmc.rjags(cjs.c.c)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan
#-------------------------------------------------------------------------------

# 7.3.1 Inclusion of Information about Latent State Variable
#-------------------------------------------------------------------------------
#To improve computation and convergence speed, we can add what we know about the latent state z (whenever we observe a marked individual we know its latent state is z=1)
#We also know that the latent state is z=1 for every occasion between the first and the last observation of an individual even if it was not seen at every occasion
#To include this, we need to create a matrix that has a value of of 1 at all occasions when the individual is alive and NAs everywhere else
#This model is conditional on the first capture, so the latent state is dependent on the first capture and we need NAs too

#Function to create a matrix with information about known latent state z
known.state.cjs <- function(ch) {
  state <- ch
  for (i in 1:dim(ch)[1]) {
    n1 <- min(which(ch[i,]==1))
    n2 <- max(which(ch[i,]==1))
    state[i, n1:n2] <- 1
    state[i, n1] <- NA
  }
  state[state== 0] <- NA
  return(state)
} #When an individual was observed OR an individual was not observed but was alive for sure: z=1

#Bundle data
jags.data <- list(y= CH, f= f, nind= dim(CH)[1], n.occasions= dim(CH)[2], z= known.state.cjs(CH))

#Function to create a matrix of initial values for latent state z
cjs.init.z <- function(ch,f){                    # initialize the state matrix.
  for (i in 1:dim(ch)[1]){
    if (sum(ch[i,])==1) next
    n2 <- max(which(ch[i,]==1)) 
    ch[i,f[i]:n2] <- NA
  }
  for (i in 1:dim(ch)[1]){ 
    ch[i,1:f[i]] <- NA
  }
  return(ch)
}

#Initial values
inits <- function() {list(z= cjs.init.z(CH,f), mean.phi= runif(1, 0, 1), mean.p= runif(1, 0, 1))}

#Parameters monitored
params <- c("mean.phi", "mean.p")

#MCMC settings
ni <- 10000
nt <- 6
nb <- 5000
nc <- 3

#Call JAGS from R
cjs.c.c <- jags(data  = jags.data,
                inits = inits,
                parameters.to.save = params,
                model.file = jags.model.txt,
                n.chains = nc,
                n.thin= nt,
                n.iter = ni,
                n.burnin = nb)

print(cjs.c.c, dig = 3) #summarize Posteriors

k<-mcmcplots::as.mcmc.rjags(cjs.c.c)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan
#The model now runs faster! its best to provide all latent state information 
#-------------------------------------------------------------------------------

# 7.4 Models with Time-Variation
# 7.4.1 Fixed Time Effects
#-------------------------------------------------------------------------------
#For this model survival and recapture vary independently over time (fixed-effects factor)

#Priors and Constraints to fit the fixed effect assumptions
for (i in 1:nind) {
  for (t in f[i]:(n.occasions-1)) {
    phi[i,t] <- alpha[t]
    p[i,t] <- beta[t]
  }
}
for (t in 1:n.occasions-1){
  alpha[t] ~ dunif(0,1) #Priors for time-specific survival
  beta[t] ~ dunif(0,1) #Priors for time-specific recapture
}
#-------------------------------------------------------------------------------

# 7.4.2 Random Time Effects 
#-------------------------------------------------------------------------------
#Define parameter values
n.occasions <- 20 #number of capture occasions
marked <- rep(30, n.occasions-1) #Annual variance of survival
mean.phi <- 0.65
var.phi <- 1 #Temporal variance of survival 
p <- rep(0.4, n.occasions-1)

#Determine annual survival probabilities
logit.phi <- rnorm(n.occasions-1, qlogis(mean.phi), var.phi^0.5)
phi <- plogis(logit.phi)

#Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol= n.occasions-1, nrow= sum(marked), byrow= TRUE)
P <- matrix(p, ncol= n.occasions-1, nrow= sum(marked))

#Simulate capture-histories
CH <- simul.cjs(PHI, P, marked)

#Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors and Constraints
  for (i in 1:nind) {
    for (t in f[i]:(n.occasions-1)) {
      logit(phi[i,t]) <- mu + epsilon[t]
      p[i,t] <- mean.p 
    }
  }
  for (t in 1:(n.occasions-1)) {
    epsilon[t] ~ dnorm(0, tau)
  }
  #mu ~ dnorm(0, 0.001) #Prior for logit of mean survival
  #mean.phi <- 1/(1+exp(-mu)) #Logit transformation
  mean.phi ~ dunif(0,1) #Prior for mean survival
  mu <- log(mean.phi/(1-mean.phi)) #Logit transformation
  sigma ~ dunif(0, 10) #Prior for standard deviation
  tau <- 1/(sigma*sigma)
  sigma2 <- sigma*sigma #Temporal variance
  mean.p ~ dunif(0,1)
  
  #Likelihood
  for (i in 1:nind) {
    #Define latent state at first capture
    z[i,f[i]] <- 1
    for (t in (f[i] + 1):n.occasions) {
      #State Process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- p[i, t-1]*z[i,t-1]
      #Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i, t-1]*z[i,t]
    }
  }
}

#Bundle data
jags.data <- list(y= CH, f= f, nind= dim(CH)[1], n.occasions= dim(CH)[2], z= known.state.cjs(CH))

#Initial values
inits <- function() {list(z= cjs.init.z(CH,f), mean.phi= runif(1, 0, 1), sigma= runif(1, 0, 10), mean.p= runif(1, 0, 1))}

#Parameters monitored
params <- c("mean.phi", "mean.p", "sigma2")

#MCMC settings
ni <- 10000
nt <- 6
nb <- 5000
nc <- 3

#Call JAGS from R
cjs.ran <- jags(data  = jags.data,
                inits = inits,
                parameters.to.save = params,
                model.file = jags.model.txt,
                n.chains = nc,
                n.thin= nt,
                n.iter = ni,
                n.burnin = nb)

print(cjs.ran, dig = 3) #summarize Posteriors

k<-mcmcplots::as.mcmc.rjags(cjs.ran)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#Create histogram
hist(cjs.ran$BUGSoutput$sims.list$sigma2, col= "gray", nclass= 35, las= 1, xlab= expression(sigma^2), main= "Figure 7.3")
abline(v= var.phi, col= "red", lwd= 2)
#-------------------------------------------------------------------------------

# 7.4.3 Temporal Covariates 
#-------------------------------------------------------------------------------
#Define parameter values
n.occasions <- 20 # Number of capture occasions
marked <- rep(15, n.occasions-1) # Annual number of newly marked individuals
mean.phi <- 0.65
p <- rep(0.4, n.occasions-1)
beta <- -0.3 # Slope of survival: winter relationship
r.var <- 0.2 # Residual temporal variance

#Draw annual survival probabilities
winter <- rnorm(n.occasions-1, 0, 1^0.5)
logit.phi <- qlogis(mean.phi) + beta*winter + rnorm(n.occasions-1, 0, r.var^0.5)
phi <- plogis(logit.phi)

#Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked), byrow = TRUE)
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))

#Simulate capture-histories
CH <- simul.cjs(PHI, P, marked)

#Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  #Priors and constraints
  for (i in 1:nind){
    for (t in f[i]:(n.occasions-1)){
      logit(phi[i,t]) <- mu + beta*x[t] + epsilon[t]
      p[i,t] <- mean.p
    } 
  } 
  for (t in 1:(n.occasions-1)){
    epsilon[t] ~ dnorm(0, tau)
    phi.est[t] <- 1 / (1+exp(-mu-beta*x[t]-epsilon[t])) #Yearly survival
  }
  
  mu ~ dnorm(0, 0.001) #Prior for logit of mean survival
  mean.phi <- 1 / (1+exp(-mu)) #Logit transformation
  beta ~ dnorm(0, 0.001); T(-10, 10) #Prior for slope parameter
  sigma ~ dunif(0, 10) #Prior on standard deviation
  tau <- 1/(sigma*sigma)
  sigma2 <- sigma*sigma#Residual temporal variance
  mean.p ~ dunif(0, 1) #Prior for mean recapture
  
  #Likelihood
  for (i in 1:nind){
    #Define latent state at first capture
    z[i,f[i]] <- 1
    for (t in (f[i]+1):n.occasions){
      #State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      #Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
    } 
  } 
}

#Bundle data
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH), x = winter)

#Initial values
inits <- function(){list(z = cjs.init.z(CH, f), mu = rnorm(1), sigma = runif(1, 0, 5), beta = runif(1, -5, 5), mean.p = runif(1, 0, 1))}

#Parameters monitored
params <- c("mean.phi", "mean.p", "phi.est", "sigma2", "beta")

# MCMC settings
ni <- 20000
nt <- 6
nb <- 10000
nc <- 3

#Call JAGS from R
cjs.ran <- jags(data  = jags.data,
                inits = inits,
                parameters.to.save = params,
                model.file = jags.model.txt,
                n.chains = nc,
                n.thin= nt,
                n.iter = ni,
                n.burnin = nb)

print(cjs.ran, dig = 3) #summarize Posteriors

k<-mcmcplots::as.mcmc.rjags(cjs.ran)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#Create graph
par(mfrow = c(1, 2), las = 1)
hist(cjs.ran$BUGSoutput$sims.list$beta, nclass = 25, col = "gray", main = "",
     xlab = expression(beta), ylab = "Frequency")
abline(v = -0.3, col = "red", lwd = 2)
hist(cjs.ran$BUGSoutput$sims.list$sigma2, nclass = 50, col = "gray", main = "",
     xlab = expression(sigma^2), ylab = "Frequency", xlim=c(0, 3))
abline(v = 0.2, col = "red", lwd = 2)
mtext("Figure 7.4", side= 3, line= -1.5, outer= T) #adding main title to multiplot
#-------------------------------------------------------------------------------

# 7.5 Models with Individual Variation
# 7.5.1 Fixed Group Effects
#-------------------------------------------------------------------------------
#Define parameter values
n.occasions <- 12 #Number of capture occasions
marked <- rep(30, n.occasions-1) #Annual number of newly marked individuals
phi.f <- rep(0.65, n.occasions-1) #Survival of females
p.f <- rep(0.6, n.occasions-1) #Recapture probability of females
phi.m <- rep(0.8, n.occasions-1) #Survival of males
p.m <- rep(0.3, n.occasions-1) #Recapture prob of males

#Define matrices with survival and recapture probabilities
PHI.F <- matrix(phi.f, ncol= n.occasions-1, nrow= sum(marked))
P.F <- matrix(p.f, ncol= n.occasions-1, nrow= sum(marked))
PHI.M <- matrix(phi.m, ncol= n.occasions-1, nrow= sum(marked))
P.M <- matrix(p.m, ncol= n.occasions-1, nrow= sum(marked))

#Simulate capture-histories
CH.F <- simul.cjs(PHI.F, P.F, marked)
CH.M <- simul.cjs(PHI.M, P.M, marked)

#Merge capture-histories by row
CH <- rbind(CH.F, CH.M)

#Create group variable
group <- c(rep(1,dim(CH.F)[1]), rep(2, dim(CH.M)[1]))

#Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors and Constraints
  for (i in 1:nind) {
    for (t in f[i]:(n.occasions-1)) {
      phi[i,t] <- phi.g[group[i]]
      p[i,t] <- p.g[group[i]]
    }
  }
  for (u in 1:g) {
    phi.g[u] ~ dunif(0,1) #Priors for group-specific survival
    p.g[u] ~ dunif(0,1) #Priors for group-specific recapture
  }
  
  #Likelihood
  for (i in 1:nind) {
    #Define latent state at first capture
    z[i, f[i]] <- 1
    for (t in (f[i]+1):n.occasions) {
      #State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1]*z[i,t-1]
      #Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i, t-1]*z[i,t]
    }
  }
}

#Bundle data
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH), g= length(unique(group)), group= group)

#Initial values
inits <- function(){list(z = cjs.init.z(CH, f), phi.g= runif(length(unique(group)), 0, 1), p.g= runif(length(unique(group)), 0, 1))}

#Parameters monitored
params <- c("phi.g", "p.g")

# MCMC settings
ni <- 5000
nt <- 3
nb <- 2000
nc <- 3

#Call JAGS from R
cjs.group <- jags(data  = jags.data,
                inits = inits,
                parameters.to.save = params,
                model.file = jags.model.txt,
                n.chains = nc,
                n.thin= nt,
                n.iter = ni,
                n.burnin = nb)

print(cjs.group, dig = 3) #summarize Posteriors

k<-mcmcplots::as.mcmc.rjags(cjs.group)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan
#-------------------------------------------------------------------------------

# 7.5.2 Random Group Effects
#-------------------------------------------------------------------------------
#Priors and Constraints to fit the random group effects
for (i in 1:nind) {
  for (t in f[i]: (n.occasions-1)) {
    logit(phi[i,t]) <- beta[group[i]]
    p[i,t] <- mean.p
  }
}
for (u in 1:g) {
  beta[u] ~ dnorm(mean.beta, tau)
  phi.g[u] <- 1/(1+exp(-beta[u])) #Back-transformed group-specific survival
}
mean.beta ~ dnorm(0, 0.001) #Prior for logit of mean survival
mean.phi <- 1/(1+exp(-mean.beta)) #Back-transformed for mean survival
sigma ~ dunif(0, 10) #Prior for sd of logit of survival variability
tau <- sigma*sigma 
mean.p ~ dunif(0,1) #Prior for mean recapture
#-------------------------------------------------------------------------------

# 7.5.3 Individual Random Effects
#-------------------------------------------------------------------------------
#A case where each individual belongs to its own group and is identified through treatment as a fixed effect

#Define parameter values
n.occasions <- 20 #Number of capture occasions
marked <- rep(30, n.occasions-1) #Annual number of newly marked individuals
mean.phi <- 0.65
p <- rep(0.4, n.occasions-1)
v.ind <- 0.5

#Draw annual survival probabilities
logit.phi <- rnorm(sum(marked), qlogis(mean.phi), v.ind^0.5)
phi <- plogis(logit.phi)

#Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol= n.occasions-1, nrow= sum(marked), byrow= FALSE)
P <- matrix(p, ncol= n.occasions-1, nrow= sum(marked))

#Simulate capture-histories
CH <- simul.cjs(PHI, P, marked)

#Create vector with occasion of marking 
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

#Specify model in JAGS language 
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors and Constraints
  for (i in 1:nind) {
    for (t in f[i]:(n.occasions-1)) {
      logit(phi[i,t]) <- mu + epsilon[i]
      p[i,t] <- mean.p
    }
  }
  for (i in 1:nind) {
    epsilon[i] ~ dnorm(0, tau)
  }
  mean.phi ~ dunif(0,1) #Prior for mean survival 
  mu <- log(mean.phi/(1-mean.phi)) #Logit transformation
  sigma ~ dunif(0,5) #Prior for standard deviation
  tau <- 1/(sigma*sigma)
  sigma2 <- sigma*sigma
  mean.p ~ dunif(0,1) #Prior for mean recapture
  
  #Likelihood
  for (i in 1:nind) {
    #Define latent state at first capture
    z[i, f[i]] <- 1
    for (t in (f[i]+1):n.occasions) {
      #State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1]*z[i,t-1]
      #Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i, t-1]*z[i,t]
    }
  }
}

#Bundle data
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH))

#Initial values
inits <- function(){list(z = cjs.init.z(CH, f), mean.phi= runif(1, 0, 1), mean.p= runif(1, 0, 1), sigma= runif(1, 0, 2))}

#Parameters monitored
params <- c("mean.phi", "mean.p", "sigma")

# MCMC settings
ni <- 5000
nt <- 3
nb <- 2000
nc <- 3

#Call JAGS from R
cjs.ind <- jags(data  = jags.data,
                  inits = inits,
                  parameters.to.save = params,
                  model.file = jags.model.txt,
                  n.chains = nc,
                  n.thin= nt,
                  n.iter = ni,
                  n.burnin = nb)

print(cjs.ind, dig = 3) #summarize Posteriors

k<-mcmcplots::as.mcmc.rjags(cjs.ind)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#Create graph
par(mfrow= c(1,2), las= 1)
hist(cjs.ind$BUGSoutput$sims.list$mean.phi, nclass= 25, col= "gray", main= "", xlab= expression(bar(phi)), ylab= "Frequency")
abline(v= mean.phi, col= "red", lwd= 2)
hist(cjs.ind$BUGSoutput$sims.list$sigma, nclass= 15, col= "gray", main="", xlab= expression(sigma^2), ylab= "Frequency", xlim= c(0,3))
abline(v= v.ind, col= "red", lwd= 2)
mtext("Figure 7.5", side= 3, line= -1.5, outer= T) #adding main title to multiplot
#Posterior distributions of mean survival and of the individual variance in survival. Red= values used for the simulation

#Priors and Constraints to estimate survival as a function of an individal covariate x
for (i in 1:nind) {
  for (t in f[i]: (n.occasions-1)) {
    logit(phi[i,t]) <- mu + beta*x[i] + epsilon[i]
    p[i,t] <- mean.p
  }
}
for (i in 1:nind) {
  epsilon[i] ~ dnrom(0, tau)
}
mean.phi ~ dunif(0,1) #Prior for mean survival
mu <- log(mean.phi/(1-mean.phi)) #logit transformation
beta ~ dnorm(0,0.001) #Prior for covariate slope
sigma ~ dunif(0,5) #Prior for standard deviation
tau <- 1/(sigma*sigma)
sigma2 <- sigma*sigma
mean.p ~ dunif(0,1) #Prior for mean recapture
#-------------------------------------------------------------------------------

# 7.6 Models with Time and Group Effects
# 7.6.1 Fixed Group and Time Effects
#-------------------------------------------------------------------------------
#Define parameter values
n.occasions <- 12 #Number of capture occasions
marked <- rep(50, n.occasions-1) #Annual number of newly marked individuals
phi.f <- c(0.6, 0.5, 0.55, 0.6, 0.5, 0.4, 0.6, 0.5, 0.55, 0.6, 0.7)
p.f <- rep(0.6, n.occasions-1)
diff <- 0.5 #Difference between male and female survival on logit scale
phi.m <- plogis(qlogis(phi.f) + diff)
p.m <- rep(0.3, n.occasions-1)

#Define matrices with survival and recapture probabilities
PHI.F <- matrix(rep(phi.f, sum(marked)), ncol= n.occasions-1, nrow= sum(marked), byrow= TRUE)
P.F <- matrix(rep(p.f, sum(marked)), ncol= n.occasions-1, nrow= sum(marked), byrow= TRUE)
PHI.M <- matrix(rep(phi.m, sum(marked)), ncol= n.occasions-1, nrow= sum(marked), byrow= TRUE)
P.M <- matrix(rep(p.m, sum(marked)), ncol= n.occasions-1, nrow= sum(marked), byrow= TRUE)

#Simulate capture-histories
CH.F <- simul.cjs(PHI.F, P.F, marked)
CH.M <- simul.cjs(PHI.M, P.M, marked)

#Merge capture-histories
CH <- rbind(CH.F, CH.M)

#Create group variable
group <- c(rep(1, dim(CH.F)[1]), rep(2, dim(CH.M)[1]))

#Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

#Specify model in JAGS language 
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors and Constraints
  for (i in 1:nind) {
    for (t in f[i]:(n.occasions-1)) {
      logit(phi[i,t]) <- beta[group[i]] + gamma[t]
      p[i,t] <- p.g[group[i]]
    }
  }
  #for survival parameters
  for (t in 1:(n.occasions-1)) {
    gamma[t] ~ dnorm(0, 0.01); T(-10,10) #Priors for time effects
    phi.g1[t] <- 1/(1+exp(-gamma[t])) #Back transformed survival of males
    phi.g2[t] <- 1/(1+exp(-gamma[t]-beta[2])) #Back transformed survival of females
  }
  beta[1] <- 0 #Corner constraint
  beta[2] ~ dnorm(0, 0.01); T(-10,10) #Prior for difference in male and female survival
  
  #For recapture parameters
  for(u in 1:g) {
    p.g[u] ~ dunif(0,1) #Priors for group-spec. recapture
  }
  
  #Likelihood
  for (i in 1:nind) {
    #Define latent state at first capture
    z[i, f[i]] <- 1
    for (t in (f[i]+1):n.occasions) {
      #State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1]*z[i,t-1]
      #Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i, t-1]*z[i,t]
    }
  }
}

#Bundle data
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH), g= length(unique(group)), group= group)

#Initial values
inits <- function(){list(z = cjs.init.z(CH, f), gamma= rnorm(n.occasions-1), beta= c(NA, rnorm(1)), p.g= runif(length(unique(group)), 0, 1))}

#Parameters monitored
params <- c("phi.g1", "phi.g2", "p.g", "beta")

# MCMC settings
ni <- 5000
nt <- 3
nb <- 2000
nc <- 3

#Call JAGS from R
cjs.add <- jags(data  = jags.data,
                inits = inits,
                parameters.to.save = params,
                model.file = jags.model.txt,
                n.chains = nc,
                n.thin= nt,
                n.iter = ni,
                n.burnin = nb)

print(cjs.add, dig = 3) #summarize Posteriors

k<-mcmcplots::as.mcmc.rjags(cjs.ind)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#Create graph
lower.f <- upper.f <- lower.m <- upper.m <- numeric() #Figure of male and female survival
for (t in 1:(n.occasions-1)){
  lower.f[t] <- quantile(cjs.add$BUGSoutput$sims.list$phi.g1[,t], 0.025)
  upper.f[t] <- quantile(cjs.add$BUGSoutput$sims.list$phi.g1[,t], 0.975)
  lower.m[t] <- quantile(cjs.add$BUGSoutput$sims.list$phi.g2[,t], 0.025)
  upper.m[t] <- quantile(cjs.add$BUGSoutput$sims.list$phi.g2[,t], 0.975)
}
plot(x=(1:(n.occasions-1))-0.1, y = cjs.add$BUGSoutput$mean$phi.g1, type = "b", pch = 16, ylim = c(0.2, 1), main= "Figure 7.6", ylab = "Survival probability", xlab = "Year", bty = "n", cex = 1.5, axes = FALSE)
axis(1, at = 1:11, labels = rep(NA,11), tcl = -0.25)
axis(1, at = seq(2,10,2), labels = c("2","4","6","8","10"))
axis(2, at = seq(0.2, 1, 0.1), labels = c("0.2", NA, "0.4", NA, "0.6", NA, "0.8", NA, "1.0"), las = 1)
segments((1:(n.occasions-1))-0.1, lower.f, (1:(n.occasions-1))-0.1, upper.f)
points(x = (1:(n.occasions-1))+0.1, y = cjs.add$BUGSoutput$mean$phi.g2, type = "b", pch = 1, lty = 2, cex = 1.5)
segments((1:(n.occasions-1))+0.1, lower.m, (1:(n.occasions-1))+0.1, upper.m)

#Priors and constraints to convey the interaction between sex and time
for (i in 1:nind) {
  for (t in f[i]: (n.occasions-1)) {
    phi[i,t] <- eta.phi[group[i], t]
    p[i,t] <- p.g[group[i]]
  }
}
#for survival parameters
for (u in 1:g) {
  for (t in 1: (n.occasions-1)) {
    eta.phi[u,t] ~ dunif(0,1) #Prior for time and group-spec. survival
  }
}
#for recapture parameters
for (u in 1:g) {
  p.g[u] ~ dunif(0,1) #Priors for group-spec. recapture
}
#-------------------------------------------------------------------------------

# 7.6.2 Fixed Group and Random Time Effects
#-------------------------------------------------------------------------------
#Used to estimate temporal variability of survival or recapture in each group separately
#Priors and constraints to conduct fixed group and random time effects
for (i in 1: nind) {
  for (t in f[i]:(n.occasions-1)) {
    logit(phi[i,t]) <- eta.phi[group[i], t]
    p[i,t] <- p.g[group[i]]
  }
}
#for survival parameters
for (u in 1:g) {
  for (t in 1:(n.occasions-1)) {
    eta.phi[u,t] <- mu.phi[u] + epsilon[u,t]
    epsilon[u,t] ~ dnorm(0, tau[u])
  }
  mean.phi[u] ~ dunif(0,1) #Priors on mean group-spec. survival
  mu.phi[u] <- log(mean.phi[u]/(1-mean.phi[u]))
  sigma[u] ~ dunif(0,10) #priors on mean group-spec.
  tau[u] <- 1/(sigma*sigma)
  sigma2[u] <- sigma*sigma
}
#for recapture parameters
for (u in 1:g) {
  p.g[u] ~ dunif(0,1) #Priors for group-spec. recapture
}

#Specify model in JAGS language 
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors and Constraints
  for (i in 1:nind) {
    for (t in f[i]:(n.occasions-1)) {
      logit(phi[i,t]) <- eta.phi[t, group[i]]
      p[i,t] <- p.g[group[i]]
    }
  }
  #for survival parameters
  for (t in 1:(n.occasions-1)) {
   eta.phi[t, 1:g] ~ dmnorm(mu.phi[], Omega[,])
  }
  for (u in 1:g) {
    mean.phi[u] ~ dunif(0,1) #Priors on mean group-spec. survival
    mean.phi[u] <- log(mean.phi[u]/(1-mean.phi[u]))
  }
  Omgega[1:g, 1:g] ~ dwish(R[,], df) #Priors for variance-covariance matrix
  Sigma[1:g, 1:g] <- inverse(Omega[,])
  
  #For recapture parameters
  for(u in 1:g) {
    p.g[u] ~ dunif(0,1) #Priors for group-spec. recapture
  }
}

#Likelihood
for (i in 1:nind) {
  #Define latent state at first capture
  z[i, f[i]] <- 1
  for (t in (f[i]+1):n.occasions) {
    #State process
    z[i,t] ~ dbern(mu1[i,t])
    mu1[i,t] <- phi[i,t-1]*z[i,t-1]
    #Observation process
    y[i,t] <- dbern(mu2[i,t])
    mu2[i,t] <- p[i, t-1]*z[i,t]
  }
}

#Bundle data
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH), g= length(unique(group)), group= group, R= matrix(c(1, 0, 0, 1), ncol= 2))

#Initial values
inits <- function(){list(z = cjs.init.z(CH, f), gamma= rnorm(n.occasions-1), beta= c(NA, rnorm(1)), p.g= runif(length(unique(group)), 0, 1))}

#Parameters monitored
params <- c("eta.phi", "p.g", "Sigma", "mean.phu")

# MCMC settings
ni <- 5000
nt <- 3
nb <- 2000
nc <- 3

#Call JAGS from R
cjs.corr <- jags(data  = jags.data,
                inits = inits,
                parameters.to.save = params,
                model.file = jags.model.txt,
                n.chains = nc,
                n.thin= nt,
                n.iter = ni,
                n.burnin = nb)

print(cjs.corr, dig = 3) #summarize Posteriors

k<-mcmcplots::as.mcmc.rjags(cjs.corr)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#To compute the temporal correlation of male and female survival
corr.coef <- cjs.corr$BUGSoutput$sims.list$Sigma[,1,2] / sqrt(cjs.corr$BUGSoutput$sims.list$Sigma[,1,1] * cjs.corr$BUGSoutput$sims.list$Sigma[,2,2])
#-------------------------------------------------------------------------------

# 7.7 Models with Age Effects
#-------------------------------------------------------------------------------
#Define parameter values
n.occasions <- 10 #Number of capture occasions
marked.j <- rep(200, n.occasions-1) #Annual number of newly marked juveniles
marked.a <- rep(30, n.occasions-1) #Annual number of newly marked adults
phi.juv <- 0.3 #Juvenile annual survival
phi.ad <- 0.65 # Adult annual survival
p <- rep(0.5, n.occasions-1) #Recapture
phi.j <- c(phi.juv, rep(phi.ad, n.occasions-2))
phi.a <- rep(phi.ad, n.occasions-1)

#Define matrices with survival and recapture probabilities
PHI.J <- matrix(0, ncol = n.occasions-1, nrow = sum(marked.j))
for (i in 1:length(marked.j)){
  PHI.J[(sum(marked.j[1:i])-marked.j[i]+1):sum(marked.j[1:i]),
        i:(n.occasions-1)] <- matrix(rep(phi.j[1:(n.occasions-i)],
                                         marked.j[i]), ncol = n.occasions-i, byrow = TRUE)
}
P.J <- matrix(rep(p, sum(marked.j)), ncol = n.occasions-1,
              nrow = sum(marked.j), byrow = TRUE)
PHI.A <- matrix(rep(phi.a, sum(marked.a)), ncol = n.occasions-1,
                nrow = sum(marked.a), byrow = TRUE)
P.A <- matrix(rep(p, sum(marked.a)), ncol = n.occasions-1,
              nrow = sum(marked.a), byrow = TRUE)
# Apply simulation function
CH.J <- simul.cjs(PHI.J, P.J, marked.j)
CH.A <- simul.cjs(PHI.A, P.A, marked.a)
# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f.j <- apply(CH.J, 1, get.first)
f.a <- apply(CH.A, 1, get.first)
# Create matrices X indicating age classes
x.j <- matrix(NA, ncol = dim(CH.J)[2]-1, nrow = dim(CH.J)[1])
x.a <- matrix(NA, ncol = dim(CH.A)[2]-1, nrow = dim(CH.A)[1])
for (i in 1:dim(CH.J)[1]){
  for (t in f.j[i]:(dim(CH.J)[2]-1)){
    x.j[i,t] <- 2
    x.j[i,f.j[i]] <- 1
  } #t
} #i
for (i in 1:dim(CH.A)[1]){
  for (t in f.a[i]:(dim(CH.A)[2]-1)){
    x.a[i,t] <- 2
  } 
} 

#Combined data sets
CH <- rbind(CH.J, CH.A)
f <- c(f.j, f.a)
x <- rbind(x.j, x.a)

#Specify model in JAGS language 
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors and Constraints
  for (i in 1:nind) {
    for (t in f[i]:(n.occasions-1)) {
      phi[i,t] <- beta[x[i,t]]
      p[i,t] <- mean.p
    }
  }
  for (u in 1:2) {
    beta[u] ~ dunif(0,1) #Priors for age specific survival
  }
  mean.p ~ dunif(0,1) #Prior for mean recapture
  
  #Likelihood
  for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- 1
    for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
    } 
  } 
}

#Bundle data
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions= dim(CH)[2], z = known.state.cjs(CH), x = x)

#Initial values
inits <- function(){list(z = cjs.init.z(CH, f), beta = runif(2, 0, 1), mean.p = runif(1, 0, 1))}

#Parameters monitored
params <- c("beta", "mean.p")

#MCMC settings
ni <- 2000
nt <- 3
nb <- 1000
nc <- 3

#Call JAGS from R
cjs.age <- jags(data  = jags.data,
                 inits = inits,
                 parameters.to.save = params,
                 model.file = jags.model.txt,
                 n.chains = nc,
                 n.thin= nt,
                 n.iter = ni,
                 n.burnin = nb)

print(cjs.age, dig = 3) #summarize Posteriors

k<-mcmcplots::as.mcmc.rjags(cjs.age)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#Create matrix X indicating age classes
x <- matrix(NA, ncol = dim(CH)[2]-1, nrow = dim(CH)[1])
for (i in 1:dim(CH)[1]){
  for (t in f[i]:(dim(CH)[2]-1)){
    x[i,t] <- t-f[i]+1
  } 
} 

#Priors and constraints for modeling survival as a linear function of age x
for (i in 1:nind){
  for (t in f[i]:(n.occasions-1)){
    logit(phi[i,t]) <- mu + beta*x[i,t]
    p[i,t] <- mean.p
  } #t
} #i
mu ~ dnorm(0, 0.01) #Prior for mean of logit survival
beta ~ dnorm(0, 0.01) #Prior for slope parameter
for (i in 1:(n.occasions-1)){
  phi.age[i] <- 1 / (1+exp(-mu-beta*i)) #Logit back-transformation
}
mean.p ~ dunif(0, 1) #Prior for mean recapture
#-------------------------------------------------------------------------------

# 7.8 Immediate Trap Response in Recapture Probability
#-------------------------------------------------------------------------------
#Import data
CH <- as.matrix(read.table(file= "/Users/shelbie.ishimaru/Documents/GitHub/BayesianPopulationAnalysis_Learning/trap.txt", sep= " "))

#Compute vector with occasion of first capture
get.first <- function(x) min(which(x!=0))
f<- apply(CH, 1, get.first)

#Create matrix m indication when an individual was captured 
m <- CH[,1:(dim(CH)[2]-1)]
u <- which(m==0)
m[u] <- 2
#in matrix: 1 denotes an individual was captured at the preceding occasion and 2 denotes that it was not captured at the preceding location

#Specify model in JAGS language 
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors and Constraints
  for (i in 1:nind) {
    for (t in f[i]:(n.occasions-1)) {
      phi[i,t] <- mean.phi
      p[i,t] <- beta[m[i,t]]
    }
  }
  mean.phi ~ dunif(0,1) #Prior for mean survival
  for (u in 1:2) {
    beta[u] ~ dunif(0,1) #Prior for recapture
  }
  
  #Likelihood
  for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- 1
    for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
    } 
  } 
}

#Bundle data
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions= dim(CH)[2], z = known.state.cjs(CH), m = m)

#Initial values
inits <- function(){list(z = cjs.init.z(CH, f), mean.phi = runif(1, 0, 1), beta = runif(2, 0, 1))}

#Parameters monitored
params <- c("mean.phi", "beta")

#MCMC settings
ni <- 20000
nt <- 3
nb <- 10000
nc <- 3

#Call JAGS from R
cjs.trap <- jags(data  = jags.data,
                inits = inits,
                parameters.to.save = params,
                model.file = jags.model.txt,
                n.chains = nc,
                n.thin= nt,
                n.iter = ni,
                n.burnin = nb)

print(cjs.trap, dig = 3) #summarize Posteriors

k<-mcmcplots::as.mcmc.rjags(cjs.trap)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan
#-------------------------------------------------------------------------------

# 7.9 Parameter Identifiable
#-------------------------------------------------------------------------------
# Define parameter values
n.occasions <- 12 #Number of capture occasions
marked <- rep(30, n.occasions-1) #Annual number of newly marked individuals
phi <- c(0.6, 0.5, 0.55, 0.6, 0.5, 0.4, 0.6, 0.5, 0.55, 0.6, 0.7)
p <- c(0.4, 0.65, 0.4, 0.45, 0.55, 0.68, 0.66, 0.28, 0.55, 0.45, 0.35)

#Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked), byrow = TRUE)
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked), byrow = TRUE)

#Simulate capture-histories
CH <- simul.cjs(PHI, P, marked)

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

#Specify model in JAGS language 
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors and Constraints
  for (i in 1:nind) {
    for (t in f[i]:(n.occasions-1)) {
      phi[i,t] <- phi.t[t]
      p[i,t] <- p.t[t]
    }
  }
  for (t in 1:(n.occasions-1)) {
    phi.t[t] ~ dunif(0,1) #Prior for the time specific survival
    p.t[t] ~ dunif(0,1) #prior for time specific recapture
  }
  
  #Likelihood
  for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- 1
    for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
    } 
  } 
}

#Bundle data
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions= dim(CH)[2], z = known.state.cjs(CH))

#Initial values
inits <- function(){list(z = cjs.init.z(CH, f), phi.t = runif((dim(CH)[2]-1),0,1), p.t = runif((dim(CH)[2]-1),0,1))}

#Parameters monitored
params <- c("phi.t", "p.t")

#MCMC settings
ni <- 25000
nt <- 3
nb <- 20000
nc <- 3

#Call JAGS from R
cjs.t.t <- jags(data  = jags.data,
                 inits = inits,
                 parameters.to.save = params,
                 model.file = jags.model.txt,
                 n.chains = nc,
                 n.thin= nt,
                 n.iter = ni,
                 n.burnin = nb)

print(cjs.t.t, dig = 3) #summarize Posteriors

k<-mcmcplots::as.mcmc.rjags(cjs.t.t)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#Plot posterior distributions of some phi and p
par(mfrow = c(2, 2), cex = 1.2, las = 1, mar=c(5, 4, 2, 1))
plot(density(cjs.t.t$BUGSoutput$sims.list$phi.t[,6]), xlim = c(0, 1), ylim = c(0, 5),
     main = "", xlab = expression(phi[6]), ylab = "Density", frame = FALSE,
     lwd = 2)
abline(h = 1, lty = 2, lwd = 2)
par(mar=c(5, 3, 2, 2))
plot(density(cjs.t.t$BUGSoutput$sims.list$phi.t[,11]), xlim = c(0, 1),
     ylim =c(0, 5), main = "", xlab = expression(phi[11]), ylab ="",
     frame = FALSE, lwd = 2)
abline(h = 1, lty = 2, lwd = 2)
par(mar=c(5, 4, 2, 1))
plot(density(cjs.t.t$BUGSoutput$sims.list$p.t[,6]), xlim = c(0, 1), ylim = c(0, 5),
     main = "", xlab = expression(p[6]), ylab = "Density", frame = FALSE,
     lwd = 2)
abline(h = 1, lty = 2, lwd = 2)
par(mar=c(5, 3, 2, 2))
plot(density(cjs.t.t$BUGSoutput$sims.list$p.t[,11]), xlim = c(0, 1), ylim =
       c(0, 5), main = "", xlab = expression(p[11]), ylab ="", frame = FALSE,
     lwd = 2)
abline(h = 1, lty = 2, lwd = 2)
mtext("Figure 7.9", side= 3, line= -1.5, outer= T) #adding main title to multiplot
#-------------------------------------------------------------------------------

# 7.10 Fitting the CJS to Data in the M-Array Format: The Multinomial Likelihood
# 7.10.1 Introduction
#-------------------------------------------------------------------------------
#Function to create a m-array based on capture-histories (CH)
marray <- function(CH){
  nind <- dim(CH)[1]
  n.occasions <- dim(CH)[2]
  m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
  
  #Calculate the number of released individuals at each time period
  for (t in 1:n.occasions){
    m.array[t,1] <- sum(CH[,t])
  }
  for (i in 1:nind){
    pos <- which(CH[i,]!=0)
    g <- length(pos)
    for (z in 1:(g-1)){
      m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
    } 
  }
  
  #Calculate the number of individuals that is never recaptured
  for (t in 1:n.occasions){
    m.array[t,n.occasions+1] <- m.array[t,1] -
      sum(m.array[t,2:n.occasions])
  }
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}
#-------------------------------------------------------------------------------

# 7.10.2 Time-Dependent Models
#-------------------------------------------------------------------------------
#Specify model in JAGS language 
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors and Constraints
  for (t in 1:(n.occasions-1)) {
    phi.t[t] ~ dunif(0,1) #Prior for the time specific survival
    p[t] ~ dunif(0,1) #prior for time specific recapture
  }
  
  #Define the multimodal likelihood
  for (t in 1:(n.occasions-1)) {
    marr[t,1:n.occasions] ~ dmulti(pr[t,], r[t])
  }
  
  #Calculate the number of birds released each year
  for (t in 1:(n.occasions-1)){
    r[t] <- sum(marr.a[t, ])
  }
  
  #Define the cell probabilities of the m-array
  #Main diagonal
  for (t in 1:(n.occasions-1)){
    q[t] <- 1-p[t] #Probability of non-recapture
    pr[t,t] <- phi[t]*p[t]
    #Above main diagonal
    for (j in (t+1):(n.occasions-1)){
      pr[t,j] <- prod(phi[t:j])*prod(q[t:(j-1)])*p[j]
    } 
    
    #Below main diagonal
    for (j in 1:(t-1)){
      pr[t,j] <- 0
    } 
  } 
  #Last column: probability of non-recapture
  for (t in 1:(n.occasions-1)){
    pr[t,n.occasions] <- 1-sum(pr[t,1:(n.occasions-1)])
  } 
  
  #Assess model fit using Freeman-Tukey statistic
  #Compute fit statistics for observed data
  for (t in 1:(n.occasions-1)){
    for (j in 1:n.occasions){
      expmarr[t,j] <- r[t]*pr[t,j]
      E.org[t,j] <- (sqrt(marr[t,j]) - sqrt(expmarr[t,j]))*(sqrt(marr[t,j]) - sqrt(expmarr[t,j]))
    } 
  } 
  
  #Generate replicate data and compute fit stats from them
  for (t in 1:(n.occasions-1)){
    marr.new[t,1:n.occasions] ~ dmulti(pr[t, ], r[t])
    for (j in 1:n.occasions){
      E.new[t,j] <- (sqrt(marr.new[t,j]) - sqrt(expmarr[t,j]))*(sqrt(marr.new[t,j]) - sqrt(expmarr[t,j]))
    } 
  } 
  fit <- sum(E.org[,])
  fit.new <- sum(E.new[,])
  }
  
#Create the m-array from the capture-histories
marr <- marray(CH)

#Bundle data
jags.data <- list(marr = marr, n.occasions = dim(marr)[2])

#Initial values
inits <- function(){list(phi = runif(dim(marr)[2]-1, 0, 1), p = runif(dim(marr)[2]-1, 0, 1))}

# Parameters monitored
params <- c("phi", "p", "fit", "fit.new")

# MCMC settings
ni <- 10000
nt <- 3
nb <- 5000
nc <- 3

#Call JAGS from R
cjs <- jags(data  = jags.data,
                inits = inits,
                parameters.to.save = params,
                model.file = jags.model.txt,
                n.chains = nc,
                n.thin= nt,
                n.iter = ni,
                n.burnin = nb)

print(cjs, dig = 3) #summarize Posteriors

k<-mcmcplots::as.mcmc.rjags(cjs)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#Evaluation of fit
plot(cjs$sims.list$fit, cjs$sims.list$fit.new, main= "Figure 7.9", xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", las = 1, ylim = c(5, 25), xlim = c(5, 25), bty ="n")
abline(0, 1, col = "black", lwd = 2)
mean(cjs$sims.list$fit.new > cjs$sims.list$fit)
#We found that the observed and simulated data are fairly similar, so the model is adequate for the data set. Bayesian p-value= 0.75
#-------------------------------------------------------------------------------

# 7.10.3 Age-Dependent Models
#-------------------------------------------------------------------------------
#Define parameter values
n.occasions <- 12 #Number of capture occasions
marked.j <- rep(200, n.occasions-1) #Annual number of newly marked juveniles
marked.a <- rep(30, n.occasions-1) #Annual number of newly marked adults
phi.juv <- 0.3 #Juvenile annual survival
phi.ad <- 0.65 #Adult annual survival
p <- rep(0.5, n.occasions-1) # Recapture
phi.j <- c(phi.juv, rep(phi.ad,n.occasions-2))
phi.a <- rep(phi.ad, n.occasions-1)

#Define matrices with survival and recapture probabilities
PHI.J <- matrix(0, ncol = n.occasions-1, nrow = sum(marked.j))
for (i in 1:(length(marked.j)-1)){
  PHI.J[(sum(marked.j[1:i])-marked.j[i]+1):sum(marked.j[1:i]),i:(n.occasions-1)] <- matrix(rep(phi.j[1:(n.occasions-i)], marked.j[i]), ncol = n.occasions-i, byrow = TRUE)
}
P.J <- matrix(rep(p, n.occasions*sum(marked.j)), ncol = n.occasions-1, nrow = sum(marked.j), byrow = TRUE)
PHI.A <- matrix(rep(phi.a, sum(marked.a)), ncol = n.occasions-1, nrow = sum(marked.a), byrow = TRUE)
P.A <- matrix(rep(p, sum(marked.a)), ncol = n.occasions-1, nrow = sum(marked.a), byrow = TRUE)

#Apply simulation function
CH.J <- simul.cjs(PHI.J, P.J, marked.j)
CH.A <- simul.cjs(PHI.A, P.A, marked.a)
cap <- apply(CH.J, 1, sum)
ind <- which(cap >= 2)
CH.J.R <- CH.J[ind,] #Juvenile CH recaptured at least once
CH.J.N <- CH.J[-ind,] #Juvenile CH never recaptured

#Remove first capture
first <- numeric()
for (i in 1:dim(CH.J.R)[1]){
  first[i] <- min(which(CH.J.R[i,]==1))
}
CH.J.R1 <- CH.J.R
for (i in 1:dim(CH.J.R)[1]){
  CH.J.R1[i,first[i]] <- 0
}

#Add grown-up juveniles to adults and create m-array
CH.A.m <- rbind(CH.A, CH.J.R1)
CH.A.marray <- marray(CH.A.m)
marr.a.r <- CH.A.marray #work around for DAG error, CHANGED BY SKI

#Create CH matrix for juveniles, ignoring subsequent recaptures
second <- numeric()
for (i in 1:dim(CH.J.R1)[1]){
  second[i] <- min(which(CH.J.R1[i,]==1))
}
CH.J.R2 <- matrix(0, nrow = dim(CH.J.R)[1], ncol = dim(CH.J.R)[2])
for (i in 1:dim(CH.J.R)[1]){
  CH.J.R2[i,first[i]] <- 1
  CH.J.R2[i,second[i]] <- 1
}

#Create m-array for these
CH.J.R.marray <- marray(CH.J.R2)

#The last column ought to show the number of juveniles not recaptured again and should all be zeros, since all of them are released as adults
CH.J.R.marray[,dim(CH.J)[2]] <- 0

#Create the m-array for juveniles never recaptured and add it to the previous m-array
CH.J.N.marray <- marray(CH.J.N)
CH.J.marray <- CH.J.R.marray + CH.J.N.marray
marr.j.r <- CH.J.marray #work around for DAG error, CHANGED BY SKI

#Specify model in JAGS language 
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors and constraints
  for (t in 1:(n.occasions-1)){
    phi.juv[t] <- mean.phijuv
    phi.ad[t] <- mean.phiad
    p[t] <- mean.p
  }
  mean.phijuv ~ dunif(0, 1) # Prior for mean juv. survival
  mean.phiad ~ dunif(0, 1) # Prior for mean ad. survival
  mean.p ~ dunif(0, 1) # Prior for mean recapture
  
  #Define the multinomial likelihood
  for (t in 1:(n.occasions-1)){
    marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,], r.j[t])
    marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,], r.a[t])
  }
  
  #Calculate the number of birds released each year
  for (t in 1:(n.occasions-1)){
    r.j[t] <- sum(marr.j.r[t,]) #CHANGED BY SKI
    r.a[t] <- sum(marr.a.r[t,]) #CHANGED BY SKI
  }
  
  #Define the cell probabilities of the m-arrays
  #Main diagonal
  for (t in 1:(n.occasions-1)){
    q[t] <- 1-p[t] #Probability of non-recapture
    pr.j[t,t] <- phi.juv[t]*p[t]
    pr.a[t,t] <- phi.ad[t]*p[t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
      pr.j[t,j] <- phi.juv[t]*prod(phi.ad[(t+1):j])*prod(q[t:(j-1)])*p[j]
      pr.a[t,j] <- prod(phi.ad[t:j])*prod(q[t:(j-1)])*p[j]
    } 
    # Below main diagonal
    for (j in 1:(t-1)){
      pr.j[t,j] <- 0
      pr.a[t,j] <- 0
    } 
  } 
  # Last column: probability of non-recapture
  for (t in 1:(n.occasions-1)){
    pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])
    pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
  } 
}

#Bundle data
jags.data <- list(marr.j = CH.J.marray, marr.j.r = marr.j.r, marr.a = CH.A.marray, marr.a.r = marr.a.r, n.occasions = dim(CH.J.marray)[2]) #CHANGED BY SKI

#Initial values
inits <- function(){list(mean.phijuv = runif(1, 0, 1), mean.phiad = runif(1, 0, 1), mean.p = runif(1, 0, 1))}

#Parameters monitored
params <- c("mean.phijuv", "mean.phiad", "mean.p")

#MCMC settings
ni <- 3000
nt <- 3
nb <- 1000
nc <- 3

cjs.2 <- jags(data  = jags.data,
            inits = inits,
            parameters.to.save = params,
            model.file = jags.model.txt,
            n.chains = nc,
            n.thin= nt,
            n.iter = ni,
            n.burnin = nb)

print(cjs.2, dig = 3) #summarize Posteriors

k<-mcmcplots::as.mcmc.rjags(cjs.2)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#Create graph
par(mfrow = c(1, 3), las = 1)
hist(cjs.2$BUGSoutput$sims.list$mean.phijuv, nclass = 30, col = "gray", main = "", xlab = "Juvenile survival", ylab = "Frequency")
abline(v = phi.juv, col = "red", lwd = 2)
hist(cjs.2$BUGSoutput$sims.list$mean.phiad, nclass = 30, col = "gray", main = "", xlab = "Adult survival", ylab = "")
abline(v = phi.ad, col = "red", lwd = 2)
hist(cjs.2$BUGSoutput$sims.list$mean.p, nclass = 30, col = "gray", main = "", xlab = "Recapture", ylab = "")
abline(v = p[1], col = "red", lwd = 2)
mtext("Figure 7.10", side= 3, line= -1.5, outer= T) #adding main title to multiplot
#-------------------------------------------------------------------------------

# 7.11 Analysis of a Real Data Set: Survival of Female Leisler's Bats
#-------------------------------------------------------------------------------
m.leisleri <- matrix(c(4,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
                       0,5,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
                       0,0,9,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
                       0,0,0,10,2,0,0,0,0,0,0,0,0,0,0,0,0,0,5,
                       0,0,0,0,10,2,1,0,0,0,0,0,0,0,0,0,0,0,6,
                       0,0,0,0,0,15,0,0,0,0,0,0,0,0,0,0,0,0,6,
                       0,0,0,0,0,0,11,2,0,1,0,0,0,0,0,0,0,0,19,
                       0,0,0,0,0,0,0,12,1,1,0,0,0,0,0,0,0,0,6,
                       0,0,0,0,0,0,0,0,13,2,0,0,0,0,0,0,0,0,4,
                       0,0,0,0,0,0,0,0,0,14,0,0,0,0,0,0,0,0,6,
                       0,0,0,0,0,0,0,0,0,0,13,1,0,0,0,1,0,0,8,
                       0,0,0,0,0,0,0,0,0,0,0,15,3,1,0,0,0,0,12,
                       0,0,0,0,0,0,0,0,0,0,0,0,12,4,0,1,0,0,7,
                       0,0,0,0,0,0,0,0,0,0,0,0,0,19,2,0,0,0,3,
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,28,1,0,0,4,
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,22,7,2,21,
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,12,2,21,
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,14,18), ncol = 19, nrow = 18, byrow = TRUE)

#Specify model in JAGS language 
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors and constraints
  for (t in 1:(n.occasions-1)){
    logit(phi[t]) <- mu + epsilon[t]
    epsilon[t] ~ dnorm(0, tau)
    p[t] <- mean.p
  }
  mean.phi ~ dunif(0, 1) # Prior for mean survival
  mu <- log(mean.phi / (1-mean.phi)) # Logit transformation
  sigma ~ dunif(0, 5) # Prior for standard deviation
  tau <- 1/(sigma*sigma)
  sigma2 <- sigma*sigma
  # Temporal variance on real scale
  sigma2.real <- sigma2 * (mean.phi*mean.phi) * ((1-mean.phi)*(1-mean.phi))
  mean.p ~ dunif(0, 1) # Prior for mean recapture
  
  #Define the multinomial likelihood
  for (t in 1:(n.occasions-1)) {
    marr[t,1:n.occasions] ~ dmulti(pr[t,], r[t])
  }
  
  #Calculate the number of birds released each year
  for (t in 1:(n.occasions-1)){
    r[t] <- sum(marr[t,])
  }
  
  #Define the cell probabilities of the m-array:
  #Main diagonal
  for (t in 1:(n.occasions-1)){
    q[t] <- 1-p[t]
    pr[t,t] <- phi[t]*p[t]
    #Above main diagonal
    for (j in (t+1):(n.occasions-1)){
      pr[t,j] <- prod(phi[t:j])*prod(q[t:(j-1)])*p[j]
    } 
    #Below main diagonal
    for (j in 1:(t-1)){
      pr[t,j]<-0
    } 
  } 
  #Last column: probability of non-recapture
  for (t in 1:(n.occasions-1)){
    pr[t,n.occasions] <- 1-sum(pr[t,1:(n.occasions-1)])
  }
  
  #Assess model fit using Freeman-Tukey statistic
  #Compute fit statistics for observed data
  for (t in 1:(n.occasions-1)){
    for (j in 1:n.occasions){
      expmarr[t,j] <- r[t]*pr[t,j]
      E.org[t,j] <- (sqrt(marr[t,j]) - sqrt(expmarr[t,j]))*(sqrt(marr[t,j]) - sqrt(expmarr[t,j]))
    }
  }
  
  #Generate replicate data and compute fit stats from them
  for (t in 1:(n.occasions-1)){
    marr.new[t,1:n.occasions] ~ dmulti(pr[t,], r[t])
    for (j in 1:n.occasions){
      E.new[t,j] <- (sqrt(marr.new[t,j]) - sqrt(expmarr[t,j]))*(sqrt(marr.new[t,j]) - sqrt(expmarr[t,j]))
    }
  }
  fit <- sum(E.org)
  fit.new <- sum(E.new)
}

#Bundle data
jags.data <- list(marr = m.leisleri, n.occasions = dim(m.leisleri)[2])

#Initial values
inits <- function(){list(mean.phi = runif(1, 0, 1), sigma = runif(1, 0, 5), mean.p = runif(1, 0, 1))}

# Parameters monitored
params <- c("phi", "mean.p", "mean.phi", "sigma2", "sigma2.real", "fit", "fit.new")

# MCMC settings
ni <- 5000
nt <- 3
nb <- 1000
nc <- 3

leis.result <- jags(data  = jags.data,
              inits = inits,
              parameters.to.save = params,
              model.file = jags.model.txt,
              n.chains = nc,
              n.thin= nt,
              n.iter = ni,
              n.burnin = nb)

print(leis.result, dig = 3) #summarize Posteriors

k<-mcmcplots::as.mcmc.rjags(leis.result)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#Produce figure of female survival probabilities
par(mfrow = c(1, 2), las = 1, mar=c(4, 4, 2, 2), mgp = c(3, 1, 0))
lower <- upper <- numeric()
T <- dim(m.leisleri)[2]-1
for (t in 1:T){
  lower[t] <- quantile(leis.result$sims.list$phi[,t], 0.025)
  upper[t] <- quantile(leis.result$sims.list$phi[,t], 0.975)
}
plot(y = leis.result$mean$phi, x = (1:T)+0.5, type = "b", pch = 16, ylim =
       c(0.3, 1), ylab = "Annual survival probability", xlab = "", axes = F)
axis(1, at = seq(1,(T+1),2), labels = seq(1990,2008,2))
axis(1, at = 1:(T+1), labels = rep("", T+1), tcl = -0.25)
axis(2, las = 1)
mtext("Year", 1, line = 2.25)
segments((1:T)+0.5, lower, (1:T)+0.5, upper)
segments(1, leis.result$mean$mean.phi, T+1, leis.result$mean$mean.phi,
         lty = 2, col = "red", lwd = 2)
segments(1, quantile(leis.result$sims.list$mean.phi,0.025), T+1,
         quantile(leis.result$sims.list$mean.phi, 0.025), lty = 2, col = "red")
segments(1, quantile(leis.result$sims.list$mean.phi, 0.975), T+1,
         quantile(leis.result$sims.list$mean.phi, 0.975), lty = 2, col = "red")
hist(leis.result$sims.list$sigma2.real, nclass = 45, col = "gray",
     main = "", las = 1, xlab = "")
mtext("Temporal variance of survival", 1, line = 2.25)
mtext("Figure 7.12", side= 3, line= -1.5, outer= T) #adding main title to multiplot

#Evaluation of fit
plot(leis.result$sims.list$fit, leis.result$sims.list$fit.new,
     main = "", xlab = "Discrepancy actual data", ylab = "Discrepancy
replicate data", las = 1, ylim = c(10, 35), xlim = c(10, 35), frame = FALSE)
abline(0, 1, col = "black")
#-------------------------------------------------------------------------------