# Bayesian Population Analysis using WinBUGS
# Chapter 5: State-Space Models for Populations Counts

library(R2jags) #to run JAGS
library(shinystan) #to run shiny stan
library(tidyverse) #to utilize pipe operators

# 5.2 A Simple Model
#-------------------------------------------------------------------------------
#Population of Ibex has an initial size of 30 individuals and was monitored over 25 year. The population grew an average of 2% annually. Following snow melt, observers use a telescope to count the population. The survey design is not perfect and we assume the variance of the observation error is 20
n.years <- 25 #Number of years monitored
N1 <- 30 #Initial population size
mean.lambda <- 1.02 #Mean annual population growth rate
sigma2.lambda <- 0.02 #Process (temporal) variation of the growth rate
sigma2.y <- 20 #Variance of the observation error

#Simulating population size under exponential growth
y <- N <- numeric(n.years)
N[1] <- N1
lambda <- rnorm(n.years-1, mean.lambda, sqrt(sigma2.lambda))
for (t in 1:(n.years-1)){
  N[t+1] <- N[t]*lambda[t]
}

#Generating the observed data conditional on the true population sizes
for (t in 1:n.years){
  y[t] <- rnorm(1, N[t], sqrt(sigma2.y))
}

#Specify model in JAGS
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors and Contraints
  N.est[1] ~ dunif(0, 500) #Prior for initial population size
  mean.lambda ~ dunif(0, 10) #Prior for mean growth rate
  sigma.proc ~ dunif(0, 10) #Prior for sd of state process
  sigma2.proc <- (sigma.proc*sigma.proc) #equivalent to pow(sigma.proc, 2) from book
  tau.proc <- 1/(sigma.proc*sigma.proc) #equivalent to pow(sigma.proc, -2) from book
  sigma.obs ~ dunif(0, 100) #Prior for sd of observation process
  sigma2.obs <- (sigma.obs*sigma.obs) #equivalent to pow(sigma.obs, 2) from book
  tau.obs <- 1/(sigma.obs*sigma.obs) #equivalent to pow(sigma.obs, -2) from book
  
  #Likelihood/State Process
  for (t in 1:(T-1)){
    lambda[t] ~ dnorm(mean.lambda, tau.proc)
    N.est[t+1] <- N.est[t]*lambda[t]
  }
  
  #Observation process
  for (t in 1:T){
    y[t] ~dnorm(N.est[t], tau.obs)
  }
}

jags.data <- list(y= y, T= n.years) #Bundle data
inits <- function(){list(sigma.proc= runif(1, 0, 5), mean.lambda= runif(1, 0.1, 2), sigma.obs= runif(1, 0, 10), N.est= c(runif(1, 20, 40), rep(NA, (n.years-1))))} #initial values
params <- c("lambda", "mean.lambda", "sigma2.obs", "sigma2.proc", "N.est") #parameters monitored

#MCMC settings
ni <- 25000 #number of iterations
nt <- 3 #thinning rate
nb <- 10000 #burn-in length
nc <- 3 #number of chains, we run multiple chains to check for convergence 

#Call JAGS from R
ssm <- jags(data  = jags.data,
             inits = inits,
             parameters.to.save = params,
             model.file = jags.model.txt,
             n.chains = nc,
             n.thin= nt,
             n.iter = ni,
             n.burnin = nb)

print(ssm, dig = 3) #summarize Posteriors

k<-mcmcplots::as.mcmc.rjags(ssm)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#Creating Figure 5.3: comparing estimated, true, and observed population sizes
#Define function to draw a graph to summarize results
graph.ssm <- function(ssm, N, y, title){
  fitted <- lower <- upper <- numeric()
  n.years <- length(y)
  for (i in 1:n.years){
    fitted[i] <- mean(ssm$BUGSoutput$sims.list$N.est[,i])
    lower[i] <- quantile(ssm$BUGSoutput$sims.list$N.est[,i], 0.025)
    upper[i] <- quantile(ssm$BUGSoutput$sims.list$N.est[,i], 0.975)
  }
  m1 <- min(c(y, fitted, N, lower))
  m2 <- max(c(y, fitted, N, upper))
  par(mar= c(4.5, 4, 1, 1), cex= 1.2)
  plot(0, 0, ylim= c(m1, m2), xlim= c(0.5, n.years), ylab= "Population Size", xlab= "Year", main= title, las= 1, col= "black", type= "l",lwd= 2, frame= F, axes= F)
  axis(2, las= 1)
  axis(1, at= seq(0, n.years, 5), labels= seq(0, n.years, 5))
  axis(1, at= 0:n.years, label= rep("", n.years+1), tcl= -0.25)
  polygon(x= c(1:n.years, n.years:1), y= c(lower, upper[n.years:1]), col= "gray90", border= "gray90")
  points(N, type= "l", col= "red", lwd= 2) #plotting true population size
  points(y, type= "l", col= "black", lwd= 2) #plotting observed counts
  points(fitted, type= "l", col= "blue", lwd= 2) #plotting estimated population size
  legend("topleft", legend= c("True", "Observed", "Estimated"), lty= c(1, 1, 1), lwd= c(2, 2, 2), col= c("red", "black", "blue"), bty= "n", cex= 1)
}
graph.ssm(ssm, N, y, "Figure 5.3") #Execute function: Produce Figure
#-------------------------------------------------------------------------------

# 5.3 Systematic Bias in the Observation Process
#-------------------------------------------------------------------------------
#The ibex population remains stable at 50 individuals over 25 years. We assume the detection probability is constant at 0.7 so we count 70% of the population on average. Each individual can be seen or missed so the annual count is a binomial random variable written as Yt ~ Bin(Nt, p)
#Defining true population sizes
n.years <- 25 #Number of years monitored
N <- rep(50, n.years)

#Simulating the counts using the binomial distribution and inspecting the numbers
p <- 0.7
y <- numeric(n.years)
for(t in 1:n.years){
  y[t] <- rbinom(1, N[t], p)
}
y #observed counts are smaller than the true population size, and counts vary despite true population size and detection probability remaining constant

#Analyze the model
jags.data <- list(y= y, T= n.years) #Bundle data
inits <- function(){list(sigma.proc= runif(1, 0, 5), mean.lambda= runif(1, 0.1, 2), sigma.obs= runif(1, 0, 10), N.est= c(runif(1, 30, 60), rep(NA, (n.years-1))))} #initial values
params <- c("lambda", "mean.lambda", "sigma2.obs", "sigma2.proc", "N.est") #parameters monitored

#MCMC settings
ni <- 25000 #number of iterations
nt <- 3 #thinning rate
nb <- 10000 #burn-in length
nc <- 3 #number of chains, we run multiple chains to check for convergence 

#Call JAGS from R
ssm <- jags(data  = jags.data,
            inits = inits,
            parameters.to.save = params,
            model.file = jags.model.txt,
            n.chains = nc,
            n.thin= nt,
            n.iter = ni,
            n.burnin = nb)

print(ssm, dig = 3) #summarize Posteriors

k<-mcmcplots::as.mcmc.rjags(ssm)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#Comparing the estimated, true population sizes, and counts by creating figure 5.4
graph.ssm(ssm, N, y, "Figure 5.4")

#Simulating where detection probability is not constant overtime but is a positive trend that is linear on the logit scale 
#Defining true population sizes
n.years <- 25 #Number of years monitored
N <- rep(50, n.years)

#Simulating the counts using the binomial distribution
lp <- -0.5 + 0.1*(1:n.years) #Increasing trend of logit p
p <- plogis(lp)
y <- numeric(n.years)
for (t in 1:n.years){
  y[t] <- rbinom(1, N[t], p[t])
}

#Analyze the model
jags.data <- list(y= y, T= n.years) #Bundle data

#Call JAGS from R
ssm <- jags(data  = jags.data,
            inits = inits,
            parameters.to.save = params,
            model.file = jags.model.txt,
            n.chains = nc,
            n.thin= nt,
            n.iter = ni,
            n.burnin = nb)

print(ssm, dig = 3) #summarize Posteriors

k<-mcmcplots::as.mcmc.rjags(ssm)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#Produce Figure 5.5
graph.ssm(ssm, N, y, "Figure 5.5")
points(N*p, col= "black", type= "l", lwd= 2, lty= 2)
legend(x= -0.5, y= 55.5, legend= "Np", lwd= 2, col= "black", lty= 2, bty= "n")
#-------------------------------------------------------------------------------

# 5.4 Real Example: House Martin Population Counts in the Village of Magden
#-------------------------------------------------------------------------------
jags.model.txt <- function(){ 
  #Priors and Constraints
  logN.est[1] ~ dnorm(5.6, 0.01) #Prior for intial population size
  mean.r ~ dnorm(1, 0.001) #Prior for mean growth rate
  sigma.proc ~ dunif(0, 1) #Prior for sd of state process
  sigma2.proc <- (sigma.proc*sigma.proc) #equivalent to pow(sigma.proc, 2) from book
  tau.proc <- 1/(sigma.proc*sigma.proc) #equivalent to pow(sigma.proc, -2) from book
  sigma.obs ~ dunif(0, 1) #Prior for sd of observation process
  sigma2.obs <- (sigma.obs*sigma.obs) #equivalent to pow(sigma.obs, 2) from book
  tau.obs <- 1/(sigma.obs*sigma.obs) #equivalent to pow(sigma.obs, -2) from book
  
  #Likelihood/State Process
  for (t in 1:(T-1)){
    r[t] ~ dnorm(mean.r, tau.proc)
    logN.est[t+1] <- logN.est[t] + r[t]
  }
  
  #Observation Process
  for (t in 1:T){
    y[t] ~ dnorm(logN.est[t], tau.obs)
  }
  
  #Population sizes on real scale
  for (t in 1:T){
    N.est[t] <- exp(logN.est[t])
  }
}

#House Martin Population Data from Magden
pyears <- 6 #Number of future years with predictions (to get us to 2015)
hm <- c(271, 261, 309, 318, 231, 216, 208, 226, 195, 226, 233, 209, 226, 192, 191, 225, 245, 205, 191, 174, rep(NA, pyears)) #number of observed breeding pairs with NAs to support prediction into future years
year <- 1990:(2009+pyears)

#Analyze the model
jags.data <- list(y= log(hm), T= length(year)) #Bundle data
inits <- function(){list(sigma.proc= runif(1, 0, 1), mean.r= rnorm(1), sigma.obs= runif(1, 0, 1), logN.est= c(rnorm(1, 5.6, 0.1), rep(NA, length(year)-1)))} #initial values
params <- c("r", "mean.r", "sigma2.obs", "sigma2.proc", "N.est") #parameters monitored

#MCMC settings
ni <- 200000 #number of iterations
nt <- 6 #thinning rate
nb <- 100000 #burn-in length
nc <- 3 #number of chains, we run multiple chains to check for convergence 

#Call JAGS from R
hm.ssm <- jags(data  = jags.data,
            inits = inits,
            parameters.to.save = params,
            model.file = jags.model.txt,
            n.chains = nc,
            n.thin= nt,
            n.iter = ni,
            n.burnin = nb)

print(hm.ssm, dig = 2) #summarize Posteriors

k<-mcmcplots::as.mcmc.rjags(hm.ssm)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#Plotting Figure 5.7
fitted <- lower <- upper <- numeric()
year <- 1990:2015
n.years <- length(hm)
for (i in 1:n.years){
  fitted[i] <- mean(hm.ssm$BUGSoutput$sims.list$N.est[,i], 0.025)
  lower[i] <- quantile(hm.ssm$BUGSoutput$sims.list$N.est[,i], 0.975)
}
m1 <- min(c(fitted, hm, lower), na.rm= T)
m2 <- max(c(fitted, hm, upper), na.rm= T)
par(mar= c(4.5, 4, 1, 1))
plot(0, 0, ylim= c(m1, m2), xlim= c(1, n.years), ylab= "Population Size", xlab= "Year", main= "Figure 5.7", col= "black", type= "l", lwd= 2, axes= F, frame= F)
axis(2, las= 1)
axis(1, at= 1:n.years, labels= year)
polygon(x= c(1:n.years, n.years:1), y= c(lower, upper[n.years:1]), col= "gray90", border= "gray90")
points(hm, type= "l", col= "black", lwd= 2)
points(fitted, type= "l", col= "blue", lwd= 2)
legend("bottomleft", legend= c("Counts", "Estimates"), lty= c(1,1), lwd= c(2,2), col= c("black", "blue"), bty= "n", cex= 1)

#Probability of N(2015) < N(2009)
mean(hm.ssm$BUGSoutput$sims.list$N.est[,26] < hm.ssm$BUGSoutput$mean$N.est[20])
#-------------------------------------------------------------------------------