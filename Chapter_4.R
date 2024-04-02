# Bayesian Population Analysis using WinBUGS
# Chapter 4: Introduction to Random Effects: Conventional Poisson GLMM for Count Data

library(R2jags) #to run JAGS
library(shinystan) #to run shiny stan
library(tidyverse) #to utilize pipe operators
library(lme4) #for REML method (section 4.1)

# 4.1 Introduction
# 4.1.1 An Example
#-------------------------------------------------------------------------------
#We examined nine individuals asp vipers in three different populations to study the mass-length relationship
#Plotting figure 4.1
mass <- c(25, 14, 68, 79, 64, 139, 49, 119, 111)
pop <- factor(c(1, 1, 1, 2, 2, 2, 3, 3, 3))
length <- c(1, 14, 22, 2, 9, 20, 2, 13, 22)
plot(length, mass, col= c(rep("red", 3), rep("blue", 3), rep("green", 3)), xlim= c(-1, 25), ylim= c(0, 140), cex= 1.5, lwd= 2, frame.plot= F, las= 1, pch= 16, xlab= "Length", ylab= "Mass", main= "Figure 4.1")
#The plot suggests a linear relationship between mass and length with a different baseline in each population

#Fit fixed-effect model: print regression parameter estimates and plot regression lines
#Dotted Lines= Estimated regressions for each population under an ANCOVA model with fixed effects
summary(lm <- lm(mass~ pop-1 + length))
abline(lm$coef[1], lm$coef[4], col= "red", lwd= 3, lty= 2)
abline(lm$coef[2], lm$coef[4], col="blue", lwd= 3, lty= 2)
abline(lm$coef[3], lm$coef[4], col= "green", lwd= 3, lty= 2)

#Fit mixed model: print random effects and plot regression lines
# Solid lines= Estimated regressions under a mixed ANCOVA model with random effects (intercepts)
summary(lmm <- lmer(mass ~ length + (1|pop)))
abline((lmm@beta[1] + ranef(lmm)$pop)[1,], lmm@beta[2], col= "red", lwd= 3) ##CHANGED lmm@fixef to lmm@beta to correctly plot
abline((lmm@beta[1] + ranef(lmm)$pop)[2,], lmm@beta[2], col= "blue", lwd= 3) ##CHANGED lmm@fixef to lmm@beta to correctly plot
abline((lmm@beta[1] + ranef(lmm)$pop)[3,], lmm@beta[2], col= "green", lwd= 3) ##CHANGED lmm@fixef to lmm@beta to correctly plot
#-------------------------------------------------------------------------------

# 4.2 Accounting for Overdispersion by Random Effects: Modeling in R and JAGS
# 4.2.1 Generation and Analysis of Simulated Data
#-------------------------------------------------------------------------------
data.fn <- function(n= 40, alpha= 3.5576, beta1= -0.0912, beta2= 0.0091, beta3= -0.00014, sd= 0.1){
  #n: Number of years
  #alpha, beta1, beta2, beta3: coefficients of a cubic polynomial of count on year
  #sd: standard deviation of normal distribution assumed for year effects
  
  #Generate values of time covariate
  year <- 1:n
  
  #First level of noise: generate random year effects
  eps <- rnorm(n= n, mean= 0, sd= sd)
  
  #Signal (plus first level of noise): build up systematic part of the GLM and add the random year effects
  log.expected.count <- alpha + beta1*year + beta2*year^2 + beta3*year^3 + eps #Expected count now includes random variables so the plot will be stochastic and not smooth!
  expected.count <- exp(log.expected.count)
  
  #Second level of noise: generate random part of the GLM: Poisson noise around expected counts
  C <- rpois(n= n, lambda= expected.count)
  
  #Plot simulated data
  plot(year, C, lwd= 2, main="Figure 4.2", las= 1, ylab= "Population size", xlab= "Year", ylim= c(0, 1.1*max(C)))
  lines(year, expected.count, wd= 3, col= "red")
  return(list(n= n, alpha= alpha, beta1= beta1, beta2= beta2, beta3= beta3, year= year, sd= sd, expected.count= expected.count, C= C))
}

data <- data.fn() #With the added random variables the response now includes 2 random components: one from the year and the other a common poisson residual

#DATA ANALYSIS 1: Frequentist Mode, using the lme4 package
yr <- factor(data$year) #create a factor year
mny <- mean(data$year)
sdy <- sd(data$year)
cov1 <- (data$year - mny)/sdy
cov2 <- cov1*cov1
cov3 <- cov1*cov1*cov1
glmm.fit <- glmer(C ~ (1|yr) + cov1 + cov2 + cov3, family= poisson, data= data)
glmm.fit

R.predictions <- exp(fixef(glmm.fit)[1] + fixef(glmm.fit)[2]*cov1 + fixef(glmm.fit)[3]*cov2 +fixef(glmm.fit)[4]*cov3 + unlist(ranef(glmm.fit)))
lines(data$year, R.predictions, col= "green", lwd= 2)

#DATA ANALYSIS 2: Bayesian Mode using JAGS
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors
  alpha ~ dunif(-20, 20) 
  beta1 ~ dunif(-10, 10)
  beta2 ~ dunif(-10, 10)
  beta3 ~ dunif(-10, 10)
  tau <- 1/(sd*sd)
  sd ~ dunif(0,5)
  
  
  #Likelihood: Note key components of a GLM on one line each
  for (i in 1:n){
    C[i] ~ dpois(lambda[i]) # 1. Distribution for random part 
    log(lambda[i]) <- alpha + beta1 * year[i] + beta2 * (year[i]*year[i]) + beta3 * (year[i]*year[i]*year[i]) + eps[i] # 3. Linear predictor, CHANGED FROM pow() to year*year because pow() does not work in JAGS
    eps[i] ~ dnorm(0, tau) # 4. Definition of random effects dist
  }
}

jags.data <- list(C= data$C, n= length(data$C), year= cov1) #Bundle data
inits <- function() list(alpha= runif(1,-2,2), beta1= runif(1,-3,3), sd= runif(1, 0, 1)) #initial values
params <- c("alpha", "beta1", "beta2", "beta3", "lambda", "sd", "eps") #parameters monitored

#MCMC settings
ni <- 30000 #number of iterations
nt <- 10 #thinning rate
nb <- 20000 #burn-in length
nc <- 3 #number of chains, we run multiple chains to check for convergence 

#Call JAGS from R
out <- jags(data  = jags.data,
            inits = inits,
            parameters.to.save = params,
            model.file = jags.model.txt,
            n.chains = nc,
            n.thin= nt,
            n.iter = ni,
            n.burnin = nb)

print(out, dig = 3)

k<-mcmcplots::as.mcmc.rjags(out)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan
JAGS.predictions <- out$BUGSoutput$mean$lambda #JAGS predicted values
lines(data$year, JAGS.predictions, col= "blue", lwd= 2, lty= 2)

#Comparing standard error of the regression estimates under GLM and GLMM
glm.fit <- glm(C ~ cov1 + cov2 + cov3, family= poisson, data= data)
summary(glm.fit)
summary(glmm.fit)
#The standard error for all estimates is slightly higher for the GLMM, so the GLMM propagates the additional uncertainty in the modeled system into the regression estimates
#-------------------------------------------------------------------------------

# 4.2.2 Analysis of Real Data
#-------------------------------------------------------------------------------
#Analyzing the true peregrine population breeding in Jura from 1964 to 2003
peregrine <- read.table("/Users/shelbie.ishimaru/Documents/GitHub/BayesianPopulationAnalysis_Learning/falcons.txt", header= T) #read in falcon data

#DATA ANALYSIS 1: Frequentist Mode, using the lme4 package
#Standardize year to avoid convergence issues
yr <- factor(peregrine$Year)
mny <- mean(peregrine$Year)
sdy <- sd(peregrine$Year)
cov1 <- (peregrine$Year - mny)/sdy
cov2 <- cov1*cov1
cov3 <- cov1*cov1*cov1
glmm <- glmer(peregrine$Pairs ~ (1|yr) + cov1 + cov2 + cov3, family= poisson, data= peregrine)
glmm #Standard Deviation = 0.10

#DATA ANALYSIS 2: Bayesian Mode using JAGS
jags.data <- list(C= peregrine$Pairs, n= length(peregrine$Pairs), year= cov1)

#MCMC settings
ni <- 30000 #number of iterations
nt <- 10 #thinning rate
nb <- 20000 #burn-in length
nc <- 3 #number of chains, we run multiple chains to check for convergence 

#Call JAGS from R
out <- jags(data  = jags.data,
            inits = inits,
            parameters.to.save = params,
            model.file = jags.model.txt,
            n.chains = nc,
            n.thin= nt,
            n.iter = ni,
            n.burnin = nb)

print(out, dig = 3) #Standard Deviation= 0.12
k<-mcmcplots::as.mcmc.rjags(out)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan
#-------------------------------------------------------------------------------

# 4.3 Mixed Models with Random Effects for Variability in a Group (site and year effects)
# 4.3.1 Generation and Analysis of Simulated Data
#-------------------------------------------------------------------------------
data.fn <- function(nsite= 5, nyear= 40, alpha= 4.18456,beta1= 1.90672, beta2= 0.10852, beta3= -1.17121, sd.site= 0.5, sd.year= 0.2){
  #nsite= Number of populations
  #nyear= number of years
  #alpha, beta, beta2, beta3: cubic polynomial coefficients or year
  #sd.site: standard deviation of the normal distribution assumed for the population intercepts alpha
  #sd.year: standard deviation of the normal distribution assumed for the year effects
  #We standardize the year covariate so that it runs from about -1 to 1
  
  #Generate data structure to hold counts and log(lambda)
  C <- log.expected.count <- array(NA, dim= c(nyear, nsite))
  
  #Generate covariate values
  year <- 1:nyear
  yr <- (year-20)/20 #Standardize
  site <- 1:nsite
  
  #Draw two sets of random effects from their respective distribution
  alpha.site <- rnorm(n= nsite, mean= alpha, sd= sd.site)
  eps.year <- rnorm(n= nyear, mean= 0, sd= sd.year)
  
  #Loop over populations
  for (j in 1:nsite){
    #Signal (plus first level of noise): build up systematic part of the GLM including random site and year effects
    log.expected.count[,j] <- alpha.site[j] + beta1*yr + beta2*yr^2 + beta3*yr^3 + eps.year
    expected.count <- exp(log.expected.count[,j])
    
    #Second level of noise: generate random part of the GLM: Poission noise around expected counts
    C[,j] <- rpois(n= nyear, lambda= expected.count)
  }
  matplot(year, C, type= "l", lty= 1, lwd= 2, main= "Figure 4.3", las= 1, ylab= "Population size", xlab= "Year")
  return(list(nsite= nsite, nyear= nyear, alpha.site= alpha.site, beta1= beta1, beta2= beta2, beta3= beta3, year= year, sd.site= sd.site, sd.year= sd.year, expected.count= expected.count, C= C))
}

data <- data.fn(nsite= 100, nyear= 40, sd.site= 0.3, sd.year= 0.2)

#Bayesian Mode using JAGS
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors
  for (j in 1:nsite){
    alpha[j] ~ dnorm(mu, tau.alpha) # 4. Random site effects
  }
  mu ~ dnorm(0, 0.01) #Hyperparameter 1
  tau.alpha <- 1/(sd.alpha*sd.alpha) #Hyperparameter 2
  sd.alpha ~ dunif(0,2)
  for (p in 1:3){
    beta[p] ~ dnorm(0, 0.01)
  }
  tau.year <- 1/(sd.year*sd.year)
  sd.year ~ dunif(0, 1) #Hyperparameter 3
  
  #Likelihood: Note key components of a GLM on one line each
  for (i in 1:nyear){
    eps[i] ~ dnorm(0, tau.year) # 4.Random year effects
    for (j in 1:nsite){
      C[i,j] ~ dpois(lambda[i,j]) # 1.Distribution for random part
      lambda[i,j] <- exp(log.lambda[i,j]) #2. Link function
      log.lambda[i,j] <- alpha[j] + beta[1]*year[i] + beta[2]*(year[i]*year[i]) + beta[3]*(year[i]*year[i]*year[i]) + eps[i] # 3.Linear predictor including random site and random year effects
    }
  }
}

jags.data <- list(C= data$C, nsite= length(data$C), nyear= nrow(data$C), year= (data$year-20)/20) #Bundle data, standardized
inits <- function() list(mu= runif(1,0,2), alpha= runif(data$nsite, -1, 1), beta= runif(3,-1,1), sd.alpha= runif(1, 0, 0.1), sd.year= runif(1, 0, 0.1)) #initial values
params <- c("mu", "alpha", "beta", "sd.alpha", "sd.year", "lambda") #parameters monitored

#MCMC settings
ni <- 100000 #number of iterations
nt <- 50 #thinning rate
nb <- 50000 #burn-in length
nc <- 3 #number of chains, we run multiple chains to check for convergence 

#Call JAGS from R
out <- jags(data  = jags.data,
            inits = inits,
            parameters.to.save = params,
            model.file = jags.model.txt,
            n.chains = nc,
            n.thin= nt,
            n.iter = ni,
            n.burnin = nb)

print(out, dig = 3)

k<-mcmcplots::as.mcmc.rjags(out)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#-------------------------------------------------------------------------------

# 4.3.2 Analysis of a Real Data Set
#-------------------------------------------------------------------------------
tits <- read.table("/Users/shelbie.ishimaru/Documents/GitHub/BayesianPopulationAnalysis_Learning/tits.txt", header= T) #read in the data
str(tits)

C <- as.matrix(tits[5:13])
obs <- as.matrix(tits[14:22])
first <- as.matrix(tits[23:31])

matplot(1999:2007, t(C), lty= 1, lwd= 2, main= "Figure 4.5", las= 1, ylab= "Territory counts", xlab= "Year", ylim= c(0, 80), frame= F)

table(obs)
length(table(obs))
apply(first, 2, sum, na.rm= T)

a <- as.numeric(levels(factor(obs))) #all the levels, numeric
newobs <- ob #get ObsID from 1:271
for (j in 1:length(a)){newobs[which(obs==a[j])] <- j}

newobs[is.na(newobs)] <- 272
table(newobs)
first[is.na(first)] <- 0
table(first)

#MODEL 1: Null or Intercept only
#A constant expected count throughtout space and time
#Specify model in JAGS
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors
  alpha ~ dnorm(0, 0.01) #log(mean count)
  
  #Likelihood:
  for (i in 1:nyear){
    for (j in 1:nsite){
      C[i,j] ~ dpois(lambda[i,j]) 
      lambda[i,j] <- exp(log.lambda[i,j])
      log.lambda[i,j] <- alpha
    }
  }
}

jags.data <- list(C= t(C), nsite= nrow(C), nyear= ncol(C)) #Bundle data
inits <- function() list(alpha= runif(1, -10, 10)) #initial values
params <- c("alpha") #parameters monitored

#MCMC settings
ni <- 1200 #number of iterations
nt <- 2 #thinning rate
nb <- 200 #burn-in length
nc <- 3 #number of chains, we run multiple chains to check for convergence 

#Call JAGS from R
out0 <- jags(data  = jags.data,
            inits = inits,
            parameters.to.save = params,
            model.file = jags.model.txt,
            n.chains = nc,
            n.thin= nt,
            n.iter = ni,
            n.burnin = nb)

print(out0, dig = 3)

k<-mcmcplots::as.mcmc.rjags(out0)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan
#Good start but not a great model

#MODEL 2: Fixed Site Effects
#Essentially a one-way ANOVA, except that its a Poisson rather than a normal response
#Specify model in JAGS
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors
  for (j in 1:nsite){
    alpha[j] ~ dnorm(0,0.01) #site effects
  }
 
  #Likelihood:
  for (i in 1:nyear){
    for (j in 1:nsite){
      C[i,j] ~ dpois(lambda[i,j]) 
      lambda[i,j] <- exp(log.lambda[i,j]) 
      log.lambda[i,j] <- alpha[j]
    }
  }
}

jags.data <- list(C= t(C), nsite= nrow(C), nyear= nrow(C)) #Bundle data, standardized
inits <- function() list(alpha= runif(235, -1, 1)) #initial values
params <- c("alpha") #parameters monitored

#MCMC settings
ni <- 1200 #number of iterations
nt <- 2 #thinning rate
nb <- 200 #burn-in length
nc <- 3 #number of chains, we run multiple chains to check for convergence 

#Call JAGS from R
out1 <- jags(data  = jags.data,
            inits = inits,
            parameters.to.save = params,
            model.file = jags.model.txt,
            n.chains = nc,
            n.thin= nt,
            n.iter = ni,
            n.burnin = nb)

print(out1, dig = 2)

k<-mcmcplots::as.mcmc.rjags(out1)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan
#Adding the site effects greatly improved the fit of the model and deviance and DIC has gone down. So now lets add fixed year effects

#MODEL 3: Fixed Site and Fixed Year Effects
#A two-way main-effects ANOVA for a Poisson response
#Specify model in JAGS
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors
  for (j in 1:nsite){
    alpha[j] ~ dnorm(0, 0.01) #site effects
  }
  for (i in 2:nyear){ #nyear-1 year effects
    esp[i] ~ dnorm(0, 0.01)
  }
  eps[1] <- 0 #Aliased
  
  #Likelihood: Note key components of a GLM on one line each
  for (i in 1:nyear){
    for (j in 1:nsite){
      C[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- exp(log.lambda[i,j]) 
      log.lambda[i,j] <- alpha[j] + eps[i]
    }
  }
}

jags.data <- list(C= t(C), nsite= nrow(C), nyear= col(C)) #Bundle data
inits <- function() list(alpha= runif(235, -1, 1), eps= c(NA, runif(8, -1, 1))) #initial values
params <- c("alpha", "eps") #parameters monitored

#MCMC settings
ni <- 1200 #number of iterations
nt <- 2 #thinning rate
nb <- 200 #burn-in length
nc <- 3 #number of chains, we run multiple chains to check for convergence 

#Call JAGS from R
out2 <- jags(data  = jags.data,
            inits = inits,
            parameters.to.save = params,
            model.file = jags.model.txt,
            n.chains = nc,
            n.thin= nt,
            n.iter = ni,
            n.burnin = nb)

print(out2, dig = 2)

k<-mcmcplots::as.mcmc.rjags(out2)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan
#Again deviance and DIC has gone down and we can see annual population fluctuations. Now lets continue to build on the model by adding random instead of fixed effects

#MODEL 4: Random Site Effects (No Year Effects)
#The linear predictor of this model is that of a one-way, random-effects ANOVA
#Specify model in JAGS
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors
  for (j in 1:nsite){
    alpha[j] ~ dnorm(mu.alpha, tau.alpha) #Random site effects
  }
  mu ~ dnorm(0, 0.01) #Hyperparameter 1
  tau.alpha <- 1/(sd.alpha*sd.alpha) #Hyperparameter 2
  sd.alpha ~ dunif(0, 5)
  
  #Likelihood:
  for (i in 1:nyear){
    for (j in 1:nsite){
      C[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- exp(log.lambda[i,j])
      log.lambda[i,j] <- alpha[j]
      }
  }
}

jags.data <- list(C= t(C), nsite= nrow(C), nyear= nrow(C)) #Bundle data
inits <- function() list(mu= runif(1,2,3)) #initial values
params <- c("alpha", "mu.alpha", "sd.alpha") #parameters monitored

#MCMC settings
ni <- 1200 #number of iterations
nt <- 2 #thinning rate
nb <- 200 #burn-in length
nc <- 3 #number of chains, we run multiple chains to check for convergence 

#Call JAGS from R
out3 <- jags(data  = jags.data,
            inits = inits,
            parameters.to.save = params,
            model.file = jags.model.txt,
            n.chains = nc,
            n.thin= nt,
            n.iter = ni,
            n.burnin = nb)

print(out3, dig = 2)

k<-mcmcplots::as.mcmc.rjags(out3)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan
#The mass of the posterior distribution of the standard deviation of the random site effects is concentrated away from zero, which confirms our conclusion that sites differ substantially in expected counts
#Next let's add random year effects and re-parameterize the model to have a single grand mean to help convergence

#MODEL 5: Random Site and Random Year Effects
#A linear model that corresponds to a main-effect ANOVA with two random factors. Now we no longer need to constrain one effect of one factor to zero to avoid overparameterization. The borrowing strength among parameters within the same random-effects factor ensures that all can be estimated
#Specify model in JAGS
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors
  mu ~ dnorm(0, 0.01)
  for (j in 1:nsite){
    alpha[j] ~ dnorm(0, tau.alpha) #Random site effects
  }
  tau.alpha <- 1/(sd.alpha*sd.alpha)
  sd.alpha ~ dunif(0, 5)
  
  for (i in 1:nyear){
    eps[i] ~ dnorm(0, tau.eps) #random year effects
  }
  tau.eps <- 1/(sd.eps*sd.eps)
  sd.eps ~ dunif(0,3)
  
  #Likelihood:
  for (i in 1:nyear){
    for (j in 1:nsite){
      C[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- exp(log.lambda[i,j])
      log.lambda[i,j] <- mu + alpha[j] + eps[i]
    }
  }
}

jags.data <- list(C= t(C), nsite= nrow(C), nyear= nrow(C)) #Bundle data
inits <- function() list(mu= runif(1, 2, 3), alpha= runif(235, -2, 2), eps= runif(9, -1, 1)) #initial values
params <- c("mu", "alpha", "eps", "sd.alpha", "sd.eps") #parameters monitored

#MCMC settings
ni <- 6000 #number of iterations
nt <- 5 #thinning rate
nb <- 1000 #burn-in length
nc <- 3 #number of chains, we run multiple chains to check for convergence 

#Call JAGS from R
out4 <- jags(data  = jags.data,
             inits = inits,
             parameters.to.save = params,
             model.file = jags.model.txt,
             n.chains = nc,
             n.thin= nt,
             n.iter = ni,
             n.burnin = nb)

print(out4, dig = 2)

k<-mcmcplots::as.mcmc.rjags(out4)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan
#The year effects do not seem to be very different at any rate, they (sd.eps) are much smaller than the variation of counts among sites (sd.alpha). Next we add a fixed first-year observer effect

#MODEL 6: Random Site and Random ear Effects and First-Year Fixed Observer Effect
#Analogous to a three-way ANOVA with two factors random and one fixed and all of them acting in an additive way
#Specify model in JAGS
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors
  mu ~ dnorm(0, 0.01) #overall mean
  beta2 ~ dnorm(0, 0.01) #first-year observer effect
  
  for (j in 1:nsite){
    alpha[j] ~ dnorm(0, tau.alpha) #Random site effects
  }
  tau.alpha <- 1/(sd.alpha*sd.alpha)
  sd.alpha ~ dunif(0, 5)
  
  for (i in 1:nyear){
    eps[i] ~ dnorm(0, tau.eps) #random year effects
  }
  tau.eps <- 1/(sd.eps*sd.eps)
  sd.eps ~ dunif(0,5)
  
  #Likelihood:
  for (i in 1:nyear){
    for (j in 1:nsite){
      C[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- exp(log.lambda[i,j])
      log.lambda[i,j] <- mu + beta2*first[i,j] + alpha[j] + eps[i]
    }
  }
}

jags.data <- list(C= t(C), nsite= nrow(C), nyear= nrow(C), first= t(first)) #Bundle data
inits <- function() list(mu= runif(1, 0, 4), beta2= runif(1, -1, 1), alpha= runif(235, -2, 2), eps= runif(9, -1, 1)) #initial values
params <- c("mu", "beta", "alpha", "eps", "sd.alpha", "sd.eps") #parameters monitored

#MCMC settings
ni <- 6000 #number of iterations
nt <- 5 #thinning rate
nb <- 1000 #burn-in length
nc <- 3 #number of chains, we run multiple chains to check for convergence 

#Call JAGS from R
out5 <- jags(data  = jags.data,
             inits = inits,
             parameters.to.save = params,
             model.file = jags.model.txt,
             n.chains = nc,
             n.thin= nt,
             n.iter = ni,
             n.burnin = nb)

print(out5, dig = 2)

k<-mcmcplots::as.mcmc.rjags(out5)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan
#Beta2=0, so there is no evidence for a first year effect. Next we'll add an overall linear time trend

#Model 7: Random Site and Random Year Effects, First-Year Fixed Observer Effect, and Overall Linear Time Trend
#The linear part of the model corresponds to an analysis of covariance (ANCOVA. In addition to the three factors of the previous model, an added effect of one continuous covariate
#Specify the model in JAGS
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors
  mu ~ dnorm(0, 0.01) #overall intercept
  beta1 ~ dnorm(0, 0,01) #overall trend
  beta2 ~ dnrom(0, 0.01) #first-year observer effect
  
  for (j in 1:nsite){
    alpha[j] ~ dnorm(0, tau.alpha) #Random site effects
  }
  tau.alpha <- 1/(sd.alpha*sd.alpha)
  sd.alpha ~ dunif(0, 5)
  
  for (i in 1:nyear){
    eps[i] ~ dnorm(0, tau.eps) #random year effects
  }
  tau.eps <- 1/(sd.eps*sd.eps)
  sd.eps ~ dunif(0,3)
  
  #Likelihood:
  for (i in 1:nyear){
    for (j in 1:nsite){
      C[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- exp(log.lambda[i,j])
      log.lambda[i,j] <- mu + beta1*year[i] + beta2*first[i,j] + alpha[j] + eps[i]
    }
  }
}

jags.data <- list(C= t(C), nsite= nrow(C), nyear= nrow(C), first= t(first), year= ((1:9)-5)/4) #Bundle data
inits <- function() list(mu= runif(1, 0, 4), beta1= runif(1, -1, 1), beta2= runif(1, -1, 1), alpha= runif(235, -2, 2), eps= runif(9, -1, 1)) #initial values
params <- c("mu", "beta1", "beta2", "alpha", "eps", "sd.alpha", "sd.eps") #parameters monitored

#MCMC settings
ni <- 12000 #number of iterations
nt <- 6 #thinning rate
nb <- 6000 #burn-in length
nc <- 3 #number of chains, we run multiple chains to check for convergence 

#Call JAGS from R
out6 <- jags(data  = jags.data,
             inits = inits,
             parameters.to.save = params,
             model.file = jags.model.txt,
             n.chains = nc,
             n.thin= nt,
             n.iter = ni,
             n.burnin = nb)

print(out6, dig = 2)

k<-mcmcplots::as.mcmc.rjags(out6)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#MODEL 8: The Full Model
#Final model with random site effects, an overall linear time trend, random observer effects, random year effects, and a fixed first-year observer effect. We constrain the random effects standard deviations sufficiently so that JAGS doesn't get lost numerically and increase chain length
#Specify the model in JAGS
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors
  mu ~ dnorm(0, 0.01) #overall intercept
  beta1 ~ dnorm(0, 0,01) #overall trend
  beta2 ~ dnrom(0, 0.01) #first-year observer effect
  
  for (j in 1:nsite){
    alpha[j] ~ dnorm(0, tau.alpha) #Random site effects
  }
  tau.alpha <- 1/(sd.alpha*sd.alpha)
  sd.alpha ~ dunif(0, 3)
  
  for (i in 1:nyear){
    eps[i] ~ dnorm(0, tau.eps) #random year effects
  }
  tau.eps <- 1/(sd.eps*sd.eps)
  sd.eps ~ dunif(0,1)
  
  for (k in 1:nobs){
    gamma[k] ~ dnorm(0, tau.gamma) #random observer effect
  }
  tau.gamma <- 1/(sd.gamma*sd.gamma)
  sd.gamma ~dunif(0,1)
  
  #Likelihood:
  for (i in 1:nyear){
    for (j in 1:nsite){
      C[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- exp(log.lambda[i,j])
      log.lambda[i,j] <- mu + beta1*year[i] + beta2*first[i,j] + alpha[j] + gamma[newobs[i,j]] + eps[i]
    }
  }
}

jags.data <- list(C= t(C), nsite= nrow(C), nyear= nrow(C), nobs= 272, newobs= t(newobs), first= t(first), year= ((1:9)-5)/4) #Bundle data
inits <- function() list(mu= runif(1, 0, 4), beta1= runif(1, -1, 1), beta2= runif(1, -1, 1), alpha= runif(235, -1, 1), gamma= runif(272, -1, 1), eps= runif(9, -1, 1)) #initial values
params <- c("mu", "beta1", "beta2", "alpha", "gamma", "eps", "sd.alpha", "sd.eps") #parameters monitored

#MCMC settings
ni <- 12000 #number of iterations
nt <- 6 #thinning rate
nb <- 6000 #burn-in length
nc <- 3 #number of chains, we run multiple chains to check for convergence 

#Call JAGS from R
out7 <- jags(data  = jags.data,
             inits = inits,
             parameters.to.save = params,
             model.file = jags.model.txt,
             n.chains = nc,
             n.thin= nt,
             n.iter = ni,
             n.burnin = nb)

print(out7, dig = 2)

k<-mcmcplots::as.mcmc.rjags(out7)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan
#The most variation in coal tit counts are due to differences among site and differences among observers. There is also a small variance component due to years. 
#-------------------------------------------------------------------------------