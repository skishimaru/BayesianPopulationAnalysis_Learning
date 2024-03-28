# Bayesian Population Analysis using WinBUGS
# Chapter 3: Introduction to the Generalized Linear Model: The Simplest Model for Count Data

library(R2jags) #to run JAGS
library(shinystan) #to run shiny stan
library(tidyverse) #to utilize pipe operators

# 3.2 Statistical Models: Response = Signal + Noise
# 3.2.1 The Noise Component
#-------------------------------------------------------------------------------
#See how a beta distribution with parameters 2 and 4 looks like
#Both lines generate a random sample of n from the specified beta distribution and produce a plot
plot(density(rbeta(n= 10^6, shape1= 2, shape2= 4)))
hist(rbeta(10^6, 2, 4), nclass= 100, col= "gray")
#-------------------------------------------------------------------------------

# 3.2.2 The Signal Component
#-------------------------------------------------------------------------------
#Example design matrix of a model using ANCOVA linear model
y <- c(25, 14, 68, 79, 64, 139, 49, 119, 111)
A <- factor(c(1, 1, 1, 2, 2, 2, 3, 3, 3))
X <- c(1, 14, 22, 2, 9, 20, 2, 13, 22)
plot(X,y,col= c(rep("red", 3), rep("blue", 3), rep("green", 3)),  xlim= c(-1, 25), ylim= c(0,140))

summary(fm <- lm(y ~ A-1+X)) #Running ANCOVA

#Further understanding the design matrix of the model
model.matrix(~A+X) #Effects or treatment contrast parameterization: Specifies in terms of baseline response
model.matrix(~A-1+X) #Means paramerization
model.matrix(~A+X) %*% fm$coefficients #Obtaining the linear predictor or the expected typical response
#-------------------------------------------------------------------------------

# 3.3 Poisson GLM in R and JAGS for Modeling Time Series Counts
# 3.3.1 Generation and Analysis of Simulated Data
#-------------------------------------------------------------------------------
#Generates Poisson counts of peregrine falcon population in the French Jura mountains over n years
#In this example, the linear predictor will be a cubic polynomial function of time
data.fn <- function(n= 40, alpha= 3.5576, beta1= -0.0912, beta2= 0.0091, beta3= -0.00014){
  #n= number of years
  #alpha, beta1, beta2, beta3= coefficients of a cubic polynomial of count on year
  year <- 1:n #Generates values of time covariate
  log.expected.count <- alpha+beta1*year+beta2*year^2+beta3*year^3 #Signal: build up systematic part of GLM
  expected.count <- exp(log.expected.count)
  C <- rpois(n=n, lambda= expected.count) #Noise: generate random part of the GLM, Poisson noise around expected counts
  
  #Simulated data
  plot(year, C, type= "b", lwd= 2, col= "black", main= "Figure 3.2A", las= 1, ylab= "Population Size", xlab= "Year", cex.lab= 1.2, cex.axis= 1.2)
  lines(year, expected.count, lwd= 3, col= "red")
  return(list(n=n, alpha=alpha, beta1=beta1, beta2=beta2, beta3=beta3, year=year, expected.count=expected.count, C=C))
}
data <- data.fn() #Output: population counts over 40 years and the trajectory overtime, realization of the stochastic process
#Black line= realized population size
#Red line= Expected population size

#DATA ANALYSIS 1: Frequentist Mode, using maximum likelihood
fm <- glm(C~ year + I(year^2) + I(year^3),family= poisson, data= data)
summary(fm)

#DATA ANALYSIS 2: Bayesian Mode, using JAGS
#Specify model in JAGS language
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #Priors
  alpha ~ dunif(-20, 20) 
  beta1 ~ dunif(-10, 10)
  beta2 ~ dunif(-10, 10)
  beta3 ~ dunif(-10, 10)
  
  #Likelihood: Note key components of a GLM on one line each
  for (i in 1:n){
    C[i] ~ dpois(lambda[i]) # 1. Distribution for random part 
    log(lambda[i]) <- alpha + beta1 * year[i] + beta2 * (year[i]*year[i]) + beta3 * (year[i]*year[i]*year[i]) # 3. Linear predictor, CHANGED FROM pow() to year*year because pow() does not work in JAGS
  }
  
} #model function

jags.data <- list(C= data$C, n= length(data$C), year= data$year) #Bundle data
inits <- function() list(alpha= runif(1,-2,2), beta1= runif(1,-3,3)) #initial values
params <- c("alpha", "beta1", "beta2", "beta3", "lambda") #parameters monitored

#MCMC settings
ni <- 2000 #number of iterations
nt <- 2 #thinning rate
nb <- 1000 #burn-in length
nc <- 3 #number of chains

#Call JAGS from R
out <- jags(data  = jags.data,
            inits = inits,
            parameters.to.save = params,
            model.file = jags.model.txt,
            n.chains = nc,
            n.iter = ni,
            n.burnin = nb)

print(out, dig = 3)

k<-mcmcplots::as.mcmc.rjags(out)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan


#DATA ANALYSIS 3:  
#Bundle data
mean.year <- mean(data$year) #mean of year covariate
sd.year <- sd(data$year) #SD of year covariate
jags.data <- list(C= data$C, n= length(data$C), year= (data$year - mean.year)/sd.year)

#Call JAGS from R 
out <- jags(data  = jags.data,
            inits = inits,
            parameters.to.save = params,
            model.file = jags.model.txt,
            n.chains = nc,
            n.iter = ni,
            n.burnin = nb)

#Summarize posteriors
print(out, dig=3)

k<-mcmcplots::as.mcmc.rjags(out)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#DATA ANALYSIS 4: Non-Convergence
#Repeat analysis with non convergence
#New MCMC settings with essentially no burn-in
ni <- 100 #number of draws per chain
nt <- 1 #thinning rate
nb <- 1 #burn-in length

#Call JAGS from R
tmp <- jags(data  = jags.data,
            inits = inits,
            parameters.to.save = params,
            model.file = jags.model.txt,
            n.chains = nc,
            n.iter = ni,
            n.burnin = nb)

print(tmp, dig=3)

k<-mcmcplots::as.mcmc.rjags(tmp)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#Plotting figure 3.2B
plot(1:40, data$C, lwd= 2, col= "black", main= "Figure 3.2B", las= 1, ylab= "Population Size", xlab= "Year")
R.predictions <- predict(glm(C ~ year + I(year^2) + I(year^3), family= poisson, data= data), type= "response") #R predicted values
lines(1:40, R.predictions, lwd= 3, col= "green") #graphing R predictions
JAGS.predictions <- out$BUGSoutput$mean$lambda #JAGS predicted values
lines(1:40, JAGS.predictions, lwd= 3, col= "blue", lty= 2) #graphing JAGS predictions
#this shows that the inferences from the R (frequentist mode) and JAGS (Bayesian mode) are numerically similar

cbind(R.predictions, JAGS.predictions) #comparing how similar they are
#-------------------------------------------------------------------------------

# 3.3.2 Analysis of a Real Data Set
#-------------------------------------------------------------------------------
#Analyzing the true peregrine population breeding in Jura from 1964 to 2003
peregrine <- read.table("falcons.txt", header= T) #WHERE DO WE GET THIS TXT?
attach(peregrine) #attached data set so we can directly write variable names 
plot(Year, Pairs, lwd= 2, main= "3.4A", las= 1, ylab= "Pair Count", xlab= "Year", ylim= c(0,200), pch= 16) #plotting the variable Pairs

#Fitting the model in JAGS
#Bundle data
mean.year <- mean(1:length(Year)) #mean of year covariate
sd.year <- sd(1:length(Year)) #SD of year covariate
jags.data <- list(C= Pairs, n= length(Pairs), year= (1:length(Year)- mean.year)/sd.year)
inits <- function()list(alpha= runif(1,-2,2), beta1= runif(1,-3,3)) #initial values
params <- c("alpha", "beta1", "beta2", "beta3", "lambda") #parameters monitored

#MCMC settings
ni= 2500
nt= 2
nb= 500
nc= 3

#Call JAGS from R
out1 <- jags.model(data= jags.data, inits= inits, parameters.to.save= params, model.file= "GLM_Poisson.txt", n.chains= nc, n.thin= nt,
             n.inter= ni, n.burnin= nb, debug= T, jags.directory= jags.dir, working.directory= getwd())

print(out1, dig=3) #summarize posteriors
#convergence looks good so now we plot the predicted population trajectory
JAGS.predictions <- out1$mean$lambda #Figure 3.4A
lines(Year, JAGS.predictions, lwd= 3, col= "blue", lty= 2)
#-------------------------------------------------------------------------------

# 3.4 Poisson GLM for Modeling Fecundity
#-------------------------------------------------------------------------------
#Modeling fecundity of the same peregrine population
plot(Year, Eyasses, lwd= 2, main= "3.4B", las= 1, ylab= "Nesting Count", xlab= "Year", ylim= c(0,260), pch= 16)

#Bundle data
mean.year= mean(1:length(Year)) #mean of year covariate
sd.year= sd(1:length(Year)) #SD of year covariate
jags.data <- list( C= Eyasses, n= length(Eyasses), year= 1:length(Year) - mean.year/sd.year)

#Call JAGS from R
out2 <- jags.model(data= jags.data, inits= inits, parameters.to.save= params, model.file= "GLM_Poisson.txt", n.chains= nc, n.thin= nt,
             n.oter= ni, n.burnin= nb, debug= T, jags.directory= jags.dir, working.directory= getwd())
#convergence looks good so now we plot figure 3.4B estimates
JAGS.predictions <- out2$mean$lambda
lines(Year, JAGS.predictions, lwd= 3, col= "blue")


#-------------------------------------------------------------------------------

# 3.5 Binomial GLM for Modeling Bounded Counts or Proportions
# 3.5.1 Generation and Analysis of Simulated Data
#-------------------------------------------------------------------------------
data.fn <- function(nyears= 40, alpha= 0, beta1= -0.1, beta2= -0.9) {
  #nyears= number of years 
  #coefficients= alpha, beta1, beta2
  
  #Generate un-transformed and transformed values of time covariate
  year <- 1:nyears
  YR <- (year-round(nyears/2))/(nyears/2)
  
  #Generate values of binomial totals (N)
  N <- round(runif(nyears, min= 20, max= 100))
  
  #Signal: build up systematic part of the GLM
  exp.p <- plogis(alpha + beta1*YR + beta2*(YR^2))
  
  #Noise: generate random part of the GLM, Binomial noise around expected counts (which is N)
  C <- rbinom(n= nyears, size= N, prob= exp.p)
  
  #Plot simulated data
  plot(year, C/N, lwd= 2, col= "black", main="", las= 1, ylab= "Proportion Successful Pairs", xlab= "Year", ylim= c(0,1))
  points(year, exp.p, lwd= 3, col= "red")
  return(list(nyears=nyears, alpha=alpha, beta1=beta1, beta2=beta2, year=year, YR=YR, exp.p=exp.p, C=C, N=N))
}

data <- data.fn(nyears= 40, alpha= 1, beta1= 0.03, beta2= -0.9) #data set inspired by section 3.5.2

#model was set up in textfile "3.5.1_JAGS_Code"
#Bundle data
jags.data <- list(C= data$C, N= data$N, nyears= length(data$C), year= data$YR)

inits <- function() list(alpha= runif(1,-1,1), beta1= runif(1,-1,1), beta2= runif(1,-1,1)) #initial values
params <- c("alpha", "beta1", "beta2", "p")

#MCMC settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

#Call JAGS from R
out <- jags.model(data=jags.data, inits=inits, parameters.to.save=params, model.file="GLM_Binomial.txt", n.chains= nc, n.thin= nt, n.inter= ni, n.burnin= nb, debug= T, jags.directory= jags.dir, working.directory= getwd())

#convergence looks good so now we plot the predicted proportion of successful pairs
JAGS.predictions <- out$mean$p
lines(1:length(data$C), JAGS.predictions, lwd= 3, col= "blue", lty= 2)
#-------------------------------------------------------------------------------

# 3.5.2 Analysis of a Real Data Set
#-------------------------------------------------------------------------------
#Read data and attach them
peregrine <- read.table("falcons.txt", header= T)
attach(peregrine)

#Bundle data (note another standardization for year)
jags.data <- list(C= R.Pairs, N= Pairs, nyears= length(Pairs), year= (Year-1985)/20)

#Initial values
inits <- function() list(alpha= runif(1,-1,1), beta1= runif(1,-1,1), beta2= runif(1,-1,1))

#Parameters monitored
params <- c("alpha", "beta1", "beta2", "p")

#MCMC settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

#Call JAGS from R
out3 <- jags.model(data= jags.data, inits= inits, parameters.to.save= params, modle.file= "GLM_Binomial.txt", n.chains= nc, n.thin= nt, n.inter= ni, n.burnin= nb, debug= T, jags.directory= jags.dir, working.directory= getwd())

#Summarize posteriors and plot estimates
print(out3, dig= 3)
plot(Year, R.Pairs/Pairs, lwd= 2, col= "black", main= "Figre 3.4C", las= 1, ylab= "Proportion Successful Pairs", xlab= "Year", ylim= c(0,1))
lines(Year, out3$mean$p, lwd= 3, col= "blue")
#-------------------------------------------------------------------------------
