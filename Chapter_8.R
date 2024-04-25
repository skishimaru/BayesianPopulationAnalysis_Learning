# Bayesian Population Analysis using WinBUGS
# Chapter 8: Estimation of Survival Using Mark-Recovery Data

library(R2jags) #to run JAGS
library(shinystan) #to run shiny stan
library(tidyverse) #to utilize pipe operators

# 8.2 The Mark-Recovery Model as a State-Space Model
# 8.2.1 Simulation of Mark-Recovery Data
#-------------------------------------------------------------------------------
#Also referred to as "dead-recovery" models because its data found tagged dead animals

#The code below simulates a mark-recovery matrix
#Define parameter values
n.occasions <- 14 #Number of release occasions
marked <- rep(50, n.occasions) #Annual number of marked individuals 
s <- rep(0.8, n.occasions)
r <- rep(0.2, n.occasions)

#Define matrices with survival and recoveru probabilities
S <- matrix(s, ncol= n.occasions, nrow= sum(marked))
R <- matrix(r, ncol= n.occasions, nrow= sum(marked))

#Define function to simulate mark-recovery data
simul.mr <- function(S, R, marked){
  n.occasions <- dim(S)[2]
  MR <- matrix(NA, ncol = n.occasions+1, nrow = sum(marked))
  #Define a vector with the occasion of marking
  mark.occ <- rep(1:n.occasions, marked)
  #Fill the CH matrix
  for (i in 1:sum(marked)){
    MR[i, mark.occ[i]] <- 1 #Write an 1 at the release occasion
    for (t in mark.occ[i]:n.occasions){
      #Bernoulli trial: has individual survived occasion?
      sur <- rbinom(1, 1, S[i,t])
      if (sur==1) next #If still alive, move to next occasion
      
      #Bernoulli trial: has dead individual been recovered?
      rp <- rbinom(1, 1, R[i,t])
      if (rp==0){
        MR[i,t+1] <- 0
        break
      }
      if (rp==1){
        MR[i,t+1] <- 1
        break
      }
    } 
  }
  
  #Replace the NA in the file by 0
  MR[which(is.na(MR))] <- 0
  return(MR)
}

#Execute function
MR <- simul.mr(S, R, marked)
#-------------------------------------------------------------------------------

# 8.2.2 Analysis of a Model with Constant Parameters
#-------------------------------------------------------------------------------
#To fit a dead-recovery model we need a vector indicating the occasion of marking each individual

#Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(MR, 1, get.first)

#Specify model in JAGS
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  #Priors and constraints
  for (i in 1:nind){
    for (t in f[i]:(n.occasions-1)){
      s[i,t] <- mean.s
      r[i,t] <- mean.r
    } #t
  } #i
  mean.s ~ dunif(0, 1) # Prior for mean survival
  mean.r ~ dunif(0, 1) # Prior for mean recapture
  # Likelihood
  for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- 1
    for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- s[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- r[i,t-1] * (z[i,t-1] - z[i,t])
    } 
  } 
}

#Define function to create a matrix with information about known latent state z
known.state.mr <- function(mr){
  state <- matrix(NA, nrow = dim(mr)[1], ncol = dim(mr)[2])
  rec <- which(rowSums(mr)==2)
  for (i in 1:length(rec)){
    n1 <- min(which(mr[rec[i],]==1))
    n2 <- max(which(mr[rec[i],]==1))
    state[rec[i],n1:n2] <- 1
    state[rec[i],n1] <- NA
    state[rec[i],n2:dim(mr)[2]] <- 0
  }
  return(state)
}

#Bundle data
jags.data <- list(y = MR, f = f, nind = dim(MR)[1], n.occasions = dim(MR)[2], z = known.state.mr(MR))

#Define function to create a matrix of initial values for latent state z
mr.init.z <- function(mr){
  ch <- matrix(NA, nrow = dim(mr)[1], ncol = dim(mr)[2])
  rec <- which(rowSums(mr)==1)
  for (i in 1:length(rec)){
    n1 <- which(mr[rec[i],]==1)
    ch[rec[i],n1:dim(mr)[2]] <- 0
    ch[rec[i],n1] <- NA
  }
  return(ch)
}

#Initial values
inits <- function(){list(z = mr.init.z(MR), mean.s = runif(1, 0, 1), mean.r = runif(1, 0, 1))}

#Parameters monitored
params <- c("mean.s", "mean.r")

#MCMC settings
ni <- 5000
nt <- 6
nb <- 2000
nc <- 3

#Call JAGS from R
mr.ss <- jags(data  = jags.data,
                  inits = inits,
                  parameters.to.save = params,
                  model.file = jags.model.txt,
                  n.chains = nc,
                  n.thin= nt,
                  n.iter = ni,
                  n.burnin = nb)

print(mr.ss, dig = 3) #summarize Posteriors
#Posterior means obtained from a state space likelihood are very similar to posterior means from a multinomial likelihood. The main difference is that a state-space likelihood is more demanding computationally and longer chains are needed to reach convergence

k<-mcmcplots::as.mcmc.rjags(mr.ss)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#-------------------------------------------------------------------------------

# 8.3 The Mark-Recovery Model Fitted with the Multinomial Likelihood
# 8.3.1 Constant Parameters
#-------------------------------------------------------------------------------
#Code below creates a m-array for mark-recovery (M-R) data

#Define function to create an m-array based for mark-recovery (MR) data
marray.dead <- function(MR){
  nind <- dim(MR)[1]
  n.occasions <- dim(MR)[2]
  m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
  #Create vector with occasion of marking
  get.first <- function(x) min(which(x!=0))
  f <- apply(MR, 1, get.first)
  #Calculate the number of released individuals at each time period
  first <- as.numeric(table(f))
  for (t in 1:n.occasions){
    m.array[t,1] <- first[t]
  }
  #Fill m-array with recovered individuals
  rec.ind <- which(apply(MR, 1, sum)==2)
  rec <- numeric()
  for (i in 1:length(rec.ind)){
    d <- which(MR[rec.ind[i],(f[rec.ind[i]]+1):n.occasions]==1)
    rec[i] <- d + f[rec.ind[i]]
    next_node <- m.array[f[rec.ind[i]], rec[i]] + 1 #SKI CHANGED!
    m.array[f[rec.ind[i]],rec[i]] <- next_node #SKI CHANGED!
  }
  # Calculate the number of individuals that are never recovered
  for (t in 1:n.occasions){
    m.array[t,n.occasions+1] <- m.array[t,1]-sum(m.array[t,2:n.occasions])
  }
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}

marr <- marray.dead(MR) #produce m-array

#Specify model in JAGS
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  #Priors and constraints
  for (t in 1:n.occasions){
    s[t] <- mean.s
    r[t] <- mean.r
  }
  mean.s ~ dunif(0, 1) # Prior for mean survival
  mean.r ~ dunif(0, 1) # Prior for mean recovery
  #Define the multinomial likelihood
  for (t in 1:n.occasions){
    marr[t,1:(n.occasions+1)] ~ dmulti(pr[t,], rel[t])
  }
  # Calculate the number of birds released each year
  for (t in 1:n.occasions){
    rel[t] <- sum(marr[t,])
  }
  #Define the cell probabilities of the m-array
  #Main diagonal
  for (t in 1:n.occasions){
    pr[t,t] <- (1-s[t])*r[t]
    #Above main diagonal
    for (j in (t+1):n.occasions){
      pr[t,j] <- prod(s[t:(j-1)])*(1-s[j])*r[j]
    } 
    #Below main diagonal
    for (j in 1:(t-1)){
      pr[t,j] <- 0
    }
  } 
  #Last column: probability of non-recovery
  for (t in 1:n.occasions){
    pr[t,n.occasions+1] <- 1-sum(pr[t,1:n.occasions])
  } 
}

#Bundle data
jags.data <- list(marr = marr, n.occasions = dim(marr)[2]-1)

#Initial values
inits <- function(){list(mean.s = runif(1, 0, 1), mean.r = runif(1, 0, 1))}

#Parameters monitored
params <- c("mean.s", "mean.r")

#MCMC settings
ni <- 5000
nt <- 6
nb <- 2000
nc <- 3

#Call JAGS from R
mr <- jags(data  = jags.data,
              inits = inits,
              parameters.to.save = params,
              model.file = jags.model.txt,
              n.chains = nc,
              n.thin= nt,
              n.iter = ni,
              n.burnin = nb)

print(mr, dig = 3) #summarize Posteriors

k<-mcmcplots::as.mcmc.rjags(mr)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#Create graph
par(mfrow = c(1, 2), las = 1)
hist(mr$sims.list$mean.s, nclass = 25, col = "gray", main = "", ylab =
       "Frequency", xlab = "Survival probability")
abline(v = 0.8, col = "red", lwd = 2)
hist(mr$sims.list$mean.r, nclass = 25, col = "gray", main = "", ylab =
       "", xlab = "Recovery probability")
abline(v = 0.2, col = "red", lwd = 2)
mtext("Figure 8.3", side= 3, line= -1.5, outer= T) #adding main title to multiplot
#-------------------------------------------------------------------------------

# 8.3.2 Age-Dependent Parameters
#-------------------------------------------------------------------------------
#Survival and recovery probabilities often depend on age!

n.occasions <- 15 # Number of occasions
marked.j <- rep(200, n.occasions) # Annual number of newly marked young
marked.a <- rep(20, n.occasions) # Annual number of newly marked adults
sjuv <- 0.3 # Juvenile survival probability
sad <- 0.8 # Adult survival probability
rjuv <- 0.25 # Juvenile recovery probability
rad <- 0.15 # Adult recovery probability
sj <- c(sjuv, rep(sad, n.occasions-1))
rj <- c(rjuv, rep(rad, n.occasions-1))

#Define matrices with survival and recovery probabilities
SJ <- matrix(0, ncol = n.occasions, nrow = sum(marked.j))
for (i in 1:length(marked.j)){
  SJ[(sum(marked.j[1:i])-marked.j[i]+1):sum(marked.j[1:i]),
     i:n.occasions] <- matrix(rep(sj[1:(n.occasions-i+1)],
                                  marked.j[i]), ncol = n.occasions-i+1, byrow = TRUE)
}
SA <- matrix(sad, ncol = n.occasions, nrow = sum(marked.a))
RJ <- matrix(0, ncol = n.occasions, nrow = sum(marked.j))
for (i in 1:length(marked.j)){
  RJ[(sum(marked.j[1:i])-marked.j[i]+1):sum(marked.j[1:i]),
     i:n.occasions] <- matrix(rep(rj[1:(n.occasions-i+1)],
                                  marked.j[i]), ncol = n.occasions-i+1, byrow = TRUE)
}
RA <- matrix(rad, ncol = n.occasions, nrow = sum(marked.a))

#Execute simulation function
MRj <- simul.mr(SJ, RJ, marked.j)
MRa <- simul.mr(SA, RA, marked.a)

#Summarize data in m-arrays
marr.j <- marray.dead(MRj)
marr.a <- marray.dead(MRa)

#Specify model in JAGS
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  #Priors and constraints
  for (t in 1:n.occasions){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    rj[t] <- mean.rj
    ra[t] <- mean.ra
  }
  mean.sj ~ dunif(0, 1) # Prior for mean juv.survival
  mean.sa ~ dunif(0, 1) # Prior for mean ad.survival
  mean.rj ~ dunif(0, 1) # Prior for mean juv.recovery
  mean.ra ~ dunif(0, 1) # Prior for mean ad.recovery
  #Define the multinomial likelihoods
  for (t in 1:n.occasions){
    marr.j[t,1:(n.occasions+1)] ~ dmulti(pr.j[t,], rel.j[t])
    marr.a[t,1:(n.occasions+1)] ~ dmulti(pr.a[t,], rel.a[t])
  }
  #Calculate the number of birds released each year
  for (t in 1:n.occasions){
    rel.j[t] <- sum(marr.j[t,])
    rel.a[t] <- sum(marr.a[t,])
  }
  #Define the cell probabilities of the juvenile m-array
  #Main diagonal
  for (t in 1:n.occasions){
    pr.j[t,t] <- (1-sj[t])*rj[t]
    #Further above main diagonal
    for (j in (t+2):n.occasions){
      pr.j[t,j] <- sj[t]*prod(sa[(t+1):(j-1)])*(1-sa[j])*ra[j]
    } 
    #Below main diagonal
    for (j in 1:(t-1)){
      pr.j[t,j] <- 0
    } 
  } 
  for (t in 1:(n.occasions-1)){
    #One above main diagonal
    pr.j[t,t+1] <- sj[t]*(1-sa[t+1])*ra[t+1]
  } 
  #Last column: probability of non-recovery
  for (t in 1:n.occasions){
    pr.j[t,n.occasions+1] <- 1-sum(pr.j[t,1:n.occasions])
  } 
  #Define the cell probabilities of the adult m-array
  #Main diagonal
  for (t in 1:n.occasions){
    pr.a[t,t] <- (1-sa[t])*ra[t]
    #Above main diagonal
    for (j in (t+1):n.occasions){
      pr.a[t,j] <- prod(sa[t:(j-1)])*(1-sa[j])*ra[j]
    }
    #Below main diagonal
    for (j in 1:(t-1)){
      pr.a[t,j] <- 0
    } 
  }
  #Last column: probability of non-recovery
  for (t in 1:n.occasions){
    pr.a[t,n.occasions+1] <- 1-sum(pr.a[t,1:n.occasions])
  }
}

#Bundle data
jags.data <- list(marr.j = marr.j, marr.a = marr.a, n.occasions = dim(marr.j)[2]-1)

#Initial values
inits <- function(){list(mean.sj = runif(1, 0, 1), mean.sa = runif(1, 0, 1), mean.rj = runif(1, 0, 1), mean.ra = runif(1, 0, 1))}

#Parameters monitored
params <- c("mean.sj", "mean.rj", "mean.sa", "mean.ra")

#MCMC settings
ni <- 5000
nt <- 6
nb <- 2000
nc <- 3

#Call JAGS from R
mr.age <- jags(data  = jags.data,
           inits = inits,
           parameters.to.save = params,
           model.file = jags.model.txt,
           n.chains = nc,
           n.thin= nt,
           n.iter = ni,
           n.burnin = nb)

print(mr.age, dig = 3) #summarize Posteriors

k<-mcmcplots::as.mcmc.rjags(mr.age)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan
#-------------------------------------------------------------------------------

# 8.4 Real-Data Example: Age-Dependent Survival in Swiss Red Kites
#-------------------------------------------------------------------------------
marray.juv <- c(42, 18, 5, 7, 4, 3, 2, 1, 2, 2, 1, 0, 1, 3, 0, 0, 1, 1388) #Juvenile "m-array" column
marray.ad <- c(3, 1, 1, 3, 0, 2, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 137) #Adult "m-array" column

#Specify model in JAGS
jags.model.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  #Priors and constraints
  sjuv ~ dbeta(4.2, 2.8) # Informative prior for juv. survival: Analysis A
  #sjuv ~ dunif(0, 1) # Non-informative for juv. survival prior: Analysis B
  ssub ~ dunif(0, 1) # Prior for subad. survival
  sad ~ dunif(0, 1) # Prior for ad. survival
  rjuv ~ dunif(0, 1) # Prior for juv. recovery
  rad ~ dunif(0, 1) # Prior for ad. recovery
  # Define the multinomial likelihoods
  marr.j[1:(n.age+1)] ~ dmulti(pr.j[], rel.j)
  marr.a[1:(n.age+1)] ~ dmulti(pr.a[], rel.a)
  #Calculate the number of birds released each year
  rel.j <- sum(marr.j[])
  rel.a <- sum(marr.a[])
  #Define the cell probabilities of the juvenile m-array
  #First element
  pr.j[1] <- (1-sjuv)*rjuv
  #Second element
  pr.j[2] <- sjuv*(1-ssub)*rad
  #Third and further elements
  for (t in 3:n.age){
    pr.j[t] <- sjuv*ssub*pow(sad,(t-3))*(1-sad)*rad
    
  }
  #Probability of non-recovery
  pr.j[n.age+1] <- 1-sum(pr.j[1:n.age])
  #Define the cell probabilities of the adult m-array
  #All elements
  for (t in 1:n.age){
    pr.a[t] <- pow(sad,(t-1))*(1-sad)*rad
  }
  #Probability of non-recovery
  pr.a[n.age+1] <- 1-sum(pr.a[1:n.age])
}

#Bundle data
jags.data <- list(marr.j = marray.juv, marr.a = marray.ad, n.age = length(marray.juv)-1)

#Initial values
inits <- function(){list(sjuv = runif(1, 0, 1), ssub = runif(1, 0, 1), sad = runif(1, 0, 1), rjuv = runif(1, 0, 1), rad = runif(1, 0, 1))}

#Parameters monitored
params <- c("sjuv", "ssub", "sad", "rjuv", "rad")

#MCMC settings
ni <- 30000
nt <- 10
nb <- 10000
nc <- 3


#Call JAGS from R
rk.ageA <- jags(data  = jags.data,
           inits = inits,
           parameters.to.save = params,
           model.file = jags.model.txt,
           n.chains = nc,
           n.thin= nt,
           n.iter = ni,
           n.burnin = nb)

print(rk.ageA, dig = 3) #summarize Posteriors

k<-mcmcplots::as.mcmc.rjags(rk.ageA)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#Call JAGS from R
rk.ageB <- jags(data  = jags.data,
                inits = inits,
                parameters.to.save = params,
                model.file = jags.model.txt,
                n.chains = nc,
                n.thin= nt,
                n.iter = ni,
                n.burnin = nb)

print(rk.ageB, dig = 3) #summarize Posteriors

k<-mcmcplots::as.mcmc.rjags(rk.ageB)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

#Plot posterior distributions
par(mfrow = c(2, 3), las = 1)
hist(rk.ageB$sims.list$sjuv, breaks = 20, col = "gray", main = "", xlab = "Juvenile survival")
hist(rk.ageB$sims.list$ssub, breaks = 20, col = "gray", main = "", xlab = "Subadult survival")
hist(rk.ageB$sims.list$sad, breaks = 20, col = "gray", main = "", xlab = "Adult survival")
hist(rk.ageB$sims.list$rjuv, breaks = 20, col = "gray", main = "", xlab = "Juvenile recovery", xlim = c(0, 0.2))
hist(rk.ageB$sims.list$rad, breaks = 20, col = "gray", main = "", xlab = "Adult recovery")
mtext("Figure 8.5", side= 3, line= -1.5, outer= T) #adding main title to multiplot

#Plot comparison posterior distributions of juvenile survival under two priors
plot(density(rk.ageA$sims.list$sjuv), ylim = c(0, 5), lwd = 2, main = "Figure 8.6", xlab = "Juvenile survival", las = 1)
points(density(rk.ageB$sims.list$sjuv), col = "red", type = "l", lwd = 2)
text(x = 0.5, y = 4.8, "Prior distributions", pos = 4, font = 3)
legend(x = 0.6, y = 4.7, legend = c("U(0,1)", "beta(4.2,2.8)"), lwd = c(2, 2), col = c("black", "red"), bty = "n")

quantile(rk.ageA$sims.list$ssub-rk.ageA$sims.list$sjuv, prob = c(0.025, 0.975))
quantile(rk.ageA$sims.list$sad-rk.ageA$sims.list$ssub, prob = c(0.025, 0.975))
#-------------------------------------------------------------------------------