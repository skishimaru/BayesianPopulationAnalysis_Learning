# Bayesian Population Analysis using WinBUGS
# Chapter 1: Introduction

# 1.3 Binomial Distribution: Sparrow Example
#-------------------------------------------------------------------------------
# An observation process that is dominated by non-detection error, where misclassification and double counts are absent
N <- 16 #population size of sparrows in the yard
p <- 0.4 #individual detection probability

rbinom(n= 1, size= N, prob= p) #simulating a single count, like field observation this changes every run

C <-rbinom(n= 10^6, size= N, prob= p) #simulating a one million counts to convey variation
mean(C)
var(C)
sd(C)
hist(C, breaks= 50, col= "Blue", main= "Figure 1.3", xlab= "Sparrow Count", ylab= "Density", las= 1, freq= F) #frequency distributions of one million counts

#The typical count C is smaller than the actual population size N. Indeed, the mean of a binomial random variable and hence the expected count of sparrows, equals the product of N and p.
#However, the counts vary quite a lot, even under totally identical conditions.Hence, there is nothing intrinsically wrong or inferior with smaller counts. Smaller counts may simply result from the random nature of the counting process in the presence of imperfect detection. Thus, any count with p< 1 will automatically tend to vary from trial to trial. Unless p= 0 or p= 1, it is impossible to eliminate that variation by the sampling design or standardization (though other components of variation may be eliminated).
#Actually, not only the mean count but also the magnitude of the variation of counts is known from statistical theory. The variance is equal to the product of N, p, and 1−p.
#-------------------------------------------------------------------------------

# 1.5 Bias and Precision: Asp Viper Example
#-------------------------------------------------------------------------------
# Mature Asp Vipers in Jura mountain average a length of 65 cm with an sd of 5
mu <- 65 #mean length of Jura mountain's mature asp viper population 
sigma <- 5 #sd  of Jura mountain's mature asp viper population 

#Coding figure 1.5A
x <- rnorm(n= 10, mean= mu, sd= sigma) #simulating a single sample (measuring 10 vipers)

#Coding figure 1.5B
reps <- 10^6 #repeating one million times
sample.means <- rep(NA, reps) #one million blank spots
for (i in 1:reps) { 
  sample.means[i] <- mean(rnorm(n= 10, mean= mu, sd= sigma))
} #fills each rep spot with a random sample, ultimately simulating a million samples

#Plotting Figure 1.5A and 1.5B
par(mfrow= c(1,2), las= 1) #formatting a 2x1 plot

hist(x, col= "gray", main= "Figure 1.5A", xlab= "Body length (cm)", ylab= "Frequency", las= 1) #creates figure 1.5A plot
abline(v= mu, lwd= 3, col= "red") #plots population mean on figure 1.5A
abline(v= mean(x), lwd= 3, col= "blue") #plots the sample mean on figure 1.5A

hist(sample.means, col= "gray", main= "Figure 1.5B", xlab= "Body length (cm)", ylab= "Density", nclass= 50, freq= F, las= 1) #creates figure 1.5B
abline(v= mu, lwd= 5, col= "red") #plots the population mean on figure 1.5B
abline(v= mean(sample.means), lwd= 5, col= "blue", lty= 2) #plots sample mean on figure 1.5B
#This plot shows that in figure 1.5B our distribution of sample means is more centered around the true population average length

sd(sample.means) #sd of the distribution of the estimates
#-------------------------------------------------------------------------------