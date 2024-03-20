# Bayesian Population Analysis using WinBUGS
# Chapter 1: Introduction

# 1.3 Sparrow Example: Binomial Distribution
# An observation process that is dominated by non-detection error, where mis-classification and double counts are absent
N <- 16 #population size of sparrows in the yard
p <- 0.4 #individual detection probability

rbinom(n= 1, size= N, prob= p) #simulating a single count, like field observation this changes every run

C <-rbinom(n= 10^6, size= N, prob= p) #simulating a one million counts to convey variation
mean(C)
var(C)
sd(C)
hist(C, breaks= 50, col= "Blue", main= "", xlab= "Sparrow Count", ylab= "Density", las= 1, freq= F) #frequency distributions of one million counts

#The typical count C is smaller than the actual population size N. Indeed, the mean of a binomial random variable and hence the expected count of sparrows, equals the product of N and p.
#However, the counts vary quite a lot, even under totally identical conditions.Hence, there is nothing intrinsically wrong or inferior with smaller counts. Smaller counts may simply result from the random nature of the counting process in the presence of imperfect detection. Thus, any count with p< 1 will automatically tend to vary from trial to trial. Unless p= 0 or p= 1, it is impossible to eliminate that variation by the sampling design or standardization (though other components of variation may be eliminated).
#Actually, not only the mean count but also the magnitude of the variation of counts is known from statistical theory. The variance is equal to the product of N, p, and 1−p.
