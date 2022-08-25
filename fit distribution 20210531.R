### Two COMPLEMENTARY methods

par(mfrow=c(2,2))
set.seed(333)
breaks = 20
normal.example = rnorm(n = 1000, mean =16, sd = 1)
hist(normal.example, breaks = breaks)


beta.example = rbeta(n=1000, shape1 = 2, shape2 = 8)
hist(beta.example, breaks = breaks)

gamma.example = rgamma(n=1000,shape = 1,rate = 2 )
hist(gamma.example, breaks = breaks)

weibull.example = rweibull(n=1000, shape=2, scale = 2)
hist(weibull.example, breaks = breaks)

par(mfrow=c(1,1))


## METHOD # 1 ==> fitdistrplus

## scenario 1
x= normal.example

library(fitdistrplus)
library(logspline)

## Step (1) Plot 
descdist( data = x , discrete = FALSE)
descdist(data = x, discrete = FALSE, boot=1000)



## Step (2) Fit

fitdist(x,"beta")
max(x)
#   values must be in [0-1] to fit a beta distribution

normal_ = fitdist(x, "norm")
weibull_ = fitdist(x, "weibull")
gamma_ = fitdist(x, "gamma")

plot(normal_)
plot(weibull_)
plot(gamma_)

## Step (3) Estimate parameters

print(normal_)
print(weibull_)
print(gamma_)

summary(normal_)
summary(weibull_)
summary(gamma_)

## scenario 2
x= beta.example

library(fitdistrplus)

## Step (1) Plot 
descdist( data = x , discrete = FALSE)
descdist(data = x, discrete = FALSE, boot=1000)



## Step (2) Fit

fitdist(x,"beta")
max(x)
#   values must be in [0-1] to fit a beta distribution
beta_ = fitdist(x, "beta")
normal_ = fitdist(x, "norm")
weibull_ = fitdist(x, "weibull")
gamma_ = fitdist(x, "gamma")

plot(beta_)
plot(normal_)
plot(weibull_)
plot(gamma_)

## Step (3) Estimate parameters

print(beta_)
print(normal_)
print(weibull_)
print(gamma_)

summary(beta_)
summary(normal_)
summary(weibull_)
summary(gamma_)
 
## scenario 3
x= rpois(1000, lambda = 2)
hist(x)

library(fitdistrplus)

## Step (1) Plot 
descdist( data = x , discrete = T)
descdist(data = x, discrete = T, boot=1000)



## Step (2) Fit

fitdist(x,"pois")
pois__ = fitdist(x,"pois")
plot(pois__)
 

## Step (3) Estimate parameters


print(pois__)

summary(pois__)


### METHOD 2: library(fitur)

library(fitur)
library(actuar)
fitur::fit_dist_addin()
