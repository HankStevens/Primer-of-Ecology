Here is a function that will generate a colored time series.
```{r}
############
source("color-funs.R")


```{r}
## Petchey 2000

```


Here we use the function to generate a time series based on our sparrow data. It will have the same mean and variance as the log-transformed data we used previously. 
```{r}
#########

This generates a time series where sparows are likely to experience several good or bad years in a row. 

Interested parties may want to try another new function.
```{r}
myForLoop.f <- function(mu, sigma, years, initial.N, gamma=1) {
  # select all R at random from 
  lrR <- one_over_f(N=years, gamma=gamma)
  rR <- exp(lrR[1:years])
  # create a vector to hold N
  N <- numeric(years+1)
  # give it an initial population size
  N[1] <- initial.N
  # Do the for-loop
  for( t in 1:years ) {
    # project the population one time step
    N[t+1] <-  N[t] * rR[t]
  }
  # return the vector of N
  N
} 
```

Our new simulations.
```{r}
sims=1e3; years=100
set.seed(2)
outmat.f <- replicate(sims,   
                      expr=myForLoop.f(mu=mu, sigma=sigma, years=years, initial.N=43)
)
matplot(0:years, outmat.f[, 1:10], type='l', log='y')
Nf.2053 <- outmat.f[51,]
quantile(Nf.2053, prob=c(0.05, .95) )
quantile(N.2053, prob=c(0.05, .95) )

#probability of extinction within 'years' years.
# extinction = N <- minimum population size.
min.pop.size <- 0.0001

pr.e.100 <- mean( apply(outmat.f, 2, function(x) {any(x<min.pop.size)}) )

### time to extinction, given extinction
# which.max(x) finds the first occurrence of a value that is the maximum
# here maximum = TRUE i.e. 1
a <- apply(outmat.f, 2, function(x) {which.max(x < min.pop.size)})
t.e <- ifelse(a==1,NA, a)
t.e
hist(t.e)
```