Here is a function that will generate a colored time series.
```{r}
############
one_over_f <- function(N, alpha = 1, my.mu=0, my.sigma=1){ 
  ## This function generates a time series of 1/f noise 
  ## with a user defined mean and SD.
  ## alpha is the color of the noise, and 
  ## alpha=0,1,2 are white, pink, and red noise respectively. 
  
  ## function to generate 1/f noise is by Simon C. on Stack Overflow
  ## code to standardize the series to a defined mean and SD is
  ## by MHHS
  
  f <- seq(from=0, to=pi, length.out=(N/2+1))[-c(1,(N/2+1))] # Fourier frequencies
  f_ <- 1 / f^alpha # Power law
  RW <- sqrt(0.5*f_) * rnorm(N/2-1) # for the real part
  IW <- sqrt(0.5*f_) * rnorm(N/2-1) # for the imaginary part
  fR <- complex(real = c(rnorm(1), RW, rnorm(1), RW[(N/2-1):1]), 
                imaginary = c(0, IW, 0, -IW[(N/2-1):1]), length.out=N)
  ## Those complex numbers that are to be back transformed for 
  ## Fourier Frequencies 0, 2pi/N, 2*2pi/N, ..., pi, ..., 2pi-1/N 
  ## Choose in a way that frequencies are complex-conjugated and symmetric around pi 
  ## 0 and pi do not need an imaginary part
  reihe <- Re(fft(fR, inverse=TRUE)) # go back into time domain
  a <- reihe/( sd(reihe)/my.sigma )
  reihe2 <- (a - mean(a)) + my.mu
  return(Re(reihe2)) # imaginary part is 0
}
one_over_f <- function(N, alpha = 1, my.mu=0, my.sigma=1){ 
  ## This function generates a time series of 1/f noise 
  ## with a user defined mean and SD.
  ## alpha is the color of the noise, and 
  ## alpha=0,1,2 are white, pink, and red noise respectively. 
  
  ## function to generate 1/f noise is by Simon C. on Stack Overflow
  ## code to standardize the series to a defined mean and SD is
  ## by MHHS
  
  f <- seq(from=0, to=pi, length.out=(N/2+1))[-c(1,(N/2+1))] # Fourier frequencies
  f_ <- 1 / f^alpha # Power law
  RW <- sqrt(0.5*f_) * rnorm(N/2-1) # for the real part
  IW <- sqrt(0.5*f_) * rnorm(N/2-1) # for the imaginary part
  fR <- complex(real = c(rnorm(1), RW, rnorm(1), RW[(N/2-1):1]), 
                imaginary = c(0, IW, 0, -IW[(N/2-1):1]), length.out=N)
  ## Those complex numbers that are to be back transformed for 
  ## Fourier Frequencies 0, 2pi/N, 2*2pi/N, ..., pi, ..., 2pi-1/N 
  ## Choose in a way that frequencies are complex-conjugated and symmetric around pi 
  ## 0 and pi do not need an imaginary part
  reihe <- Re(fft(fR, inverse=TRUE)) # go back into time domain
  return(Re(reihe)) # imaginary part is 0
}
```


```{r}
## Petchey 2000

```


Here we use the function to generate a time series based on our sparrow data. It will have the same mean and variance as the log-transformed data we used previously. 
```{r}
#########
sp <- stats::spectrum(ldeaths)
str(sp)
y <- log(sp$spec^2)
plot(sp$freq, y)
## Generate a series with known mean and variance
years <- 50
a <- one_over_f(N=years, alpha=2, my.mu=0.007, my.sigma=0.34 )
as <- stats::spectrum(a)
lsa <- log(as$spec^2)
plot(as$freq, lsa)
m <- lm(lsa~as$freq)
summary( m )
mean(a);sd(a)
plot(a, type='l')
R1 <- exp(a)
qplot(1:years, R1, geom=c('line','point'))
```

This generates a time series where sparows are likely to experience several good or bad years in a row. 

Interested parties may want to try another new function.
```{r}
myForLoop.f <- function(mu, sigma, years, initial.N, alpha=1) {
  # select all R at random from 
  lrR <- one_over_f(N=years*2, alpha = 1, my.mu=mu, my.sigma=sigma)
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
set.seed(3)
outmat.f <- replicate(sims,   
                      expr=myForLoop.f(mu=mu, sigma=sigma, years=years, initial.N=43)
)
matplot(0:years, outmat.f[, 1:10], type='l', log='y')
Nf.2053 <- outmat.f[51,]
quantile(Nf.2053, prob=c(0.05, .95) )
quantile(N.2053, prob=c(0.05, .95) )
```