### Functions that perform spectral mimicry (Cohen et al. 1999, Cohen et al. 1998, Petchey et al. 1997)

one_over_f <- function(gamma=1, N=200){
  ## Generate 1/f noise with power = gamma
  ## after Petchey SAS code, etc.
  N.2 <- N/2
  sine.waves <- matrix(NA, nr=N, nc=N.2)
  steps=2*pi*(1:N)/N 
  phase <- runif(N.2, 0, 2*pi)
  for(i in 1:N.2) {
    freq <- i
    weight <- 1/(freq^gamma)
    y <- weight*sin(freq*steps+phase[i])
    sine.waves[,i] <- y
    
  }
  my.series <- rowSums(sine.waves)
  return(my.series)
}

spec_mimic <- function(X, Y=NULL, gamma=1){
  ## Based on J. E. Cohen, C. M. Newman, A. E. Cohen, O. L. Petchey, and A. Gonzalez. Spectral mimicry: a method of synthesizing matching time series with different {F}ourier spectra. Circuits, Systems and Signal Processing, 18:431–442, 1999.
  
  ## X is the raw data that we want to rearrange to match some 
  ## particular spectrum
  
  ## Y is the time series whose spectum we would like to mimic.
  
  ## gamma is the power of a simulated spectrum.
  ## If we do not specify Y, then this function will generate a random
  ## time series with with 1/f noise with this power.
  
  if(!is.null(Y)) {
    Z <- sort(X)[rank(Y)]
  } else {
    n <- length(X)
    Y2 <- one_over_f(N=n, gamma = gamma)
      Z <- sort(X)[rank(Y2)]
  }
  return(Z)
  
}

plot_f <- function(z){
  ## plotting function, which also estimates of the slope of
  ## log(amplitude) - log(freq)
  ## 
  n <- length(z)
  z1 <- z[1:(n/2)][-1]
  spz <- stats::spectrum(z1)
  lsa <- log(spz$spec)
  lf <- log(spz$freq)
  m <- lm(lsa~lf)
  par(mfrow=c(1,2))
  plot(1:n, z, type='l')
  plot(lf, lsa)
  abline(m)
  confint(m)[2,]
}


## Example
N <- 1000
gamma <- 1
X <- sample(1:N, size=N, replace=TRUE)
Y <- one_over_f(N, gamma = gamma)
Z <- spec_mimic(X, Y)
plot_f(Z)

Z2 <- spec_mimic(X, gamma=1)
plot_f(Z2)

################################################

## Pulled off of StackOverflow
  one_over_f_simon <- function(N, gamma = 0){ 
    ## This function generates a time series of 1/f noise 
    ## gamma is the color of the noise, and 
    ## gamma=0,1,2 are white, pink, and red noise respectively. 
    
    ## function to generate 1/f noise is by Simon C. on Stack Overflow
    
    f <- seq(from=0, to=pi, length.out=(N/2+1))[-c(1,(N/2+1))] # Fourier frequencies
    f_ <- 1 / f^gamma # Power law
    RW <- sqrt(0.5*f_) * rnorm(N/2-1) # for the real part
    IW <- sqrt(0.5*f_) * rnorm(N/2-1) # for the imaginary part
    fR <- complex(real = c(rnorm(1), RW, rnorm(1), RW[(N/2-1):1]), 
                  imaginary = c(0, IW, 0, -IW[(N/2-1):1]), length.out=N)
    ## Those complex numbers that are to be back transformed for 
    ## Fourier Frequencies 0, 2pi/N, 2*2pi/N, ..., pi, ..., 2pi-1/N 
    ## Choose in a way that frequencies are complex-conjugated and symmetric around pi 
    ## 0 and pi do not need an imaginary part
    reihe <- Re(fft(fR, inverse=TRUE)) # go back into time domain
    return(Re(reihe)) # clean up because imaginary part is 0
  }
  
  plot.f <- function(z){
    n <- length(z)
    z1 <- z[1:(n/2)][-1]
    spz <- stats::spectrum(z1)
    lsa <- log(spz$spec)
    lf <- log(spz$freq)
    m <- lm(lsa~lf)
    par(mfrow=c(1,2))
    plot(1:n, z, type='l')
    plot(lf, lsa)
    abline(m)
    summary( m )$coefficients
  }
)
