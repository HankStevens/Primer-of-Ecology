library(tidyverse)
library(lubridate)
butler <- read.csv("data/2255083.csv")
summary(butler)
head(butler)
butler <- butler %>% mutate(date = ymd(DATE),
                            year = year(DATE))

b.g <- butler %>% group_by(year )

b.s <- summarize(b.g,
                 PRCP.T = mean(PRCP, na.rm=TRUE),
                 TAVG.m = mean(TAVG, na.rm=TRUE),
                 TMID.m = mean( (TMAX+TMIN)/2, na.rm=TRUE),
                 TMAX.m = mean(TMAX, na.rm=TRUE),
                 TMIN.m = mean(TMIN, na.rm=TRUE)
                 )

ggplot(b.s, aes(year, TMIN.m)) + geom_line()

b.s2 <- subset(b.s, year > 1966)

ggplot(b.s2, aes(year, TMIN.m)) + geom_line()

data("sparrows")
s2 <- subset(sparrows, Year > 1966 & Year < 2003)
sb <- left_join(s2, b.s2, by=c("Year"="year"))

l <- sb$Count[-1]/sb$Count[-nrow(sb)]
i <- sb$Year[-1] - sb$Year[-nrow(sb)]

r <- log(l)

N <- 30
PSD <- rgamma(N, 2)^2
ACS <- fft(PSD,inverse = TRUE)
PSDeven <- c(PSD,PSD[N:2]) 
ACS <- fft(PSDeven,inverse = TRUE)/N
ACS <- Re(ACS)
plot(ACS, type ='b')

############
one_over_f <- function(N, alpha = 1, my.mean=0, my.SD=1){ 
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
  a <- reihe/( sd(reihe)/my.sd )
  reihe2 <- (a - mean(a)) + my.mean
  return(Re(reihe2)) # imaginary part is 0
}


#########
## Generate a series with known mean and variance

a <- one_over_f(100,2)
plot(a, type='l')
R1 <- exp(a4)
plot(R1, type='l')

