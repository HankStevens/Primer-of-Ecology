###################################################
### chunk number 1: 
###################################################
setwd("~/Documents/Projects/Book/draft/Chap09")
rm(list=ls())
badfonts = FALSE
.detach.stuff = function() {
 s1 = grep("package|Autoloads",search())
 nattach = length(search())
 xattach = search()[c(-1,-s1)]
 for (i in xattach)
   eval(substitute(detach(i),list(i=i)))
}
.detach.stuff()
## ps.options(family="NimbusSan")
ps.options(colormodel="cmyk")
options(width=75, SweaveHooks=list(fig=function() par(mar=c(5,4,.5,1.5))) ) 

## palette(gray((0:8)/8))
library(primer)
#source("DMBpplane.R")
require(deSolve);


###################################################
### chunk number 2: BSsuccFig
###################################################
trellis.par.set(canonical.theme(color=FALSE))
#buellsmallsucc <- read.xls("Stevens.xls")
#BSsucc <- melt(buellsmallsucc, id="AGE", measure=c("Annual","Perennial","Woody"))
data(BSsucc)
print(xyplot(value~ AGE, groups=variable, data=BSsucc, 
             type='l', ylab="Percent Cover",
             xlab="Age Since Abandonment (y)",
             auto.key=list(columns=3, space="top", lines=TRUE, points=FALSE)))


###################################################
### chunk number 3: 
###################################################
compcol <- function(t, y, params) {
  p1 <- y[1]; p2 <- y[2]
  with( as.list(params), {
    dp1.dt <- c1*p1*(1-p1) - m1*p1
    dp2.dt <- c2*p2*(1-p1-p2) -m2*p2 - c1*p1*p2
    return( list( c(dp1.dt, dp2.dt) ) )
    })
}


###################################################
### chunk number 4: 
###################################################
cd1 <- .5; cd2 <- 0.95
md1 <- .1; md2 <- .1


###################################################
### chunk number 5: ccnoint
###################################################
t <- 20; ps <- matrix(0, nrow=t+1, ncol=2); 
for(i in 1:t) ps[i+1,] <- {p1 <- ps[i,1] + cd1*(1-ps[i,1]) - md1*ps[i,1]
                          p2 <- ps[i,2] + cd2*(1-ps[i,2]) - md2*ps[i,2]
                          c(p1,p2)}
matplot(0:t+1, ps, type='b', ylab="Proportion of Sites", xlab="Time", xlim=c(0,t+1), ylim=c(0,1))


###################################################
### chunk number 6: ccint
###################################################
ps2  <- matrix(0, nrow=t+1, ncol=2); ps2[1,] <- ps[2,];
for(i in 1:t) ps2[i+1,] <- {p1 <- ps2[i,1] + cd1*ps2[i,1]*(1-ps2[i,1]) - md1*ps2[i,1]
                          p2 <- ps2[i,2] + cd2*ps2[i,2]*(1-ps2[i,2]) - md2*ps2[i,2] - cd1*ps2[i,2]
                          c(p1,p2)}
matplot(1:t+1, ps2[-(t+1),], type='b', ylab="Proportion of Sites", xlab="Time",xlim=c(0,t+1), ylim=c(0,1))


###################################################
### chunk number 7: 
###################################################
compcolM <- function(t, y, params) {
## This models S species with a competition--colonization tradeoff
  ## This requires 'params' is a list with named elements
  S <- params[["S"]]; D <- params[["D"]]
  with(params,
list( dpi.dt <- sapply(1:S, function(i) {
params[["ci"]][i] * y[i] * ( 1 - D - sum(y[1:i]) ) - 
   params[["m"]][i] * y[i] - sum( params[["ci"]][0:(i-1)] * y[0:(i-1)] * y[i] )  }
                       )
     )
       )
}


###################################################
### chunk number 8: 
###################################################
S <- 10; ranks <- 1:S
d <- .20; m=0.04;
geo.p <- expression(d*(1-d)^(ranks-1))
ci <- expression(m/(1-d)^(2*ranks-1) )


###################################################
### chunk number 9: ccM1
###################################################
par(mar=c(5,4,1,4), mgp=c(2,.75,0))
plot(ranks, eval(geo.p), type="b", ylab="Proportional Abundance", xlab="Rank", xlim=c(1,S))
par(new=TRUE); 
plot(ranks, eval(ci), type="b", axes=FALSE, ann=FALSE, lty=2)
axis(4)
mtext("Colonization Rates", side=4, line=2)


###################################################
### chunk number 10: 
###################################################
params <- list(ci=eval(ci), m=rep(m, S), S=S, D=0)
init.N <- rep(0.01, S); t=seq(1, 200, .1)
cc.out <- ode(init.N, t, compcolM, params)


###################################################
### chunk number 11: ccM2
###################################################
par(mgp=c(2,.75,0))
matplot(t, cc.out[, -1], type="l", ylab="Proportion of Habitat", xlab="Years", col=1)


###################################################
### chunk number 12: MultiHD
###################################################
params["D"] <- 0.25
init.N <- eval(geo.p)
cchd.out <- ode(init.N, t, compcolM, params)
matplot(t, cchd.out[,-1], type="l", lty=1:5, lwd=rep(c(3,1), each=5), col=1,
        ylab="Proportion of Habitat", xlab="Years")


###################################################
### chunk number 13: 
###################################################
succniche <- function(t, y, params) {
  S <- y[1]; E <- y[2]; M <- y[3]; R <- y[4]
    F <- max(c(0,1 - params["D"] - S - E - M - R))
  with(as.list(params), {
    dS = c*(S+R+M)*F - a*c*(M+E)*S - g*S - m*S
    dR = g*(S+M) - m*R
    dM = a*c*(M+E)*S + c*(S+R+M)*E - g*M - m*M
    dE = a*c*(M+E)*F - c*(S+R+M)*E - m*E
    return( list( c(dS, dE, dM, dR)))
})
}


###################################################
### chunk number 14: 
###################################################
params.suc <- c(a=7, c=0.2, g=5, m=0.04, D=0)


###################################################
### chunk number 15: 
###################################################
t=seq(0,50,.1)
init.suc <- c(S=0.01, E=0.03, M=0.0, R=0.00)
ccg.out <- data.frame(ode(init.suc, t, succniche, params.suc))


###################################################
### chunk number 16: SNPR1
###################################################

matplot(t, ccg.out[,-1], type="l", ylab="Relative Frequency", xlab="Time",
        ylim=c(0,1), col=1)
legend("right", colnames(ccg.out)[5:2], lty=4:1,  bty="n")


###################################################
### chunk number 17: 
###################################################
params.suc["g"] <- .1; exp(-.1)
ccg.out <- data.frame(ode(init.suc, t, succniche, params.suc))


###################################################
### chunk number 18: SNPR2
###################################################
matplot(t, ccg.out[,-1], type="l", ylab="Relative Frequency", xlab="Time",
        ylim=c(0,1), col=1)


###################################################
### chunk number 19: SNPR3
###################################################
params.suc["g"] <- 5;
params.suc["m"] <- abs(log(.9))
ccg.out <- data.frame(ode(init.suc, t, succniche, params.suc))
matplot(t, ccg.out[,-1], type="l", ylab="Relative Frequency", xlab="Time",
        ylim=c(0,1), col=1)


###################################################
### chunk number 20: SNPR4
###################################################
params.suc <- c(a=1, c=.7, g=.1, m=0.04, D=0)
ccg.out <- data.frame(ode(init.suc, t, succniche, params.suc))
matplot(t, ccg.out[,-1], type="l", ylab="Relative Frequency", xlab="Time",
        ylim=c(0,1), col=1)


###################################################
### chunk number 21: 
###################################################
params.suc1 <- c(a=10, c=.1, g=10, m=0.04, D=0)
Xstar1 <- ode(init.suc, 1:500, succniche, params.suc1)[500, -1]


###################################################
### chunk number 22: ED1
###################################################
params.suc1D <- c(a=10, c=.1, g=10, m=0.04, D=as.numeric(Xstar1["R"]))
t=1:100
ccg.out1 <- data.frame(ode(Xstar1, t, succniche, params.suc1D))
matplot(t, ccg.out1[,-1], type="l", col=1, ylab="Relative Frequency", xlab="Time")


###################################################
### chunk number 23: 
###################################################
params.suc2 <- c(a=1, c=1, g=.1, m=0.04, D=0)
Xstar2 <- ode(init.suc, 1:500, succniche, params.suc2)[500, -1]


###################################################
### chunk number 24: ED2
###################################################
params.suc2D <- c(a=1, c=.7, g=.1, m=0.04, D=as.numeric(Xstar1["R"]))
ccg.out2 <- data.frame(ode(Xstar2, t, succniche, params.suc2D))
matplot(t, ccg.out2[,-1], type="l", ylab="Relative Frequency", xlab="Time", col=1)
legend("topright", colnames(ccg.out2[5:2]), lty=4:1, bty='n')


###################################################
### chunk number 25: env
###################################################
years <- 100
t <- 1:years
variability = 4
env <- rnorm(years, m=0, sd=variability)
plot(t, env, type='l')


###################################################
### chunk number 26: 
###################################################
w.rare <- .5
w.comm <- 1


###################################################
### chunk number 27: 
###################################################
rho <- sd(env)


###################################################
### chunk number 28: rho
###################################################
hist.env <- hist(env,col='lightgray',main="Histogram of Environment")
abline(v=c(c(-rho, rho)/2), lty=3)
arrows(x0=-rho/2, y0=mean(hist.env[["counts"]]), 
       x1=rho/2, y1=mean(hist.env[["counts"]]), code=3, length=.1)
text(0, mean(hist.env[["counts"]]), quote(italic(rho)), adj=c(1.5,0), cex=1.5)
text(min(hist.env[["breaks"]]),mean(hist.env[["counts"]]),
     "Common sp.\ngrows best", adj=c(0,0))
text(max(hist.env[["breaks"]]),mean(hist.env[["counts"]]),
     "Rare sp.\ngrows best", adj=c(1,0))


###################################################
### chunk number 29: diff
###################################################
a.rare <-  (env+rho/2)*w.rare
a.comm <- -(env-rho/2)*w.comm


###################################################
### chunk number 30: Etime
###################################################
Es <- matrix(NA, nrow=years, ncol=2)
Es[,1] <- ifelse(a.rare > 0, a.rare, 0)
Es[,2] <- ifelse(a.comm > 0, a.comm, 0)
matplot(t, Es, type='l', col=1)


###################################################
### chunk number 31: Eenv
###################################################
matplot(env,Es, col=1)


###################################################
### chunk number 32: 
###################################################
d <- 0.1


###################################################
### chunk number 33: alpha
###################################################
alpha <- c(2*1e-5,1e-5)


###################################################
### chunk number 34: 
###################################################
Ns <- matrix(NA, nrow=years+1, ncol=2)
Cs <- matrix(NA, nrow=years, ncol=2)
Rs <- matrix(NA, nrow=years, ncol=2)


###################################################
### chunk number 35: 
###################################################
Ns[1,] <- c(1e3,1e5)


###################################################
### chunk number 36: 
###################################################
for(i in 1:years) Ns[i+1,] <- { 
  juveniles <- sum( exp(Es[i,]) * Ns[i,] )
  Cs[i,]  <-  alpha*juveniles
  #Cs[i,] <- log(juveniles/sum(d*Ns[i,]))
  Rs[i,] <- exp(Es[i,]-Cs[i,])
  (1-d) * Ns[i,] + Rs[i,]*Ns[i,]
                              }


###################################################
### chunk number 37: storepops
###################################################
matplot(c(0,t), Ns, type='b', log='y')


###################################################
### chunk number 38: chessonsim
###################################################
layout(matrix(1:3, nc=1))
par(mar=c(2,4,1,1) )
plot(t, env, type='l', ylab="Environment") 
matplot(t, Es, type='l',  ylab="E", col=1)
par(mar=c(5,4,1,1) )
matplot(c(t,years+1), Ns,  log='y', type='l', col=1, xlab="Time",
        ylab="Population Size (N)")


###################################################
### chunk number 39: 
###################################################
Nt1 <- rowSums(Ns)[1:years]
R.obs <- Ns[-1,]/Ns[-(years+1),]
Rmax <- apply(R.obs,2,max)


###################################################
### chunk number 40: covar
###################################################
nu <- log( t(Rmax / t(R.obs) ) )
colnames(nu) <- c("nu.rare", "nu.comm")
var(Nt1,nu)


###################################################
### chunk number 41: cv
###################################################
apply(Ns[round(years/2):years,], 2, function(x) sd(x)/mean(x) * 100)


###################################################
### chunk number 42: loadchesson
###################################################
source("chesson.r")
#help.start()


###################################################
### chunk number 43: 
###################################################
outA <- chesson(years=500, specialization=1, spread=.1)
outB <- chesson(years=500, specialization=5, spread=.67)
outA$overlap
outB$overlap


###################################################
### chunk number 44: chessonfunc1
###################################################
matplot(outB[['env']], outB[['Bs']], pch=1:2, xlim=c(-.6,.6),
        ylab="Fitness-independent Response", xlab="Environment")
matplot(outA[['env']], outA[['Bs']], pch=c(19,17), add=TRUE)
bs <- outB[["Bs"]][order(outB[["env"]]),]
rho.y <- apply(bs, 1, min )
envs <- sort(outB[["env"]])
polygon(envs, rho.y, col='grey', border=FALSE)


###################################################
### chunk number 45: chessonfunc2
###################################################
matplot(outB[['env']], outB[['Es']], pch=1:2, xlim=c(-.6,.6),
        ylab="Density-independent Reproduction", xlab="Environment")
matplot(outA[['env']], outA[['Es']], pch=c(19,17), add=TRUE)



