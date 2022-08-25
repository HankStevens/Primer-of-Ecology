## ----web, echo=FALSE, fig.cap="A food web is a map of feeding relations. Here, snails (S) and invertebrates (I) feed on algae, and fish, F, feed on invertebrates and algae. We can represent this as a graph or a matrix.", out.width="40%", fig.show='hold', fig.ncol=2, fig.height=5, fig.width=5----
gid <- graph(c('S','A', 'A','S', 
               'I','A', 'A','I', 
               'I','F', 'F','I',
               'A','F', 'F','A'))
signs <- c("-",'+', "-",'+', "-",'+', "-",'+')
plot(gid, layout=matrix(c(0, .5, 1, .6, 0.5, 0.1, 0.5, .85), nc=2),
     edge.curved=.5,edge.label=signs, edge.label.font=2,  
     edge.label.cex=1.5, vertex.size=60, 
     margin=c(0,0,0,0), rescale=TRUE)


## ----motifs, echo=FALSE, fig.cap="There are 13 different possible motifs for three species. These links record only consumption. For instance, S1 is a simple linear food chain and S4 is a single species that is consumed by two other species. motifs D1-D8 all include mutual consumption. This is likely to occur between species in the same guild, and when adults of one feed on juveniles of another."----
include_graphics("figs/motifs.png")


## ----PLABE, echo=FALSE, fig.cap="Three of the six webs investigated by Pimm and Lawton (1977). Left to right, these correspond to Pimm and Lawton (1977) Figs. 1A, E, and B. All basal species exhibit negative density dependence.", out.width=c("50%","50%","50%"), fig.show='hold'----
# fig.ncol=3, fig.height=10, fig.width=5
resc <- TRUE
gA <- graph(c(1,1, 1,2,2,1, 2,3,3,2, 3,4,4,3))
signs <- c("-", '+',"-", '+',"-", '+',"-")
plot(gA, layout=matrix(c(.5,.5,.5,.5, 0,.33,.67,1), nc=2), edge.curved=.5,      edge.label=signs, edge.label.font=1,  edge.label.cex=1, 
     vertex.size=20, edge.loop.angle=pi/4,
     margin=c(0,0,0,0), rescale=resc)

gE <- graph(c(1,1, 2,2, 3,3, 1,4,4,1, 2,4,4,2, 3,4,4,3))
signs <- c("-","-","-", '+',"-", '+',"-", '+',"-")
plot(gE, layout=matrix(c(0,.5,1,.5, 0,0,0,.5), nc=2), edge.curved=.5,      edge.label=signs, edge.label.font=1,  edge.label.cex=1, 
     vertex.size=30, edge.loop.angle=pi/4,
     margin=c(0,0,0,0), rescale=resc)

gB <- graph(c(1,1, 1,2,2,1, 
              2,3,3,2, 
              3,4,4,3, 
              2,4,4,2))
signs <- c("-", '+',"-", 
           '+',"-", 
           '+',"-", 
           '+',"-")
Bl <- matrix(c(.3,.5,.7,.3, 0,.33,.67,1), nc=2)
Bl <- matrix(c(.5,.5,.51,.49, 0,.33,.67,1), nc=2)
plot(gB, layout=Bl, 
     edge.curved=.5,edge.label=signs, edge.label.font=1,  edge.label.cex=1, 
     vertex.size=30, edge.loop.angle=pi/4,
     margin=c(0,0,0,0), rescale=resc)


## ----may1, fig.cap="Relation between the average interaction strength and the number of species able to coexist (here directed connectance is $C_D = 0.3$). The line represents the maximum number of species that are predicted to be able to coexist at equilibrium. Fewer species could coexist, but, on average, more species cannot coexist at equilibrium.", fig.width=5, fig.height=5, results='hide', echo=FALSE, out.width="60%"----
S <- 40; C=.3
par(mar=c(3,3,0,0), mgp=c(2,.5, 0))
curve(C*x^ -2, .1, 1, ylab=expression("Number of Species ("*italic(S)*")"),
      xlab = expression("Interaction Strength ("*italic(I)*")") )
text(.6, 20, "Unstable Webs", srt=45)
text(.2, 3.75, "Stable Webs" , srt=-45, cex=.8)


## ---------------------------------------------------------------------------------------
Aq = matrix(c(
  -1,   -10,   0,     0,
  0.1, 0,   -10,  0,
  0,    0.1,   0,    -10,
  0,      0,    0.1,   0),
  nrow=4, byrow=TRUE)


## ---------------------------------------------------------------------------------------
library(diagram)
plotmat(Aq,relsize=.8)


## ---------------------------------------------------------------------------------------
S <- nrow(Aq)


## ---------------------------------------------------------------------------------------
M <- Aq * runif(S^2)


## ---------------------------------------------------------------------------------------
eM <- eigen(M)[["values"]]


## ---------------------------------------------------------------------------------------
deM <- max( Re(eM) )


## ---------------------------------------------------------------------------------------
intraNDD <- sqrt(sum(diag(M)^2)/S)


## ---------------------------------------------------------------------------------------
diag(M) <- 0
IS <- sqrt( sum(M^2)/(S*(S-1)) )


## ----PL2, echo=FALSE, results='hide'----------------------------------------------------
pimmlawton <- function(mat, N=1, omni.i=NA, omni.j=NA, omega=NULL){
  S <- nrow(mat)
  if(is.na(omni.i)) {
  out <- matrix(NA, nrow=N, ncol=4)
    colnames(out) <- c("DomEig", "Im", "IntraDD", "I")
  for(n in 1:N) out[n,] <- {
  M = runif(S^2) * mat
  ## Do eigenanalysis
  eigs <- eigen(M)[["values"]]
  mx <- which.max(Re(eigs))
  deM = Re(eigs)[mx]
  deMi = Im(eigs)[mx]
  intraNDD <- sqrt(sum(diag(M)^2)/S)
  diag(M) <- 0
IS = sqrt( sum(M^2)/(S*(S-1)) )
 c(deM, deMi, intraNDD, IS)
  } } else {
  out <- matrix(NA, nrow=N, ncol=5)
    colnames(out) <- c("DomEig", "Im", "IntraDD", "I", "I.omni")
  for(n in 1:N) out[n,] <- {
  M = runif(S^2) * mat
  ## Adjust for McCann type omnivory
  if(!is.null(omega))  {M[omni.i,omni.j] <- omega*M[omni.i+1,omni.j]
                      M[omni.i+1,omni.j] <- (1-omega)*M[omni.i+1,omni.j]
                      M[omni.j,omni.i] <- omega*M[omni.j,omni.i+1]
                      M[omni.j,omni.i+1] <- (1-omega)*M[omni.j,omni.i+1]}
  ## Do eigenanalysis
  eigs <- eigen(M)[["values"]]
  mx <- which.max(Re(eigs))
  deM = Re(eigs)[mx]
  deMi = Im(eigs)[mx]
  intraNDD <- sqrt(sum(diag(M)^2)/S)
  diag(M) <- 0
IS = sqrt( sum(M^2)/(S*(S-1)) )
    omnivory <- sqrt(mean(c(M[omni.i,omni.j],M[omni.j,omni.i])^2))
 c(deM, deMi,intraNDD, IS, omnivory)
  }
}
    return(as.data.frame(out))
}


## ---------------------------------------------------------------------------------------
args(pimmlawton)


## ---------------------------------------------------------------------------------------
set.seed(1)
pimmlawton(Aq)


## ---------------------------------------------------------------------------------------
out.A <- pimmlawton(Aq, N=2000) 


## ---------------------------------------------------------------------------------------
summary(out.A)


## ----pairsA, fig.cap="Perturbations at the equilibrium tend to dissipate more rapidly (more negative dominant eigenvalues) with greater intraspecific negative density dependence (IntraDD) and greater interspecific interaction strength (I). This graph also demonstrates the independence of IntraDD and I in these simulations."----
out.A$I.rel <- with(out.A, log(I/IntraDD) )
pairs(out.A)
# out.A$Idd <-  cut(out.A$IntraDD, 
#          quantile(out.A$IntraDD), include.lowest=TRUE )
# ggplot(data=out.A, aes(I, DomEig)) + geom_point() +
#   geom_smooth() + facet_wrap(~ Idd)


## ---------------------------------------------------------------------------------------
  RT.A <- -1/out.A[["DomEig"]]
  summary(RT.A)


## ---------------------------------------------------------------------------------------
  sum( RT.A > 150 ) / 2000


## ----histA, fig.cap="*A linear four-species food chain has is only somewhat likely to have reasonble return times.*"----
A.fast <-RT.A[RT.A < 150]
histA <- hist(A.fast, breaks=seq(0,150, by=5), main=NULL)


## ---------------------------------------------------------------------------------------
Eq = matrix(c(
  -1,   0,   0, -10,
  0,   -1,   0, -10,
  0,    0,   -1, -10,
  0.1, 0.1, 0.1,  0),
  nrow=4, byrow=TRUE)
# to confirm we got it right:
# plotmat(Eq, relsize=.8)


## ---------------------------------------------------------------------------------------
out.E <- pimmlawton(Eq, N=2000)
summary(out.E)


## ----E1, fig.cap="For a food chain with two levels, and three basal species, perturbation growth rate (lambda_1) declines with increasing intraspecific negative density dependence (IntraDD) and is unrelated to predator-prey interaction strengths.", out.width=9----
layout(matrix(1:2,nr=1))
plot(DomEig ~ IntraDD, data=out.E)
plot(DomEig ~ I, data=out.E)


## ----histE, fig.cap="*One predator and three prey, each with negative density dependence results in very stable webs.*"----
RT.E <- -1/out.E[["DomEig"]]
E.fast <-RT.E[RT.E < 150]
histE <- hist(E.fast, breaks=seq(0,150, by=5), main=NULL)


## ---------------------------------------------------------------------------------------
Bq <- matrix(c(
  -1, -10,   0,     0,
  0.1,  0,   -10,   -2,
  0,    0.02,   0,    -10,
  0,    0.1,    0.1,   0),
  nrow=4, byrow=TRUE)
#plotmat(Bq, relsize=.8)
# recall we have to tell the function which links are omnivorous
out.B <-  pimmlawton(Bq, N=2000, omni.i=2, omni.j=4)
summary(out.B)


## ---------------------------------------------------------------------------------------
out.B.stable <- out.B[out.B[,"DomEig"]<0,]
pairs(out.B.stable)


## ----histB, fig.cap="*Of the relatively few omnivorous webs that were qualitatively stable, a moderate fraction had low return times.*"----
RT.B <- -1/out.B[["DomEig"]]
B.fast <- RT.B[RT.B < 150 & RT.B > 0 ]
out.B.fast <- out.B[RT.B < 150 & RT.B > 0,]
histB <- hist(B.fast, breaks=seq(0,150, by=5),
              main=NULL)


## ---------------------------------------------------------------------------------------
ret.time <- data.frame(web=rep(LETTERS[1:2], each=2000),
                       return.time = c(RT.A, RT.B)) %>% 
  filter(return.time>0)


## ----histsAB, fig.cap="*Comparing the distributions of return times for chain **A** and **B**. 'Density' is probability density.*"----
ggplot(ret.time, aes(return.time, fill=web)) + 
  geom_histogram(aes(y = ..density..),
                 alpha=.5, position="identity") + 
  geom_density(aes(fill=NULL, colour=web)) + 
  scale_x_log10() 


## ----eval=FALSE-------------------------------------------------------------------------
## pkgs <- c("bipartite", "igraph", "rARPACK", "lavaan", "semPlot")
## # install.packages(pkgs)
## invisible( # to hide messages
##   lapply(pkgs, library, character.only=TRUE)
## )


## ----TandF, echo=FALSE, fig.cap="*Possible effects of species richness, connectance, modularity, and nestedness on network resilience.*", out.width="50%"----
tf <-
"digraph{
  rankdir = BT;
//node[shape=box];
//{rank=max; R;}
//{rank=min; S; C;}
// oval is the default node shape
S -> M; 
S -> N;
S -> R
C -> M;
C -> N;
C -> R;
M -> R;
N -> R;
}"

grViz(tf) %>%
  export_svg %>% charToRaw %>% rsvg_png("TandF.png")

include_graphics("TandF.png")


## ----cars-------------------------------------------------------------------------------
# The number of species of plants
Sp <- 5
# The number of species of pollinators
Sa <- 12
# the total number of species
S <- Sa+Sp

# The connectance in this web. 
contc <- .4


## ---------------------------------------------------------------------------------------
# size is the number of coins to flip at once.
rbinom(n=10, size = 1, prob=0.5 )


## ---------------------------------------------------------------------------------------
# size is the number of coins to flip at once.
rbinom(n=1, size = 10, prob=0.5 )


## ---------------------------------------------------------------------------------------
# Create interactions at random, drawing 0's and 1's with a probability of contc=0.5
interactions <- rbinom( n = Sa * Sp, size=1, prob=contc)


## ---------------------------------------------------------------------------------------
bip.network <- matrix(interactions, nrow=Sp, ncol=Sa)
rownames(bip.network) <- paste("plant", LETTERS[1:Sp], sep="")
colnames(bip.network) <- paste("bug", LETTERS[1:Sa], sep="")
#...and look at it
bip.network


## ---------------------------------------------------------------------------------------
# Exponentially distributed of interaction strengths
# Random numbers from and exponential distribution
xp <- rexp(Sa*Sp)

# express these as a probability distribution (divide by the sum)
xp <- xp/sum(xp)

# assign those exponential probabilities to our interactions,
# creating the Quantitative network
bip.networkQ <- bip.network * xp


## ----toys, fig.cap="*Illustrations of a hypothetical plant-pollinator network. The two matrices organize the species to highlight either nestedness (left) or modularity (right).*"----
{
  # library(bipartite) if needed
  # visweb and plotweb are in the bipartite package
  layout(matrix(c(1,2,3,4), nc=2, byrow=TRUE) )
  hist(bip.networkQ, main="", xlab="interaction strength")
  bipartite::plotweb(bip.networkQ, method='normal') 
  bipartite::visweb(bip.networkQ, type='nested')
  bipartite::visweb(bip.networkQ, type='diagonal')
  }


## ----network_properties-----------------------------------------------------------------
# Nestedness:
# quant. version of NODF  (Almeida-Neto et al. 2008, 2011)
networklevel(bip.networkQ, index=c("weighted NODF"))
metaComputeModules(bip.networkQ, N=5, method="Beckett") @ likelihood


## ----network_stab-----------------------------------------------------------------------
n <- matrix(0, nr=S, nc=S) 
# upper triangular part
n[1:Sp, (Sp+1):S] <- bip.networkQ
# lower triangular part
n[(Sp+1):S, 1:Sp] <- t(bip.networkQ) 
# Add negative density dependence that is the 
# same magnitude as other species interactions
diag(n) <- -1*max(rowSums(n))

#Calculate stability
Re(eigen(n)$values[1])


## ----bip_stability----------------------------------------------------------------------
#library(primer)
b <- bip_stability(Spp.a = 16:60, Spp.p = 8:30, C=c(.05, 0.5), reps=5000, quant=FALSE)


## ---------------------------------------------------------------------------------------
summary(b)


## ---------------------------------------------------------------------------------------
pairs(b)


## ---------------------------------------------------------------------------------------
out.tr <- (b)^.75


## ---------------------------------------------------------------------------------------
dat <- data.frame( scale(out.tr) ) 


## ---------------------------------------------------------------------------------------
# library(lavvan)
# Step 1: Specify model
model.m <- 'resilience.m ~ nestedness + modularity + S + C
          nestedness ~  S + C
          modularity ~ S + C '

# Step 2: Estimate model of mutualism
mod.m.fit <- sem(model.m, data=dat)
coef(mod.m.fit)


## ---------------------------------------------------------------------------------------
# Step 1: Specify model
model.h <- 'resilience.h ~ nestedness + modularity + S + C
          nestedness ~ S + C
          modularity ~ S + C'

# Step 2: Estimate model of herbivory
mod.h.fit <- sem(model.h, data=dat)
#coef(mod.h.fit)  


## ----paths, fig.cap="*Given the same bipartite network structure and the identical interaction strengths, mutualistic (left) and antagonistic webs (right) show slightly different relations among emergent properties. Recall also that resilience is about *", fig.show='hold', out.width="50%"----
# library(semPlot)
semPaths(mod.m.fit, what='path', whatLabels="std", 
         curve=1, curvature=1, sizeMan=10,
         edge.label.cex = 1.5, nCharNodes=4,
         edge.label.position=.7, layout='tree2')
semPaths(mod.h.fit, what='path', whatLabels="std", 
         curve=1, curvature=1, sizeMan=10,
         edge.label.cex = 1.5, nCharNodes=4,
         edge.label.position=.7, layout='tree2')



## ----plot1, fig.cap="*Example of a random graph using the same number of vertices and links as in plot 1 (Beiler et al. (2015)).*", out.width="70%"----
# node clustering coefficient: 0.84
# for Plot 1
V <- 26
E <- 258

gr <- erdos.renyi.game(V, E, type = "gnm")
# graph "order" and "size"
gorder(gr); gsize(gr)
plot(gr, vertex.label= NA, edge.arrow.size=0.02,vertex.size = 0.5, xlab = "Random Network: G(V,E) model")


## ---------------------------------------------------------------------------------------
mean(transitivity(gr, type = "local"))


## ----nccHist, fig.cap="*Observed node clustering appears unuusally high, relatively to that of random networks.*"----
ncc.sims <- replicate(100, { 
  gr <- erdos.renyi.game(V, E, type = "gnm")
  mean(transitivity(gr, type = "local"))
  } )
qplot(x=ncc.sims, geom="histogram") + 
  geom_vline(xintercept = 0.84, linetype=2)


## ----plot1sf, fig.cap="*Emergent relation between node clustering and node degree in a scale-free random network, based on the number of nodes and links in plot 1.*"----
# node clustering coefficients (per node)
ncc <- transitivity(gr, type = "local")
# degrees of each node
d <- degree(gr)
qplot(d, ncc, geom="jitter")


## ---------------------------------------------------------------------------------------
# vertices
V <- c(26, 13, 31, 27, 27, 41)
# egdes
E <- c(258, 66, 283, 254, 177, 434)
# proportion of trees colonized
pV <- c(.846, .5, .875, .846, .778, .889)
# mean node degree
mean.deg <- c(18, 695, 17.15, 17.52, 12.21, 20.19)
# mean node clustering coef.
mean.ncc <- c(.84, .63, .84, .85, .83, .8)
plot <- 1:6
df <- data.frame(plot=plot, V=V, E=E, pV=pV, mean.deg=mean.deg, mean.ncc=mean.ncc)
df


## ----echo=FALSE, eval=FALSE-------------------------------------------------------------
## my_ER_Stats <- function(ve){
##   # ve is a two-element vector of vertices and edges totals
## 
##   V <- ve[1];  E <- ve[2]
##   gr <- erdos.renyi.game(V, E, type = "gnm")
## 
##   # betweenness centrality
##   bc <- igraph::betweenness(gr, directed = FALSE, weights = NULL,
##   nobigint = TRUE, normalized = FALSE)
##   bcbar <- mean(bc)
## 
##   # Network clustering coefficient
##   gcc <- transitivity(gr) # global
## 
##   ncc <- transitivity(gr, type = "local")
## # mean node clustering coefficient
##  nccbar.o <- mean(ncc) # cf 0.84
## # degrees of each node
## d <- igraph::degree(gr)
## # mean degree
## dbar.o <- mean(d) # cf 18
## return( c(bcbar.op=bcbar, gcc.o=gcc, nccbar.o=nccbar.o, dbar.o=dbar.o))
## }


## ---------------------------------------------------------------------------------------
exp.out <- 2.5


## ----allplotsf, fig.cap="*Emergent relations between node clustering and node degree in scale-free random networks, based on the number of nodes and links in each plot.*"----
df2 <- apply(df[,1:3], 1, function(ve){
  #make a scale-free random graph, given V and E
  gsf <- sample_fitness_pl(ve[2], ve[3], exponent.out= exp.out,
                           exponent.in = -1, loops = FALSE,
                           multiple = FALSE, 
                           finite.size.correction = TRUE)
  # calculate prpoerties
  ncc <- transitivity(gsf, type = "local")
  deg <- igraph::degree(gsf)
  
  #combine in a matrix
  m <- matrix(c(plot=rep(ve[1], length(deg)), ncc, deg), ncol=3)
  m
 })

# bind together both sites and label columns
df.sf <- as.data.frame( do.call(rbind, df2) )
names(df.sf) <- c("plot", "ncc", "deg")

# plot the relation between node clustering coefficient
# and node degree
ggplot(df.sf, aes(deg, ncc, colour=as.factor(plot)) ) +
  geom_point()  + geom_smooth( method="lm", se=FALSE) + 
  scale_y_log10() 

