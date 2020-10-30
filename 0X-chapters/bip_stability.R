### Weighted NODF
#Almeida-Neto, M., Loyola, R.D., Ulrich, W., Guimaraes, P., Guimaraes, Jr., P.R. 2008. A consistent metric for nestedness analysis in ecological systems: reconciling concept and measurement. Oikos 117, 1227–1239

## Modularity
# Clauset, A.; Newman, M. E. J. & Moore, C. Finding community structure in very large networks, Phyisical Review E 2004, 70, 066111


bip_stability <- function( Spp.a = 10:50, Spp.p = 5:25,  reps=10, C=c(.05,.3), quant=FALSE){
  
  
  out.m <- data.frame(S=numeric(reps), C = numeric(reps),
                      modularity = numeric(reps),
                      nestedness = numeric(reps),
                      resilience.m=numeric(reps),
                      resilience.h=numeric(reps)
  )
  
  for(i in 1:reps) {
    Sp <- sample(Spp.p, size=1)
    Sa <- sample(Spp.a, size=1)
    S = Sp + Sa
    
    #random draw for C
    contc <- runif(1, min=C[1], max=C[2])
    
    # random draw for bipartite interactions, each = C
    interactions <- rbinom( n = Sa * Sp, size=1, prob=contc)
    
    # Bipartite Network
    bn <- matrix(interactions, nr=Sp, nc=Sa)
    rownames(bn) <- paste("plant", LETTERS[1:Sp], sep="")
    colnames(bn) <- paste("bug", LETTERS[1:Sa], sep="")
    
    # determine final web, based on "observed" interactions
    # remove rows and cols with no interactions
    rs <- rowSums(bn)
    cs <- colSums(bn)
    bn <- bn[rs>0, cs>0]
    
    # Find "observed" S.p and S.a
    Sp <- dim(bn)[1]
    Sa <- dim(bn)[2]
    S <- Sa+Sp
    
    # exponentially distributed interactions
    # what is average interaction strength (IS)?
    # expected IS is 1/rate
    x <- rexp(Sa*Sp, rate=1)
    
    # Does this step will create a negative relation between
    # S and average IS or make it independent?
    xp <- x/sum(x)
    
    # Bipartite Network with Exponentially distributed IS
    bne <- bn*xp
    n <- matrix(0, nr=S, nc=S)
    
    n[1:Sp, (Sp+1):S] <- bne
    n[(Sp+1):S, 1:Sp] <- t(bne)
    n2 <- n
    n2[1:Sp, (Sp+1):S] <- bne * -1
    
    
    # create negative density dep. for all species
    diag(n) <- -1*max(rowSums(n))
    diag(n2) <- -1*max(rowSums(n2))
    
    # quant. version of NODF  Almeida-Neto et al. (2008, 2011)
    WNODF <- bipartite::networklevel(bne, index=c("weighted NODF"))
    
    #connectance
    L <- sum( abs(bne) > 0 )
    connectance <- L/prod(dim(bne))
    
    # modularity: quantitative or based on adjacency matrix
    # The latter is much faster and is the default.
    if(!(quant)){
      m <- (n > 0) + 1 - 1
      g <- igraph::graph_from_adjacency_matrix(m)
      wtc <- igraph::cluster_walktrap(g)
      mod <- igraph::modularity(wtc)
    } else {
      mod <- bipartite::computeModules(bne, method="Beckett") @ likelihood
    }
    
    # resilience
    l1 <- rARPACK::eigs(n, k=2)$values[1]
    l1n <- rARPACK::eigs(n2, k=2)$values[1]
    
    # store measures
    out.m[i,] <- c(S, connectance, mod,
                   WNODF, -Re(l1), -Re(l1n) )
    
  }
  
  return( out.m )
}
