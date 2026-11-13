DemoInfo <- function(mat) {
  ## Right eigenvector: stable stage distribution ----
  eR <- eigen(mat)
  iR <- which.max(Re(eR$values))   # Perron root: real, largest real part
  L1 <- Re(eR$values[iR])
  w  <- Re(eR$vectors[, iR]); w <- w / sum(w)
  
  ## Left eigenvector: reproductive value ----
  eL <- eigen(t(mat))
  iL <- which.max(Re(eL$values))   # index separately; ordering can differ
  v  <- Re(eL$vectors[, iL])
  v  <- if (abs(v[1]) > 1e-12) v / v[1] else v / sum(v)
  
  ## Sensitivity and elasticity ----
  sens <- outer(v, w) / sum(v * w)
  elas <- (mat / L1) * sens
  
  list(lambda = L1, SSD = w, RV = v, Sensitivities = sens,
       Elasticities = elas, PPM = mat)
}
