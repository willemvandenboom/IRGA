VAMP <- function(y, X, lambda, psi, sigma.sq = NULL, rho = 1, max.iter = 100, min.rho = 1, verbose = FALSE, a.0 = 1, b.0 = 1) {
  # Vector approximate message passing (VAMP)
  # This implements Algorithm S1 from the Supplementary Material to `Approximating high-dimensional posteriors with nuisance parameters via integrated roated Gaussian approximation`.
  
  # The prior implemented here is:
  # Pi(beta) = lambda * N(0, psi) + (1-lambda) delta(0)
  # which is a spike-and-slab prior.
  
  # The parameter `rho` controls the dampening in the VAMP updates.
  
  y <- as.vector(y)
  n <- length(y) # Number of observations
  X <- as.matrix(X)
  p <- ncol(X) # Number of predictors
  
  
  # Check whether sigma.sq is known
  s2unknown <- is.null(sigma.sq0 <- sigma.sq)
  
  
  ## Step 1
  # Compute the SVD of X = UDV'
  tmp <- svd(X)
  U <- tmp$u
  diag_D <- tmp$d # Diagonal of D
  V <- tmp$v
  rm(tmp)
  
  
  ## Precompute quantities that are used repeatedly in Step 3c
  UtY <- crossprod(U, y) # U^T y
  D_Vt <- diag_D * t(V) # D V^T
  
  
  ## Step 2
  # Initialize r_0 and t_0^2
  r <- numeric(p)
  t.sq <- psi
  
  
  ## Step 3
  # The VAMP algorithm
  if(verbose) pb <- txtProgressBar(max = max.iter, style = 3)
  for(k in 0:max.iter) {
    
    ## Step 3a
    incProb <- 1 / (1 + exp( log(1-lambda) - log(lambda) +
      dnorm(x = r, sd = sqrt(t.sq), log = TRUE) - dnorm(x = r, sd = sqrt(psi + t.sq), log = TRUE)
    )) # Inclusion probabilities
    
    nonzero_mean <- r / (t.sq / psi + 1) # Mean of beta given that beta is nonzero
    
    beta.new <- incProb * nonzero_mean
    s.sq.new <- mean(incProb / (1/psi + 1/t.sq) + nonzero_mean^2 * incProb * (1-incProb)) # Follows from the law of total variance
    
    # Dampened updates
    if(k == 0) {
      beta <- beta.new
      s.sq <- s.sq.new
    } else {
      beta <- rho * beta.new + (1-rho) * beta
      s.sq <- rho * s.sq.new + (1-rho) * s.sq
    }
    
    
    ## Step 3b
    t.tilde.sq <- 1/(1/s.sq - 1/t.sq)
    
    r.tilde <- (t.sq * beta - s.sq * r) / (t.sq - s.sq)
    
    
    # Update sigma.sq if sigma.sq is unknown
    if(s2unknown) sigma.sq <- (2*b.0+sum((y-X%*%beta)^2))/(2*a.0+n)
    
    
    ## Step 3c
    tmp <- 1/(sigma.sq / diag_D / t.tilde.sq + diag_D) # (\sigma^2 D^{-1} / t.tilde^2 + D)^{-1}
    
    beta.tilde.new <- r.tilde + V %*% (tmp * (UtY - D_Vt %*% r.tilde))
    s.tilde.sq.new <- t.tilde.sq * ( 1 - sum(diag_D * tmp)/p )
    
    # Dampened updates
    if(k == 0) {
      beta.tilde <- beta.tilde.new
      s.tilde.sq <- s.tilde.sq.new
    } else {
      beta.tilde <- rho * beta.tilde.new + (1-rho) * beta.tilde
      s.tilde.sq <- rho * s.tilde.sq.new + (1-rho) * s.tilde.sq
    }
    
    
    ## Step 3d
    t.sq <- 1 / (1/s.tilde.sq - 1/t.tilde.sq)
    t.sq <- min(max(t.sq, 1e-11), 1e11)
    
    r.old <- r
    r <- (t.tilde.sq * beta.tilde - s.tilde.sq * r.tilde) / (t.tilde.sq - s.tilde.sq)
    
    
    ## Check for convergence
    Delta <- r-r.old
    Delta <- sum(Delta*Delta)/sum(r*r) < 1E-8
    
    if(is.na(Delta)) {
      if(rho <= min.rho) {
        warning("\nVAMP failed to converge: Stopping condition evaluated as 'NA'.\n")
        break
      }
      
      if(verbose) cat("\nVAMP failed to converge: Stopping condition evaluated as 'NA'. Setting rho to", rho/2, "\n")
      return(VAMP(y, X, lambda, psi, sigma.sq = sigma.sq0, rho/2, max.iter, min.rho, verbose, a.0, b.0)) # If VAMP fails to converge, add dampenin
    }
    
    # Stop iterating if VAMP has sufficiently converged
    if(Delta) {
      if(verbose) cat("\nVAMP has converged with rho equal to", rho)
      break
    }
    
    if(verbose) setTxtProgressBar(pb, k)
  }
  if(verbose) close(pb)
  
  if(k == max.iter) {
    if(rho > min.rho) {
      if(verbose) cat("\nVAMP failed to converge: Maximum number of iterations reached. Setting rho to", rho/2, "\n")
      return(VAMP(y, X, lambda, psi, sigma.sq = sigma.sq0, rho/2, max.iter, min.rho, verbose, a.0, b.0)) # If VAMP fails to converge, add dampening
    }
    
    warning("VAMP failed to converge: Maximum number of iterations reached")
  }
  
  return(list(mean=as.vector(beta), variance=s.sq, sigma.sq = sigma.sq))
}