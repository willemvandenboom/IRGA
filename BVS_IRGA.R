# For intToBin()
if(!'R.utils' %in% rownames(installed.packages())) install.packages('R.utils', dependencies = T)

# For parallel computing
if(!'doParallel' %in% rownames(installed.packages())) install.packages('doParallel', dependencies = T)


## Exact computation of posterior inclusion probabilities (PIPs)
# with non-iid errors with covariance matrix Sigma
PIP_exact <- function(y, X, lambda, psi, Sigma) {
  log.sum.exp <- function(a, b) {
    # Computes log(exp(a)+exp(b)), vectorizing over b argument.
    # Uses offset trick to avoid numeric overflow.
    min.ab <- pmin(a, b)
    max.ab <- pmax(a, b)
    
    return(log1p(exp(min.ab-max.ab)) + max.ab)
  }
  
  X <- as.matrix(X)
  p <- ncol(X) # Number of candidate variables
  y <- as.vector(y)
  
  # Log-density evaluation for a zero-mean multivariate normal
  # This is a simplified and faster version of `dmvnorm` from the `mvtnorm` R package.
  # The normalization constant involving pi has been left out as
  # we only need the density up to a proportionality constant.
  ldmvnorm <- function(Sigma) {
    dec <- chol(Sigma)
    tmp <- backsolve(dec, y, transpose = TRUE)
    return(-sum(log(diag(dec))) - 0.5 * sum(tmp*tmp))
  }
  
  # The marginal log-likelihood of the model without predictors
  LL <- p*log(1-lambda) + ldmvnorm(Sigma)
  
  # Object to store the marginal log-likelihood that a beta_j is *not* included
  LLj <- rep(LL, p)
  
  for(k in 1:(2^p-1)) {
    gam <- integer(p)
    gam[min(p-ceiling(log2(k+1)-1),p):p] <- as.numeric(strsplit(R.utils::intToBin(k), split = "")[[1]])
    index <- which(gam==1)
    X.gam <- X[,index, drop=FALSE] # `drop=FALSE` such that X.gam does not reduce to a vector if s.gam=1.
    s.gam <- length(index)
    LL.gam <- s.gam*log(lambda) + (p-s.gam)*log(1-lambda) + ldmvnorm(Sigma + psi*tcrossprod(X.gam))
    
    LL <- log.sum.exp(a=LL.gam, b=LL)
    LLj[-index] <- log.sum.exp(a=LL.gam, b=LLj[-index])
  }
  
  return(-expm1(LLj - LL)) # 1 - exp(LLj)/exp(LL)
}




## Implementation of IRGA

BVS_IRGA_12 <- function(
  y, X, Z, lambda, psi, sigma.sq = NULL, max.iter = 100,
  min.rho = 1, rho = 1, estimator = "VAMP", verbose = FALSE
) {
  # Steps 1 and 2 of
  # Integrated rotated Gaussian approximation (IRGA, Algorithm 1)
  # using spike-and-slab priors as described in Section 3.1,
  # with a spike-and-slab prior on alpha as well:
  # Pi(alpha) = lambda * N(0, psi) + (1-lambda) delta(0).
  
  p <- NCOL(X) # Length of beta
  n <- length(y) # Number of observations
  
  
  ## Step 1
  # Compute the QR decomposition of X
  # to obtain the rotation matrix Q.
  Q <- qr.Q(qr(X), complete=TRUE)
  
  # Compute the rotated quantities
  QtY <- as.vector(crossprod(Q, y))
  RtY <- QtY[1:p]
  
  RtX <- crossprod(Q[,1:p, drop=FALSE], X)
  
  QtZ <- crossprod(Q, Z)
  RtZ <- QtZ[1:p,, drop=FALSE]
  
  # If p equals the number of observations,
  # then the nuisance parameter eta follows its prior
  # since StZ is empty then.
  if(p == n) {
    if(is.null(sigma.sq)) stop(
      'sigma.sq is unknown while IRGA is called with p = n such that sigma.sq cannot be inferred.'
    )
    
    # The prior variance of alpha is lambda * psi * I_q
    covariance <- lambda * psi * tcrossprod(RtZ)
    
    return(list(RtY=RtY, RtX=RtX, mean=numeric(p), covariance=covariance, sigma.sq=sigma.sq))
  }
  
  StY <- QtY[-(1:p)]
  StZ <- QtZ[-(1:p),, drop=FALSE]
  rm(y, QtY, QtZ)
  
  
  ## Step 2
  # Estimate the mean and variance of the nuisance parameter alpha
  # using vector approximate message passing (VAMP) or lasso.
  
  if(estimator == "VAMP") {
    
    # Load the function for vector approximate message passing (VAMP)
    source('VAMP.R')
    
    VAMPresult <- VAMP(
      y=StY, X=StZ, lambda, psi, sigma.sq,
      rho = rho, max.iter = max.iter,
      min.rho = min.rho, verbose = verbose
    )
    mean <- RtZ %*% VAMPresult$mean
    covariance <- VAMPresult$variance * tcrossprod(RtZ, RtZ)
    
    # If sigma.sq was unknown, it has been estimated by VAMP.
    sigma.sq <- VAMPresult$sigma.sq
    
  } else {
    
    # Use the debiased lasso (Javanmard & Montanari, 2013)
    q <- NCOL(StZ)
    
    # Keep the largest lasso estimate which has less than lambda*q+1 nonzero coefficients.
    lars.result <- coef(lars::lars(x=StZ, y=StY, type = "lasso", intercept = FALSE))
    alpha.lasso <- lars.result[max(which(rowSums(lars.result != 0) < lambda*q+1)),]
    
    lasso.u <- as.vector(alpha.lasso+t(StZ)%*%(StY-StZ%*%alpha.lasso)/(n-p))
    Sigma.hat <- t(StZ)%*%StZ/(n-p)
    mean <- as.vector((diag(q)-Sigma.hat)%*%alpha.lasso)
    
    a.0 <- 1
    b.0 <- 1
    sigma.sq <- ( b.0+sum((StY - StZ %*% alpha.lasso)^2)/2 ) / ( a.0 + (n-p)/2 )
    
    covariance <- sigma.sq*solve(Sigma.hat)/(n-p)
    
    mean <- RtZ %*% mean
    covariance <- tcrossprod(RtZ %*% covariance, RtZ)
    
  }
  
  return(list(RtY=RtY, RtX=RtX, mean=mean, covariance=covariance, sigma.sq=sigma.sq))
}


BVS_IRGA <- function(y, X, Z, lambda, psi, sigma.sq = NULL, max.iter = 100, estimator = "VAMP") {
  # Integrated rotated Gaussian approximation (IRGA, Algorithm 1)
  # using spike-and-slab priors as described in Section 3.1,
  # with a spike-and-slab prior on beta as well:
  # Pi(beta) = lambda * N(0, psi) + (1-lambda) delta(0) and
  # Pi(alpha) = lambda * N(0, psi) + (1-lambda) delta(0).
  # `max.iter` is the maximum number of VAMP iterations before adding dampening.
  # `estimator` specifies whether the mean and variance of (4)
  # are estimated using VAMP or lasso.
  
  p <- NCOL(X) # Length of beta
  
  ## Steps 1 and 2 of Algorithm 1:
  BVS_IRGA_12result <- BVS_IRGA_12(
    y, X, Z, lambda, psi, sigma.sq, max.iter, estimator = estimator
  )
  
  
  ## Step 3 of Algorithm 1:
  # Compute the posterior inclusion probabilities (PIPs).
  return(PIP_exact(
    y = BVS_IRGA_12result$RtY - BVS_IRGA_12result$mean,
    X = BVS_IRGA_12result$RtX,
    lambda, psi,
    # If sigma.sq was unknown, it has been estimated in Step 2.
    Sigma <- BVS_IRGA_12result$sigma.sq * diag(p) + BVS_IRGA_12result$covariance
  ))
}


# Computation of posterior inclusion probabilities (PIPs) using IRGA
# with iid erros with variance `sigma.sq`, potentially unknown.
PIP_IRGA <- function(
  y, A, lambda, psi, sigma.sq = NULL, p = NULL,
  max.iter = 100, estimator = "VAMP"
) {
  # `p` is the (maximum) length of beta in the executions of IRGA
  # `max.iter` is the maximum number of VAMP iterations before adding dampening.
  
  A <- as.matrix(A)
  r <- ncol(A) # Number of candidate variables
  y <- as.vector(y)
  
  # `n_splits` specifies the number of executions of IRGA
  if(is.null(p)) {
    
    # p = O(log(r)) yields a favorable computational complexity.
    p = floor(log(r))
    n_splits <- ceiling(r/p)
    
    # To make best use of compute power,
    # the number of splits should not be less than the number of CPU cores.
    n_splits <- max(n_cores, n_splits)
    
    # The number of splits cannot be so small that p > n,
    # because then IRGA does not apply.
    n_splits <- max(r %/% length(y), n_splits)
    
  } else {
    if(p > n) stop('p>n such that IRGA does not apply.')
    n_splits <- ceiling(r/p)
  }
  
  # Compute the start and end indices of the splits
  split_ind <- 0:n_splits * ceiling(r/n_splits)
  
  # If r %% n_splits != 0, we need to make the splits of unequal size.
  diff <- split_ind[n_splits+1] - r
  if(diff > 0) for(d in 1:diff) {
      split_ind[n_splits+2-d] <- split_ind[n_splits+2-d] - diff + d - 1
  }
  
  s <- 1
    beta_ind <- (split_ind[s]+1):split_ind[s+1]
  
  BVS_IRGA(
    y,
    X = A[, beta_ind],
    Z = A[, -beta_ind],
    lambda, psi, sigma.sq, max.iter,
    estimator = estimator
  )
  
  # Setup for parallel computing
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  
  # Compute the PIPs using IRGA and parallel computing.
  PIP <- foreach::`%dopar%`(
    foreach::foreach(
      s=1:n_splits, .combine='c',
      .export = c('BVS_IRGA', 'BVS_IRGA_12', 'PIP_exact')
    ),
    {
      beta_ind <- (split_ind[s]+1):split_ind[s+1]
      
      BVS_IRGA(
        y,
        X = A[, beta_ind],
        Z = A[, -beta_ind],
        lambda, psi, sigma.sq, max.iter,
        estimator = estimator
      )
    }
  )
  
  parallel::stopCluster(cl)
  
  return(PIP)
}