# For the variational Bayes method from Carbonetto & Stephens (2012, doi:10.1214/12-BA703)
# The version on CRAN contains a bug. We therefore obtain a more recent version from Github.
if(!'varbvs' %in% rownames(installed.packages()) || packageDescription("varbvs")$Version == "2.4-0"){
  if(!'remotes' %in% rownames(installed.packages())) install.packages('remotes', dependencies = T)
  remotes::install_github('pcarbo/varbvs', subdir = 'varbvs-R')
}


# Computation of posterior inclusion probabilities (PIPs)
# using the variational Bayes (VB) method from Carbonetto & Stephens (2012, doi:10.1214/12-BA703).
PIP_VB <- function(y, X, lambda, psi, sigma.sq = NULL) {
  # Determine whether sigma.sq is given.
  s2unknown <- is.null(sigma.sq)
  
  if(s2unknown) {
    sigma <- sd(y)
    
    # The prior variance in this variational Bayes algorithm
    # is given by psi/sigma.sq.
    # Since sigma.sq is unkown, we estimate psi/sigma.sq as follows.
    sa <- psi/sd(y)
  } else {
    sigma <- sigma.sq
    sa <- psi/sigma.sq
  }
  
  # `varbvs` automatically includes an intercept,
  # which is not in the model that IRGA fits.
  # However, any differences as a result should be small.
  return(as.vector(varbvs::varbvs(
    X, Z = NULL, y = as.vector(y),
    family = "gaussian",
    sigma = sigma, sa = sa,
    logodds = log10(lambda) - log10(1-lambda),
    update.sigma = s2unknown, update.sa = FALSE
  )$alpha))
}


# Computation of posterior inclusion probabilities (PIPs)
# using the variational Bayes (VB) method from
# Ormerod et al. (2012, Algorithm 1, doi:10.1214/17-EJS1332).
PIP_VB_Ormerod <- function(y, X,
  rho, # prior inclusion probability, same as `lambda`
  sigma.sq.beta = 1, # slab variance, same as `psi`
  sigma.sq = NULL
) {
  X <- as.matrix(X)
  # Number of candidate predictors
  p <- ncol(X)
  
  # Number of observations
  n <- length(y)
  
  # Determine whether sigma.sq is given.
  s2unknown <- is.null(sigma.sq)
  
  if(s2unknown) {
    tau <- 1000
    
    # Shape parameter of inverse-gamma prior on sigma.sq
    A <- a.0
    
    # Rate parameter of inverse-gamma prior on sigma.sq
    B <- b.0
  } else tau <- 1/sigma.sq
  
  
  # Initialize the posterior inclusion probabilities at their prior
  w <- rep(rho, p)
  
  # Precompute
  XtX <- crossprod(X)
  XtY <- crossprod(X, y)
  YtY <- as.vector(crossprod(y))
  
  # Line 3
  lambda <- log(rho)-log(1-rho) # logit(rho)
  
  for(t in 1:1e3) {
    # Line 5
    W <- diag(w)
    Omega <- tcrossprod(w) + diag(w * (1-w))
    
    # Line 6
    Sigma <- solve(
      tau * XtX * Omega + diag(p)/sigma.sq.beta
    )
    
    mu <- tau * Sigma %*% W %*% XtY
    
    if(s2unknown) {
      # Line 7
      # `as.vector` ensures `s` does no become a 1x1 matrix
      s <- B + .5 * as.vector((
        YtY - 2 * crossprod(XtY, W) %*% mu + sum(diag(
          (XtX * Omega) %*% (tcrossprod(mu) + Sigma)
        ))
      ))
      
      # Line 8
      # `as.vector` ensures `tau` does no become a 1x1 matrix
      tau <- as.vector( (A + n/2)/s )
    }
    
    # Line 9
    w.star <- w
    
    # Line 10
    for(j in 1:p) {
      
      # Line 11
      # `as.vector` ensures `eta.j` does no become a 1x1 matrix
      eta.j <- as.vector(
        lambda - .5 * tau * (mu[j]^2 + Sigma[j,j]) * crossprod(X[,j])
        + tau * crossprod(
          X[,j],
          y*mu[j] - X[,-j] %*% diag(w.star[-j]) %*% (mu[-j]*mu[j] + Sigma[-j,j])
        )
      )
      
      # Line 12
      w.star[j] <- 1/(1+exp(-eta.j)) # expit(eta.j)
    }
    
    # Line 13
    w <- w.star
    
    # Compute the log of the evidence lower bound
    # We drop any parts that do not change with t
    ELBO <- ifelse(s2unknown, -(A + n/2)*log(s), 0) + .5*log(det(Sigma)) - sum(diag(
      tcrossprod(mu) + Sigma
    )) / (2*sigma.sq.beta) + sum(
      ifelse(
        w == 0,
        0,
        w * (log(rho) - log(w))
      ) + ifelse(
        w == 1,
        0,
        (1-w) * (log(1-rho) - log(1-w))
      )
    )
    
    if(t > 1) {
      # Line 14
      if(is.na(abs(ELBO - ELBO.old) < 1e-12)) {
        warning("Stopping condition returned 'na'")
        break
      }
      if(abs(ELBO - ELBO.old) < 1e-12) {
        cat("Converged after", t, "iterations\n")
        return(w)
      }
    }
    
    ELBO.old <- ELBO
  }
  
  warning("VBVS did not converge.")
  return(w)
}



# Computation of posterior inclusion probabilities (PIPs)
# using the expectation propagation (EP) method from Hernandez-Lobato et al. (2015, doi:10.1007/s10994-014-5475-7).
PIP_EP <- function(Y, X,
  p0 = 0.5, # prior inclusion probability, same as `lambda`
  v = 1, # slab variance, same as `psi`
  sigma.sq = NULL) {
  
  # Load the EP code by Jose Miguel Hernandez Lobato
  source("epBVS.R")
  
  # Determine whether sigma.sq is given.
  s2unknown <- is.null(sigma.sq)
  
  if(!s2unknown) {
    ret <- epBVSinternal(X, Y, beta = 1/sigma.sq, p0, v)
    return(as.vector(1/(1+exp(-ret$phi))))
  }
  
  
  # If `sigma.sq` is not given, we need to learn it:
  
  # Initialize the noise precision beta
  beta <- 1
  
  # We find the optimal configuration for the hyper-parameters
  
  target <- function(param) {
    
    beta <- exp(param[ 1 ])
    #p0 <- logistic(param[ 2 ])
    #v <- exp(param[ 3 ])
    
    -epBVSinternal(X, Y, beta, p0, v)$evidence
  }
  
  # We call the optimization method
  
  startPoint <- c(log(beta)) #, invLogistic(p0), log(v))
  ret <- optim(startPoint, target, method = "Nelder-Mead", control = list(trace = T, maxit = 40))$par
  
  time <- system.time(ret <- epBVSinternal(X, Y, exp(ret[ 1 ]), p0, v))
  
  ret$timeInternal <- time[[ 1 ]] + time[[ 2 ]]
  
  return(as.vector(1/(1+exp(-ret$phi))))
}



# Computation of posterior inclusion probabilities (PIPs) using a Gibbs sampler
PIP_Gibbs <- function(y, X, lambda, psi, sigma.sq = NULL, beta.0 = NULL, iter = 1e4, MCerror = FALSE) {
  # `beta.0` is the initialization of the sampler.
  # `iter` is the number of MCMC iterations,
  # the first 10% of which are discarded as burnin iterations.
  # If `MCerror` is true, an estimate of the Monte Carlo standard
  # error is produced using overlapping batch means (Flegal & Jones, 2010, doi:10.1214/09-AOS735, Sec. 3)
  
  if(MCerror) gibbsStart <- Sys.time()
  
  X <- as.matrix(X)
  p <- ncol(X) # Number of predictors
  n <- length(y) # Number of observations
  
  # Matrix to store whether the sample beta is zero or not
  gamma.gibbs <- matrix(NA, nrow = iter, ncol = p)
  
  # Current value for beta
  if(is.null(beta.0)) {
    beta <- numeric(p)
  } else beta <- beta.0
  
  s2unknown <- is.null(sigma.sq) # Record whether sigma.sq is given.
  
  # Initialize sigma.sq
  if(s2unknown) sigma.sq <- b.0/a.0
  
  l.lam <- log(1-lambda)-log(lambda)
  column.norm.sq <- colSums(X*X)
  
  pb <- txtProgressBar(max = iter, style = 3)
  for(k in 1:iter) {
    
    for(j in 1:p) {
      Xj <- X[,j]
      
      # Residual vector, recomputed in full every 100th time to avoid compounding numerical error
      if(j %% 1e2 == 1) {
        r <- y-X[,-j]%*%beta[-j]
      } else r <- r+Xj*beta[j]
      
      inner.prod <- sum(r*Xj)
      
      if(
        gamma.gibbs[k,j] <- runif(1) < 1 / ( 1 + exp(
          l.lam + .5 * log1p(column.norm.sq[j] * psi/sigma.sq)
          - .5 / sigma.sq * inner.prod*inner.prod / (sigma.sq/psi + column.norm.sq[j])
        ))
      ) {
        # If beta is nonzero, resample it.
        post.var <- 1/(1/psi+column.norm.sq[j]/sigma.sq)
        beta[j] <- rnorm(1, mean = max(min(post.var*inner.prod/sigma.sq, 1e12), -1e12), sd = sqrt(post.var))
        r <- r-Xj*beta[j]
      } else beta[j] <- 0
    }
    
    # Update sigma.sq
    if(s2unknown) sigma.sq <- 1/rgamma(n = 1, shape = a.0+n/2, rate = b.0+sum(r*r))
    
    setTxtProgressBar(pb, k)
  }
  close(pb)
  
  if(MCerror) {
    cat("Gibbs computation time:")
    print(Sys.time() - gibbsStart)
    
    cat("Computing Monte Carlo standard error\n")
    
    # For `mcse`
    if(!'mcmcse' %in% rownames(installed.packages())) install.packages('mcmcse', dependencies = T)
    
    mSE <- 0
    pb <- txtProgressBar(max = p, style = 3)
    for(j in 1:p) {
      mSE <- mSE + mcmcse::mcse(x = gamma.gibbs[,j], method = "obm", warn = TRUE)$se
      setTxtProgressBar(pb, j)
    }
    close(pb)
    cat("Average Monte Carlo standard error:", mSE/p)
  }
  
  # The PIP estimates follow as the mean of gamma,
  # excluding the first 10% of iterations as burnin.
  return(colMeans(gamma.gibbs[-(1:(iter %/% 10)), ]))
}