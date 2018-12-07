# The leukemia data is available in the 'varbvs' package
library(varbvs)

data(leukemia)
A <- leukemia$x

n <- nrow(A)
r <- ncol(A)

predictor.selected <- 647
X <- A[, predictor.selected]
Z <- A[, -predictor.selected]


## Step 1 of Algorithm 1
# Compute the QR decomposition of X
# to obtain the rotation matrix Q.
Q <- qr.Q(qr(X), complete=TRUE)

# Compute the rotated quantities
QtZ <- crossprod(Q, Z)
p <- 1
RtZ <- QtZ[1:p,, drop=FALSE]
StZ <- QtZ[-(1:p),, drop=FALSE]


# Specify the spike-and-slab prior
psi <- 1
lambda <- 40/r



# Set seed for reproducibility
set.seed(1)

# Generate theta according to the prior
theta <- rep(0, r)
theta[sample.int(n = r, size = p)] <- rnorm(n = p, sd = sqrt(psi))

# Number of Gibbs samples used
n.iter <- 1e5


# Drawing posterior samples of beta using a Gibbs sampler,
# where beta is the parameter vector in the standard linear model y ~ N(X beta, sigma.sq I)
beta_Gibbs <- function(y, sigma.sq, X = StZ, iter = n.iter) {
  # `beta.0` is the initialization of the sampler.
  # `iter` is the number of MCMC iterations,
  # the first 10% of which are discarded as burnin iterations.
  
  X <- as.matrix(X)
  p <- ncol(X) # Number of predictors
  n <- length(y) # Number of observations
  
  # Matrix to store the samples of beta
  beta.gibbs <- matrix(NA_real_, nrow = iter, ncol = p)
  
  # Initialize beta the truth
  beta <- theta[-predictor.selected]
  
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
        runif(1) < 1 / ( 1 + exp(
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
    
    # Record the current beta sample
    beta.gibbs[k,] <- beta
    
    setTxtProgressBar(pb, k)
  }
  close(pb)
  
  
  # Return the posterior samples of beta,
  # excluding the first 10% as burnin iterations.
  return(beta.gibbs[-(1:(iter %/% 10)), ])
}



# Set sigma.sq equal to a half
sigma.sq <- .5

# Generate data
y <- rnorm(n = n, mean = A%*%theta, sd = sqrt(sigma.sq))

# Compute the rotated quantities
QtY <- as.vector(crossprod(Q, y))
StY <- QtY[-(1:p)]

RtZ.alpha.half <- as.vector(tcrossprod(RtZ, beta_Gibbs(y = StY, sigma.sq = sigma.sq)))


# Repeat with sigma.sq equal to 1:
sigma.sq <- 1

# Generate data
y <- rnorm(n = n, mean = A%*%theta, sd = sqrt(sigma.sq))

# Compute the rotated quantities
QtY <- as.vector(crossprod(Q, y))
StY <- QtY[-(1:p)]

RtZ.alpha.unit <- as.vector(tcrossprod(RtZ, beta_Gibbs(y = StY, sigma.sq = sigma.sq)))


pl.samp <- function(x, main) {
  qqnorm(x, pch = 16, main = main)
  qqline(x)
}

pdf("fig_1_check_normal.pdf", width = 6, height = 2, pointsize = 11)
par(mfcol = c(1, 2), cex = .5)
pl.samp(RtZ.alpha.half[(1:900) * (n.iter %/% 1e3)], main = expression(sigma^2*" = 1/2"))
pl.samp(RtZ.alpha.unit[(1:900) * (n.iter %/% 1e3)], main = expression(sigma^2*" = 1"))
dev.off()