# Setup the R environment:
source("setup.R")


# Number of repetitions
n.rep <- 20

# Number of predictors r
r <- 60

psi <- 1
sigma.sq <- .5
n <- 100
n_nonzero <- 40
lambda <- n_nonzero / r


sequential_PIP <- array(NA_real_, dim = c(n.rep, r))
random_PIP <- array(NA_real_, dim = c(n.rep, r))
belsley_PIP <- array(NA_real_, dim = c(n.rep, r))
spectral_PIP <- array(NA_real_, dim = c(n.rep, r))

Gibbs_PIP <- array(NA_real_, dim = c(n.rep, r))

# Set seed for reproducibility
set.seed(1)

pb <- txtProgressBar(max = n.rep, style = 3)
for (i in 1:n.rep) {
  
  A <- matrix(rnorm(n * r), nrow = n, ncol = r)
  A <- scale(A)
  
  theta <- rep(0, r)
  theta[sample.int(n = r, size = n_nonzero)] <- rnorm(n = n_nonzero, sd = sqrt(psi))
  y <- rnorm(n = n, mean = A %*% theta, sd = sqrt(sigma.sq))
  
  
  Gibbs_PIP[i, ] <- PIP_Gibbs(y, X = A, lambda, psi, sigma.sq, beta.0 = theta, iter = 1e5)
  
  sequential_PIP[i, ]  <- PIP_IRGA(y, A, lambda, psi, sigma.sq, theta_split = "sequential")
  
  tmp <- rep(0, r)
  for (j in 1:10) tmp <- tmp + PIP_IRGA(y, A, lambda, psi, sigma.sq, theta_split = "random")
  random_PIP[i, ]  <- tmp/10
  
  belsley_PIP[i, ]  <- PIP_IRGA(y, A, lambda, psi, sigma.sq, theta_split = "belsley")
  
  spectral_PIP[i, ]  <- PIP_IRGA(y, A, lambda, psi, sigma.sq, theta_split = "spectral")
  
  setTxtProgressBar(pb, i)
  
}
close(pb)



## Compute absolute errors
sequential_error <- abs(Gibbs_PIP - sequential_PIP)
random_error <- abs(Gibbs_PIP - random_PIP)
belsley_error <- abs(Gibbs_PIP - belsley_PIP)
spectral_error <- abs(Gibbs_PIP - spectral_PIP)

summary(as.vector(sequential_error))
summary(as.vector(random_error))
summary(as.vector(belsley_error))
summary(as.vector(spectral_error))