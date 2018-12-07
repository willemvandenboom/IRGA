# Setup the R environment:
source("setup.R")


# Number of repetitions
n.rep <- 10

# The leukemia data is available in the 'varbvs' package
library(varbvs)

data(leukemia)
A <- leukemia$x

n <- nrow(A)
r <- ncol(A)


psi <- 1
sigma.sq <- .5
p <- 40 # `p` is the number of nonzeros in `theta`
lambda <- p/r


IRGA_PIP <- array(NA_real_, dim = c(n.rep, r))
IRGA_time <- array(NA_real_, dim = c(n.rep))

VB_PIP <- array(NA_real_, dim = c(n.rep, r))
VB_time <- array(NA_real_, dim = c(n.rep))

EP_PIP <- array(NA_real_, dim = c(n.rep, r))
EP_time <- array(NA_real_, dim = c(n.rep))

Gibbs_PIP <- array(NA_real_, dim = c(n.rep, r))

# Set seed for reproducibility
set.seed(1)
for(j in 1:n.rep) {
  cat("j:", j, "\n")
  
  theta <- rep(0, r)
  theta[sample.int(n = r, size = p)] <- rnorm(n = p, sd = sqrt(psi))
  y <- rnorm(n = n, mean = A%*%theta, sd = sqrt(sigma.sq))
  
  
  Gibbs_PIP[j,] <- PIP_Gibbs(y, X = A, lambda, psi, sigma.sq, beta.0 = theta, iter = 1e5)
  
  
  ## IRGA
  irgaStart <- Sys.time()
  IRGA_PIP[j,]  <- PIP_IRGA(y, A, lambda, psi, sigma.sq)
  
  # Record computation time
  IRGA_time[j] <- as.numeric(Sys.time() - irgaStart, unit = 'secs')
  
  
  ## EP
  epStart <- Sys.time()
  EP_PIP[j,] <- PIP_EP(y, X = A, lambda, psi, sigma.sq)

  # Record computation time
  EP_time[j] <- as.numeric(Sys.time() - epStart, unit = 'secs')


  ## VB
  vbStart <- Sys.time()
  VB_PIP[j,] <- PIP_VB(y, X = A, lambda, psi, sigma.sq)

  # Record computation time
  VB_time[j] <- as.numeric(Sys.time() - vbStart, unit = 'secs')
}



IRGA_error <- as.vector(abs(Gibbs_PIP - IRGA_PIP))

cat("IRGA absolute difference in posterior inclusion probability:")
summary(IRGA_error)
summary(IRGA_time)


VB_error <- as.vector(abs(Gibbs_PIP - VB_PIP))

cat("VB absolute difference in posterior inclusion probability:")
summary(VB_error)
summary(VB_time)


EP_error <- as.vector(abs(Gibbs_PIP - EP_PIP))

cat("EP absolute difference in posterior inclusion probability:")
summary(EP_error)
summary(EP_time)