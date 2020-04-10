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



sequential_PIP <- array(NA_real_, dim = c(n.rep, r))
random_PIP <- array(NA_real_, dim = c(n.rep, r))
spectral_PIP <- array(NA_real_, dim = c(n.rep, r))


Gibbs_PIP <- array(NA_real_, dim = c(n.rep, r))

# Set seed for reproducibility
set.seed(1)

for (j in 1:n.rep) {
  cat("j:", j, "\n")
  
  theta <- rep(0, r)
  theta[sample.int(n = r, size = p)] <- rnorm(n = p, sd = sqrt(psi))
  y <- rnorm(n = n, mean = A%*%theta, sd = sqrt(sigma.sq))
  
  
  Gibbs_PIP[j, ] <- PIP_Gibbs(y, X = A, lambda, psi, sigma.sq, beta.0 = theta, iter = 1e5)
  
  sequential_PIP[j, ]  <- PIP_IRGA(y, A, lambda, psi, sigma.sq, theta_split = "sequential")
  
  tmp <- rep(0, r)
  for (i in 1:10) tmp <- tmp + PIP_IRGA(y, A, lambda, psi, sigma.sq, theta_split = "random")
  random_PIP[j, ]  <- tmp/10
  
  spectral_PIP[j, ]  <- PIP_IRGA(y, A, lambda, psi, sigma.sq, theta_split = "spectral")
  
}


# We cap the log odds since the approximate PIP can equal 1.
log_odds <- function(PIP) pmin(as.vector(log(PIP) - log(1-PIP)), 50)
Gibbs_log_odds <- log_odds(Gibbs_PIP)
sequential_log_odds <- log_odds(sequential_PIP)
random_log_odds <- log_odds(random_PIP)
spectral_log_odds <- log_odds(spectral_PIP)


summary(as.vector(abs(Gibbs_log_odds - sequential_log_odds)))
summary(as.vector(abs(Gibbs_log_odds - random_log_odds)))
summary(as.vector(abs(Gibbs_log_odds - spectral_log_odds)))