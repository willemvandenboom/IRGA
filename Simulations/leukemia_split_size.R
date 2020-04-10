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

split.size <- c(1, 2, 4, 8, 16)
n.split.size <- length(split.size)

sequential_PIP <- array(NA_real_, dim = c(n.split.size, n.rep, r))


Gibbs_PIP <- array(NA_real_, dim = c(n.rep, r))

# Set seed for reproducibility
set.seed(1)

for (j in 1:n.rep) {
  cat("j:", j, "\n")
  
  theta <- rep(0, r)
  theta[sample.int(n = r, size = p)] <- rnorm(n = p, sd = sqrt(psi))
  y <- rnorm(n = n, mean = A%*%theta, sd = sqrt(sigma.sq))
  
  
  Gibbs_PIP[j, ] <- PIP_Gibbs(y, X = A, lambda, psi, sigma.sq, beta.0 = theta, iter = 1e5)
  
  for (k in 1:n.split.size) sequential_PIP[k, j, ]  <- PIP_IRGA(y, A, lambda, psi, sigma.sq, p = split.size[k])
  
}


# We cap the log odds since the approximate PIP can equal 1.
log_odds <- function(PIP) pmin(as.vector(log(PIP) - log(1-PIP)), 50)
Gibbs_log_odds <- log_odds(Gibbs_PIP)

for (k in 1:n.split.size) {
  print(split.size[k])
  print(summary(as.vector(abs(Gibbs_log_odds - log_odds(sequential_PIP[k, , ])))))
}