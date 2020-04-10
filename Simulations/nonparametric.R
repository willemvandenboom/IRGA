# Setup the R environment:
source("setup.R")


# For rmvnorm(), random number generation from a multivariate normal distribution
if(!'mvtnorm' %in% rownames(installed.packages())) install.packages('mvtnorm', dependencies = T)


set.seed(1) # For reproducibility

# Set up problem parameters and generate data
n <- 100
p <- 3
q <- 2
beta <- 2*c(2, -2, 2)
g <- function(f) f^2

# The derivative of g
grad_g <- Vectorize(function(f) 2*f)

phi <- .9

Phi <- toeplitz(phi^(0:(p-1)))
X <- mvtnorm::rmvnorm(n = n, sigma = Phi)

Z <- X[,1]


# Compute the prior covariance `Sigma` of F
# implied by the Gaussian process prior on f
Sigma <- matrix(NA_real_, nrow = n, ncol = n)
for(i in 1:n) for(j in 1:n) Sigma[i,j] <- .1*sum((Z[i]-Z[j])^2)
Sigma <- exp(-Sigma)
Sigma_inv <- solve(Sigma + 1e-13*diag(n))

# The prior mean `mu` of F equals zero.
mu <- numeric(n)

# Set f equal to a draw from the prior
F <- as.vector(mvtnorm::rmvnorm(n = 1, mean = mu, sigma = Sigma))

# Visualize f
ord.Z <- order(Z)
plot(x = Z[ord.Z], y = F[ord.Z], type = 'l')

# Error variance
sigma.sq <- 1

y <- rnorm(n = n, mean = X%*%beta + g(F), sd = sqrt(sigma.sq))


# Prior variance for beta
psi <- 16



### The IRGA Algorithm

irgaStart <- Sys.time()

## Step 1
# Compute the QR decomposition of X
# to obtain the rotation matrix Q.
Q <- qr.Q(qr(X), complete=TRUE)

# Compute the rotated quantities
QtY <- as.vector(crossprod(Q, y))
MtY <- QtY[1:p]

MtX <- crossprod(Q[,1:p, drop=FALSE], X)

M <- Q[, 1:p]
S <- Q[, -(1:p)]
StY <- QtY[-(1:p)]



## Step 2
# Estimate the mean and variance of the nuisance parameter alpha
# using a Gauss-Newton algorithm derived from iteratively linearizing `g`.

# Initialize the estimate `m` of E[F | StY]
m <- as.vector(mvtnorm::rmvnorm(n = 1, mean = mu, sigma = Sigma))

SSt <- tcrossprod(S)
Sigma_inv_mu <- Sigma_inv %*% mu

for(t in 1:1e3) {
  m.old <- m
  
  # Set the learning rate `rho`
  rho <- 1/(1 + t^.55)
  
  # Compute the Jacobian at the current value of `m`
  J_m <- diag(grad_g(m))
  
  # Compute the estimate of Cov[F | StY] using the current linearization.
  Sigma.star <- solve( crossprod(J_m, SSt) %*% J_m + Sigma_inv )
  m.new <- as.vector( Sigma.star %*% (
    1/sigma.sq * crossprod(J_m, SSt) %*% (y - g(m) + J_m %*% m)
    + Sigma_inv_mu
  ) )
  
  m <- (1 - rho) * m.old + rho * m.new
  
  if(sum((m.new - m.old)^2) < 1e-8) break
}
cat('Stopped iterative scheme after', t, 'iterations.\n')


# Estimate the mean and covariance for the Gaussian approximation via simple Monte Carlo
# from the approximate P(eta | StY) = P(G(F) | StY).
eta.samples <- g(mvtnorm::rmvnorm(n = 1e3, mean = m, sigma = (Sigma.star + t(Sigma.star))/2))
Mt_eta.samples <- eta.samples %*% M


Gauss_mean <- colMeans(Mt_eta.samples)
Gauss_covariance <- cov(Mt_eta.samples)


## Step 3
# Compute the posterior mean and covariance for beta
tmp <- solve(Gauss_covariance + sigma.sq*diag(p))
IRGA_cov <- solve(crossprod(MtX, tmp%*%MtX) + diag(p)/psi)
IRGA_mean <- IRGA_cov%*%crossprod(MtX, tmp%*%(MtY - Gauss_mean))


irgaEnd <- Sys.time()






### Inference using MCMC

MCMCstart <- Sys.time()

# Loglikelihood evaluation for the samples of f, with beta integrated out
ll.precision <- solve(sigma.sq*diag(n) + psi * tcrossprod(X))

ll.F <- function(FF) {
  tmp <- g(as.vector(FF)) - y
  # We leave out the normalization constant involving pi as
  # we only need the density up to a proportionality constant.
  return(-.5 * crossprod(tmp, ll.precision) %*% tmp)
}

# Number of Monte Carlo samples
n.samp <- 1e5

# Matrix to store the posterior sample of F
F.samp <- matrix(NA_real_, nrow = n.samp, ncol = n)

# Initialize the sampler at the f with which y was generated
F.current <- numeric(n)#f(z1 = Z[,1], z2 = Z[,2])

# Compute the current loglikelihood
ll.current <- ll.F(F.current)

# Set the radom walk step size
step_size <- 1

# Count the number of accepts
n.accept <- 0

pb <- txtProgressBar(max = n.samp, style = 3)
for(s in 1:n.samp) {
  # Generate a proposal that is reversible with respect to the prior
  F.proposal <- sqrt(1 - step_size^2) * F.current + step_size * mvtnorm::rmvnorm(n = 1, mean = mu, sigma = Sigma)
  
  # Compute the loglikelihood for the proposal
  ll.proposal <- ll.F(F.proposal)
  
  if(log(runif(1)) < ll.proposal - ll.current) {
    F.current <- F.proposal
    n.accept <- n.accept+1
  }
  
  F.samp[s,] <- F.current
  
  setTxtProgressBar(pb, s)
}
close(pb)

MCMCend <- Sys.time()

cat("The acceptance rate is", n.accept/n.samp, "\n")


# Plot samples of F to get a sense of MCMC convergence

plot(x = Z[ord.Z], y = g(F[ord.Z]), type = 'l')

plot(x = Z[ord.Z], y = g(F.samp[n.samp %/% 100,ord.Z]), type = 'l', ylim = c(0, 8))
for(k in 2:100) lines(x = Z[ord.Z], y = g(F.samp[k * n.samp %/% 100,ord.Z]))




## Plot the results
# Plot marginal posteriors for beta
grid.size <- 1e3
x.range <- 7
grid <- seq(from = -x.range, to = x.range, length.out = grid.size)


# Compute the marginal densities from IRGA
IRGA_dens <- matrix(NA_real_, nrow = 3, ncol = grid.size)
for(j in 1:p) IRGA_dens[j,] <- dnorm(x = grid, mean = IRGA_mean[j], sd = sqrt(IRGA_cov[j,j]))


# Compute the marginal densities from MCMC
MCMC_dens <- matrix(0, nrow = 3, ncol = grid.size)

MCMC_cov <- solve(crossprod(X)/sigma.sq + diag(p)/psi)

ind <- (n.samp %/% 10 + 1):n.samp

pb <- txtProgressBar(min = min(ind), max = n.samp, style = 3)
s.ind <- unique(round(ind, -2))
for(s in s.ind) {
  temp.mean <- MCMC_cov %*% crossprod(X, y - g(F.samp[s,])) / sigma.sq
  
  for(j in 1:p) MCMC_dens[j,] <- MCMC_dens[j,] + dnorm(x = grid, mean = temp.mean[j], sd = sqrt(MCMC_cov[j,j]))
  
  setTxtProgressBar(pb, s)
}
close(pb)

MCMC_dens <- MCMC_dens/length(s.ind)


# Compute the marginal densities resulting from setting eta = 0
zero_dens <- matrix(NA_real_, nrow = 3, ncol = grid.size)
temp.mean <- MCMC_cov %*% crossprod(X, y) / sigma.sq
for(j in 1:p) zero_dens[j,] <- dnorm(x = grid, mean = temp.mean[j], sd = sqrt(MCMC_cov[j,j]))


x.min <- c(3, -6, 3.5)
x.max <- c(5, -3, 5.5)

ylabs <- c(
  expression(paste(pi, "(", beta[1], " | y) and ", hat(pi), "(", beta[1], " | y)")),
  expression(paste(pi, "(", beta[2], " | y) and ", hat(pi), "(", beta[2], " | y)")),
  expression(paste(pi, "(", beta[3], " | y) and ", hat(pi), "(", beta[3], " | y)"))
)

xlabs <- c(
  expression(beta[1]),
  expression(beta[2]),
  expression(beta[3])
)

pdf("simulation_nonparametric.pdf", width = 6, height = 2)
par(mfrow = c(1,3), mar = c(4, 5, 2, 1)+.1)
for(j in 1:p) {
  plot(x = grid, y = MCMC_dens[j,], type = 'l', xlim = c(x.min[j], x.max[j]), ylab = ylabs[j], xlab = xlabs[j], ylim = c(0, 1.15*max(MCMC_dens[j,], IRGA_dens[j,], zero_dens[j,])))
  lines(x = grid, y = IRGA_dens[j,], lty = "11", lwd = 3)
  lines(x = grid, y = zero_dens[j,], lty = 3, lwd = 1)
}
dev.off()


print(irgaEnd - irgaStart)
print(MCMCend - MCMCstart)