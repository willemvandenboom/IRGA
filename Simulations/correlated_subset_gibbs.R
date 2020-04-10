# Setup the R environment:
source("setup.R")


# Number of different q
n.q <- 7
q.list <- 15* 2^(0:(n.q-1))

# Number of repetitions at for each q
n.rep <- 20

# Sample size
n <- 100


p <- 2

psi <- 1
sigma.sq <- .5


IRGA_error <- array(NA_real_, dim = c(n.q, n.rep, 2))
IRGA_time <- array(NA_real_, dim = c(n.q, n.rep))

EP_error <- array(NA_real_, dim = c(n.q, n.rep, 2))
EP_time <- array(NA_real_, dim = c(n.q, n.rep))

VB_error <- array(NA_real_, dim = c(n.q, n.rep, 2))
VB_time <- array(NA_real_, dim = c(n.q, n.rep))

Ormerod_error <- array(NA_real_, dim = c(n.q, n.rep, 2))

Gibbs_PIP <- array(NA_real_, dim = c(n.q, n.rep, 2))

set.seed(1)
for(i in 1:n.q) for(j in 1:n.rep) {
  cat("\ni:", i, "j:", j, "\n")
  
  q <- q.list[i]
  
  r <- p+q
  lambda <- p/r
  
  theta <- numeric(r)
  theta[1:2] <- 1:2
  
  A <- matrix(rnorm(n*r), nrow = n, ncol = r)
  A[,2] <- .01*A[,2]+.99*A[,1]
  A <- scale(A)
  
  y <- rnorm(n = n, mean = A%*%theta, sd = sqrt(sigma.sq))
  
  PIP <- PIP_Gibbs(y, X = A, lambda, psi, sigma.sq, beta.0 = theta, iter = 1e5)[1:2]
  Gibbs_PIP[i,j,] <- PIP
  
  
  ## IRGA
  irgaStart <- Sys.time()
  IRGA_PIP <- BVS_IRGA(y, X = A[,1:2], Z = A[,-(1:2)], lambda, psi, sigma.sq)
  
  # Record computation time
  IRGA_time[i,j] <- as.numeric(Sys.time() - irgaStart, unit = 'secs')
  
  # Compute absolute difference
  IRGA_error[i,j,] <- abs(PIP - IRGA_PIP)
  
  
  ## EP
  epStart <- Sys.time()
  EP_PIP <- PIP_EP(y, X = A, lambda, psi, sigma.sq)[1:2]
  
  # Record computation time
  EP_time[i,j] <- as.numeric(Sys.time() - epStart, unit = 'secs')
  
  # Compute absolute difference
  EP_error[i,j,] <- abs(PIP - EP_PIP)
  
  
  ## VB
  vbStart <- Sys.time()
  VB_PIP <- PIP_VB(y, X = A, lambda, psi, sigma.sq)[1:2]
  
  # Record computation time
  VB_time[i,j] <- as.numeric(Sys.time() - vbStart, unit = 'secs')
  
  # Compute absolute difference
  VB_error[i,j,] <- abs(PIP - VB_PIP)
  
  
  ## Ormerod's VB method
  Ormerod_PIP <- PIP_VB_Ormerod(y, X = A, lambda, psi, sigma.sq)[1:2]
  
  # Compute absolute difference
  Ormerod_error[i,j,] <- abs(PIP - Ormerod_PIP)
}



## Plot the results

# Compute summary statistics of errors
IRGA_error_median <- apply(IRGA_error, 1, function(x) log(median(x)))
IRGA_error_min <- apply(IRGA_error, 1, function(x) log(quantile(x, prob = .25)))
IRGA_error_max <- apply(IRGA_error, 1, function(x) log(quantile(x, prob = .75)))

EP_error_median <- apply(EP_error, 1, function(x) log(median(x)))
EP_error_min <- apply(EP_error, 1, function(x) log(quantile(x, prob = .25)))
EP_error_max <- apply(EP_error, 1, function(x) log(quantile(x, prob = .75)))

VB_error_median <- apply(VB_error, 1, function(x) log(median(x)))
VB_error_min <- apply(VB_error, 1, function(x) log(quantile(x, prob = .25)))
VB_error_max <- apply(VB_error, 1, function(x) log(quantile(x, prob = .75)))

Ormerod_error_median <- apply(Ormerod_error, 1, function(x) log(median(x)))
Ormerod_error_min <- apply(Ormerod_error, 1, function(x) log(quantile(x, prob = .25)))
Ormerod_error_max <- apply(Ormerod_error, 1, function(x) log(quantile(x, prob = .75)))


ymax = 1.05*max(IRGA_error_max, VB_error_max, EP_error_max, Ormerod_error_max)
mid.symbol <- 20
q.symbol <- 4


pdf("simulation_correlated_subset_gibbs.pdf", width = 5.5, height = 2)
par(mfrow = c(1,1), bty = 'l',
    cex = .5)

plot(IRGA_error_median, ylim = c(min(IRGA_error_min), ymax), type = 'p', pch = mid.symbol,
     ylab = "Log of the absolute error in PIP",
     xlab = "q",
     xaxt = "n"
)
axis(side = 1, at = 1:n.q, labels = q.list)
points(IRGA_error_min, pch = q.symbol)
points(IRGA_error_max, pch = q.symbol)
for(i in 1:n.q) lines(x = rep(i, 2), y = c(IRGA_error_min[i], IRGA_error_max[i]))

points(x = (1:n.q)+.1, y = EP_error_median, col = 'blue', pch = mid.symbol)
points(x = (1:n.q)+.1, y = EP_error_min, col = 'blue', pch = q.symbol)
points(x = (1:n.q)+.1, y = EP_error_max, col = 'blue', pch = q.symbol)
for(i in 1:n.q) lines(x = rep(i, 2)+.1, y = c(EP_error_min[i], EP_error_max[i]), col = 'blue')

points(x = (1:n.q)+.2, y = VB_error_median, col = 'red', pch = mid.symbol)
points(x = (1:n.q)+.2, y = VB_error_min, col = 'red', pch = q.symbol)
points(x = (1:n.q)+.2, y = VB_error_max, col = 'red', pch = q.symbol)
for(i in 1:n.q) lines(x = rep(i, 2)+.2, y = c(VB_error_min[i], VB_error_max[i]), col = 'red')

# points(Ormerod_error_median, col = 'gray', pch = mid.symbol)
# points(Ormerod_error_min, col = 'gray', pch = q.symbol)
# points(Ormerod_error_max, col = 'gray', pch = q.symbol)

dev.off()



# Compute summary statistics of time
IRGA_time_median <- apply(IRGA_time, 1, function(x) log(median(x)))
IRGA_time_min <- apply(IRGA_time, 1, function(x) log(quantile(x, prob = .25)))
IRGA_time_max <- apply(IRGA_time, 1, function(x) log(quantile(x, prob = .75)))

EP_time_median <- apply(EP_time, 1, function(x) log(median(x)))
EP_time_min <- apply(EP_time, 1, function(x) log(quantile(x, prob = .25)))
EP_time_max <- apply(EP_time, 1, function(x) log(quantile(x, prob = .75)))

VB_time_median <- apply(VB_time, 1, function(x) log(median(x)))
VB_time_min <- apply(VB_time, 1, function(x) log(quantile(x, prob = .25)))
VB_time_max <- apply(VB_time, 1, function(x) log(quantile(x, prob = .75)))


ymax <- 1.05*max(IRGA_time_max, VB_time_max, EP_time_max)
ymin <- min(IRGA_time_min, VB_time_min, EP_time_min)


pdf("simulation_correlated_subset_gibbs_time.pdf", width = 5.5, height = 2)
par(mfrow = c(1,1), bty = 'l',
    cex = .5)

plot(IRGA_time_median, ylim = c(ymin, ymax), type = 'p', pch = mid.symbol,
     ylab = "Log of computation time in seconds",
     xlab = "q",
     xaxt = "n"
)
axis(side = 1, at = 1:n.q, labels = q.list)
points(IRGA_time_min, pch = q.symbol)
points(IRGA_time_max, pch = q.symbol)
for(i in 1:n.q) lines(x = rep(i, 2), y = c(IRGA_time_min[i], IRGA_time_max[i]))

points(EP_time_median, col = 'blue', pch = mid.symbol)
points(EP_time_min, col = 'blue', pch = q.symbol)
points(EP_time_max, col = 'blue', pch = q.symbol)
for(i in 1:n.q) lines(x = rep(i, 2), y = c(EP_time_min[i], EP_time_max[i]), col = 'blue')

points(VB_time_median, col = 'red', pch = mid.symbol)
points(VB_time_min, col = 'red', pch = q.symbol)
points(VB_time_max, col = 'red', pch = q.symbol)
for(i in 1:n.q) lines(x = rep(i, 2), y = c(VB_time_min[i], VB_time_max[i]), col = 'red')

dev.off()