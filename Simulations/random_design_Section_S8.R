# Setup the R environment:
source("setup.R")


# Number of repetitions
n.rep <- 20

# Number of different r considered
n.r <- 5
r.list <- 60* 2^(0:(n.r-1))
r.max <- max(r.list)

psi <- 1
sigma.sq <- .5
n <- 100
n_nonzero <- 40


IRGA_PIP <- array(NA_real_, dim = c(n.r, n.rep, r.max))
IRGA_time <- array(NA_real_, dim = c(n.r, n.rep))

VB_PIP <- array(NA_real_, dim = c(n.r, n.rep, r.max))
VB_time <- array(NA_real_, dim = c(n.r, n.rep))

EP_PIP <- array(NA_real_, dim = c(n.r, n.rep, r.max))
EP_time <- array(NA_real_, dim = c(n.r, n.rep))


Gibbs_PIP <- array(NA_real_, dim = c(n.r, n.rep, r.max))

# Set seed for reproducibility
set.seed(1)
for(i in 1:n.r) for(j in 1:n.rep) {
  cat("\ni:", i, "j:", j, "\n")
  
  r <- r.list[i]
  
  lambda <- n_nonzero/r
  
  A <- matrix(rnorm(n*r), nrow = n, ncol = r)
  A <- scale(A)
  
  theta <- rep(0, r)
  theta[sample.int(n = r, size = n_nonzero)] <- rnorm(n = n_nonzero, sd = sqrt(psi))
  y <- rnorm(n = n, mean = A%*%theta, sd = sqrt(sigma.sq))
  
  
  Gibbs_PIP[i,j,1:r] <- PIP_Gibbs(y, X = A, lambda, psi, sigma.sq, beta.0 = theta, iter = 1e5)
  
  
  ## IRGA
  irgaStart <- Sys.time()
  IRGA_PIP[i,j,1:r]  <- PIP_IRGA(y, A, lambda, psi, sigma.sq)
  
  # Record computation time
  IRGA_time[i,j] <- as.numeric(Sys.time() - irgaStart, unit = 'secs')
  
  
  ## EP
  epStart <- Sys.time()
  EP_PIP[i,j,1:r] <- PIP_EP(y, X = A, lambda, psi, sigma.sq)
  
  # Record computation time
  EP_time[i,j] <- as.numeric(Sys.time() - epStart, unit = 'secs')
  
  
  ## VB
  vbStart <- Sys.time()
  VB_PIP[i,j,1:r] <- PIP_VB(y, X = A, lambda, psi, sigma.sq)
  
  # Record computation time
  VB_time[i,j] <- as.numeric(Sys.time() - vbStart, unit = 'secs')
  
}



## Compute absolute errors
IRGA_error <- abs(Gibbs_PIP - IRGA_PIP)
EP_error <- abs(Gibbs_PIP - EP_PIP)
VB_error <- abs(Gibbs_PIP - VB_PIP)



## Plot the results

# Compute summary statistics of errors
IRGA_error_median <- apply(IRGA_error, 1, function(x) log(median(x, na.rm = TRUE)))
IRGA_error_min <- apply(IRGA_error, 1, function(x) log(quantile(x, prob = .25, na.rm = TRUE)))
IRGA_error_max <- apply(IRGA_error, 1, function(x) log(quantile(x, prob = .75, na.rm = TRUE)))

EP_error_median <- apply(EP_error, 1, function(x) log(median(x, na.rm = TRUE)))
EP_error_min <- apply(EP_error, 1, function(x) log(quantile(x, prob = .25, na.rm = TRUE)))
EP_error_max <- apply(EP_error, 1, function(x) log(quantile(x, prob = .75, na.rm = TRUE)))

VB_error_median <- apply(VB_error, 1, function(x) log(median(x, na.rm = TRUE)))
VB_error_min <- apply(VB_error, 1, function(x) log(quantile(x, prob = .25, na.rm = TRUE)))
VB_error_max <- apply(VB_error, 1, function(x) log(quantile(x, prob = .75, na.rm = TRUE)))

VB_error_min[!is.finite(VB_error_min)] <- -50
EP_error_min[!is.finite(EP_error_min)] <- -50

ymax = 1.05*max(IRGA_error_max, VB_error_max, EP_error_max)
mid.symbol <- 20
q.symbol <- 4


pdf("simulation_random_design.pdf", width = 5.5, height = 2)
par(mfrow = c(1,1), bty = 'l',
    cex = .5)

plot(IRGA_error_median, ylim = c(-10, ymax), xlim = c(1, n.r+.2), type = 'p', pch = mid.symbol,
     ylab = "Log of the absolute error in PIP",
     xlab = "r",
     xaxt = "n"
)
axis(side = 1, at = 1:n.r, labels = r.list)
points(IRGA_error_min, pch = q.symbol)
points(IRGA_error_max, pch = q.symbol)
for(i in 1:n.r) lines(x = rep(i, 2), y = c(IRGA_error_min[i], IRGA_error_max[i]))

points(x = (1:n.r)+.1, y = EP_error_median, col = 'blue', pch = mid.symbol)
points(x = (1:n.r)+.1, y = EP_error_min, col = 'blue', pch = q.symbol)
points(x = (1:n.r)+.1, y = EP_error_max, col = 'blue', pch = q.symbol)
for(i in 1:n.r) lines(x = rep(i, 2)+.1, y = c(EP_error_min[i], EP_error_max[i]), col = 'blue')

points(x = (1:n.r)+.2, y = VB_error_median, col = 'red', pch = mid.symbol)
points(x = (1:n.r)+.2, y = VB_error_min, col = 'red', pch = q.symbol)
points(x = (1:n.r)+.2, y = VB_error_max, col = 'red', pch = q.symbol)
for(i in 1:n.r) lines(x = rep(i, 2)+.2, y = c(VB_error_min[i], VB_error_max[i]), col = 'red')

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


pdf("simulation_random_design_time.pdf", width = 5.5, height = 2)
par(mfrow = c(1,1), bty = 'l',
    cex = .5)

plot(IRGA_time_median, ylim = c(ymin, ymax), type = 'p', pch = mid.symbol,
     ylab = "Log of computation time in seconds",
     xlab = "q",
     xaxt = "n"
)
axis(side = 1, at = 1:n.r, labels = r.list)
points(IRGA_time_min, pch = q.symbol)
points(IRGA_time_max, pch = q.symbol)
for(i in 1:n.r) lines(x = rep(i, 2), y = c(IRGA_time_min[i], IRGA_time_max[i]))

points(EP_time_median, col = 'blue', pch = mid.symbol)
points(EP_time_min, col = 'blue', pch = q.symbol)
points(EP_time_max, col = 'blue', pch = q.symbol)
for(i in 1:n.r) lines(x = rep(i, 2), y = c(EP_time_min[i], EP_time_max[i]), col = 'blue')

points(VB_time_median, col = 'red', pch = mid.symbol)
points(VB_time_min, col = 'red', pch = q.symbol)
points(VB_time_max, col = 'red', pch = q.symbol)
for(i in 1:n.r) lines(x = rep(i, 2), y = c(VB_time_min[i], VB_time_max[i]), col = 'red')

dev.off()