# Setup the R environment:
source("setup.R")

# Set seed for reproduciblity
set.seed(1)

# The diabetes data is available in the 'lars' package
if(!'lars' %in% rownames(installed.packages())) install.packages('lars', dependencies = F)
library(lars)
data(diabetes)


A <- diabetes$x2

n <- nrow(A)
r <- ncol(A)

# Standardize the outcome
y <- diabetes$y-mean(diabetes$y)
y <- y/sd(y)

psi <- 1
lambda <- .5


# Set prior for sigma.sq
a.0 <- 1
b.0 <- 1



irgaStart <- Sys.time()
IRGA_PIP <- PIP_IRGA(y, A, lambda, psi, sigma.sq = NULL, p = NULL)
irgaEnd <- Sys.time()

irga_lassoStart <- Sys.time()
IRGA_lasso_PIP <- PIP_IRGA(y, A, lambda, psi, sigma.sq = NULL, p = NULL, estimator = "lasso")
irga_lassoEnd <- Sys.time()

epStart <- Sys.time()
EP_PIP <- PIP_EP(Y = y, X = A, lambda, psi, sigma.sq = NULL)
epEnd <- Sys.time()

vbStart <- Sys.time()
VB_PIP <- PIP_VB(y, X = A, lambda, psi, sigma.sq = NULL)
vbEnd <- Sys.time()

Ormerod_PIP <- PIP_VB_Ormerod(y, X = A, lambda, psi, sigma.sq = NULL)

Gibbs_PIP <- PIP_Gibbs(y, X = A, lambda, psi, sigma.sq = NULL, iter = 1e5, MCerror = TRUE)


log_odds <- function(PIP) log(PIP) - log(1-PIP)
Gibbs_log_odds <- log_odds(Gibbs_PIP)
IRGA_log_odds <- log_odds(IRGA_PIP)
IRGA_lasso_log_odds <- log_odds(IRGA_lasso_PIP)
EP_log_odds <- log_odds(EP_PIP)
VB_log_odds <- log_odds(VB_PIP)
Ormerod_log_odds <- log_odds(Ormerod_PIP)


cat("IRGA absolute difference in posterior log odds of inclusion:")
summary(abs(Gibbs_log_odds - IRGA_log_odds))
print(irgaEnd - irgaStart)

cat("IRGA with lasso absolute difference in posterior log odds of inclusion:")
summary(abs(Gibbs_log_odds - IRGA_lasso_log_odds))
print(irga_lassoEnd - irga_lassoStart)

cat("EP absolute difference in posterior log odds of inclusion:")
summary(abs(Gibbs_log_odds - EP_log_odds))
print(epEnd - epStart)

cat("VB absolute difference in posterior log odds of inclusion:")
summary(abs(Gibbs_log_odds - VB_log_odds))
print(vbEnd - vbStart)

summary(abs(Gibbs_log_odds - Ormerod_log_odds))



# Comparison of different methods to split theta
sequential_PIP  <- IRGA_PIP

tmp <- rep(0, r)
for (j in 1:10) tmp <- tmp + PIP_IRGA(y, A, lambda, psi, sigma.sq = NULL, theta_split = "random")
random_PIP  <- tmp/10

belsley_PIP  <- PIP_IRGA(y, A, lambda, psi, sigma.sq = NULL, theta_split = "belsley")

spectral_PIP  <- PIP_IRGA(y, A, lambda, psi, sigma.sq = NULL, theta_split = "spectral")


sequential_log_odds <- log_odds(sequential_PIP)
random_log_odds <- log_odds(random_PIP)
belsley_log_odds <- log_odds(belsley_PIP)
spectral_log_odds <- log_odds(spectral_PIP)


summary(as.vector(abs(Gibbs_log_odds - sequential_log_odds)))
summary(as.vector(abs(Gibbs_log_odds - random_log_odds)))
summary(as.vector(abs(Gibbs_log_odds - belsley_log_odds)))
summary(as.vector(abs(Gibbs_log_odds - spectral_log_odds)))