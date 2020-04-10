# Setup the R environment:
source("setup.R")


## Loading the data

demographics <- read.delim("ftp://ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/technical/working/20130606_sample_info/20130606_g1k.ped", colClasses = c(
  "NULL", "character", "NULL", "NULL", "factor", "NULL", "factor", rep("NULL", 5)
))


# The following file can be downloaded from:
# "https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz"
expression <- read.delim("GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz", colClasses = c(
  "NULL", "character", "NULL", "NULL", rep("numeric", 462)
))

E2F2expression <- expression[which(expression$Gene_Symbol == "ENSG00000007968.6"), -1]
rm(expression)

matchID <- match(colnames(E2F2expression), demographics$Individual.ID)
data <- demographics[matchID,]
rm(demographics, matchID)

data$expression <- as.numeric(E2F2expression[1,])

# Remove cases that have no expression level measured
data <- data[!is.na(data$expression),]

y <- as.vector(scale(data$expression))
rm(E2F2expression)


# Package to read in .vcf files
if(!require(vcfR)) {
  install.packages("vcfR")
  library(vcfR)
}


# Maximum number `q` of SNPs that we would select using sure independence screening (SIS)
q.max <- 1e4

# Count total number of variants
total.n.variants <- 0

for(chr in 1:22) {
  cat("Working on chromosome", chr, "\n")

  gc(reset = TRUE)

  tmp <- read.vcfR(
    file = paste0(
      "https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/GEUVADIS.chr",
      chr,
      ".PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz"
    ),
    limit = 1e11
  )

  tmp <- extract.gt(tmp, IDtoRowNames = FALSE)

  tmp[tmp %in% c("0|0", "0/0")] <- 0
  tmp[tmp %in% c("0|1", "0/1", "1|0", "1/0")] <- 1
  tmp[tmp %in% c("1|1", "1/1")] <- 2
  mode(tmp) <- "numeric"

  # Delete constant rows from tmp
  n.SNPs <- nrow(tmp)
  del <- logical(n.SNPs)
  
  pb <- txtProgressBar(max = n.SNPs, style = 3)
  for(i in 1:n.SNPs) {
    del[i] <- any(is.na(tmp[i,])) || all(tmp[i,] == 0) || all(tmp[i,] == 1) || all(tmp[i,] == 2)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  tmp <- tmp[!del,]
  
  # Standardize before computing the SIS score
  tmp <- scale(t(tmp))
  tmp <- data.frame(tmp, Individual.ID = rownames(tmp))
  data <- merge(x = data, y = tmp, by = "Individual.ID", all.x = TRUE, all.y = FALSE)
  
  
  print("Only keep q.max most important SNPs as determined by SIS")
  col.data <- ncol(data)
  col.tmp <- ncol(tmp)-1
  rm(tmp)
  
  
  # Removing those with missing expressions created variants constant accross all cases
  # After scaling, these appear as NaN
  data <- data[, !is.na(data[1,])]
  
  col.tmp <- col.tmp - (col.data - ncol(data))
  col.data <- ncol(data)
  
  total.n.variants <- total.n.variants + col.tmp
  
  if(col.tmp > q.max) {
    print("More variants remaining than q.max.")
    SIS <- abs( crossprod(as.matrix(data[,(col.data-col.tmp+1):col.data]), y) )
  
    SIS_cutoff = sort(SIS, decreasing = TRUE)[q.max]
    SIS_index = which(SIS >= SIS_cutoff)
    
    data <- cbind(data[,1:(col.data-col.tmp)], data[,((col.data-col.tmp+1):col.data)[SIS_index]])
  }
}


## Analyzing the data

data$Population <- droplevels(data$Population)

# Relevel to make 'GBR', the most common population type the reference level
data$Population <- relevel(data$Population, ref = "GBR")

X <- model.matrix(expression ~ Gender + Population - 1, data = data)

# Remove the most common gender "Gender2" from X to obtain a design matrix without intercept
X <- X[,-2]

# Standardize X
X <- scale(X)


# Z_all contains all candidate SNPs from which we select a subset using SIS.
Z_all <- data[,-(1:4)]
rm(data)


Z_all <- scale(Z_all)


# Removing those with missing expressions created variants constant accross all cases
# After scaling, these appear as NaN
Z_all <- Z_all[, !is.nan(Z_all[1,])]

n <- length(y)
p <- ncol(X)


## Select q predictors from Z_all using SIS
q <- q.max

SIS <- abs( crossprod(Z_all, y) )

SIS_cutoff = sort(SIS, decreasing = TRUE)[q]
SIS_index = which(SIS >= SIS_cutoff)

Z <- Z_all[, SIS_index]


psi <- 1/n
lambda <- n/q/10


# Set prior for sigma.sq
a.0 <- 1
b.0 <- 1

IRGA_12_result <- BVS_IRGA_12(y, X, Z, lambda, psi, sigma.sq = NULL, max.iter = 1e3, min.rho = 1e-12, verbose = TRUE, rho = 1/2^4)



# Density evaluation for multivariate normal
# This is a simplified version of 'dmvnorm' from the 'mvtnorm' R package.
dmvnorm <- function (x, mean = numeric(p), sigma = diag(p), log = FALSE) {
  p <- length(x)
  dec <- chol(sigma)
  tmp <- backsolve(dec, x - mean, transpose = TRUE)
  logretval <- -sum(log(diag(dec))) - 0.5 * p * log(2 * pi) - 0.5 * sum(tmp*tmp)
  return(ifelse(log, logretval, exp(logretval)))
}


# Computation of posterior means and posterior model probabilities
# for Bayesian variable selection with error covariance Sigma.
full.posterior <- function(y, X, lambda, psi, Sigma) {
  Sigma.inv <- solve(Sigma)
  
  p = NCOL(X) # Number of predictors
  n = length(y) # Sample size
  
  probs <- numeric(2^p)
  means <- matrix(0, nrow = 2^p, ncol = p)
  
  probs[1] <- (1-lambda)^p * dmvnorm(x = y, sigma = Sigma) # The marginal likelihood of the model without predictors
  
  for(k in 1:(2^p-1)) {
    gam <- as.logical(intToBits(k)[1:p])
    
    gam.index <- which(gam)
    
    X.gam <- X[, gam.index]
    s.gam <- length(gam.index)
    
    probs[k+1] <- lambda^s.gam*(1-lambda)^(p-s.gam)* dmvnorm(x = y, sigma = Sigma + psi*X.gam%*%t(X.gam))
    
    means[k+1, gam.index] <- solve(crossprod(X.gam, Sigma.inv %*% X.gam) + diag(s.gam)/psi) %*% crossprod(X.gam, Sigma.inv%*%y)
  }
  
  return(list(
    probs = probs/sum(probs),
    means = means
  ))
}


psi.X <- 1
lambda.X = .5
fullResult <- full.posterior(y = IRGA_12_result$MtY - IRGA_12_result$mean, X = IRGA_12_result$MtX, lambda = lambda.X, psi = psi.X, Sigma = IRGA_12_result$sigma.sq*diag(p) + IRGA_12_result$covariance)

# Posterior mean of beta from IRGA
post_mean <- crossprod(fullResult$means, fullResult$probs)

inc_prob <- numeric(p)
for(j in 1:p) inc_prob[j] <- sum(fullResult$probs[fullResult$means[,j] == 0])


# The same as the above but with the effect of the SNPs set to zero
fullResult2 <- full.posterior(y = IRGA_12_result$MtY, X = IRGA_12_result$MtX, lambda = lambda.X, psi = psi.X, Sigma = IRGA_12_result$sigma.sq*diag(p))

# Posterior means of beta when the SNPs are ignored
post_mean2 <- crossprod(fullResult2$means, fullResult2$probs)

inc_prob2 <- numeric(p)
for(j in 1:p) inc_prob2[j] <- sum(fullResult2$probs[fullResult2$means[,j] != 0])


# Print the results for inspection.
post_mean
post_mean2
inc_prob
inc_prob2



# Computation of posterior inclusion probabilities (PIPs)
# using the expectation propagation (EP) method from Hernandez-Lobato et al. (2015, doi:10.1007/s10994-014-5475-7).
PIP_EP_geuvadis <- function(Y, X,
                            p0 = rep(0.5, ncol(X)), # prior inclusion probability, same as `lambda`
                            v = rep(1, ncol(X)), # slab variance, same as `psi`
                            sigma.sq = NULL) {
  
  # Load the EP code by Jose Miguel Hernandez Lobato modified to deal with vector-valued p0 and v.
  source("epBVS_vec.R")
  
  ret <- epBVSinternal_vec(X, Y, beta = 1/sigma.sq, p0, v)
  
  return(list(
    post_mean = ret$m,
    incProb = as.vector(1/(1+exp(-ret$phi)))
  ))
}


EPstart <- Sys.time()

EPresult <- PIP_EP_geuvadis(
  Y = y,
  X = cbind(X, Z),
  p0 = c(rep(lambda.X, p), rep(lambda, q)),
  v = c(rep(psi.X, p), rep(psi, q)),
  sigma.sq = IRGA_12_result$sigma.sq
)

EPend <- Sys.time()
print(paste("EP computation time:", EPend - EPstart))

EPresult$incProb[1:p]
EPresult$post_mean[1:p]



# Computation of posterior inclusion probabilities (PIPs) using a Gibbs sampler
PIP_Gibbs_geuvadis <- function(
  sigma.sq = NULL,
  beta.0 = post_mean,
  iter = 1e4,
  MCerror = FALSE
) {
  # `beta.0` is the initialization of the sampler.
  # `iter` is the number of MCMC iterations,
  # the first 10% of which are discarded as burnin iterations.
  # If `MCerror` is true, an estimate of the Monte Carlo standard
  # error is produced using overlapping batch means (Flegal & Jones, 2010, doi:10.1214/09-AOS735, Sec. 3)
  
  if (MCerror) gibbsStart <- Sys.time()
  
  X <- as.matrix(X)
  p <- ncol(X) # Number of predictors
  q <- ncol(Z) # Number of covariates
  n <- length(y) # Number of observations
  
  # Matrix to store the Gibbs samples of beta.
  beta.gibbs <- matrix(NA_real_, nrow = iter, ncol = p)
  
  # Current value for alpha
  alpha <- numeric(q)
  
  # Current value for beta
  if (is.null(beta.0)) {
    beta <- numeric(p)
  } else beta <- beta.0
  
  s2unknown <- is.null(sigma.sq) # Record whether sigma.sq is given.
  
  # Initialize sigma.sq
  if (s2unknown) sigma.sq <- b.0 / a.0
  
  l.lam <- log(1 - lambda) - log(lambda)
  l.lam.X <- log(1 - lambda.X) - log(lambda.X)
  column.norm.sq <- colSums(Z * Z)
  column.norm.sq.X <- colSums(X * X)
  
  pb <- txtProgressBar(max = iter, style = 3)
  for (k in 1:iter) {
    
    # Update beta
    y_tmp <- y - Z %*% alpha
    
    for (j in 1:p) {
      Xj <- X[,j]
      
      # Residual vector, recomputed in full every 100th time to avoid compounding numerical error
      if(j %% 1e2 == 1) {
        res <- y_tmp - X[,-j] %*% beta[-j]
      } else res <- res + Xj * beta[j]
      
      inner.prod <- sum(res * Xj)
      
      if (
        runif(1) < 1 / ( 1 + exp(
          l.lam.X + .5 * log1p(column.norm.sq.X[j] * psi.X / sigma.sq)
          - .5 / sigma.sq * inner.prod * inner.prod / (sigma.sq / psi.X + column.norm.sq.X[j])
        ))
      ) {
        # If beta is nonzero, resample it.
        post.var <- 1 / (1 / psi.X + column.norm.sq.X[j] / sigma.sq)
        beta[j] <- rnorm(1, mean = max(min(post.var * inner.prod / sigma.sq, 1e12), -1e12), sd = sqrt(post.var))
        res <- res - Xj * beta[j]
      } else beta[j] <- 0
    }
    
    beta.gibbs[k, ] <- beta
    
    # Update alpha
    y_tmp <- y - X %*% beta
    
    for (j in 1:q) {
      Zj <- Z[,j]

      # Residual vector, recomputed in full every 100th time to avoid compounding numerical error
      if(j %% 1e2 == 1) {
        res <- y_tmp - Z[,-j] %*% alpha[-j]
      } else res <- res + Zj * alpha[j]

      inner.prod <- sum(res * Zj)

      if (
        runif(1) < 1 / ( 1 + exp(
          l.lam + .5 * log1p(column.norm.sq[j] * psi / sigma.sq)
          - .5 / sigma.sq * inner.prod * inner.prod / (sigma.sq / psi + column.norm.sq[j])
        ))
      ) {
        # If alpha is nonzero, resample it.
        post.var <- 1 / (1 / psi + column.norm.sq[j] / sigma.sq)
        alpha[j] <- rnorm(1, mean = max(min(post.var * inner.prod / sigma.sq, 1e12), -1e12), sd = sqrt(post.var))
        res <- res - Zj * alpha[j]
      } else alpha[j] <- 0
    }
    
    # Update sigma.sq
    if (s2unknown) sigma.sq <- 1 / rgamma(n = 1, shape = a.0+n/2, rate = b.0 + sum(res * res))
    
    setTxtProgressBar(pb, k)
  }
  close(pb)
  
  gamma.gibbs <- beta.gibbs != 0
  
  if (MCerror) {
    cat("Gibbs computation time:")
    print(Sys.time() - gibbsStart)
    
    cat("Computing Monte Carlo standard error\n")
    
    # For `mcse`
    if (!'mcmcse' %in% rownames(installed.packages())) install.packages('mcmcse', dependencies = T)
    
    mSE <- 0
    pb <- txtProgressBar(max = p, style = 3)
    for (j in 1:p) {
      mSE <- mSE + mcmcse::mcse(x = gamma.gibbs[,j], method = "obm", warn = TRUE)$se
      setTxtProgressBar(pb, j)
    }
    close(pb)
    cat("Average Monte Carlo standard error:", mSE/p)
  }
  
  # The PIP estimates follow as the mean of gamma,
  # excluding the first 10% of iterations as burnin.
  return(list(
    post_mean = colMeans(beta.gibbs[-(1:(iter %/% 10)), ]),
    incProb = colMeans(gamma.gibbs[-(1:(iter %/% 10)), ])
  ))
}

set.seed(1)
GibbsResult <- PIP_Gibbs_geuvadis(sigma.sq = IRGA_12_result$sigma.sq, iter = 1e5, MCerror = TRUE)
GibbsResult



# Estimate beta using a mixed effects model.
if (!'lme4' %in% rownames(installed.packages())) install.packages('lme4', dependencies = T)

# Transform Z to a logical matrix: Either there is a variant or there isn't.
Z_logical <- apply(Z, MARGIN = 2, FUN = function(col) col != min(col))

# A too large number of  SNPs fails to converge.
# Let us focus on a subset using sure independence screening.
Z_logical_ordered <- Z_logical[, order(
  abs( crossprod(scale(Z_logical), y) ),
  decreasing = TRUE
)]

q2 <- 2e3
tmp <- Z_logical_ordered[, 1:q2]
n.groups <- sum(!duplicated(tmp))

groups <- numeric(n)

for (g in 1:n.groups) {
  current.row <- tmp[which(groups == 0)[1], ]
  for (i in 1:n) if (all(tmp[i, ] == current.row)) groups[i] <- g
}


LMMstart <- Sys.time()
lmm <- lme4::lmer(y ~ X + (1 | groups))
LMMend <- Sys.time()

print(paste("LMM computation time:", LMMend - LMMstart))
lmm