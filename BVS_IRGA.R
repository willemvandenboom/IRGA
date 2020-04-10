# For intToBin()
if(!'R.utils' %in% rownames(installed.packages())) install.packages('R.utils', dependencies = T)

# For parallel computing
if(!'doParallel' %in% rownames(installed.packages())) install.packages('doParallel', dependencies = T)


## Exact computation of posterior inclusion probabilities (PIPs)
# with non-iid errors with covariance matrix Sigma
PIP_exact <- function(y, X, lambda, psi, Sigma) {
  log.sum.exp <- function(a, b) {
    # Computes log(exp(a)+exp(b)), vectorizing over b argument.
    # Uses offset trick to avoid numeric overflow.
    min.ab <- pmin(a, b)
    max.ab <- pmax(a, b)
    
    return(log1p(exp(min.ab-max.ab)) + max.ab)
  }
  
  X <- as.matrix(X)
  p <- ncol(X) # Number of candidate variables
  y <- as.vector(y)
  
  # Log-density evaluation for a zero-mean multivariate normal
  # This is a simplified and faster version of `dmvnorm` from the `mvtnorm` R package.
  # The normalization constant involving pi has been left out as
  # we only need the density up to a proportionality constant.
  ldmvnorm <- function(Sigma) {
    dec <- chol(Sigma)
    tmp <- backsolve(dec, y, transpose = TRUE)
    return(-sum(log(diag(dec))) - 0.5 * sum(tmp*tmp))
  }
  
  # The marginal log-likelihood of the model without predictors
  LL <- p*log(1-lambda) + ldmvnorm(Sigma)
  
  # Object to store the marginal log-likelihood that a beta_j is *not* included
  LLj <- rep(LL, p)
  
  for(k in 1:(2^p-1)) {
    gam <- integer(p)
    gam[min(p-ceiling(log2(k+1)-1),p):p] <- as.numeric(strsplit(R.utils::intToBin(k), split = "")[[1]])
    index <- which(gam==1)
    X.gam <- X[,index, drop=FALSE] # `drop=FALSE` such that X.gam does not reduce to a vector if s.gam=1.
    s.gam <- length(index)
    LL.gam <- s.gam*log(lambda) + (p-s.gam)*log(1-lambda) + ldmvnorm(Sigma + psi*tcrossprod(X.gam))
    
    LL <- log.sum.exp(a=LL.gam, b=LL)
    LLj[-index] <- log.sum.exp(a=LL.gam, b=LLj[-index])
  }
  
  return(-expm1(LLj - LL)) # 1 - exp(LLj)/exp(LL)
}




## Implementation of IRGA

BVS_IRGA_12 <- function(
  y, X, Z, lambda, psi, sigma.sq = NULL, max.iter = 100,
  min.rho = 1, rho = 1, estimator = "VAMP", verbose = FALSE
) {
  # Steps 1 and 2 of
  # Integrated rotated Gaussian approximation (IRGA, Algorithm 1)
  # using spike-and-slab priors as described in Section 3.1,
  # with a spike-and-slab prior on alpha as well:
  # Pi(alpha) = lambda * N(0, psi) + (1-lambda) delta(0).
  
  p <- NCOL(X) # Length of beta
  n <- length(y) # Number of observations
  
  
  ## Step 1
  # Compute the QR decomposition of X
  # to obtain the rotation matrix Q.
  Q <- qr.Q(qr(X), complete=TRUE)
  
  # Compute the rotated quantities
  QtY <- as.vector(crossprod(Q, y))
  MtY <- QtY[1:p]
  
  MtX <- crossprod(Q[,1:p, drop=FALSE], X)
  
  QtZ <- crossprod(Q, Z)
  MtZ <- QtZ[1:p,, drop=FALSE]
  
  # If p equals the number of observations,
  # then the nuisance parameter eta follows its prior
  # since StZ is empty then.
  if(p == n) {
    if(is.null(sigma.sq)) stop(
      'sigma.sq is unknown while IRGA is called with p = n such that sigma.sq cannot be inferred.'
    )
    
    # The prior variance of alpha is lambda * psi * I_q
    covariance <- lambda * psi * tcrossprod(MtZ)
    
    return(list(MtY=MtY, MtX=MtX, mean=numeric(p), covariance=covariance, sigma.sq=sigma.sq))
  }
  
  StY <- QtY[-(1:p)]
  StZ <- QtZ[-(1:p),, drop=FALSE]
  rm(y, QtY, QtZ)
  
  
  ## Step 2
  # Estimate the mean and variance of the nuisance parameter alpha
  # using vector approximate message passing (VAMP) or lasso.
  
  if(estimator == "VAMP") {
    
    # Load the function for vector approximate message passing (VAMP)
    source('VAMP.R')
    
    VAMPresult <- VAMP(
      y=StY, X=StZ, lambda, psi, sigma.sq,
      rho = rho, max.iter = max.iter,
      min.rho = min.rho, verbose = verbose
    )
    mean <- MtZ %*% VAMPresult$mean
    covariance <- VAMPresult$variance * tcrossprod(MtZ, MtZ)
    
    # If sigma.sq was unknown, it has been estimated by VAMP.
    sigma.sq <- VAMPresult$sigma.sq
    
  } else {
    
    # Use the debiased lasso (Javanmard & Montanari, 2013)
    q <- NCOL(StZ)
    
    # Keep the largest lasso estimate which has less than lambda*q+1 nonzero coefficients.
    lars.result <- coef(lars::lars(x=StZ, y=StY, type = "lasso", intercept = FALSE, use.Gram = q <= n))
    alpha.lasso <- lars.result[max(which(rowSums(lars.result != 0) < lambda*q+1)),]
    
    lasso.u <- as.vector(alpha.lasso+t(StZ)%*%(StY-StZ%*%alpha.lasso)/(n-p))
    Sigma.hat <- t(StZ)%*%StZ/(n-p)
    mean <- as.vector((diag(q)-Sigma.hat)%*%alpha.lasso)
    
    a.0 <- 1
    b.0 <- 1
    sigma.sq <- ( b.0+sum((StY - StZ %*% alpha.lasso)^2)/2 ) / ( a.0 + (n-p)/2 )
    
    covariance <- sigma.sq*solve(Sigma.hat)/(n-p)
    
    mean <- MtZ %*% mean
    covariance <- tcrossprod(MtZ %*% covariance, MtZ)
    
  }
  
  return(list(MtY=MtY, MtX=MtX, mean=mean, covariance=covariance, sigma.sq=sigma.sq))
}


BVS_IRGA <- function(y, X, Z, lambda, psi, sigma.sq = NULL, max.iter = 100, estimator = "VAMP") {
  # Integrated rotated Gaussian approximation (IRGA, Algorithm 1)
  # using spike-and-slab priors as described in Section 3.1,
  # with a spike-and-slab prior on beta as well:
  # Pi(beta) = lambda * N(0, psi) + (1-lambda) delta(0) and
  # Pi(alpha) = lambda * N(0, psi) + (1-lambda) delta(0).
  # `max.iter` is the maximum number of VAMP iterations before adding dampening.
  # `estimator` specifies whether the mean and variance of (4)
  # are estimated using VAMP or lasso.
  
  p <- NCOL(X) # Length of beta
  
  ## Steps 1 and 2 of Algorithm 1:
  BVS_IRGA_12result <- BVS_IRGA_12(
    y, X, Z, lambda, psi, sigma.sq, max.iter, estimator = estimator
  )
  
  
  ## Step 3 of Algorithm 1:
  # Compute the posterior inclusion probabilities (PIPs).
  return(PIP_exact(
    y = BVS_IRGA_12result$MtY - BVS_IRGA_12result$mean,
    X = BVS_IRGA_12result$MtX,
    lambda, psi,
    # If sigma.sq was unknown, it has been estimated in Step 2.
    Sigma <- BVS_IRGA_12result$sigma.sq * diag(p) + BVS_IRGA_12result$covariance
  ))
}


# Computation of posterior inclusion probabilities (PIPs) using IRGA
# with iid errors with variance `sigma.sq`, potentially unknown.
PIP_IRGA <- function(
  y, A, lambda, psi, sigma.sq = NULL, p = NULL,
  max.iter = 100, estimator = "VAMP", theta_split = "sequential"
) {
  # `p` is the (maximum) length of beta in the executions of IRGA
  # `max.iter` is the maximum number of VAMP iterations before adding dampening.
  # `theta_split` specifies how A is split into X and Z.
  # The options are "sequential", "random", "spectral" and "belsley".
  # The latter two group correlated predictors in X.
  # "belsley" cannot be used if r > n.
  
  A <- as.matrix(A)
  r <- ncol(A) # Number of candidate variables
  y <- as.vector(y)
  
  # `n_splits` specifies the number of executions of IRGA
  if (is.null(p)) {
    
    # p = O(log(r)) yields a favorable computational complexity.
    p = floor(log(r))
    n_splits <- ceiling(r / p)
    
    # To make best use of compute power,
    # the number of splits should not be less than the number of CPU cores.
    n_splits <- max(n_cores, n_splits)
    
    # The number of splits cannot be so small that p > n,
    # because then IRGA does not apply.
    n_splits <- max(r %/% length(y), n_splits)
    
  } else {
    if (p > n) stop("p>n such that IRGA does not apply.")
    n_splits <- ceiling(r / p)
  }
  
  # Compute the size of each split.
  split_size <- rep(ceiling(r / n_splits), n_splits)
  
  # If r %% n_splits != 0, we need to make the splits of unequal size.
  diff <- sum(split_size) - r
  if (diff > 0) split_size[n_splits - 0:(diff - 1)] <- split_size[n_splits - 0:(diff - 1)] - 1
  
  theta_split <- match.arg(tolower(theta_split), c("sequential", "random", "spectral", "belsley"))
  
  # List to store the column numbers of A constituting each beta
  beta_ind <- list()
  
  if (theta_split == "sequential") {
    
    for (s in 1:n_splits) beta_ind[[s]] <- sum(split_size[seq_len(s - 1)]) + 1:split_size[s]
    
  } else if (theta_split == "random") {
    
    tmp <- sample.int(n = r)
    for (s in 1:n_splits) beta_ind[[s]] <- tmp[sum(split_size[seq_len(s - 1)]) + 1:split_size[s]]
    
  } else if (theta_split == "spectral") {
    
    # Spectral clustering
    
    # Create the similarity matrix.
    W <- diag(r)
    for (i in 2:r) for (j in 1:(i - 1)) W[i, j] <- cor(A[, i], A[, j])
    W <- abs(W)
    W <- W + t(W)
    
    # Compute the Laplacian matrix.
    L <- diag(colSums(W)) - W
    
    U <- eigen(x = L, symmetric = TRUE)$vectors[, (r - n_splits):(r - 1)]
    
    # Need to ensure that the clusters have equal size.
    # We first run k-means:
    tmp <- kmeans(x = U, centers = n_splits, iter.max = 100)
    cluster_centers <- tmp$centers
    
    # Then, we assign "surplus" of the largest cluster to the closest clusters:
    
    # `theta_cluster` is a factor to ensure also empty clusters are included in the output of `table`.
    theta_cluster <- as.factor(tmp$cluster)
    
    # Vector to keep track of which clusters have already been done
    cluster_done <- logical(n_splits)
    
    for (s in 1:n_splits) {
      
      if (s == 1) {
        cluster_table <- table(theta_cluster)
      } else {
        # We do not consider the clusters that we already reduced to the correct size.
        cluster_table <- table(theta_cluster[-unlist(beta_ind)])
      }
      
      cluster_index <- which.max(cluster_table)
      theta_ind <- which(theta_cluster == cluster_index)
      cluster_size <- cluster_table[cluster_index]
      
      # If the largest cluster is small enough, then we are done.
      if (cluster_size <= split_size[1]) {
        for (cl in which(!cluster_done)) {
          beta_ind[[s]] <- which(theta_cluster == cl)
          s <- s+1
        }
        break
      }
      
      # Compute the distance to the center for each point in the largest cluster.
      distance <- numeric(cluster_size)
      for (i in 1:cluster_size) distance[i] <- sum((U[theta_ind[i], ] - cluster_centers[cluster_index, ])^2)
      
      # Find the `split_size[1]` smallest distances.
      tmp <- order(distance)[1:split_size[1]]
      
      # The points corresponding to these smallest distances remain in the cluster.
      beta_ind[[s]] <- theta_ind[tmp]
      cluster_done[cluster_index] <- TRUE
      
      # The other points are assigned to the closest remaining clusters.
      theta_ind <- theta_ind[-tmp]
      
      # Vector to keep track of minimum distance
      min_distance <- rep(Inf, cluster_size - split_size[1])
      
      for (cl in which(!cluster_done)) for (i in 1:(cluster_size - split_size[1])) {
        
        tmp <- sum((U[theta_ind[i], ] - cluster_centers[cl, ])^2)
        
        if (tmp < min_distance[i]) {
          min_distance[i] <- tmp
          theta_cluster[theta_ind[i]] <- cl
        }
        
      }
      
    }
    
  } else { # theta_split == "belsley"
    
    if (r > n) stop("`theta_split = 'belsley'` cannot be used if `r > n`.")
    
    # Cluster collinear predictors based on Section 3.2 from
    # "Regression Diagnostics: Identifying Influential Data and Sources of Collinearity"
    # by David A. Belsley, Edwin Kuh, and Roy E. Welsch (1980).
    tmp <- svd(x = A, nu = 0, nv = r)
    phi <- tmp$v^2 %*% diag(1 / tmp$d^2)
    pi_mat <- t(diag(1 / rowSums(phi)) %*% phi)
    
    for (s in 1:n_splits) {
      beta_ind[[s]] <- order(pi_mat[s, ], decreasing = TRUE)[1:split_size[s]]
      pi_mat[, beta_ind[[s]]] <- NA
    }
    
  }
  
  # Setup for parallel computing
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  
  # Compute the PIPs using IRGA and parallel computing.
  PIP_list <- foreach::`%dopar%`(
    foreach::foreach(
      s=1:n_splits,
      .export = c('BVS_IRGA', 'BVS_IRGA_12', 'PIP_exact')
    ),
    {
      BVS_IRGA(
        y,
        X = A[, beta_ind[[s]]],
        Z = A[, -beta_ind[[s]]],
        lambda, psi, sigma.sq, max.iter,
        estimator = estimator
      )
    }
  )
  
  parallel::stopCluster(cl)
  
  PIP <- numeric(r)
  for (s in 1:n_splits) PIP[beta_ind[[s]]] <- PIP_list[[s]]
  
  return(PIP)
}