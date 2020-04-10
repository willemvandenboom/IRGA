## This is an edited version of epBVS.R that can deal with vector-valued p0 and v.

#########################################################################################################
# Module: epBVS.R
# Date  : January 2010
# Author: Jose Miguel Hernandez Lobato
# email : josemiguel.hernandez@uam.es
#
# Module that approximates the posterior of the parameters of the the Bayesian
# variable selection model for the linear regression problem. The procedure used to perform
# the approximation is Expectation Propagation.
#
#########################################################################################################
#
# EXPORTED FUNCTIONS: epBVS
#
#########################################################################################################
#
# The main function that should be used in this module is "epBVS". You have to call it with the arguments:
#
# 	X -> Design matrix for the regression problem.
#       Y -> Target vector for the regression problem.
#	beta -> Noise precision.
#	p0 -> Prior probability that a feature is relevant for solving the regression problem.
#	v -> Variance of the slab.
#
# "epBVS" returns the approximate distribution for the posterior as a list with components:
#
#	m -> Mean vector for the marginals of the posterior.
#	v -> Variance vector for the marginals of the posterior.
#	phi -> Vector with the marginal prbobabilities of activation of the latent variables.
#
#	t1Hat -> A list with the approximation for the likelihood.
#		mHat -> Mean vector of the factorized Gaussian approximation.
#		vHat -> Variance vector for the factorized Gaussian approximation.
#
#	t2Hat -> A list with the approximation for the spike and slab prior.
#		mHat -> Mean vector of the factorized Gaussian approximation.
#		vHat -> Variance vector for the factorized Gaussian approximation.
#		phiHat -> Parameter vector for the Bernoulli approximation.
#
#	t3Hat -> A list with the approximation for Bernoulli prior.
#		phiHat -> Parameter vector for the Bernoulli approximation.
#
#	evidence -> The approximation for the evidence given by EP.
#

invLogistic <- function(x) log(x / (1 - x))

epBVS <- function(X, Y, beta = 1, p0 = 0.5, v = 1) {

	# We find the optimal configuration for the hyper-parameters

	target <- function(param) {
		
		beta <- exp(param[ 1 ])
		p0 <- logistic(param[ 2 ])
		v <- exp(param[ 3 ])

		-epBVSinternal(X, Y, beta, p0, v)$evidence
	}

	# We call the optimization method

	startPoint <- c(log(beta), invLogistic(p0), log(v))
	ret <- optim(startPoint, target, method = "Nelder-Mead", control = list(trace = T, maxit = 40))$par

	time <- system.time(ret <- epBVSinternal(X, Y, exp(ret[ 1 ]), logistic(ret[ 2 ]), exp(ret[ 3 ])))

	ret$timeInternal <- time[[ 1 ]] + time[[ 2 ]]

	ret
}

epBVSinternal_vec <- function(X, Y, beta = 1, p0 = rep(0.5, ncol(X)), v = rep(1, ncol(X))) {

	d <- ncol(X)
	n <- nrow(X)

	# Precomputation

	tXX <- t(X) %*% X
	tXY <- t(X) %*% Y

	# We initialize the approximation

	t1Hat <- list(mHat = rep(0, d), vHat = rep(Inf, d))
	t2Hat <- list(mHat = rep(0, d), vHat = rep(Inf, d), phiHat = rep(0, d))
	t3Hat <- list(phiHat = rep(0, d))
	a <- list(m = rep(0, d), v = rep(Inf, d), phi = rep(0, d), p = rep(NA, d), t1Hat = t1Hat, t2Hat = t2Hat, t3Hat = t3Hat, indexNegative = c())

	# We process the approximate term for the Bernoulli prior

	a$t3Hat$phiHat <- log(p0 / (1 - p0))

	# We process the approximate term for the spike and slab prior

	a$t2Hat$vHat <- p0 * v

	# Main loop of ep. Repeated till the algorithm reaches convegence

	i <- 1
	damping <- 0.99
	convergence <- FALSE
	while(!convergence && i < 1000) {

		aOld <- a

		# We refine the approximate term for the likelihood

		if (d > n) {
			inverseWoodbury <- solve(diag(rep(beta^-1, n)) + (X * matrix(a$t2Hat$vHat, n, d, byrow = T)) %*% t(X))
			vectorAux <- a$t2Hat$vHat^-1 * a$t2Hat$mHat + beta * tXY
			a$m <- as.double(a$t2Hat$vHat * (vectorAux - (t(X) %*% (inverseWoodbury %*% (X %*% (a$t2Hat$vHat * vectorAux))))))
			a$v <- as.double(a$t2Hat$vHat - a$t2Hat$vHat^2 * (rep(1, n) %*% (X * (inverseWoodbury %*% X))))
		} else {
			Sigma <- solve(diag(as.double(a$t2Hat$vHat)^-1) + beta * tXX)
			a$m <- Sigma %*% (a$t2Hat$vHat^-1 * a$t2Hat$mHat + beta * tXY)
			a$v <- diag(Sigma)
		}

		a$t1Hat$mHat <- (damping * (a$m * a$v^-1 - a$t2Hat$mHat * a$t2Hat$vHat^-1) + (1 - damping) * a$t1Hat$mHat * a$t1Hat$vHat^-1)
		a$t1Hat$vHat <- 1 / (damping * (1 / a$v - 1 / a$t2Hat$vHat) + (1 - damping) / a$t1Hat$vHat)
		a$t1Hat$mHat <- a$t1Hat$vHat * a$t1Hat$mHat

		damping <- damping * 0.99

		# We refine the approximate term for the spike and slab prior

		phiHatNew <- 0.5 * log(a$t1Hat$vHat) - 0.5 * log(a$t1Hat$vHat + v) + 0.5 * a$t1Hat$mHat^2 * (a$t1Hat$vHat^-1 - (a$t1Hat$vHat + v)^-1)
		aa <- logistic(phiHatNew + a$t3Hat$phiHat) * a$t1Hat$mHat * (a$t1Hat$vHat + v)^-1 +
			logistic(-phiHatNew - a$t3Hat$phiHat) * a$t1Hat$mHat * a$t1Hat$vHat^-1
		bb <- logistic(phiHatNew + a$t3Hat$phiHat) * (a$t1Hat$mHat^2 - a$t1Hat$vHat - v) * (a$t1Hat$vHat + v)^-2 +
			logistic(-phiHatNew - a$t3Hat$phiHat) * (a$t1Hat$mHat^2 * a$t1Hat$vHat^-2 - a$t1Hat$vHat^-1)

		vHatNew <- (aa^2 - bb)^-1 - a$t1Hat$vHat
		mHatNew <- a$t1Hat$mHat - aa * (vHatNew + a$t1Hat$vHat)

		a$indexNegative <- which(vHatNew < 0)

		# We minimize the KL divergence with vHatNew constrained to be positive.

		vHatNew[ a$indexNegative ] <- 100
		mHatNew[ a$indexNegative ] <- a$t1Hat$mHat[ a$indexNegative ] - aa[ a$indexNegative ] *
			(vHatNew[ a$indexNegative ] + a$t1Hat$vHat[ a$indexNegative ])

		a$t2Hat$phiHat <- phiHatNew * damping + a$t2Hat$phiHat * (1 - damping)
		a$t2Hat$mHat <- damping * mHatNew * vHatNew^-1 + (1 - damping) * a$t2Hat$mHat * a$t2Hat$vHat^-1
		a$t2Hat$vHat <- 1 / (damping / vHatNew + (1 - damping) / a$t2Hat$vHat)
		a$t2Hat$mHat <- a$t2Hat$mHat * a$t2Hat$vHat

		# We compute the posterior approximation from the approximate terms

		a$v <- 1 / (1 / a$t1Hat$vHat + 1 / a$t2Hat$vHat)
		a$m <- a$v * (a$t1Hat$mHat / a$t1Hat$vHat + a$t2Hat$mHat / a$t2Hat$vHat)
		a$phi <- a$t2Hat$phiHat + a$t3Hat$phiHat
		a$p <- logistic(a$phi)

		convergence <- checkConvergence(a, aOld)

		i <- i + 1
	}

	# We compute the evidence

	a$evidence <- computeEvidence(a, Y, X, beta, v)

	a$beta <- beta
	a$p0 <- p0
	a$vSlab <- v

	a$Sigma <- solve(diag(as.double(a$t2Hat$vHat)^-1) + beta * tXX)

	# We return the current approximation

	a
}

##
# The logistic function
#

logistic <- function(x) {

	1 / (1 + exp(-x))
}

##
# Checks convergence of the EP algorithm.
#
# Input:
# 	aOld -> The previous approximation.
# 	aNew -> The new approximation.
# Output:
# 	TRUE if the values in aOld are differ from those in aNew by less than a small constant.
#

checkConvergence <- function(aNew, aOld) {

	tol <- 1e-4

	convergence <- max(max(abs(aNew$m - aOld$m)))
	convergence <- max(convergence, max(abs(aNew$v - aOld$v)))

	print(convergence)

	if (convergence < tol)
		TRUE
	else
		FALSE
}

##
# Function that computes the log evidence
#

computeEvidence <- function(a, Y, X, beta, v) {

	n <- nrow(X)
	d <- ncol(X)

	vPrior <- a$tHatPriorSS$vHat
	mPrior <- a$tHatPriorSS$mHat

	# We compute the logarithm of s1 and s2

	if (n > d)
		alpha <- det(diag(as.double(a$t2Hat$vHat)) %*% t(X) %*% X * beta + diag(rep(1, d)))
	else
		alpha <- det((matrix(a$t2Hat$vHat, n, d, byrow = T) * X) %*% t(X) * beta + diag(rep(1, n)))
	
	logs1 <- -n / 2 * log(2 * pi / beta) - 0.5 * beta * sum(Y^2) +
		0.5 * sum((a$t2Hat$vHat^-1 * a$t2Hat$mHat + t(X) %*% Y * beta) * a$m) -0.5 * sum(a$t2Hat$mHat^2 * a$t2Hat$vHat^-1) - 0.5 * log(alpha) +
		1 / 2 * sum(log(1 + a$t2Hat$vHat * a$t1Hat$vHat^-1)) +
		1 / 2 * sum(a$t2Hat$mHat^2 * a$t2Hat$vHat^-1 + a$t1Hat$mHat^2 * a$t1Hat$vHat^-1 - a$m^2 * a$v^-1)

	c <- logistic(a$t3Hat$phiHat) * dnorm(0, a$t1Hat$mHat, sqrt(a$t1Hat$vHat + v)) + logistic(-a$t3Hat$phiHat) *
		dnorm(0, a$t1Hat$mHat, sqrt(a$t1Hat$vHat))

	logs2 <- sum(log(c) + 1 / 2 * log(1 + a$t1Hat$vHat * a$t2Hat$vHat^-1) +
		 1 / 2 * (a$t2Hat$mHat^2 * a$t2Hat$vHat^-1 + a$t1Hat$mHat^2 * a$t1Hat$vHat^-1 - a$m^2 * a$v^-1) +
		 log(logistic(a$phi) / logistic(a$t3Hat$phiHat) + logistic(-a$phi) / logistic(-a$t3Hat$phiHat)))

	aux <- d / 2 * log(2 * pi) + 0.5 * sum(log(a$v)) - 0.5 * sum(a$t1Hat$mHat^2 / a$t1Hat$vHat) -
		0.5 * sum(a$t2Hat$mHat^2 / a$t2Hat$vHat) + 0.5 * sum(a$m^2 / a$v)

	logs1 + logs2 + aux + sum(log(logistic(a$t2Hat$phiHat) * logistic(a$t3Hat$phiHat) + logistic(-a$t2Hat$phiHat) * logistic(-a$t3Hat$phiHat)))
}
