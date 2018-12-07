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
p <- ncol(A)

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

vbStart <- Sys.time()
VB_PIP <- PIP_VB(y, X = A, lambda, psi, sigma.sq = NULL)
vbEnd <- Sys.time()

Ormerod_PIP <- PIP_VB_Ormerod(y, X = A, lambda, psi, sigma.sq = NULL)

epStart <- Sys.time()
EP_PIP <- PIP_EP(Y = y, X = A, lambda, psi, sigma.sq = NULL)
epEnd <- Sys.time()

Gibbs_PIP <- PIP_Gibbs(y, X = A, lambda, psi, sigma.sq = NULL, iter = 1e5, MCerror = TRUE)



cat("IRGA absolute difference in posterior inclusion probability:")
summary(abs(Gibbs_PIP - IRGA_PIP))
print(irgaEnd - irgaStart)

cat("VB absolute difference in posterior inclusion probability:")
summary(abs(Gibbs_PIP - VB_PIP))
print(vbEnd - vbStart)

summary(abs(Gibbs_PIP - Ormerod_PIP))

cat("EP absolute difference in posterior inclusion probability:")
summary(abs(Gibbs_PIP - EP_PIP))
print(epEnd - epStart)


## Plot the results
ymax <- .4
colGibbs= rgb(0,0,0,alpha=.5) 
colVB= rgb(0,0,0,alpha=.7)

pdf("application_diabetes.pdf", width = 5.5, height = 2)
par(mfrow = c(1,1), bty = 'l',
    cex = .5)
ylab <- "Error in PIP"
plot(x = 1:64-.25, y = abs(IRGA_PIP - Gibbs_PIP), type = 'h', cex = .6, ylim = c(0, ymax), xlim = c(1, p), ylab = "Difference in PIP", xlab = "Predictor")
lines(x = 1:64, y = abs(EP_PIP - Gibbs_PIP), type = 'h', col = colGibbs, cex = .7)
lines(x = 1:64+.25, y = abs(VB_PIP - Gibbs_PIP), type = 'h', col = colVB, cex = .7)

points(x = 1:64-.25, y = abs(IRGA_PIP - Gibbs_PIP), pch = 4, cex = .3)
points(x = 1:64, y = abs(EP_PIP - Gibbs_PIP), pch = 16, col = colGibbs, cex = .5)
points(x = 1:64+.25, y = abs(VB_PIP - Gibbs_PIP), pch = 17, col = colVB, cex = .5)

dev.off()