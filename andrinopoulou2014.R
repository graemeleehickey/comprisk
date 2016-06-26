##*********************************************************
## Data
##*********************************************************

library(joineR)
library(splines)
library(nlme)
library(R2WinBUGS)
library(arm)

data(epileptic)
head(epileptic)

epileptic$with.time <- epileptic$with.time / 365.25
epileptic$time <- epileptic$time / 365.25
epileptic$treat <- as.numeric(epileptic$treat == "LTG")
epileptic$event <- with(epileptic, 
  as.numeric(with.status.uae | with.status.isc))

# Time-to-event data (one row per subject)
survdat <- epileptic
survdat <- survdat[!duplicated(survdat$id), ]

# Number of subjects
N <- nrow(survdat) 

# Number of measurements per subject
ni <- as.vector(as.numeric(table(epileptic$id))) # Num
offset <- c(1, 1 + cumsum(ni))

# Longitudinal outcome
y <- epileptic$dose

# Design matrix for longitudinal fixed effects
X <- with(epileptic, data.frame(
  rep(1, nrow(epileptic)),
  time,
  treat,
  time*treat
))
X <- as.matrix(X)
colnames(X) <- NULL

# Design matrix for longitudinal random effects
Z <- with(epileptic, data.frame(
  rep(1, nrow(epileptic)),
  time
))
Z <- as.matrix(Z)
colnames(Z) <- NULL

# Number of random effects
nb <- ncol(Z)

# Event indicators for UAE and ISC
eventR <- survdat$with.status.uae
eventD <- survdat$with.status.isc

# Vector of zeroes
zeros <- rep(0, N)

# Cause-specific design matrices for time-to-event fixed effects
WR <- survdat$treat
WD <- WR

# Gauss-Kronrod points
gaussKronrod <- JMbayes:::gaussKronrod
wk <- gaussKronrod()$wk
sk <- gaussKronrod()$sk
ordsk <- order(sk)
sk <- sk[ordsk]
wk <- wk[ordsk]
K <- length(sk)

# Failure-times (divided by 2)
P <- survdat$with.time / 2

# Number of baseline hazard splines
Q <- 5

# Design matrices for baseline hazard functions
st <- outer(P, sk + 1)
qs <- with(survdat, quantile(with.time[event == 1], seq(0, 1, len = Q + 1), names = FALSE)[-c(1, Q + 1)])
qs <- c(0, qs, max(survdat$with.time) + 1)
W2R <- splineDesign(qs, survdat$with.time, ord = 1)
W2D <- W2R
W2sR <- splineDesign(qs, c(t(st)), ord = 1)
W2sD <- W2sR

# Design matrices for baseline hazard functions (old code)
# qs <- with(survdat, quantile(with.time[event == 1], seq(0, 1, len = Q + 1), names = FALSE)[-c(1, Q + 1)])
# qs <- c(0, qs, max(survdat$with.time) + 1)
# ind <- findInterval(survdat$with.time, qs, rightmost.closed = TRUE)
# D <- matrix(0, length(ind), Q)
# D[cbind(seq_along(ind), ind)] <- 1
# W2R <- D
# W2D <- D
# W2sR <- W2R[rep(1:N, each = 15), ]
# W2sD <- W2D[rep(1:N, each = 15), ]

# Column lengths of matrices
ncW2R <- ncol(W2R)
ncW2D <- ncol(W2D)
ncX <- ncol(X)
ncZ <- ncol(Z)
#ncWR <- ncol(WR)
#ncWD <- ncol(WD)

# Constant for zeroes-trick
C <- 5000

model.data <- list(
  # Main model data (observed)
  N = N, K = K,
  offset = offset,
  X = X, Z = Z, y = y,
  eventR = eventR, eventD = eventD,
  zeros = zeros,
  WR = WR, WD = WD,
  W2R = W2R, W2D = W2D, W2sD = W2sD, W2sR = W2sR,
  ncX = ncX, ncZ = ncZ, 
  ncW2R = ncW2R, ncW2D = ncW2D,
  C = C, P = P,
  wk = wk,
  nb = nb,
  # Prior hyper-parameters below
  mu0 = rep(0, nb),
  priorMean.betas = rep(0, ncX),
  priorTau.betas = diag(rep(0.01, ncX)),
  priorA.tau = 0.001,
  priorB.tau = 0.001,
  priorMean.gammas = 0,
  priorTau.gammas = 0.01,
  priorMean.alphas = rep(0, nb),
  priorTau.alphas = diag(rep(0.01, nb)),
  priorMean.Bs.gammas = rep(0, Q),
  priorTau.Bs.gammas = diag(rep(0.01, Q)),
  priorR.D = diag(rep(1, nb)),
  priorK.D = 2
)

##*********************************************************
## Fit model
##*********************************************************

inits <- function() {
  list(
    betas = rnorm(4, 0, 0.2) + c(1.9, 0.15, -0.1, 0.2),
    tau = runif(1, min = 4.1, max = 6.1),
    gammasR = rnorm(1, -0.9, 1),
    gammasD = rnorm(1, 0, 1),
    alphasR = rnorm(2),#c(-1.3, 1),
    alphasD = rnorm(2),#c(-0.2, 3),
    inv.D = matrix(c(1.4, -0.4, -0.4, 5), nr = 2, byrow = TRUE)
  )
}

#bugs.data(model.data)
#R2WinBUGS:::bugs.inits(inits, 1, digits = 6)

ptm <- proc.time()
joint.mcmc <- bugs(data = model.data, 
                   inits = inits,
                   model.file = "andrinopoulous2014_model.txt", 
                   parameters = c("betas", "sigma", "gammasR", "gammasD", "alphasR", "alphasD", "D"),
                   n.chains = 2, 
                   n.iter = 350000,
                   n.burnin = 250000,
                   n.thin = 10,
                   debug = FALSE,
                   working.directory = "C:/Users/graeme/Dropbox/Research/Biostatistics/Competing Risks/")
proc.time() - ptm # CPU time

##*********************************************************
## Diagnostics
##*********************************************************

pnorm(abs(geweke.diag(as.mcmc.list(joint.mcmc))[[1]][[1]]), 
      lower.tail = FALSE)*2

plot(as.mcmc.list(joint.mcmc))

library(runjags)
joint.mcmc10 <- combine.mcmc(as.mcmc.list(joint.mcmc), thin = 10,
                             collapse.chains = FALSE)
pnorm(abs(geweke.diag(joint.mcmc10)[[1]][[1]]), 
      lower.tail = FALSE)*2
summary(joint.mcmc10)
traceplot(joint.model10)
autocorr.plot(joint.mcmc10)

##*********************************************************
## JAGS model
##*********************************************************

# library(rjags)
# library(runjags)
# 
# joint.mcmc <- run.jags(
#   model = "C:/Users/graeme/Dropbox/Research/Biostatistics/Competing Risks/andrinopoulous2014_model.txt", 
#   monitor = c("betas", "sigma", "gammasR", "gammasD", "alphasR", "alphasD", "D"),
#   data = model.data, 
#   inits = inits,
#   adapt = 15000,
#   burnin = 30000,
#   sample = 2000,
#   n.chains = 3,
#   thin = 100,
#   modules = "glm",
#   method = "parallel"
# )
# 
# joint.mcmc10 <- combine.mcmc(as.mcmc.list(joint.mcmc), thin = 30,
#                              collapse.chains = FALSE)



