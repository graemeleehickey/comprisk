##*********************************************************
## Data
##*********************************************************

library(joineR)
library(nlme)
library(JM)

epileptic <- read.table("epileptic.txt", header = TRUE)
head(epileptic)

epileptic$time <- epileptic$time / 365.25
epileptic$with.time <- epileptic$with.time / 365.25
epileptic$treat <- as.numeric(epileptic$treat == "LTG")

# Time-to-event data (one row per subject)
survdat <- epileptic
survdat <- survdat[!duplicated(survdat$id), ]
survdat$status <- rep("alive", nrow(survdat))
survdat$status[survdat$with.status.uae == 1] <- "UAE"
survdat$status[survdat$with.status.isc == 1] <- "ISC"
survdat <- survdat[ , c(1, 4, 8, 12)]
survdat.CR <- crLong(survdat, statusVar = "status",
                     censLevel = "alive", nameStrata = "CR")

##*********************************************************
## Custom summary function
##*********************************************************

# JM uses contrasts, therefore we need to extract cov-var terms
# in order to get the SEs and effects we are interested in

jm.summary <- function(model.in, D = 0) {
  
  out <- confint(model.in)
  
  vcv <- solve(model.in$Hessian)[6:11, 6:11] # take up to 11 just in case extra alphas
  SE.beta <- sqrt(vcv[1,1] + vcv[2,2] + 2*vcv[1,2])
  SE.alpha <- sqrt(vcv[3,3] + vcv[4,4] + 2*vcv[3,4])
  
  sigma <- c(NA, model.in$coefficients$sigma, NA)
  Sigma11 <- c(NA, model.in$coefficients$D[1, 1], NA)
  Sigma22 <- c(NA, model.in$coefficients$D[2, 2], NA)
  Sigma12 <- c(NA, model.in$coefficients$D[2, 1], NA)  
  
  out[6, 2] <- sum(model.in$coefficients$gammas)
  out[6, 1] <- out[6, 2] - qnorm(0.975)*SE.beta
  out[6, 3] <- out[6, 2] + qnorm(0.975)*SE.beta
  if (D == 0 | D == 2) out[8, 2] <- sum(model.in$coefficients$alpha)
  if (D == 1) out[8, 2] <- sum(model.in$coefficients$Dalpha)
  out[8, 1] <- out[8, 2] - qnorm(0.975)*SE.alpha
  out[8, 3] <- out[8, 2] + qnorm(0.975)*SE.alpha
  
  if (D == 2) {
    SE.alpha2 <- sqrt(vcv[5,5] + vcv[6,6] + 2*vcv[5,6])
    out[10, 2] <- sum(model.in$coefficients$Dalpha)
    out[10, 1] <- out[10, 2] - qnorm(0.975)*SE.alpha2
    out[10, 3] <- out[10, 2] + qnorm(0.975)*SE.alpha2
  }
  
  out <- t(out)
  out <- cbind(out, sigma, Sigma11, Sigma22, Sigma12)
  out <- round(out, 3)
  ci <- paste("(", paste(out[1, ], out[3, ], sep = ", "), ")", sep = "")
  out <- rbind(out[2, ], ci)
  
  if (ncol(out) == 14) out <- out[ , c(1, 3, 2, 4, 11:14, 5:10)]
  else out <- out[ , c(1, 3, 2, 4, 9:12, 5:8)]
  
  return(out)
  
}

##*********************************************************
## Fit models
##*********************************************************

## Separate models

# Longitudinal submodel
lmeFit <- lme(dose ~ treat * time,
              random = ~ time | id,
              data = epileptic)
summary(lmeFit)

# Time-to-event submodel
coxFit <- coxph(Surv(with.time, status2) ~ treat*CR + strata(CR),
                data = survdat.CR, 
                x = TRUE)

summary(coxFit)

# Current value parameterization
ptm <- proc.time()
jointFit1 <- jointModel(lmeFit, coxFit,
                        timeVar = "time",
                        method = "spline-PH-aGH",
                        CompRisk = TRUE,
                        interFact = list(value = ~ CR, data = survdat.CR),
                        #GHk = 9,
                        #GKk = 15,
                        #iter.qN = 500,
                        #numeriDeriv = "cd",
                        verbose = FALSE)

summary(jointFit1)
proc.time() - ptm # CPU time

# Time-dependent slopes parameterization
ptm <- proc.time()
dform <- list(
  fixed = ~ 1 + treat, indFixed = 3:4,
  random = ~ 1, indRandom = 2)

jointFit2 <- jointModel(lmeFit, coxFit,
                        timeVar = "time",
                        method = "spline-PH-aGH",
                        CompRisk = TRUE,
                        interFact = list(value = ~ CR, slope = ~ CR, 
                                         data = survdat.CR),
                        #GHk = 9,
                        #GKk = 15,
                        #iter.qN = 500,
                        #numeriDeriv = "cd",
                        parameterization = "both",
                        derivForm = dform,
                        verbose = FALSE)
summary(jointFit2)
proc.time() - ptm # CPU time

anova(jointFit1, jointFit2)

# Lagged effects (6-months)
ptm <- proc.time()
jointFit3 <- update(jointFit1, lag = 0.5)
summary(jointFit3)
proc.time() - ptm # CPU time

# Cummulative effects parameterization
ptm <- proc.time()
iform <- list(
  fixed = ~ -1 + time + I(treat * time) + 
    I(time^2 / 2) + I(treat * time^2 / 2), 
  indFixed = 1:4,
  random = ~ -1 + time + I(time^2 / 2), 
  indRandom = 1:2)

jointFit4 <- jointModel(lmeFit, coxFit,
                        timeVar = "time",
                        method = "spline-PH-aGH",
                        CompRisk = TRUE,
                        interFact = list(slope = ~ CR, 
                                         data = survdat.CR),
                        #GHk = 9,
                        #GKk = 15,
                        #iter.qN = 500,
                        #numeriDeriv = "cd",
                        parameterization = "slope",
                        derivForm = iform,
                        verbose = FALSE)
summary(jointFit4)
proc.time() - ptm # CPU time

# Weighted-cummulative effects parameterization
ptm <- proc.time()
g <- function(u, pow = 0) {
  f <- function(t) integrate(function(s) s^pow * dnorm(t - s), 0, t)$value
  sapply(u, f)
}

iformW <- list(
  fixed = ~ -1 + I(g(time)) + I(treat * g(time)) + 
    I(g(time, 1)) + I(treat * g(time, 1)), 
  indFixed = 1:4,
  random = ~ -1 + I(g(time)) + I(g(time, 1)), 
  indRandom = 1:2)

jointFit5 <- jointModel(lmeFit, coxFit,
                        timeVar = "time",
                        method = "spline-PH-aGH",
                        CompRisk = TRUE,
                        interFact = list(slope = ~ CR, 
                                         data = survdat.CR),
                        #GHk = 9,
                        #GKk = 15,
                        #iter.qN = 500,
                        #numeriDeriv = "cd",
                        parameterization = "slope",
                        derivForm = iformW,
                        verbose = FALSE)
summary(jointFit5)
proc.time() - ptm # CPU time

# Random effects parameterization (fixed + random slopes only)

ptm <- proc.time()

dform2 <- list(
  fixed = ~ 1, indFixed = 3,
  random = ~ 1, indRandom = 2)

# The re-scaling by 0.45 is the SD from model without the scales
# dform2 <- list(
#   fixed = ~ -1 + I(rep(1/0.45, length(treat))), indFixed = 3,
#   random = ~ -1 +  I(rep(1/0.45, length(treat))), indRandom = 2)

jointFit6 <- jointModel(lmeFit, coxFit,
                        timeVar = "time",
                        method = "spline-PH-aGH",
                        CompRisk = TRUE,
                        interFact = list(slope = ~ CR, 
                                         data = survdat.CR),
                        #GHk = 9,
                        #GKk = 15,
                        #iter.qN = 500,
                        #numeriDeriv = "cd",
                        parameterization = "slope",
                        derivForm = dform2,
                        verbose = FALSE)
summary(jointFit6)
proc.time() - ptm # CPU time

ests <- rbind(
  jm.summary(jointFit1, D = 0),
  jm.summary(jointFit3, D = 0),
  jm.summary(jointFit4, D = 1),
  jm.summary(jointFit5, D = 1),
  jm.summary(jointFit6, D = 1)
)
write.csv(ests, "Riz2012_models13456.csv", row.names = FALSE)
write.csv(jm.summary(jointFit2, D = 2), "Riz2012_models2.csv", row.names = FALSE)

##*********************************************************
## Diagnostics
##*********************************************************

## Longitudinal

pdf("resid_plot_JM.pdf", width = 9, height = 5)
#png("resid_plot.png", width = 9, height = 5, units = "in", res = 300)
par(mfrow = c(1, 2))
plot(jointFit1)
dev.off()

## Event time

par(mfrow = c(1, 2))

# This doesn't work with 'process = "Event' (p. 151 of JM book) -- suspect not Martingale resids
martRes <- resid(jointFit1, type = "Martingale") 
# With 'type = "EventTime"', the lengths of mi.t and martRes are different!
mi.t <- fitted(jointFit1, process = "Longitudinal")
plot(mi.t, martRes)
lines(lowess(y = martRes, x = mi.t, iter = 0), col = 2)

plot(epileptic$dose, martRes)
lines(lowess(y = martRes, x = epileptic$dose, iter = 0), col = 2)

