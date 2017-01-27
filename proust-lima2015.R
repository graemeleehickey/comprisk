##*********************************************************
## Data
##*********************************************************

library(joineR)
library(nlme)
library(lcmm)

data(epileptic)
head(epileptic)

epileptic$time <- epileptic$time / 365.25
epileptic$with.time <- epileptic$with.time / 365.25
epileptic$treat <- as.numeric(epileptic$treat == "LTG")

epileptic$status <- rep(0, nrow(epileptic))
epileptic$status[epileptic$with.status.uae == 1] <- 1
epileptic$status[epileptic$with.status.isc == 1] <- 2

##*********************************************************
## Weibull baseline hazard
##*********************************************************

jlcmFit1 <- Jointlcmm(
  fixed = dose ~ treat * time,
  random = ~ time,
  subject = "id",
  survival = Surv(with.time, status) ~ cause(treat),
  hazard = "Weibull",
  hazardtype = "PH",
  ng = 1,
  logscale = TRUE,
  data = epileptic)

jlcmFit2 <- Jointlcmm(
  fixed = dose ~ treat * time,
  random = ~ time,
  mixture = ~ time + treat:time,
  subject = "id",
  survival = Surv(with.time, status) ~ cause(treat),
  hazard = "Weibull",
  hazardtype = "PH",
  ng = 2,
  #partialH = FALSE,
  logscale = TRUE,  
  data = epileptic)

jlcmFit3 <- Jointlcmm(
  fixed = dose ~ treat * time,
  random = ~ time,
  mixture = ~ time + treat:time,
  subject = "id",
  survival = Surv(with.time, status) ~ cause(treat),
  hazard = "Weibull",
  hazardtype = "PH",
  ng = 3,
  #partialH = FALSE,
  logscale = TRUE,  
  data = epileptic)

jlcmFit4 <- Jointlcmm(
  fixed = dose ~ treat * time,
  random = ~ time,
  mixture = ~ time + treat:time,
  subject = "id",
  survival = Surv(with.time, status) ~ cause(treat),
  hazard = "Weibull",
  hazardtype = "PH",
  ng = 4,
  #partialH = FALSE,
  logscale = TRUE,  
  data = epileptic)

ptm <- proc.time()

jlcmFit4.gs <- gridsearch(
  m = Jointlcmm(
    fixed = dose ~ treat * time,
    random = ~ time,
    mixture = ~ time + treat:time,
    subject = "id",
    survival = Surv(with.time, status) ~ cause(treat),
    hazard = "Weibull",
    hazardtype = "PH",
    ng = 4,
    logscale = TRUE,    
    data = epileptic),
  rep = 20,
  maxiter = 25,
  minit = jlcmFit1
)

time4 <- proc.time() - ptm # CPU time
time4

jlcmFit5 <- Jointlcmm(
  fixed = dose ~ treat * time,
  random = ~ time,
  mixture = ~ time + treat:time,
  subject = "id",
  survival = Surv(with.time, status) ~ cause(treat),
  hazard = "Weibull",
  hazardtype = "PH",
  ng = 5,
  #partialH = FALSE,
  logscale = TRUE,  
  data = epileptic)

ptm <- proc.time()

jlcmFit5.gs <- gridsearch(
  m = Jointlcmm(
    fixed = dose ~ treat * time,
    random = ~ time,
    mixture = ~ time + treat:time,
    subject = "id",
    survival = Surv(with.time, status) ~ cause(treat),
    hazard = "Weibull",
    hazardtype = "PH",
    ng = 5,
    logscale = TRUE,    
    data = epileptic),
  rep = 40,
  maxiter = 35,
  minit = jlcmFit1
)

time5 <- proc.time() - ptm # CPU time
time5

jlcmFit6 <- Jointlcmm(
  fixed = dose ~ treat * time,
  random = ~ time,
  mixture = ~ time + treat:time,
  subject = "id",
  survival = Surv(with.time, status) ~ cause(treat),
  hazard = "Weibull",
  hazardtype = "PH",
  ng = 6,
  #partialH = TRUE,
  logscale = TRUE,    
  data = epileptic)

ptm <- proc.time()

jlcmFit6.gs <- gridsearch(
  m = Jointlcmm(
    fixed = dose ~ treat * time,
    random = ~ time,
    mixture = ~ time + treat:time,
    subject = "id",
    survival = Surv(with.time, status) ~ cause(treat),
    hazard = "Weibull",
    hazardtype = "PH",
    ng = 6,
    logscale = TRUE,    
    data = epileptic),
  rep = 20,
  maxiter = 25,
  minit = jlcmFit1
)

time6 <- proc.time() - ptm # CPU time
time6

##*********************************************************
## Model assessments
##*********************************************************

summarytable(jlcmFit1, jlcmFit2, jlcmFit3, jlcmFit4, jlcmFit4.gs, 
             jlcmFit5, jlcmFit5.gs, jlcmFit6, jlcmFit6.gs)

summary(jlcmFit5) # Best on BIC
confint(jlcmFit5)

bhat <- jlcmFit5$best
id <- 1:length(bhat)
round(bhat - qnorm(0.975)*sqrt(jlcmFit5$V[id*(id+1)/2]), 3)
round(bhat + qnorm(0.975)*sqrt(jlcmFit5$V[id*(id+1)/2]), 3)

postprob(jlcmFit5)

# Confidence interval on the cholesky parameters -> variances
mat <- matrix(c(jlcmFit5$cholesky[1:2], 0, jlcmFit5$cholesky[3]), nr = 2)
mat %*% t(mat)
mat.low <- matrix(c(-0.58735216, -0.11084451, 0, 0.12906991), nr = 2)
mat.upp <- matrix(c(-0.478004253, -0.007790368, 0, 0.199695505), nr = 2)
mat.low %*% t(mat.low)
mat.upp %*% t(mat.upp)

summary(jlcmFit6.gs) # First to pass C.I. test
confint(jlcmFit6.gs)

## Residual plots

pdf("resid_plot_lcmm.pdf", width = 11, height = 11)
#png("resid_plot_lcmm.png", width = 11, height = 11, units = "in", res = 300)
plot(jlcmFit5)
dev.off()

## Cumulative incidence plot
pdf(file = "JLCM (time-to-event).pdf", width = 8, height = 8)
par(mfrow = c(2, 1))
plot(cuminc(jlcmFit5, time = seq(0, 5, 0.1), treat = c(0, 1)), profil = 1, event = 1,
     main = "UAE", bty = "n", ylim = c(0, 1), xlim = c(0, 6), lwd = 2,
     legend.loc = "topright",
     xlab = "Time from randomization (years)",
     ylab = "Cumulative incidence of UAE")
plot(cuminc(jlcmFit5, time = seq(0, 5, 0.1), treat = c(0, 1)), profil = 2, event = 1, add = TRUE, lty = 2,
     legend = NULL, bty = "l", ylim = c(0, 1), lwd = 2)
legend("bottomright", c("CBZ", "LTG"), lty = 1:2, lwd = 2, bty = "n", inset = c(0.02, 0.02))
plot(cuminc(jlcmFit5, time = seq(0, 5, 0.1), treat = c(0, 1)), profil = 1, event = 2,
     main = "ISC", bty = "l", ylim = c(0, 1), xlim = c(0, 6), lwd = 2,
     legend.loc = "topright",
     xlab = "Time from randomization (years)",
     ylab = "Cumulative incidence of ISC")
plot(cuminc(jlcmFit5, time = seq(0, 5, 0.1), treat = c(0, 1)), profil = 2, event = 2, add = TRUE, lty = 2,
     legend = NULL, bty = "n", ylim = c(0, 1), lwd = 2)
legend("bottomright", c("CBZ", "LTG"), lty = 1:2, lwd = 2, bty = "n", inset = c(0.02, 0.02))
dev.off()

## Hazard plot
pdf(file = "JLCM (time-to-event) - hazards.pdf", width = 9, height = 5)
par(mfrow = c(1, 2))
plot(jlcmFit5, which = "hazard",
     bty = "n", lwd = 2, event = 1, legend = NULL,
     main = "Baseline hazards: UAE")
plot(jlcmFit5, which = "hazard",
     bty = "n", lwd = 2, event = 2,
     main = "Baseline hazards: ISC")
dev.off()

## Longitudinal plot
pdf(file = "JLCM (longitudinal).pdf", width = 6, height = 6)
newdata0 <- data.frame("time" = seq(0, 5, 0.1), 2)
newdata0$treat <- 0
newdata1 <- data.frame("time" = seq(0, 5, 0.1), 2)
newdata1$treat <- 1
jlcmFit5.fit0 <- predictY(jlcmFit5, 
                          newdata = newdata0, 
                          var.time = "time")
jlcmFit5.fit1 <- predictY(jlcmFit5, 
                           newdata = newdata1, 
                           var.time = "time")
plot(jlcmFit5.fit0, bty = "l", col = 1:5, lty = 1, lwd = 2,
     legend.loc = "topleft",
     xlab = "Time from randomization (years)",
     ylab = "Calibrated dose")
for (i in 1:5) lines(jlcmFit5.fit1$times$time, jlcmFit5.fit1$pred[ , i], col = i, lty = 2, lwd = 2)
legend(x = 1.5, y = 14, c("CBZ", "LTG"), lwd = 2, lty = c(1, 2), bty = "n")
dev.off()

## Investigation of negative slope for class 2
i <- jlcmFit5$pprob$id[jlcmFit5$pprob$class == "2"]
subset(epileptic, id %in% i)

##*********************************************************
## Piecewise constant baseline hazard (5-knots; quantiles)
##*********************************************************

jlcmFit5pw <- Jointlcmm(
  fixed = dose ~ treat * time,
  random = ~ time,
  mixture = ~ time + treat:time,
  subject = "id",
  survival = Surv(with.time, status) ~ cause(treat),
  hazard = "5-quant-piecewise",
  hazardtype = "PH",
  ng = 5,
  logscale = TRUE,    
  data = epileptic)

summary(jlcmFit5pw)

# ptm <- proc.time()
# 
# jlcmFit1pw <- Jointlcmm(
#   fixed = dose ~ treat * time,
#   random = ~ time,
#   subject = "id",
#   survival = Surv(with.time, status) ~ cause(treat),
#   hazard = "5-quant-piecewise",
#   hazardtype = "PH",
#   ng = 1,
#   #partialH = FALSE,
#   logscale = TRUE,  
#   data = epileptic)
# 
# jlcmFit5pw.gs <- gridsearch(
#   m = Jointlcmm(
#     fixed = dose ~ treat * time,
#     random = ~ time,
#     mixture = ~ time + treat:time,
#     subject = "id",
#     survival = Surv(with.time, status) ~ cause(treat),
#     hazard = "5-quant-piecewise",
#     hazardtype = "PH",
#     ng = 5,
#     logscale = TRUE,    
#     data = epileptic),
#   rep = 20,
#   maxiter = 25,
#   minit = jlcmFit1pw
# )
# 
# time5pw <- proc.time() - ptm # CPU time
# time5pw
# 
# summary(jlcmFit5pw.gs)

bhat <- jlcmFit5pw$best
id <- 1:length(bhat)
round(bhat - qnorm(0.975)*sqrt(jlcmFit5pw$V[id*(id+1)/2]), 3)
round(bhat + qnorm(0.975)*sqrt(jlcmFit5pw$V[id*(id+1)/2]), 3)

##*********************************************************
## Cubic spline baseline hazard
##*********************************************************

# Could not achieve convergence!

jlcmFit5s <- Jointlcmm(
  fixed = dose ~ treat * time,
  random = ~ time,
  mixture = ~ time + treat:time,
  subject = "id",
  survival = Surv(with.time, status) ~ cause(treat),
  hazard = "3-quant-splines",
  hazardtype = "PH",
  ng = 5,
  partialH = FALSE,
  logscale = TRUE,  
  data = epileptic)

ptm <- proc.time()

jlcmFit1s <- Jointlcmm(
  fixed = dose ~ treat * time,
  random = ~ time,
  subject = "id",
  survival = Surv(with.time, status) ~ cause(treat),
  hazard = "3-quant-splines",
  hazardtype = "PH",
  ng = 1,
  partialH = FALSE,
  logscale = TRUE,  
  data = epileptic)

outs <- gridsearch(
    m = Jointlcmm(
      fixed = dose ~ treat * time,
      random = ~ time,
      mixture = ~ time + treat:time,
      subject = "id",
      survival = Surv(with.time, status) ~ cause(treat),
      hazard = "3-quant-splines",
      hazardtype = "PH",
      ng = 3,
      logscale = TRUE,    
      data = epileptic),
    rep = 20,
    maxiter = 50,
    minit = jlcmFit1s
  )

time3s <- proc.time() - ptm # CPU time
time3s


