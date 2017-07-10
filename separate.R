##*********************************************************
## Data
##*********************************************************

library(joineR)
library(ggplot2)

epileptic <- read.table("epileptic.txt", header = TRUE)
head(epileptic)

epileptic$time <- epileptic$time / 365.25
epileptic$with.time <- epileptic$with.time / 365.25
epileptic$status <- rep("Censored", nrow(epileptic))
epileptic$status[epileptic$with.status.uae == 1] <- "UAE"
epileptic$status[epileptic$with.status.isc == 1] <- "ISC"

# epileptic$time2 <- do.call("c", 
#   by(epileptic, epileptic$id, function(u) u$time - max(u$time)))

## Profile plots

ggplot(aes(x = time - with.time, y = dose), data = epileptic) +
  geom_line(aes(group = id), colour = "grey", size = 0.8) +
  facet_grid(status ~ treat) +
  #theme_bw() +
  labs(
    x = "Time before treatment failure or censoring (years)", 
    y = "Calibrated dose"
  ) +
  geom_smooth(aes(group = 1), colour = "red", se = FALSE, size = 1.5,
              method = "loess", span = 0.7) +
  theme(text = element_text(size = 16))
ggsave("figure1.pdf", width = 7, height = 7)
ggsave("figure1.png", width = 7, height = 7)

##*********************************************************
## Seperate models
##*********************************************************

## Model data

# Time-to-event data (one row per subject)
survdat <- epileptic[ , c(1, 4, 6, 7, 8)]
nobs <- table(survdat$id)
survdat <- survdat[!duplicated(survdat$id), ]

# Logitudinal data
longdat <- epileptic[ , c(1:3, 8)]

ptm <- proc.time()

## Longitudinal submodel

lmeFit <- lme(dose ~ treat * time,
              random = ~ time | id,
              data = epileptic)
summary(lmeFit)
getVarCov(lmeFit)

fixef(lmeFit) - qnorm(0.975)*sqrt(diag(lmeFit$varFix))
fixef(lmeFit) + qnorm(0.975)*sqrt(diag(lmeFit$varFix))

## Time-to-event submodel

coxFit.uae <- coxph(Surv(with.time, with.status.uae) ~ treat,
                    data = survdat)
summary(coxFit.uae)
confint(coxFit.uae)

coxFit.isc <- coxph(Surv(with.time, with.status.isc) ~ treat,
                    data = survdat)
summary(coxFit.isc)
confint(coxFit.isc)

proc.time() - ptm # CPU time

## Random effects parameterization w/o joint modelling

survdat$b1 <- ranef(lmeFit)[,1] + fixef(lmeFit)[1]
survdat$b2 <- ranef(lmeFit)[,2] + fixef(lmeFit)[3]

coxFit1 <- coxph(Surv(with.time, with.status.isc) ~ treat + b1 + b2,
                    data = survdat)
coxFit2 <- coxph(Surv(with.time, with.status.uae) ~ treat + b1 + b2,
                 data = survdat)
summary(coxFit1)
summary(coxFit2)





