##*********************************************************
## Data
##*********************************************************

library(joineR)

data(epileptic)
head(epileptic)

epileptic$time <- epileptic$time / 365.25
epileptic$with.time <- epileptic$with.time / 365.25
epileptic$treat <- as.numeric(epileptic$treat == "LTG")
epileptic$status <- rep(0, nrow(epileptic))
epileptic$status[epileptic$with.status.uae == 1] <- 1
epileptic$status[epileptic$with.status.isc == 1] <- 2


y <- with(epileptic, data.frame(
  dose,
  rep(1, nrow(epileptic)),
  time,
  rep(1, nrow(epileptic)),
  time,
  treat,
  time*treat
))

c <- with(epileptic, data.frame(
  with.time,
  status,
  treat
))
c <- c[!duplicated(epileptic$id), ]

m <- data.frame(as.numeric(table(epileptic$id)))

write.table(y, "y.txt", sep = "    ", row.names = FALSE, col.names = FALSE)
write.table(c, "c.txt", sep = "    ", row.names = FALSE, col.names = FALSE)
write.table(m, "m.txt", sep = "    ", row.names = FALSE, col.names = FALSE)

##*********************************************************
## Results
##*********************************************************

beta <- c(1.926670, 0.148609, -0.141785, 0.122646)
beta.se <- c(0.050951, 0.029018, 0.068613, 0.035016)

beta - qnorm(0.975)*beta.se
beta + qnorm(0.975)*beta.se

gamma <- c(-0.542853, -0.306389)
gamma.se <- c(0.231458, 0.223083)

gamma - qnorm(0.975)*gamma.se
gamma + qnorm(0.975)*gamma.se

v <- c(1, -1.501823)
v.se <- c(0, 0.224280)

v - qnorm(0.975)*v.se
v + qnorm(0.975)*v.se

sigma2 <- 0.184340
sigma2.se <- 0.002213

sigma2 + c(-1, 1)*qnorm(0.975)*sigma2.se

# In the order of Sigma_11, Sigma_12, Sigma_13, Sigma_22, Sigma_23, Sigma_33
Sigma <- c(0.660305, 0.119664, -0.092769, 0.222904, -0.402497, 0.790844)
# In the order of Sigma_11, Sigma_22, Sigma_33, Sigma_12, Sigma_23, Sigma_13
Sigma.se <- c(0.041392, 0.024250, 0.111955, 0.021178, 0.046978, 0.047834)

Sigma - qnorm(0.975)*Sigma.se[c(1, 4, 6, 2, 5, 3)]
Sigma + qnorm(0.975)*Sigma.se[c(1, 4, 6, 2, 5, 3)]
