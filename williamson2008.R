##*********************************************************
## Estimation functions
##*********************************************************

# FIRST PROGRAM TO FIT EXTENDED W&T TO COMPETING RISK SCENARIO
# ONLY CONSIDER 2 RISKS INITIALLY

# longst <- function(longdat) {
#   
#   list(
#     b1vec = data.frame(fixef(long.start)), 
#     sigma.z = data.frame(noise), 
#     sigma.u = sig, 
#     corr = data.frame(corr), 
#     log.like = data.frame(l1)
#   )
# 
# }

longst <- function(longdat) {

  # Separate longitudinal model to extract initial estimates
  long.start <- lme(dose ~ time + treat + interaction,
                    data = longdat,
                    random = ~ time | id)
  
  #summary(long.start)
  noise <- long.start$sigma^2
  sig <- as.matrix(getVarCov(long.start))
  corr <- as.numeric(VarCorr(long.start)[2, 3])
  l1 <- as.numeric(logLik(long.start))
  
  list(
    b1vec = data.frame(fixef(long.start)), 
    sigma.z = data.frame(noise), 
    sigma.u = sig, 
    corr = data.frame(corr), 
    log.like = data.frame(l1)
  )
  
}

# NOW WANT TWO OF THESE, ONE FOR EACH RISK, GET 2 HAZARDS AND 2 ID
# VECS AS WELL AS 2 VECTORS OF DISTINCT SURVIVAL TIMES. CALL THEM
# survst1 AND survst2

survsta <- function(survdat) {

  n = length(survdat[ , 2])
  surv.time = survdat[ , 2]
  cen = survdat[ , 3]
  #if(cen[1] == 0) {cen[1] = 1}
  p2 = dim(survdat)[2] - 5
  if (p2 == 0) {surv.start = coxph(Surv(surv.time, cen) ~ 0)}
  if(p2 > 0) {
    X2 = as.matrix(survdat[ , 6:dim(survdat)[2]])
    surv.start = coxph(Surv(surv.time, cen) ~ X2)
  }
  alpha.0 = basehaz(surv.start, FALSE)
  l = length(alpha.0[ , 2])
  haz = vector("numeric", l)
  s.dist = vector("numeric", l)
  haz[1] = alpha.0[1, 1]
  haz[2:l] = diff(alpha.0[ , 1])
  s.dist = alpha.0[ , 2]
  id.1 = match(surv.time[1:n], s.dist)
  dummy = match(s.dist[1:l], surv.time)
  d.dummy = diff(dummy)
  id.2 = rep(id.1[dummy[1:(l-1)]], d.dummy[1:(l-1)])
  id.3 = rep(id.1[dummy[l]], (n - dummy[l] + 1))
  id.4 = vector("numeric", n)
  id.4 = c(id.2, id.3)
  id.5 = c(rep(0, match(1, id.1) - 1), id.4)
  l2 = surv.start$loglik
  
  list(
    b2vec.a = data.frame(c(coef(surv.start), 0)), 
    a.0a = data.frame(haz), 
    iden.a = data.frame(id.5), 
    surv.dista = data.frame(s.dist), 
    log.likea = data.frame(l2)
  )

}

survstb <- function(survdat) {

  n = length(survdat[ , 2])
  surv.time = survdat[ , 2]
  cen = survdat[ , 4]
  #if(cen[1] == 0) {cen[1] = 1}
  p2 = dim(survdat)[2] - 5
  if (p2 == 0) {surv.start = coxph(Surv(surv.time, cen) ~ 0)}
  if(p2 > 0) {
    X2 = as.matrix(survdat[ , 6:dim(survdat)[2]])
    surv.start = coxph(Surv(surv.time, cen) ~ X2)
  }
  alpha.0 = basehaz(surv.start, FALSE)
  l = length(alpha.0[ , 2])
  haz = vector("numeric", l)
  s.dist = vector("numeric", l)
  haz[1] = alpha.0[1, 1]
  haz[2:l] = diff(alpha.0[ , 1])
  s.dist = alpha.0[ , 2]
  id.1 = match(surv.time[1:n], s.dist)
  dummy = match(s.dist[1:l], surv.time)
  d.dummy = diff(dummy)
  id.2 = rep(id.1[dummy[1:(l-1)]], d.dummy[1:(l-1)])
  id.3 = rep(id.1[dummy[l]], (n - dummy[l] + 1))
  id.4 = vector("numeric", n)
  id.4 = c(id.2, id.3)
  id.5 = c(rep(0, match(1, id.1) - 1), id.4)
  l2 = surv.start$loglik
  
  list(
    b2vec.b = data.frame(c(coef(surv.start), 0)), 
    a.0b = data.frame(haz), 
    iden.b = data.frame(id.5), 
    surv.distb = data.frame(s.dist), 
    log.likeb = data.frame(l2)
  )

}

em.alg.cr <- function(longdat, survdat, paraests) {
  
  Y = longdat[ , 2]
  lda.time = longdat[ , 3]
  X1 = as.matrix(longdat[ , 4:dim(longdat)[2]])
  n = length(survdat[ , 2])
  surv.time = survdat[ , 2]
  cen.a = survdat[ , 3]
  cen.b = survdat[ , 4]
  n.obs = survdat[ , 5]
  p1 = dim(X1)[2]
  p2 = dim(survdat)[2] - 5
  X2 = 0
  if(p2 > 0) {
    X2 = as.matrix(survdat[ , 6:dim(survdat)[2]])
  }
  b1 = paraests$b1vec[ , 1]
  sig = matrix(0, 2, 2)
  sig[1, 1] = paraests$sigma.u[1, 1]
  sig[2, 2] = paraests$sigma.u[2, 2]
  sig[1, 2] = paraests$sigma.u[1, 2]
  sig[2, 1] = sig[1, 2]
  rho = paraests$corr[ , 1]
  varz = paraests$sigma.z[ , 1]
  
  # NOW NEED SURVIVAL INITIAL HAZARDS, DISTINCT FAILURE VECTORS
  haz.a = paraests$a.0a[ , 1]
  s.dista = paraests$surv.dista[ , 1]
  id.a = paraests$iden.a[ , 1]
  haz.b = paraests$a.0b[ , 1]
  s.distb = paraests$surv.distb[ , 1]
  id.b = paraests$iden.b[ , 1]
  N = sum(n.obs)
  maxn = max(n.obs)
  gpt = 3
  ran = 2
  lat = 1
  b2.a = vector("numeric", p2 + lat)
  b2.b = vector("numeric", p2 + lat)
  b2.a = paraests$b2vec.a[ , 1]
  b2.b = paraests$b2vec.b[ , 1]
  ab = vector("numeric", gpt)
  w = vector("numeric", gpt)
  
  # ENTER ABSCISSAE & WEIGHTS
  ab[1] = 1.22474487
  ab[2] = 0.0
  ab[3] = -ab[1]
  w[1] = 0.295408975
  w[2] = 1.1816359
  w[3] = w[1]
  gammat = matrix(0, gpt^2, ran)
  gammat[ , 1] = rep(ab, each = gpt)
  gammat[ , 2] = rep(ab, gpt)
  wvec = as.vector(w %*% t(w))
  EU = matrix(0, n, 2)
  EUU = matrix(0, n, 3)
  EexpU.a = matrix(0, n, length(haz.a))
  EU0expU.a = matrix(0, n, length(haz.a))
  EU1expU.a = matrix(0, n, length(haz.a))
  EU0U0expU.a = matrix(0, n, length(haz.a))
  EU0U1expU.a = matrix(0, n, length(haz.a))
  EU1U1expU.a = matrix(0, n, length(haz.a))
  EexpU.b = matrix(0, n, length(haz.b))
  EU0expU.b = matrix(0, n, length(haz.b))
  EU1expU.b = matrix(0, n, length(haz.b))
  EU0U0expU.b = matrix(0, n, length(haz.b))
  EU0U1expU.b = matrix(0, n, length(haz.b))
  EU1U1expU.b = matrix(0, n, length(haz.b))
  W11 = diag(maxn)
  W3 = matrix(0, maxn, 2)
  cvar = matrix(0, ran, ran)
  cvarch = matrix(0, ran, ran)
  
  # LOOP PART BEGINS HERE...
  major.it = 25
  minor.it = 10
  iter = 0
  for (it in 1:major.it) {
    for (it.2 in 1:minor.it) {
      iter = iter + 1
      W22 = sig
      b2x = matrix(0, n, 1)
      if(p2 > 0) {
        b2temp.a = c(b2.a[1:p2])
        b2temp.b = c(b2.b[1:p2])
      }
      if(p2 > 0) {
        b2x.a = X2 %*% b2temp.a
        b2x.b = X2 %*% b2temp.b
        b2x = b2x.a + b2x.b
      }
      rlong = Y - (X1 %*% b1)
      count = 1
      for (i in 1:n) {
        W21 = matrix(0, 2, n.obs[i])
        rvec = rlong[count:(count + n.obs[i]-1)]
        W11 = (sig[1, 2] + sig[2, 2] * lda.time[count:(count + n.obs[i]-1)]) %*% t(lda.time[count:(count + n.obs[i]-1)]) + sig[1, 1] + sig[1, 2] * lda.time[count:(count + n.obs[i]-1)]
        W21[1, 1:n.obs[i]] = sig[1, 1] + sig[1, 2] * lda.time[count:(count + n.obs[i]-1)]
        W21[2, 1:n.obs[i]] = sig[1, 2] + sig[2, 2] * lda.time[count:(count + n.obs[i]-1)]
        W11 = W11 + (varz * diag(n.obs[i]))
        count = count + n.obs[i]
        W3 = solve(W11, t(W21))
        cvar = W22-W21 %*% W3
        
        # TRANSFORM TO INDEPENDENT VARIABLES, WHICH WE'LL CALL GAMMA
        cvar = cvar * 2
        cvarch = chol(cvar)
        cvarch = t(cvarch)
        cm = t(W3) %*% rvec
        cmmat = matrix(0, gpt^2, ran)
        cmmat[ , 1] = rep(cm[1], gpt^2)
        cmmat[ , 2] = rep(cm[2], gpt^2)
        newumat = cvarch %*% t(gammat) + t(cmmat)
        fvec = exp(cen.a[i] * b2.a[p2 + lat] * (newumat[1, ] + newumat[2, ] * surv.time[i]) + (cen.b[i] * b2.b[p2 + lat] * (newumat[1, ] + newumat[2, ] * surv.time[i])))
        ssvec.a = 1
        #}
        
        if(id.a[i] > 0) {
          ssvec.a = exp(b2.a[p2 + lat] * (newumat[1, ] + newumat[2, ] %*% t(s.dista[1:id.a[i]]))) %*% haz.a[1:id.a[i]]
        }
        ssvec.b = 1
        if(id.b[i] > 0) {
          ssvec.b = exp(b2.b[p2 + lat] * (newumat[1, ] + newumat[2, ] %*% t(s.distb[1:id.b[i]]))) %*% haz.b[1:id.b[i]]
        }
        ssvec = ssvec.a * ssvec.b * exp(b2x[i, ])
        fvec = fvec * wvec * exp(-ssvec)
        den = sum(fvec)
        EU[i, 1:2] = t(newumat %*% fvec / den)
        EUU[i, 1:2] = t((newumat^2) %*% fvec / den)
        EUU[i, 3] = (newumat[1, ] * newumat[2, ]) %*% fvec / den
        const.a = exp(b2.a[p2 + lat] * (newumat[1, ] + newumat[2, ] %*% t(s.dista[1:id.a[i]])))
        if(id.a[i] == 0) {const.a = const.a^0}
        EexpU.a[i, 1:id.a[i]] = t(fvec) %*% const.a / den
        EU0expU.a[i, 1:id.a[i]] = t(fvec) %*% (newumat[1, ] * const.a) / den
        EU1expU.a[i, 1:id.a[i]] = t(fvec) %*% (newumat[2, ] * const.a) / den
        EU0U0expU.a[i, 1:id.a[i]] = t(fvec) %*% (newumat[1, ]^2 * const.a) / den
        EU0U1expU.a[i, 1:id.a[i]] = t(fvec) %*% (newumat[1, ] * newumat[2, ] * const.a) / den
        EU1U1expU.a[i, 1:id.a[i]] = t(fvec) %*% (newumat[2, ]^2 * const.a) / den
        const.b = exp(b2.b[p2 + lat] * (newumat[1, ] + newumat[2, ] %*% t(s.distb[1:id.b[i]])))
        
        if(id.b[i] == 0) {const.b = const.b^0}
        EexpU.b[i, 1:id.b[i]] = t(fvec) %*% const.b / den
        EU0expU.b[i, 1:id.b[i]] = t(fvec) %*% (newumat[1, ] * const.b) / den
        EU1expU.b[i, 1:id.b[i]] = t(fvec) %*% (newumat[2, ] * const.b) / den
        EU0U0expU.b[i, 1:id.b[i]] = t(fvec) %*% (newumat[1, ]^2 * const.b) / den
        
        EU0U1expU.b[i, 1:id.b[i]] = t(fvec) %*% (newumat[1, ] * newumat[2, ] * const.b) / den
        EU1U1expU.b[i, 1:id.b[i]] = t(fvec) %*% (newumat[2, ]^2 * const.b) / den
      }
      
      # NOW M-STEP...
      paracopy.em <- data.frame(c(b1, b2.a, b2.b, varz, sig, rho))
      # UPDATE FOR BASELINE HAZARDS
      # THERE WILL NOW BE A PAIR OF THESE
      ndist.a = max(id.a)
      ndist.b = max(id.b)
      for (i in 1:(ndist.a-1)) {
        sum3 = sum(exp(b2x.a[match(i, id.a):n]) * (EexpU.a[match(i, id.a):n, i]))
        nfail = sum(cen.a[match(i, id.a):(match(i + 1, id.a)-1)])
        haz.a[i] = nfail / sum3
      }
      sum3 = sum(exp(b2x.a[match(ndist.a, id.a):n]) * (EexpU.a[match(ndist.a, id.a):n, i]))
      nfail = sum(cen.a[match(ndist.a, id.a):n])
      haz.a[ndist.a] = nfail / sum3
      for (i in 1:(ndist.b-1)) {
        sum3 = sum(exp(b2x.b[match(i, id.b):n]) * (EexpU.b[match(i, id.b):n, i]))
        nfail = sum(cen.b[match(i, id.b):(match(i + 1, id.b)-1)])
        haz.b[i] = nfail / sum3
      }
      sum3 = sum(exp(b2x.b[match(ndist.b, id.b):n]) * (EexpU.b[match(ndist.b, id.b):n, i]))
      nfail = sum(cen.b[match(ndist.b, id.b):n])
      haz.b[ndist.b] = nfail / sum3
      
      # NEED SOME MATRICES SET-UP FOR ESTIMATING BETA_1 & LATER BETA_2
      EUmat = matrix(0, N, 2)
      EUUmat = matrix(0, N, 3)
      EUmat[ , 1] = rep(EU[ , 1], n.obs)
      EUmat[ , 2] = rep(EU[ , 2], n.obs)
      EUUmat[ , 1] = rep(EUU[ , 1], n.obs)
      EUUmat[ , 2] = rep(EUU[ , 2], n.obs)
      EUUmat[ , 3] = rep(EUU[ , 3], n.obs)
      summat.a = matrix(0, n, 1)
      summat2.a = matrix(0, n, 2)
      summat3.a = matrix(0, n, 3)
      summat.b = matrix(0, n, 1)
      summat2.b = matrix(0, n, 2)
      summat3.b = matrix(0, n, 3)
      
      for (i in 1:n) {
        if (id.a[i]>0) {
          summat.a[i, 1] = sum(EexpU.a[i, 1:id.a[i]] * haz.a[1:id.a[i]])
          summat2.a[i, 1] = sum(EU0expU.a[i, 1:id.a[i]] * haz.a[1:id.a[i]])
          summat2.a[i, 2] = sum(EU1expU.a[i, 1:id.a[i]] * s.dista[1:id.a[i]] * haz.a[1:id.a[i]])
          summat3.a[i, 1] = sum(EU0U0expU.a[i, 1:id.a[i]] * haz.a[1:id.a[i]])
          summat3.a[i, 2] = sum(EU1U1expU.a[i, 1:id.a[i]] * (s.dista[1:id.a[i]]^2) * haz.a[1:id.a[i]])
          summat3.a[i, 3] = sum(EU0U1expU.a[i, 1:id.a[i]] * s.dista[1:id.a[i]] * haz.a[1:id.a[i]])}
        if (id.b[i]>0) {
          summat.b[i, 1] = sum(EexpU.b[i, 1:id.b[i]] * haz.b[1:id.b[i]])
          summat2.b[i, 1] = sum(EU0expU.b[i, 1:id.b[i]] * haz.b[1:id.b[i]])
          summat2.b[i, 2] = sum(EU1expU.b[i, 1:id.b[i]] * s.distb[1:id.b[i]] * haz.b[1:id.b[i]])
          summat3.b[i, 1] = sum(EU0U0expU.b[i, 1:id.b[i]] * haz.b[1:id.b[i]])
          summat3.b[i, 2] = sum(EU1U1expU.b[i, 1:id.b[i]] * (s.distb[1:id.b[i]]^2) * haz.b[1:id.b[i]])
          summat3.b[i, 3] = sum(EU0U1expU.b[i, 1:id.b[i]] * s.distb[1:id.b[i]] * haz.b[1:id.b[i]])}
      }
      
      # UPDATE BETA_1 PARAMETER VECTOR
      tEUmat = EUmat[ , 2] * lda.time
      sum = EUmat[ , 1] + tEUmat
      Ystar = Y-sum
      XTX = t(X1) %*% X1 
      XTY = t(X1) %*% Ystar
      b1 = solve(XTX, XTY)
      
      # GET UPDATED NOISE TERM
      bx = X1 %*% b1
      r = Y-bx
      sum2 = r^2-2 * r * (EUmat[ , 1] + tEUmat) + EUUmat[ , 1] + (EUUmat[ , 2] * (lda.time^2))
      sum2 = sum2 + 2 * EUUmat[ , 3] * lda.time
      varz = sum(sum2) / N
      
      # RANDOM EFFECTS COVARIANCE MATRIX
      sig[1, 1] = sum(EUU[ , 1]) / n
      sig[2, 2] = sum(EUU[ , 2]) / n
      sig[1, 2] = sum(EUU[ , 3]) / n
      sig[2, 1] = sig[1, 2]
      rho = sig[1, 2] / sqrt(sig[1, 1] * sig[2, 2])
      
      # CALCULATE DERIVATIVES REQUIRED FOR BETA_2 NEWTON-RAPHSON
      fd.a = vector("numeric", p2 + lat)
      sd.a = matrix(0, p2 + lat, p2 + lat)
      eb2x.a = exp(b2x.a)
      fd.b = vector("numeric", p2 + lat)
      sd.b = matrix(0, p2 + lat, p2 + lat)
      eb2x.b = exp(b2x.b)
      if(p2>0) {
        for (i in 1:p2) {
          fd.a[i] = sum(cen.a * X2[ , i])-sum(X2[ , i] * eb2x.a * summat.a[ , 1])
          sd.a[i, p2 + lat] = (-sum(X2[ , i] * eb2x.a * (summat2.a[ , 1] + summat2.a[ , 2])))
          fd.b[i] = sum(cen.b * X2[ , i])-sum(X2[ , i] * eb2x.b * summat.b[ , 1])
          sd.b[i, p2 + lat] = (-sum(X2[ , i] * eb2x.b * (summat2.b[ , 1] + summat2.b[ , 2])))}
      }
      fd.a[p2 + lat] = sum(cen.a * (EU[ , 1] + EU[ , 2] * surv.time))-
        sum(eb2x.a * (summat2.a[ , 1] + summat2.a[ , 2]))
      fd.b[p2 + lat] = sum(cen.b * (EU[ , 1] + EU[ , 2] * surv.time))-
        sum(eb2x.b * (summat2.b[ , 1] + summat2.b[ , 2]))
      sd.a = sd.a + t(sd.a)
      sd.b = sd.b + t(sd.b)
      if(p2>0) {
        for (i in 1:p2) {
          for (j in 1:p2) {
            sd.a[i, j] = (-sum(X2[ , i] * X2[ , j] * eb2x.a * summat.a[ , 1]))
            sd.b[i, j] = (-sum(X2[ , i] * X2[ , j] * eb2x.b * summat.b[ , 1]))}}}
      sd.a[p2 + lat, p2 + lat] = (-sum(eb2x.a * (summat3.a[ , 1] + 2 * summat3.a[ , 3] + summat3.a[ , 2])))
      sd.b[p2 + lat, p2 + lat] = (-sum(eb2x.b * (summat3.b[ , 1] + 2 * summat3.b[ , 3] + summat3.b[ , 2])))
      
      # PERFORM NEWTON-RAPHSON STEP
      b2.a = b2.a-solve(sd.a, fd.a)
      b2.b = b2.b-solve(sd.b, fd.b)
    }
    
    para.em <- data.frame(c(b1, b2.a, b2.b, varz, sig, rho))
    check = sum(abs(paracopy.em-para.em) > 0.001)
    if(check == 0) {break}
    
  }
  
  list(
    b1vec = data.frame(b1), 
    b2vec.a = data.frame(b2.a), 
    b2vec.b = data.frame(b2.b), 
    sigma.z = data.frame(varz), 
    sigma.u = data.frame(sig), 
    corr = data.frame(rho), 
    conv = data.frame(c(iter, check))
  )

}

sortcr.dat <- function(longdat, survdat) {
  # Sorting survival & longitudinal data by survival time
  sort.long = matrix(0, dim(longdat)[1], dim(longdat)[2])
  index = rep(survdat[ , 2], survdat[ , 5])
  sort.long = longdat[order(index), ]
  sort.surv = matrix(0, dim(survdat)[1], dim(survdat)[2])
  sort.surv[ , c(1, 3:(dim(survdat)[2]))] = as.matrix(survdat[order(survdat[ , 2]), c(1, 3:(dim(survdat)[2]))])
  sort.surv[ , 2] = sort(survdat[ , 2])
  list(long.s = data.frame(sort.long), surv.s = data.frame(sort.surv))
}


fitWT.cr <- function(longdat, survdat) {
  # Check data are OK - still to do
  # Eg order(unique(longdat[ , 1])) =  = order(survdat[ , 1])
  # dat.check = sum(as.integer(order(unique(long.dat[ , 1])) =  = order(surv.dat[ , 1])))
  # ifelse(dat.check =  = n,  RUN ALG HERE,  break)
  sort = sortcr.dat(longdat, survdat)
  ldaests = longst(longdat)
  longdat = as.matrix(sort$long.s)
  survdat = as.matrix(sort$surv.s)
  survests.a = survsta(survdat)
  survests.b = survstb(survdat)
  paraests = c(ldaests, survests.a, survests.b)
  em.alg.cr(longdat, survdat, paraests)
} 

##*********************************************************
## Data
##*********************************************************

library(joineR)
data(epileptic)

head(epileptic)
summary(epileptic)

epileptic$time <- epileptic$time / 365.25
epileptic$with.time <- epileptic$with.time / 365.25
epileptic$treat <- as.numeric(epileptic$treat == "LTG")

# Time-to-event data (one row per subject)
survdat <- epileptic[ , c(1, 4, 6, 7, 8)]
nobs <- table(survdat$id)
survdat <- survdat[!duplicated(survdat$id), ]
survdat$nobs <- nobs
survdat <- survdat[ , c(1, 2, 3, 4, 6, 5)]

# Logitudinal data
longdat <- epileptic[ , c(1:3, 8)]
longdat$interaction <- with(longdat, treat*time)
longdat$time.fix <- longdat$time
longdat$intercept <- rep(1, nrow(longdat))
longdat <- longdat[ , c(1:3, 7, 6, 4, 5)]

##*********************************************************
## Fit model
##*********************************************************

library(nlme)

# Separate longitudinal model to extract initial estimates
long.start <- lme(dose ~ time + treat + interaction,
                  data = longdat,
                  random = ~ time | id)
# summary(long.start)
# noise <- long.start$sigma^2
# sig <- as.matrix(getVarCov(long.start))
# corr <- as.numeric(VarCorr(long.start)[2, 3])
# l1 <- as.numeric(logLik(long.start))

ptm <- proc.time()
fit <- fitWT.cr(longdat, survdat)
proc.time() - ptm # CPU time

##*********************************************************
## Bootstrap SEs / 95% CIs
##*********************************************************

joint.stats <- function(d, i) {

  # Generate cluster-sampled epileptic dataset
  N <- length(i)
  out <- lapply(i, function(u) subset(d, d$id == u))
  m <- do.call("c", lapply(out, nrow))
  id <- rep(1:N, m)
  epileptic.boot <- do.call("rbind", out)
  epileptic.boot$id <- id
 
  # Cluster-sampled time-to-event data
  survdat <- epileptic.boot[ , c(1, 4, 6, 7, 8)]
  nobs <- table(survdat$id)
  survdat <- survdat[!duplicated(survdat$id), ]
  survdat$nobs <- nobs
  survdat <- survdat[ , c(1, 2, 3, 4, 6, 5)]

  # Cluster-sampled logitudinal data
  longdat <- epileptic.boot[ , c(1:3, 8)]
  longdat$interaction <- with(longdat, treat*time)
  longdat$time.fix <- longdat$time
  longdat$intercept <- rep(1, nrow(longdat))
  longdat <- longdat[ , c(1:3, 7, 6, 4, 5)]
  
   fit <- fitWT.cr(longdat, survdat)
   return(unlist(fit))

}

library(snow)

ptm <- proc.time()
set.seed(1)
cl <- makeCluster(4, type = "SOCK")
clusterCall(cl, fun = function() {
  library("nlme")
  library("stats")
  library("survival")
})
clusterExport(cl, c("epileptic", "joint.stats", "fitWT.cr", "sortcr.dat", 
                    "em.alg.cr", "longst", "survsta", "survstb")) 
joiner.boot <- boot(epileptic, joint.stats, R = 500, parallel = "snow", ncpus = 4, cl = cl)
stopCluster(cl)
proc.time() - ptm # CPU time

do.call("rbind", lapply(1:16, function(i) boot.ci(joiner.boot, index = i, type = "perc")$percent[1 , 4:5]))

# ## Using clusterApply()
# set.seed(1)
# cl <- makeCluster(4, type = "SOCK")
# clusterCall(cl, fun = function() {
#   library("nlme")
#   library("stats")
#   library("survival")
# })
# clusterExport(cl, c("epileptic", "joint.stats", "fitWT.cr", "sortcr.dat", 
#                     "em.alg.cr", "longst", "survsta", "survstb")) 
# system.time(clusterApply(cl, 1:10, function(i) joint.stats(epileptic, sample(1:605, replace=T))))
# stopCluster(cl)

# ## Mac OSX equivalent
# system.time(mclapply(1:10, function(i) joint.stats(epileptic, sample(1:605, replace = TRUE))))


