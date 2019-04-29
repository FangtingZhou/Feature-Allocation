library(snow)
library(e1071)
library(doSNOW)
library(truncnorm)

load('~/Desktop/Advanced Applied Statistics/project/Microbial-Data.RData')

X0 = as.matrix(seedLev2Counts)
Y0 = t(t(pmax(X0, 1)) / colSums(pmax(X0, 1)))

UpdateA = function(Z, A, B, W, zeta, rho, sigmaw = 0.1, K = 20, M = 1) {
  Azero = Aone = Anew = A
  for(i in 1 : nrow(Anew)) {
    index = NULL
    for(k in 1 : ncol(Anew)) {
      if(sum(Anew[-i, k]) != 0) {
        Aone[i, k] = 1
        Azero[i, k] = 0
        postone = log(sum(Anew[-i, k]) + M / K)
        postzero = log(nrow(Anew) - sum(Anew[-i, k]))
        for(l in 1 : nrow(B)) {
          oddsone = Aone[i, ] %*% diag(W[l, ], ncol(W), ncol(W)) %*% B[l, ] + zeta[l]
          postone = oddsone * (Z[i, l] == 1) - log(exp(oddsone) + 1) + postone
          oddszero = Azero[i, ] %*% diag(W[l, ], ncol(W), ncol(W)) %*% B[l, ] + zeta[l]
          postzero = oddszero * (Z[i, l] == 1) - log(exp(oddszero) + 1) + postzero
        }
        if(abs(postone - postzero) > 1e2) Anew[i, k] = as.numeric(postzero < postone) else {
          post = (postone + postzero) / 2
          Anew[i, k] = sample(0 : 1, 1, prob = c(exp(postzero - post), exp(postone-post)))
        }
        Azero = Aone = Anew
      } else {
        index = c(index, k)
      }
    }
    if((length(index) != 0) & (ncol(Anew) - length(index) > 1)) {
      Anew = as.matrix(Anew[, -index])
      Azero = Aone = Anew
      B = as.matrix(B[, -index])
      W = as.matrix(W[, -index])
    }
    if(ncol(Anew) < K) {
      kstar = rpois(1, M / nrow(Anew))
      if(kstar == 0) Azero = Aone = Anew
      if(kstar != 0) {
        anew = matrix(0, nrow(Anew), kstar)
        anew[i, ] = rep(1, kstar)
        anew = cbind(Anew, anew)
        bnew = matrix(runif(nrow(B) * kstar) > rho, nrow(B), kstar)
        bnew = cbind(B, bnew)
        wnew = matrix(rgamma(nrow(W) * kstar, shape = 1, sigmaw), nrow(W), kstar)
        wnew = cbind(W, wnew)
        num = den = 0
        for(l in 1 : nrow(B)) {
          oddsnum = anew[i, ] %*% diag(wnew[l, ], ncol(wnew), ncol(wnew)) %*% bnew[l, ] + zeta[l]
          num = oddsnum * (Z[i, l] == 1) - log(exp(oddsnum) + 1) + num
          oddsden = Anew[i, ] %*% diag(W[l, ], ncol(W), ncol(W)) %*% B[l, ] + zeta[l]
          den = oddsden * (Z[i, l] == 1) - log(exp(oddsden) + 1) + den
        }
        if(runif(1) <= min(1, exp(num - den))) {
          Anew = anew
          Azero = Aone = Anew
          B = bnew
          W = wnew
        } else Azero = Aone = Anew
      }
    }
  }
  return(list(A = Anew, B = B, W = W))
}

UpdateZ = function(Ga, A, B, W, zeta, theta) {
  Z = matrix(NA, nrow(Ga), ncol(Ga))
  for(i in 1 : nrow(Ga)) {
    for(l in 1 : ncol(Ga)) {
      odds = A[i, ] %*% diag(W[l, ], ncol(W), ncol(W)) %*% B[l, ] + zeta[l]
      postone = odds + dgamma(Ga[i, l], theta[i] + 1, log = TRUE)
      postzero = dgamma(Ga[i, l], 1 / (theta[i] + 1), log = TRUE)
      if(abs(postone - postzero) > 1e2) Z[i, l] = as.numeric(postzero < postone) else {
        post = (postone + postzero) / 2
        Z[i, l] = sample(0 : 1, 1, prob = c(exp(postzero - post), exp(postone - post)))
      }
    }
  }
  return(Z)
}

UpdateB = function(Z, A, B, W, zeta, rho) {
  Bzero = Bone = Bnew = B
  for(l in 1 : nrow(B)) {
    for (k in 1 : ncol(B)) {
      Bone[l, k] = 1
      Bzero[l, k] = 0
      oddsone = A %*% diag(W[l, ], ncol(W), ncol(W)) %*% Bone[l, ] + zeta[l]
      postone = log(rho) + sum(oddsone * (Z[, l] == 1)) - sum(log(exp(oddsone) + 1))
      oddszero = A %*% diag(W[l, ], ncol(W), ncol(W)) %*% Bzero[l, ] + zeta[l]
      postzero = log(1 - rho) + sum(oddszero * (Z[, l] == 1)) - sum(log(exp(oddszero) + 1))
      if(abs(postone - postzero) > 1e2) Bnew[l, k] = as.numeric(postzero < postone) else {
        post = (postzero + postone) / 2
        Bnew[l, k] = sample(0 : 1, 1, prob = c(exp(postzero - post), exp(postone - post)))
      }
      Bzero = Bone = Bnew
    }
  }
  return(Bnew)
}

UpdateW = function(Z, A, B, W, zeta, sigma = 1, sigmaw = 0.1) {
  Wnew = W
  for(l in 1: nrow(W)) {
    for(k in 1 : ncol(W)) {
      Wnew[l, k] = rtruncnorm(1, 0, Inf, W[l, k], sigma)
      oddsold = A %*% diag(W[l, ], ncol(W), ncol(W)) %*% B[l, ] + zeta[l]
      den = sum(oddsold * (Z[, l] == 1)) - sum(log(exp(oddsold) + 1)) + dgamma(W[l, k], 1, sigmaw, log = TRUE) + 
        log(dtruncnorm(Wnew[l, k], 0, Inf, W[l, k], sigma))
      oddsnew = A %*% diag(Wnew[l, ], ncol(Wnew), ncol(Wnew)) %*% B[l, ] + zeta[l]
      num = sum(oddsnew * (Z[, l] == 1)) - sum(log(exp(oddsnew) + 1)) + dgamma(Wnew[l, k], 1, sigmaw, log = TRUE) + 
        log(dtruncnorm(W[l, k], 0, Inf, Wnew[l, k], sigma))
      if(runif(1) <= min(1, exp(num - den))) W[l, k] = Wnew[l, k] else
        Wnew[l, k] = W[l, k]
    }
  }
  return(Wnew)
}

UpdateZeta = function(Z, A, B, W, zeta, sigma = 1, sigmaz = 100) {
  zetanew = rep(NA, length(zeta))
  for(l in 1: length(zeta)) {
    z = rnorm(1, zeta[l], sigma)
    oddsold = A %*% diag(W[l, ], ncol(W), ncol(W)) %*% B[l, ] + zeta[l]
    den = sum(oddsold * (Z[, l] == 1)) - sum(log(exp(oddsold) + 1)) + dnorm(zeta[l], 0, sigmaz, log = TRUE)
    oddsnew = A %*% diag(W[l, ], ncol(W), ncol(W)) %*% B[l, ] + z
    num = sum(oddsnew * (Z[, l] == 1)) - sum(log(exp(oddsnew) + 1)) + dnorm(z, 0, sigmaz, log = TRUE)
    if(runif(1) <= min(1, exp(num - den))) zetanew[l] = z else
      zetanew[l] = zeta[l]
  }
  return(zetanew)
}

UpdateRho = function(B, a = 1, b = 1) {
  anew = a + sum(B == 1)
  bnew = b + sum(B == 0)
  rho = rbeta(1, anew, bnew)
  return(rho)
}

UpdateTheta = function(Z, theta, Ga, a = 1, b = 0.1, sigma = 1) {
  for(i in 1 : length(theta)) {
    thetanew = rtruncnorm(1, 0, Inf, theta[i], sigma)
    den = log(dtruncnorm(thetanew, 0, Inf, theta[i], sigma)) + dgamma(theta[i], a, b, log = TRUE) +
      sum(dgamma(Ga[i, Z[i, ] == 0], 1 / (theta[i] + 1), log = TRUE)) + 
      sum(dgamma(Ga[i, Z[i, ] == 1], theta[i] + 1, log = TRUE))
    num = log(dtruncnorm(theta[i], 0, Inf, thetanew, sigma)) + dgamma(thetanew, a, b, log = TRUE) +
      sum(dgamma(Ga[i, Z[i, ] == 0], 1 / (thetanew + 1), log = TRUE)) + 
      sum(dgamma(Ga[i, Z[i, ] == 1], thetanew + 1, log = TRUE))
    if(runif(1) <= min(1, exp(num - den))) theta[i] = thetanew else theta[i] = theta[i]
  }
  return(theta)
}

UpdateGa = function(X, Z, theta, Ga, sigma = 1) {
  Ganew = Ga
  for(i in 1 : nrow(X)) {
    for(l in 1 : ncol(X)) {
      Ganew[i, l] = rtruncnorm(1, 0, Inf, Ga[i, l], sigma)
      para = ifelse(Z[i, l] == 0, 1 / (1 + theta[i]), 1 + theta[i])
      den = dmultinom(X[, l], prob = Ga[, l] / sum(Ga[, l]), log = TRUE) + dgamma(Ga[i, l], para, log = TRUE) +
        log(dtruncnorm(Ganew[i, l], 0, Inf, Ga[i, l], sigma))
      num = dmultinom(X[, l], prob = Ganew[, l] / sum(Ganew[, l]), log = TRUE) + dgamma(Ganew[i, l], para, log = TRUE) +
        log(dtruncnorm(Ga[i, l], 0, Inf, Ganew[i, l], sigma))
      if(runif(1) <= min(1, exp(num - den))) Ga[i, l] = Ganew[i, l] else Ganew[i, l] = Ga[i, l]
    }
  }
  return(Ga)
}

UpdateMCMC = function(maxit, X0, Y0, seed) {
  set.seed(seed)
  ## Initialize parameters
  rho = runif(1)
  A = matrix(as.numeric(runif(nrow(X0)) > 0.5), nrow(X0), 1)
  B = matrix(as.numeric(runif(ncol(X0)) > rho), ncol(X0), 1)
  W = matrix(rgamma(ncol(X0), shape = 1), ncol(X0), 1)
  zeta = rnorm(ncol(X0))
  Z = matrix(NA, nrow(Y0), ncol(Y0))
  theta = rep(NA, nrow(Y0))
  for(i in 1 : nrow(Y0)) {
    Z[i, ] = as.numeric(Y0[i, ] > median(Y0[i, ]))
    theta[i] = mean(Y0[i, Y0[i, ] > median(Y0[i, ])]) / mean(Y0[i, Y0[i, ] <= median(Y0[i, ])])
    theta[i] = theta[i] - 1
  }
  Ga = matrix(NA, nrow(Y0), ncol(Y0))
  for(i in 1 : nrow(Y0)) {
    Ga[i, Z[i, ] == 1] = rgamma(sum(Z[i, ] == 1), theta[i] + 1)
    Ga[i, Z[i, ] == 0] = rgamma(sum(Z[i, ] == 0), 1 / (theta[i] + 1))
  }
  ## record the trace
  recordA = list(A)
  recordB = list(B)
  recordZ = list(Z)
  ## MCMC iteration
  iter = 1
  while(iter < maxit) {
    Z = UpdateZ(Ga, A, B, W, zeta, theta)
    Ga = UpdateGa(X0, Z, theta, Ga)
    theta = UpdateTheta(Z, theta, Ga)
    B = UpdateB(Z, A, B, W, zeta, rho)
    rho = UpdateRho(B)
    W = UpdateW(Z, A, B, W, zeta)
    zeta = UpdateZeta(Z, A, B, W, zeta)
    result = UpdateA(Z, A, B, W, zeta, rho)
    A = result$A
    B = result$B
    W = result$W
    ## record the result
    recordA = c(recordA, list(A))
    recordB = c(recordB, list(B))
    recordZ = c(recordZ, list(Z))
    iter = iter + 1
  }
  return(list(recordA = recordA, recordB = recordB, recordZ = recordZ))
}

Result = UpdateMCMC(20000, X0, Y0, seed = 123)

A = Result$recordA[seq(10000, 20000, 5)]
B = Result$recordB[seq(10000, 20000, 5)]
Z = Result$recordZ[seq(10000, 20000, 5)]

K = NULL
for(i in 1 : length(A)) {
  index = NULL
  for(j in 1 : ncol(A[[i]])) {
    if(sum(A[[i]][, j]) == 1) index = c(index, j)
  }
  if(length(index) != 0) {
    A[[i]] = A[[i]][, -index]
    B[[i]] = B[[i]][, -index]
  }
  K = c(K, ncol(A[[i]]))
}

prob = table(K) / length(K)

ggplot() + geom_histogram(aes(x = K), binwidth = 1, fill = "white", color = "black")


ResultA = NULL
for(i in 1 : length(A)) {
  if(ncol(A[[i]]) == Khat) ResultA = list(ResultA, A[[i]])
}

ResultA = UpdateMCMC(20000, X0, Y0, seed = 456)
ResultB = UpdateMCMC(20000, X0, Y0, seed = 789)

FeatureNumA = ResultA$Num
FeatureNumB = ResultB$Num

gelman.diag(list(as.mcmc(FeatureNumA), as.mcmc(FeatureNumB)))

calHamming = function(MA, MB) {
  A = as.vector(MA)
  ## list all permutations
  permutation = permutations(ncol(MA))
  
  DHamming = rep(NA, nrow(permutation))
  for(l in 1 : nrow(permutation)) {
    B = as.vector(MB[, permutation[l, ]])
    DHamming[l] = hamming.distance(A, B)
  }
  
  result = min(DHamming)
  return(result)
}

dist = matrix(0, length(ResultA), length(ResultA))

numcore = parallel::detectCores()
cluster = makeCluster(numcore)
registerDoSNOW(cluster)
for(i in 1 : length(ResultA)) {
  x = foreach(j = i : length(ResultA), .combine = "c") %dopar% {library(e1071); 
    calHamming(ResultA[[i]], ResultA[[j]])}
  dist[i, i : length(ResultA)] = x
  dist[i : length(ResultA), i] = x
}
stop(cluster)

Ahat = ResultA[[which.min(rowMeans(dist))]]
Baht = ResultB[[which.min(rowMeans(dist))]]

Ahat = melt(Ahat, varnames = c("OTU", "LatentFeature"), 
                 value.name = "Indicator")

ggplot(Ahat, aes(x = LatentFeature, y = OTU)) +
  geom_tile(aes(fill = Indicator), color = "white") +
  scale_fill_gradient(low = "white", high = "seagreen") +
  scale_x_discrete(expand = c(0, 0), breaks = 1 : Khat) +
  scale_y_discrete(expand = c(0, 0), breaks = 1 : length(OTU)) +
  theme(legend.position = "none")
