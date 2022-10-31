

#****************************************************************
#
# Multivariate Test under High Dimensions
# Adjusting covariate matrix X must include intercept 
#
#****************************************************************

mvt_ud <- function(Y, X, G, Sigma2 = NULL, w = NULL){
  
  n <- nrow(X)
  d <- ncol(X)
  p <- ncol(G)
  K <- ncol(Y)
  
  # use the estimator under the null if sigma2 not provided
  if(is.null(Sigma2)){
    resM <- matrix(NA, n, K)
    for(k in 1:K){
      fit <- lm(Y[, k] ~ X - 1)
      resM[, k] <- fit$residuals
    }
    Sigma2 <- cov(resM) * (n - 1) / (n - d)
  }
  
  Sigma22 <- Sigma2 ^ 2
  
  Px <- X %*% solve(t(X) %*% X) %*% t(X)
  IPx <- diag(1, n) - Px
  halfA <- t(G) %*% IPx
  vecD <- 1 / diag(halfA %*% t(halfA))
  halfA2 <- sweep(halfA, 1, vecD, FUN = "*")
  A <- t(halfA) %*% halfA2
  A1 <- A - p / (n - d) * IPx
  
  Q_star <- diag(t(Y) %*% A %*% Y)
  mu_star <- p * diag(Sigma2)
  
  sa2 <- sum(A ^ 2)
  sa12 <- sum(A1 ^ 2)
  
  ts1 <- t(Q_star - mu_star) %*% solve(Sigma22) %*% (Q_star - mu_star)
  Tl1 <- c(ts1) / (2 * sa2)
  
  Th1 <- c(ts1) / (2 * sa12)
  
  pl1 <- 1 - pchisq(Tl1, df = K)
  ph1 <- 1 - pchisq(Th1, df = K)
  
  # NOT Sigma22 Here.
  ts2 <- sum(Q_star - mu_star) / sqrt(2 * sum(Sigma2 ^ 2))
  
  Tl2 <- ts2 / sqrt(sa2)
  Th2 <- ts2 / sqrt(sa12)

  pl2 <- 2 * pnorm(abs(Tl2), lower.tail = FALSE) 
  ph2 <- 2 * pnorm(abs(Th2), lower.tail = FALSE) 


  res <- list(Tl1 = Tl1, pl1 = pl1, 
              Th1 = Th1, ph1 = ph1,
              Tl2 = Tl2, pl2 = pl2,
              Th2 = Th2, ph2 = ph2)
  

  res
}

