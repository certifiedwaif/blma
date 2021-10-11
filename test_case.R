R -d gdb
break all_correlations_main
y
run
library(blma); library(MASS)
dat <- UScrime
dat[,-c(2,ncol(UScrime))] <- log(dat[,-c(2,ncol(UScrime))])
vy <- dat$y
mX <- data.matrix(cbind(dat[1:15]))
colnames(mX) <- c("log(AGE)","S","log(ED)","log(Ex0)","log(Ex1)", "log(LF)","log(M)","log(N)","log(NW)","log(U1)","log(U2)","log(W)",
"log(X)","log(prison)","log(time)")
blma_result <- blma(vy, mX, prior="ZE")
blma_result
normalize <- function(y, X)
{
  n <- length(y)
  p <- ncol(X)

  mu.y <- mean(y)
  sigma2.y <- (n - 1) * var(y) / n
  vy <- (y - mu.y) / sqrt(sigma2.y)

  # Normalise covariates
  mX <- matrix(0, n, p)
  mu.x <- c()
  sigma2.x <- c()
  for (j in 1:p)
  {
    mu.x[j] <- mean(X[, j])
    sigma2.x[j] <- (n - 1) * var(X[, j]) / n
    mX[, j] <- (X[, j] - mu.x[j]) / sqrt(sigma2.x[j])
  }

  return(list(vy = vy, mX = mX, mu.y = mu.y, sigma2.y = sigma2.y, mu.x = mu.x, sigma2.x = sigma2.x))
}
norm_UScrime <- normalize(vy, mX)
normed_blma_result <- blma(norm_UScrime$vy, norm_UScrime$mX, prior="ZE")
normed_blma_result
