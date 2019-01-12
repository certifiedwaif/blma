context("comCrime")

library(tidyverse)

get_comCrime  <- function()
{
    # Fill in
    data(comData)
    Y <- comData[, 1:18]
    X <- comData[, 19:142]
    # Data preparation

    sum.na <- function(x) {  sum(is.na(x)); }
    inds <- which(apply(X,2,sum.na)==0)
    X2 <- X[,inds]
    X3 <- X2[,!colnames(X2)%in%c("ownHousQrange","rentQrange")]

    y <- Y %>% pull(18)
    inds <- which(is.na(y))

    vy <- y[-inds]
    mX <- X3[-inds,] %>% as.matrix
    mX.til <- cbind(1,mX)

    n <- length(vy)
    p <- ncol(mX)

    mult <- sqrt(n/(n-1))
    mX <- mX
    for (j in 1:p) {
      mX[,j] = mult*(mX[,j] - mean(mX[,j]))/sd(mX[,j])
    }
    vy <- mult*(vy - mean(vy))/sd(vy)
    return(list(vy=vy, mX=mX))
}
