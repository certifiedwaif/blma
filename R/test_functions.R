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

#' get_Kakadu
#'
#' @return A list containing vy and mX for the Kakadu data set
#' @export
get_Kakadu <- function()
{
	data(Kakadu)
	y.t <- as.vector(Kakadu$income)
	X.f <- Kakadu[,c(2:21,23)]
	X.f <- model.matrix(~.,data=X.f)[,-1]
	res <- normalize(y.t, X.f)
	vy <- res$vy
	mX <- res$mX
	return(list(vy=vy, mX=mX))
}

#' get_UScrime
#'
#' @return A list containing vy and mX for the UScrime data set
#' @export
get_UScrime <- function()
{
	mD <- MASS::UScrime

	notlog <- c(2,ncol(MASS::UScrime))

	mD[,-notlog] <- log(mD[,-notlog])

	for (j in 1:ncol(mD)) {
		mD[,j] <- (mD[,j] - mean(mD[,j]))/sd(mD[,j])
	}
	vy <- mD$y
	mX <- data.matrix(cbind(mD[1:15]))
	norm <- normalize(vy, mX)
	vy <- norm$vy
	mX <- norm$mX
	return(list(vy=vy, mX=mX))
}

#' get_eyeData
#'
#' @return A list containing vy and mX for the eyeData data set
#' @export
get_eyeData <- function()
{
    data(eyeData)
    vy <- y
    mX <- x
    norm <- normalize(vy, mX)
    vy <- norm$vy
    mX <- norm$mX
    return(list(vy=vy, mX=mX))
}

#' get_comCrime
#'
#' @return A list containing vy and mX for the comCrime data set
#' @export
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

    y <- Y[, 18]
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
