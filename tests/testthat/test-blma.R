context("blma")

gen_UScrime <- function()
{
  library(MASS)

  mD <- UScrime
  notlog <- c(2,ncol(UScrime))
  mD[,-notlog] <- log(mD[,-notlog])

  for (j in 1:ncol(mD)) {
    mD[,j] <- (mD[,j] - mean(mD[,j]))/sd(mD[,j])
  }

  varnames <- c(
    "log(AGE)",
    "S",
    "log(ED)",
    "log(Ex0)",
    "log(Ex1)",
    "log(LF)",
    "log(M)",
    "log(N)",
    "log(NW)",
    "log(U1)",
    "log(U2)",
    "log(W)",
    "log(X)",
    "log(prison)",
    "log(time)")

  y.t <- mD$y
  X.f <- data.matrix(cbind(mD[1:15]))
  colnames(X.f) <- varnames 
  return(list(y.t=y.t, X.f=X.f))
}


test_that("blma uniform modelprior", {
  UScrime_data <- gen_UScrime()
  y.t <- UScrime_data$y.t
  X.f <- UScrime_data$X.f
  blma_result <- blma(y.t, X.f, "ZE")
  expect_equal(blma_result$vinclusion_prob, c(0.6551317, 0.2287509,0.8691121,0.6964901,0.4236469,0.2017927,
  																						0.3243432,0.5691166,0.358086,0.2434952,0.5019137,0.2157295,
  																						0.9969491,0.849205,0.2255266), tolerance=1e-5)
})


test_that("blma beta-binomial modelprior", {
  UScrime_data <- gen_UScrime()
  y.t <- UScrime_data$y.t
  X.f <- UScrime_data$X.f
  blma_result <- blma(y.t, X.f, "BIC", "beta-binomial", c(2, 1))
  expect_equal(blma_result$vinclusion_prob, c(0.730, 0.232,0.911,0.738,0.398,0.207,
  																						0.312,0.637,0.419,0.286,0.600,0.228,
  																						0.999,0.898,0.222), tolerance=1e-3)
}) 


test_that("blma bernoulli modelprior", {
  UScrime_data <- gen_UScrime()
  y.t <- UScrime_data$y.t
  X.f <- UScrime_data$X.f
  blma_result <- blma(y.t, X.f, "BIC", "bernoulli", rep(1/15, 15))
  expect_equal(blma_result$vinclusion_prob, c(0.0472, 0.01270852,0.06138092, 0.0483536,
  																						0.02467466, 0.01054849, 0.01803831,
  																						0.04042873, 0.02461602, 0.0146143,
  																						0.03722526, 0.01159106, 0.06661513,
  																						0.06018038, 0.01175592), tolerance=1e-3)
})
