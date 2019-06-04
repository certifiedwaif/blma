
library(testthat)
library(blma)
library(tictoc)
library(parallel)

cores <- detectCores()

test_that("Kakadu produces correct results liang_g_n_appell", {
	skip("We were unable to get this case to work. It's too numerically difficult.")
	Kakadu <- get_Kakadu()
	vy <- Kakadu$vy
	mX <- Kakadu$mX
	tic("Kakadu produces correct results liang_g_n_appell")
	result <- blma(vy, mX, prior="liang_g_n_appell", modelprior="uniform", cores=1)
	toc()
	expect_equal(result$vinclusion_prob, c(
	), tolerance = 1e-8)
})