
library(testthat)
library(blma)
library(tictoc)
library(parallel)

cores <- detectCores()

test_that("kakadu produces correct results robust_bayarri2", {
	Kakadu <- get_Kakadu()
	vy <- Kakadu$vy
	mX <- Kakadu$mX
	tic("kakadu produces correct results robust_bayarri2")
	toc()
	result <- blma(vy, mX, prior="robust_bayarri2", modelprior="uniform", cores=cores)
	expect_equal(result$vinclusion_prob, c(
		0.2645609801933496,0.4872708350999790,0.1119549057642774,0.4428310083492054,
		0.8859403686976074,0.3327165328329659,0.1381842637776008,0.1633673877358905,
		0.1104117927905364,0.8591578715885575,0.9435211879405069,0.9996313610316295,
		0.1012638529130122,0.2837032904835066,0.2586857176888843,0.7998311650014662,
		0.1284988708795221,0.8194918851458997,0.9999999999992317,0.5357989689844272,
		0.9999987783204403,0.1987407247563582
	), tolerance = 1e-5)
})