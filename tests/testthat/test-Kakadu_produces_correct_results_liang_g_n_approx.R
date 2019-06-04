
library(testthat)
library(blma)
library(tictoc)
library(parallel)

cores <- detectCores()

test_that("Kakadu produces correct results liang_g_n_approx", {
	Kakadu <- get_Kakadu()
	vy <- Kakadu$vy
	mX <- Kakadu$mX
	tic("Kakadu produces correct results liang_g_n_approx")
	result <- blma(vy, mX, prior="liang_g_n_approx", modelprior="uniform", cores=cores)
	toc()
	expect_equal(result$vinclusion_prob, c(
		0.3296136691294210,0.4997609497180022,0.1579846499019243,0.4631387408540441,
		0.9006995511063154,0.4005475083511888,0.1982925439416997,0.2209076900143787,
		0.1575545973403185,0.8998253808608440,0.9453194636814199,0.9996900361221982,
		0.1456947072049949,0.3654752994812290,0.3329298710932005,0.8268178993825527,
		0.1806657458627090,0.8573731512895385,0.9999999999992192,0.6083132163536824,
		0.9999983508359184,0.2714475290242928
	), tolerance = 1e-8)
})