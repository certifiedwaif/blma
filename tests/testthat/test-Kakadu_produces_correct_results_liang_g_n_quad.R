
library(testthat)
library(blma)
library(tictoc)
library(parallel)

cores <- detectCores()

test_that("Kakadu produces correct results liang_g_n_quad", {
	Kakadu <- get_Kakadu()
	vy <- Kakadu$vy
	mX <- Kakadu$mX
	tic("Kakadu produces correct results liang_g_n_quad")
	result <- blma(vy, mX, prior="liang_g_n_quad", modelprior="uniform", cores=cores)
	toc()
	expect_equal(result$vinclusion_prob, c(
		0.3203591831238068,0.4977774306550474,0.1512859110442455,0.4601315123732897,
		0.8984723576023463,0.3910244742565688,0.1894701351919556,0.2126145072480971,
		0.1506655077406368,0.8944304522063428,0.9448935925720320,0.9996801335893382,
		0.1391831130834718,0.3538609740786210,0.3223705659773634,0.8229190449064608,
		0.1730543020195790,0.8521810234828430,0.9999999999993161,0.5982779020650832,
		0.9999983946351332,0.2609239705895227
	), tolerance = 1e-8)
})