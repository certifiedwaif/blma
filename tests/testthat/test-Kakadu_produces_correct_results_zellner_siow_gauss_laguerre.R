
library(testthat)
library(blma)
library(tictoc)
library(parallel)

cores <- detectCores()

test_that("Kakadu produces correct results zellner_siow_gauss_laguerre", {
	Kakadu <- get_Kakadu()
	vy <- Kakadu$vy
	mX <- Kakadu$mX
	tic("Kakadu produces correct results zellner_siow_gauss_laguerre")
	result <- blma(vy, mX, prior="zellner_siow_gauss_laguerre", modelprior="uniform", cores=1)
	toc()
	expect_equal(result$vinclusion_prob, 
c(
0.20578507196599863,0.47307551042974938,0.07612923670551421,
0.42112901060440383,0.86215067202350637,0.26921889275215943,
0.09053504955485733,0.11279139046662487,0.07323093390440316,
0.78214745397319363,0.93756983435089636,0.99939292585211370,
0.06720265691059692,0.20221362639019935,0.18783534784339340,
0.75520354356412289,0.08676136168658173,0.75285716863824470,
0.99999999999929501,0.44505775064088876,0.99999873218671564,
0.13452936396921064
)
	, tolerance = 1e-8)
})