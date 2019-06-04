
library(testthat)
library(blma)
library(tictoc)
library(parallel)

cores <- detectCores()

test_that("UScrime produces correct results robust_bayarri1", {
	UScrime <- get_UScrime()
	vy <- UScrime$vy
	mX <- UScrime$mX
	tic("UScrime produces correct results robust_bayarri1")
	result <- blma(vy, mX, prior="robust_bayarri1", modelprior="uniform", cores=cores)
	toc()
	expect_equal(result$vinclusion_prob, c(
	  0.6473883232855026,0.2450615600524966,0.8559561722433557,0.6902233614648986,
	  0.4408054425401610,0.2203974917215597,0.3408084056726517,0.5647124777936052,
	  0.3635295995251649,0.2577731178654963,0.4966009918811810,0.2339193088407784,
	  0.9954515035508796,0.8345346895344056,0.2451791177061960
	), tolerance = 1e-8)
})