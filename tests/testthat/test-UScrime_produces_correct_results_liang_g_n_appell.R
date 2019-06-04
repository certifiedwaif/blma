
library(testthat)
library(blma)
library(tictoc)
library(parallel)

cores <- detectCores()

test_that("UScrime produces correct results liang_g_n_appell", {
	UScrime <- get_UScrime()
	vy <- UScrime$vy
	mX <- UScrime$mX
	tic("UScrime produces correct results liang_g_n_appell")
	result <- blma(vy, mX, prior="liang_g_n_appell", modelprior="uniform", cores=1)
	toc()
	expect_equal(result$vinclusion_prob, c(
	  0.6510317307247586,0.2291245644515895,0.8651014024703148,0.6951167568300936,
	  0.4252167906333386,0.2025875318978383,0.3258937597030980,0.5663482122162754,
	  0.3560815607057162,0.2429056433304466,0.4974513113404140,0.2163153102625932,
	  0.9966071300164177,0.8454731057912949,0.2265273203490445
	), tolerance = 1e-8)
})