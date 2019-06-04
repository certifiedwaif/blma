
library(testthat)
library(blma)
library(tictoc)
library(parallel)

cores <- detectCores()

test_that("UScrime produces correct results liang_g_n_approx", {
	UScrime <- get_UScrime()
	vy <- UScrime$vy
	mX <- UScrime$mX
	tic("UScrime produces correct results liang_g_n_approx")
	result <- blma(vy, mX, prior="liang_g_n_approx", modelprior="uniform", cores=cores)
	toc()
	expect_equal(result$vinclusion_prob, c(
	  0.6572101998512460,0.2246852625428866,0.8724241652237846,0.6988559837509296,
	  0.4187705059074147,0.1972665987685162,0.3200139348087251,0.5707013245228147,
	  0.3570504620645473,0.2399912284225987,0.5038302968886970,0.2112239440381476,
	  0.9971640443211224,0.8531712185011106,0.2205035229588246
	), tolerance = 1e-8)
})