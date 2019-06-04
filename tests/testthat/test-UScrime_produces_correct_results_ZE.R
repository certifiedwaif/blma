
library(testthat)
library(blma)
library(tictoc)
library(parallel)

cores <- detectCores()

test_that("UScrime produces correct results ZE", {
	UScrime <- get_UScrime()
	vy <- UScrime$vy
	mX <- UScrime$mX
	tic("UScrime produces correct results ZE")
	result <- blma(vy, mX, prior="ZE", modelprior="uniform", cores=cores)
	toc()
	expect_equal(result$vinclusion_prob, c(
		0.6551129989168902,0.2287509395743590,0.8691121295594134,0.6964900562039598,
		0.4236469337382322,0.2017927025100451,0.3243432277162467,0.5691166496984528,
		0.3580859703810375,0.2434951699745767,0.5019137481840252,0.2157295344341948,
		0.9969490852909634,0.8492050249433519,0.2255266067521766
	), tolerance = 1e-8)
})