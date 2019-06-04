
library(testthat)
library(blma)
library(tictoc)
library(parallel)

cores <- detectCores()

test_that("UScrime produces correct results liang_g2", {
	UScrime <- get_UScrime()
	vy <- UScrime$vy
	mX <- UScrime$mX
	tic("UScrime produces correct results liang_g2")
	result <- blma(vy, mX, prior="liang_g2", modelprior="uniform", cores=cores)
	toc()
	expect_equal(result$vinclusion_prob, c(
		0.6592570082818349,0.2551890195947963,0.8623073076292518,0.6920362397961510,
		0.4460647456912886,0.2306155050885771,0.3455094697420375,0.5734193328978877,
		0.3765879959676907,0.2706204512952343,0.5124771409697698,0.2445691137122302,
		0.9950499664783639,0.8387376358181811,0.2549263051824416
	), tolerance = 1e-8)
})