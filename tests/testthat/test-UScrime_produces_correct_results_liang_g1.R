
library(testthat)
library(blma)
library(tictoc)
library(parallel)

cores <- detectCores()

test_that("UScrime produces correct results liang_g1", {
	UScrime <- get_UScrime()
	vy <- UScrime$vy
	mX <- UScrime$mX
	tic("UScrime produces correct results liang_g1")
	result <- blma(vy, mX, prior="liang_g1", modelprior="uniform", cores=cores)
	toc()
	expect_equal(result$vinclusion_prob, c(
		0.6592570082818343,0.2551890195947963,0.8623073076292519,0.6920362397961507,
		0.4460647456912883,0.2306155050885771,0.3455094697420374,0.5734193328978876,
		0.3765879959676905,0.2706204512952343,0.5124771409697690,0.2445691137122302,
		0.9950499664783635,0.8387376358181808,0.2549263051824415
	), tolerance = 1e-8)
})