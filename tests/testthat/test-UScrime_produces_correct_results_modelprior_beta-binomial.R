
library(testthat)
library(blma)
library(tictoc)
library(parallel)

cores <- detectCores()

test_that("UScrime produces correct results modelprior beta-binomial", {
	UScrime <- get_UScrime()
	vy <- UScrime$vy
	mX <- UScrime$mX
	p <- ncol(mX)
	tic("UScrime produces correct results modelprior beta-binomial")
	result <- blma(vy, mX, prior="BIC", modelprior="beta-binomial", modelpriorvec = c(1, p), cores=cores)
	toc()
	expect_equal(result$vinclusion_prob, c(
		0.28382041872785302,0.05401025213196230,0.52472573681024470,
		0.67899545738306066,0.34415858271523636,0.05257683401220613,
		0.31709892405155543,0.32163128053123047,0.09511945604467670,
		0.04703321836137517,0.13984669234378327,0.05271970967353467,
		0.99021404499476029,0.54819562466110150,0.06534213149560644
	), tolerance = 1e-8)
})