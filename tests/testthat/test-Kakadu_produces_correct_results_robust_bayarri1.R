
library(testthat)
library(blma)
library(tictoc)
library(parallel)

cores <- detectCores()

test_that("Kakadu produces correct results robust_bayarri1", {
	Kakadu <- get_Kakadu()
	vy <- Kakadu$vy
	mX <- Kakadu$mX
	tic("Kakadu produces correct results robust_bayarri1")
	result <- blma(vy, mX, prior="robust_bayarri1", modelprior="uniform", cores=cores)
	toc()
	expect_equal(result$vinclusion_prob, c(
		0.2645609801935807,0.4872708351002377,0.1119549057643435,0.4428310083494666,
		0.8859403686981394,0.3327165328332443,0.1381842637777420,0.1633673877360185,
		0.1104117927906037,0.8591578715890693,0.9435211879410100,0.9996313610321645,
		0.1012638529130835,0.2837032904837962,0.2586857176889374,0.7998311650018346,
		0.1284988708797292,0.8194918851464008,1.0000000000003251,0.5357989689853954,
		0.9999987783206028,0.1987407247561971
	), tolerance = 1e-8)
})