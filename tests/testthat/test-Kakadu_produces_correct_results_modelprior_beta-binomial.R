
library(testthat)
library(blma)
library(tictoc)
library(parallel)

cores <- detectCores()

test_that("Kakadu produces correct results modelprior beta-binomial", {
	Kakadu <- get_Kakadu()
	vy <- Kakadu$vy
	mX <- Kakadu$mX
	p <- ncol(mX)
	tic("Kakadu produces correct results modelprior beta-binomial")
	result <- blma(vy, mX, prior="BIC", modelprior="beta-binomial", modelpriorvec = c(1, p), cores=cores)
	toc()
	expect_equal(result$vinclusion_prob, c(
		0.056281843799996970,0.306314199561247891,0.007552489928519051,
		0.247616597500152053,0.742709650039889202,0.095631352679238363,
		0.006351015631741534,0.009101941725939111,0.005241750496263586,
		0.153817084200635196,0.907193762314498886,0.984666525605746967,
		0.005113387083452965,0.017754736585571605,0.021245104548298487,
		0.371729547784183711,0.006733387430678084,0.216129983894335492,
		0.999999999998894218,0.110633351931119187,0.999999030703775826,
		0.009702560359548688
	), tolerance = 1e-8)
})