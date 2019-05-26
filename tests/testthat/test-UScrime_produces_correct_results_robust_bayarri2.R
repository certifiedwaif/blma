test_that("UScrime produces correct results robust_bayarri2", {
	UScrime <- get_UScrime()
	vy <- UScrime$vy
	mX <- UScrime$mX
	tic("UScrime produces correct results robust_bayarri2")
	result <- blma(vy, mX, prior="robust_bayarri2", modelprior="uniform", cores=cores)
	toc()
	expect_equal(result$vinclusion_prob, c(
	  0.6473659645368532,0.2450935725931694,0.8559277873673219,0.6902103956460465,
	  0.4408397852569301,0.2204342991741803,0.3408415449185910,0.5647004484926239,
	  0.3635349998539846,0.2577975787131784,0.4965810378824649,0.2339548960414573,
	  0.9954498578038506,0.8345042702863100,0.2452186439348585
	), tolerance = 1e-8)
})