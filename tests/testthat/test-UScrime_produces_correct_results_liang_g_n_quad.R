test_that("UScrime produces correct results liang_g_n_quad", {
	UScrime <- get_UScrime()
	vy <- UScrime$vy
	mX <- UScrime$mX
	tic("UScrime produces correct results liang_g_n_quad")
	result <- blma(vy, mX, prior="liang_g_n_quad", modelprior="uniform", cores=cores)
	toc()
	expect_equal(result$vinclusion_prob, c(
	  0.6510319616402197,0.2291246515020669,0.8651016578499309,0.6951167646184457,
	  0.4252168286795515,0.2025876028333554,0.3258937799141847,0.5663483825834370,
	  0.3560816955345109,0.2429057369231485,0.4974515067645519,0.2163153896529443,
	  0.9966071426484949,0.8454733725259110,0.2265274029341707
	), tolerance = 1e-8)
})