test_that("UScrime produces correct results BIC", {
	UScrime <- get_UScrime()
	vy <- UScrime$vy
	mX <- UScrime$mX
	tic("UScrime produces correct results BIC")
	result <- blma(vy, mX, prior="BIC", modelprior="uniform", cores=cores)
	toc()
	expect_equal(result$vinclusion_prob, c(
		0.7086661153970550,0.1906278707126815,0.9207138264989706,0.7253039423829655,
		0.3701199526247380,0.1582273002726970,0.2705746685548234,0.6064309892368921,
		0.3692402458531630,0.2192145020585733,0.5583788608712796,0.1738658933193491,
		0.9992269387743786,0.9027057249464372,0.1763387397295991
	), tolerance = 1e-8)
})