test_that("UScrime produces correct results zellner_siow_gauss_laguerre", {
	UScrime <- get_UScrime()
	vy <- UScrime$vy
	mX <- UScrime$mX
	tic("UScrime produces correct results zellner_siow_gauss_laguerre")
	result <- blma(vy, mX, prior="zellner_siow_gauss_laguerre", modelprior="uniform", cores=1)
	toc()
	expect_equal(result$vinclusion_prob, 
c(
0.6529769401658134,0.2267683042922373,0.8679006779519631,0.6958693064404039,
0.4227195813764409,0.1997499011462964,0.3232070948480962,0.5673433815457091,
0.3556018636146127,0.2408930357991375,0.4988904773520176,0.2135955436467171,
0.9969385837823601,0.8482357484625307,0.2236271676907234
)

	, tolerance = 1e-8)
})