context("exact UScrime")

library(parallel)

cores <- detectCores()

normalize <- function(y, X)
{
  n <- length(y)
  p <- ncol(X)

  mu.y <- mean(y)
  sigma2.y <- (n - 1) * var(y) / n
  vy <- (y - mu.y) / sqrt(sigma2.y)

  # Normalise covariates
  mX <- matrix(0, n, p)
  mu.x <- c()
  sigma2.x <- c()
  for (j in 1:p)
  {
    mu.x[j] <- mean(X[, j])
    sigma2.x[j] <- (n - 1) * var(X[, j]) / n
    mX[, j] <- (X[, j] - mu.x[j]) / sqrt(sigma2.x[j])
  }

  return(list(vy = vy, mX = mX, mu.y = mu.y, sigma2.y = sigma2.y, mu.x = mu.x, sigma2.x = sigma2.x))
}

get_UScrime <- function()
{
	mD <- MASS::UScrime

	notlog <- c(2,ncol(MASS::UScrime))

	mD[,-notlog] <- log(mD[,-notlog])

	for (j in 1:ncol(mD)) {
		mD[,j] <- (mD[,j] - mean(mD[,j]))/sd(mD[,j])
	}
	vy <- mD$y
	mX <- data.matrix(cbind(mD[1:15]))
	norm <- normalize(vy, mX)
	vy <- norm$vy
	mX <- norm$mX
	return(list(vy=vy, mX=mX))
}

# modelprior = uniform
test_that("UScrime produces correct results BIC", {
	skip("I don't care about this right now")
	UScrime <- get_UScrime()
	vy <- UScrime$vy
	mX <- UScrime$mX
	result <- blma(vy, mX, prior="BIC", modelprior="uniform", cores=cores)
	expect_equal(result$vinclusion_prob, c(
		0.7086661153970550,0.1906278707126815,0.9207138264989706,0.7253039423829655,
		0.3701199526247380,0.1582273002726970,0.2705746685548234,0.6064309892368921,
		0.3692402458531630,0.2192145020585733,0.5583788608712796,0.1738658933193491,
		0.9992269387743786,0.9027057249464372,0.1763387397295991
	), tolerance = 1e-8)
})

test_that("UScrime produces correct results ZE", {
	skip("I don't care about this right now")
	UScrime <- get_UScrime()
	vy <- UScrime$vy
	mX <- UScrime$mX
	result <- blma(vy, mX, prior="ZE", modelprior="uniform", cores=cores)
	expect_equal(result$vinclusion_prob, c(
		0.6551129989168902,0.2287509395743590,0.8691121295594134,0.6964900562039598,
		0.4236469337382322,0.2017927025100451,0.3243432277162467,0.5691166496984528,
		0.3580859703810375,0.2434951699745767,0.5019137481840252,0.2157295344341948,
		0.9969490852909634,0.8492050249433519,0.2255266067521766
	), tolerance = 1e-8)
})

test_that("UScrime produces correct results liang_g1", {
	skip("I don't care about this right now")
	UScrime <- get_UScrime()
	vy <- UScrime$vy
	mX <- UScrime$mX
	result <- blma(vy, mX, prior="liang_g1", modelprior="uniform", cores=cores)
	expect_equal(result$vinclusion_prob, c(
		0.6592570082818343,0.2551890195947963,0.8623073076292519,0.6920362397961507,
		0.4460647456912883,0.2306155050885771,0.3455094697420374,0.5734193328978876,
		0.3765879959676905,0.2706204512952343,0.5124771409697690,0.2445691137122302,
		0.9950499664783635,0.8387376358181808,0.2549263051824415
	), tolerance = 1e-8)
})

test_that("UScrime produces correct results liang_g2", {
	skip("I don't care about this right now")
	UScrime <- get_UScrime()
	vy <- UScrime$vy
	mX <- UScrime$mX
	result <- blma(vy, mX, prior="liang_g2", modelprior="uniform", cores=cores)
	expect_equal(result$vinclusion_prob, c(
		0.6592570082818349,0.2551890195947963,0.8623073076292518,0.6920362397961510,
		0.4460647456912886,0.2306155050885771,0.3455094697420375,0.5734193328978877,
		0.3765879959676907,0.2706204512952343,0.5124771409697698,0.2445691137122302,
		0.9950499664783639,0.8387376358181811,0.2549263051824416
	), tolerance = 1e-8)
})

test_that("UScrime produces correct results liang_g_n_appell", {
	skip("I don't care about this right now")
	UScrime <- get_UScrime()
	vy <- UScrime$vy
	mX <- UScrime$mX
	result <- blma(vy, mX, prior="liang_g_n_appell", modelprior="uniform", cores=cores)
	expect_equal(result$vinclusion_prob, c(
	  0.6510317307247586,0.2291245644515895,0.8651014024703148,0.6951167568300936,
	  0.4252167906333386,0.2025875318978383,0.3258937597030980,0.5663482122162754,
	  0.3560815607057162,0.2429056433304466,0.4974513113404140,0.2163153102625932,
	  0.9966071300164177,0.8454731057912949,0.2265273203490445
	), tolerance = 1e-8)
})

test_that("UScrime produces correct results liang_g_n_approx", {
	skip("I don't care about this right now")
	UScrime <- get_UScrime()
	vy <- UScrime$vy
	mX <- UScrime$mX
	result <- blma(vy, mX, prior="liang_g_n_approx", modelprior="uniform", cores=cores)
	expect_equal(result$vinclusion_prob, c(
	  0.6572101998512460,0.2246852625428866,0.8724241652237846,0.6988559837509296,
	  0.4187705059074147,0.1972665987685162,0.3200139348087251,0.5707013245228147,
	  0.3570504620645473,0.2399912284225987,0.5038302968886970,0.2112239440381476,
	  0.9971640443211224,0.8531712185011106,0.2205035229588246
	), tolerance = 1e-8)
})

test_that("UScrime produces correct results liang_g_n_quad", {
	skip("I don't care about this right now")
	UScrime <- get_UScrime()
	vy <- UScrime$vy
	mX <- UScrime$mX
	result <- blma(vy, mX, prior="liang_g_n_quad", modelprior="uniform", cores=cores)
	expect_equal(result$vinclusion_prob, c(
	  0.6510319616402197,0.2291246515020669,0.8651016578499309,0.6951167646184457,
	  0.4252168286795515,0.2025876028333554,0.3258937799141847,0.5663483825834370,
	  0.3560816955345109,0.2429057369231485,0.4974515067645519,0.2163153896529443,
	  0.9966071426484949,0.8454733725259110,0.2265274029341707
	), tolerance = 1e-8)
})

test_that("UScrime produces correct results robust_bayarri1", {
	skip("I don't care about this right now")
	UScrime <- get_UScrime()
	vy <- UScrime$vy
	mX <- UScrime$mX
	result <- blma(vy, mX, prior="robust_bayarri1", modelprior="uniform", cores=cores)
	expect_equal(result$vinclusion_prob, c(
	  0.6473883232855026,0.2450615600524966,0.8559561722433557,0.6902233614648986,
	  0.4408054425401610,0.2203974917215597,0.3408084056726517,0.5647124777936052,
	  0.3635295995251649,0.2577731178654963,0.4966009918811810,0.2339193088407784,
	  0.9954515035508796,0.8345346895344056,0.2451791177061960
	), tolerance = 1e-8)
})

test_that("UScrime produces correct results robust_bayarri2", {
	skip("I don't care about this right now")
	UScrime <- get_UScrime()
	vy <- UScrime$vy
	mX <- UScrime$mX
	result <- blma(vy, mX, prior="robust_bayarri2", modelprior="uniform", cores=cores)
	expect_equal(result$vinclusion_prob, c(
	  0.6473659645368532,0.2450935725931694,0.8559277873673219,0.6902103956460465,
	  0.4408397852569301,0.2204342991741803,0.3408415449185910,0.5647004484926239,
	  0.3635349998539846,0.2577975787131784,0.4965810378824649,0.2339548960414573,
	  0.9954498578038506,0.8345042702863100,0.2452186439348585
	), tolerance = 1e-8)
})

test_that("UScrime produces correct results zellner_siow_gauss_laguerre", {
	browser()
	UScrime <- get_UScrime()
	vy <- UScrime$vy
	mX <- UScrime$mX
	result <- blma(vy, mX, prior="zellner_siow_gauss_laguerre", modelprior="uniform", cores=cores)
	expect_equal(result$vinclusion_prob, c(
		0.6529769401655269,0.2267683042922874,0.8679006779518879,0.6958693064404211,
		0.4227195813764282,0.1997499011463117,0.3232070948482000,0.5673433815457235,
		0.3556018636146491,0.2408930357989732,0.4988904773514005,0.2135955436462610,
		0.9969385837823613,0.8482357484619247,0.2236271676913612
		
	), tolerance = 1e-8)
})

# prior = BIC, modelprior = beta-binomial
test_that("UScrime produces correct results modelprior beta-binomial", {
	UScrime <- get_UScrime()
	vy <- UScrime$vy
	mX <- UScrime$mX
	p <- ncol(mX)
	result <- blma(vy, mX, prior="BIC", modelprior="beta-binomial", modelpriorvec = c(1, p), cores=cores)
	expect_equal(result$vinclusion_prob, c(
		0.28382041872785302,0.05401025213196230,0.52472573681024470,
		0.67899545738306066,0.34415858271523636,0.05257683401220613,
		0.31709892405155543,0.32163128053123047,0.09511945604467670,
		0.04703321836137517,0.13984669234378327,0.05271970967353467,
		0.99021404499476029,0.54819562466110150,0.06534213149560644
	), tolerance = 1e-8)
})
