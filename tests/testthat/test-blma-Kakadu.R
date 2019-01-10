context("blma Kakadu")

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

get_Kakadu <- function()
{
	data(Kakadu)
	y.t <- as.vector(Kakadu$income)
	X.f <- Kakadu[,c(2:21,23)]
	X.f <- model.matrix(~.,data=X.f)[,-1]
	res <- normalize(y.t, X.f)
	vy <- res$vy
	mX <- res$mX
	return(list(vy=vy, mX=mX))
}

# modelprior = uniform
test_that("Kakadu produces correct results BIC", {
	Kakadu <- get_Kakadu()
	vy <- Kakadu$vy
	mX <- Kakadu$mX
	result <- blma(vy, mX, prior="BIC", modelprior="uniform", cores=144)
	expect_equal(result$vinclusion_prob, c(
		0.11956375866122319,0.43598166701152036,0.03004602756871903,
		0.37141951236536819,0.81868011537400154,0.16833752435869168,
		0.03221165979213939,0.04298271555297178,0.02617847083497719,
		0.52531028106054500,0.92513375271585307,0.99815130365593652,
		0.02454673388614505,0.08100948178569978,0.08169622092074036,
		0.62985862997075337,0.03271752259200703,0.54749994288839865,
		0.99999999999905875,0.26634597085271527,0.99999886673686156,
		0.04953096133476929
	), tolerance = 1e-8)
})

test_that("Kakadu produces correct results ZE", {
	Kakadu <- get_Kakadu()
	vy <- Kakadu$vy
	mX <- Kakadu$mX
	result <- blma(vy, mX, prior="ZE", modelprior="uniform", cores=144)
	expect_equal(result$vinclusion_prob, c(
		0.20355684150991538,0.47244398575065949,0.07485696352142855,
		0.42017464791967363,0.86114935243204360,0.26671059509879097,
		0.08889221213490464,0.11092303936894550,0.07192086903966989,
		0.77777906132245778,0.93733237674906777,0.99937746373248582,
		0.06601135712668889,0.19907890141491377,0.18511189624094807,
		0.75297589180094671,0.08527867659504840,0.74927591041652553,
		0.99999999999929812,0.44106509651296133,0.99999873524470562,
		0.13221495773128050
	), tolerance = 1e-8)
})

test_that("Kakadu produces correct results liang_g1", {
	Kakadu <- get_Kakadu()
	vy <- Kakadu$vy
	mX <- Kakadu$mX
	result <- blma(vy, mX, prior="liang_g1", modelprior="uniform", cores=144)
	expect_equal(result$vinclusion_prob, c(
		0.3469442223401997,0.5033763713275147,0.1710391998737254,0.4686716902393477,
		0.9041282693203453,0.4182787535192226,0.2152694787243199,0.2365794398160191,
		0.1709195101828257,0.9081285313534126,0.9457607536048467,0.9996987392264256,
		0.1584196650685204,0.3866126303418119,0.3524291472817372,0.8329753936885851,
		0.1953702526810372,0.8655055789962616,0.9999999999992009,0.6255640528176040,
		0.9999982103316476,0.2911572010612115
	), tolerance = 1e-8)
})


test_that("Kakadu produces correct results liang_g2", {
	Kakadu <- get_Kakadu()
	vy <- Kakadu$vy
	mX <- Kakadu$mX
	result <- blma(vy, mX, prior="liang_g2", modelprior="uniform", cores=144)
	expect_equal(result$vinclusion_prob, c(
		0.3469442223402000,0.5033763713275141,0.1710391998737261,0.4686716902393479,
		0.9041282693203457,0.4182787535192236,0.2152694787243207,0.2365794398160196,
		0.1709195101828260,0.9081285313534141,0.9457607536048475,0.9996987392264269,
		0.1584196650685208,0.3866126303418126,0.3524291472817367,0.8329753936885853,
		0.1953702526810377,0.8655055789962616,0.9999999999992023,0.6255640528176055,
		0.9999982103316490,0.2911572010612116
	), tolerance = 1e-8)
})

#test_that("Kakadu produces correct results liang_g_n_appell", {
	#Kakadu <- get_Kakadu()
	#vy <- Kakadu$vy
	#mX <- Kakadu$mX
	#result <- blma(vy, mX, prior="liang_g_n_appell", modelprior="uniform")
	#expect_equal(result$vinclusion_prob, c(
	#), tolerance = 1e-8)
#})

#test_that("Kakadu produces correct results liang_g_n_approx", {
	#Kakadu <- get_Kakadu()
	#vy <- Kakadu$vy
	#mX <- Kakadu$mX
	#result <- blma(vy, mX, prior="liang_g_n_approx", modelprior="uniform")
	#expect_equal(result$vinclusion_prob, c(
	#), tolerance = 1e-8)
#})

#test_that("Kakadu produces correct results liang_g_n_quad", {
	#Kakadu <- get_Kakadu()
	#vy <- Kakadu$vy
	#mX <- Kakadu$mX
	#result <- blma(vy, mX, prior="liang_g_n_quad", modelprior="uniform")
	#expect_equal(result$vinclusion_prob, c(
	#), tolerance = 1e-8)
#})

#test_that("Kakadu produces correct results robust_bayarri1", {
	#Kakadu <- get_Kakadu()
	#vy <- Kakadu$vy
	#mX <- Kakadu$mX
	#result <- blma(vy, mX, prior="robust_bayarri1", modelprior="uniform")
	#expect_equal(result$vinclusion_prob, c(
	#), tolerance = 1e-8)
#})

#test_that("Kakadu produces correct results robust_bayarri2", {
	#Kakadu <- get_Kakadu()
	#vy <- Kakadu$vy
	#mX <- Kakadu$mX
	#result <- blma(vy, mX, prior="robust_bayarri2", modelprior="uniform")
	#expect_equal(result$vinclusion_prob, c(
	#), tolerance = 1e-8)
#})

#test_that("Kakadu produces correct results zellner_siow_gauss_laguerre", {
	#Kakadu <- get_Kakadu()
	#vy <- Kakadu$vy
	#mX <- Kakadu$mX
	#result <- blma(vy, mX, prior="zellner_siow_gauss_laguerre", modelprior="uniform")
	#expect_equal(result$vinclusion_prob, c(
		
	#), tolerance = 1e-8)
#})

## prior = BIC, modelprior = beta-binomial
#test_that("Kakadu produces correct results modelprior beta-binomial", {
	#Kakadu <- get_Kakadu()
	#vy <- Kakadu$vy
	#mX <- Kakadu$mX
	#p <- ncol(mX)
	#result <- blma(vy, mX, prior="BIC", modelprior="beta-binomial", modelpriorvec = c(1, p))
	#expect_equal(result$vinclusion_prob, c(
	#), tolerance = 1e-8)
#})
