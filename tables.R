library(blma)
library(parallel)
library(tictoc)

save_table_data <- function()
{
	data_sets <- c("Kakadu", "UScrime", "comCrime", "eyeData")
	priors <- c("BIC", "ZE", "liang_g1", "liang_g2", "liang_g_n_approx", "liang_g_n_quad", "robust_bayarri1", "robust_bayarri2", "zellner_siow_gauss_laguerre")
	results <- list()
	for (data_set in data_sets) {
		results[data_set] <- list()
		for (prior in priors) {
			cat(data_set, prior)
			ds <- NULL
			if (data_set == "Kakadu") ds <- get_Kakadu()
			if (data_set == "UScrime") ds <- get_UScrime()
			if (data_set == "comCrime") ds <- get_comCrime()
			if (data_set == "eyeData") ds <- get_eyeData()
			if (is.null(ds)) break
			vy <- ds$vy
			mX <- ds$mX
			p <- ncol(mX)
 	 		cores <- detect_cores()
 	 		# Our package is currently not thread-safe for these cases
			if (data_set %in% c("Kakadu", "UScrime") && prior %in% c("liang_g_n_appell", "zellner_siow_gauss_laguerre")) { 
				cores <- 1
			}
			tic()
			if (data_set %in% c("Kakadu", "UScrime")) {
				cat("exact")
				modelprior <- "uniform"
				results[data_set][prior] = blma(vy, mX, prior = prior, modelprior = modelprior, cores = cores)
			} else {
				cat("sampler")
				modelprior <- "beta-binomial"
				modelpriorvec <- c(1, p)
				results[data_set][prior] = sampler(
												vy,
												mX,
												prior = prior,
												modelprior = modelprior,
												modelpriorvec=modelpriorvec,
												cores = cores
											)
			}
			results[data_set][prior]["tictoc"] <- toc()
			save(results, file="results.RData")
			cat("\n")
		}
	}
}

save_table_data()
