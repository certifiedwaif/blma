library(blma)
library(parallel)
library(tictoc)
library(tidyverse)
library(glue)

RESULTS_FILE_NAME <- "results.RData"
data_sets <- c("Kakadu", "UScrime", "comCrime", "eyeData")
priors <- c("BIC", "ZE", "liang_g1", "liang_g2", "liang_g_n_approx", "liang_g_n_quad", "robust_bayarri1", "robust_bayarri2", "zellner_siow_gauss_laguerre")

save_table_data <- function()
{
	results <- list()
	for (data_set in data_sets) {
		results[[data_set]] <- list()
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
 	 		cores <- detectCores()
 	 		# Our package is currently not thread-safe for these cases
			if (data_set %in% c("Kakadu", "UScrime") && prior %in% c("liang_g_n_appell", "zellner_siow_gauss_laguerre")) { 
				cores <- 1
			}
			tic()
			if (data_set %in% c("Kakadu", "UScrime")) {
				cat(" exact ")
				modelprior <- "uniform"
				vinclusion_prob <- blma(vy, mX, prior = prior, modelprior = modelprior, cores = cores)$vinclusion_prob
				results[[data_set]][[prior]] <- list(vinclusion_prob=vinclusion_prob)
			} else {
				cat(" sampler ")
				modelprior <- "beta-binomial"
				modelpriorvec <- c(1, p)
				vinclusion_prob <- sampler(
												100000,
												vy,
												mX,
												prior = prior,
												modelprior = modelprior,
												modelpriorvec=modelpriorvec,
												cores = cores
											)$vinclusion_prob
				results[[data_set]][[prior]] <- list(vinclusion_prob=vinclusion_prob)
			}
			results[[data_set]][[prior]][["tictoc"]] <- toc()
			save(results, file=RESULTS_FILE_NAME)
		}
	}
	return(results)
}

produce_tables <- function(results)
{
	for (data_set in names(results)) {
		produce_table(results[[data_set]])
	}
}

produce_table <- function(result)
{
	table_cols <- length(result)
	table_rows <- length(result[[1]]$vinclusion_prob)
	# Produce start of table
	cat("\\begin{sidewaystable}[h!]\n\\begin{center}\n{\\tiny\n\\tabular{\n")
	cat("c|")
	# TODO(Mark): Write code to put | at the column splits
	for (col in 1:table_cols) {
		cat("r")
	}
	cat("}\n")
	# Produce header rows
	cat("Package")
	for (col in 1:table_cols) {
		cat("&", "BLMA")
	}
	# Produce each table row for inclusion probabilities
	for (row in 1:table_rows) {
		cat(row, "&")
		for (col in 1:table_cols) {
			cat(round(result[[col]]$vinclusion_prob[row] * 100.0, 2))
			if (col < table_cols) cat("&")
		}
		cat("\\\\")
	}
	# Produce end of table
	cat("\\hline\n")
	cat("\\end{tabular}\n}\n\\end{center}")
}

if (file.exists(RESULTS_FILE_NAME)) {
	load(file=RESULTS_FILE_NAME)
} else {
	results <- save_table_data()
}
produce_tables(results)
