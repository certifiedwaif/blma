library(tidyverse)
library(blma)
library(parallel)

cores <- detectCores()
set.seed(2019)

data_sets <- c('eyeData', 'comCrime')
priors <- c('BIC', 'liang_g1', 'robust_bayarri1', 'zellner_siow_gauss_laguerre')

call_sampler <- function(data_set, prior)
{
  # Load data set
  if (data_set == 'eyeData') {
    eyeData <- get_eyeData()
    vy <- eyeData$vy
    mX <- eyeData$mX
  }
  if (data_set == 'comCrime') {
    comCrime <- get_comCrime()
    vy <- comCrime$vy
    mX <- comCrime$mX
  }

  n <- nrow(mX)
  p <- ncol(mX)

  if (n < p) {
    modelprior <- 'beta-binomial' 
    modelpriorvec <- c(1, p)
  } else {
    modelprior <-'uniform'
    modelpriorvec <- NULL
  }

  sampler(100000, vy, mX, prior=prior, modelprior=modelprior, modelpriorvec=modelpriorvec, cores=cores)
}

sampler_calls <- cross2(data_sets, priors) %>%
                   transpose %>%
                   pmap(call_sampler)
save(sampler_calls, file='sampler_calls.rda')
sampler_calls_df <- cross_df(list(data_sets=data_sets, priors=priors)) %>%
                      mutate(sampler_calls=sampler_calls)
save(sampler_calls_df, file='sampler_calls_df.rda')

get_top_inclusions <- function(result) {
  ord <- order(result$vinclusion_prob, decreasing=TRUE)[1:10]
  vinclusion_probs <- result$vinclusion_prob[ord]
  colnames(vinclusion_probs) <- colnames(eyeData$mX)[ord]
  vinclusion_probs
}
