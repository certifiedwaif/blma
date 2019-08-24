library(tidyverse)
library(blma)
library(parallel)

cores <- detectCores()
set.seed(2019)

data_sets <- c('eyeData', 'comCrime')
priors <- c('BIC', 'hyper-g', 'robust', 'beta-prime')

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
  p <- ncol(mX)

  modelprior <- 'beta-binomial' 
  modelpriorvec <- c(1, p)
  sampler(100000, vy, mX, prior=prior, modelprior=modelprior, modelpriorvec=modelpriorvec, cores=cores)
}

sampler_calls <- cross2(data_sets, priors) %>%
                   transpose %>%
                   pmap(call_sampler)
