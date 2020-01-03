# Monte Carlo Simulation
#
# The script investigates resampling methods, in particular, bootstrap and 
# jacknife procedure through Monte Carlo simlulation.
#
# Author: Ming-Chen Lu (mingchlu@umich.edu)
# Updated: Jan 3, 2020
#80: ---------------------------------------------------------------------------

# Libraries: -------------------------------------------------------------------
library(tidyverse)

# Part a: ----------------------------------------------------------------------
get_mc_jacknife_ci = function(xmat, ymat, alpha = .05) {
  # Inputs: x, y - two matrices representing Monte Carlo replicates
  # Outputs: a list of string representing 100*(1-alpha)% C.I. for theta_hat = 
  #          x_bar/y_bar using the jackknife estimate for the standard error
  
  # Test the input
  stopifnot( is.matrix(xmat) )
  stopifnot( is.matrix(ymat) )
  
  # Compute mean ratio of two matrices within each row
  theta_hat = rowMeans(xmat) / rowMeans(ymat)
  
  # Compute the standard error using the jackknife estimate
  xmat_i = (rowSums(xmat) - xmat) / (ncol(xmat) - 1)
  ymat_i = (rowSums(ymat) - ymat) / (ncol(ymat) - 1)
  theta_hat_i = cbind(xmat_i / rowMeans(ymat), rowMeans(xmat) / ymat_i)
  theta_bar = rowMeans(theta_hat_i)
  se_jack = sqrt((ncol(theta_hat_i) - 1) * rowMeans((theta_hat_i - theta_bar)^2))
  
  # Upper and lower bound of confidence interval
  lcb = theta_hat - qnorm(1 - alpha/2) * se_jack
  ucb = theta_hat + qnorm(1 - alpha/2) * se_jack
  
  return(
    tibble(
      measure = 'jk', 
      lwr = lcb, 
      pnt_est = theta_hat, 
      upr = ucb)
    )
}

# Part b: ----------------------------------------------------------------------
get_mc_bootstrap_ci = function(xmat, ymat, alpha = .05, n_boot = 1e3) {
  # Inputs: 
  #     x, y - two matrices representing Monte Carlo replicates
  #     alpha - a numeric representing level of confidence interval
  #     n_boot - number of bootstrap samples
  # Outputs: a numeric vector representing 100*(1-alpha)% C.I. 
  
  # Test the input
  stopifnot( is.matrix(xmat) )
  stopifnot( is.matrix(ymat) )
  
  # Original mean ratio of the observed data
  theta_hat = rowSums(xmat) / rowSums(ymat)

  # Each row represents a Monte Carlo replicate
  # Bootstrap samples from column indices n_boot times
  # In this case, since all replicates are independent, we use the same indices
  # for each replicate.
  b_xmat = xmat[,sample(col(xmat), n_boot*ncol(xmat), replace = TRUE)]
  dim(b_xmat) = c(nrow(xmat), ncol(xmat), n_boot)
  
  b_ymat = ymat[,sample(col(ymat), n_boot*ncol(ymat), replace = TRUE)]
  dim(b_ymat) = c(nrow(ymat), ncol(ymat), n_boot)
  
  # Transpose an array by permuting and resizing its dimensions and compute
  # each bootstraping sample mean
  n_theta_hat = 
    colMeans(aperm(b_xmat, c(2,1,3))) / colMeans(aperm(b_ymat, c(2,1,3)))
  
  # alpha-quantile of the boortrap distribution
  qr = apply(n_theta_hat, 1, 
                 function(b) quantile(b, c(alpha/2, .5, 1-alpha/2)))
  
  # Standard deviation of the bootstrap distribution
  mat_c = n_theta_hat - rowMeans(n_theta_hat)
  sd = sqrt(rowSums((mat_c)^2) / {dim(n_theta_hat)[2] - 1})

  # return point estimates and 100(1 - alpha)% C.I for each method
  lcb = qr[1,]
  qr_est = qr[2,]
  ucb = qr[3,]
  
  return(
    tibble(
      measure = rep(c('qt', 'bs', 'nbs'), each = nrow(xmat)),
      lwr = c(lcb, 2*theta_hat - ucb, theta_hat - qnorm(1 - alpha/2)*sd),
      pnt_est = c(qr_est, theta_hat, theta_hat),
      upr = c(ucb, 2*theta_hat - lcb, theta_hat + qnorm(1 - alpha/2)*sd)
    )
  )
}

# Part c: ----------------------------------------------------------------------
# Comparing below 3 quantities for each of four CIs: ---------------------------
# i. coverage probability: the percentage of samples that contain the true value
# ii. the average length of the CIs produced by each method
# iii. the average shape of CIs produced by each method: e.g. the ratio of 
#      lengths on either side of the point estimate

# Monte Carlo simulations
mcrep = 1e4
n_x = 30
n_y = 30

# Generate x and y from normal distribution
# Each row is a dataset
xmat = rpois(n_x * mcrep, lambda = 4)
dim(xmat) = c(mcrep, n_x) 
ymat = rpois(n_y * mcrep, lambda = 2)
dim(ymat) = c(mcrep, n_y)

# Combine four confidence intervals
tb = 
  rbind(
    get_mc_jacknife_ci(xmat, ymat, alpha = .05), 
    get_mc_bootstrap_ci(xmat, ymat, alpha = .05, n_boot = 1e3)
  )

# To compute coverage probability, the true value of E[X]/E[Y] = 2 since 
# X~Poi(lambda = 4) and Y~Poi(lambda = 2)
target = 2

result = 
  tb %>%
  group_by(measure) %>%
  summarise( 
    cvrg_prob = mean( {lwr < target} & {upr > target} ),
    avg_len_ci = mean(upr - lwr),
    avg_sp_ci = mean( (upr - pnt_est) / (pnt_est - lwr) )
  ) %>%
  arrange(cvrg_prob, avg_len_ci, avg_sp_ci)

# Produce a table: -------------------------------------------------------------
result$measure = 
  result$measure %>% 
  recode(bs = 'basic bootstrap', 
         qt = 'percentile', 
         nbs = 'normal approx w/ bootstrap std error', 
         jk = 'jackknife estimate for std error')

cn = c('measure', 'Coverage Probability', 'Average Length', 'Average Shape')
caps = paste(
  '*Confidence interval comparisons.* This table compares the quantities for ',
  'each of the four confidence interval types. ',
  'The quantities being computed include the coverage probability, the average',
  ' length and the average shape of the confidence intervals.'
)

# knitr::kable(result, caption = caps, col.names = cn)

