#' @title NMGS posterior summaries
#' @description: Compute summary statistics of the sampled MCMC parameters,
#'        excluding the burn-in samples and thinning. 
#' @param samples Posterior samples of the NMGS model (output from the read_nmgs function)
#' @param burnin Burn-in sample 
#' @param thinning Thinning
#' @return median, mean, standard deviation and upper and lower 95\% confidence 
#'         intervals for the biodiversity parameter theta and the immigration rates. 
#'         Parameters x summaries matrix.
#' @details The confidence intervals are calculated as the 2.5\% and 97.5\% quantiles.
#' @export
#' @references See citation("NMGS") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities
nmgs_posterior_summaries <- function (samples, burnin = 0, thinning = 1) {

  # Exclude the burn-in samples	and return median for all parameters
  # plus lower and upper quantiles
  s <- seq(burnin + 1, nrow(samples), thinning) 

  # Samples medians
  med <- apply(samples[s, -1], 2, "median")

  # Samples means
  ave <- apply(samples[s, -1], 2, "mean")

  # Samples stds
  sds <- apply(samples[s, -1], 2, "sd")
  
  # Lower and upper quantiles
  lc <- apply(samples[s, -1], 2, function (x) {quantile(x, 0.025)})
  uc <- apply(samples[s, -1], 2, function (x) {quantile(x, 0.975)})

  cbind(lower = lc, median = med, mean = ave, sd = sds, upper = uc)

}



