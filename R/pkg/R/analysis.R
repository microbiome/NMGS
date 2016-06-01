#' Description: Compute summary statistics of the sampled MCMC parameters,
#' excluding the burn-in samples and thinning. 
#'
#' Arguments:
#'   @param samples Posterior samples of the NMGS model (output from the read_nmgs function)
#'   @param burnin Burn-in sample 
#'   @param thinning Thinning
#'
#' Returns:
#' @return median, mean, standard deviation and upper and lower 95\% confidence 
#'         intervals for the biodiversity parameter theta and the immigration rates. 
#'         Parameters x summaries matrix.
#'
#' @details The confidence intervals are calculated as the 2.5\% and 97.5\% quantiles.
#'
#' @export
#'
#' @references See citation("NMGS") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities
#'
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

#' Determine if community appears neutral. Calculate the
#' proportion of sample A that exceeds the median of sample B.
#'
#' Arguments:
#'   @param stats Statistics of the samples communities as in 
#'          the read_nmgs_s output.
#'   @param model Model to test: "local" or "full"
#'
#' Returns:
#'   @return A vector: med1 med2 n1 nT p, where med1 is the median for the first value, med2 second, n1 is the number of samples of the first variable that exceed the median of the second, nT is the total number of samples and p the proportion that exceed i.e. what we take as the pseudo p-value.
#'
#' @export
#'
#' @references See citation("NMGS") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

nmgs_neutrality <- function (stats, model) {

  if (model == "local") {i <- "LL"}
  if (model == "full") {i <- "LN"}
  j <- "LO"

  # Medians for the two samples
  med1 <- median(stats[, i])
  med2 <- median(stats[, j])

  # The number of samples of the first variable that exceed the second
  n1 <- sum(stats[, i] > med2)

  # Total number of samples
  nT <- nrow(stats)

  # Proportion that exceed i.e. what we take as the pseudo p-value.
  p <- n1/nT

  data.frame(list(median1 = med1, median2 = med2, n1 = n1, n = nT, pseudo.pvalue = p, model = I(model)))

}


#' Average over the metapopulations.
#'
#' Arguments:
#'   @param metacommunity
#'
#' Returns:
#'   @return The average
#'
#' @export
#'
#' @references See citation("NMGS") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

nmgs_metapopulation_average <- function (metacommunity) {

  # Average over species proportions across the local assembly samples			    
  pmeans <- colMeans(metacommunity$p)

  # TODO for the full neutral model
  # metacommunity$q

  list(local = pmeans)

}




