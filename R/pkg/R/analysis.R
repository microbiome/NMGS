#' Description: Compute medians from the sampled MCMC parameters,
#' excluding the burn-in samples. For the biodiversity parameter theta
#' and the immigration rates. 
#'
#' Arguments:
#'   @param samples Posterior samples of the NMGS model. As in the output from the read_nmgs function.
#'   @param burnin Burn-in sample 
#'
#' Returns:
#'   @return median values with upper and lower 95\% confidence intervals (lower / median / upper)
#'
#' @details The cofidence intervals are calculated as the 2.5\% and 97.5\% quantiles.
#'
#' @export
#'
#' @references See citation("NMGS") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

nmgs_median <- function (samples, burnin) {

  # Exclude the burn-in samples	and return median for all parameters
  med <- apply(samples[samples$iteration > burnin, -1], 2, median)
  lc <- quantile(samples[samples$iteration > burnin, 2], 0.025)
  uc <- quantile(samples[samples$iteration > burnin, 2], 0.975)

  rbind(lower = lc, median = med, upper = uc)

}



#' Determine if community appears neutral. Calculate the
#' proportion of sample A that exceeds the median of sample B.
#'
#' Arguments:
#'   @param stats Statistics of the samples communities as in the read_nmgs_s output.
#'   @param i index of sample A
#'   @param j index of sample B
#'
#' Returns:
#'   @return A vector: med1 med2 n1 nT p, where med1 is the median for the first value, med2 second, n1 is the number of samples of the first variable that exceed the median of the second, nT is the total number of samples and p the proportion that exceed i.e. what we take as the pseudo p-value.
#'
#' @export
#'
#' @references See citation("NMGS") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

nmgs_neutrality <- function (stats, i, j) {

  # Medians for the two samples
  med1 <- median(stats[, i])
  med2 <- median(stats[, j])

  # The number of samples of the first variable that exceed the second
  n1 <- sum(stats[, i] > med2)

  # Total number of samples
  nT <- nrow(stats)

  # Proportion that exceed i.e. what we take as the pseudo p-value.
  p <- n1/nT

  c(median1 = med1, median2 = med2, n1 = n1, n = nT, pseudo.pvalue = p)

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

  # TODO also for the full neutral model
  # metacommunity$q

  list(local = pmeans)

}




