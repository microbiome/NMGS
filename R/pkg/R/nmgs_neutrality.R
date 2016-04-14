#' @title Determine if community appears neutral
#' @description Determine if community appears neutral
#'
#' @param stats Statistics of the samples communities as in 
#'          the read_nmgs_s output.
#' @param model Model to test: "local" or "full"
#'
#' @return A vector: med1 med2 n1 nT p, where med1 is the median for the first value, med2 second, n1 is the number of samples of the first variable that exceed the median of the second, nT is the total number of samples and p the proportion that exceed i.e. what we take as the pseudo p-value.
#'
#' @export
#'
#' @references See citation("NMGS") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities
nmgs_neutrality <- function (stats, model) {

  if (model == "local") {i <- "LL"} # local community sample
  if (model == "full")  {i <- "LN"} # full neutral sample

  j <- "LO" # Local observed sample

  # Medians for the two samples
  med1 <- median(stats[, i])
  med2 <- median(stats[, j])

  # The number of samples of the neutral model sample likelihood exceed observed sample likelihood
  n1 <- sum(stats[, i] > med2)

  # Total number of samples
  nT <- nrow(stats)

  # Proportion that exceed i.e. what we take as the pseudo p-value.
  p <- n1/nT

  data.frame(list(median1 = med1, median2 = med2, n1 = n1, n = nT, pseudo.pvalue = p, model = I(model)))

}
