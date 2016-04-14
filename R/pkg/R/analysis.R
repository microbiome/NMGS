#' @title Average over the metapopulations
#' @description Average over the metapopulations
#' @param metacommunity metacommunity
#' @return The average
#' @export
#' @references See citation("NMGS") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities
nmgs_metapopulation_average <- function (metacommunity) {

  # Average over species proportions across the local assembly samples			    
  pmeans <- colMeans(metacommunity$p)

  # TODO for the full neutral model
  # metacommunity$q

  list(local = pmeans, full = "TODO")

}




