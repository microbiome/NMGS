#' @title Read in Simulation_out_s.csv
#' @description Gives statistics on the sampled communities.
#' @details Format:
#' SampleN,LN,LL,LO,HN,HL,HO,SN,SL,SO
#' where SampleN is nth sample generated under the neutral model with
#' fitted parameters taken from the corresponding MCMC sample.
#' LN,LL,LO are the log-likelihoods of the full neutral sample, the local
#' community sample and the observed sample respectively. Then HN,HL,HO
#' and SN,SL,SO are the species diversities (entropies) and richnesses.
#' @param file NMGS output file name
#' @return the statistics matrix (see above)
#' @export
#' @references See citation("NMGS") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities
read_nmgs_stats <- function (file) {

  stats <- read.csv(file, header = FALSE)
  colnames(stats) <- c("SampleN","LN","LL","LO","HN","HL","HO","SN","SL","SO")	    
  stats	    

}
