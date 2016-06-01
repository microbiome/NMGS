#' @title Read in NMGS generated metacommunities file
#' @description Simulation_out_m.csv gives the generated metacommunities with format:
#'   SampleN,SL,SN
#'   p1,...,pSL
#'   q1,...,qSN
#' where SampleN is nth sample generated from the neutral model with
#' parameters from the corresponding MCMC sample. SL and SN are the
#' number of species in the local assembly sample (always S) and the
#' sample generated under the full neutral model respectively, p and q
#' are then the metacommunity distributions.
#'
#' @param file NMGS output file name
#' @return Generated metacommunities from the NMGS model. A list with the following elements:
#'     	     N: A matrix. Each row corresponds to a sample. SampleN,SL,SN for each sample (see above)
#'	     p: metacommunity distribution (local assemply; see above)
#'	     q: metacommunity distribution (sample generated under full neutral model; see above)
#' @export
#' @references See citation("NMGS") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities
read_nmgs_metacommunity <- function ( file ) {

  dat <- readLines(file)

  # Separate the information to distinct, matched matrices
  N <- array(NA, dim = c(length(dat)/3, 3))
  colnames(N) <- c("SampleN", "SL", "SN")

  p <- list() # array(NA, dim = c(length(dat)/3, 3))
  q <- list() # array(NA, dim = c(length(dat)/3, 3))
  
  for (i in seq(1, length(dat)/3)) {

    k <- (i-1) * 3

    # SampleN,SL,SN
    N[i, ] <- as.numeric(unlist(strsplit(dat[k + 1], ",")))
    
    #' p1,...,pSL
    p[[i]] <- as.numeric(unlist(strsplit(dat[k + 2], ",")))
    
    #' q1,...,qSN
    q[[i]] <- as.numeric(unlist(strsplit(dat[k + 3], ",")))
  
  }

  # Convert p into a matrix (as it has always the same number of samples)
  p <- t(sapply(p, identity))

  list(N = N, p = p, q = q)

}


