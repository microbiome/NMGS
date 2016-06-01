#An example generate Gibbs samples from fitting to the simulated data set
#./NMGS -in Simulation.csv -out Simulation_out -v -s
#-in specifies the input abundances across samples
#-out the output file stub
#-v prints progress to screen
#-s generates samples to determine whether community appears neutral and output sampled metacommunity relative abundances

# ------------------------------------------------------------------------------------

#Produces three output files

#Simulation_out.csv
#Simulation_out_m.csv
#Simulation_out_s.csv

# -------------------------------------------------------------------------------------

#' Read in NMGS posterior samples
#' 
#' Simulation_out.csv: Contains the MCMC sampled parameters in format:
#' iter1,theta,i1,..,iN
#' ..
#' iterM,theta,i1,..,iN
#' where theta is the biodiversity parameter and i1,...,iN are the immigration rates to each sample
#'
#' Arguments:
#'   @param file NMGS output file name
#'
#' Returns:
#'   @return Posterior samples of the NMGS parameters
#'
#' @export
#'
#' @references See citation("NMGS") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

read_nmgs <- function (file) {

  # iteration: MCMC iteration
  # theta: the biodiversity parameter
  #  i1,...,iN are the immigration rates to each sample
  parameter.samples <- read.csv(file)
  colnames(parameter.samples) <- c("iteration", "theta", paste("i", 1:(ncol(parameter.samples)-2), sep = ""))

  parameter.samples

}


#' Read in NMGS generated metacommunities file
#' Simulation_out_m.csv: Gives the generated metacommunities with format:
#'   SampleN,SL,SN
#'   p1,...,pSL
#'   q1,...,qSN
#'
#' where SampleN is nth sample generated from the neutral model with
#' parameters from the corresponding MCMC sample. SL and SN are the
#' number of species in the local assembly sample (always S) and the
#' sample generated under the full neutral model respectively, p and q
#' are then the metacommunity distributions
#'
#' Arguments:
#'   @param file NMGS output file name
#'
#' Returns:
#'   @return Generated metacommunities from the NMGS model. A list with the following elements:
#'     	     N: A matrix. Each row corresponds to a sample. SampleN,SL,SN for each sample (see above)
#'	     p: metacommunity distribution (local assemply; see above)
#'	     q: metacommunity distribution (sample generated under full neutral model; see above)
#'
#' @export
#'
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

    k <- (i-1)*3

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


#' Read in Simulation_out_s.csv: 
#' Gives statistics on the sampled communities in format:
#' SampleN,LN,LL,LO,HN,HL,HO,SN,SL,SO
#'
#' where SampleN is nth sample generated under the neutral model with
#' fitted parameters taken from the corresponding MCMC sample. LN,LL,LO
#' are the log-likelihoods of the full neutral sample, the local
#' community sample and the observed sample respectively. Then HN,HL,HO
#' and SN,SL,SO are the species entropies and richness's in the same
#' order.
#'
#' Arguments:
#'   @param file NMGS output file name
#'
#' Returns:
#'   @return the statistics matrix (see above)
#'
#' @export
#'
#' @references See citation("NMGS") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

read_nmgs_stats <- function (file) {

  stats <- read.csv(file, header = FALSE)
  colnames(stats) <- c("SampleN","LN","LL","LO","HN","HL","HO","SN","SL","SO")	    
  stats	    

}


