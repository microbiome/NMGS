#' @title Read in NMGS posterior samples
#' @description Simulation_out.csv: Contains the MCMC sampled parameters in format:
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
read_nmgs <- function (file, burnin = 0) {

  #An example generate Gibbs samples from fitting to the simulated data set
  #./NMGS -in Simulation.csv -out Simulation_out -v -s
  #-in specifies the input abundances across samples
  #-out the output file stub
  #-v prints progress to screen
  #-s generates samples to determine whether community appears neutral and output sampled metacommunity relative abundances

  # ------------------------------------------------------------------------------------

  # Produces three output files

  # Simulation_out.csv
  # Simulation_out_m.csv
  # Simulation_out_s.csv

# -------------------------------------------------------------------------------------

  # iteration: MCMC iteration
  # theta: the biodiversity parameter
  #  i1,...,iN are the immigration rates to each sample
  parameter.samples <- read.csv(file, header = FALSE)
  colnames(parameter.samples) <- c("iteration", "theta", paste("i", 1:(ncol(parameter.samples)-2), sep = ""))

  # Remove burn in samples
  parameter.samples = parameter.samples[(burnin+1):nrow(parameter.samples),]

  parameter.samples

}


