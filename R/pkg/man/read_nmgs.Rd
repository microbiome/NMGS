\name{read_nmgs}
\alias{read_nmgs}
\title{Read in NMGS posterior samples}
\usage{
  read_nmgs(file)
}
\arguments{
  \item{file}{NMGS output file name

  Returns:}
}
\value{
  Posterior samples of the NMGS parameters
}
\description{
  Simulation_out.csv: Contains the MCMC sampled parameters
  in format: iter1,theta,i1,..,iN .. iterM,theta,i1,..,iN
  where theta is the biodiversity parameter and i1,...,iN
  are the immigration rates to each sample
}
\details{
  Arguments:
}
\author{
  Contact: Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("NMGS")
}
\keyword{utilities}

