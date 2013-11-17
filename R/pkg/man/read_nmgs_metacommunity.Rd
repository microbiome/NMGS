\name{read_nmgs_metacommunity}
\alias{read_nmgs_metacommunity}
\title{Read in NMGS generated metacommunities file
Simulation_out_m.csv: Gives the generated metacommunities with format:
  SampleN,SL,SN
  p1,...,pSL
  q1,...,qSN}
\usage{
  read_nmgs_metacommunity(file)
}
\arguments{
  \item{file}{NMGS output file name

  Returns:}
}
\value{
  Generated metacommunities from the NMGS model. A list
  with the following elements: N: A matrix. Each row
  corresponds to a sample. SampleN,SL,SN for each sample
  (see above) p: metacommunity distribution (local
  assemply; see above) q: metacommunity distribution
  (sample generated under full neutral model; see above)
}
\description{
  where SampleN is nth sample generated from the neutral
  model with parameters from the corresponding MCMC sample.
  SL and SN are the number of species in the local assembly
  sample (always S) and the sample generated under the full
  neutral model respectively, p and q are then the
  metacommunity distributions
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

