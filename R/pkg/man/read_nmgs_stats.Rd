\name{read_nmgs_stats}
\alias{read_nmgs_stats}
\title{Read in Simulation_out_s.csv:
Gives statistics on the sampled communities in format:
SampleN,LN,LL,LO,HN,HL,HO,SN,SL,SO}
\usage{
read_nmgs_stats(file)
}
\arguments{
  \item{file}{NMGS output file name

  Returns:}
}
\value{
the statistics matrix (see above)
}
\description{
where SampleN is nth sample generated under the neutral
model with fitted parameters taken from the corresponding
MCMC sample. LN,LL,LO are the log-likelihoods of the full
neutral sample, the local community sample and the observed
sample respectively. Then HN,HL,HO and SN,SL,SO are the
species entropies and richness's in the same order.
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

