\name{nmgs_median}
\alias{nmgs_median}
\title{Description: Compute medians from the sampled MCMC parameters,
excluding the burn-in samples. For the biodiversity parameter theta
and the immigration rates.}
\usage{
  nmgs_median(samples, burnin)
}
\arguments{
  \item{samples}{Posterior samples of the NMGS model. As in
  the output from the read_nmgs function.}

  \item{burnin}{Burn-in sample

  Returns:}
}
\value{
  median values with upper and lower 95\% confidence
  intervals (lower / median / upper)
}
\description{
  Arguments:
}
\details{
  The cofidence intervals are calculated as the 2.5\% and
  97.5\% quantiles.
}
\author{
  Contact: Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("NMGS")
}
\keyword{utilities}

