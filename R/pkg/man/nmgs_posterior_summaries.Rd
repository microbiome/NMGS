\name{nmgs_posterior_summaries}
\alias{nmgs_posterior_summaries}
\title{Description: Compute summary statistics of the sampled MCMC parameters,
excluding the burn-in samples and thinning.}
\usage{
  nmgs_posterior_summaries(samples, burnin = 0,
    thinning = 1)
}
\arguments{
  \item{samples}{Posterior samples of the NMGS model
  (output from the read_nmgs function)}

  \item{burnin}{Burn-in sample}

  \item{thinning}{Thinning

  Returns:}
}
\value{
  median, mean, standard deviation and upper and lower 95\%
  confidence intervals for the biodiversity parameter theta
  and the immigration rates.  Parameters x summaries
  matrix.
}
\description{
  Arguments:
}
\details{
  The confidence intervals are calculated as the 2.5\% and
  97.5\% quantiles.
}
\author{
  Contact: Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("NMGS")
}
\keyword{utilities}

