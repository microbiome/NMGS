\name{nmgs_neutrality}
\alias{nmgs_neutrality}
\title{Determine if community appears neutral. Calculate the
proportion of sample A that exceeds the median of sample B.}
\usage{
  nmgs_neutrality(stats, model)
}
\arguments{
  \item{stats}{Statistics of the samples communities as in
  the read_nmgs_s output.}

  \item{model}{Model to test: "local" or "full"

  Returns:}
}
\value{
  A vector: med1 med2 n1 nT p, where med1 is the median for
  the first value, med2 second, n1 is the number of samples
  of the first variable that exceed the median of the
  second, nT is the total number of samples and p the
  proportion that exceed i.e. what we take as the pseudo
  p-value.
}
\description{
  Arguments:
}
\author{
  Contact: Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("NMGS")
}
\keyword{utilities}

