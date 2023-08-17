\name{mean_abs_change_scalar_ou}
\alias{mean_abs_change_scalar_ou}
\title{
Compute the expected absolute change of an Ornstein-Uhlenbeck process.
}
\description{
Given a scalar Ornstein-Uhlenbeck process at stationarity, compute the expected absolute net change over a specific time interval. In other words, if \eqn{X(t)} is the process at time \eqn{t}, compute the conditional expectation of \eqn{|X(t)-X(0)|} given that \eqn{X(0)} is randomly drawn from the stationary distribution. This quantity may be used as a measure for the speed at which a continuous trait evolves over time.
}
\usage{
mean_abs_change_scalar_ou(stationary_mean,
                          stationary_std,
                          decay_rate,
                          delta,
                          rel_error = 0.001,
                          Nsamples = NULL)
}
\arguments{
\item{stationary_mean}{
Numeric, the stationary mean of the OU process, i.e., its equilibrium (\eqn{\mu}).
}
\item{stationary_std}{
Positive numeric, the stationary standard deviation of the OU process. Note that this is \eqn{\sigma/\sqrt{2\lambda}}, where \eqn{\sigma} is the volatility.
}
\item{decay_rate}{
Positive numeric, the decay rate or "rubber band" parameter of the OU process (\eqn{\lambda}).
}
\item{delta}{
Positive numeric, the time step for which to compute the expected absolute change.
}
\item{rel_error}{
Positive numeric, the relative tolerable standard estimation error (relative to the true mean absolute displacement).
}
\item{Nsamples}{
Integer, number of Monte Carlo samples to use for estimation. If \code{NULL}, this is determined automatically based on the desired accuracy (\code{rel_error}).
}
}
\details{
The scalar OU process is a continuous-time stochastic process that satisfies the following stochastic differential equation:
\deqn{
dX = \lambda\cdot(\mu-X)\ dt + \sigma\ dW,
}
where \eqn{W} is a Wiener process (aka. "standard Brownian Motion"), \eqn{\mu} is the equilibrium, \eqn{\sigma} is the volatility and \eqn{\lambda} is the decay rate. The OU process is commonly used to model the evolution of a continuous trait over time. The decay rate \eqn{\lambda} alone is not a proper measure for how fast a trait changes over time (despite being erroneously used for this purpose in some sudies), as it only measures how fast the trait tends to revert to \eqn{\mu} when it is far away from \eqn{\mu}. Similarly, the volatility \eqn{\sigma} alone is not a proper measure of evolutionary rate, because it only describes how fast a trait changes when it is very close to the equilibrium \eqn{\mu}, where the tendency to revert is negligible and the process behaves approximately as a Brownian Motion.

This function uses Monte Carlo integration to estimate the expected absolute change, by repeatedly sampling start values \eqn{X(0)} from the OU's stationary distribution, computing the conditional expected absolute change given the sampled start value, and then averaging those conditional expectations.
}

\value{
A non-negative numeric, specifying the expected absolute change of the OU process.
}


\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{simulate_ou_model}}
}

\examples{
# compute the expected absolute change of an OU process after 10 time units
expected_abs_change = mean_abs_change_scalar_ou(stationary_mean=5, 
                                                stationary_std=1,
                                                decay_rate=0.1,
                                                delta=10)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Ornstein-Uhlenbeck}
