\name{simulate_deterministic_hbd}
\alias{simulate_deterministic_hbd}
\title{
Simulate a deterministic homogenous birth-death model.
}
\description{
Given a homogenous birth-death (HBD) model, i.e., with speciation rate \eqn{\lambda}, extinction rate \eqn{\mu} and sampling fraction \eqn{\rho}, calculate various deterministic features of the model backwards in time, such as the total diversity over time. The speciation and extinction rates may be time-dependent. ``Homogenous'' refers to the assumption that, at any given moment in time, all lineages exhibit the same speciation/extinction rates (in the literature this is sometimes referred to simply as ``birth-death model''; Morlon et al. 2011). ``Deterministic'' refers to the fact that all calculated properties are completely determined by the model's parameters (i.e. non-random), as if an infinitely large tree was generated (aka. ``continuum limit'').
Alternatively to \eqn{\lambda}, one may provide the pulled diversification rate (PDR; Louca et al. 2018) and the present-day speciation rate \eqn{\lambda(0)}. In either case, the time-profiles of \eqn{\lambda}, \eqn{\mu} or the PDR are specified as piecewise polynomial functions (splines), defined on a discrete grid of ages.
}
\usage{
simulate_deterministic_hbd(Ntips,
                           oldest_age,
                           rho            = 1,
                           age_grid       = NULL,
                           lambda         = NULL,
                           mu             = NULL,
                           PDR            = NULL,
                           lambda0        = NULL,
                           splines_degree = 1,
                           relative_dt    = 1e-3)
}
\arguments{
\item{Ntips}{
The assumed number of sampled extant species at present, defining the initial condition for the simulation. If the HBD model is supposed to describe a specific timetree, then Ntips should correspond to the number of tips in the tree.
}
\item{oldest_age}{
Strictly positive numeric, specifying the oldest time before present (``age'') to include in the simulation.
}
\item{rho}{
Numeric between 0 (exclusive) and 1 (inclusive), specifying the sampling fraction \eqn{\rho} of the tree, i.e. the fraction of extant species included (aka. "rarefaction"). Note that if \eqn{rho<1}, species are assumed to have been sampled randomly at equal probabilities. Can also be \code{NULL}, in which case \eqn{\rho=1} is assumed.
}
\item{age_grid}{
Numeric vector, listing discrete ages (time before present) on which either \eqn{\lambda} and \eqn{\mu}, or the PDR and \eqn{\mu}, are specified. Listed ages must be strictly increasing, and must cover at least the full considered age interval (from 0 to \code{oldest_age}). Can also be \code{NULL} or a vector of size 1, in which case the speciation rate, extinction rate and PDR are assumed to be time-independent.
}
\item{lambda}{
Numeric vector, of the same size as \code{age_grid} (or size 1 if \code{age_grid==NULL}), listing speciation rates (\eqn{\lambda}, in units 1/time) at the ages listed in \code{age_grid}. Speciation rates should be non-negative, and are assumed to vary polynomially between grid points (see argument \code{splines_degree}). If \code{NULL}, then \code{PDR} and \code{lambda0} must be provided.
}
\item{mu}{
Numeric vector, of the same size as \code{age_grid} (or size 1 if \code{age_grid==NULL}), listing extinction rates (\eqn{\mu}, in units 1/time) at the ages listed in \code{age_grid}. Extinction rates should be non-negative, and are assumed to vary polynomially between grid points (see argument \code{splines_degree}). Note that \code{mu} must always be provided.
}
\item{PDR}{
Numeric vector, of the same size as \code{age_grid} (or size 1 if \code{age_grid==NULL}), listing pulled diversification rates (in units 1/time) at the ages listed in \code{age_grid}. PDRs can be negative or positive, and are assumed to vary polynomially between grid points (see argument \code{splines_degree}). If \code{NULL}, then \code{lambda} must be provided.
}
\item{lambda0}{
Non-negative numeric, specifying the present-day extinction rate (in units 1/time). Either \code{lambda0} or \code{lambda} must be provided, but not both.
}
\item{splines_degree}{
Integer, either 0,1,2 or 3, specifying the polynomial degree of the provided \code{lambda}, \code{mu} and \code{PDR} between grid points in \code{age_grid}. For example, if \code{splines_degree==1}, then the provided \code{lambda}, \code{mu} and \code{PDR} are interpreted as piecewise-linear curves; if \code{splines_degree==2} they are interpreted as quadratic splines; if \code{splines_degree==3} they are interpreted as cubic splines. The \code{splines_degree} influences the analytical properties of the curve, e.g. \code{splines_degree==1} guarantees a continuous curve, \code{splines_degree==2} guarantees a continuous curve and continuous derivative, and so on.
}
\item{relative_dt}{
Strictly positive numeric (unitless), specifying the maximum relative time step allowed for integration over time. Smaller values increase integration accuracy but increase computation time. Typical values are 0.0001-0.001. The default is usually sufficient.
}
}
\details{
This function supports two alternative parameterizations of HBD models, either using the speciation and extinction rates (\eqn{\lambda}, \eqn{\mu}), or using the extinction rate, the pulled diversification rate (PDR) and the present-day speciation rate \eqn{\lambda(0)}. The PDR is defined as \eqn{PDR = \lambda-\mu+\lambda^{-1}d\lambda/d\tau}, where \eqn{\tau} is age. To avoid ambiguities, only one set of parameters is accepted, i.e. you must either provide \code{lambda} and \code{mu}, or provide \code{mu}, \code{PDR} and \code{lambda0}. The sampling fraction \code{rho} should be provided in either case; setting it to \code{NULL} is equivalent to setting it to 1.
}

\value{
A named list with the following elements:
\item{success}{
Logical, indicating whether the calculation was successful. If \code{FALSE}, then the returned list includes an additional `\code{error}' element (character) providing a description of the error; all other return variables may be undefined.
}
\item{ages}{Numerical vector of size NG, listing discrete ages (time before present) on which all returned time-curves are specified. Listed ages will be in ascending order, will cover exactly the range 0 - \code{oldest_age}, may be irregularly spaced, and may be finer than the original provided \code{age_grid}. Note that \code{ages[1]} corresponds to the present, while \code{ages[NG]} corresponds to the oldest time point (\code{oldest_age}).}
\item{total_diversity}{Numerical vector of size NG, listing the predicted (deterministic) total diversity (number of extant species, denoted \eqn{N}) at the ages given in \code{ages[]}.}
\item{shadow_diversity}{Numerical vector of size NG, listing the predicted (deterministic) ``shadow diversity'' at the ages given in \code{ages[]}. The shadow diversity is defined as \eqn{N_s=N\cdot \rho\lambda(0)/\lambda}.}
\item{Pmissing}{Numeric vector of size NG, listing the probability that a lineage, extant at a given age, will be absent from the tree either due to extinction or due to incomplete sampling.}
\item{Pextinct}{Numeric vector of size NG, listing the probability that a lineage, extant at a given age, will be fully extinct at present. Note that always \code{Pextinct<=Pmissing}.}
\item{LTT}{Numeric vector of size NG, listing the number of lineages represented in the tree at any given age, also known as ``lineages-through-time'' (LTT) curve. Note that \code{LTT[1]} will be equal to Ntips, and that values in \code{LTT} will be non-increasing with age.}
\item{lambda}{Numeric vector of size NG, listing the speciation rate (in units 1/time) at the ages given in \code{ages[]}.}
\item{mu}{Numeric vector of size NG, listing the extinction rate (in units 1/time) at the ages given in \code{ages[]}.}
\item{diversification_rate}{Numeric vector of size NG, listing the net diversification rate (\eqn{\lambda-\mu}) at the ages given in \code{ages[]}.}
\item{PDR}{Numeric vector of size NG, listing the pulled diversification rate (PDR, in units 1/time) at the ages given in \code{ages[]}.}
\item{PND}{Numeric vector of size NG, listing the pulled normalized diversity (PND, in units 1/time) at the ages given in \code{ages[]}. The PND is defined as \eqn{PND=(N/N(0))\cdot\lambda(0)/\lambda}.}
\item{SER}{Numeric vector of size NG, listing the ``shadow extinction rate'' (SER, in units 1/time) at the ages given in \code{ages[]}. The SER is defined as \eqn{SER=\rho\lambda(0)-PDR}.}
\item{PER}{Numeric vector of size NG, listing the ``pulled extinction rate'' (PER, in units 1/time) at the ages given in \code{ages[]}. The PER is defined as \eqn{SER=\lambda(0)-PDR} (Louca et al. 2018).}
\item{rholambda0}{Non-negative numeric, specifying the product of the sampling fraction and the present-day speciation rate, \eqn{\rho\cdot\lambda(0)}.}
}


\author{Stilianos Louca}

\references{
H. Morlon, T. L. Parsons, J. B. Plotkin (2011). Reconciling molecular phylogenies with the fossil record. Proceedings of the National Academy of Sciences. 108:16327-16332.

S. Louca et al. (2018). Bacterial diversification through geological time. Nature Ecology & Evolution. 2:1458-1467.
}

\seealso{
\code{\link{loglikelihood_hbd}}
}

\examples{
# define an HBD model with exponentially decreasing speciation/extinction rates
Ntips     = 1000
beta      = 0.01 # exponential decay rate of lambda over time
oldest_age= 10
age_grid  = seq(from=0,to=oldest_age,by=0.1) # choose a sufficiently fine age grid
lambda    = 1*exp(beta*age_grid) # define lambda on the age grid
mu        = 0.2*lambda # assume similarly shaped but smaller mu

# simulate deterministic HBD model
simulation = simulate_deterministic_hbd(Ntips, 
                                        oldest_age, 
                                        rho       = 0.5,
                                        age_grid  = age_grid,
                                        lambda    = lambda,
                                        mu        = mu)

# plot deterministic LTT
plot( x = simulation$ages, y = simulation$LTT, type='l',
      main='dLTT', xlab='age', ylab='lineages')
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{tree model}
\keyword{simulation}
\keyword{homogenous birth-death model}