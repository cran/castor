\name{simulate_diversification_model}
\alias{simulate_diversification_model}
\title{
Simulate a speciation/extinction model.
}
\description{
Simulate a speciation/extinction model for diversity over time, corresponding to a tree generation model. Speciation (birth) and extinction (death) rates can each be constant or power-law functions of the number of extant species. For example,
\deqn{
B = I + F\cdot N^E,
}
where \eqn{B} is the birth rate, \eqn{I} is the intercept, \eqn{F} is the power-law factor, \eqn{N} is the current number of extant species and \eqn{E} is the power-law exponent. Note that currently only the deterministic version of the model can be simulated; for the stochastic version (Poisson process) see \code{\link{generate_random_tree}}.
}
\usage{
simulate_diversification_model( times,
                                parameters          = list(),
                                start_time          = NULL,
                                start_diversity	    = 1,
                                coalescent          = FALSE,
                                reverse             = FALSE,
                                include_event_rates = FALSE,
                                include_Nevents     = FALSE,
                                max_runtime         = NULL)
}
\arguments{
\item{times}{
Numeric vector, listing the times for which to calculate diversities, as predicted by the model. Values must be in ascending order.
}
\item{parameters}{
A named list specifying the birth-death model parameters, with one or more of the following entries:

\code{birth_rate_intercept}: Non-negative number. The intercept of the Poissonian rate at which new species (tips) are added. In units 1/time.

\code{birth_rate_factor}:
Non-negative number. The power-law factor of the Poissonian rate at which new species (tips) are added. In units 1/time.

\code{birth_rate_exponent}:
Numeric. The power-law exponent of the Poissonian rate at which new species (tips) are added. Unitless.

\code{death_rate_intercept}:
Non-negative number. The intercept of the Poissonian rate at which extant species (tips) go extinct. In units 1/time.

\code{death_rate_factor}:
Non-negative number. The power-law factor of the Poissonian rate at which extant species (tips) go extinct. In units 1/time.

\code{death_rate_exponent}:
Numeric. The power-law exponent of the Poissonian rate at which extant species (tips) go extinct. Unitless.

}
\item{start_time}{
Numeric. If \code{reverse==FALSE} (see below), then this is the start time of the tree (<=\code{times[1]}). If \code{reverse==TRUE}, then this is the final time of the tree (>=\code{max(times)}). Can also be \code{NULL}, in which case it is set to the first or last value in \code{times} (depending on \code{reverse}).
}
\item{start_diversity}{
Numeric. True diversity at \code{start_time}. For example, if \code{reverse==TRUE}, then this should be the final diversity of the tree.

}
\item{coalescent}{
Logical, specifying whether the diversity corresponding to a coalescent tree (i.e. the tree spanning only extant tips) should be calculated. If \code{coalescent==TRUE} and the death rate is non-zero, then the returned diversities will generally be lower than the number of extant species calculated by the model at each time point.
}
\item{reverse}{
Logical. If \code{TRUE}, then the tree model is simulated in backward time direction. In that case, \code{start_diversity} is interpreted as the known diversity at the last time point, and all diversities at previous time points are calculated based on the model. If \code{FALSE}, then the model is simulated in forward-time.
}
\item{include_event_rates}{
Logical. If \code{TRUE}, then the birth (speciation) and death (extinction) rates (for each time point) are included as returned values. This comes at a moderate computational overhead.
}
\item{include_Nevents}{
Logical. If \code{TRUE}, then the cumulative birth (speciation) and death (extinction) events (for each time point) are included as returned values. This comes at a moderate computational overhead.
}
\item{max_runtime}{
Numeric. Maximum runtime (in seconds) allowed for the simulation. If this time is surpassed, the simulation aborts.
}
}


\details{
In the special case where per-capita birth and death rates are constant (i.e. \eqn{I=0} and \eqn{E=1} for birth and death rates), this function uses an explicit analytical solution to the underlying differential equations, and is thus much faster than in the general case.
}


\value{
A named list with the following elements:
\item{success}{
Logical, indicating whether the simulation was successful. If the simulation aborted due to runtime constraints (option \code{max_runtime}), \code{success} will be \code{FALSE}.
}
\item{diversities}{
Numeric vector of the same size as \code{times}, listing the diversity (if \code{coalescent==FALSE}) or the subset of diversity seen in the coalescent tree (if \code{coalescent==TRUE}), for each time point in \code{times}.
}
\item{birth_rates}{
Numeric vector of the same size as \code{times}, listing the speciation (birth) rate at each time point. Only relevant if \code{include_event_rates==TRUE}.
}
\item{death_rates}{
Numeric vector of the same size as \code{times}, listing the extinction (death) rate at each time point. Only relevant if \code{include_event_rates==TRUE}.
}
\item{Nbirths}{
Numeric vector of the same size as \code{times}, listing the cumulative number of speciation (birth) events up to each time point. Only relevant if \code{include_Nevents==TRUE}.
}
\item{Ndeaths}{
Numeric vector of the same size as \code{times}, listing the cumulative number of extinction (death) events up to each time point. Only relevant if \code{include_Nevents==TRUE}.
}
}

\author{Stilianos Louca}


\seealso{
\code{\link{generate_random_tree}},
%\code{\link{fit_tree_model}}
\code{\link{count_clades_over_time}}
}

\examples{
# Generate a tree
max_time = 100
parameters = list(birth_rate_intercept  = 1, 
                  birth_rate_factor     = 0,
                  birth_rate_exponent   = 0,
                  death_rate_intercept  = 0,
                  death_rate_factor     = 0,
                  death_rate_exponent   = 0)
tree = generate_random_tree(parameters,max_time=max_time)$tree

# Calculate diversity-vs-time curve for the tree
times = seq(from=0,to=0.99*max_time,length.out=10)
tree_diversities = count_clades_over_time(tree, times=times)$diversities

# simulate diversity curve based on deterministic model
model_diversities = simulate_diversification_model(times,parameters)$diversities

# compare diversities in the tree to the simulated ones
plot(tree_diversities,model_diversities,xlab="tree diversities",ylab="simulated diversities")
abline(a=0,b=1,col="#A0A0A0") # show diagonal for reference
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% Only 1 keyword per line
\keyword{tree model}
