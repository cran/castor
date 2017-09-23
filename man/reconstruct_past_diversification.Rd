\name{reconstruct_past_diversification}
\alias{reconstruct_past_diversification}
\title{
Reconstruct past diversification dynamics from a diversity time series.
}
\description{
Given a time series of past diversities (coalescent or not), this function estimates instantaneous birth (speciation) and death (extinction) rates that would lead to the observed diversity time series. The function is based on a deterministic model, in which instantaneous birth and death rates lead to a predictable growth of a tree (one new species per birth event). The reconstruction is non-parametric. Of particular importance is the case where the time series is coalescent, i.e. represents the diversity that would be visible from a coalescent tree generated throuhg a Poissonian birth-death process (see \code{\link{generate_random_tree}}).
}
\usage{
reconstruct_past_diversification( times,
                                  diversities,
                                  birth_rates_pc  = NULL,
                                  rarefaction     = 1.0,
                                  max_age         = 0,
                                  coalescent      = FALSE,
                                  smoothing_span  = 0,
                                  smoothing_order = 1)
}
\arguments{
\item{times}{
Numeric vector, listing the times at which diversities are given. Values must be in ascending order.
}
\item{diversities}{
Numeric vector of the same size as \code{times}, listing diversities (coalescent or not) at each time point.
}
\item{birth_rates_pc}{
Numeric vector of the same size as \code{times}, listing known or assumed per-capita birth rates. Can also be of size 1, in which case the same per-capita birth rate is assumed throughout. Alternatively if \code{coalescent==TRUE}, then this vector can also be empty, in which case a constant per-capita birth rate is assumed and estimated from the slope of the coalescent diversities at the last time point. The last alternative is not available when \code{coalescent==FALSE}.
}
\item{rarefaction}{
Numeric between 0 and 1.  Optional rarefaction fraction assumed for the diversities at the very end. Set to 1 to assume no rarefaction was performed.
}
\item{max_age}{
Numeric. Optional maximum distance from the end time to be considered. If <=0, all provided time points are considered.
}
\item{coalescent}{
Logical, indicating whether the provided diversities are from a coalescent tree (only including clades with extant representatives) or true diversities (extant species at each time point).
}
\item{smoothing_span}{
Non-negative integer. Optional sliding window size (number of time points) for smoothening the diversities time series via Savitzky-Golay-filter. If <=2, no smoothing is done. Smoothening the time series can reduce the effects of noise on the reconstructed diversity dynamics.
}
\item{smoothing_order}{
Integer between 1 and 4. Polynomial order of the Savitzky-Golay smoothing filter to be applied. Only relevant if \code{smoothing_span>2}. A value of 1 or 2 is typically recommended.
}
}


\details{
For a coalescent diversity time series \eqn{N_c(\tau)} at various ages \eqn{\tau} (distance from the end of the time series), reconstruction of the true diversities \eqn{N(\tau)} is based on the following formulas:
\deqn{
P(\tau)=1+\frac{\nu(\tau)}{\beta(\tau)},\\
}
\deqn{
N(\tau)=\frac{N_c}{1-P(\tau)},
}
\deqn{
\nu(\tau)=\frac{1}{N_c(\tau)}\frac{dN_c(\tau)}{d\tau}
}
where \eqn{P(\tau)} is the probability that a clade of size 1 at age \eqn{\tau} went extinct by the end of the time series and \eqn{\beta} is the per-capita birth rate. If the per-capita birth rate is not explicitly provided for each time point (see argument \code{birth_rate_pc}), the function assumes that the per-capita birth rate is constant at all times. If \code{birth_rates_pc==NULL} and \code{coalescent==TRUE}, the constant per-capita birth rate is estimated as
\deqn{
\beta = \frac{\nu(0)\rho}{1-2\rho-(1-\rho)^2},
}
where \eqn{\rho} is the fraction of species kept after rarefaction (see argument \code{rarefaction}).

The reconstruction is only accurate in the deterministic limit, i.e. for high diversities where the stochastic nature of Poissonian birth/death processes diminishes.
}


\value{
A named list with the following elements:
\item{success}{
Logical, specifying whether the reconstruction was successful. If \code{FALSE}, the remaining elements may not be defined.
}
\item{Ntimes}{
Integer. Number of time points for which reconstruction is returned.
}
\item{true_diversities}{
Numeric vector of the same size as \code{times}, listing the true diversities at each time point (number of extant species at each time point). If \code{coalescent==FALSE}, then these are the same as the \code{diversities} passed to the function.
}
\item{coalescent_diversities}{
Numeric vector of the same size as \code{times}, listing the coalescent diversities at each time point (number of species with at least one extant descendant at the last time point). If \code{coalescent==TRUE}, then these are the same as the \code{diversities} passed to the function.
}
\item{birth_rates}{
Numeric vector of the same size as \code{times}, listing the estimated birth rates (speciation events per time unit).
}
\item{death_rates}{
Numeric vector of the same size as \code{times}, listing the estimated death rates (extinction events per time unit).
}
\item{survival_chances}{
Numeric vector of the same size as \code{times}, listing the estimated survival probabilities of size-1 clades at each time point. \code{survival_chances[i]} is the probability that a clade of size 1 at time \code{times[i]} will be extant by the end of the time series.
}
\item{total_births}{
Numeric, giving the estimated total number of birth events that occurred between times \code{T-max_age} and \code{T}, where \code{T} is the last time point of the time series.
}
\item{total_deaths}{
Numeric, giving the estimated total number of death events that occurred between times \code{T-max_age} and \code{T}, where \code{T} is the last time point of the time series.
}
\item{last_birth_rate_pc}{
The provided or estimated (if not provided) per-capita birth rate at the last time point. This corresponds to the birth rate divided by the estimated true diversity (prior to rarefaction) at the last time point.
}
\item{last_death_rate_pc}{
The provided or estimated (if not provided) per-capita death rate at the last time point. This corresponds to the death rate divided by the estimated true diversity (prior to rarefaction) at the last time point.
}
}

\author{Stilianos Louca}

%\references{
%Louca et al (2017). In review
%}

\seealso{
\code{\link{generate_random_tree}},
%\code{\link{fit_tree_model}}
\code{\link{count_clades_over_time}}
}

\examples{
#####################################################
# EXAMPLE 1

# Generate a coalescent tree
params = list(birth_rate_intercept  = 0, 
              birth_rate_factor     = 1,
              birth_rate_exponent   = 1,
              death_rate_intercept  = 0,
              death_rate_factor     = 0.05,
              death_rate_exponent   = 1.3,
              rarefaction           = 1)
simulation = generate_random_tree(params,max_time_eq=1,coalescent=TRUE)
tree = simulation$tree
time_span = simulation$final_time - simulation$root_time
cat(sprintf("Generated tree has \%d tips, spans \%g time units\n",length(tree$tip.label),time_span))

# Calculate diversity time series from the tree
counter = count_clades_over_time(tree, times=seq(0,0.99*time_span,length.out=100))

# print coalescent diversities
print(counter$diversities)

# reconstruct diversification dynamics based on diversity time series
results = reconstruct_past_diversification( counter$times,
                                            counter$diversities,
                                            coalescent      = TRUE,
                                            smoothing_span  = 3,
                                            smoothing_order = 1)
                                            
# print reconstructed true diversities
print(results$true_diversities)
                                                  
# plot coalescent and reconstructed true diversities
matplot(x     = counter$times, 
        y     = matrix(c(counter$diversities,results$true_diversities), ncol=2, byrow=FALSE),
        type  = "b", 
        xlab  = "time", 
        ylab  = "# clades",
        lty   = c(1,2), pch = c(1,0), col = c("red","blue"))
legend( "topleft", 
        legend  = c("coalescent (simulated)","true (reconstructed)"), 
        col     = c("red","blue"), lty = c(1,2), pch = c(1,0));
        
        
        
#####################################################
# EXAMPLE 2

# Generate a non-coalescent tree
params = list(birth_rate_intercept  = 0, 
              birth_rate_factor     = 1,
              birth_rate_exponent   = 1,
              death_rate_intercept  = 0,
              death_rate_factor     = 0.05,
              death_rate_exponent   = 1.3,
              rarefaction           = 1)
simulation = generate_random_tree(params,max_time_eq=1,coalescent=FALSE)
tree = simulation$tree
time_span = simulation$final_time - simulation$root_time
cat(sprintf("Generated tree has \%d tips, spans \%g time units\n",length(tree$tip.label),time_span))

# Calculate diversity time series from the tree
counter = count_clades_over_time(tree, times=seq(0,0.99*time_span,length.out=100))

# print true diversities
print(counter$diversities)

# reconstruct diversification dynamics based on diversity time series
results = reconstruct_past_diversification( counter$times,
                                            counter$diversities,
                                            birth_rates_pc  = params$birth_rate_factor,
                                            coalescent      = FALSE,
                                            smoothing_span  = 3,
                                            smoothing_order = 1)
                                            
# print coalescent diversities
print(results$coalescent_diversities)
                                                  
# plot coalescent and reconstructed true diversities
matplot(x     = counter$times, 
        y     = matrix(c(results$coalescent_diversities,counter$diversities), ncol=2, byrow=FALSE),
        type  = "b", 
        xlab  = "time", 
        ylab  = "# clades",
        lty   = c(1,2), pch = c(1,0), col = c("red","blue"))
legend( "topleft", 
        legend  = c("coalescent (reconstructed)","true (simulated)"), 
        col     = c("red","blue"), lty = c(1,2), pch = c(1,0));
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% Only 1 keyword per line
\keyword{tree model}
\keyword{fitting}
