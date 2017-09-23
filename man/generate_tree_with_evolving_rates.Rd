\name{generate_tree_with_evolving_rates}
\alias{generate_tree_with_evolving_rates}
\title{
Generate a random tree with birth/death rates evolving under Brownian motion.
}
\description{
Generate a random phylogenetic tree via simulation of a Poissonian speciation/extinction (birth/death) process. New species are added (born) by splitting of a randomly chosen extant tip. Per-capita birth and death rates evolve under a Brownian motion model constrained to a finite interval (via reflection). Thus, the probability rate of a tip splitting or going extinct depends on the tip, with closely related tips having more similar per-capita birth and death rates.
}
\usage{
generate_tree_with_evolving_rates(parameters           = list(),
                                  max_tips             = NULL, 
                                  max_time             = NULL,
                                  max_time_eq          = NULL,
                                  coalescent           = TRUE,
                                  as_generations       = FALSE,
                                  Nsplits              = 2, 
                                  tip_basename         = "", 
                                  node_basename        = NULL,
                                  include_birth_times  = FALSE,
                                  include_death_times  = FALSE,
                                  include_rates        = FALSE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{parameters}{
A named list specifying the birth-death model parameters, with one or more of the following entries:

\code{birth_rate_diffusivity}: Non-negative number. Diffusivity constant for the Brownian motion model of the evolving per-capita birth rate. In units 1/time^3. See \code{\link{simulate_bm_model}} for an explanation of the diffusivity parameter.

\code{min_birth_rate_pc}:
Non-negative number. The minimum allowed per-capita birth rate of a clade. In units 1/time. By default this is 0.

\code{max_birth_rate_pc}:
Non-negative number. The maximum allowed per-capita birth rate of a clade. In units 1/time. By default this is 1.

\code{death_rate_diffusivity}: Non-negative number. Diffusivity constant for the Brownian motion model of the evolving per-capita death rate. In units 1/time^3. See \code{\link{simulate_bm_model}} for an explanation of the diffusivity parameter.

\code{min_death_rate_pc}:
Non-negative number. The minimum allowed per-capita death rate of a clade. In units 1/time. By default this is 0.

\code{max_death_rate_pc}:
Non-negative number. The maximum allowed per-capita death rate of a clade. In units 1/time. By default this is 1.

\code{rarefaction}:
Numeric between 0 and 1. Rarefaction to be applied at the end of the simulation (fraction of random tips kept in the tree). 
Note that if \code{coalescent==FALSE}, rarefaction may remove both extant as well as extinct clades. Set \code{rarefaction=1} to not perform any rarefaction.
}
\item{max_tips}{
Maximum number of tips of the tree to be generated. If \code{coalescent=TRUE}, this refers to the number of extant tips. Otherwise, it refers to the number of extinct + extant tips. If \code{NULL} or <=0, the number of tips is unlimited (so be careful).
}
\item{max_time}{
Maximum duration of the simulation. If \code{NULL} or <=0, this constraint is ignored.
}
\item{max_time_eq}{
Maximum duration of the simulation, counting from the first point at which speciation/extinction equilibrium is reached, i.e. when (birth rate - death rate) changed sign for the first time. If \code{NULL} or <0, this constraint is ignored.
}
\item{coalescent}{
Logical, specifying whether only the coalescent tree (i.e. the tree spanning the extant tips) should be returned. If \code{coalescent==FALSE} and the death rate is non-zero, then the tree may include non-extant tips (i.e. tips whose distance from the root is less than the total time of evolution). In that case, the tree will not be ultrametric.
}
\item{as_generations}{
Logical, specifying whether edge lengths should correspond to generations. If FALSE, then edge lengths correspond to time.
}
\item{Nsplits}{
Integer greater than 1. Number of child-tips to generate at each diversification event. If set to 2, the generated tree will be bifurcating. If >2, the tree will be multifurcating.
}
\item{tip_basename}{
Character. Prefix to be used for tip labels (e.g. "tip."). If empty (""), then tip labels will be integers "1", "2" and so on.
}
\item{node_basename}{
Character. Prefix to be used for node labels (e.g. "node."). If \code{NULL}, no node labels will be included in the tree.
}
\item{include_birth_times}{
Logical. If \code{TRUE}, then the times of speciation events (in order of occurrence) will also be returned.
}
\item{include_death_times}{
Logical. If \code{TRUE}, then the times of extinction events (in order of occurrence) will also be returned.
}
\item{include_rates}{
Logical. If \code{TRUE}, then the bper-capita birth & death rates of all tips and nodes will also be returned.
}
}


\details{
If \code{max_time==NULL}, then the returned tree will always contain \code{max_tips} tips. In particular, if at any moment during the simulation the tree only includes a single extant tip, the death rate is temporarily set to zero to prevent the complete extinction of the tree. If \code{max_tips==NULL}, then the simulation is ran as long as specified by \code{max_time}. If neither \code{max_time} nor \code{max_tips} is \code{NULL}, then the simulation halts as soon as the time exceeds \code{max_time} or the number of tips (extant tips if \code{coalescent} is \code{TRUE}) exceeds \code{max_tips}. If \code{max_tips!=NULL} and \code{Nsplits>2}, then the last diversification even may generate fewer than \code{Nsplits} children, in order to keep the total number of tips within the specified limit.

The per-capita birth and death rates of the root are chosen randomly, independently and uniformly from their respective allowed intervals.
}


\value{
A named list with the following elements:
\item{success}{
Logical, indicating whether the simulation was successful. If \code{FALSE}, an additional element \code{error} (of type character) is included containing an explanation of the error; in that case the value of any of the other elements is undetermined.
}
\item{tree}{
A rooted bifurcating (if \code{Nsplits==2}) or multifurcating (if \code{Nsplits>2}) tree of class "phylo", generated according to the specified birth/death model. 

If \code{coalescent==TRUE} or if all death rates are zero, and only if \code{as_generations==FALSE}, then the tree will be ultrametric. If \code{as_generations==TRUE} and \code{coalescent==FALSE}, all edges will have unit length.
}
\item{root_time}{
Numeric, giving the time at which the tree's root was first split during the simulation. 
Note that if \code{coalescent==TRUE}, this may be later than the first speciation event during the simulation.
}
\item{final_time}{
Numeric, giving the final time at the end of the simulation. If \code{coalescent==TRUE}, then this may be greater than the total time span of the tree (since the root of the coalescent tree need not correspond to the first speciation event).
}
\item{equilibrium_time}{
Numeric, giving the first time where the sign of (death rate - birth rate) changed from the beginning of the simulation, i.e. when speciation/extinction equilibrium was reached. May be infinite if the simulation stoped before reaching this point.
}
\item{Nbirths}{
Total number of birth events (speciations) that occurred during tree growth. This may be lower than the total number of tips in the tree if death rates were non-zero and \code{coalescent==TRUE}, or if \code{Nsplits>2}.
}
\item{Ndeaths}{
Total number of deaths (extinctions) that occurred during tree growth.
}
\item{birth_times}{
Numeric vector, listing the times of speciation events during tree growth, in order of occurrence. Note that if \code{coalescent==TRUE}, then \code{speciation_times} may be greater than the phylogenetic distance to the coalescent root.
}
\item{death_times}{
Numeric vector, listing the times of extinction events during tree growth, in order of occurrence. Note that if \code{coalescent==TRUE}, then \code{speciation_times} may be greater than the phylogenetic distance to the coalescent root.
}
\item{birth_rates_pc}{
Numeric vector, listing the per-capita birth rate of each tip and node in the tree. The length of an edge in the tree was thus drawn from an exponential distribution with rate equal to the per-capita birth rate of the child tip or node.
}
\item{death_rates_pc}{
Numeric vector, listing the per-capita death rate of each tip and node in the tree.
}
}

\author{Stilianos Louca}

\references{
D. J. Aldous (2001). Stochastic models and descriptive statistics for phylogenetic trees, from Yule to today. Statistical Science. 16:23-34.
}

%\seealso{
%}

\examples{
# Generate tree
parameters = list(birth_rate_diffusivity  = 1,
                  min_birth_rate_pc       = 1,
                  max_birth_rate_pc       = 2,
                  death_rate_diffusivity  = 0.5,
                  min_death_rate_pc       = 0,
                  max_death_rate_pc       = 1)
simulation = generate_tree_with_evolving_rates(parameters,max_tips=1000,include_rates=TRUE)
tree       = simulation$tree
Ntips      = length(tree$tip.label)

# plot per-capita birth & death rates of tips
plot( x=simulation$birth_rates_pc[1:Ntips], 
      y=simulation$death_rates_pc[1:Ntips], 
      type='p', 
      xlab="pc birth rate", 
      ylab="pc death rate", 
      main="Per-capita birth & death rates across tips",
      las=1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% Only 1 keyword per line
\keyword{random}
\keyword{speciation}
\keyword{tree model}
