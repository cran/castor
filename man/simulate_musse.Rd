\name{simulate_musse}
\alias{simulate_musse}
\title{
Simulate a Multiple State Speciation and Extinction (MuSSE) model.
}
\description{
Generate a random phylogenetic tree via simulation of a Poissonian speciation/extinction (birth/death) process, whereby birth and death rates are determined by a co-evolving discrete trait. New species are added (born) by splitting of a randomly chosen extant tip. The discrete trait, whose values determine birth/eath rates, evolves according to a discrete-space continuous-time Markov process with fixed transition rates between states. This is the Multiple State Speciation and Extinction (MuSSE) model described by FitzJohn et al. (2009), as an extension to the Binary State Speciation and Extinction (BiSSE) model by Maddison et al. (2007). Optionally, the model can be turned into a Hidden State Speciation and Extinction model (Beaulieu and O'meara, 2016), by replacing the simulated tip/node states with "proxy" states, thus hiding the original states actually influencing speciation/extinction rates.
}
\usage{
simulate_musse(Nstates,
               NPstates             = NULL,
               proxy_map            = NULL,
               parameters           = list(),
               root_state           = NULL,
               max_tips             = NULL, 
               max_time             = NULL,
               max_time_eq          = NULL,
               sampling_fractions   = NULL,
               reveal_fractions     = NULL,
               coalescent           = TRUE,
               as_generations       = FALSE,
               all_Mk_transitions   = TRUE,
               no_full_extinction   = TRUE,
               Nsplits              = 2, 
               tip_basename         = "", 
               node_basename        = NULL,
               include_birth_times  = FALSE,
               include_death_times  = FALSE,
               include_rates        = FALSE,
               include_labels       = TRUE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{Nstates}{
Integer, specifying the number of possible discrete states a tip can have, influencing speciation/extinction rates. For example, if \code{Nstates==2} then this corresponds to the common Binary State Speciation and Extinction (BiSSE) model (Maddison et al., 2007). In the case of a HiSSE model, \code{Nstates} refers to the total number of diversification rate categories, as described by Beaulieu and O'meara (2016).
}
\item{NPstates}{
Integer, optionally specifying a number of "proxy-states" that are observed instead of the underlying speciation/extinction-modulating states. To simulate a HiSSE model, this should be smaller than \code{Nstates}. Each state corresponds to a different proxy-state, as defined using the variable \code{proxy_map} (see below). For BiSSE/MuSSE with no hidden states, \code{NPstates} can be set to either \code{NULL} or equal to \code{Nstates}, and proxy-states are equivalent to states.
}
\item{proxy_map}{
Integer vector of size \code{Nstates} and with values in 1,..\code{NPstates}, specifying the correspondence between states (i.e. diversification-rate categories) and (observed) proxy-states, in a HiSSE model. Specifically, \code{proxy_map[s]} indicates which proxy-state the state s is represented by. Each proxy-state can represent multiple states (i.e. proxies are ambiguous), but each state must be represented by exactly one proxy-state. For non-HiSSE models, set this to \code{NULL}. See below for more details.
}
\item{parameters}{
A named list specifying the MuSSE/HiSSE model parameters, such as the transition rates between states and the state-dependent birth/death rates (see details below).
}
\item{root_state}{
Integer within 1,..,\code{Nstates}, specifying the initial state of the root. If left unspecified, this is chosen randomly and uniformly among all possible states.
}
\item{max_tips}{
Maximum number of tips in the generated tree, prior to any subsampling. If \code{coalescent=TRUE}, this refers to the number of extant tips, prior to subsampling. Otherwise, it refers to the number of extinct + extant tips, prior to subsampling. If \code{NULL} or <=0, the number of tips is not limited, so you should use \code{max_time} and/or \code{max_time_eq} to stop the simulation.
}
\item{max_time}{
Maximum duration of the simulation. If \code{NULL} or <=0, this constraint is ignored.
}
\item{max_time_eq}{
Maximum duration of the simulation, counting from the first point at which speciation/extinction equilibrium is reached, i.e. when (birth rate - death rate) changed sign for the first time. If \code{NULL} or <0, this constraint is ignored.
}
\item{sampling_fractions}{
A single number, or a numeric vector of size \code{NPstates}, listing tip sub-sampling fractions, depending on proxy-state. \code{sampling_fractions[p]} is the probability of including a tip in the final tree, if its proxy-state is p. If \code{NULL}, all tips (or all extant tips, if \code{coalescent==TRUE}) are included in the tree. If a single number, all tips are included with the same probability, i.e. regardless of their proxy-state.
}
\item{reveal_fractions}{
Numeric vector of size \code{NPstates}, listing reveal fractions of tip proxy-states, depending on proxy state. \code{reveal_fractions[p]} is the probability of knowing a tip's proxy-state, if its proxy state is p. Can also be NULL, in which case all tip proxy states will be known.
}
\item{coalescent}{
Logical, specifying whether only the coalescent tree (i.e. the tree spanning the extant tips) should be returned. If \code{coalescent==FALSE} and the death rate is non-zero, then the tree may include non-extant tips (i.e. tips whose distance from the root is less than the total time of evolution). In that case, the tree will not be ultrametric.
}
\item{as_generations}{
Logical, specifying whether edge lengths should correspond to generations. If FALSE, then edge lengths correspond to time.
}
\item{all_Mk_transitions}{
Logical, specifying whether to simulate all state transitions along edges. If \code{FALSE}, then state transitions are only simulated (i.e. states are only updated) during speciation events. This should typically be \code{TRUE}, however it may increase computation time if transition rates are comparable to or greater than speciation rates.
}
\item{no_full_extinction}{
Logical, specifying whether to prevent complete extinction of the tree. Full extinction is prevented by temporarily disabling extinctions whenever the number of extant tips is 1. if \code{no_full_extinction==FALSE} and death rates are non-zero, the tree may go extinct during the simulation; if \code{coalescent==TRUE}, then the returned tree would be empty, hence the function will return unsuccessfully (i.e. \code{success} will be \code{FALSE}). By default \code{no_full_extinction} is \code{TRUE}, however in some special cases it may be desirable to allow full extinctions to ensure that the generated trees are statistically distributed exactly according to the underlying cladogenic model.
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
Logical. If \code{TRUE}, then the per-capita birth & death rates of all tips and nodes will also be returned.
}
\item{include_labels}{
Logical, specifying whether to include tip-labels and node-labels (if available) as names in the returned state vectors (e.g. \code{tip_states} and \code{node_states}). In any case, returned states are always listed in the same order as tips and nodes in the tree. Setting this to \code{FALSE} may increase computational efficiency for situations where labels are not required.
}
}


\details{
Birth/death (speciation/extinction) rates are determined by a tip's current "state", which evolves according to a continuous-time discrete-state Markov chain (Mk model) with constant transition rates. The argument \code{parameters} should be a named list including one or more of the following elements:
\itemize{
\item\code{birth_rates}: Numeric vector of size Nstates, listing the per-capita birth rate (speciation rate) at each state. Can also be a single number (all states have the same birth rate).
\item\code{death_rates}: Numeric vector of size Nstates, listing the per-capita death rate (extinction rate) at each state. Can also be a single number (all states have the same death rate).
\item\code{transition_matrix}: 2D numeric matrix of size Nstates x Nstates. Transition rate matrix for the Markov chain model of birth/death rate evolution.
}
If \code{max_time==NULL} and \code{max_time_eq==NULL}, then the returned tree will always contain \code{max_tips} tips. In particular, if at any moment during the simulation the tree only includes a single extant tip, the death rate is temporarily set to zero to prevent the complete extinction of the tree. If \code{max_tips==NULL}, then the simulation is ran as long as specified by \code{max_time} and/or \code{max_time_eq}. If neither \code{max_time}, \code{max_time_eq} nor \code{max_tips} is \code{NULL}, then the simulation halts as soon as the time exceeds \code{max_time}, or the time since equilibration exceeds \code{max_time_eq}, or the number of tips (extant tips if \code{coalescent} is \code{TRUE}) exceeds \code{max_tips}, whichever occurs first. If \code{max_tips!=NULL} and \code{Nsplits>2}, then the last diversification even may generate fewer than \code{Nsplits} children, in order to keep the total number of tips within the specified limit.

HiSSE models (Beaulieu and O'meara, 2016) are closely related to BiSSE/MuSSE models, the main difference being the fact that the actual diversification-modulating states are not directly observed. Hence, this function is also able to simulate HiSSE models, with appropriate choice of the input variables \code{Nstates}, \code{NPstates} and \code{proxy_map}. For example, \code{Nstates=4}, \code{NPstates=2} and \code{proxy_map=c(1,2,1,2)} specifies that states 1 and 3 are represented by proxy-state 1, and states 2 and 4 are represented by proxy-state 2. This is the original case described by Beaulieu and O'meara (2016); in their terminology, there would be 2 "hidden"" states ("0" and "1") and 2 "observed" (proxy) states ("A" and "B"), and the 4 diversification rate categories (\code{Nstates=4}) would be called "0A", "1A", "0B" and "1B", respectively. The somewhat different terminology used here allows for easier generalization to an arbitrary number of diversification-modulating states and an arbitrary number of proxy states. For example, if there are 6 diversification modulating states, represented by 3 proxy-states as 1->A, 2->A, 3->B, 4->C, 5->C, 6->C, then one would set \code{Nstates=6}, \code{NPstates=3} and \code{proxy_map=c(1,1,2,3,3,3)}.
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
Numeric, giving the first time where the sign of (death rate - birth rate) changed from the beginning of the simulation, i.e. when speciation/extinction equilibrium was reached. May be infinite if the simulation stopped before reaching this point.
}
\item{Nbirths}{
Numeric vector of size Nstates, listing the total number of birth events (speciations) that occurred at each state. The sum of all entries in \code{Nbirths} may be lower than the total number of tips in the tree if death rates were non-zero and \code{coalescent==TRUE}, or if \code{Nsplits>2}.
}
\item{Ndeaths}{
Numeric vector of size Nstates, listing the total number of death events (extinctions) that occurred at each state.
}
\item{Ntransitions}{
2D numeric matrix of size Nstates x Nstates, listing the total number of transition events that occurred between each pair of states. For example, \code{Ntransitions[1,2]} is the number of transitions that occured from state 1 to state 2.
}
\item{tip_states}{
Integer vector of size Ntips and with values in 1,..,Nstates, listing the state of each tip in the tree.
}
\item{node_states}{
Integer vector of size Nnodes and with values in 1,..,Nstates, listing the state of each node in the tree.
}
\item{tip_proxy_states}{
Integer vector of size Ntips and with values in 1,..,NPstates, listing the proxy state of each tip in the tree. Only included in the case of HiSSE models.
}
\item{node_proxy_states}{
Integer vector of size Nnodes and with values in 1,..,NPstates, listing the proxy state of each node in the tree. Only included in the case of HiSSE models.
}
\item{root_state}{
Integer, specifying the initial state of the root (either provided during the function call, or generated randomly).
}
\item{birth_times}{
Numeric vector, listing the times of speciation events during tree growth, in order of occurrence. Note that if \code{coalescent==TRUE}, then \code{speciation_times} may be greater than the phylogenetic distance to the coalescent root.
}
\item{death_times}{
Numeric vector, listing the times of extinction events during tree growth, in order of occurrence. Note that if \code{coalescent==TRUE}, then \code{speciation_times} may be greater than the phylogenetic distance to the coalescent root.
}
\item{clade_birth_rates}{
Numeric vector of size Ntips+Nnodes, listing the per-capita birth rate of each tip and node in the tree. Only included if \code{include_rates==TRUE}.
}
\item{clade_death_rates}{
Numeric vector of size Ntips+Nnodes, listing the per-capita death rate of each tip and node in the tree. Only included if \code{include_rates==TRUE}.
}
}

\author{Stilianos Louca}

\references{
W. P. Maddison, P. E. Midford, S. P. Otto (2007). Estimating a binary character's effect on speciation and extinction. Systematic Biology. 56:701-710.

R. G. FitzJohn, W. P. Maddison, S. P. Otto (2009). Estimating trait-dependent speciation and extinction rates from incompletely resolved phylogenies. Systematic Biology. 58:595-611

R. G. FitzJohn (2012). Diversitree: comparative phylogenetic analyses of diversification in R. Methods in Ecology and Evolution. 3:1084-1092

J. M. Beaulieu and B. C. O'Meara (2016). Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. Systematic Biology. 65:583-601.
}

\seealso{
\code{\link{fit_musse}}
}

\examples{
# Simulate a tree under a BiSSE model
Q = get_random_mk_transition_matrix(Nstates=2, rate_model="ER", max_rate=0.1)
parameters = list(birth_rates       = c(1,1.5),
                  death_rates       = 0.5,
                  transition_matrix = Q)
simulation = simulate_musse(Nstates       = 2, 
                            parameters    = parameters, 
                            max_tips      = 1000, 
                            include_rates = TRUE)
tree       = simulation$tree
Ntips      = length(tree$tip.label)

# plot distribution of per-capita birth rates of tips
rates = simulation$clade_birth_rates[1:Ntips]
barplot(table(rates)/length(rates), 
        xlab="rate", 
        main="Distribution of pc birth rates across tips (BiSSE model)")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% Only 1 keyword per line
\keyword{random}
\keyword{tree model}
\keyword{MuSSE model}
\keyword{BiSSE model}
