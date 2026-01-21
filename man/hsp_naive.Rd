\name{hsp_naive}
\alias{hsp_naive}
\title{
Naive hidden state prediction.
}
\description{
Reconstruct ancestral discrete states of nodes and predict unknown (hidden) states of tips on a tree based on empirical state proportions across tips on the entire tree. This method completely ignores any phylogenetic relationships, and is not designed for making reliable state predictions. It is only meant for testing and debugging purposes, e.g. as a "null" model!}
\usage{
hsp_naive(tree, 
          tip_states, 
          Nstates     = NULL,
          check_input = TRUE)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". Note that the tree structure is not actually used internally.
}

\item{tip_states}{
An integer vector of size Ntips, specifying the state of each tip in the tree as an integer from 1 to Nstates, where Nstates is the possible number of states (see below). \code{tip_states} can include \code{NA} to indicate an unknown tip state that is to be predicted.
}
\item{Nstates}{
Either \code{NULL}, or an integer specifying the number of possible states of the trait. If \code{NULL}, then it will be computed based on the maximum non-\code{NA} value encountered in \code{tip_states}
}
\item{check_input}{
Logical, specifying whether to perform some basic checks on the validity of the input data. If you are certain that your input data are valid, you can set this to \code{FALSE} to reduce computation.
}
}


\details{
For this function, the trait's states must be represented by integers within 1,..,Nstates, where Nstates is the total number of possible states. If the states are originally in some other format (e.g. characters or factors), you should map them to a set of integers 1,..,Nstates. You can easily map any set of discrete states to integers using the function \code{\link{map_to_state_space}}. Any \code{NA} entries in \code{tip_states} are interpreted as unknown states. 

The function calculates the "global" empirical proportions of known states, i.e., across the entire tree, and then sets the state likelihoods of tips with unknown state and nodes to this distribution. States are "predicted" for each tip with unknown state and each node by randomly drawing states according to the same global empirical distribution. This function has asymptotic time complexity O(Ntips x Nstates).

Tips must be represented in \code{tip_states} in the same order as in \code{tree$tip.label}. The vector \code{tip_states} need not include names; if it does, however, they are checked for consistency (if \code{check_input==TRUE}).

This function is meant for reconstructing ancestral states in all nodes of a tree as well as predicting the states of tips with an a priory unknown state, according to the "null" model where phylogenetic relationships are ignored but empirical frequencies of known states are still accounted for.
}


\value{
A list with the following elements:
\item{success}{
Logical, indicating whether HSP was successful. If \code{FALSE}, some return values may be \code{NULL}.
}
\item{likelihoods}{
A 2D numeric matrix, listing the probability of each tip and node being in each state. This matrix will have (Ntips+Nnodes) rows and Nstates columns, where Nstates was either explicitly provided as an argument or set to the maximum value found in \code{tip_states} (if \code{Nstates} was passed as NULL). The rows in this matrix will be in the order in which tips and nodes are indexed in the tree, i.e. the rows 1,..,Ntips store the probabilities for tips, while rows (Ntips+1),..,(Ntips+Nnodes) store the probabilities for nodes. Each row in this matrix will sum up to 1. Note that the return value is named this way for compatibility with other HSP functions.
}
\item{states}{
Integer vector of length Ntips+Nnodes, with values in \{1,..,Nstates\}, specifying the "predicted" state for each tip & node. For nodes and tips with a priori unknown state, the predicted state is simply a random draw from the empirical distribution.
}
}


\author{Stilianos Louca}

\references{
J. R. Zaneveld and R. L. V. Thurber (2014). Hidden state prediction: A modification of classic ancestral state reconstruction algorithms helps unravel complex symbioses. Frontiers in Microbiology. 5:431.
}

\seealso{
\code{\link{hsp_max_parsimony}}, 
\code{\link{hsp_mk_model}}, 
\code{\link{hsp_empirical_probabilities}}
\code{\link{hsp_mk_model}}, 
}

\examples{\dontrun{
# generate random tree
Ntips = 100
tree = generate_random_tree(list(birth_rate_intercept=1),max_tips=Ntips)$tree

# simulate a discrete trait
Nstates = 5
Q = get_random_mk_transition_matrix(Nstates, rate_model="ER", max_rate=0.1)
tip_states = simulate_mk_model(tree, Q)$tip_states

# print states of first 20 tips
print(tip_states[1:20])

# set half of the tips to unknown state
tip_states[sample.int(Ntips,size=as.integer(Ntips/2),replace=FALSE)] = NA

# reconstruct all tip states using the naive method
HSP = hsp_naive(tree, tip_states=tip_states, Nstates=Nstates)

# print estimated states of first 20 tips
print(HSP$states[1:20])
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{maximum parsimony}
\keyword{ancestral state reconstruction}
\keyword{Sankoff's algorithm}
\keyword{hidden state prediction}
