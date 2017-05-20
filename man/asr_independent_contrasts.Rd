\name{asr_independent_contrasts}
\alias{asr_independent_contrasts}
\title{
Ancestral state reconstruction via phylogenetic independent contrasts.
}
\description{
Reconstruct ancestral states for a continuous (numeric) trait using phylogenetic independent contrasts (PIC; Felsenstein, 1985).
}
\usage{
asr_independent_contrasts(tree, tip_states,  weighted=TRUE, check_input=TRUE)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{tip_states}{
A numeric vector of size Ntips, specifying the known state of each tip in the tree.
}
\item{weighted}{
Logical, specifying whether to weight tips and nodes by the inverse length of their incoming edge, as in the original method by Felsenstein (1985). If \code{FALSE}, edge lengths are treated as if they were 1.
}
\item{check_input}{
Logical, specifying whether to perform some basic checks on the validity of the input data. If you are certain that your input data are valid, you can set this to \code{FALSE} to reduce computation.
}
}


\details{
The function traverses the tree in postorder (tips-->root) and estimates the state of each node as a convex combination of the estimated states of its chilren. These estimates are the intermediate "X" variables introduced by Felsenstein (1985) in his phylogenetic independent contrasts method. For the root, this yields the same globally parsimonious state as the squared-changes parsimony algorithm implemented in \code{asr_squared_change_parsimony} (Maddison 1991). For any other node, PIC only yields locally parsimonious reconstructions, i.e. reconstructed states only depend on the subtree descending from the node (see discussion by Maddison 1991).

If \code{weighted==TRUE}, then this function yields the same ancestral state reconstructions as

\code{ape::ace(phy=tree, x=tip_states, type="continuous", method="pic", model="BM", CI=FALSE)}

in the \code{ape} package (v. 0.5-64).

If \code{tree$edge.length} is missing, each edge in the tree is assumed to have length 1. This is the same as setting \code{weighted=FALSE}. The tree may include multi-furcations (i.e. nodes with more than 2 children) as well as mono-furcations (i.e. nodes with only one child). Edges with length 0 will be adjusted internally to some tiny length if needed (if \code{weighted==TRUE}).

Tips must be represented in \code{tip_states} in the same order as in \code{tree$tip.label}. The vector \code{tip_states} need not include item names; if it does, however, they are checked for consistency (if \code{check_input==TRUE}). All tip states must be non-NA; otherwise, consider using one of the functions for hidden-state-prediction (e.g., \code{\link{hsp_independent_contrasts}}).

The function has asymptotic time complexity O(Nedges).
}


\value{
A list with the following elements:
\item{ancestral_states}{
A numeric vector of size Nnodes, listing the reconstructed state of each node. The entries in this vector will be in the order in which nodes are indexed in the tree.
}
\item{total_sum_of_squared_changes}{
The total sum of squared changes in tree, minimized by the (optionally weighted) squared-change parsimony algorithm. This is equation 7 in (Maddison, 1991).
}
}


\author{Stilianos Louca}

\references{
J. Felsenstein (1985). Phylogenies and the Comparative Method. The American Naturalist. 125:1-15.

W. P. Maddison (1991). Squared-change parsimony reconstructions of ancestral states for continuous-valued characters on a phylogenetic tree. Systematic Zoology. 40:304-314.
}

\seealso{
\code{\link{asr_squared_change_parsimony}},
\code{\link{asr_max_parsimony}},
\code{\link{asr_mk_model}}
}

\examples{
# generate random tree
Ntips = 100
tree = generate_random_tree(max_tips=Ntips, birth_rate_intercept=1);

# simulate a continuous trait
tip_states = simulate_ou_model(tree, stationary_mean=0, spread=1, decay_rate=0.001)$tip_states

# reconstruct node states via weighted PIC
node_states = asr_independent_contrasts(tree, tip_states, weighted=TRUE)$ancestral_states
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{maximum parsimony}
\keyword{ancestral state reconstruction}
\keyword{BM model}
