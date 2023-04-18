\name{get_independent_sister_tips}
\alias{get_independent_sister_tips}
\title{
Extract disjoint tip pairs with independent relationships.
}
\description{
Given a rooted tree, extract disjoint pairs of sister tips with disjoint connecting paths. These tip pairs can be used to compute independent contrasts of numerical traits evolving according to an arbitrary time-reversible process.
}
\usage{
get_independent_sister_tips(tree)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo".
}
}

\details{
If the input tree only contains monofurcations and bifurcations (recommended), it is guaranteed that at most one unpaired tip will be left (i.e., if Ntips was odd).
}


\value{
An integer matrix of size NP x 2, where NP is the number of returned tip pairs. Each row in this matrix lists the indices (from 1 to Ntips) of a tip pair that can be used to compute a phylogenetic independent contrast.
}


\author{Stilianos Louca}

\references{
J. Felsenstein (1985). Phylogenies and the Comparative Method. The American Naturalist. 125:1-15.

T. Garland Jr., P. E. Midford, A. R. Ives (1999). An introduction to phylogenetically based statistical methods, with a new method for confidence intervals on ancestral values. American Zoologist. 39:374-388. 
}

\seealso{
\code{\link{get_independent_contrasts}}
}

\examples{
# generate random tree
Ntips = 100
tree = generate_random_tree(list(birth_rate_intercept=1),Ntips)$tree

# simulate a continuous trait on the tree
tip_states = simulate_bm_model(tree, diffusivity=0.1, include_nodes=FALSE)$tip_states

# get independent tip pairs
tip_pairs = get_independent_sister_tips(tree)

# calculate non-scaled PICs for the trait using the independent tip pairs
PICs = tip_states[tip_pairs]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{BM model}
