\name{get_ancestral_nodes}
\alias{get_ancestral_nodes}
\title{
Compute ancestral nodes.
}
\description{
Given a rooted phylogenetic tree and a set of tips and/or nodes ("descendants"), determine their ancestral node indices, traveling a specific number of splits back in time (i.e., towards the root). Optionally, one may restrict the back-traversal to a specific subset of nodes.
}
\usage{
get_ancestral_nodes(tree,
                    descendants,
                    Nsplits     = 1,
                    only_nodes  = NULL,
                    omit_nodes  = NULL)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo".
}
\item{descendants}{
An integer vector or character vector, specifying the tips/nodes for each of which to determine the ancestral node. If an integer vector, it must list indices of tips (from 1 to Ntips) and/or nodes (from Ntips+1 to Ntips+Nnodes), where Ntips and Nnodes is the number of tips and nodes in the tree, respectively. If a character vector, it must list tip and/or node names. In this case \code{tree} must include \code{tip.label}, as well as \code{node.label} if nodes are included in \code{descendants}.
}
\item{Nsplits}{
Either a single integer or an integer vector of the same length as \code{descendants}, with values >=1, specifying how many splits to travel backward at most. For example, Nsplits=1 will yield the parent node of each tip/node in \code{descendants} whenever possible. Note that this is an upper bound on the allowed number of backward steps; if fewer backward steps are possible (e.g., if a descendant is close to the root), then the number of backward steps taken will be less than \code{Nsplits}.
}
\item{only_nodes}{
Optional integer vector with values in 1,..,Nnodes, specifying nodes to restrict the back-traversal to. Other nodes will be skipped and not counted towards \code{Nsplits}. If \code{NULL} or empty, this filter is ignored.
}
\item{omit_nodes}{
Optional integer vector with values in 1,..,Nnodes, specifying nodes to exlude from the back-traversal. These nodes will be skipped and not counted towards \code{Nsplits}.
}
}
\details{
The tree may include multi-furcations (i.e. nodes with more than 2 children) as well as mono-furcations (i.e. nodes with only one child).
}
\value{
An integer vector of the same length as \code{descendants}, with values in 1,..,Nnodes, listing the node indices representing the ancestors of \code{descendants} traveling backward up to \code{Nsplits}.
}


\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{get_pairwise_mrcas}}, \code{\link{get_mrca_of_set}}
}

\examples{
# generate a random tree
tree = generate_random_tree(list(birth_rate_intercept=1),
                            max_tips = 50,
                            tip_basename = "tip.",
                            node_basename = "node.")$tree

# pick 3 tips
descendants=c("tip.5", "tip.7","tip.10")

# determine the immediate parent node of each tip
ancestors = castor::get_ancestral_nodes(tree, descendants, Nsplits=1)
print(tree$node.label[ancestors])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ancestor}
