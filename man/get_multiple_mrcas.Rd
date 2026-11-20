\name{get_multiple_mrcas}
\alias{get_multiple_mrcas}
\title{
Get most recent common ancestors of one or more sets of tips/nodes.
}
\description{
Given a rooted phylogenetic tree and one or more sets of tips and/or nodes, for each set of tips/nodes find the most recent common ancestor (MRCA).
This is the multi-set version of \code{get_mrca_of_set}, and tends to be more efficient than calling \code{get_mrca_of_set} many times on the same tree.
}
\usage{
get_multiple_mrcas(tree, descendants)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{descendants}{
List of vectors, each specifying a set of tips/nodes to compute the MRCA of. Each \code{descendants[s]} must either be a character vector (listing tip/node names) or an integer vector (listing indices of tips/nodes). If an integer vector, it must list indices of tips (from 1 to Ntips) and/or nodes (from Ntips+1 to Ntips+Nnodes).
}
}
\details{
The tree may include multi-furcations as well as mono-furcations (i.e. nodes with only one child).
If \code{tree$edge.length} is missing, then each edge is assumed to be of length 1. 
Note that in some cases the MRCA of two tips may be a tip, namely when both tips are the same.

If some of the elements of \code{descendants} are character vectors, then \code{tree$tip.label} must exist. If node names are included in \code{descendants}, then \code{tree$node.label} must also exist.
}
\value{
An integer vector of size Nsets with values in 1,..,Ntips (tips) and/or in (Ntips+1),..,(Ntips+Nnodes) (nodes), with the s-th element being the index of the MRCA of tips/nodes listed in \code{descendants[s]}.
}


\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{get_mrca_of_set}},
\code{\link{get_pairwise_mrcas}},
\code{\link{get_tips_for_mrcas}}
}

\examples{
# generate a random tree
Ntips = 100
tree = generate_random_tree(list(birth_rate_intercept=1),Ntips)$tree

# pick 3 random sets of tips or nodes
Nsets = 3
descendants = list(sample.int(n=(Ntips+tree$Nnode), size=3, replace=FALSE),
                sample.int(n=(Ntips+tree$Nnode), size=5, replace=FALSE),
                sample.int(n=(Ntips+tree$Nnode), size=2, replace=FALSE))

# calculate the MRCA for each set
MRCAs = get_multiple_mrcas(tree, descendants)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{MRCA}
