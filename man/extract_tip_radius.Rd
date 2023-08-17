\name{extract_tip_radius}
\alias{extract_tip_radius}
\title{
Extract a subtree spanning tips within a certain radius.
}
\description{
Given a rooted tree, a focal tip and a phylogenetic (patristic) distance radius, extract a subtree comprising all tips located within the provided radius from the focal tip.
}
\usage{
extract_tip_radius(tree, 
                   focal_tip,
                   radius,
                   include_subtree = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{tree}{
A rooted tree of class "phylo".
}
\item{focal_tip}{
Either a character, specifying the name of the focal tip, or an integer between 1 and Ntips, specifying the focal tip's index.
}
\item{radius}{
Non-negative numeric, specifying the patristic distance radius to consider around the focal tip.
}
\item{include_subtree}{
Logical, whether to actually extract the subtree, rather than just returning the inferred within-radius tips.
}
}


\details{
The "patristic distance" between two tips and/or nodes is the shortest cumulative branch length that must be traversed along the tree in order to reach one tip/node from the other. If \code{tree$edge.length} is missing, then each edge is assumed to be of length 1.

The tree may include multi-furcations as well as mono-furcations (i.e. nodes with only one child). The input tree must be rooted at some node for technical reasons (see function \code{\link{root_at_node}}), but the choice of the root node does not influence which tips are extracted.
}


\value{
A named list with the following elements:
\item{radius_tips}{
Integer vector with values in 1,..,Ntips, specifying which tips were found to be within the specified radius of the focal tip.
}
\item{subtree}{
A new tree of class "phylo", containing only the tips within the specified radius from the focal tip, and the nodes & edges connecting those tips to the root. Only returned if \code{include_subtree=TRUE}.
}
\item{new2old_tip}{
Integer vector of length Ntips_extracted (=number of tips in the extracted subtree, i.e., within the specified distance radius from the focal tip) with values in 1,..,Ntips, mapping tip indices of the extracted subtree to tip indices in the original tree. In particular, \code{tree$tip.label[new2old_tip]} will be equal to \code{subtree$tip.label}. Only returned if \code{include_subtree=TRUE}.
}
}

\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{get_all_distances_to_tip}}
}

\examples{
# generate a random tree
Ntips = 50
tree  = generate_random_tree(list(birth_rate_factor=0.1),
                             max_tips = Ntips,
                             tip_basename="tip.")$tree

# extract all tips within radius 50 from a focal tip
subtree = extract_tip_radius(tree, 
                             focal_tip="tip.39",
                             radius=50)$subtree
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% Only 1 keyword per line
\keyword{pruning}
\keyword{subtree}
\keyword{phylogenetic distance}
