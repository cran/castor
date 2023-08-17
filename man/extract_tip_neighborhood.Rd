\name{extract_tip_neighborhood}
\alias{extract_tip_neighborhood}
\title{
Extract a subtree spanning tips within a certain neighborhood.
}
\description{
Given a rooted tree and a focal tip, extract a subtree comprising various representative nodes and tips in the vicinity of the focal tip using a heuristic algorithm. This may be used for example to display closely related taxa from a reference tree.
}
\usage{
extract_tip_neighborhood(tree, 
                   focal_tip,
                   Nbackward,
                   Nforward,
                   force_tips 		= NULL,
                   include_subtree 	= TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{tree}{
A rooted tree of class "phylo".
}
\item{focal_tip}{
Either a character, specifying the name of the focal tip, or an integer between 1 and Ntips, specifying the focal tip's index.
}
\item{Nbackward}{
Integer >=1, specifying how many splits backward (towards the root) to explore. A larger value of \code{Nbackward} will generally lead to a larger extracted subtree (i.e., including deeper splits).
}
\item{Nforward}{
Non-negative integer, specifying how many splits forward (towards the tips) to explore. A larger value of \code{Nforward} will generally lead to a larger extracted subtree (i.e., including more representative tips from each sister branch).
}
\item{force_tips}{
Optional integer or character list, specifying indices or names of tips to force-include in any case.
}
\item{include_subtree}{
Logical, whether to actually extract the subtree, rather than just returning the inferred neighbor tips.
}
}

\details{
The tree may include multi-furcations as well as mono-furcations (i.e. nodes with only one child). The input tree must be rooted at some node for technical reasons (see function \code{\link{root_at_node}}), but the choice of the root node does not influence which tips are extracted.
}


\value{
A named list with the following elements:
\item{neighbor_tips}{
Integer vector with values in 1,..,Ntips, specifying which tips were found to be neighbors of the focal tip.
}
\item{subtree}{
A new tree of class "phylo", containing a subset of the tips and nodes in the vicinity of the focal tip. Only returned if \code{include_subtree=TRUE}.
}
\item{new2old_tip}{
Integer vector of length Ntips_extracted (=number of tips in the extracted subtree) with values in 1,..,Ntips, mapping tip indices of the extracted subtree to tip indices in the original tree. In particular, \code{tree$tip.label[new2old_tip]} will be equal to \code{subtree$tip.label}. Only returned if \code{include_subtree=TRUE}.
}
}

\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{get_subtree_with_tips}}
}

\examples{
# generate a random tree
Ntips = 50
tree  = generate_random_tree(list(birth_rate_factor=0.1),
                             max_tips = Ntips,
                             tip_basename="tip.")$tree

# extract a subtree in the vicinity of a focal tip
subtree = extract_tip_neighborhood(tree, 
                                   focal_tip="tip.39",
                                   Nbackward=5,
                                   Nforward=2)$subtree
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% Only 1 keyword per line
\keyword{pruning}
\keyword{subtree}
