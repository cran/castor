\name{extract_deep_frame}
\alias{extract_deep_frame}
\title{
Extract tips representing a tree's deep splits.
}
\description{
Given a rooted phylogenetic tree, extract a subset of tips representing the tree's deepest splits (nodes), thus obtaining a rough "frame" of the tree. For example, if \code{Nsplits=1} and the tree is bifurcating, then two tips will be extracted representing the two clades splitting at the root.
}
\usage{
extract_deep_frame(tree,
                   Nsplits    = 1,
                   only_tips  = FALSE)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo".
}
\item{Nsplits}{
Strictly positive integer, specifying the maximum number of splits to descend from the root. A larger value generally implies that more tips will be extracted, representing a larger number of splits.
}
\item{only_tips}{
Boolean, specifying whether to only return the subset of extracted tips, rather than the subtree spanned by those tips. If \code{FALSE}, a subtree is returned in addition to the tips (this comes at a computational cost).
}
}
\details{
The tree may include multi-furcations (i.e. nodes with more than 2 children) as well as mono-furcations (i.e. nodes with only one child). No guarantee is made as to the precise subset of tips extracted.
}
\value{
A named list with the following elements:
\item{tips}{
Integer vector with values from 1 to Ntips, listing the indices of the extracted tips.
}
\item{subtree}{
A new tree of class "phylo", the subtree spanned by the extracted tips. Only included if \code{only_tips==FALSE}.
}
}


\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{get_pairwise_mrcas}}, \code{\link{get_tips_for_mrcas}}
}

\examples{
# generate a random tree
Ntips = 100
tree = generate_random_tree(list(birth_rate_factor=0.1),Ntips)$tree

# extract a subtree representing the deep splits
subtree = extract_deep_frame(tree, Nsplits=3)$subtree

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{MRCA}
