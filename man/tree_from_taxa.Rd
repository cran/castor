\name{tree_from_taxa}
\alias{tree_from_taxa}
\title{
Construct a rooted tree from lists of taxa.
}
\description{
Given a collection of taxon lists, construct a rooted taxonomic tree. Each taxon list is defined by a parent name and the names of its children (i.e., immediate descendants).
}
\usage{
tree_from_taxa(taxa)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{taxa}{
Named list, whose elements are character vectors, each representing a parent and its children. The element names of \code{taxa} are parents. Each element \code{taxa[n]} is a character vector listing an arbitrary number of taxon names (the children), which immediately descend from taxon \code{names(taxa)[n]}.
}
}


\details{
The following rules apply:
\itemize{
\item Each taxon must appear at most once as a parent and at most once as a child.
\item Any taxa found as parents but not as children, will be assumed to descend directly from the root. If only one such taxon was found, it will become the root itself.
\item Any taxa found as a child but not as a parent, will become tips.
\item Any parents without children will be considered tips.
\item Empty parent names (i.e., "") are not allowed.
\item Taxa can be specified in any arbitrary order, including breadth-first, depth-first etc.
}
Since the returned tree is a taxonomy, it will contain no edge lengths.
}


\value{
A rooted tree of class "phylo".
}

\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{consensus_taxonomies}},
\code{\link{place_tips_taxonomically}}
}


\examples{
# define a list of taxa, with taxon "A" being the root
# Taxa G, H, I, J, K, L, M, N and O will be tips
taxa = list(A = c("B","C","D"), 
            B = c("E","I"), 
            C = c("J","F"), 
            D = c("M", "N", "O"), 
            E = c("G", "H"),
            F = c("K", "L"))
tree = castor::tree_from_taxa(taxa)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% Only 1 keyword per line
\keyword{taxonomy}
