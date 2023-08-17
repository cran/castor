\name{consensus_taxonomies}
\alias{consensus_taxonomies}
\title{
Compute consensus taxonomies across a tree.
}
\description{
Given a rooted phylogenetic tree and taxonomies for all tips, compute corresponding consensus taxonomies for all nodes in the tree based on their descending tips. The consensus taxonomy of a given node is the longest possible taxonomic path (i.e., to the lowest possible level) such that the taxonomies of all descending tips are nested within that taxonomic path. Some tip taxonomies may be incomplete, i.e., truncated at higher taxonomic levels. In that case, consensus taxonomy building will be conservative, i.e., no assumptions will be made about the missing taxonomic levels.
}
\usage{
consensus_taxonomies(tree, tip_taxonomies = NULL, delimiter = ";")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{tree}{
A rooted tree of class "phylo".
}
\item{tip_taxonomies}{
Optional, character vector of length Ntips, listing taxonomic paths for the tips. If \code{NULL}, then tip labels in the tree are assumed to be tip taxonomies.
}
\item{delimiter}{
Character, the delimiter between taxonomic levels (e.g., ";" for SILVA and Greengenes taxonomies).
}
}


\details{
Examples:
\itemize{
\item If the descending tips of a node have taxonomies "A;B;C" and "A;B;C;D" and "A;B;C;E", then their consensus taxonomy is "A;B;C".
\item If the descending tips of a node have taxonomies "A;B" and "A;B;C;D" and "A;B;C;E", then their consensus taxonomy is "A;B".
\item If the descending tips of a node have taxonomies "X;B;C" and "Y;B;C", then their consensus taxonomy is "" (i.e., empty).
}
}


\value{
A character vector of length Nnodes, listing the inferred consensus taxonomies for all nodes in the tree.
}


\author{Stilianos Louca}

%\references{
%}

%\seealso{
%\code{\link{find_nearest_tips}}
%}

\examples{
# generate a random tree
tree = generate_random_tree(list(birth_rate_factor=0.1), max_tips=7)$tree

# define a character vector storing hypothetical tip taxonomies
tip_taxonomies = c("R;BB;C", "AA;BB;C;DD", "AA;BB;C;E", 
                   "AA;BB", "AA;BB;D", "AA;BB;C;F", "AA;BB;C;X")

# compute consensus taxonomies and store them in the tree as node labels
tree$node.label = castor::consensus_taxonomies(tree, 
                                               tip_taxonomies = tip_taxonomies, 
                                               delimiter = ";")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ Pruning }
