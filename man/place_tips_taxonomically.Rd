\name{place_tips_taxonomically}
\alias{place_tips_taxonomically}
\title{
Place queries on a tree based on taxonomic identities.
}
\description{
Given a rooted tree with associated tip & node taxonomies, as well as a list of query taxonomies, place the queries on nodes of the tree based on taxonomic identity. Each query is placed at the deepest possible node (furthest from the root in terms of splits) for which it is certain that the query is a descendant of.
}
\usage{
place_tips_taxonomically( tree,
                          query_labels,
                          query_taxonomies      = NULL,
                          tip_taxonomies        = NULL,
                          node_taxonomies       = NULL,
                          tree_taxon_delimiter  = ";",
                          query_taxon_delimiter = ";",
                          include_expanded_tree = TRUE)
}
\arguments{
\item{tree}{
Rooted tree of class "phylo".
}
\item{query_labels}{
Character vector of length Nqueries, listing labels for the newly placed tips.
}
\item{query_taxonomies}{
Optional character vector of length Nqueries, listing the taxonomic paths of the queries. If \code{NULL}, it is assumed that \code{query_labels} are taxonomies.
}
\item{tip_taxonomies}{
Optional character vector of length Ntips, listing taxonomic paths for the tree's tips. If \code{NULL}, then tip labels are assumed to be tip taxonomies.
}
\item{node_taxonomies}{
Optional character vector of length Nnodes, listing taxonomic paths for the tree's nodes. If \code{NULL}, then node labels are assumed to be node taxonomies.
}
\item{tree_taxon_delimiter}{
Character, the delimiter between taxonomic levels in the tree's tip & node taxonomies (e.g., ";" for SILVA taxonomies).
}
\item{query_taxon_delimiter}{
Character, the delimiter between taxonomic levels in \code{query_taxonomies}.
}
\item{include_expanded_tree}{
If \code{TRUE}, the expanded tree (i.e., including the placements) is returned as well, at some computational cost. If \code{FALSE}, only the placement info is returned, but no tree expansion is done. 
}
}
\details{
This function assumes that the tip & node taxonomies are somewhat consistent with each other and with the tree's topology.
}
\value{
A named list with the following elements:
\item{placement_nodes}{
Integer vector of length Nqueries, with values in 1,..,Nnodes, specifying for each query the node on which it was placed. For queries that could not be placed on the tree, the value 0 is used.
}
If \code{include_expanded_tree} was \code{TRUE}, the following additional elements are included:
\item{tree}{
Object of class "phylo", the extended tree constructed by adding the placements on the original tree.
}
\item{placed_tips}{
Integer vector of length Nqueries, specifying which tips in the returned tree correspond to placements. For queries that could not be placed on the tree, the value 0 is used.
}
}


\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{expanded_tree_from_jplace}}
}

%\examples{
%\dontrun{
%}
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{placement}
\keyword{taxonomy}
