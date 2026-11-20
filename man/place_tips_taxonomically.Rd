\name{place_tips_taxonomically}
\alias{place_tips_taxonomically}
\title{
Place queries on a tree based on taxonomic identities.
}
\description{
Given a rooted tree with associated tip & node taxonomies, as well as a list of query taxonomies, place the queries on nodes of the tree based on taxonomic identity. 
Loosely speaking, each query is placed at the deepest possible node (furthest from the root in terms of splits) for which it is certain that the query is a descendant of (see below for details).
}
\usage{
place_tips_taxonomically( tree,
                          query_labels,
                          query_taxonomies      = NULL,
                          tip_taxonomies        = NULL,
                          node_taxonomies       = NULL,
                          tree_taxon_delimiter  = ";",
                          query_taxon_delimiter = ";",
                          include_expanded_tree = TRUE,
                          placement_length      = "mean_depth",
                          strict                = TRUE,
                          omit_placeholders     = NULL)
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
Optional character vector of length Ntips, listing taxonomic paths for the tree's tips.
If \code{NULL}, then tip labels are assumed to be tip taxonomies.
Only used if \code{strict=TRUE}.
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
\item{placement_length}{
Character or numeric, specifying how to determine the edge length leading into each placed query. Options are "mean_depth" (mean distance of the placement node to any of its descending tips), "max_depth" (maximum distance of the placement node to any of its descending tips), "min_depth" (minimum distance of the placement node to any of its descending tips) or any non-negative numeric value (fixed edge length). Only relevant if the tree has edge lengths and if \code{include_expanded_tree==TRUE}.
}
\item{strict}{
Logical, whether to strictly follow taxonomic paths, i.e., without skipping any mismatched taxon names.
If there is a possibility that the same taxon name might refer to different clades depending on its ancestors, leave this at \code{TRUE}. 
If each taxon name defines a unique clade, this can be set to \code{FALSE}, which is more forgiving to cases where the tree's taxonomies and query taxonomies are not fully consistent.
}
\item{omit_placeholders}{
Optional character, case-insensitive regular expression for \code{grepl} (POSIX extended) for detecting and omitting placeholder taxa (such as "\emph{Incertae Sedis}") from the tip, node and query taxonomies.
This is most critical when \code{strict=FALSE}. If \code{NULL} or "", no such removal is done. 
For example, you might set this to \code{"^(incertae sedis|unknown( [a-z]+)?|family [ivxl]+|uncultured.*|unidentified.*|unclassified.*|unassigned|metagenome.*|environmental.*)$"}.
}
}
\details{
If \code{strict=TRUE}, each query is placed at the most recent common ancestor (MRCA) of all possible ancestral nodes with maximal height (PAMH nodes) and maximum taxonomic precision.
A node is considered a PAMH node if and only if the query could in principle descend from that node (based on the provided node and query taxonomies) but surely cannot descend from any of the node's child nodes.
The taxonomic precision of a node is the number of meaningful available taxonomic levels; only PAMH nodes with the highest taxonomic precision are considered for the MRCA.
This function only makes sense if the tip & node taxonomies are somewhat consistent with each other and with the tree's topology.

If \code{strict=FALSE}, each taxon name is assumed to define a unique clade, regardless of the taxonomic path leading to it. 
Each node's taxon is the lowest level of its taxonomy, ignoring placeholder taxa matching \code{omit_placeholders}.
Each query is placed at the MRCA of all nodes whose taxon matches the lowest possible level of the query's taxonomy.
Queries without any match are placed at the root if the root's taxonomy is empty, and are otherwise not placed.
This mode tolerates inconsistencies between query and tree taxonomies, such as extra or missing levels (e.g., kingdoms) or renamed higher taxa.
Only node taxonomies are used (\code{tip_taxonomies} is ignored), so queries are never placed onto tips.
Since recurring placeholder names (e.g., \dQuote{\emph{Incertae Sedis}} or \dQuote{\emph{Family XI}}) may otherwise cause misplacements, these should be omitted via \code{omit_placeholders}.

Taxon matching is case-insensitive throughout, and any whitespace flanking taxonomic levels is ignored.
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
Integer vector of length Nqueries, specifying the newly placed tips in the expanded tree, in the same order as the original provided queries. For queries that could not be placed on the tree, the value 0 is used.
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
