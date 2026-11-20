\name{tree_from_taxonomies}
\alias{tree_from_taxonomies}
\title{
Construct a rooted tree from a collection of taxonomies.
}
\description{
Given a collection of taxonomies, construct a compatible rooted taxonomic tree in which each tip corresponds to an input taxonomy.
Each taxon list is defined by a parent name and the names of its children (i.e., immediate descendants).
}
\usage{
tree_from_taxonomies(taxonomies,
                     tip_names           = NULL,
                     delimiter           = ";",
                     undefined_taxa      = c(),
                     label_nodes         = FALSE,
                     max_level           = NULL,
                     keep_monofurcations = FALSE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{taxonomies}{
Character vector of length Ntips, specifying a taxonomic path for each tip.
}
\item{tip_names}{
Optional character vector of length Ntips, specifying tip names to use regardless of their taxonomies.
For example, these may bey OTU IDs. If \code{NULL}, tip names will be ``tip.1'', ``tip.2'' and so on.
}
\item{delimiter}{
Character, the delimiter used between taxonomic levels.
}
\item{undefined_taxa}{
Optional character vector specifying taxon names or wildcards (case-insensitive) to be considered as practically undefined/unknown, for example \code{c("Incertae sedis", "*unknown*")}.
The empty string ("") is always considered undefined, so there's no need to include it here.
}
\item{label_nodes}{
Logical. If \code{TRUE}, then partial taxonomies will be used as node labels. If \code{FALSE}, no node labels are stored.
}
\item{max_level}{
Optional non-negative integer, maximum number of taxonomic levels to consider. 
Any splits past this threshold will be collapsed to multifurcations. 
This threshold applies to the raw taxonomies, i.e., including any undefined/unknown levels. 
For example, if max_level=3, then the taxonomies "A;B;C;D;E", "A;B;C;X;Y", "A;B;C;S;T" and "A;B;C;;W" will all branch off from the same "A;B;C" node.
}
\item{keep_monofurcations}{
Logical, whether to keep monofurcations in the final tree. If \code{FALSE}, monofurcating nodes are collapsed.
}
}


\details{
The following rules apply:
\itemize{
\item Each input taxonomy will correspond to a distinct tip in the tree. Tips are indexed in the same order as the input taxonomies.
\item Any tips sharing the same taxonomic path up until a given level, will descend from a node corresponding to that shared path.
\item If all taxonomies share the same top-level taxon name, this taxon name will correspond to the tree's root.
\item If multiple top-level taxon names are encountered, each of them will correspond to a separate node, all of which will descend from a common root.
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
\code{\link{tree_from_taxa}},
\code{\link{place_tips_taxonomically}}
}


\examples{
taxonomies = c("Bacillales;Bacillaceae;Priestia;Bacillus megaterium",
               "Paenibacillales;Paenibacillaceae;Paenibacillus;P. taiwanensis",
               "Bacillales;Bacillaceae;Bacillus;Bacillus amyloliquefaciens",
               "Bacillales;Bacillaceae;Bacillus;",
               "Lactobacillales;Enterococcaceae;Enterococcus;E. mundtii",
               "Bacillales;Bacillaceae;Ureibacillus;U. thermosphaericus")
tree = castor::tree_from_taxonomies(taxonomies, label_nodes = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% Only 1 keyword per line
\keyword{taxonomy}
