\name{inherit_taxonomies}
\alias{inherit_taxonomies}
\title{
Inherit ancestral taxonomies from a phylogeny.
}
\description{
Given a rooted tree and taxon names associated with some nodes and/or tips, infer the inherited taxonomic paths for all tips and/or nodes by concatenating the ancestral taxon names.
Tip/node taxon names can be specified as separate vectors or taken from the tip/node labels stored in the tree itself.
}
\usage{
inherit_taxonomies(tree,
                   tip_taxa      = NULL,
                   node_taxa     = NULL,
                   delimiter     = ";",
                   only_proper   = FALSE,
                   halt          = "last",
                   strict_proper = TRUE,
                   focal_tips    = NULL,
                   focal_nodes   = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{tree}{
A rooted tree of class "phylo".
}
\item{tip_taxa}{
Optional character vector of length Ntips, listing taxon names for the tips. These will be appended to the inherited taxonomies of the tips. If \code{NULL}, then tip labels are used as taxon names. If you don't want any taxon names associated with tips, pass an empty vector. Some elements in \code{tip_taxa} may be ""/\code{NULL}/\code{NA}.
}
\item{node_taxa}{
Optional character vector of length Nnodes, listing taxon names for the nodes. If NULL, then node labels are assumed to be taxon names. Some elements in \code{node_taxa} may be ""/\code{NULL}/\code{NA}.
}
\item{delimiter}{
Character, the delimiter to use between taxonomic levels (e.g., ";" for SILVA taxonomies)
}
\item{only_proper}{
Logical, whether to only include proper taxa in the inferred taxonomies. See option \code{strict_proper} on how proper taxa are defined.
}
\item{halt}{
Character, specifying when to halt accumulating taxon levels, moving from root to tips. Options are \code{"last"} (halt after the last non-empty taxon, i.e., don't skip any levels), \code{"last_proper"} (halt at the last properly name taxon, including the latter), \code{"first_improper"} (halt at the first improperly named non-empty taxon, not including the latter) and \code{"first_empty"} (halt at the first empty taxon).
}
\item{strict_proper}{
Logical, specifying whether to define taxa as proper in the strict sense, i.e., only if the taxon name is a properly encoded taxon name (i.e., contains no numeric digits, dashes or "candidatus"). If \code{FALSE}, then a taxon is considered proper as long as any of the descendants are proper in the strict sense; for example, the taxon "WOR-3" would be considered proper if one of its descending clades has taxon name "Hydrothermia".
}
\item{focal_tips}{
Optional specification of tips to focus on. Either a character vector listing tip names, or an integer vector listing tip indices. If \code{NULL}, inherited taxonomies will be returned for all tips.
}
\item{focal_nodes}{
Optional specification of nodes to focus on. Either a character vector listing node names, or an integer vector listing node indices. If \code{NULL}, inherited taxonomies will be returned for all nodes.
}
}


\details{
This function may be used for example to infer taxonomic paths of new genomes placed on a phylogeny/taxonomy such as the GTDB.
Some taxa may be empty (""), i.e., not all tips/nodes are expected to have taxon names; empty taxa are never included in the returned taxonomies, not even as placeholders.
}


\value{
A named list with the following elements:
\itemize{
\item \code{tip_taxonomies}: Character vector of length Ntips (if \code{focal_tips=NULL}) or of the same length as \code{focal_tips}, listing the inferred "inheritted" taxonomies at the tips.
\item \code{node_taxonomies}: Character vector of length Nnodes (if \code{focal_nodes=NULL}) or of the same length as \code{focal_nodes}, listing the inferred "inheritted" taxonomies at each node.
\item \code{tip_levels}: Integer vector of the same length as \code{tip_taxonomies}, listing the number of concatenated taxonomic levels for each inferred tip taxonomy. Note that each taxon name is counted once, even if it consists of multiple levels; for example, if a node has taxon name "WOR-3; Hydrothermia" then any descending taxonomies will count this taxon as one level.
\item \code{node_levels}: Integer vector of the same length as \code{node_taxonomies}, listing the number of concatenated taxonomic levels for each inferred node taxonomy, similarly to \code{tip_taxonomies}.
}
}


\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{place_tips_taxonomically}}, \code{\link{consensus_taxonomies}}
}

\examples{
# generate a random tree
tree = generate_random_tree(list(birth_rate_factor=0.1), max_tips=7)$tree

# define a character vector storing hypothetical node taxa
node_taxa = paste("taxon.",seq_len(tree$Nnode),sep="")

# infer inheritted taxonomies for all tips
inheritted = castor::inherit_taxonomies(tree, node_taxa = node_taxa)
print(inheritted$tip_taxonomies)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{taxonomy}
