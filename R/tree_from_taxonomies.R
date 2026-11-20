# Given a collection of taxonomic paths (for example, taxonomies assigned to OTUs), construct a taxonomic tree based on those taxonomies.
# Rules:
#	Each taxonomy will correspond to a distinct tip in the tree.
#	Any tips sharing the same taxonomic path up until a given level L, will descend from a node corresponding to that shared path.
#	If all taxonomies share the same top-level taxon name, this taxon name will correspond to the tree's root.
#	If multiple top-level taxon names are encountered, each of them will correspond to a separate node, all of which will descend from a common root.
# Since the returned tree is a taxonomy, it will contain no edge lengths.
# Tip indices will match the order in which the taxonomies were provided
tree_from_taxonomies = function(taxonomies,						# character vector of length Ntips, specifying a taxonomic path for each tip
								tip_names 		 	= NULL, 	# optional character vector of length Ntips, specifying tip labels to use regardless of their taxonomies. For example, these may bey OTU IDs. If NULL, tip names will be tip.1, tip.2 etc.
								delimiter		 	= ";",		# character at which to split taxonomies into multiple levels.
								undefined_taxa	 	= c(),		# optional character vector specifying taxon names (or shell wildcards, case-insensitive) to be considered practically undefined/unknown, for example c("Incertae sedis", "unknown"). The empty string ("") is always considered unknown, so there's no need to include it here.
								label_nodes		 	= FALSE,	# if TRUE, then partial taxonomies will be used as node labels. If FALSE, no node labels are stored.
								max_level			= NULL,		# Optional integer, maximum number of taxonomic levels to consider. Any splits past this threshold will be collapsed to multifurcations. This threshold applies to the raw taxonomies, i.e., including any undefined/unknown levels. For example, if max_level=3, then the taxonomies "A;B;C;D;E", "A;B;C;X;Y", "A;B;C;S;T" and "A;B;C;;W" will all branch off from the same "A;B;C" node.
								keep_monofurcations = FALSE){
	# basic input checks
	Ntips = length(taxonomies)
	if(Ntips<2) stop("Need at least 2 taxonomies to make a tree")
	if(!is.character(taxonomies)) stop("taxonomies must be a character vector")
	if((!is.null(tip_names)) && (length(tip_names)!=Ntips)) stop(sprintf("tip_names must have length equal to Ntips (%d)",Ntips))
	
	# split taxonomies into their various levels
	taxonomies = strsplit(taxonomies, split=delimiter, fixed=TRUE)
	
	# truncate taxonomies to max_level
	if(!is.null(max_level)) taxonomies = lapply(taxonomies, FUN=function(taxonomy) taxonomy[1:min(length(taxonomy),max_level)])

	# omit empty levels and undefined levels
	taxonomies = lapply(taxonomies, FUN=function(taxonomy) taxonomy[nzchar(taxonomy)])
	for(undefined in undefined_taxa){
		undefined_regex = utils::glob2rx(undefined)
		taxonomies = lapply(taxonomies, FUN=function(taxonomy) taxonomy[!grepl(undefined_regex, taxonomy, ignore.case=TRUE)])
	}

	# determine the root name
	root_name = unique(sapply(taxonomies, FUN=function(taxonomy) (if(length(taxonomy)==0) "" else taxonomy[1])))
	if((length(root_name)!=1) || (root_name=="")){
		# either there are multiple level-1 taxa, or all taxonomies are essentially empty
		# Hence, define another taxonomic level on top of all taxonomies, representing the root
		root_name = "root"
		taxonomies = lapply(taxonomies, FUN=function(taxonomy) c(root_name,taxonomy))
	}

	# compute series of subtaxonomies leading to each tip
	# each subtaxonomies[tip] will be a character vector of length <= max_length
	subtaxonomies = lapply(taxonomies, FUN=function(taxonomy){ sapply(seq_along(taxonomy), FUN=function(level) paste(taxonomy[1:level], collapse=delimiter))})
	
	# define nodes corresponding to unique subtaxonomies
	# It is guaranteed that the first element in unique_subtaxonomies is the root name, since unique() preserves order of occurrence
	unique_subtaxonomies = unique(unlist(subtaxonomies))
	Nnodes 				 = length(unique_subtaxonomies)
	subtaxonomy2node 	 = setNames(seq_len(Nnodes), unique_subtaxonomies)
	
	# define all lineages, with each lineage[[tip]] being a sequence of integers (from 1 to Nclades) traversing the tree
	lineages = lapply(seq_len(Ntips), FUN=function(tip) unname(c(Ntips+subtaxonomy2node[subtaxonomies[[tip]]], tip)))
	
	# define edges based on the unique pairs of successive clades encountered in the lineages
	edges_to_each_tip = lapply(lineages, function(lineage) { cbind(head(lineage, -1L), tail(lineage, -1L)) })
	edges = unique(do.call(rbind, edges_to_each_tip), MARGIN=1)

	
	# construct the tree object
	tree = list(Nnode 		= Nnodes,
				tip.label 	= (if(is.null(tip_names)) sprintf("tip.%d",seq_len(Ntips)) else tip_names),
				edge 		= edges,
				root 		= Ntips+1)
	if(label_nodes) tree$node.label = unique_subtaxonomies
	class(tree) = "phylo"
	attr(tree,"order") = NULL
	
	# eliminate monofurcating nodes if needed
	if(!keep_monofurcations) tree = collapse_monofurcations(tree, force_keep_root=FALSE)$tree
	tree$root.edge = 0 # add a dummy root.edge, so that ape knows the tree is rooted. Otherwise, if the root is multifurcating, ape will think the tree is not rooted.

	return(tree)
}
