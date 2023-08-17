# Given a rooted tree and taxonomies for all tips, figure out consensus taxonomies for each node in the tree
# The consensus taxonomy of a given node is the longest possible taxonomic path (i.e., to the lowest possible level) such that the taxonomies of all descending tips are nested within that taxonomic path.
# Some tip taxonomies may be incomplete, i.e., truncated at higher taxonomic levels. In that case, consensus taxonomy building will be conservative, i.e., no assumptions will be made about the missing taxonomic levels.
# Examples: 
#	If the descending tips of a node have taxonomies "A;B;C" and "A;B;C;D" and "A;B;C;E", then their consensus taxonomy is "A;B;C".
#	If the descending tips of a node have taxonomies "A;B" and "A;B;C;D" and "A;B;C;E", then their consensus taxonomy is "A;B".
consensus_taxonomies = function(tree,					# rooted tree of class phylo, 
								tip_taxonomies 	= NULL, # optional character vector of length Ntips, listing taxonomic paths for the tips. If NULL, then tip labels are assumed to be tip taxonomies.
								delimiter 		= ";"){	# character, the delimiter between taxonomic levels (e.g., ";" for SILVA taxonomies)
	if(is.null(tip_taxonomies)) tip_taxonomies = tree$tip.label
	results = consensus_taxonomies_CPP(	Ntips 			= length(tree$tip.label),
										Nnodes			= tree$Nnode,
										Nedges			= nrow(tree$edge),
										tree_edge		= as.vector(t(tree$edge)) - 1,
										tip_taxonomies	= tip_taxonomies,
										delimiter		= delimiter)
	return(results$node_taxonomies)
}