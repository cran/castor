# Given a rooted tree and taxonomies for all tips, figure out consensus taxonomies for each node in the tree.
# The consensus taxonomy of a given node is the longest possible taxonomic path (i.e., to the lowest possible level) such that the taxonomies of all descending tips are nested within that taxonomic path.
# Some tip taxonomies may be incomplete, i.e., truncated at higher taxonomic levels. In that case, consensus taxonomy building will be conservative, i.e., no assumptions will be made about the missing taxonomic levels.
# Examples: 
#	If the descending tips of a node have taxonomies "A;B;C" & "A;B;C;D" & "A;B;C;E", then their consensus taxonomy is "A;B;C".
#	If the descending tips of a node have taxonomies "A;B" & "A;B;C;D" & "A;B;C;E", then their consensus taxonomy is "A;B".
consensus_taxonomies = function(tree,					# rooted tree of class phylo, 
								tip_taxonomies 	= NULL, # optional character vector of length Ntips, listing taxonomic paths for the tips. If NULL, then tip labels are assumed to be tip taxonomies.
								delimiter 		= ";",	# character, the delimiter between taxonomic levels (e.g., ";" for SILVA taxonomies)
								duplicates		= "keep_all"){	# how to handle duplicate consensus taxonomies, i.e., when multiple nodes get the same consensus taxonomies. Options are "keep_all", "keep_oldest", "keep_none".
	if(is.null(tip_taxonomies)) tip_taxonomies = tree$tip.label
	results = consensus_taxonomies_CPP(	Ntips 			= length(tree$tip.label),
										Nnodes			= tree$Nnode,
										Nedges			= nrow(tree$edge),
										tree_edge		= as.vector(t(tree$edge)) - 1,
										tip_taxonomies	= tip_taxonomies,
										delimiter		= delimiter)
	if(any(duplicated(results$node_taxonomies))){
		if(tolower(duplicates)=="keep_all"){
			# nothing to do
		}else if(tolower(duplicates)=="keep_none"){
			# eliminate (set to blank) any consensus taxonomies that are duplicated
			results$node_taxonomies[duplicated(results$node_taxonomies) | duplicated(results$node_taxonomies, fromLast=TRUE)] = ""
		}else if(tolower(duplicates)=="keep_oldest"){
			node_order = order(get_all_distances_to_root(tree))
			results$node_taxonomies[node_order[duplicated(results$node_taxonomies[node_order])]] = ""
		}else{
			stop(sprintf("Invalid choice '%s' for duplicates argument.",duplicates))
		}
	}
	return(results$node_taxonomies)
}