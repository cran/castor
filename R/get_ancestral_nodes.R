# Given a rooted tree and one or more "descendant" tips or nodes, determine their ancestral node indices, traveling N splits back in time (i.e., N splits towards the root).
get_ancestral_nodes = function(	tree,			# rooted phylogenetic tree of type "phylo"
								descendants, 	# integer or character vector of length ND, specifying tip/node indices or names
								Nsplits){		# either a single integer or an integer vector of length ND, with values >=1, specifying how many splits backward to travel. For example, Nsplits=1 will yield the immediate ancestral nodes.
	Ntips  = length(tree$tip.label);
	Nnodes = tree$Nnode;
	descendants = map_tip_or_node_names_to_indices(tree, descendants, type="both", list_title="descendants", check_input=TRUE)
	ancestors = get_ancestral_nodes_CPP(Ntips			= Ntips,
										Nnodes			= Nnodes,
										Nedges			= nrow(tree$edge),
										tree_edge		= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
										descendants		= descendants-1,
										Nsplits			= Nsplits)
	return(ancestors+1L)
}