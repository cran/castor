# Given a rooted tree and a focal tip, compute the patristic distances of all clades to the focal tip
# If the tree lacks edge lengths, then each edge is assumed to have length 1.
get_all_distances_to_tip = function(tree,		# rooted phylogenetic tree of type "phylo"
									focal_tip){	# integer or character, specifying the index or name of a focal tip
	focal_tip = map_tip_or_node_names_to_indices(tree, focal_tip, type="tip", list_title="tip", check_input=TRUE)
	distances = get_all_distances_to_tip_CPP(	Ntips		= length(tree$tip.label),
												Nnodes		= tree$Nnode,
												Nedges		= nrow(tree$edge),
												tree_edge	= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
												edge_length	= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
												focal_tip	= focal_tip-1)		
	return(distances)
}