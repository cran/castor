# Given a rooted tree and a focal tip, extract a subtree comprising a few representative nearby tips using a heuristic algorithm
# This may be used for example for displaying sister clades of some focal novel strain, pulled from a database
extract_tip_neighborhood = function(tree,		# rooted phylogenetic tree of type "phylo"
									focal_tip, 	# integer or character, specifying the index or name of a focal tip
									Nbackward,	# integer >=1, specifying how many splits backward (towards the root) to explore
									Nforward,	# integer >=0, specifying how many splits forward (towards the tips) to explore
									force_tips 		= NULL,	# optional integer or character list, specifying tips to force-include in any case
									include_subtree = TRUE){# whether to actually extract the subtree, rather than just returning the list of neighbor tips
	Nbackward = max(1,Nbackward)
	Nforward  = max(0,Nforward)
	focal_tip = map_tip_or_node_names_to_indices(tree, focal_tip, type="tip", list_title="tip", check_input=TRUE)

	# determine neighbor tips
	neighbor_tips = extract_tip_neighborhood_CPP(Ntips		= length(tree$tip.label),
												Nnodes		= tree$Nnode,
												Nedges		= nrow(tree$edge),
												tree_edge	= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
												focal_tip	= focal_tip-1,
												Nbackward	= Nbackward,
												Nforward	= Nforward) + 1L

	if((!is.null(force_tips)) && (length(force_tips)>0)){
		force_tips = map_tip_or_node_names_to_indices(tree, force_tips, type="tip", list_title="force_tips", check_input=TRUE)
		neighbor_tips = unique(c(neighbor_tips,force_tips))
	}

	if(include_subtree){
		# extract subtree, spanning the focal tip & its neighbors
		subtreeing = get_subtree_with_tips(tree, only_tips=neighbor_tips, collapse_monofurcations=TRUE, force_keep_root=FALSE)
		return(list(neighbor_tips	= neighbor_tips,
					subtree 		= subtreeing$subtree,
					new2old_tip 	= subtreeing$new2old_tip))
	}else{
		return(list(neighbor_tips = neighbor_tips))
	}
}