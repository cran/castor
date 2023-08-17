# Given a rooted tree and a focal tip, extract a subtree comprising all tips within a certain patristic-distance radius of the focal tip.
# If the tree lacks edge lengths, then each edge is assumed to have length 1.
extract_tip_radius = function(	tree,		# rooted phylogenetic tree of type "phylo"
								focal_tip, 	# integer or character, specifying the index or name of a focal tip
								radius,		# positive numeric, the patristic distance radius
								include_subtree = TRUE){ # whether to actually extract the subtree, rather than just returning the list of within-radius tips
	Ntips 		= length(tree$tip.label)
	distances 	= get_all_distances_to_tip(tree=tree, focal_tip=focal_tip)
	radius_tips = which(distances[1:Ntips]<=radius)
									
	# extract subtree comprisign radius_tips
	if(include_subtree){
		subtreeing = get_subtree_with_tips(tree, only_tips=radius_tips, collapse_monofurcations=TRUE, force_keep_root=FALSE)
		return(list(radius_tips	= radius_tips,
					subtree 	= subtreeing$subtree,
					new2old_tip = subtreeing$new2old_tip))
	}else{
		return(list(radius_tips = radius_tips))
	}
}