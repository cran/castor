# Extract a subset of tips representing the tree's deepest splits, thus obtaining a rough "frame" of the tree.
# For example, if Nsplits=1 and the tree is bifurcating, then two tips will be extracted representing the two clades splitting at the root.
extract_deep_frame = function(tree, Nsplits=1, only_tips=FALSE){
	Nsplits = max(1,Nsplits)
	frame_tips = extract_deep_frame_CPP(Ntips 		= length(tree$tip.label),
										Nnodes		= tree$Nnode,
										Nedges		= nrow(tree$edge),
										tree_edge 	= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
										Nsplits		= Nsplits) + 1L
	if(!only_tips){
		subtree = get_subtree_with_tips(tree, only_tips=frame_tips, collapse_monofurcations=TRUE, force_keep_root=TRUE)$subtree
	}
	return(list(tips 	= frame_tips,
				subtree = (if(only_tips) NULL else subtree)))
}