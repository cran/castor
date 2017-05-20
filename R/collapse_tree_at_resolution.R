# Collapse tree nodes (and their descending subtrees) into tips, whenever all descending tips have a distance from a node below a certain phylogenetic resolution threshold
# Any node whose distance to all its descending tips is <=distance_threshold, will be collapsed into a single tip
# This function can be used to get the "coarse structure" of a tree
collapse_tree_at_resolution = function(tree, resolution=0, by_edge_count=FALSE){
	Ntips  = length(tree$tip.label);
	Nnodes = tree$Nnode;
	Nedges = nrow(tree$edge);
	
	# collapse
	results = collapse_tree_at_resolution_CPP(	Ntips		= Ntips,
												Nnodes		= Nnodes,
												Nedges 		= Nedges,
												tree_edge 	= as.vector(t(tree$edge)) - 1,
												edge_length	= (if(by_edge_count || is.null(tree$edge.length)) numeric() else tree$edge.length),
												resolution	= resolution);
	
	# reformat results into a valid "phylo" object
	# note that some of the old nodes may have turned into new tips
	Ntips_new  		= results$Ntips_new
	Nnodes_new	 	= results$Nnodes_new
	Nclades_new		= Ntips_new+Nnodes_new
	new2old_clade 	= results$new2old_clade + 1; # switch to 1-based indices
	new2old_edge	= results$new2old_edge + 1;
	clade_labels	= c(tree$tip.label, tree$node.label)
	collapsed_tree = list(	Nnode 		= Nnodes_new,
							tip.label 	= clade_labels[new2old_clade[1:Ntips_new]],
							node.label 	= (if(is.null(tree$node.label)) NULL else clade_labels[new2old_clade[(Ntips_new+1):Nclades_new]]),
							edge 		= matrix(results$new_tree_edge,ncol=2,byrow=TRUE) + 1,
							edge.length = (if(is.null(tree$edge.length)) NULL else tree$edge.length[new2old_edge]),
							root 		= results$new_root+1)
	class(collapsed_tree) = "phylo";
	return(list(collapsed_tree	= collapsed_tree, 
				new2old_clade	= new2old_clade, 
				new2old_edge	= new2old_edge));
}
