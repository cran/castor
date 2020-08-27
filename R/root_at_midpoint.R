# adjust edge directions such that the new root is at the midpoint (node with smallest max_distance_to_a_tip)
# Note that the number of tips & nodes remains the same
# If update_indices==FALSE, then tip & node indices also remain the same
root_at_midpoint = function(tree, 
							update_indices	= TRUE,
							as_edge_counts 	= FALSE,	# calculate distances in terms of cumulative edge counts (as as if each edge had length 1)
							is_rooted		= FALSE){ 	# if TRUE, the caller guarantees that the input tree is rooted
	Ntips 	= length(tree$tip.label);
	Nnodes	= tree$Nnode;
	
	# root arbitrarily if needed
	if(!is_rooted){
		tree = root_at_node(tree, 1, update_indices=FALSE)
	}
	
	# figure out midpoint root node
	results = get_farthest_tip_per_clade_CPP(	Ntips					= Ntips,
												Nnodes					= Nnodes,
												Nedges					= nrow(tree$edge),
												tree_edge				= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
												edge_length				= (if(as_edge_counts || is.null(tree$edge.length)) numeric() else tree$edge.length),
												onlyToTips				= integer(),
												only_descending_tips	= FALSE,
												verbose					= FALSE,
												verbose_prefix			= "");
	new_root_node = which.min(results$farthest_distances[(Ntips+1):(Ntips+Nnodes)])
	
	# place root at midpoint node
	new_edges = root_tree_at_node_CPP(	Ntips			= Ntips,
										Nnodes			= Nnodes,
										Nedges			= nrow(tree$edge),
										tree_edge		= as.vector(t(tree$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
										new_root_node	= new_root_node-1);
	new_edges = matrix(new_edges+1, ncol=2, byrow=TRUE); # unflatten returned table and shift clade indices to 1-based
	tree$edge = new_edges;
	tree$root = Ntips + new_root_node;
	
	# update node indices if required
	correct_root_node = 1; # correct index that the root node should have
	if(update_indices && (new_root_node!=correct_root_node)){
		# swap indices with wrong node
		temp_root_index = 0;
		tree$edge[tree$edge==(Ntips+new_root_node)] 	= temp_root_index;
		tree$edge[tree$edge==(Ntips+correct_root_node)] = Ntips+new_root_node;
		tree$edge[tree$edge==temp_root_index] 			= Ntips+correct_root_node;
		tree$root 										= Ntips+correct_root_node;
		
		if(!is.null(tree$node.label)){
			root_label 							= tree$node.label[new_root_node];
			tree$node.label[new_root_node] 		= tree$node.label[correct_root_node];
			tree$node.label[correct_root_node] 	= root_label;
		}
	}
	
	return(tree);
}