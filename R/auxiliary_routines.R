get_adjacent_edges_per_edge = function(tree){
	Nedges = nrow(tree$edge)
	adjacents = get_adjacent_edges_per_edge_CPP(	Ntips 		= length(tree$tip.label),
													Nnodes		= tree$Nnode,
													Nedges		= Nedges,
													tree_edge	= as.vector(t(tree$edge)) - 1);
	
	# update indices from 0-based to 1-based
	return(lapply(1:Nedges,FUN = function(edge) adjacents[[edge]]+1))
}



get_outgoing_edges_per_clade = function(tree){
	Nclades = length(tree$tip.label) + tree$Nnode;
	outgoing_edges = get_outgoing_edges_per_clade_CPP(	Ntips 		= length(tree$tip.label),
														Nnodes		= tree$Nnode,
														Nedges		= nrow(tree$edge),
														tree_edge	= as.vector(t(tree$edge)) - 1);
	
	# update indices from 0-based to 1-based
	return(lapply(1:Nclades,FUN = function(clade) outgoing_edges[[clade]]+1))
}


get_incoming_edges_per_clade = function(tree){
	Nclades = length(tree$tip.label) + tree$Nnode;
	incoming_edges = get_incoming_edges_per_clade_CPP(	Ntips 		= length(tree$tip.label),
														Nnodes		= tree$Nnode,
														Nedges		= nrow(tree$edge),
														tree_edge	= as.vector(t(tree$edge)) - 1);
	
	# update indices from 0-based to 1-based
	return(lapply(1:Nclades,FUN = function(clade) incoming_edges[[clade]]+1))
}



get_paths_root_to_tips = function(tree){
	Ntips = length(tree$tip.label)
	paths = get_paths_root_to_tips_CPP(	Ntips 		= Ntips,
										Nnodes		= tree$Nnode,
										Nedges		= nrow(tree$edge),
										tree_edge	= as.vector(t(tree$edge)) - 1);
	
	# update indices from 0-based to 1-based
	return(lapply(1:Ntips,FUN = function(tip) paths[[tip]]+1))
}



# type can be 'tip', 'node' or 'both'
map_tip_or_node_names_to_indices = function(tree, A, type, list_title, check_input=TRUE){
	if(type=='tip'){
		item_title 	= 'tip'
		Nmax_title 	= 'Ntips'
		Nmax 		= length(tree$tip.label)
		if(is.character(A)) name_pool = tree$tip.label;
	}else if(type=='node'){
		item_title 	= 'node';
		Nmax_title 	= 'Nnodes'
		Nmax 		= tree$Nnode;
		if(is.character(A)) name_pool = tree$node.label;
	}else{
		item_title 	= 'tip or node'
		Nmax_title 	= 'Ntips+Nnodes'
		Nmax 		= length(tree$tip.label)+tree$Nnode;
		if(is.character(A)) name_pool = c(tree$tip.label,tree$node.label);
	}
	if((!is.character(A)) && (!is.numeric(A))) stop(sprintf("ERROR: %s must be a character or integer vector",list_title))
	if(is.character(A)){
		name2index = 1:Nmax;
		names(name2index) = name_pool;
		Ai = name2index[A]; 
		if(check_input && any(is.na(Ai))) stop(sprintf("ERROR: Unknown %s name '%s'",item_title,A[which(is.na(Ai))[1]]));
		A = Ai;
	}else if(check_input){
		minA = min(A); maxA = max(A);
		if(minA<1 || maxA>Nmax) stop(sprintf("ERROR: %s must contain values between 1 and %s (%d); instead found values from %d to %d",list_title,Nmax_title,Nmax,minA,maxA));
	}
	return(A);
}



collapse_monofurcations = function(tree, force_keep_root=TRUE, as_edge_counts=FALSE){
	Ntips 	= length(tree$tip.label)
	Nnodes	= tree$Nnode
	Nedges	= nrow(tree$edge)
	
	results = get_tree_with_collapsed_monofurcations_CPP(	Ntips 			= length(tree$tip.label),
															Nnodes			= tree$Nnode,
															Nedges			= nrow(tree$edge),
															tree_edge		= as.vector(t(tree$edge)) - 1,
															edge_length		= (if(is.null(tree$edge.length) || as_edge_counts) numeric() else tree$edge.length),
															force_keep_root = force_keep_root);
	# reformat results into a valid "phylo" object
	# note that some of the old nodes may have turned into new tips
	Nnodes_new	 	= results$Nnodes_new
	new2old_node	= results$new2old_node + 1; # switch to 1-based indices
	collapsed_tree = list(	Nnode 		= Nnodes_new,
							tip.label 	= tree$tip.label,
							node.label 	= (if(is.null(tree$node.label)) NULL else tree$node.label[new2old_node]),
							edge 		= matrix(results$new_tree_edge,ncol=2,byrow=TRUE) + 1,
							edge.length = (if(is.null(tree$edge.length)) NULL else results$new_edge_length),
							root 		= results$new_root+1)
	class(collapsed_tree) = "phylo";
	
	return(list(tree			= collapsed_tree, 
				new2old_node	= new2old_node, 
				Nnodes_removed	= Nnodes-Nnodes_new));
}




# given a Markov transition rate matrix, calculate the transition probability matrix conditional upon a single transition occurring
# input: Q[i,j] is probability rate of transition i-->j
# output: P[i,j] will be the probability of transition i-->j, provided that a single transition i-->* occurred
get_conditional_transition_probabilities = function(Q){
	S = rowSums(Q)-diag(Q)
	S[S<=0] = 1
	P = Q/S
	diag(P) = 0
	return(P)
}



# extend the terminal edges (edges leading to tips) so that each tip has the same fixed distance (new_height) from the root
# if a tip already extends beyond the specified new_height, its incoming edge remains unchanged
# this is a quick-and-dirty way to make the tree ultrametric
# if new_height<0 or new_height==NULL, then it is set to the max_distance_to_root of the input tree
extend_tree_to_height = function(tree, new_height=NULL){
	results = extend_tree_to_height_CPP(	Ntips 			= length(tree$tip.label),
											Nnodes			= tree$Nnode,
											Nedges			= nrow(tree$edge),
											tree_edge		= as.vector(t(tree$edge)) - 1,
											edge_length		= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
											new_height 		= (if(is.null(new_height)) -1.0 else new_height));
	tree$edge.length = results$new_edge_length;
	return(list(tree=tree, max_extension=results$max_extension))
}



