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



find_edge_splitting_tree = function(tree, target_tips, is_rooted=FALSE){
	Ntips 	= length(tree$tip.label)
	Nnodes 	= tree$Nnode
	if(is.character(target_tips)){
		# target tips are given as names, not indices
		indices	= match(target_tips, tree$tip.label)
		if(any(is.na(indices))) stop(sprintf("ERROR: Some target_tips (e.g. '%s') were not found in the tree",target_tips[which(is.na(indices))[1]]))
		target_tips = indices
	}
	
	results = find_edge_splitting_tree_CPP(	Ntips				= Ntips,
											Nnodes				= tree$Nnode,
											Nedges				= nrow(tree$edge),
											tree_edge			= as.vector(t(tree$edge)) - 1,	# flatten in row-major format and adjust clade indices to 0-based
											is_rooted			= is_rooted,
											target_tips			= target_tips - 1,
											include_misplaced 	= TRUE)
													
	return(list(edge					= (if(results$edge<0) NA else as.integer(results$edge+1)),
				Nmisplaced_targets		= results$Nmisplaced_targets,
				Nmisplaced_nontargets	= results$Nmisplaced_nontargets,
				Ntargets_upstream 		= results$Ntargets_upstream,
				Ntargets_downstream 	= results$Ntargets_downstream,
				misplaced_targets		= results$misplaced_targets,
				misplaced_nontargets	= results$misplaced_nontargets));
				
}






