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




get_subtree_with_clades = function(	tree, 
									clades_to_keep = NULL, 	# integer vector listing tip/node indices to keep
									collapse_monofurcations = TRUE, 
									force_keep_root = FALSE, 
									keep_all_children_of_explicit_clades_to_keep = FALSE, 
									keep_all_tips_of_explicit_clades_to_keep = FALSE){ 
	Ntips  = length(tree$tip.label);
	Nnodes = tree$Nnode;
	
	results = get_subtree_with_specific_clades_CPP(	Ntips 					= Ntips,
													Nnodes 					= Nnodes,
													Nedges					= nrow(tree$edge),
													tree_edge				= as.vector(t(tree$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
													edge_length				= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
													clades_to_keep			= clades_to_keep-1,
													collapse_monofurcations = collapse_monofurcations,
													force_keep_root			= force_keep_root,
													keep_all_children_of_explicit_clades_to_keep 	= keep_all_children_of_explicit_clades_to_keep,
													keep_all_tips_of_explicit_clades_to_keep 		= keep_all_tips_of_explicit_clades_to_keep)
	Ntips_kept  	= results$Ntips_kept
	Nnodes_kept 	= results$Nnodes_kept
	new2old_clade 	= results$new2old_clade + 1L # switch to 1-based indices
	subtree = list(	Nnode 		= Nnodes_kept,
					tip.label 	= (if(is.null(tree$tip.label)) NULL else tree$tip.label[new2old_clade[1:Ntips_kept]]),
					node.label 	= (if(is.null(tree$node.label)) NULL else tree$node.label[new2old_clade[(Ntips_kept+1):(Ntips_kept+Nnodes_kept)]-Ntips]),
					edge 		= matrix(as.integer(results$new_tree_edge),ncol=2,byrow=TRUE) + 1L,
					edge.length = results$new_edge_length,
					root 		= results$new_root+1L,
					root.edge	= (if(force_keep_root && (!is.null(tree$root.edge))) tree$root.edge else NULL));
	class(subtree) = "phylo";
	attr(subtree,"order") = "none";
	
	return(list(tree 			= subtree,
				root_shift		= results$root_shift, # distance between old & new root (will always be non-negative)
				new2old_tip		= new2old_clade[1:Ntips_kept], 
				new2old_node	= new2old_clade[(Ntips_kept+1):(Ntips_kept+Nnodes_kept)]-Ntips));
}



# calculate the geometric placement of tips & nodes for plotting a tree as a phylogram (with the root on the left and tips on the right end, edges extend horizontally left to right)
get_phylogram_geometry = function(tree){
	results = get_phylogram_geometry_CPP(	Ntips 					= length(tree$tip.label),
											Nnodes 					= tree$Nnode,
											Nedges					= nrow(tree$edge),
											tree_edge				= as.vector(t(tree$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
											edge_length				= (if(is.null(tree$edge.length)) numeric() else tree$edge.length));
	return(list(clade_x	= results$clade_x,		# x-coordinates of tips & nodes
				clade_y	= results$clade_y,		# y-coordinates of tips & nodes
				root_y	= results$root_y,		# the y-coordinate of the root
				min_x	= results$min_x,		# the minimum x-coordinate of any clade (normally this is 0)
				max_x	= results$max_x,		# the maximum x-coordinate of any tip
				min_y	= results$min_y,		# the minimum y-coordinate of any tip
				max_y	= results$max_y));		# the maximum y-coordinate of any tip
}




# assign tips & nodes of a tree to groups, such that each group is monophyletic (a "taxon") represented by exactly one of given representative tips
# each representative tip is taken to represent a different taxon
# tip2taxon[n] or node2taxon[n] will be -1 if the tip/node could not be unambiguously assigned to a taxon (e.g., it contains multiple descending representatives)
assign_clades_to_taxa = function(tree, representative_tips){
	Ntips  = length(tree$tip.label);
	Nnodes = tree$Nnode;
	results = assign_clades_to_taxa_CPP(	Ntips 			= Ntips,
											Nnodes 			= tree$Nnode,
											Nedges			= nrow(tree$edge),
											tree_edge		= as.vector(t(tree$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
											representatives	= representative_tips-1);
	return(list(tip2taxon 	= results$clade2taxon[1:Ntips]+1L,
				node2taxon	= results$clade2taxon[(Ntips+1):(Ntips+Nnodes)]+1L));
}




# congruify trees (map nodes in the target tree to "equivalent" nodes in the reference tree)
# [Eastman et al (2013). Congruification: support for time scaling large phylogenetic trees. Methods in Ecology and Evolution. 4:688-691]
# mapping must be one of the following:
#	A 2D integer array of size NM x 2 (with NM<=TNtips), listing Ttips mapped to Rtips (mapping[m,1] --> mapping[m,2])
#	A 2D character array of size NM x 2 (with NM<=TNtips), listing Ttip names mapped to Rtip names (mapping[m,1] --> mapping[m,2])
#	A data frame of size NM x 1, whose row names are target tip labels and whose entries are either integers (R tip indices) or strings (R tip labels). This is the format used by geiger::congruify.phylo
#	A vector of size NM, whose names are target tip labels and whose entries are either integers (R tip indices) or strings (R tip labels).
congruify_trees = function(reference_tree, target_tree, mapping){
	TNtips = length(target_tree$tip.label)
	RNtips = length(reference_tree$tip.label)

	# re-format mapping if needed
	if(is.data.frame(mapping)){
		mapped_Ttips = rownames(mapping)
		mapped_Rtips = (if(is.numeric(mapping[,1])) mapping[,1] else as.character(mapping[,1]))
	}else if(is.vector(mapping)){
		mapped_Ttips = names(mapping)
		mapped_Rtips = (if(is.numeric(mapping)) mapping else as.character(mapping))
	}else{
		# assuming mapping is a 2D array of size NM x 2
		mapped_Ttips = mapping[,1]
		mapped_Rtips = mapping[,2]
	}
	if(is.character(mapped_Ttips)){
		# mapped target tips given as names, not indices
		indices = match(mapped_Ttips, target_tree$tip.label)
		if(any(is.na(indices))) stop(sprintf("ERROR: Some mapped target tips (e.g. '%s') were not found in the target tree",mapped_Ttips[which(is.na(indices))[1]]))
		mapped_Ttips = indices
	}
	if(is.character(mapped_Rtips)){
		# mapped reference tips given as names, not integers
		indices = match(mapped_Rtips, reference_tree$tip.label)
		if(any(is.na(indices))) stop(sprintf("ERROR: Some mapped reference tips (e.g. '%s') were not found in the reference tree",mapped_Rtips[which(is.na(indices))[1]]))
		mapped_Rtips = indices
	}
	mapping = matrix(c(mapped_Ttips,mapped_Rtips),ncol=2,byrow=FALSE)
		
	# congruify
	results = congruify_trees_CPP(	RNtips		= RNtips,
									RNnodes		= reference_tree$Nnode,
									RNedges		= nrow(reference_tree$edge),
									Rtree_edge	= as.vector(t(reference_tree$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
									TNtips		= TNtips,
									TNnodes		= target_tree$Nnode,
									TNedges		= nrow(target_tree$edge),
									Ttree_edge	= as.vector(t(target_tree$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
									mapping		= as.vector(t(mapping))-1) # flatten in row-major format and adjust tip indices to 0-based
	return(list(reference_nodes = results$mapped_Rnodes+1L,
				target_nodes 	= results$mapped_Tnodes+1L))
}



# map nodes in tree A to nodes in tree B, assuming both trees have the same topologies (but are potentially indexed differently)
# if tipsA2B is NULL, tips are matched by name
# This function returns an error (success=FALSE) if the trees don't have equivalent topologies, so it can also be used as a simple equivalence test
match_tree_nodes = function(treeA, treeB, tipsA2B=NULL){
	Ntips  = length(treeA$tip.label)
	Nnodes = treeA$Nnode
	if((Ntips!=length(treeB$tip.label)) || (Nnodes!=treeB$Nnode)) return(list(success=FALSE, error=sprintf("Tree sizes don't match: TreeA has %d tips and %d nodes, treeB has %d tips and %d nodes",Ntips,Nnodes,length(treeB$tip.label),treeB$Nnode)))
	if(is.null(tipsA2B)){
		tipsA2B = match(treeA$tip.label, treeB$tip.label)
		if(any(is.na(tipsA2B))) return(list(success=FALSE, error=sprintf("Tip labels in treeA don't match tip labels in treeB")))
	}
	results = match_tree_nodes_CPP(	Ntips		= Ntips,
									Nnodes		= Nnodes,
									Nedges		= nrow(treeA$edge),
									tree_edgeA	= as.vector(t(treeA$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
									tree_edgeB	= as.vector(t(treeB$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
									tipsA2B		= tipsA2B-1);
	if(!results$success) return(list(success=FALSE, error=results$error));
	return(list(success 	= TRUE,
				rootA		= results$rootA,
				rootB		= results$rootB,
				nodesA2B	= results$nodesA2B+1L));
}


get_complement = function(N, indices){
	include 		 = rep(TRUE, times=N)
	include[indices] = FALSE
	return(which(include))
}


# extract the values of independent rates from a transition_matrix, based on a provided index_matrix
# The index matrix should be as generated by get_transition_index_matrix, i.e. a 2D matrix of size Nstates x Nstates, and values in 1,..,Nrates, where Nrates is the number of independent rate variables
extract_independent_rates_from_transition_matrix = function(transition_matrix, index_matrix){
	flattened_index_matrix 	= as.vector(index_matrix)
	flattened_transition_matrix = as.vector(transition_matrix)
	independents 			= seq_along(flattened_index_matrix)[!duplicated(flattened_index_matrix)]
	independents 			= independents[flattened_index_matrix[independents]>0] # omit any zeros from index_matrix
	independent_rates 		= rep(NA,length(independents))
	independent_rates[flattened_index_matrix[independents]] = flattened_transition_matrix[independents]
	return(independent_rates)
}



# get a 1D lookup matrix of size Nstates, mapping birth-rates to indices of a subset of independent birth-rates
# model can be "ER" or "ARD" or a custom index_matrix as if it was generated by this function (in which case it is merely used to determine Nrates)
# This function is analogous to get_transition_index_matrix(..), but for 1D vectors
get_rate_index_vector = function(Nstates, rate_model){
	if (is.character(rate_model)) {
		if(rate_model == "ER"){
			Nrates = 1;
			index_vector = rep(1,Nstates)
		}else if(rate_model == "ARD"){
			Nrates = Nstates;
			index_vector = 1:Nrates;		
		}else{
			stop(sprintf("ERROR: Unknown rate_model '%s'",rate_model))
		}
	}else{
		if(length(rate_model)!=Nstates) stop(sprintf("ERROR: Wrong number of elements in rate model (expected %d, found %d)",Nstates,length(rate_model)));
		index_vector = rate_model
		Nrates = max(rate_model)
	}
	return(list(index_vector=index_vector, Nrates=Nrates))
}


extract_independent_rates_from_vector = function(rates, index_vector){
	independents 			= seq_along(index_vector)[!duplicated(index_vector)]
	independent_rates 		= rep(NA,length(independents))
	independent_rates[index_vector[independents]] = rates[independents]
	return(independent_rates)
}




# guesstimate an Mk transition matrix Q based on transitions along edges, as inferred via max-parsimony ASR
# Convention: Q[i,j] will be an estimate for the probability rate of the transition i-->j
# at least one of tip_states[] or tip_priors[] must be given; tip_states[] is given priority
guesstimate_Mk_transition_rates_via_max_parsimony_ASR = function(	tree, 
																	tip_states			= NULL,	# 1D vector of size Ntips, or NULL
																	tip_priors			= NULL,	# 2D array of size Ntips x Nstates, or NULL 
																	Nstates			 	= NULL, 
																	transition_costs 	= "all_equal"){
	# basic error checking & input formatting
	if(is.null(tip_states)){
		if(is.null(tip_priors)) return(list(success=FALSE, error="Missing tip_states or tip_priors"));
		tip_states  = max.col(tip_priors, ties.method="first")
		tip_states2 = max.col(tip_priors, ties.method="last")
		tip_states[tip_states!=tip_states2] = NA # in case of ambiguity, assign NA
	}
	if(length(tip_states)!=length(tree$tip.label)) return(list(success=FALSE, error=sprintf("Number of provided tip states (%d) does not match number of tips in the tree (%d)",length(tip_states),length(tree$tip.label))))
	
	# only consider subtree with known tip states
	known_tips = which(!is.na(tip_states));
	if(length(known_tips)<length(tip_states)){
		extraction	= get_subtree_with_tips(tree, only_tips=known_tips, omit_tips=FALSE, collapse_monofurcations=TRUE, force_keep_root=TRUE);
		tree		= extraction$subtree;
		tip_states	= tip_states[extraction$new2old_tip]
	}
	Ntips  = length(tree$tip.label)
	Nedges = nrow(tree$edge)
	
	# perform ASR max-parsimony on known subtree
	asr = asr_max_parsimony(	tree				= tree, 
								tip_states			= tip_states, 		
								Nstates				= Nstates, 
								transition_costs	= transition_costs);
	if(!asr$success) return(list(success=FALSE, error="ASR max-parsimony failed"))
	Nstates = ncol(asr$ancestral_likelihoods);
								
	# determine most likely node states
	node_states 	 	= max.col(asr$ancestral_likelihoods, ties.method="first")
	clade_states 		= c(tip_states, node_states)
	state_transitions 	= cbind(clade_states[tree$edge[,1]],clade_states[tree$edge[,2]]) # state_transitions[e,1]-->state_transitions[e,2] is the transition along edge e.
	
	# determine entries Q[i,j] based on transitions i-->j among all edges
	Q 			 = matrix(0,nrow=Nstates,ncol=Nstates)
	Ntransitions = matrix(0,nrow=Nstates,ncol=Nstates)
	edge_lengths = (if(is.null(tree$edge.length)) rep(1.0,Nedges) else tree$edge.length)
	for(i in 1:Nstates){
		for(j in 1:Nstates){
			transitions = which((state_transitions[,1]==i) & (state_transitions[,2]==j));
			Ntransitions[i,j] = length(transitions)
			if(i!=j){ # keep diagonal zero, for correct normalization afterwards
				if(length(transitions)>0){
					mean_transition_time = mean(edge_lengths[transitions])
					Q[i,j] = 1.0/mean_transition_time
				}else{
					Q[i,j] = 0;
				}
			}
		}
	}
	
	# make sure Q has zero sum in each row
	diag(Q) = -rowSums(Q, na.rm = TRUE);
	return(list(success=TRUE, Q=Q, Ntransitions=Ntransitions));
}



# return the ages of all branching events in an ultrametric timetree, i.e. all nodes accounting for multifurcations
# if the tree is purely bifurcating, this is the same as getting all node ages
# However, if the tree includes multifurcations, these are counted multiple times (since they represent multiple nearby bifurcations)
# Monofurcations are not returned
# Assumes that the tree is ultrametric
get_all_branching_ages = function(tree){
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode;
	depths = get_mean_depth_per_node_CPP(	Ntips			= Ntips,
											Nnodes			= Nnodes,
											Nedges			= nrow(tree$edge),
											tree_edge		= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
											edge_length		= (if(is.null(tree$edge.length)) numeric() else tree$edge.length));
	Nchildren = get_child_count_per_node_CPP(	Ntips			= Ntips,
												Nnodes			= Nnodes,
												Nedges			= nrow(tree$edge),
												tree_edge		= as.vector(t(tree$edge))-1);
	branch_ages = rep(depths,times=Nchildren-1);
	return(branch_ages);
}



# given a piecewise linear function f(x), defined on some x-grid, calculate its antiderivative A(x):=\int_{Xstart}^x f(u) du for an arbitrary number of target x values
# this function is most efficient when the requested target x-values are monotonically increasing or decreasing
get_antiderivative_of_piecewise_linear_function = function(	Xgrid,		# numeric vector of size NG, listing x-values in ascending order
															Xstart,		# numeric, lower end of the integration, i.e. x-value where antiderivative is set to zero
															Ygrid,		# numeric vector of size NG, listing y-values in ascending order
															splines_degree,	# integer, either 0,1,2 or 3, specifying the splines degree assumed for Y between grid points
															Xtarget){	# numeric vector of size N, specifying the target x values on which to evaluate the antiderivative
	A = get_antiderivative_CPP(	Xgrid, Xstart, Ygrid, splines_degree, Xtarget);
	return(A);
}



# given a lineages-through-time curve, defined as a time series on some discrete age grid, extract the branching ages that would have generated that LTT
# ages[] should be a 1D vector of ages (time before present) in ascending order
# LTT[] should be a 1D vector of the same size as ages[], listing the number of lineages at each age
# the LTT is assumed to be linear between adjacent age grid points
# branching points will be associated with those times where the LTT passes through an integer value
get_branching_ages_from_LTT = function(ages, LTT){
	results = get_branching_ages_from_LTT_CPP(ages, LTT);
	if(!results$success){
		return(list(success=FALSE, error=results$error));
	}else{
		return(list(success=TRUE, branching_ages = results$branching_ages));
	}
}





