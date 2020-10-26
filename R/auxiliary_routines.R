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
	Ntips = length(tree$tip.label)
	if(type=='tip'){
		item_title 	= 'tip'
		Nmax_title 	= 'Ntips'
		Nmax 		= Ntips
		if(is.character(A)) name_pool = tree$tip.label;
	}else if(type=='node'){
		item_title 	= 'node';
		Nmax_title 	= 'Nnodes'
		Nmax 		= tree$Nnode;
		if(is.character(A)) name_pool = tree$node.label;
	}else{
		item_title 	= 'tip or node'
		Nmax_title 	= 'Ntips+Nnodes'
		Nmax 		= Ntips+tree$Nnode;
		if(is.character(A)) name_pool = c(tree$tip.label,tree$node.label);
	}
	if((!is.character(A)) && (!is.numeric(A)) && (!is.integer(A))) stop(sprintf("ERROR: %s must be a character or integer vector",list_title))
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

get_intersection = function(N, indicesA, indicesB){
	includeA 		 	= rep(FALSE, times=N)
	includeB 		 	= rep(FALSE, times=N)
	includeA[indicesA] 	= TRUE
	includeB[indicesB] 	= TRUE
	return(which(includeA & includeB))
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
																	transition_costs 	= "all_equal",
																	allow_ties			= FALSE){	# logical, specifying whether to include tips with ties in the tip priors, i.e. tips that have the same tip_prior for various states. If TRUE, then ties are broken at random. If FALSE, then tips with ties are omitted from the analysis.
	# basic error checking & input formatting
	if(is.null(tip_states)){
		if(is.null(tip_priors)) return(list(success=FALSE, error="Missing tip_states or tip_priors"));
		if(allow_ties){
			tip_states = max.col(tip_priors, ties.method="random")
		}else{
			tip_states  = max.col(tip_priors, ties.method="first")
			tip_states2 = max.col(tip_priors, ties.method="last")
			tip_states[tip_states!=tip_states2] = NA # in case of ambiguity, assign NA
		}
	}
	if(length(tip_states)!=length(tree$tip.label)) return(list(success=FALSE, error=sprintf("Number of provided tip states (%d) does not match number of tips in the tree (%d)",length(tip_states),length(tree$tip.label))))
	
	# only consider subtree with known tip states
	known_tips = which(!is.na(tip_states));
	if(length(known_tips)<=1) return(list(success=FALSE, error=sprintf("All or almost all tips have unknown or ambiguous state")))
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



# guesstimate the scaled of Mk transition rates based on transitions along independent contrasts
# Attention: Strictly speaking, this function can only estimate the single rate of an "ER" Mk model
# at least one of tip_states[] or tip_priors[] must be given; tip_states[] is given priority
guesstimate_Mk_rate_via_independent_contrasts = function(	tree, 					# rooted phylogenetic tree of class "phylo"
															Nstates,
															tip_states	= NULL,		# 1D vector of size Ntips, or NULL
															tip_priors	= NULL,		# 2D array of size Ntips x Nstates, or NULL 
															allow_ties	= FALSE){	# logical, specifying whether to include tips with ties in the tip priors, i.e. tips that have the same tip_prior for various states. If TRUE, then ties are broken at random. If FALSE, then tips with ties are omitted from the analysis.
	# basic error checking & input formatting
	if(is.null(tip_states)){
		if(is.null(tip_priors)) return(list(success=FALSE, error="Missing tip_states or tip_priors"));
		if(allow_ties){
			tip_states = max.col(tip_priors, ties.method="random")
		}else{
			tip_states  = max.col(tip_priors, ties.method="first")
			tip_states2 = max.col(tip_priors, ties.method="last")
			tip_states[tip_states!=tip_states2] = NA # in case of ambiguity, assign NA
		}
	}
	if(length(tip_states)!=length(tree$tip.label)) return(list(success=FALSE, error=sprintf("Number of provided tip states (%d) does not match number of tips in the tree (%d)",length(tip_states),length(tree$tip.label))))
	
	# only consider subtree with known tip states
	known_tips = which(!is.na(tip_states));
	if(length(known_tips)<=1) return(list(success=FALSE, error=sprintf("All or almost all tips have unknown or ambiguous state")))
	if(length(known_tips)<length(tip_states)){
		extraction	= get_subtree_with_tips(tree, only_tips=known_tips, omit_tips=FALSE, collapse_monofurcations=TRUE, force_keep_root=TRUE);
		tree		= extraction$subtree;
		tip_states	= tip_states[extraction$new2old_tip]
	}
	Ntips  = length(tree$tip.label)
	Nedges = nrow(tree$edge)
	
	# fit a symmetric Mk model to the data
	fit = fit_symmetric_mk(	tree, Nstates = Nstates, tip_states	= tip_states, rate_model = "ER", Ntrials = 1)
	if(fit$success){	
		rate = fit$transition_matrix[1,2]
	}else{
		# fitting failed, so use first-fuess rate
		rate = fit$guess_rate
	}
	
	return(list(success=TRUE, rate=rate));
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
	branching_ages = rep(depths,times=Nchildren-1);
	return(branching_ages);
}



# given a piecewise polynomial (splines) function f(x), defined as a time series on some x-grid, calculate its antiderivative A(x):=\int_{Xstart}^x f(u) du for an arbitrary number of target x values
# this function is most efficient when the requested target x-values are monotonically increasing or decreasing
get_antiderivative_of_splines_function = function(	Xgrid,			# numeric vector of size NG, listing x-values in ascending order
													Xstart,			# numeric, lower end of the integration, i.e. x-value where antiderivative is set to zero
													Ygrid,			# numeric vector of size NG, listing y-values along Xgrid
													splines_degree,	# integer, either 0,1,2 or 3, specifying the splines degree assumed for Y between grid points
													Xtarget){		# numeric vector of size N, specifying the target x values on which to evaluate the antiderivative. The function is most efficient if Xtarget are in ascending or descending order.
	A = get_antiderivative_CPP(	Xgrid, Xstart, Ygrid, splines_degree, Xtarget);
	return(A);
}


# given a piecewise polynomial (splines) function f(x), defined as a time series on some x-grid, calculate its first derivative A(x):=df(x)/dx at an arbitrary number of target x values
# this function is most efficient when the requested target x-values are monotonically increasing or decreasing
get_derivative_of_splines_function = function(	Xgrid,			# numeric vector of size NG, listing x-values in ascending order
												Ygrid,			# numeric vector of size NG, listing y-values along Xgrid
												splines_degree,	# integer, either 0,1,2 or 3, specifying the splines degree assumed for Y between grid points
												Xtarget){		# numeric vector of size N, specifying the target x values on which to evaluate the derivative. The function is most efficient if Xtarget are in ascending or descending order.
	D = get_derivative_CPP(Xgrid, Ygrid, splines_degree, Xtarget);
	return(D);
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


# given some density curve on an X-interval, define a non-uniform X-grid on that interval so that the density of grid points reflects the requested density
# this can be used for example to define an age grid, with the grid density reflecting the number of lineages in a timetree at any given age, e.g. for fitting purposes
# the density curve is specified as a piecewise linear function. The density must be non-negative, and have non-zero total area under the curve.
get_inhomogeneous_grid_1D = function(	Xstart,
										Xend, 
										Ngrid, 		# integer, number of grid points to return, including the edges Xstart & Xend
										densityX, 	# numeric vector of size ND, listing X-values for defining the density, in ascending order
										densityY,	# numeric vector of size ND, listing density values at densityX. 
										extrapolate = FALSE){	# extrapolate density grid as needed, to cover Xstart & Xend. The density will be extrapolated as a constant.
	if(Ngrid<2){
		stop(sprintf("Ngrid must be at least 2"));
	}else if(densityX[1]>=tail(densityX,1)){
		stop(sprintf("Values in densityX must be strictly increasing"))
	}
	if(Xstart<densityX[1]){
		if(extrapolate){
			densityX = c(Xstart,densityX);
			densityY = c(densityY[1], densityY);
		}else{
			stop(sprintf("Xstart (%g) is not covered by the density grid (which starts at %g). Consider setting extrapolate=TRUE.",Xstart,densityX[1]));
		}
	}else if(Xend>tail(densityX,1)){
		if(extrapolate){
			densityX = c(densityX,Xend);
			densityY = c(densityY,tail(densityY,1));
		}else{
			stop(sprintf("Xend (%g) is not covered by the density grid (which ends at %g). Consider setting extrapolate=TRUE.",Xend,tail(densityX,1)));
		}
	}
	return(get_inhomogeneous_grid_1D_CPP(	Xstart		= Xstart, 
											Xend		= Xend, 
											Ngrid		= Ngrid, 
											densityX	= densityX, 
											densityY	= densityY, 
											xepsilon	= 0.0000001*(Xend-Xstart)/Ngrid));
}



# calculate the pulled speciation rate (PSR) of an HBD congruence class for a given pulled diversification rate (PDR) and product rho*lambda(0)
get_PSR_from_PDR_HBD = function(oldest_age,
								age_grid,				# numeric vector of size NG, listing grid ages in ascending order. Must cover at least age0 and oldest_age.
								PDR				= 0,	# numeric vector of size NG, listing PDRs on the corresponding grid points
								age0			= 0,	# non-negative numeric, specifying the age at which rholambda0 is given, i.e. rholambda0=rho(age0)*lambda(age0)
								rholambda0		= 1,	# positive numeric, product rho(age0)*lambda(age0), where rho is the sampling fraction and lambda is the speciation rate
								splines_degree	= 1,	# either 1, 2 or 3, specifying the degree of the splines defined by the PDR on the age grid.
								relative_dt		= 1e-3,	# numeric, maximum relative time step allowed for integration. Smaller values increase integration accuracy. Typical values are 0.0001-0.001.
								include_nLTT0	= FALSE){	# (logical) whether to also calculate the ratio nLTT0:=LTT(age0)/LTT(present-day)
	# basic error checking
	if(is.null(PDR)) stop("Missing PDR")
	if(is.null(age_grid) || (length(age_grid)<=1)){
		if((!is.null(PDR)) && (length(PDR)!=1)) return(list(success=FALSE, error=sprintf("Invalid number of PDR values; since no age grid was provided, you must provide a single (constant) PDR")))
		# create dummy age grid
		NG 			= 2;
		age_grid	= seq(from=0,to=oldest_age,length.out=NG)
		if(!is.null(PDR)) PDR = rep(PDR,times=NG);
	}else{
		NG = length(age_grid);
		if((age_grid[1]>oldest_age) || (age_grid[NG]<oldest_age)) return(list(success=FALSE, error=sprintf("Age grid must cover the entire requested age interval, including oldest_age (%g)",oldest_age)))
		if((age_grid[1]>age0) || (age_grid[NG]<age0)) return(list(success=FALSE, error=sprintf("Age grid must cover the entire requested age interval, including age0 (%g)",age0)))
		if(include_nLTT0 && (age_grid[1]>0) || (age_grid[NG]<0)) return(list(success=FALSE, error=sprintf("Age grid must cover the present-day age (0) in order to calculate nLTT0")))
		if((!is.null(PDR)) && (length(PDR)!=1) && (length(PDR)!=NG)) return(list(success=FALSE, error=sprintf("Invalid number of PDR values; since an age grid of size %d was provided, you must either provide one or %d PDR",NG,NG)))
		if((!is.null(PDR)) && (length(PDR)==1)) PDR = rep(PDR,times=NG);
	}
	if(rholambda0<=0) return(list(success=FALSE, error=sprintf("rholambda0 must be strictly positive; instead, got %g",rholambda0)))
	if(!(splines_degree %in% c(0,1,2,3))) return(list(success=FALSE, error=sprintf("Invalid splines_degree (%d): Expected one of 0,1,2,3.",splines_degree)))
	if(age_grid[1]>tail(age_grid,1)) return(list(success=FALSE, error=sprintf("Values in age_grid must be strictly increasing")))

	# calculate PSR
	results = get_PSR_from_PDR_HBD_CPP(	age0			= age0,
										oldest_age 		= oldest_age,
										age_grid		= age_grid,
										PDR				= PDR,
										rholambda0		= rholambda0,
										splines_degree	= splines_degree,
										relative_dt		= relative_dt,
										include_nLTT0	= include_nLTT0)
	if(results$success){
		return(list(success	= TRUE, 
					ages	= results$refined_age_grid,	# numeric vector listing (potentially refined) grid ages, spanning [max(0,age_grid[1]), oldest_age]
					PSR		= results$PSR, 				# numeric vector of the same size as ages[], listing the PSR on the refined grid
					nLTT0	= (if(include_nLTT0) results$nLTT0 else NULL)))
	}else{
		return(list(success=FALSE, error=results$error))
	}
}





# calculate the pulled speciation rate (PSR) of an HBD model for a given speciation rate (lambda), extinction rate (mu) and the sampling fraction (rho0) at some age0>=0
get_PSR_of_HBD_model = function(oldest_age,					# oldest age until which to calculate the PSR
								age_grid		= NULL,		# numeric vector of size NG, listing grid ages in ascending order. Must cover at least age0 and oldest_age. Can also be NULL, in which case the same lambda & mu are assumed everywhere.
								lambda			= 0,		# numeric vector of size NG, listing speciation rates on the corresponding grid points. Can also be a single constant.
								mu				= 0,		# numeric vector of size NG, listing extinction rates on the corresponding grid points. Can also be a single constant.
								age0			= 0,		# numeric, age (time before present) at which the sampling fraction (rho) is specified
								rho0			= 1,		# positive numeric, sampling fraction at age0
								splines_degree	= 1,		# either 1, 2 or 3, specifying the degree of the splines defined by the PDR on the age grid.
								relative_dt		= 1e-3){	# numeric, maximum relative time step allowed for integration. Smaller values increase integration accuracy. Typical values are 0.0001-0.001.
	# basic error checking
	if(is.null(age_grid) || (length(age_grid)<=1)){
		if((!is.null(lambda)) && (length(lambda)!=1)) return(list(success = FALSE, error = sprintf("Invalid number of lambda values (%d); since no age grid was provided, you must either provide a single (constant) lambda or none",length(lambda))))
		if((!is.null(mu)) && (length(mu)!=1)) return(list(success = FALSE, error = sprintf("Invalid number of mu values (%d); since no age grid was provided, you must provide a single (constant) mu",length(mu))))
		# create dummy age grid
		NG 		 = 2;
		age_grid = seq(from=0,to=1.01*oldest_age,length.out=NG)
		if(!is.null(lambda)) lambda = rep(lambda,times=NG);
		if(!is.null(mu)) mu = rep(mu,times=NG);
	}else{
		NG = length(age_grid);
		if((age_grid[1]>oldest_age) || (age_grid[NG]<oldest_age)) return(list(success = FALSE, error = sprintf("Age grid must cover the entire requested age interval, including oldest_age (%g)",oldest_age)))
		if((age_grid[1]>age0) || (age_grid[NG]<age0)) return(list(success = FALSE, error = sprintf("Age grid must cover the entire requested age interval, including age0 (%g)",age0)))
		if((!is.null(lambda)) && (length(lambda)!=1) && (length(lambda)!=NG)) return(list(success = FALSE, error = sprintf("Invalid number of lambda values (%d); since an age grid of size %d was provided, you must either provide zero, one or %d lambdas",length(lambda),NG,NG)))
		if((!is.null(mu)) && (length(mu)!=1) && (length(mu)!=NG)) return(list(success = FALSE, error = sprintf("Invalid number of mu values (%d); since an age grid of size %d was provided, you must either provide one or %d mus",length(mu),NG,NG)))
		if((!is.null(lambda)) && (length(lambda)==1)) lambda = rep(lambda,times=NG);
		if((!is.null(mu)) && (length(mu)==1)) mu = rep(mu,times=NG);
	}
	if(rho0<=0) return(list(success=FALSE, error=sprintf("rho must be strictly positive; instead, got %g",rho0)))
	if(!(splines_degree %in% c(0,1,2,3))) return(list(success=FALSE, error=sprintf("Invalid splines_degree (%d): Expected one of 0,1,2,3.",splines_degree)))
	if(age_grid[1]>tail(age_grid,1)) return(list(success=FALSE, error=sprintf("Values in age_grid must be strictly increasing")))

	# calculate PSR from lambda, mu & rho0=rho(age0)
	results = get_PSR_of_HBD_model_CPP(	age0			= age0,
										oldest_age 		= oldest_age,
										age_grid		= age_grid,
										lambda			= lambda,
										mu				= mu,
										rho0			= rho0,
										splines_degree	= splines_degree,
										relative_dt		= relative_dt)
	if(results$success){
		return(list(success	= TRUE, 
					ages	= results$refined_age_grid, # numeric vector listing (potentially refined) grid ages in ascending order, spanning [max(0,age_grid[1]), oldest_age]
					PSR		= results$PSR))	# numeric vector of the same size as ages[], listing the PSR on the refined grid
	}else{
		return(list(success=FALSE, error=results$error))
	}
}




# generate a random variable bounded from below but not from above
# used to pick random start params for fitting trials
random_semiboxed_left = function(lower_bound, default, typical_scale, orders_of_magnitude){
	if((default==0) && (lower_bound==0)){
		if(typical_scale==0){
			return(runif(n=1,min=0,max=1))
		}else{
			return(abs(typical_scale) * 10**runif(n=1,min=0,max=orders_of_magnitude))
		}
	}else if((default>lower_bound) && (default>0)){
		return(if(rbinom(n=1,size=1,prob=0.5)==1) (default * 10**runif(n=1,min=0,max=orders_of_magnitude)) else (lower_bound + (default-lower_bound)*runif(n=1,min=0,max=1)))
	}else if(default>lower_bound){
		return(lower_bound + (default-lower_bound) * 10**runif(n=1,min=-orders_of_magnitude/2,max=orders_of_magnitude/2))
	}else{
		return(lower_bound + (0-lower_bound) * 10**runif(n=1,min=-orders_of_magnitude/2,max=orders_of_magnitude/2))
	}
}


# choose random parameter values within boxed constraints
# each lower & upper bound may be Inf
# defaults[], lower_bounds[], upper_bounds[] and scales[] must be 1D numeric vectors of the same length, and must not include NaN or NA 
# lower_bounds[] and upper_bounds[] may include +Inf and -Inf
get_random_params = function(defaults, lower_bounds, upper_bounds, scales, orders_of_magnitude){
	start_values = defaults
	boxed_left	 = which((!is.infinite(lower_bounds)) & is.infinite(upper_bounds))
	boxed_right	 = which((!is.infinite(upper_bounds)) & is.infinite(lower_bounds))
	boxed_dual   = which(!(is.infinite(lower_bounds) | is.infinite(upper_bounds))); # determine fitted params that are boxed, i.e. constrained to within finite lower & upper bounds
	unboxed 	 = which(is.infinite(lower_bounds) & is.infinite(upper_bounds))
	if(length(boxed_dual)>0) 	start_values[boxed_dual] = lower_bounds[boxed_dual] + (upper_bounds[boxed_dual]-lower_bounds[boxed_dual]) * runif(n=length(boxed_dual),min=0,max=1)
	if(length(unboxed)>0) 	 	start_values[unboxed]	 = 10**runif(n=length(unboxed), min=-orders_of_magnitude/2.0, max=orders_of_magnitude/2.0) * start_values[unboxed]
	if(length(boxed_left)>0) 	start_values[boxed_left] = sapply(boxed_left, FUN=function(fp) random_semiboxed_left(lower_bound=lower_bounds[fp], default=start_values[fp], typical_scale=scales[fp], orders_of_magnitude=orders_of_magnitude))
	if(length(boxed_right)>0) 	start_values[boxed_right]= sapply(boxed_right, FUN=function(fp) -random_semiboxed_left(lower_bound=-upper_bounds[fp], default=-start_values[fp], typical_scale=scales[fp], orders_of_magnitude=orders_of_magnitude))
	start_values = pmax(lower_bounds,pmin(upper_bounds,start_values))
	return(start_values)
}


# check validity and sanitize (standardize format) of model parameters for fitting
# used by various fitting routines
sanitize_parameters_for_fitting = function(	param_values,				# numeric vector of size NP, specifying fixed values for a some or all parameters. For fitted (i.e. non-fixed) parameters, use NaN or NA.
											param_guess		= NULL,		# numeric vector of size NP, listing an initial guess for each parameter. For fixed parameters, guess values are ignored.
											param_min		= -Inf,		# numeric vector of size NP, specifying lower bounds for the model parameters. For fixed parameters, bounds are ignored. May also be a single scalar, in which case the same lower bound is assumed for all params.
											param_max		= +Inf,		# numeric vector of size NP, specifying upper bounds for the model parameters. For fixed parameters, bounds are ignored. May also be a single scalar, in which case the same upper bound is assumed for all params.
											param_scale		= NULL){	# numeric vector of size NP, specifying typical scales for the model parameters. For fixed parameters, scales are ignored. If NULL, scales are automatically estimated from other information (such as provided guess and bounds). May also be a single scalar, in which case the same scale is assumed for all params.
	NP 			= length(param_values);
	param_names	= names(param_values);
	if(is.null(param_guess)){
		if(any(!is.finite(param_values))){
			return(list(success=FALSE, error=sprintf("Missing guessed parameter values")))
		}else{
			param_guess = rep(NA, times=NP);
		}
	}
	if(length(param_guess)!=NP){
		return(list(success=FALSE, error=sprintf("Number of guessed parameters (%d) differs from number of model parameters (%d)",length(param_guess),NP)))
	}else if(!is.null(param_names)){
		names(param_guess) = param_names;
	}
	if((!is.null(param_names)) && (length(param_names)!=NP)){
		return(list(success=FALSE, error=sprintf("Number of parameter names (%d) differs from number of model parameters (%d)",length(param_names),NP)))
	}
	if(is.null(param_min)){
		param_min = rep(-Inf,times=NP);
	}else if(length(param_min)==1){
		param_min = rep(param_min,times=NP);
	}else if(length(param_min)!=NP){
		return(list(success=FALSE, error=sprintf("Length of param_min[] (%d) differs from number of model parameters (%d)",length(param_min),NP)))
	}
	if(is.null(param_max)){
		param_max = rep(+Inf,times=NP);
	}else if(length(param_max)==1){
		param_max = rep(param_max,times=NP);
	}else if(length(param_max)!=NP){
		return(list(success=FALSE, error=sprintf("Length of param_max[] (%d) differs from number of model parameters (%d)",length(param_max),NP)))
	}
	if(is.null(param_scale)){
		param_scale = rep(NA,times=NP);
	}else if(length(param_scale)==1){
		param_scale = rep(param_scale,times=NP);
	}else if(length(param_scale)!=NP){
		return(list(success=FALSE, error=sprintf("Length of param_scale[] (%d) differs from number of model parameters (%d)",length(param_scale),NP)))
	}
	param_values[is.nan(param_values)] = NA # standardize representation of non-fixed params
	param_scale[is.nan(param_scale)] = NA	# standardize representation of unknown param scales		
		
	# determine which parameters are to be fitted
	fitted_params	= which(is.na(param_values))
	fixed_params	= which(!is.na(param_values))
	NFP				= length(fitted_params);
	param_guess[fixed_params] = param_values[fixed_params] # make sure guessed param values are consistent with fixed param values

	if(any(!is.finite(param_guess))) return(list(success=FALSE, error=sprintf("Some guessed parameter values are NA or NaN or Inf; you must specify a valid guess for each non-fixed model parameter")));
	if(any((!is.na(param_scale)) & (param_scale==0))) return(list(success=FALSE, error=sprintf("Some provided parameter scales are zero; expecting non-zero scale for each parameter")));
	
	# determine typical parameter scales, whenever these are not provided
	for(p in fitted_params){
		if(is.na(param_scale[p])){
			if(param_guess[p]!=0){
				param_scale[p] = abs(param_guess[p]);
			}else if((is.finite(param_min[p]) && (param_min[p]!=0)) || (is.finite(param_max[p]) && (param_max[p]!=0))){
				param_scale[p] = mean(abs(c((if(is.finite(param_min[p]) && (param_min[p]!=0)) param_min[p] else NULL), (if(is.finite(param_max[p]) && (param_max[p]!=0)) param_max[p] else NULL))));
			}else{
				param_scale[p] = 1;
			}
		}
	}
	
	return(list(success			= TRUE,
				NP				= NP,
				NFP				= NFP,
				param_names		= param_names,
				param_values	= param_values,
				param_guess		= param_guess,
				param_min		= param_min,
				param_max		= param_max,
				param_scale		= param_scale,
				fitted_params	= fitted_params,
				fixed_params	= fixed_params))
}
	



# function for reformatting/sanitizing generic fit parameters provided by the user
# used by various fitting routines
# prepare_generic_fit_params = function(	param_values,					# numeric vector of size NP, specifying fixed values for a some or all parameters. For fitted (i.e. non-fixed) parameters, use NaN or NA.
# 										param_guess			= NULL,		# numeric vector of size NP, listing an initial guess for each parameter. For fixed parameters, guess values are ignored.
# 										param_min			= -Inf,		# numeric vector of size NP, specifying lower bounds for the model parameters. For fixed parameters, bounds are ignored. May also be a single scalar, in which case the same lower bound is assumed for all params.
# 										param_max			= +Inf,		# numeric vector of size NP, specifying upper bounds for the model parameters. For fixed parameters, bounds are ignored. May also be a single scalar, in which case the same upper bound is assumed for all params.
# 										param_scale			= NULL){	# numeric vector of size NP, specifying typical scales for the model parameters. For fixed parameters, scales are ignored. If NULL, scales are automatically estimated from other information (such as provided guess and bounds). May also be a single scalar, in which case the same scale is assumed for all params.
# 	NP = length(param_values)
# 	param_names = names(param_values);
# 	if(is.null(param_guess)){
# 		if(any(is.finite(param_values))){
# 			return(list(success=FALSE, error=sprintf("Missing guessed parameter values")))
# 		}else{
# 			param_guess = rep(NA, times=NP);
# 		}
# 	}
# 	if(length(param_guess)!=NP){
# 		return(list(success=FALSE, error=sprintf("Number of guessed parameters (%d) differs from number of model parameters (%d)",length(param_guess),NP)))
# 	}else if(!is.null(param_names)){
# 		names(param_guess) = param_names;
# 	}
# 	if((!is.null(param_names)) && (length(param_names)!=NP)){
# 		return(list(success=FALSE, error=sprintf("Number of parameter names (%d) differs from number of model parameters (%d)",length(param_names),NP)))
# 	}
# 	if(is.null(param_min)){
# 		param_min = rep(-Inf,times=NP);
# 	}else if(length(param_min)==1){
# 		param_min = rep(param_min,times=NP);
# 	}else if(length(param_min)!=NP){
# 		return(list(success=FALSE, error=sprintf("Length of param_min[] (%d) differs from number of model parameters (%d)",length(param_min),NP)))
# 	}
# 	if(is.null(param_max)){
# 		param_max = rep(+Inf,times=NP);
# 	}else if(length(param_max)==1){
# 		param_max = rep(param_max,times=NP);
# 	}else if(length(param_max)!=NP){
# 		return(list(success=FALSE, error=sprintf("Length of param_max[] (%d) differs from number of model parameters (%d)",length(param_max),NP)))
# 	}
# 	if(is.null(param_scale)){
# 		param_scale = rep(NA,times=NP);
# 	}else if(length(param_scale)==1){
# 		param_scale = rep(param_scale,times=NP);
# 	}else if(length(param_scale)!=NP){
# 		return(list(success=FALSE, error=sprintf("Length of param_scale[] (%d) differs from number of model parameters (%d)",length(param_scale),NP)))
# 	}
# 	return(list(success=TRUE, param_values=param_values, param_guess=param_guess, param_min=param_min, param_max=param_max, param_scale))
# }


# given an undirected graph (nodes,edges), find its maximal connected subgraphs
# any two nodes may be connected by zero, one or multiple edges
# edges[] should be a 2D array of size Nedges x 2, listing source & target nodes of the graph
get_connected_subgraphs = function(Nnodes, edges){
	results = split_undirected_graph_CPP(Nnodes=Nnodes, Nedges=nrow(edges), edges = as.vector(t(edges))-1);
	return(list(Nsubgraphs 		= results$Nsubgraphs,
				subgraph2nodes	= lapply(1:results$Nsubgraphs, FUN=function(n) results$subgraph2nodes[[n]]+1),
				subgraph2edges	= lapply(1:results$Nsubgraphs, FUN=function(n) results$subgraph2edges[[n]]+1),
				node2subgraph	= results$node2subgraph+1,
				edge2subgraph	= results$edge2subgraph+1));
}



# Sarting from a mapping pool-->group, calculate the reverse mapping group-->member_list
# this function returns a list of length Ngroups, where each element is a 1D vector listing member indices
# values in pool2group[] must be between 1 and Ngroups, or NA or <=0 (in which case the item is not affiliated with any group)
get_member_lists_from_group_assignments = function(Ngroups, pool2group){
	if(Ngroups==0) return(list())
	pool2group[is.na(pool2group)] = -1;
	results = get_member_lists_from_group_assignments_CPP(Ngroups=Ngroups, pool2group=pool2group-1)
	return(lapply(1:Ngroups,FUN=function(g) results$group2members[[g]]+1))
}



# evaluate a mathematical expression (univariate function of X) for various X-values
# the input X[] can either be a 1D vector or a 2D matrix
evaluate_univariate_expression = function(expression, Xname="x", X){
	if(is.vector(X)){
		results = evaluate_univariate_expression_CPP(expression=expression, Xname=Xname, X=X);
		return(list(success=results$success, error=results$error, Y=results$Y))
	}else if(is.matrix(X)){
		results = evaluate_univariate_expression_CPP(expression=expression, Xname=Xname, X=as.vector(t(X)));
		if(!results$success) return(list(success=FALSE, error=results$error))
		return(list(success=TRUE, Y=matrix(results$Y,ncol=ncol(X),byrow=TRUE)))
	}else{
		return(list(success=FALSE, error="Unknown data format X: Expecting either a vector or a matrix"))
	}
}


# Split tree into pairs of sister-tips, such that the paths within distinct pairs do not overlap
# If the input tree only contains monofurcations and bifurcations (recommended), it is guaranteed that at most one unpaired tip will be left (i.e., if Ntips was odd)
extract_independent_sister_tips = function(tree){
	results = extract_independent_sister_tips_CPP(	Ntips		= length(tree$tip.label),
													Nnodes		= tree$Nnode,
													Nedges		= nrow(tree$edge),
													tree_edge	= as.vector(t(tree$edge))-1);
	tip_pairs = matrix(as.integer(results$tip_pairs),ncol=2,byrow=TRUE) + 1L;
	return(tip_pairs);
}


# calculate geodesic angle (aka. central angle, in radians) between two geographical locations (assuming the Earth is a sphere)
# based on the Vincenty formula with equal major and minor axis
# geographic coordinates should be given in decimal degrees
geodesic_angle = function(latitude1, longitude1, latitude2, longitude2){
	theta1	= pi*latitude1/180;
	theta2	= pi*latitude2/180;
	phi1 	= pi*longitude1/180.0;
	phi2 	= pi*longitude2/180.0;
	delta	= abs(phi1-phi2);
	angle 	= abs(atan2(sqrt((cos(theta2)*sin(delta))^2 + (cos(theta1)*sin(theta2)-sin(theta1)*cos(theta2)*cos(delta))^2), (sin(theta1)*sin(theta2)+cos(theta1)*cos(theta2)*cos(delta))));
	return(angle);
}

# calculate geodesic angles (in radians) between pairs of coordinates
# this function returns a 1D vector of size equal to the size of the input lists latitudes1[], latitudes2[] etc
geodesic_angles = function(latitudes1, longitudes1, latitudes2, longitudes2){
	return(geodesic_angles_CPP(latitudes1, longitudes1, latitudes2, longitudes2));
}

# calculate all geodesic angles (in radians) between two sets of coordinates
# This function returns a 2D matrix of size N1 x N2
all_pairwise_geodesic_angles = function(latitudes1, longitudes1, latitudes2, longitudes2){
	angles = get_all_pairwise_geodesic_angles_CPP(latitudes1, longitudes1, latitudes2, longitudes2);
	return(matrix(angles, ncol=length(latitudes2), byrow=TRUE));
}


# calculate expectation of sin^2(omega) across multiple independent contrasts with non-equal time steps, where omega is the central transition angle, for the spherical Brownian Motion
get_expected_SBM_sinsquare = function(time_steps, diffusivity, radius){
	return((2/3)*(1-mean(exp(-6*(diffusivity/(radius^2))*time_steps))))
}

get_expected_SBM_transition_angle = function(time_step, diffusivity, radius){
	tD = time_step * diffusivity/radius^2;
	return(SBM_get_average_transition_angle_CPP(tD = tD, max_error = 1e-7, max_Legendre_terms = 200))
}


# autocorrelation function ACF(t) of Spherical Brownian Motion with constant diffusivity D
# ACF(t) := E <n(t),n(0)> = E cos(omega(t))
# where n(t) is the unit vector pointing to the particle's random location on the unit phere at time t, and <,> is the scalar product, and omega is the transition angle.
get_SBM_ACF = function(times, diffusivity, radius){
	return(exp(-2*(diffusivity/radius^2)*times))
}



# given some rooted timetree, extract information about various sampling and branching events, as defined by the homogenous birth-death-sampling model with continuous sampling and concentrated sampling efforts
# monofurcating nodes are interpreted as sampled nodes, i.e. sampled lineages with additional subsequently sampled descendants
# bifurcating and multifurcating nodes are interpreted as branching events, with multifurcations counted multiple times (i.e., as if they are first split into bifurcations)
extract_HBDS_events_from_tree = function(	tree,
											root_age = NULL, # optional numeric, the age of the root. Can be used to define a time offset, e.g. if the last tip was not actually sampled at the present. If NULL, this will be calculated from the treem and it will be assumed that the last tip was sampled at the present
											CSA_ages = numeric(0), # 1D vector listing ages of concentrated sampling efforts, in ascending order
											age_epsilon = 0){
	Ntips  	= length(tree$tip.label)
	Nnodes 	= tree$Nnode
	NCSA		= length(CSA_ages)
	if((NCSA>=2) && (CSA_ages[1]>tail(CSA_ages,1))) return(list(success=FALSE, error="CSA_ages must be in ascending order"))
	
	# determine branching ages & tip ages & sampled node ages
	if(is.null(root_age)) root_age = get_tree_span(tree)$max_distance
	node2Nchildren = get_child_count_per_node_CPP(	Ntips		= Ntips,
													Nnodes		= Nnodes,
													Nedges		= nrow(tree$edge),
													tree_edge	= as.vector(t(tree$edge))-1);
	clade_heights = get_distances_from_root_CPP(Ntips		= Ntips,
												Nnodes		= Nnodes,
												Nedges		= nrow(tree$edge),
												tree_edge	= as.vector(t(tree$edge))-1,
												edge_length	= (if(is.null(tree$edge.length)) numeric() else tree$edge.length))
	sampled_nodes		= which(node2Nchildren==1)
	sampled_node_ages 	= pmax(0, root_age - clade_heights[Ntips + sampled_nodes])
	tip_ages 			= pmax(0, root_age - clade_heights[1:Ntips])
	branchings			= which(node2Nchildren>=2) # determine which nodes are branchings
	branching_ages		= pmax(0, root_age - clade_heights[Ntips + branchings])
	branching_ages		= rep(branching_ages,times=node2Nchildren[branchings]-1) # count multifurcations multiple times
	
	if(length(CSA_ages)>0){
		# determine which tips were sampled during concentrated sampling efforts versus due to Poissonian sampling
		sorted2original = order(tip_ages)
		CSAbinning		= place_sorted_values_into_bins_CPP(items=tip_ages[sorted2original], bin_mins=CSA_ages-age_epsilon, bin_maxs=CSA_ages+age_epsilon)
		tip2CSA 		= 1L + CSAbinning$item2bin
		Ptip_ages		= tip_ages[sorted2original[which(CSAbinning$item2bin<0)]]
		CSA2tips		= lapply(1:NCSA, FUN=function(ce) sorted2original[1L + CSAbinning$bin2items[[ce]]]);
		CSA_tip_counts 	= sapply(1:NCSA, FUN=function(ce) length(CSA2tips[[ce]]))

		# determine which sampled nodes were sampled during concentrated sampling efforts versus due to Poissonian sampling
		sorted2original = order(sampled_node_ages)
		CSAbinning 		= place_sorted_values_into_bins_CPP(items=sampled_node_ages[sorted2original], bin_mins=CSA_ages-age_epsilon, bin_maxs=CSA_ages+age_epsilon)
		Pnode_ages		= sampled_node_ages[sorted2original[which(CSAbinning$item2bin<0)]]
		CSA2nodes		= lapply(1:NCSA, FUN=function(ce) sampled_nodes[sorted2original[1L + CSAbinning$bin2items[[ce]]]]);
		CSA_node_counts = sapply(1:NCSA, FUN=function(ce) length(CSA2nodes[[ce]]))
	}else{
		CSA2tips 		= integer(0)
		CSA_tip_counts 	= integer(0)
		CSA2nodes 		= integer(0)
		CSA_node_counts = integer(0)
		Ptip_ages 		= sort(tip_ages)
		Pnode_ages 		= sort(sampled_node_ages)
	}
	
	return(list(concentrated_tips 			= CSA2tips,
				concentrated_tip_counts 	= CSA_tip_counts,
				concentrated_nodes 			= CSA2nodes,
				concentrated_node_counts 	= CSA_node_counts,
				Ptip_ages					= Ptip_ages, 	# ages of Poissonian tip ages, i.e. of Poissonian terminal sampling events, in ascending order
				Pnode_ages					= Pnode_ages,	# ages of Poissonian node ages, i.e. of Poissonian non-terminal sampling events, in ascending order
				branching_ages				= sort(branching_ages)))
}



# Given some HBDS congruence class (in terms of PDR, PSR lambda_psi and CSA_pulled_probs), as well as some desired profile for a specific variable (e.g., the sampling rate psi), calculate the corresponding model in the congruence class
# The model in the congruence class must be specified through exactly one of the following variable choices:
# 	lambda (speciation rate)
#	mu (extinction rate) and lambda0 (present-day speciation rate)
#	psi (sampling rate)
#	Reff (effective reproduction ratio) and lambda0
#	removal_rate (aka. become-uninfectious rate, mu+psi) and lambda0
# ATTENTION: This function is currently only implemented for HBDS models without retention, i.e. where kappa=0
# ATTENTION: Not all psi are allowed. In particular, psi(CSA_ages[k]) must be the same for all models in the congruence class whenever CSA_ages[k]>0.
get_congruent_hbds = function(	age_grid,						# numeric vector of size NG, listing discrete ages in ascending order
								PSR,							# numeric vector of size NG, listing the PSR on the age_grid
								PDR,							# numeric vector of size NG, listing the PDR on the age_grid
								lambda_psi,						# numeric vector of size NG, listing lambda*psi on the age_grid
								lambda				= NULL,		# numeric vector of size NG. Exactly one of lambda, psi, mu, Reff or removal_rate must be provided.
								mu					= NULL,		# numeric vector of size NG. Exactly one of psi, mu, Reff or removal_rate must be provided.
								psi					= NULL,		# numeric vector of size NG. Exactly one of psi, mu, Reff or removal_rate must be provided.
								Reff				= NULL,		# numeric vector of size NG. Exactly one of psi, mu, Reff or removal_rate must be provided.
								removal_rate		= NULL,		# numeric vector of size NG. Exactly one of psi, mu, Reff or removal_rate must be provided.
								lambda0				= NULL,		# numeric, specifying lambda at age 0 (present-day). Only relevant if mu or Reff or removal_rate is provided.
								CSA_ages			= NULL,		# numeric vector of size NCSA, listing the ages of concentrated sampling attempts, in ascending order
								CSA_pulled_probs 	= NULL,		# numeric vector of size NCSA, listing the pulled probabilities of concentrated sampling attempts, in ascending order
								CSA_PSRs			= NULL, 	# numeric vector of size NCSA, listing the PSR at the concentrated sampling attempts, in ascending order
								splines_degree		= 1, 		# integer between 1 and 3. 0-splines are not allowed, because intervally the 1st derivatives are needed.
								ODE_relative_dt		= 0.001,
								ODE_relative_dy		= 1e-4){
	# basic error checking
	NCSA = (if(is.null(CSA_ages)) 0 else length(CSA_ages))
	if((NCSA==0) && (!is.null(CSA_pulled_probs)) && (length(CSA_pulled_probs)>0)) return(list(success=FALSE, error="No CE ages were provided, but CSA_pulled_probs were"))
	if((NCSA>0) && is.null(CSA_pulled_probs)) return(list(success=FALSE, error="Missing CSA_pulled_probs"))
	if((NCSA>0) && (length(CSA_pulled_probs)!=NCSA)) return(list(success=FALSE, error=sprintf("Expected %d CSA_pulled_probs, but instead got %d",NCSA,length(CSA_pulled_probs))))
	if((NCSA>0) && is.null(CSA_PSRs)) return(list(success=FALSE, error="Missing CSA_PSRs"))
	if((NCSA>0) && (length(CSA_PSRs)!=NCSA)) return(list(success=FALSE, error=sprintf("Expected %d CSA_PSRs, but instead got %d",NCSA,length(CSA_PSRs))))
	if(is.null(CSA_ages)){
		CSA_ages 			= numeric(0)
		CSA_pulled_probs 	= numeric(0)
		CSA_PSRs 			= numeric(0)
	}
	NG = length(age_grid)
	if(is.null(PSR)){
		PSR = rep(0,times=NG)
	}else if(length(PSR)==1){
		PSR = rep(PSR,times=NG)
	}else if(length(PSR)!=NG){
		return(list(success=FALSE, error=sprintf("Expected %d PSR values, but instead got %d",NG,length(PSR))))
	}
	if(is.null(PDR)){
		PDR = rep(0,times=NG)
	}else if(length(PDR)==1){
		PDR = rep(PDR,times=NG)
	}else if(length(PDR)!=NG){
		return(list(success=FALSE, error=sprintf("Expected %d PDR values, but instead got %d",NG,length(PDR))))
	}
	if(is.null(lambda_psi)){
		lambda_psi = rep(0,times=NG)
	}else if(length(psi)==1){
		lambda_psi = rep(lambda_psi,times=NG)
	}else if(length(lambda_psi)!=NG){
		return(list(success=FALSE, error=sprintf("Expected %d lambda_psi values, but instead got %d",NG,length(lambda_psi))))
	}
	if(is.null(lambda) && is.null(mu) && is.null(psi) && is.null(Reff) && is.null(removal_rate)) return(list(success=FALSE, error=sprintf("Expecting either lambda, mu, psi, Reff or removal_rate")))
	if(sum(!c(is.null(lambda),is.null(mu),is.null(psi),is.null(Reff),is.null(removal_rate)))>1) return(list(success=FALSE, error=sprintf("Only one of lambda, mu, psi, Reff or removal_rate must be provided")))
	if((!is.null(lambda)) && (!is.null(lambda0))) return(list(success=FALSE, error=sprintf("lambda0 must not be provided if lambda is provided")))
	if((!is.null(mu)) && is.null(lambda0)) return(list(success=FALSE, error=sprintf("lambda0 must be provided when mu is provided")))
	if((!is.null(psi)) && (!is.null(lambda0))) return(list(success=FALSE, error=sprintf("lambda0 must not be provided if psi is provided")))
	if((!is.null(Reff)) && is.null(lambda0)) return(list(success=FALSE, error=sprintf("lambda0 must be provided when Reff is provided")))
	if((!is.null(removal_rate)) && is.null(lambda0)) return(list(success=FALSE, error=sprintf("lambda0 must be provided when removal_rate is provided")))
	if(!is.null(lambda)){
		if(length(lambda)==1){
			lambda = rep(lambda,times=NG)
		}else if(length(lambda)!=NG){
			return(list(success=FALSE, error=sprintf("Expected %d lambda values, but instead got %d",NG,length(lambda))))
		}
		if(NCSA>0) return(list(success=FALSE, error=sprintf("Providing lambda to define a model is only available in the absence of CSAs")))
	}
	if(!is.null(mu)){
		if(length(mu)==1){
			mu = rep(mu,times=NG)
		}else if(length(mu)!=NG){
			return(list(success=FALSE, error=sprintf("Expected %d mu values, but instead got %d",NG,length(mu))))
		}
		if(NCSA>0) return(list(success=FALSE, error=sprintf("Providing mu to define a model is only available in the absence of CSAs")))
	}
	if(!is.null(psi)){
		if(length(psi)==1){
			psi = rep(psi,times=NG)
		}else if(length(psi)!=NG){
			return(list(success=FALSE, error=sprintf("Expected %d psi values, but instead got %d",NG,length(psi))))
		}
	}
	if(!is.null(Reff)){
		if(length(Reff)==1){
			Reff = rep(Reff,times=NG)
		}else if(length(Reff)!=NG){
			return(list(success=FALSE, error=sprintf("Expected %d Reff values, but instead got %d",NG,length(Reff))))
		}
		if(NCSA>0) return(list(success=FALSE, error=sprintf("Providing Reff to define a model is only available in the absence of CSAs")))
	}
	if(!is.null(removal_rate)){
		if(length(removal_rate)==1){
			removal_rate = rep(removal_rate,times=NG)
		}else if(length(removal_rate)!=NG){
			return(list(success=FALSE, error=sprintf("Expected %d removal_rate values, but instead got %d",NG,length(removal_rate))))
		}
		if(NCSA>0) return(list(success=FALSE, error=sprintf("Providing removal_rate to define a model is only available in the absence of CSAs")))
	}
	if(!(splines_degree %in% c(1,2,3))) return(list(success = FALSE, error = sprintf("Invalid splines_degree (%d): Expected one of 1,2,3.",splines_degree)))

	results = get_congruent_HBDS_CPP(	CSA_ages			= CSA_ages,
										CSA_pulled_probs	= CSA_pulled_probs,
										CSA_PSRs			= CSA_PSRs,
										age_grid			= age_grid,
										PSRs				= PSR,
										PDRs				= PDR,
										lambda_psis			= lambda_psi,
										lambdas				= (if(is.null(lambda)) numeric(0) else lambda),
										mus					= (if(is.null(mu)) numeric(0) else mu),
										psis				= (if(is.null(psi)) numeric(0) else psi),
										Reffs				= (if(is.null(Reff)) numeric(0) else Reff),
										removal_rates		= (if(is.null(removal_rate)) numeric(0) else removal_rate),
										lambda0				= (if(is.null(lambda0)) 0 else lambda0),
										splines_degree		= splines_degree,
										ODE_relative_dt		= ODE_relative_dt,
										ODE_relative_dy		= ODE_relative_dy,
										runtime_out_seconds	= -1)

	if(!results$success) return(list(success=FALSE, error=results$error))
	
	return(list(success			= TRUE,
				valid			= results$valid,
				ages			= age_grid,
				lambda			= results$lambdas,
				mu				= results$mus,
				psi				= results$psis,
				lambda_psi		= results$lambda_psis,
				Pmissing		= results$Pmissings,
				CSA_probs		= results$CSA_probs,
				CSA_Pmissings	= results$CSA_Pmissings,
				Reff			= results$Reffs,
				removal_rate	= results$removal_rates))
}


generate_OU_time_series = function(	times,					# numeric vector of size NT
									start_value,			# optional numeric. If NA or NaN or NULL, it will be chosen randomly from the stationary distribution
									stationary_mean,		# numeric
									stationary_std,			# non-negative numeric
									decay_rate,				# strictly positive numeric, in units 1/time
									constrain_min=-Inf,		# optional lower bound for the returned values
									constrain_max=+Inf){	# optional upper bound for the returned values
	values = get_Ornstein_Uhlenbeck_time_series_CPP(times			= times,
													start_value		= (if((!is.null(start_value)) && is.finite(start_value)) start_value else NaN),
													stationary_mean = stationary_mean,
													stationary_std	= stationary_std,
													decay_rate		= decay_rate)$values
	return(list(times=times, values=pmin(constrain_max,pmax(constrain_min,values))))
}




# fit a homogenous birth-death model on a grid to a given extant timetree, choosing the "best" grid size according to AIC or BIC
fit_hbd_model_on_best_grid_size = function(	tree, 
											oldest_age			= NULL,		# numeric, the oldest age to consider in the evaluation of the likelihood as well as for defining the age grid. Typically set to the stem age or root age. Can be NULL (equivalent to the root age).
											age0				= 0,		# non-negative numeric, youngest age (time before present) to consider when fitting and with respect to which rho is defined (rho(age0) is the fraction of lineages extant at age0 that are included in the tree)
											grid_sizes			= c(1,10),	# integer vector, listing the grid sizes to consider
											uniform_grid		= TRUE,		# logical, specifying whether the age grid should be uniform (equidistant age intervals). If FALSE, then the grid point density is chosen proportional to the square root of the LTT, hence resulting in higher resolution grid near the present.
											criterion			= "AIC",	# character, how to choose the optimal grid point. Options are "AIC" or "BIC".
											exhaustive			= TRUE,		# logical, whether to try all grid sizes for choosing the "best" one. If FALSE, the grid size is gradually increased until the selectin criterio (e.g., AIC) starts becoming worse, at which point the search is halted.
											min_lambda			= 0,		# numeric, lower bound for the fitted lambdas (applying to all grid points).
											max_lambda			= +Inf,		# numeric, upper bound for the fitted lambdas (applying to all grid points).
											min_mu				= 0,		# numeric, lower bound for the fitted mus (applying to all grid points).
											max_mu				= +Inf,		# numeric, upper bound for the fitted mus (applying to all grid points).
											min_rho0			= 1e-10,	# numeric, lower bound for the fitted rho. Note that rho is always within (0,1]
											max_rho0			= 1,		# numeric, upper bound for the fitted rho.
											guess_lambda		= NULL,		# initial guess for the lambda. Either NULL (an initial guess will be computed automatically), or a single numeric (guessing a constant lambda at all ages).
											guess_mu			= NULL,		# initial guess for the mu. Either NULL (an initial guess will be computed automatically), or a single numeric (guessing a constant mu at all ages).
											guess_rho0			= 1,		# initial guess for rho. Either NULL (an initial guess will be computed automatically) or a single strictly-positive numeric.
											fixed_lambda		= NULL,		# optional fixed lambda value. Either NULL (none of the lambdas are fixed), or a single scalar (all lambdas are fixed).
											fixed_mu			= NULL,		# optional fixed mu value. Either NULL (none of the mus are fixed), or a single scalar (all mus are fixed).
											fixed_rho0			= NULL,		# optional fixed value for rho. If non-NULL and non-NA, then rho is not fitted. 
											const_lambda		= FALSE,	# logical, whether to enforce a constant (time-independent) fitted speciation rate. Only relevant if lambdas are non-fixed.
											const_mu			= FALSE,	# logical, whether to enforce a constant (time-independent) fitted extinction rate. Only relevant if mus are non-fixed.
											splines_degree		= 1,		# integer, either 1 or 2 or 3, specifying the degree for the splines defined by lambda and mu on the age grid.
											condition			= "auto",	# one of "crown" or "stem" or "none" or "auto", specifying whether to condition the likelihood on the survival of the stem group or the crown group. It is recommended to use "stem" when oldest_age!=root_age, and "crown" when oldest_age==root_age. This argument is similar to the "cond" argument in the R function RPANDA::likelihood_bd. Note that "crown" really only makes sense when oldest_age==root_age.
											relative_dt			= 1e-3,		# maximum relative time step allowed for integration. Smaller values increase the accuracy of the computed likelihoods, but increase computation time. Typical values are 0.0001-0.001. The default is usually sufficient.
											Ntrials				= 1,
											Nthreads			= 1,
											max_model_runtime	= NULL,		# maximum time (in seconds) to allocate for each likelihood evaluation. Use this to escape from badly parameterized models during fitting (this will likely cause the affected fitting trial to fail). If NULL or <=0, this option is ignored.
											fit_control			= list(),
											verbose				= FALSE,
											verbose_prefix		= ""){
	# basic error checking
	if(verbose) cat(sprintf("%sChecking input parameters..\n",verbose_prefix))
	if((!is.null(guess_lambda))	&& (length(guess_lambda)!=1)) return(list(success=FALSE, error="Expecting either exactly one guess_lambda, or NULL"))									
	if((!is.null(guess_mu))	&& (length(guess_mu)!=1)) return(list(success=FALSE, error="Expecting either exactly one guess_mu, or NULL"))									
	if((!is.null(fixed_lambda))	&& (length(fixed_lambda)!=1)) return(list(success=FALSE, error="Expecting either exactly one fixed_lambda, or NULL"))									
	if((!is.null(fixed_mu))	&& (length(fixed_mu)!=1)) return(list(success=FALSE, error="Expecting either exactly one fixed_mu, or NULL"))									
	if(length(min_lambda)!=1) return(list(success=FALSE, error=sprintf("Expecting exactly one min_lambda; instead, received %d",length(min_lambda))))
	if(length(max_lambda)!=1) return(list(success=FALSE, error=sprintf("Expecting exactly one max_lambda; instead, received %d",length(max_lambda))))
	if(!(criterion %in% c("AIC", "BIC"))) return(list(success=FALSE, error=sprintf("Invalid model selection criterion '%s'. Expected 'AIC' or 'BIC'",criterion)))
	root_age = get_tree_span(tree)$max_distance
	Nmodels  = length(grid_sizes)
	
	# calculate tree LTT if needed
	if(!uniform_grid){
		LTT = count_lineages_through_time(tree=tree, Ntimes = max(100,10*max(grid_sizes)), regular_grid = TRUE)
		LTT$ages = root_age - LTT$times
	}
	
	# determine order in which to examine models
	if(exhaustive){
		model_order = c(1:Nmodels)
	}else{
		# examine models in the order of increasing grid sizes
		model_order = order(grid_sizes)
	}
	
	# fit HBD model on various grid sizes, keeping track of the "best" Ngrid
	if(verbose) cat(sprintf("%sFitting models with %s%d different grid sizes..\n",verbose_prefix,(if(exhaustive) "" else "up to "),Nmodels))
	AICs 		= rep(NA, times=Nmodels)
	BICs 		= rep(NA, times=Nmodels)
	best_fit	= NULL
	for(m in model_order){
		Ngrid = grid_sizes[m]
		if(uniform_grid){
			age_grid = seq(from=age0, to=oldest_age, length.out=Ngrid)
		}else{
			age_grid = get_inhomogeneous_grid_1D(Xstart = age0, Xend = oldest_age, Ngrid = Ngrid, densityX = rev(LTT$ages), densityY=sqrt(rev(LTT$lineages)), extrapolate=TRUE)
		}
		if(verbose) cat(sprintf("%s  Fitting model with grid size %d..\n",verbose_prefix,Ngrid))
		fit = fit_hbd_model_on_grid(tree				= tree, 
									oldest_age			= oldest_age,
									age0				= age0,
									age_grid			= age_grid,
									min_lambda			= min_lambda,
									max_lambda			= max_lambda,
									min_mu				= min_mu,
									max_mu				= max_mu,
									min_rho0			= min_rho0,
									max_rho0			= max_rho0,
									guess_lambda		= guess_lambda,
									guess_mu			= guess_mu,
									guess_rho0			= guess_rho0,
									fixed_lambda		= fixed_lambda,
									fixed_mu			= fixed_mu,
									fixed_rho0			= fixed_rho0,
									const_lambda		= const_lambda,
									const_mu			= const_mu,
									splines_degree		= splines_degree,
									condition			= condition,
									relative_dt			= relative_dt,
									Ntrials				= Ntrials,
									Nthreads			= Nthreads,
									max_model_runtime	= max_model_runtime,
									fit_control			= fit_control)
		if(!fit$success) return(list(success=FALSE, error=sprintf("Fitting model with grid size %d failed: %s",Ngrid,fit$error)))
		criterion_value = fit[[criterion]]
		if(is.null(best_fit)){
			best_fit = fit
			worsened = FALSE
		}else if(criterion_value<best_fit[[criterion]]){
			best_fit = fit
			worsened = FALSE
		}else{
			worsened = TRUE
		}
		AICs[m] = fit$AIC
		BICs[m] = fit$BIC
		if(verbose) cat(sprintf("%s  --> %s=%.10g. Best grid size so far: %d\n",verbose_prefix,criterion,criterion_value,length(best_fit$age_grid)))
		if((!exhaustive) && worsened) break; # model selection criterion became worse compared to the previous grid size, so stop search and keep best model found so far
	}
	
	return(list(success	 	= TRUE,
				best_fit 	= best_fit,
				grid_sizes	= grid_sizes,
				AICs		= AICs,
				BICs		= BICs))
}


# convert a list to a vector using unlist, after converting NULLs to NAs
# this function may be needed if your list contains NULLs, in which case the standard unlist removes those entries from the returned vector
unlist_with_nulls = function(L){
	L[sapply(L, is.null)] = NA
	return(unlist(L))
}


# given a list of objects, unpack every object that is a pure list (non recursively)
# "pure list" means that class(object) = c("list")
# For example, list(a=3, b=list(c=4,d=5), e=list(f=list(6,7), g=8, h=phylo_object, i=list())) will become list(a=3, c=4, d=5, f=list(6,7), g=8, h=phylo_object)
# This function may be needed when returning a large list of lists from an Rcpp function (which may necessitate list nesting due to limitations of Rcpp) and you want to unpack the first-level lists.
flatten_list_first_level = function(L){
	N  = length(L)
	NF = sum(sapply(seq_len(N), FUN=function(i) (if(all(class(L[[i]])=="list")) length(L[[i]]) else 1))) # predicted length of the flattened list
	LF = vector(mode="list", NF)
	Lnames = names(L)
	LFnames = character(NF)
	f  = 1 # next index in flattened list at which to place item
	for(i in seq_len(N)){
		if(all(class(L[[i]])=="list")){
			if(length(L[[i]])>0){
				LF[f:(f+length(L[[i]])-1)] = L[[i]]
				if((!is.null(Lnames)) && (!is.null(names(L[[i]])))) LFnames[f:(f+length(L[[i]])-1)] = names(L[[i]])
			}
			f = f + length(L[[i]])
		}else{
			if(!is.null(Lnames)) LFnames[f] = Lnames[i]
			LF[[f]] = L[[i]]
			f = f + 1
		}
	}
	names(LF) = LFnames
	return(LF)
}



# given a list of birth events and death events, calculate the corresponding LTT
# Note that birth & death events need not be consistent with a tree topology, i.e. the LTT simply decreases (over increasing time) with every death event and increases with every birth event
count_lineages_through_time_BD = function(	birth_times, 	# numeric vector of arbitrary length, listing times of birth events. Need not necessarily be sorted.
											death_times, 	# numeric vector of arbitrary length, listing times of death events. Need not necessarily be sorted.
											time_grid, 		# times at which lineage counts are requested. Need not necessarily be sorted.
											initial_count){ # integer, the number of lineages assumed at -infinity
	if(length(time_grid)==0) return(numeric())
	if(any(diff(time_grid)<0)){
		# time grid is not sorted in ascending order, so figure out proper order
		time_grid_order = order(time_grid)
	}else{
		time_grid_order = c(1:length(time_grid))
	}
	lineages = numeric(length(time_grid))
	lineages[time_grid_order] = initial_count + get_LTT_BD_CPP(birth_times = sort(birth_times), death_times = sort(death_times), time_grid=time_grid[time_grid_order]);
}



# extract independent contrasts from a bifurcating tree, for fitting SBM diffusion models
get_SBM_independent_contrasts = function(	tree,
											tip_latitudes, 						# numeric vector of size Ntips, listing geographical latitudes of the tips (in decimal degrees)
											tip_longitudes, 					# numeric vector of size Ntips, listing geographical longitudes of the tips (in decimal degrees)
											radius,								# numeric, radius to assume for the sphere (e.g. Earth). Use this e.g. if you want to hange the units in which diffusivity is estimated. Earth's mean radius is about 6371e3 m.
											clade_states			= NULL,		# optional, either an integer vector of length Ntips+Nnodes (if trees[] is a single tree) or a list of 1D vectors (if trees[] is a list of trees), specifying the discrete "state" of each tip and node in each tree. This can be used to limit independent contrasts to tip pairs whose total number of state-transitions (along their shortest path) is zero.
											planar_approximation	= FALSE,	# logical, specifying whether the estimation formula should be based on a planar approximation of Earth's surface, i.e. geodesic angles are converted to distances and then those are treated as if they were Euclideanon a 2D plane. This approximation substantially increases the speed of computations.
											only_basal_tip_pairs	= FALSE,	# logical, specifying whether only immediate sister tips should be considered, i.e. tip pairs with at most 2 edges between the two tips
											only_distant_tip_pairs	= FALSE,	# logical, whether to only consider tip pairs located at distinct geographic locations
											min_MRCA_time			= 0,		# numeric, specifying the minimum allowed height (distance from root) of the MRCA of sister tips considered in the fitting. In other words, an independent contrast is only considered if the two sister tips' MRCA has at least this distance from the root. Set min_MRCA_time=0 to disable this filter.
											max_MRCA_age			= Inf,		# numeric, specifying the maximum allowed age (distance from youngest tip) of the MRCA of sister tips considered in the fitting. In other words, an independent contrast is only considered if the two sister tips' MRCA has at most this age (time to present). Set max_MRCA_age=Inf to disable this filter.
											no_state_transitions	= FALSE,	# if TRUE, only tip pairs without state transitions along their shortest paths are considered. In particular, only tips in the same state are considered. Requires that clade_states[] is provided.
											only_state				= NULL){	# optional integer, specifying the state in which tip pairs (and their connecting ancestors) must be in order to be considered. Requires that clade_states[] is provided.
	Ntips = length(tree$tip.label)
	if(("list" %in% class(tip_latitudes)) && (length(tip_latitudes)==Ntips)){
		tip_latitudes = unlist(tip_latitudes)
	}
	if(("list" %in% class(tip_latitudes)) && (length(tip_longitudes)==Ntips)){
		tip_longitudes = unlist(tip_longitudes)
	}
	if((!is.null(clade_states)) && ("list" %in% class(clade_states)) && (length(clade_states)==Ntips+tree$Nnode)){
		clade_states = unlist(clade_states)
	}
	if((!is.null(only_state)) && is.null(clade_states)) return(list(success=FALSE, error="Missing clade_states[], needed when only_state is specified"))
	if(no_state_transitions && is.null(clade_states)) return(list(success=FALSE, error="Missing clade_states[], needed when no_state_transitions=TRUE"))
		
	# make sure tree does not have multifurcations
	# note that this does not change the tip & node indices, so the returned tip_pairs[] will still refer to the correct tips in the original tree
	if(tree_has_multifurcations_CPP(Ntips = Ntips, Nnodes = tree$Nnode, Nedges = nrow(tree$edge), tree_edge = as.vector(t(tree$edge)) - 1)){
		tree = multifurcations_to_bifurcations(tree)$tree
	}

	# extract independent pairs of sister tips
	tip_pairs = extract_independent_sister_tips(tree)
	if(only_basal_tip_pairs){
		# calculate number of nodes between tip pairs
		edge_counts = get_pairwise_distances(tree, A=tip_pairs[,1], B=tip_pairs[,2], as_edge_counts=TRUE, check_input=FALSE)
		# only keep tip pairs with at most 2 edges connecting them
		keep_pairs 	= which(edge_counts<=2)
		tip_pairs 	= tip_pairs[keep_pairs,,drop=FALSE]
	}
	
	# calculate MRCAs
	MRCAs = get_pairwise_mrcas(tree, tip_pairs[,1], tip_pairs[,2], check_input=FALSE)
	
	# calculate clade times
	clade_times = castor::get_all_distances_to_root(tree)
	
	# filter tip pairs if needed
	if((min_MRCA_time>0) || (max_MRCA_age<Inf)){
		tree_span 	= max(clade_times)
		keep_pairs	= which((clade_times[MRCAs]>=min_MRCA_time) & (tree_span-clade_times[MRCAs]<=max_MRCA_age))
		tip_pairs	= tip_pairs[keep_pairs,,drop=FALSE]
		MRCAs		= MRCAs[keep_pairs]
	}
	if(nrow(tip_pairs)==0) return(list(success=FALSE, error="No valid tip pairs left for extracting independent contrasts"))
		
	# calculate phylogenetic divergences and geodesic distances between sister tips
	phylodistances 	= get_pairwise_distances(tree, A=tip_pairs[,1], B=tip_pairs[,2], check_input=FALSE)
	geodistances 	= radius * sapply(1:nrow(tip_pairs), FUN=function(p) geodesic_angle(tip_latitudes[tip_pairs[p,1]],tip_longitudes[tip_pairs[p,1]],tip_latitudes[tip_pairs[p,2]],tip_longitudes[tip_pairs[p,2]]))

	# omit tip pairs with zero phylogenetic distance, because in that case the likelihood density is pathological
	# also omit tip pairs located at the same geographic location, if requested
	keep_pair	 	= (phylodistances>0)
	if(only_distant_tip_pairs) keep_pair = keep_pair & (geodistances>0)
	if(no_state_transitions || (!is.null(only_state))){
		Ntransitions = count_transitions_between_clades(tree=tree, A=tip_pairs[,1], B=tip_pairs[,2], states=clade_states, check_input=TRUE)
		if(no_state_transitions) keep_pair = keep_pair & (Ntransitions==0)
		if(!is.null(only_state)) keep_pair = keep_pair & (Ntransitions==0) & (clade_states[tip_pairs[,1]]==only_state) & (clade_states[tip_pairs[,2]]==only_state)
	}
	tip_pairs		= tip_pairs[keep_pair,,drop=FALSE]
	MRCAs			= MRCAs[keep_pair]
	phylodistances 	= phylodistances[keep_pair]
	geodistances 	= geodistances[keep_pair]
	NC 				= length(phylodistances)
	if(NC==0) return(list(success=FALSE, error="No valid tip pairs left for extracting independent contrasts"))

	# determine MRCA and tip times
	MRCA_times = clade_times[MRCAs]
	tip_times1 = clade_times[tip_pairs[,1]]
	tip_times2 = clade_times[tip_pairs[,2]]

	return(list(success			= TRUE,
				NC				= NC,
				tip_pairs		= tip_pairs,		# numeric matrix of size NC x 2, listing the tip indices used to define the indepent contrasts
				phylodistances	= phylodistances,
				geodistances	= geodistances,
				MRCA_times		= MRCA_times,
				child_times1	= tip_times1,
				child_times2	= tip_times2))
}



# fit a homogenous birth-death model on a grid to a given extant timetree, choosing the "best" grid size according to AIC or BIC
fit_sbm_on_best_grid_size = function(	tree, 
										tip_latitudes, 						# numeric vector of size Ntips, listing geographical latitudes of the tips (in decimal degrees)
										tip_longitudes, 					# numeric vector of size Ntips, listing geographical longitudes of the tips (in decimal degrees)
										radius,								# numeric, radius to assume for the sphere (e.g. Earth). Use this e.g. if you want to hange the units in which diffusivity is estimated. Earth's mean radius is about 6371e3 m.
										clade_states			= NULL,		# optional, either an integer vector of length Ntips+Nnodes (if trees[] is a single tree) or a list of 1D vectors (if trees[] is a list of trees), specifying the discrete "state" of each tip and node in each tree. This can be used to limit independent contrasts to tip pairs whose total number of state-transitions (along their shortest path) is zero.
										planar_approximation	= FALSE,	# logical, specifying whether the estimation formula should be based on a planar approximation of Earth's surface, i.e. geodesic angles are converted to distances and then those are treated as if they were Euclideanon a 2D plane. This approximation substantially increases the speed of computations.
										only_basal_tip_pairs	= FALSE,	# logical, specifying whether only immediate sister tips should be considered, i.e. tip pairs with at most 2 edges between the two tips
										only_distant_tip_pairs	= FALSE,	# logical, whether to only consider tip pairs located at distinct geographic locations
										min_MRCA_time			= 0,		# numeric, specifying the minimum allowed height (distance from root) of the MRCA of sister tips considered in the fitting. In other words, an independent contrast is only considered if the two sister tips' MRCA has at least this distance from the root. Set min_MRCA_time=0 to disable this filter.
										max_MRCA_age			= Inf,		# numeric, specifying the maximum allowed age (distance from youngest tip) of the MRCA of sister tips considered in the fitting. In other words, an independent contrast is only considered if the two sister tips' MRCA has at most this age (time to present). Set max_MRCA_age=Inf to disable this filter.
										no_state_transitions	= FALSE,	# if TRUE, only tip pairs without state transitions along their shortest paths are considered. In particular, only tips in the same state are considered. Requires that clade_states[] is provided.
										only_state				= NULL,		# optional integer, specifying the state in which tip pairs (and their connecting ancestors) must be in order to be considered. Requires that clade_states[] is provided.
										grid_sizes				= c(1,10),	# integer vector, listing the grid sizes to consider
										uniform_grid			= TRUE,		# logical, specifying whether the age grid should be uniform (equidistant age intervals). If FALSE, then the grid point density is chosen proportional to the square root of the independent contrasts density, hence resulting in higher resolution grid in areas where there is more data available.
										guess_diffusivity		= NULL,		# optional numeric, first guess for the diffusivity (at all grid points). If NULL, this will be automatically chosen.
										min_diffusivity			= NULL,		# numeric, lower bound for the fitted diffusivity (applying to all grid points). If NULL, this is automatically chosen.
										max_diffusivity			= Inf,		# numeric, upper bound for the fitted diffusivity (applying to all grid points).
										criterion				= "AIC",	# character, how to choose the optimal grid point. Options are "AIC" or "BIC".
										exhaustive				= TRUE,		# logical, whether to try all grid sizes for choosing the "best" one. If FALSE, the grid size is gradually increased until the selectin criterio (e.g., AIC) starts becoming worse, at which point the search is halted.
										Ntrials					= 1,
										Nthreads				= 1,
										Nbootstraps				= 0,		# (integer) optional number of parametric-bootstrap samples for estimating confidence intervals of fitted parameters. If 0, no parametric bootstrapping is performed. Typical values are 10-100.
										Ntrials_per_bootstrap	= NULL,		# (integer) optional number of fitting trials for each bootstrap sampling. If NULL, this is set equal to Ntrials. A smaller Ntrials_per_bootstrap will reduce computation, at the expense of increasing the estimated confidence intervals (i.e. yielding more conservative estimates of confidence).
										NQQ						= 0,		# (integer) optional number of simulations to perform for creating Q-Q plots of the theoretically expected distribution of geodistances vs the empirical distribution of geodistances (across independent contrasts). The resolution of the returned QQ plot will be equal to the number of independent contrasts used for fitting.
										fit_control				= list(),
										SBM_PD_functor			= NULL,		# internally used SBM probability density functor
										verbose					= FALSE,
										verbose_prefix			= ""){
	# basic error checking
	if(verbose) cat(sprintf("%sChecking input parameters..\n",verbose_prefix))
	if((!is.null(guess_diffusivity)) && (length(guess_diffusivity)!=1)) return(list(success=FALSE, error="Expecting either exactly one guess_diffusivity, or NULL"))									
	if((!is.null(min_diffusivity)) && (length(min_diffusivity)!=1)) return(list(success=FALSE, error="Expecting either exactly one min_diffusivity, or NULL"))									
	if((!is.null(max_diffusivity)) && (length(max_diffusivity)!=1)) return(list(success=FALSE, error="Expecting either exactly one min_diffusivity, or NULL"))									
	if(!(criterion %in% c("AIC", "BIC"))) return(list(success=FALSE, error=sprintf("Invalid model selection criterion '%s'. Expected 'AIC' or 'BIC'",criterion)))
	root_age = get_tree_span(tree)$max_distance
	Nmodels  = length(grid_sizes)
		
	# determine order in which to examine models
	if(exhaustive){
		model_order = c(1:Nmodels)
	}else{
		# examine models in the order of increasing grid sizes
		model_order = order(grid_sizes)
	}
	
	if(!uniform_grid){
		# get independent contrasts density through time, needed for defining non-uniform grid
		if(verbose) cat(sprintf("%sDetermining independent contrasts density..\n",verbose_prefix))
		ICs = get_SBM_independent_contrasts(tree					= tree,
											tip_latitudes			= tip_latitudes,
											tip_longitudes			= tip_longitudes,
											radius					= radius,
											clade_states			= clade_states,
											planar_approximation	= planar_approximation,
											only_basal_tip_pairs	= only_basal_tip_pairs,
											only_distant_tip_pairs	= only_distant_tip_pairs,
											min_MRCA_time			= min_MRCA_time,
											max_MRCA_age			= max_MRCA_age,
											no_state_transitions	= no_state_transitions,
											only_state				= only_state)
		ICdensity = list(times = seq(from=root_age/1000, to=root_age, length.out=200))
		ICdensity$lineages = count_lineages_through_time_BD(birth_times		= c(ICs$MRCA_times,ICs$MRCA_times),
															death_times		= c(ICs$child_times1,ICs$child_times2),
															time_grid		= ICdensity$times, 
															initial_count	= 0)
		if(!ICs$success) return(list(success=FALSE, error=sprintf("Failed to acquire independent contrasts: %s",ICs$error)))
	}
	
	# fit SBM model on various grid sizes, keeping track of the "best" Ngrid
	if(verbose) cat(sprintf("%sFitting SBM models with %s%d different grid sizes..\n",verbose_prefix,(if(exhaustive) "" else "up to "),Nmodels))
	AICs 		= rep(NA, times=Nmodels)
	BICs 		= rep(NA, times=Nmodels)
	best_fit	= NULL
	for(m in model_order){
		Ngrid = grid_sizes[m]
		if(uniform_grid || (Ngrid==1)){
			time_grid = seq(from=0, to=root_age, length.out=Ngrid)
		}else{
			time_grid = get_inhomogeneous_grid_1D(Xstart = 0, Xend = root_age, Ngrid = Ngrid, densityX = ICdensity$times, densityY=sqrt(ICdensity$lineages), extrapolate=TRUE)
		}
		if(verbose) cat(sprintf("%s  Fitting SBM model with grid size %d..\n",verbose_prefix,Ngrid))
		fit = fit_sbm_on_grid(	tree					= tree, 
								tip_latitudes			= tip_latitudes,
								tip_longitudes			= tip_longitudes,
								radius					= radius,
								clade_states			= clade_states,
								planar_approximation	= planar_approximation,
								only_basal_tip_pairs	= only_basal_tip_pairs,
								only_distant_tip_pairs	= only_distant_tip_pairs,
								min_MRCA_time			= min_MRCA_time,
								max_MRCA_age			= max_MRCA_age,
								no_state_transitions	= no_state_transitions,
								only_state				= only_state,
								time_grid				= time_grid,
								guess_diffusivity		= guess_diffusivity,
								min_diffusivity			= min_diffusivity,
								max_diffusivity			= max_diffusivity,
								Ntrials					= Ntrials,
								Nthreads				= Nthreads,
								Nbootstraps				= 0,
								Ntrials_per_bootstrap	= Ntrials_per_bootstrap,
								NQQ						= NQQ,
								fit_control				= fit_control,
								SBM_PD_functor			= (if(!is.null(best_fit)) best_fit$SBM_PD_functor else SBM_PD_functor),
								verbose					= verbose,
								verbose_prefix			= paste0(verbose_prefix,"    "))
		if(!fit$success) return(list(success=FALSE, error=sprintf("Fitting SBM model with grid size %d failed: %s",Ngrid,fit$error)))
		criterion_value = fit[[criterion]]
		if(is.null(best_fit)){
			best_fit = fit
			worsened = FALSE
		}else if(criterion_value<best_fit[[criterion]]){
			best_fit = fit
			worsened = FALSE
		}else{
			worsened = TRUE
		}
		AICs[m] = fit$AIC
		BICs[m] = fit$BIC
		if(verbose) cat(sprintf("%s  --> %s=%.10g. Best grid size so far: %d\n",verbose_prefix,criterion,criterion_value,length(best_fit$time_grid)))
		if((!exhaustive) && worsened) break; # model selection criterion became worse compared to the previous grid size, so stop search and keep best model found so far
	}
	
	if(Nbootstraps>0){
		if(verbose) cat(sprintf("%sRepeating fit of model with grid size %d, for bootstrapping..\n",verbose_prefix,length(best_fit$time_grid)))
		best_fit = fit_sbm_on_grid(	tree					= tree, 
									tip_latitudes			= tip_latitudes,
									tip_longitudes			= tip_longitudes,
									radius					= radius,
									clade_states			= clade_states,
									planar_approximation	= planar_approximation,
									only_basal_tip_pairs	= only_basal_tip_pairs,
									only_distant_tip_pairs	= only_distant_tip_pairs,
									min_MRCA_time			= min_MRCA_time,
									max_MRCA_age			= max_MRCA_age,
									no_state_transitions	= no_state_transitions,
									only_state				= only_state,
									time_grid				= best_fit$time_grid,
									guess_diffusivity		= guess_diffusivity,
									min_diffusivity			= min_diffusivity,
									max_diffusivity			= max_diffusivity,
									Ntrials					= Ntrials,
									Nthreads				= Nthreads,
									Nbootstraps				= Nbootstraps,
									Ntrials_per_bootstrap	= Ntrials_per_bootstrap,
									NQQ						= NQQ,
									fit_control				= fit_control,
									SBM_PD_functor			= best_fit$SBM_PD_functor,
									verbose					= verbose,
									verbose_prefix			= paste0(verbose_prefix,"  "))

	}
	
	return(list(success	 	= TRUE,
				best_fit 	= best_fit,
				grid_sizes	= grid_sizes,
				AICs		= AICs,
				BICs		= BICs))
}


read_fasta = function(	file,
						include_headers		= TRUE,
						include_sequences	= TRUE,
						truncate_headers_at	= NULL){ # optional needle string, at which to truncate headers (i.e. remove everything at and after the first instance of the needle)
	results = read_fasta_from_file_CPP(	fasta_path			= file,
										include_headers		= include_headers,
										include_sequences	= include_sequences)
	if(!results$success) return(list("success"=FALSE, error=results$error))
	if(include_headers && (!is.null(truncate_headers_at))){
		results$headers = sapply(seq_len(length(results$headers)), FUN=function(h) strsplit(results$headers[h],split=truncate_headers_at,fixed=TRUE)[[1]][1])
	}
	return(list(headers		= (if(include_headers) results$headers else NULL),
				sequences	= (if(include_sequences) results$sequences else NULL),
				Nlines		= results$Nlines,
				Nsequences	= results$Nsequences))
}



monotonize_time_series = function(	values,				# 1D array of size N, listing the scalar time series values. May include NaN.
									increasing, 		# logical, specifying whether the resulting time series should be monotonically increasing (rather than decreasing)
									prefer_later_data){	# logical, specifying whether later data (rather than earlier data) should be kept when resolving monotonicity conflicts
	results = monotonize_time_series_CPP(values = values, increasing = increasing, prefer_later_data = prefer_later_data)
	return(results)
}

