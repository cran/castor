# Perform random simulation of a Brownian Motion model of continuous multivariate trait evolution, moving from root to tips.
# The root's state is explicitly specified at each simulation.
# Optionally, multiple independent simulations can be performed using the same model (e.g. as part of some Monte Carlo integration)
# The diffusivity matrix D is a non-negative definite symmetric matrix, such that exp(-X^T*D^{-1}*X/(4*L))/sqrt(det(2*pi*D)) is the probability density for the multidimensional trait vector X after phylogenetic distance L, if initially located at the origin.
simulate_bm_model = function(	tree, 
								diffusivity,				# either a single number, or a 2D array of size Ntraits, diffusivity matrix of the model (in units trait_units^2/edge_length)
								root_states 	= NULL, 	# 2D numeric matrix of size NR x Ntrait (where NR can be arbitrary), specifying states for the root. If NR is smaller than Nsimulations, then values are recycled. If NULL, zero is used as root state for all traits.
								include_tips	= TRUE, 
								include_nodes	= TRUE, 
								Nsimulations	= 1,
								drop_dims	= TRUE){
	Ntips  	 	= length(tree$tip.label);
	Nnodes  	= tree$Nnode;
	Ntraits 	= (if(is.vector(diffusivity)) 1 else nrow(diffusivity))
	if((Ntraits>1) && (nrow(diffusivity)!=ncol(diffusivity))) stop(sprintf("ERROR: Diffusivity matrix must be quadratic (instead, it has dimensions %d x %d)",nrow(diffusivity),ncol(diffusivity)))
	if(is.null(root_states)) root_states = numeric()
	
	if(Ntraits==1){
		results = simulate_scalar_Brownian_motion_model_CPP(Ntips				= Ntips,
															Nnodes				= Nnodes,
															Nedges				= nrow(tree$edge),
															tree_edge 			= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based,
															edge_length		 	= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
															root_states			= root_states,
															diffusivity			= diffusivity,
															include_tips		= include_tips,
															include_nodes		= include_nodes,
															Nsimulations		= Nsimulations);	
	}else{
		cholesky = t(chol(diffusivity)); # lower-triangular part of Chlesky decomposition, i.e. such that diffusivity = cholesky * cholesky^T
		results = simulate_multivariate_Brownian_motion_model_CPP(	Ntips				= Ntips,
																	Nnodes				= Nnodes,
																	Nedges				= nrow(tree$edge),
																	Ntraits				= Ntraits,
																	tree_edge 			= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based,
																	edge_length		 	= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
																	root_states			= as.vector(t(root_states)),
																	diffusivity			= as.vector(t(diffusivity)), 	# flatten in row-major format
																	cholesky			= as.vector(t(cholesky)), 		# flatten in row-major format
																	include_tips		= include_tips,
																	include_nodes		= include_nodes,
																	Nsimulations		= Nsimulations);
	}

	# unflatten returned arrays
	tip_states  = NULL
	node_states = NULL
	if(include_tips){
		if(drop_dims && (Nsimulations==1) && (Ntraits==1)){ 
			tip_states = results$tip_states;
		}else if(drop_dims && (Ntraits==1)){
			tip_states = matrix(results$tip_states, ncol=Ntips, byrow=TRUE)
		}else if(drop_dims && (Nsimulations==1)){
			tip_states = matrix(results$tip_states, ncol=Ntraits, byrow=TRUE)
		}else{ 
			tip_states = aperm(array(results$tip_states,dim=c(Ntraits,Ntips,Nsimulations)),c(3,2,1))
		}
	}
	if(include_nodes){
		if(drop_dims && (Nsimulations==1) && (Ntraits==1)){ 
			node_states = results$node_states;
		}else if(drop_dims && (Ntraits==1)){
			node_states = matrix(results$node_states, ncol=Nnodes, byrow=TRUE)
		}else if(drop_dims && (Nsimulations==1)){
			node_states = matrix(results$node_states, ncol=Ntraits, byrow=TRUE)
		}else{ 
			node_states = aperm(array(results$node_states,dim=c(Ntraits,Nnodes,Nsimulations)),c(3,2,1))
		}
	}	
	return(list(tip_states	= tip_states, 	# 3D matrix of size Nsimulations x Ntips x Ntraits
				node_states	= node_states));# 3D matrix of size Nsimulations x Nnodes x Ntraits
}
