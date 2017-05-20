# Perform random simulation of a Brownian Motion model of continuous trait evolution, moving from root to tips.
# The root's state is explicitly specified at each simulation.
# Optionally, multiple independent simulations can be performed using the same model (e.g. as part of some Monte Carlo integration)
simulate_bm_model = function(	tree, 
								diffusivity,				# diffusion rate of the model (in units trait_units^2/edge_length)
								root_states 	= 0, 		# numeric vector of arbitrary size, specifying states for the root. If smaller than Nsimulations, then values are recycled. If empty, zero is used as root state.
								include_tips	= TRUE, 
								include_nodes	= TRUE, 
								Nsimulations	= 1,
								drop			= TRUE){
	Ntips  	 	= length(tree$tip.label);
	Nnodes  	= tree$Nnode;
	results = simulate_Brownian_motion_model_CPP(	Ntips				= Ntips,
													Nnodes				= Nnodes,
													Nedges				= nrow(tree$edge),
													tree_edge 			= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based,
													edge_length		 	= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
													root_states			= root_states,
													diffusivity			= diffusivity,
													include_tips		= include_tips,
													include_nodes		= include_nodes,
													Nsimulations		= Nsimulations);

	tip_states  = NULL
	node_states = NULL
	if(include_tips) tip_states = (if(drop && Nsimulations==1) results$tip_states else matrix(results$tip_states, ncol=Ntips, byrow=TRUE));
	if(include_nodes) node_states = (if(drop && Nsimulations==1) results$node_states else matrix(results$node_states, ncol=Nnodes, byrow=TRUE));
	return(list(tip_states=tip_states, node_states=node_states));
}