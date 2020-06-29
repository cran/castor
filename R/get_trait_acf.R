get_trait_acf = function(	tree, 
							tip_states,
							Npairs=10000,
							Nbins=10){
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode;
	Nedges = nrow(tree$edge);
	if(!is.numeric(tip_states)) stop("ERROR: tip_states must be numeric")
	results = ACF_continuous_trait_CPP(	Ntips 				= Ntips,
										Nnodes 				= Nnodes,
										Nedges 				= Nedges,
										tree_edge			= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
										edge_length		 	= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
										state_per_tip 		= tip_states,
										Npairs				= Npairs,
										Nbins				= Nbins,
										verbose 			= FALSE,
										verbose_prefix 		= "");
	return(list(distances=results$distance_grid,
				autocorrelations=results$autocorrelations,
				mean_abs_differences=results$mean_abs_deviations,
				mean_rel_differences=results$mean_rel_deviations,
				Npairs_per_distance=results$N_per_grid_point));
}