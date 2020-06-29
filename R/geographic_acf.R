geographic_acf = function(	tree, 
							tip_latitudes,
							tip_longitudes,
							Npairs	= 10000,
							Nbins	= 10){
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode;
	Nedges = nrow(tree$edge);
	
	# basic option checking
	if(!is.numeric(tip_latitudes)) return(list(success=FALSE, error="tip_latitudes must be numeric"))
	if(!is.numeric(tip_longitudes)) return(list(success=FALSE, error="tip_longitudes must be numeric"))
	if(length(tip_latitudes)!=Ntips) return(list(success=FALSE, error=sprintf("tip_latitudes has length %d, but expected exactly %d (=Ntips) entries",length(tip_latitudes),Ntips)))
	if(length(tip_longitudes)!=Ntips) return(list(success=FALSE, error=sprintf("tip_longitudes has length %d, but expected exactly %d (=Ntips) entries",length(tip_longitudes),Ntips)))
	
	results = ACF_spherical_CPP(	Ntips 				= Ntips,
									Nnodes 				= Nnodes,
									Nedges 				= Nedges,
									tree_edge			= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
									edge_length		 	= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
									tip_latitudes 		= tip_latitudes,
									tip_longitudes 		= tip_longitudes,
									Npairs				= Npairs,
									Nbins				= Nbins,
									verbose 			= FALSE,
									verbose_prefix 		= "");
	return(list(success 			= TRUE,
				distances			= results$distance_grid,
				autocorrelations	= results$autocorrelations,
				mean_geodistances	= results$mean_geodistances,
				Npairs_per_distance	= results$N_per_grid_point));
}