# count the number of clades that existed on each of multiple equidistant time points, or at a specific set of time points
# time = distance from root
# Ntimes = number of time points to consider, spanning 0 to the maximum time
# if tree$edge.length is missing, edges are assumed to have length 1.
count_clades_over_time = function(	tree, 
									Ntimes			= NULL, 		# number of equidistant time points for which to calculate diversities
									times 			= NULL, 	# 1D array of time points in increasing order, for which to calculate diversities
									include_slopes	= FALSE){
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode;
	if((!is.null(Ntimes)) && (!is.null(times))) stop("ERROR: Either Ntimes or times must be non-NULL, but not both")
	if(is.null(Ntimes) && is.null(times)) stop("ERROR: Both Ntimes and times are NULL; please specify one of the two")
	
	if(is.null(times)){
		results = count_clades_at_regular_times_CPP(Ntips			= Ntips,
													Nnodes			= Nnodes,
													Nedges			= nrow(tree$edge),
													tree_edge		= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
													edge_length		= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
													Ntimes			= Ntimes,
													include_slopes 	= include_slopes);
		return(list(times			= results$time_points, 
					diversities		= results$diversities, 
					slopes			= (if(include_slopes) results$slopes else NULL),
					relative_slopes	= (if(include_slopes) results$relative_slopes else NULL)));
	}else{
		Ntimes = length(times)
		diversities = count_clades_at_times_CPP(	Ntips 		= Ntips,
													Nnodes 		= Nnodes,
													Nedges 		= nrow(tree$edge),
													tree_edge	= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
													edge_length	= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
													times		= times);
		if(include_slopes){
			slopes = c((diversities[2]-diversities[1])/(times[2]-times[1]), (diversities[3:Ntimes]-diversities[1:(Ntimes-2)])/(times[3:Ntimes]-times[1:(Ntimes-2)]), (diversities[Ntimes]-diversities[Ntimes-1])/(times[Ntimes]-times[Ntimes-1]))
			CC = c(0.5*(diversities[2]+diversities[1]), (1/3.0)*(diversities[1:(Ntimes-2)]+diversities[2:(Ntimes-1)]+diversities[3:Ntimes]), 0.5*(diversities[Ntimes]+diversities[Ntimes-1]))
			relative_slopes = slopes/CC;
		}else{
			slopes = NULL
			relative_slopes = NULL
		}
		return(list(Ntimes			= length(times),
					times			= times,
					diversities 	= diversities,
					slopes			= slopes,
					relative_slopes	= relative_slopes));
	}
}