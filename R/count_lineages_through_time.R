# count the number of clades (lineages) in the tree at each of multiple equidistant time points, or at a specific set of time points
# This is the classical lineages-through-time (LTT) curve
# time = distance from root
# Ntimes = number of time points to consider, spanning 0 to the maximum time
# if tree$edge.length is missing, edges are assumed to have length 1.
count_lineages_through_time = function(	tree, 
										Ntimes			= NULL, 	# number of equidistant time points at which to calculate lineages
										min_time		= NULL,		# minimum time to consider. If NULL, will be set to the minimum possible
										max_time		= NULL,		# maximum time to consider. If NULL, will be set to the maximum possible
										times 			= NULL, 	# 1D array of time points in increasing order, for which to calculate lineages
										include_slopes	= FALSE){
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode;
	if((!is.null(Ntimes)) && (!is.null(times))) stop("ERROR: Either Ntimes or times must be non-NULL, but not both")
	if(is.null(Ntimes) && is.null(times)) stop("ERROR: Both Ntimes and times are NULL; please specify one of the two")
	if((!is.null(min_time)) && (!is.null(times))) stop("ERROR: min_time and times cannot both be non-NULL; choose one method to specify time points")
	if((!is.null(max_time)) && (!is.null(times))) stop("ERROR: max_time and times cannot both be non-NULL; choose one method to specify time points")
	
	if(is.null(times)){
		results = count_clades_at_regular_times_CPP(Ntips			= Ntips,
													Nnodes			= Nnodes,
													Nedges			= nrow(tree$edge),
													tree_edge		= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
													edge_length		= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
													Ntimes			= Ntimes,
													min_time		= (if(is.null(min_time)) 0 else min_time),
													max_time		= (if(is.null(max_time)) Inf else max_time),
													include_slopes 	= include_slopes);
		return(list(Ntimes			= length(results$time_points),
					times			= results$time_points, 
					lineages		= results$diversities, 
					slopes			= (if(include_slopes) results$slopes else NULL),
					relative_slopes	= (if(include_slopes) results$relative_slopes else NULL)));
		
	}else{
		Ntimes = length(times)
		lineages = count_clades_at_times_CPP(	Ntips 		= Ntips,
												Nnodes 		= Nnodes,
												Nedges 		= nrow(tree$edge),
												tree_edge	= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
												edge_length	= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
												times		= times);
		if(include_slopes){
			slopes = c((lineages[2]-lineages[1])/(times[2]-times[1]), (lineages[3:Ntimes]-lineages[1:(Ntimes-2)])/(times[3:Ntimes]-times[1:(Ntimes-2)]), (lineages[Ntimes]-lineages[Ntimes-1])/(times[Ntimes]-times[Ntimes-1]))
			CC = c(0.5*(lineages[2]+lineages[1]), (1/3.0)*(lineages[1:(Ntimes-2)]+lineages[2:(Ntimes-1)]+lineages[3:Ntimes]), 0.5*(lineages[Ntimes]+lineages[Ntimes-1]))
			relative_slopes = slopes/CC;
		}else{
			slopes = NULL
			relative_slopes = NULL
		}
		return(list(Ntimes			= length(times),
					times			= times,
					lineages 		= lineages,
					slopes			= slopes,
					relative_slopes	= relative_slopes));
	}
}