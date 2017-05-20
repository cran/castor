# count the number of clades that existed on each of multiple equidistant time points
# time = distance from root
# resolution = number of time points to consider, spanning 0 to the maximum time
# if tree$edge.length is missing, edges are assumed to have length 1.
count_clades_over_time = function(tree, Ntimes, include_slopes=FALSE){
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode;
	results = count_clades_per_time_point_CPP(	Ntips			= Ntips,
												Nnodes			= Nnodes,
												Nedges			= nrow(tree$edge),
												tree_edge		= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
												edge_length		= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
												Ntimes			= Ntimes,
												include_slopes = include_slopes);
	return(list(time_points		= results$time_points, 
				clade_counts	= as.numeric(results$clade_counts), 
				slopes			= (if(include_slopes) results$slopes else NULL),
				relative_slopes	= (if(include_slopes) results$relative_slopes else NULL)));
}