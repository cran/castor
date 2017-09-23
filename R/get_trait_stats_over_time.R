# calculate statistics (mean & standard deviation) of a numeric trait on a tree over time
# if tree$edge.length is missing, edges are assumed to have length 1.
get_trait_stats_over_time = function(	tree, 
										states,					# 1D numeric vector of size Nclades, listing the trait's value for each tip & node
										Ntimes		= NULL, 	# number of equidistant time points for which to calculate trait stats
										times		= NULL, 	# 1D array of time points in increasing order, for which to calculate trait stats
										check_input	= TRUE){
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode;
	
	# basic error checking
	if((!is.null(Ntimes)) && (!is.null(times))) stop("ERROR: Either Ntimes or times must be non-NULL, but not both")
	if(is.null(Ntimes) && is.null(times)) stop("ERROR: Both Ntimes and times are NULL; please specify one of the two")
	if(length(states)!=(Ntips+Nnodes)) stop("ERROR: Number of provided states does not much number of tips+nodes in the tree")
	if(check_input){
		if((!is.null(names(states))) && any(names(states[1:Ntips])!=tree$tip.label)) stop("ERROR: Names in states and tip labels in tree don't match (must be in the same order).")
		if((!is.null(names(states))) && (!is.null(tree$node.label)) && any(names(states[(Ntips+1):(Ntips+Nnodes)])!=tree$node.label)) stop("ERROR: Names in states and node labels in tree don't match (must be in the same order).")
	}
	
	if(is.null(times)){
		tree_span = castor::get_tree_span(tree)
		times = seq(from=0, to=(1-1e-8)*tree_span$max_distance, length.out=Ntimes)
	}else{
		Ntimes = length(times)
	}
	trait_stats = get_trait_stats_at_times_CPP(	Ntips 			= Ntips,
												Nnodes 			= Nnodes,
												Nedges 			= nrow(tree$edge),
												tree_edge		= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
												edge_length		= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
												times			= times,
												states			= states)
	
	return(list(Ntimes			= length(times),
				times			= times,
				clade_counts	= trait_stats$clade_counts,
				means 			= trait_stats$means,	# arithmetic means of trait
				stds 			= trait_stats$stds))	# populaton standard deviations of trait
}