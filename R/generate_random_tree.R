# generate a random phylogenetic tree, by randomly splitting tips at a certain rate and ranodmly killing tips at a certain rate
# birth & death rates (=speciation & extintion rates) can be arbitrary power-law functions of extant_tip_counts
# For example: birth_rate = intercept + factor * extant_tip_count^exponent
# The simulation is halted as soon as Ntips>=max_tips (if max_tips>0) and/or time>=max_time (if max_time>0)
generate_random_tree = function( max_tips, 
								 max_time				= NULL,
								 birth_rate_intercept	= 0, 
								 birth_rate_factor 		= 0,
								 birth_rate_exponent 	= 1,
								 death_rate_intercept 	= 0,
								 death_rate_factor		= 0,
								 death_rate_exponent	= 1,
								 coalescent 			= TRUE){
	if(is.null(max_tips) && is.null(max_time)) stop("ERROR: At least one of max_tips and/or max_time must be non-NULL")
	
	results = generate_random_tree_CPP(	max_tips				= (if(is.null(max_tips)) -1 else max_tips),
										max_time				= (if(is.null(max_time)) -1 else max_time),
										birth_rate_intercept 	= birth_rate_intercept, 
										birth_rate_factor 		= birth_rate_factor,
										birth_rate_exponent 	= birth_rate_exponent, 
										death_rate_intercept 	= death_rate_intercept,
										death_rate_factor		= death_rate_factor,
										death_rate_exponent		= death_rate_exponent,
										coalescent				= coalescent);
	tree = list(Nnode 		= results$Nnodes,
				tip.label 	= as.character(1:results$Ntips),
				node.label 	= NULL,
				edge 		= matrix(results$tree_edge,ncol=2,byrow=TRUE) + 1,
				edge.length = results$edge_length,
				root 		= results$root+1)
	class(tree) = "phylo";
	return(tree);
	
}