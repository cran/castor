# Calculate most-recent-common-ancestors (MRCA) for one or more sets of clades
# There's some initial overhead involved with this function, but for large number of sets this becomes more efficient
# Returns an IntegerVector of size Nsets, with each entry being the clade index of the MRCA of one set
# Requirements:
#   The input tree must be rooted (root will be determined automatically, as the node that has no incoming edge)
#   The input tree can be multifurcating and/or monofurcating
get_multiple_mrcas = function(tree,
							  descendants){ # list of vectors, each specifying a set of tips/nodes to compute the MRCA of. Each descendants[s] must either be an integer vector (specifying clade indices) or a character vector (specifying tip/node names)
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode
	if(any(sapply(descendants,is.character))){
		# for at least some sets, the tips/nodes are given as names, not indices
		name2clade 			= 1:(Ntips+Nnodes)
		names(name2clade) 	= c(tree$tip.label,tree$node.label)
		descendants	= lapply(descendants, FUN=function(clades){ 
									if(is.character(clades)){ 
										clades = name2clade[clades]; 
										clades = clades[!is.na(clades)];
									}
									return(clades)})
	}
	# flatten descendants into a single integer vector, for passing to C++
	flat_descendants = flatten_list_of_lists(descendants)
	
	mrcas = get_multiple_MRCAs_CPP(	Ntips		= Ntips,
									Nnodes		= Nnodes,
									Nedges		= nrow(tree$edge),
									tree_edge	= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
									descendants	= flat_descendants$flattened-1,
									set_starts	= flat_descendants$starts-1);
	mrcas[mrcas<0] = NA # use NA to indicate failure
	mrcas = mrcas+1L
	return(mrcas);
}