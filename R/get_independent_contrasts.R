# Calculate Phylogenetic Independent Contrasts (PIC) for a single continuous trait [Felsenstein 1985, page 10]
# PICs are often used to calculate correlations between multiple traits, while accounting for shared evolutionary history at the tips.
# If the tree is bifurcating, then one PIC is returned for each node.
# If only_bifurcations==FALSE and multifurcations are present, these are internally expanded to bifurcations and an additional PIC is returned for each such bifurcation.
# Hence, the number of returned PICs is the number of bifurcations in the tree, (potentially after multifurcations have been expanded to bifurcations, if only_bifurcations==FALSE).
# References:
#    Felsenstein (1985). Phylogenies and the Comparative Method. The American Naturalist. 125:1-15.
#  Requirements:
# 	Tree can include monofurcations and multifurcations. Multifurcations are internally expanded to bifurcations. One additional PIC is returned for each such bifurcation.
# 	Tree can also include edges with length zero.
#	Tree must be rooted.
get_independent_contrasts = function(	tree, 									# a phylogenetic tree of class "phylo"
										tip_states, 							# numeric vector of size Ntips, or 2D matrix of size Ntips x Ntraits
										scaled			 			= TRUE,		# rescale PICs by the square root of their corresponding distances (typically done to standardize their variances)
										only_bifurcations 			= FALSE,	# if TRUE, then only existing bifurcating nodes are considered. Multifurcations will not be expanded.
										include_zero_phylodistances	= FALSE,	# if TRUE, then returned PICs may include cases where the phylodistance is zero (this can only happen if the tree has edges with length 0).
										check_input					= TRUE){
	Ntips  	= length(tree$tip.label)
	scalar 	= is.vector(tip_states)
	Ntraits = (if(scalar) 1 else ncol(tip_states))

	# basic input checking
	if(scalar && (length(tip_states)!=Ntips)) stop(sprintf("ERROR: Length of tip_states (%d) is not the same as the number of tips in the tree (%d)",length(tip_states),Ntips));
	if((!scalar) && (nrow(tip_states)!=Ntips)) stop(sprintf("ERROR: Number of rows in tip_states (%d) is not the same as the number of tips in the tree (%d)",nrow(tip_states),Ntips));
	if(!is.matrix(tip_states)) tip_states = as.matrix(tip_states)
	if(!is.numeric(tip_states)) stop(sprintf("ERROR: tip_states must be numeric"))
	if(check_input){
		if(scalar && (!is.null(names(tip_states))) && any(names(tip_states)!=tree$tip.label)) stop("ERROR: Names in tip_states and tip labels in tree don't match (must be in the same order).")
		if((!scalar) && (!is.null(rownames(tip_states))) && any(rownames(tip_states)!=tree$tip.label)) stop("ERROR: Row names in tip_states and tip labels in tree don't match (must be in the same order).")
	}

	results = get_phylogenetic_independent_contrasts_CPP(	Ntips						= Ntips,
															Nnodes						= tree$Nnode,
															Nedges						= nrow(tree$edge),
															Ntraits						= Ntraits,
															tree_edge					= as.vector(t(tree$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
															edge_length					= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
															tip_states					= as.vector(t(tip_states)), # flatten in row-major format
															scaled						= scaled,
															only_bifurcations			= only_bifurcations,
															include_zero_phylodistances	= include_zero_phylodistances)
	PICs = (if(scalar) results$PICs else matrix(results$PICs,ncol=Ntraits,byrow=TRUE)) # unflatten PICs if needed
	if((!scalar) && (!is.null(colnames(tip_states)))) colnames(PICs) = colnames(tip_states);
	nodes = 1 + results$nodes; # shift indices C++ to R
	if(!only_bifurcations) nodes[nodes<=0] = NA; # May include NA, indicating temporary nodes created during expansion of multifurcations
		
	return(list(PICs				= PICs,					# either a vector of size Npics, or a 2D matrix of size Npics x Ntraits
				distances			= results$distances,	# phylogenetic distances corresponding to the PICs (in units of edge lengths), analogous to the edge length spanning each PIC
				nodes				= nodes,				# indices of nodes for which PICs are returned
				root_state			= results$root_state,	# (1D array of size Ntraits) globally estimated state for the root
				root_standard_error	= results$root_standard_error,	# (1D array of size Ntraits) standard error for the root (under a Brownian motion model)
				root_CI95			= stats::qt(0.975, df=length(PICs))*results$root_standard_error));
}