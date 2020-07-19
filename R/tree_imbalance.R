# Given a rooted tree, calculate various measures of tree imbalance
# Supported imbalance types:
#	Sackin (sum of tip distances to root, counted in number of edges). [Sackin 1972. "Good" and "Bad" Phenograms. Systematic Biology. 21:225-226]
#	Colless (non-normalized, differences of tip_counts between a node's children summed across nodes). [Shao 1990. Tree Balance. Systematic Biology. 39:266-276]
#	Colless_normalized (normalized by the maximum possible value for bifurcating trees)
# Requirements:
#	The input tree must be rooted
#   The input tree can be multifurcating and/or monofurcating
tree_imbalance = function(tree, type){
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode;

	if(type=="Sackin"){
		imbalance = sum(get_all_distances_to_root(tree, as_edge_count=TRUE))
	}else if(type=="Colless"){
		imbalance = get_Colless_Imbalance_CPP(Ntips, Nnodes, Nedges = nrow(tree$edge), tree_edge = as.vector(t(tree$edge)) - 1, normalized=FALSE)
	}else if(type=="Colless_normalized"){
		imbalance = get_Colless_Imbalance_CPP(Ntips, Nnodes, Nedges = nrow(tree$edge), tree_edge = as.vector(t(tree$edge)) - 1, normalized=TRUE)
	}else{
		imbalance = NA
	}				
	return(imbalance)
}