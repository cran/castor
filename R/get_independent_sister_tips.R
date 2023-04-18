# Split tree into pairs of sister-tips, such that the paths within distinct pairs do not overlap
# If the input tree only contains monofurcations and bifurcations (recommended), it is guaranteed that at most one unpaired tip will be left (i.e., if Ntips was odd)
get_independent_sister_tips = function(tree){
	results = get_independent_sister_tips_CPP(	Ntips		= length(tree$tip.label),
													Nnodes		= tree$Nnode,
													Nedges		= nrow(tree$edge),
													tree_edge	= as.vector(t(tree$edge))-1);
	tip_pairs = matrix(as.integer(results$tip_pairs),ncol=2,byrow=TRUE) + 1L;
	return(tip_pairs);
}
