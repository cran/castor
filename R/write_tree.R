# Write tree in Newick format into a string or to a file
# The tree does not need to be rooted. In that case, it is rooted temporarily at the first node.
write_tree = function(tree, file="", append=FALSE, digits=10){
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode;
	result = tree_to_Newick_string_CPP(	Ntips,
										Nnodes,
										Nedges				= nrow(tree$edge),
										tree_edge			= as.vector(t(tree$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
										edge_length			= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
										tip_labels			= (if(is.null(tree$tip.label)) character() else tree$tip.label),
										node_labels			= (if(is.null(tree$node.label)) character() else tree$node.label),
										digits				= digits,
										root_edge_length	= (if(is.null(tree$root.edge)) -1 else tree$root.edge));
	if(file!=""){
		cat(result, file=file, append=append);
	}else{
		return(result);
	}
}