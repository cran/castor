# Read tree from a Newick-formatted string or file
read_tree = function(	string="", 
						file="", 
						include_edge_lengths 	= TRUE, 
						include_node_labels 	= TRUE, 
						underscores_as_blanks 	= FALSE, 
						check_label_uniqueness 	= FALSE){
	if(file!=""){
		if(string!="") stop("ERROR: Either string or file must be specified, but not both")
		string = readChar(file, file.info(file)$size)
	}

	results = read_Newick_string_CPP(	input = string, 
										underscores_as_blanks = underscores_as_blanks)								
	if(!results$success) stop(sprintf("ERROR: Could not parse Newick string: %s",results$error))
	if(check_label_uniqueness){
		duplicates = which(duplicated(results$tip_names))
		if(length(duplicates)>0) stop(sprintf("ERROR: Duplicate tip labels (e.g. '%s') found in input tree",results$tip_names[duplicates[1]]))
	}
	
	tree = list(Nnode 		= results$Nnodes,
				tip.label 	= results$tip_names,
				node.label 	= (if((!include_node_labels) || all(results$node_names=="")) NULL else results$node_names),
				edge 		= matrix(results$tree_edge+1L, ncol=2, byrow=TRUE), # unflatten row-major array
				edge.length = (if((!include_edge_lengths) || all(is.nan(results$edge_length))) NULL else results$edge_length),
				root 		= results$root+1L,
				root.edge	= (if(is.nan(results$root_edge)) NULL else results$root_edge))
	class(tree) = "phylo";
	
	return(tree)
}