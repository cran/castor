# Given a jplace file (e.g., as generated by pplacer or EPA-NG), construct an expanded tree consisting of the original reference tree and additional tips representing the placed query sequences. 
# The reference tree and placements are loaded from the jplace file.
# If multiple placements are listed for a query, this function can either add the best (maximum-likelihood) placement or all listed placements.
expanded_tree_from_jplace = function(file_path,						# path to the jplace file
									only_best_placements	= TRUE,	# only keep the best placement of each query, i.e., the placement with maximum likelihood
									max_names_per_query		= 1){	# maximum number of sequence names to keep from each query
	J = load_jplace(file_path, only_best_placements=only_best_placements, max_names_per_query=max_names_per_query)
	results = place_tips_on_edges(tree				= J$tree,
								edges				= J$edges,
								distal_lengths		= J$distal_lengths,
								pendant_lengths		= J$pendant_lengths,
								placed_tip_labels	= J$sequence_names)
	results$reference_tree = J$tree
	return(results)
}