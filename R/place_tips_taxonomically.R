# Given a rooted tree and taxonomies for all or some tips & nodes (typically these will be consensus taxonomies), as well as a list of additional query taxonomies (e.g., of new strains), place the queries onto the tree based on taxonomic identity.
# Each placement defines a new tip, descending directly from some node in the input tree.
# Each query is placed at the deepest possible node (furthest from the root in terms of splits) for which it is certain that the query is a descendant of.
# If include_expanded_tree, each placed query's distance from its placement node will be equal to that node's average distance to all descending tips (i.e., mean node depth).
# Examples: 
#	If a parental node has taxonomy N1="A;B;C;D" and its two children have taxonomies N2="A;B;C;D;E1" & N2="A;B;C;D;E2", then the query taxonomy "A;B;C;D;F" will be placed on the parental node N1, while the taxonomy "A;B;C;D;E1;F" would be placed at or below the node N2.
place_tips_taxonomically = function(tree,					# rooted tree of class phylo, 
									query_labels,						# character vector of length Nqueries, listing labels for the newly placed tips 
									query_taxonomies		= NULL,		# optional character vector of length Nqueries, listing query taxonomies to be placed on the tree. If NULL, it is assumed that query_labels are taxonomies.
									tip_taxonomies 			= NULL,		# optional character vector of length Ntips, listing taxonomic paths for the tips. If NULL, then tip labels are assumed to be tip taxonomies.
									node_taxonomies 		= NULL,		# optional character vector of length Nnodes, listing taxonomic paths for the nodes. If NULL, then node labels are assumed to be node taxonomies.
									tree_taxon_delimiter	= ";",		# character, the delimiter between taxonomic levels in the tip & node taxonomies (e.g., ";" for SILVA taxonomies)
									query_taxon_delimiter	= ";",		# character, the delimiter between taxonomic levels in query taxonomies
									include_expanded_tree	= TRUE){	# if TRUE, the expanded tree (i.e., including the placements) is returned as well, at some computational cost. If FALSE, only the placement info is returned, but no tree expansion is done. 
	Ntips  	 = length(tree$tip.label)
	Nnodes 	 = tree$Nnode
	Nqueries = length(query_labels)
	
	if(is.null(tip_taxonomies)) tip_taxonomies 		= tree$tip.label
	if(is.null(node_taxonomies)) node_taxonomies 	= tree$node.label
	if(is.null(query_taxonomies)) query_taxonomies 	= query_labels
	clade_taxonomies = c(tip_taxonomies, node_taxonomies)
	
	# determine where to place each query
	placements = place_tips_taxonomically_CPP(	Ntips					= Ntips,
												Nnodes					= Nnodes,
												Nedges					= nrow(tree$edge),
												tree_edge				= as.vector(t(tree$edge))-1,
												clade_taxonomies		= clade_taxonomies,
												query_taxonomies		= query_taxonomies,
												tree_taxon_delimiter	= tree_taxon_delimiter,
												query_taxon_delimiter	= query_taxon_delimiter,
												allow_placement_at_tips	= FALSE)
	placement_nodes = placements$placement_clades + 1 - Ntips
	placement_nodes[placement_nodes<=0] = 0 # use 0 for queries that could not be placed on the tree
	results = list(placement_nodes = placement_nodes)
	
	if(include_expanded_tree){
		placed_queries = which(placement_nodes>0)
		if(!is.null(tree$edge.length)){
			mean_node_depths =  get_mean_depth_per_node_CPP(Ntips		= Ntips,
															Nnodes		= Nnodes,
															Nedges		= nrow(tree$edge),
															tree_edge	= as.vector(t(tree$edge))-1,
															edge_length	= tree$edge.length)
		}
		expansion = place_tips_on_nodes(tree			  = tree,
										placement_nodes	  = placement_nodes[placed_queries],
										placed_tip_labels = query_labels[placed_queries],
										placement_lengths = (if(is.null(tree$edge.length)) NULL else mean_node_depths[placement_nodes[placed_queries]]))
		results$tree = expansion$tree
		# construct mapping placed_tips[], as an integer vector of length Nqueries, specifying the newly placed tips corresponding to the input queries. For unplaced queries, a value 0 is used.
		results$placed_tips = integer(Nqueries)
		results$placed_tips[placed_queries] = expansion$placed_tips
	}
	return(results)
}