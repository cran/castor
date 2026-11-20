# Given a rooted tree and taxonomies for all or some tips & nodes (typically these will be consensus taxonomies), as well as a list of additional query taxonomies (e.g., of new strains), place the queries onto the tree based on taxonomic identity.
# Each placement defines a new tip, descending directly from some node in the input tree.
# Each query is placed at the deepest possible node (furthest from the root in terms of splits) for which it is certain that the query is a descendant of.
# Examples: 
#	If a parental node has taxonomy N1="A;B;C;D" and its two children have taxonomies N2="A;B;C;D;E1" & N3="A;B;C;D;E2", then the query taxonomy "A;B;C;D;F" will be placed on the parental node N1, while the taxonomy "A;B;C;D;E1;F" would be placed at or below the node N2.
place_tips_taxonomically = function(tree,					# rooted tree of class phylo, 
									query_labels,						# character vector of length Nqueries, listing labels for the newly placed tips 
									query_taxonomies		= NULL,		# optional character vector of length Nqueries, listing query taxonomies to be placed on the tree. If NULL, it is assumed that query_labels are taxonomies.
									tip_taxonomies 			= NULL,		# optional character vector of length Ntips, listing taxonomic paths for the tips. If NULL, then tip labels are assumed to be tip taxonomies.
									node_taxonomies 		= NULL,		# optional character vector of length Nnodes, listing taxonomic paths for the nodes. If NULL, then node labels are assumed to be node taxonomies.
									tree_taxon_delimiter	= ";",		# character, the delimiter between taxonomic levels in the tip & node taxonomies (e.g., ";" for SILVA taxonomies)
									query_taxon_delimiter	= ";",		# character, the delimiter between taxonomic levels in query taxonomies
									include_expanded_tree	= TRUE,		# if TRUE, the expanded tree (i.e., including the placements) is returned as well, at some computational cost. If FALSE, only the placement info is returned, but no tree expansion is done. 
									placement_length		= "mean_depth",	# how to determine the edge length leading into each placed query. Options are "mean_depth", "max_depth", "min_depth" or any non-negative numeric value. Only relevant if the tree has edge lengths and if include_expanded_tree==TRUE.
									strict					= TRUE,		# whether to strictly follow taxonomic paths, i.e., without skipping any mismatched taxon names. If there is a possibility that the same taxon name might refer to different clades depending on its ancestors, leave this at TRUE. If each taxon name defines a unique clade, this can be set to FALSE, which is more forgiving to cases where the tree's taxonomies and query taxonomies are not fully consistent.
									omit_placeholders		= NULL){	# optional regular expression for grepl for detecting and omit placeholder taxa (such as "Incertae Sedis") from the tip, node and query taxonomies. Node that this is most critical when strict=FALSE. If NULL or "", no such removal is done. For example, "^(incertae sedis|unknown( [a-z]+)?|family [ivxl]+|uncultured.*|unidentified.*|unclassified.*|unassigned|metagenome.*|environmental.*)$".
	Ntips  	 = length(tree$tip.label)
	Nnodes 	 = tree$Nnode
	Nqueries = length(query_labels)
	
	if(is.null(tip_taxonomies)) tip_taxonomies 		= tree$tip.label
	if(is.null(node_taxonomies)) node_taxonomies 	= tree$node.label
	if(is.null(query_taxonomies)) query_taxonomies 	= query_labels
	erase_placeholders = ((!is.null(omit_placeholders)) && (omit_placeholders!=""))

	# split taxonomies into their levels, in flattened form for efficiency
	# Taxon names are stripped of flanking whitespace and converted to lowercase, so that all subsequent taxon matching is case-insensitive.
	# Empty levels are omitted if drop_empty==TRUE. Placeholder taxa are omitted if requested (see omit_placeholders).
	# Returns a named list with elements taxa (character vector, listing all retained taxonomic levels of all taxonomies, in their original order) and owners (integer vector, specifying for each element in taxa which taxonomy it came from)
	split_taxonomies = function(taxonomies, delimiter, drop_empty=TRUE){
		taxonomies[is.na(taxonomies)] = ""
		levels_per_taxonomy	= strsplit(taxonomies, split=delimiter, fixed=TRUE)
		taxa				= tolower(trimws(unlist(levels_per_taxonomy, use.names=FALSE)))
		owners				= rep(seq_along(levels_per_taxonomy), times=lengths(levels_per_taxonomy))
		keep				= (if(drop_empty) (taxa!="") else rep(TRUE,length(taxa)))
		if(erase_placeholders) keep = keep & (!grepl(omit_placeholders, taxa, ignore.case=TRUE))
		return(list(taxa=taxa[keep], owners=owners[keep]))
	}

	# determine where to place each query
	if(strict){
		# normalize the tip, node & query taxonomies (strip flanking whitespace from each level, convert to lowercase, and omit placeholder taxa if requested), by splitting them into levels and joining the retained levels back together
		# Retained levels keep their original order, so the taxonomic paths can still be compared level by level. Empty levels are only omitted when placeholders are omitted, so that otherwise the positions of all levels are preserved.
		join_taxonomies = function(taxonomy_levels, Ntaxonomies, delimiter){
			joined = rep("", Ntaxonomies)
			if(length(taxonomy_levels$taxa)>0){
				pasted = vapply(split(taxonomy_levels$taxa, taxonomy_levels$owners), FUN.VALUE=character(1), FUN=paste, collapse=delimiter)
				joined[as.integer(names(pasted))] = pasted
			}
			return(joined)
		}
		clade_taxonomies = c(tip_taxonomies, node_taxonomies)
		clade_taxonomies = join_taxonomies(split_taxonomies(clade_taxonomies, tree_taxon_delimiter, drop_empty=erase_placeholders), length(clade_taxonomies), tree_taxon_delimiter)
		query_taxonomies = join_taxonomies(split_taxonomies(query_taxonomies, query_taxon_delimiter, drop_empty=erase_placeholders), Nqueries, query_taxon_delimiter)
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
	}else{
		# place each query onto the deepest tree node whose taxon matches a taxon in the query taxonomy
		# Here, a node's "taxon" is the last (lowest) level of its taxonomy. Each query is matched via its lowest-level taxon that is also the taxon of some node,
		# and is then placed at the MRCA of all nodes with that taxon, i.e., the root of that taxon's clade.
		# Placeholder taxa (e.g., "Incertae Sedis" or "Family XI", see omit_placeholders) are erased from all node & query taxonomies beforehand, since they do not define unique clades.
		# All other taxon names are assumed to be unique, i.e., to each define a single clade.
		if(is.null(node_taxonomies)) stop("Missing node taxonomies, needed for non-strict taxonomic placement")

		# determine each node's taxon, i.e. the last retained level of its taxonomy (NA if none)
		node_levels	= split_taxonomies(node_taxonomies, tree_taxon_delimiter)
		last_levels	= which(!duplicated(node_levels$owners, fromLast=TRUE))
		node_taxa	= rep(NA_character_, Nnodes)
		node_taxa[node_levels$owners[last_levels]] = node_levels$taxa[last_levels]

		# for each node taxon, find the MRCA of all nodes with that taxon
		# If the taxon forms a single clade in the tree, this is simply the root of that clade. If the taxon is split into multiple separate clades (non-monophyletic in the tree), this is their common ancestor, consistent with strict mode.
		valid_nodes		= which(!is.na(node_taxa))
		taxon2nodes		= split(Ntips+valid_nodes, node_taxa[valid_nodes]) # clade indices of all nodes with each taxon
		taxon2node		= setNames(get_multiple_mrcas(tree, descendants=taxon2nodes) - Ntips, names(taxon2nodes))

		# for the q'th query, find the lowest-level taxon in query_taxonomies[q] that matches some node's taxon, and place the query at that taxon's MRCA node
		# Queries without any matching taxon are placed at the root if the root's taxonomy is empty (consistent with strict mode), and are otherwise not placed (placement node 0)
		root_node			= find_root(tree) - Ntips
		unmatched_placement = (if(is.na(node_taxa[root_node])) as.integer(root_node) else 0L)
		query_levels	= split_taxonomies(query_taxonomies, query_taxon_delimiter)
		level2taxon		= match(query_levels$taxa, names(taxon2node)) # for each query level, the index of the matching node taxon (NA if none)
		matched_levels	= which(!is.na(level2taxon))
		last_matched	= matched_levels[!duplicated(query_levels$owners[matched_levels], fromLast=TRUE)] # the lowest matched level of each query (if any)
		placement_nodes	= rep(unmatched_placement, Nqueries)
		placement_nodes[query_levels$owners[last_matched]] = as.integer(taxon2node[level2taxon[last_matched]])
		results = list(placement_nodes = placement_nodes)
	}
	
	if(include_expanded_tree){
		placed_queries = which(placement_nodes>0)
		if(is.null(tree$edge.length)){
			lengths = NULL
		}else if(placement_length=="mean_depth"){
			mean_node_depths =  get_mean_depth_per_node_CPP(Ntips		= Ntips,
															Nnodes		= Nnodes,
															Nedges		= nrow(tree$edge),
															tree_edge	= as.vector(t(tree$edge))-1,
															edge_length	= tree$edge.length)
			lengths = mean_node_depths[placement_nodes[placed_queries]]
		}else if(placement_length %in% c("min_depth","max_depth")){
			min_max_node_depths = get_min_max_depth_per_node_CPP(Ntips		= Ntips,
																Nnodes		= Nnodes,
																Nedges		= nrow(tree$edge),
																tree_edge	= as.vector(t(tree$edge))-1,
																edge_length	= tree$edge.length)
			if(placement_length=="min_depth"){
				lengths = min_max_node_depths$min_depths[placement_nodes[placed_queries]]
			}else if(placement_length=="max_depth"){
				lengths = min_max_node_depths$max_depths[placement_nodes[placed_queries]]
			}
		}else if(is.numeric(placement_length)){
			lengths = rep(placement_length, length(placed_queries))
		}else{
			stop(sprintf("Invalid choice '%s' for placement_length",placement_length))
		}
		expansion = place_tips_on_nodes(tree			  = tree,
										placement_nodes	  = placement_nodes[placed_queries],
										placed_tip_labels = query_labels[placed_queries],
										placement_lengths = lengths)
		results$tree = expansion$tree
		# construct mapping placed_tips[], as an integer vector of length Nqueries, specifying the newly placed tips corresponding to the input queries. For unplaced queries, a value 0 is used.
		results$placed_tips = integer(Nqueries)
		results$placed_tips[placed_queries] = expansion$placed_tips
	}
	return(results)
}