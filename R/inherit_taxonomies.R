# Given a rooted tree and taxon names associated with some nodes and/or tips, infer the inherited taxonomic paths for all tips and/or nodes.
# Some taxa may be empty (""), i.e., not all tips/nodes are expected to have taxon names; these are not included in the returned taxonomies, not even as empty placeholders. Tips/node taxa can be specified as a separate vector or taken from the tip/node labels.
inherit_taxonomies = function(tree,						# rooted tree of class phylo, 
							  tip_taxa 		= NULL, 	# optional character vector of length Ntips, listing taxon names for the tips. These will be appended to the inherited taxonomies of the tips. If NULL, then tip labels are assumed to be taxon names. If you don't want any taxon names associated with tips, pass an empty vector. Some elements may be ""/NULL/NA.
							  node_taxa 	= NULL, 	# optional character vector of length Nnodes, listing taxon names for the nodes. If NULL, then node labels are assumed to be taxon names. Some entries may be empty/NULL/NA.
							  delimiter 	= ";",		# character, the delimiter to use between taxonomic levels (e.g., ";" for SILVA taxonomies)
							  only_proper	= FALSE,	# whether to only include proper taxon names in the inferred taxonomies
							  halt			= "last",	# when to halt accumulating taxon levels, moving from root to tips. Options are "last" (halt after the last non-empty taxon, i.e., don't skip any levels), "last_proper" (halt at the last properly name taxon, including the latter), "first_improper" (halt at the first improperly named non-empty taxon, not including the latter) or "first_empty" (halt at the first empty taxon)
							  strict_proper	= TRUE,		# whether to define taxa as proper in the strict sense, i.e., only if the taxon name is a proper taxon name. If FALSE, then a taxon name is considered proper as long as any of the descendants are proper in the strict sense; for example, the taxon "WOR-3" would be considered proper if one of its descending clades has taxon name "Hydrothermia".
							  focal_tips	= NULL,		# optional specification of tips to focus on. Either a character vector listing tip names, or an integer vector listing tip indices. If NULL, inherited taxonomies will be returned for all tips.
							  focal_nodes	= NULL){	# optional specification of nodes to focus on. Either a character vector listing node names, or an integer vector listing node indices. If empty, inherited taxonomies will be returned for all nodes.
	# basic option checks
	Ntips = length(tree$tip.label)
	halt = tolower(halt)
	if(!(halt %in% c('last', 'last_proper', 'first_improper', 'first_empty'))) stop(sprintf("Invalid choice '%s' for argument 'halt'.",halt))
	if(is.null(tip_taxa)){
		tip_taxa = tree$tip.label
	}else if((length(tip_taxa)>0) && (length(tip_taxa)!=Ntips)){
		stop(sprintf("Number of passed tip taxa (%d) does not match the number of tips in the tree (%d)",length(tip_taxa),Ntips))
	}else{
		tip_taxa[is.null(tip_taxa) | is.na(tip_taxa)] = ""	
	}
	if(is.null(node_taxa)){
		if(is.null(tree$node.label)){
			stop("Missing node taxa, since node_taxa was not provided and the tree lacks node labels.")
		}else{
			node_taxa = tree$node.label
		}
	}else if(length(node_taxa)!=tree$Nnode){
		stop(sprintf("Number of passed node taxa (%d) does not match the number of nodes in the tree (%d)",length(node_taxa),tree$Nnode))
	}else{
		node_taxa[is.null(node_taxa) | is.na(node_taxa)] = ""
	}
	
	focal_tips  = map_tip_or_node_names_to_indices(tree, focal_tips, type="tip", list_title="focal_tips", check_input=TRUE)
	focal_nodes = map_tip_or_node_names_to_indices(tree, focal_nodes, type="node", list_title="focal_nodes", check_input=TRUE)

	# determine the inherited taxonomic path for each focal clade
	results = inherit_taxonomies_CPP(Ntips 			= length(tree$tip.label),
									 Nnodes			= tree$Nnode,
									 Nedges			= nrow(tree$edge),
									 tree_edge		= as.vector(t(tree$edge)) - 1,
									 tip_taxa		= tip_taxa,
									 node_taxa		= node_taxa,
									 focal_clades	= c(focal_tips,focal_nodes)-1,
									 only_proper	= only_proper,
									 halt			= halt,
									 strict_proper	= strict_proper,
									 delimiter		= delimiter)

	return(list(tip_taxonomies	= results$taxonomies[seq_len(length(focal_tips))],
				node_taxonomies	= results$taxonomies[upward_seq(from=(length(focal_tips)+1), to=(length(focal_tips)+length(focal_nodes)))],
				tip_levels		= results$levels[seq_len(length(focal_tips))],
				node_levels		= results$levels[upward_seq(from=(length(focal_tips)+1), to=(length(focal_tips)+length(focal_nodes)))]))
}