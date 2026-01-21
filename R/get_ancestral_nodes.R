# Given a rooted tree and one or more "descendant" tips or nodes, determine their ancestral node indices, traveling N splits back in time (i.e., N splits towards the root).
get_ancestral_nodes = function(	tree,			# rooted phylogenetic tree of type "phylo"
								descendants, 	# integer or character vector of length ND, specifying tip/node indices or names
								Nsplits	   = 1,		# either a single integer or an integer vector of length ND, with values >=1, specifying how many splits backward to travel. For example, Nsplits=1 will yield the immediate ancestral (i.e., parent) nodes.
								only_nodes = NULL,	# optional integer vector with values in 1,..,Nnodes, specifying nodes to restrict the back-traversal to. Other nodes will be skipped and not counted towards Nsplits. If NULL or empty, this filter is ignored.
								omit_nodes = NULL){	# optional integer vector with values in 1,..,Nnodes, specifying nodes to exlude from the back-traversal. These nodes will be skipped and not counted towards Nsplits.
	Ntips  = length(tree$tip.label);
	Nnodes = tree$Nnode;
	descendants = map_tip_or_node_names_to_indices(tree, descendants, type="both", list_title="descendants", check_input=TRUE)
	ancestors = get_ancestral_nodes_CPP(Ntips			= Ntips,
										Nnodes			= Nnodes,
										Nedges			= nrow(tree$edge),
										tree_edge		= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
										descendants		= descendants-1,
										Nsplits			= Nsplits,
										only_nodes		= (if(is.null(only_nodes)) integer(0) else (only_nodes-1)),
										omit_nodes		= (if(is.null(omit_nodes)) integer(0) else (omit_nodes-1)))
	return(ancestors+1L)
}