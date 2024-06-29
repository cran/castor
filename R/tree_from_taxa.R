# Given a collection of taxon lists, construct a rooted taxonomic tree.
# Each taxon list is defined by an parent name and the names of its children (i.e., immediate descendants).
# Rules:
#	Each taxon must appear at most once as an parent and at most once as a child.
# 	Any taxa found as parents but not as children, will be assumed to descend directly from the root. If only one such taxon was found, it will become the root itself.
# 	Any taxa found as a child but not as an parent, will become tips.
# 	Any parents without children will be considered tips.
# Since the returned tree is a taxonomy, it will contain no edge lengths.
tree_from_taxa = function(taxa){
	# basic input checking
	NTL = length(taxa)
	if(!(class(taxa) %in% c("list"))) stop("taxa must be a list of vectors")
	parent_names = names(taxa)
	if(any(parent_names=="")) stop("One of the parent names is empty, which is not allowed")
	duplicate_parents = unique(parent_names[duplicated(parent_names)])
	if(length(duplicate_parents)>0) stop(sprintf("%d %s multiple times as parent: %s",length(duplicate_parents),ifelse(length(duplicate_parents)==1,"taxon appears","taxa appear"),paste(duplicate_parents,collapse=", ")))
	child_names	= unlist(taxa)
	duplicate_children = unique(child_names[duplicated(child_names)])
	if(length(duplicate_children)>0) stop(sprintf("%d %s multiple times as a child: %s",length(duplicate_children),ifelse(length(duplicate_children)==1,"taxon appears","taxa appear"),paste(duplicate_children,collapse=", ")))
	
	# determine tips & nodes
	taxon_sizes = sapply(taxa, FUN=function(taxon_list) length(taxon_list))
	node_names 	= parent_names[taxon_sizes>=1]
	tip_names	= c(setdiff(child_names, node_names), parent_names[taxon_sizes==0])
	root_names	= setdiff(parent_names, child_names)
	if(length(root_names)==0) stop("Missing root: All parents also appear as children")

	# determine tree dimensions
	Ntips   = length(tip_names)
	Nnodes  = ifelse(length(root_names)==1, length(node_names), length(node_names)+1)
	Nclades = Ntips + Nnodes
	Nedges  = sum(taxon_sizes) + ifelse(length(root_names)==1, 0, length(root_names))
	
	# define edges
	edges = matrix(0L, nrow=Nedges, ncol=2)
	if(length(root_names)==1){
		# found a unique root. By convention, this should be the first node in the tree.
		node_names = c(root_names, setdiff(node_names,root_names))
		clade_name2index = setNames(c(1:Nclades), nm=c(tip_names,node_names))
		next_edge = 1
	}else{
		# found multiple roots, so combine them as children of a single true root with empty name ("")
		# By convention, the root will be the first node in the tree
		node_names = c("", node_names)
		clade_name2index = setNames(c(1:Nclades), nm=c(tip_names,node_names))
		edges[1:length(root_names),1] = Ntips+1
		edges[1:length(root_names),2] = clade_name2index[root_names]
		next_edge = 1 + length(root_names)
	}
	for(tl in seq_len(NTL)){
		# define edges connecting the parent of this taxon list to all of its children
		parent_clade = clade_name2index[parent_names[tl]]
		child_clades = clade_name2index[taxa[[tl]]]
		edges[next_edge:(next_edge+taxon_sizes[tl]-1),1] = parent_clade
		edges[next_edge:(next_edge+taxon_sizes[tl]-1),2] = child_clades
		next_edge = next_edge + taxon_sizes[tl]
	}
	
	# construct the tree object
	tree = list(Nnode 		= Nnodes,
				tip.label 	= tip_names,
				node.label 	= node_names,
				edge 		= edges,
				root 		= Ntips+1)
	class(tree) = "phylo"
	attr(tree,"order") = NULL
	return(tree)
}



# construct a rooted tree from taxon lists loaded from a text file
# Each taxon list in the text file begins with the parent name, followed by the children names in separate lines.
# Taxon lists are separated by at least one empty line.
# Comments are allowed, and must be preceded by the character #.
# For example:
#	Enterobacteriaceae # family
#	Enterobacter	# genus
#	Escherichia		# genus
#	Klebsielle		# genus
#
#	Erwiniaceae		# family
#	Erwinia			# genus
#	Buchnera		# genus
#	Pantoea			# genus
#
#	Escherichia		# genus
#	E. albertii		# species
#	E. coli			# species
#	E. hermannii	# species
tree_from_taxa_file = function(	file_path,			 # path to an input text file, containing taxon lists. This file may be gzipped.
								prefix_sep = NULL){  # optional separator character, any part prior to which in taxon names will be omitted. For example, if ":", then "genus:Escherichia" will become "Escherichia", while "Staphylococcus" will remain "Staphylococcus".
	fin   = open_file(file_path, "rt")
	lines = suppressWarnings(readLines(fin, file.info(file_path)$size))
	close(fin)

	new_taxon = TRUE
	children  = character(0)
	parent 	  = NULL
	taxa 	  = list()
	for(line in lines){
		line = gsub("#.*","",line,fixed=FALSE) # remove any comments
		line = trimws(line) # remove any flanking whitespace
		
		# remove prefix is applicable
		if(!is.null(prefix_sep)){
			parts = strsplit(line, prefix_sep, fixed=TRUE)[[1]]
			if(length(parts)>1){
				line = paste(parts[2:length(parts)], collapse=prefix_sep)
			}
		}
		
		if(line==""){
			if(!is.null(parent)){
				taxa[[parent]] = children
				children = character(0)
				parent = NULL
			}
		}else if(is.null(parent)){
			parent = line
		}else{
			children = c(children, line)
		}
	}
	
	# add remaining taxon
	if(!is.null(parent)){
		taxa[[parent]] = children
		children = character(0)
		parent = NULL
	}
	
	return(list(taxa=taxa, tree=tree_from_taxa(taxa)))
	
}