# Given a timetree, estimate the tip densities and node densities (tips or nodes per time interval) on an age grid.
# Optionally, the densities can be normalized by the local number of lineages.
# If the tree is full (includes extinct + extant clades), then the normalized tip (node) density is an estimate for the per-capita extinction (speciation) rate.
# age = distance from youngest tip
# Nages = number of time points to consider, spanning the youngest tip to the root.
clade_densities = function(	tree, 
							Nbins			= NULL, 	# number of equidistant age bins for which to calculate densities
							min_age			= NULL,		# minimum age to consider. If NULL, will be set to the minimum possible
							max_age			= NULL,		# maximum age to consider. If NULL, will be set to the maximum possible
							normalize		= TRUE,		# logical, whether to normalize densities by the local number of lineages (in addition to dividing by the age interval)
							ultrametric		= FALSE){	# logical, whether the input tree is guaranteed to be ultrametric (even in the presence of numerical inaccuracies)			
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode
	
	# determine the age bins (left edges)
	root_age = castor::get_tree_span(tree)$max_distance
	if(is.null(min_age)) min_age = 0
	if(is.null(max_age)) max_age = root_age
	age_mins = seq(from=min_age, to=max_age, length.out=(Nbins+1))[1:Nbins]
	age_maxs = c(age_mins[2:Nbins], max_age)
	ages	 = 0.5*(age_mins + age_maxs)
	
	# compute the ages of tips & nodes
	clade_ages = root_age - get_all_distances_to_root(tree)
	tip_ages   = clade_ages[1:Ntips]
	node_ages  = clade_ages[(Ntips+1):(Ntips+Nnodes)]
	
	# bin tip ages
	tip_binning = place_sorted_values_into_bins_CPP(items = sort(tip_ages), bin_mins = age_mins, bin_maxs = age_maxs)
	tip_densities = sapply(seq_len(Nbins), FUN=function(b) length(tip_binning$bin2items[[b]]))/(age_maxs-age_mins)

	# bin node ages
	node_binning = place_sorted_values_into_bins_CPP(items = sort(node_ages), bin_mins = age_mins, bin_maxs = age_maxs)
	node_densities = sapply(seq_len(Nbins), FUN=function(b) length(node_binning$bin2items[[b]]))/(age_maxs-age_mins)
	
	if(normalize){
		# normalize densities by the local number of lineages
		LTT 		   = count_lineages_through_time(tree, times=root_age-ages, ultrametric=ultrametric)
		tip_densities  = tip_densities/LTT$lineages
		node_densities = node_densities/LTT$lineages
	}

	return(list(Nbins			= Nbins,
				ages			= ages,
				tip_densities	= tip_densities,
				node_densities	= node_densities))
}