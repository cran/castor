# Calculate distance (single number) between two trees with identical tips.
# The trees may include multifurcations and monofurcations.
# If normalized==TRUE, the distance will always be within 0 and 1.
# tipsA2B[] should be a 1D array mapping A-tip indices to B-tip indices (one-to-one mapping). If NULL, it is determined by matching tip labels.
tree_distance = function(	treeA, 
							treeB, 
							tipsA2B		= NULL,
							metric		= "RFrooted",	# which distance metric to use
							normalized	= FALSE){
	Ntips  = length(treeA$tip.label)
	if(Ntips!=length(treeB$tip.label)) stop(sprintf("Tip counts don't match: TreeA has %d tips, treeB has %d tips",Ntips,length(treeB$tip.label)))
	if(is.null(tipsA2B)){
		tipsA2B = match(treeA$tip.label, treeB$tip.label)
		if(any(is.na(tipsA2B))) stop(sprintf("Tip labels in treeA don't match tip labels in treeB"))
	}else{
		if(Ntips!=length(tipsA2B)) stop(sprintf("The length of tipsA2B (%d) differs from the number of tips (%d)",length(tipsA2B),Ntips))
	}
	if(metric=="RFrooted"){
		# Robinson-Foulds distance between rooted trees (i.e., outcome depends on rooting)
		results = get_Robinson_Foulds_distance_CPP(	Ntips		= Ntips,
													NnodesA		= treeA$Nnode,
													NedgesA		= nrow(treeA$edge),
													tree_edgeA	= as.vector(t(treeA$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
													NnodesB		= treeB$Nnode,
													NedgesB		= nrow(treeB$edge),
													tree_edgeB	= as.vector(t(treeB$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
													tipsA2B		= tipsA2B-1);
		distance = treeA$Nnode+treeB$Nnode - 2*results$Nmatches
		if(normalized) distance = distance/(treeA$Nnode+treeB$Nnode)
	}else if(metric=="MeanPathLengthDifference"){
		# path-difference metric between unrooted trees (outcome does not depend on rooting)
		# [Steel and Penny, 1993, Systematic Biology. 42:126-141.]
		tip_distancesA = as.vector(get_all_pairwise_distances(treeA, only_clades = c(1:Ntips), as_edge_counts = FALSE, check_input = FALSE))
		tip_distancesB = as.vector(get_all_pairwise_distances(treeB, only_clades = tipsA2B, as_edge_counts	= FALSE, check_input = FALSE))
		if(normalized){
			# square root of mean squared difference of normalized patristic distances between all tip pairs
			# normalization of distances is done by dividing by the maximum of the two distances being compared
			both_zero = ((tip_distancesA==0) & (tip_distancesB==0))
			tip_distancesA[both_zero] = 1 # avoid comparison when distance is zero in both trees, due to division by zero
			tip_distancesB[both_zero] = 1 # avoid comparison when distance is zero in both trees, due to division by zero
			distance = sqrt(sum(((tip_distancesA-tip_distancesB)/pmax(tip_distancesA,tip_distancesB))**2)/length(tip_distancesA))
		}else{
			# square root of mean squared difference of patristic distances between all tip pairs
			distance = sqrt(sum((tip_distancesA-tip_distancesB)**2)/length(tip_distancesA))
		}
	}else{
		stop(sprintf("Unknown metric '%s'",metric))
	}
	return(distance)
}
