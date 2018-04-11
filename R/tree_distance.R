# Calculate distance (single number) between two trees.
# The trees may include multifurcations and monofurcations.
# If normalized==TRUE, the distance will always be within 0 and 1.
# tipsA2B[] should be a 1D array mapping A-tip indices to B-tip indices (one-to-one mapping). If NULL, it is determined by matching tip labels.
tree_distance = function(treeA, treeB, tipsA2B=NULL, metric="RF", normalized=FALSE){
	Ntips  = length(treeA$tip.label)
	if(Ntips!=length(treeB$tip.label)) return(list(success=FALSE, error=sprintf("Tip counts don't match: TreeA has %d tips, treeB has %d tips",Ntips,length(treeB$tip.label))))
	if(is.null(tipsA2B)){
		tipsA2B = match(treeA$tip.label, treeB$tip.label)
		if(any(is.na(tipsA2B))) return(list(success=FALSE, error=sprintf("Tip labels in treeA don't match tip labels in treeB")))
	}
	if(metric=="RF"){
		# Robinson-Foulds
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
	}else{
		stop(sprintf("Unknown metric '%s'",metric))
	}
	return(distance)
}
