# Calculate the correlation between pairwise phylogenetic & geographic distances
# Optionally, the P-value is estimated by randomly permutting tip coordinates (essentially a phylogenetic Mantel test).
correlate_phylo_geodistances = function(tree, 						# phylogenetic tree of class "phylo"
										tip_latitudes, 				# numeric vector of length Ntips, decimal latitudes of tips
										tip_longitudes,				# numeric vector of length Ntips, decimal longitudes of tips
										correlation_method,			# "pearson" or "spearman" or "kendall"
										max_phylodistance 	= Inf,	# optional upper bound for the phylodistances to consider. This may be used to examine the correlation between only closely related tips
										Npermutations 		= 1000,	# number of random tip permutations to perform, for estimating the significance (P-value) of the correlation. If 0, the P value will not be computed
										alternative 		= "right",	# which sides of the null-model-correlation's distribution to take as P-value. Either "two_sided", "left" or "right"
										radius				= 1){	# optional radius to assume for the sphere. If 1, all geodistances are measured in multiples of the sphere radius. This does not affect the correlation or P-value, but it affects the returned geodistances.
	Ntips = length(tree$tip.label)
	if(is.infinite(max_phylodistance)){
		# consider all tip pairs
		phylodistances_flat	= as.vector(get_all_pairwise_distances(tree, only_clades=seq_len(Ntips)))
		geodistances   		= radius * all_pairwise_geodesic_angles(tip_latitudes, tip_longitudes, tip_latitudes, tip_longitudes)
		geodistances_flat 	= as.vector(geodistances)
	}else{
		# consider only tip pairs with phylodistance below threshold
		phylodistances  	= get_all_pairwise_distances(tree, only_clades=seq_len(Ntips))
		only_tip_pairs		= which(phylodistances<=max_phylodistance, arr.ind=TRUE)
		only_tips1 			= only_tip_pairs[,1]
		only_tips2			= only_tip_pairs[,2]
		phylodistances_flat	= phylodistances[only_tip_pairs]
		geodistances_flat	= radius * geodesic_angles(tip_latitudes[only_tips1], tip_longitudes[only_tips1], tip_latitudes[only_tips2], tip_longitudes[only_tips2])
		rm(only_tip_pairs) # free some RAM
		invisible(gc())
	}
	correlation = stats::cor(x=phylodistances_flat, y=geodistances_flat, method=correlation_method)
	results = list(correlation=correlation, Npairs=length(phylodistances_flat))
	if(Npermutations>0){
		Pvalue = 0
		mean_random_correlation = 0
		for(r in seq_len(Npermutations)){
			# randomly shuffle tip coordinates and recompute pairwise geodistances
			random_order = sample.int(Ntips, replace=FALSE)
			if(is.infinite(max_phylodistance)){
				random_geodistances_flat = as.vector(geodistances[random_order,random_order])
			}else{
				random_geodistances_flat = radius*geodesic_angles(tip_latitudes[random_order[only_tips1]], tip_longitudes[random_order[only_tips1]], tip_latitudes[random_order[only_tips2]], tip_longitudes[random_order[only_tips2]])
			}
			# recompute correlation with randomized data and compare to true one
			random_correlation = stats::cor(x=phylodistances_flat, y=random_geodistances_flat, method=correlation_method)
			if((alternative=="two_sided") && (abs(random_correlation)>=abs(correlation))){
				Pvalue = Pvalue + 1
			}else if((alternative=="left") && (random_correlation<=correlation)){
				Pvalue = Pvalue + 1
			}else if((alternative=="right") && (random_correlation>=correlation)){
				Pvalue = Pvalue + 1
			}
			mean_random_correlation = mean_random_correlation + random_correlation
		}
		Pvalue = Pvalue / Npermutations 
		mean_random_correlation = mean_random_correlation / Npermutations
		results$Pvalue = Pvalue
		results$mean_random_correlation = mean_random_correlation
	}
	results$phylodistances	= phylodistances_flat
	results$geodistances	= geodistances_flat
	return(results)
}

