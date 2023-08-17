# compute the expected geodesic distances traversed by a Spherical Brownian Motion process after specific time intervals (deltas)
expected_distances_sbm = function(diffusivity, 	# non-negative numeric, diffusivity of the SBM
								  radius,		# positive numeric, radius of the sphere
								  deltas){		# numeric vector, listing time intervals for which to compute the expected geodesic distances
	return(sapply(deltas, FUN=function(delta) expected_SBM_distance(delta, diffusivity=diffusivity, radius=radius)))
}
