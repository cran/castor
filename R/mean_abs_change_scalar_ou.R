# compute the expected absolute change (jump) of a scalar (i.e., 1-dimensional) Ornstein-Uhlenbeck process after a specific time step delta, at stationarity
# Hence, compute: E{|X(t)-X(0)|} assuming that X(0) is at stationarity, where X is the Ornstein-Uhlenbeck process
mean_abs_change_scalar_ou = function(stationary_mean, 		# stationary mean of the OU process, i.e., the deterministic equilibrium
									stationary_std,			# stationary standard deviation of the OU process
									decay_rate,				# decay rate (aka. lambda) of the OU process
									delta,					# time step for the change
									rel_error = 0.001,		# relative tolerable standard estimation error (relative to the true mean abs change)
									Nsamples = NULL){		# number of random samples to use for estimation. If NULL, this is determined automatically based on the desired accuracy (rel_error)
	if(delta==0) return(0)
	if(is.null(Nsamples)){
		# start with a decent but small sample size, then increase it if needed
		Nsamples_here = 1000
	}else{
		Nsamples_here = Nsamples
	}
	W = (stationary_std^2)*(1 - exp(-2*decay_rate*delta)) # conditional variance of the changes. Note that this does not depend on X(0).
	# Method 1: Monte Carlo integration
	# Generate many random X(0), then compute the corresponding conditional expected absolute changes, and average those
	X0s = rnorm(n=Nsamples_here, mean=stationary_mean, sd=stationary_std)
	Ms = (stationary_mean-X0s)*(1-exp(-decay_rate*delta)) # conditional means of the changes, i.e., conditioned on the X(0)
	conditional_mean_abs_changes = sqrt(2*W/pi) * exp(-(Ms**2)/(2*W)) + Ms*(1-2*pnorm(-Ms/sqrt(W)))
	mean_abs_change = mean(conditional_mean_abs_changes)
	
	if(is.null(Nsamples)){
		# this was actually just a small trial to estimate the necessary sample size for the desired accuracy
		Nsamples = min(1000000000,(sd(conditional_mean_abs_changes)/(mean_abs_change*rel_error))**2)
		# repeat estimationg using the appropriate sample size
		if(Nsamples>Nsamples_here){
			mean_abs_change = mean_abs_change_scalar_ou(stationary_mean=stationary_mean, decay_rate=decay_rate, stationary_std=stationary_std, delta=delta, rel_error=rel_error, Nsamples=Nsamples)
		}
	}
	return(mean_abs_change)
}
