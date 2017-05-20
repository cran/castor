# Fit a model of random phylogenetic tree generation, to the clades-over-time curve of another observed tree
# In the model, tips (species) are born and killed at Poissonian rates, each of which is a power-lar function of the number of extant tips
# The tip split at each speciation event, as well as the tip killed, are chosen randomly uniformly among all extant tips
# The model is parameterized via the various power-lar parameters for the birth/death rates
# Any parameter set to NULL will be fitted (within reasonable ranges). To fix a parameter, set its value explicitly.
fit_random_tree_model = function(tree, 
								 birth_rate_intercept	= NULL, 
								 birth_rate_factor 		= NULL,
								 birth_rate_exponent 	= NULL,
								 death_rate_intercept 	= NULL,
								 death_rate_factor		= NULL,
								 death_rate_exponent	= NULL,
								 Ntime_points			= 100, 	# number of time-points for clades-over-time curve to be fitted
								 Nrandom_trees			= 10,	# number of random trees to generate for each evaluation of the objective function during fitting
								 coalescent 			= TRUE,
								 Ntrials				= 1,
								 Nthreads				= 1){
	Ntips = length(tree$tip.label)
	
	
	# figure out which parameters remain to be fitted
	fitted_parameter_names = c(	(if(is.null(birth_rate_intercept)) "birth_rate_intercept" else NULL),
								(if(is.null(birth_rate_factor)) "birth_rate_factor" else NULL),
								(if(is.null(birth_rate_exponent)) "birth_rate_exponent" else NULL),
								(if(is.null(death_rate_intercept)) "death_rate_intercept" else NULL),
								(if(is.null(death_rate_factor)) "death_rate_factor" else NULL),
								(if(is.null(death_rate_exponent)) "death_rate_exponent" else NULL));
	if(is.null(fitted_parameter_names) || (length(fitted_parameter_names)==0)) stop("ERROR: All model parameters are fixed")
	NFP = length(fitted_parameter_names);
	return_value_on_failure = rep(NA, NFP)
	

	# get clade-over-time (COT) curve of observed tree
	results 		= castor::count_clades_over_time(tree, Ntimes=Ntime_points, include_slopes=FALSE);
	clade_counts 	= results$clade_counts;
	time_points 	= results$time_points;
	max_time		= max(time_points)

	# determine typical parameter values
	typical_params = unlist(list(	birth_rate_intercept 	= Ntips/max_time, 
									birth_rate_factor 		= log(Ntips)/max_time, 
									birth_rate_exponent 	= 1, 
									death_rate_intercept	= birth_rate_intercept/10, 
									death_rate_factor		= birth_rate_factor/10, 
									death_rate_exponent		= 1))

	# define some auxiliary functions
	expand_fitted_model_parameters = function(fitted_params){
		expanded_params = list(birth_rate_intercept, birth_rate_factor, birth_rate_exponent, death_rate_intercept, death_rate_factor, death_rate_exponent);
		if(length(fitted_params)==0) return(expanded_params)
		for(param_name in fitted_parameter_names) expanded_params[[param_name]] = fitted_params[param_name];
		return(expanded_params)
	}
	
	get_random_cot_curve = function(params){
		random_tree = generate_random_tree( max_tips				= NULL, 
											max_time				= max_time,
											birth_rate_intercept	= params$birth_rate_intercept, 
											birth_rate_factor 		= params$birth_rate_factor,
											birth_rate_exponent 	= params$birth_rate_exponent,
											death_rate_intercept 	= params$death_rate_intercept,
											death_rate_factor		= params$death_rate_factor,
											death_rate_exponent		= params$death_rate_exponent,
											coalescent 				= coalescent);
		cot_curve = castor::count_clades_over_time(random_tree, Ntimes=Ntime_points, include_slopes=FALSE);
		return(cot_curve)	
	}

	objective_function = function(fitted_params){
		params = expand_fitted_model_parameters(fitted_params);
		mean_random_clade_counts = rowMeans(vapply(1:Nrandom_trees, function(r) get_random_cot_curve(params)$clade_counts, FUN.VALUE=clade_counts))
		return(sum(((clade_counts-mean_random_clade_counts)/clade_counts)**2))
	}
	
	# fit with various starting points
	fit_single_trial = function(trial){
		initial_fitted_params = runif(n=NFP, min=0.1, max = 10) * typical_params[[fitted_parameter_names]];
		fit = stats::nlminb(initial_fitted_params, objective_function, lower=rep(0, NFP), upper=rep(1e+50, NFP))
		return(list(SSE=fit$objective, fit=fit));
	}
	
	# run one or more independent fitting trials
	if(Nthreads>1){
		# run trials in parallel using multiple forks
		# Note: Forks (and hence shared memory) are not available on Windows
		fits = parallel::mclapply(	1:Ntrials, 
									FUN = function(trial) fit_single_trial(trial), 
									mc.cores = min(Nthreads, Ntrials), 
									mc.preschedule = FALSE, 
									mc.cleanup = TRUE);
	}else{
		# run in serial mode
		fits = sapply(1:Ntrials,function(x) NULL)
		for(trial in 1:Ntrials){
			fits[[trial]] = fit_single_trial(trial)
		}
	}

	# extract information from best fit (note that some fits may have LL=NaN or NA)
	SSEs 				= sapply(1:Ntrials, function(trial) fits[[trial]]$SSE)
	valids				= which((!is.na(SSEs)) & (!is.nan(SSEs)) & (!is.null(SSEs)))
	if(length(valids)==0) return(return_value_on_failure); # fitting failed for all trials
	best 				= valids[which.max(sapply(valids, function(i) SSEs[i]))]
	SSE					= fits[[best]]$SSE;
	fitted_params 		= fits[[best]]$fit$par;
	params			 	= expand_fitted_model_parameters(fitted_params);			
	if(is.null(SSE) || any(is.na(fitted_params)) || any(is.nan(fitted_params))) return(return_value_on_failure); # fitting failed

	return(params);	
}
