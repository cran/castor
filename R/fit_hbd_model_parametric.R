# Fit a homogenous-birth-death cladogenic model to an ultrametric timetree, by estimating parameters for functional forms of lambda & mu
# An HBD model is defined by a time-dependent speciation rate (lambda), a time-dependent extinction rate (mu) and a rarefaction (rho, subsampling fraction)
#
# References:
#	Morlon et al. (2011). Reconciling molecular phylogenies with the fossil record. PNAS 108:16327-16332
fit_hbd_model_parametric = function(tree, 
									param_values,					# numeric vector of size NP, specifying fixed values for a some or all parameters. For fitted (i.e. non-fixed) parameters, use NaN or NA.
									param_guess			= NULL,		# numeric vector of size NP, listing an initial guess for each parameter. For fixed parameters, guess values are ignored.
									param_min			= -Inf,		# numeric vector of size NP, specifying lower bounds for the model parameters. For fixed parameters, bounds are ignored. May also be a single scalar, in which case the same lower bound is assumed for all params.
									param_max			= +Inf,		# numeric vector of size NP, specifying upper bounds for the model parameters. For fixed parameters, bounds are ignored. May also be a single scalar, in which case the same upper bound is assumed for all params.
									param_scale			= NULL,		# numeric vector of size NP, specifying typical scales for the model parameters. For fixed parameters, scales are ignored. If NULL, scales are automatically estimated from other information (such as provided guess and bounds). May also be a single scalar, in which case the same scale is assumed for all params.
									oldest_age			= NULL,		# either a numeric specifying the stem age or NULL (equivalent to the root age). This is similar to the "tot_time" option in the R function RPANDA::likelihood_bd
									lambda,							# function handle, mapping age & model_parameters to the current speciation rate, (age,param_values) --> lambda. Must be defined for all ages in [0:oldest_age] and for all parameters within the imposed bounds. Must be vectorized in the age argument, i.e. return a vector the same size as age[].
									mu,								# function handle, mapping age & model_parameters to the current extinction rate, (age,param_values) --> mu. Must be defined for all ages in [0:oldest_age] and for all parameters within the imposed bounds. Must be vectorized in the age argument, i.e. return a vector the same size as age[].
									rho,							# function handle, mapping model_parameters to the sampling fraction (aka. rarefaction), (param_values) --> rho. Must be defined for all parameters within the imposed bounds.
									age_grid			= NULL,		# numeric vector of size NG>=1, listing ages in ascending order, on which the lambda and mu functionals should be evaluated. This age grid must be fine enough to capture the possible variation in lambda() and mu() over time. If NULL or of length 1, then lambda & mu are assumed to be time-independent.
									condition			= "stem",	# one of "crown" or "stem", specifying whether to condition the likelihood on the survival of the stem group or the crown group. It is recommended to use "stem" when oldest_age>root_age, and "crown" when oldest_age==root_age. This argument is similar to the "cond" argument in the R function RPANDA::likelihood_bd. Note that "crown" really only makes sense when oldest_age==root_age.
									relative_dt			= 1e-3,		# maximum relative time step allowed for integration. Smaller values increase the accuracy of the computed likelihoods, but increase computation time. Typical values are 0.0001-0.001. The default is usually sufficient.
									Ntrials				= 1,
									Nthreads			= 1,
									max_model_runtime	= NULL,		# maximum time (in seconds) to allocate for each likelihood evaluation. Use this to escape from badly parameterized models during fitting (this will likely cause the affected fitting trial to fail). If NULL or <=0, this option is ignored.
									fit_control			= list()){	# a named list containing options for the nlminb fitting routine (e.g. iter.max and rel.tol)
	Ntips	= length(tree$tip.label);
	Nnodes	= tree$Nnode;
	NP		= length(param_values);

	# pre-compute some tree stats
	lineage_counter  = count_lineages_through_time(tree, Ntimes=log2(Ntips), include_slopes=TRUE);
	sorted_node_ages = sort(get_all_branching_ages(tree));
	root_age 		 = tail(sorted_node_ages,1);
	age_epsilon		 = 1e-4*mean(tree$edge.length);

	# basic input control
	if((Ntips<2) || (Nnodes<2)) return(list(success = FALSE, error="Tree is too small"));
	if(is.null(oldest_age)) oldest_age = root_age;
	if(Ntrials<1) return(list(success = FALSE, error = sprintf("Ntrials must be at least 1")))
	if(!(condition %in% c("stem","crown","none"))) return(list(success = FALSE, error = sprintf("Invalid condition option '%s'; expected either 'crown', or 'stem' or 'none' (not recommended)",condition)))
	if(is.null(age_grid)) age_grid = 0;
	param_names = names(param_values);
	if(is.null(param_guess)){
		if(any(is.finite(param_values))){
			return(list(success=FALSE, error=sprintf("Missing guessed parameter values")))
		}else{
			param_guess = rep(NA, times=NP);
		}
	}
	if(length(param_guess)!=NP){
		return(list(success=FALSE, error=sprintf("Number of guessed parameters (%d) differs from number of model parameters (%d)",length(param_guess),NP)))
	}else if(!is.null(param_names)){
		names(param_guess) = param_names;
	}
	if((!is.null(param_names)) && (length(param_names)!=NP)){
		return(list(success=FALSE, error=sprintf("Number of parameter names (%d) differs from number of model parameters (%d)",length(param_names),NP)))
	}
	if(is.null(param_min)){
		param_min = rep(-Inf,times=NP);
	}else if(length(param_min)==1){
		param_min = rep(param_min,times=NP);
	}else if(length(param_min)!=NP){
		return(list(success=FALSE, error=sprintf("Length of param_min[] (%d) differs from number of model parameters (%d)",length(param_min),NP)))
	}
	if(is.null(param_max)){
		param_max = rep(+Inf,times=NP);
	}else if(length(param_max)==1){
		param_max = rep(param_max,times=NP);
	}else if(length(param_max)!=NP){
		return(list(success=FALSE, error=sprintf("Length of param_max[] (%d) differs from number of model parameters (%d)",length(param_max),NP)))
	}
	if(is.null(param_scale)){
		param_scale = rep(NA,times=NP);
	}else if(length(param_scale)==1){
		param_scale = rep(param_scale,times=NP);
	}else if(length(param_scale)!=NP){
		return(list(success=FALSE, error=sprintf("Length of param_scale[] (%d) differs from number of model parameters (%d)",length(param_scale),NP)))
	}
	if(is.null(max_model_runtime)) max_model_runtime = 0;
	if(any(is.nan(param_guess) | is.na(param_guess))) return(list(success=FALSE, error=sprintf("Some guessed parameter values are NA or NaN; you must specify a valied guess for each model parameter")));
	param_values[is.nan(param_values)] = NA # standardize representation of non-fixed params
	param_scale[is.nan(param_scale)] = NA	# standardize representation of unknown param scales
	if(any((!is.na(param_scale)) & (param_scale==0))) return(list(success=FALSE, error=sprintf("Some provided parameter scales are zero; expecting non-zero scale for each parameter")));
	
	# check if functionals are valid at least on the initial guess
	lambda_guess 	= lambda(0,param_guess)
	mu_guess 		= mu(0,param_guess)
	rho_guess 		= rho(param_guess)
	if(!is.finite(lambda_guess)) return(list(success=FALSE, error=sprintf("lambda is not a valid number for guessed parameters, at age 0")));
	if(!is.finite(mu_guess)) return(list(success=FALSE, error=sprintf("mu is not a valid number for guessed parameters, at age 0")));
	if(!is.finite(rho_guess)) return(list(success=FALSE, error=sprintf("rho is not a valid number for guessed parameters")));
						
	#################################
	# PREPARE PARAMETERS TO BE FITTED
	
		
	# determine which parameters are to be fitted
	# convention: parameters are indexed as follows: [lambda[], mu[], rho]
	fitted_params	= which(is.na(param_values))
	fixed_params	= which(!is.na(param_values))
	NFP				= length(fitted_params);
	param_guess[fixed_params] = param_values[fixed_params] # make sure guessed param values are consistent with fixed param values
	
	# determine typical parameter scales
	for(p in fitted_params){
		if(is.na(param_scale[p])){
			if(param_guess[p]!=0){
				param_scale[p] = abs(param_guess[p]);
			}else if((is.finite(param_min[p]) && (param_min[p]!=0)) || (is.finite(param_max[p]) && (param_max[p]!=0))){
				param_scale[p] = mean(abs(c((if(is.finite(param_min[p]) && (param_min[p]!=0)) param_min[p] else NULL), (if(is.finite(param_max[p]) && (param_max[p]!=0)) param_max[p] else NULL))));
			}else{
				param_scale[p] = 1;
			}
		}
	}
	

	################################
	# FITTING
	
	# objective function: negated log-likelihood
	# input argument is the subset of fitted parameters, rescaled according to param_scale
	objective_function = function(fparam_values){
		params = param_values; params[fitted_params] = fparam_values * param_scale[fitted_params];
		if(any(is.nan(params)) || any(is.infinite(params))) return(Inf); # catch weird cases where params become NaN
		if(!is.null(param_names)) names(params) = param_names;
		lambdas = lambda(age_grid,params)
		mus 	= mu(age_grid,params)
		rho 	= rho(params)
		if(!(all(is.finite(lambdas)) && all(is.finite(mus)) && is.finite(rho))) return(Inf); # catch weird cases where lambda/mu/rho become NaN
		if(length(age_grid)==1){
			# while age-grid has only one point (i.e., lambda & mu are constant over time), we need to provide a least 2 grid points to the loglikelihood calculator, spanning the interval [0,oldest_age]
			input_age_grid 	= c(0,oldest_age);
			input_lambdas	= c(lambdas, lambdas);
			input_mus		= c(mus, mus);
		}else{
			input_age_grid 	= age_grid;
			input_lambdas	= lambdas
			input_mus 		= mus
		}
		results = get_HBD_model_loglikelihood_CPP(	branching_ages		= sorted_node_ages,
													oldest_age			= oldest_age,
													rarefaction			= rho,
													age_grid 			= input_age_grid,
													lambdas 			= input_lambdas,
													mus 				= input_mus,
													splines_degree		= 1,
													condition			= condition,
													relative_dt			= relative_dt,
													runtime_out_seconds	= max_model_runtime);
		if(!results$success) return(Inf);
		LL = results$loglikelihood;
		if(is.na(LL) || is.nan(LL) || is.infinite(LL)) return(Inf);
		return(-LL);
	}
	

	# fit with various starting points
	fit_single_trial = function(trial){
		scales		 = param_scale[fitted_params]
		lower_bounds = param_min[fitted_params]
		upper_bounds = param_max[fitted_params]
		# randomly choose start values for fitted params
		start_values = param_guess[fitted_params]
		if(trial>1){
			boxed   = which(!(is.infinite(lower_bounds) || is.infinite(upper_bounds))); # determine fitted params that are boxed, i.e. constrained to within finite lower & upper bounds
			unboxed = completement(NFP, boxed);
			if(length(boxed)>0) start_values[boxed] = lower_bounds[boxed] + (upper_bounds[boxed]-lower_bounds[boxed]) * runif(n=length(boxed),min=0,max=1)
			if(length(unboxed)>0) start_values[unboxed]	= 10**runif(n=length(unboxed), min=-2, max=2) * start_values[unboxed]
		}
		# make sure start fparams are within bounds
		start_values = pmax(lower_bounds,pmin(upper_bounds,start_values))
		# run fit
		fit = stats::nlminb(start_values/scales, 
							objective	= objective_function, 
							lower		= lower_bounds/scales, 
							upper		= upper_bounds/scales, 
							control		= fit_control)
		return(list(objective_value=fit$objective, fparam_values = fit$par*scales, converged=(fit$convergence==0), Niterations=fit$iterations, Nevaluations=fit$evaluations[1]));
	}
	
	################################

	# run one or more independent fitting trials
    if((Ntrials>1) && (Nthreads>1) && (.Platform$OS.type!="windows")){
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
	objective_values	= sapply(1:Ntrials, function(trial) fits[[trial]]$objective_value);
	valids				= which((!is.na(objective_values)) & (!is.nan(objective_values)) & (!is.null(objective_values)) & (!is.infinite(objective_values)));
	if(length(valids)==0) return(list(success=FALSE, error=sprintf("Fitting failed for all trials")));
	best 				= valids[which.min(sapply(valids, function(i) objective_values[i]))]
	objective_value		= -fits[[best]]$objective_value;
	loglikelihood		= objective_value;
	fitted_param_values = param_values; fitted_param_values[fitted_params] = fits[[best]]$fparam_values;
	if(is.null(objective_value) || any(is.na(fitted_param_values)) || any(is.nan(fitted_param_values))) return(list(success=FALSE, error=sprintf("Some fitted parameters are NaN")));
	if(!is.null(param_names)) names(fitted_param_values) = param_names;
		
	# return results
	return(list(success					= TRUE,
				objective_value			= objective_value,
				objective_name			= "loglikelihood",
				loglikelihood			= loglikelihood,
				param_fitted			= fitted_param_values,
				param_guess				= param_guess,
				NFP						= NFP,
				AIC						= 2*NFP - 2*loglikelihood,
				converged				= fits[[best]]$converged,
				Niterations				= fits[[best]]$Niterations,
				Nevaluations			= fits[[best]]$Nevaluations));
}



completement = function(N, indices){
	pool = rep(TRUE,N);
	pool[indices] = FALSE;
	return(which(pool));
}


