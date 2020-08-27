# Fit a homogenous-birth-death cladogenic model to an ultrametric timetree, by estimating parameters for functional forms of lambda & mu
# An HBD model is defined by a time-dependent speciation rate (lambda), a time-dependent extinction rate (mu) and a rarefaction (rho0, sampling fraction)
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
									age0				= 0,		# non-negative numeric, youngest age (time before present) to consider when fitting and with respect to which rho is defined (rho0:=rho(age0) is the fraction of lineages extant at age0 that are included in the tree)
									lambda,							# function handle, mapping age & model_parameters to the current speciation rate, (age,param_values) --> lambda. Must be defined for all ages in [0:oldest_age] and for all parameters within the imposed bounds. Must be vectorized in the age argument, i.e. return a vector the same size as age[].
									mu					= 0,		# function handle, mapping age & model_parameters to the current extinction rate, (age,param_values) --> mu. Must be defined for all ages in [0:oldest_age] and for all parameters within the imposed bounds. Must be vectorized in the age argument, i.e. return a vector the same size as age[]. Can also be a single numeric.
									rho0				= 1,		# function handle, mapping model_parameters to the sampling fraction at age0 (aka. rarefaction), (param_values) --> rho. Must be defined for all parameters within the imposed bounds. Can also be a single numeric.
									age_grid			= NULL,		# numeric vector of size NG>=1, listing ages in ascending order, on which the lambda and mu functionals should be evaluated. This age grid must be fine enough to capture the possible variation in lambda() and mu() over time. If NULL or of length 1, then lambda & mu are assumed to be time-independent.
									condition			= "auto",	# one of "crown" or "stem" or "none" or "auto", specifying whether to condition the likelihood on the survival of the stem group or the crown group. It is recommended to use "stem" when oldest_age>root_age, and "crown" when oldest_age==root_age. This argument is similar to the "cond" argument in the R function RPANDA::likelihood_bd. Note that "crown" really only makes sense when oldest_age==root_age.
									relative_dt			= 1e-3,		# maximum relative time step allowed for integration. Smaller values increase the accuracy of the computed likelihoods, but increase computation time. Typical values are 0.0001-0.001. The default is usually sufficient.
									Ntrials				= 1,		# number of fitting trials to perform, each time starting with random parameter values
									max_start_attempts	= 1,		# number of times to attempt finding a valid start point (per trial) before giving up. Randomly choosen start parameters may result in Inf/undefined objective, so this option allows the algorithm to to keep looking for valid starting points.
									Nthreads			= 1,
									max_model_runtime	= NULL,		# maximum time (in seconds) to allocate for each likelihood evaluation. Use this to escape from badly parameterized models during fitting (this will likely cause the affected fitting trial to fail). If NULL or <=0, this option is ignored.
									fit_control			= list()){	# a named list containing options for the nlminb fitting routine (e.g. iter.max and rel.tol)
	# basic input error checking
	if(tree$Nnode<2) return(list(success = FALSE, error="Tree is too small"));
	if(age0<0) return(list(success = FALSE, error="age0 must be non-negative"));
	if((!is.null(oldest_age)) && (!is.null(age_grid)) && (tail(age_grid,1)<oldest_age)) return(list(success=FALSE, error=sprintf("Provided age grid must cover oldest_age (%g)",oldest_age)))
	if(is.null(max_model_runtime)) max_model_runtime = 0;
	if(is.null(age_grid)){
		age_grid = 0
	}else{
		if(any(diff(age_grid)<0)) age_grid = sort(age_grid) # avoid common errors where age_grid is in reverse order
	}
	max_start_attempts 	= max(1,max_start_attempts)
	Ntrials 			= max(1,Ntrials)
	Nthreads 			= max(1,Nthreads)

	# check if some of the functionals are actually fixed numbers
	model_fixed = TRUE
	if(is.numeric(lambda) && (length(lambda)==1)){
		# the provided lambda is actually a single number, so convert to a function
		lambda_value = lambda
		lambda = function(ages,params){ rep(lambda_value, times=length(ages)) }
	}else{
		model_fixed = FALSE
	}
	if(is.numeric(mu) && (length(mu)==1)){
		# the provided mu is actually a single number, so convert to a function
		mu_value = mu
		mu = function(ages,params){ rep(mu_value, times=length(ages)) }
	}else{
		model_fixed = FALSE
	}
	if(is.numeric(rho0) && (length(rho0)==1)){
		# the provided rho0 is actually a single number, so convert to a function
		rho0_value = rho0
		rho0 = function(params){ rho0_value }
	}else{
		model_fixed = FALSE
	}

	# trim tree at age0 if needed, while shifting time for the subsequent analyses (i.e. new ages will start counting at age0)
	if(age0>0){
		root_age = get_tree_span(tree)$max_distance
		if(root_age<age0) return(list(success=FALSE, error=sprintf("age0 (%g) is older than the root age (%g)",age0,root_age)));
		if((!is.null(oldest_age)) && (oldest_age<age0)) return(list(success=FALSE, error=sprintf("age0 (%g) is older than the oldest considered age (%g)",age0,oldest_age)));
		tree = trim_tree_at_height(tree,height=root_age-age0)$tree
		if(tree$Nnode<2) return(list(success = FALSE, error=sprintf("Tree is too small after trimming at age0 (%g)",age0)));
		if(!is.null(oldest_age)) oldest_age	= oldest_age - age0	
		if(!is.null(age_grid)) age_grid 	= age_grid - age0
		root_age = root_age - age0
	}

	# pre-compute some tree stats
	sorted_node_ages = sort(get_all_branching_ages(tree));
	root_age 		 = tail(sorted_node_ages,1);
	age_epsilon		 = 1e-4*mean(tree$edge.length);

	# more input error checking
	if(is.null(oldest_age)) oldest_age = root_age;
	if(!(condition %in% c("crown","stem","auto","none"))) return(list(success = FALSE, error = sprintf("Invalid condition '%s': Extected 'stem', 'crown', 'none' or 'auto'.",condition)));
	if(condition=="auto") condition = (if(abs(oldest_age-root_age)<=1e-10*root_age) "crown" else "stem")
	if((!is.null(oldest_age)) && (!is.null(age_grid)) && (tail(age_grid,1)<oldest_age)) return(list(success=FALSE, error=sprintf("Provided age grid must cover oldest_age (%g)",oldest_age)))

	# sanitize model parameters
	sanitized_params = sanitize_parameters_for_fitting(param_values, param_guess = param_guess, param_min = param_min, param_max = param_max, param_scale = param_scale)
	if(!sanitized_params$success) return(list(success=FALSE, error=sanitized_params$error))
	NP				= sanitized_params$NP
	NFP				= sanitized_params$NFP
	param_names		= sanitized_params$param_names
	param_values	= sanitized_params$param_values
	param_guess		= sanitized_params$param_guess
	param_min		= sanitized_params$param_min
	param_max		= sanitized_params$param_max
	param_scale		= sanitized_params$param_scale
	fitted_params	= sanitized_params$fitted_params
	fixed_params	= sanitized_params$fixed_params
	
	if((NFP>0) && model_fixed) return(list(success=FALSE, error="At least one model parameter is fitted, however all model aspects (lambda, mu, rho, kappa etc) are fixed"))
	
	# check if functionals are valid at least on the initial guess
	lambda_guess 	= lambda(age_grid+age0,param_guess)
	mu_guess 		= mu(age_grid+age0,param_guess)
	rho0_guess 		= rho0(param_guess)
	if(!all(is.finite(lambda_guess))) return(list(success=FALSE, error=sprintf("lambda is not a valid number for guessed parameters, at some ages")));
	if(!all(is.finite(mu_guess))) return(list(success=FALSE, error=sprintf("mu is not a valid number for guessed parameters, at some ages")));
	if(!is.finite(rho0_guess)) return(list(success=FALSE, error=sprintf("rho0 is not a valid number for guessed parameters")));
	if(length(lambda_guess)!=length(age_grid)) return(list(success=FALSE, error=sprintf("lambda function must return vectors of the same length as the input ages")));
	if(length(mu_guess)!=length(age_grid)) return(list(success=FALSE, error=sprintf("mu function must return vectors of the same length as the input ages")));
	
	# set fit-control options, unless provided by the caller
	if(is.null(fit_control)) fit_control = list()
	if(is.null(fit_control$step.min)) fit_control$step.min = 0.001
	if(is.null(fit_control$x.tol)) fit_control$x.tol = 1e-8
	if(is.null(fit_control$iter.max)) fit_control$iter.max = 1000
	if(is.null(fit_control$eval.max)) fit_control$eval.max = 2 * fit_control$iter.max * NFP


	################################
	# FITTING
	
	# objective function: negated log-likelihood
	# input argument is the subset of fitted parameters, rescaled according to param_scale
	objective_function = function(fparam_values){
		params = param_values; params[fitted_params] = fparam_values * param_scale[fitted_params];
		if(any(is.nan(params)) || any(is.infinite(params))) return(Inf); # catch weird cases where params become NaN
		if(!is.null(param_names)) names(params) = param_names;
		lambdas 	= lambda(age_grid+age0,params)
		mus 		= mu(age_grid+age0,params)
		input_rho0 	= rho0(params)
		if(!(all(is.finite(lambdas)) && all(is.finite(mus)) && is.finite(input_rho0))) return(Inf); # catch weird cases where lambda/mu/rho become NaN
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
													rarefaction			= input_rho0,
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
		
	# calculate loglikelihood for initial guess
	guess_loglikelihood = -objective_function(param_guess[fitted_params]/param_scale[fitted_params])	

	# fit with various starting points
	fit_single_trial = function(trial){
		scales		 = param_scale[fitted_params]
		lower_bounds = param_min[fitted_params]
		upper_bounds = param_max[fitted_params]
		# randomly choose start values for fitted params (keep trying up to max_start_attempts times)
		Nstart_attempts = 0
		while(Nstart_attempts<max_start_attempts){
			start_values = param_guess[fitted_params]
			if(trial>1){
				boxed_left	= which((!is.infinite(lower_bounds)) & is.infinite(upper_bounds))
				boxed_right	= which((!is.infinite(upper_bounds)) & is.infinite(lower_bounds))
				boxed_dual  = which(!(is.infinite(lower_bounds) | is.infinite(upper_bounds))); # determine fitted params that are boxed, i.e. constrained to within finite lower & upper bounds
				unboxed 	= which(is.infinite(lower_bounds) & is.infinite(upper_bounds))
				if(length(boxed_dual)>0) 	start_values[boxed_dual] = lower_bounds[boxed_dual] + (upper_bounds[boxed_dual]-lower_bounds[boxed_dual]) * runif(n=length(boxed_dual),min=0,max=1)
				if(length(unboxed)>0) 	 	start_values[unboxed]	 = 10**runif(n=length(unboxed), min=-2, max=2) * start_values[unboxed]
				if(length(boxed_left)>0) 	start_values[boxed_left] = sapply(boxed_left, FUN=function(fp) random_semiboxed_left(lower_bound=lower_bounds[fp], default=start_values[fp], typical_scale=scales[fp], orders_of_magnitude=4))
				if(length(boxed_right)>0) 	start_values[boxed_right]= sapply(boxed_right, FUN=function(fp) -random_semiboxed_left(lower_bound=-upper_bounds[fp], default=-start_values[fp], typical_scale=scales[fp], orders_of_magnitude=4))
			}
			# make sure start fparams are within bounds
			start_values 	= pmax(lower_bounds,pmin(upper_bounds,start_values))
			start_objective = objective_function(start_values/scales);
			Nstart_attempts = Nstart_attempts + 1
			if(is.finite(start_objective)) break;
		}
		# run fit
		if(is.finite(start_objective)){
			fit = tryCatch({ stats::nlminb(start_values/scales, 
								objective	= objective_function, 
								lower		= lower_bounds/scales, 
								upper		= upper_bounds/scales, 
								control		= fit_control)
							}, error = function(e){ list(objective=NaN, par=NA, convergence=1) })
			return(list(objective_value=fit$objective, fparam_values = fit$par*scales, converged=(fit$convergence==0), Niterations=fit$iterations, Nevaluations=fit$evaluations[[1]], Nstart_attempts=Nstart_attempts, start_values=start_values, start_objective=start_objective));
		}else{
			return(list(objective_value=NA, fparam_values = NA, converged=FALSE, Niterations=0, Nevaluations=0, Nstart_attempts=Nstart_attempts, start_values=start_values, start_objective=start_objective));
		}
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
	objective_values	= unlist_with_nulls(sapply(1:Ntrials, function(trial) fits[[trial]]$objective_value))
	valids				= which((!is.na(objective_values)) & (!is.nan(objective_values)) & (!is.null(objective_values)) & (!is.infinite(objective_values)));
	if(length(valids)==0) return(list(success=FALSE, error=sprintf("Fitting failed for all trials")));
	best 				= valids[which.min(sapply(valids, function(i) objective_values[i]))]
	objective_value		= -fits[[best]]$objective_value;
	loglikelihood		= objective_value;
	fitted_param_values = param_values; fitted_param_values[fitted_params] = fits[[best]]$fparam_values;
	if(is.null(objective_value) || any(is.na(fitted_param_values)) || any(is.nan(fitted_param_values))) return(list(success=FALSE, error=sprintf("Some fitted parameters are NaN")));
	if(!is.null(param_names)) names(fitted_param_values) = param_names;
		
	# reverse any time shift due to earlier tree trimming
	age_grid = age_grid + age0

	# return results
	return(list(success					= TRUE,
				objective_value			= objective_value,
				objective_name			= "loglikelihood",
				loglikelihood			= loglikelihood,
				param_fitted			= fitted_param_values,
				param_guess				= param_guess,
				guess_loglikelihood		= guess_loglikelihood,
				NFP						= NFP,
				AIC						= 2*NFP - 2*loglikelihood,
				BIC						= log(sum((sorted_node_ages<=oldest_age) & (sorted_node_ages>=age0)))*NFP - 2*loglikelihood,
				condition				= condition,
				converged				= fits[[best]]$converged,
				Niterations				= fits[[best]]$Niterations,
				Nevaluations			= fits[[best]]$Nevaluations,
				trial_start_objectives	= -unlist_with_nulls(sapply(1:Ntrials, function(trial) fits[[trial]]$start_objective)),
				trial_objective_values	= -objective_values,
				trial_Nstart_attempts	= unlist_with_nulls(sapply(1:Ntrials, function(trial) fits[[trial]]$Nstart_attempts)),
				trial_Niterations		= unlist_with_nulls(sapply(1:Ntrials, function(trial) fits[[trial]]$Niterations)),
				trial_Nevaluations		= unlist_with_nulls(sapply(1:Ntrials, function(trial) fits[[trial]]$Nevaluations))));
}

