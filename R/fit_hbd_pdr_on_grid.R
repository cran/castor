# Fit a homogenous-birth-death cladogenic model-congruence-class to an ultrametric timetree, by fitting the pulled diversification rate (PDR)
# An HBD model is defined by a time-dependent speciation rate (lambda), a time-dependent extinction rate (mu) and a rarefaction (rho, subsampling fraction)
# However, for each specific model and a given timetree there exists a continuum of alternative models that would all generate the same deterministic lineages-through-time (LTT) curve (when calculated backward in time), and all of these models actually have the same likelihood.
# Hence, each model is part of an "equivalence class" of models, and likelihood-based approaches can only discern between model classes, but not between the individual model members in a class
# It turns out that each HBD model-class is uniquely defined by its "pulled diversification rate" (PDR) and the product rho*lambda(0)=:rholambda0.
# This function thus fits model-classes, rather than models, by fitting the PDR and the parameter rholambda0.
#
# References:
#	Morlon et al. (2011). Reconciling molecular phylogenies with the fossil record. PNAS 108:16327-16332
fit_hbd_pdr_on_grid = function(	tree, 
									oldest_age			= NULL,		# either a numeric specifying the stem age or NULL (equivalent to the root age). This is similar to the "tot_time" option in the R function RPANDA::likelihood_bd
									age_grid			= NULL,		# either NULL, or a numeric vector of size NG, listing ages in ascending order, on which the PDR is defined as a piecewise linear curve. If NULL, the PDR is assumed to be time-independent.
									min_PDR				= -Inf,		# optional lower bound for the fitted PDRs. Either a single numeric (applying to all age-grid-points) or a numeric vector of size NG, specifying the lower bound at each age-grid point.
									max_PDR				= +Inf,		# optional upper bound for the fitted PDRs. Either a single numeric (applying to all age-grid-points) or a numeric vector of size NG, specifying the upper bound at each age-grid point.
									min_rholambda0		= 1e-10,	# optional lower bound for the fitted rholambda0. Note that rholambda0 is always non-negative
									max_rholambda0		= +Inf,		# optional upper bound for the fitted rholambda0
									guess_PDR			= NULL,		# initial guess for the PDR. Either NULL (an initial guess will be computed automatically), or a single numeric (guessing a constant PDR at all ages) or a numeric vector of size NG specifying an initial guess for the PDR at each age-grid point (can include NAs)
									guess_rholambda0	= NULL,		# initial guess for the product rho*lambda(0). Either NULL (an initial guess will be computed automatically) or a single strictly-positive numeric.
									fixed_PDR			= NULL,		# optional fixed PDR values, on one or more of the age grid points. Either NULL (none of the PDRs are fixed), or a single scalar (all PDRs are fixed) or a numeric vector of size NG (some or all PDRs are fixed, can include NAs).
									fixed_rholambda0	= NULL,		# optional fixed value for rholambda0. If non-NULL and non-NA, then rholambda0 is not fitted. 
									splines_degree		= 1,		# integer, either 1 or 2 or 3, specifying the degree for the splines defined by the PDR on the age grid.
									condition			= "stem",	# one of "crown" or "stem", specifying whether to condition the likelihood on the survival of the stem group or the crown group. It is recommended to use "stem" when oldest_age>root_age, and "crown" when oldest_age==root_age. This argument is similar to the "cond" argument in the R function RPANDA::likelihood_bd. Note that "crown" really only makes sense when oldest_age==root_age.
									relative_dt			= 1e-3,		# maximum relative time step allowed for integration. Smaller values increase the accuracy of the computed likelihoods, but increase computation time. Typical values are 0.0001-0.001. The default is usually sufficient.
									Ntrials				= 1,
									Nthreads			= 1,
									max_model_runtime	= NULL,		# maximum time (in seconds) to allocate for each likelihood evaluation. Use this to escape from badly parameterized models during fitting (this will likely cause the affected fitting trial to fail). If NULL or <=0, this option is ignored.
									fit_control			= list()){	# a named list containing options for the nlminb fitting routine (e.g. iter.max and rel.tol)
	Ntips	= length(tree$tip.label);
	Nnodes	= tree$Nnode;

	# pre-compute some tree stats
	lineage_counter 	= count_lineages_through_time(tree, Ntimes=log2(Ntips), include_slopes=TRUE);
	sorted_node_ages	= sort(get_all_branching_ages(tree));
	root_age 		 	= tail(sorted_node_ages,1);
	age_epsilon		 	= 1e-4*mean(tree$edge.length);

	# basic error checking
	if((Ntips<2) || (Nnodes<2)){
		# tree is trivial (~empty)
		return(list(success = FALSE, error="Tree is too small"));
	}
	if(Ntrials<1) return(list(success = FALSE, error = sprintf("Ntrials must be at least 1")))
	if(is.null(oldest_age)) oldest_age = root_age;
	if(is.null(age_grid)){
		if((!is.null(guess_PDR)) && (length(guess_PDR)>1)) return(list(success = FALSE, error = sprintf("Invalid number of guessed PDRs; since no age grid was provided, you must provide a single (constant) guess_PDR or none at all")));
		age_grid = 0 # single-point grid, means that PDRs are assumed time-independent
		NG = 1
	}else{
		NG = length(age_grid)
		if((!is.null(guess_PDR)) && (length(guess_PDR)!=1) && (length(guess_PDR)!=NG)) return(list(success = FALSE, error = sprintf("Invalid number of guessed PDRs (%d); since an age grid of size %d was provided, you must either provide one or %d PDRs",length(guess_PDR),NG)));
		if((length(age_grid)>1) && (age_grid[NG]>oldest_age-1e-5*(age_grid[NG]-age_grid[NG-1]))) age_grid[NG] = max(age_grid[NG],oldest_age); # if age_grid "almost" covers oldest_age (i.e. up to rounding errors), then fix the remaining difference
		if((length(age_grid)>1) && ((age_grid[1]>0) || (age_grid[NG]<oldest_age))) return(list(success = FALSE, error=sprintf("Provided age-grid range (%g - %g) does not cover entire required age range (0 - %g)",age_grid[1],tail(age_grid,1),oldest_age)));
	}
	if(is.null(max_model_runtime)) max_model_runtime = 0;
	if(!(splines_degree %in% c(0,1,2,3))) return(list(success = FALSE, error = sprintf("Invalid splines_degree: Extected one of 0,1,2,3.")));
	if(NG==1) splines_degree = 1; # no point in using splines since PDR is assumed to be time-independent
	
	# reformat shape of input params to an internally standardized format
	min_rholambda0 = max(0,min_rholambda0);
	max_rholambda0 = max(0,max_rholambda0);
	if(length(min_PDR)==1) min_PDR = rep(min_PDR,times=NG);
	if(length(max_PDR)==1) max_PDR = rep(max_PDR,times=NG);
	if(is.null(guess_rholambda0)) guess_rholambda0 = NA;
	if(is.null(fixed_rholambda0)) fixed_rholambda0 = NA;
	if(is.null(guess_PDR)){
		guess_PDR = rep(NA,times=NG);
	}else if(length(guess_PDR)==1){
		guess_PDR = rep(guess_PDR,times=NG);
	}
	if(is.null(fixed_PDR)){
		fixed_PDR = rep(NA,times=NG);
	}else if(length(fixed_PDR)==1){
		fixed_PDR = rep(fixed_PDR,times=NG);
	}

	# verify that fixed params are within the imposed bounds
	if((!is.na(fixed_rholambda0)) && ((fixed_rholambda0<min_rholambda0) || (fixed_rholambda0>max_rholambda0))){
		return(list(success = FALSE, error=sprintf("Fixed rholambda0 (%g) is outside of the requested bounds (%g - %g)",fixed_rholambda0,min_rholambda0,max_rholambda0)));
	}
	if(any(fixed_PDR[!is.na(fixed_PDR)]<min_PDR[!is.na(fixed_PDR)]) || any(fixed_PDR[!is.na(fixed_PDR)]>max_PDR[!is.na(fixed_PDR)])){
		return(list(success = FALSE, error=sprintf("Some fixed PDRs are outside of the requested bounds")));
	}
						
	#################################
	# PREPARE PARAMETERS TO BE FITTED
	
	# guess reasonable start params, if not provided
	default_guess_PDR = mean(lineage_counter$relative_slopes); # a reasonable guesstimate for the average PDR is the average of the relative LTT-slope
	guess_PDR[is.na(guess_PDR)] = default_guess_PDR;
	if(is.na(guess_rholambda0)) guess_rholambda0 = tail(lineage_counter$relative_slopes,1);
	
	# make sure initial guess is within the imposed bounds
	guess_PDR = pmin(max_PDR, pmax(min_PDR, guess_PDR));
	guess_rholambda0 = min(max_rholambda0, max(min_rholambda0, guess_rholambda0))
	
	# determine which parameters are to be fitted
	# convention: parameters are indexed as follows: [PDR[], rholambda0]
	fixed_param_values 	= c(fixed_PDR, fixed_rholambda0); # may contain NAs, corresponding to non-fixed parameters
	fitted_params		= which(is.na(fixed_param_values))
	fixed_params		= which(!is.na(fixed_param_values))
	guess_param_values 	= c(guess_PDR, guess_rholambda0); # should contain a valid numeric for each parameter, even if the parameter is fixed
	guess_param_values[fixed_params] = fixed_param_values[fixed_params] # make sure guessed param values are consistent with fixed param values
	min_param_values	= c(min_PDR,min_rholambda0);
	max_param_values	= c(max_PDR,max_rholambda0);
	NFP					= length(fitted_params);
	
	# determine typical parameter scales
	scale_PDR = abs(guess_PDR); scale_PDR[scale_PDR==0] = mean(scale_PDR);
	scale_rholambda0 = abs(guess_rholambda0);
	if(scale_rholambda0==0) scale_rholambda0 = log2(Ntips)/root_age;
	param_scales = c(rep(scale_PDR,times=NG),scale_rholambda0);


	################################
	# FITTING
	
	# objective function: negated log-likelihood
	# input argument is the subset of fitted parameters, rescaled according to param_scales
	objective_function = function(fparam_values){
		param_values = fixed_param_values; param_values[fitted_params] = fparam_values * param_scales[fitted_params];
		if(any(is.nan(param_values)) || any(is.infinite(param_values))) return(Inf); # catch weird cases where params become NaN
		PDRs = param_values[1:NG]; 
		rholambda0 = param_values[NG+1];
		if(length(age_grid)==1){
			# while age-grid has only one point (i.e., PDRs are constant over time), we need to provide a least 2 grid points to the loglikelihood calculator, spanning the interval [0,oldest_age]
			input_age_grid 	= c(0,oldest_age);
			input_PDRs 		= c(PDRs, PDRs);
		}else{
			input_age_grid 	= age_grid;
			input_PDRs 		= PDRs
		}
		results = get_HBD_PDR_loglikelihood_CPP(branching_ages		= sorted_node_ages,
												oldest_age			= oldest_age,
												rholambda0 			= rholambda0,
												age_grid 			= input_age_grid,
												PDRs 				= input_PDRs,
												splines_degree		= splines_degree,
												condition			= condition,
												relative_dt			= relative_dt,
												runtime_out_seconds	= max_model_runtime);
		if(!results$success) return(Inf);
		LL = results$loglikelihood;
		if(is.na(LL) || is.nan(LL)) return(Inf);
		return(-LL);
	}
	

	# fit with various starting points
	fit_single_trial = function(trial){
		scales		 = param_scales[fitted_params]
		lower_bounds = min_param_values[fitted_params]
		upper_bounds = max_param_values[fitted_params]
		# randomly choose start values for fitted params
		start_values = guess_param_values[fitted_params]
		if(trial>1){
			boxed   = which(!(is.infinite(lower_bounds) || is.infinite(upper_bounds))); # determine fitted params that are boxed, i.e. constrained to within finite lower & upper bounds
			unboxed = completement(NFP, boxed);
			if(length(boxed)>0) start_values[boxed] = lower_bounds[boxed] + (upper_bounds[boxed]-lower_bounds[boxed]) * runif(n=length(boxed),min=0,max=1)
			if(length(unboxed)>0) start_values[unboxed]	= 10**runif(n=length(unboxed), min=-2, max=2) * start_values[unboxed]
		}
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
	loglikelihood		= objective_value
	fitted_param_values = fixed_param_values; fitted_param_values[fitted_params] = fits[[best]]$fparam_values;
	if(is.null(objective_value) || any(is.na(fitted_param_values)) || any(is.nan(fitted_param_values))) return(list(success=FALSE, error=sprintf("Some fitted parameters are NaN")));
		
	# return results
	return(list(success						= TRUE,
				objective_value				= objective_value,
				objective_name				= "loglikelihood",
				loglikelihood				= loglikelihood,
				fitted_PDR					= fitted_param_values[1:NG],
				fitted_rholambda0			= fitted_param_values[NG+1], 
				guess_PDR					= guess_param_values[1:NG],
				guess_rholambda0			= guess_param_values[NG+1],
				age_grid					= age_grid,
				NFP							= NFP,
				AIC							= 2*NFP - 2*loglikelihood,
				converged					= fits[[best]]$converged,
				Niterations					= fits[[best]]$Niterations,
				Nevaluations				= fits[[best]]$Nevaluations));
}



completement = function(N, indices){
	pool = rep(TRUE,N);
	pool[indices] = FALSE;
	return(which(pool));
}


