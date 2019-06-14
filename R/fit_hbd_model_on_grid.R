# Fit a homogenous-birth-death cladogenic model to an ultrametric timetree, by estimating lambda & mu at discrete ages
# An HBD model is defined by a time-dependent speciation rate (lambda), a time-dependent extinction rate (mu) and a rarefaction (rho, subsampling fraction)
# The speciation rate lambda and extinction rate mu are specified on a discrete age-grid, and assumed to vary linearly (or polynomially, as splines) between grid points (see "degree" argument).
#
# Note that for each specific model and a given timetree there exists a continuum of alternative models that would all generate the same deterministic lineages-through-time (LTT) curve (when calculated backward in time), and all of these models actually have the same likelihood.
#   Hence, each model is part of an "equivalence class" of models, and likelihood-based approaches can only discern between model classes, but not between the individual model members in a class
#   It turns out that each HBD model-class is uniquely defined by its "pulled diversification rate" (lambda) and the product rho*lambda(0)=:rho.
#   You should thus seriously consider using the function fit_HBD_class() instead, unless you can a priori strongly constrain the speciation or extinction rates, as well as the rarefaction.
#
# References:
#	Morlon et al. (2011). Reconciling molecular phylogenies with the fossil record. PNAS 108:16327-16332
fit_hbd_model_on_grid = function(	tree, 
									oldest_age			= NULL,		# either a numeric specifying the stem age or NULL (equivalent to the root age). This is similar to the "tot_time" option in the R function RPANDA::likelihood_bd
									age_grid			= NULL,		# either NULL, or a numeric vector of size NG, listing ages in ascending order, on which the lambda is defined as a piecewise linear curve. If NULL, the lambda is assumed to be time-independent.
									min_lambda			= 0,		# optional lower bound for the fitted lambdas. Either a single numeric (applying to all age-grid-points) or a numeric vector of size NG, specifying the lower bound at each age-grid point.
									max_lambda			= +Inf,		# optional upper bound for the fitted lambdas. Either a single numeric (applying to all age-grid-points) or a numeric vector of size NG, specifying the upper bound at each age-grid point.
									min_mu				= 0,		# optional lower bound for the fitted mus. Either a single numeric (applying to all age-grid-points) or a numeric vector of size NG, specifying the lower bound at each age-grid point.
									max_mu				= +Inf,		# optional upper bound for the fitted mus. Either a single numeric (applying to all age-grid-points) or a numeric vector of size NG, specifying the upper bound at each age-grid point.
									min_rho				= 1e-10,	# optional lower bound for the fitted rho. Note that rho is always within (0,1]
									max_rho				= 1,		# optional upper bound for the fitted rho
									guess_lambda		= NULL,		# initial guess for the lambda. Either NULL (an initial guess will be computed automatically), or a single numeric (guessing a constant lambda at all ages) or a numeric vector of size NG specifying an initial guess for the lambda at each age-grid point (can include NAs)
									guess_mu			= NULL,		# initial guess for the mu. Either NULL (an initial guess will be computed automatically), or a single numeric (guessing a constant mu at all ages) or a numeric vector of size NG specifying an initial guess for the mu at each age-grid point (can include NAs)
									guess_rho			= 1,		# initial guess for the product rho*lambda(0). Either NULL (an initial guess will be computed automatically) or a single strictly-positive numeric.
									fixed_lambda		= NULL,		# optional fixed lambda values, on one or more of the age grid points. Either NULL (none of the lambdas are fixed), or a single scalar (all lambdas are fixed) or a numeric vector of size NG (some or all lambdas are fixed, can include NAs).
									fixed_mu			= NULL,		# optional fixed mu values, on one or more of the age grid points. Either NULL (none of the mus are fixed), or a single scalar (all mus are fixed) or a numeric vector of size NG (some or all mus are fixed, can include NAs).
									fixed_rho			= NULL,		# optional fixed value for rho. If non-NULL and non-NA, then rho is not fitted. 
									const_lambda		= FALSE,	# logical, whether to enforce a constant (time-independent) fitted speciation rate. Only relevant for those lambdas that are fitted (i.e. fixed lambda values are kept as is).
									const_mu			= FALSE,	# logical, whether to enforce a constant (time-independent) fitted extinction rate. Only relevant for those lambdas that are fitted (i.e. fixed lambda values are kept as is).
									splines_degree		= 1,		# integer, either 1 or 2 or 3, specifying the degree for the splines defined by lambda and mu on the age grid.
									condition			= "stem",	# one of "crown" or "stem", specifying whether to condition the likelihood on the survival of the stem group or the crown group. It is recommended to use "stem" when oldest_age>root_age, and "crown" when oldest_age==root_age. This argument is similar to the "cond" argument in the R function RPANDA::likelihood_bd. Note that "crown" really only makes sense when oldest_age==root_age.
									relative_dt			= 1e-3,		# maximum relative time step allowed for integration. Smaller values increase the accuracy of the computed likelihoods, but increase computation time. Typical values are 0.0001-0.001. The default is usually sufficient.
									Ntrials				= 1,
									Nthreads			= 1,
									max_model_runtime	= NULL,		# maximum time (in seconds) to allocate for each likelihood evaluation. Use this to escape from badly parameterized models during fitting (this will likely cause the affected fitting trial to fail). If NULL or <=0, this option is ignored.
									fit_control			= list()){	# a named list containing options for the nlminb fitting routine (e.g. iter.max and rel.tol)
	Ntips	= length(tree$tip.label);
	Nnodes	= tree$Nnode;

	# pre-compute some tree stats
	lineage_counter = count_lineages_through_time(tree, Ntimes=log2(Ntips), include_slopes=TRUE);
	sorted_node_ages = sort(get_all_branching_ages(tree));
	root_age 		 = tail(sorted_node_ages,1)
	age_epsilon		 = 1e-4*mean(tree$edge.length);

	# basic input checking
	if((Ntips<2) || (Nnodes<2)){
		# tree is trivial (~empty)
		return(list(success = FALSE, error="Tree is too small"));
	}
	if(!(condition %in% c("stem","crown","none"))) return(list(success = FALSE, error = sprintf("Invalid condition option '%s'; expected either 'crown', or 'stem' or 'none' (not recommended)",condition)))
	if(is.null(fixed_rho)) fixed_rho = NA;
	if((!is.na(fixed_rho)) && ((fixed_rho<=0) || (fixed_rho>1))) return(list(success = FALSE, error=sprintf("Fixed rho (%g) is outside of the accepted range (0,1].",fixed_rho)));
	if(is.null(oldest_age)) oldest_age = root_age;
	if(Ntrials<1) return(list(success = FALSE, error = sprintf("Ntrials must be at least 1")))
	if(is.null(age_grid) || (length(age_grid)<=1)){
		if((!is.null(guess_lambda)) && (length(guess_lambda)>1)) return(list(success = FALSE, error = sprintf("Invalid number of guessed lambdas; since no age grid was provided, you must provide a single (constant) guess_lambda or none at all")));
		if((!is.null(guess_mu)) && (length(guess_mu)>1)) return(list(success = FALSE, error = sprintf("Invalid number of guessed mus; since no age grid was provided, you must provide a single (constant) guess_mu or none at all")));
		age_grid = 0 # single-point grid, means that lambdas are assumed time-independent
		NG = 1
	}else{
		NG = length(age_grid)
		if((!is.null(guess_lambda)) && (length(guess_lambda)!=1) && (length(guess_lambda)!=NG)) return(list(success = FALSE, error = sprintf("Invalid number of guessed lambdas (%d); since an age grid of size %d was provided, you must either provide one or %d lambdas",length(guess_lambda),NG)));
		if((!is.null(guess_mu)) && (length(guess_mu)!=1) && (length(guess_mu)!=NG)) return(list(success = FALSE, error = sprintf("Invalid number of guessed mus (%d); since an age grid of size %d was provided, you must either provide one or %d mus",length(guess_mu),NG)));
		if((length(age_grid)>1) && (age_grid[NG]>oldest_age-1e-5*(age_grid[NG]-age_grid[NG-1]))) age_grid[NG] = max(age_grid[NG],oldest_age); # if age_grid "almost" covers oldest_age (i.e. up to rounding errors), then fix the remaining difference
		if((length(age_grid)>1) && ((age_grid[1]>0) || (age_grid[NG]<oldest_age))) return(list(success = FALSE, error=sprintf("Provided age-grid range (%g - %g) does not cover entire required age range (0 - %g)",age_grid[1],tail(age_grid,1),oldest_age)));
	}
	if(is.null(max_model_runtime)) max_model_runtime = 0;
	if(!(splines_degree %in% c(0,1,2,3))) return(list(success = FALSE, error = sprintf("Invalid splines_degree: Extected one of 0,1,2,3.")));
	if(NG==1) splines_degree = 1; # no point in using splines since lambda & mu are assumed to be time-independent
	
	# reformat shape of input params to an internally standardized format
	min_rho = max(0,min_rho);
	max_rho = max(0,max_rho);
	if(length(min_lambda)==1) min_lambda = rep(min_lambda,times=NG);
	if(length(max_lambda)==1) max_lambda = rep(max_lambda,times=NG);
	if(length(min_mu)==1) min_mu = rep(min_mu,times=NG);
	if(length(max_mu)==1) max_mu = rep(max_mu,times=NG);
	if(is.null(guess_rho)) guess_rho = NA;
	if(is.null(guess_lambda)){
		guess_lambda = rep(NA,times=NG);
	}else if(length(guess_lambda)==1){
		guess_lambda = rep(guess_lambda,times=NG);
	}
	if(is.null(fixed_lambda)){
		fixed_lambda = rep(NA,times=NG);
	}else if(length(fixed_lambda)==1){
		fixed_lambda = rep(fixed_lambda,times=NG);
	}
	if(is.null(guess_mu)){
		guess_mu = rep(NA,times=NG);
	}else if(length(guess_mu)==1){
		guess_mu = rep(guess_mu,times=NG);
	}
	if(is.null(fixed_mu)){
		fixed_mu = rep(NA,times=NG);
	}else if(length(fixed_mu)==1){
		fixed_mu = rep(fixed_mu,times=NG);
	}

	# verify that fixed params are within the imposed bounds
# 	if(any(fixed_lambda[!is.na(fixed_lambda)]<min_lambda[!is.na(fixed_lambda)]) || any(fixed_lambda[!is.na(fixed_lambda)]>max_lambda[!is.na(fixed_lambda)])){
# 		return(list(success = FALSE, error=sprintf("Some fixed lambdas are outside of the requested bounds")));
# 	}
# 	if(any(fixed_mu[!is.na(fixed_mu)]<min_mu[!is.na(fixed_mu)]) || any(fixed_mu[!is.na(fixed_mu)]>max_mu[!is.na(fixed_mu)])){
# 		return(list(success = FALSE, error=sprintf("Some fixed mus are outside of the requested bounds")));
# 	}
# 	if((!is.na(fixed_rho)) && ((fixed_rho<min_rho) || (fixed_rho>max_rho))){
# 		return(list(success = FALSE, error=sprintf("Fixed rho (%g) is outside of the requested bounds (%g - %g)",fixed_rho,min_rho,max_rho)));
# 	}
	
						
	#################################
	# PREPARE PARAMETERS TO BE FITTED
	
	# guess reasonable start params, if not provided
	if(is.na(guess_rho)) guess_rho = 1;
	default_guess_PDR 	 				= mean(lineage_counter$relative_slopes); # a reasonable guesstimate for the average PDR is the average of the relative LTT-slope
	default_guess_lambda 				= tail(lineage_counter$relative_slopes,1); # a reasonable guesstimate for the present-day lambda is the relative LTT-slope at age=0
	guess_lambda[is.na(guess_lambda)] 	= default_guess_lambda;
	default_guess_mu 					= max(0,default_guess_lambda - default_guess_PDR);
	guess_mu[is.na(guess_mu)] 			= default_guess_mu;
		
	# make sure initial guess is within the imposed bounds
	guess_lambda = pmin(max_lambda, pmax(min_lambda, guess_lambda));
	guess_mu	 = pmin(max_mu, pmax(min_mu, guess_mu));
	guess_rho 	 = min(max_rho, max(min_rho, guess_rho))
	
	# determine which parameters are to be fitted
	# convention: parameters are indexed as follows: [lambda[], mu[], rho]
	fixed_param_values 	= c(fixed_lambda, fixed_mu, fixed_rho); # may contain NAs, corresponding to non-fixed parameters
	fitted_params		= which(is.na(fixed_param_values))
	fixed_params		= which(!is.na(fixed_param_values))
	guess_param_values 	= c(guess_lambda, guess_mu, guess_rho); # should contain a valid numeric for each parameter, even if the parameter is fixed
	guess_param_values[fixed_params] = fixed_param_values[fixed_params] # make sure guessed param values are consistent with fixed param values
	min_param_values	= c(min_lambda,min_mu,min_rho);
	max_param_values	= c(max_lambda,max_mu,max_rho);

	# determine free (i.e. independent) fitted parameters
	# for example, if lambda is enforced to be time-independent, this reduces the number of free parameters
	# free2fitted[frp] (where frp=1,..,Nfree) will be a list of fitted parameter indices represented by the frp-th free parameter
	# fitted2free[fp] will be the index of the free parameter representing the fp-th fitted parameter
	NFlambda		= sum(is.na(fixed_lambda)) 	# number of non-fixed lambda
	NFmu			= sum(is.na(fixed_mu))		# number of non-fixed mu
	fitted2free		= (if(NFlambda==0) c() else (if(const_lambda) rep(1,times=NFlambda) else c(1:NFlambda)))
	fitted2free		= c(fitted2free, (if(NFmu==0) c() else tail(fitted2free,1)+(if(const_mu) rep(1,times=NFmu) else c(1:NFmu))))
	fitted2free		= c(fitted2free, (if(is.na(fixed_rho)) tail(fitted2free,1)+1 else c()))
	Nfree			= length(unique(fitted2free)); # number of free (i.e. independently) fitted parameters
	free2fitted		= lapply(1:Nfree, FUN=function(frp) which(fitted2free==frp))
	
	# determine typical parameter scales
	scale_lambda = abs(guess_lambda); scale_lambda[scale_lambda==0] = mean(scale_lambda);
	scale_mu 	 = abs(guess_mu);
	if(all(scale_mu==0)){ 
		scale_mu[] = scale_lambda; 
	}else{ 
		scale_mu[scale_mu==0] = mean(scale_mu); 
	}
	scale_rho 	 = abs(guess_rho);
	if(scale_rho==0) scale_rho = 1;
	param_scales = c(rep(scale_lambda,times=NG),rep(scale_mu,times=NG),scale_rho);
	
	# define auxiliary function for obtaining full parameter list from rescaled free fitted parameters
	# input: fparam_values[] is a 1D vector of length NFP, listing rescaled values for the free fitted parameters
	# output: param_values[] will be a 1D vector of length NP, listing all model parameter values
	fparam_scales = sapply(1:Nfree, FUN = function(frp) mean(param_scales[fitted_params[free2fitted[[frp]]]]))
	expand_free_fitted_params = function(fparam_values){
		fparam_values = fparam_values * fparam_scales;
		param_values = fixed_param_values; 
		param_values[fitted_params] = pmax(min_param_values[fitted_params], pmin(max_param_values[fitted_params], fparam_values[fitted2free]))
		return(param_values)
	}

	

	################################
	# FITTING
	
	# objective function: negated log-likelihood
	# input argument is the subset of fitted parameters, rescaled according to param_scales
	objective_function = function(fparam_values){
		param_values = expand_free_fitted_params(fparam_values);
		if(any(is.nan(param_values)) || any(is.infinite(param_values))) return(Inf); # catch weird cases where params become NaN
		lambdas = param_values[1:NG]; 
		mus 	= param_values[(NG+1):(NG+NG)]; 
		rho 	= param_values[2*NG+1];
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
													splines_degree		= splines_degree,
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
		lower_bounds = sapply(1:Nfree, FUN = function(frp) max(min_param_values[fitted_params[free2fitted[[frp]]]]))
		upper_bounds = sapply(1:Nfree, FUN = function(frp) min(max_param_values[fitted_params[free2fitted[[frp]]]]))
		# randomly choose start values for fitted params
		start_values = sapply(1:Nfree, FUN = function(frp) guess_param_values[fitted_params[free2fitted[[frp]][1]]])
		if(trial>1){
			boxed   = which(!(is.infinite(lower_bounds) || is.infinite(upper_bounds))); # determine fitted params that are boxed, i.e. constrained to within finite lower & upper bounds
			unboxed = complement(Nfree, boxed);
			if(length(boxed)>0) start_values[boxed] = lower_bounds[boxed] + (upper_bounds[boxed]-lower_bounds[boxed]) * runif(n=length(boxed),min=0,max=1)
			if(length(unboxed)>0) start_values[unboxed]	= 10**runif(n=length(unboxed), min=-2, max=2) * start_values[unboxed]
		}
		start_values = pmax(lower_bounds,pmin(upper_bounds,start_values))
		# run fit
		fit = stats::nlminb(start_values/fparam_scales, 
							objective	= objective_function, 
							lower		= lower_bounds/fparam_scales, 
							upper		= upper_bounds/fparam_scales, 
							control		= fit_control)
		return(list(objective_value=fit$objective, fparam_values = fit$par, converged=(fit$convergence==0), Niterations=fit$iterations, Nevaluations=fit$evaluations[1]));
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
	fitted_param_values = expand_free_fitted_params(fits[[best]]$fparam_values)
	if(is.null(objective_value) || any(is.na(fitted_param_values)) || any(is.nan(fitted_param_values))) return(list(success=FALSE, error=sprintf("Some fitted parameters are NaN")));
		
	# return results
	return(list(success					= TRUE,
				objective_value			= objective_value,
				objective_name			= "loglikelihood",
				loglikelihood			= loglikelihood,
				fitted_lambda			= fitted_param_values[1:NG],
				fitted_mu				= fitted_param_values[(NG+1):(NG+NG)],
				fitted_rho				= fitted_param_values[2*NG+1], 
				guess_lambda			= guess_param_values[1:NG],
				guess_mu				= guess_param_values[(NG+1):(NG+NG)],
				guess_rho				= guess_param_values[2*NG+1],
				age_grid				= age_grid,
				NFP						= Nfree,
				AIC						= 2*Nfree - 2*loglikelihood,
				converged				= fits[[best]]$converged,
				Niterations				= fits[[best]]$Niterations,
				Nevaluations			= fits[[best]]$Nevaluations));
}



complement = function(N, indices){
	pool = rep(TRUE,N);
	pool[indices] = FALSE;
	return(which(pool));
}

