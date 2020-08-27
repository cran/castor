# Fit a homogenous-birth-death cladogenic model-congruence-class to an ultrametric timetree, 
#   by fitting the pulled diversification rate (PDR) on a discrete age grid, where the grid point locations are themselves fitted
#
# This function can optionally fit the congruence class by considering only the age interval [age0, oldest_age].
# oldest_age is enforced by (mathematically) cutting the tree into multiple subtrees stemming at age oldest_age (omitting anything older) and treating each subtree as an independent realization of the same stochastic BD process
# age0 is enforced by (actually) trimming the tree at age0 (omitting anything younger than age0) and fitting the HBD class to the new (shorter) timetree. 
# In that case, the fitted rholambda0 is actually the product Phi(age0)*lambda(age0), where Phi(age0) is the probability that a lineage extant at age tau would be eventually included in the original timetree; Phi(age0) is essentially the sampling fraction at age0.
#
# References:
#	Morlon et al. (2011). Reconciling molecular phylogenies with the fossil record. PNAS 108:16327-16332
fit_hbd_pdr_and_grid = function(tree, 
								oldest_age				= NULL,		# either a numeric specifying the stem age or NULL (equivalent to the root age). This is similar to the "tot_time" option in the R function RPANDA::likelihood_bd
								age0					= 0,		# non-negative numeric, youngest age (time before present) to consider when fitting and with respect to which rholambda0 is defined (i.e. rholambda0 = rho(age0)*lambda(age0))
								Ngrid					= 1,		# integer, the number of age grid points to consider. The grid will always span from oldest_age to age0. If Ngrid=1, this means a single (time-independent) PDR will be fitted throughout.
								min_age_grid			= NULL,		# optional lower bound for the fitted grid ages. Either NULL, or a numeric vector of length Ngrid. Only revelant of Ngrid>=3.
								max_age_grid			= NULL,		# optional upper bound for the fitted grid ages. Either NULL, or a numeric vector of length Ngrid. Only revelant of Ngrid>=3.
								min_PDR					= -Inf,		# optional lower bound for the fitted PDRs. Either a single numeric (applying to all age-grid-points) or a numeric vector of size Ngrid, specifying the lower bound at each age-grid point.
								max_PDR					= +Inf,		# optional upper bound for the fitted PDRs. Either a single numeric (applying to all age-grid-points) or a numeric vector of size Ngrid, specifying the upper bound at each age-grid point.
								min_rholambda0			= 1e-10,	# optional lower bound for the fitted rholambda0. Note that rholambda0 is always non-negative.
								max_rholambda0			= +Inf,		# optional upper bound for the fitted rholambda0
								guess_age_grid			= NULL,		# initial guess for the age grid. Either NULL (an initial guess will be computed automatically) or a numeric vector of length Ngrid (if Ngrid>=2), listing discrete ages in ascending order. Only relevant if Ngrid>=3.
								guess_PDR				= NULL,		# initial guess for the PDR. Either NULL (an initial guess will be computed automatically), or a single numeric (guessing a constant PDR at all ages) or a numeric vector of size Ngridspecifying an initial guess for the PDR at each age-grid point (can include NAs)
								guess_rholambda0		= NULL,		# initial guess for the product rho*lambda(0). Either NULL (an initial guess will be computed automatically) or a single strictly-positive numeric.
								fixed_age_grid			= NULL,		# optional fixed grid ages, for one or more of the grid points. Either NULL (none of the grid ages are fixed), or a numeric vector of size Ngrid (some or all grid ages are fixed, can include NAs). The first and last elements are ignored, since these are always fixed at age0 and oldest_age, respectively. Only relevant if Ngrid>=3. fixed_age_grid[n] corresponds to the n-th grid age, when sorted in increasing order.
								fixed_PDR				= NULL,		# optional fixed PDR values, on one or more of the age grid points. Either NULL (none of the PDRs are fixed), or a single scalar (all PDRs are fixed) or a numeric vector of size Ngrid(some or all PDRs are fixed, can include NAs).
								fixed_rholambda0		= NULL,		# optional fixed value for rholambda0. If non-NULL and non-NA, then rholambda0 is not fitted. 
								min_age_step			= 0,		# optional minimum allowed distance between adjacent grid ages.
								splines_degree			= 1,		# integer, either 1 or 2 or 3, specifying the degree for the splines defined by the PDR on the age grid.
								condition				= "auto",	# one of "crown" or "stem" or "auto", specifying whether to condition the likelihood on the survival of the stem group or the crown group. It is recommended to use "stem" when oldest_age<root_age, and "crown" when oldest_age==root_age. This argument is similar to the "cond" argument in the R function RPANDA::likelihood_bd. Note that "crown" really only makes sense when oldest_age==root_age.
								relative_dt				= 1e-3,		# maximum relative time step allowed for integration. Smaller values increase the accuracy of the computed likelihoods, but increase computation time. Typical values are 0.0001-0.001. The default is usually sufficient.
								Ntrials					= 1,
								Nbootstraps				= 0,		# (integer) optional number of parametric-bootstrap samples (random trees generated using the fitted PDR & rholambda0) for estimating confidence intervals of fitted parameters. If 0, no parametric bootstrapping is performed. Typical values are 10-100.
								Ntrials_per_bootstrap	= NULL,		# (integer) optional number of fitting trials for each bootstrap sampling. If NULL, this is set equal to Ntrials. A smaller Ntrials_per_bootstrap will reduce computation, at the expense of increasing the estimated confidence intervals (i.e. yielding more conservative estimates of confidence).
								Nthreads				= 1,
								max_model_runtime		= NULL,		# maximum time (in seconds) to allocate for each likelihood evaluation. Use this to escape from badly parameterized models during fitting (this will likely cause the affected fitting trial to fail). If NULL or <=0, this option is ignored.
								fit_control				= list(),	# a named list containing options for the optim and constrOptim fitting routines (e.g. maxit, reltol and max_outer_iterations). Note that every outer iteration is a full optim() call, and maxit refers to the total number of iterations withi the optim() algorithm, counted across all outer iterations.
								verbose					= FALSE,	# boolean, specifying whether to print informative messages
								verbose_prefix			= ""){		# string, specifying the line prefix when printing messages. Only relevant if verbose==TRUE.	
	# basic error checking
	if(verbose) cat(sprintf("%sChecking input variables..\n",verbose_prefix))
	original_Ntips = length(tree$tip.label)
	if(tree$Nnode<2) return(list(success = FALSE, error="Input tree is too small"));
	if(age0<0) return(list(success = FALSE, error="age0 must be non-negative"));
	min_age_step = max(0.0,min_age_step)
	if(Ntrials<1) return(list(success = FALSE, error = sprintf("Ntrials must be at least 1")))
	if(!(splines_degree %in% c(0,1,2,3))) return(list(success = FALSE, error = sprintf("Invalid splines_degree: Extected one of 0,1,2,3.")));
	if(is.null(max_model_runtime)) max_model_runtime = 0;
	root_age = get_tree_span(tree)$max_distance
	if(is.null(oldest_age)) oldest_age = root_age;
	if(!(condition %in% c("crown","stem","auto"))) return(list(success = FALSE, error = sprintf("Invalid condition '%s': Extected 'stem', 'crown' or 'auto'.",condition)));
	if(condition=="auto") condition = (if(abs(oldest_age-root_age)<=1e-10*root_age) "crown" else "stem")
	if(root_age<age0) return(list(success=FALSE, error=sprintf("age0 (%g) is older than the root age (%g)",age0,root_age)));
	if(oldest_age<age0) return(list(success=FALSE, error=sprintf("age0 (%g) is older than the oldest considered age (%g)",age0,oldest_age)));
	if(Ngrid<1) return(list(success=FALSE, error=sprintf("Ngrid is %d, but must be >=1",Ngrid)));
	if(Ngrid==1){
		splines_degree 	= 1; # no point in using splines since PDR is assumed to be time-independent
		guess_age_grid 	= age0
		min_age_grid 	= age0
		max_age_grid 	= age0
		fixed_age_grid 	= age0
	}
	if((!is.null(guess_age_grid)) && (length(guess_age_grid)!=Ngrid)) return(list(success=FALSE, error=sprintf("guess_age_grid must either be NULL or a numeric vector of length Ngrid (%d)",Ngrid)));
	if((!is.null(min_age_grid)) && (length(min_age_grid)!=Ngrid)) return(list(success=FALSE, error=sprintf("min_age_grid must either be NULL or a numeric vector of length Ngrid (%d)",Ngrid)));
	if((!is.null(max_age_grid)) && (length(max_age_grid)!=Ngrid)) return(list(success=FALSE, error=sprintf("max_age_grid must either be NULL or a numeric vector of length Ngrid (%d)",Ngrid)));
	if((!is.null(fixed_age_grid)) && (length(fixed_age_grid)!=Ngrid)) return(list(success=FALSE, error=sprintf("fixed_age_grid must either be NULL or a numeric vector of length Ngrid (%d)",Ngrid)));
	if((!is.null(guess_PDR)) && (length(guess_PDR)!=1) && (length(guess_PDR)!=Ngrid)) return(list(success = FALSE, error = sprintf("Invalid number of guessed PDRs (%d); since an age grid of size %d was provided, you must either provide one or %d PDRs",length(guess_PDR),Ngrid)));
	if((!is.null(fixed_age_grid)) && (Ngrid>1) && (!is.na(fixed_age_grid[1])) && (fixed_age_grid[1]!=age0)) return(list(success = FALSE, error=sprintf("The first grid age must be age0 (%g); instead, it was set to %g",age0,fixed_age_grid[1])))
	if((!is.null(fixed_age_grid)) && (Ngrid>1) && (!is.na(fixed_age_grid[Ngrid])) && (fixed_age_grid[Ngrid]!=oldest_age)) return(list(success = FALSE, error=sprintf("The last grid age must be oldest_age (%g); instead, it was set to %g",oldest_age,fixed_age_grid[Ngrid])))
	if((Ngrid>1) && (min_age_step>(oldest_age-age0)/(Ngrid-1))) return(list(success = FALSE, error = sprintf("min_age_step (%g) is too large for this number of grid points (%d) and age span (%g)",min_age_step,Ngrid,oldest_age-age0)))

	# save some of the input arguments for potential bootstrapping
	original_min_age_grid		= min_age_grid
	original_max_age_grid		= max_age_grid
	original_min_PDR			= min_PDR
	original_max_PDR			= max_PDR
	original_min_rholambda0		= min_rholambda0
	original_max_rholambda0		= max_rholambda0
	original_guess_age_grid 	= guess_age_grid
	original_guess_PDR			= guess_PDR
	original_guess_rholambda0 	= guess_rholambda0

	# trim tree at age0 if needed, while shifting time for the subsequent analyses (i.e. new ages will start counting at age0)
	if(age0>0){
		if(verbose) cat(sprintf("%sTrimming tree at age0=%g..\n",verbose_prefix,age0))
		tree = trim_tree_at_height(tree,height=root_age-age0)$tree
		if(tree$Nnode<2) return(list(success = FALSE, error=sprintf("Tree is too small after trimming at age0 (%g)",age0)));
		if(!is.null(oldest_age)) 		oldest_age		= oldest_age - age0	
		if(!is.null(min_age_grid)) 		min_age_grid 	= min_age_grid - age0
		if(!is.null(max_age_grid)) 		max_age_grid 	= max_age_grid - age0
		if(!is.null(guess_age_grid)) 	guess_age_grid 	= guess_age_grid - age0
		if(!is.null(fixed_age_grid)) 	fixed_age_grid 	= fixed_age_grid - age0
		root_age = root_age - age0
	}
								
	# pre-compute some tree stats
	if(verbose) cat(sprintf("%sPrecomputing some stats about the tree..\n",verbose_prefix))
	LTT0				= length(tree$tip.label);
	sorted_node_ages	= sort(get_all_branching_ages(tree));
	root_age 		 	= tail(sorted_node_ages,1);
	age_epsilon		 	= 1e-4*mean(tree$edge.length);
	lineage_counter 	= count_lineages_through_time(tree, Ntimes=log2(LTT0), include_slopes=TRUE);
	lineage_counter$ages = root_age - lineage_counter$times
	
	# reformat shape of input params to an internally standardized format
	if(verbose) cat(sprintf("%sPreparing for fitting..\n",verbose_prefix))
	if(is.null(min_age_grid)) min_age_grid 		= rep(NA, times=Ngrid)
	if(is.null(max_age_grid)) max_age_grid 		= rep(NA, times=Ngrid)
	if(is.null(guess_age_grid)) guess_age_grid 	= rep(NA, times=Ngrid)
	if(is.null(fixed_age_grid)) fixed_age_grid 	= rep(NA, times=Ngrid)
	min_rholambda0 = max(0,min_rholambda0);
	max_rholambda0 = max(0,max_rholambda0);
	if(length(min_PDR)==1) min_PDR = rep(min_PDR,times=Ngrid);
	if(length(max_PDR)==1) max_PDR = rep(max_PDR,times=Ngrid);
	if(is.null(guess_rholambda0)) guess_rholambda0 = NA;
	if(is.null(fixed_rholambda0)) fixed_rholambda0 = NA;
	if(is.null(guess_PDR)){
		guess_PDR = rep(NA,times=Ngrid);
	}else if(length(guess_PDR)==1){
		guess_PDR = rep(guess_PDR,times=Ngrid);
	}
	if(is.null(fixed_PDR)){
		fixed_PDR = rep(NA,times=Ngrid);
	}else if(length(fixed_PDR)==1){
		fixed_PDR = rep(fixed_PDR,times=Ngrid);
	}

	# verify that min-max make sense, and that fixed params are within the imposed bounds
	if((!is.na(fixed_rholambda0)) && ((fixed_rholambda0<min_rholambda0) || (fixed_rholambda0>max_rholambda0))){
		return(list(success = FALSE, error=sprintf("Fixed rholambda0 (%g) is outside of the requested bounds (%g - %g)",fixed_rholambda0,min_rholambda0,max_rholambda0)));
	}
	if(any(fixed_PDR[!is.na(fixed_PDR)]<min_PDR[!is.na(fixed_PDR)]) || any(fixed_PDR[!is.na(fixed_PDR)]>max_PDR[!is.na(fixed_PDR)])){
		return(list(success = FALSE, error=sprintf("Some fixed PDRs are outside of the requested bounds")));
	}
	if(any(fixed_age_grid[!is.na(fixed_age_grid)]<age0) || any(fixed_age_grid[!is.na(fixed_age_grid)]>oldest_age)){
		return(list(success = FALSE, error=sprintf("Some fixed grid ages are outside of the requested bounds")));
	}else if(sum(!is.na(fixed_age_grid))>=2){
		if(any(diff(fixed_age_grid[!is.na(fixed_age_grid)])<=0)){
			return(list(success = FALSE, error=sprintf("Fixed grid ages must be strictly increasing")));
		}else if(any(diff(fixed_age_grid[!is.na(fixed_age_grid)])<min_age_step)){
			return(list(success = FALSE, error=sprintf("Some adjacent fixed grid ages are below the minimum allowed age step (%g)",min_age_step)));
		}
	}
	if((sum(!is.na(min_age_grid))>=2) && any(diff(min_age_grid[!is.na(min_age_grid)])<0)){
		return(list(success = FALSE, error=sprintf("Min grid ages must be non-decreasing")));
	}
	if((sum(!is.na(max_age_grid))>=2) && any(diff(max_age_grid[!is.na(max_age_grid)])<0)){
		return(list(success = FALSE, error=sprintf("Max grid ages must be non-decreasing")));
	}
	if(any((!is.na(min_age_grid)) & (!is.na(max_age_grid)) & (max_age_grid<min_age_grid))){
		return(list(success = FALSE, error=sprintf("Some grid age minima are greater than their maxima")));
	}
	if(any((!is.na(min_PDR)) & (!is.na(max_PDR)) & (max_PDR<min_PDR))){
		return(list(success = FALSE, error=sprintf("Some PDR minima are greater than their maxima")));
	}
	if((!is.na(min_rholambda0)) && (!is.na(max_rholambda0)) && (max_rholambda0<min_rholambda0)){
		return(list(success = FALSE, error=sprintf("rholambda0 minimum is greater than the rholambda0 maximum")));
	}
	
	# make sure guess_age_grid[] and fixed_age_grid[] are consistent in terms of order
	composite_age_grid = guess_age_grid; composite_age_grid[!is.na(fixed_age_grid)] = fixed_age_grid[!is.na(fixed_age_grid)];
	if(!all(is.na(composite_age_grid))){
		if(any(diff(composite_age_grid[!is.na(composite_age_grid)])<0)){
			return(list(success = FALSE, error=sprintf("Provided guess grid ages and fixed grid ages are inconsistent; their combination must yield a list of grid ages in ascending order")));
		}
	}
		
	# first and last grid points are always fixed
	fixed_age_grid[1] 		= 0; 
	fixed_age_grid[Ngrid] 	= oldest_age;
	
	# set min & max grid ages
	max_age_grid[!is.na(max_age_grid)] = pmin(oldest_age, max_age_grid[!is.na(max_age_grid)])
	min_age_grid[!is.na(min_age_grid)] = pmax(0, min_age_grid[!is.na(min_age_grid)])
	min_age_grid[!is.na(fixed_age_grid)] = fixed_age_grid[!is.na(fixed_age_grid)]
	max_age_grid[!is.na(fixed_age_grid)] = fixed_age_grid[!is.na(fixed_age_grid)]
	if(Ngrid>1){
		lower = 0
		for(g in 1:Ngrid){
			if(is.na(min_age_grid[g])){
				min_age_grid[g] = lower
			}else{
				min_age_grid[g] = max(lower,min_age_grid[g])
			}
			lower = max(lower,min_age_grid[g]+min_age_step)
		}
		upper = oldest_age
		for(g in Ngrid:1){
			if(is.na(max_age_grid[g])){
				max_age_grid[g] = upper
			}else{
				max_age_grid[g] = min(upper,max_age_grid[g])
			}
			upper = min(upper,max_age_grid[g]-min_age_step)
		}
	}
	
	# if some min-max constrain a parameter to a single value, fix that parameter (because the optimization algorithm fails otherwise)
	fixed_age_grid 		= fix_constrained_params_if_needed(min_age_grid, max_age_grid, fixed_age_grid)
	fixed_PDR 			= fix_constrained_params_if_needed(min_PDR, max_PDR, fixed_PDR)
	fixed_rholambda0 	= fix_constrained_params_if_needed(min_rholambda0, max_rholambda0, fixed_rholambda0)
			
						
	#################################
	# PREPARE PARAMETERS TO BE FITTED
	
	# guess reasonable start params, if not provided
	if(Ngrid>1){
		guess_age_grid[!is.na(fixed_age_grid)] = fixed_age_grid[!is.na(fixed_age_grid)]
		if(all(is.na(guess_age_grid)) && all(is.na(min_age_grid)) && all(is.na(max_age_grid))){
			guess_age_grid = get_inhomogeneous_grid_1D(Xstart=age0, Xend=oldest_age, Ngrid=Ngrid, densityX=rev(lineage_counter$ages), densityY=rev(sqrt(lineage_counter$lineages)), extrapolate=TRUE)
		}else if(any(is.na(guess_age_grid))){
			guess_age_grid[1] = 0
			guess_age_grid[Ngrid] = oldest_age
			last_nonna = 1
			for(g in 2:Ngrid){
				if(!is.na(guess_age_grid[g])){
					if(last_nonna<g-1){
						# we passed by one or more NAs, so set default values for those NAs as regular sub-grid
						guess_age_grid[(last_nonna+1):(g-1)] = seq(from=guess_age_grid[last_nonna],to=guess_age_grid[g],length.out=g-last_nonna+1)[2:(g-last_nonna)]
					}
					last_nonna = g;
				}
			}
		}
	}
	default_guess_PDR = mean(lineage_counter$relative_slopes); # a reasonable guesstimate for the average PDR is the average of the relative LTT-slope
	guess_PDR[is.na(guess_PDR)] = default_guess_PDR;
	if(is.na(guess_rholambda0)) guess_rholambda0 = tail(lineage_counter$relative_slopes,1);
	
	# make sure initial guess is in the interior of the imposed bounds
	# being at the edges is not good enough, because the constrOptim algorithm below requires starting points to be in the *interior* of the feasible region
	guess_age_grid		= enforce_min_grid_step(min_values=min_age_grid, max_values=max_age_grid, min_step=min_age_step, guess_values=sort(place_guessed_params_into_interior(min_values=min_age_grid, max_values=max_age_grid, guess_values=guess_age_grid)))
	if(is.null(guess_age_grid)) return(list(success = FALSE, error=sprintf("Could not figure out a reasonable guessed age grid; perhaps the min/max ages are too restrictive, and/or the min_age_step (%g) is too large?",min_age_step)));
	guess_PDR			= place_guessed_params_into_interior(min_values=min_PDR, max_values=max_PDR, guess_values=guess_PDR)
	guess_rholambda0	= place_guessed_params_into_interior(min_values=min_rholambda0, max_values=max_rholambda0, guess_values=guess_rholambda0)
	
	# determine which parameters are to be fitted
	# convention: parameters are indexed as follows: [PDR[], rholambda0]
	fixed_param_values 	= c(fixed_age_grid, fixed_PDR, fixed_rholambda0); # may contain NAs, corresponding to non-fixed parameters
	fitted_params		= which(is.na(fixed_param_values))
	fixed_params		= which(!is.na(fixed_param_values))
	guess_param_values 	= c(guess_age_grid, guess_PDR, guess_rholambda0); # should contain a valid numeric for each parameter, even if the parameter is fixed
	guess_param_values[fixed_params] = fixed_param_values[fixed_params] # make sure guessed param values are consistent with fixed param values
	min_param_values	= c(min_age_grid, min_PDR, min_rholambda0);
	max_param_values	= c(max_age_grid, max_PDR, max_rholambda0);
	NP					= length(fixed_param_values)
	NFP					= length(fitted_params);
	NFG 				= sum(is.na(fixed_age_grid)) # number of fitted grid ages
	
	# determine typical parameter scales
	scale_age_grid 		= seq(from=0, to=oldest_age, length.out=Ngrid); scale_age_grid[1] = oldest_age/Ngrid;
	scale_PDR 			= abs(guess_PDR); scale_PDR[scale_PDR==0] = mean(scale_PDR);
	scale_rholambda0 	= abs(guess_rholambda0);
	if(scale_rholambda0==0) scale_rholambda0 = log2(LTT0)/root_age;
	param_scales = c(scale_age_grid,scale_PDR,scale_rholambda0);
	
	# prepare inequality matrix for all parameters, ensuring that fitted grid ages are always in ascending order and all parameters are within their box constraints
	# inequalities are encoded as: A * x >= B, where x are the fitted parameters
	# every inequality will be a separate row in the matrix A
	ineq_B 	= c()
	ineq_A 	= matrix(0,nrow=0,ncol=NFP)
	zeros 	= rep(0,times=NFP)
	for(fp in 1:NFP){
		p = fitted_params[fp]
		if(p<=Ngrid){
			g = p
			# encode inequalities for g-th grid age
			if((g>1) && (is.na(fixed_age_grid[g-1]))){
				# inequality: age_grid[g] - age_grid[g-1] >= min_age_step
				# age_grid[g-1] is fitted as well, thus corresponding to fp-1
				ineq_A = rbind(ineq_A,zeros); ineq_A[nrow(ineq_A),fp] = 1
				ineq_A[nrow(ineq_A),fp-1] = -1
				ineq_B = c(ineq_B,min_age_step)
			}
			#   inequality: age_grid[g]>=min_age_grid[g]
			if(is.finite(min_age_grid[g])){
				ineq_A = rbind(ineq_A,zeros); ineq_A[nrow(ineq_A),fp] = 1
				ineq_B = c(ineq_B,min_age_grid[g])
			}	
			#   inequality: age_grid[g]<=max_age_grid[g]
			if(is.finite(max_age_grid[g])){
				ineq_A = rbind(ineq_A,zeros); ineq_A[nrow(ineq_A),fp] = -1
				ineq_B = c(ineq_B,-max_age_grid[g])
			}
		}else if(p<=2*Ngrid){
			# encode box-constraint inequalities for g-th PDR
			g = p-Ngrid
			if(is.finite(min_PDR[g])){
				ineq_A = rbind(ineq_A,zeros); ineq_A[nrow(ineq_A),fp] = 1
				ineq_B = c(ineq_B,min_PDR[g])
			}	
			if(is.finite(max_PDR[g])){
				ineq_A = rbind(ineq_A,zeros); ineq_A[nrow(ineq_A),fp] = -1
				ineq_B = c(ineq_B,-max_PDR[g])
			}
		}else{
			# encode inequalities for rholambda0
			# min_rholambda0 <= rholambda0 & rholambda0 <= max_rholambda0
			if(is.finite(min_rholambda0)){
				ineq_A = rbind(ineq_A,zeros); ineq_A[nrow(ineq_A),fp] = 1
				ineq_B = c(ineq_B,min_rholambda0)
			}	
			if(is.finite(max_rholambda0)){
				ineq_A = rbind(ineq_A,zeros); ineq_A[nrow(ineq_A),fp] = -1
				ineq_B = c(ineq_B,-max_rholambda0)
			}
		}
	}


	################################
	# FITTING

	# ensure some fit_control options are named correctly for the constrOptim algorithm used here
	if(is.null(fit_control)) fit_control = list()
	if((!is.null(fit_control$iter.max)) && (is.null(fit_control$maxit))) fit_control$maxit = fit_control$iter.max
	if((!is.null(fit_control$rel.tol)) && (is.null(fit_control$reltol))) fit_control$reltol = fit_control$rel.tol
	if((!is.null(fit_control$abs.tol)) && (is.null(fit_control$abstol))) fit_control$abstol = fit_control$abs.tol
	outer_maxit 	= (if(is.null(fit_control$outer_maxit)) 100 else fit_control$outer_maxit) # max number of optim() calls by the barrier algorithm
	outer_reltol 	= fit_control$outer_reltol
	if(is.null(fit_control$maxit)) fit_control$maxit = outer_maxit*100*NFP # max cumulative number of iterations across all optim calls
	fit_control$outer_maxit = NULL
	fit_control$outer_reltol = NULL
	

	fitted_grid_params = fitted_params[fitted_params!=(2*Ngrid+1)]
	evaluate_model = function(fparam_values, include_gradient){
		param_values = fixed_param_values; param_values[fitted_params] = fparam_values * param_scales[fitted_params];
		if(any(is.nan(param_values)) || any(is.infinite(param_values))) return(list(success=FALSE)); # catch weird cases where params become NaN
		age_grid	= param_values[1:Ngrid];
		PDRs 		= param_values[(Ngrid+1):(2*Ngrid)];
		rholambda0 	= param_values[2*Ngrid+1];
		if(length(age_grid)==1){
			# while age-grid has only one point (i.e., PDRs are constant over time), we need to provide a least 2 grid points to the loglikelihood calculator, spanning the interval [0,oldest_age]
			input_age_grid 	= c(0,oldest_age);
			input_PDRs 		= c(PDRs, PDRs);
		}else{
			input_age_grid 	= age_grid;
			input_PDRs 		= PDRs
		}
		if(include_gradient){
			# calculate differentials of PDR in the directions of the various fitted parameters (fitted grid ages and PDRs)
			diff_PDR_degree = 1
			if(Ngrid==1){
				diff_PDR = (if(is.na(fixed_PDR)) c(1,0,1,0) else numeric())
			}else{
				diff_PDR_all = derivatives_of_grid_curve_CPP(Xgrid=age_grid, Ygrid=PDRs); # yields a 3D array of size (2*Ngrid)*Ngrid*2, flattened in layer-row-major format
				diff_PDR = unlist(lapply(fitted_grid_params, FUN=function(p) diff_PDR_all[((p-1)*Ngrid*(diff_PDR_degree+1) + 1):((p-1)*Ngrid*(diff_PDR_degree+1) + Ngrid*(diff_PDR_degree+1))])) # extract only differentials along fitted parameters
			}
		}else{
			diff_PDR = numeric()
			diff_PDR_degree = 0
		}
		results = get_HBD_PDR_loglikelihood_CPP(branching_ages		= sorted_node_ages,
												oldest_age			= oldest_age,
												rholambda0 			= rholambda0,
												age_grid 			= input_age_grid,
												PDRs 				= input_PDRs,
												splines_degree		= splines_degree,
												condition			= condition,
												relative_dt			= relative_dt,
												runtime_out_seconds	= max_model_runtime,
												diff_PDR			= diff_PDR,
												diff_PDR_degree		= diff_PDR_degree);
		if(results$success && (is.na(results$loglikelihood) || is.nan(results$loglikelihood) || is.infinite(results$loglikelihood))) results$success=FALSE
		if(!results$success) return(results)
		if(include_gradient){
			gradient_full = rep(NA, times=NP)
			gradient_full[2*Ngrid+1] = results$dLL_drholambda0
			gradient_full[fitted_grid_params] = results$dLL_dPDR
			results$gradient = gradient_full[fitted_params] * param_scales[fitted_params]
		}
		return(results)
	}


	# objective function: negated log-likelihood
	# input argument is the subset of fitted parameters, rescaled according to param_scales
	objective_function = function(fparam_values){
		results = evaluate_model(fparam_values, include_gradient=FALSE)
		if(!results$success) return(Inf);
		return(-results$loglikelihood);
	}
	
	
	gradient_function = function(fparam_values){
		if(splines_degree!=1) return(NaN); # only implemented for splines_degree=1
		results = evaluate_model(fparam_values, include_gradient=TRUE)
		if(!results$success) return(rep(Inf,times=NFP));
		return(-results$gradient);
	}
			

	# fit with various starting points
	fit_single_trial = function(trial){
		scales		 = param_scales[fitted_params]
		lower_bounds = min_param_values[fitted_params]
		upper_bounds = max_param_values[fitted_params]
		# randomly choose start values for fitted params
		start_values = guess_param_values[fitted_params]
		if(trial>1){
			boxed_left	= which((!is.infinite(lower_bounds)) & is.infinite(upper_bounds))
			boxed_right	= which((!is.infinite(upper_bounds)) & is.infinite(lower_bounds))
			boxed_dual  = which(!(is.infinite(lower_bounds) | is.infinite(upper_bounds))); # determine fitted params that are boxed, i.e. constrained to within finite lower & upper bounds
			unboxed 	= which(is.infinite(lower_bounds) & is.infinite(upper_bounds))
			if(length(boxed_dual)>0) 	start_values[boxed_dual] = lower_bounds[boxed_dual] + (upper_bounds[boxed_dual]-lower_bounds[boxed_dual]) * runif(n=length(boxed_dual),min=0,max=1)
			if(length(unboxed)>0) 	 	start_values[unboxed]	 = 10**runif(n=length(unboxed), min=-2, max=2) * start_values[unboxed]
			if(length(boxed_left)>0) 	start_values[boxed_left] = sapply(boxed_left, FUN=function(fp) random_semiboxed_left(lower_bound=lower_bounds[fp], default=start_values[fp], typical_scale=scales[fp], orders_of_magnitude=4))
			if(length(boxed_right)>0) 	start_values[boxed_right]= sapply(boxed_right, FUN=function(fp) -random_semiboxed_left(lower_bound=-upper_bounds[fp], default=-start_values[fp], typical_scale=scales[fp], orders_of_magnitude=4))
			
			# make sure randomly chosen start grid ages are in ascending order and within the feasible region
			if(NFG>0){
				start_age_grid = fixed_param_values[1:Ngrid]; start_age_grid[fitted_params[1:NFG]] = sort(start_values[1:NFG]);
				start_age_grid = enforce_min_grid_step(min_values=min_age_grid,max_values=max_age_grid,min_step=min_age_step,guess_values=start_age_grid)
				if(is.null(start_age_grid)) return(list(objective_value=NA, fparam_values = rep(NA,times=NFP), converged=FALSE, loglikelihood=NA, runtime=NA, Niterations=0, Nevaluations=0, start_loglikelihood=NA, start_runtime=NA)); # something went wrong while adjusting the start_age_grid
				start_values[1:NFG] = start_age_grid[fitted_params[1:NFG]]
			}
		}
		# check if start values yield NaN
		start_results = evaluate_model(start_values/scales, include_gradient=(splines_degree==1))
		if(!start_results$success) return(list(	objective_value=NA, fparam_values = rep(NA,times=NFP), converged=FALSE, 
												loglikelihood=NA, runtime=NA,
												Niterations=0, Nevaluations=1, 
												start_loglikelihood=NA, start_runtime=start_results$runtime));
												
		# fine-tune some fitting controls based on the initial model evaluation
		fit_control$fnscale = abs(start_results$loglikelihood)
		if(is.null(fit_control$reltol)) fit_control$reltol = max(1e-30,min(1e-5,0.001/abs(start_results$loglikelihood)))
		if(is.null(outer_reltol)) outer_reltol = max(1e-30,min(1e-5,0.001/abs(start_results$loglikelihood))) # relative tolerance for the outer algorithm
		
		# debug
		start_ages = fixed_age_grid
		start_ages[fitted_params[1:NFG]] = start_values[1:NFG]
		if(any(diff(start_ages)<min_age_step)){
			cat(sprintf("OOPS: age grid steps are at or below minimum:\n"))
			print(start_ages)
		}
		if(any(start_values<=lower_bounds)){
			cat(sprintf("OOPS: some start params are below minimum:\n"))
			print(lower_bounds)
			print(start_values)
		}
		if(any(start_values>=upper_bounds)){
			cat(sprintf("OOPS: some start params are above maximum:\n"))
			print(upper_bounds)
			print(start_values)
		}
			

		# old code. Using gradient to improve starting point. Does not seem to improve the fitting.
		#if(splines_degree==1){
		#	# improve the starting point by first running a preliminary BFGS search
		#	fit = stats::constrOptim(	theta	= start_values/scales, 
		#								f		= objective_function, 
		#								grad	= gradient_function,
		#								ui		= (if(NFP==1) ineq_A*scales else ineq_A%*%diag(scales)), 
		#								ci		= ineq_B, 
		#								mu 				 = 1e-4, 
		#								control 		 = fit_control,
		#								method 			 = "BFGS",
		#								outer.iterations = outer_maxit, 
		#								outer.eps 		 = outer_reltol,
		#								hessian 		 = FALSE)
		#	if(is.finite(fit$value) && all(is.finite(fit$par)) && (fit$value<-start_results$loglikelihood)) start_values = fit$par*scales
		#	start_results = evaluate_model(fit$par, include_gradient=(splines_degree==1))
		#}
		
		fit = stats::constrOptim(	theta	= start_values/scales, 
									f		= objective_function, 
									ui		= (if(NFP==1) ineq_A*scales else ineq_A%*%diag(scales)), 
									ci		= ineq_B, 
									mu 				 = 1e-4, 
									control 		 = fit_control,
									method 			 = "Nelder-Mead",
									outer.iterations = outer_maxit, 
									outer.eps 		 = outer_reltol,
									hessian 		 = FALSE)
		if((!is.finite(fit$value)) || (!all(is.finite(fit$par)))){
			final_results = list()
		}else{
			final_results = evaluate_model(fit$par, include_gradient=FALSE)
		}
		

		# debug
		if(is.finite(fit$value)){
			final_values = fit$par * scales
			final_ages = fixed_age_grid
			final_ages[fitted_params[1:NFG]] = final_values[1:NFG]
			D = ineq_A %*% final_values - ineq_B
			if(any(diff(final_ages)<min_age_step)){
				cat(sprintf("OOPS %d: fitted age grid steps are at or below minimum:\n",trial))
				print(final_ages)
				cat(sprintf("     %d: min(A*X-B)=%g\n",trial,min(D)))
			}
			if(any(final_values<=lower_bounds)){
				cat(sprintf("OOPS %d: some fitted params are below minimum:\n",trial))
				print(lower_bounds)
				print(final_values)
			}
			if(any(final_values>=upper_bounds)){
				cat(sprintf("OOPS %d: some fitted params are above maximum:\n",trial))
				print(upper_bounds)
				print(final_values)
			}
			if(any(D<0)){
				cat(sprintf("OOPS %d: D=A*X-B<0, min(D)=%g\n",trial,min(D)))
				print(D)
			}
		}
		
		return(list(objective_value=fit$value, fparam_values = fit$par*scales, converged=(fit$convergence==0), 
				loglikelihood=final_results$loglikelihood, runtime=final_results$runtime,
				Niterations=fit$outer.iterations, Nevaluations=fit$counts[1], 
				start_loglikelihood=start_results$loglikelihood, start_runtime=start_results$runtime));
		
		# alternative optimization algorithm
		#fit = tryCatch({
		#		NlcOptim::solnl(X 		= start_values/scales, 
		#						objfun 	= objective_function,
		#						A 		= -ineq_A%*%diag(scales),
		#						B 		= -ineq_B,
		#						tolX 	= fit_control$reltol,
		#						tolFun 	= 1e-06,
		#						tolCon 	= 1e-06,
		#						maxnFun = fit_control$maxit*NFP,
		#						maxIter = fit_control$maxit)
		#		}, error = function(e){ list(fn=NaN, par=rep(NaN,times=NFP), counts=c(NaN,NaN)) })
		#if((!is.finite(fit$fn)) || (!all(is.finite(fit$par)))){
		#	final_results = list()
		#}else{
		#	final_results = evaluate_model(fit$par, include_gradient=(splines_degree==1))
		#}
		return(list(objective_value=fit$fn, fparam_values = fit$par*scales, converged=TRUE, Niterations=NULL, Nevaluations=fit$counts[1],
					loglikelihood=final_results$loglikelihood, runtime=final_results$runtime,
					start_loglikelihood=start_results$loglikelihood, start_runtime=start_results$runtime));
	}
	
	################################

	# run one or more independent fitting trials
    if((Ntrials>1) && (Nthreads>1) && (.Platform$OS.type!="windows")){
		# run trials in parallel using multiple forks
		# Note: Forks (and hence shared memory) are not available on Windows
		if(verbose) cat(sprintf("%sFitting %d model parameters (%d trials, parallelized)..\n",verbose_prefix,NFP,Ntrials))
		fits = parallel::mclapply(	1:Ntrials, 
									FUN = function(trial) fit_single_trial(trial), 
									mc.cores = min(Nthreads, Ntrials), 
									mc.preschedule = FALSE, 
									mc.cleanup = TRUE);
	}else{
		# run in serial mode
		if(verbose) cat(sprintf("%sFitting %d model parameters (%s)..\n",verbose_prefix,NFP,(if(Ntrials==1) "1 trial" else sprintf("%d trials",Ntrials))))
		fits = sapply(1:Ntrials,function(x) NULL)
		for(trial in 1:Ntrials){
			fits[[trial]] = fit_single_trial(trial)
		}
	}

	# extract information from best fit (note that some fits may have LL=NaN or NA)
	objective_values	= unlist_with_nulls(sapply(1:Ntrials, function(trial) fits[[trial]]$objective_value))
	loglikelihoods		= unlist_with_nulls(sapply(1:Ntrials, function(trial) fits[[trial]]$loglikelihood))
	runtimes			= unlist_with_nulls(sapply(1:Ntrials, function(trial) fits[[trial]]$runtime))
	start_loglikelihoods= unlist_with_nulls(sapply(1:Ntrials, function(trial) fits[[trial]]$start_loglikelihood))
	start_runtimes		= unlist_with_nulls(sapply(1:Ntrials, function(trial) fits[[trial]]$start_runtime))
	valids				= which((!is.na(objective_values)) & (!is.nan(objective_values)) & (!is.null(objective_values)) & (!is.infinite(objective_values)));
	if(length(valids)==0) return(list(success=FALSE, error=sprintf("Fitting failed for all trials")));
	best 				= valids[which.min(sapply(valids, function(i) objective_values[i]))]
	objective_value		= -fits[[best]]$objective_value;
	loglikelihood		= fits[[best]]$loglikelihood;
	fitted_param_values = fixed_param_values; fitted_param_values[fitted_params] = fits[[best]]$fparam_values;
	fitted_age_grid		= fitted_param_values[1:Ngrid]
	fitted_PDR			= fitted_param_values[(Ngrid+1):(2*Ngrid)]
	fitted_rholambda0	= fitted_param_values[2*Ngrid+1] 
	if(is.null(objective_value) || any(is.na(fitted_param_values)) || any(is.nan(fitted_param_values))) return(list(success=FALSE, error=sprintf("Some fitted parameters are NaN")));
	
	# reverse any time shift due to earlier tree trimming
	oldest_age 		= oldest_age + age0
	root_age 		= root_age + age0
	fitted_age_grid = fitted_age_grid + age0
	
	
	#######################################################################
	# estimate confidence intervals if needed, via parametric bootstrapping
	
	if(Nbootstraps>0){
		if(verbose) cat(sprintf("%sEstimating confidence intervals using %d parametric bootstraps..\n",verbose_prefix,Nbootstraps))
		if(verbose) cat(sprintf("%s  Calculating pulled speciation rate from PDR, for simulating trees..\n",verbose_prefix))
		# first calculate the PSR from the PDR & rholambda0
		#   include a dummy age grid point at the end of age_grid if needed (extrapolating the fitted PDR as a constant), so as to cover the root age
		#   also include a dummy age grid point at the beginning if necessary (extrapolating PDR as a constant), to cover present-day (age 0)
		sim_age_grid = fitted_age_grid
		sim_PDR		 = fitted_PDR
		if(tail(sim_age_grid,1)<root_age){
			# extrapolate sim_PDR (as a constant) all the way to the root (a bit more to avoid rounding errors)
			sim_age_grid = c(sim_age_grid, root_age*1.01)
			sim_PDR 	 = c(sim_PDR, tail(sim_PDR,1));
		}
		if(sim_age_grid[1]>0){
			# extrapolate sim_PDR (as a constant) all the way to the present-day (age 0)
			sim_age_grid = c(0,sim_age_grid)
			sim_PDR		 = c(sim_PDR[1],sim_PDR)
		}
		sim = get_PSR_from_PDR_HBD(	age0			= age0,
									oldest_age 		= tail(sim_age_grid,1),
									age_grid		= sim_age_grid,
									PDR				= sim_PDR,
									rholambda0		= fitted_rholambda0,
									splines_degree	= splines_degree,
									relative_dt		= relative_dt,
									include_nLTT0	= TRUE);
		if(!sim$success) return(list(success=FALSE, error=sprintf("Bootstrapping failed: Could not calculate PSR corresponding to fitted PDR: %s",sim$error), fitted_age_grid=fitted_age_grid, fitted_PDR=fitted_PDR, fitted_rholambda0=fitted_rholambda0, loglikelihood=loglikelihood));
		if(is.null(Ntrials_per_bootstrap)) Ntrials_per_bootstrap = max(1,Ntrials)
		NBA 			 = max(Ngrid*10,1000) # number of ages on whcih to evaluate bootstrap-fitted PDR
		bootstrap_ages	 = seq(from=age0, to=oldest_age, length.out=NBA) # age grid on which bootstrap-fitted PDR will be evaluated
		bootstrap_params = matrix(NA,nrow=Nbootstraps,ncol=length(bootstrap_ages)+1)
		NBsucceeded		 = 0
		for(b in 1:Nbootstraps){
			# simulate model with fitted parameters
			if(verbose) cat(sprintf("%s  Bootstrap #%d..\n",verbose_prefix,b))			
			bootstrap = castor::generate_tree_hbd_reverse(	Ntips			= LTT0/sim$nLTT0,
															crown_age		= root_age,
															age_grid		= sim$ages, 
															PSR				= sim$PSR,
															splines_degree	= 1,
															relative_dt		= relative_dt)
			if(!bootstrap$success) return(list(success=FALSE, error=sprintf("Bootstrapping failed: Could not generate tree for the fitted PDR: %s",bootstrap$error), fitted_age_grid=fitted_age_grid, fitted_PDR=fitted_PDR, fitted_rholambda0=fitted_rholambda0, loglikelihood=loglikelihood));
			bootstrap_tree = bootstrap$trees[[1]]
			
			# fit PSR to simulated tree
			fit = fit_hbd_pdr_and_grid(	tree				= bootstrap_tree, 
										oldest_age			= oldest_age,
										age0				= age0,
										Ngrid				= Ngrid,
										min_age_grid		= original_min_age_grid,
										max_age_grid		= original_max_age_grid,
										min_PDR				= original_min_PDR,
										max_PDR				= original_max_PDR,
										min_rholambda0		= original_min_rholambda0,
										max_rholambda0		= original_max_rholambda0,
										guess_age_grid		= original_guess_age_grid,
										guess_PDR			= original_guess_PDR,
										guess_rholambda0	= original_guess_rholambda0,
										fixed_age_grid		= fixed_age_grid,
										fixed_PDR			= fixed_PDR,
										fixed_rholambda0	= fixed_rholambda0,
										min_age_step		= min_age_step,
										splines_degree		= splines_degree,
										condition			= condition,
										relative_dt			= relative_dt,
										Ntrials				= Ntrials_per_bootstrap,
										Nbootstraps			= 0,
										Nthreads			= Nthreads,
										max_model_runtime	= max_model_runtime,
										fit_control			= fit_control,
										verbose				= verbose,
										verbose_prefix		= paste0(verbose_prefix,"    "))
			if(!fit$success){
				if(verbose) cat(sprintf("%s  WARNING: Fitting failed for this bootstrap: %s\n",verbose_prefix,fit$error))
			}else{
				bootstrap_params[b,1:NBA] = evaluate_spline(Xgrid=fit$fitted_age_grid,Ygrid=fit$fitted_PDR, splines_degree=splines_degree, Xtarget=bootstrap_ages)
				bootstrap_params[b,NBA+1] = fit$fitted_rholambda0
				NBsucceeded 			  = NBsucceeded + 1
			}			
		}
		# calculate standard errors and confidence intervals from distribution of bootstrapped parameters
		standard_errors_flat = sqrt(pmax(0, colMeans(bootstrap_params^2, na.rm=TRUE) - colMeans(bootstrap_params, na.rm=TRUE)^2))
		standard_errors = list(PDR=standard_errors_flat[1:NBA], rholambda0=standard_errors_flat[NBA+1])
		quantiles = sapply(1:ncol(bootstrap_params), FUN=function(p) quantile(bootstrap_params[,p], probs=c(0.25, 0.75, 0.025, 0.975), na.rm=TRUE, type=8))
		CI50lower = list(PDR=quantiles[1,1:NBA], rholambda0=quantiles[1,NBA+1])
		CI50upper = list(PDR=quantiles[2,1:NBA], rholambda0=quantiles[2,NBA+1])
		CI95lower = list(PDR=quantiles[3,1:NBA], rholambda0=quantiles[3,NBA+1])
		CI95upper = list(PDR=quantiles[4,1:NBA], rholambda0=quantiles[4,NBA+1])
	}

	# return results
	return(list(success					= TRUE,
				objective_value			= objective_value,
				objective_name			= "loglikelihood",
				loglikelihood			= loglikelihood,
				runtime					= fits[[best]]$runtime,
				fitted_age_grid			= fitted_age_grid,
				fitted_PDR				= fitted_PDR,
				fitted_rholambda0		= fitted_rholambda0,
				guess_age_grid			= guess_param_values[1:Ngrid],
				guess_PDR				= guess_param_values[(Ngrid+1):(2*Ngrid)],
				guess_rholambda0		= guess_param_values[2*Ngrid+1],
				NFP						= NFP,
				AIC						= 2*NFP - 2*loglikelihood,
				BIC						= log(sum((sorted_node_ages<=oldest_age) & (sorted_node_ages>=age0)))*NFP - 2*loglikelihood,
				converged				= fits[[best]]$converged,
				Niterations				= fits[[best]]$Niterations,
				Nevaluations			= fits[[best]]$Nevaluations,
				all_start_runtimes		= start_runtimes, # start runtimes for all trials
				all_start_loglikelihoods= start_loglikelihoods, # start loglikelihoods for all trials
				all_final_runtimes		= runtimes,
				all_final_loglikelihoods= loglikelihoods,
				bootstrap_ages			= (if(Nbootstraps>0) bootstrap_ages else NULL),
				standard_errors			= (if(Nbootstraps>0) standard_errors else NULL),
				CI50lower				= (if(Nbootstraps>0) CI50lower else NULL),
				CI50upper				= (if(Nbootstraps>0) CI50upper else NULL),
				CI95lower				= (if(Nbootstraps>0) CI95lower else NULL),
				CI95upper				= (if(Nbootstraps>0) CI95upper else NULL)))
}



completement = function(N, indices){
	pool = rep(TRUE,N);
	pool[indices] = FALSE;
	return(which(pool));
}


# given an age grid (in ascending order), calculate the corresponding parameters zeta\in[0,1]^Ngrid
ages2zetas = function(age0,oldest_age,ages){
	Ngrid 	= length(ages)
	lower	= age0
	upper 	= oldest_age
	zetas	= rep(NA,Ngrid)
	for(g in 1:Ngrid){
		zetas[g] = (ages[g]-lower)/(upper-lower);
		lower = ages[g]
	}
	return(zetas)
}


zetas2ages = function(age0,oldest_age,zetas){
	Ngrid 	= length(zetas)
	ages	= rep(NA,Ngrid)
	lower	= age0
	upper 	= oldest_age
	for(g in 1:Ngrid){
		ages[g] = lower + zetas[g] * (upper-lower);
		lower = ages[g]
	}
	return(ages)
}


# fix any parameters for which min & max are specified and identical
fix_constrained_params_if_needed = function(min_values, max_values, fixed_values){
	NP = length(fixed_values)
	free_params_with_min_max = which((!is.na(min_values)) & (!is.na(max_values)) & (is.na(fixed_values)))
	free_params_with_min_equal_max = free_params_with_min_max[which(min_values[free_params_with_min_max]==max_values[free_params_with_min_max])]
	fixed_values[free_params_with_min_equal_max] = min_values[free_params_with_min_equal_max]
	return(fixed_values)
}

# make sure guessed parameter values are in the interior of their box constraints
place_guessed_params_into_interior = function(min_values,max_values,guess_values){
	left_outliers  = which(guess_values<=min_values)
	right_outliers = which(guess_values>=max_values)
	guess_values[left_outliers]  = min_values[left_outliers] + 0.1 * (max_values[left_outliers] - min_values[left_outliers])
	guess_values[right_outliers] = max_values[right_outliers] - 0.1 * (max_values[right_outliers] - min_values[right_outliers])
	return(guess_values)
}


# ensure that adjacent grid ages are strictly greater than min_step apart
# returns NULL upon failure
enforce_min_grid_step = function(min_values,max_values,min_step,guess_values){
	NG = length(guess_values)
	if(NG<=1) return(guess_values)
	if(min_step>0) guess_values = sort(guess_values)
	buffer_left  = guess_values[1] - min_values[1]
	for(g in 2:NG){
		difference = guess_values[g] - guess_values[g-1]
		correction = min_step-difference
		if(correction>=0){
			# push [g-1] as far left as possible until min_step is satisfied
			if(buffer_left>correction){
				# there's enough space on the left
				correction = (if(correction<1e-6*buffer_left) correction+0.01*(buffer_left-correction) else min(1.01*correction,correction+0.5*(buffer_left-correction)))
				guess_values[g-1] = guess_values[g-1] - correction;
			}else{
				# not enough space on the left, so push only as far as it goes
				buffer_right = max_values[g]-guess_values[g]
				if(buffer_right+buffer_left<=correction) return(NULL) # unable to make enough space
				correction = min(1.01*correction, correction+0.5*(buffer_right+buffer_left-correction))
				# there's enough space on left+right to accommodate this correction, but we need to figure out how to allocate it; try to use as much of the left buffer as possible
				correction_right = correction - 0.99*buffer_left
				while(correction_right>=buffer_right){
					correction_right = correction_right/2.0
				}
				guess_values[g-1] = guess_values[g-1] - (correction-correction_right);
				guess_values[g]   = guess_values[g] + correction_right;
			}
		}
		buffer_left = min(guess_values[g]-guess_values[g-1]-min_step, guess_values[g]-min_values[g]);
	}
	return(guess_values)
}


