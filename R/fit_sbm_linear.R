# Fit a Spherical Brownian Motion (SBM) model with a diffusivity that varies linearly over time
#
fit_sbm_linear = function(	tree, 
							tip_latitudes, 						# numeric vector of size Ntips, listing geographical latitudes of the tips (in decimal degrees)
							tip_longitudes, 					# numeric vector of size Ntips, listing geographical longitudes of the tips (in decimal degrees)
							radius,								# numeric, radius to assume for the sphere (e.g. Earth). Use this e.g. if you want to hange the units in which diffusivity is estimated. Earth's mean radius is about 6371e3 m.
							planar_approximation	= FALSE,	# logical, specifying whether the estimation formula should be based on a planar approximation of Earth's surface, i.e. geodesic angles are converted to distances and then those are treated as if they were Euclideanon a 2D plane. This approximation substantially increases the speed of computations.
							only_basal_tip_pairs	= FALSE,	# logical, specifying whether only immediate sister tips should be considered, i.e. tip pairs with at most 2 edges between the two tips
							min_MRCA_time			= 0,		# numeric, specifying the minimum allowed height (distance from root) of the MRCA of sister tips considered in the fitting. In other words, an independent contrast is only considered if the two sister tips' MRCA has at least this distance from the root. Set min_MRCA_time=0 to disable this filter.
							Ntrials					= 1,		# number of fitting trials to perform, each time starting with random parameter values
							Nthreads				= 1,
							Nbootstraps				= 0,		# (integer) optional number of parametric-bootstrap samples for estimating confidence intervals of fitted parameters. If 0, no parametric bootstrapping is performed. Typical values are 10-100.
							Ntrials_per_bootstrap	= NULL,		# (integer) optional number of fitting trials for each bootstrap sampling. If NULL, this is set equal to Ntrials. A smaller Ntrials_per_bootstrap will reduce computation, at the expense of increasing the estimated confidence intervals (i.e. yielding more conservative estimates of confidence).
							Nsignificance			= 0,		# (integer) optional number of simulations to perform (under a const-diffusivity model) for testing the statistical significance of the fitted slope. Set to 0 to not calculate the significance of the slope.
							fit_control				= list(),	# a named list containing options for the nlminb fitting routine (e.g. iter.max and rel.tol)
							verbose					= FALSE,	# boolean, specifying whether to print informative messages
							verbose_prefix			= ""){		# string, specifying the line prefix when printing messages. Only relevant if verbose==TRUE.
	# basic input error checking
	if(verbose) cat(sprintf("%sChecking input variables..\n",verbose_prefix))
	if(tree$Nnode<2) return(list(success = FALSE, error="Input tree is too small"));
	root_age	= get_tree_span(tree)$max_distance
	Ntrials		= max(1,Ntrials)
	Nthreads	= max(1,Nthreads)
	if(min_MRCA_time<0) min_MRCA_time = Inf;
	Ntrials = pmax(1,Ntrials)
	if((class(tip_latitudes)=="list") && (length(tip_latitudes)==length(tree$tip.label))){
		tip_latitudes = unlist(tip_latitudes)
	}
	if((class(tip_latitudes)=="list") && (length(tip_longitudes)==length(tree$tip.label))){
		tip_longitudes = unlist(tip_longitudes)
	}
	if(is.null(Nbootstraps) || is.na(Nbootstraps) || (Nbootstraps<0)) Nbootstraps = 0;
					
					
	####################################
	# fit time-independent SBM to get a rough estimate

	if(verbose) cat(sprintf("%sFitting const-diffusivity model..\n",verbose_prefix))
	fit_const = fit_sbm_const(	tree,
								tip_latitudes		= tip_latitudes,
								tip_longitudes		= tip_longitudes,
								radius				= radius,
								only_basal_tip_pairs= only_basal_tip_pairs,
								min_MRCA_time		= min_MRCA_time,
								Nbootstraps			= 0)
	if(!fit_const$success) return(list(success=FALSE, error=sprintf("Failed to fit constant-diffusivity model: %s",fit_const$error)))
		

	####################################
	# Fit linear diffusivity model
	
	if(verbose) cat(sprintf("%sPre-computing SBM probability density functor..\n",verbose_prefix))
	SBM_PD_functor = SBM_get_SBM_PD_functor_CPP(max_error = 1e-7, max_Legendre_terms = 200)

	
	if(verbose) cat(sprintf("%sFitting linear-diffusivity model..\n",verbose_prefix))
	diffusivity_functor = function(times, params){
		return(params[1] + (times/root_age)*(params[2]-params[1]))
	}
	time_grid = seq(0,root_age,length.out=2)
	fit_linear = fit_sbm_parametric(tree,
									tip_latitudes			= tip_latitudes,
									tip_longitudes			= tip_longitudes,
									radius 					= radius,
									param_values 			= c(NA,NA),
									param_guess 			= c(fit_const$diffusivity,fit_const$diffusivity),
									diffusivity 			= diffusivity_functor,
									time_grid 				= time_grid,
									planar_approximation	= planar_approximation,
									only_basal_tip_pairs	= only_basal_tip_pairs,
									min_MRCA_time			= min_MRCA_time,
									param_min				= 0.00001*c(fit_const$diffusivity,fit_const$diffusivity),
									param_max				= c(Inf,Inf),
									param_scale				= c(fit_const$diffusivity,fit_const$diffusivity),
									Ntrials 				= Ntrials,
									Nthreads				= Nthreads,
									Nbootstraps				= Nbootstraps,
									Ntrials_per_bootstrap	= Ntrials_per_bootstrap,
									fit_control				= fit_control,
									SBM_PD_functor			= SBM_PD_functor,
									verbose					= verbose,
									verbose_prefix			= sprintf("%s  ",verbose_prefix))
	if(!fit_linear$success) return(list(success=FALSE, error=sprintf("Failed to fit linear model: %s",fit_linear$error)))
	
	# Calculate statistical significance of the slope
	if(Nsignificance>0){
		if(verbose) cat(sprintf("%sEstimating statistical significance of slope using %d simulations..\n",verbose_prefix,Nsignificance))
		abs_diff = abs(fit_linear$param_fitted[1] - fit_linear$param_fitted[2])
		random_abs_diffs = rep(NA,times=Nsignificance)
		for(r in 1:Nsignificance){
			# simulate a const-diffusivity model
			if(verbose) cat(sprintf("%s  Simulation #%d..\n",verbose_prefix,r))
			simulation = simulate_sbm(tree, radius = radius, diffusivity = fit_const$diffusivity)
			if(!simulation$success) return(list(success=FALSE, error=sprintf("Simulation #%d failed: Could not simulate SBM with constant diffusivity: %s",r,simulation$error), param_fitted=fit_linear$param_fitted, loglikelihood=fit_linear$loglikelihood, const_diffusivity=fit_const$diffusivity));
			fit = fit_sbm_parametric(tree,
									tip_latitudes			= simulation$tip_latitudes,
									tip_longitudes			= simulation$tip_longitudes,
									radius 					= radius,
									param_values 			= c(NA,NA),
									param_guess 			= c(fit_const$diffusivity,fit_const$diffusivity),
									diffusivity 			= diffusivity_functor,
									time_grid 				= time_grid,
									planar_approximation	= planar_approximation,
									only_basal_tip_pairs	= only_basal_tip_pairs,
									min_MRCA_time			= min_MRCA_time,
									param_min				= 0.00001*c(fit_const$diffusivity,fit_const$diffusivity),
									param_max				= c(Inf,Inf),
									param_scale				= c(fit_const$diffusivity,fit_const$diffusivity),
									Ntrials 				= Ntrials,
									Nthreads				= Nthreads,
									Nbootstraps				= 0,
									fit_control				= fit_control,
									SBM_PD_functor			= SBM_PD_functor,
									verbose					= verbose,
									verbose_prefix			= sprintf("%s    ",verbose_prefix))
			if(!fit$success){
				if(verbose) cat(sprintf("%s  WARNING: Fitting failed for this simulation: %s\n",verbose_prefix,fit$error))
			}else{
				random_abs_diffs[r] = abs(fit$param_fitted[2]-fit$param_fitted[1])
			}
		}
		significance = sum(random_abs_diffs>=abs_diff,na.rm=TRUE)/sum(!is.nan(random_abs_diffs))
	}
	
	####################################
		
	# return results
	return(list(success					= TRUE,
				objective_value			= fit_linear$objective_value,
				objective_name			= fit_linear$objective_name,
				diffusivity				= setNames(fit_linear$param_fitted,c("root","present")),
				loglikelihood			= fit_linear$loglikelihood,
				NFP						= fit_linear$NFP,
				AIC						= fit_linear$AIC,
				BIC						= fit_linear$BIC,
				converged				= fit_linear$converged,
				Niterations				= fit_linear$Niterations,
				Nevaluations			= fit_linear$Nevaluations,
				trial_start_objectives	= fit_linear$trial_start_objectives,
				trial_objective_values	= fit_linear$trial_objective_values,
				trial_Nstart_attempts	= fit_linear$trial_Nstart_attempts,
				trial_Niterations		= fit_linear$trial_Niterations,
				trial_Nevaluations		= fit_linear$trial_Nevaluations,
				standard_errors			= fit_linear$standard_errors,
				CI50lower				= fit_linear$CI50lower,
				CI50upper				= fit_linear$CI50upper,
				CI95lower				= fit_linear$CI95lower,
				CI95upper				= fit_linear$CI95upper,
				consistency				= fit_linear$consistency,
				const_diffusivity		= fit_const$diffusivity,
				significance			= (if(Nsignificance>0) significance else NULL)))

}



