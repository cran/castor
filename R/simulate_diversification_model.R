# Simulate a deterministic speciation/extinction model
# The model is parameterized similarly to tree models, but the returned value is a time series of diversity (number of clades) over time
simulate_diversification_model = function(	times,									# times at which to calculate diversities, in ascending order
											parameters				= list(), 		# named list of model parameters. For entries and default values see the main function body below
											start_time				= NULL,			# If reverse==TRUE, then this is the end time of the tree (>=times[NT]), otherwise this is the beginning time of the tree (<=times[1]). If NULL, this is set to times[1] or times[NT], depending on reverse.
											start_diversity			= 1,			# diversity of extant clades at start_time. If reverse==TRUE, then this is the visible diversity after rarefaction.
											coalescent				= FALSE,		# (boolean) if true, the diversities at any time are replaced by the diversities that would remain in a coalescent tree (i.e. only including ancestors of extant tips)
											reverse					= FALSE,		# (boolean) if true, then the tree model is integrated in backward time direction. In that case, start_diversity is interpreted as the true diversity at times.back()
											include_event_rates		= FALSE,		# (boolean) include birth & death rates in returned values
											include_Nevents			= FALSE,		# (boolean) include an estimate of the total birth (speciation) and death (extinction) events during the simulation
											max_runtime				= NULL){		# (numeric) max allowed runtime in seconds. If NULL or <=0, this option is ignored
	NT = length(times)
	
	# basic error checking
	if((!is.null(start_time)) && (!reverse) && (start_time>times[1])) stop(sprintf("start_time must be equal to or smaller than the first requested time point (got start_time=%g, first requested time point = %g)",start_time,times[1]))
	if((!is.null(start_time)) && reverse && (start_time<times[NT])) stop(sprintf("start_time must be equal to or larger than the last requested time point (got start_time=%g, last requested time point = %g)",start_time,times[NT]))
	if(is.null(start_time)) start_time = (if(reverse) times[NT] else times[1]);
							
	# set default model parameters
	if(is.null(parameters$birth_rate_intercept)) 	parameters$birth_rate_intercept = 0;
	if(is.null(parameters$birth_rate_factor)) 		parameters$birth_rate_factor = 0;
	if(is.null(parameters$birth_rate_exponent)) 	parameters$birth_rate_exponent = 1;
	if(is.null(parameters$death_rate_intercept)) 	parameters$death_rate_intercept = 0;
	if(is.null(parameters$death_rate_factor))		parameters$death_rate_factor = 0;
	if(is.null(parameters$death_rate_exponent)) 	parameters$death_rate_exponent = 1;
	if(is.null(parameters$rarefaction)) 			parameters$rarefaction = 1;
						
	# run simulation		
	simulation = simulate_deterministic_diversity_growth_CPP(	birth_rate_intercept 		= parameters$birth_rate_intercept,
																birth_rate_factor 			= parameters$birth_rate_factor,
																birth_rate_exponent 		= parameters$birth_rate_exponent,
																death_rate_intercept 		= parameters$death_rate_intercept,
																death_rate_factor 			= parameters$death_rate_factor,
																death_rate_exponent 		= parameters$death_rate_exponent,
																rarefaction					= parameters$rarefaction,
																Nsplits 					= 2,
																times 						= times,
																start_time					= start_time,
																start_diversity				= start_diversity,
																reverse						= reverse,
																coalescent					= coalescent,
																include_probabilities		= TRUE,
																include_birth_rates			= include_event_rates,
																include_death_rates			= include_event_rates,
																include_Nevents				= include_Nevents,
																runtime_out_seconds			= (if(is.null(max_runtime)) 0 else max_runtime));

	# get probabilities of survival/discovery/representation
	# note that the simulation may return empty vectors (e.g. if this calculation was not requested) or NULL (calculation not available)
	Psurvival 		= simulation$Psurvival	 # probability of survival to the present (regardless of discovery)
	Pdiscovery		= simulation$Pdiscovery
	Prepresentation	= simulation$Prepresentation # probability of survival to the present and discovery
	if((!is.null(Psurvival)) && (length(Psurvival)==0)) Psurvival = NULL
	if((!is.null(Pdiscovery)) && (length(Pdiscovery)==0)) Pdiscovery = NULL
	if((!is.null(Prepresentation)) && (length(Prepresentation)==0)) Prepresentation = NULL
	if(is.null(Psurvival) && (!is.null(Prepresentation)) && (!is.null(Pdiscovery))){
		Psurvival = Prepresentation/Pdiscovery
	}else if(is.null(Prepresentation) && (!is.null(Psurvival)) && (!is.null(Pdiscovery))){
		Prepresentation = Psurvival * Pdiscovery
	}else if(is.null(Pdiscovery) && (!is.null(Psurvival)) && (!is.null(Prepresentation))){
		Pdiscovery = Prepresentation/Psurvival
	}

	return(list(success				= simulation$success,
				diversities			= simulation$diversities,
				Psurvival			= Psurvival,
				Pdiscovery			= Pdiscovery,
				Prepresentation		= Prepresentation,
				birth_rates			= (if(include_event_rates) simulation$birth_rates else NULL),
				death_rates			= (if(include_event_rates) simulation$death_rates else NULL),
				Nbirths				= simulation$Nbirths,
				Ndeaths				= simulation$Ndeaths,
				error				= (if(!simulation$success) simulation$error else NULL)));
}
