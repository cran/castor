# Predict various deterministic features of a homogenous birth-death cladogenic model, backward in time.
# The speciation rate lambda and extinction rate mu are specified on a discrete age-grid, and assumed to vary linearly (or polynomially, as splines) between grid points (see "degree" argument).
# This function calculates, among others, the following features over time:
#	Deterministic LTT curve
#	Deterministic total diversity
#	Deterministic shadow diversity
#   Pulled diversification rate (PDR)
# Alternatively, the lambda may be omitted and instead the PDR may be provided together with only the present birth_rate (lambda0).
#
# References:
#	Morlon et al. (2011). Reconciling molecular phylogenies with the fossil record. PNAS 108:16327-16332
#
simulate_deterministic_hbd = function(	Ntips, 						# number of extant species represented in the tree, i.e. after rarefaction. This is equal to the value of the LTT at present.
										oldest_age,					# numeric, specifying how far back (time before present) to simulate the model
										rho				= 1,		# numeric within (0,1], specifying the fraction of extant diversity represented in the tree. Can also be NULL, which is equivalent to setting rarefaction=1.
										age_grid		= NULL,		# either NULL, or empty, or a numeric vector of size NG, listing ages in ascending order, on which birth/mu are specified. If NULL or empty, then lambda and mu mut be a single scalar. The returned time series will be defined on an age-grid that may be finer than this grid.
										lambda			= NULL,		# either NULL, or a single numeric (constant speciation rate over time), or a numeric vector of size NG (listing speciation rates at each age in grid_ages[]).
										mu				= NULL,		# either a single numeric (constant extinction rate over time) or a numeric vector of size NG (listing extinction rates at each age in grid_ages[]).
										PDR				= NULL,		# either NULL, or a single numeric (constant PDR over time), or a numeric vector of size NG (listing PDR at each age in grid_ages[]). Only needed if lambda is NULL.
										lambda0			= NULL,		# either NULL, or a single numeric specifying the present-day speciation rate (i.e. at age 0). Only needed if lambda is NULL.
										splines_degree	= 1,		# integer, either 1 or 2 or 3, specifying the degree for the splines defined by lambda, mu and PDR on the age grid.
										relative_dt		= 1e-3){	# maximum relative time step allowed for integration. Smaller values increase integration accuracy but increase computation time. Typical values are 0.0001-0.001. The default is usually sufficient.	
	# check validity of input variables
	if(is.null(rho)) rho = 1;
	if(is.null(mu)) return(list(success = FALSE, error = sprintf("Missing mu")))
	if(is.null(lambda)){
		if(is.null(PDR)) return(list(success = FALSE, error = sprintf("PDR must be provided when lambda are omitted")))
		if(is.null(lambda0)) return(list(success = FALSE, error = sprintf("lambda0 must be provided when lambda are omitted")))
	}else{
		if(!is.null(lambda0)) return(list(success = FALSE, error = sprintf("lambda0 must not be explicitly provided when lambda are provided (due to potential ambiguity)")))
		if(!is.null(PDR)) return(list(success = FALSE, error = sprintf("PDR must not be explicitly provided when lambda are provided (due to potential ambiguity)")))
	}
	if(is.null(age_grid) || (length(age_grid)<=1)){
		if((!is.null(lambda)) && (length(lambda)!=1)) return(list(success = FALSE, error = sprintf("Invalid number of lambda; since no age grid was provided, you must either provide a single (constant) lambda or none")))
		if(length(mu)!=1) return(list(success = FALSE, error = sprintf("Invalid number of mu; since no age grid was provided, you must provide a single (constant) mu")))
		# create dummy age grid
		NG 			= 2;
		age_grid	= seq(from=0,to=oldest_age,length.out=NG)
		if(!is.null(lambda)) lambda = rep(lambda,times=NG);
		if(!is.null(mu)) mu = rep(mu,times=NG);
	}else{
		NG = length(age_grid);
		if((age_grid[1]>0) || (age_grid[NG]<oldest_age)) return(list(success = FALSE, error = sprintf("Age grid must cover the entire requested age interval, including age 0 and oldest_age (%g)",oldest_age)))
		if((!is.null(lambda)) && (length(lambda)!=1) && (length(lambda)!=NG)) return(list(success = FALSE, error = sprintf("Invalid number of lambda; since an age grid of size %d was provided, you must either provide zero, one or %d lambda",NG,NG)))
		if((length(mu)!=1) && (length(mu)!=NG)) return(list(success = FALSE, error = sprintf("Invalid number of mu; since an age grid of size %d was provided, you must either provide one or %d mu",NG,NG)))
		if((!is.null(lambda)) && (length(lambda)==1)) lambda = rep(lambda,times=NG);
		if((!is.null(mu)) && (length(mu)==1)) mu = rep(mu,times=NG);
	}
	if(!(splines_degree %in% c(0,1,2,3))) return(list(success = FALSE, error = sprintf("Invalid splines_degree: Expected one of 0,1,2,3.")))
		
	# simulate model backward in time
	simulation = simulate_deterministic_HBD_model_CPP(	Ntips 			= Ntips,
														oldest_age		= oldest_age,
														rarefaction 	= rho,
														age_grid 		= age_grid,
														lambdas 		= (if(is.null(lambda)) numeric() else lambda),
														mus 			= mu,
														PDRs	 		= (if(is.null(PDR)) numeric() else PDR),
														lambda0			= (if(is.null(lambda0)) NaN else lambda0),
														splines_degree	= splines_degree,
														relative_dt		= relative_dt);
	if(!simulation$success) return(list(success = FALSE, error = sprintf("Could not simulate model: %s",simulation$error)))
	rholambda0 = rho*simulation$lambda0;

	return(list(success							= TRUE,
				ages							= simulation$refined_age_grid, # potentially refined ages grid, on which all returned variables are defined
				total_diversity					= simulation$total_diversity,
				shadow_diversity				= simulation$shadow_diversity,
				Pmissing						= simulation$Pmissing,
				Pextinct						= simulation$Pextinct,
				LTT								= simulation$LTT,
				lambda							= simulation$lambda,
				mu								= simulation$mu,
				diversification_rate			= simulation$diversification_rate,
				PDR								= simulation$PDR,				# pulled diversification rate
				PND								= simulation$PND,				# pulled normalized diversity
				SER								= rholambda0 - simulation$PDR, 	# shadow extinction rate
				PER								= simulation$lambda0 - simulation$PDR, # pulled extinction rate
				rholambda0						= rholambda0));
}

