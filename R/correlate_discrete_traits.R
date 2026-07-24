# Examines the association between two discrete traits via correlation coefficients and (optionally) a likelihood-ratio test
# Can account for optional error & detection probabilities supplied by the user for each trait, assuming that these are independent for the two traits.
correlate_discrete_traits = function(trees, 						# either a single tree in phylo format, or a list of trees
									tip_states1,					# either a 1D integer vector of size Ntips (if only a single tree is provided) or a list of 1D integer vectors (if one or more trees are provided), listing the state of each tip in each tree, for trait 1. Values should be integers between 1 and Nstates1. May also contain NA, if a tip state has not been observed.
									tip_states2,					# Similar to tip_states1, but for trait 2.
									value_probs1		= NULL,		# optional, either a 3D numeric array of size Ntips x Nstates1 x Nstates  (if only a single tree is provided) or a list of such 3D numeric arrays (if one or more trees are provided), listing the observation probabilities for all tips (aka. emission matrixes), for trait 1. Hence, value_probs1[tp,s1,s2] is the probability of observing state s2 at tip tp conditioned on the true state being s1, and conditioned on a recording being made in the first place. For any given tp & s1, the entries in the tube value_probs1[tp,s1,] must sum to 1. value_probs1 can also be NULL, in which case observations are assumed to be error-free.
									value_probs2		= NULL,		# Similar to value_probs1, but for trait 2.
									inclusion_probs1	= NULL,		# optional, either a numeric vectors of size Nstates1 (applying equally to all trees) or a list of Ntrees such numeric vectors, listing the inclusion probabilities for tips depending on state, for trait 1. Hence, inclusion_probs1[[tr]][s1] is the probability of recording the trait conditioned on the true state being s1, for any given tip in tree tr. inclusion_probs1 can also be NULL, in which case the probability of inclusion/recording of a tip's state is assumed to be the same for all tips and independent of state. If inclusion probabilities don't depend on the state of trait 1, i.e., there are no inclusion biases, there is no benefit in providing inclusion_probs1. For any given tree, only the relative values in inclusion_probs1 matter, not their absolute values; they will be normalized to match the number of non-NAs in the tree under the null's stationarity.
									inclusion_probs2	= NULL,		# similar to inclusion_probs1, but for trait 2.
									null_rate_model1 	= "ARD",	# character, rate model to assume for trait 1 under the null hypothesis. Either "ER" or "SYM" or "ARD" or "SUEDE". The format and interpretation is the same as for fit_mk.
									null_rate_model2 	= "ARD",	# character, rate model to assume for trait 1 under the null hypothesis. Either "ER" or "SYM" or "ARD" or "SUEDE". The format and interpretation is the same as for fit_mk.
									root_prior			= "stationary", # character, prior probability distribution of the root's states, used to calculate the each model's likelihood during fitting. For details see the documentation of fit_mk()
									dependent_trait		= 0,		# integer, specifying which trait to consider as the "dependent variable" during full model fitting. Options are 0 (both traits influence each others rates), 1 (trait 2 influences the rates of trait 1, but not vise versa) and 2 (trait 1 influences the rates of trait 2, but not vise versa). Only relevant if include_lrt=TRUE. This is analogous to the dev.var argument in phytools::fitPagel.
									Ntrials 			= 10,		# (int) number of trials (starting points) for fitting each transition matrix to the real data. A greater Ntrials tends to improve fitting accuracy at the expense of more computing time.
									Ntrials_sim 		= 1,		# (int) number of trials (starting points) for fitting each transition matrix to the simulated data. Only relevant if include_lrt==TRUE.
									Nthreads 			= -1,		# integer, number of parallel threads to use whenever applicable. If negative, this defines the number of desired threads minus the number of available cores. Only works on Linux/UNIX/MacOS; for Windows, all computation is done in a single thread.
									Nsimulations		= 1000,		# integer, number of simulations to perform of the null, for estimating statistical significances. A larger Nsims improves the accuracy of the estimated P-values (especially if P values are small), at the expense of longer computation.
									max_model_runtime	= NULL,		# optional maximum time (in seconds) to allocate for each likelihood evaluation per tree, during fitting. Use this to escape from badly parameterized models during fitting (this will likely cause the affected fitting trial to fail). If NULL or <=0, this option is ignored.
									include_lrt			= FALSE,	# whether to include a likelihood-ratio test. This requires additional computation due to many additional model fits. If only the correlation analysis is required, set this to FALSE to save computing resources.
									verbose				= FALSE,	# boolean, specifying whether to print informative messages
									verbose_prefix		= ""){		# string, specifying the line prefix when printing messages. Only relevant if verbose==TRUE.
					
	# basic input checking
	if(verbose) cat(sprintf("%sChecking input variables..\n",verbose_prefix))
	if(Nthreads<0) Nthreads = max(1,parallel::detectCores() + Nthreads)
	if("phylo" %in% class(trees)){
		# trees[] is actually a single tree
		trees = list(trees)
		Ntrees = 1
	}else if("list" %in% class(trees)){
		# trees[] is a list of trees
		Ntrees = length(trees)
	}else{
		stop(sprintf("Unknown data format '%s' for input trees[]: Expected a list of phylo trees or a single phylo tree",class(trees)[1]))
	}
	Ntips_total = sum(sapply(trees, FUN=function(tree) length(tree$tip.label)))
	
	# sanitize input tip states & probabilities
	aux_sanitize_states_and_probs = function(tip_states, value_probs, inclusion_probs, trait){	
		if(!(("list" %in% class(tip_states)) && (length(tip_states)==Ntrees))){
			# something is wrong with tip_states, perhaps we can fix it
			if((Ntrees==1) && (length(tip_states)==length(trees[[1]]$tip.label))){
				tip_states = list(unlist(tip_states))
			}else{
				stop(sprintf("Invalid input format for tip_states%d: Expected a list of vectors, each listing the tip states of a specific tree",trait))
			}
		}
		# check that the correct number of tip states was given for each tree
		for(tr in seq_len(Ntrees)){
			if(length(tip_states[[tr]])!=length(trees[[tr]]$tip.label)) stop(sprintf("Wrong number of trait %d tip states (%d) provided for tree %d, which has %d tips.",trait,length(tip_states[[tr]]),tr,length(trees[[tr]]$tip.label)))
			if(any(tip_states[[tr]][!is.na(tip_states[[tr]])]<=0)) stop(sprintf("Trait %d tip states for tree %d seem problematic - found zero or negative values.",trait,tr))
		}
		tip_states_flat = unlist(tip_states)
		if(diff(range(tip_states_flat, na.rm=TRUE)) == 0) return(list(success=FALSE, error=sprintf("All non-NA tip_states for trait %d are equal, so correlations between the two traits are undefined.",trait)))
		Nstates = length(unique(tip_states_flat[!is.na(tip_states_flat)]))
		# compute empirical probabilities of tip states
		empirical_probs = sapply(seq_len(Nstates), FUN=function(st) mean(tip_states_flat==st, na.rm=TRUE))
		# check and sanitize value_probs
		if(!is.null(value_probs)){
			if(!(("list" %in% class(value_probs)) && (length(value_probs)==Ntrees))){
				# something is wrong with value_probs, perhaps we can fix it
				if((Ntrees==1) && (nrow(value_probs)==length(trees[[1]]$tip.label)) && (ncol(value_probs)==Nstates) && (dim(value_probs)[3]==Nstates)){
					value_probs = list(value_probs)
				}else{
					stop(sprintf("Invalid input format for value_probs%d: Expected a list of %d 3D arrays, each listing the trait %d value probabilities for tips on a specific tree",trait,Ntrees,trait))
				}
			}
			for(tr in seq_len(Ntrees)){
				value_probs0 = value_probs[[tr]]
				Ntips0 		 = length(tip_states[[tr]])
				if(nrow(value_probs0)!=Ntips0) stop(sprintf("Wrong number of trait %d value probs (%d) provided for tree %d, which has %d tips.",trait,nrow(value_probs0),tr,Ntips0))
				if((ncol(value_probs0)!=Nstates) || (dim(value_probs0)[3]!=Nstates)) stop(sprintf("value_probs%d[[%d]] has wrong shape %s: Expected shape %d x %d x %d.",trait,tr,sprint_shape(value_probs0),Ntips0,Nstates,Nstates))
				if(any(abs(apply(value_probs0, MARGIN = c(1, 2), FUN = sum)-1)>0.00001)) stop(sprintf("value_probs%d for tree %d are problematic - Entries in some tubes [t,s,] do not sum to 1",trait,tr))
				if(any(is.na(value_probs0) | is.nan(value_probs0))) stop(sprintf("value_probs0%d for tree %d contain NA or NaN, which is not allowed.",trait,tr))
				if(any((value_probs0<0) | (value_probs0>1))) stop(sprintf("value_probs0%d for tree %d are problematic - Entries must be between 0 and 1.",trait,tr))
			}
		}
		# check and sanitize inclusion_probs
		if(!is.null(inclusion_probs)){
			if("list" %in% class(inclusion_probs)){
				if(length(inclusion_probs)!=Ntrees) stop(sprintf("Invalid number of inclusion_probs given for trait %d: Expected as many vectors as there are trees (%d), but instead got %d",trait,Ntrees,length(inclusion_probs)))
			}else{
				# use the same inclusion probs for all trees
				inclusion_probs = lapply(seq_len(Ntrees), FUN=function(tr) inclusion_probs)
			}
				
			for(tr in seq_len(Ntrees)){
				inclusion_probs0 = inclusion_probs[[tr]]
				if(length(inclusion_probs0)!=Nstates) stop(sprintf("Wrong number of trait %d inclusion probs provided for tree %d: Expected %d probabilities, but got %d",trait,tr,Nstates,length(inclusion_probs0)))
				if(any((inclusion_probs0<0) | (inclusion_probs0>1))) stop(sprintf("inclusion_probs%d for tree %d are problematic - Entries must be between 0 and 1.",trait,tr))
			}
# 		}else{
# 			# estimate inclusion probs based on the fraction of non-NA tip states, separately for each tree
# 			if(verbose) cat(sprintf("%sNote: Inclusion probs were not provided for trait %d, so estimating them based on the fraction of non-NA tip states\n",verbose_prefix,trait))
# 			inclusion_probs = vector(mode="list", length=Ntrees)
# 			for(tr in seq_len(Ntrees)){
# 				tip_states0 = tip_states[[tr]]
# 				Ntips0		= length(tip_states0)
# 				inclusion_probs[[tr]] = matrix((1-mean(is.na(tip_states0))), nrow=Ntips0, nrow=2)
# 			}
		}
		return(list(success=TRUE, tip_states=tip_states, value_probs=value_probs, inclusion_probs=inclusion_probs, Nstates=Nstates, empirical_probs=empirical_probs))
	}
	
	# auxiliary function for constructing the tip_priors matrix for a specific trait, given observed tip states and optionally given value_probs & inclusion_probs
	# tip_priors[[tr]][i,st] will be the probability of the observed state at the i'th tip of tree tr, conditioned on the true state being st
	aux_get_tip_priors = function(tip_states, value_probs, inclusion_probs, Nstates){
		tip_priors = vector(mode="list", length=Ntrees)
		for(tr in seq_len(Ntrees)){
			tip_states0 = tip_states[[tr]]
			Ntips0 		= length(tip_states0)
			if(is.null(value_probs) || is.null(value_probs[[tr]])){
				# no value_probs given, so construct tip_priors assuming that measurements are error-free
				tip_priors0 = matrix(.Machine$double.eps, nrow=Ntips0, ncol=Nstates)
				tip_priors0[cbind(seq_len(Ntips0), tip_states0)] = 1.0
			}else{
				value_probs0 	 = value_probs[[tr]]
				tip_priors0 	 = value_probs0[cbind(seq_len(Ntips0), rep(seq_len(Nstates), each=Ntips0), rep(tip_states0, times=Nstates))]
				dim(tip_priors0) = c(Ntips0,Nstates)

			}
			# apply inclusion probs
			if(is.null(inclusion_probs) || is.null(inclusion_probs[[tr]])){
				# no inclusion probs given for this tree, so estimate those from the fraction of NA tip states
				inclusion_probs0 = 1-mean(is.na(tip_states0))
			}else{
				inclusion_probs0 = inclusion_probs[[tr]]
			}
			tip_priors0[is.na(tip_states0),]  = 1-inclusion_probs0 # inclusion_probs0 is either a scalar or a vector of length Nstates, in which case it is recycled row-wise
			tip_priors0[!is.na(tip_states0),] = tip_priors0[!is.na(tip_states0),] * inclusion_probs0
			tip_priors[[tr]] = tip_priors0
		}

# 
# 		if(is.null(value_probs)){
# 			# no value_probs given, so construct tip_priors assuming that measurements are error-free
# 			for(tr in seq_len(Ntrees)){
# 				tip_states0 	= tip_states[[tr]]
# 				Ntips0 			= length(tip_states0)
# 				tip_priors0 	= matrix(.Machine$double.eps, nrow=Ntips0, ncol=Nstates)
# 				tip_priors0[cbind(seq_len(Ntips0), tip_states0)] = 1.0
# 				tip_priors0[is.na(tip_states0),]  = 1-inclusion_probs0
# 				tip_priors0[!is.na(tip_states0),] = tip_priors0[!is.na(tip_states0),] * inclusion_probs0
# 				tip_priors[[tr]] = tip_priors0
# 			}
# 		}else{
# 			for(tr in seq_len(Ntrees)){
# 				tip_states0 	= tip_states[[tr]]
# 				value_probs0 	= value_probs[[tr]]
# 				Ntips0 			= length(tip_states0)
# 				tip_priors0 	= value_probs0[cbind(seq_len(Ntips0), rep(seq_len(Nstates), each=Ntips0), rep(tip_states0, times=Nstates))]
# 				dim(tip_priors0) = c(Ntips0,Nstates)
# 				if(is.null(inclusion_probs) || is.null(inclusion_probs[[tr]])){
# 					# no inclusion probs given, so estimate those from the fraction of NA tip states
# 					inclusion_probs0 = 1 - mean(is.na(tip_states0))
# 				}else{
# 					inclusion_probs0 = inclusion_probs[[tr]]
# 				}
# 				tip_priors0 = tip_priors0 * inclusion_probs0
# 				if(is.null(dim(inclusion_probs0))){
# 					# inclusion_probs0 is a scalar
# 					tip_priors0[is.na(tip_states0),] = 1-inclusion_probs0
# 				}else{
# 					# inclusion_probs0 is a 2D matrix of size Ntips0 x Nstates
# 					tip_priors0[is.na(tip_states0),] = 1-inclusion_probs0[is.na(tip_states0),]
# 				}
# # 				tip_states0_extended = tip_states0
# # 				tip_states0_extended[is.na(tip_states0)] = Nstates+1
# # 				tip_priors0 	 = value_probs0[cbind(seq_len(Ntips0), rep(seq_len(Nstates), each=Ntips0), rep(tip_states0_extended, times=Nstates))]
# # 				dim(tip_priors0) = c(Ntips0,Nstates) ## reshape into a 2D matrix
# 				tip_priors[[tr]] = tip_priors0
# 			}
# 		}
		return(tip_priors)
	}
	
	# sanitize and merge tip states & probs into lists for easier iteration
	sp1 = aux_sanitize_states_and_probs(tip_states1, value_probs1, inclusion_probs1, trait=1)
	sp2 = aux_sanitize_states_and_probs(tip_states2, value_probs2, inclusion_probs2, trait=2)
	if(!sp1$success) return(list(success=FALSE, error=sp1$error))
	if(!sp2$success) return(list(success=FALSE, error=sp2$error))
	tip_states 		= list(sp1$tip_states, sp2$tip_states)
	value_probs 	= list(sp1$value_probs, sp2$value_probs)
	inclusion_probs	= list(sp1$inclusion_probs, sp2$inclusion_probs)
	Nstates	   		= c(sp1$Nstates, sp2$Nstates)
	empirical_probs	= list(sp1$empirical_probs, sp2$empirical_probs)
	rm(tip_states1, tip_states2, value_probs1, value_probs2, inclusion_probs1, inclusion_probs2, sp1, sp2) # remove outdated variables to avoid bugs below
	
	# construct the tip priors for fitting the null models
	# tip_priors[[trait]][[tr]][i,st] will be the probability of the observed state at the i'th tip of tree tr, conditioned on the true state being st, for a given trait
	tip_priors	= list(	aux_get_tip_priors(tip_states=tip_states[[1]], value_probs=value_probs[[1]], inclusion_probs=inclusion_probs[[1]], Nstates=Nstates[[1]]),
						aux_get_tip_priors(tip_states=tip_states[[2]], value_probs=value_probs[[2]], inclusion_probs=inclusion_probs[[2]], Nstates=Nstates[[2]]))
	
	# Function for computing various test statistics for the input data
	aux_summary_stats = function(tip_states1, tip_states2){
		tip_states1 = unlist(tip_states1)
		tip_states2 = unlist(tip_states2)
		if((diff(range(tip_states1, na.rm=TRUE)) == 0) || (diff(range(tip_states2, na.rm=TRUE)) == 0)){
			# All tip_states1 and/or all tip_states2 are equal, so correlations are not defined
			Pearson = NA
			Spearman= NA
		}else{
			Pearson  = stats::cor(x=tip_states1, y=tip_states2, use="na.or.complete", method="pearson")
			Spearman = stats::cor(x=tip_states1, y=tip_states2, use="na.or.complete", method="spearman")
		}
		NMI = get_mutual_information(X=tip_states1, Y=tip_states2, normalize=TRUE, na.rm=TRUE) # normalized mutual information
		return(list(Pearson=Pearson, Spearman=Spearman, NMI=NMI))
	}
	stats = aux_summary_stats(tip_states[[1]], tip_states[[2]])
	

	# Function for constructing the tip priors for the full (i.e., combined) Mk model
	# Each full_tip_priors[[tr]] will be a 2D matrix of size Ntips x NCstates, with full_tip_priors[[tr]][i,s1+Nstates1*(s2-1)] being equal to the product of probabilities tip_priors[[1]][[tr]][i,s1] * tip_priors[[2]][[tr]][i,s2]
	# Only needed for the likelihood ratio test
	aux_get_full_tip_priors = function(tip_priors1, tip_priors2){
		# form the outer product of the last axes of the two input tip_prior matrixes and immediately store the result in a 2D matrix of size Ntips0 x NCstates in column-major format
		lapply(seq_len(Ntrees), FUN=function(tr){ tip_priors1[[tr]][, rep(seq_len(Nstates[[1]]), times=Nstates[[2]])] * tip_priors2[[tr]][, rep(seq_len(Nstates[[2]]), each=Nstates[[1]])]})
	}
	
	# fit and simulate the null model for each trait
	null_fits 	= vector(mode="list", length=2)
	null_sims 	= vector(mode="list", length=2) # null_sims[[trait]][[r]][[tr]] will be an integer vector of size Ntips[[tr]] with values in 1,..,Nstates[trait], specifying the simulated tip states for the tr'th tree during the r'th simulation for the given trait
	null_Nrates = integer(2)
	for(trait in c(1,2)){
		if(Nstates[[trait]]==1) return(list(success=FALSE, error=sprintf("Trait %d needs to comprise at least 2 states for a meaningful correlation analysis.",trait)))
		if(verbose) cat(sprintf("%sFitting Markov chain null model to trait %d..\n",verbose_prefix,trait))
		Nstates0 = Nstates[[trait]]
		fit = fit_mk(trees					= trees,
					Nstates					= Nstates[[trait]],
					tip_priors 				= tip_priors[[trait]],
					rate_model 				= (if(trait==1) null_rate_model1 else null_rate_model2),
					root_prior 				= root_prior,
					Ntrials 				= Ntrials,
					max_model_runtime		= max_model_runtime,
					optim_max_iterations	= 1000,
					check_input 			= TRUE,
					Nthreads 				= Nthreads)
		if(!fit$success) return(list(success=FALSE, error=sprintf("Could not fit Markov chain null model to trait %d: %s",trait,fit$error)))
		null_fits[[trait]] 	 = fit
		null_Nrates[[trait]] = fit$Nrates
		stationary_probs0 	 = get_stationary_distribution(fit$transition_matrix)

		if(Nsimulations>0){
			if(verbose) cat(sprintf("%sSimulating Markov chain null model for trait %d..\n",verbose_prefix,trait))
			sims = vector(mode="list", length=Ntrees)
			for(tr in seq_len(Ntrees)){
				Ntips0 = length(trees[[tr]]$tip.label)
				# simulate the true tip states
				states0 = simulate_mk_model(tree				= trees[[tr]], 
											Q					= fit$transition_matrix, 
											root_probabilities	= "stationary",
											include_tips		= TRUE, 
											include_nodes		= FALSE, 
											Nsimulations		= Nsimulations,
											drop_dims			= FALSE)$tip_states
				# map the true tip states to "measured" tip states or NAs
				if(is.null(value_probs[[trait]])){
					sims0 = states0
				}else{
					# simulate errors in the 'observed' tip states based on the simulated true tip states and the provided value probs
					# Hence, each state states0[sm,tp] is mapped to a randomly chosen observation sims0[sm,tp] according to the probabilities value_probs[[trait]][[tr]][tp,states0[sm,tp],]
					probs0 = matrix(value_probs[[trait]][[tr]][cbind(rep(seq_len(Ntips0), each=Nsimulations, times=Nstates0), rep(as.vector(states0),times=Nstates0), rep(seq_len(Nstates0),each=Nsimulations*Ntips0))], nrow=Nsimulations*Ntips0, ncol=Nstates0) # probs0[i,] are the probabilities of observed states for each unique (r,tp) pair, where (r,tp) enumerates all simulation-tip combos in column-major format				
					sims0  = matrix(apply(probs0, 1, FUN=function(p){ sample(Nstates0, 1, prob=p) }), nrow=Nsimulations, ncol=Ntips0, byrow=FALSE)				
				}
				if(is.null(inclusion_probs[[trait]]) || is.null(inclusion_probs[[trait]][[tr]])){
					# set the same tips to NA as in the input data
					sims0[,is.na(tip_states[[trait]][[tr]])] = NA
				}else{
					# rescale this tree's inclusion_probs for simulation purposes, to be consistent with the fraction of observed non-NAs in the tree under stationarity
					inclusion_probs0 = inclusion_probs[[trait]][[tr]]
					inclusion_probs0 = inclusion_probs0 * (1-mean(is.na(tip_states[[trait]])))/sum(stationary_probs0*inclusion_probs0)
					inclusion_probs0 = inclusion_probs0 * 1/max(inclusion_probs0) # ensure inclusion_probs0 don't exceed 1
					# randomly choose some simulated observations to be NA, according to state-dependent inclusion probabilities
					sims0[runif(length(states0)) < (1-inclusion_probs0[states0])] = NA
# 					na_probs0 = 1-matrix(inclusion_probs[[trait]][[tr]][cbind(as.vector(col(states0)), as.vector(states0))], nrow=nrow(states0), ncol=ncol(states0)) # na_probs0[sm,i] is the probability of a NA observation at tip i during simulation sm
# 					sims0[runif(length(states0)) < na_probs0] = NA
				}
				sims[[tr]] = sims0
			}
			null_sims[[trait]] = lapply(seq_len(Nsimulations), FUN=function(r){ lapply(seq_len(Ntrees), FUN=function(tr){ sims[[tr]][r,] })})
		}
	}
	null_loglikelihood = null_fits[[1]]$loglikelihood + null_fits[[2]]$loglikelihood # combined/joint likelihood of the two null models
	
	# compute correlation coefficients for simulated data
	# Eventually, sim_stats will be a data.frame with size Nsimulations x Ncorrelation_coefficients
	if(verbose && (Nsimulations>0)) cat(sprintf("%sComputing correlations for simulated data..\n",verbose_prefix))
	sim_stats = parallel::mclapply(seq_len(Nsimulations), FUN=function(r){ aux_summary_stats(null_sims[[1]][[r]], null_sims[[2]][[r]]) })
	sim_stats = data.frame(	"Pearson"  = vapply(sim_stats, '[[', 0.0, "Pearson"),
							"Spearman" = vapply(sim_stats, '[[', 0.0, "Spearman"),
							"NMI"      = vapply(sim_stats, '[[', 0.0, "NMI"))
					
	# include naive P-values for the various correlations, i.e., under a simple permutation null model ignoring phylogenetic relationships
	# Eventually, naive_stats will be a data.frame with size Nsimulations x Ncorrelation_coefficients
	if(verbose && (Nsimulations>0)) cat(sprintf("%sComputing correlations for naively permutted data..\n",verbose_prefix))
	flat_tip_states = list(unlist(tip_states[[1]]), unlist(tip_states[[2]])) # flat_tip_states[[k]] will be an integer vector listing the states of all tips across all trees for the k'th trait, i.e., concatenated across all trees
	naive_stats = parallel::mclapply(seq_len(Nsimulations), FUN=function(r){ aux_summary_stats(sample(flat_tip_states[[1]],size=Ntips_total,replace=FALSE), flat_tip_states[[2]]) })
	naive_stats = data.frame("Pearson"  = vapply(naive_stats, '[[', 0.0, "Pearson"),
							"Spearman" = vapply(naive_stats, '[[', 0.0, "Spearman"),
							"NMI"      = vapply(naive_stats, '[[', 0.0, "NMI"))

	# compute the loglikelihood difference for the LRT, if needed
	# This involves fitting a full Mk model whose state space spans all state combos between the two traits
	# By convention, we enumerate state combos (st1,st2) in column-major format, i.e., with trait1 being the faster-changing index.
	# Simultaneous transitions of states in both traits are not allowed, i.e., they are fixed at zero.
	if(include_lrt){
		if(verbose) cat(sprintf("%sFitting combined Markov chain model to both traits to get LRT..\n",verbose_prefix))
		NCstates = Nstates[[1]]*Nstates[[2]] # number of possible state combos (st1,st2)
		# construct tip_priors for the combined traits
		full_tip_priors = aux_get_full_tip_priors(tip_priors[[1]], tip_priors[[2]])
		# construct the combined rate model
		full_rate_model = combine_Mk_rate_models(Nstates1			= Nstates[[1]], 
												 Nstates2			= Nstates[[2]],
												 rate_model1		= null_rate_model1,
												 rate_model2		= null_rate_model2,
												 transition_matrix1 = null_fits[[1]]$transition_matrix,
												 transition_matrix2 = null_fits[[2]]$transition_matrix,
												 dependent_trait	= dependent_trait)
		# fit the full model
		full_fit = fit_mk(trees					= trees,
						Nstates					= full_rate_model$Nstates,
						tip_priors 				= full_tip_priors,
						rate_model 				= full_rate_model$index_matrix,
						root_prior 				= root_prior,
						guess_transition_matrix	= full_rate_model$transition_matrix,
						Ntrials 				= Ntrials,
						max_model_runtime		= max_model_runtime,
						optim_max_iterations	= 1000,
						check_input 			= TRUE,
						Nthreads 				= Nthreads)
		if(!full_fit$success) return(list(success=FALSE, error=sprintf("Could not fit full Markov chain model to the combined traits: %s",full_fit$error)))
		LRT 		= 2*(full_fit$loglikelihood - null_loglikelihood) # the likelihood ratio test statistic
		LRT_df 		= full_rate_model$Nrates-sum(null_Nrates)
		LRT_P1_chi2 = pchisq(LRT, df=LRT_df, lower.tail=FALSE) # chi2-approximation of the significance of the LRT statistic, valid asymptotically when the number of trait transitions in the tree is very large
	}
	
	# Function for computing the LRT test statistic for a single simulated dataset
	# Used for paralellization purposes
	aux_lrt_single_simulation = function(r, # integer in 1,..,Nsimulations, the simulation index to focus on
										 paralellized=TRUE){ # whether this function is being called within a parallelized loop, i.e., running concurrently with other similar function calls
		# fit null models to this simulated dataset
		sim_null_loglikelihood = 0
		sim_tip_priors = vector(mode="list", length=2)
		for(trait in c(1,2)){
			sim_tip_priors[[trait]] = aux_get_tip_priors(tip_states=null_sims[[trait]][[r]], value_probs=value_probs[[trait]], inclusion_probs=inclusion_probs[[trait]], Nstates=Nstates[[trait]])
			fit = fit_mk(trees					= trees,
						Nstates					= Nstates[[trait]],
						tip_priors 				= sim_tip_priors[[trait]],
						rate_model 				= (if(trait==1) null_rate_model1 else null_rate_model2),
						root_prior 				= "max_likelihood",
						guess_transition_matrix	= null_fits[[trait]]$transition_matrix,
						Ntrials 				= Ntrials_sim,
						max_model_runtime		= max_model_runtime,
						optim_max_iterations	= 1000,
						check_input 			= TRUE,
						Nthreads 				= (if(paralellized) 1 else Nthreads))
			if(!fit$success) return(NA)
			sim_null_loglikelihood = sim_null_loglikelihood + fit$loglikelihood
		}
		
		# fit the full model to this simulated dataset
		sim_full_tip_priors = aux_get_full_tip_priors(sim_tip_priors[[1]], sim_tip_priors[[2]])
		sim_full_fit = fit_mk(	trees					= trees,
								Nstates					= full_rate_model$Nstates,
								tip_priors 				= sim_full_tip_priors,
								rate_model 				= full_rate_model$index_matrix,
								root_prior 				= "max_likelihood",
								guess_transition_matrix	= full_fit$transition_matrix,
								Ntrials 				= Ntrials_sim,
								max_model_runtime		= max_model_runtime,
								optim_max_iterations	= 1000,
								check_input 			= TRUE,
								Nthreads 				= (if(paralellized) 1 else Nthreads))
		if(!fit$success) return(NA)
		sim_LRT = 2*(sim_full_fit$loglikelihood - sim_null_loglikelihood) # the likelihood ratio test statistic for this simulated dataset
		return(sim_LRT)
	}	

	# estimate the statistical significance of the LRT statistic using simulations
	if(include_lrt && (Nsimulations>0)){
		if(verbose) cat(sprintf("%sEstimating the significance of the LRT by re-fitting models to simulated data..\n",verbose_prefix))
		if((Nsimulations>1) && (Nthreads>1) && (.Platform$OS.type!="windows")){
			# run trials in parallel using multiple forks
			# Note: Forks (and hence shared memory) are not available on Windows
			LRT_sims = parallel::mclapply(	seq_len(Nsimulations), 
											FUN 			= aux_lrt_single_simulation,
											mc.cores 		= min(Nthreads, Nsimulations), 
											mc.preschedule 	= FALSE, 
											mc.cleanup 		= TRUE)
			LRT_sims = unlist(LRT_sims)
		}else{
			# run in serial mode
			LRT_sims = rep(NA, Nsimulations)
			for(r in seq_len(Nsimulations)){
				LRT_sims[r] = aux_lrt_single_simulation(r, paralellized=FALSE)
			}
		}
	}else{
		LRT_sims = numeric(0)
	}


	results = list(	success				= TRUE, 
					Nstates1			= Nstates[[1]],
					Nstates2			= Nstates[[2]],
					null_fit1			= null_fits[[1]],
					null_fit2			= null_fits[[2]],
					null_loglikelihood	= null_loglikelihood,
					null_AIC			= (null_fits[[1]]$AIC+null_fits[[2]]$AIC))
	for(stat in colnames(sim_stats)){
		results[[stat]] = stats[[stat]]
		results[[sprintf("%s_P1",stat)]] 		= (if(stats[[stat]]>=0) mean(sim_stats[,stat]>=stats[[stat]], na.rm=TRUE) else mean(sim_stats[,stat]<=stats[[stat]], na.rm=TRUE)) # one-sided significance
		results[[sprintf("%s_mean_null",stat)]] = mean(sim_stats[,stat], na.rm=TRUE)
		results[[sprintf("%s_std_null",stat)]]  = sd(sim_stats[,stat], na.rm=TRUE)
		results[[sprintf("%s_sims",stat)]] 		= sim_stats[,stat]
		if(stat!="NMI") results[[sprintf("%s_P2",stat)]] = mean(abs(sim_stats[,stat])>=abs(stats[[stat]]), na.rm=TRUE) # two-sided significance
	}
	for(stat in colnames(naive_stats)){
		results[[sprintf("%s_P1_naive",stat)]] = (if(stats[[stat]]>=0) mean(naive_stats[,stat]>=stats[[stat]], na.rm=TRUE) else mean(naive_stats[,stat]<=stats[[stat]], na.rm=TRUE)) # one-sided significance
		if(stat!="NMI") results[[sprintf("%s_P2_naive",stat)]] = mean(abs(naive_stats[,stat])>=abs(stats[[stat]]), na.rm=TRUE) # two-sided significance
	}
	if(include_lrt){
		results$LRT 				= LRT # likelihood ratio test statistic
		results$LRT_P1_chi2 		= LRT_P1_chi2 # chi2-approximation of the one-sided significance of the likelihood ratio test statistic
		results$LRT_df				= LRT_df # degrees of freedom of asymptotic chi2 approximation of the LRT statistic
		results$LRT_P1				= mean(LRT_sims>=LRT, na.rm=TRUE) # one-sided significance of the likelihood ratio test statistic estimated from the simulations
		results$LRT_mean_null 		= mean(LRT_sims, na.rm=TRUE)
		results$LRT_std_null 		= sd(LRT_sims, na.rm=TRUE)
		results$LRT_sims			= LRT_sims
		results$full_fit 			= full_fit
		results$full_loglikelihood	= full_fit$loglikelihood
		results$full_AIC			= full_fit$AIC
	}
	return(results)
}





