# generate a random phylogenetic tree, by randomly splitting tips at a certain rate and randomly killing tips at a certain rate
# birth & death rates (=speciation & extintion rates) are themselves evolving under one of the following models:
#	'BM': Constrained Brownian Motion model (constrained into an interval via reflection)
#	'Mk': Discrete-state continuous-time Markov model, as defined by the transition rate matrix
# The simulation is halted as soon as Ntips>=max_tips (if max_tips>0) and/or time>=max_time (if max_time>0) and/or time>=max_time_eq+equilibrium_time (if max_time_eq>=0)
generate_tree_with_evolving_rates = function(parameters				= list(), 	# named list of model parameters. For entries and default values (different for each rate_model) see the main function body below
											 rate_model				= 'BM',		# how speciation/extinction rates are evolving alogn lineages over time. Options are 'BM' (brownian motion) or 'Mk'.
											 max_tips				= NULL, 
											 max_time				= NULL,
											 max_time_eq			= NULL,
											 coalescent 			= TRUE,
											 as_generations			= FALSE,	# if FALSE, then edge lengths correspond to time. If TRUE, then edge lengths correspond to generations (hence if coalescent==false, all edges will have unit length).
											 all_Mk_transitions		= TRUE,		# simulate all Mk transition events along edges (not just at branching points). This should typically be TRUE.
											 Nsplits				= 2,	 	# number of children generated at each diversification event. If set to 2, a bifurcating tree is generated. If >2, the tree will be multifurcating.
											 tip_basename			= "",		# basename for tips (e.g. "tip."). 
											 node_basename			= NULL,		# basename for nodes (e.g. "node."). If NULL, then nodes will not have any labels.
											 include_birth_times	= FALSE,
											 include_death_times	= FALSE,
											 include_rates			= FALSE){
	if(is.null(max_tips) && is.null(max_time) && is.null(max_time_eq)) stop("ERROR: At least one of max_tips and/or max_time and/or max_time_eq must be non-NULL")
	
	# set default model parameters
	if(rate_model=="BM"){
		if(is.null(parameters$birth_rate_diffusivity)) 	parameters$birth_rate_diffusivity = 1;
		if(is.null(parameters$min_birth_rate_pc)) 		parameters$min_birth_rate_pc = 0;
		if(is.null(parameters$max_birth_rate_pc)) 		parameters$max_birth_rate_pc = 1;
		if(is.null(parameters$death_rate_diffusivity)) 	parameters$death_rate_diffusivity = 1;
		if(is.null(parameters$min_death_rate_pc))		parameters$min_death_rate_pc = 0;
		if(is.null(parameters$max_death_rate_pc)) 		parameters$max_death_rate_pc = 1;
	}else if(rate_model=="Mk"){
		if(is.null(parameters$Nstates)) 			parameters$Nstates = 1;
		if(is.null(parameters$transition_matrix))	parameters$transition_matrix = matrix(0,nrow=parameters$Nstates, ncol=parameters$Nstates);
		# pc birth rates corresponding to each state
		if(is.null(parameters$state_birth_rates)){
			parameters$state_birth_rates = rep(1, times=parameters$Nstates);
		}else if(length(parameters$state_birth_rates)==1){
			parameters$state_birth_rates = rep(parameters$state_birth_rates, times=parameters$Nstates);
		}
		# pc death rates corresponding to each state
		if(is.null(parameters$state_death_rates)){
			parameters$state_death_rates = rep(1, times=parameters$Nstates);
		}else if(length(parameters$state_death_rates)==1){
			parameters$state_death_rates = rep(parameters$state_death_rates, times=parameters$Nstates);
		}
	}else{
		stop(sprintf("ERROR: Invalid rate_model '%s'",rate_model))
	}
	if(is.null(parameters$rarefaction)) parameters$rarefaction = 1;
	
	if(parameters$rarefaction<=0 || parameters$rarefaction>1) stop("ERROR: rarefaction parameter must be between 0 (non-inclusive) and 1 (inclusive).")
	
	if(rate_model=="BM"){
		results = generate_random_tree_BM_rates_CPP(max_tips					= (if(is.null(max_tips)) -1 else max_tips),
													max_time					= (if(is.null(max_time)) -1 else max_time),
													max_time_since_equilibrium	= (if(is.null(max_time_eq)) -1 else max_time_eq),
													birth_rate_diffusivity 		= parameters$birth_rate_diffusivity, 
													min_birth_rate_pc 			= parameters$min_birth_rate_pc,
													max_birth_rate_pc 			= parameters$max_birth_rate_pc, 
													death_rate_diffusivity 		= parameters$death_rate_diffusivity,
													min_death_rate_pc			= parameters$min_death_rate_pc,
													max_death_rate_pc			= parameters$max_death_rate_pc,
													coalescent					= coalescent,
													Nsplits						= Nsplits,
													as_generations				= as_generations,
													include_birth_times			= include_birth_times,
													include_death_times			= include_death_times,
													include_rates				= include_rates);
	}else{
		results = generate_random_tree_Mk_rates_CPP(max_tips					= (if(is.null(max_tips)) -1 else max_tips),
													max_time					= (if(is.null(max_time)) -1 else max_time),
													max_time_since_equilibrium	= (if(is.null(max_time_eq)) -1 else max_time_eq),
													Nstates						= parameters$Nstates,
													state_birth_rates			= parameters$state_birth_rates, 
													state_death_rates			= parameters$state_death_rates,
													root_state					= sample.int(n=parameters$Nstates,size=1)-1,
													transition_matrix 			= as.vector(t(parameters$transition_matrix)), # flatten in row-major format
													coalescent					= coalescent,
													Nsplits						= Nsplits,
													as_generations				= as_generations,
													all_transitions				= all_Mk_transitions,
													include_birth_times			= include_birth_times,
													include_death_times			= include_death_times,
													include_rates				= include_rates);
	
	}
	if(!results$success) return(list(success=FALSE, error=results$error)); # something went wrong
	Ntips	= results$Ntips
	Nnodes 	= results$Nnodes
	tree = list(Nnode 		= Nnodes,
				tip.label 	= paste(tip_basename, 1:Ntips, sep=""),
				node.label 	= (if(is.null(node_basename)) NULL else paste(node_basename, 1:Nnodes, sep="")),
				edge 		= matrix(results$tree_edge,ncol=2,byrow=TRUE) + 1L,
				edge.length = results$edge_length,
				root 		= results$root+1L)
	class(tree) = "phylo";
	
	
	# rarefy if needed
	Nrarefied = 0;
	if(parameters$rarefaction<1){
		rarefaction_depth = parameters$rarefaction*Ntips;
		if(rarefaction_depth<2){
			return(list(success=FALSE, error=sprintf("Rarefaction (%g) is too low for the generated tree (%d tips)", parameters$rarefaction,Ntips)))
		}
		keep_tips 	= sample.int(n=Ntips, size=rarefaction_depth, replace=FALSE)
		rarefaction = castor::get_subtree_with_tips(tree, only_tips=keep_tips, omit_tips=FALSE, collapse_monofurcations=TRUE)
		tree 		= rarefaction$subtree
		Nrarefied 	= Ntips - length(tree$tip.label)
		results$root_time = results$root_time + rarefaction$root_shift; # update root time, in case root has changed
		if(include_rates){
			results$birth_rates_pc	= c(results$birth_rates_pc[rarefaction$new2old_tip],results$birth_rates_pc[Ntips+rarefaction$new2old_node])
			results$death_rates_pc	= c(results$death_rates_pc[rarefaction$new2old_tip],results$death_rates_pc[Ntips+rarefaction$new2old_node])
		}
		Ntips 	= length(tree$tip.label)
		Nnodes 	= results$Nnodes
		if(!is.null(results$clade_states)) results$clade_states = results$clade_states[rarefaction$new2old_tip]
	}
	
	return(list(success				= TRUE,
				tree				= tree,
				root_time			= results$root_time,
				final_time			= results$final_time,
				equilibrium_time	= results$equilibrium_time,
				Nbirths		 		= results$Nbirths,
				Ndeaths				= results$Ndeaths,
				Nrarefied			= Nrarefied, # number of tips removed via rarefaction at the end
				states				= results$clade_states+1L, # only relevant for discrete-state rate models
				birth_times			= (if(include_birth_times) results$birth_times else NULL),
				death_times			= (if(include_death_times) results$death_times else NULL),
				birth_rates_pc		= (if(include_rates) results$birth_rates_pc else NULL),
				death_rates_pc		= (if(include_rates) results$death_rates_pc else NULL)));
	
}