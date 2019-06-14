# simulate a Discrete-State Speciation and Extinction (dSSE) model, whereby tree generation is coupled with the evolution of a discrete trait
# birth & death rates (=speciation & extintion rates) are evolving according to a discrete-state continuous-time Markov model.
# Transitions between states can occur either along an edge (anagenetically) and/or during speciation events (cladogenetically).
# Anagenetic transition rates are specified through a transition_matrix_A, while cladogenic transition proabilities are specified through a transition_matrix_C
# The simulation is halted as soon as Ntips>=max_tips (if max_tips>0) and/or time>=max_time (if max_time>0) and/or time>=max_time_eq+equilibrium_time (if max_time_eq>=0)
simulate_dsse = function(	Nstates,							# number of discrete possible states for the trait
							NPstates				= NULL,		# optional number of proxy states, for hiding the original states (i.e. according to a Hidden State Speciation Extinction model)
							proxy_map				= NULL,		# optional 1D integer vector of size Nstates, mapping states to proxy-states, in the case of a HiSSE model. Hence, proxy_map[s] is an integer in 1:NPstates, specifying which proxy-state the state s belongs to. Only relevant if NPstates!=NULL and NPstates!=Nstates
							parameters				= list(), 	# named list of dSSE model parameters. For names and default values see the main function body below.
							start_state				= NULL,		# integer between 1 and Nstates, specifying the state of the first lineage. If NULL, the root state is chosen randomly.
							max_tips				= NULL, 	# integer, specifying the max number of tips in the simulated tree (prior to any subsampling)
							max_time				= NULL,
							max_time_eq				= NULL,
							sampling_fractions		= NULL,		# numeric vector of size NPstates, listing sampling fractions depending on state. sampling_fractions[p] = probability of including a species in the tree, if its proxy state is p
							reveal_fractions		= NULL,		# numeric vector of size NPstates, listing reveal fractions depending on state. reveal_fractions[p] = probability of knowing a tip's proxy state, if its proxy state is p
							coalescent 				= TRUE,
							as_generations			= FALSE,	# if FALSE, then edge lengths correspond to time. If TRUE, then edge lengths correspond to generations (hence if coalescent==false, all edges will have unit length).
							no_full_extinction		= TRUE,		# if true, then extinction of the entire tree is prevented. This is done by temporarily disabling extinctions when the number of extant tips is 1.
							Nsplits					= 2,	 	# number of children generated at each diversification event. If set to 2, a bifurcating tree is generated. If >2, the tree will be multifurcating.
							tip_basename			= "",		# basename for tips (e.g. "tip."). 
							node_basename			= NULL,		# basename for nodes (e.g. "node."). If NULL, then nodes will not have any labels.
							include_birth_times		= FALSE,
							include_death_times		= FALSE,
							include_rates			= FALSE,
							include_labels			= TRUE){	# whether to include tip-labels and node-labels as names in the returned state vectors (e.g. tip_states and node_states). Setting this to FALSE may slightly increase computational time/memory efficiency.
	# basic input checking
	if(is.null(max_tips) && is.null(max_time) && is.null(max_time_eq)) stop("ERROR: At least one of max_tips and/or max_time and/or max_time_eq must be non-NULL")
	is_hisse_model = !(is.null(NPstates) || (NPstates==0) || (NPstates==Nstates))
	if(is_hisse_model && is.null(proxy_map)) stop("ERROR: Missing proxy_map, needed for HiSSE model")
	if(is_hisse_model && (length(proxy_map)!=Nstates)) stop("ERROR: proxy_map has length %d, but should have length %d (Nstates)",length(proxy_map),Nstates)
	if(is_hisse_model && (length(unique(proxy_map))!=NPstates)) stop("ERROR: Not all %d proxy states are represented in proxy_map",NPstates)
	if((!is_hisse_model) && (!is.null(proxy_map)) & ((length(proxy_map)!=Nstates) || (any(proxy_map!=(1:Nstates))))) stop("ERROR: Non-trivial proxy_map contradicts non-HiSSE model")
	if(!is_hisse_model) NPstates = Nstates
	
	# set default model parameters
	parameter_names = c("transition_matrix_A", "transition_matrix_C", "birth_rates", "death_rates")
	if(is.null(parameters$transition_matrix_A)) parameters$transition_matrix_A = matrix(0,nrow=Nstates, ncol=Nstates);
	if(is.null(parameters$transition_matrix_C)) parameters$transition_matrix_C = diag(Nstates);
	# pc birth rates corresponding to each state
	if(is.null(parameters$birth_rates)){
		parameters$birth_rates = rep(1, times=Nstates);
	}else if(length(parameters$birth_rates)==1){
		parameters$birth_rates = rep(parameters$birth_rates, times=Nstates);
	}else if(length(parameters$birth_rates)!=Nstates){
		stop(sprintf("ERROR: Invalid number of birth_rates; expected %d, but got %d",Nstates,length(parameters$birth_rates)))
	}
	# pc death rates corresponding to each state
	if(is.null(parameters$death_rates)){
		parameters$death_rates = rep(1, times=Nstates);
	}else if(length(parameters$death_rates)==1){
		parameters$death_rates = rep(parameters$death_rates, times=Nstates);
	}else if(length(parameters$death_rates)!=Nstates){
		stop(sprintf("ERROR: Invalid number of death_rates; expected %d, but got %d",Nstates,length(parameters$death_rates)))
	}
	# start state
	if(is.null(start_state)){
		start_state = sample.int(n=Nstates,size=1)
	}else if(!(start_state %in% (1:Nstates))){
		stop(sprintf("ERROR: Invalid start_state (%d): Must be an integer between 1 and %d",start_state,Nstates))
	}
	# prepare sampling fractions
	if(is.null(sampling_fractions) || (length(sampling_fractions)==0)){
		sampling_fractions = rep(1,NPstates);
	}else if(length(sampling_fractions)==1){
		sampling_fractions = rep(sampling_fractions,NPstates);
	}else if(length(sampling_fractions)!=NPstates){
		stop(sprintf("ERROR: Invalid number of sampling fractions (%d), expected either 0, 1 or %d (NPstates)",length(sampling_fractions),NPstates))
	}
	# prepare reveal fractions  = probability of knowing a tip's state, depending on its actual state
	if(is.null(reveal_fractions) || (length(reveal_fractions)==0)){
		reveal_fractions = rep(1,NPstates);
	}else if(length(reveal_fractions)==1){
		reveal_fractions = rep(reveal_fractions,NPstates);
	}else if(length(reveal_fractions)!=NPstates){
		stop(sprintf("ERROR: Invalid number of reveal fractions (%d), expected either 0, 1 or %d (NPstates)",length(reveal_fractions),NPstates))
	}

	# check biological/physical validity of model parameters	
	if(any(sampling_fractions<=0 | sampling_fractions>1)) stop("ERROR: sampling_fractions must be between 0 (non-inclusive) and 1 (inclusive).")
	if(any(abs(rowSums(parameters$transition_matrix_A))>1e-6*max(abs(parameters$transition_matrix_A)))) stop("ERROR: Anagenetic transition rate matrix does not seem to be valid; some row sums are not zero.")
	if(any(parameters$transition_matrix_C<0)) stop("ERROR: Cladogenic transition probability matrix does not seem to be valid; some entries are negative.")
	if(any(abs(rowSums(parameters$transition_matrix_C)-1)>1e-6)) stop("ERROR: Cladogenic transition probability matrix does not seem to be valid; some row sums differ from 1.")


	# check if some passed parameters are not recognized
	invalids = setdiff(names(parameters),parameter_names)
	if(length(invalids)>0) stop(sprintf("ERROR: Unknown parameter '%s'",invalids[1]))

	
	results = generate_random_tree_Mk_rates_CPP(max_tips					= (if(is.null(max_tips)) -1 else max_tips),
												max_time					= (if(is.null(max_time)) -1 else max_time),
												max_time_since_equilibrium	= (if(is.null(max_time_eq)) -1 else max_time_eq),
												Nstates						= Nstates,
												state_birth_rates			= parameters$birth_rates, 
												state_death_rates			= parameters$death_rates,
												start_state					= max(1,min(Nstates, start_state)) - 1,
												transition_matrix_A			= as.vector(t(parameters$transition_matrix_A)), # flatten in row-major format
												transition_matrix_C			= as.vector(t(parameters$transition_matrix_C)), # flatten in row-major format
												coalescent					= coalescent,
												Nsplits						= Nsplits,
												as_generations				= as_generations,
												no_full_extinction			= no_full_extinction,
												include_birth_times			= include_birth_times,
												include_death_times			= include_death_times,
												include_rates				= include_rates);

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
	attr(tree,"order") = "none";
	
	
	# sub-sample tips (rarefy) if needed
	Nrarefied = 0;
	if(any(sampling_fractions<1)){
		if(length(unique(sampling_fractions))==1){
			keep_tips = sample.int(n=Ntips, size=sampling_fractions[1]*Ntips, replace=FALSE);
		}else{
			tip_pstates = results$clade_states[1:Ntips]+1L
			if(is_hisse_model) tip_pstates = proxy_map[tip_pstates]
			keep_tip = logical(Ntips)
			for(pstate in 1:NPstates){
				tips_with_pstate = which(tip_pstates==pstate)
				keep_tip[tips_with_pstate] = as.logical(rbinom(n=length(tips_with_pstate), size=1, prob=sampling_fractions[pstate])) # probability of keeping a tip in the tree is sampling_fractions[pstate]
			}
			keep_tips = which(keep_tip)
		}
		if(length(keep_tips)<2) return(list(success=FALSE, error=sprintf("Sampling fractions are too low for the generated tree (%d tips)",Ntips)))
		rarefaction = castor::get_subtree_with_tips(tree, only_tips=keep_tips, collapse_monofurcations=TRUE, force_keep_root=FALSE)
		tree 		= rarefaction$subtree
		Nrarefied 	= Ntips - length(tree$tip.label)
		results$root_time = results$root_time + rarefaction$root_shift; # update root time, in case root has changed
		if(include_rates){
			results$birth_rates_pc	= results$birth_rates_pc[rarefaction$new2old_clade]
			results$death_rates_pc	= results$death_rates_pc[rarefaction$new2old_clade]
		}
		Ntips 	= length(tree$tip.label)
		Nnodes 	= tree$Nnode
		if(!is.null(results$clade_states)) results$clade_states = results$clade_states[rarefaction$new2old_clade]
	}

	tip_states  = (if(is.null(results$clade_states)) NULL else results$clade_states[1:Ntips]+1L)
	node_states = (if(is.null(results$clade_states)) NULL else results$clade_states[(Ntips+1):(Ntips+Nnodes)]+1L)
	tip_proxy_states = NULL; node_proxy_states = NULL;
	if(is_hisse_model){
		if(!is.null(tip_states)) tip_proxy_states = proxy_map[tip_states]
		if(!is.null(node_states)) node_proxy_states = proxy_map[node_states]
	}

	# extract tip/node states & proxy states (if is_hisse_model)
	tip_states  = (if(is.null(results$clade_states)) NULL else results$clade_states[1:Ntips]+1L)
	node_states = (if(is.null(results$clade_states)) NULL else results$clade_states[(Ntips+1):(Ntips+Nnodes)]+1L)
	tip_proxy_states = NULL; node_proxy_states = NULL;
	if(is_hisse_model){
		if(!is.null(tip_states)) tip_proxy_states = proxy_map[tip_states]
		if(!is.null(node_states)) node_proxy_states = proxy_map[node_states]
	}
	
	# make some tip states (or proxy states, if is_hisse_model) unknown (i.e. assign state=NA)
	if((!is.null(tip_states)) && any(reveal_fractions<1)){
		tip_known = logical(Ntips)
		for(state in 1:NPstates){
			if(is_hisse_model){
				tips_with_state = which(tip_proxy_states==state)
			}else{
				tips_with_state = which(tip_states==state)
			}
			tip_known[tips_with_state] = as.logical(rbinom(n=length(tips_with_state), size=1, prob=reveal_fractions[state]))
		}
		if(is_hisse_model){
			tip_proxy_states[!tip_known] = NA
		}else{
			tip_states[!tip_known] = NA
		}		
	}
	
	# add labels to tip & node states
	if(include_labels){
		if(!is.null(tip_states)) names(tip_states) = tree$tip.label
		if((!is.null(node_states)) && (!is.null(tree$node.label))) names(node_states) = tree$node.label
		if(!is.null(tip_proxy_states)) names(tip_proxy_states) = tree$tip.label
		if((!is.null(node_proxy_states)) && (!is.null(tree$node.label))) names(node_proxy_states) = tree$node.label
	}
	
	return(list(success				= TRUE,
				tree				= tree,
				root_time			= results$root_time,
				final_time			= results$final_time,
				equilibrium_time	= results$equilibrium_time,
				Nbirths		 		= results$Nbirths,  	# number of birth events in each state
				Ndeaths				= results$Ndeaths,  	# number of death events in each state
				Ntransitions_A		= matrix(results$Ntransitions_A,ncol=Nstates,byrow=TRUE), # number of anagenetic transition events between each pair of states
				Ntransitions_C		= matrix(results$Ntransitions_C,ncol=Nstates,byrow=TRUE), # number of cladogenic transition events between each pair of states
				Nrarefied			= Nrarefied, # number of tips removed via rarefaction at the end
				tip_states			= tip_states,
				node_states			= node_states,
				tip_proxy_states	= tip_proxy_states,
				node_proxy_states	= node_proxy_states,
				start_state			= start_state,
				birth_times			= (if(include_birth_times) results$birth_times else NULL),
				death_times			= (if(include_death_times) results$death_times else NULL),
				clade_birth_rates	= (if(include_rates) results$birth_rates_pc else NULL),
				clade_death_rates	= (if(include_rates) results$death_rates_pc else NULL)));
	
}