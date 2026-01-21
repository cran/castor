# Hidden state prediction of discrete states for tips & nodes on a tree, based solely on the tree-wide ("global") proportions of known states but ignoring any phylogenetic relationships
# Used mostly as a neutral/null model for comparisons.
# The function returns likelihoods for all tips & nodes reflecting the global empirical state proportions, as well as randomly drawn states for tips & nodes drawn from the same empirical distribution
hsp_naive = function(tree, 
					tip_states, 	# integer vector of size Ntips. Unknown states are encoded as NA.
					Nstates		= NULL, 
					check_input	= TRUE){
	Ntips   = length(tree$tip.label)
	Nclades = Ntips + tree$Nnode
	
	# basic error checking
	if(length(tip_states)!=Ntips) stop(sprintf("ERROR: Length of tip_states (%d) is not the same as the number of tips in the tree (%d)",length(tip_states),Ntips));
	if(!is.numeric(tip_states)) stop(sprintf("ERROR: tip_states must be integers"))
	if(check_input){
		if((!is.null(names(tip_states))) && any(names(tip_states)!=tree$tip.label)) stop("ERROR: Names in tip_states and tip labels in tree don't match (must be in the same order).")
	}
	
	# find known_tips and determine known state proportions
	known_tips = which(!is.na(tip_states))
	if(length(known_tips)==0) stop("ERROR: All tip states are hidden")
	known_tip_states = tip_states[known_tips]
	if(is.null(Nstates)) Nstates = max(known_tip_states)
	if(check_input){
		min_tip_state = min(known_tip_states)
		max_tip_state = max(known_tip_states)
		if((min_tip_state<1) || (max_tip_state>Nstates)) stop(sprintf("ERROR: Non-NA tip_states must be integers between 1 and %d, but found values between %d and %d",Nstates,min_tip_state,max_tip_state))
	}
	if(Nstates==1) return(list(	success		= TRUE, 
								likelihoods	= rep(1, Nclades),
								states		= rep(1,Nclades)))
	proportions = sapply(seq_len(Nstates), FUN=function(state){mean(known_tip_states==state)})
	
	# infer likelihoods & draw random states for all tips & nodes based on the empirical proportions
	likelihoods = matrix(rep(proportions, times=Nclades), nrow=Nclades, byrow = TRUE)
	likelihoods[known_tips, ] = 0.0
	likelihoods[cbind(known_tips, known_tip_states)] = 1.0
	states = sapply(seq_len(nrow(likelihoods)), FUN=function(r){ sample.int(ncol(likelihoods), size=1, prob=likelihoods[r,]) })
		
	return(list(success 	= TRUE,
				likelihoods	= likelihoods,
				states		= states))
}


