# Fit a Brownian Motion model of continuous trait evolution for one or more traits
# If more than one trait is modelled, then this function estimates the full diffusivity matrix D (=2*sigma^2)
# 	D is a non-negative definite symmetric matrix of size Ntraits x Ntraits, 
# 	such that exp(-X^T*D^{-1}*X/(4*L))/sqrt(det(2*pi*D)) is the probability density for the multidimensional trait vector X after phylogenetic distance L, if initially located at the origin.
fit_bm_model = function(tree, 
						tip_states,				# numeric vector of size Ntips
						check_input	= TRUE){
	Ntips  	= length(tree$tip.label);
	Nnodes	= tree$Nnode;
	scalar  = is.vector(tip_states)
	Ntraits	= (if(scalar) 1 else ncol(tip_states))
	
	# calculate independent contrasts (as independent increments of a Brownian motion process)
	pic_results = get_independent_contrasts(tree, tip_states, scaled = TRUE, only_bifurcations = FALSE, check_input = check_input);
	Npics = (if(scalar) length(pic_results$PICs) else nrow(pic_results$PICs))
	X = pic_results$PICs
	
	# estimate diffusivity based on vectorial increments
	# maximum-likelihood estimator on the intrinsic geometry of positive-definite matrices
	if(scalar){
		diffusivity 	= 0.5 * mean(X**2)
		loglikelihood 	= -0.25 * sum((X**2)/diffusivity) - 0.5*Npics*log(2*pi) - 0.5*sum(log(2*pic_results$distances)) - 0.5*Npics*log(diffusivity)
	}else{
		diffusivity = matrix(0, ncol=Ntraits, nrow=Ntraits);
		for(t1 in 1:Ntraits){
			for(t2 in t1:Ntraits){
				diffusivity[t1,t2] = 0.5 * mean(X[,t1]*X[,t2])
				diffusivity[t2,t1] = diffusivity[t1,t2]; 
			}
		}
		inverse_diffusivity = solve(diffusivity)
		loglikelihood = -0.25 * sum(sapply(1:Npics, function(p) X[p,,drop=FALSE] %*% inverse_diffusivity %*% t(X[p,,drop=FALSE]))) - 0.5*Ntraits*Npics*log(2*pi) - 0.5*Ntraits*sum(log(2*pic_results$distances)) - 0.5*Npics*determinant(diffusivity,logarithm=TRUE)$modulus
	}
	
	return(list(diffusivity		= diffusivity, 
				loglikelihood	= loglikelihood))
}