# reroot a non-dated tree based on tip sampling times, by optimizing the root-to-tip (RTT) regression goodness of fit (regression of tip times vs phylodistances-from-root)
# The precise objective optimized can be "correlation", "R2" or "SSR" (sum of squared residuals)
# The input tree's edge lengths should be measured in substitutions per site, and tip sampling times should be measured in forward time (i.e., younger tips have a greater time).
# Adjusted (and improved) from ape::rtt v5.7-1.
root_via_rtt = function(tree, 								# tree of class "phylo". The tree may be rooted or unrooted, and the root placement does not matter.
						tip_times,							# numeric vector of length Ntips, listing the sampling times of all tips.
						objective 			= "R2",
						force_positive_rate = FALSE, 		# (logical) whether to force the implied mutation rate to be >=0, by constraining the root placements.
						Nthreads 			= 1,			# (integer) number of threads to use in parallel, where applicable
						optim_algorithm		= "nlminb",	# (character) either "optimize" or "nlminb". What algorithm to use for fitting.
						relative_error 		= 1e-9){
	if(objective == "correlation")
		objective_function = function(x, y){
			return(stats::cor(y, x))
		}
	else if(objective == "R2")
		objective_function = function(x, y){
			fit = stats::lm.fit(cbind(rep(1,Ntips),x), y)
			if(force_positive_rate && (fit$coefficients[2]<0)) return(-1e50)
			R2 = 1 - mean(fit$residuals^2)/stats::var(y)
			return(R2)
		}
	else if(objective == "SSR"){
		objective_function = function(x, y){
			fit = stats::lm.fit(cbind(rep(1,Ntips),x), y)
			if(force_positive_rate && (fit$coefficients[2]<0)) return(-1e50)
			return(-sum(fit$residuals^2) / fit$df.residual)
		}
	}

	Ntips 			= length(tree$tip.label)
	tree  			= ape::unroot(tree)
	phylodistances 	= castor::get_all_pairwise_distances(tree=tree, check_input=FALSE)[, 1:Ntips] # get all pairwise phylodistances between clades & tips

	aux_objective_function_at_edge = function(x, parent, child){
		if(!is.finite(x)) return(-1e50)
		edge_distaces = x * phylodistances[parent, ] + (1 - x) * phylodistances[child,]
		return(objective_function(x=tip_times, y=edge_distaces))
	}
	
	aux_fit_single_edge = function(e){
		if(optim_algorithm == "optimize"){
			fit = stats::optimize(f		 = function(x) aux_objective_function_at_edge(x, tree$edge[e,1], tree$edge[e,2]), 
								interval = c(0, 1), 
								maximum  = TRUE, 
								tol 	 = relative_error)
		}else if(optim_algorithm == "nlminb"){
			fit = stats::nlminb(start		= 0.5, 
								objective	= function(x) -aux_objective_function_at_edge(x, tree$edge[e,1], tree$edge[e,2]),
								lower 		= 0,
								upper 		= 1,
								control 	= list(iter.max = 10000, eval.max = 10000, rel.tol = relative_error))
			fit$objective = -fit$objective
		}
		return(fit$objective)
	}

	if(Nthreads>1){
		obj_edge = unlist(parallel::mclapply(1:nrow(tree$edge), FUN=aux_fit_single_edge, mc.cores=Nthreads, mc.preschedule = TRUE, mc.cleanup = TRUE))
	}else{
		obj_edge = sapply(1:nrow(tree$edge), FUN=aux_fit_single_edge)
	}
		
	valid_edges = which(is.finite(obj_edge))
	best.edge = valid_edges[which.max(obj_edge[valid_edges])]

	best.edge.parent = tree$edge[best.edge, 1]
	best.edge.child  = tree$edge[best.edge, 2]
	best.edge.length = tree$edge.length[best.edge]

	opt.fun  = function(x) aux_objective_function_at_edge(x, best.edge.parent, best.edge.child)
	best.pos = optimize(opt.fun, c(0, 1), maximum = TRUE, tol = relative_error)$maximum

	new_root = list(edge = matrix(c(2L, 1L), 1, 2), tip.label = "new_root", edge.length = 1, Nnode = 1L, root.edge = 1)
	class(new_root) = "phylo"
	tree = ape::bind.tree(tree, new_root, where = best.edge.child, position = best.pos*best.edge.length)
	tree = ape::collapse.singles(tree)
	tree = ape::root(tree, "new_root")
	tree = ape::drop.tip(tree, "new_root")
	return(list(tree=tree))
}