# date a tree (make ultrametric) using penalized likelihood
date_tree_PL = function(	tree,		# rooted tree
							lambda,		# smoothing parameter
							min_age			= 1,
							max_age			= NULL,
							node			= "root",
							Nsites			= 1,		# number of sites in the sequences. Set to 1 if treelogenetic distances are in mean number of subtitotions.
							tolerance		= 1e-8,
							optim_control	= list(),	# a named list containing options for the optim fitting routine (e.g. maxit)
							Ntrials			= 1,
							Nthreads		= 1){
	Ntips 	= length(tree$tip.label)
	Nnodes	= tree$Nnode
	Nedges	= nrow(tree$edge)
	
	ROOT = Ntips + 1L
	if (identical(node, "root")) node = ROOT
	if (any(node <= Ntips)) stop("node numbers should be greater than the number of tips")
	zerobl = which(tree$edge.length <= 0)
	if (length(zerobl)) {
		if (any(tree$edge[zerobl, 2] <= Ntips)) {
			stop("at least one terminal branch is of length zero:\n  you should remove it to have a meaningful estimation.")
		}else {
			warning("at least one internal branch is of length zero:\n  it was collapsed and some nodes have been deleted.")
			if (length(node) == 1 && node == ROOT) {
				tree = castor::merge_short_edges(tree)
			}else {
				tmp = FALSE
				if (is.null(tree$node.label)) {
				  tmp = !tmp
				  tree$node.label = paste("node", 1:tree$Nnode)
				}
				node.lab = tree$node.label[node - Ntips]
				tree = castor::merge_short_edges(tree)
				node = match(node.lab, tree$node.label) + Ntips
				if (tmp) tree$node.label = NULL
			}
		}
	}
	el 			= tree$edge.length
	e1 			= tree$edge[, 1L]
	e2 			= tree$edge[, 2L]
	Nedges		= length(e1)
	TIPS 		= 1:Ntips
	EDGES 		= 1:Nedges
	ini.rate 	= pmax(el,tolerance) # make sure ini.rate is larger than the tolerance threshold
	el 			= el/Nsites
	basal 		= which(e1 == ROOT)
	Nbasal 		= length(basal)
	ind 		= matrix(0L, Nedges - Nbasal, 2)
	ind[, 1] 	= EDGES[-basal]
	ind[, 2] 	= match(e1[EDGES[-basal]], e2)
	age 		= numeric(Ntips + Nnodes)
    seq.nod		= get_paths_root_to_tips(tree)
	ini.time 	= age
	ini.time[ROOT:(Ntips + Nnodes)] = NA
	ini.time[node] = if (is.null(max_age)) 
		min_age
	else (min_age + max_age)/2
	if (is.na(ini.time[ROOT])) {
		ini.time[ROOT] = if (is.null(max_age)) 3 * max(min_age) else 3 * max(max_age)
	}
	ISnotNA.ALL = unlist(lapply(seq.nod, function(x) sum(!is.na(ini.time[x]))))
	o = order(ISnotNA.ALL, decreasing = TRUE)
	for (y in seq.nod[o]) {
		ISNA = is.na(ini.time[y])
		if (any(ISNA)) {
			i = 2L
			while (i <= length(y)) {
				if (ISNA[i]) {
				  j = i + 1L
				  while (ISNA[j]) j = j + 1L
				  nb.val = j - i
				  by = (ini.time[y[i - 1L]] - ini.time[y[j]])/(nb.val + 1)
				  ini.time[y[i:(j - 1L)]] = ini.time[y[i - 1L]] - by * seq_len(nb.val)
				  i = j + 1L
				}
				else i = i + 1L
			}
		}
	}
	real.edge.length = ini.time[e1] - ini.time[e2]
	if (any(real.edge.length <= 0))  stop("some initial branch lengths are zero or negative;\n  maybe you need to adjust the given dates -- see '?chronopl' for details")
	node.bak = node
	unknown.ages = Ntips + 1:Nnodes
	lower = rep(tolerance, length(unknown.ages))
	upper = rep(1/tolerance, length(unknown.ages))
	if (!is.null(max_age)) {
		lower[node - Ntips] = min_age
		upper[node - Ntips] = max_age
		interv = which(min_age != max_age)
        if(length(interv)>0){
			node = node[-interv]
			if (length(node)) age[node] = min_age[-interv]
		}else{
			age[node] = min_age
		}
	}else{
		age[node] = min_age
	}
	if (length(node)) {
		unknown.ages = unknown.ages[Ntips - node]
		lower = lower[Ntips - node]
		upper = upper[Ntips - node]
	}
	known.ages = c(TIPS, node)
	lower = c(rep(tolerance, Nedges), lower)
	upper = c(rep(1 - tolerance, Nedges), upper)
	
    # prepare auxiliary lookup lists
    adjacent_edges = get_adjacent_edges_per_edge(tree)
	
	minusploglik.gr = function(rate, node.time) {
		grad = numeric(Nedges + length(unknown.ages))
		age[unknown.ages] = node.time
		real.edge.length = age[e1] - age[e2]
		if (any(real.edge.length < 0)) {
			grad[] = 0
			return(grad)
		}
		grad[EDGES] = real.edge.length - el/rate
		if (Nbasal == 2) {
			grad[basal[1]] = grad[basal[1]] + lambda * (rate[basal[1]] - rate[basal[2]])
			grad[basal[2]] = grad[basal[2]] + lambda * (rate[basal[2]] - rate[basal[1]])
		} else {
			for (i in 1:Nbasal) grad[basal[i]] = grad[basal[i]] + lambda * (2 * rate[basal[i]] * (1 - 1/Nbasal) - 2 * sum(rate[basal[-i]])/Nbasal)/(Nbasal - 1)
		}
		for (i in EDGES) {
			ii = c(which(e2 == e1[i]), which(e1 == e2[i]))
			if (!length(ii)) next
			grad[i] = grad[i] + lambda * (2 * length(ii) * rate[i] - 2 * sum(rate[ii]))
		}
		for (i in 1:length(unknown.ages)) {
			nd = unknown.ages[i]
			ii = which(e1 == nd)
			grad[i + Nedges] = sum(rate[ii] - el[ii]/real.edge.length[ii])
			if (nd != ROOT) {
				ii = which(e2 == nd)
				grad[i + Nedges] = grad[i + Nedges] - rate[ii] + el[ii]/real.edge.length[ii]
			}
		}
		grad
	}
	
	minusploglik = function(rate, node.time) {
		age[unknown.ages] = node.time
		real.edge.length = age[e1] - age[e2]
		if (any(real.edge.length < 0)) return(1e+50)
		B = rate * real.edge.length
		loglik = sum(-B + el * log(B) - lfactorial(el))
		value = -(loglik - lambda * (sum((rate[ind[, 1]] - rate[ind[, 2]])^2) + stats::var(rate[basal])))
		return(value)
	}
	
#	 out = nlminb(start		= c(ini.rate, ini.time[unknown.ages]), 
#	 			objective 	= function(p) minusploglik(p[EDGES], p[-EDGES]), 
#	 			gradient 	= function(p) minusploglik.gr(p[EDGES], p[-EDGES]), 
# 				control 	= list(eval.max = eval.max, iter.max = iter.max, ...), 
# 				lower 		= lower, 
# 				upper 		= upper)
#	 attr(tree, "ploglik") = -out$objective
#	 attr(tree, "rates") = out$par[EDGES]
#	 attr(tree, "message") = out$message
#	 age[unknown.ages] = out$par[-EDGES]
	out = optim(par		= c(ini.rate, ini.time[unknown.ages]), 
				fn		= function(p) minusploglik(p[EDGES], p[-EDGES]),
				gr		= function(p) minusploglik.gr(p[EDGES], p[-EDGES]),
				method	= "L-BFGS-B",
				control	= optim_control,
				lower	= lower,
				upper	= upper)
	attr(tree, "ploglik") = -out$value
	attr(tree, "rates") = out$par[EDGES]
	attr(tree, "message") = out$message
	age[unknown.ages] = out$par[-EDGES]
	
	tree$edge.length = age[e1] - age[e2]
	return(list(dated_tree=tree))
}
