# Plot a rooted phylogenetic tree in radial layout
plot_tree_radial = function(tree,							# rooted phylogenetic tree of class "phylo"
							file_path 		  	= NULL,		# optional path to save the plot to. Any existing file is silently replaced.
							make_ultrametric  	= FALSE,	# whether to make the tree ultrametric first, using a quick-and-dirty algorithm, so that all tips appear at the same distance from the center.
							symmetric_arcs		= TRUE,		# If TRUE, then nodes appear in the middle of their arc. If FALSE, nodes are radially placed in the middle of their descending tips, which could result in nodes not being in the center of their arc.
							opening				= 0,		# optional opening angle (in degrees), i.e., separating the two flanks of the tree. If 0, the tips will uniformly cover the whole circle. Must be between 0 and 360.
							rotate				= 0,		# optional angle (degrees) by which to rotate the whole plot, counter-clockwise
							plot_width		  	= 5,		# plot width & height in inches. Only relevant if file_path is specified.
							plot_title			= NULL,		# optional plot title
							show_scale_bar		= FALSE,
							align_tip_labels	= FALSE,	# whether to align all tip labels on a single circle, rather than placing each tip label immediately adjacent to each tip. Only relevant for non-ultrametric trees and only if tip_label_cex!=0. If colored rings are shown, then tip labels are always aligned.
							tip_cex				= 0,		# point size rescale for tips. If 0, no tip circles are shown.
							node_cex			= 0,		# point size rescale for nodes. If 0, no tip circles are shown.
							tip_label_cex		= NULL,		# size rescale for tip labels. If 0, no tip labels are shown. If NULL, this is chosen automatically to try to make tips fit.
							scale_edge_widths 	= FALSE,	# whether to scale the width (thickness) of each edge and node arc proportionally to the square root of the number of tips descending from it. If FALSE, all edges will have the same width.
							base_edge_width		= 1,		# base (smallest) width of edges, e.g. at the tips (and everywhere else if scale_edge_widths==FALSE)
							tip_color			= "black",	# either a single character or a character vector of length Ntips or an RGB tuple or an RGBA tuple, specifying tip colors & alpha.
							node_color			= "black",	# either a single character or a character vector of length Nnodes or an RGB tuple or an RGBA tuple, specifying node colors & alpha.
							edge_color			= "#00000080",	# either a single character or a character vector of length Nedges or an RGB tuple or an RGBA tuple, specifying edge colors & alpha.
							node_arc_color		= "#00000080",	# either a single character or a character vector of length Nodes or an RGB tuple or an RGBA tuple or a list of RGB/RGBA tuples, specifying node arc colors & alpha
							root_edge_color		= "#00000080",	# color of the root edge, if applicable
							tip_label_color		= "black",	# either a single character or a character vector of length Ntips or an RGB tuple or an RGBA tuple or a list of RGB/RGBA tuples, specifying tip label colors & alpha.
							ring_colors			= NULL,		# optional data.frame or matrix of size Ntips x NC, specifying one or more multi-colored rings to show around the tree. The entries of ring_colors[] must be valid HEX color characters, so ring_colors[t,r] is the color for tip t on ring r. A value of NULL or NA means no color segment is drawn for a tip, effectively keeping the background color. If ring_colors[] has row names, these will be matched to the tip labels, otherwise rows are assumed to already be synchronized with the tree's tips.
							ring_width			= 5,		# numeric, line width to use for colored rings around the tree. Only relevant if ring_colors is provided.
							ring_border_width	= 1,		# line thickness for the ring borders
							ring_border_color	= "#FFFFFF",# color to use for the ring borders, i.e., separating the rings from each other and from the tree
							legend_colors		= NULL,		# optional vector specifying colors to be defined in a legend. These can be for example colors in ring_colors[], or colors used for the tips & nodes & edges.
							legend_labels		= NULL){	# optional vector specifying labels to be given in a legend, synchronized with legend_colors[].
	# basic input checks
	if(is.null(legend_labels)!=is.null(legend_colors)) stop("legend_colors[] and legend_labels[] must either both be NULL, or both non-NULL.")
	if(is.null(plot_title)) plot_title = ""
	Nrings = (if(is.null(ring_colors)) 0 else ncol(ring_colors))
	
	# convert all user-supplied angles from degrees to radians
	opening = opening * (pi/180)
	rotate	= rotate * (pi/180)

	# compute the main geometry of the radial tree, e.g. polar coordinates of clades & edges
	if(make_ultrametric) tree = date_tree_red(tree, anchor_age=NULL)$tree
	Ntips	= length(tree$tip.label)
	Nnodes	= tree$Nnode
	Nclades	= Ntips+Nnodes
	geometry = get_radial_tree_geometry_CPP(Ntips			= Ntips,
											Nnodes			= Nnodes,
											Nedges			= nrow(tree$edge),
											tree_edge		= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
											edge_length		= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
											root_edge		= (if(is.null(tree$root.edge)) 0 else tree$root.edge),
											symmetric_arcs	= symmetric_arcs,
											opening			= opening)
	geometry$root = geometry$root + 1 # convert to 1-based indexing
	
	# rotate the entire tree if needed
	if(rotate!=0){
		geometry$edge_phi 		= geometry$edge_phi + rotate
		geometry$clade_phi 		= geometry$clade_phi + rotate
		geometry$node_arc_phi0 	= geometry$node_arc_phi0 + rotate
		geometry$node_arc_phi1 	= geometry$node_arc_phi1 + rotate
		geometry$clade_phi0 	= geometry$clade_phi0 + rotate
		geometry$clade_phi1 	= geometry$clade_phi1 + rotate
	}
	
	# instantiate plot file with proper dimensions and margins
	if(is.null(tip_label_cex)) tip_label_cex = min(1, 3*plot_width*max(0,2*pi-opening)*mean(geometry$clade_rho[1:Ntips])/(Ntips*geometry$max_rho))
	margins_chars = c(1+2*show_scale_bar, # bottom
				1, # left
				2+(length(strsplit(plot_title, "\n")[[1]])), # top
				1+(if(is.null(legend_labels)) 0 else 4+0.4*max(sapply(legend_labels, FUN=nchar)))) # right
	if(tip_label_cex!=0) margins_chars = margins_chars + tip_label_cex*0.5*max(sapply(tree$tip.label, FUN=nchar)) # increase margin to accommodate tip labels
	graphics_params = get_dummy_graphics_params()
	margins_inches  = margins_chars * graphics_params$csi
	if(!is.null(file_path)) initiate_plot_file(file_path=file_path, width=plot_width+margins_inches[2]+margins_inches[4], height=plot_width+margins_inches[1]+margins_inches[3], dpi=200)
	graphics::par(mar = margins_chars, xpd = NA)

	# establish basic plot region
	graphics::plot(	x 		= 0,
					y 		= 0,
					type 	= "n",
					main 	= plot_title,
					xlab	= "",
					ylab	= "",
					asp		= 1,
					axes	= FALSE,
					xlim	= c(-1.1*geometry$max_rho, 1.1*geometry$max_rho),
					ylim	= c(-1.1*geometry$max_rho, 1.1*geometry$max_rho))
	lw2coordinates = (graphics::grconvertX(1, "inches", "user") - graphics::grconvertX(0, "inches", "user"))/96 # conversion factor between line widths (e.g., used for lwd) and user coordinates

	# show edges as lines
	if(scale_edge_widths) mass2width = function(mass){ base_edge_width*(1+4*sqrt((mass-1)/(Ntips-1))); }
	graphics::segments(	x0	= geometry$edge_rho0*cos(geometry$edge_phi),
						y0	= geometry$edge_rho0*sin(geometry$edge_phi),
						x1	= geometry$edge_rho1*cos(geometry$edge_phi),
						y1	= geometry$edge_rho1*sin(geometry$edge_phi),
						col	= edge_color,
						lwd	= (if(scale_edge_widths) mass2width(geometry$edge_mass) else base_edge_width),
						lend= "butt")

	# show node arcs, i.e., connecting all child edges of a node
	for(node in seq_len(Nnodes)){
		plot_arc(phi0	= geometry$node_arc_phi0[node], 
				 phi1	= geometry$node_arc_phi1[node],
				 rho	= geometry$clade_rho[Ntips+node],
				 col	= (if(length(node_arc_color)==1) node_arc_color else node_arc_color[[node]]),
				 lwd	= (if(scale_edge_widths) mass2width(geometry$node_arc_mass[node]) else base_edge_width),
				 lend	= "square")
	}

	# show root edge
	if(!is.null(tree$root.edge)){
		graphics::segments(	x0	= 0,
							y0	= 0,
							x1	= geometry$clade_rho[geometry$root]*cos(geometry$clade_phi[geometry$root]),
							y1	= geometry$clade_rho[geometry$root]*sin(geometry$clade_phi[geometry$root]),
							col	= root_edge_color,
							lwd	= (if(scale_edge_widths) mass2width(geometry$node_arc_mass[geometry$root-Ntips]) else base_edge_width),
							lend= "butt")
	}

	# show tips as points
	if(tip_cex>0){
		graphics::points(x 		= geometry$clade_rho[1:Ntips]*cos(geometry$clade_phi[1:Ntips]),
						y 		= geometry$clade_rho[1:Ntips]*sin(geometry$clade_phi[1:Ntips]),
						type 	= "p",
						col 	= tip_color,
						pch		= 19,
						cex		= tip_cex)
	}

	# show nodes as points
	if(node_cex>0){
		graphics::points(x 		= geometry$clade_rho[(Ntips+1):Nclades]*cos(geometry$clade_phi[(Ntips+1):Nclades]),
						y 		= geometry$clade_rho[(Ntips+1):Nclades]*sin(geometry$clade_phi[(Ntips+1):Nclades]),
						type 	= "p",
						col 	= node_color,
						pch		= 19,
						cex		= node_cex)
	}

	# show ring colors
	if(Nrings>0){
		if(!is.null(rownames(ring_colors))){
			tip2row = match(tree$tip.label,rownames(ring_colors))
			if(any(is.na(tip2row))) stop(sprintf("%d out of %d tree tips did not match any row names in ring_colors[], for example '%s'.",sum(is.na(tip2row)),Ntips,tree$tip.label[is.na(tip2row)][[1]]))
			ring_colors = ring_colors[tip2row,,drop=FALSE]
		}
		ring_border_rhos = c(geometry$max_rho+lw2coordinates*ring_border_width*2)
		for(ring in seq_len(Nrings)){
			ring_rho = geometry$max_rho + lw2coordinates*(ring_border_width*2 + ring_width/2 + (ring-1)*(ring_width+ring_border_width))
			for(tip in seq_len(Ntips)){
				color = ring_colors[tip,ring]
				if(is.null(color) || is.na(color)) next # skip this arc
				plot_arc(phi0	= geometry$clade_phi0[tip], 
						 phi1	= geometry$clade_phi1[tip],
						 rho	= ring_rho,
						 col	= color,
						 lwd	= ring_width,
						 lend	= "butt")
	
			}
			ring_border_rhos = c(ring_border_rhos, ring_rho + lw2coordinates*(ring_width/2+ring_border_width/2))
		}
		# add ring borders. We do this at the end to avoid covering borders up with subsequent rings
		for(ring_border_rho in ring_border_rhos){
			plot_arc(phi0	= geometry$clade_phi0[geometry$root],
					 phi1	= geometry$clade_phi1[geometry$root],
					 rho	= ring_border_rho,
					 col	= ring_border_color,
					 lwd	= ring_border_width)
		}
		geometry$max_rho = max(ring_rho + lw2coordinates*(ring_width/2 + ring_border_width), geometry$max_rho) # update max_rho to the end of the outer-most ring
	}

	#	show tip labels
	if(tip_label_cex!=0){
		max_label_rho_outer = 0 # keep track of the maximum radius increase due to a tip label
		for(tip in seq_len(Ntips)){
			label_srt = (180/pi)*geometry$clade_phi[tip]
			label_adj = c(0, 0.5)
			if((label_srt>90) && (label_srt<270)){
				label_srt = label_srt-180
				label_adj = c(1,0.5)
			}
			label_rho = (if(align_tip_labels || (Nrings>0)) geometry$max_rho else geometry$clade_rho[tip]) + tip_label_cex*graphics::par("cxy")[1]*0.7
			graphics::text(	x 		= label_rho * cos(geometry$clade_phi[tip]),
							y 		= label_rho * sin(geometry$clade_phi[tip]),
							labels 	= tree$tip.label[tip],
							col		= (if(length(tip_label_color)==1) tip_label_color else tip_label_color[[tip]]),
							cex		= tip_label_cex,
							srt 	= label_srt,
							adj 	= label_adj)
			max_label_rho_outer = max(max_label_rho_outer, label_rho + tip_label_cex*0.6*graphics::par("cxy")[1]*nchar(tree$tip.label[tip]))
		}
		geometry$max_rho = max_label_rho_outer # update max_rho to account for tip labels
	}

	# show scale bar
	if(show_scale_bar){
		bar_length = get_nice_scale_bar(geometry$max_rho*0.5)
		scale_bar_y = -1.05*geometry$max_rho
		graphics::segments(	x0	= geometry$max_rho - bar_length,
							y0	= scale_bar_y,
							x1	= geometry$max_rho,
							y1	= scale_bar_y,
							lwd= 2)
		graphics::text(	x 		= geometry$max_rho - (bar_length/2), 
						y		= scale_bar_y - 5*lw2coordinates,
						adj		= c(0.5,1),
						labels	= format(bar_length, digits=4, scientific=FALSE))
	}
	
	# show legend
	if(!is.null(legend_labels)){
		graphics::legend(x		= 1.05 * geometry$max_rho,
						y		= geometry$max_rho,
						legend	= legend_labels,
						col		= legend_colors,
						xjust	= 0,
						yjust	= 1,
						lwd 	= 4,
						xpd		= TRUE)
	}

	if(!is.null(file_path)) invisible(grDevices::dev.off())
}


#######################################################
# Auxiliary functions

# get a "nice" scale bar length, i.e., rounded to the nearest value from the sequence ..., 0.05, 0.1, 0.5, 1, 5, 10, 50, ...
get_nice_scale_bar = function(target_length){
	base = 10^floor(log10(target_length))
	candidates = c(1, 5, 10) * base
	return(candidates[which.min(abs(candidates - target_length))])
}


plot_arc = function(phi0, phi1, rho, col="black", lwd=1, lend="round"){
	phis = seq(phi0, phi1, length.out = max(10, as.integer((phi1-phi0)/0.02)))
	graphics::lines(x	= rho * cos(phis), 
					y	= rho * sin(phis), 
					col = col, 
					lwd	= lwd,
					lend= lend)
}


get_dummy_graphics_params = function(){
	grDevices::pdf(NULL)
	params = list(csi=graphics::par("csi"))
	invisible(grDevices::dev.off())
	return(params)
}

