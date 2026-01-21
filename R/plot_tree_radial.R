# Plot a rooted phylogenetic tree in radial layout
plot_tree_radial = function(tree,							# rooted phylogenetic tree of class "phylo"
							file_path 		  		= NULL,		# optional path to save the plot to. Any existing file is silently replaced.
							make_ultrametric  		= FALSE,	# whether to make the tree ultrametric first, using a quick-and-dirty algorithm, so that all tips appear at the same distance from the center.
							symmetric_arcs			= TRUE,		# If TRUE, then nodes appear in the middle of their arc. If FALSE, nodes are radially placed in the middle of their descending tips, which could result in nodes not being in the center of their arc.
							opening					= 0,		# optional opening angle (in degrees), i.e., separating the two flanks of the tree. If 0, the tips will uniformly cover the whole circle. Must be between 0 and 360.
							rotate					= 0,		# optional angle (degrees) by which to rotate the whole plot, counter-clockwise
							plot_width		  		= 5,		# plot width & height in inches. Only relevant if file_path is specified.
							plot_title				= NULL,		# optional plot title
							show_scale_bar			= FALSE,
							align_tip_labels		= FALSE,	# whether to align all tip labels on a single circle, rather than placing each tip label immediately adjacent to each tip. Only relevant for non-ultrametric trees and only if tip_label_cex!=0. If colored rings are shown, then tip labels are always aligned.
							tip_cex					= 0,		# point size rescale for tips. If 0, no tip circles are shown.
							node_cex				= 0,		# point size rescale for nodes. If 0, no tip circles are shown.
							tip_label_cex			= NULL,		# size rescale for tip labels. If 0, no tip labels are shown. If NULL, this is chosen automatically to try to make tips fit.
							scale_edge_widths 		= FALSE,	# whether to scale the width (thickness) of each edge and node arc proportionally to the square root of the number of tips descending from it. If FALSE, all edges will have the same width.
							base_edge_width			= 1,		# base (smallest) width of edges, e.g. at the tips (and everywhere else if scale_edge_widths==FALSE)
							tip_color				= "black",	# either a single character or a character vector of length Ntips or an RGB tuple or an RGBA tuple, specifying tip colors & alpha.
							node_color				= "black",	# either a single character or a character vector of length Nnodes or an RGB tuple or an RGBA tuple, specifying node colors & alpha.
							edge_color				= "#00000080",	# either a single character or a character vector of length Nedges or an RGB tuple or an RGBA tuple, specifying edge colors & alpha.
							node_arc_color			= "#00000080",	# either a single character or a character vector of length Nodes or an RGB tuple or an RGBA tuple or a list of RGB/RGBA tuples, specifying node arc colors & alpha
							root_edge_color			= "#00000080",	# color of the root edge, if applicable
							tip_label_color			= "black",	# either a single character or a character vector of length Ntips or an RGB tuple or an RGBA tuple or a list of RGB/RGBA tuples, specifying tip label colors & alpha.
							ring_colors				= NULL,		# optional data.frame or matrix of size Ntips x NC, specifying one or more multi-colored rings to show around the tree. The entries of ring_colors[] must be valid HEX color characters, so ring_colors[t,r] is the color for tip t on ring r. A value of NULL or NA means no color segment is drawn for a tip, effectively keeping the background color. If ring_colors[] has row names, these will be matched to the tip labels, otherwise rows are assumed to already be synchronized with the tree's tips.
							ring_width				= 5,		# numeric, line width to use for colored rings around the tree. Only relevant if ring_colors is provided.
							ring_border_width		= 1,		# line thickness for the ring borders
							ring_border_color		= "#FFFFFF",# color to use for the ring borders, i.e., separating the rings from each other and from the tree
							dotstack_heights		= NULL,		# optional data.frame or matrix of size Ntips x ND, specifying one or more dot-stacks to show around the tree, i.e. for visualizing discrete quantitative data. The entries of dotstack_heights[] must be non-negative integers, so dotstack_heights[t,d] is the number of dots to stack for tip t in the d'th dot-stack. If dotstack_heights[] has row names, these will be matched to the tip labels, otherwise rows are assumed to already be synchronized with the tree's tips.
							dotstack_color			= "black",	# either a single color or a vector colors, specifying the dot colors to use for the dot-stacks. Colors are recycled if needed.
							dotstack_border_width	= 1,		# line thickness for the dotstack borders. If NULL, this is set equal to ring_border_width.
							dotstack_border_color	= "black",	# color to use for the dotstack borders, i.e., separating the dotstacks from each other. If NULL, this is set to ring_border_color.
							tip_sectors				= NULL,		# optional vector of length Ntips, specifying a sector ID for each tip. IDs may be characters, integers, or some other reasonable scalar. Tips will be grouped by sector, and each such group will be visually highlighted using an outer arc and label. This may be used e.g. to highlight major phyla in a tree. tip_sectors[] may also contain NA or NULL, representing sector-less tips not to be highlighted. It is expected that all or most sectors are monophyletic (but see option fragmented_sectors).
							sector_color			= "black",	# color to use for sector arcs. Only relevant if tip_sectors is not NULL.
							sector_width			= 2,		# line thickness for the sector arcs
							sector_gap				= 15,		# gap between the sector arcs and the inner-next element (typically tip labels), measured in the same units as linewidths
							min_sector_size			= 2,		# minimum number of contiguous tips per displayed sector (after splitting, if fragmented_sectors=="split").
							sector_label_cex		= 1,		# size rescale for sector labels. If 0, no sector labels are shown.
							sector_label_color		= "black",	# color specification for sector labels
							fragmented_sectors		= "error",	# how to handle fragmented (i.e., non-monophyletic) tip sectors. Options include "error", "omit" or "split".
							legend_colors			= NULL,		# optional vector specifying colors to be defined in a legend. These can be for example colors in ring_colors[], or colors used for the tips & nodes & edges.
							legend_labels			= NULL){	# optional vector specifying labels to be given in a legend, synchronized with legend_colors[].
	# basic input checks
	if(is.null(legend_labels)!=is.null(legend_colors)) stop("ERROR: legend_colors[] and legend_labels[] must either both be NULL, or both non-NULL.")
	if((!is.null(legend_labels)) && (length(legend_labels)!=length(legend_colors))) stop("ERROR: legend_colors[] and legend_labels[] have different lengths")
	if(is.null(plot_title)) plot_title = ""
	if(is.null(dotstack_border_width)) dotstack_border_width = ring_border_width
	if(is.null(dotstack_border_color)) dotstack_border_color = dotstack_border_color
	Nrings = (if(is.null(ring_colors)) 0 else ncol(ring_colors))
	NDS	   = (if(is.null(dotstack_heights)) 0 else ncol(dotstack_heights))
	
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
	angular_tip_order = order(geometry$clade_phi[1:Ntips]) # integer vector of length Ntips, specifying the geometric order of tip placement along the tree's perimeter, i.e., in order of increasing angle. The first & last tips mark the tree's 'opening'. This is NOT the order in which tips are indexed in the tree or in the geometry object.
	
	# rotate the entire tree if needed
	if(rotate!=0){
		geometry$edge_phi 		= geometry$edge_phi + rotate
		geometry$clade_phi 		= geometry$clade_phi + rotate
		geometry$node_arc_phi0 	= geometry$node_arc_phi0 + rotate
		geometry$node_arc_phi1 	= geometry$node_arc_phi1 + rotate
		geometry$clade_phi0 	= geometry$clade_phi0 + rotate
		geometry$clade_phi1 	= geometry$clade_phi1 + rotate
	}
	
	# parse and check validity of tip sectors
	if(!is.null(tip_sectors)){
		if(!is.null(names(tip_sectors))){
			tip2i = match(tree$tip.label,names(tip_sectors))
			if(any(is.na(tip2i))) stop(sprintf("%d out of %d tree tips did not match any names in tip_sectors[], for example '%s'.",sum(is.na(tip2i)),Ntips,tree$tip.label[is.na(tip2i)][[1]]))
			tip_sectors = tip_sectors[tip2i]
		}
		Nsectors		= 0
		sector_labels 	= character(0) # name of each sector. A name may in principle appear multiple times here, if fragmented_sectors=="split".
		sector_firsts 	= integer(0)   # index of the first tip in each sector. Note that not all tips between sector_firsts[s] and sector_lasts[s] need be in sector s, since the order in which tips are stored in the tree does not necessarily reflect their visual order on the tree's perimeter.
		sector_lasts	= integer(0)   # index of the last tip in each sector
		sector_sizes	= integer(0)   # number of contiguous tips covered per sector
		previous_name	= NULL
		for(tip in angular_tip_order){
			if(is.null(tip_sectors[tip]) || is.na(tip_sectors[tip])) next
			name = as.character(tip_sectors[tip])
			if(is.null(previous_name) || is.na(previous_name) || (previous_name!=name)){
				# start a new sector
				Nsectors 				= Nsectors + 1
				sector_labels[Nsectors] 	= name
				sector_firsts[Nsectors]	= tip
				sector_lasts[Nsectors]	= tip
				sector_sizes[Nsectors]	= 1
			}else{
				# continue the previous sector
				sector_lasts[Nsectors] = tip
				sector_sizes[Nsectors] = sector_sizes[Nsectors] + 1
			}
			previous_name = name
		}
		# check for fragmented sectors, i.e., whose name appears multiple times in sector_labels[]
		fragmented = (duplicated(sector_labels) | duplicated(sector_labels, fromLast=TRUE))
		if(any(fragmented)){
			fragment_names = unique(sector_labels[fragmented])
			if(fragmented_sectors=="error"){
				stop(sprintf("ERROR: Encountered %d fragmented sectors, for example '%s'.",length(fragment_names),paste0(fragment_names[1:min(3,length(fragment_names))],collapse="', '")))
			}else if(fragmented_sectors=="omit"){
				sector_labels  = sector_labels[!fragmented]
				sector_firsts = sector_firsts[!fragmented]
				sector_lasts  = sector_lasts[!fragmented]
				sector_sizes  = sector_sizes[!fragmented]
			}else if(fragmented_sectors=="split"){
				# fragmented sectors are already split, so no further action needed
			}else{
				stop(sprintf("ERROR: Unknown choice '%s' for fragmented_sectors",fragmented_sectors))
			}
		}
		# filter out short sectors
		long_sectors = which(sector_sizes>=min_sector_size)
		if(length(long_sectors)<length(sector_sizes)){
			sector_labels  = sector_labels[long_sectors]
			sector_firsts = sector_firsts[long_sectors]
			sector_lasts  = sector_lasts[long_sectors]
			sector_sizes  = sector_sizes[long_sectors]
		}
		Nsectors = length(sector_labels)
		# print(data.frame(name=sector_labels,first=sector_firsts,last=sector_lasts,start_phi=180*geometry$clade_phi[sector_firsts]/pi,end_phi=180*geometry$clade_phi[sector_lasts]/pi)) # debug
	}else{
		Nsectors = 0
	}
	
	# instantiate plot file with proper dimensions and margins
	if(is.null(tip_label_cex)) tip_label_cex = min(1, 3*plot_width*max(0,2*pi-opening)*mean(geometry$clade_rho[1:Ntips])/(Ntips*geometry$max_rho))
	legend_margin = (if(is.null(legend_labels) || (length(legend_labels)==0)) 0 else 4+0.4*max(sapply(legend_labels, FUN=nchar))) # margin that should be reserved around the figure to fit the tip labels
	margins_chars = c(1+2*show_scale_bar, # bottom
					1, # left
					2+(1+length(strsplit(plot_title, "\n", fixed=TRUE)[[1]])), # top
					1+legend_margin) # right
	if(tip_label_cex>0) margins_chars = margins_chars + tip_label_cex*0.8*max(sapply(tree$tip.label, FUN=nchar)) # increase margin to accommodate tip labels
	if((Nsectors>0) && (sector_label_cex>0)) margins_chars = margins_chars + sector_label_cex*0.7*(1+max(sapply(sector_labels, FUN=nchar))) # increase margin to accommodate sector labels
	graphics_params = get_dummy_graphics_params()
	margins_inches  = (margins_chars * graphics_params$csi 
					   + (Nrings*(ring_width+ring_border_width)+ring_border_width)/72
					   + NDS*dotstack_border_width/72
					   + (if(NDS==0) 0 else sum(apply(dotstack_heights,2,max)+2)*(1.1*plot_width/2)*(2*pi-opening)/Ntips)
					   + (if(Nsectors==0) 0 else (sector_gap+sector_width)/72))
	margins_chars	= margins_inches/graphics_params$csi # update margin sizes in char units
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
	
	# show dot stacks
	if(NDS>0){
		if(!is.null(rownames(dotstack_heights))){
			tip2row = match(tree$tip.label,rownames(dotstack_heights))
			if(any(is.na(tip2row))) stop(sprintf("%d out of %d tree tips did not match any row names in dotstack_heights[], for example '%s'.",sum(is.na(tip2row)),Ntips,tree$tip.label[is.na(tip2row)][[1]]))
			dotstack_heights = dotstack_heights[tip2row,,drop=FALSE]
		}
		dot_diameter = 0.9*geometry$max_rho*(2*pi-opening)/Ntips # select dot size to correspond to the space allocated to each tip
		for(ds in seq_len(NDS)){
			lowest_dot_rho = geometry$max_rho + dot_diameter
			color	  = pick_recycled(dotstack_color,ds)
			dot_phis  = rep(geometry$clade_phi[1:Ntips], times=dotstack_heights[,ds])
			dot_rhos  = lowest_dot_rho + unlist(lapply(dotstack_heights[,ds], FUN=function(height){ dot_diameter*(seq_len(height)-1) }))
			highest_dot_rho = (if(length(dot_rhos)==0) lowest_dot_rho else max(dot_rhos))
			graphics::symbols(x 	= dot_rhos*cos(dot_phis),
							y 		= dot_rhos*sin(dot_phis),
							circles	= rep(dot_diameter/2.25,length(dot_rhos)),
							col 	= color,
							bg		= color,
							lwd		= 0.00001, # make circle border width negligible so that dot_diameter fully determines the true circle diameter
							add		= TRUE,
							inches	= FALSE)
			# add border outside of this dot stack
			plot_arc(phi0	= geometry$clade_phi0[geometry$root],
					 phi1	= geometry$clade_phi1[geometry$root],
					 rho	= highest_dot_rho + dot_diameter + lw2coordinates*dotstack_border_width/2,
					 col	= dotstack_border_color,
					 lwd	= dotstack_border_width)
			# update max_rho to the end of this dot stack
			geometry$max_rho = highest_dot_rho + dot_diameter + dotstack_border_width*lw2coordinates
		}
	}

	# show tip labels
	if(tip_label_cex!=0){
		max_label_rho_outer = 0 # keep track of the maximum radius increase due to a tip label
		for(tip in seq_len(Ntips)){
			label_srt = (180/pi)*geometry$clade_phi[tip]
			label_adj = c(0, 0.5)
			if((label_srt>90) && (label_srt<270)){
				label_srt = label_srt-180
				label_adj = c(1,0.5)
			}
			label_rho = (if(align_tip_labels || (Nrings>0) || (NDS>0)) geometry$max_rho else geometry$clade_rho[tip]) + tip_label_cex*graphics::par("cxy")[1]*0.7
			graphics::text(	x 		= label_rho * cos(geometry$clade_phi[tip]),
							y 		= label_rho * sin(geometry$clade_phi[tip]),
							labels 	= tree$tip.label[tip],
							col		= (if(length(tip_label_color)==1) tip_label_color else tip_label_color[[tip]]),
							cex		= tip_label_cex,
							srt 	= label_srt,
							adj 	= label_adj)
			max_label_rho_outer = max(max_label_rho_outer, label_rho + tip_label_cex*0.7*graphics::par("cxy")[1]*nchar(tree$tip.label[tip]))
		}
		geometry$max_rho = max_label_rho_outer # update max_rho to account for tip labels
	}
	
	# highlight sectors
	if(Nsectors>0){
		sector_rho 		 = geometry$max_rho + lw2coordinates*(sector_gap+sector_width*0.5)
		cap_len			 = lw2coordinates*3*sector_width
		label_rho 		 = sector_rho + sector_label_cex*graphics::par("cxy")[1]
		geometry$max_rho = max(label_rho, sector_rho + lw2coordinates*sector_width*0.5 + cap_len)
		cap_len 		 = lw2coordinates*sector_width*3
		for(k in seq_along(sector_labels)){
			# draw sector arc
			sector_phi0 = geometry$clade_phi[sector_firsts[k]]
			sector_phi1 = geometry$clade_phi[sector_lasts[k]]
			plot_arc(phi0	= sector_phi0,
					 phi1	= sector_phi1,
					 rho	= sector_rho,
					 col	= sector_color,
					 lwd	= sector_width,
					 lend	= "cap",
					 cap_len= cap_len)
			# add sector label
			label_phi = 0.5*(sector_phi0+sector_phi1)
			label_srt = (180/pi)*label_phi
			label_adj = c(0, 0.5)
			if((label_srt>90) && (label_srt<270)){
				label_srt = label_srt-180
				label_adj = c(1,0.5)
			}
			graphics::text(	x 		= label_rho * cos(label_phi),
							y 		= label_rho * sin(label_phi),
							labels 	= sector_labels[k],
							col		= sector_label_color,
							cex		= sector_label_cex,
							srt 	= label_srt,
							adj 	= label_adj)
			geometry$max_rho = max(geometry$max_rho, label_rho + sector_label_cex*0.7*graphics::par("cxy")[1]*nchar(sector_labels[k]))
		}
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
	if((!is.null(legend_labels)) && (length(legend_labels)>0)){
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


plot_arc = function(phi0, 	# start angle in radians
					phi1,	# end angle in radians
					rho,	# radius, in plot coordinates
					col="black", # color
					lwd=1,	# line width
					lend="round", # line ending. Either "round", "butt", "square" or "cap" (perpendicular caps). 
					cap_len = NULL){ # length of the caps in plot coordinates in each direction. If NULL, this is chosen heuristically
	phis = seq(phi0, phi1, length.out = max(10, as.integer((phi1-phi0)/0.02)))
	if(lend=="cap"){
		# draw perpendicular caps on each side
		if(is.null(cap_len)) cap_len = 0.01*rho
		graphics::segments(	x0	= (rho-cap_len)*cos(phis[c(1,length(phis))]),
							y0	= (rho-cap_len)*sin(phis[c(1,length(phis))]),
							x1	= (rho+cap_len)*cos(phis[c(1,length(phis))]),
							y1	= (rho+cap_len)*sin(phis[c(1,length(phis))]),
							col = col,
							lwd	= lwd,
							lend= "butt")
		lend = "square"
	}
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



