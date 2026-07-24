\name{plot_tree_radial}
\alias{plot_tree_radial}
\title{
Plot a tree in radial layout.
}
\description{
Given a rooted phylogenetic tree, visualize the tree in radial form, i.e., with the root placed in the center and descending clades radiating outwards.
Optionally, categorial tip data may be displayed in the form of colored rings and discrete-quantitative tip data may be displayed in the form of dot-stacks around the tree.
}
\usage{
plot_tree_radial(tree,
                file_path             = NULL,
                make_ultrametric      = FALSE,
                symmetric_arcs        = TRUE,
                opening               = 0,
                rotate                = 0,
                plot_width            = 5,
                plot_title            = NULL,
                show_scale_bar        = FALSE,
                align_tip_labels      = FALSE,
                tip_cex               = 0,
                node_cex              = 0,
                tip_label_cex         = NULL,
                scale_edge_widths     = FALSE,
                base_edge_width       = 1,
                tip_color             = "black",
                node_color            = "black",
                edge_color            = "#00000080",
                node_arc_color        = "#00000080",
                root_edge_color       = "#00000080",
                tip_label_color       = "black",
                ring_colors           = NULL,
                ring_width            = 5,
                ring_border_width     = 1,
                ring_border_color     = "#FFFFFF",
                dotstack_heights      = NULL,
                dotstack_color        = "black",
                dotstack_border_width = 1,
                dotstack_border_color = "black",
                tip_sectors           = NULL,
                sector_color          = "black",
                sector_width          = 2,
                sector_gap            = 15,
                min_sector_size       = 2,
                sector_label_cex      = 1,
                sector_label_color    = "black",
                fragmented_sectors    = "error",
                legend_colors         = NULL,
                legend_labels         = NULL)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{file_path}{
Optional character, file path to save the plot to. Any existing file is silently replaced. The image type is determined based on the file name extension (pdf, jpg/jpeg, png, tif/tiff, bmp).
}
\item{make_ultrametric}{
Logical, whether to make the tree ultrametric first, so that all tips appear at the same distance from the center.
The transformation is done based on relative evolutionary divergences (Parks et al. 2018).
}
\item{symmetric_arcs}{
Logical, whether nodes should appear in the middle of their arc. 
If \code{FALSE}, nodes are radially placed in the middle of their descending tips, which could result in nodes not being in the center of their arc.
}
\item{opening}{
Optional numeric, opening angle of the tree in degrees, i.e., separating the two flanks of the tree.
If 0, the tips will uniformly cover the whole circle. Must be between 0 and 360.
}
\item{rotate}{
Optional numeric, angle (degrees) by which to rotate the whole plot, counter-clockwise. By default the tree is oriented such that the outer-most tip is placed on the right-most point of the circle, i.e., at "3'o clock".
}
\item{plot_width}{
Plot width & height in inches. Only refers to the tree portion, i.e., excluding rings, dot-stacks, sector arcs/labels and margins, which are allocated space automatically. Only relevant if \code{file_path} is specified.
}
\item{plot_title}{
Optional character, plot title.
}
\item{show_scale_bar}{
Logical, whether to show a scale bar in the corner of the plot.
}
\item{align_tip_labels}{
Logical, whether to align all tip labels on a single circle, rather than placing each tip label immediately adjacent to each tip. Only relevant for non-ultrametric trees, only if \code{show_tip_labels=TRUE} and only if no rings are shown (see \code{ring_colors}). If colored rings are shown, then tip labels are always aligned.
}
\item{tip_cex}{
Numeric, circle rescaling for tips. If >0, tips are shown as circles.
}
\item{node_cex}{
Numeric, circle rescaling for nodes. If >0, nodes are shown as circles.
}
\item{tip_label_cex}{
Numeric, tip label rescaling. If 0, tip labels are not shown. If \code{NULL}, this is adjusted automatically to try to make tip labels fit.
}
\item{scale_edge_widths}{
Logical, whether to scale the width (thickness) of each edge and node arc proportionally to the square root of the number of tips descending from it.
If \code{FALSE}, all edges will have the same width.
}
\item{base_edge_width}{
Numeric, base (smallest) width of edges, e.g. at the tips (and everywhere else if \code{scale_edge_widths==FALSE}).
}
\item{tip_color}{
Either a single character or a character vector of length \code{Ntips} or an RGB tuple or an RGBA tuple, specifying tip colors & alpha. This vector must be synchronized with the tree tips.
}
\item{node_color}{
Either a single character or a character vector of length \code{Nnodes} or an RGB tuple or an RGBA tuple, specifying node colors & alpha. This vector must be synchronized with the tree nodes.
}
\item{edge_color}{
Either a single character or a character vector of length \code{Nedges} or an RGB tuple or an RGBA tuple, specifying edge colors & alpha. This vector must be synchronized with the tree edges, i.e., as specified in \code{tree$edge} or \code{tree$edge.length}.
}
\item{node_arc_color}{
Either a single character or a character vector of length \code{Nodes} or an RGB tuple or an RGBA tuple, specifying node arc colors & alpha. This vector must be synchronized with the tree nodes.
}
\item{root_edge_color}{
Color of the root edge, if applicable.
}
\item{tip_label_color}{
Either a single character or a character vector of length \code{Ntips} or an RGB tuple or an RGBA tuple, specifying tip label colors & alpha. This vector must be synchronized with the tree tips.
}
\item{ring_colors}{
Optional specification of categorial data to be visualized as multi-colored rings around the tree. Must be either a data.frame or a matrix of size Ntips x NC, specifying one or more multi-colored rings to show around the tree.
The entries of \code{ring_colors} must be valid HEX color characters, with \code{ring_colors[t,r]} being the color for tip t on ring r. An entry of \code{NULL} or \code{NA} means no color segment is drawn for a tip, effectively keeping the background color.
If \code{ring_colors} has row names, these will be matched to the tip labels, otherwise rows are assumed to already be synchronized with the tree's tips.
}
\item{ring_width}{
Numeric, line width to use for colored rings around the tree. Only relevant if \code{ring_colors} is provided.
}
\item{ring_border_width}{
Numeric, line width to use for the ring borders. A value of 0 means zero space between successive rings.
}
\item{ring_border_color}{
Color to use for the ring borders, i.e., separating the rings from each other and from the tree.
}
\item{dotstack_heights}{
Optional specification of discrete quantitative data to be visualized as stacked dots around the tree. Must be either a data.frame or a matrix of size Ntips x ND, specifying one or more dot stacks to show around the tree.
The entries of \code{dotstack_heights} must be non-negative integers, with \code{dotstack_heights[t,d]} being the dot-stack height (number of dots) for tip t on dot stack d.
If \code{dotstack_heights} has row names, these will be matched to the tip labels, otherwise rows are assumed to already be synchronized with the tree's tips.
An example use case would be to show the number of copies found for a specific gene in each genome on a genome phylogeny.
}
\item{dotstack_color}{
Either a single color or a vector colors, specifying the dot colors to use for the dot-stacks. For example, colors may be specified as HEX codes or color names. Colors are recycled as needed.
}
\item{dotstack_border_width}{
Numeric, line thickness for the dot-stack borders. If \code{NULL}, this is set equal to \code{ring_border_width}.
}
\item{dotstack_border_color}{
Color to use for the dot-stack borders, i.e., separating the dot-stacks from each other. If \code{NULL}, this is set to \code{ring_border_color}.
}
\item{tip_sectors}{
Optional vector of length Ntips, specifying a sector ID for each tip. IDs may be characters, integers, or some other reasonable scalar. Tips will be grouped by sector, and each such group will be visually highlighted using an outer arc and label. This may be used e.g. to highlight major phyla in a tree. \code{tip_sectors} may also contain \code{NA} or \code{NULL}, representing sector-less tips not to be highlighted. It is expected that all or most sectors are monophyletic (but see option \code{fragmented_sectors}). If \code{tip_sectors} has names, these will be matched to the tip labels, otherwise elements in \code{tip_sectors} are assumed to already be synchronized with the tree's tips.
}
\item{sector_color}{
Color to use for sector arcs, if applicable.
}
\item{sector_width}{
Numeric, line thickness for the sector arcs.
}
\item{sector_gap}{
Numeric, gap between the sector arcs and the inner-next element (typically tip labels), measured in the same units as line widths.
}
\item{min_sector_size}{
Integer, minimum number of contiguous tips per displayed sector (after splitting, if \code{fragmented_sectors=="split"}).
}
\item{sector_label_cex}{
Numeric, sector label rescaling if applicable. If 0, sector labels are not shown.
}
\item{sector_label_color}{
Color for the sector labels, if applicable.
}
\item{fragmented_sectors}{
Character, specifying how to handle fragmented (i.e., non-monophyletic) tip sectors. Options include "\code{error}" (throw an error), "\code{omit}" (omit the affected sectors) and "\code{split}" (split the affected sectors into multiple parts where needed).
}
\item{legend_colors}{
Vector specifying colors to be defined in a legend. These can be for example colors in \code{ring_colors[]}, or colors used for the tips & nodes & edges.
}
\item{legend_labels}{
Optional character vector specifying labels to be given in a legend, synchronized with \code{legend_colors[]}.
}
}


\details{
The input tree may include multi-furcations (i.e. nodes with more than 2 children) as well as mono-furcations (i.e. nodes with only one child). 
If \code{tree$edge.length} is \code{NULL}, then all edges in the input tree are assumed to have length 1.
}


\value{
This function returns no value. It either displays a plot (if \code{file_path==NULL}) or saves a plot to a file.
}

\author{Stilianos Louca}

\references{
D. H. Parks, M. Chuvochina et al. (2018). A proposal for a standardized bacterial taxonomy based on genome phylogeny. bioRxiv 256800. DOI:10.1101/256800
}

\seealso{
\code{\link{write_tree}}
}

\examples{
# generate a random ultrametric tree
tree = generate_random_tree(list(birth_rate_intercept=1), max_tips=50, tip_basename="tip.")$tree

# plot in radial layout
castor::plot_tree_radial(tree           = tree,
                        show_scale_bar  = TRUE,
                        tip_cex         = 0.5,
                        node_cex        = 0)
}
\keyword{plotting}
