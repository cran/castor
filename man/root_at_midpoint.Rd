\name{root_at_midpoint}
\alias{root_at_midpoint}
\title{
Root or re-root a tree at the midpoint node.
}
\description{
Given a tree (rooted or unrooted), this function changes the direction of edges (\code{tree$edge}) such that the midpoint node becomes the new root (i.e. has no incoming edges and all other tips and nodes descend from it). The number of tips and the number of nodes remain unchanged. The midpoint node is the node whose maximum distance to any tip is smallest.
}
\usage{
root_at_midpoint( tree, 
                  update_indices  = TRUE,
                  as_edge_counts  = FALSE,
                  is_rooted       = FALSE)
}
\arguments{
\item{tree}{
A tree object of class "phylo". Can be unrooted or rooted (but see option \code{is_rooted}).
}
\item{update_indices}{
Logical, specifying whether to update the node indices such that the new root is the first node in the list (as is common convention). This will modify \code{tree$node.label} (if it exists) and also the node indices listed in \code{tree$edge}.
}
\item{as_edge_counts}{
Logical, specifying whether phylogenetic distances should be measured as cumulative edge counts. This is the same if all edges had length 1.
}
\item{is_rooted}{
Logical, specifying whether the input tree can be assumed to be rooted. If you are not certain that the tree is rooted, set this to \code{FALSE}.
}
}

\details{
The input tree may include an arbitrary number of incoming and outgoing edges per node (but only one edge per tip), and the direction of these edges can be arbitrary. Of course, the undirected graph defined by all edges must still be a valid tree. Only set \code{is_rooted=TRUE} if you are sure that the input tree is rooted.

If \code{update_indices==FALSE}, then node indices remain unchanged.
If \code{update_indices==TRUE} (default), then node indices are modified such that the new root is the first node (i.e. with index Ntips+1 in \code{edge} and with index 1 in \code{node.label}), as is common convention. Setting \code{update_indices=FALSE} reduces the computation required for rerooting. Tip indices always remain unchanged.

The asymptotic time complexity of this function is O(Nedges).
}

\value{
A tree object of class "phylo", with the \code{edge} element modified such that the maximum distance of the root to any tip is minimized. The elements \code{tip.label}, \code{edge.length} and \code{root.edge} (if they exist) are the same as for the input tree. If \code{update_indices==FALSE}, then the element \code{node.label} will also remain the same.
}

\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{root_via_outgroup}},
\code{\link{root_at_node}},
\code{\link{root_in_edge}}
}

\examples{
# generate a random tree
Ntips = 100
tree = generate_random_tree(list(birth_rate_intercept=1),max_tips=Ntips)$tree

# reroot the tree at its midpoint node
tree = root_at_midpoint(tree)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ Rerooting }
