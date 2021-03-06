\name{simulate_bm_model}
\alias{simulate_bm_model}
\title{
Simulate a Brownian motion model for multivariate trait co-evolution.
}
\description{
Given a rooted phylogenetic tree and a Brownian motion (BM) model for the co-evolution of one or more continuous (numeric) unbounded traits, simulate random outcomes of the model on all nodes and/or tips of the tree. The function traverses nodes from root to tips and randomly assigns a multivariate state to each node or tip based on its parent's previously assigned state and the specified model parameters. The generated states have joint distributions consistent with the multivariate BM model. Optionally, multiple independent simulations can be performed using the same model.
}
\usage{
simulate_bm_model(tree, diffusivity=NULL, sigma=NULL,
                  include_tips=TRUE, include_nodes=TRUE, 
                  root_states=NULL, Nsimulations=1, drop_dims=TRUE)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{diffusivity}{
Either NULL, or a single number, or a 2D quadratic positive definite symmetric matrix of size Ntraits x Ntraits. Diffusivity matrix ("\eqn{D}") of the multivariate Brownian motion model (in units trait^2/edge_length). The convention is that if the root's state is fixed, then the covariance matrix of a node's state at distance \eqn{L} from the root will be \eqn{2LD} (see mathematical details below).
}
\item{sigma}{
Either NULL, or a single number, or a 2D matrix of size Ntraits x Ndegrees, where Ndegrees refers to the degrees of freedom of the model. Noise-amplitude coefficients of the multivariate Brownian motion model (in units trait/sqrt(edge_length)). This can be used as an alternative way to specify the Brownian motion model instead of through the diffusivity \eqn{D}. Note that \eqn{sigma\cdot\sigma^T=2D} (see mathematical details below).
}
\item{include_tips}{
Include random states for the tips. If \code{FALSE}, no states will be returned for tips.
}
\item{include_nodes}{
Include random states for the nodes. If \code{FALSE}, no states will be returned for nodes.
}
\item{root_states}{
Numeric matrix of size NR x Ntraits (where NR can be arbitrary), specifying the state of the root for each simulation. If NR is smaller than \code{Nsimulations}, values in \code{root_states} are recycled in rotation. If \code{root_states} is \code{NULL} or empty, then the root state is set to 0 for all traits in all simulations.
}
\item{Nsimulations}{
Number of random independent simulations to perform. For each node and/or tip, there will be \code{Nsimulations} random states generated.
}
\item{drop_dims}{
Logical, specifying whether singleton dimensions should be dropped from \code{tip_states} and \code{node_states}, if \code{Nsimulations==1} and/or Ntraits==1. If \code{drop_dims==FALSE}, then \code{tip_states} and \code{tip_nodes} will always be 3D matrices.
}
}


\details{
The BM model for Ntraits co-evolving traits is defined by the stochastic differential equation
\deqn{
dX = \sigma \cdot dW
}
where \eqn{W} is a multidimensional Wiener process with Ndegrees independent components and \eqn{\sigma} is a matrix of size Ntraits x Ndegrees. Alternatively, the same model can be defined as a Fokker-Planck equation for the probability density \eqn{\rho}:
\deqn{
\frac{\partial \rho}{\partial t} = \sum_{i,j}D_{ij}\frac{\partial^2\rho}{\partial x_i\partial x_j}.
}
The matrix \eqn{D} is referred to as the diffusivity matrix (or diffusion tensor), and \eqn{2D=\sigma\cdot\sigma^T}. Either \code{diffusivity} (\eqn{D}) or \code{sigma} (\eqn{\sigma}) may be used to specify the BM model, but not both.

If \code{tree$edge.length} is missing, each edge in the tree is assumed to have length 1. The tree may include multi-furcations (i.e. nodes with more than 2 children) as well as mono-furcations (i.e. nodes with only one child). The asymptotic time complexity of this function is O(Nedges*Nsimulations*Ntraits).
}


\value{
A list with the following elements:
\item{tip_states}{
Either \code{NULL} (if \code{include_tips==FALSE}), or a 3D numeric matrix of size Nsimulations x Ntips x Ntraits. The [r,c,i]-th entry of this matrix will be the state of trait i at tip c generated by the r-th simulation. If \code{drop_dims==TRUE} and \code{Nsimulations==1} and Ntraits==1, then \code{tip_states} will be a vector.
}
\item{node_states}{
Either \code{NULL} (if \code{include_nodes==FALSE}), or a 3D numeric matrix of size Nsimulations x Nnodes x Ntraits. The [r,c,i]-th entry of this matrix will be the state of trait i at node c generated by the r-th simulation. If \code{drop_dims==TRUE} and \code{Nsimulations==1} and Ntraits==1, then \code{node_states} will be a vector.
}
}

\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{simulate_ou_model}}, 
\code{\link{simulate_rou_model}},
\code{\link{simulate_mk_model}},
\code{\link{fit_bm_model}} 
}

\examples{
# generate a random tree
tree = generate_random_tree(list(birth_rate_intercept=1),max_tips=10000)$tree

# Example 1: Scalar case
# - - - - - - - - - - - - - - -
# simulate scalar continuous trait evolution on the tree
tip_states = simulate_bm_model(tree, diffusivity=1)$tip_states

# plot histogram of simulated tip states
hist(tip_states, breaks=20, xlab="state", main="Trait probability distribution", prob=TRUE)

# Example 2: Multivariate case
# - - - - - - - - - - - - - - -
# simulate co-evolution of 2 traits with 3 degrees of freedom
Ntraits  = 2
Ndegrees = 3
sigma    = matrix(stats::rnorm(n=Ntraits*Ndegrees, mean=0, sd=1), ncol=Ndegrees)
tip_states = simulate_bm_model(tree, sigma=sigma, drop_dims=TRUE)$tip_states

# generate scatterplot of traits across tips
plot(tip_states[,1],tip_states[,2],xlab="trait 1",ylab="trait 2",cex=0.5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{BM model}
\keyword{random}
\keyword{simulation}
