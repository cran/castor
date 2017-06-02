\name{fit_bm_model}
\alias{fit_bm_model}
\title{
Fit a Brownian motion model for multivariate trait evolution.
}
\description{
Given a rooted phylogenetic tree and states of one or more continuous (numeric) traits on the tree's tips, fit a multivariate Brownian motion model of correlated co-evolution of these traits. This estimates a single diffusivity matrix, which describes the variance-covariance structure of each trait's random walk. The model assumes a fixed diffusivity matrix on the entire tree.
}
\usage{
fit_bm_model(tree, tip_states, check_input=TRUE)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{tip_states}{
A numeric vector of size Ntips, or a 2D numeric matrix of size Ntips x Ntraits, specifying the numeric state of each trait at each tip in the tree.
}
\item{check_input}{
Logical, specifying whether to perform some basic checks on the validity of the input data. If you are certain that your input data are valid, you can set this to \code{FALSE} to reduce computation.
}
}


\details{
The BM model is described by the stochastic differential equation
\deqn{
dX = \sigma dB
}
where \eqn{B} is a multivariate standard Brownian motion, or equivalently, as a diffusion equation for the probability density
\deqn{
\frac{\partial \rho}{\partial t} = D\frac{\partial^2\rho}{\partial x^2}.
}
The matrix \eqn{D} is referred to as the diffusion coefficient matrix (or diffusivity matrix), and \eqn{D=\sigma\cdot\sigma^T/2}.

The function uses phylogenetic independent contrasts (Felsenstein, 1985) to retrieve independent increments of the multivariate random walk. The diffusivity matrix \eqn{D} is then fitted using maximum-likelihood on the intrinsic geometry of positive-definite matrices.

If \code{tree$edge.length} is missing, each edge in the tree is assumed to have length 1. The tree may include multifurcations (i.e. nodes with more than 2 children) as well as monofurcations (i.e. nodes with only one child). Note that multifurcations are internally expanded to bifurcations, prior to model fitting. The asymptotic time complexity of this function is O(Nedges*Nsimulations*Ntraits).
}


\value{
A list with the following elements:
\item{diffusivity}{
Either a single non-negative number (if \code{tip_states} was a vector) or a 2D quadratic non-negative-definite matrix (if \code{tip_states} was a 2D matrix). The fitted diffusivity matrix of the multivariate Brownian motion model.
}
\item{loglikelihood}{
The log-likelihood of the fitted model, given the provided tip states data.
}
}

\author{Stilianos Louca}

\references{
J. Felsenstein (1985). Phylogenies and the Comparative Method. The American Naturalist. 125:1-15.
}

\seealso{
\code{\link{simulate_bm_model}},
\code{\link{get_independent_contrasts}}
}

\examples{
# generate a random tree
tree = generate_random_tree(max_tips=10000, birth_rate_intercept=1);


# Example 1: Univariate case
# - - - - - - - - - - - - - -
# simulate scalar continuous trait on the tree
D = 1
tip_states = simulate_bm_model(tree, diffusivity=D)$tip_states

# estimate original diffusivity from the generated data
fit = fit_bm_model(tree, tip_states)
cat(sprintf("True D=\%g, fitted D=\%g\n",D,fit$diffusivity))


# Example 2: Multivariate case
# - - - - - - - - - - - - - - -
# simulate vector-valued continuous trait on the tree
D = get_random_diffusivity_matrix(Ntraits=5)
tip_states = simulate_bm_model(tree, diffusivity=D)$tip_states

# estimate original diffusivity matrix from the generated data
fit = fit_bm_model(tree, tip_states)

# compare true and fitted diffusivity matrices
cat("True D:\n"); print(D)
cat("Fitted D:\n"); print(fit$diffusivity)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{BM model}
