\name{fit_and_compare_bm_models}
\alias{fit_and_compare_bm_models}
\title{
Fit and compare Brownian Motion models for multivariate trait evolution between two data sets.
}
\description{
Given two rooted phylogenetic trees and states of one or more continuous (numeric) traits on the trees' tips, fit a multivariate Brownian motion model of correlated evolution to each data set and compare the fitted models. This function estimates the diffusivity matrix for each data set (i.e., each tree/tip-states set) via maximum-likelihood and assesses whether the log-difference between the two fitted diffusivity matrixes is statistically significant, under the null hypothesis that the two data sets exhibit the same diffusivity. Optionally, multiple trees can be used as input for each data set, under the assumption that the trait evolved on each tree according to the same BM model. For more details on how BM is fitted to each data set see the function \code{\link{fit_bm_model}}.
}
\usage{
fit_and_compare_bm_models(  trees1, 
                            tip_states1,
                            trees2,
                            tip_states2,
                            Nbootstraps     = 0,
                            Nsignificance   = 0,
                            check_input     = TRUE,
                            verbose         = FALSE,
                            verbose_prefix  = "")
}
\arguments{
\item{trees1}{
Either a single rooted tree or a list of rooted trees, of class "phylo", corresponding to the first data set on which a BM model is to be fitted. Edge lengths are assumed to represent time intervals or a similarly interpretable phylogenetic distance.
}
\item{tip_states1}{
Numeric state of each trait at each tip in each tree in the first data set. If \code{trees1} is a single tree, then \code{tip_states1} must either be a numeric vector of size Ntips or a 2D numeric matrix of size Ntips x Ntraits, listing the trait states for each tip in the tree.
If \code{trees1} is a list of Ntrees trees, then \code{tip_states1} must be a list of length Ntrees, each element of which lists the trait states for the corresponding tree (as a vector or 2D matrix, similarly to the single-tree case).
}
\item{trees2}{
Either a single rooted tree or a list of rooted trees, of class "phylo", corresponding to the second data set on which a BM model is to be fitted. Edge lengths are assumed to represent time intervals or a similarly interpretable phylogenetic distance.
}
\item{tip_states2}{
Numeric state of each trait at each tip in each tree in the second data set, similarly to \code{tip_states1}.
}
\item{Nbootstraps}{
Integer, specifying the number of parametric bootstraps to perform for calculating the confidence intervals of BM diffusivities fitted to each data set. If <=0, no bootstrapping is performed.
}
\item{Nsignificance}{
Integer, specifying the number of simulations to perform for assessing the statistical significance of the log-transformed difference between the diffusivities fitted to the two data sets, i.e. of \eqn{|\log(D_1)-\log(D_2)|}. Set to 0 to not calculate the statistical significance. See below for additional details.
}
\item{check_input}{
Logical, specifying whether to perform some basic checks on the validity of the input data. If you are certain that your input data are valid, you can set this to \code{FALSE} to reduce computation.
}
\item{verbose}{
Logical, specifying whether to print progress report messages to the screen.
}
\item{verbose_prefix}{
Character, specifying a prefix to include in front of progress report messages on each line. Only relevant if \code{verbose==TRUE}.
}
}


\details{
For details on the Brownian Motion model see \code{\link{fit_bm_model}} and \code{\link{simulate_bm_model}}. This function separately fits a single-variate or multi-variate BM model with constant diffusivity (diffusivity matrix, in the multivariate case) to each data set; internally, this function applies \code{fit_bm_model} to each data set. 

If \code{Nsignificance>0}, the statistical significance of the log-transformed difference of the two fitted diffusivity matrixes, \eqn{|\log(D_1)-\log(D_2)|}, is assessed, under the null hypothesis that both data sets were generated by the same common BM model. The diffusivity of this common BM model is estimated by fitting to both datasets at once, i.e. after merging the two datasets into a single dataset of trees and tip states (see return variable \code{fit_common} below). For each of the \code{Nsignificance} random simulations of the common BM model on the two tree sets, the diffusivities are again separately fitted on the two simulated sets and the resulting log-difference is compared to the one of the original data sets. The returned \code{significance} is the probability that the diffusivities would have a log-difference larger than the observed one, if the two data sets had been generated under the common BM model.

If \code{tree$edge.length} is missing, each edge in the tree is assumed to have length 1. The tree may include multifurcations (i.e. nodes with more than 2 children) as well as monofurcations (i.e. nodes with only one child). Note that multifurcations are internally expanded to bifurcations, prior to model fitting.
}


\value{
A list with the following elements:
\item{success}{
Logical, indicating whether the fitting was successful for both data sets. If \code{FALSE}, then an additional return variable, \code{error}, will contain a description of the error; in that case all other return variables may be undefined.
}
\item{fit1}{
A named list containing the fitting results for the first data set, in the same format as returned by \code{\link{fit_bm_model}}. In particular, the diffusivity fitted to the first data set will be stored in \code{fit1$diffusivity}.
}
\item{fit2}{
A named list containing the fitting results for the second data set, in the same format as returned by \code{\link{fit_bm_model}}. In particular, the diffusivity fitted to the second data set will be stored in \code{fit2$diffusivity}.
}
\item{log_difference}{
The absolute difference between the log-transformed diffusivities, i.e. \eqn{|\log(D_1)-\log(D_2)|}. In the multivariate case, this will be a matrix of size Ntraits x Ntraits.
}
\item{significance}{
Numeric, statistical significance of the observed log-difference under the null hypothesis that the two data sets were generated by a common BM model. Only returned if \code{Nsignificance>0}.
}
\item{fit_common}{
A named list containing the fitting results for the two data sets combined, in the same format as returned by \code{\link{fit_bm_model}}. The common diffusivity, \code{fit_common$diffusivity} is used for the random simulations when assessing the statistical significance of the log-difference of the separately fitted diffusivities. Only returned if \code{Nsignificance>0}.
}
}

\author{Stilianos Louca}

\references{
J. Felsenstein (1985). Phylogenies and the Comparative Method. The American Naturalist. 125:1-15.
}

\seealso{
\code{\link{simulate_bm_model}},
\code{\link{fit_bm_model}},
\code{\link{get_independent_contrasts}}
}

\examples{
# simulate distinct BM models on two random trees
D1 = 1
D2 = 2
tree1 = generate_random_tree(list(birth_rate_factor=1),max_tips=100)$tree
tree2 = generate_random_tree(list(birth_rate_factor=1),max_tips=100)$tree
tip_states1 = simulate_bm_model(tree1, diffusivity = D1)$tip_states
tip_states2 = simulate_bm_model(tree2, diffusivity = D2)$tip_states

# fit and compare BM models between the two data sets
fit = fit_and_compare_bm_models(trees1          = tree1, 
                                tip_states1     = tip_states1, 
                                trees2          = tree2,
                                tip_states2     = tip_states2,
                                Nbootstraps     = 100,
                                Nsignificance   = 100)

# print summary of results
cat(sprintf("Fitted D1 = \%g, D2 = \%g, significance of log-diff. = \%g\n",
            fit$fit1$diffusivity, fit$fit2$diffusivity, fit$significance))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{BM model}
