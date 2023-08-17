\name{root_via_rtt}
\alias{root_via_rtt}
\title{
Root a tree via root-to-tip regression.
}
\description{
Root a non-dated tree based on tip sampling times, by optimizing the goodness of fit of a linear root-to-tip (RTT) regression (regression of tip times vs phylogenetic distances from root). The precise objective optimized can be chosen by the user, typical choices being \code{R2} or \code{SSR} (sum of squared residuals). This method is only suitable for clades that are "measurably evolving". The input tree's edge lengths should be measured in substitutions per site.
}
\usage{
root_via_rtt(tree,
             tip_times,
             objective           = "R2",
             force_positive_rate = FALSE,
             Nthreads            = 1,
             optim_algorithm     = "nlminb",
             relative_error      = 1e-9)
}
\arguments{
\item{tree}{
A tree object of class "phylo". Can be unrooted or rooted (the root placement does not matter). Edge lengths should be measured in expected substitutions per site.
}
\item{tip_times}{
Numeric vector of length Ntips, listing the sampling times of all tips. Time is measured in forward direction, i.e., younger tips have a greater time value. Note that if you originally have tip sampling dates, you will first need to convert these to numeric values (for example decimal years or number of days since the start of the experiment).
}
\item{objective}{
Character, specifying the goodness-of-fit measure to consider for the root-to-tip regression. Must be one of \code{correlation}, \code{R2} (fraction of explained variance) or \code{SSR} (sum of squared residuals).
}
\item{force_positive_rate}{
Logical, whether to force the mutation rate implied by the root placement to be positive (>=0).
}
\item{Nthreads}{
Integer, number of parallel threads to use where applicable.
}
\item{optim_algorithm}{
Character, the optimization algorithm to use. Must be either \code{nlminb} or \code{optimize}.
}
\item{relative_error}{
Positive numeric, specifying the acceptable relative error when optimizing the goodness of fit. The precise interpretation depends on the optimization algorithm used. Smaller values may increase accuracy but also computing time.
}
}

%\details{
%}

\value{
A named list with the following elements (more may be added in the future):
\item{tree}{
The rooted tree. A tree object of class "phylo", with the same tips as the original tree (not necessarily in the original order).
}
}

\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{root_via_outgroup}},
\code{\link{root_at_node}},
\code{\link{root_in_edge}}
\code{\link{root_at_midpoint}}
}

\examples{
# generate a random tree
Ntips = 10
tree = generate_random_tree(list(birth_rate_intercept=1),max_tips=Ntips)$tree

# construct a vector with hypothetical tip sampling times
tip_times = c(2010.5, 2010.7, 2011.3, 2008.7,
              2009.1, 2013.9, 2013.8, 2011.4,
              2011.7, 2005.2)

# reroot the tree via root-to-tip regression
tree = root_via_rtt(tree, tip_times=tip_times)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{rooting}
