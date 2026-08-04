\name{correlate_categorical_traits}
\alias{correlate_categorical_traits}
\title{
Correlations between categorical traits on trees.
}
\description{
Given one or more rooted phylogenetic trees and states for two categorical traits on the tree tips, examine the correlation between the two traits while accounting for phylogenetic relationships between tips. Regular correlation tests treating each tip as an independent sample fail to account for statistical dependencies between closely related tips, and thus misestimate statistical significances; this function avoids this fallacy by using a stochastic null model for trait evolution on phylogenies.

This function can handle missing states encoded as \code{NA}, erroneous state measurements provided that error rates are a priori known, as well as sampling biases provided that the biases are a priori known. The function can handle binary traits (e.g., presence/absence of wings) as well as multi-state traits (e.g., representing multiple levels of wing complexity). Multiple test statistics are examined, including Pearson correlation, Spearman rank correlation, normalized mutual information (Strehl & Ghosh 2002) and likelihood ratio test (Pagel 1994). The function supportsconventional homogeneous models sensu Pagel (1994) as well as heterogeneous hidden Markov models sensu Boyko & Beaulieu (2023).
}
\usage{
correlate_categorical_traits(trees,
                             tip_states1,
                             tip_states2,
                             value_probs1      = NULL,
                             value_probs2      = NULL,
                             inclusion_probs1  = NULL,
                             inclusion_probs2  = NULL,
                             null_rate_model1  = "ARD",
                             null_rate_model2  = "ARD",
                             root_prior        = "stationary",
                             dependent_trait   = 0,
                             Ntrials           = 10,
                             Ntrials_sim       = 1,
                             Nthreads          = -1,
                             Nsimulations      = 1000,
                             max_model_runtime = NULL,
                             include_lrt       = FALSE,
                             heterogeneous     = FALSE,
                             verbose           = FALSE,
                             verbose_prefix    = "")
}
\arguments{
\item{trees}{
A single rooted tree of class "phylo", or a list of NT such trees. Typically, edge lengths will be in units of time (i.e., the trees are dated), but may also be in other units such as nucleotide substitutions per site. The edge length units will influence the precise meaning of the null models used.
}
\item{tip_states1}{
A list of NT integer vectors, specifying the state of each tip in each tree, for trait 1. Hence, \code{tip_states1[[i]][j]} is the state of trait 1 at tip j of tree i. If only a single tree is provided (NT=1), then \code{tip_states1} may also be given as a single integer vector. Values must be integers between 1 and Nstates1, where Nstates1 is the total number of possible states for trait 1. Some states may be unknown, i.e., represented by \code{NA}. 

Note that tip states must be encoded as integers >=1; for example, binary traits must be encoded using the values 1 and 2, not 0 and 1. This convention allows a direct correspondence between tip states and the indexing of other auxiliary inputs such as \code{value_probs1} and \code{inclusion_probs1} (see below). Further, this convention ensures that for ordinal traits the states are enumerated in the right order whenever applicable; for example, if a discretized quantitative trait can have states "low", "medium", "high", these should be encoded as 1, 2 and 3. 

For non-ordinal traits, i.e., lacking any intrinsic ordering (e.g., plumage color), states may be encoded in any arbitrary order, as long as this is done consistently. Note that for such traits some computed correlation metrics such as Pearson and Spearman may not be meaningful.
}
\item{tip_states2}{
Similar to \code{tip_states1}, but for trait 2.
}
\item{value_probs1}{
Optional list of NT numeric 3-dimensional arrays, each specifying conditional probabilities for measured states of trait 1 on each tree. This argument allows incorporating custom error probabilities, separately for each tip on each tree and depending on a tip's true state. Each \code{value_probs1[[i]]} is a 3D array of size \code{Ntips[i] x Nstates1 x Nstates1}, where \code{Ntips[i]} is the number of tips in the \code{i}'th tree, and where \code{value_probs1[[i]][j,x,y]} is the conditional probability of measuring state \code{y} at the \code{j}'th tip of the \code{i}'th tree, given that its true state is \code{x} and given that a state measurement took place. For any given \code{i,j,x}, the values in \code{value_probs1[[i]][j,x,]} must sum to 1. For any given tip, the matrix \code{value_probs[[i]][j,,]} is also known as "emission matrix" in the context of hidden Markov models. These emission matrixes influence the model fitting and simulations used to assess statistical significances. If only a single tree is provided (i.e., NT=1), then \code{value_probs1} may also be provided as a single 3D array.

For example, suppose trait 1 represents the absence (state 1) or presence (state 2) of a gene family in a bacterial strain, and a given tip \code{j} in tree \code{i} represents an incompletely sequenced genome with completeness \eqn{0\leq C\leq 1}, in which the gene is detected using standard bioinformatics tools. A reasonable assumption could then be that the gene is never erroneously detected if it is absent in the strain, and is detected with probability \eqn{C} if present in the strain. To encode these probabilities, set \code{value_probs1[[i]][j,1,1]=1}, \code{value_probs1[[i]][j,1,2]=0}, \code{value_probs1[[i]][j,2,1]=1-C} and \code{value_probs1[[i]][j,2,2]=C}.
}
\item{value_probs2}{
Similar to \code{value_probs1}, but for trait 2.
}
\item{inclusion_probs1}{
Optional list of NT numeric vectors, each specifying inclusion probabilities for trait 1 on a tree. This argument allows specifying the probability of any state being measured on a tip (as opposed to being \code{NA}), separately for each tree and as a function of a tip's true state, thus accounting for potential sampling biases. Hence, each \code{inclusion_probs1[[i]]} is a numeric vector of size \code{Nstates1}, with \code{inclusion_probs1[[i]][x]} being for any given tip on tree \code{i} the conditional probability of measuring any state (regardless of right or wrong), given that the tip's true state was \code{x}. 

While each \code{inclusion_probs1[[i]][x]} should be a number between 0 and 1, the entries in \code{inclusion_probs1[[i]]} don't need to sum to 1. In fact, for any given tree \code{i} only the relative values in \code{inclusion_probs1[[i]]} matter, i.e., the ratios between \code{inclusion_probs1[[i]][1]}, \code{inclusion_probs1[[i]][2]} and so on, but not their absolute values. During simulations of the null, each \code{inclusion_probs1[[i]]} is rescaled to be consistent with the fraction of non-NA tip states in that tree under stationarity.

For convenience, \code{inclusion_probs1} may also be provided as a single numeric vector of size \code{Nstates1}, in which case it applies equally to all trees. If \code{inclusion_probs1} is \code{NULL}, then a state-independent (i.e., unbiased) inclusion probability is estimated from the fraction of non-\code{NA} tip states, separately for each tree.

As an example, suppose binary trait 1 represents the ability of microbial species to perform methanogenesis (1:no, 2:yes), and suppose non-methanogens are known to be twice as likely than methanogens to be studied in sufficient detail to assess the trait. Then \code{inclusion_probs1[[i]][1]} should be set to twice as high as \code{inclusion_probs1[[i]][2]}.
}
\item{inclusion_probs2}{
Similar to \code{inclusion_probs1}, but for trait 2.
}
\item{null_rate_model1}{
Character, specifying the null Markov transition rate model to use for trait 1. Choices are "ARD" (all rates different), "ER" (all rates equal), "SYM" (equal rates in each direction), "SUEDE" and "SRD". For details see the documentation of \code{\link{fit_mk}}. Also used to construct the full rate model for the likelihood ratio test (see option \code{include_lrt}).
}
\item{null_rate_model2}{
Similar to \code{null_rate_model1}, but for trait 2.
}
\item{root_prior}{
Character, prior probability distribution of the root's states, used to calculate each model's likelihood during fitting. For options and explanations see the documentation of \code{\link{fit_mk}}. The default is recommended, but other options are available for comparison with other software that might use a different prior. This argument does NOT affect the model simulations, where the root state is always drawn from the stationary distribution.
}
\item{dependent_trait}{
Integer, specifying which trait to consider as the dependent (influenced) trait during full model fitting. Options are 0 (both traits influence each other's transition rates), 1 (trait 2 influences the rates of trait 1, but not vise versa) and 2 (trait 1 influences the rates of trait 2, but not vise versa). Only relevant if \code{include_lrt=TRUE}. This is analogous to the \code{dev.var} argument in \code{phytools::fitPagel}.
}
\item{Ntrials}{
Strictly positive integer, number of independent fitting trials to attempt when fitting models to the input data. A larger number increases the chances of an optimal model fit (i.e., avoiding non-global maxima in the likelihood function), at the expense of increased computation. This argument is passed directly to \code{\link{fit_mk}}.
}
\item{Ntrials_sim}{
Similar to \code{Ntrials}, but specifically for fitting models to simulated data. Only relevant for likelihood ratio tests (see option \code{include_lrt}), which involve refitting models to simulated data. This argument is passed directly to \code{\link{fit_mk}}.
}
\item{Nthreads}{
Integer, number of parallel threads to use whenever applicable. If negative, this defines the number of desired threads minus the number of available cores; for example -1 means use all cores except one. Only relevant for Linux/UNIX/MacOS; for Windows, computation is always performed in a single thread.
}
\item{Nsimulations}{
Strictly positive integer, number of simulations to perform of the null model, for estimating statistical significances (aka. P-values). A larger \code{Nsimulations} improves the accuracy of the estimated P-values (especially if P values are small), at the expense of increased computing time and memory usage. If \code{Nsimulations<=0}, then no simulations are performed and the only P-value that is possibly computed is \code{LRT_P1_chi2}, if \code{include_lrt=TRUE}.
}
\item{max_model_runtime}{
Optional numeric, maximum time (in seconds) to allocate for each likelihood evaluation per tree, during fitting. Use this to escape from badly parameterized models during fitting (this will likely cause the affected fitting trial to fail). This argument is passed directly to \code{\link{fit_mk}}.
}
\item{include_lrt}{
Logical, whether to also perform a likelihood ratio test, i.e., fitting a "full"" Markov model where both traits may influence each other's transition rates, and comparing the maximized likelihood of the full model to the maximized likelihood of the null model (where each trait evolves independently). Setting this to \code{TRUE} will yield one more test statistic (LRT, defined as twice the difference between the model log-likelihood) that may be used to assess mutual trait dependencies, but at the expense of substantially increased computing. This is because, in contrast to the other test statistics (for example, the Pearson correlation), evaluating the LRT under the null requires repeatedly fitting the null and full models to simulated data.
}
\item{heterogeneous}{
Logical, whether to use hidden Markov models with rate heterogeneity sensu Boyko & Beaulieu (2023), i.e. where the transition rates between states depend on a hidden binary character that evolves independently according to a separate Markov model. If \code{FALSE}, the classical homogeneous models sensu Pagel (1994) are fitted. A heterogeneous model may be useful for reducing false support for correlations in the presence of rare ``bursts'', at the cost of increased computation time, reduced parameter identifiability and possibly increased type II (false negative) errors.
}
\item{verbose}{
Logical, whether to print out progress messages.
}
\item{verbose_prefix}{
Character, optional line prefix for printed messages, for example to achieve a certain indentation.
}
}
\details{
To compute the statistical significance (P value) of any given correlation coefficient, this function repeatedly simulates hypothetical tip states using the fitted null model, each time recomputing the corresponding "random" correlation coefficients, and then examines the distribution of those random correlation coefficients.

Under the null, each trait evolves along the phylogenetic trees according to a continuous-time Markov chain model ("Mk model"), independently of the other trait. In the homogeneous case (\code{heterogeneous==FALSE}), the transition rates between states are constant in time, and the null model consists of two decoupled Mk models (henceforth ``partial'' null models), one for each trait. For binary traits, for example, this implies 4 fitted rate parameters if both \code{null_rate_model1} and \code{null_rate_model2} as "ARD" (default). In the heterogeneous case sensu Boyko & Beaulieu (2023), transition rates are further influenced by a hidden, independently evolving binary character; for binary traits, this implies up to 10 fitted rate parameters (4+4+2).

The fitted null model is subsequently simulated multiple times to generate hypothetical trait data on the same trees and re-evaluating the various considered test statistics (e.g., Pearson correlations and LRT). If \code{value_probs} and/or \code{inclusion_probs} are given for a trait, simulated tip states are further modified and/or replaced by \code{NA} according to those probabilities.

If a likelihood ratio test is included (option \code{include_lrt}), then a "full" model is also fitted, and its maximized likelihood is compared to the maximized likelihood of the null. In the full model, both traits evolve according to a continuous-time Markov chain, and the state of one trait can impact the transition rates of the other trait depending on the option \code{dependent_trait}. Simultaneous transitions of both traits are not permitted in this model. In the homogeneous case (\code{heterogeneous==FALSE}), this full model thus has Nstates1 x Nstates2 states, corresponding to all possible combinations of the two traits. In the heterogeneous case, the full model has Nstates1 x Nstates2 x 2 states, due to the additional hidden binary character.
For any given state of trait 1 (and hidden binary character, if applicable), the transition rates of trait 2 will be consistent with the provided \code{null_rate_model1}, and vise versa. The statistical significance of the LRT is estimated by re-fitting the null and full model to each simulated dataset, and re-computing the corresponding LRT. Simulation-based LRT is thus the most computationally heavy mode in which this function may be used. The function also returns an asymptotic chi2-approximation for the P-value, which should generally not be used for small datasets nor for the heterogeneous case, and is provided mostly for quick preliminary checks and comparisons with other software. Nevertheless, if you only want the chi2 approximation and not the simulation-based significance, you can set \code{Nsimulations=0} for efficiency.
}
\value{
A named list with the following elements:
\item{success}{
Logical, indicating whether the analyses was performed successfully.
}
\item{error}{
Character, an error message if \code{success=FALSE}. Not relevant if \code{success=TRUE}.
}
\item{Nstates1}{
Integer, the inferred number of categorical states for trait 1.
}
\item{Nstates2}{
Integer, the inferred number of categorical states for trait 2.
}
\item{null_fit1}{
Named list, containing the raw partial null model fit results for trait 1, as returned by \code{\link{fit_mk}}. Mostly for diagnosis purposes. Only included if \code{heterogeneous==FALSE}.
}
\item{null_fit2}{
Similar to \code{null_fit1}, but for trait 2. Only included if \code{heterogeneous==FALSE}.
}
\item{null_fit}{
Named list, containing the raw heterogeneous null model fit results for both traits, as returned by \code{\link{fit_mk}}. Mostly for diagnosis purposes. Only included if \code{heterogeneous==TRUE}.
}
\item{full_fit}{
Named list, containing the raw full model fit results, as returned by \code{\link{fit_mk}}. Mostly for diagnosis purposes. Only included if \code{include_lrt} was \code{TRUE}.
}
\item{null_transition_matrix1}{
Numeric matrix of size Nstates1 x Nstates1, listing the estimated transition rates for trait 1, under the null model. The entry \code{null_transition_matrix1[r,k]} is the transition rate from state r to state k. Only included if \code{heterogeneous==FALSE}.
}
\item{null_transition_matrix2}{
Similar to \code{null_transition_matrix2}, but for trait 2. Only included if \code{heterogeneous==FALSE}.
}
\item{null_transition_matrix}{
Numeric 6-dimensional array of size Nstates1 x Nstates2 x 2 x Nstates1 x Nstates2 x 2, listing the estimated joint transition rates between multi-states (s1,s2,h) under the null model, where s1 enumerates states of trait 1, s2 enumerates states of trait 2, and h enumerates the two states of the hidden rate-modulating character. For example, \code{null_transition_matrix[s1,s2,1,s1,s2,2]} is the transition rate from state 1 to 2 in the hidden character, when trait 1 is in state s1 and trait 2 is in state s2. Only included if \code{heterogeneous==TRUE}.
}
\item{null_Nrates}{
Integer, the total number of free/fitted rate parameters in the null model.
}
\item{full_Nrates}{
Integer, the total number of free/fitted rate parameters in the full model.
}
\item{null_loglikelihood}{
The maximized log-likelihood of the null model.
}
\item{full_loglikelihood}{
The maximized log-likelihood of the full model. Only included if \code{include_lrt} was \code{TRUE}.
}
\item{null_AIC}{
The Akaike Information Criterion of the null model.
}
\item{full_AIC}{
The Akaike Information Criterion of the full model. Only included if \code{include_lrt} was \code{TRUE}.
}
\item{AIC_weights}{
Numeric vector of length 2, listing the AIC weights for the null and full model, respectively (Burnham & Anderson, 2002).
}
\item{Pearson}{
Numeric, the Pearson correlation coefficient between the two traits across all tips and across all trees. For binary traits this is also known as "Matthew's correlation coefficient" or "Phi coefficient". Note that for non-binary traits without a clear quantitative interpretation the Pearson correlation may not be meaningful, and the Spearman correlation or normalized mutual information may be more suitable.
}
\item{Spearman}{
Numeric, the Spearman rank correlation coefficient between the two traits across all tips and across all trees. Note that for non-ordinal traits, i.e., without a clear ranking order, the Spearman correlation may not be meaningful, and the normalized mutual information may be more suitable.
}
\item{NMI}{
Numeric, the mutual information between the two traits across all tips and across all trees, normalized by the square root of the product of entropies. In contrast to the raw mutual information, this normalized version is guaranteed to be between 0 (statistical independence) and 1 (full functional dependence). In contrast to the Pearson and Spearman correlations, the NMI is also suitable for comparing nominal categorical traits, i.e., that lack a quantitative interpretation and ranking order (e.g., beak color or cell shape).
}
\item{LRT}{
The likelihood ratio test statistic, defined as \eqn{2\cdot[\ln(L_{\rm full}) - \ln(L_{\rm null})]}, where \eqn{L_{\rm full}} and \eqn{L_{\rm null}} are the maximized likelihoods of the full and null model, respectively. Only included if \code{include_lrt} was \code{TRUE}.
}
\item{Pearson_P1}{
Numeric, the one-sided statistical significance of the Pearson correlation, i.e., the estimated probability under the null model that the Pearson correlation would be smaller than the value obtained for the input data (if the latter is negative) or larger than the value obtained for the input data (if the latter is positive). Only available if \code{Nsimulations>0}.
}
\item{Pearson_P2}{
Numeric, the two-sided statistical significance of the Pearson correlation, i.e., the estimated probability under the null model that the modulus of the Pearson correlation would be larger than observed for the input data. Only available if \code{Nsimulations>0}.
}
\item{Pearson_mean_null}{
Numeric, the average Pearson correlation expected under the null model. Only available if \code{Nsimulations>0}. Only available if \code{Nsimulations>0}.
}
\item{Pearson_std_null}{
Numeric, the standard deviation of the Pearson correlation under the null model. Only available if \code{Nsimulations>0}.
}
\item{Pearson_sims}{
Numeric vector of length \code{Nsimulations}, Pearson correlation coefficients obtained from the simulations of the null model.
}
\item{Pearson_P1_naive}{
Numeric, the one-sided statistical significance of the Pearson correlation under a naive null hypothesis that ignores evolutionary relationships, i.e., under random permutation of the tips regardless of their location in the trees. This should NOT be used for any serious inference of trait dependencies; it is merely provided for comparison purposes, for example to check how much evolutionary relationships impact significance estimates. Only available if \code{Nsimulations>0}.
}
\item{Pearson_P2_naive}{
Numeric, the two-sided statistical significance of the Pearson correlation under a naive null hypothesis that ignores evolutionary relationships, i.e., under random permutation of the tips regardless of their location in the trees. This should NOT be used for any serious inference of trait dependencies; it is merely provided for comparison purposes, for example to check how much evolutionary relationships impact significance estimates. Only available if \code{Nsimulations>0}.
}
\item{Spearman_P1, Spearman_mean_null, Spearman_std_null, Spearman_sims, Spearman_P1_naive, Spearman_P2_naive}{
Similar to the values returned for the Pearson correlation, but referring to the Spearman rank correlation.
}
\item{NMI_P1, NMI_mean_null, NMI_std_null, NMI_sims, NMI_P1_naive}{
Similar to the values returned for the Pearson correlation, but referring to the normalized mutual information.
}
\item{LRT_P1, LRT_mean_null, LRT_std_null, LRT_sims}{
Similar to the values returned for the Pearson correlation, but referring to the likelihood ratio test statistic.
}
\item{LRT_P1_chi2}{
Numeric, the chi2-approximation analytical of the one-sided significance of the LRT. This approximation is only accurate asymptotically in the limit of infinite data, and only for the homogeneous models (\code{heterogeneous==FALSE}). It is mostly provided for comparison purposes. In contrast to the other P-values, this one is returned even if \code{Nsimulations<=0}, since it does not rely on simulated data. Only included if \code{heterogeneous==FALSE}.
}
\item{LRT_df}{
Integer, the degrees of freedom of the asymptotic chi2-approximation of the LRT. This is the difference in the number of free parameters between the full and null model.
}
}


\author{Stilianos Louca}

\references{
J. D. Boyko, J. M. Beaulieu (2023). Reducing the biases in false correlations between discrete characters. Systematic Biology. 72:476-488.

K. P. Burnham & D. R. Anderson (2002). Model Selection and Multimodel Inference. Springer.

M. Pagel (1994). Detecting correlated evolution on phylogenies: a general method for the comparative analysis of discrete characters. Proceedings of the Royal Society of London B. 255:37-45.

A. Strehl & J. Ghosh (2002). Cluster ensembles - a knowledge reuse framework for combining multiple partitions. Journal of Machine Learning Research. 3:583-617.
}

\seealso{
\code{\link{fit_mk}},
\code{\link{simulate_mk_model}},
\code{\link{correlate_phylo_geodistances}}
}

\examples{
# EXAMPLE 1: Independent traits on one tree, homogeneous models

# generate a tree
Ntips = 100
tree = generate_random_tree(list(birth_rate_factor=0.1),max_tips=Ntips)$tree

# generate states for two independent categorical traits on all tips
# by simulating Mk models
# One trait has 3 states, the other has two states
Q1 = get_random_mk_transition_matrix(Nstates=3, rate_model="ARD", max_rate=1)
Q2 = get_random_mk_transition_matrix(Nstates=2, rate_model="ARD", max_rate=1)
tip_states1 = simulate_mk_model(tree, Q1, Nsimulations=1)$tip_states
tip_states2 = simulate_mk_model(tree, Q2, Nsimulations=1)$tip_states

# analyze correlations between the two traits
correlations = correlate_categorical_traits(trees       = tree,
                                         tip_states1 = tip_states1,
                                         tip_states2 = tip_states2,
                                         Nthreads    = 1)

# print normalized mutual information and its one-sided P value
print(correlations$NMI)
print(correlations$NMI_P1)

######################################################################
# EXAMPLE 2: Co-dependent traits on multiple trees, homogeneous models

# generate multiple trees
Ntips = 100
Ntrees= 3
trees = lapply(seq_len(Ntrees), FUN=function(tr)
          { generate_random_tree(list(birth_rate_factor=0.1),max_tips=Ntips)$tree})

# simulate two binary trates influencing each other's transition rates
Q = get_random_mk_transition_matrix(Nstates=4, rate_model="ARD", max_rate=0.02)
state_combos = lapply(seq_len(Ntrees), FUN=function(tr) 
                  simulate_mk_model(trees[[tr]], Q, Nsimulations=1)$tip_states)
tip_states1 = lapply(seq_len(Ntrees), FUN=function(tr)
                  { 1L + ((state_combos[[tr]]==2) | (state_combos[[tr]]==4))} )
tip_states2 = lapply(seq_len(Ntrees), FUN=function(tr)
                  { 1L + ((state_combos[[tr]]==3) | (state_combos[[tr]]==4))} )

# set some states in the to NA to simulate missing trait observations
tip_states1[[1]][c(1,3)] = NA
tip_states2[[2]][c(5,6)] = NA

# analyze correlations between the two traits
correlations = correlate_categorical_traits(trees        = trees,
                                         tip_states1  = tip_states1,
                                         tip_states2  = tip_states2,
                                         Nthreads     = 1)

# print Pearson correlation and its two-sided P value
print(correlations$Pearson)
print(correlations$Pearson_P2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{phylogeography}
\keyword{dispersal}
