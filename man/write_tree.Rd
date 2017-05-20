\name{write_tree}
\alias{write_tree}
\title{
Write a tree in Newick (parenthetic) format.
}
\description{
Write a phylogenetic tree to a file or a string, in Newick (parenthetic) format. If the tree is unrooted, it is first rooted internally at the first node.
}
\usage{
write_tree(tree, file="", append=FALSE, digits=10)
}

\arguments{
\item{tree}{
A tree of class "phylo".
}
\item{file}{
An optional path to a file, to which the tree should be written. The file may be overwritten without warning. If left empty (default), then a string is returned representing the tree.
}
\item{append}{
Logical, specifying whether the tree should be appended at the end of the file, rather than replacing the entire file (if it exists).
}
\item{digits}{
Number of significant digits for writing edge lengths.
}
}


\details{
This function is comparable to (but typically much faster than) the \code{ape} function \code{write.tree}.
}

\value{
If \code{file==""}, then a string is returned containing the Newick representation of the tree. Otherwise, the tree is directly written to the file and no value is returned.
}


\author{Stilianos Louca}

%\references{
%}

%\seealso{
%}

\examples{
# generate a random tree
tree = generate_random_tree(max_tips=100, birth_rate_intercept=1)

# obtain a string representation of the tree in Newick format
Newick_string = write_tree(tree)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ Prunning }
