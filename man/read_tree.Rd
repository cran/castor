\name{read_tree}
\alias{read_tree}
\title{
Load a tree from a string or file in Newick (parenthetic) format.
}
\description{
Load a phylogenetic tree from a file or a string, in Newick (parenthetic) format. Any valid Newick format is acceptable.
}
\usage{
read_tree(  string  = "", 
            file    = "", 
            include_edge_lengths    = TRUE, 
            include_node_labels     = TRUE, 
            underscores_as_blanks   = FALSE, 
            check_label_uniqueness  = FALSE)
}

\arguments{
\item{string}{
A string containing a single tree in Newick format. Can be used alternatively to \code{file}.
}
\item{file}{
A path to an input text file containing a single tree in Newick format. Can be used alternatively to \code{string}.
}
\item{include_edge_lengths}{
Logical, specifying whether edge lengths (if available) should be included in the generated tree.
}
\item{include_node_labels}{
Logical, specifying whether node labels (if available) should be included in the generated tree.
}
\item{underscores_as_blanks}{
Logical, specifying whether underscores ("_") in tip and node labels should be replaced by spaces (" "). This is common behavior in other tree parsers. In any case, tip and node labels are also allowed to contain explicit whitespace (other than newlines).
}
\item{check_label_uniqueness}{
Logical, specifying whether to check if all tip labels are unique.
}
}


\details{
This function is comparable to (but typically much faster than) the \code{ape} function \code{read.tree}. The function supports trees with monofurcations and multifurcations, trees with or without tip/node labels, and trees with or without edge lengths. The time complexity is linear in the number of edges in the tree.

Either \code{file} or \code{string} must be specified, but not both. The tree may be arbitrarily split across multiple lines, but no other non-whitespace text is permitted in \code{string} or in the input file. Flanking whitespace (space, tab, newlines) is ignored.
}

\value{
A single rooted phylogenetic tree in ``phylo'' format.
}


\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{write_tree}}
}

\examples{
# generate a random tree
tree = generate_random_tree(list(birth_rate_intercept=1),max_tips=100)$tree

# obtain a string representation of the tree in Newick format
Newick_string = write_tree(tree)

# re-parse tree from string
parsed_tree = read_tree(Newick_string)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ Pruning }