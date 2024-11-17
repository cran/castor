\name{read_fasta}
\alias{read_fasta}
\title{
Load a fasta file.
}
\description{
Efficiently load headers & sequences from a fasta file.
}
\usage{
read_fasta(file,
		   include_headers		= TRUE,
		   include_sequences	= TRUE,
		   truncate_headers_at	= NULL)
}

\arguments{
\item{file}{
A character, path to the input fasta file. This may be gzipped (with extension .gz).
}
\item{include_headers}{
Logical, whether to load the headers. If you don't need the headers you can set this to \code{FALSE} for efficiency.
}
\item{include_sequences}{
Logical, whether to load the sequences. If you don't need the sequences you can set this to \code{FALSE} for efficiency.
}
\item{truncate_headers_at}{
Optional character, needle at which to truncate headers. Everything at and after the first instance of the needle will be removed from the headers.
}
}


\details{
This function is a fast and simple fasta loader. Note that all sequences and headers are loaded into memory at once.
}

\value{
A named list with the following elements:
\item{headers}{
Character vector, listing the loaded headers in the order encountered. Only included if \code{include_headers} was \code{TRUE}.
}
\item{sequences}{
Character vector, listing the loaded sequences in the order encountered. Only included if \code{include_sequences} was \code{TRUE}.
}
\item{Nlines}{
Integer, number of lines encountered.
}
\item{Nsequences}{
Integer, number of sequences encountered.
}
}


\author{Stilianos Louca}

\seealso{
\code{\link{read_tree}}
}

\examples{\dontrun{
# load a gzipped fasta file
fasta = read_faste(file="myfasta.fasta.gz")

# print the first sequence
cat(fasta$sequences[1])
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{fasta}
