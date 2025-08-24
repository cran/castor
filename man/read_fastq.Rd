\name{read_fastq}
\alias{read_fastq}
\title{
Load a fastq file.
}
\description{
Efficiently load headers, sequences & qualities from a fastq file.
}
\usage{
read_fastq(file,
		   include_headers      = TRUE,
		   include_sequences    = TRUE,
		   include_qualities    = TRUE,
		   include_phred_scores = FALSE,
		   include_error_probs  = FALSE,
		   truncate_headers_at  = NULL,
		   phred_offset         = NULL,
		   max_sequences        = Inf,
		   max_lines            = Inf)
}

\arguments{
\item{file}{
A character, path to the input fastq file. This file may be gzipped with extension ".gz".
}
\item{include_headers}{
Logical, whether to load the headers. If you don't need the headers you can set this to \code{FALSE} for efficiency.
}
\item{include_sequences}{
Logical, whether to load the sequences. If you don't need the sequences you can set this to \code{FALSE} for efficiency.
}
\item{include_qualities}{
Logical, whether to load the raw qualities, encoded as ASCII characters. If you don't need the raw qualities you can set this to \code{FALSE} for efficiency.
}
\item{include_phred_scores}{
Logical, whether to compute and return the Phred quality scores, in the form of integers. These contain the same information as the raw qualities, but converted from character representation to integer scores. Also see option \code{phred_offset}.
}
\item{include_error_probs}{
Logical, whether to compute and return the nominal error probabilities, based on the qualities. The nominal error probability of each nucleobase is computed as \eqn{10^{-Q/10}}, where \eqn{Q} is the Phred score.
}
\item{truncate_headers_at}{
Optional character, needle at which to truncate headers. Everything at and after the first instance of the needle will be removed from the headers.
}
\item{phred_offset}{
Optional integer, Phred offset to assume for converting raw quality characters to Phred scores. If \code{NULL}, this is automatically chosen among either 33 or 64.
}
\item{max_sequences}{
Optional integer, maximum number of sequences to load. Note that in the case of a gzipped input file the whole file is temporarily decompressed (up to \code{max_lines} lines) regardless of \code{max_sequences}.
}
\item{max_lines}{
Optional integer, maximum number of lines to load. Any trailing sequence truncated due to this limit will be discarded.
In contrast to \code{max_sequences}, for gzipped inputs this limit is already applied at the decompression stage, so it is more effective at reducing computing time.
Keep in mind that typically each fastq record (header+sequence+qualities) spans 4 lines, however in some rare cases sequences and/or qualities may be split across multiple lines.
}
}


\details{
This function is a fast and simple fastq loader. It can be used to load entire files into memory, or to only sample a small portion of sequences without reading the entire file (using \code{max_lines}).
}

\value{
A named list with the following elements:
\item{success}{
Logical, indicating whether the file was loaded successfully. If FALSE, then an error message will be specified by the element \code{error}, and all other elements may be undefined.
}
\item{headers}{
Character vector, listing the loaded headers in the order encountered. Only included if \code{include_headers} was \code{TRUE}.
}
\item{sequences}{
Character vector, listing the loaded sequences in the order encountered. Only included if \code{include_sequences} was \code{TRUE}.
}
\item{qualities}{
Character vector, listing the loaded raw qualities in the order encountered. Only included if \code{include_qualities} was \code{TRUE}.
}
\item{phred_scores}{
List of integer vectors, listing the loaded Phred scores in the order encountered. Hence, \code{phred_scores[[k]]} is an integer vector specifying the Phred scores for the k-th sequence. Only included if \code{include_phred_scores} was \code{TRUE}.
}
\item{error_probs}{
List of numeric vectors, listing the loaded error probabilities in the order encountered. Hence, \code{error_probs[[k]]} is a numeric vector specifying the error probabilities for the k-th loaded sequence. Only included if \code{include_error_probs} was \code{TRUE}.
}
\item{Nlines}{
Integer, number of lines encountered.
}
\item{Nsequences}{
Integer, number of sequences loaded.
}
}


\author{Stilianos Louca}

\seealso{
\code{\link{read_fasta}}, 
\code{\link{read_tree}}
}

\examples{\dontrun{
# load a gzipped fastq file, considering only the first 1000 lines
fastq = read_fastq(file="mysequences.fastq.gz", max_lines=1000)

# print the first sequence and its error probabilities
cat(fastq$sequences[1])
print(fastq$error_probs[[1]])
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{fastq}
