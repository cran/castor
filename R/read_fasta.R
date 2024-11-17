# Load the contents of a fasta file
read_fasta = function(	file,						# character, path to the input fasta file. This may be gzipped (with extension .gz).
						include_headers		= TRUE,
						include_sequences	= TRUE,
						truncate_headers_at	= NULL){ # optional needle string, at which to truncate headers (i.e. remove everything at and after the first instance of the needle)
	uncompressed_file = ensure_uncompressed(file)
	results = read_fasta_from_file_CPP(	fasta_path			= uncompressed_file$file_path,
										include_headers		= include_headers,
										include_sequences	= include_sequences)
	if(uncompressed_file$was_compressed) unlink(uncompressed_file$file_path) # delete temporary uncompressed input fasta
	if(!results$success) return(list("success"=FALSE, error=results$error))
	if(include_headers && (!is.null(truncate_headers_at))){
		results$headers = sapply(seq_len(length(results$headers)), FUN=function(h){ strsplit(results$headers[h],split=truncate_headers_at,fixed=TRUE)[[1]][1] })
	}
	return(list(headers		= (if(include_headers) results$headers else NULL),
				sequences	= (if(include_sequences) results$sequences else NULL),
				Nlines		= results$Nlines,
				Nsequences	= results$Nsequences))
}
