# Load the contents of a fastq file
read_fastq = function(	file,						# character, path to the input fasta file. This may be gzipped (with extension .gz).
						include_headers		= TRUE,
						include_sequences	= TRUE,
						include_qualities	= TRUE, 	# whether to also load and return the quality scores in their raw format, i.e. characters as written in the fastq file
						include_phred_scores= FALSE, 	# whether to also return the Phred integer scores corresponding to the loaded quality characters. These are typically integers in the range 20-40.
						include_error_probs	= FALSE,	# whether to also return the nominal error probability at each nucleotide based on the associated quality/Phred scores.
						truncate_headers_at	= NULL, 	# optional needle string, at which to truncate headers (i.e. remove everything at and after the first instance of the needle)
						phred_offset		= NULL,		# optional integer, Phred offset to assume for quality scores. If NULL, this is automatically chosen among either 33 or 64.
						max_sequences		= Inf,		# optional maximum number of sequences to load.
						max_lines			= Inf){		# optional maximum number of lines to read from the input file. Any truncated trailing sequence will be discarded. In contrast to max_sequences, this option is already applied at the decompression stage (for gzipped inputs), so it is more effective at reducing computing time.
	uncompressed_file = ensure_uncompressed(file, max_lines=(max_lines+1))
	results = read_fastq_from_file_CPP(	fastq_path			= uncompressed_file$file_path,
										include_headers		= include_headers,
										include_sequences	= include_sequences,
										include_qualities	= (include_qualities || include_phred_scores || include_error_probs),
										max_sequences		= max_sequences,
										max_lines			= max_lines)
	if(uncompressed_file$was_compressed) unlink(uncompressed_file$file_path) # delete temporary uncompressed input fasta
	if(!results$success) return(list(success=FALSE, error=results$error))
	if(include_headers && (!is.null(truncate_headers_at))){
		results$headers = sapply(seq_len(length(results$headers)), FUN=function(h){ strsplit(results$headers[h],split=truncate_headers_at,fixed=TRUE)[[1]][1] })
	}
	if(include_phred_scores || include_error_probs){
		qualities_int = lapply(results$qualities, FUN=function(qualities0){ as.integer(charToRaw(paste0(strsplit(qualities0, NULL)[[1]], collapse = ""))) }) # convert characters to ASCIIs
		if(is.null(phred_offset)){
			# guestimate the proper PHRED offset
			min_quality_int  = min(sapply(seq_len(min(1000,length(qualities_int))), FUN=function(k){ min(0,qualities_int[[k]]) }))
			mean_quality_int = mean(sapply(seq_len(min(1000,length(qualities_int))), FUN=function(k){ mean(qualities_int[[k]]) }))
			if(min_quality_int<64) phred_offset = 33
			else if(mean_quality_int>83) phred_offset = 64
			else phred_offset = 33 # ambiguous, so just pick the most common one
		}
		phred_scores = lapply(qualities_int, FUN=function(qualities_int0){ qualities_int0-phred_offset })
		if(include_error_probs) error_probs = lapply(phred_scores, FUN=function(phred_scores0){ 10^(-phred_scores0/10) })
	}
	return(list(success		= TRUE,
				headers		= (if(include_headers) results$headers else NULL),
				sequences	= (if(include_sequences) results$sequences else NULL),
				qualities	= (if(include_qualities) results$qualities else NULL),
				phred_scores= (if(include_phred_scores) phred_scores else NULL),
				error_probs	= (if(include_error_probs) error_probs else NULL),
				Nlines		= results$Nlines,
				Nsequences	= results$Nsequences))
}

