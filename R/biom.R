# Functions for working with BIOM tables

# recursively convert a named or non-named list/vector into a biom string
list_to_biom_string = function(L){
	if(is.null(L)){
		string = "null"
	}else if(is.list(L) || (is.vector(L) && (length(L)!=1))){
		if(is.null(names(L))){
			string = paste0("[",paste(sapply(seq_len(length(L)),function(i) list_to_biom_string(L[[i]])), collapse=", "),"]");
		}else{
			Lnames = names(L);
			string = paste0("{",paste(sapply(seq_len(length(L)),function(i) sprintf('"%s": %s',Lnames[i],list_to_biom_string(L[[i]]))), collapse=", "),"}");
		}
	}else if(is.vector(L) && (length(L)==1)){
		# this seems like a scalar
		if(is.integer(L) || is.numeric(L)){
			if(is.nan(L)){
				string = "NaN";
			}else if(is.na(L)){
				string = "null";
			}else if(is.finite(L)){
				string = as.character(L);
			}else{
				string = paste0('"',as.character(L),'"');
			}
		}else{
			string = paste0('"',as.character(L),'"')
		}
	}
	return(string);
}




# turn a matrix into a sparse matrix
# this function can handle huge matrices with more than 2^31 entries (contrary to vanilla R)
# if the input is already a sparse matrix, it is returned as-is (transposed if needed)
make_matrix_sparse = function(input_matrix, transpose=FALSE){
	if(is(input_matrix, "sparseMatrix")){
		if(transpose){
			sparse_matrix = t(input_matrix);
		}else{
			sparse_matrix = input_matrix;
		}
	}else{
		if(is.integer(input_matrix)){
			non_zeros = find_non_zeros_int_CPP(nrow(input_matrix), ncol(input_matrix), input_matrix, transpose=transpose)
		}else{
			non_zeros = find_non_zeros_float_CPP(nrow(input_matrix), ncol(input_matrix), input_matrix, transpose=transpose)		
		}
		if(transpose){
			dims 	 = c(ncol(input_matrix),nrow(input_matrix));
			dimnames = list(colnames(input_matrix), rownames(input_matrix));
		}else{
			dims 	 = c(nrow(input_matrix),ncol(input_matrix));
			dimnames = list(rownames(input_matrix), colnames(input_matrix));
		}
		if(length(non_zeros$rows)==0){
			# treat case of empty matrix in special way (to help Matrix::sparseMatrix determine data type)
			sparse_matrix = Matrix::sparseMatrix(	i			= integer(), 
													j			= integer(), 
													x			= (if(is.integer(input_matrix)) integer() else numeric()),
													dims		= dims,
													dimnames	= dimnames,
													index1		= TRUE,
													repr		= "C")
		}else{
			sparse_matrix = Matrix::sparseMatrix(	i			= unlist(non_zeros$rows+1L), 
													j			= unlist(non_zeros$columns+1L), 
													x			= unlist(non_zeros$values), 
													dims		= dims,
													dimnames	= dimnames,
													index1		= TRUE,
													repr		= "C")
		}
	}
	return(sparse_matrix)
}


# # Create a BIOM table object (of class "biom", BIOM format v1.0.0) in sparse matrix format, based on a provided matrix & row/column metadata
# # This function creates a similar object as biomformat::make_biom, with the following advantages:
# #	It can handle sparse-matrices as input, and generates a sparse biom-format table (biomformat::make_biom automatically creates a dense matrix internally)
# #	It automatically determines the numerical type (float or int)
# make_sparse_biom_table = function(	input_matrix, 			# Either a dense matrix of class matrix or class data.frame, or a sparse matrix of class Matrix::CsparseMatrix. Row & column names must be set.
# 									row_metadata	= NULL,	# Either a data.frame of size NR x NRM (listing NRM metadatas), or a list of size NR, each element of which is a non-named or named list of metadata values for a specific row
# 									column_metadata	= NULL,	# Either a data.frame of size NC x NCM (listing NCM metadatas), or a list of size NC, each element of which is a non-named or named list of metadata values for a specific column
# 									id				= NULL,
# 									type			= NULL,
# 									generated_by	= "castor"){
# 	if(is.data.frame(input_matrix)) input_matrix = as.matrix(input_matrix)
# 	if(!is(input_matrix, "sparseMatrix")) input_matrix = make_matrix_sparse(input_matrix) # ensure matrix is sparse
# 
# 	NR  		 = nrow(input_matrix)
# 	NC  		 = ncol(input_matrix)
# 	row_names 	 = rownames(input_matrix)
# 	column_names = colnames(input_matrix)
# 	
# 	non_zeros 	  = Matrix::summary(input_matrix); # non_zeros will be a data.frame of size NNZ x 3, listing [row,col,nz-value] in each row
# 	NNZ 	 	  = nrow(non_zeros);
# 	non_zeros[,1] = non_zeros[,1] - 1; # make NZ row indices 0-based, since indices are 0-based in internal data format
# 	non_zeros[,2] = non_zeros[,2] - 1; # make NZ column indices 0-based
# 	
# 	# sort non-zero entries in row-major format
# 	non_zeros = non_zeros[order(non_zeros[,1],non_zeros[,2]),]; # sort non-zeros by row, then by column
# 	print(class(non_zeros)) # debug
# 	
# 	# add basic data
# 	print("making list..") # debug
# 	biom_table = list(	id 				= id, 
# 						format 			= "Biological Observation Matrix 1.0.0",
# 						format_url 		= "http://biom-format.org",
# 						generated_by 	= generated_by,
# 						date 			= format(Sys.time(), format = "%Y-%m-%dT%H:%M:%OS6%z"),
# 						matrix_element_type = ifelse((NNZ>0) && is.integer(non_zeros[1,3]),"int","float"),
# 						shape 			= c(NR,NC),
# 						type 			= type,
# 						matrix_type 	= "sparse",
# 						data 			= lapply(seq_len(NNZ), function(n) non_zeros[n,])); # indices are 0-based in internal data format
# 	print("done list..") # debug
# 						
# 	# add row info (names and metadata)
# 	if(is.null(row_names)) row_names = as.character(c(1:NR)-1);
# 	if(is.null(row_metadata)){
# 		biom_table$rows = lapply(1:NR, function(i) list(id=row_names[i], metadata=NULL));
# 	}else if(is.data.frame(row_metadata)){
# 		biom_table$rows = lapply(1:NR, function(i) list(id=row_names[i], metadata=setNames(row_metadata[i,,drop=TRUE],colnames(row_metadata))));
# 	}else{
# 		biom_table$rows = lapply(1:NR, function(i) list(id=row_names[i], metadata=row_metadata[[i]]));
# 	}
# 	
# 	# add column info (names and metadata)
# 	if(is.null(column_names)) column_names = as.character(c(1:NC)-1);
# 	if(is.null(column_metadata)){
# 		biom_table$columns = lapply(1:NC, function(i) list(id=column_names[i], metadata=NULL));
# 	}else if(is.data.frame(column_metadata)){
# 		biom_table$columns = lapply(1:NC, function(i) list(id=column_names[i], metadata=setNames(column_metadata[i,,drop=TRUE],colnames(column_metadata))));
# 	}else{
# 		biom_table$columns = lapply(1:NC, function(i) list(id=column_names[i], metadata=column_metadata[[i]]));
# 	}
# 	class(biom_table) = "biom";
# 	return(biom_table);
# }


# Create a BIOM table object (of class "biom", BIOM format v1.0.0), based on a provided matrix & row/column metadata
# This function creates a similar object as biomformat::make_biom, with the following advantages:
#	It can handle sparse matrices as input, and can generate a sparse biom-format table (biomformat::make_biom automatically creates a dense matrix internally)
#	It automatically determines the numerical type (float or int)
make_biom_table = function(	input_matrix, 			# Either a dense matrix of class matrix or class data.frame, or a sparse matrix of class Matrix::CsparseMatrix. Row & column names must be set.
							row_metadata	= NULL,	# Either a data.frame of size NR x NRM (listing NRM metadatas), or a list of size NR, each element of which is a non-named or named list of metadata values for a specific row
							column_metadata	= NULL,	# Either a data.frame of size NC x NCM (listing NCM metadatas), or a list of size NC, each element of which is a non-named or named list of metadata values for a specific column
							id				= NULL,
							type			= NULL,
							matrix_type		= "auto", # either "dense", "sparse" or "auto", specifying the desired type of the generated biom table. If "auto", the original format of the input matrix is followed.
							generated_by	= NULL){
	NR  		 = nrow(input_matrix)
	NC  		 = ncol(input_matrix)
	row_names 	 = rownames(input_matrix)
	column_names = colnames(input_matrix)

	# construct biom table object's headers
	biom_table = list(	id 					= id, 
						format 				= "Biological Observation Matrix 1.0.0",
						format_url 			= "http://biom-format.org",
						generated_by 		= generated_by,
						date 				= format(Sys.time(), format = "%Y-%m-%dT%H:%M:%OS6%z"),
						matrix_element_type = ifelse(is.integer(input_matrix),"int","float"),
						shape 				= c(NR,NC),
						type 				= type,
						matrix_type 		= ifelse(matrix_type=="auto",ifelse(is(input_matrix, "sparseMatrix"),"sparse","dense"),matrix_type))
# 						data 				= lapply(seq_len(NNZ), function(n) non_zeros[n,])); # indices are 0-based in internal data format

	# construct numerical data part				
	if(is.data.frame(input_matrix)) input_matrix = as.matrix(input_matrix)
	if((matrix_type %in% c("sparse","auto")) && is(input_matrix, "sparseMatrix")){
		# keep matrix sparse
		non_zeros 	  	= Matrix::summary(input_matrix) # non_zeros will be a data.frame of size NNZ x 3, listing [row,col,nz-value] in each row
		non_zeros[,1] 	= non_zeros[,1] - 1L # make NZ row indices 0-based, since indices are 0-based in internal data format
		non_zeros[,2] 	= non_zeros[,2] - 1L # make NZ column indices 0-based
		non_zeros 	  	= non_zeros[order(non_zeros[,1],non_zeros[,2]),] # sort non-zeros by row, then by column
		biom_table$data	= lapply(seq_len(nrow(non_zeros)), function(n) unlist(non_zeros[n,]))
	}else if((matrix_type=="sparse") && (!is(input_matrix, "sparseMatrix"))){
		# turn a dense matrix into sparse
		if(is.integer(input_matrix)){
			non_zeros = find_non_zeros_int_CPP(nrow(input_matrix), NC, input_matrix, transpose=FALSE)
		}else{
			non_zeros = find_non_zeros_float_CPP(nrow(input_matrix), NC, input_matrix, transpose=FALSE)		
		}
		biom_table$data	= lapply(seq_len(length(non_zeros$values)), function(n) c(non_zeros$rows[n], non_zeros$columns[n], non_zeros$values[n]))
	}else if((matrix_type %in% c("dense","auto")) && (!is(input_matrix, "sparseMatrix"))){
		# keep matrix dense
		biom_table$data	= lapply(seq_len(nrow(input_matrix)), function(r) input_matrix[r,,drop=TRUE])
	}else if((matrix_type=="dense") && is(input_matrix, "sparseMatrix")){
		# turn a sparse matrix into dense
		dense_input_matrix = as.matrix(input_matrix)
		biom_table$data = lapply(seq_len(nrow(dense_input_matrix)), function(r) dense_input_matrix[r,,drop=TRUE])
	}

	# add row info (names and metadata)
	if(is.null(row_names)) row_names = as.character(c(1:NR)-1);
	if(is.null(row_metadata)){
		biom_table$rows = lapply(1:NR, function(i) list(id=row_names[i], metadata=NULL));
	}else if(is.data.frame(row_metadata)){
		biom_table$rows = lapply(1:NR, function(i) list(id=row_names[i], metadata=setNames(row_metadata[i,,drop=TRUE],colnames(row_metadata))));
	}else{
		biom_table$rows = lapply(1:NR, function(i) list(id=row_names[i], metadata=row_metadata[[i]]));
	}
	
	# add column info (names and metadata)
	if(is.null(column_names)) column_names = as.character(c(1:NC)-1);
	if(is.null(column_metadata)){
		biom_table$columns = lapply(1:NC, function(i) list(id=column_names[i], metadata=NULL));
	}else if(is.data.frame(column_metadata)){
		biom_table$columns = lapply(1:NC, function(i) list(id=column_names[i], metadata=setNames(column_metadata[i,,drop=TRUE],colnames(column_metadata))));
	}else{
		biom_table$columns = lapply(1:NC, function(i) list(id=column_names[i], metadata=column_metadata[[i]]));
	}
	class(biom_table) = "biom";
	return(biom_table);
}


# write a BIOM table to a file in JSON format (v1.0)
# For a specification of the file format (JSON BIOM v1.0) see: http://biom-format.org/documentation/format_versions/biom-1.0.html
# This function is similar to biomformat::write_biom, with the following advantages:
#	Supports gzipped output files
#	Supports BIOM JSON files >2 GB
#	Is about 10 times faster
#	Correctly saves sample metadata (biomformat::write_biom has a bug) 
#	Can save BIOM tables in sparse format, even if the passed object is dense
# if force_sparce==TRUE, the table is first turned into truly sparse format (all zeros omitted), even if it already was in sparse format (structure-wise). Note that this comes at a computational cost.
# Ideally biom_table is an object of class BIOM as defined by the biomformat package.
write_biom_table = function(biom_table, file_path, force_sparse=FALSE){
	dir.create(dirname(file_path), showWarnings = FALSE, recursive=TRUE);
	if(endsWith(tolower(file_path),".gz")){
		fout = gzfile(file_path, open="w");
	}else{
		fout = file(file_path, open="w");
	}
	
	NR = biom_table$shape[1]
	NC = biom_table$shape[2]
	
	# check matrix type (don't trust $matrix_type attribute) based on $data; this is only possible if NC!=3 (otherwise it's harder to distinguish sparse from dense based on $data alone)
	if(NC!=3){
		if(length(biom_table$data[[1]])==NC){
			biom_table$matrix_type = "dense";
		}else{
			biom_table$matrix_type = "sparse";
		}
	}
	
	# make table sparse if needed
	# even if the table already is sparse, find any spurious zeros in the stored data and remove those as well
	# Note that row & column indices should be 0-based in $data if sparse, consistent with the BIOM JSON file format
	if(force_sparse){
		if(biom_table$matrix_type=="dense"){
			row2nonzero_columns = sapply(seq_len(NR), function(r) which(biom_table$data[[r]]!=0));
			biom_table$data = unlist(lapply(seq_len(NR), function(r) lapply(row2nonzero_columns[[r]], function(k) c(r-1,k-1,biom_table$data[[r]][k]))), recursive=FALSE);
		}else{
			biom_table$data = biom_table$data[sapply(seq_len(length(biom_table$data)), function(i) biom_table$data[[i]][3]!=0)];
		}
		biom_table$matrix_type = "sparse";
	}

	# print header
	cat(sprintf('{"id": %s, "format": "Biological Observation Matrix 1.0.0", "format_url": "http://biom-format.org", "generated_by": %s, "date": "%s", "matrix_element_type": "%s", "shape": [%d, %d], "type": %s, "matrix_type": "%s", ',
				ifelse(is.null(biom_table$id),"null",paste0('"',biom_table$id,'"')),
				ifelse(is.null(biom_table$generated_by),"null",paste0('"',biom_table$generated_by,'"')),
				ifelse(is.null(biom_table$date),format(Sys.time(), format="%Y-%m-%dT%H:%M:%OS6%z"),biom_table$date),
				biom_table$matrix_element_type,biom_table$shape[1],biom_table$shape[2],
				(if(is.null(biom_table$type)) "null" else paste0('"',biom_table$type,'"')),
				biom_table$matrix_type), file=fout)

	# print data in chunks, so as not to exceed R's length limit for strings
	if(biom_table$matrix_type=="dense"){
		# $data is a list of size NR, each element of which is a vector of size NC
		chunk_size = max(1,min(10000,as.integer(1e8/NC))); # number of rows to render and print in each round
		NR_printed = 0;
		while(NR_printed<NR){
			last_r = min(NR_printed+chunk_size, NR); # last row to be printed in this round
			cat(paste0(ifelse(NR_printed>0,",",'"data": ['),paste(sapply(c((NR_printed+1):last_r), function(r) paste0("[",paste(sprintf("%g",biom_table$data[[r]]),collapse=","),"]")),collapse=",")), file=fout);
			NR_printed = last_r;
		}
	}else{
		# $data is a list of size Nnz, each element of which is a triple of format [row,col,non-zero-value]
		# Note that row & column indices are 0-based in $data, consistent with the BIOM JSON file format
		Nnz 			= length(biom_table$data); # number of values explicitly stored in the table, typically equal to the number of non-zero entries (but could in prinziple include zeros)
		chunk_size 		= 100000; # number of non-zero entries ([row,col,value]-tuples) to render and print in each round
		Nnz_printed 	= 0;
		while(Nnz_printed<Nnz){
			last_nz = min(Nnz_printed+chunk_size, Nnz); # last nz-entry to be printed in this round
			cat(paste0(ifelse(Nnz_printed>0,",",'"data": ['),paste(sapply(c((Nnz_printed+1):last_nz), function(i) sprintf("[%d,%d,%g]",biom_table$data[[i]][[1]],biom_table$data[[i]][[2]],biom_table$data[[i]][[3]])),collapse=",")), file=fout);
			Nnz_printed = last_nz;
		}
	}
	cat("],\n", file=fout);

	# print row names and metadata, in chunks
	chunk_size = 10000;
	NR_printed = 0;
	row_metadata = lapply(1:NR, function(r) biom_table$rows[[r]]$metadata);
	while(NR_printed<NR){
		last_r = min(NR_printed+chunk_size, NR); # last row to be printed in this round
		cat(paste0(ifelse(NR_printed>0,",",'"rows": ['),paste(sapply(c((NR_printed+1):last_r), function(r) sprintf('{"id": "%s", "metadata": %s}',biom_table$rows[[r]]$id,ifelse(is.null(row_metadata[[r]]) || (length(row_metadata[[r]])==1 && is.na(row_metadata[[r]])),"null",list_to_biom_string(biom_table$rows[[r]]$metadata)))),collapse=",")), file=fout);
		NR_printed = last_r;
	}
	cat("],\n", file=fout);

	# print column names and metadata, in chunks
	chunk_size = 10000;
	NC_printed = 0;
	column_metadata = lapply(1:NC, function(k) biom_table$columns[[k]]$metadata);
	while(NC_printed<NC){
		last_c = min(NC_printed+chunk_size, NC); # last column to be printed in this round
		cat(paste0(ifelse(NC_printed>0,",",'"columns": ['),paste(sapply(c((NC_printed+1):last_c), function(k) sprintf('{"id": "%s", "metadata": %s}',biom_table$columns[[k]]$id,ifelse(is.null(column_metadata[[k]]) || (length(column_metadata[[k]])==1 && is.na(column_metadata[[k]])),"null",list_to_biom_string(biom_table$columns[[k]]$metadata)))),collapse=",")), file=fout);
		NC_printed = last_c;
	}
	cat("]}\n", file=fout)
	close(fout)
}


# check the basic validity of a BIOM JSON object
check_biom_json_validity = function(json_object){
	expected_keys = c("id", "format", "format_url", "rows", "columns", "matrix_type", "matrix_element_type", "shape", "data")
  
  	missings = which(!(expected_keys %in% names(json_object)))
	if(length(missings)>0) return(list(valid=FALSE, error=sprintf("Missing top-level keys in putative biom file: %s",paste(expected_keys[missings],collapse=", "))))

	if(length(json_object$matrix_element_type)!=1L){
		return(list(valid=FALSE, error=sprintf("Expected exactly 1 value in matrix_element_type, but found %d",length(json_object$matrix_element_type))))
	}else if(!(json_object$matrix_element_type %in% c("int", "float", "unicode"))){
		return(list(valid=FALSE, error=sprintf("matrix_element_type '%s' is unsupported.",json_object$matrix_element_type)))
	}
	
	if(length(json_object$shape)!=2){
		return(list(valid=FALSE, error=sprintf("Invalid biom shape length %d, expected a length of 2",length(json_object$shape))))
	}	
	if(json_object$shape[1]!=length(json_object$rows)) return(list(valid=FALSE, error=sprintf("Number of rows specified by shape element (%d) differs from the number of defined rows (%d)",json_object$shape[1],length(json_object$rows))))
	if(json_object$shape[2]!=length(json_object$columns)) return(list(valid=FALSE, error=sprintf("Number of columns specified by shape element (%d) differs from the number of defined columns (%d)",json_object$shape[2],length(json_object$columns))))
	if(!json_object$matrix_type %in% c("sparse", "dense")) return(list(valid=FALSE, error=sprintf("Invalid matrix type '%s'",json_object$matrix_type)))
	if(!json_object$matrix_element_type %in% c("int", "float", "unicode")) return(list(valid=FALSE, error=sprintf("Invalid matrix element type '%s'",json_object$matrix_element_type)))
	
	lengths = sapply(json_object$data, length)
	if((json_object$matrix_type=="sparse") && (!all(lengths==3L))) return(list(valid=FALSE, error=sprintf("Found data fields of lengths %d - %d, but expected all to have length 3",min(lengths),max(lengths))))
	else if((json_object$matrix_type=="dense") && (!all(lengths==json_object$shape[2]))) return(list(valid=FALSE, error=sprintf("Found data fields of lengths %d - %d, but expected all to have length %d (number of columns)",min(lengths),max(lengths),json_object$shape[2])))
	
	return(list(valid=TRUE))
}


# Read a BIOM table from a JSON string or file
read_biom_table = function(	string		= "", 	# optional character, encoding a BIOM table in JSON format
							file_path	= ""){	# optional character, path to an input BIOM file in JSON format
	# load BIOM file as a single string, if needed
	if(file_path!=""){
		if(string!="") stop("ERROR: Either string or file must be specified, but not both")
		fin  	= open_file(file_path, "rt")
		string 	= suppressWarnings(paste(readLines(fin, file.info(file_path)$size), collapse="\n"))
		close(fin)
	}
	# parse as JSON
	biom_table = jsonlite::fromJSON(string, simplifyVector=FALSE)
	# check if the loaded JSON object has all the basic elements of a BIOM table
	validity = check_biom_json_validity(biom_table)
	if(!validity$valid) stop(sprintf("ERROR: Invalid BIOM table: %s",validity$error))
	biom_table$shape  = unlist(biom_table$shape)
	biom_table$data   = lapply(biom_table$data, FUN=function(record) unlist(record))
	class(biom_table) = "biom"
	return(biom_table)
}


# extract the values of a numerical biom table, in the form of a sparse or dense matrix
# if matrix_type=="auto", then the returned matrix type will be chosen automatically to match the internal data of the biom table (dense or sparse)
# if matrix_type=="sparse" or "dense", the matrix will be converted to (or kept in) sparse or dense format, respectively
# Note that only the numerical data & row names & column names are included in the returned object - row & column metadata are not represented
matrix_from_biom_table = function(	biom_table,					# a "biom" object, for example as generated by read_biom_table() or make_sparse_biom_table()
									include_dimnames = TRUE,	# if TRUE, the returned matrix will have row & column names, as specified in the biom table. If FALSE, the returned matrix will not have row or column names.
									matrix_type		 = "auto"){ # either "dense", "sparse" or "auto", specifying the desired type of the generated matrix. If "auto", the original format of the biom table's data is kept.
	if(include_dimnames){
		dimnames = list(sapply(biom_table$rows, function(row) row$id), sapply(biom_table$columns, function(column) column$id))
	}else{
		dimnames = NULL
	}
	if((matrix_type %in% c("auto","sparse")) && (biom_table$matrix_type=="sparse")){
		# keep sparse matrix format
		if((biom_table$shape[1]==0) || (biom_table$shape[2]==0)){
			# treat case of empty matrix in special way, to help Matrix::sparseMatrix determine data type
			M = Matrix::sparseMatrix(i			= integer(), 
									j			= integer(), 
									x			= (if(biom_table$matrix_element_type=="int") integer() else numeric()),
									dims		= biom_table$shape,
									dimnames	= dimnames,
									index1		= TRUE,
									repr		= "C")
		}else{
			M = Matrix::sparseMatrix(i			= sapply(biom_table$data, FUN=function(d) d[1]+1), # row indices of non-zero entries. Note that in the BIOM format row & column indices are 0-based
									j			= sapply(biom_table$data, FUN=function(d) d[2]+1), # column indices of non-zero entries. Note that in the BIOM format row & column indices are 0-based
									x			= sapply(biom_table$data, FUN=function(d) d[3]),
									dims		= biom_table$shape,
									dimnames	= dimnames,
									index1		= TRUE,
									repr		= "C")
		}
	}else if((matrix_type %in% c("auto","dense")) && (biom_table$matrix_type=="dense")){
		# keep dense matrix format
		M = matrix(unlist(biom_table$data), nrow=biom_table$shape[1], ncol=biom_table$shape[2], dimnames=dimnames)
	}else if(biom_table$matrix_type=="dense"){
		# create dense matrix from sparse representation
		M = matrix((if(biom_table$matrix_element_type=="int") 0L else  0.0), nrow=biom_table$shape[1], ncol=biom_table$shape[2], dimnames=dimnames)
		M[cbind(sapply(biom_table$data, FUN=function(d) d[1]+1), sapply(biom_table$data, FUN=function(d) d[2]+1))] = sapply(biom_table$data, FUN=function(d) d[3])
	}else{
		# create sparse matrix from dense representation
		M = make_matrix_sparse(input_matrix=matrix(unlist(biom_table$data), nrow=biom_table$shape[1], ncol=biom_table$shape[2], dimnames=dimnames), transpose=FALSE)
	}
	return(M)
}



# extract row or column metadata from a biom table, returning it as a list of named lists.
# Hence, metadata[i] is a named list containing the metadata available for the i-th row (if axis==1) or i-th column (if axis==2).
metadata_from_biom_table = function(biom_table,	# a "biom" object, for example as generated by read_biom_table() or make_sparse_biom_table()
									axis,		# either 1 (rows) or 2 (columns)
									only_indices 		= NULL, 	# optional integer vector, with values in 1,..,biom_table$shape[axis], specifying the indices of rows or columns for which to get metadata. If NULL, metadata are extracted for all rows or columns.
									only_metadata_names	= NULL){	# optional character vector, specifying the names of metadata fields to extract. If NULL, all available metadata fields are included. If metadata are missing for some of these field names, NULLs are returned.
	if(is.null(only_indices)) only_indices = seq_len(biom_table$shape[axis])
	if(is.null(only_metadata_names)){
		if(axis==1){
			metadata = lapply(only_indices, function(r) biom_table$rows[[r]]$metadata)
		}else{
			metadata = lapply(only_indices, function(k) biom_table$columns[[k]]$metadata)
		}
	}else{
		if(axis==1){
			metadata = lapply(only_indices, function(r) biom_table$rows[[r]]$metadata[only_metadata_names])
		}else{
			metadata = lapply(only_indices, function(k) biom_table$columns[[k]]$metadata[only_metadata_names])
		}
	}
	return(metadata)
}

