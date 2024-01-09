# given a list of character states, find the unique states and map them to integers (1,..,Nstates)
map_to_state_space = function(raw_states, 	# vector of state values, not necessarily as characters. May include NAs, which are not included in the constructed state space.
							  fill_gaps=FALSE,
							  sort_order="natural"){
	# determine valid states
	valids = which(!is.na(raw_states))
	
	# make sure values are strings at this point (assumed below), not numbers
	raw_states = as.character(raw_states)

	# determine possible states
	# Make sure that state_names are sorted (e.g. alphabetically or by numeric value)
	# Sorting is needed for algorithms (e.g. "MPR") that assume trait states are subject to an order relation.
	state_names = unique(raw_states[valids])
	if(sort_order=="natural"){
		state_names = naturalsort::naturalsort(state_names)
	}else if(sort_order=="alphabetical"){
		state_names = sort(state_names)
	}else{
		stop(sprintf("ERROR: Unknown sort_order '%s'",sort_order))
	}

	# expand state domain if needed
	if(fill_gaps){
		state_names = fill_gaps_in_integer_state_list(state_names)
	}

	# transform trait states from strings to indices (indexing state space)
	# construct lookup table for mapping state names --> index in 1:Nstates. Caution: This assumes that state_names[] contains strings
	Nstates 				= length(state_names);
	name2index				= setNames(1:Nstates, state_names)
	mapped_states			= rep(NA, length(raw_states))
	mapped_states[valids] 	= name2index[raw_states[valids]] # mapped_states[i] is an integer in 1:Nstates, and state_names[mapped_states[i]] is the state_name of raw_states[i]
	names(mapped_states)	= names(raw_states)
	
	results = list(	Nstates 		= Nstates,			# number of unique states
					state_names 	= state_names,		# original names (strings) of the mapped states
					mapped_states	= mapped_states,	# the transformed states (i.e. with values in 1,..,Nstates). May also contain NAs, where raw_states was NA.
					name2index		= name2index,
					state_values 	= as.numeric(state_names)) # numerical value of the various states. This may not always be meaningful, depending on how the states are defined, e.g., "red", "blue" don't have a numerical interpretation.
	return(results)
}





# non-integer values are kept as is
# note that here a string x is defined as integer if as.character(as.integer(x)) == x
# e.g. c("2", "2.5", "a", "5", "7.0") will become c("2", "3", "4", "5", "2.5", "a", "7.0")
fill_gaps_in_integer_state_list = function(state_names){
	values = suppressWarnings(as.integer(state_names));
	is_int = ((!is.na(values)) & (as.character(values)==state_names))
	int_values = values[is_int]
	str_values = state_names[!is_int]
	min_value = min(int_values);
	max_value = max(int_values);
	return(c(as.character(min_value:max_value), str_values))
}