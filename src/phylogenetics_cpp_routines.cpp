/*
This C++ code library includes routines for efficiently working with huge phylogenetic trees in R (using the Rcpp interface).

Most code supports multifurcating trees, as well as trees containing monofurcations (i.e. some nodes having only one child).
In most cases, the input tree must be rooted.
The computational complexity of most routines is O(Nedges).

For example, this library includes a Weighted Maximum Parsimony (WMPR) ancestral state reconstruction (ASR) algorithm for discrete traits.
The routine performs WMPR ASR on a phylogenetic tree.
It is meant for large (>100,000 tips) trees, structured in the conventional "phylo" format in R.

Similarly, the library contains functions for efficient maximum-likelihood ancestral state reconstruction with fixed-rate continuous-time Markov models (Mk models) via rerooting.

In the R "phylo" format, the tree topology (ignoring branch lengths) is encoded in a 2D array edge[,] of size Nedges x 2, where tree$edge[e,:] encodes the e-th edge, tree$edge[e,1] --> tree$edge[e,2], with tree$edge[e,1] and tree$edge[e,2] being indices in 1:(Ntips+Nnodes).
Note that in C++ (including this code) indices are zero-based (in contrast to R). 
Hence, tips are indexed 0..(Ntips-1) and nodes are indexed 0..(Nnodes-1), and edge[,] has values in 0,..,(Ntips+Nnodes-1)
All arrays passed to this code must be flattened and stored in row-major-format.

General terminology and indexing conventions for trees, as used below:
	A 'node' is always an internal node.
	Any index called 'node' will be within 0:(Nnodes-1)
	Any index called 'tip' will be within 0:(Ntips-1)
	Any index called 'parent' and 'child' will be within 0:(Ntips+Nnodes-1)
	Any index called 'edge' will be within 0:(Nedges-1)
	Any index called 'root' will be within Ntips:(Ntips+Nnodes-1)
	tree_edge[] will always be of size Nedge x 2 (flattened in row-major-format), with entries in 0:(Ntips+Nnodes-1)
	edge_length[] will always be of length Nedges


Stilianos Louca
February 13, 2017
*/

#include <new>
#include <limits>
#include <cmath>
#include <string>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <complex>
#include <Rcpp.h>
#include <time.h>
#include <sys/time.h>

#ifndef INFTY_D
#define INFTY_D numeric_limits<double>::infinity()
#endif

#ifndef NAN_D
#define NAN_D std::numeric_limits<double>::quiet_NaN()
#endif

#ifndef RELATIVE_EPSILON
#define RELATIVE_EPSILON 1e-8
#endif

#ifndef STRANDOM_EPSILON 
#define STRANDOM_EPSILON 1e-30
#endif


#ifdef __MACH__ 
	// headers specific to Mac OS X to substitute for stuff available in Linux, for timing functions
	#include <sys/types.h> 
	#include <mach/mach_time.h> // needed for high-resolution monotonic timer
	
	#include <mach/mach_init.h>
	#include <mach/thread_act.h>
	#include <mach/mach_port.h>
#endif

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
	#ifndef IS_WINDOWS
		#define IS_WINDOWS
	#endif
#endif


typedef std::complex<double> cdouble;

using namespace Rcpp;
using namespace std;

// defines outcome of request for vector field during simulations
typedef enum {
	RequestedDynamicsRateOfChange,
	RequestedDynamicsForceJumpToState,
	RequestedDynamicsInvalidState
} RequestedDynamics;

// define answers to question whether simulation crossed the domain boundary
typedef enum {
	CrossedBoundaryYesButFixedByReducingTimeStep,
	CrossedBoundaryYesButFixedBruteForce,
	CrossedBoundaryNo
} CrossedBoundary;


typedef enum _MathError{
	MathErrorNone = 0,
	MathErrorNoData = -1,
	MathErrorInvalidData = -2,
	MathErrorUndefined = -2,
	MathErrorNoSolution = -3,
	MathErrorDivisionByZero = -4,
	MathErrorUnknownError = -5
} MathError;


// ****************************************************** //
// BASIC AUXILIARY FUNCTIONS

#pragma mark -
#pragma mark Auxiliary functions
#pragma mark 


inline bool XOR(bool a, bool b){
	return ((!a) && b) || (a && (!b));
}

template<class TYPE>
inline TYPE SQ(TYPE value){
	return value * value;
}

template<class TYPE>
inline TYPE Qube(TYPE value){
	return value * value * value;
}

template<class TYPE>
inline TYPE QTR(TYPE value){
	return SQ(SQ(value));
}

template<class TYPE>
inline int sgn(TYPE value){
	return ((value)>0 ? 1 : ((value)<0 ? -1 : 0));
}


void cast_ComplexVector_to_CPP(	const ComplexVector &A,		// (INPUT)
								vector<cdouble>		&B){	// (OUTPUT)
	B.resize(A.size());
	for(long i=0; i<A.size(); ++i){
		B[i] = cdouble(A[i].r, A[i].i);
	}
}


void set_array_to_value(long X[], long N, long value){
	for(long n=0; n<N; ++n) X[n] = value;
}

void set_array_to_value(double X[], long N, double value){
	for(long n=0; n<N; ++n) X[n] = value;
}

template<class TYPE1, class TYPE2>
long find_first_occurrence(const std::vector<TYPE1> &haystack, const TYPE2 needle){
	for(long n=0; n<haystack.size(); ++n){
		if(haystack[n]==needle) return n;
	}
	return -1;
}

template<class TYPE1, class TYPE2>
long find_first_non_occurrence(const std::vector<TYPE1> &haystack, const TYPE2 antineedle){
	for(long n=0; n<haystack.size(); ++n){
		if(haystack[n]!=antineedle) return n;
	}
	return -1;
}


template<class ARRAY_TYPE>
void multiply_array_with_scalar(ARRAY_TYPE &X, const double scalar){
	for(long i=0; i<X.size(); ++i) X[i] *= scalar;
}




template<class ARRAY_TYPE>
void reverse_array(ARRAY_TYPE &X){
	const long N = X.size();
	long scratch;
	for(long n=0; n<(N/2.0); ++n){
		scratch = X[n];
		X[n] = X[N-n-1];
		X[N-n-1] = scratch;
	}
}



template<class ARRAY_TYPE>
double get_array_min(const ARRAY_TYPE &X){
	const long N = X.size();
	if(N==0) return NAN_D;
	double minX = X[0];
	for(long n=0; n<N; ++n){
		if(X[n]<minX) minX = X[n];
	}
	return minX;
}


template<class ARRAY_TYPE>
double get_array_nonzero_min(const ARRAY_TYPE &X){
	const long N = X.size();
	double minX = NAN_D;
	for(long n=0; n<N; ++n){
		if((X[n]!=0) && (isnan(minX) || (X[n]<minX))) minX = X[n];
	}
	return minX;
}


template<class ARRAY_TYPE>
double get_array_max(const ARRAY_TYPE &X){
	const long N = X.size();
	if(N==0) return NAN_D;
	double maxX = X[0];
	for(long n=0; n<N; ++n){
		if(X[n]>maxX) maxX = X[n];
	}
	return maxX;
}



double get_array_min(const std::vector<double> &X, long start_index, long end_index){
	if(end_index<start_index) return NAN_D;
	double minX = X[start_index];
	for(long n=start_index; n<=end_index; ++n){
		if(X[n]<minX) minX = X[n];
	}
	return minX;
}


long vector_sum(const std::vector<long> &values){
	long S = 0;
	for(long i=0; i<values.size(); ++i) S += values[i];
	return S;
}


long vector_sum(const std::vector<char> &values){
	long S = 0;
	for(long i=0; i<values.size(); ++i) S += values[i];
	return S;
}


long vector_count_zeros(const std::vector<long> &values){
	long S = 0;
	for(long i=0; i<values.size(); ++i) S += (values[i]==0 ? 1 : 0);
	return S;
}


template<class ARRAY_TYPE>
long count_values_below_threshold(const ARRAY_TYPE &values, double threshold){
	long N = 0;
	for(long i=0; i<values.size(); ++i){
		if(values[i]<=threshold) ++N;
	}
	return N;
}


template<class TIME_ARRAY>
long integrate1D(const TIME_ARRAY &times, const std::vector<double> &values, const long start, const long end, const bool ignore_inf){
	double S = 0;
	long last_valid_t = -1;
	for(long t=max(0l,start); t<=end; ++t){
		if(std::isnan(values[t]) || (ignore_inf && abs(values[t])==INFTY_D)){
			continue;
		}else if(last_valid_t<0){ 
			last_valid_t = t;
			continue;
		}
		S += (times[t]-times[last_valid_t])*0.5*(values[t]+values[last_valid_t]);
		last_valid_t = t;
	}
	return S;
}


inline double SQR(double x){
	return x*x;
}


// get index of item in haystack[] closest to needle
// assumes that haystack is in ascending order
template<class TYPE>
long get_nearest_index(const std::vector<TYPE> &haystack, const TYPE needle){
	for(long i=0; i<(haystack.size()-1); ++i){
		if(abs(needle-haystack[i+1])>abs(needle-haystack[i])) return i;
	}
	return (haystack.size()-1);
}


long get_nearest_index(const IntegerVector &haystack, const long needle){
	for(long i=0; i<(haystack.size()-1); ++i){
		if(abs(needle-haystack[i+1])>abs(needle-haystack[i])) return i;
	}
	return (haystack.size()-1);
}



// Convert a dense binary matrix from row-major format (i.e. provided as a list of lists of non-zero column entries)
// to column-major format, i.e. return a list of lists of non-zero row entries
// Note that since the matrix is assumed to be binary (0/1), the actual values don't need to be provided and neither are they returned; we only need to keep track of which cells are non-zero
// [[Rcpp::export]]
Rcpp::List dense_binary_matrix_row2column_major_CPP(const long 			NR,
													const long 			NC,
													const Rcpp::List	&dense_rows, 	// (INPUT) list of size NR, with the r-th entry being an integer std::vector of arbitrary size containing the columns with non-zero entries in row r
													const bool			Rindexing){		// (INPUT) if true, then indices listed in dense_rows are 1-based (as opposed to 0-based, as expected in C++), and the indices listed in the returned dense_columns should also be 1-based. If false, then 0-based indexing is used in input and output
	const long index_shift = (Rindexing ? -1 : 0);
	std::vector< std::vector<long> > dense_columns(NC); // dense_columns[c][] will be a list of row indices, so that dense_columns[c][i] is the i-th non-zero entry in column c (matrix[i,c]!=0)
	vector<long> dense_row;
	for(long r=0; r<NR; ++r){
		dense_row = Rcpp::as<vector<long> >(dense_rows[r]);
		for(long i=0; i<dense_row.size(); ++i){
			dense_columns[dense_row[i]+index_shift].push_back(r-index_shift);
		}
	}
	return Rcpp::wrap(dense_columns);
}



template<typename TYPE>
long aux_qsortIndices_partition(const std::vector<TYPE> &values, std::vector<long> &indices, long start, long end){
	const long pi = start+(end-start)/2; //middle
	const TYPE &pv = values[indices[pi]];
	long temp;
	
	// swap
	temp = indices[pi]; indices[pi] = indices[end]; indices[end] = temp;
	
	long si = start;
	for(long i=start; i<end; ++i){
		if(values[indices[i]] < pv){
			// swap
			temp = indices[i]; indices[i] = indices[si]; indices[si] = temp;
			++si;
		}
	}
	// swap
	temp = indices[si]; indices[si] = indices[end]; indices[end] = temp;
	
	return si;
}



template<typename TYPE>
void aux_qsortIndices(const std::vector<TYPE> &values, std::vector<long> &indices, long start, long end){
	if(start<end){
		long p = aux_qsortIndices_partition(values, indices, start, end);
		aux_qsortIndices(values, indices, start, p-1);
		aux_qsortIndices(values, indices, p+1, end);
	}
}



//quick-sort (average order n*log(n)) of indices pointing to values
//sortedIndices[] will contain the original positions of the sorted values in values[], that is
//  values[sortedIndices[0]] <= values[sortedIndices[1]] <= ... 
template<typename TYPE>
void qsortIndices(const std::vector<TYPE> &values, std::vector<long> &sortedIndices){
	sortedIndices.resize(values.size());
	for(long n=0; n<sortedIndices.size(); ++n) sortedIndices[n] = n;
	aux_qsortIndices(values, sortedIndices, 0, sortedIndices.size()-1);
}


// same as qsortIndices(..), but only considering a subset of the values[]
template<typename TYPE>
void qsortIndices(const std::vector<TYPE> &values, const std::vector<long> &onlyIndices, std::vector<long> &sortedIndices){
	sortedIndices = onlyIndices;
	aux_qsortIndices(values, sortedIndices, 0, sortedIndices.size()-1);
}


// same as qsortIndices(..), but only considering a subset of the values[], namely those for which includeIndex[i]==true
template<typename TYPE>
void qsortIndices(const std::vector<TYPE> &values, const std::vector<char> &includeIndex, std::vector<long> &sortedIndices){
	sortedIndices.clear();
	sortedIndices.reserve(values.size());
	for(long i=0; i<values.size(); ++i){ if(includeIndex[i]) sortedIndices.push_back(i); }
	aux_qsortIndices(values, sortedIndices, 0, sortedIndices.size()-1);
}



//returns the positive (negative) modulo M of division numerator/denominator if denominator>0 (denominator<0)
//that is: 	numerator = D * denominator + M, 
//			whereas M in [0,denominator) if denominator>0 or M in (denominator, 0] if denominator<0
//			and D integer
//For example: 	-5 mod +3 = +1 
//				-5 mod -3 = -2
//				+5 mod -3 = -1
double modulo(double numerator, double denominator){
	if((numerator>=0) && (denominator >= 0)){
		return numerator - denominator * floor(numerator/denominator);
	}else if((numerator < 0) && (denominator < 0)){
		return numerator - denominator * floor(numerator/denominator);
	}else if((numerator < 0) && (denominator > 0)){
		return numerator + denominator * ceil(abs(numerator/denominator));
	}else if((numerator >0) && (denominator < 0)){
		return numerator + denominator * ceil(numerator/abs(denominator));
	}
	return 0; //should not have arrived here
}


//Calculates equivalence class of value modulo L=(intervalMax - intervalMin)
//Returns equivalence class using representants from within the interval [intervalMin, intervalMax]
//Use for example as moduloInterval(angle, 0, 2*PI) or moduloInterval(angle, -PI, +PI)
inline double moduloInterval(double value, double intervalMin, double intervalMax){
	return modulo(value - intervalMin, intervalMax - intervalMin) + intervalMin;
}


double smallest_nonzero_time_step(const std::vector<double> &times){
	double S = INFTY_D;
	for(long i=0; i<times.size(); ++i){
		if(times[i+1]>times[i]){
			S = min(S, times[i+1]-times[i]);
		}
	}
	return S;
}





#pragma mark -
#pragma mark Auxiliary functions for debugging
#pragma mark 

// print a 2D matrix (in row-major format) as a python list-of-lists
// Mainly used for debugging
template<class TYPE>
void print_as_python_matrix(const long NR, const long NC, const TYPE A[]){
	Rcout << "[";
	for(long r=0; r<NR; ++r){
		Rcout << (r>0 ? ", [" : "[");
		for(long c=0; c<NC; ++c){
			Rcout << (c>0 ? ", " : "") << A[r*NC+c];
		}
		Rcout << "]";
	}
	Rcout << "]\n";
}


template<class ARRAY_TYPE>
void print_as_python_matrix(const long NR, const long NC, const ARRAY_TYPE &A){
	Rcout << "[";
	for(long r=0; r<NR; ++r){
		Rcout << (r>0 ? ", [" : "[");
		for(long c=0; c<NC; ++c){
			Rcout << (c>0 ? ", " : "") << A[r*NC+c];
		}
		Rcout << "]";
	}
	Rcout << "]\n";
}

template<class ARRAY_TYPE>
void print_as_matrix(const long NR, const long NC, const ARRAY_TYPE &A){
	for(long r=0; r<NR; ++r){
		for(long c=0; c<NC; ++c){
			Rcout << (c>0 ? ", " : "") << A[r*NC+c];
		}
		Rcout << "\n";
	}
	Rcout << "\n";
}


template<class TYPE>
void print_as_matrix(const long NR, const long NC, const TYPE A[]){
	for(long r=0; r<NR; ++r){
		for(long c=0; c<NC; ++c){
			Rcout << (c>0 ? ", " : "") << A[r*NC+c];
		}
		Rcout << "\n";
	}
	Rcout << "\n";
}



#pragma mark -
#pragma mark Random numbers
#pragma mark 


inline double uniformWithinInclusiveRight(double minimum, double maximum){
	return minimum + STRANDOM_EPSILON + (maximum - minimum - STRANDOM_EPSILON)*R::runif(0.0,1.0);
}

inline double uniformWithinInclusiveLeft(double minimum, double maximum){
	return minimum + (maximum - minimum - STRANDOM_EPSILON)*R::runif(0.0,1.0);
}


inline long uniformIntWithin(long minimum, long maximum){
	//return min(maximum, (long) floor(minimum + (maximum-minimum+1)*(double(rand())/RAND_MAX))); // rand() is discouraged by R package builders
	return min(maximum, (long) floor(minimum + (maximum-minimum+1) * R::runif(0.0,1.0)));
}


// draw a standard-normal random variable
double random_standard_normal(){
	return sqrt(-2.0*log(uniformWithinInclusiveRight(0, 1)))*cos(2.0*PI*uniformWithinInclusiveRight(0,1));
}



// generate a sample from an OU process, conditional upon a previous sample
inline double get_next_OU_sample(	double mean,
									double decay_rate,
									double stationary_std,
									double dt,
									double previous){
	// use transition probability density (Gaussian, known analytically)
	const double std 		 = stationary_std * sqrt(1-exp(-2*dt*decay_rate));
	const double expectation = mean + (previous-mean) * exp(-dt*decay_rate);
	return expectation + std*random_standard_normal();
	// Alternative: Use correlation structure, and the fact that OU is Gaussian
	//const double rho = ((1.0/decay_rate)<dt*STRANDOM_EPSILON ? 0 : exp(-dt*decay_rate));
	//return mean*(1-rho) + rho*previous + sqrt(1-rho*rho)*random_standard_normal()*stationary_std;								
}


// generate a sample from a Brownian motion process, conditional upon a previous sample
inline double get_next_BM_sample(	double	diffusivity,
									double	dt,
									double 	previous){
	return previous + sqrt(2 * diffusivity * dt) * random_standard_normal();									
}


// generate a sample from a bounded Brownian motion process (constrained in an interval via reflection), conditional upon a previous sample
// in practice, this is done by first simulating an unbounded BM, and then reflecting at the boundaries as much as needed
double get_next_bounded_BM_sample(	double	diffusivity,
									double	min_state,
									double	max_state,
									double	dt,
									double 	previous){
	// first perform calculations assuming min_state=0, then add it at the end
	double step  = max_state - min_state;
	if(step<=0) return min_state;
	double state = (previous-min_state) + sqrt(2 * diffusivity * dt) * random_standard_normal();
	state = abs(state); // reflection at the origin
	long wrap_count = floor(state/step);
	if(wrap_count%2==0){
		state = state - (wrap_count*step);
	}else if(wrap_count%2==1){
		state = step - (state - (wrap_count*step));
	}
	state += min_state;
	return state;
}




template<class TYPE>
long random_int_from_distribution(const TYPE probabilities[], const long N){
	double p = R::runif(0.0,1.0);
	for(long i=0; i<N; ++i){
		if(p<=probabilities[i]) return i;
		p -= probabilities[i];
	}
	return N-1;
}


// pick an index within 0:(N-1) at probability proportional to weights[index_pool[i]]
// for efficiency, the caller guarantees that total_weight = sum_i weights[index_pool[i]]
long random_int_from_distribution(const std::vector<long> &index_pool, const std::vector<double> &weights, const double total_weight){
	const long N = index_pool.size();
	double p = R::runif(0.0,1.0);
	for(long i=0; i<N; ++i){
		if(p<=weights[index_pool[i]]/total_weight) return i;
		p -= weights[index_pool[i]]/total_weight;
	}
	return N-1;
}



// generate exponentially distributed random variable, with PDF f(x) = lambda*exp(-lambda*x)
double random_exponential_distribution(double lambda){
	return -log(R::runif(0.0,1.0))/lambda;
}


inline bool random_bernoulli(double p){
	//return ((double(rand())/RAND_MAX) <= p); // rand() is discouraged by R package builders
	return (R::runif(0.0,1.0)<=p);
}



// returns the number of seconds since an arbitrary (but consistent) point in the past
// The timer is monotonic (i.e. not affected by manual changes in the system time), and is specific to the calling thread
// Note that for Windows this function is not thread-specific 
double get_thread_monotonic_walltime_seconds(){
	#if __MACH__ 
		mach_timebase_info_data_t info;
		int error_code = mach_timebase_info(&info);
		if (error_code != KERN_SUCCESS) return 0.0;
		return 1e-9 * mach_absolute_time() * double(info.numer)/double(info.denom);
		/*
		// wrong: non-monotonic, i.e. affected by manual changes in system time
		struct timeval T;
		gettimeofday(&T, NULL);
		return double(T.tv_sec) + 1e-6*T.tv_usec;
		*/
	#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__linux__)
		// POSIX code
		// For details on clock_gettime() see: http://www.tin.org/bin/man.cgi?section=3&topic=clock_gettime
		struct timespec T;
		clock_gettime(CLOCK_MONOTONIC, &T); // Note that CLOCK_MONOTONIC_RAW is not available on all Linux distros
		return double(T.tv_sec) + 1e-9*T.tv_nsec;
	#elifdef IS_WINDOWS
		//return GetTickCount()*1e-6; // note that this is not thread-specific, but it's the best you can get on Windows
		return clock()/double(CLOCKS_PER_SEC);
	#else
		return 0; // not implemented for other systems
	#endif
}



// ##########################################################
// MATRIX ALGEBRA

#pragma mark -
#pragma mark Matrix algebra
#pragma mark 





// LU decomposition of square matrix
// Crout’s method with partial pivoting
// matrix[] (input/output) will store the LU decomposition as a return value (in row-major format). This output can be used as input for solving linear equations (Ax=b) or inverting the original matrix
// pivoting_indices[] (output) will store the pivoting indices, i.e. the row permutation effected by the partial pivoting. 
// pivoting_indices[] should already be allocated before calling (size at least N).
// Psign (output) returns the signature of the permutation pivoting_indices, i.e. will be +1 (-1) if the number of row permutations is even (uneven)
template<class TYPE>
bool LUDecomposition(	TYPE			matrix[], 	// (INPUT/OUTPUT) in row major format
						unsigned long 	N, 			// (INPUT) matrix row count = column count
						unsigned long 	pivoting_indices[], 
						int 			&Psign){
	if(N==0) return false;
	long i,imax,j,k;
	double big,dummy,temp;
	TYPE sum,dummyT;
	double *W = new double[N];
	Psign = 1;
	
	for(i=0; i<N; ++i){
		for(j=0, big=0; j<N; ++j){
			if((temp = abs(matrix[i*N+j])) > big){ big = temp; }
		}
		if(big == 0){
			delete[] W;
			return false; //singular matrix
		}
		W[i] = 1.0/big;
	}
	for(j=0; j<N; ++j){
		for(i=0; i<j; ++i){
			sum = matrix[i*N+j];
			for(k=0; k<i; ++k){ sum -= matrix[i*N+k] * matrix[k*N+j]; }
			matrix[i*N+j] = sum;
		}
		big=0.0;
		for(i=j; i<N; ++i){
			sum = matrix[i*N+j];
			for (k=0; k<j; ++k){ sum -= matrix[i*N+k] * matrix[k*N+j]; }
			matrix[i*N+j] = sum;
			if((dummy=W[i]*abs(sum)) >= big){
				big	 = dummy;
				imax = i;
			}
		}
		if(j != imax){
			for(k=0; k<N; ++k){
				dummyT 			= matrix[imax*N+k];
				matrix[imax*N+k] = matrix[j*N+k]; 
				matrix[j*N+k] 	= dummyT;
			}
			Psign = -Psign;
			W[imax] = W[j];
		}
		pivoting_indices[j] = imax;
		if(abs(matrix[j*N+j]) == 0.0){
			matrix[j*N+j] = RELATIVE_EPSILON; //If the pivot element is zero the matrix is singular (at least to the precision of the algorithm). For some applications on singular matrices, it is desirable to substitute EPSILON for zero. 
		}
		if(j < N-1){
			dummyT = 1.0/(matrix[j*N+j]);
			for(i=j+1; i<N; ++i){ matrix[i*N+j] *= dummyT; }
		}
	}
	delete[] W;
	return true;
}


// Solve the set of n linear equations Ax=b for the vector x and a rectangular (non-singular) matrix A.
// Uses forward substitution and back substitution.
// Matrix needs to be in LU decomposition row-major format, as returned by the routine LUDecomposition(..) defined above.
// pivoting_indices[] should store the pivoting indices from the LU decomposition, as returned by LUDecomposition(..)
// b[] (input/output) will return the solution vector x upon completion.
// This routine takes into account the possibility that b[] will begin with many zero elements, so it is efficient for use in matrix inversion.
template<class TYPE>
void LUSolveLinearSystem(	const TYPE			LUmatrix[],  // (INPUT)
							unsigned long 		N, 
							const unsigned long pivoting_indices[], 
							TYPE				b[]){
	long i, ii=-1, ip, j; 
	TYPE sum;
	for(i=0; i<N; ++i){
		ip = pivoting_indices[i];
		sum = b[ip];
		b[ip] = b[i];
		if(ii>=0){
			for(j=ii; j<=i-1; ++j){ sum -= LUmatrix[i*N+j] * b[j]; }
		}else if(abs(sum)>RELATIVE_EPSILON){
			ii = i; 
		}
		b[i] = sum;
	}
	for(i=N-1; i>=0; --i){
		sum = b[i];
		for(j=i+1; j<N; ++j){ sum -= LUmatrix[i*N+j]*b[j]; }
		b[i] = sum/LUmatrix[i*N+i];
	}
}



template<class TYPE>
double errorInLinearSolution(	const TYPE			matrix[], // (INPUT) in row-major format
								unsigned long 		N, 
								const TYPE 			b[], 
								const TYPE 			x[]){
	double error = 0;
	TYPE Ax;
	unsigned long i,j;
	for(i=0; i<N; ++i){
		for(j=0, Ax=0; j<N; ++j){ Ax += matrix[i*N+j]*x[j]; }
		error += SQ(abs(Ax-b[i]));
	}
	return std::sqrt(error);
}



// Iteratively improve solution x[] to linear system Ax=b
// Uses LU-decomposed matrix as well as original matrix
// Only iterates once, so repeat this step as much as needed to achieve the wanted accuracy
// x[] (input/output) should initially store the approximate (but inaccurate) solution. Upon return, it will store the improved solution.
template<class TYPE>
void LUImproveSolutionToLinearSystem(	const TYPE			matrix[], 	//original square matrix (in row-major format)
										const TYPE			LUmatrix[], //LU decomposed matrix in row-major format, as returned by LUDecomposition()
										unsigned long 		N, 
										const unsigned long pivoting_indices[], 	//as returned by LUDecomposition()
										const TYPE 			b[], 
										TYPE 				x[]){ 		//approximate solution we want to improve
	long i,j;
	TYPE s;
	TYPE *r = new TYPE[N];
	for(i=0; i<N; ++i){
		s = -b[i];
		for(j=0; j<N; ++j){ s += matrix[i*N+j] * x[j]; }
		r[i] = s;
	}
	LUSolveLinearSystem(LUmatrix, N, pivoting_indices, r);
	for(i=0; i<N; ++i){ x[i] -= r[i]; }
	delete[] r;
}



// Solve linear system Ax=b for vector x
// scratchSpace[] is only needed for the internal calculations, and should be allocated prior to calling (size at least N*N).
// maxError defines the maximum acceptable L2 norm E(X):=|AX-b|, where X is the approximate solution.
// As long as E(X)>maxError, the solution X is improved iteratively.
// If maxError<=0 or maxImprovements==0, no iterative improvement is performed, i.e. X is set to the solution obtained from regular LU decomposition.
// x[] will contain the solution vector upon return.
// x[] should be allocated prior to calling (size at least N).
// Returns false on error (e.g. when matrix is singular).
template<class TYPE>
bool solveLinearSystem(	const TYPE 			matrix[], // in row major format
						TYPE				scratchSpace[],
						unsigned long 		N,
						const TYPE 			b[],
						double 				maxError, 			//max allowed L2 error |Ax-b|
						unsigned int 		maxImprovements, 	//max allowed number of iterative improvements of solution. Use this to cap computing time if your matrices are unpredictably *pathological*
						TYPE 				x[]){
	if(N==0) return false;
	int Psign;
	long i,j;
	unsigned long *pivoting_indices = new unsigned long[N];
	for(i=0; i<N; ++i){
		for(j=0; j<N; ++j){
			scratchSpace[i*N+j] = matrix[i*N+j];
		}
		x[i] = b[i];
	}
	
	if(!LUDecomposition(scratchSpace, N, pivoting_indices, Psign)){
		delete[] pivoting_indices;
		return false; //error on LU decomposition
	}
	
	LUSolveLinearSystem(scratchSpace, N, pivoting_indices, x);
	
	if((maxError>0) && (maxImprovements>0)){
		int impCount = 0;
		while((errorInLinearSolution(matrix, N, b, x) > maxError) && (impCount<maxImprovements)){
			LUImproveSolutionToLinearSystem(matrix, scratchSpace, N, pivoting_indices, b, x);
			++impCount;
		}
	}
	
	
	delete[] pivoting_indices;
	return true;
}





void get_identity_matrix(	const long 			NR,
							std::vector<double>	&A){	// (INPUT/OUTPUT) will be the identity matrix of size NR x NR, in row-major format
	A.assign(NR*NR,0);
	for(long r=0; r<NR; ++r){
		A[r*NR + r] = 1;
	}
}


// calculate the Hilbert-Schmidt norm (L2-norm) of a matrix: ||A||_2 = SQRT(sum_{r,c}|A_rc|^2)
double get_matrix_norm_L2(	const long 	 NR,	// (INPUT) number of rows & columns of the matrix A
							const double A[]){ 	// (INPUT) matrix of size NR x NR, in row-major format
	double S = 0;
	for(long r=0; r<NR; ++r){
		for(long c=0; c<NR; ++c){
			S += SQR(A[r*NR+c]);
		}
	}
	return sqrt(S);
}


template<class ARRAY_TYPE>
double get_norm_L2_of_vector(const ARRAY_TYPE &X){
	double S=0;
	for(long i=0; i<X.size(); ++i) S += SQR(X[i]);
	return sqrt(S);
}


template<class ARRAY_TYPE>
double get_norm_L2_of_inverted_vector(const ARRAY_TYPE &X){
	double S=0;
	for(long i=0; i<X.size(); ++i) S += SQR(1.0/X[i]);
	return sqrt(S);
}


double get_column_norm_L2(	const long 					NR,		// (INPUT) number of rows & columns in the matrix A
							const long 					c,		// (INPUT) column index for which to calculate the L2 norm
							const std::vector<double>	&A){	// (INPUT) quadratic matrix of size NR x NR, in row-major format
	double S = 0;
	for(long r=0; r<NR; ++r) S += SQR(A[r*NR + c]);
	return sqrt(S);
}


double get_row_norm_L2(	const long 					NR,		// (INPUT) number of rows & columns in the matrix A
						const long 					r,		// (INPUT) row index for which to calculate the L2 norm
						const std::vector<double>	&A){	// (INPUT) quadratic matrix of size NR x NR, in row-major format
	double S = 0;
	for(long c=0; c<NR; ++c) S += SQR(A[r*NR + c]);
	return sqrt(S);
}


// transform a quadratic matrix A using another diagonal matrix D: A --> D*A*D^{-1}  or  D^{-1}*A*D
void diagonally_transform_matrix(	const long	 				NR,				// number of rows & columns in A & D
									const std::vector<double> 	&D,				// array of size NR, diagonal elements of the diagonal transformation matrix D
									const bool					inverse_right,	// (INPUT) if true, then the transformed matrix will be D*A*D^{-1}, otherwise it will be D^{-1}*A*D.
									double 						A[]){			// (INPUT/OUTPUT) 2D array of size NR x NR, in row-major format. Matrix to be transformed in-situ.
	for(long r=0; r<NR; ++r){
		for(long c=0; c<NR; ++c){
			A[r*NR + c] *= (inverse_right ? D[r]/D[c] : D[c]/D[r]);
		}
	}
}


// replace a quadratic matrix A in-situ by its square, A*A
void get_square_matrix(	const long		NR,				// (INPUT) number of rows & columns in the matrix A
						const double	A[],			// (INPUT) 2D matrix of size NR x NR, in row-major format. Matrix to be squared.
						double			Asquared[]){	// (OUTPUT) 2D matrix of size NR x NR, in row-major format, equal to A*A. Must be pre-allocated.
	for(long r=0; r<NR; ++r){
		for(long c=0; c<NR; ++c){
			Asquared[r*NR + c] = 0;
			for(long k=0; k<NR; ++k){
				Asquared[r*NR + c] += A[r*NR + k] * A[k*NR + c];
			}
		}
	}
}






// Calculate the product Y = A * X
template<class TYPE1,class TYPE2,class TYPE3>
void multiply_matrix_with_vector(	const long			NR,
									const long			NC,
									TYPE1				A[],	// (INPUT) array of size NR*NC, in row-major format
									TYPE2				X[],	// (INPUT) pre-allocated array of size NC or greater
									std::vector<TYPE3>	&Y){	// (OUTPUT) product A*X, of size NR
	Y.assign(NR,0);
	for(long r=0; r<NR; ++r){
		for(long c=0; c<NC; ++c){
			Y[r] += A[r*NC+c] * X[c];
		}
	}
}

// Calculate the product Y = X^T * A
template<class TYPE1,class TYPE2,class TYPE3>
void multiply_vector_with_matrix(	const long			NR,
									const long			NC,
									TYPE1				A[],	// (INPUT) array of size NR*NC, in row-major format
									TYPE2				X[],	// (INPUT) pre-allocated array of size NC or greater
									std::vector<TYPE3>	&Y){	// (OUTPUT) product A*X, of size NR
	Y.assign(NR,0);
	for(long r=0; r<NR; ++r){
		for(long c=0; c<NC; ++c){
			Y[r] += X[c] * A[c*NC+r];
		}
	}
}




// multiply matrix A with column-vector X
// matrix is assumed in row-major format
template<class TYPE1, class TYPE2, class TYPE3>
void multiply_matrix_with_vector(	const long					NR,
									const long					NC,
									const std::vector<TYPE1>	&A,		// (INPUT) matrix of size NR*NC, in row-major format
									const std::vector<TYPE2>	&X,		// (INPUT) std::vector of size NC
									std::vector<TYPE3>			&Y){	// (OUTPUT) product A*X, of size NR
	Y.assign(NR,0);
	for(long r=0; r<NR; ++r){
		for(long c=0; c<NC; ++c){
			Y[r] += A[r*NC+c] * X[c];
		}
	}
}



// multiply matrix A with another matrix B
// matrices are assumed to be in row-major format
template<class ARRAY_TYPE1, class ARRAY_TYPE2, class ARRAY_TYPE3>
void multiply_matrices(	const long			NRa,	// (INPUT) number of rows in A
						const long			NCa,	// (INPUT) number of columns in A. Must be equal to the number of rows in B.
						const long			NCb,	// (INPUT) number of columns in B
						const ARRAY_TYPE1	&A,		// (INPUT) 2D matrix of size NRa x NCa, in row-major format
						const ARRAY_TYPE2	&B,		// (INPUT) 2D matrix of size NCa x NCb, in row-major format
						ARRAY_TYPE3			&AB){	// (OUTPUT) 2D matrix of size NRa x NCb, in row-major format, containing the product A*B.
	const long NR = NRa;
	const long NC = NCb;
	AB.resize(NR*NC);
	for(long r=0; r<NR; ++r){
		for(long c=0; c<NC; ++c){
			AB[r*NC + c] = 0;
			for(long k=0; k<NCa; ++k){
				AB[r*NC + c] += A[r*NCa + k] * B[k*NCb + c];
			}
		}
	}
}



// multiply matrix A with another matrix B
// matrices are assumed to be in row-major format
// Special case: cast the output to real valued, by omitting the imaginary part
template<>
void multiply_matrices(	const long					NRa,	// (INPUT) number of rows in A
						const long					NCa,	// (INPUT) number of columns in A. Must be equal to the number of rows in B.
						const long					NCb,	// (INPUT) number of columns in B
						const std::vector<cdouble>	&A,		// (INPUT) 2D matrix of size NRa x NCa, in row-major format
						const std::vector<cdouble>	&B,		// (INPUT) 2D matrix of size NCa x NCb, in row-major format
						std::vector<double>			&AB){	// (OUTPUT) 2D matrix of size NRa x NCb, in row-major format, containing the product A*B.
	const long NR = NRa;
	const long NC = NCb;
	AB.resize(NR*NC);
	for(long r=0; r<NR; ++r){
		for(long c=0; c<NC; ++c){
			cdouble S = 0;
			for(long k=0; k<NCa; ++k){
				S += A[r*NCa + k] * B[k*NCb + c];
			}
			AB[r*NC + c] = S.real();			
		}
	}
}



// multiply matrix A with another diagonal matrix D
// this is the same as multiplying the c-th column in A with the c-th diagonal entry in D
template<class ARRAY_TYPE1, class ARRAY_TYPE2, class ARRAY_TYPE3>
void multiply_matrix_with_diagonal_matrix(	const long			NR,		// (INPUT) number of rows in A
											const long			NC,		// (INPUT) number of columns in A. Must be equal to the number of rows & columns in D.
											const ARRAY_TYPE1	&A,		// (INPUT) 2D matrix of size NR x NC, in row-major format
											const ARRAY_TYPE2	&D,		// (INPUT) 1D array of size NC, containing the diagonal elements of the diagonal matrix D.
											ARRAY_TYPE3			&AD){	// (OUTPUT) matrix of size NR x NC, in row-major format, containing the product A*D.
	AD.resize(NR*NC);
	for(long r=0; r<NR; ++r){
		for(long c=0; c<NC; ++c){
			AD[r*NC + c] = A[r*NC + c] * D[c];
		}
	}
}






// calculate the appropriate base-2 scaling power for the scaling-and-squaring in matrix exponentiation
// and automatically rescale the matrix
void get_scaling_power_for_matrix_exponentiation(	const long 			&NR,				// (INPUT) number of rows & column in the matrix A
													std::vector<double> &A,					// (INPUT/OUTPUT) 2D matrix of size NR x NR, in row-major format. Matrix to be rescaled and replaced by A/2^scaling_power.
													long				&scaling_power,		// (OUTPUT) base-2 power of the scaling applied, i.e. the rescaled matrix is A/2^scaling_power
													double				&rescaled_norm){	// (OUTPUT) norm of the rescaled matrix, i.e. ||A/2^scaling_power||_2
	scaling_power = 0;
	const double norm = get_matrix_norm_L2(NR, &A[0]);
	if(norm<=1.0){
		scaling_power = 0.0; // no rescaling needed
		rescaled_norm = norm;
	}else{
		scaling_power = (long)ceil(log(norm)/log(2.0));
		double factor = pow(0.5, scaling_power);
		multiply_array_with_scalar(A,factor);
		rescaled_norm = norm * factor;
	}
}



// balance a quadratic matrix A using a diagonal transformation matrix D ("balances matrix"): A --> A'=D^{-1}*A*D
// so that the L2-norms of the columns of the transformed A' are roughly the same as the L2 norms of the corresponding rows
// Uses the balancing algorithm proposed by: Rodney James et al (2014). On matrix balancing and eigenvector computation. Algorithm 3.
void balance_matrix_with_diagonal_transformation(	const long			NR,		// (INPUT) number of rows & columns in the matrix A
													std::vector<double>	&A,		// (INPUT/OUTPUT) 2D matrix of size NR x NR, in row-major format. Input matrix to be balanced, i.e. replaced by D^{-1}*A*D
													std::vector<double>	&D){	// (OUTPUT) 1D array of size NR, storing the diagonal elements of the balances matrix D.
	D.assign(NR,1.0);
	const double beta = 2;
	bool converged = false;
	double C,R,S,f;
	while(!converged){
		converged = true;
		for(long i=0; i<NR; ++i){
			C = get_column_norm_L2(NR,i,A);
			R = get_row_norm_L2(NR,i,A);
			if((C<=1e-32*R) || (R<=1e-32*C)){
				// row and/or column is practically zero, so no point in rescaling. Move on to next row & column.
				D[i] = 1;
				continue;
			}
			S = SQR(C) + SQR(R);
			f = 1;
			while(C<R/beta){
				C *= beta;
				R /= beta;
				f *= beta;
			}
			while(C>=R*beta){
				C /= beta;
				R *= beta;
				f /= beta;
			}
			if((SQR(C)+SQR(R))<0.95*S){
				converged = false;
				D[i] *= f;
				for(long j=0; j<NR; ++j){
					A[j*NR + i] *= f;
					A[i*NR + j] /= f;
				}
			}
		}
	}
}



// estimate the minimum number of polynomials A^p/p! to include for exponentiating an irreducible matrix A via Taylor expansion, such that the approximated exp(A) does not contain any zeros.
// For an irreducible matrix A, exp(A) cannot include zeros.
// In Mk-model ancestral state reconstruction algorithms with transition matrix Q, false zeros in exp(Q) can mess up the algorithm
// This function determines the minimum number P (P<=NR+1) such that sum_{p=0}^{P-1} abs(A^p) has positive entries everywhere
// Note that if the matrix is reducible, the returned value would still be <=NR+1.
template<class ARRAY_TYPE>
long min_polynomials_for_positive_exponential_of_irreducible_matrix(const long			NR,		// (INPUT) number of rows & columns of the matrix A
																	const ARRAY_TYPE	&A){	// (INPUT) 2D array of size NR x NR, in row-major format
	std::vector<double> scratch_power1, scratch_power2;
	get_identity_matrix(NR,scratch_power1);
	std::vector<double> sum_of_powers = scratch_power1;
	std::vector<double> &last_power = scratch_power1;
	for(long p=1; p<=NR; ++p){
		if((p%2)==1){
			multiply_matrices(NR, NR, NR, scratch_power1, A, scratch_power2);
			last_power = scratch_power2;
		}else{
			multiply_matrices(NR, NR, NR, scratch_power2, A, scratch_power1);
			last_power = scratch_power1;
		}
		for(long i=0; i<NR*NR; ++i) sum_of_powers[i] += abs(last_power[i]);
		if(find_first_occurrence(sum_of_powers, 0.0)<0) return (p+1);
	}
	return (NR+1); // this may be returned e.g. if the matrix is reducible
}



// calculate polynomials C_p=A^p/p! of some quadratic matrix A and for p=0,..,NP-1
// The number of polynomials NP is determined based on the requested accuracy, i.e. such that ||exp(A)-sum_{p=0}^{NP-1}C_p||<epsilon
// Optionally, a scalar scaling factor T may be provided, so that T*A is considered in all calculations instead of A.
// This function is meant to be used in conjunction with the function get_matrix_exponential_using_polynomials(..).
void calculate_matrix_polynomials(	const long					NR,					// (INPUT) number of rows & columns of the matrix A
									const std::vector<double>	&A,					// (INPUT) 2D array of size NR x NR, in row-major format
									const double				T,					// (INPUT) scalar scaling factor for input matrix (i.e. use T*A instead of A in all calculations). Set to 1.0 for no rescaling.
									const double				epsilon,			// (INPUT) norm threshold for calculated polynomials C_p=A^p/p!, i.e. stop calculating polynomials as soon as ||exp(A)-sum_{p=0}^{NP-1}C_p||<epsilon. Norm refers to the Hilbert-Schmidt L2 norm.
									const long					NPmin,				// (INPUT) minimum number of polynomials to calculate if possible (including A^0), regardless of the pursued accuracy epsilon. For sparse Markov transition matrix it is recommended to set this to NR+1, so that the matrix exponential does not contain zeros that it shouldn't contain (assuming A is irreducible). The presence of false zeros in exp(A) can mess up ancestral state reconstruction algorithms.
									const long					NPmax,				// (INPUT) maximum possible number of polynomials to calculate (including A^0), regardless of the pursued accuracy epsilon. Used as safety vault, but may break the guaranteed accuracy.
									std::vector<double>			&polynomials,		// (OUTPUT) array of size NP x NR x NR, containing the polynomials C_0,C_1,..,C_{NP-1} in layer-row-major format, i.e. such that polynomials[p*NR*NR + r*NR + c] is the (r,c)-th element of the p-th polynomial
									std::vector<double>			&polynomial_norms,	// (OUTPUT) Hilbert-Schmidt (L2) norms ||C_p||_2 of the polynomials
									long						&NP){				// (OUTPUT) number of polynomials calculated
	const double base_norm = abs(T)*get_matrix_norm_L2(NR,&A[0]);
	
	// set first polynomial to identity matrix
	get_identity_matrix(NR,polynomials);
	polynomial_norms.assign(1, get_matrix_norm_L2(NR,&polynomials[0]));
	NP = 1;
	double suzuki_error = base_norm*exp(base_norm);
	for(long p=1; p<NPmax; ++p){
		NP = p+1;
		// calculate polynomial C_p, resize storage array if needed
		polynomials.resize(NP*NR*NR,0);
		for(long r=0; r<NR; ++r){
			for(long c=0; c<NR; ++c){
				const long pentry = p*NR*NR + r*NR + c;
				polynomials[pentry] = 0;
				for(long k=0; k<NR; ++k){
					polynomials[pentry] += polynomials[(p-1)*NR*NR + r*NR + k] * (T/p) * A[k*NR + c];
				}
			}
		}
		// calculate norm
		polynomial_norms.resize(NP);
		polynomial_norms[p] = get_matrix_norm_L2(NR, &polynomials[p*NR*NR]);
		// check accuracy
		suzuki_error *= base_norm/NP; // suzuki_error = ||A||^NP * exp(||A||)/(NP!)
		if(p>=NPmin-1){
			const double R = base_norm/NP; // it is guaranteed that the norm of any polynomial p>=NP is not greater than polynomial_norms[NP-1] * R^(p-NP+1)
			if(R>=0.999) continue; // it is guaranteed that eventually R will fall below 1 (i.e. for sufficiently large NP)
			const double max_global_error = polynomial_norms[NP-1] * R/(1-R); // it is guaranteed that sum_{p=NP}^INF polynomial_norm[p] is less than max_global_error
			if(min(suzuki_error,max_global_error)<epsilon) break;
		}
	}
}




// calculate exp(tau*A) for some quadratic matrix A and some scalar tau
// provided that the polynomials C0=A^0/0!, C1=A^1/1!, C2=A^2/2!, ... have been pre-computed (up to sufficient power)
// these polynomials can be calculated e.g. using the function calculate_matrix_polynomials(..)
// For example, if the polynomials were calculated using calculate_matrix_polynomials(T*A,epsilon), 
//		then get_matrix_exponential_using_polynomials(tau,epsilon) will yield exp(tau*T*A) with an error of less than epsilon (in terms of the Hilbert-Schmidt L2-norm as well as in terms of the entrywise maximum norm).
// It is recommended to only call this function with |tau|<=1 and prepare the polynomials for some rescaled T*A such that tau*T covers the original desired range of values.
// This way the accuracy threshold (epsilon) can be guaranteed.
// Note that in principle you can provide any epsilon here (i.e. regardless of the epsilon used to calculate the polynomials), but then there may not be enough polynomials available to achieve this accuracy.
void get_matrix_exponential_using_polynomials(	const long 					NR,					// (INPUT) number of rows & columns of the matrix to be exponentiated
												const long					NP,					// (INPUT) number of available pre-computed polynomials of the matrix				
												const std::vector<double> 	&polynomials,		// (INPUT) array of size NP x NR x NR, containing the pre-computed polynomials of matrix: Cp:=A^p/p! for p=0,..,NP-1. Polynomials are stored in layer-row-major format, polynomials[p*NR*NR + r*NR + c] is (r,c)-th-element of A^p/p!
												const std::vector<double>	&polynomial_norms,	// (INPUT) array of size NP, containing the Hilbert-Schmidt L2 norm of each polynomial Cp, ||Cp||_2. Used to reduce the number of incorporated polynomials to the only the necessary number (for optimization reasons).
												const double 				tau,				// (INPUT) scaling factor in exponent
												const double				epsilon,			// (INPUT) maximum allowed error in terms of the L2 norm. Inclusion of polynomials stops as soon as the error (per component) can be guaranteed to be below epsilon.
												const long					NPmin,				// (INPUT) minimum number of polynomials to include if possible (including A^0), regardless of the pursued accuracy epsilon. For sparse Markov transition matrix it is recommended to set this to NR+1, so that the matrix exponential does not contain zeros that it shouldn't contain (assuming A is irreducible). The presence of false zeros in exp(A) can mess up ancestral state reconstruction algorithms.
												std::vector<double>			&exponential){		// (OUTPUT) array of size NR x NR, containing the exponentiated matrix exp(tau*A), in row-major format.
	const double base_norm = abs(tau)*polynomial_norms[1];

	exponential.assign(NR*NR,0);
	double ptau = 1.0;
	double suzuki_error = base_norm*exp(base_norm);
	for(long p=0; p<NP; ++p, ptau*=tau){
		for(long r=0; r<NR; ++r){
			for(long c=0; c<NR; ++c){
				exponential[r*NR + c] += ptau * polynomials[p*NR*NR + r*NR + c];
			}
		}

		// check accuracy
		suzuki_error *= base_norm/(p+1.0); // suzuki_error = ||A||^NP * exp(||A||)/(NP!)
		if(p>=NPmin-1){
			const double R = base_norm/(p+1.0);
			if(R>=0.999) continue; // it is guaranteed that eventually R will fall below 1 (i.e. for sufficiently large p)
			const double max_global_error = abs(ptau)*polynomial_norms[p] * R/(1-R); // it is guaranteed that sum_{k=p+1}^INF ||tau^k*polynomial[k]|| is less than max_global_error
			if(min(suzuki_error,max_global_error)<epsilon) break;
		}
	}
}






// Calculate matrix polynomials similarly to above, but after applying a "balancing transformation" A-->D^{-1}*A*D with some appropriately chosen diagonal matrix D ("balances matrix")
// The matrix D is chosen so that the column L2-norms in D^{-1}*A*D are roughly equal to its row L2-norms
// The returned polynomials are then actually the polynomials of D^{-1}*(T*A)*D
// The error threshold epsilon is automatically adjusted internally to account for the transformation D
// This function is meant to be used in conjunction with the function get_matrix_exponential_using_balanced_polynomials(..)
// For an example of its use see the routine exponentiate_matrix_for_multiple_scalings_CPP(..) below
void calculate_balanced_matrix_polynomials(	const long			NR,					// (INPUT) number of rows & columns of the matrix A
											std::vector<double>	A,					// (INPUT) 2D array of size NR x NR, in row-major format
											const double		T,					// (INPUT) scalar scaling factor for input matrix (i.e. use T*A instead of A in all calculations). Set to 1.0 for no rescaling.
											double				epsilon,			// (INPUT) norm threshold for calculated polynomials C_p=A^p/p!, i.e. stop calculating polynomials as soon as ||exp(A)-sum_{p=0}^{NP-1}C_p||<epsilon. Norm refers to the Hilbert-Schmidt L2 norm.
											const long			NPmin,				// (INPUT) minimum number of polynomials to calculate if possible (including A^0), regardless of the pursued accuracy epsilon. For sparse Markov transition matrix it is recommended to set this to NR+1, so that the matrix exponential does not contain zeros that it shouldn't contain (assuming A is irreducible). The presence of false zeros in exp(A) can mess up ancestral state reconstruction algorithms.
											const long			NPmax,				// (INPUT) maximum possible number of polynomials to calculate, regardless of the pursued accuracy epsilon. Used as safety vault, but may break the guaranteed accuracy.
											std::vector<double>	&polynomials,		// (OUTPUT) array of size NP x NR x NR, containing the polynomials C_0,C_1,..,C_{NP-1} in layer-row-major format, i.e. such that polynomials[p*NR*NR + r*NR + c] is the (r,c)-th element of the p-th polynomial
											std::vector<double>	&polynomial_norms,	// (OUTPUT) Hilbert-Schmidt (L2) norms ||C_p||_2 of the polynomials
											long				&NP,				// (OUTPUT) number of polynomials calculated
											std::vector<double>	&balances,			// (OUTPUT) 1D array of size NR, storing the diagonal elements of a calculated diagonal matrix D that was applied for balancing A prior to exponentiation. That is, the returned polynomials will be those for D^{-1}*A*D instead of A.	
											long				&scaling_power){	// (OUTPUT) base-2 scaling power used to rescale the matrix prior to calculating the polynomials (according to the scaling-and-squaring method). This rescaling will need to be corrected for (via repeated squaring) when calculating exponentials with the polynomials
	multiply_array_with_scalar(A, T);	
	balance_matrix_with_diagonal_transformation(NR, A, balances);
	
	// adjust error threshold (make more strict) to account for balancing (A -> D^{-1}*A*D)
	// in the worst-case scenario the error can be emplified by a factor ||D||_2 * ||D^{-1}||_2
	epsilon /= get_norm_L2_of_vector(balances) * get_norm_L2_of_inverted_vector(balances);
	
	double rescaled_norm;
	get_scaling_power_for_matrix_exponentiation(NR, A, scaling_power, rescaled_norm);
	
	// adjust error threshold (make more strict) to account for rescaling
	// in the worst-case scenario, the error can be amplified by (roughly) 2^scaling_power * ||exp(rescaled_A)||_2 <= 2^scaling_power * exp(rescaled_norm)
	epsilon /= pow(2.0,scaling_power) * exp(rescaled_norm);

	calculate_matrix_polynomials(	NR,
									A,
									1.0,
									epsilon,
									NPmin,
									NPmax,
									polynomials,
									polynomial_norms,
									NP);
}




// This function is meant to be used in conjunction with the function calculate_balanced_matrix_polynomials(..)
// For an example of its use see the routine exponentiate_matrix_for_multiple_scalings_CPP(..) below
void get_matrix_exponential_using_balanced_polynomials(	const long 					NR,					// (INPUT) number of rows & columns of the matrix to be exponentiated
														const long					NP,					// (INPUT) number of available pre-computed polynomials of the matrix				
														const std::vector<double> 	&polynomials,		// (INPUT) array of size NP x NR x NR, containing the pre-computed polynomials of matrix: Cp:=A^p/p! for p=0,..,NP-1. Polynomials are stored in layer-row-major format, polynomials[p*NR*NR + r*NR + c] is (r,c)-th-element of A^p/p!
														const std::vector<double>	&polynomial_norms,	// (INPUT) array of size NP, containing the Hilbert-Schmidt L2 norm of each polynomial Cp, ||Cp||_2. Used to reduce the number of incorporated polynomials to the only the necessary number (for optimization reasons).
														double		 				tau,				// (INPUT) scaling factor in exponent
														double						epsilon,			// (INPUT) maximum allowed error in terms of the L2 norm. Inclusion of polynomials stops as soon as the error (per component) can be guaranteed to be below epsilon.
														const long					NPmin,				// (INPUT) minimum number of polynomials to include if possible (including A^0), regardless of the pursued accuracy epsilon. For sparse Markov transition matrix it is recommended to set this to NR+1, so that the matrix exponential does not contain zeros that it shouldn't contain (assuming A is irreducible). The presence of false zeros in exp(A) can mess up ancestral state reconstruction algorithms.
														const std::vector<double>	&balances,			// (INPUT) 1D array of size NR, storing the diagonal elements of a diagonal matrix D that was applied for balancing A prior to polynomial calculation. This transformation will now be reversed after exponentiation.
														long						scaling_power,		// (INPUT) base-2 scaling power that was applied to the matrix prior to calculating the polynomials. This scaling will now be reversed after exponentiation, via repeated squaring.
														std::vector<double>			&exponential){		// (OUTPUT) array of size NR x NR, containing the exponentiated matrix exp(tau*A), in row-major format.
	
	// reverse some of the scaling beforehand as permitted (this operation is faster than repeated squaring afterwards)
	// make sure tau stays below 1.0 and we only rescale by at most 2^scaling_power
	if((tau<1.0) && (scaling_power>0)){
		const long power_decrease = min(scaling_power,(long)floor(log(1.0/tau)/log(2.0)));
		scaling_power -= power_decrease;
		tau *= pow(2.0, power_decrease);
	}
	
	// correct error threshold to account for balancing and rescaling
	// in the worst-case scenario the error can be emplified due to balancing by a factor ||D||_2 * ||D^{-1}||_2
	// in the worst-case scenario, the error can be amplified due to scaling by (roughly) 2^scaling_power * ||exp(rescaled_A)||_2 <= 2^scaling_power * exp(rescaled_norm)
	epsilon /= (get_norm_L2_of_vector(balances) * get_norm_L2_of_inverted_vector(balances));
	epsilon /= pow(2.0,scaling_power) * exp(tau*polynomial_norms[1]);
	
	// calculate exponential of balanced & rescaled matrix
	get_matrix_exponential_using_polynomials(	NR,
												NP,
												polynomials,
												polynomial_norms,
												tau,
												epsilon,
												NPmin,
												exponential);
												
	// reverse scaling by repeated squaring (move last version between exponential[] and scratch[], since squaring can't be done in-situ).
	std::vector<double> scratch(NR*NR);
	for(long i=0; i<scaling_power; ++i){
		if((i%2)==0) get_square_matrix(NR, &exponential[0], &scratch[0]);
		else get_square_matrix(NR, &scratch[0], &exponential[0]);
	}
	if((scaling_power%2)==1) exponential = scratch;
	
	// reverse balancing in-situ
	diagonally_transform_matrix(NR, balances, true,	&exponential[0]);
}





// Calculate the exponential of a real-valued matrix A using its eigendecomposition.
// EVmatrix contains eigenvectors of A in its columns, hence inverse_EVmatrix * A * EVmatrix is a diagonal matrix D with eigenvalues of A in the diagonal
void get_matrix_exponential_using_eigendecomposition(	const long 					NR,					// (INPUT) number of rows & columns of A
														const std::vector<cdouble>	&eigenvalues,		// (INPUT) 1D array of size NR, listing the complex eigenvalues of A
														const std::vector<cdouble>	&EVmatrix,			// (INPUT) 2D array of size NR x NR, in row-major format, whose columns are the eigenvectors of A
														const std::vector<cdouble>	&inverse_EVmatrix,	// (INPUT) 2D array of size NR x NR, in row-major format, the inverse of EVmatrix
														double		 				tau,				// (INPUT) scaling factor in exponent, i.e. calculate exp(tau*A).
														std::vector<cdouble> 		&scratch,			// (SCRATCH SPACE)
														vector<double>				&exponential){		// (OUTPUT) 2D array of size NR x NR, containing the exponentiated matrix exp(tau*A), in row-major format.
	// Note that inverse_EVmatrix * A * EVmatrix is a diagonal matrix D with entries eigenvalues in the diagonal
	// Hence, first exponentiate D and then reverse diagonalization
	std::vector<cdouble> exponentiated_eigenvalues(NR);
	for(long i=0; i<NR; ++i) exponentiated_eigenvalues[i] = exp(tau * eigenvalues[i]);

	// reverse diagonalization
	scratch.resize(NR*NR);
	exponential.resize(NR*NR);
	multiply_matrix_with_diagonal_matrix(NR, NR, EVmatrix,exponentiated_eigenvalues, scratch);
	multiply_matrices(NR, NR, NR, scratch, inverse_EVmatrix, exponential);

}





// calculate the exponential exp(A) of some quadratic matrix A
// Returns the exponential as a flattened array of size NR x NR in row-major format, i.e. with exponential[r*NR + c] being the (r,c)-th entry of exp(A)
void exponentiate_matrix(	const long	 				NR,				// (INPUT) number of rows & columns of the matrix A
							const std::vector<double>	&A,				// (INPUT) 2D array of size NR x NR, in row-major format
							const double				epsilon,		// (INPUT) Absolute error threashold for the approximation of exp(T*A), in terms of the Hilbert-Schmidt L2 norm.
							const long					NPmin,			// (INPUT) minimum number of polynomials to include (including A^0), regardless of the pursued accuracy epsilon. For sparse Markov transition matrix it is recommended to set this to NR+1, so that the matrix exponential does not contain zeros that it shouldn't contain (assuming A is irreducible). The presence of false zeros in exp(A) can mess up ancestral state reconstruction algorithms.
							const long					NPmax,			// (INPUT) maximum possible number of polynomials to calculate (RAM required will scale linearly with NPmax * NR^2). Used as safety vault, but may break the guaranteed accuracy.
							const bool					enforce_probability_matrix, // (INPUT) if true, then the sum along each column is enforced to be 1
							std::vector<double>			&exponential){	// (OUTPUT) exponentiated matrix exp(A); 2D array of size NR x NR, in row-major format	
	
	// prepare data structures for exponentiations of transition matrix
	std::vector<double> polynomials, polynomial_norms, balances;
	long Npolynomials, scaling_power;
	calculate_balanced_matrix_polynomials(	NR,
											std::vector<double>(A.begin(), A.end()),
											1.0,
											epsilon,
											NPmin,
											NPmax,
											polynomials,
											polynomial_norms,
											Npolynomials,
											balances,
											scaling_power);
																
	// calculate exponential using the pre-computed polynomials
	get_matrix_exponential_using_balanced_polynomials(	NR,
														Npolynomials,
														polynomials,
														polynomial_norms,
														1.0,
														epsilon,
														NPmin,
														balances,
														scaling_power,
														exponential);
	
	if(enforce_probability_matrix){
		for(long c=0; c<NR; ++c){
			double col_sum = 0;
			for(long r=0; r<NR; ++r){
				exponential[r*NR + c] = max(0.0, exponential[r*NR + c]);
				if(r!=c) col_sum += exponential[r*NR + c];
			}
			exponential[c*NR + c] = 1 - col_sum;
		}
	}
}


// calculate the exponential exp(T*A) for some quadratic matrix A and for a large number of scalar scaling factors T
// Returns the exponentials (one per scaling) as a flattened array of size NS x NR x NR in layer-row-major format, i.e. with exponentials[s*NR*NR + r*NR + c] being the (r,c)-th entry of exp(scalings[s]*A)
// This routine is most efficient when T is very large.
// [[Rcpp::export]]
NumericVector exponentiate_matrix_for_multiple_scalings_CPP(const long	 			NR,			// (INPUT) number of rows & columns of the matrix A
															const NumericVector		&A,			// (INPUT) 2D array of size NR x NR, in row-major format
															const NumericVector		&scalings,	// (INPUT) 1D numeric array of size NS containing scalar scaling factors
															const double			epsilon,	// (INPUT) Absolute error threashold for the approximation of exp(T*A), in terms of the Hilbert-Schmidt L2 norm.
															const long				NPmin,		// (INPUT) minimum number of polynomials to include (including A^0), regardless of the pursued accuracy epsilon. For sparse Markov transition matrix it is recommended to set this to NR+1, so that the matrix exponential does not contain zeros that it shouldn't contain (assuming A is irreducible). The presence of false zeros in exp(A) can mess up ancestral state reconstruction algorithms.
															const long				NPmax,		// (INPUT) maximum possible number of polynomials to calculate (RAM required will scale linearly with NPmax * NR^2). Used as safety vault, but may break the guaranteed accuracy.
															const bool				enforce_probability_matrix){ // (INPUT) if true, then the sum along each column is enforced to be 1
	const double max_scaling = get_array_max(scalings);
	const long NS = scalings.size();
	
	// prepare data structures for exponentiations of transition matrix
	std::vector<double> polynomials, polynomial_norms, balances;
	long Npolynomials, scaling_power;
	calculate_balanced_matrix_polynomials(	NR,
											std::vector<double>(A.begin(), A.end()),
											max_scaling,
											epsilon,
											NPmin,
											NPmax,
											polynomials,
											polynomial_norms,
											Npolynomials,
											balances,
											scaling_power);
																
	// calculate exponentials using the pre-computed polynomials
	NumericVector exponentials(NS*NR*NR);
	std::vector<double> scratch_exponential;
	for(long s=0; s<NS; ++s){
		get_matrix_exponential_using_balanced_polynomials(	NR,
															Npolynomials,
															polynomials,
															polynomial_norms,
															scalings[s],
															epsilon,
															NPmin,
															balances,
															scaling_power,
															scratch_exponential);
		for(long r=0; r<NR; ++r){
			for(long c=0; c<NR; ++c){
				exponentials[s*NR*NR + r*NR + c] = scratch_exponential[r*NR + c];
			}
		}
		if(enforce_probability_matrix){
			for(long c=0; c<NR; ++c){
				double col_sum = 0;
				for(long r=0; r<NR; ++r){
					exponentials[s*NR*NR + r*NR + c] = max(0.0, exponentials[s*NR*NR + r*NR + c]);
					if(r!=c) col_sum += exponentials[s*NR*NR + r*NR + c];
				}
				exponentials[s*NR*NR + c*NR + c] = 1 - col_sum;
			}
		}
	}
	return exponentials;
}





#pragma mark -
#pragma mark Time series analysis
#pragma mark 



//Returns linear regression fit (y=trafo*x + offset) for points passed
//Algorithm tries to minimize \sum_i (y_i-y_{fit})^2
MathError fitLinearRegression(double pointsX[], double pointsY[], long count, double &slope, double &offset){
	if(count == 0){
		slope = offset = NAN_D;
		return MathErrorNoData;
	}else if(count == 1){
		slope = offset = NAN_D;
		return MathErrorUndefined;
	}
	
	double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;
	for(long i=0; i<count; ++i){
		sumX  += pointsX[i];
		sumY  += pointsY[i];
		sumXY += pointsX[i]*pointsY[i];
		sumX2 += pointsX[i]*pointsX[i];
	}
	slope = (count*sumXY - sumX*sumY)/(count*sumX2 - SQ(sumX));
	offset = sumY/count - slope*sumX/count;
	return MathErrorNone;
}



MathError fitLinearRegressionNANSensitive(double pointsX[], double pointsY[], long count, double &slope, double &offset){
	double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;
	long counted = 0;
	for(long i=0; i<count; ++i){
		if(std::isnan(pointsX[i]) || std::isnan(pointsY[i])) continue;
		++counted;
		sumX  += pointsX[i];
		sumY  += pointsY[i];
		sumXY += pointsX[i]*pointsY[i];
		sumX2 += pointsX[i]*pointsX[i];
	}
	if(counted == 0){
		slope = offset = NAN_D;
		return MathErrorNoData;
	}else if(counted == 1){
		slope = offset = NAN_D;
		return MathErrorUndefined;
	}
	slope = (count*sumXY - sumX*sumY)/(count*sumX2 - SQ(sumX));
	offset = sumY/count - slope*sumX/count;
	return MathErrorNone;
}



template<class TYPE_X, class TYPE_Y>
double interpolateBetween2PointsLinear(const TYPE_X &x1, const TYPE_Y &y1, const TYPE_X &x2, const TYPE_Y &y2, const TYPE_X &x){
	if(x1 == x2){
		return TYPE_Y(NAN_D);
	}else{
		return y1 + (x-x1)*(y2-y1)/(x2-x1);
	}
}




// interpolates an 'old' time series at new time points
template<class TYPE, class TIME_ARRAY_TYPE2>
bool interpolateTimeSeriesAtTimes(	const std::vector<double> 	&oldTimes,					// (INPUT)
									const std::vector<TYPE> 	&valuesAtOldTimes,			// (INPUT)
									long						oldStart,
									long						oldEnd,
									const TIME_ARRAY_TYPE2		&newTimes,					// (INPUT)
									long						newStart,
									long						newEnd,
									long						&includedNewTimesStart, 	// (OUTPUT)
									long						&includedNewTimesEnd,	 	// (OUTPUT)
									std::vector<TYPE>			&valuesAtNewTimes){		 	// (OUTPUT)
	valuesAtNewTimes.clear();
	if((newStart<0) || (newStart>newEnd) || (oldStart<0) || (oldStart>oldEnd)) return true;
	if((oldTimes[oldEnd]<newTimes[newStart]) || (newTimes[newEnd]<oldTimes[oldStart])) return true;

	// find first time point to interpolate at
	long tpNew = newStart;
	while((tpNew<=newEnd) && (newTimes[tpNew]<oldTimes[oldStart])) ++tpNew;
	if((tpNew>newEnd) || (newTimes[tpNew]>oldTimes[oldEnd])) return true; // new time domain does not overlap with old time domain
	includedNewTimesStart = tpNew;
	
	// find last time point to interpolate at (mainly for space pre-allocation purposes)
	long tpLastNew = newEnd;
	while((tpLastNew>=tpNew) && (newTimes[tpLastNew]>oldTimes[oldEnd])) --tpLastNew;
	if((tpLastNew<tpNew) || (newTimes[tpLastNew]<oldTimes[oldStart])) return true; // new time domain does not overlap with old time domain
	includedNewTimesEnd = tpLastNew;
	
	valuesAtNewTimes.reserve(includedNewTimesEnd - includedNewTimesStart + 1);
	
	long tpOld = oldStart;
	double time;
	for(; tpNew<=tpLastNew; ++tpNew){
		time = newTimes[tpNew];
		
		// at this point it is guaranteed that   oldTimes[tpOld]<=time   and that oldTimes[oldEnd]>=time
		
		//find left-closest oldTime to time
		while((tpOld+1<=oldEnd) && (oldTimes[tpOld+1]<=time)){ ++tpOld; }
		
		// at this point it is guaranteed that either (tpOld==oldEnd and oldTimes[tpOld]==time), or (tpOld<oldEnd and oldTimes[tpOld]<=time and oldTimes[tpOld+1]>time)
				
		// interpolate old series at focal time
		valuesAtNewTimes.push_back(tpOld==oldEnd ? valuesAtOldTimes[oldEnd] : interpolateBetween2PointsLinear(oldTimes[tpOld], valuesAtOldTimes[tpOld], oldTimes[tpOld+1], valuesAtOldTimes[tpOld+1], time));
	}
	return true;
}






// fits the curve y = A + B*x + C*x^2 to the data set defined by (x,y) using least squares
// returns false if an error occured, in which case the coefficients will be set to NAN_D.
// epsilon specifies the aimed-for fitting accuracy (in case of iterative refinement), relative to the mean{|y|}. epsilon=1e-6 would typically be reasonable.
// This function is NAN-sensitive, i.e. it avoids NAN in the x[] and y[] data.
bool fitLeastSquares_Quadratic(const std::vector<double> &x, const std::vector<double> &y, double &A, double &B, double &C, double epsilon){
	// calculate problem coefficients
	std::vector<double> b(3,0), sol(3), scratchSpace(3*3), matrix(3*3,0);
	const long N = x.size();
	
	// get scales
	double scaleX = 0, scaleY = 0;
	long counted = 0;
	for(long i=0; i<N; ++i){
		if(std::isnan(x[i]) || std::isnan(y[i])) continue;
		++counted;
		scaleX += abs(x[i]);
		scaleY += abs(y[i]);
	}
	if(counted<3){ A=B=C=NAN_D; return false; }
	scaleX = (counted==0 ? 1.0 : scaleX/counted);
	scaleY = (counted==0 ? 1.0 : scaleY/counted);
	
	for(long i=0; i<N; ++i){
		if(std::isnan(x[i]) || std::isnan(y[i])) continue;
		b[0] += y[i]/scaleY;
		b[1] += (x[i]/scaleX) * (y[i]/scaleY);
		b[2] += SQ(x[i]/scaleX) * (y[i]/scaleY);

		// Vandermonde matrix
		matrix[0*3+0] += 1;
		matrix[0*3+1] += x[i]/scaleX;
		matrix[0*3+2] += SQ(x[i]/scaleX);
		
		matrix[1*3+0] += x[i]/scaleX;
		matrix[1*3+1] += SQ(x[i]/scaleX);
		matrix[1*3+2] += Qube(x[i]/scaleX);
		
		matrix[2*3+0] += SQ(x[i]/scaleX);
		matrix[2*3+1] += Qube(x[i]/scaleX);
		matrix[2*3+2] += QTR(x[i]/scaleX);
	}
	
	if(!solveLinearSystem(	&matrix[0],
							&scratchSpace[0],
							3,
							&b[0],
							epsilon*3,
							5,
							&sol[0])){
		A=B=C=NAN_D;
		return false;
	}
	A = scaleY*sol[0]; 
	B = scaleY*sol[1]/scaleX; 
	C = scaleY*sol[2]/SQ(scaleX); 
	return true;
}






// fits the curve y = A + B*x + C*x^2 + D*x^3 to the data set defined by (x,y) using least squares
// returns false if an error occured, in which case the coefficients will be set to NAN_D.
// epsilon specifies the aimed-for fitting accuracy (in case of iterative refinement), relative to the mean{|y|}. epsilon=1e-6 would typically be reasonable.
// This function is NAN-sensitive, i.e. it avoids NAN in the x[] and y[] data
bool fitLeastSquares_Qubic(const std::vector<double> &x, const std::vector<double> &y, double &A, double &B, double &C, double &D, double epsilon){
	// calculate problem coefficients
	std::vector<double> b(4,0), sol(4), scratchSpace(4*4), matrix(4*4,0);
	const long N = x.size();
	
	// get scales
	double scaleX = 0, scaleY = 0;
	long counted = 0;
	for(long i=0; i<N; ++i){
		if(std::isnan(x[i]) || std::isnan(y[i])) continue;
		++counted;
		scaleX += abs(x[i]);
		scaleY += abs(y[i]);
	}
	if(counted<4){ A=B=C=D=NAN_D; return false; }
	scaleX = (counted==0 ? 1.0 : scaleX/counted);
	scaleY = (counted==0 ? 1.0 : scaleY/counted);
	
	for(long i=0; i<N; ++i){
		if(std::isnan(x[i]) || std::isnan(y[i])) continue;
		b[0] += y[i]/scaleY;
		b[1] += (x[i]/scaleX) * (y[i]/scaleY);
		b[2] += SQ(x[i]/scaleX) * (y[i]/scaleY);
		b[3] += Qube(x[i]/scaleX) * (y[i]/scaleY);
		
		// Vandermonde matrix
		matrix[0*4+0] += 1;
		matrix[0*4+1] += x[i]/scaleX;
		matrix[0*4+2] += SQ(x[i]/scaleX);
		matrix[0*4+3] += Qube(x[i]/scaleX);
		
		matrix[1*4+0] += x[i]/scaleX;
		matrix[1*4+1] += SQ(x[i]/scaleX);
		matrix[1*4+2] += Qube(x[i]/scaleX);
		matrix[1*4+3] += QTR(x[i]/scaleX);
		
		matrix[2*4+0] += SQ(x[i]/scaleX);
		matrix[2*4+1] += Qube(x[i]/scaleX);
		matrix[2*4+2] += QTR(x[i]/scaleX);
		matrix[2*4+3] += (x[i]/scaleX)*QTR(x[i]/scaleX);
		
		matrix[3*4+0] += Qube(x[i]/scaleX);
		matrix[3*4+1] += QTR(x[i]/scaleX);
		matrix[3*4+2] += (x[i]/scaleX)*QTR(x[i]/scaleX);
		matrix[3*4+3] += SQ(x[i]/scaleX)*QTR(x[i]/scaleX);
	}
	
	if(!solveLinearSystem(	&matrix[0],
							&scratchSpace[0],
							4,
							&b[0],
							epsilon*4,
							5,
							&sol[0])){
		A=B=C=D=NAN_D;
		return false;
	}
	A = scaleY*sol[0]; 
	B = scaleY*sol[1]/scaleX; 
	C = scaleY*sol[2]/SQ(scaleX); 
	D = scaleY*sol[3]/Qube(scaleX);
	return true;
}	







// fits the curve y = A + B*x + C*x^2 + D*x^3 + E*x^4 to the data set defined by (x,y) using least squares
// returns false if an error occured, in which case the coefficients will be set to NAN_D.
// epsilon specifies the aimed-for fitting accuracy (in case of iterative refinement), relative to the mean{|y|}. epsilon=1e-6 would typically be reasonable.
// This function is NAN-sensitive, i.e. it avoids NAN in the x[] and y[] data
bool fitLeastSquares_Quartic(const std::vector<double> &x, const std::vector<double> &y, double &A, double &B, double &C, double &D, double &E, double epsilon){
	// calculate problem coefficients
	std::vector<double> b(5,0), sol(5), scratchSpace(5*5), matrix(5*5,0);
	const long N = x.size();
	
	// get scales
	double scaleX = 0, scaleY = 0;
	long counted = 0;
	for(long i=0; i<N; ++i){
		if(std::isnan(x[i]) || std::isnan(y[i])) continue;
		++counted;
		scaleX += abs(x[i]);
		scaleY += abs(y[i]);
	}
	if(counted<5){ A=B=C=D=E=NAN_D; return false; }
	scaleX = (counted==0 ? 1.0 : scaleX/counted);
	scaleY = (counted==0 ? 1.0 : scaleY/counted);
	
	for(long i=0; i<N; ++i){
		if(std::isnan(x[i]) || std::isnan(y[i])) continue;
		b[0] += y[i]/scaleY;
		b[1] += (x[i]/scaleX) * (y[i]/scaleY);
		b[2] += SQ(x[i]/scaleX) * (y[i]/scaleY);
		b[3] += Qube(x[i]/scaleX) * (y[i]/scaleY);
		b[4] += QTR(x[i]/scaleX) * (y[i]/scaleY);
		
		// Vandermonde matrix
		matrix[0*5+0] += 1;
		matrix[0*5+1] += x[i]/scaleX;
		matrix[0*5+2] += SQ(x[i]/scaleX);
		matrix[0*5+3] += Qube(x[i]/scaleX);
		matrix[0*5+4] += QTR(x[i]/scaleX);
		
		matrix[1*5+0] += x[i]/scaleX;
		matrix[1*5+1] += SQ(x[i]/scaleX);
		matrix[1*5+2] += Qube(x[i]/scaleX);
		matrix[1*5+3] += QTR(x[i]/scaleX);
		matrix[1*5+4] += (x[i]/scaleX)*QTR(x[i]/scaleX);
		
		matrix[2*5+0] += SQ(x[i]/scaleX);
		matrix[2*5+1] += Qube(x[i]/scaleX);
		matrix[2*5+2] += QTR(x[i]/scaleX);
		matrix[2*5+3] += (x[i]/scaleX)*QTR(x[i]/scaleX);
		matrix[2*5+4] += SQ(x[i]/scaleX)*QTR(x[i]/scaleX);
		
		matrix[3*5+0] += Qube(x[i]/scaleX);
		matrix[3*5+1] += QTR(x[i]/scaleX);
		matrix[3*5+2] += (x[i]/scaleX)*QTR(x[i]/scaleX);
		matrix[3*5+3] += SQ(x[i]/scaleX)*QTR(x[i]/scaleX);
		matrix[3*5+4] += Qube(x[i]/scaleX)*QTR(x[i]/scaleX);
		
		matrix[4*5+0] += QTR(x[i]/scaleX);
		matrix[4*5+1] += (x[i]/scaleX)*QTR(x[i]/scaleX);
		matrix[4*5+2] += SQ(x[i]/scaleX)*QTR(x[i]/scaleX);
		matrix[4*5+3] += Qube(x[i]/scaleX)*QTR(x[i]/scaleX);
		matrix[4*5+4] += QTR(x[i]/scaleX)*QTR(x[i]/scaleX);
	}
	
	if(!solveLinearSystem(	&matrix[0],
							&scratchSpace[0],
							5,
							&b[0],
							epsilon*5,
							5,
							&sol[0])){
		A=B=C=D=E=NAN_D;
		return false;
	}
	A = scaleY*sol[0]; 
	B = scaleY*sol[1]/scaleX; 
	C = scaleY*sol[2]/SQ(scaleX); 
	D = scaleY*sol[3]/Qube(scaleX);
	E = scaleY*sol[4]/QTR(scaleX);
	return true;
}




// Performs a least residual sum of squares fitting of the real scalar linear function y(x) = slope*x to the data (x[], y[]), by choice of slope.
// Here, x is the predictor and y the response variable.
// returns false on error
template<class REAL_TYPE>
bool fitLeastSquares_linear_real_scalar(const vector<REAL_TYPE> 	&x,		// (input)	Measured covariates
										const vector<REAL_TYPE> 	&y,		// (input)	Measured response variables
										unsigned long				start,	// (input) first considered data index (ignore any entries prior)
										unsigned long				end,	// (input) last considered data index (ignore any entries after)
										REAL_TYPE 					&slope,	// (output) Fitted slope
										REAL_TYPE					&RSS){	// (output) The sum of squared residuals
	REAL_TYPE sumX2(0), sumXY(0), sumY2(0);
	if((x.size()<=end) || (y.size()<=end)) return false;
	if(start>end) return false; //need at least one data point
	unsigned long k;
	for(k=start; k<=end; ++k){
		sumX2 += SQ(x[k]);
		sumY2 += SQ(y[k]);
		sumXY += x[k]*y[k];
	}
	slope 	= sumXY/sumX2;
	RSS		= sumY2 - SQ(sumXY)/sumX2;
	return true;
}




template<class REAL_TYPE>
REAL_TYPE affine_real_scalar_RSS(REAL_TYPE slope, REAL_TYPE intercept, long count, REAL_TYPE meanX, REAL_TYPE meanY, REAL_TYPE meanX2, REAL_TYPE meanY2, REAL_TYPE sumXY){
	return count * SQ(intercept)
			+ 2*slope*intercept*meanX*count 
			- 2*intercept*meanY*count 
			+ SQ(slope)*meanX2*count 
			+ meanY2*count 
			- 2*slope*sumXY;
}


// Performs a least residual sum of squares fitting of the real scalar affine function y(x) = intercept + slope*x to the data (x[], y[]), by choice of intercept & slope.
// Here, x is the covariate and y the response variable.
// returns false on error
template<class REAL_TYPE>
bool fitLeastSquares_affine_real_scalar(const vector<REAL_TYPE> 	&x,			// (input) Measured covariates
										const vector<REAL_TYPE> 	&y,			// (input) Measured response variables
										unsigned long				start,		// (input) first considered data index (ignore any entries prior)
										unsigned long				end,		// (input) last considered data index (ignore any entries after)
										bool						forcePositiveSlope,		// (input)
										bool						forcePositiveIntercept,	// (input)
										REAL_TYPE					&intercept, // (output) Fitted intercept
										REAL_TYPE 					&slope,		// (output) Fitted slope
										REAL_TYPE					&RSS){		// (output) The sum of squared residuals
	REAL_TYPE meanX(0), meanY(0), meanX2(0), meanY2(0), sumXY(0);
	if((x.size()<=end) || (y.size()<=end)) return false;
	if(start>=end) return false; //need at least two data points
	long k, count=end-start+1;
	
	for(k=start; k<=end; ++k){
		meanX	+= x[k];
		meanY	+= y[k];
		meanX2 	+= SQ(x[k]);
		meanY2 	+= SQ(y[k]);
		sumXY 	+= x[k]*y[k];
	}
	meanX	/= count;
	meanY	/= count;
	meanX2	/= count;
	meanY2	/= count;
	
	slope 		= 	(sumXY/count - meanX*meanY)/(meanX2 - SQ(meanX));
	intercept	= 	(meanX2*meanY - meanX*sumXY/count)/(meanX2 - SQ(meanX));
	
	if(forcePositiveSlope && (!forcePositiveIntercept)){
		if(slope<0){
			slope = 0;
			intercept = meanY;
		}
	}else if((!forcePositiveSlope) && forcePositiveIntercept){
		if(intercept<0){
			intercept = 0;
			slope = sumXY/(count*meanX2);
		}
	}else if(forcePositiveSlope && forcePositiveIntercept){
		if(slope<0 || intercept<0){
			REAL_TYPE RSS00 = affine_real_scalar_RSS(0.0, 0.0, count, meanX, meanY, meanX2, meanY2, sumXY);
			REAL_TYPE RSS01 = affine_real_scalar_RSS(0.0, meanY, count, meanX, meanY, meanX2, meanY2, sumXY);
			REAL_TYPE RSS10 = affine_real_scalar_RSS(sumXY/(count*meanX2), 0.0, count, meanX, meanY, meanX2, meanY2, sumXY);
			if(meanY<0){ 
				intercept = 0;
				slope = ((RSS00<RSS10) || (sumXY<0) ? 0 : sumXY/(count*meanX2));
			}else{
				intercept 	= ((RSS01<RSS10) || (sumXY<0) ? meanY 	: 0);
				slope		= ((RSS01<RSS10) || (sumXY<0) ? 0 		: sumXY/(count*meanX2));
			}
		}
	}
	RSS = affine_real_scalar_RSS(slope, intercept, count, meanX, meanY, meanX2, meanY2, sumXY);

	if(isnan(slope) || isnan(intercept) || isnan(RSS)) return false;
	return true;
}




// Fits the real scalar exponential function y(x) = intercept*exp(x/scale) to the data (x[], y[]), by choice of scale and (possibly) intercept.
// Here, x is the covariate and y the response variable.
// Fitting minimizes sum of squared residuals of log-transformed responses to log-transformed predictions.
// returns false on error
template<class ARRAY_TYPE1, class ARRAY_TYPE2>
bool fitLeastLogSquares_exponential_real_scalar(const ARRAY_TYPE1	&x,					// (input) Measured covariates
												const ARRAY_TYPE2 	&y,					// (input) Measured response variables.
												long				start,				// (input) first considered data index (ignore any entries prior)
												long				end,				// (input) last considered data index (ignore any entries after)
												bool				intercept_fixed,	// (input) if true, intercept is fixed to the value given. Otherwise, intercept is fitted as well.
												double				&intercept, 		// (input/output) Fitted y-intercept if intercept_fixed==false, otherwise left unchanged.
												double 				&scale,				// (output) Fitted caracteristic scale (=inverse rate)
												double				&RSS){				// (output) The sum of squared residuals (on y-logarithmic scale)
	long k;
	if((end<=start) || (end>=x.size()) || (end>=y.size()) || (end<0) || (start<0)) return false;
	vector<double> X = vector<double>(x.begin()+start,x.begin()+end+1);
	vector<double> Y = vector<double>(y.begin()+start,y.begin()+end+1);
	start = 0; end = X.size()-1;
	for(k=start; k<=end; ++k){
		Y[k] = log(Y[k]);
	}
	if(intercept_fixed){
		double logintercept = log(intercept);
		for(k=start; k<=end; ++k){
			Y[k] -= logintercept;
		}
		if(!fitLeastSquares_linear_real_scalar(X, Y, start, end, scale, RSS)){
			return false;
		}
	}else{
		if(!fitLeastSquares_affine_real_scalar(X, Y, start, end, false, false, intercept, scale, RSS)){
			return false;
		}
		intercept = exp(intercept);
	}
	scale = 1.0/scale;
	return true;
}



// smoothen time series using Savitzky-Golay filtering (unweighted local least-squares polynomial fitting)
// for details see: http://www.mathworks.com/help/curvefit/smoothing-data.html
// Note: If the sliding window does not cover enough time points, a lower-order filter may be applied (locally)
// This function supports non-evenly spaced data, which means fitting a polynomial at each time point, making it computationally expensive
// The time series is assumed to be sorted in time
// This function is NAN-sensitive, i.e. NANs in data[] are avoided.
template<class TIME_ARRAY, class VALUE_ARRAY>
bool smoothenTimeSeriesSavitzkyGolay(	const TIME_ARRAY	&times,				// time points in ascending order
										const VALUE_ARRAY	&data, 	
										double				windowTimeSpan,		// span of fitting window in time units. Ignored if <=0.
										long				windowIndexSpan,	// span of fitting window in terms of the number of data points included (will be rounded up to nearest odd number). Can be used as an alternative to windowTimeSpan. Ignored if <=1.
										int					order,			
										vector<double> 		&smoothenedData){
	if(times.size() != data.size()) return false;
	if((order<0) || (order>4)) return false;
	smoothenedData.clear();
	const long N = times.size();
	if(N<1) return true;	// nothing to do
	double coeff_A, coeff_B, coeff_C, coeff_D, coeff_E;
	long localOrder, distinctTimePoints;
	vector<double> localTimes, localData;
	if(windowTimeSpan>0) windowTimeSpan = max(0.0, min(windowTimeSpan, times[N-1]-times[0]));
	else if(windowIndexSpan>1) windowIndexSpan = min(windowIndexSpan,N-1);

	smoothenedData.resize(N);
	for(long n=0, nl,nr,M; n<N; ++n){
		const double t = times[n];
		// determine data points to include
		if(windowIndexSpan>1){
			nl = n-windowIndexSpan/2;
			nr = n+windowIndexSpan/2;
			if(nl<0) nr += (0-nl); // shift to the right since we're at the left edge
			else if(nr>=N) nl -= (nr-N+1); // shift to the left since we're at the right edge
		}else{
			const double tl = t-windowTimeSpan/2;
			const double tr = t+windowTimeSpan/2;
			if(times[0]>=tl){
				nl=0;
			}else{
				for(nl=n; times[nl-1]>=tl; --nl){}
			}
			if(times[N-1]<=tr){
				nr=N-1;
			}else{
				for(nr=n; times[nr+1]<=tr; ++nr){}
			}
		}
		nl = max(0l,nl);
		nr = min((long)(N-1),nr);
		M  = nr-nl+1; // number of local data points covered by window
		
		// count distinct time points covered by window
		distinctTimePoints = 1;
		for(long i=nl+1; i<=nr; ++i){ if(times[i]>times[i-1]) ++distinctTimePoints; }
				
		if(M==1){
			smoothenedData[n] = data[n]; // only one data point. 
		}else if((distinctTimePoints==1) || (order==0)){
			// all local data are at the same time point (or polynomial order==0), so take arithmetic average
			double sum = 0; 
			long counted = 0;
			for(long i=nl; i<=nr; ++i){ if(!isnan(data[i])){ ++counted; sum += data[i]; } }
			smoothenedData[n] = (counted>0 ? sum/counted : NAN_D);
		
		}else{			
			// prepare local subset of data for polynomial fitting
			localTimes.resize(M);
			localData.resize(M);
			for(long i=nl; i<=nr; ++i){
				localTimes[i-nl] = times[i];
				localData[i-nl]	 = data[i];
			}
		
			// fit polynomial using least squares and evaluate at centre point
			localOrder = (int)min(distinctTimePoints-1,(long)order); // if too few data points available, use lower order
			switch(localOrder){
			case 1: 
				fitLinearRegressionNANSensitive(&localTimes[0], &localData[0], M, coeff_B, coeff_A); 
				smoothenedData[n] = coeff_A + t*coeff_B;
				break;			
			case 2:
				fitLeastSquares_Quadratic(localTimes, localData, coeff_A, coeff_B, coeff_C, 1e-6);	
				smoothenedData[n] = coeff_A + t*coeff_B + SQ(t)*coeff_C;
				break;			
			case 3:
				fitLeastSquares_Qubic(localTimes, localData, coeff_A, coeff_B, coeff_C, coeff_D, 1e-6);	
				smoothenedData[n] = coeff_A + t*coeff_B + SQ(t)*coeff_C + Qube(t)*coeff_D;
				break;	
			case 4:
				fitLeastSquares_Quartic(localTimes, localData, coeff_A, coeff_B, coeff_C, coeff_D, coeff_E, 1e-6);	
				smoothenedData[n] = coeff_A + t*coeff_B + SQ(t)*coeff_C + Qube(t)*coeff_D + QTR(t)*coeff_E;
				break;	
			}
		}
	}
	return true;												
}


// Rcpp wrapper for the homonymous function above
// [[Rcpp::export]]
Rcpp::List smoothenTimeSeriesSavitzkyGolay_CPP(	const NumericVector	&times,				// time points in ascending order
												const NumericVector	&data, 	
												double				windowTimeSpan,		// span of fitting window in time units. Ignored if <=0.
												long				windowIndexSpan,	// span of fitting window in terms of the number of data points included (will be rounded up to nearest odd number). Can be used as an alternative to windowTimeSpan. Ignored if <=1.
												int					order){
	std::vector<double> smoothened_data;
	const bool success = smoothenTimeSeriesSavitzkyGolay(times, data, windowTimeSpan, windowIndexSpan, order, smoothened_data);
	return Rcpp::List::create(	Rcpp::Named("success") 			= success,
								Rcpp::Named("smoothened_data") 	= smoothened_data);
}



// Use Savitzky-Golay filter to estimate nth-derivative of a time series
// Elaboration: At each data point, SG filter fits a polynomial to local (nearby) data points using least squares
//				The coefficients of that polynomial are then used to estimate the derivative at that points
//	Note: 	This method DOES NOT ensure that the antiderivative of the estimated derivative (i.e. its integral) equals the original signal or the SG-smoothened signal
//			If you need consistency between the integrated derivative and the smoothened signal, consider SG-smoothening the signal and then taking the derivative using some finite differences scheme
// Currently only 1st and 2nd derivatives are supported (but it shouldn't be hard to extend this function to higher orders)
// Only fitting polynomials of orders 2-4 are supported
// Time points are assumed to be sorted in ascending order
// If the sliding window contains too few data points, the order of the fitted polynomials may be locally reduced
// Note: Time points for which the derivative cannot be estimated (e.g. too few time points available in sliding window) will have derivate==NAN_D
// This function is NAN-sensitive, i.e. NANs in data[] are avoided.
// The set of time points included for each fitting is either specified via windowTimeSpan or windowIndexSpan
template<class TIME_ARRAY, class VALUE_ARRAY>
bool derivativeOfTimeSeries_SavitzkyGolay(	const TIME_ARRAY	&times,				// time points in ascending order
											const VALUE_ARRAY	&data, 	
											double				windowTimeSpan,		// span of fitting window in time units. Ignored if <=0.
											long				windowIndexSpan,	// span of fitting window in terms of the number of data points included (will be rounded up to nearest odd number). Can be used as an alternative to windowTimeSpan. Ignored if <=1.
											int					polynomialOrder,
											int					derivativeOrder,
											vector<double> 		&derivative){
	if(times.size() != data.size()) return false;
	if((polynomialOrder<0) || (polynomialOrder>4) || (derivativeOrder<1) || (derivativeOrder>2) || (polynomialOrder<=derivativeOrder)) return false;
	derivative.clear();
	const long N = times.size();
	if(N<1) return true;	// nothing to do
	double coeff_A, coeff_B, coeff_C, coeff_D, coeff_E;
	long M, localPolynomialOrder, distinctTimePoints;
	vector<double> localTimes, localData;
	if(windowTimeSpan>0) windowTimeSpan = max(0.0, min(windowTimeSpan, times[N-1]-times[0]));
	else if(windowIndexSpan>1) windowIndexSpan = min(windowIndexSpan,N-1);

	derivative.resize(N);
	for(long n=0, nl,nr; n<N; ++n){
		const double t = times[n];
		// determine data points to include
		if(windowIndexSpan>1){
			nl = n-windowIndexSpan/2;
			nr = n+windowIndexSpan/2;
			if(nl<0) nr += (0-nl); // shift to the right since we're at the left edge
			else if(nr>=N) nl -= (nr-N+1); // shift to the left since we're at the right edge
		}else{
			const double tl = t-windowTimeSpan/2;
			const double tr = t+windowTimeSpan/2;
			if(times[0]>=tl){
				nl=0;
			}else{
				for(nl=n; times[nl-1]>=tl; --nl){}
			}
			if(times[N-1]<=tr){
				nr=N-1;
			}else{
				for(nr=n; times[nr+1]<=tr; ++nr){}
			}
		}
		nl = max(0l,nl);
		nr = min((long)(N-1),nr);
		M  = nr-nl+1; // number of local data points covered by window
						
		if(M<derivativeOrder+2){
			// definitely too few data points for estimating derivative 
			derivative[n] = NAN_D;
		
		}else{
			// count distinct time points covered by window
			distinctTimePoints = 1;
			for(long i=nl+1; i<=nr; ++i){ if(times[i]>times[i-1]) ++distinctTimePoints; }
			
			// use lower polynomial order if too few distinct time points
			localPolynomialOrder = (int)min(distinctTimePoints-1,(long)polynomialOrder);
			
			if(localPolynomialOrder<derivativeOrder+1){
				// polynomial order is too low for this derivative order
				derivative[n] = NAN_D;
			}else{
				// prepare local subset of data for polynomial fitting
				localTimes.resize(M);
				localData.resize(M);
				for(long i=nl; i<=nr; ++i){
					localTimes[i-nl] = times[i];
					localData[i-nl]	 = data[i];
				}

				// fit polynomial
				coeff_A = coeff_B = coeff_C = coeff_D = coeff_E = 0;
				switch(localPolynomialOrder){
				case 1: 
					fitLinearRegressionNANSensitive(&localTimes[0], &localData[0], M, coeff_B, coeff_A); 
					break;			
				case 2:
					fitLeastSquares_Quadratic(localTimes, localData, coeff_A, coeff_B, coeff_C, 1e-6);	
					break;			
				case 3:
					fitLeastSquares_Qubic(localTimes, localData, coeff_A, coeff_B, coeff_C, coeff_D, 1e-6);	
					break;	
				case 4:
					fitLeastSquares_Quartic(localTimes, localData, coeff_A, coeff_B, coeff_C, coeff_D, coeff_E, 1e-6);	
					break;	
				}
				
				// calculate derivative using polynomial coefficients
				switch(derivativeOrder){
				case 1: derivative[n] = coeff_B + 2*coeff_C*t + 3*coeff_D*SQ(t) + 4*coeff_E*Qube(t); break;
				case 2: derivative[n] = 2*coeff_C + 6*coeff_D*t + 12*coeff_E*SQ(t); break;
				}
			}
		}
	}
	return true;														
}









#pragma mark -
#pragma mark Class: LinearInterpolationFunctor
#pragma mark -


//Class for defining piecewiese linear function using a reference grid
template<class VALUE_TYPE>
class LinearInterpolationFunctor{
private:
	std::vector<double> 	referencePoints;
	std::vector<VALUE_TYPE>	referenceValues;
	
	double		domain_min, domain_max;
	double		domainStep;
	double		lengthScale;
	
	bool		periodic;
	VALUE_TYPE	out_of_scope_value;
	bool		is_zero, is_constant;
	
	void 	set_to_regular_grid_values(	long 				referenceCount, 
										const double		domain_min, 
										const double		domain_max, 
										VALUE_TYPE 			referenceValues[],
										bool				periodic,
										const VALUE_TYPE	&out_of_scope_value);		// value to use for extending time series if needed, if periodic==false. Irrelevant if periodic==true.
public:

	LinearInterpolationFunctor(): domain_min(0), domain_max(0), domainStep(0), lengthScale(1), periodic(false), is_zero(true), is_constant(true) {}
	
	//define piecewise linear function
	//domain is assumed to span from lowest to highest referencePoints
	//if coordinate is periodic, the last reference point is identified with the first one (its reference value is taken to be (referenceValues[0]+referenceValues[last])/2)
	LinearInterpolationFunctor(	const std::vector<double> 		&referencePoints, // domain grid points, in ascending order
								const std::vector<VALUE_TYPE> 	&referenceValues,
								bool							periodic,	// if true, time series are extended periodically if needed (e.g. if the reference domain is too small). Otherwise they are extended with zeros.
								const VALUE_TYPE				&out_of_scope_value,
								bool							preInterpolateOnRegularGrid); // if true, then the function is internally pre-interpolated on a regular domain grid. In that case, value retrieval later on will be more efficient.
	
	//define piecewise linear function
	//referenceValues[i] should correspond to domain_min + i * intervalWidth(domain)/(referenceCount-1)
	//in particular, referenceValue[referenceCount-1] should correspond to domain_max
	//if coordinate is periodic, the last reference point is identified with the first one (its reference value is taken to be (referenceValues[0]+referenceValues[last])/2)
	//since the domain grid is assumed to be regular, interpolations defined this way are much for efficient
	LinearInterpolationFunctor(	long 				referenceCount, 
								const double		domain_min, 
								const double		domain_max, 
								VALUE_TYPE 			referenceValues[],
								bool				periodic,
								const VALUE_TYPE	&out_of_scope_value);
							
	VALUE_TYPE operator()(double x) const;
	
	bool isNULL() const{ return referenceValues.empty(); }
	bool isPeriodic() const{ return periodic; }
	void getDomain(double &_domain_min, double &_domain_max) const{ _domain_min=domain_min; _domain_max=domain_max; }
	long getReferenceCount() const{ return referenceValues.size(); }
	bool isConstant() const{ return is_constant; }
	bool isZero() const{ return is_zero; }
	
	void setTypicalLengthScale(double _lengthScale){ lengthScale = _lengthScale; } 
	double typicalLengthScale() const{ return lengthScale; }
	
	const std::vector<double> &getReferencePoints() const{ return referencePoints; }
	const std::vector<double> &getReferenceValues() const{ return referenceValues; }
};





template<class VALUE_TYPE>
LinearInterpolationFunctor<VALUE_TYPE>::LinearInterpolationFunctor(	const std::vector<double> 		&_referencePoints, 
																	const std::vector<VALUE_TYPE> 	&_referenceValues,
																	bool							_periodic,
																	const VALUE_TYPE				&_out_of_scope_value,
																	bool							preInterpolateOnRegularGrid){
	periodic 			= _periodic;
	is_zero 			= true;
	is_constant 		= true;
	out_of_scope_value 	= _out_of_scope_value;
	referencePoints.clear(); 
	referenceValues.clear();
	if(_referencePoints.empty()) return;

	if(preInterpolateOnRegularGrid){
		// pre-interpolate on regular domain grid, then setup interpolator functor
		const double step 			= smallest_nonzero_time_step(_referencePoints);
		const long NR 				= 1 + (_referencePoints.back() - _referencePoints.front())/step;
		const double _domain_min 	= _referencePoints.front();
		std::vector<double> regular_reference_points(NR), regular_reference_values;
		for(long i=0; i<NR; ++i) regular_reference_points[i] = _domain_min + i*step;
		long includedNewTimesStart, includedNewTimesEnd;
		interpolateTimeSeriesAtTimes(	_referencePoints,
										_referenceValues,
										0,
										_referencePoints.size()-1,
										regular_reference_points,
										0,
										NR-1,
										includedNewTimesStart,
										includedNewTimesEnd,
										regular_reference_values);
		set_to_regular_grid_values(	1+includedNewTimesEnd-includedNewTimesStart, 
									regular_reference_points[includedNewTimesStart], 
									regular_reference_points[includedNewTimesEnd], 
									&regular_reference_values[includedNewTimesStart], 
									_periodic,
									_out_of_scope_value);		
		
	}else{
		// use provided irregular grid
		referencePoints = _referencePoints;
		referenceValues = _referenceValues;
		const long referenceCount = referencePoints.size();
	
		// check if all values are zero
		is_zero = true;
		for(long i=0; i<referencePoints.size(); ++i){
			if(referenceValues[i]!=0) is_zero = false;
		}
	
		// figure out domain
		domain_min	= referencePoints.front();
		domain_max	= referencePoints.back();
		lengthScale	= domain_max - domain_min;
	
		// check constantness
		is_constant = true;
		for(long i=1; i<referenceCount; ++i){
			if(referenceValues[i]!=referenceValues[i-1]){ is_constant = false; break; }
		}
	
		// force periodic boundaries if necessary
		if(periodic){
			double boundaryValue = (referenceValues[0] + referenceValues[referenceCount-1])/2.0;
			referenceValues[0] = referenceValues[referenceCount-1] = boundaryValue;
		}
	}
}


template<class VALUE_TYPE>
void LinearInterpolationFunctor<VALUE_TYPE>::set_to_regular_grid_values(long 				referenceCount, 
																		const double 		_domain_min, 
																		const double 		_domain_max, 
																		VALUE_TYPE 			*_referenceValues,
																		bool				_periodic,
																		const VALUE_TYPE	&_out_of_scope_value){

	periodic			= _periodic;
	domain_min			= _domain_min;
	domain_max			= _domain_max;
	lengthScale			= (domain_max-domain_min);
	domainStep			= (domain_max-domain_min)/max(1.0, referenceCount-1.0);
	is_zero 			= true;
	is_constant 		= true;
	out_of_scope_value 	= _out_of_scope_value;
	referencePoints.clear(); 
	referenceValues.clear();
	if(referenceCount == 0) return;
	
	referenceValues.resize(referenceCount);
	for(long i=0; i<referenceCount; ++i){
		referenceValues[i] = _referenceValues[i];
		if(referenceValues[i] != 0) is_zero = false;
	}
	
	// check constantness
	is_constant = true;
	for(long i=1; i<referenceCount; ++i){
		if(referenceValues[i]!=referenceValues[i-1]){ is_constant = false; break; }
	}
	
	// enforce periodic boundary values if necessary
	if(periodic){
		double boundaryValue 	= (referenceValues[0] + referenceValues[referenceCount-1])/2.0;
		referenceValues[0] 		= referenceValues[referenceCount-1] = boundaryValue;
	}															
}


template<class VALUE_TYPE>
LinearInterpolationFunctor<VALUE_TYPE>::LinearInterpolationFunctor(	long 				referenceCount, 
																	const double 		_domain_min, 
																	const double 		_domain_max, 
																	VALUE_TYPE 			_referenceValues[],
																	bool				_periodic,
																	const VALUE_TYPE	&_out_of_scope_value){
	set_to_regular_grid_values(referenceCount, _domain_min, _domain_max, _referenceValues, _periodic, _out_of_scope_value);
}




template<class VALUE_TYPE>
VALUE_TYPE LinearInterpolationFunctor<VALUE_TYPE>::operator()(double x) const{
	if(referenceValues.empty()){
		//nothing available
		return 0;
	}
	if(is_zero) return 0;
	const long referenceCount = referenceValues.size();
	if(periodic){
		x = moduloInterval(x, domain_min, domain_max);
	}else if((x <= domain_min) || ( x >= domain_max)){
		// requested x is outside of the reference domain
		return out_of_scope_value;
	}
	if(referenceCount == 1){
		return referenceValues[0];
	}
	
	long j;
	if(referencePoints.empty()){
		j = floor((x - domain_min)/domainStep); //left border of relevant domain-interval
		j = min(j, referenceCount - 1); //avoid out-of-bounds errors caused by numerical inaccuracies
		if(j == referenceCount - 1){
			return referenceValues[referenceCount - 1];
		}else{
			return referenceValues[j] + (x - (domain_min + j*domainStep)) * (referenceValues[j+1] - referenceValues[j]) / domainStep; //linear interpolation between referenceValues
		}
	}else{
		for(j=1; j<referenceCount; ++j){
			if(referencePoints[j] > x){
				return referenceValues[j-1] + (x - referencePoints[j-1]) * (referenceValues[j] - referenceValues[j-1]) / (referencePoints[j] - referencePoints[j-1]);
			}
		}
		return referenceValues[referenceCount - 1];
	}
}







#pragma mark -
#pragma mark Deterministic diversity dynamics
#pragma mark



class ProgressReporter{
private:
	long reportCount;
	string prefix, suffix;
	bool asPercentage;
	mutable long lastReportedCase;
	bool silent;
public:
	ProgressReporter(const bool _silent){ reportCount = 10; lastReportedCase=INFTY_D; silent = _silent; }
	ProgressReporter(long _reportCount){ reportCount = _reportCount; lastReportedCase=INFTY_D; }
	ProgressReporter(long _reportCount, const string &_prefix, const string &_suffix, bool _asPercentage){ reportCount = _reportCount; prefix = _prefix; suffix = _suffix; asPercentage = _asPercentage; lastReportedCase = INFTY_D; }
	
	void setReportCount(long count){ reportCount = count; }
	void setPrefix(const string &_prefix){ prefix = _prefix; }
	void setSuffix(const string &_suffix){ suffix = _suffix; }
	void setAsPercentage(bool _asPercentage){ asPercentage = _asPercentage; }
	void reset(){ lastReportedCase = INFTY_D; }
	
	void operator()(long casesFinished, long totalCases, double exactFraction) const{
		if(reportCount<=0) return;
		if(silent) return;
		const double p = (1.0*casesFinished)/totalCases;
		const double last_p = (1.0*lastReportedCase)/totalCases;
		const double dp = 1.0/reportCount;
		if((lastReportedCase<=casesFinished) && (floor(p/dp)<=floor(last_p/dp))) return; //don't report this case
		if(floor(p/dp)==0) return;
		lastReportedCase = casesFinished;
		const long rounder = pow(10.0,1+log10(max(reportCount,1l)));
		Rcout << prefix;
		if(asPercentage){ Rcout << long(rounder*100.0*exactFraction)/rounder << " %"; }
		else{ Rcout << long(floor(exactFraction*totalCases)) << " out of " << totalCases; }
		Rcout << suffix;
	}
};



// Class for representing simulated diversity at any time point, including cumulative speciation & extinction events
class TreeStateHistory{
public:
	double diversity;		// current number of extant species
	double Nbirths;			// cumulative number of speciation events since the beginning of a simulation
	double Ndeaths;			// cumulative number of extinction events since the beginning of a simulation
	double Pextinction;		// extinction probability of a size-1 clade until max time. Only makes sense for simulations performed in reverse.
	TreeStateHistory(){ diversity = Nbirths = Ndeaths = Pextinction = 0; }
	TreeStateHistory(const double start_diversity){ diversity = start_diversity; Nbirths = Ndeaths = 0; Pextinction = 0; }
	bool isnan() const{ return (std::isnan(diversity) || std::isnan(Nbirths) || std::isnan(Ndeaths) || std::isnan(Pextinction)); }
};


inline TreeStateHistory operator*(TreeStateHistory x, double scalar){
	x.diversity 	*= scalar;
	x.Nbirths 		*= scalar;
	x.Ndeaths 		*= scalar;
	x.Pextinction 	*= scalar;
	return x;
}


TreeStateHistory& operator*=(TreeStateHistory &x, double scalar) {
	x.diversity 	*= scalar;
	x.Nbirths 	*= scalar;
	x.Ndeaths 	*= scalar;
	x.Pextinction	*= scalar;
	return x;
}


inline TreeStateHistory operator+(TreeStateHistory x, const TreeStateHistory &y){
	x.diversity 	+= y.diversity;
	x.Nbirths 		+= y.Nbirths;
	x.Ndeaths 		+= y.Ndeaths;
	x.Pextinction	+= y.Pextinction;
	return x;
}


inline TreeStateHistory &operator+=(TreeStateHistory &x, const TreeStateHistory &y){
	x.diversity 	+= y.diversity;
	x.Nbirths 	+= y.Nbirths;
	x.Ndeaths 	+= y.Ndeaths;
	x.Pextinction	+= y.Pextinction;
	return x;
}


inline TreeStateHistory operator-(TreeStateHistory x, const TreeStateHistory &y){
	x.diversity 	-= y.diversity;
	x.Nbirths 		-= y.Nbirths;
	x.Ndeaths 		-= y.Ndeaths;
	x.Pextinction	-= y.Pextinction;
	return x;
}


inline TreeStateHistory &operator-=(TreeStateHistory &x, const TreeStateHistory &y){
	x.diversity 	-= y.diversity;
	x.Nbirths 		-= y.Nbirths;
	x.Ndeaths 		-= y.Ndeaths;
	x.Pextinction	-= y.Pextinction;
	return x;
}


inline TreeStateHistory operator/(TreeStateHistory x, const TreeStateHistory &y){
	x.diversity 	/= y.diversity;
	x.Nbirths 		/= y.Nbirths;
	x.Ndeaths 		/= y.Ndeaths;
	x.Pextinction	/= y.Pextinction;
	return x;
}


inline TreeStateHistory operator/(TreeStateHistory x, double scalar){
	x.diversity 	/= scalar;
	x.Nbirths 		/= scalar;
	x.Ndeaths 		/= scalar;
	x.Pextinction	/= scalar;
	return x;
}


inline TreeStateHistory &operator/=(TreeStateHistory &x, double scalar){
	x.diversity 	/= scalar;
	x.Nbirths 		/= scalar;
	x.Ndeaths 		/= scalar;
	x.Pextinction	/= scalar;
	return x;
}



// class for species speciation-extinction model for tree generation
// includes integration of cumulative birth & death events
// can be integrated in reverse (time counted backwards) as well as in forward direction
// When integrated in reverse, the extinction probability P(t,T) is also simulated
class TreeSpeciationExtinctionModel{
	double 	min_valid_diversity;
	bool 	reverse;
	double 	reflection_time;
public:
	std::vector<TreeStateHistory> trajectory;
	std::vector<double> times;
	std::vector<double> survival_chances;
	double initial_diversity;
	double rarefaction;

	TreeSpeciationExtinctionModel(){ min_valid_diversity = 0; reverse = false; rarefaction = 1; }
	double birth_rate_intercept, birth_rate_factor, birth_rate_exponent;
	double death_rate_intercept, death_rate_factor, death_rate_exponent;
	long Nsplits;
		
	// model is warned that a time series of size count will need to be stored
	void reserveSpaceForTimeSeries(long Ntimes){ 
		trajectory.clear();
		trajectory.reserve(Ntimes); 
		times.clear();
		times.reserve(Ntimes); 
	}
	
	// reverse the dynamics of the model, i.e. return -V(_reflection_time-t) whenever V(t) is requested.
	void set_reverse(double _reflection_time){
		reflection_time = _reflection_time;
		reverse = true;
	}
	
	
	// calculate extinction probabilities, by integrating in reverse
	// this should be done after the end of a simulation, since the ODE's has initial condition and rate of change depends on the final trajectory of diversity
	// assumes that times[] are in ascending order
	void calculate_extinction_probabilities(){
		const long NMT = times.size();
		trajectory[NMT-1].Pextinction = (1.0-rarefaction); // probability of extinction begins at (1-rarefacton_fraction)
		for(long t=NMT-2; t>=0; --t){
			const double birth_rate_pc = get_speciation_rate_at_state(times[t+1],trajectory[t+1].diversity)/trajectory[t+1].diversity;
			const double death_rate_pc = get_extinction_rate_at_state(times[t+1],trajectory[t+1].diversity)/trajectory[t+1].diversity;
			const double dt = times[t+1]-times[t];
			trajectory[t].Pextinction = trajectory[t+1].Pextinction + dt*(death_rate_pc - trajectory[t+1].Pextinction*(birth_rate_pc+death_rate_pc) + pow(trajectory[t+1].Pextinction,Nsplits)*birth_rate_pc);
		}
	}
	
	
	// adjust (reduce) diversities time series to reflect a coalescent tree, based on the diversities & Pextinctions
	void make_coalescent(){
		if(!reverse) calculate_extinction_probabilities();  // if not integrated in reverse, then Pextinctions still need to be calculated
		for(long t=0; t<times.size(); ++t){
			trajectory[t].diversity *= (1-trajectory[t].Pextinction);
		}
	}

	bool getInitialState(double time, TreeStateHistory &state) const{ 
		state = TreeStateHistory(initial_diversity); 
		state.Pextinction = rarefaction;
		if(reverse) state.Pextinction = (1.0-rarefaction);
		else state.Pextinction = 0; // not implemented for forward integration; Pextinction can be calculated afterwards (in reverse) using calculate_extinction_probabilities()
		return true; 
	}

	// record time series point
	void registerState(double time, const TreeStateHistory &state){
		trajectory.push_back(state); 
		times.push_back(time); 
		if(reverse) trajectory.back().Pextinction = min(1.0,max(0.0,trajectory.back().Pextinction));
	}

	double get_speciation_rate_at_state(double time, const double diversity) const{
		return birth_rate_intercept + birth_rate_factor*pow(diversity,birth_rate_exponent);
	}
	
	double get_extinction_rate_at_state(double time, const double diversity) const{
		return death_rate_intercept + death_rate_factor*pow(diversity,death_rate_exponent);
	}

	RequestedDynamics getRateOfChangeAtState(double time, const TreeStateHistory &current_state, TreeStateHistory &rate_of_change, TreeStateHistory &jump_state){
		if(reverse) time = reflection_time-time;
		const double N = current_state.diversity;
		const double P = current_state.Pextinction;
		const double B = get_speciation_rate_at_state(time,current_state.diversity);
		const double D = get_extinction_rate_at_state(time,current_state.diversity);
		rate_of_change.Nbirths = B;
		rate_of_change.Ndeaths = D;
		rate_of_change.diversity 	= (Nsplits-1)*B - D;
		if(reverse) rate_of_change.Pextinction	= D/N - P*(B+D)/N + pow(P,Nsplits)*B/N;
		else rate_of_change.Pextinction = 0; // not implemented for forward integration
		if(reverse) rate_of_change.diversity = -rate_of_change.diversity;
		return RequestedDynamicsRateOfChange; 
	}

	// returns true if for some reason the time step should be refined, e.g. if the domain boundary is crossed
	bool checkShouldRefineTimeStep(double time, const TreeStateHistory &current_state, double dt, const TreeStateHistory &candidate_state) const{ 
		if(reverse){
			return (candidate_state.diversity<min_valid_diversity) || (candidate_state.Pextinction<0) || (candidate_state.Pextinction>1); 
		}else{
			return (candidate_state.diversity<min_valid_diversity);
		}
	}
	
	// check if candidateState is outside of the domain boundaries.
	// In that case, tries to correct the candidateState to be the "last" valid state on the linear trajectory
	// If this is not possible (e.g. because previousState was already on the boundary), the problematic components are brute-force adjusted to be within the domain. In this case, the routine returns CrossedBoundaryYesButFixedBruteForce.
	CrossedBoundary checkCrossedDomainBoundaryAndFix(	double					previousTime,
														const TreeStateHistory	&previousState,					// previous state (assumed to be valid!)
														double 					&dt,							// (INPUT/OUTPUT) will be adjusted (reduced) if candidateState crossed the boundary. The modified value is guaranteed to be within (0,dt]
														TreeStateHistory		&candidateState) const{			// (INPUT/OUTPUT) if candidateState is outside of domain boundaries, then this may become the "last" state (on the linear trajectory from previousState to candidateState) within the domain (if possible).
		if(candidateState.diversity>=min_valid_diversity){
			return CrossedBoundaryNo;
		}else if(previousState.diversity>min_valid_diversity){
			double lambda = 1;
			if(candidateState.diversity<min_valid_diversity) lambda = min(lambda,(previousState.diversity-min_valid_diversity)/(previousState.diversity-candidateState.diversity));
			candidateState = previousState*(1-lambda) + candidateState*lambda;
			dt *= lambda;
			return CrossedBoundaryYesButFixedByReducingTimeStep;
		}else{
			candidateState.diversity 	= max(min_valid_diversity,candidateState.diversity);
			candidateState.Pextinction 	= min(1.0,max(0.0,candidateState.Pextinction));
			return CrossedBoundaryYesButFixedBruteForce;
		}
		/*
		if((candidateState.diversity>=min_valid_diversity) && ((!reverse) || ((candidateState.Pextinction>=0) && (candidateState.Pextinction<=1)))){
			return CrossedBoundaryNo;
		}else if((previousState.diversity>min_valid_diversity) && ((!reverse) || ((previousState.Pextinction>=0) && (previousState.Pextinction<=1)))){
			double lambda = 1;
			if(candidateState.diversity<min_valid_diversity) lambda = min(lambda,(previousState.diversity-min_valid_diversity)/(previousState.diversity-candidateState.diversity));
			if(reverse && (candidateState.Pextinction<0)) lambda = min(lambda,(previousState.Pextinction-0)/(previousState.Pextinction-candidateState.Pextinction));
			if(reverse && (candidateState.Pextinction>1)) lambda = min(lambda,(1-previousState.Pextinction)/(candidateState.Pextinction-previousState.Pextinction));
			candidateState = previousState*(1-lambda) + candidateState*lambda;
			dt *= lambda;
			return CrossedBoundaryYesButFixedByReducingTimeStep;
		}else{
			candidateState.diversity 	= max(min_valid_diversity,candidateState.diversity);
			candidateState.Pextinction 	= min(1.0,max(0.0,candidateState.Pextinction));
			return CrossedBoundaryYesButFixedBruteForce;
		}
		*/
	}
	
	bool stateIsNaN(const TreeStateHistory &state) const{
		return state.isnan();
	}
	
	// estimate the maximum rate of change of diversity, for diversities within the interval [0, max_diversity]
	// may be useful for choosing the time step for simulations
	double estimate_max_rate_of_change(const double max_time, const double max_diversity, bool per_capita) const{
		double zenith_diversity;
		const double B3 = (per_capita ? birth_rate_exponent-1 : birth_rate_exponent); 	// adjust exponent for per-capita rates
		const double D3 = (per_capita ? death_rate_exponent-1 : death_rate_exponent);	// adjust exponent for per-capita rates
		if(death_rate_exponent==birth_rate_exponent) zenith_diversity = max_diversity;
		else if(((death_rate_factor==0) || (D3==0)) && (D3>B3)) zenith_diversity = max_diversity; // zenith is not defined
		else if(((birth_rate_factor==0) || (B3==0)) && (B3>D3)) zenith_diversity = max_diversity; // zenith is not defined
		else zenith_diversity = max(0.0,min(max_diversity, pow(((Nsplits-1)*birth_rate_factor*B3)/(death_rate_factor*D3), 1.0/(D3-B3))));
		const double R1 = abs((Nsplits-1)*get_speciation_rate_at_state(max_time,max_diversity) - get_extinction_rate_at_state(max_time,max_diversity))/(per_capita ? max_diversity : 1.0);
		const double R2 = abs((Nsplits-1)*get_speciation_rate_at_state(max_time,zenith_diversity) - get_extinction_rate_at_state(max_time,zenith_diversity))/(per_capita ? zenith_diversity : 1.0);
		const double R3 = abs((Nsplits-1)*get_speciation_rate_at_state(max_time,0) - get_extinction_rate_at_state(max_time,0));
		return max(R1, max(R2, R3));
	}
	
};



// RungeKutta2 integrator, handling:
//    1. Suggestions by the model for temporarily refined time steps
//    2. Situations where the domain boundary is crossed during an iteration (in which case the time step is temporarily decreased as needed in order to not overshoot)
//    3. Requests by the model to jump discontinuously to another state
// requested times might reverse (i.e. times at which model dynamics are requested might not be strictly increasing), but recorded time series will be strictly forward in time
template<class COORDINATE, class MODEL, class PROGRESS_REPORTER>
bool RungeKutta2(	double						startTime,						// (INPUT) simulation start time
					double 						endTime, 						// (INPUT) simulation end time
					double 						dt, 							// (INPUT) default integration time step
					MODEL 						&model,							// (INPUT/OUTPUT) object defining model dynamics, also handling time series storage
					long						maxRecordedPoints, 				// (INPUT) if small, some intermediate trajectory points are skipped (i.e. not recorded). This might be useful if accurate integration requires a small time step, but would produce to large of a time series.
					long						maxTimeStepRefinements,			// (INPUT) max allowed number of refinements of local time step when encountering invalid states. Only relevant if abortOnInvalidState==false.
					const PROGRESS_REPORTER		&reporter,						// (INPUT) callback functor to handle progress reporting during simulation
					const double				runtime_out_seconds,			// (INPUT) max allowed runtime in seconds. If <=0, this is ignored.
					string						&warningMessage){				// (OUTPUT) will be non-empty in case of error, or if non-fatal problems occurred
	COORDINATE currentPoint, candidatePoint, point2, jumpPoint;
	COORDINATE k1, k2, kConsensus;
	double t=startTime, t2, current_dt1, current_dt2;
	long recorded, refinements;
	RequestedDynamics dynamics;
	CrossedBoundary crossedBoundary;
	const double simulationTime = endTime-startTime; // target duration of simulation
	const double start_runtime  = get_thread_monotonic_walltime_seconds();
	
	// keep track of warnings
	bool crossedBoundaryButFixedByReducingTimeStep = false;
	bool crossedBoundaryYesButFixedBruteForce = false;
	bool crossedBoundaryYesButFixedUsingEuler = false;
	bool simulationExceededMaxProcessingTime = false;

	//preliminary error checking
	if(dt<simulationTime*RELATIVE_EPSILON){
		warningMessage = "Time step too small";
		return false;
	}
	if(simulationTime < dt){
		warningMessage = "Time step exceeds simulation time";
		return false;
	}
	if(maxRecordedPoints < 1){
		warningMessage = "Requested zero recorded points";
		return false;
	}
	
	const double recordingTimeStep = simulationTime/maxRecordedPoints;
	model.reserveSpaceForTimeSeries(maxRecordedPoints);
	bool k1AlreadyCalculated = false;
	
	//initialization
	if(!model.getInitialState(t, currentPoint)){
		warningMessage = "Failed to get initial state";
		return false;
	}
	model.registerState(t, currentPoint);
	double lastRecordedTime = t;	
	
	//run simulation
	//default time step is dt, but temporarily dt might be reduced to a smaller value
	for(recorded=1, current_dt1=current_dt2=dt, refinements=0; t<endTime; /* increment of t handled in loop */ ){
		// at this stage currentPoint is guaranteed to be a valid state
		// t should always correspond to currentPoint
		// t2 should always correspond to point2
		// current_dt1 is used for moving from currentPoint-->point2
		// current_dt2 is used for moving from currentPoint-->candidatePoint (= potential next currentPoint)
		
		// check runtime-out
		if((runtime_out_seconds>0) && (get_thread_monotonic_walltime_seconds()-start_runtime>=runtime_out_seconds)){
			warningMessage = "Simulation timed out";
			return false;
		}
				
		// get dynamics (k1) at currentPoint if needed
		if(!k1AlreadyCalculated){
			dynamics = model.getRateOfChangeAtState(t, currentPoint, k1, jumpPoint);
			if(dynamics==RequestedDynamicsForceJumpToState){
				currentPoint = jumpPoint;
				t = t + current_dt1;
				goto RK2_REGISTER_STATE; 
			}
			k1AlreadyCalculated = true;
		}
		
		// get point2 (in forward k1-direction)
		point2 = currentPoint + k1*current_dt1;
		t2 = t + current_dt1; // t2 should always correspond to point2
		
		// check if point2 crossed the boundary and move point2 backward as much as needed
		// point2 and current_dt1 will be adjusted if this routine returns true
		crossedBoundary = model.checkCrossedDomainBoundaryAndFix(t, currentPoint, current_dt1, point2);
		if(crossedBoundary == CrossedBoundaryYesButFixedBruteForce){ crossedBoundaryYesButFixedBruteForce = true; }
		else if(crossedBoundary == CrossedBoundaryYesButFixedByReducingTimeStep){ crossedBoundaryButFixedByReducingTimeStep = true; }
		
		// check if the time step should be reduced (as suggested by the model)
		while(model.checkShouldRefineTimeStep(t, currentPoint, current_dt1, point2) && (refinements<maxTimeStepRefinements)){
			current_dt1 /= 2;
			++refinements;
			point2 = currentPoint + k1*current_dt1;
			t2 = t + current_dt1; // t2 should always correspond to point2
		}		
		
		// get dynamics (k2) at point2
		dynamics = model.getRateOfChangeAtState(t2, point2, k2, jumpPoint);
		if(dynamics==RequestedDynamicsForceJumpToState){
			currentPoint = jumpPoint;
			t = t2;
			goto RK2_REGISTER_STATE; 
		}

		// use k1 & k2 to move the currentPoint forward, potentially after refining the time step multiple times
		kConsensus = (k1+k2)*0.5;
		candidatePoint = currentPoint + kConsensus*current_dt2;
		crossedBoundary = model.checkCrossedDomainBoundaryAndFix(t, currentPoint, current_dt2, candidatePoint);
		if(crossedBoundary==CrossedBoundaryYesButFixedBruteForce){
			if(refinements<maxTimeStepRefinements){
				++refinements; 
				current_dt2 /= 2;
				crossedBoundaryButFixedByReducingTimeStep = true;
				continue; // repeat the whole iteration with a smaller time step
			}else{
				// use simple Euler scheme for this time step, since kConsensus throws us out of the boundary at even arbitrarily small time steps
				currentPoint = point2;
				t = t2;
				crossedBoundaryYesButFixedUsingEuler = true;
				goto RK2_REGISTER_STATE;
			}
		}
				
		// check if the time step should be reduced (as suggested by the model)
		if(model.checkShouldRefineTimeStep(t, currentPoint, current_dt2, candidatePoint) && (refinements<maxTimeStepRefinements)){
			++refinements; 
			current_dt2 /= 2; 
			continue; // repeat the last part with a smaller time step
		}
		
		// seems everything worked out as normal, so set the new currentPoint
		currentPoint = candidatePoint;
		t = t + current_dt2; // t should always correspond to currentPoint
		goto RK2_REGISTER_STATE;
		
		// check and register new state if needed	
		RK2_REGISTER_STATE:
		// check sanity (some ODEs may explode)		
		if(std::isnan(t) || model.stateIsNaN(currentPoint)){
			warningMessage = "Solution is NaN";
			return false;
		}
		k1AlreadyCalculated = false;
		if(t-lastRecordedTime > recordingTimeStep){
			model.registerState(t, currentPoint);
			++recorded;
			lastRecordedTime = t;
			reporter(recorded, maxRecordedPoints, 1-(endTime-t)/simulationTime);
		}
		
		// initialize counters for next iteration (note: this step may be skipped by a 'continue' statement)
		current_dt1  = dt;
		current_dt2  = dt;
		refinements = 0;
	}
	
	// register final state if needed
	if(t>lastRecordedTime){
		model.registerState(t, currentPoint);
	}
	
	// construct warnings message if needed
	if(crossedBoundaryYesButFixedUsingEuler) warningMessage += string(warningMessage=="" ? "" : "\n") + "Simulation crossed domain boundary and was fixed by locally using forward Euler.";
	if(crossedBoundaryYesButFixedBruteForce) warningMessage += string(warningMessage=="" ? "" : "\n") + "Simulation crossed domain boundary and was fixed by brute-force adjusting the trajectory.";
	if(crossedBoundaryButFixedByReducingTimeStep) warningMessage += string(warningMessage=="" ? "" : "\n") + "Simulation crossed domain boundary and was fixed by locally reducing the time step.";
	if(simulationExceededMaxProcessingTime) warningMessage += string(warningMessage=="" ? "" : "\n") + "Simulation was aborted prematurely because it reached the maximum allowed processing time.";

	return true;
}




// simulate the deterministic trajectory of a speciation-extinction model for tree growth, calculating the diversities at specific times
// Optionally, the model can be integrated backwards in time, so that the "initial condition" (start_diversity) refers to the final true diversity of the tree.
// Typically, if 'coalescent' is true then 'reverse' should be true, unless the true diversity (i.e. of extant clades) is known for the past in the coalescent tree.
// [[Rcpp::export]]
Rcpp::List simulate_deterministic_diversity_growth_CPP(	const double 		birth_rate_intercept,	// (INPUT) intercept of Poissonian rate at which new tips are added to the tree
														const double 		birth_rate_factor,		// (INPUT) power-law factor of Poissonian rate at which new tips are added to the tree
														const double 		birth_rate_exponent,	// (INPUT) power-law exponent of Poissonian rate at which new tips are added to the tree
														const double 		death_rate_intercept,	// (INPUT) intercept of Poissonian rate at which extant tips are removed from the tree
														const double 		death_rate_factor,		// (INPUT) power-law factor of Poissonian rate at which extant tips are removed from the tree
														const double 		death_rate_exponent,	// (INPUT) power-law exponent of Poissonian rate at which extant tips are removed from the tree
														const double		rarefaction,			// (INPUT) optional rarefaction fraction to apply, i.e. fraction of tips remaining after rarefaction. Set this to 1 for no rarefaction. Note that rarefaction affects the full coalescent diversity curve in a non-uniform manner. If coalescent==false, rarefaction has no effect.
														const long			Nsplits,				// (INPUT) number of children to create at each diversification event. Must be at least 2. For a bifurcating tree this should be set to 2. If >2, then the tree will be multifurcating.
														const NumericVector	&times,					// (INPUT) times at which to calculate deterministic diversities, in ascending order
														const double		start_time,				// (INPUT) Time at which to start simulation. If reverse==true, then this should be the ending time of the tree (>= times.back()), otherwise this should be the beginning of the tree (<=times[0])
														const double		start_diversity,		// (INPUT) diversity of extant clades at start_time. If reverse==true, this is the visible diversity after rarefaction.
														const bool			reverse,				// (OUTPUT) if true, then the tree model is integrated in backward time direction. In that case, start_diversity is interpreted as the true diversity at times.back()
														const bool			coalescent,				// (INPUT) if true, the diversities at any time are replaced by the diversities that would remain in a coalescent tree (i.e. only including ancestors of extant tips)
														const bool			include_survival_chances,	// (INPUT) if true, the survival probabilities of clades (for each time point) will also be returned
														const bool			include_birth_rates,		// (INPUT) if true, the speciation rates corresponding to times[] will also be returned
														const bool			include_Nbirths,			// (INPUT) if true, then the total predicted birth events (up to each time point) will also be calculated and returned
														const double		runtime_out_seconds){		// (INPUT) max allowed runtime in seconds. If <=0, this option is ignored.
	const long NT = times.size();
	const double max_time = (reverse ? start_time : times[NT-1]);

	// initialize model for tree growth
	TreeSpeciationExtinctionModel model;
	model.birth_rate_intercept 	= birth_rate_intercept;
	model.birth_rate_factor 	= birth_rate_factor;
	model.birth_rate_exponent 	= birth_rate_exponent;
	model.death_rate_intercept 	= death_rate_intercept;
	model.death_rate_factor 	= death_rate_factor;
	model.death_rate_exponent 	= death_rate_exponent;
	model.Nsplits				= Nsplits;
	model.initial_diversity		= (reverse ? 1.0/rarefaction : 1.0) * start_diversity; // for curves simulated in reverse, apply the rarefaction to the initial condition (i.e. increase diversity from tau=0 to tau=epsilon)
	model.rarefaction			= (coalescent ? rarefaction : 1.0); // if not coalescent, ignore effects of rarefaction (i.e. simulate true diversities)
	if(reverse) model.set_reverse(max_time);
	
	// run simulation
	string warningMessage;
	const double default_dt = min((max_time-times[0])/NT, (reverse ? 1e-2*rarefaction*start_diversity/model.estimate_max_rate_of_change(max_time, start_diversity,false) : (max_time-times[0])/NT)); // guess appropriate simulation time step
	const bool success = RungeKutta2<TreeStateHistory,TreeSpeciationExtinctionModel,ProgressReporter>
								((reverse ? 0 : times[0]),
								default_dt+(reverse ? max_time-times[0] : max_time),
								default_dt,
								model,
								10*times.size(),
								2l,
								ProgressReporter(true),
								runtime_out_seconds,
								warningMessage);
	const long NMT = model.times.size();
	if(!success) return Rcpp::List::create(Rcpp::Named("success") = false, Rcpp::Named("error") = warningMessage ); // simulation failed
		
	if(reverse){
		// reverse time course of trajectory
		TreeStateHistory state;
		double time;
		for(long t=0, s; t<NMT/2; ++t){
			s = NMT-t-1;
			state = model.trajectory[t];
			model.trajectory[t] = model.trajectory[s];
			model.trajectory[s] = state;
			time = model.times[t];
			model.times[t] = model.times[s];
			model.times[s] = time;
		}
		// correct direction of accumulation of time, Nbirths & Ndeaths
		const double max_Nbirths = model.trajectory[0].Nbirths;
		const double max_Ndeaths = model.trajectory[0].Ndeaths;
		for(long t=0; t<NMT; ++t){
			model.times[t] = max_time - model.times[t];
			model.trajectory[t].Nbirths = max_Nbirths - model.trajectory[t].Nbirths;
			model.trajectory[t].Ndeaths = max_Ndeaths - model.trajectory[t].Ndeaths;
		}
	}
					
	// adjust (reduce) past diversities to reflect a coalescent tree
	// note that the effects of rarefaction are taken into account here
	if(coalescent) model.make_coalescent();
	
	// extract diversities from simulated trajectory
	std::vector<double> model_diversities(NMT);
	for(long t=0; t<NMT; ++t){
		model_diversities[t] = model.trajectory[t].diversity;
	}
				
	// interpolate simulated trajectory on requested time points
	long includedNewTimesStart, includedNewTimesEnd;
	std::vector<double> diversities;
	interpolateTimeSeriesAtTimes(	model.times,
									model_diversities,
									0,
									NMT-1,
									times,
									0,
									NT-1,
									includedNewTimesStart,
									includedNewTimesEnd,
									diversities);
									
	// calculate clade survival chances for requested time points
	// this only makes sense for reverse integrations
	std::vector<double> survival_chances;
	if(include_survival_chances){
 		// extract from model trajectory
 		std::vector<double> model_survival_chances(NMT);
		for(long t=0; t<NMT; ++t){
			model_survival_chances[t] = 1-model.trajectory[t].Pextinction;
		}
		// interpolate onto requested time points
		interpolateTimeSeriesAtTimes(	model.times,
										model_survival_chances,
										0,
										NMT-1,
										times,
										0,
										NT-1,
										includedNewTimesStart,
										includedNewTimesEnd,
										survival_chances);
	}
								
	// calculate speciation rates on returned time series, if needed
	std::vector<double> birth_rates;
	if(include_birth_rates){
		birth_rates.resize(NT);
		for(long t=0; t<NT; ++t){
			birth_rates[t] = model.get_speciation_rate_at_state(times[t],diversities[t]);
		}
	}
	
	// extract cumulative speciation & extinction events, if needed
	// double Nbirths = NAN_D;
	std::vector<double> Nbirths, Ndeaths;
 	if(include_Nbirths){
 		// extract from model trajectory
 		std::vector<double> model_Nbirths(NMT), model_Ndeaths(NMT);
		for(long t=0; t<NMT; ++t){
			model_Nbirths[t] = model.trajectory[t].Nbirths;
			model_Ndeaths[t] = model.trajectory[t].Ndeaths;
		}
		// interpolate onto requested time points
		interpolateTimeSeriesAtTimes(model.times, model_Nbirths, 0, NMT-1, times, 0, NT-1, includedNewTimesStart, includedNewTimesEnd, Nbirths);
		interpolateTimeSeriesAtTimes(model.times, model_Ndeaths, 0, NMT-1, times, 0, NT-1, includedNewTimesStart, includedNewTimesEnd, Ndeaths);
 	}
 		
	return Rcpp::List::create(	Rcpp::Named("success")			= true,
								Rcpp::Named("diversities") 		= Rcpp::wrap(diversities),
								Rcpp::Named("survival_chances")	= Rcpp::wrap(survival_chances),
								Rcpp::Named("birth_rates")		= Rcpp::wrap(birth_rates),
								Rcpp::Named("Nbirths")			= Nbirths,
								Rcpp::Named("Ndeaths")			= Ndeaths);
}





// Estimate past diversity, birth & death rates, based on a time series of coalescent diversities
// This reconstruction is non-parametric, i.e. no particular model is assumed (except for knowledge or constancy of the per-capita birth rate)
// Note: This function is currently only implemented for bifurcating trees. 
// It's possible to extend to multifurcating trees (Nsplits>2), but it requires solving (e.g. numerically) a higher-order polynomial at each age.
// Input:
//	A time series of coalescent diversities, i.e. as would be visible in a coalescent phylogenetic tree.
// 	Corresponding assumed per-capita birth rates. Alternatively, these can be assumed to be constant and be estimated directly from the coalescent_diversity time series.
// Output:
//  The estimated true past diversities (N(t))
//  The estimated corresponding death (extinction) rates (delta(t))
//  The probability of a size-1 clade surviving from each time point to the present (P(t))
//  The estimated total number of speciation & extinction events, over the considered age interval
// [[Rcpp::export]]
Rcpp::List reconstruct_past_diversity_from_coalescent_CPP(	const std::vector<double>	&times,						// (INPUT) 1D numeric array of size NT, listing times in ascending order. The last time point corresponds to the "present" (age=0)
															const std::vector<double>	&raw_coalescent_diversities,// (INPUT) 1D numeric array of size NT, listing coalescent diversities as visible after rarefaction. Should be unsmoothened.
															const std::vector<double> 	&birth_rates_pc,			// (INPUT) 1D numeric array of size NT, listing known or assumed per-capita birth rates. Can also be of size 1, in which case the same per-capita birth rate is assumed throughout. Can also be empty, in which case a constant per-capita birth rate is assumed and estimated from the last slope of the coalescent_diversity curve.
															const double				rarefaction,				// (INPUT) optional rarefaction fraction assumed to have occurred at the very end (fraction of kept tips). Set to 1.0 for no rarefaction.
															const double 				max_age,					// (INPUT) max age (distance from last time point) to consider for integrating total births & deaths. If <=0, all times are considered.
															const long					smoothing_span,				// (INPUT) Integer. Optional number of time points for smoothening the diversities time series via Savitzky-Golay-filter. If <=2, no smoothing is done. Smoothening the coalescent diversity can reduce the noise in the non-parametric reconstruction. 
															const long					smoothing_order){			// (INPUT) Integer. Optional polynomial order of the smoothing model.
	const long NT 			= times.size();
	const double max_time 	= times[NT-1];
	const double Nsplits	= 2; // currently not implemented for other Nsplits
	
	// smoothen time series if needed
	// using a Savitzky-Golay-filter of 2nd-order (1st order if smoothen==3)
	const bool smoothen = (smoothing_span>2);
	std::vector<double> smooth_coalescent_diversities(NT);
	if(smoothen){
		if(!smoothenTimeSeriesSavitzkyGolay(times, raw_coalescent_diversities, 0.0, smoothing_span, min(smoothing_span-2,smoothing_order), smooth_coalescent_diversities)){
			return Rcpp::List::create(Rcpp::Named("success") = false, Rcpp::Named("error") = "Smoothing of coalescent diversity time series failed");
		}
		smooth_coalescent_diversities[NT-1] = raw_coalescent_diversities[NT-1];
	}
	const std::vector<double> &coalescent_diversities = (smoothen ? smooth_coalescent_diversities : raw_coalescent_diversities);
	
	// determine latest (most recent) per-capita birth rate if needed
	double last_birth_rate_pc;
	const bool constant_birth_rate_pc = (birth_rates_pc.size()<=1);
	if(birth_rates_pc.size()==0){
		// estimate from slope of raw_coalescent_diversities
		// the raw (unsmoothened) data is preferred for estimation of birth_rate_pc, so as to preserve the underlying (by assumption) exponential structure
		last_birth_rate_pc = log(raw_coalescent_diversities[NT-1]/raw_coalescent_diversities[NT-2])/(times[NT-1] - times[NT-2]);
		// correct for the effects of rarefaction
		const double zeta = 1-rarefaction;
		if(rarefaction<1) last_birth_rate_pc *= rarefaction/(Nsplits-1-Nsplits*zeta+pow(zeta,Nsplits));
	}else{
		// either birth_rates_oc are provided for all time points, or as a single constant number
		last_birth_rate_pc = birth_rates_pc[birth_rates_pc.size()-1];
	}
	if(last_birth_rate_pc<0) return Rcpp::List::create(Rcpp::Named("success") = false, Rcpp::Named("error") = "Latest per-capita birth rate is negative or zero");
	
	// determine first considered time point
	long first_time_point;
	if(max_age>0){
		first_time_point = NT-1;
		for(long t=0; t<NT; ++t){
			if((max_time-times[t])<=max_age){
				first_time_point = t;
				break;
			}
		}
	}else{
		first_time_point = 0;
	}
			
	// calculate survival chances & true diversities over time
	std::vector<double> survival_chances(NT), true_diversities(NT), nu(NT);
	for(long t=NT-1, tl, tr; t>=0; --t){
		tl = (t==0 ? 0 : t-1);
		tr = (t==0 ? 1 : t); // use left-sided derivative, to avoid the discontinuity at age=0
		const double birth_rate_pc = (constant_birth_rate_pc ? last_birth_rate_pc : birth_rates_pc[t]);
		nu[t] = -log(coalescent_diversities[tr]/coalescent_diversities[tl])/(times[tr]-times[tl]); // (1/Nc)*dNc/dtau = dlog(Nc)/dtau, where Nc: coalescent_diversity
		survival_chances[t] = -nu[t]/birth_rate_pc;
		if((survival_chances[t]<=1e-2/coalescent_diversities[NT-1]) && (coalescent_diversities[t]>0)){
			// diversity not defined (survival seems spuriously low)
			true_diversities[t] = NAN_D;
		}else if(survival_chances[t]>1){
			// diversity not defined (survival chance > 1)
			true_diversities[t] = NAN_D;			
		}else{
			true_diversities[t] = coalescent_diversities[t]/survival_chances[t];
		}
	}
	true_diversities[NT-1] = coalescent_diversities[NT-1]/rarefaction; // coalescent diversity last time point is assumed to be equal to true diversity multiplied by rarefaction fraction
	
	// calculate birth & death rates & per-capita (exponential) pulled growth rate
	std::vector<double> death_rates(NT), birth_rates(NT), pulled_diversification_rates(NT);
	for(long t=NT-1, tl, tr; t>=0; --t){
		tl = (t==0 ? 0 : t-1);
		tr = (t==0 ? 1 : t); // use left-sided derivative, to avoid the discontinuity at age=0
		const double birth_rate_pc = (constant_birth_rate_pc ? last_birth_rate_pc : birth_rates_pc[t]);
		birth_rates[t] = birth_rate_pc*true_diversities[t];
		death_rates[t] = max(0.0, -(true_diversities[tr] - true_diversities[tl])/(times[tr] - times[tl]) + birth_rates[t]); // death rate must be non-negative
		pulled_diversification_rates[t] = (log(-coalescent_diversities[tr]/nu[tr]) - log(-coalescent_diversities[tl]/nu[tl]))/(times[tr] - times[tl]);
	}
	
	
	
	// calculate total number of births & deaths
	const double total_births = integrate1D(times,birth_rates,first_time_point,NT-1,true);
	const double total_deaths = integrate1D(times,death_rates,first_time_point,NT-1,true);
	
	return Rcpp::List::create(	Rcpp::Named("success")						= true,
								Rcpp::Named("true_diversities")				= Rcpp::wrap(true_diversities),
								Rcpp::Named("birth_rates") 					= Rcpp::wrap(birth_rates),
								Rcpp::Named("death_rates") 					= Rcpp::wrap(death_rates),
								Rcpp::Named("survival_chances")				= Rcpp::wrap(survival_chances),
								Rcpp::Named("total_births")					= total_births,
								Rcpp::Named("total_deaths")					= total_deaths,
								Rcpp::Named("last_birth_rate_pc")			= last_birth_rate_pc,
								Rcpp::Named("last_death_rate_pc")			= death_rates[NT-1]/true_diversities[NT-1],
								Rcpp::Named("pulled_diversification_rates")	= Rcpp::wrap(pulled_diversification_rates));
}



// Estimate past diversity, birth & death rates, based on a time series of true diversities
// Input:
//	A time series of diversities, i.e. as would be visible in a non-coalescent phylogenetic tree (including extinct species).
// 	Corresponding assumed per-capita birth rates.
// Output:
//  The estimated corresponding death (extinction) rates (delta(t))
//  The probability of a size-1 clade surviving from each time point to the present (P(t))
//  The estimated total number of speciation & extinction events, over the considered age interval
// [[Rcpp::export]]
Rcpp::List reconstruct_past_diversifications_CPP(	const std::vector<double>	&times,				// (INPUT) 1D numeric array of size NT, listing times in ascending order. The last time point corresponds to the "present" (age=0)
													const std::vector<double>	&raw_diversities,	// (INPUT) 1D numeric array of size NT, listing true diversities. Should be unsmoothened.
													const std::vector<double> 	&birth_rates_pc,	// (INPUT) 1D numeric array of size NT, listing known or assumed per-capita birth rates. Can also be of size 1, in which case the same per-capita birth rate is assumed throughout.
													const double				rarefaction,		// (INPUT) optional rarefaction fraction, to apply when calculating survival chances
													const long					Nsplits,			// (INPUT)
													const double 				max_age,			// (INPUT) max age (distance from last time point) to consider for integrating total births & deaths. If <=0, all times are considered.
													const long					smoothing_span,		// (INPUT) Integer. Optional number of time points for smoothening the diversities time series via Savitzky-Golay-filter. If <=2, no smoothing is done. Smoothening the coalescent diversity can reduce the noise in the non-parametric reconstruction. 
													const long					smoothing_order){	// (INPUT) Integer. Optional polynomial order of the smoothing model.
	const long NT 			= times.size();
	const double max_time 	= times[NT-1];
	const bool const_birth_rate_pc = (birth_rates_pc.size()==1);
	
	// smoothen time series if needed
	// using a Savitzky-Golay-filter of 2nd-order (1st order if smoothen==3)
	const bool smoothen = (smoothing_span>2);
	std::vector<double> smooth_diversities(NT);
	if(smoothen){
		if(!smoothenTimeSeriesSavitzkyGolay(times, raw_diversities, 0.0, smoothing_span, min(smoothing_span-2,smoothing_order), smooth_diversities)){
			return Rcpp::List::create(Rcpp::Named("success") = false, Rcpp::Named("error") = "Smoothing of coalescent diversity time series failed");
		}
		smooth_diversities[NT-1] = raw_diversities[NT-1];
	}
	const std::vector<double> &true_diversities = (smoothen ? smooth_diversities : raw_diversities);
			
	// determine first considered time point
	long first_time_point;
	if(max_age>0){
		first_time_point = NT-1;
		for(long t=0; t<NT; ++t){
			if((max_time-times[t])<=max_age){
				first_time_point = t;
				break;
			}
		}
	}else{
		first_time_point = 0;
	}
			
	// calculate birth & death rates & per-capita (exponential) growth rates
	std::vector<double> birth_rates(NT), death_rates(NT), diversification_rates(NT);
	for(long t=NT-1, tl, tr; t>=0; --t){
		tl = (t==0 ? 0 : t-1);
		tr = (t==NT-1 ? t : t+1);
		const double birth_rate_pc = (const_birth_rate_pc ? birth_rates_pc[0] : birth_rates_pc[t]);
		birth_rates[t] = birth_rate_pc*true_diversities[t];
		death_rates[t] = (true_diversities[tr] - true_diversities[tl])/(times[tr] - times[tl]) + birth_rates[t];
		diversification_rates[t] = (log(true_diversities[tr])-log(true_diversities[tl]))/(times[tr] - times[tl]);
	}
	
	// integrate backwards to calculate probability of extinction (and thus survival chances)
	// probability of extinction at age=0 begins at (1-rarefaction)
	std::vector<double> survival_chances(NT);
	std::vector<double> coalescent_diversities(NT);
	survival_chances[NT-1] = rarefaction;
	coalescent_diversities[NT-1] = true_diversities[NT-1]*survival_chances[NT-1];
	double Pextinction = (1-rarefaction); 
	for(long t=NT-2; t>=0; --t){
		const double dt = times[t+1]-times[t];
		const double birth_rate_pc = (const_birth_rate_pc ? birth_rates_pc[0] : birth_rates_pc[t]);
		const double death_rate_pc = 0.5*(death_rates[t+1]/true_diversities[t+1] + death_rates[t]/true_diversities[t]);
		Pextinction += dt*(death_rate_pc - Pextinction*(birth_rate_pc+death_rate_pc) + pow(Pextinction,Nsplits)*birth_rate_pc);
		survival_chances[t] = 1 - Pextinction;
		coalescent_diversities[t] = true_diversities[t]*survival_chances[t];
	}

	
	// calculate total number of births & deaths
	const double total_births = integrate1D(times,birth_rates,first_time_point,NT-1,true);
	const double total_deaths = integrate1D(times,death_rates,first_time_point,NT-1,true);
	
	return Rcpp::List::create(	Rcpp::Named("success")					= true,
								Rcpp::Named("birth_rates") 				= Rcpp::wrap(birth_rates),
								Rcpp::Named("death_rates") 				= Rcpp::wrap(death_rates),
								Rcpp::Named("survival_chances")			= Rcpp::wrap(survival_chances),
								Rcpp::Named("coalescent_diversities")	= Rcpp::wrap(coalescent_diversities),
								Rcpp::Named("total_births")				= total_births,
								Rcpp::Named("total_deaths")				= total_deaths,
								Rcpp::Named("diversification_rates")	= Rcpp::wrap(diversification_rates));
}



// Based on the timings of birth (speciation) & death (extinction) events, calculate the diversity-over-time curve
// diversity at time t = start_diversity + (Nsplits-1)*total_birth_counts_until[t] - total_death_counts_until[t]
// [[Rcpp::export]]
Rcpp::List get_diversities_from_birth_and_death_events_CPP(	const NumericVector &times,
															const NumericVector &birth_times,		// 1D array of size NT, listing birth times in ascending order
															const NumericVector &death_times,		// 1D array of size NT, listing death times in ascending order
															const double 		start_diversity,	// (INPUT) diversity prior to any recorded birth & death events. Will typically be 0 or 1.
															const double		Nsplits){			// (INPUT) number of new species emerging from each speciation event
	const long NT = times.size();
	const long NB = birth_times.size();
	const long ND = death_times.size();
	
	// bin birth events into slots
	// births_per_time[t] = number of births that occurred between times[t-1] and times[t]
	std::vector<double> births_per_time(NT,0);
	long time_slot = 0;
	for(long b=0; b<NB; ++b){
		const double birth_time = birth_times[b];
		while((times[time_slot]<birth_time) && (time_slot<NT)){ ++time_slot; }
		if(time_slot<NT) ++births_per_time[time_slot];
	}
	
	// bin death events into slots
	// deaths_per_time[t] = number of deaths that occurred between times[t-1] and times[t]
	std::vector<double> deaths_per_time(NT,0);
	time_slot = 0;
	for(long d=0; d<ND; ++d){
		const double death_time = death_times[d];
		while((times[time_slot]<death_time) && (time_slot<NT)){ ++time_slot; }
		if(time_slot<NT) ++deaths_per_time[time_slot];
	}
	
	// count cumulative births & deaths at each time slot
	double diversity = start_diversity;
	std::vector<double> diversities(NT);
	for(long t=0; t<NT; ++t){
		diversity += (Nsplits-1)*births_per_time[t] - deaths_per_time[t];
		diversities[t] = diversity;
	}
	
	return Rcpp::List::create(Rcpp::Named("diversities") = diversities);
}






// ##########################################################


#pragma mark -
#pragma mark Building auxiliary data structures
#pragma mark 



// retrieve the parent clade of each clade
// returns a 1D array of size (Ntips+Nnodes), with elements in (Ntips):(Ntips+Nnodes-1)
// for the root, the parent will be set to -1 (i.e. not available)
template<class ARRAY_TYPE>
void get_parent_per_clade(	const long			Ntips,
							const long 			Nnodes,
							const long			Nedges,
							const ARRAY_TYPE 	&tree_edge,			// (INPUT) 2D array (in row-major format) of size Nedges x 2
							std::vector<long>	&clade2parent){		// (OUTPUT) 1D array of size (Nnodes+Ntips), normally with values in 0:(Nnodes+Ntipes-1). Clades with no parent will have have value -1.
	clade2parent.assign(Ntips+Nnodes, -1);
	for(long edge=0; edge<Nedges; ++edge){
		clade2parent[tree_edge[edge*2+1]] = tree_edge[edge*2+0];
	}
}


// given a mapping clade-->parent (with the root's parent being <0), determine the root of a tree by traversing upwards
// this function is more efficient than just searching for the clade without parent
long get_root_from_clade2parent(const long 				first_guess,		// (INPUT) first guess for what the root should be. The better the guess (i.e. the closer to the root), the faster the true root will be found. A good first guess is usually Ntips;
								const std::vector<long> &clade2parent){		// (INPUT) 1D array of size (Nnodes+Ntips), normally with values in 0:(Nnodes+Ntipes-1). Clades with no parent will have have value -1.
	long clade = first_guess;
	while(clade2parent[clade]>=0){
		clade = clade2parent[clade];
	}
	return clade;
}


// given a mapping clade-->incoming_edge (with the root's incoming edge being <0), determine the root of a tree by traversing upwards
// this function is more efficient than just searching for the clade without parent
template<class ARRAY_TYPE>
long get_root_from_incoming_edge_per_clade(	const long 				first_guess,				// (INPUT) first guess for what the root should be. The better the guess (i.e. the closer to the root), the faster the true root will be found. A good first guess is usually Ntips;
											const ARRAY_TYPE		&tree_edge,					// (INPUT) 2D array of size Nedges x 2, in row-major format
											const std::vector<long> &incoming_edge_per_clade){	// (INPUT) 1D array of size (Nnodes+Ntips), normally with values in 0:(Nedges-1). Clades with no incoming edge will have have value -1.
	long clade = first_guess;
	while(incoming_edge_per_clade[clade]>=0){
		clade = tree_edge[2*incoming_edge_per_clade[clade]+0];
	}
	return clade;
}



// given a mapping clade-->inout edges, determine the root of a tree by traversing upwards
// this function is more efficient than just searching for the clade without parent
template<class ARRAY_TYPE>
long get_root_from_clade2inout_edges(	const long 				first_guess,				// (INPUT) first guess for what the root clade should be. The better the guess (i.e. the closer to the root), the faster the true root will be found. A good first guess is usually Ntips;
										const ARRAY_TYPE		&tree_edge,					// (INPUT) 2D array of size Nedges x 2, in row-major format
										const std::vector<long>	&clade2first_inout_edge,	// (INPUT) 1D array of size Nclades, with values in 0:(2*Nedges-1), mapping clades to their first incoming or outgoing edge.
										const std::vector<long>	&clade2last_inout_edge,		// (INPUT) 1D array of size Nclades, with values in 0:(2*Nedges-1), mapping clades to their last incoming or outgoing edge.
										const std::vector<long>	&inout_edges){				// (INPUT) 1D array of size 2*Nedges, with values in 0:(Nedges-1). Maps internal edge indices (i.e. as listed in clade2first_inout_edge[] and clade2last_inout_edge[]) to original edge indices.

	long clade = first_guess;
	bool found_parent;
	do{
		found_parent = false;
		for(long e=clade2first_inout_edge[clade], edge; e<=clade2last_inout_edge[clade]; ++e){
			edge  = inout_edges[e];
			if(tree_edge[2*edge+1]==clade){
				// found incoming edge, so transition to parent
				clade = tree_edge[2*edge+0];
				found_parent = true;
				break;
			}
		}
	}while(found_parent);
	return clade;
}






// find the incoming edge per clade
// returns a 1D array of size (Ntips+Nnodes), with elements in 0,..,Nedges-1
// for the root, the incoming edge will be set to -1 (i.e. not available)
template<class ARRAY_TYPE>
void get_incoming_edge_per_clade(	const long			Ntips,
									const long 			Nnodes,
									const long			Nedges,
									const ARRAY_TYPE	&tree_edge,					// (INPUT) 2D array of size Nedges x 2, in row-major format
									std::vector<long>	&incoming_edge_per_clade){	// (OUTPUT) 1D array of size (Nnodes+Ntips), with values in 0,..,Nedges-1. 
	incoming_edge_per_clade.assign(Ntips+Nnodes,-1);
	for(long edge=0; edge<Nedges; ++edge){
		incoming_edge_per_clade[tree_edge[edge*2+1]] = edge;
	}
}



// for a given phylogenetic tree, create list of incoming edges for each clade
// normally each clade has either 0 or 1 incoming edges
// the tree need not be rooted, and may include multifurcations and monofurcations
// [[Rcpp::export]]
std::vector<std::vector<long> > get_incoming_edges_per_clade_CPP(	const long				Ntips,
																	const long 				Nnodes,
																	const long				Nedges,
																	const std::vector<long>	&tree_edge){	// (INPUT) 2D array of size Nedges x 2, in row-major format
	std::vector<std::vector<long> > incoming_edges_per_clade(Ntips+Nnodes);
	for(long edge=0, child; edge<Nedges; ++edge){
		child = tree_edge[edge*2+1];
		incoming_edges_per_clade[child].push_back(edge);
	}
	return incoming_edges_per_clade;
}




// find the incoming edge per clade
// returns a 1D array of size (Ntips+Nnodes), with elements in 0,..,Nedges-1
// for the root, the incoming edge will be set to -1 (i.e. not available)
template<class ARRAY_TYPE>
void get_child_count_per_node(	const long			Ntips,
								const long 			Nnodes,
								const long			Nedges,
								const ARRAY_TYPE	&tree_edge,					// (INPUT) 2D array of size Nedges x 2, in row-major format
								std::vector<long>	&child_count_per_node){		// (OUTPUT) 1D array of size Nnodes
	child_count_per_node.assign(Nnodes,0);
	for(long edge=0; edge<Nedges; ++edge){
		if(tree_edge[edge*2+0]<Ntips) continue;
		child_count_per_node[tree_edge[edge*2+0]-Ntips] += 1;
	}
}




// determine which nodes are basal (i.e. all their children are tips)
template<class ARRAY_TYPE>
void determine_basal_nodes(	const long			Ntips,
							const long 			Nnodes,
							const long			Nedges,
							const ARRAY_TYPE	&tree_edge,			// (INPUT) 2D array of size Nedges x 2, in row-major format
							std::vector<char>	&node_is_basal){	// (OUTPUT) 1D array of size Nnodes
	node_is_basal.assign(Nnodes,1);
	for(long edge=0; edge<Nedges; ++edge){
		if(tree_edge[edge*2+0]<Ntips) continue;
		if(tree_edge[edge*2+1]>=Ntips) node_is_basal[tree_edge[edge*2+0]-Ntips] = 0;
	}
}





// determine root of a tree (as the node with no incoming edge)
template<class ARRAY_TYPE>
long get_root_clade(const long			Ntips,
					const long 			Nnodes,
					const long			Nedges,
					const ARRAY_TYPE	&tree_edge){			// (INPUT) 2D array (in row-major format) of size Nedges x 2
	const long Nclades = Ntips+Nnodes;
	std::vector<bool> clade_has_parent(Nclades,false);
	for(long edge=0; edge<Nedges; ++edge){
		clade_has_parent[tree_edge[edge*2+1]] = true;
	}
	for(long c=0; c<Nclades; ++c){
		if(!clade_has_parent[c]) return c;
	}
	return -1;
}



// Find the root of a tree (= node with no parent)
// This is an Rcpp wrapper for get_root_clade()
// [[Rcpp::export]]
long get_root_clade_CPP(const long			Ntips,
						const long 			Nnodes,
						const long			Nedges,
						const IntegerVector	&tree_edge){	// (INPUT) 2D array (in row-major format) of size Nedges x 2
	return get_root_clade(Ntips, Nnodes, Nedges, tree_edge);
}






// returns a mapping from nodes to their child clades (Nclades = Ntips+Nnodes)
// the tree can be multifurcating, and can also include nodes with a single child
// root must specify the root index in the tree (typically root = Ntips)
// Returned values:
//	node2first_child[p] will be an index pointing node p (p=0:(Nnodes-1)) to its first child in children[]
//	node2last_child[p] will be an index pointing node p (p=0:(Nnodes-1)) to its last child in children[]
// 	children[] will be a list of clade indices (i.e. in 0:(Nclades-1)), such that children[node2first_child[p]],...,children[node2last_child[p]] is the set of children of node p
void get_children_per_node(	const long			Ntips,
							const long 			Nnodes,
							const long			Nedges,
							const long 			root, 				// (INPUT) index of root node, i.e. an integer in 0:(Ntips+Nnodes-1)
							const IntegerVector &tree_edge, 		// (INPUT) 2D array (in row-major format) of size Nedges x 2
							std::vector<long>	&node2first_child,	// (OUTPUT) 1D array of size Nnodes, with values in 0:(Nedges-1).
							std::vector<long>	&node2last_child,	// (OUTPUT) 1D array of size Nnodes, with values in 0:(Nedges-1).
							std::vector<long>	&children){			// (OUTPUT) 1D array of size Nedges, with values in 0:(Nclades-1).
	// Terminology in this function:
	// 	'node' runs from 0:(Nnodes-1)
	// 	'tip' runs from 0:(Ntips-1)
	// 	'parent' and 'child' runs from 0:(Ntips+Nnodes-1)
	// 	'edge' runs from 0:(Nedges-1)
	// Recall that:
	// 	tree_edge[] is of size Nedge x 2 (flattened in row-major-format), with entries in 0:(Ntips+Nnodes-1)
	children.resize(Nedges);
	node2first_child.resize(Nnodes);
	node2last_child.resize(Nnodes);
	long node;

	// determine number of children per parent
	// child_count_per_node[n] will be the number of direct children of node n (n=0:(Nnodes-1))
	std::vector<long> child_count_per_node(Nnodes, 0);
	for(long e=0; e<Nedges; ++e){
		child_count_per_node[tree_edge[e*2+0] - Ntips] += 1;
	}
	
	// collect children per parent
	node2first_child[0] = 0;
	node2last_child[0]  = node2first_child[0]+child_count_per_node[0] - 1;
	if(Nnodes>1){
		for(long n=1; n<Nnodes; ++n){
			node2first_child[n] = node2last_child[n-1]+1;
			node2last_child[n]  = node2first_child[n]+child_count_per_node[n] - 1;
		}
	}
	for(long e=0; e<Nedges; ++e){
		node = tree_edge[e*2+0] - Ntips;
		children[node2first_child[node]+child_count_per_node[node]-1] = tree_edge[e*2+1];
		child_count_per_node[node] -= 1;
	}		
}



// Calculate lookup tables mapping nodes to their outgoing (children) edges
// Requirements:
//    The tree can be multifurcating, and can also include nodes with a single child
//    The tree can be rooted or unrooted (all information on edge direction is taken from the tree_edge[] table)
// Returned values:
//	  node2first_edge[p] will be an index pointing node p (p=0:(Nnodes-1)) to edges[]
//	  node2last_edge[p] will be an index pointing node p (p=0:(Nnodes-1)) to edges[]
// 	  edges[] will be a list of edge indices (i.e. in 0:(Nedges-1)), such that edges[node2first_edge[p]],...,edges[node2last_edge[p]] is the set of edges leaving node p
template<class ARRAY_TYPE>
void get_node2edge_mappings(const long			Ntips,
							const long 			Nnodes,
							const long			Nedges,
							const ARRAY_TYPE	&tree_edge, 		// (INPUT) 2D array (in row-major format) of size Nedges x 2
							std::vector<long>	&node2first_edge,	// (OUTPUT) 1D array of size Nnodes, with values in 0:(Nedges-1), mapping nodes to their first outgoing edge.
							std::vector<long>	&node2last_edge,	// (OUTPUT) 1D array of size Nnodes, with values in 0:(Nedges-1), mapping nodes to their last outgoing edge.
							std::vector<long>	&edges){			// (OUTPUT) 1D array of size Nedges, with values in 0:(Nedges-1). Maps internal edge indices (i.e. as listed in node2first_edge[] and node2last_edge[]) to original edge indices.
	// Terminology in this function:
	// 	'node' runs from 0:(Nnodes-1)
	// 	'tip' runs from 0:(Ntips-1)
	// 	'parent' and 'child' runs from 0:(Ntips+Nnodes-1)
	// 	'edge' runs from 0:(Nedges-1)
	// Recall that:
	// 	tree_edge[] is of size Nedge x 2 (flattened in row-major-format), with entries in 0:(Ntips+Nnodes-1)

	edges.resize(Nedges);
	node2first_edge.resize(Nnodes);
	node2last_edge.resize(Nnodes);

	// determine number of children/edges per parent
	// child_count_per_node[n] will be the number of direct children of node n (n=0:(Nnodes-1))
	std::vector<long> child_count_per_node(Nnodes, 0);
	for(long e=0; e<Nedges; ++e){
		child_count_per_node[tree_edge[e*2+0] - Ntips] += 1;
	}
	// collect children per parent
	node2first_edge[0] = 0;
	node2last_edge[0]  = node2first_edge[0]+child_count_per_node[0] - 1;
	if(Nnodes>1){
		for(long n=1; n<Nnodes; ++n){
			node2first_edge[n] = node2last_edge[n-1]+1;
			node2last_edge[n]  = node2first_edge[n]+child_count_per_node[n] - 1;
		}
	}
	for(long e=0, node; e<Nedges; ++e){
		node = tree_edge[e*2+0] - Ntips;
		edges[node2first_edge[node]+child_count_per_node[node]-1] = e;
		child_count_per_node[node] -= 1;
	}
}




// returns a list of all nodes (and optionally tips) of a tree, such that each node appears prior to its children, and such that nodes closer to the root (in terms of branching counts) appear first
// also returns a list mapping nodes to their outgoing (children) edges (e.g. as listed in tree_edge)
// the tree can be multifurcating, and can also include nodes with a single child
// root must specify the root index in the tree (typically root = Ntips)
// Returned values:
//	queue: A 1D array of integers in 0:(Ntips+Nnodes-1) if include_tips==true, or (Ntips):(Ntips+Nnodes-1) if include_tips==false
//	node2first_edge[p] will be an index pointing node p (p=0:(Nnodes-1)) to edges[]
//	node2last_edge[p] will be an index pointing node p (p=0:(Nnodes-1)) to edges[]
// 	edges[] will be a list of edge indices (i.e. in 0:(Nedges-1)), such that edges[node2first_edge[p]],...,edges[node2last_edge[p]] is the set of edges leaving node p
template<class ARRAY_TYPE>
void get_tree_traversal_root_to_tips(	const long			Ntips,
										const long 			Nnodes,
										const long			Nedges,
										const long 			root, 							// (INPUT) index of root node, i.e. an integer in 0:(Ntips+Nnodes-1)
										const ARRAY_TYPE	&tree_edge, 					// (INPUT) 2D array (in row-major format) of size Nedges x 2
										const bool			include_tips,					// (INPUT) if true, then tips are included in the returned queue[]. This does not affect the returned arrays node2first_edge[], node2last_edge[], edges[].
										const bool			precalculated_edge_mappings,	// (INPUT) if true, then the edge mapping tables node2first_edge[], node2last_edge[] and edges[] are taken as is. Otherwise, they are calculated from scratch.
										std::vector<long>	&queue,							// (OUTPUT) 1D array of size Nnodes if include_tips==false, or size (Ntips+Nnodes) if include_tips==true.
										std::vector<long>	&node2first_edge,				// (INPUT/OUTPUT) 1D array of size Nnodes, with values in 0:(Nedges-1), mapping nodes to their first outgoing edge. Either pre-calculated, or to be calculated by this function.
										std::vector<long>	&node2last_edge,				// (INPUT/OUTPUT) 1D array of size Nnodes, with values in 0:(Nedges-1), mapping nodes to their last outgoing edge. Either pre-calculated, or to be calculated by this function.
										std::vector<long>	&edges,							// (INPUT/OUTPUT) 1D array of size Nedges, with values in 0:(Nedges-1). Maps internal edge indices (i.e. as listed in node2first_edge[] and node2last_edge[]) to original edge indices. Either pre-calculated, or to be calculated by this function.
										const bool			verbose,
										const string		&verbose_prefix){
	// get node-->edge mappings if needed
	if(!precalculated_edge_mappings){
		get_node2edge_mappings(	Ntips,
								Nnodes,
								Nedges,
								tree_edge,
								node2first_edge,
								node2last_edge,
								edges);
	}
	
	// fill queue from root to tips
	long child,node;
	queue.clear();
	queue.reserve(include_tips ? Ntips+Nnodes : Nnodes);
	queue.push_back(root);
	long queue_pointer = 0;
	while(queue_pointer<queue.size()){
		node = queue[queue_pointer] - Ntips;
		queue_pointer += 1;
		if(node<0) continue; // queue[queue_pointer] was actually a tip, not an internal node
		if(node2first_edge[node]>node2last_edge[node]){
			// this should not happen (every node needs to have at least one child)
			if(verbose) Rcout << verbose_prefix << "WARNING: Node " << node << " has no children\n";
			continue;
		}
		for(long ei=node2first_edge[node]; ei<=node2last_edge[node]; ++ei){
			child = tree_edge[edges[ei]*2+1];
			if((!include_tips) && (child<Ntips)) continue; // this child is a tip, so skip as requested
			// append child to queue
			queue.push_back(child);
		}
	}
}


// R wrapper for the function get_tree_traversal_root_to_tips(..) from above
// [[Rcpp::export]]
Rcpp::List get_tree_traversal_root_to_tips_CPP(	const long			Ntips,
												const long 			Nnodes,
												const long			Nedges,
												const IntegerVector	&tree_edge, 		// (INPUT) 2D array (in row-major format) of size Nedges x 2
												const bool			include_tips){		// (INPUT) if true, then tips are included in the returned queue[]. This does not affect the returned arrays node2first_edge[], node2last_edge[], edges[].
	// determine root
	const long root = get_root_clade(Ntips, Nnodes, Nedges, tree_edge);
	
	// get tree traversal
	std::vector<long> queue, node2first_edge, node2last_edge, edges;
	get_tree_traversal_root_to_tips(Ntips,
									Nnodes,
									Nedges,
									root,
									tree_edge,
									include_tips,
									false,
									queue,
									node2first_edge,
									node2last_edge,
									edges,
									false,
									"");
	return Rcpp::List::create(	Rcpp::Named("queue") 			= Rcpp::wrap(queue),
								Rcpp::Named("node2first_edge") 	= Rcpp::wrap(node2first_edge),
								Rcpp::Named("node2last_edge")	= Rcpp::wrap(node2last_edge),
								Rcpp::Named("edges")			= Rcpp::wrap(edges));

}



template<class ARRAY_TYPE>
void get_tree_traversal_depth_first_search(	const long			Ntips,
											const long 			Nnodes,
											const long			Nedges,
											const long 			root, 							// (INPUT) index of root node, i.e. an integer in 0:(Ntips+Nnodes-1)
											const ARRAY_TYPE	&tree_edge, 					// (INPUT) 2D array (in row-major format) of size Nedges x 2
											const bool			include_tips,					// (INPUT) if true, then tips are included in the returned queue[]. This does not affect the returned arrays node2first_edge[], node2last_edge[], edges[].
											const bool			precalculated_edge_mappings,	// (INPUT) if true, then the edge mapping tables node2first_edge[], node2last_edge[] and edges[] are taken as is. Otherwise, they are calculated from scratch.
											std::vector<long>	&queue,							// (OUTPUT) 1D array of size Nnodes if include_tips==false, or size (Ntips+Nnodes) if include_tips==true.
											std::vector<long>	&node2first_edge,				// (INPUT/OUTPUT) 1D array of size Nnodes, with values in 0:(Nedges-1), mapping nodes to their first outgoing edge. Either pre-calculated, or to be calculated by this function.
											std::vector<long>	&node2last_edge,				// (INPUT/OUTPUT) 1D array of size Nnodes, with values in 0:(Nedges-1), mapping nodes to their last outgoing edge. Either pre-calculated, or to be calculated by this function.
											std::vector<long>	&edges){						// (INPUT/OUTPUT) 1D array of size Nedges, with values in 0:(Nedges-1). Maps internal edge indices (i.e. as listed in node2first_edge[] and node2last_edge[]) to original edge indices. Either pre-calculated, or to be calculated by this function.

	// get node-->edge mappings if needed
	if(!precalculated_edge_mappings){
		get_node2edge_mappings(	Ntips,
								Nnodes,
								Nedges,
								tree_edge,
								node2first_edge,
								node2last_edge,
								edges);
	}
	
	// fill stack in depth-first-search direction
	// use a scratch_stack for traversing nodes
	long node,clade;
	std::vector<long> scratch_stack;
	scratch_stack.reserve(floor(2*log(Ntips)/log(2.0))); // rough estimate of typical tree depth x 2. scratch_stack may be resized along the way if needed.
	scratch_stack.push_back(root);
	queue.clear();
	queue.reserve(include_tips ? Ntips+Nnodes : Nnodes);
	while(scratch_stack.size()>0){
		clade = scratch_stack.back();
		scratch_stack.pop_back();
		if(include_tips || (clade>=Ntips)) queue.push_back(clade);
		if(clade>=Ntips){
			node = clade - Ntips;
			for(long ei=node2first_edge[node]; ei<=node2last_edge[node]; ++ei){
				scratch_stack.push_back(tree_edge[edges[ei]*2+1]);
			}
		}
	}
}



// given a phylogenetic tree, create a list of outgoing edges for each clade
// the tree need not be rooted, and may include monofurcations and multifurcations
// [[Rcpp::export]]
std::vector<std::vector<long> > get_outgoing_edges_per_clade_CPP(	const long				Ntips,
																	const long 				Nnodes,
																	const long				Nedges,
																	const std::vector<long>	&tree_edge){ // (INPUT) 2D array (in row-major format) of size Nedges x 2
	const long Nclades = Ntips+Nnodes;
	std::vector<std::vector<long> > edges_per_clade(Nclades);

	// determine number of edges per clade
	// edge_count_per_clade[n] will be the number of direct children of node n (n=0:(Nnodes-1))
	std::vector<long> edge_count_per_clade(Nclades, 0);
	for(long e=0; e<Nedges; ++e){
		edge_count_per_clade[tree_edge[e*2+0]] += 1;
	}
	// collect outgoing edges per clade
	for(long clade=0; clade<Nclades; ++clade){
		edges_per_clade.reserve(edge_count_per_clade[clade]);
	}
	for(long edge=0, parent; edge<Nedges; ++edge){
		parent = tree_edge[edge*2+0];
		edges_per_clade[parent].push_back(edge);
	}
	
	return edges_per_clade;
}



// Given a phylogenetic tree, create lookup tables listing incoming & outgoing edges for each clade (tip & node)
// Requirements:
//    The tree can include monofurcations and/or multifurcations
//    The tree need not be rooted
template<class ARRAY_TYPE>
void get_inout_edges_per_clade(	const long			Ntips,
								const long 			Nnodes,
								const long			Nedges,
								const ARRAY_TYPE	&tree_edge, 		// (INPUT) 2D array (in row-major format) of size Nedges x 2
								std::vector<long>	&clade2first_edge,	// (OUTPUT) 1D array of size Nnodes, with values in 0:(Nedges-1), mapping clades to their first incoming or outgoing edge.
								std::vector<long>	&clade2last_edge,	// (OUTPUT) 1D array of size Nnodes, with values in 0:(Nedges-1), mapping clades to their last incoming or outgoing edge.
								std::vector<long>	&edges){			// (OUTPUT) 1D array of size 2*Nedges, with values in 0:(Nedges-1). Maps internal edge indices (i.e. as listed in clade2first_edge[] and clade2last_edge[]) to original edge indices. Note that each edge is listed twice in this array (once as an outgoing and once as an incoming edge).
	const long Nclades = Ntips+Nnodes;
	edges.resize(2*Nedges);
	clade2first_edge.resize(Nclades);
	clade2last_edge.resize(Nclades);

	// determine number of edges per clade
	// edge_count_per_clade[n] will be the number of direct children of node n (n=0:(Nnodes-1))
	std::vector<long> edge_count_per_clade(Nclades, 0);
	for(long e=0; e<Nedges; ++e){
		edge_count_per_clade[tree_edge[e*2+0]] += 1;	// count edge towards parent
		edge_count_per_clade[tree_edge[e*2+1]] += 1;	// count edge towards child
	}
	// collect edges per clade
	clade2first_edge[0] = 0;
	clade2last_edge[0]  = clade2first_edge[0]+edge_count_per_clade[0] - 1;
	if(Nclades>1){
		for(long c=1; c<Nclades; ++c){
			clade2first_edge[c] = clade2last_edge[c-1]+1;
			clade2last_edge[c]  = clade2first_edge[c]+edge_count_per_clade[c] - 1;
		}
	}
	for(long e=0, clade; e<Nedges; ++e){
		// add edge to parent
		clade = tree_edge[e*2+0];
		edges[clade2first_edge[clade]+edge_count_per_clade[clade]-1] = e;
		edge_count_per_clade[clade] -= 1;

		// add edge to child
		clade = tree_edge[e*2+1];
		edges[clade2first_edge[clade]+edge_count_per_clade[clade]-1] = e;
		edge_count_per_clade[clade] -= 1;
	}
}



// Given a phylogenetic tree, create lookup tables listing incoming & outgoing edges for each node
// Requirements:
//    The tree can include monofurcations and/or multifurcations
//    The tree need not be rooted
template<class ARRAY_TYPE>
void get_inout_edges_per_node(	const long			Ntips,
								const long 			Nnodes,
								const long			Nedges,
								const ARRAY_TYPE	&tree_edge, 		// (INPUT) 2D array (in row-major format) of size Nedges x 2
								std::vector<long>	&node2first_edge,	// (OUTPUT) 1D array of size Nnodes, with values in 0:(Nedges-1), mapping nodes to their first incoming or outgoing edge.
								std::vector<long>	&node2last_edge,	// (OUTPUT) 1D array of size Nnodes, with values in 0:(Nedges-1), mapping nodes to their last incoming or outgoing edge.
								std::vector<long>	&edges){			// (OUTPUT) 1D array of size 2*Nedges, with values in 0:(Nedges-1). Maps internal edge indices (i.e. as listed in node2first_edge[] and node2last_edge[]) to original edge indices. Note that some edges will be listed twice in this array (once as an outgoing and once as an incoming edge).
	edges.resize(2*Nedges-Ntips);
	node2first_edge.resize(Nnodes);
	node2last_edge.resize(Nnodes);

	// determine number of edges per node
	// edge_count_per_node[n] will be the number of direct children of node n (n=0:(Nnodes-1))
	std::vector<long> edge_count_per_node(Nnodes, 0);
	for(long e=0; e<Nedges; ++e){
		if(tree_edge[e*2+0]>=Ntips) edge_count_per_node[tree_edge[e*2+0]-Ntips] += 1;	// count edge towards parent
		if(tree_edge[e*2+1]>=Ntips) edge_count_per_node[tree_edge[e*2+1]-Ntips] += 1;	// count edge towards child
	}
	// collect edges per clade
	node2first_edge[0] = 0;
	node2last_edge[0]  = node2first_edge[0]+edge_count_per_node[0] - 1;
	if(Nnodes>1){
		for(long n=1; n<Nnodes; ++n){
			node2first_edge[n] = node2last_edge[n-1]+1;
			node2last_edge[n]  = node2first_edge[n]+edge_count_per_node[n] - 1;
		}
	}
	for(long e=0, node; e<Nedges; ++e){
		// add edge to parent
		node = tree_edge[e*2+0]-Ntips;
		if(node>=0){
			edges[node2first_edge[node]+edge_count_per_node[node]-1] = e;
			edge_count_per_node[node] -= 1;
		}

		// add edge to child
		node = tree_edge[e*2+1]-Ntips;
		if(node>=0){
			edges[node2first_edge[node]+edge_count_per_node[node]-1] = e;
			edge_count_per_node[node] -= 1;
		}
	}
}



// calculate the path from the root to each of the tree's tips
// [[Rcpp::export]]
std::vector<std::vector<long> > get_paths_root_to_tips_CPP(	const long 				Ntips,
															const long 				Nnodes,
															const long 				Nedges,
															const std::vector<long> &tree_edge){
	// determine parent clade for each clade
	std::vector<long> clade2parent;
	get_parent_per_clade(Ntips, Nnodes, Nedges, tree_edge, clade2parent);
	
	// find root using the mapping clade2parent
	const long root = get_root_from_clade2parent(Ntips, clade2parent);
	
	std::vector<std::vector<long> > paths_to_tips(Ntips);
	for(long tip=0, clade; tip<Ntips; ++tip){
		paths_to_tips[tip].reserve(floor(2*log(Ntips)/log(2.0))); // rough estimate of typical tree depth x 2
		paths_to_tips[tip].push_back(tip);
		clade = tip;
		while(clade!=root){
			clade = clade2parent[clade];
			paths_to_tips[tip].push_back(clade);
		}
		// reverse direction (make root-->tip)
		reverse_array(paths_to_tips[tip]);
	}
	
	return paths_to_tips;
}


// determine edges adjacent to each edge (i.e. attached to each edge's child & parent)
// returns a list of size Nedges, each entry of which is a list of edge indices adjacent to the focal edge (the upstream adjacent edge is listed first)
// [[Rcpp::export]]
std::vector<std::vector<long> > get_adjacent_edges_per_edge_CPP(const long 				Ntips,
																const long 				Nnodes,
																const long 				Nedges,
																const std::vector<long> &tree_edge){
	// get incoming edge for each clade
	std::vector<long> incoming_edge_per_clade;
	get_incoming_edge_per_clade(Ntips, Nnodes, Nedges, tree_edge, incoming_edge_per_clade);
	
	std::vector<vector<long> > adjacent_edges_per_edge(Nedges);
	// record all upstream adjacents first
	for(long edge=0, parent, upstream_edge; edge<Nedges; ++edge){
		parent = tree_edge[edge*2+0];
		upstream_edge = incoming_edge_per_clade[parent];
		if(upstream_edge<0) continue;
		adjacent_edges_per_edge[edge].push_back(upstream_edge);
	}
	// record all downstream adjacents
	for(long edge=0, parent, upstream_edge; edge<Nedges; ++edge){
		parent = tree_edge[edge*2+0];
		upstream_edge = incoming_edge_per_clade[parent];
		if(upstream_edge<0) continue;
		adjacent_edges_per_edge[upstream_edge].push_back(edge);
	}
	
	return adjacent_edges_per_edge;
}	




#pragma mark -
#pragma mark Generating & manipulating trees
#pragma mark



// extract a subset of a tree
// the subset of clades (noder & tips) to keep is specified explicitly; edges are kept iff both their parent and child clade is kept
// If a node is kept but none of its descendants, then that node is turned into a tip
// This function assumes that the clades to be kept define a proper tree (and not a forest); no expansion is done.
// this function guarantees that the extracted tree follows the "phylo" conventions on indexing tips, nodes and the root
template<class INT_ARRAY>
void get_arbitrary_subtree(	const long				Ntips,				// (INPUT)
							const long 				Nnodes,				// (INPUT)
							const long				Nedges,				// (INPUT)
							const INT_ARRAY			&tree_edge,			// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
							const std::vector<char>	&keep_clade,		// (INPUT) 1D array of size Nclades, specifying whether a clade is to be kept or not
							long					&Ntips_new,			// (OUTPUT) Number of tips kept
							long					&Nnodes_new,		// (OUTPUT) Number of nodes kept
							long					&Nedges_new,		// (OUTPUT) Number of edges kept
							long					&new_root,			// (OUTPUT) new root clade index. Is actually guaranteed to be = Ntips_new.
							std::vector<long>		&new2old_clade,		// (OUTPUT) 1D array of size Nclades_new, mapping new clade indices to old indices
							std::vector<long>		&old2new_clade,		// (OUTPUT) 1D array of size Nclades, mapping old clade indices to new indices (or -1 if a clade is not kept)
							std::vector<long>		&new2old_edge,		// (OUTPUT) 1D array of size Nedges_new, mapping new edge indices to old indices
							std::vector<long>		&new_tree_edge){	// (OUTPUT) 2D array of size Nedges_new x 2, in row-major format, storing the edges of the extracted tree
	const long Nclades = Ntips+Nnodes;

	// determine number of tips & nodes & edges to keep
	Nedges_new = Ntips_new = Nnodes_new = 0;
	std::vector<bool> clade_becomes_tip(Nclades, true);
	for(long edge=0, parent, child; edge<Nedges; ++edge){
		parent 	= tree_edge[edge*2+0];
		child 	= tree_edge[edge*2+1];
		Nedges_new += (keep_clade[child] && keep_clade[parent] ? 1 : 0);
		if(keep_clade[child]) clade_becomes_tip[parent] = false;
	}
	for(long clade=0; clade<Nclades; ++clade){
		if(keep_clade[clade]){
			if(clade_becomes_tip[clade]) ++Ntips_new;
			else ++Nnodes_new;
		}
	}
	const long Nclades_new = Ntips_new + Nnodes_new;

	// create new2old clade mappings
	new_tree_edge.resize(2*Nedges_new);
	new2old_clade.resize(Nclades_new);
	new2old_edge.resize(Nedges_new);
	old2new_clade.assign(Nclades,-1);
	long next_tip_new  = 0;
	long next_node_new = 0;
	for(long clade=0, clade_new; clade<Nclades; ++clade){
		if(keep_clade[clade]){
			if(clade_becomes_tip[clade]) clade_new = next_tip_new++;
			else clade_new = Ntips_new + (next_node_new++);
			new2old_clade[clade_new] = clade;
			old2new_clade[clade] = clade_new;
		}
	}
	
	// determine new root (= new clade with no incoming kept edge)
	new_root = -1;
	std::vector<bool> new_clade_is_root(Nclades_new,true);
	for(long edge=0, parent, child; edge<Nedges; ++edge){
		parent 	= tree_edge[edge*2+0];
		child 	= tree_edge[edge*2+1];
		if(keep_clade[parent] && keep_clade[child]) new_clade_is_root[old2new_clade[child]] = false;
	}
	for(long clade_new=0; clade_new<Nclades_new; ++clade_new){
		if(new_clade_is_root[clade_new]){
			new_root = clade_new;
			break;
		}
	}
		
	// enforce common convention that new_root=Ntips_new 
	// (swap indices with clade previously mapped to Ntips_new)
	long new_root_in_old_tree				= new2old_clade[new_root];
	old2new_clade[new_root_in_old_tree] 	= Ntips_new;
	old2new_clade[new2old_clade[Ntips_new]] = new_root;
	new2old_clade[new_root] 				= new2old_clade[Ntips_new];
	new2old_clade[Ntips_new] 				= new_root_in_old_tree;
	new_root 								= Ntips_new;
	
	// create new2old edge mappings, and new edges
	long next_edge_new = 0;
	for(long edge=0, parent, child; edge<Nedges; ++edge){
		parent 	= tree_edge[edge*2+0];
		child 	= tree_edge[edge*2+1];
		if(keep_clade[child] && keep_clade[parent]){
			new2old_edge[next_edge_new] = edge;
			new_tree_edge[next_edge_new*2+0] = old2new_clade[parent];
			new_tree_edge[next_edge_new*2+1] = old2new_clade[child];
			++next_edge_new;
		}
	}
}



// Re-index clades in a tree such that tips are indexed 0,..,Ntips-1 and nodes are indexed Ntips,..,Nclades-1
// The tree is a-priori only defined based on tree_edge, but no convention of tip & node indexing is assumed
// The number of tips & nodes is inferred based on the tree topology (tree_edge)
// Optionally: The root of the tree (if existent) is ensured to have the index Ntips
void reindex_clades(	const long				Nclades,			// (INPUT) number of clades (tips or nodes) in the tree
						const long				Nedges,				// (INPUT) number of edges in the tree
						const std::vector<long>	&tree_edge,			// (INPUT) 2D array of size Nedges x 2, in row-major format
						const bool				root_convention,	// (INPUT) If true, the root of the tree (if existent) is ensured to obtain the index = Ntips
						long					&Ntips,				// (OUTPUT) the inferred number of tips in the tree
						long					&Nnodes,			// (OUTPUT) the inferred number of nodes in the tree
						std::vector<long>		&old2new_clade){	// (OUTPUT) 1D array of size Nclades, mapping old-->new clade indices
	// determine tips & nodes
	std::vector<bool> clade_is_tip(Nclades,true);
	for(long edge=0; edge<Nedges; ++edge){
		clade_is_tip[tree_edge[edge*2+0]] = false;
	}
	Ntips = Nnodes = 0;
	for(long clade=0; clade<Nclades; ++clade){
		if(clade_is_tip[clade]) ++Ntips;
		else ++Nnodes;
	}
		
	// re-index clades
	old2new_clade.resize(Nclades);
	long next_tip=0, next_node=0;
	for(long clade=0; clade<Nclades; ++clade){
		if(clade_is_tip[clade]) old2new_clade[clade] = (next_tip++);
		else old2new_clade[clade] = Ntips+(next_node++);
	}
		
	// make sure root is indexed = Ntips, if requested
	if(root_convention){
		std::vector<bool> clade_is_root(Nclades,true);
		for(long edge=0; edge<Nedges; ++edge){
			clade_is_root[tree_edge[edge*2+1]] = false;
		}
		long root = -1, occupier=-1;
		for(long clade=0; clade<Nclades; ++clade){
			if(clade_is_root[clade]) root = clade;
			if(old2new_clade[clade]==Ntips) occupier = clade; // clade currently re-indexed to Ntips
		}
		if(root>=0){
			// tree has a root, so correct index
			long temp_root_new 		= old2new_clade[root];
			old2new_clade[root] 	= Ntips;
			old2new_clade[occupier] = temp_root_new;
		}
	}
}



// sort tree edges such that they are in preorder (root-->tips) traversal
template<class ARRAY_TYPE>
void sort_tree_edges_root_to_tips(	const long			Ntips,
									const long			Nnodes,
									const long			Nedges,
									const bool			depth_first_search,	// (INPUT) if false, then breadth-first-seacrh is used instead
									const bool			root_to_tips,		// (INPUT) if false, then edges will be sorted tips-->root
									const ARRAY_TYPE	&tree_edge,			// (INPUT) 2D array of size Nedges x 2, in row-major format, listing edges
									std::vector<long>	&new2old_edge){		// (OUTPUT) 1D array of size Nedges, mapping new-->old edge index
	// get incoming edge for each clade
	std::vector<long> incoming_edge_per_clade;
	get_incoming_edge_per_clade(Ntips, Nnodes, Nedges, tree_edge, incoming_edge_per_clade);
		
	const long root = get_root_from_incoming_edge_per_clade(Ntips, tree_edge, incoming_edge_per_clade);

	// get node-->edge mappings
	std::vector<long> node2first_edge, node2last_edge, edges;
	get_node2edge_mappings(	Ntips,
							Nnodes,
							Nedges,
							tree_edge,
							node2first_edge,
							node2last_edge,
							edges);
	
	long child,node,clade;
	long next_new_edge = 0;
	new2old_edge.resize(Nedges);
	if(depth_first_search){
		// renumber edges in depth-first-search direction
		// use a scratch_stack for traversing nodes
		std::vector<long> scratch_stack;
		scratch_stack.reserve(floor(2*log(Ntips)/log(2.0))); // rough estimate of typical tree depth x 2. scratch_stack may be resized along the way if needed.
		scratch_stack.push_back(root);
		while(scratch_stack.size()>0){
			clade = scratch_stack.back();
			scratch_stack.pop_back();
			if(incoming_edge_per_clade[clade]>=0) new2old_edge[next_new_edge++] = incoming_edge_per_clade[clade];
			if(clade>=Ntips){
				node = clade-Ntips;
				for(long e=node2first_edge[node]; e<=node2last_edge[node]; ++e){
					scratch_stack.push_back(tree_edge[edges[e]*2+1]);
				}
			}
		}
		
	}else{
		// renumber edges in breadth-first-search direction
		// use a scratch_queue for traversing nodes
		std::vector<long> scratch_queue;
		scratch_queue.reserve(Nnodes);
		scratch_queue.push_back(root);
		long queue_pointer = 0;
		while(queue_pointer<scratch_queue.size()){
			node = scratch_queue[queue_pointer] - Ntips;
			queue_pointer += 1;
			for(long e=node2first_edge[node]; e<=node2last_edge[node]; ++e){
				new2old_edge[next_new_edge++] = edges[e];
				child = tree_edge[edges[e]*2+1];
				if(child>=Ntips) scratch_queue.push_back(child);	// add child to queue for further exploration in the next iteration
			}
		}
	}
	
	if(!root_to_tips) reverse_array(new2old_edge);
}




// sort tree edges such that parent nodes are in preorder (root-->tips) traversal
// Returns an integer array new2old_edge mapping new-->old edge indices
// This is an Rcpp wrapper for sort_tree_edges_root_to_tips(..)
// [[Rcpp::export]]
IntegerVector sort_tree_edges_root_to_tips_CPP(	const long			Ntips,
												const long 			Nnodes,
												const long			Nedges,
												const bool			depth_first_search,	// (INPUT) if false, then breadth-first-seacrh is used instead.
												const bool			root_to_tips,		// (INPUT) if false, then edges will be sorted tips-->root
												const IntegerVector	&tree_edge){		// (INPUT) 2D array of size Nedges x 2, in row-major format, listing edges.
	std::vector<long> new2old_edge;
	sort_tree_edges_root_to_tips(	Ntips,
									Nnodes,
									Nedges,
									depth_first_search,
									root_to_tips,
									tree_edge,
									new2old_edge);
	return Rcpp::wrap(new2old_edge);
}







// Adjust edge directions so that a specific node becomes the root of the tree (i.e. has no parent and all other tips & nodes descend from it)
// Note that node & tip indices remain unchanged, which may break the common convention that node #1 is the root.
// This function operates in-situ
template<class ARRAY_TYPE>
void root_tree_at_node(	const long 	Ntips,
						const long 	Nnodes,
						const long 	Nedges,
						ARRAY_TYPE 	&tree_edge,		// (INPUT/OUTPUT) 2D array of size Nedges x 2 in row-major format, with values in 0,..,(Nclades-1). Will be modified in-situ
						const long	new_root_node){	// (INPUT) index of node to be turned into root. Must be within 1,..,Nnodes.
	const long Nclades = Ntips+Nnodes;
	
	// create lookup tables listing all edges (regardless of direction) per clade
	std::vector<long> clade2first_edge, clade2last_edge, edges;
	get_inout_edges_per_clade(	Ntips,
								Nnodes,
								Nedges,
								tree_edge,
								clade2first_edge,
								clade2last_edge,
								edges);
						
	// determine correct edge direction and keep track of which edge is incoming for each clade
	// use a scratch_stack for traversing nodes in depth-first-search order
	const long root = new_root_node+Ntips;
	std::vector<long> incoming_edge_per_clade(Nclades,-1);
	long child,clade;
	std::vector<long> scratch_stack;
	scratch_stack.reserve(floor(2*log(Ntips)/log(2.0))); // rough estimate of typical tree depth x 2. scratch_stack may be resized along the way if needed.
	scratch_stack.push_back(root);
	while(scratch_stack.size()>0){
		clade = scratch_stack.back();
		scratch_stack.pop_back();
		for(long e=clade2first_edge[clade], edge; e<=clade2last_edge[clade]; ++e){
			edge = edges[e];
			if(edge==incoming_edge_per_clade[clade]) continue; // skip incoming edge
			// at this point it is guaranteed that edge should be outgoing from clade, so adjust direction if needed
			if(tree_edge[edge*2+0]!=clade){
				// edge is in wrong direction, so reverse
				tree_edge[edge*2+1] = tree_edge[edge*2+0];
				tree_edge[edge*2+0] = clade;
			}
			// place child on stack (if a node) and note its incoming edge
			child = tree_edge[edge*2+1];
			incoming_edge_per_clade[child] = edge;
			if(child>=Ntips) scratch_stack.push_back(child);
		}
	}
}


// Adjust edge directions so that a specific node becomes the root of the tree (i.e. has no parent and all other tips & nodes descend from it)
// This function returns a 2D array of size Nedges x 2, in row-major format, listing the new tree edges.
// This function is an Rcpp wrapper for the function root_tree_at_node(..)
// [[Rcpp::export]]
IntegerVector root_tree_at_node_CPP(const long 			Ntips,
									const long 			Nnodes,
									const long 			Nedges,
									const IntegerVector &tree_edge,		// (INPUT) 2D array of size Nedges x 2 in row-major format, with values in 0,..,(Nclades-1)
									const long			new_root_node){	// (INPUT) index of node to be turned into root. Must be within 1,..,Nnodes.
	IntegerVector new_tree_edge = tree_edge;
	root_tree_at_node(	Ntips,
						Nnodes,
						Nedges,
						new_tree_edge,
						new_root_node);
	return new_tree_edge;
}
	





// collapse monofurcations in a rooted phylogenetic tree
// i.e. whenever a node only has one child, then combine its 2 adjacent edges into one and remove that node
// note that the tips are always kept
template<class ARRAY_INT, class ARRAY_DOUBLE>
void get_tree_with_collapsed_monofurcations(const long 				Ntips,
											const long 				Nnodes,
											const long 				Nedges,
											const long 				root,
											const bool				force_keep_root,	// (INPUT) if true, the root is always kept even if it only has one child
											const ARRAY_INT			&tree_edge,			// (INPUT) 2D array of size Nedges x 2 (in row-major format)
											const ARRAY_DOUBLE 		&edge_length, 		// (INPUT) 1D array of size Nedges, or an empty array (all branches have length 1)
											std::vector<long>		&new_tree_edge,		// (OUTPUT) 2D matrix (in row major format)
											std::vector<double>		&new_edge_length,	// (OUTPUT)
											std::vector<long>		&new2old_node,		// (OUTPUT)
											long					&new_root,					// (OUTPUT) index of the root in the collapsed tree. In newer implementations, this is actually guaranteed to be = Ntips (as is common convention).
											double					&root_shift){	// (OUTPUT) phylogenetic distance of new root from the old root (will be 0 if the old root is kept)
	const long Nclades = Ntips+Nnodes;
	long node, edge, incoming_edge, child, parent;
	
	// get tree traversal route (root --> tips)											
	std::vector<long> traversal_queue, traversal_node2first_edge, traversal_node2last_edge, traversal_edges;
	get_tree_traversal_root_to_tips(	Ntips,
										Nnodes,
										Nedges,
										root,
										tree_edge,
										false,
										false,
										traversal_queue,
										traversal_node2first_edge,
										traversal_node2last_edge,
										traversal_edges,
										false,
										"");
										
	// determine which nodes to keep
	long Nnodes_kept = 0;
	std::vector<bool> keep_node(Nnodes);
	for(long node=0; node<Nnodes; ++node){
		keep_node[node] = ((Ntips+node)==root && force_keep_root) || (traversal_node2first_edge[node]<traversal_node2last_edge[node]); // only keep node if it has at least 2 edges leaving from it
		if(keep_node[node]) ++Nnodes_kept;
	}
	
	// calculate the mappings new node <--> old node
	// traverse root-->tips, so that the old root is mapped to the lowest new node index (0)
	new2old_node.clear();
	new2old_node.reserve(Nnodes_kept);
	std::vector<long> old2new_node(Nnodes,-1);
	long Nedges_kept = Nedges;
	for(long q=0; q<traversal_queue.size(); ++q){
		node = traversal_queue[q]-Ntips;
		if(keep_node[node]){
			new2old_node.push_back(node);
			old2new_node[node] = new2old_node.size()-1;
		}else{ 
			Nedges_kept -= (traversal_node2last_edge[node]-traversal_node2first_edge[node]+1); // edges leaving this node will not be kept
		}
	}
	
	// determine incoming edge per clade
	std::vector<long> incoming_edge_per_clade(Nclades,-1);
	for(edge=0; edge<Nedges; ++edge){
		incoming_edge_per_clade[tree_edge[edge*2+1]] = edge;
	}
	
	// collapse edges (traverse root --> tips)
	// note that at first, tree_edge[,] will list old clade indices (will be renamed later)
	new_tree_edge.resize(Nedges_kept*2);
	new_edge_length.resize(Nedges_kept);
	std::vector<long> new_incoming_edge_per_clade = incoming_edge_per_clade; // make a copy, then modify in-situ
	long next_new_edge = 0;
	std::vector<long> old2new_edge(Nedges,-1);
	for(long q=0; q<traversal_queue.size(); ++q){
		parent = traversal_queue[q];
		node = parent-Ntips;
		for(long ei=traversal_node2first_edge[node]; ei<=traversal_node2last_edge[node]; ++ei){
			edge  = traversal_edges[ei];
			child = tree_edge[edge*2+1];
			if(keep_node[node]){
				new_tree_edge[next_new_edge*2+0] = parent;
				new_tree_edge[next_new_edge*2+1] = child;
				new_edge_length[next_new_edge] = (edge_length.size()==0 ? 1 : edge_length[edge]);
				old2new_edge[edge] = next_new_edge;
				++next_new_edge;
			}else{
				// merge this departing edge with the incoming edge
				// note that this should only occur once per node (we're not removing nodes with more than 1 departing edges)
				incoming_edge = new_incoming_edge_per_clade[parent];
				if(incoming_edge>=0){
					new_tree_edge[old2new_edge[incoming_edge]*2+1] 	= child;
					new_edge_length[old2new_edge[incoming_edge]] 	+= (edge_length.size()==0 ? 1 : edge_length[edge]); // append this edge's length to the incoming edge
				}
				new_incoming_edge_per_clade[child] = incoming_edge; // update incoming edge of child, since this child was assigned as the destination of incoming_edge. If incoming_edge was -1 (i.e. clade was root), then child becomes the de-facto new root
			}
		}
	}

	// rename nodes in new_tree_edge[]
	for(long new_edge=0; new_edge<Nedges_kept; ++new_edge){
		parent = new_tree_edge[new_edge*2+0];
		child  = new_tree_edge[new_edge*2+1];
		new_tree_edge[new_edge*2+0] = Ntips + old2new_node[parent-Ntips];
		if(child>=Ntips) new_tree_edge[new_edge*2+1] = Ntips + old2new_node[child-Ntips];
	}
	
	// determine new root (new node with no parent)
	root_shift = 0;
	if(force_keep_root){
		new_root = Ntips + old2new_node[root-Ntips];
	}else{
		new_root = get_root_clade(Ntips, Nnodes_kept, Nedges_kept, new_tree_edge);
		// traverse from new to old root and calculate cumulative distance
		// use old tree structure
		long clade = Ntips+new2old_node[new_root-Ntips];
		while(incoming_edge_per_clade[clade]>=0){
			edge = incoming_edge_per_clade[clade];
			root_shift += (edge_length.size()==0 ? 1 : edge_length[edge]);
			clade = tree_edge[edge*2 + 0];
		}
	}
}


// Rcpp wrapper for get_tree_with_collapsed_monofurcations()
// [[Rcpp::export]]
Rcpp::List get_tree_with_collapsed_monofurcations_CPP(	const long 					Ntips,
														const long 					Nnodes,
														const long 					Nedges,
														const std::vector<long>		&tree_edge,			// (INPUT) 2D array of size Nedges x 2 (in row-major format)
														const std::vector<double> 	&edge_length, 		// (INPUT) 1D array of size Nedges, or an empty array (all branches have length 1)
														const bool					force_keep_root){	// (INPUT) if true, the root is always kept even if it only has one child
	// find root
	const long root = get_root_clade(Ntips, Nnodes, Nedges, tree_edge);

	std::vector<long> new_tree_edge, new2old_node;
	std::vector<double> new_edge_length;
	long new_root;
	double root_shift;
	get_tree_with_collapsed_monofurcations(	Ntips,
											Nnodes,
											Nedges,
											root,
											force_keep_root,
											tree_edge,	
											edge_length,
											new_tree_edge,
											new_edge_length,
											new2old_node,
											new_root,
											root_shift);

	return Rcpp::List::create(	Rcpp::Named("Nnodes_new")		= long(new2old_node.size()), // explicitly cast to long, otherwise Rcpp does not know how to wrap it (on Windows)
								Rcpp::Named("new_tree_edge") 	= Rcpp::wrap(new_tree_edge),
								Rcpp::Named("new_edge_length") 	= Rcpp::wrap(new_edge_length), // if the original edge_length[] was empty, then new_edge_length[e] will be the number of combined edges making up the new edge e
								Rcpp::Named("new2old_node") 	= Rcpp::wrap(new2old_node),
								Rcpp::Named("new_root") 		= new_root,
								Rcpp::Named("root_shift") 		= root_shift);
}






// Extract subtree with subset of tips and/or nodes
// Requirements:
//   tree can include multifucations and monofurcations
//   tree must be rooted (root will be determined automatically, as the node without a parent)
// [[Rcpp::export]]
Rcpp::List get_subtree_with_specific_clades_CPP(const long 			Ntips,
												const long 			Nnodes,
												const long 			Nedges,
												const IntegerVector &tree_edge,			// 2D array of size Nedges x 2 in row-major format
												const NumericVector &edge_length, 		// 1D array of size Nedges, or an empty std::vector (all branches have length 1)
												const IntegerVector &clades_to_keep,	// 1D array with values in 0,..,(Nclades-1)
												const bool			collapse_monofurcations,	// if true, nodes that are left with only one child (after pruning) will be removed (and the adjacent edges will be combined into a single edge)
												const bool			force_keep_root,	// if true, then the root is kept even if collapse_monofurcations==true and the root is monofurcating.
												bool				keep_all_children_of_explicit_clades_to_keep,	// if true, then for each clade in clades_to_keep[], all children are also kept. Otherwise only one child is kept.
												bool				keep_all_tips_of_explicit_clades_to_keep){		// if true, then for each clade in clades_to_keep[], all descending tips are also kept. Otherwise (typically) only one tip is kept per child.
	const long Nclades = Ntips+Nnodes;
	long clade, parent, child, node;
	if(keep_all_tips_of_explicit_clades_to_keep) keep_all_children_of_explicit_clades_to_keep = false;

	// get parent of each clade
	std::vector<long> clade2parent;
	get_parent_per_clade(Ntips, Nnodes, Nedges, tree_edge, clade2parent);

	// find root using the mapping clade2parent
	const long root = get_root_from_clade2parent(Ntips, clade2parent);	
											
	// get tree traversal route (tips --> root)											
	std::vector<long> traversal_queue_root2tips, traversal_node2first_edge, traversal_node2last_edge, traversal_edges;
	get_tree_traversal_root_to_tips(Ntips,
									Nnodes,
									Nedges,
									root,
									tree_edge,
									true,
									false,
									traversal_queue_root2tips,
									traversal_node2first_edge,
									traversal_node2last_edge,
									traversal_edges,
									false,
									"");
			
	// determine which clades to keep
	std::vector<bool> keep_clade_explicitly(Nclades,false);
	std::vector<bool> keep_clade_ancestrally(Nclades,false);
	std::vector<bool> keep_clade_descending(Nclades,false);
	for(long c=0; c<clades_to_keep.size(); ++c){
		clade = clades_to_keep[c];
		keep_clade_explicitly[clade]  = true;
		keep_clade_ancestrally[clade] = true;
		keep_clade_descending[clade]  = true;
	}
	// step 1: propagate status of to-be-kept upwards (traverse tips-->root)
	// so that all ancestors of a clade to be kept are also kept
	for(long q=traversal_queue_root2tips.size()-1; q>=0; --q){
		clade = traversal_queue_root2tips[q];
		if(keep_clade_ancestrally[clade] && (clade!=root)){
			// clade is to be kept, so propagate upwards
			keep_clade_ancestrally[clade2parent[clade]] = true;
		}
	}
	// step 2: propagate status of to-be-kept downwards (traverse root-->tips)
	// so that each clade to be kept leads to at least one (or multiple) tips
	for(long q=0; q<traversal_queue_root2tips.size(); ++q){
		clade = traversal_queue_root2tips[q];
		if(keep_clade_descending[clade] && (clade>=Ntips)){
			node = clade - Ntips;
			const bool keep_all_children = (keep_all_tips_of_explicit_clades_to_keep || (keep_clade_explicitly[clade] && keep_all_children_of_explicit_clades_to_keep));
			if((!keep_all_children) && (!keep_clade_explicitly[clade]) && keep_clade_ancestrally[clade]) continue; // clade has already been included due to one of its descendants, so no need to include more of its children
			for(long e=traversal_node2first_edge[node]; e<=(keep_all_children ? traversal_node2last_edge[node] : traversal_node2first_edge[node]); ++e){
				child = tree_edge[traversal_edges[e]*2+1];
				keep_clade_descending[child] = true;
			}
		}
	}
	// step 3: merge various reasons of clade inclusion
	std::vector<bool> keep_clade(Nclades,false);
	long Nedges_kept = 0;
	for(long c=0; c<Nclades; ++c){
		keep_clade[c] = (keep_clade_explicitly[c] || keep_clade_ancestrally[c] || keep_clade_descending[c]);
		if(keep_clade[c] && (c!=root)) ++Nedges_kept; // for each kept clade, there corresponds exactly one kept edge (unless the clade is root)
	}

		
	// translate old clade indices to new (map to -1 if not included)
	long Nclades_kept = 0;
	long Ntips_kept   = 0;
	std::vector<long> old2new_clade(Nclades,-1);
	for(long c=0; c<Nclades; ++c){
		if(keep_clade[c]){
			old2new_clade[c] = (Nclades_kept++);
			if(c<Ntips) Ntips_kept += 1;
		}
	}
	long new_root 	 = old2new_clade[root];
	long Nnodes_kept = Nclades_kept-Ntips_kept;
	
	// calculate the reverse mapping (new --> old clade index)
	std::vector<long> new2old_clade(Nclades_kept,-1);
	for(long c=0; c<Nclades; ++c){
		if(keep_clade[c]) new2old_clade[old2new_clade[c]] = c;
	}
	
	// enforce common convention that new_root=Ntips_kept (swap indices with clade previously mapped to Ntips_kept)
	old2new_clade[root] 						= Ntips_kept;
	old2new_clade[new2old_clade[Ntips_kept]] 	= new_root;
	new2old_clade[new_root] 					= new2old_clade[Ntips_kept];
	new2old_clade[Ntips_kept] 					= root;
	new_root 									= Ntips_kept;
	
	// extract subset of kept edges
	IntegerVector new_tree_edge(Nedges_kept*2);
	NumericVector new_edge_length(Nedges_kept);
	long next_new_edge = 0;
	for(long edge=0; edge<Nedges; ++edge){
		parent = tree_edge[edge*2+0];
		child  = tree_edge[edge*2+1];
		if(keep_clade[parent] && keep_clade[child]){
			new_tree_edge[next_new_edge*2+0] = old2new_clade[parent];
			new_tree_edge[next_new_edge*2+1] = old2new_clade[child];
			new_edge_length[next_new_edge] = (edge_length.size()==0 ? 1 : edge_length[edge]);
			++next_new_edge;
		}
	}
	
	
	// collapse monofurcations if needed
	double root_shift = 0;
	if(collapse_monofurcations){
		long newer_root;
		std::vector<long> newer_tree_edge, newer2new_node;
		std::vector<double> newer_edge_length;
		get_tree_with_collapsed_monofurcations(	Ntips_kept, 
												(Nclades_kept-Ntips_kept), 
												Nedges_kept,
												new_root,
												force_keep_root,
												new_tree_edge,
												new_edge_length,
												newer_tree_edge,
												newer_edge_length,
												newer2new_node,
												newer_root,
												root_shift);
		new_tree_edge 	= Rcpp::wrap(newer_tree_edge);
		new_edge_length = Rcpp::wrap(newer_edge_length);
		Nnodes_kept		= newer2new_node.size();
		Nclades_kept 	= Ntips_kept+Nnodes_kept;
		std::vector<long> newer2old_node(Nnodes_kept);
		for(long n=0; n<Nnodes_kept; ++n){
			newer2old_node[n] = new2old_clade[Ntips_kept+newer2new_node[n]]-Ntips;
		}
		new2old_clade.resize(Nclades_kept);
		for(long n=0; n<Nnodes_kept; ++n){
			new2old_clade[Ntips_kept+n] = Ntips + newer2old_node[n];
		}
		new_root = newer_root;
	}
	
	return Rcpp::List::create(	Rcpp::Named("new_tree_edge") 	= Rcpp::wrap(new_tree_edge),
								Rcpp::Named("new_edge_length") 	= Rcpp::wrap(new_edge_length), // if the original edge_length[] was empty, then new_edge_length[e] will be the number of combined edges making up the new edge e
								Rcpp::Named("new2old_clade") 	= Rcpp::wrap(new2old_clade),
								Rcpp::Named("new_root") 		= new_root,
								Rcpp::Named("Ntips_kept") 		= Ntips_kept,
								Rcpp::Named("Nnodes_kept") 		= Nnodes_kept,
								Rcpp::Named("Nedges_kept") 		= Nedges_kept,
								Rcpp::Named("root_shift") 		= root_shift);
}







// extract the subtree descending from a specific node
// The new root is guaranteed to be the node with index = 0.
// [[Rcpp::export]]
Rcpp::List get_subtree_at_node_CPP(	const long 			Ntips,
									const long 			Nnodes,
									const long 			Nedges,
									const IntegerVector &tree_edge,		// (INPUT) 2D array of size Nedges x 2 in row-major format
									const long			new_root_node){	// (INPUT) node at which to extract subtree. This node will become the root of the extracted subtree. Must be within 1,..,Nnodes
	const long Nclades = Ntips + Nnodes;
	
	// get edge mapping for easier traversal
	std::vector<long> node2first_edge, node2last_edge, edges;
	get_node2edge_mappings(	Ntips,
							Nnodes,
							Nedges,
							tree_edge,
							node2first_edge,
							node2last_edge,
							edges);
							

	// determine size of subtree (count descendants of new_root_node) and generate old2new mappings
	// use a scratch_stack for traversing nodes in a depth-first-search manner
	long child,node,edge;
	std::vector<long> old2new_node(Nnodes,-1), old2new_tip(Ntips,-1), old2new_edge(Nedges,-1);
	std::vector<long> scratch_stack;
	scratch_stack.reserve(floor(2*log(Ntips)/log(2.0))); // rough estimate of typical tree depth x 2. scratch_stack may be resized along the way if needed.
	scratch_stack.push_back(new_root_node+Ntips);
	long Ntips_kept  = 0;
	long Nnodes_kept = 1;
	long Nedges_kept = 0;
	old2new_node[new_root_node] = Nnodes_kept-1;
	while(scratch_stack.size()>0){
		node = scratch_stack.back() - Ntips;
		scratch_stack.pop_back();
		for(long e=node2first_edge[node]; e<=node2last_edge[node]; ++e){
			edge  = edges[e];
			child = tree_edge[edge*2+1];
			if(child>=Ntips) scratch_stack.push_back(child); // add child to stack for further exploration in the next iteration
			// increment counters
			++Nedges_kept;
			if(child<Ntips){
				++Ntips_kept;
				old2new_tip[child] = Ntips_kept-1;
			}else{
				++Nnodes_kept;
				old2new_node[child-Ntips] = Nnodes_kept-1;
			}
			old2new_edge[edge] = Nedges_kept-1;
		}
	}
	const long Nclades_kept = Ntips_kept + Nnodes_kept;
	
					
	// Extract subtree and generate new2old mappings
	std::vector<long> new_tree_edge(2*Nedges_kept), new2old_clade(Nclades_kept), new2old_edge(Nedges_kept);
	for(long edge=0, new_edge, clade; edge<Nedges; ++edge){
		new_edge = old2new_edge[edge];
		if(new_edge<0) continue; // this edge is not to be kept
		new2old_edge[new_edge] 		= edge;
		clade = tree_edge[edge*2+0];
		new_tree_edge[new_edge*2+0] = (clade<Ntips ? old2new_tip[clade] : old2new_node[clade-Ntips]+Ntips_kept);
		clade = tree_edge[edge*2+1];
		new_tree_edge[new_edge*2+1] = (clade<Ntips ? old2new_tip[clade] : old2new_node[clade-Ntips]+Ntips_kept);
	}
	for(long clade=0, new_clade, new_node; clade<Nclades; ++clade){
		if(clade<Ntips){
			new_clade = old2new_tip[clade];
			if(new_clade<0) continue;
		}else{
			new_node = old2new_node[clade-Ntips];
			if(new_node<0) continue;
			new_clade = new_node+Ntips_kept;
		}
		new2old_clade[new_clade] = clade;
	}
	
	return Rcpp::List::create(	Rcpp::Named("new_tree_edge") 	= Rcpp::wrap(new_tree_edge),
								Rcpp::Named("new2old_clade") 	= Rcpp::wrap(new2old_clade),
								Rcpp::Named("new2old_edge") 	= Rcpp::wrap(new2old_edge),
								Rcpp::Named("new_root") 		= old2new_node[new_root_node]+Ntips_kept,
								Rcpp::Named("Ntips_kept") 		= Ntips_kept,
								Rcpp::Named("Nnodes_kept") 		= Nnodes_kept,
								Rcpp::Named("Nedges_kept") 		= Nedges_kept);
}







// Extract subtree with a specific subset of tips (and including all of their ancestors)
// If the original edge_length[] was empty, then new_edge_length[e] will be the number of combined edges making up the new edge e
// This function guarantees that the new root of the subtree will have index = Ntips_kept.
// Requirements:
//   tree can include multifucations and monofurcations
//   tree must be rooted (root will be determined automatically, as the node without a parent)
template<class ARRAY_INT,class ARRAY_DOUBLE>
void get_subtree_with_specific_tips(const long 			Ntips,
									const long 			Nnodes,
									const long 			Nedges,
									const ARRAY_INT 	&tree_edge,						// (INPUT) 2D array of size Nedges x 2 in row-major format
									const ARRAY_DOUBLE 	&edge_length, 					// (INPUT) 1D array of size Nedges, or an empty std::vector (all branches have length 1)
									const ARRAY_INT		&tips_to_keep,					// (INPUT) 1D array with values in 0,..,(Ntips-1)
									bool				collapse_monofurcations,		// (INPUT) if true, nodes that are left with only one child (after pruning) will be removed (and the adjacent edges will be combined into a single edge)
									bool				force_keep_root,				// (INPUT) alwasy keep root, even if collapse_monofurcations==true and root is monofurcating
									std::vector<long>	&new_tree_edge,					// (OUTPUT) 2D array of size Nedges x 2, in row-major format
									std::vector<double>	&new_edge_length,				// (OUTPUT) 1D array of size Nedges
									std::vector<long>	&new2old_clade,					// (OUTPUT) 1D array of size Nclades_kept
									long				&new_root,						// (OUTPUT) root index in the new tree. In newer implementations this is actually guaranteed to be Ntips_kept+1.
									long				&Ntips_kept,					// (OUTPUT)
									long				&Nnodes_kept,					// (OUTPUT)
									long				&Nedges_kept,					// (OUTPUT)
									double				&root_shift){					// (OUTPUT)
	const long Nclades = Ntips+Nnodes;
	long clade, parent, child;
	
	// get parent of each clade
	std::vector<long> clade2parent;
	get_parent_per_clade(Ntips, Nnodes, Nedges, tree_edge, clade2parent);

	// find root using the mapping clade2parent
	const long root = get_root_from_clade2parent(Ntips, clade2parent);
								
	// get tree traversal route
	std::vector<long> traversal_queue_root2tips, traversal_node2first_edge, traversal_node2last_edge, traversal_edges;
	get_tree_traversal_root_to_tips(Ntips,
									Nnodes,
									Nedges,
									root,
									tree_edge,
									true,
									false,
									traversal_queue_root2tips,
									traversal_node2first_edge,
									traversal_node2last_edge,
									traversal_edges,
									false,
									"");
		
	// determine which clades to keep
	Nedges_kept = 0;
	std::vector<bool> keep_clade(Nclades,false);
	for(long t=0; t<tips_to_keep.size(); ++t){
		keep_clade[tips_to_keep[t]] = true;
	}
	// propagate status of to-be-kept upwards (traverse tips-->root)
	// so that all ancestors of a clade to be kept are also kept
	for(long q=traversal_queue_root2tips.size()-1; q>=0; --q){
		clade = traversal_queue_root2tips[q];
		if(keep_clade[clade] && (clade!=root)){
			// clade is to be kept, so propagate upwards
			keep_clade[clade2parent[clade]] = true;
			++Nedges_kept;
		}
	}
		
	// translate old clade indices to new (map to -1 if not included)
	long Nclades_kept = 0;
	Ntips_kept = 0;
	std::vector<long> old2new_clade(Nclades,-1);
	for(long c=0; c<Nclades; ++c){
		if(keep_clade[c]){
			old2new_clade[c] = (Nclades_kept++);
			if(c<Ntips) Ntips_kept += 1;
		}
	}
	new_root 	= old2new_clade[root];
	Nnodes_kept = Nclades_kept - Ntips_kept;
	
	// calculate the reverse mapping (new --> old clade index)
	new2old_clade.assign(Nclades_kept,-1);
	for(long c=0; c<Nclades; ++c){
		if(keep_clade[c]) new2old_clade[old2new_clade[c]] = c;
	}
	
	// enforce common convention that new_root=Ntips_kept (swap indices with clade previously mapped to Ntips_kept)
	old2new_clade[root] 						= Ntips_kept;
	old2new_clade[new2old_clade[Ntips_kept]] 	= new_root;
	new2old_clade[new_root] 					= new2old_clade[Ntips_kept];
	new2old_clade[Ntips_kept] 					= root;
	new_root 									= Ntips_kept;
	
	// extract subset of kept edges
	new_tree_edge.resize(Nedges_kept*2);
	new_edge_length.resize(Nedges_kept);
	long next_new_edge = 0;
	for(long edge=0; edge<Nedges; ++edge){
		parent = tree_edge[edge*2+0];
		child  = tree_edge[edge*2+1];
		if(keep_clade[parent] && keep_clade[child]){
			new_tree_edge[next_new_edge*2+0] = old2new_clade[parent];
			new_tree_edge[next_new_edge*2+1] = old2new_clade[child];
			new_edge_length[next_new_edge] = (edge_length.size()==0 ? 1 : edge_length[edge]);
			++next_new_edge;
		}
	}
	
	// collapse monofurcations if needed
	root_shift = 0;
	if(collapse_monofurcations){
		long newer_root;
		std::vector<long> newer_tree_edge, newer2new_node;
		std::vector<double> newer_edge_length;
		get_tree_with_collapsed_monofurcations(	Ntips_kept, 
												Nnodes_kept, 
												Nedges_kept,
												new_root,
												force_keep_root,
												new_tree_edge,
												new_edge_length,
												newer_tree_edge,
												newer_edge_length,
												newer2new_node,
												newer_root,
												root_shift);
		new_tree_edge 	= newer_tree_edge;
		new_edge_length = newer_edge_length;
		Nnodes_kept		= newer2new_node.size();
		Nclades_kept 	= Ntips_kept+Nnodes_kept;
		std::vector<long> newer2old_node(Nnodes_kept);
		for(long n=0; n<Nnodes_kept; ++n){
			newer2old_node[n] = new2old_clade[Ntips_kept+newer2new_node[n]]-Ntips;
		}
		new2old_clade.resize(Nclades_kept);
		for(long n=0; n<Nnodes_kept; ++n){
			new2old_clade[Ntips_kept+n] = Ntips + newer2old_node[n];
		}
		new_root = newer_root;
	}
}




// Extract subtree with a specific subset of tips (and including all of their ancestors)
// This is an Rcpp wrapper for get_subtree_with_specific_tips()
// Requirements:
//   tree can include multifucations and monofurcations
//   tree must be rooted (root will be determined automatically, as the node without a parent)
// [[Rcpp::export]]
Rcpp::List get_subtree_with_specific_tips_CPP(	const long 			Ntips,
												const long 			Nnodes,
												const long 			Nedges,
												const IntegerVector &tree_edge,			// (INPUT) 2D array of size Nedges x 2 in row-major format
												const NumericVector &edge_length, 		// (INPUT) 1D array of size Nedges, or an empty std::vector (all branches have length 1)
												const IntegerVector &tips_to_keep,		// (INPUT) 1D array with values in 0,..,(Ntips-1)
												bool				collapse_monofurcations,	// (INPUT) if true, nodes that are left with only one child (after pruning) will be removed (and the adjacent edges will be combined into a single edge)
												bool				force_keep_root){	// (INPUT) alwasy keep root, even if collapse_monofurcations==true and root is monofurcating

	std::vector<long> new_tree_edge, new2old_clade;
	std::vector<double> new_edge_length;
	long new_root, Ntips_kept, Nnodes_kept, Nedges_kept;
	double root_shift;
	get_subtree_with_specific_tips(	Ntips,
									Nnodes,
									Nedges,
									tree_edge,
									edge_length,
									tips_to_keep,
									collapse_monofurcations,
									force_keep_root,
									new_tree_edge,
									new_edge_length,
									new2old_clade,
									new_root,
									Ntips_kept,
									Nnodes_kept,
									Nedges_kept,
									root_shift);
	
	return Rcpp::List::create(	Rcpp::Named("new_tree_edge") 	= Rcpp::wrap(new_tree_edge),
								Rcpp::Named("new_edge_length") 	= Rcpp::wrap(new_edge_length),
								Rcpp::Named("new2old_clade") 	= Rcpp::wrap(new2old_clade),
								Rcpp::Named("new_root") 		= new_root,
								Rcpp::Named("Ntips_kept") 		= Ntips_kept,
								Rcpp::Named("Nnodes_kept") 		= Nnodes_kept,
								Rcpp::Named("Nedges_kept") 		= Nedges_kept,
								Rcpp::Named("root_shift")		= root_shift);
}




// returns a random tip in the tree, by randomly traversing the tree root-->tips.
// each child of a node is equally probable to be traversed (when passed through that node)
long get_tip_by_random_uniform_traversal(	const long			Ntips,
											const long			root,				// (INPUT) root index. Typically root=Ntips.
											std::vector<long>	&node2first_child,	// (INPUT) 1D array of size Nnodes, with values in 0:(Nedges-1). Should be in the format returned by get_children_per_node(..).
											std::vector<long>	&node2last_child,	// (INPUT) 1D array of size Nnodes, with values in 0:(Nedges-1). Should be in the format returned by get_children_per_node(..).
											std::vector<long>	&children){			// (INPUT) 1D array of size Nedges, with values in 0:(Nclades-1). Should be in the format returned by get_children_per_node(..).
	long clade = root;
	while(clade>=Ntips){
		clade = children[uniformIntWithin(node2first_child[clade-Ntips],node2last_child[clade-Ntips])];
	}
	return clade;
}



// returns a random tip in the tree, by randomly traversing the tree root-->tips.
// The probability of a child of a node is given externally (via clade2probability[]). Note t
long get_tip_by_random_traversal(	const long					Ntips,
									const long					root,					// (INPUT) root index. Typically root=Ntips.
									const std::vector<long>		&node2first_child,		// (INPUT) 1D array of size Nnodes, with values in 0:(Nedges-1). Should be in the format returned by get_children_per_node(..).
									const std::vector<long>		&node2last_child,		// (INPUT) 1D array of size Nnodes, with values in 0:(Nedges-1). Should be in the format returned by get_children_per_node(..).
									const std::vector<long>		&children,				// (INPUT) 1D array of size Nedges, with values in 0:(Nclades-1). Should be in the format returned by get_children_per_node(..).
									const std::vector<double>	&clade2weight){			// (INPUT) 1D array of non-negative numbers, of size Nclades, giving probability weights to each clade. Weights are normalized among children of each node, during traversal.
	double p, total_weight;
	long node;
	long clade = root;
	while(clade>=Ntips){
		node = clade-Ntips;
		// determine normalization factor for weights
		total_weight = 0;
		for(long c=node2first_child[node]; c<=node2last_child[node]; ++c){
			total_weight += clade2weight[children[c]];
		}
		// choose random child
		if(total_weight==0){
			// this should not happen, but if it does continue anyway by choosing a random child
			clade = children[uniformIntWithin(node2first_child[node],node2last_child[node])];
			continue;
		}
		clade = children[node2last_child[node]]; // leftover child (if random p<1)
		// p = double(rand())/RAND_MAX; // random number within [0,1]. Note that R package builders discourage the use of rand(), so use R::runif() instead
		p = R::runif(0.0,1.0);
		for(long c=node2first_child[node]; c<=node2last_child[node]; ++c){
			p -= clade2weight[children[c]]/total_weight;
			if((p<=0) && (clade2weight[children[c]]>0)){
				clade = children[c];
				break;
			}
		}
	}
	return clade;
}





// ###############################################################################################
// Ancestral state reconstruction of discrete characters
// using fixed-rate continuous-time Markov models and the rerooting method by Yang et al. (1995).




// Set an internal node as a new root for a phylogenetic tree, by redefining edge directions
// The number of tips & nodes & edges (and their indices) will remain the same, but some edges will change direction i.f. tree_edge[e,0] will be swapped with tree_edge[e,1].
// The original root may become a singleton node (i.e. with only one outgoing edge) and the new root may become trifurcating (i.e. with 3 outgoing edges)
// Requirements:
//	 The original tree can include mono- and multi-furcations
//   The original tree must be rooted.
template<class ARRAY_TYPE>
void reroot_tree_at_node(	const long 	Ntips,
							const long 	Nnodes,
							const long 	Nedges,
							const long	old_root,		// (INPUT) integer in Nnodes,..,Ntips+Nnodes-1.
							const long 	new_root,		// (INPUT) integer in Nnodes,..,Ntips+Nnodes-1. Node to be set as new root.	
							ARRAY_TYPE	&tree_edge){	// (INPUT/OUTPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1). After rerooting, this array will be modified to reflect the new edge directions (i.e. tree_edge[e,0] will be swapped with tree_edge[e,1] for some e).			
	// get incoming edge for each clade
	std::vector<long> incoming_edge_per_clade;
	get_incoming_edge_per_clade(Ntips, Nnodes, Nedges, tree_edge, incoming_edge_per_clade);
	
	// traverse from new_root towards old_root
	long clade = new_root;
	long parent, edge;
	while(clade!=old_root){
		edge	= incoming_edge_per_clade[clade];
		parent 	= tree_edge[edge*2+0];
		// swap parent & child in this edge
		tree_edge[edge*2+0] = clade;
		tree_edge[edge*2+1] = parent;
		clade	= parent;
	}
}



// Set an internal node as a new root for a phylogenetic tree, by redefining edge directions
// Similar to the function above, but here incoming_edge_per_clade[] is also updated
template<class ARRAY_TYPE>
void reroot_tree_at_node(	const long 			Ntips,
							const long 			Nnodes,
							const long 			Nedges,
							const long			old_root,					// (INPUT) integer in Ntips,..,Ntips+Nnodes-1.
							const long 			new_root,					// (INPUT) integer in Ntips,..,Ntips+Nnodes-1. Node to be set as new root.	
							ARRAY_TYPE			&tree_edge,					// (INPUT/OUTPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1). After rerooting, this array will be modified to reflect the new edge directions (i.e. tree_edge[e,0] will be swapped with tree_edge[e,1] for some e).
							std::vector<long>	&incoming_edge_per_clade){	// (INPUT/OUTPUT) 1D array of size Nclades, with elements in 0,..,Nedges-1. Will be updated after the rerooting.
	
	// traverse from new_root towards old_root
	long clade = new_root;
	long parent, edge, previous_edge=-1;
	while(clade!=old_root){
		edge	= incoming_edge_per_clade[clade];
		parent 	= tree_edge[edge*2+0];
		if(previous_edge>=0) incoming_edge_per_clade[clade] = previous_edge;
		// swap parent & child in this edge
		tree_edge[edge*2+0] = clade;
		tree_edge[edge*2+1] = parent;
		previous_edge = edge; // keep a record of this edge for the next iteration, so that we can update the incoming edge after retrieving some old information
		clade = parent;
	}
	if(old_root!=new_root){
		incoming_edge_per_clade[old_root] = previous_edge;
		incoming_edge_per_clade[new_root] = -1;
	}
}




// Collapse tree nodes (and their descending subtrees) into tips, whenever all descending tips have a distance from a node below a certain threshold (or whenever the sum of edges descending from the node is below the threshold, see option criterion)
// If shorten==true:
//   Collapsed nodes will be turned into tips, while retaining the length of their incoming edge (thus the tree is shortened)
// If shorten==false:
//   Whenever all tips descending from some node have a distance from the node below a certain threshold, remove all tips and make the node into a tip, extending its incoming edge by length L, where L was the longest distance to any of the tips
//   In other words, from each subtree below a certain diversity threshold, replace it with a single tip extending to the maximum depth of the subtree
// This function can be used to get the "coarse structure" of a tree
// This function guarantees that the new_root will have index = Ntips_new
// [[Rcpp::export]]
Rcpp::List collapse_tree_at_resolution_CPP(	const long			Ntips,
											const long 			Nnodes,
											const long			Nedges,
											const IntegerVector	&tree_edge,			// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
											const NumericVector &edge_length, 		// (INPUT) 1D array of size Nedges, or an empty std::vector (all branches have length 1)
											const double		resolution,			// (INPUT) maximum (inclusive) phylogenetic distance of descending tips from the node to be collapsed
											const bool			shorten,			// (INPUT) if true, then collapsed nodes are replaced with a tip at the same location, hence potentially shortening the tree
											const std::string	&criterion){		// (INPUT) criterion by which to collapse nodes. 'sum_tip_paths': interpret resolution as a max allowed sum of tip distances from the collapsed node. 'max_tip_depth': resolution refers to the max tip distance from the collapsed node. 'max_tip_pair_dist': resolution refers to the max distance between any pair of tips descending from the collapsed node.
	const long Nclades = Ntips + Nnodes;
	long clade, parent, edge, node;
	const bool unit_edge_lengths = (edge_length.size()==0);
	const bool sum_tip_paths = (criterion=="sum_tip_paths");
	const bool max_tip_depth = (criterion=="max_tip_depth");
	const bool tip_pair_dist = (criterion=="max_tip_pair_dist");
	
	// determine incoming edge per clade
	std::vector<long> incoming_edge_per_clade(Nclades,-1);
	for(long edge=0; edge<Nedges; ++edge){
		incoming_edge_per_clade[tree_edge[edge*2+1]] = edge;
	}

	// find root using the mapping clade-->incoming_edge
	const long root = get_root_from_incoming_edge_per_clade(Ntips, tree_edge, incoming_edge_per_clade);
		
	// get tree traversal route (root --> tips)											
	std::vector<long> traversal_queue, node2first_edge, node2last_edge, edges;
	get_tree_traversal_root_to_tips(	Ntips,
										Nnodes,
										Nedges,
										root,
										tree_edge,
										true,	// include tips
										false,	// edge mapping is not precalculated
										traversal_queue,
										node2first_edge,
										node2last_edge,
										edges,
										false,
										"");
										
	// Step 1: calculate distance to furthest descending tip for each node (traverse tips-->root)
	std::vector<double> clade2max_tip_depth(Nclades,0);
	std::vector<double> clade2sum_tip_paths((sum_tip_paths ? Nclades : 0),0);
	std::vector<double> clade2max_tip_pair_dist((tip_pair_dist ? Nclades : 0),0);
	std::vector<double> clade2farthest_tip(Nclades,-1);
	for(long tip=0; tip<Ntips; ++tip){
		clade2farthest_tip[tip] = tip;
	}
	for(long q=traversal_queue.size()-1, clade; q>=1; --q){
		clade  = traversal_queue[q];
		edge   = incoming_edge_per_clade[clade];
		parent = tree_edge[edge*2+0];
		const double candidate_distance = clade2max_tip_depth[clade] + (unit_edge_lengths ? 1.0 : edge_length[edge]);
		if(tip_pair_dist){
			// keep track of max distance between any pair of descending tips
			// this must come before updating clade2max_tip_depth[parent]
			if(clade2farthest_tip[parent]>=0){
				clade2max_tip_pair_dist[parent] = max(clade2max_tip_pair_dist[parent], max(clade2max_tip_depth[parent]+candidate_distance, clade2max_tip_pair_dist[clade]));
			}else{
				clade2max_tip_pair_dist[parent] = clade2max_tip_pair_dist[clade];
			}
		}
		if((clade2farthest_tip[parent]<0) || (clade2max_tip_depth[parent]<candidate_distance)){
			// keep track of farthest descending tip
			clade2max_tip_depth[parent] = candidate_distance;
			clade2farthest_tip[parent]	= clade2farthest_tip[clade];
		}
		if(sum_tip_paths){
			// keep track of sum of paths to descending tips
			clade2sum_tip_paths[parent] += clade2sum_tip_paths[clade] + (unit_edge_lengths ? 1.0 : edge_length[edge]);
		}
	}
	if(tip_pair_dist){
		// make sure nodes with just 
	}
	const std::vector<double> &clade2criterion = (sum_tip_paths ? clade2sum_tip_paths : (max_tip_depth ? clade2max_tip_depth : clade2max_tip_pair_dist));
	
	// Step 2: count how many tips/nodes/edges to extract
	// traverse root-->tips in depth-first-search order, stoppping at nodes below the threshold
	// use a scratch_stack for traversing nodes
	long Ntips_new  = 0;
	long Nnodes_new = 0;
	long Nedges_new = 0;
	std::vector<long> scratch_stack;
	scratch_stack.reserve(floor(2*log(Ntips)/log(2.0))); // rough estimate of typical tree depth x 2. scratch_stack may be resized along the way if needed.
	scratch_stack.push_back(root);
	while(scratch_stack.size()>0){
		clade = scratch_stack.back();
		node  = clade - Ntips;
		scratch_stack.pop_back();
		if((clade2criterion[clade]<=resolution) || (clade<Ntips)){
			// this is a tip, or a node that should be collapsed into a new tip
			++Ntips_new;
			continue; // don't expand children
		}else{
			++Nnodes_new;
		}
		// traverse outgoing edges and children
		for(long e=node2first_edge[node]; e<=node2last_edge[node]; ++e){
			edge  = edges[e];
			scratch_stack.push_back(tree_edge[edge*2+1]); // add child to stack for further exploration in the next iteration
			++Nedges_new;
		}
	}
	const long Nclades_new 	= Ntips_new + Nnodes_new;
	const long Ncollapsed 	= Nnodes - Nnodes_new; // number of collapsed nodes
	
	// Step 3: Traverse again root-->tips (depth-first-search) and create mappings old-->new
	std::vector<long> old2new_clade(Nclades,-1), old2new_edge(Nedges,-1), collapsed_nodes, representative_tips;
	collapsed_nodes.reserve(Ncollapsed);
	representative_tips.reserve(Ncollapsed);
	long next_new_tip  = 0;
	long next_new_node = 0;
	long next_new_edge = 0;
	scratch_stack.clear();
	scratch_stack.push_back(root);
	while(scratch_stack.size()>0){
		clade = scratch_stack.back();
		node  = clade - Ntips;
		scratch_stack.pop_back();
		if((clade2criterion[clade]<=resolution) || (clade<Ntips)){
			// this is a tip, or a node that should be collapsed into a new tip
			old2new_clade[clade] = (next_new_tip++);
			if(clade>=Ntips){
				// keep record of collapsed node
				collapsed_nodes.push_back(clade-Ntips);
				representative_tips.push_back(clade2farthest_tip[clade]);
			}
			continue; // don't expand children
		}else{
			old2new_clade[clade] = Ntips_new + (next_new_node++);
		}
		// traverse outgoing edges and children
		for(long e=node2first_edge[node]; e<=node2last_edge[node]; ++e){
			edge  = edges[e];
			scratch_stack.push_back(tree_edge[edge*2+1]); // add child to stack for further exploration in the next iteration
			old2new_edge[edge] = (next_new_edge++);
		}
	}
	
	// Step 4: Extract collapsed tree and generate new2old mappings
	std::vector<long> new_tree_edge(2*Nedges_new), new2old_clade(Nclades_new), new2old_edge(Nedges_new);
	std::vector<double> new_edge_length(shorten ? 0 : Nedges_new);
	for(long edge=0, new_edge, child; edge<Nedges; ++edge){
		new_edge = old2new_edge[edge];
		if(new_edge<0) continue; // this edge is not to be kept
		child = tree_edge[edge*2+1];
		new2old_edge[new_edge] = edge;
		new_tree_edge[new_edge*2+0] = old2new_clade[tree_edge[edge*2+0]];
		new_tree_edge[new_edge*2+1] = old2new_clade[child];
		if(!shorten){
			new_edge_length[new_edge] = (unit_edge_lengths ? 1.0 : edge_length[edge]);
			if((child>=Ntips) && (old2new_clade[child]<Ntips_new)) new_edge_length[new_edge] += clade2max_tip_depth[child]; // this child is a collapsed node, so extend incoming edge to maximum depth of any old descendant
		}
	}
	for(long clade=0, new_clade; clade<Nclades; ++clade){
		new_clade = old2new_clade[clade];
		if(new_clade>=0) new2old_clade[new_clade] = clade;
	}

	return Rcpp::List::create(	Rcpp::Named("new_tree_edge") 		= Rcpp::wrap(new_tree_edge),
								Rcpp::Named("new_edge_length") 		= Rcpp::wrap(new_edge_length), // only relevant if (shorten==false)
								Rcpp::Named("new2old_clade") 		= Rcpp::wrap(new2old_clade),
								Rcpp::Named("new2old_edge") 		= Rcpp::wrap(new2old_edge),
								Rcpp::Named("old2new_clade") 		= Rcpp::wrap(old2new_clade),
								Rcpp::Named("collapsed_nodes") 		= Rcpp::wrap(collapsed_nodes),
								Rcpp::Named("representative_tips") 	= Rcpp::wrap(representative_tips),
								Rcpp::Named("new_root") 			= old2new_clade[root], // in newer implementations this is actually guaranteed to be = Ntips_new
								Rcpp::Named("Ntips_new") 			= Ntips_new,
								Rcpp::Named("Nnodes_new") 			= Nnodes_new,
								Rcpp::Named("Nedges_new") 			= Nedges_new);
}






// Trim a phylogenetic tree by cutting off tips and nodes, so that all remaining tips have height<=max_height.
// Note that some edge lengths will be notified (edges cut will be shortened)
// If initially all tips had height>=max_height, then the trimmed tree will be ultrametric
// [[Rcpp::export]]
Rcpp::List trim_tree_at_height_CPP(	const long			Ntips,
									const long 			Nnodes,
									const long			Nedges,
									const IntegerVector	&tree_edge,					// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
									const NumericVector &edge_length, 				// (INPUT) 1D array of size Nedges, or an empty std::vector (all branches have length 1)
									const double		max_distance_from_root){	// (INPUT) phylogenetic distance from the root, at which to trim the tree
	const long Nclades = Ntips + Nnodes;
	long parent, child, edge, node;
	const bool unit_edge_lengths = (edge_length.size()==0);
	
	// find root
	const long root = get_root_clade(Ntips, Nnodes, Nedges, tree_edge);
		
	// get tree traversal route (root --> tips)											
	std::vector<long> traversal_queue, node2first_edge, node2last_edge, edges;
	get_tree_traversal_root_to_tips(	Ntips,
										Nnodes,
										Nedges,
										root,
										tree_edge,
										true,	// include tips
										false,	// edge mapping is not precalculated
										traversal_queue,
										node2first_edge,
										node2last_edge,
										edges,
										false,
										"");
																				
	// Step 1: Determine which clades to keep (traverse root-->tips)
	// (include calculation of distances_from_root)
	std::vector<char> keep_clade(Nclades);
	std::vector<double> distances_from_root(Nclades);
	distances_from_root[root] = 0;
	keep_clade[root] = true;
	for(long q=0; q<traversal_queue.size(); ++q){
		parent = traversal_queue[q];
		if(parent<Ntips) continue;
		node = parent - Ntips;
		for(long e=node2first_edge[node]; e<=node2last_edge[node]; ++e){
			edge 	= edges[e];
			child 	= tree_edge[edge*2+1];
			distances_from_root[child] = (unit_edge_lengths ? 1.0 : edge_length[edge]) + distances_from_root[parent];
			keep_clade[child] = ((distances_from_root[parent]<max_distance_from_root) || (distances_from_root[child]<=max_distance_from_root)); // keep a child iff the child is not above the max_distance_from_root, or the parent is strictly below max_distance_from_root
		}
	}
	
	// Step 2: Extract subtree
	// edges to be kept is determined based on clades to be kept
	long Ntips_new, Nnodes_new, Nedges_new, new_root;
	std::vector<long> new2old_clade, old2new_clade, new2old_edge, new_tree_edge;
	get_arbitrary_subtree(	Ntips,
							Nnodes,
							Nedges,
							tree_edge,
							keep_clade,
							Ntips_new,
							Nnodes_new,
							Nedges_new,
							new_root,
							new2old_clade,
							old2new_clade,
							new2old_edge,
							new_tree_edge);
							
	// figure out which nodes became tips
	long Ntips_ex_nodes = 0;
	for(long tip_new=0; tip_new<Ntips_new; ++tip_new){
		if(new2old_clade[tip_new]>=Ntips) ++Ntips_ex_nodes;
	}
	std::vector<long> new_tips_ex_nodes(Ntips_ex_nodes);
	for(long tip_new=0, counter=0; tip_new<Ntips_new; ++tip_new){
		if(new2old_clade[tip_new]>=Ntips) new_tips_ex_nodes[counter++] = tip_new;
	}
	

							
	// Step 3: Calculate edge lengths for extracted subtree
	// (trim terminal edges if needed)
	std::vector<long> new_edge_length(Nedges_new);
	long Nedges_trimmed = 0;
	for(long edge_new=0, edge; edge_new<Nedges_new; ++edge_new){
		edge 	= new2old_edge[edge_new];
		parent 	= tree_edge[edge*2+0];
		child 	= tree_edge[edge*2+1];
		const double old_length = (unit_edge_lengths ? 1.0 : edge_length[edge]);
		if(old_length>max_distance_from_root-distances_from_root[parent]){
			new_edge_length[edge_new] = max_distance_from_root-distances_from_root[parent];
			++Nedges_trimmed;
		}else{
			new_edge_length[edge_new] = old_length;
		}
	}
					
	return Rcpp::List::create(	Rcpp::Named("Ntips_new")		= Ntips_new,
								Rcpp::Named("Nnodes_new")		= Nnodes_new,
								Rcpp::Named("Nedges_new")		= Nedges_new,		// number of edges kept
								Rcpp::Named("Nedges_trimmed")	= Nedges_trimmed,	// number of kept edges that were trimmed (i.e. reduced in length)
								Rcpp::Named("new_root")			= new_root,			// new root clade index. This is actually guaranteed to be = Ntips_new.
								Rcpp::Named("new2old_clade")	= Rcpp::wrap(new2old_clade),
								Rcpp::Named("new2old_edge")		= Rcpp::wrap(new2old_edge),
								Rcpp::Named("new_tree_edge")	= Rcpp::wrap(new_tree_edge),
								Rcpp::Named("new_edge_length")	= Rcpp::wrap(new_edge_length),
								Rcpp::Named("new_tips_ex_nodes")= Rcpp::wrap(new_tips_ex_nodes)); // new tips that used to be nodes
}


// Eliminate multifurcations in a tree by replacing them with multiple bifurcations
// Tips indices remain the same, but edge indices and the total number of nodes/edges may increase (if the tree includes multifurcations)
// Old nodes retain their index, and new nodes will have indices Nnodes,...,(Nnew_nodes-1)
// The tree need not be rooted
template<class INTEGER_ARRAY, class NUMERIC_ARRAY>
void multifurcations_to_bifurcations(	const long			Ntips,
										const long 			Nnodes,
										const long			Nedges,
										const INTEGER_ARRAY	&tree_edge,			// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
										const NUMERIC_ARRAY &edge_length,		// (INPUT) 1D array of size Nedges listing edge lengths, or an empty vector (all edges have length 1)
										const double		dummy_edge_length,	// (INPUT) length to be used for new auxiliary edges (connecting old to new nodes) when splitting multifurcations. This will typically be zero, or a small number if zero-length edges are not desired.
										long				&Nnew_nodes,		// (OUTPUT) number of nodes in the new tree
										long				&Nnew_edges,		// (OUTPUT) number of edges in the new tree
										std::vector<long>	&new_tree_edge,		// (OUTPUT) 2D array of size Nnew_edges x 2, in row-major format, with elements in 0,..,(Nclades-1)
										std::vector<double>	&new_edge_length,	// (OUTPUT) 1D array of size Nnew_edges, listing the lengths of edges in the new tree
										std::vector<long>	&old2new_edge){		// (OUTPUT) 1D array of size Nedges, with values in 0,..,(Nnew_edges-1) mapping old edges to new edges
	long edge, child, Nchildren;
	
	// get edge mappings
	std::vector<long> node2first_edge, node2last_edge, edge_mapping;
	get_node2edge_mappings(	Ntips,
							Nnodes,
							Nedges,
							tree_edge,
							node2first_edge,
							node2last_edge,
							edge_mapping);
							
	// determine number of nodes/edges in the new tree (based on the number & size of multifurcations)
	Nnew_edges = Nedges;
	for(long node=0; node<Nnodes; ++node){
		Nchildren = node2last_edge[node] - node2first_edge[node] + 1;
		if(Nchildren>2) Nnew_edges += (Nchildren-2);
	}
	Nnew_nodes = Nnodes + (Nnew_edges - Nedges);
	
	if(Nnew_edges==Nedges){
		// return tree without modification
		new_tree_edge 	= Rcpp::as< std::vector<long> >(tree_edge);
		new_edge_length	= (edge_length.size()==0 ? std::vector<double>(Nedges,1) : Rcpp::as< std::vector<double> >(edge_length));
		old2new_edge.resize(Nedges);
		for(edge=0; edge<Nedges; ++edge) old2new_edge[edge] = edge;
		return;
	}
							
	// iterate through nodes and expand any multifurcations
	new_tree_edge.clear();
	new_edge_length.clear();
	new_tree_edge.reserve(2*Nnew_edges);
	new_edge_length.reserve(Nnew_edges);
	old2new_edge.resize(Nedges);
	long next_new_node = Nnodes;
	for(long node=0, clade, next_parent; node<Nnodes; ++node){
		clade = Ntips + node;
		Nchildren = node2last_edge[node] - node2first_edge[node] + 1;
		if(Nchildren>0){
			// the first child is always preserved
			edge = edge_mapping[node2first_edge[node]];
			new_tree_edge.push_back(clade);
			new_tree_edge.push_back(tree_edge[2*edge+1]);
			new_edge_length.push_back(edge_length.size()==0 ? 1.0 : edge_length[edge]);
			old2new_edge[edge] = new_edge_length.size()-1;
		}
		if(Nchildren<=2){
			// node does not multifurcate, so also preserve 2nd child (if available) and move on to the next node
			if(Nchildren>1){
				edge = edge_mapping[node2first_edge[node]+1];
				new_tree_edge.push_back(clade);
				new_tree_edge.push_back(tree_edge[2*edge+1]);
				new_edge_length.push_back(edge_length.size()==0 ? 1.0 : edge_length[edge]);
				old2new_edge[edge] = new_edge_length.size()-1;
			}
			continue;
		}
		// for all children but the first and last create a new node
		next_parent = clade;
		for(long e=1+node2first_edge[node]; e<node2last_edge[node]; ++e){
			edge  = edge_mapping[e];
			child = tree_edge[2*edge+1];
			// create new edge next_parent --> next_new_node
			new_tree_edge.push_back(next_parent);
			new_tree_edge.push_back(Ntips+next_new_node);
			new_edge_length.push_back(dummy_edge_length);
			// add child to the new node
			new_tree_edge.push_back(Ntips+next_new_node);
			new_tree_edge.push_back(child);
			new_edge_length.push_back(edge_length.size()==0 ? 1.0 : edge_length[edge]);
			old2new_edge[edge] = new_edge_length.size()-1;
			next_parent = Ntips + next_new_node;
			++next_new_node;
		}
		// add last child to the last new node (next_parent)
		edge  = edge_mapping[node2last_edge[node]];
		child = tree_edge[2*edge+1];
		new_tree_edge.push_back(next_parent);
		new_tree_edge.push_back(child);
		new_edge_length.push_back(edge_length.size()==0 ? 1.0 : edge_length[edge]);
		old2new_edge[edge] = new_edge_length.size()-1;
	}
}



// Replace multifurcations with multiple bifurcations
// The tree need not be rooted
// Rcpp wrapper function
// [[Rcpp::export]]
Rcpp::List multifurcations_to_bifurcations_CPP(	const long			Ntips,
												const long 			Nnodes,
												const long			Nedges,
												const IntegerVector	&tree_edge,			// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
												const NumericVector &edge_length,		// (INPUT) 1D array of size Nedges listing edge lengths, or an empty vector (all edges have length 1)
												const double		dummy_edge_length){	// (INPUT) length to be used for new auxiliary edges (connecting old to new nodes) when splitting multifurcations. This will typically be zero, or a small number if zero-length edges are not desired.
	long Nnew_nodes, Nnew_edges;
	std::vector<long> new_tree_edge, old2new_edge;
	std::vector<double> new_edge_length;
	multifurcations_to_bifurcations(Ntips,
									Nnodes,
									Nedges,
									tree_edge,
									edge_length,
									dummy_edge_length,
									Nnew_nodes,
									Nnew_edges,
									new_tree_edge,
									new_edge_length,
									old2new_edge);
											
	return Rcpp::List::create(	Rcpp::Named("Nnew_nodes") 		= Nnew_nodes,
								Rcpp::Named("Nnew_edges") 		= Nnew_edges,
								Rcpp::Named("new_tree_edge") 	= Rcpp::wrap(new_tree_edge),
								Rcpp::Named("new_edge_length") 	= Rcpp::wrap(new_edge_length),
								Rcpp::Named("old2new_edge") 	= Rcpp::wrap(old2new_edge));
}


// Eliminate short edges by merging affected nodes/tips into multifurcations
// The tree must be rooted
// [[Rcpp::export]]
Rcpp::List merge_short_edges_CPP(	const long					Ntips,
									const long 					Nnodes,
									const long					Nedges,
									const std::vector<long>		&tree_edge,				// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
									const std::vector<double> 	&edge_length,			// (INPUT) 1D array of size Nedges listing edge lengths, or an empty vector (all edges have length 1)
									const double				edge_length_epsilon,	// (INPUT) non-negative number. Max edge length for an edge to be considered zero-length. Typically 0 or some small positive number.
									const bool					force_keep_tips){		// (INPUT) if true, all tips are kept even if their incoming edges are short. Note that tip indices may still change.
	const long Nclades = Ntips + Nnodes;
	const bool unit_edge_lengths = (edge_length.size()==0);

	// get edge mappings
	std::vector<long> node2first_edge, node2last_edge, edge_mapping;
	get_node2edge_mappings(	Ntips,
							Nnodes,
							Nedges,
							tree_edge,
							node2first_edge,
							node2last_edge,
							edge_mapping);
							
	// get root
	const long root = get_root_clade(Ntips, Nnodes, Nedges, tree_edge);
							
	// get a lower bound for the number of edges & clades in the new tree
	// note that some short edges may be kept if they are extended by ancestral short edges
	long Nnew_edges_min = 0;
	for(long edge=0, child; edge<Nedges; ++edge){
		child = tree_edge[edge*2+1];
		const double length = (unit_edge_lengths ? 1.0 : edge_length[edge]);
		if((length>edge_length_epsilon) || (force_keep_tips && (child<Ntips))) ++Nnew_edges_min;
	}
	
	if(Nnew_edges_min==Nedges){
		// return tree without modification
		std::vector<long> new2old_clade(Nclades), new2old_edge(Nedges);
		for(long clade=0; clade<Nclades; ++clade) new2old_clade[clade] = clade;
		for(long edge=0; edge<Nedges; ++edge) new2old_edge[edge] = edge;
		return Rcpp::List::create(	Rcpp::Named("Nnew_tips") 		= Ntips,
									Rcpp::Named("Nnew_nodes") 		= Nnodes,
									Rcpp::Named("Nnew_edges") 		= Nedges,
									Rcpp::Named("new_tree_edge") 	= Rcpp::wrap(tree_edge),
									Rcpp::Named("new_edge_length") 	= Rcpp::wrap(unit_edge_lengths ? std::vector<double>(Nedges,1.0) : edge_length),
									Rcpp::Named("new2old_clade") 	= Rcpp::wrap(new2old_clade),
									Rcpp::Named("new2old_edge") 	= Rcpp::wrap(new2old_edge),
									Rcpp::Named("root")				= root);
	}
	
	// iterate through edges (root-->tips) and merge any short edges
	// use a scratch_stack to traverse root-->tips via depth-first-search
	std::vector<long> new2old_clade, new2old_edge;
	std::vector<long> 	new_tree_edge; // first populate with old node/tip indices, then update
	std::vector<double> new_edge_length;
	new2old_clade.reserve(Nnew_edges_min+1);
	new2old_edge.reserve(Nnew_edges_min);
	new_tree_edge.reserve(Nnew_edges_min*2);
	new_edge_length.reserve(Nnew_edges_min);
	std::vector<long> 	current_tree_edge = tree_edge; 	// parents will be updated here, as short edges are being eliminated (parents of sub-edges are moved closer to the root). Hence, edges must be traversed root-->tip.
	std::vector<double> current_edge_length; 			// edge lengths will be updated here, as short edges are being eliminated and their lengths appended to sub-edges
	if(!unit_edge_lengths) current_edge_length = edge_length;
	else current_edge_length.assign(Nedges,1);

	new2old_clade.push_back(root); // always keep root
	std::vector<long> scratch_stack; // will list nodes whose edges still need to be checked
	scratch_stack.reserve(floor(2*log(Ntips)/log(2.0))); // rough estimate of typical tree depth x 2. scratch_stack may be resized along the way if needed.
	scratch_stack.push_back(root); // start traversal at root
	while(scratch_stack.size()>0){
		const long clade = scratch_stack.back();
		scratch_stack.pop_back();
		const long node = clade - Ntips;
		for(long e=node2first_edge[node], edge, parent, child; e<=node2last_edge[node]; ++e){
			edge 	= edge_mapping[e];
			parent	= current_tree_edge[edge*2+0]; // may be different from clade, if this edge was extended due to elimination of its parent edge
			child 	= current_tree_edge[edge*2+1];
			if(child>=Ntips) scratch_stack.push_back(child); // add child node to stack for further exploration later on
	
			const double length = current_edge_length[edge];
			if((length>edge_length_epsilon) || (force_keep_tips && (child<Ntips))){
				// keep this edge and its child
				new2old_clade.push_back(child);
				const long next_new_edge = new2old_edge.size();
				new2old_edge.push_back(edge);
				new_tree_edge.resize(new_tree_edge.size()+2);
				new_tree_edge[next_new_edge*2+0] = parent;
				new_tree_edge[next_new_edge*2+1] = child;
				new_edge_length.push_back(length);
			}else{
				// eliminate this short edge by merging child into parent
				// this eliminates the child tip/node
				// the child's children become the parent's children
				if(child<Ntips){
					// child is a tip, so just eliminate this tip (i.e. do nothing)
				}else{
					// child is a node, so attach its children to the parent node
					const long cnode = child-Ntips;
					for(long se=node2first_edge[cnode], sub_edge; se<=node2last_edge[cnode]; ++se){
						sub_edge = edge_mapping[se];
						current_tree_edge[sub_edge*2+0] = parent;
						current_edge_length[sub_edge] 	+= length;
					}
				}
			
			}
		}
	}
	const long Nnew_edges  = new2old_edge.size();
	const long Nnew_clades = 1+Nnew_edges;
	
	// update clade indices in new_tree_edge
	std::vector<long> old2new_clade(Nclades,-1);
	for(long new_clade=0; new_clade<Nnew_clades; ++new_clade){
		old2new_clade[new2old_clade[new_clade]] = new_clade;
	}
	for(long new_edge=0; new_edge<Nnew_edges; ++new_edge){
		new_tree_edge[new_edge*2+0] = old2new_clade[new_tree_edge[new_edge*2+0]];
		new_tree_edge[new_edge*2+1] = old2new_clade[new_tree_edge[new_edge*2+1]];
	}
	
	// correct new clade indices so that tips are indexed 0,..,Nnew_tips-1 and nodes are indexed Nnew_tips,..,Nnew_clades-1
	long Nnew_tips, Nnew_nodes;
	std::vector<long> new2newer_clade;
	reindex_clades(	Nnew_clades,
					Nnew_edges,
					new_tree_edge,
					true, // ensure root is re-indexed to Nnew_tips
					Nnew_tips,
					Nnew_nodes,
					new2newer_clade);
	for(long clade=0; clade<Nclades; ++clade){
		if(old2new_clade[clade]>=0){
			old2new_clade[clade] = new2newer_clade[old2new_clade[clade]];
			new2old_clade[old2new_clade[clade]] = clade;
		}
	}
	for(long new_edge=0; new_edge<Nnew_edges; ++new_edge){
		new_tree_edge[new_edge*2+0] = new2newer_clade[new_tree_edge[new_edge*2+0]];
		new_tree_edge[new_edge*2+1] = new2newer_clade[new_tree_edge[new_edge*2+1]];
	}

	return Rcpp::List::create(	Rcpp::Named("Nnew_tips") 		= Nnew_tips,
								Rcpp::Named("Nnew_nodes") 		= Nnew_nodes,
								Rcpp::Named("Nnew_edges") 		= Nnew_edges,
								Rcpp::Named("new_tree_edge") 	= Rcpp::wrap(new_tree_edge),
								Rcpp::Named("new_edge_length") 	= Rcpp::wrap(new_edge_length),
								Rcpp::Named("new2old_clade") 	= Rcpp::wrap(new2old_clade),
								Rcpp::Named("new2old_edge") 	= Rcpp::wrap(new2old_edge),
								Rcpp::Named("root")				= Nnew_tips);
}


// Generate a random phylogenetic tree under a simple speciation model, where species are born or go extinct as a Poissonian process
// New species are added by splitting one of the currently extant tips (chosen randomly) into Nsplits new tips
// Special case is the Yule model: New species appear as a Poisson process with a constant per-capita birth rate, and without extinctions
// More generally, the species birth rate can be a power-law function of extant tips count: birth_rate = intercept + factor*number_of_extant_tips^exponent
// Similarly, the death rate of tips can be a power-law function of extant tip count: death_rate = intercept + factor*number_of_extant_tips^exponent
// The resulting tree will be bifurcating (if Nsplits=2) or multifurcating (if Nsplits>2).
// The simulation is halted as soon as Ntips>=max_tips (if max_tips>0) and/or time>=max_time (if max_time>0) and/or time>=max_time_since_equilibrium+equilibrium_time (if max_time_since_equilibrium>0)
// Reference:
//   Steel and McKenzie (2001). Properties of phylogenetic trees generated by Yule-type speciation models. Mathematical Biosciences. 170:91-112.
// [[Rcpp::export]]
Rcpp::List generate_random_tree_CPP(const long 	 	max_tips,					// (INPUT) max number of tips (extant tips, if coalescent==true). If <=0, this constraint is ignored.
									const double	max_time,					// (INPUT) max simulation time. If <=0, this constraint is ignored.
									const double	max_time_since_equilibrium,	// (INPUT) max simulation time, counting from the first time point where death_rate-birth_rate changed sign. This may be used as an alternative to (or in conjunction with) max_time to ensure the tree has reached speciation/extinction equilibrium. If <0, this constraint is ignored.
									const double 	birth_rate_intercept,		// (INPUT) intercept of Poissonian rate at which new tips are added to the tree
									const double 	birth_rate_factor,			// (INPUT) power-law factor of Poissonian rate at which new tips are added to the tree
									const double 	birth_rate_exponent,		// (INPUT) power-law exponent of Poissonian rate at which new tips are added to the tree
									const double 	death_rate_intercept,		// (INPUT) intercept of Poissonian rate at which extant tips are removed from the tree
									const double 	death_rate_factor,			// (INPUT) power-law factor of Poissonian rate at which extant tips are removed from the tree
									const double 	death_rate_exponent,		// (INPUT) power-law exponent of Poissonian rate at which extant tips are removed from the tree
									const std::vector<double>	&additional_rates_times,	// (INPUT) optional 1D array of size NAR, listing time points (in ascending order) for custom additional birth and/or death rates. Can be empty. The time series is repeated periodically if needed.
									const std::vector<double>	&additional_birth_rates_pc,	// (INPUT) optional 1D array of size NAR, listing custom per-capita birth rates (additive to the power law). Can be empty.
									const std::vector<double>	&additional_death_rates_pc,	// (INPUT) optional 1D array of size NAR, listing custom per-capita birth rates (additive to the power law). Can be empty.
									const bool		additional_periodic,		// (INPUT) if true, additional pc birth & death rates are extended periodically if needed. Otherwise they are extended with zeros.
									const bool		coalescent,					// (INPUT) whether to return only the coalescent tree (i.e. including only extant tips)
									const long		Nsplits,					// (INPUT) number of children to create at each diversification event. Must be at least 2. For a bifurcating tree this should be set to 2. If >2, then the tree will be multifurcating.
									const bool		as_generations,				// (INPUT) if false, then edge lengths correspond to time. If true, then edge lengths correspond to generations (hence if coalescent==false, all edges will have unit length).
									const bool		include_birth_times,		// (INPUT) if true, then the times of speciations (in order of occurrence) will also be returned
									const bool		include_death_times){		// (INPUT) if true, then the times of extinctions (in order of occurrence) will also be returned
	const long expected_Nclades = (max_tips<0 ? 2l : max_tips);
	long next_Nsplits = max(2l, Nsplits);
	std::vector<long> tree_edge;
	std::vector<long> extant_tips;
	std::vector<double> clade2end_time;
	std::vector<double> birth_times, death_times;
	tree_edge.reserve(expected_Nclades*2);
	extant_tips.reserve(ceil(expected_Nclades/2.0)); // keep track of which clades are extant tips, as the tree is built
	clade2end_time.reserve(expected_Nclades); // keep track of time at which each clade split or went extinct (negative if clade is an extant tip)
	
	// prepare interpolators for additional birth & death rates
	LinearInterpolationFunctor<double> added_birth_rates_pc, added_death_rates_pc;
	const bool has_added_birth_rates = ((additional_rates_times.size()>0) && (additional_birth_rates_pc.size()>0));
	const bool has_added_death_rates = ((additional_rates_times.size()>0) && (additional_death_rates_pc.size()>0));
	if(has_added_birth_rates) added_birth_rates_pc = LinearInterpolationFunctor<double>(additional_rates_times, additional_birth_rates_pc,additional_periodic,0.0,true);
	if(has_added_death_rates) added_death_rates_pc = LinearInterpolationFunctor<double>(additional_rates_times, additional_death_rates_pc,additional_periodic,0.0,true);
	
	// create the first tip (which is also the root)
	long Ntips = 0; 	// current number of extant + extinct tips
	long Nclades = 0;	// current number of clades
	long root = 0;
	extant_tips.push_back(Nclades++);
	clade2end_time.push_back(-1);
	++Ntips;
	
	// create additional tips, by splitting existing tips at each step (turning the split parent tip into a node)
	long Nedges 		= 0;
	long Nbirths		= 0;
	long Ndeaths 		= 0;
	double time 		= 0;
	double total_rate	= INFTY_D;
	double equilibrium_time = INFTY_D;
	double initial_growth_rate = NAN_D;
	while(((max_tips<=0) || ((coalescent ? extant_tips.size() : Ntips)<max_tips)) && ((max_time<=0) || (time+1/total_rate<max_time)) && ((max_time_since_equilibrium<0) || (time-equilibrium_time+1/total_rate<max_time_since_equilibrium))){
		// determine time of next speciation or extinction event
		// prevent deaths if only one tip is left
		const double NEtips 	= extant_tips.size();
		const double birth_rate = birth_rate_intercept + birth_rate_factor * pow(NEtips, birth_rate_exponent) + (has_added_birth_rates ? max(0.0,added_birth_rates_pc(time))*NEtips : 0.0);		
		const double death_rate = (NEtips<=1 ? 0 : (death_rate_intercept + death_rate_factor * pow(NEtips, death_rate_exponent)) + (has_added_death_rates ? max(0.0,added_death_rates_pc(time))*NEtips : 0.0));
		if(std::isnan(initial_growth_rate)) initial_growth_rate = birth_rate-death_rate;
		if(((birth_rate-death_rate)*initial_growth_rate<0) && (equilibrium_time>time)) equilibrium_time = time; // first crossing of equilibrium state, so keep record
		total_rate = birth_rate+death_rate;
		time += random_exponential_distribution(total_rate);
		const bool birth = random_bernoulli(birth_rate/(birth_rate+death_rate));
				
		// randomly pick an existing tip to split or kill
		long tip   = uniformIntWithin(0,NEtips-1);
		long clade = extant_tips[tip];
		clade2end_time[clade] = time;
		
		if(birth){
			// split chosen tip into Nsplits daughter-tips & create new edges leading into those tips
			if(max_tips>0) next_Nsplits = min(1+max_tips-long(coalescent ? NEtips : Ntips), max(2l, Nsplits)); // temporarily reduce Ntips to stay within limits
			// child 1:
			++Nedges;
			++Nbirths;
			tree_edge.push_back(clade);
			tree_edge.push_back(Nclades);
			extant_tips[tip] = (Nclades++); // replace the old tip with one of the new ones
			clade2end_time.push_back(-1);
			if(include_birth_times) birth_times.push_back(time);
			// remaining children:
			for(long c=1; c<next_Nsplits; ++c){
				++Nedges;
				tree_edge.push_back(clade);
				tree_edge.push_back(Nclades);
				extant_tips.push_back(Nclades++);
				clade2end_time.push_back(-1);
				++Ntips;
			}
		}else{
			// kill chosen tip (remove from pool of extant tips); note that it still remains a tip, but it can't diversify anymore
			extant_tips[tip] = extant_tips.back();
			extant_tips.pop_back();
			++Ndeaths;
			if(include_death_times) death_times.push_back(time);
		}
		// abort if the user has interrupted the calling R program
		Rcpp::checkUserInterrupt();
	}
	
	// add a small dt at the end to make all edges non-zero length
	time += random_exponential_distribution(total_rate);
	if(max_time>0) time = min(time, max_time); // prevent going past max_time
	if(max_time_since_equilibrium>=0) time = min(time, max_time_since_equilibrium+equilibrium_time); // prevent going past max_time_since_equilibrium

	std::vector<double> edge_length(Nedges);
	if(as_generations){
		// set all edge lengths to 1
		edge_length.assign(Nedges,1);
	}else{
		// calculate edge lengths based on end times
		for(long edge=0, child; edge<Nedges; ++edge){
			child = tree_edge[edge*2+1];
			if(clade2end_time[child]>=0) edge_length[edge] = clade2end_time[child] - clade2end_time[tree_edge[edge*2+0]];
			else edge_length[edge] = time - clade2end_time[tree_edge[edge*2+0]];
		}	
	}
	double root_time = clade2end_time[root]; // root_time = time at which the root split
	
	// identify tips as the clades with no outgoing edge
	std::vector<bool> clade_is_tip(Nclades,true);
	for(long edge=0; edge<Nedges; ++edge){
		clade_is_tip[tree_edge[edge*2+0]] = false;
	}
	
	// re-number tip & node indices to conform with the phylo format, where tips are indexed first (0,..,Ntips-1) and nodes last (Ntips,..,Ntips+Nnodes-1).
	long Nnodes = Nclades - Ntips;
	std::vector<long> old2new_clade(Nclades,-1);
	long next_new_tip  = 0;
	long next_new_node = 0;
	for(long clade=0; clade<Nclades; ++clade){
		if(clade_is_tip[clade]) old2new_clade[clade] = (next_new_tip++);
		else old2new_clade[clade] = Ntips + (next_new_node++);
	}
	for(long edge=0; edge<Nedges; ++edge){
		tree_edge[edge*2+0] = old2new_clade[tree_edge[edge*2+0]];
		tree_edge[edge*2+1] = old2new_clade[tree_edge[edge*2+1]];
	}
	for(long tip=0; tip<extant_tips.size(); ++tip){
		extant_tips[tip] = old2new_clade[extant_tips[tip]];
	}
	root = old2new_clade[root];
		
	// remove extinct tips if needed (make coalescent)
	if(coalescent && ((death_rate_intercept!=0) || (death_rate_factor!=0))){
		std::vector<long> pruning_new_tree_edge, pruning_new2old_clade;
		std::vector<double> pruning_new_edge_length;
		long pruning_new_root, pruning_Ntips_kept, pruning_Nnodes_kept, pruning_Nedges_kept;
		double root_shift;
		get_subtree_with_specific_tips(	Ntips,
										Nnodes,
										Nedges,
										tree_edge,
										edge_length,
										extant_tips,
										true, // collapse monofurcations
										false,
										pruning_new_tree_edge,
										pruning_new_edge_length,
										pruning_new2old_clade,
										pruning_new_root,
										pruning_Ntips_kept,
										pruning_Nnodes_kept,
										pruning_Nedges_kept,
										root_shift);
		tree_edge 	= pruning_new_tree_edge;
		edge_length = pruning_new_edge_length;
		Ntips		= pruning_Ntips_kept;
		Nnodes		= pruning_Nnodes_kept;
		Nclades 	= Ntips + Nnodes;
		Nedges		= pruning_Nedges_kept;
		root 		= pruning_new_root;
		root_time 	+= root_shift; // alternative: = clade2end_time[new2old_clade[pruning_new2old_clade[pruning_new_root]]];
	}
	
	return Rcpp::List::create(	Rcpp::Named("tree_edge") 		= Rcpp::wrap(tree_edge),
								Rcpp::Named("edge_length") 		= Rcpp::wrap(edge_length),
								Rcpp::Named("Nnodes") 			= Nnodes,
								Rcpp::Named("Ntips") 			= Ntips,
								Rcpp::Named("Nedges") 			= Nedges,
								Rcpp::Named("root")				= root, // this is actually guaranteed to be = Ntips
								Rcpp::Named("Nbirths")			= Nbirths,
								Rcpp::Named("Ndeaths")			= Ndeaths,
								Rcpp::Named("root_time")		= root_time,
								Rcpp::Named("final_time")		= time,
								Rcpp::Named("equilibrium_time")	= equilibrium_time,
								Rcpp::Named("birth_times")		= Rcpp::wrap(birth_times),
								Rcpp::Named("death_times")		= Rcpp::wrap(death_times));
}




// Generate a random phylogenetic tree under a speciation/extinction model, where species are born or go extinct as a Poissonian process
// New species are added by splitting one of the currently extant tips (chosen randomly) into Nsplits new tips
// This function is similar to generate_random_tree_CPP(..) above, with one important difference:
//    Per-capita speciation and extinction rates are modelled as Brownian motions evolving along the tree edges (constrained in a finite interval via reflection)
//	  Hence per-capita speciation/extinction rates are scalar traits specific to each node & tip.
// [[Rcpp::export]]
Rcpp::List generate_random_tree_BM_rates_CPP(	const long 	 	max_tips,					// (INPUT) max number of tips (extant tips, if coalescent==true). If <=0, this constraint is ignored.
												const double	max_time,					// (INPUT) max simulation time. If <=0, this constraint is ignored.
												const double	max_time_since_equilibrium,	// (INPUT) max simulation time, counting from the first time point where death_rate-birth_rate changed sign. This may be used as an alternative to (or in conjunction with) max_time to ensure the tree has reached speciation/extinction equilibrium. If <0, this constraint is ignored.
												const double	birth_rate_diffusivity,		// (INPUT) diffusivity of the evolving per-capita birth rate
												const double 	min_birth_rate_pc,			// (INPUT) minimum allowed per-capita birth rate
												const double 	max_birth_rate_pc,			// (INPUT) maximum allowed per-capita birth rate
												const double	death_rate_diffusivity,		// (INPUT) diffusivity of the evolving per-capita death rate
												const double 	min_death_rate_pc,			// (INPUT) minimum allowed per-capita death rate
												const double 	max_death_rate_pc,			// (INPUT) maximum allowed per-capita death rate
												const bool		coalescent,					// (INPUT) whether to return only the coalescent tree (i.e. including only extant tips)
												const long		Nsplits,					// (INPUT) number of children to create at each diversification event. Must be at least 2. For a bifurcating tree this should be set to 2. If >2, then the tree will be multifurcating.
												const bool		as_generations,				// (INPUT) if false, then edge lengths correspond to time. If true, then edge lengths correspond to generations (hence if coalescent==false, all edges will have unit length).
												const bool		include_birth_times,		// (INPUT) if true, then the times of speciations (in order of occurrence) will also be returned
												const bool		include_death_times,		// (INPUT) if true, then the times of extinctions (in order of occurrence) will also be returned
												const bool		include_rates){				// (INPUT) if true, then the per-capita birth & death rates for each clade will also be returned
	const long expected_Nclades = (max_tips<0 ? 2l : max_tips);
	long next_Nsplits = max(2l, Nsplits);
	std::vector<long> tree_edge;
	std::vector<long> extant_tips;
	std::vector<long> clade2parent;
	std::vector<double> clade2end_time;
	std::vector<double> clade2birth_rate_pc; // per-capita birth rate for each node/tip (i.e. determining the waiting time until splitting)
	std::vector<double> clade2death_rate_pc; // per-capita death rate for each node/tip (i.e. determining the waiting time until extinction)
	std::vector<double> birth_times, death_times;
	tree_edge.reserve(expected_Nclades*2);
	extant_tips.reserve(ceil(expected_Nclades/2.0)); 	// keep track of which clades are extant tips, as the tree is built
	clade2parent.reserve(expected_Nclades); 			// keep track of parent of each clade
	clade2end_time.reserve(expected_Nclades); 			// keep track of time at which each clade split or went extinct (negative if clade is an extant tip)
	clade2birth_rate_pc.reserve(expected_Nclades);
	clade2death_rate_pc.reserve(expected_Nclades);
	
	// create the first tip (which is also the root)
	long Ntips 		= 0; 	// current number of extant + extinct tips
	long Nclades 	= 0;	// current number of clades
	long root 		= 0;
	extant_tips.push_back(Nclades++);
	clade2parent.push_back(-1); // root has no parent
	clade2end_time.push_back(-1);
	clade2birth_rate_pc.push_back(uniformWithinInclusiveRight(min_birth_rate_pc,max_birth_rate_pc)); // randomly choose the root's per-capita birth rate
	clade2death_rate_pc.push_back(uniformWithinInclusiveRight(min_death_rate_pc,max_death_rate_pc)); // randomly choose the root's per-capita death rate
	++Ntips;
	
	// create additional tips, by splitting existing tips at each step (turning the split parent tip into a node)
	long Nedges 		= 0;
	long Nbirths		= 0;
	long Ndeaths 		= 0;
	double time 		= 0;
	double total_rate	= INFTY_D;
	double total_birth_rate		= clade2birth_rate_pc[root];
	double total_death_rate		= clade2death_rate_pc[root];
	double equilibrium_time 	= INFTY_D;
	double initial_growth_rate 	= NAN_D; // keep track of the net growth rate (birth rate - death rate) at the beginning of the simulation
	while(((max_tips<=0) || ((coalescent ? extant_tips.size() : Ntips)<max_tips)) && ((max_time<=0) || (time+1/total_rate<max_time)) && ((max_time_since_equilibrium<0) || (time-equilibrium_time+1/total_rate<max_time_since_equilibrium))){
		// determine time of next speciation or extinction event
		// prevent deaths if only one tip is left
		const double restricted_death_rate = (extant_tips.size()<=1 ? 0 : total_death_rate);
		if(std::isnan(initial_growth_rate)) initial_growth_rate = total_birth_rate - restricted_death_rate;
		if(((total_birth_rate - restricted_death_rate)*initial_growth_rate<0) && (equilibrium_time>time)) equilibrium_time = time; // first crossing of equilibrium state, so keep record
		total_rate = total_birth_rate + restricted_death_rate;
		time += random_exponential_distribution(total_rate);
		const bool birth = random_bernoulli(total_birth_rate/(total_birth_rate+restricted_death_rate));
				
		// randomly pick an existing tip to split or kill
		long tip   = random_int_from_distribution(extant_tips, (birth ? clade2birth_rate_pc : clade2death_rate_pc), (birth ? total_birth_rate : total_death_rate));
		long clade = extant_tips[tip];
		clade2end_time[clade] = time;
		const double edge_length = clade2end_time[clade]-(clade==root ? 0 : clade2end_time[clade2parent[clade]]);
		
		// update total birth & death rates for the removal of the chosen clade from the pool of tips
		total_birth_rate -= clade2birth_rate_pc[clade];
		total_death_rate -= clade2death_rate_pc[clade];
				
		if(birth){
			// split chosen tip into Nsplits daughter-tips & create new edges leading into those tips
			if(max_tips>0) next_Nsplits = min(1+max_tips-long(coalescent ? extant_tips.size() : Ntips), max(2l, Nsplits)); // temporarily reduce Ntips to stay within limits
			// child 1:
			++Nedges;
			++Nbirths;
			tree_edge.push_back(clade);
			tree_edge.push_back(Nclades);
			extant_tips[tip] = (Nclades++); // replace the old tip with one of the new ones
			clade2parent.push_back(clade);
			clade2end_time.push_back(-1);
			if(include_birth_times) birth_times.push_back(time);
			// pick per-capita birth & death rates for this splitting point (node) and update total birth & death rates
			clade2birth_rate_pc.push_back(get_next_bounded_BM_sample(birth_rate_diffusivity,min_birth_rate_pc,max_birth_rate_pc,edge_length,clade2birth_rate_pc[clade]));
			clade2death_rate_pc.push_back(get_next_bounded_BM_sample(death_rate_diffusivity,min_death_rate_pc,max_death_rate_pc,edge_length,clade2death_rate_pc[clade]));
			total_birth_rate += clade2birth_rate_pc.back();
			total_death_rate += clade2death_rate_pc.back();
			
			// remaining children:
			for(long c=1; c<next_Nsplits; ++c){
				++Nedges;
				++Ntips;
				tree_edge.push_back(clade);
				tree_edge.push_back(Nclades);
				extant_tips.push_back(Nclades++);
				clade2parent.push_back(clade);
				clade2end_time.push_back(-1);
				// pick per-capita birth & death rates for this splitting point (node) and update total birth & death rates
				clade2birth_rate_pc.push_back(get_next_bounded_BM_sample(birth_rate_diffusivity,min_birth_rate_pc,max_birth_rate_pc,edge_length,clade2birth_rate_pc[clade]));
				clade2death_rate_pc.push_back(get_next_bounded_BM_sample(death_rate_diffusivity,min_death_rate_pc,max_death_rate_pc,edge_length,clade2death_rate_pc[clade]));
				total_birth_rate += clade2birth_rate_pc.back();
				total_death_rate += clade2death_rate_pc.back();
			}
		}else{
			// kill chosen tip (remove from pool of extant tips); note that it still remains a tip, but it can't diversify anymore
			extant_tips[tip] = extant_tips.back();
			extant_tips.pop_back();
			++Ndeaths;
			if(include_death_times) death_times.push_back(time);
		}
		// abort if the user has interrupted the calling R program
		Rcpp::checkUserInterrupt();
	}
	
	// add a small dt at the end to make all edges non-zero length
	time += random_exponential_distribution(total_rate);
	if(max_time>0) time = min(time, max_time); // prevent going past max_time
	if(max_time_since_equilibrium>=0) time = min(time, max_time_since_equilibrium+equilibrium_time); // prevent going past max_time_since_equilibrium

	std::vector<double> edge_length(Nedges);
	if(as_generations){
		// set all edge lengths to 1
		edge_length.assign(Nedges,1);
	}else{
		// calculate edge lengths based on end times
		for(long edge=0, child; edge<Nedges; ++edge){
			child = tree_edge[edge*2+1];
			if(clade2end_time[child]>=0) edge_length[edge] = clade2end_time[child] - clade2end_time[tree_edge[edge*2+0]];
			else edge_length[edge] = time - clade2end_time[tree_edge[edge*2+0]];
		}	
	}
	double root_time = clade2end_time[root]; // root_time = time at which the root split
	
	// identify tips as the clades with no outgoing edge
	std::vector<bool> clade_is_tip(Nclades,true);
	for(long edge=0; edge<Nedges; ++edge){
		clade_is_tip[tree_edge[edge*2+0]] = false;
	}
	
	// re-number tip & node indices to conform with the phylo format, where tips are indexed first (0,..,Ntips-1) and nodes last (Ntips,..,Ntips+Nnodes-1).
	long Nnodes = Nclades - Ntips;
	std::vector<long> old2new_clade(Nclades,-1);
	long next_new_tip  = 0;
	long next_new_node = 0;
	for(long clade=0; clade<Nclades; ++clade){
		if(clade_is_tip[clade]) old2new_clade[clade] = (next_new_tip++);
		else old2new_clade[clade] = Ntips + (next_new_node++);
	}
	for(long edge=0; edge<Nedges; ++edge){
		tree_edge[edge*2+0] = old2new_clade[tree_edge[edge*2+0]];
		tree_edge[edge*2+1] = old2new_clade[tree_edge[edge*2+1]];
	}
	for(long tip=0; tip<extant_tips.size(); ++tip){
		extant_tips[tip] = old2new_clade[extant_tips[tip]];
	}
	root = old2new_clade[root];
	
	// also update pc birth & death rates to new indices if needed
	std::vector<double> birth_rates_pc, death_rates_pc;
	if(include_rates){
		birth_rates_pc.resize(Nclades);
		death_rates_pc.resize(Nclades);
		for(long c=0; c<Nclades; ++c){
			birth_rates_pc[old2new_clade[c]] = clade2birth_rate_pc[c];
			death_rates_pc[old2new_clade[c]] = clade2death_rate_pc[c];
		}
	}
		
	// remove extinct tips if needed (make coalescent)
	if(coalescent && (max_death_rate_pc>0)){
		std::vector<long> pruning_new_tree_edge, pruning_new2old_clade;
		std::vector<double> pruning_new_edge_length;
		long pruning_new_root, pruning_Ntips_kept, pruning_Nnodes_kept, pruning_Nedges_kept;
		double root_shift;
		get_subtree_with_specific_tips(	Ntips,
										Nnodes,
										Nedges,
										tree_edge,
										edge_length,
										extant_tips,
										true, // collapse monofurcations
										false,
										pruning_new_tree_edge,
										pruning_new_edge_length,
										pruning_new2old_clade,
										pruning_new_root,
										pruning_Ntips_kept,
										pruning_Nnodes_kept,
										pruning_Nedges_kept,
										root_shift);
		tree_edge 	= pruning_new_tree_edge;
		edge_length = pruning_new_edge_length;
		Ntips		= pruning_Ntips_kept;
		Nnodes		= pruning_Nnodes_kept;
		Nclades 	= Ntips + Nnodes;
		Nedges		= pruning_Nedges_kept;
		root 		= pruning_new_root;
		root_time 	+= root_shift; // alternative: = clade2end_time[new2old_clade[pruning_new2old_clade[pruning_new_root]]];
		
		// also update pc birth & death rates to new indices if needed
		if(include_rates){
			std::vector<double> old_birth_rates_pc(birth_rates_pc);
			std::vector<double> old_death_rates_pc(death_rates_pc);
			birth_rates_pc.resize(Nclades);
			death_rates_pc.resize(Nclades);
			for(long c=0; c<Nclades; ++c){
				birth_rates_pc[c] = old_birth_rates_pc[pruning_new2old_clade[c]];
				death_rates_pc[c] = old_death_rates_pc[pruning_new2old_clade[c]];
			}
		}
	}
	
	return Rcpp::List::create(	Rcpp::Named("tree_edge") 		= Rcpp::wrap(tree_edge),
								Rcpp::Named("edge_length") 		= Rcpp::wrap(edge_length),
								Rcpp::Named("Nnodes") 			= Nnodes,
								Rcpp::Named("Ntips") 			= Ntips,
								Rcpp::Named("Nedges") 			= Nedges,
								Rcpp::Named("root")				= root, // this is actually guaranteed to be = Ntips
								Rcpp::Named("Nbirths")			= Nbirths,
								Rcpp::Named("Ndeaths")			= Ndeaths,
								Rcpp::Named("root_time")		= root_time,
								Rcpp::Named("final_time")		= time,
								Rcpp::Named("equilibrium_time")	= equilibrium_time,
								Rcpp::Named("birth_times")		= Rcpp::wrap(birth_times),
								Rcpp::Named("death_times")		= Rcpp::wrap(death_times),
								Rcpp::Named("birth_rates_pc")	= Rcpp::wrap(birth_rates_pc),	// birth_rates_pc[c] will be the per-capita birth rate of clade c (prior to its extinction or splitting)
								Rcpp::Named("death_rates_pc")	= Rcpp::wrap(death_rates_pc));	// death_rates_pc[c] will be the per-capita death rate of clade c (prior to its extinction or splitting)
}




#pragma mark -
#pragma mark Tree dating
#pragma mark 


// given a phylogenetic tree and a vector in [0,1]^Nnodes ("relative ages"), map relative ages (R) to absolute ages (A = distance from present) for each node
// The mapping proceeds from root to tips
// The mapping is defined as follows:
//		For any clade, let C1,C2,..,Cn be the sequence of clades from the root to that clade
//		Let L[C_n] := min(max_abs_node_ages[C_n], A[C_{n-1}])
//		Then A[C_n] = L[C_n] + R[C_n]*(min_abs_node_ages[C_n]-L[C_n])
// Tips are implicitly assumed to have age 0.
// This mapping may be useful for fitting absolute node ages while searching within classical box-constraints
// [[Rcpp::export]]
std::vector<double> relative_to_absolute_node_ages_CPP(	const long 					Ntips,
														const long 					Nnodes,
														const long 					Nedges,
														const std::vector<long>		&tree_edge,				// (INPUT) 2D array of size Nedges x 2, in row-major format
														const std::vector<long>		&traversal_queue,		// (INPUT) traversal queue from root-->tips (listing clade indices), not including tips
														const std::vector<double>	&relative_node_ages,	// (INPUT) 1D array of size Nnodes, listing relative ages of each node. Values should be within [0,1].
														const std::vector<double>	&min_abs_node_ages,		// (INPUT) 1D array of size Nnodes
														const std::vector<double>	&max_abs_node_ages){	// (INPUT) 1D array of size Nnodes
	// determine parent clade for each clade
	std::vector<long> clade2parent;
	get_parent_per_clade(Ntips, Nnodes, Nedges, tree_edge, clade2parent);
	
	const long root_node = traversal_queue[0] - Ntips;
	
	std::vector<double> abs_node_ages(Nnodes);
	abs_node_ages[root_node] = max_abs_node_ages[root_node] + relative_node_ages[root_node]*(min_abs_node_ages[root_node] - min_abs_node_ages[root_node]);
	for(long q=1, clade, node, pnode; q<traversal_queue.size(); ++q){
		clade 	= traversal_queue[q];
		node 	= clade-Ntips;
		pnode 	= clade2parent[clade]-Ntips;
		const double L = min(max_abs_node_ages[node], abs_node_ages[pnode]);
		abs_node_ages[node] = L + relative_node_ages[node]*(min_abs_node_ages[node] - L);
	}
	return abs_node_ages;
}


// [[Rcpp::export]]
std::vector<double> propagate_min_ages_upstream_CPP(const long 					Ntips,
													const long 					Nnodes,
													const long 					Nedges,
													const std::vector<long>		&tree_edge,
													const std::vector<long>		&traversal_queue,		// (INPUT) traversal queue from root-->tips (listing clade indices), not including tips
													const std::vector<long>		&anchor_nodes,
													const std::vector<long>		&anchor_min_ages){
	// determine parent clade for each clade
	std::vector<long> clade2parent;
	get_parent_per_clade(Ntips, Nnodes, Nedges, tree_edge, clade2parent);

	std::vector<double> min_node_ages(Nnodes,0);
	for(long a=0; a<anchor_nodes.size(); ++a){
		min_node_ages[anchor_nodes[a]] = anchor_min_ages[a];
	}
	
	// propagate min_ages upstream
	for(long q=traversal_queue.size()-1, clade, node, pnode; q>=1; --q){
		clade	= traversal_queue[q];
		node	= clade - Ntips;
		pnode 	= clade2parent[clade] - Ntips;
		min_node_ages[pnode] = max(min_node_ages[pnode], min_node_ages[node]);
	}
	
	return min_node_ages;
}


// [[Rcpp::export]]
std::vector<double> propagate_max_ages_downstream_CPP(	const long 					Ntips,
														const long 					Nnodes,
														const long 					Nedges,
														const std::vector<long>		&tree_edge,
														const std::vector<long>		&traversal_queue,		// (INPUT) traversal queue from root-->tips (listing clade indices), not including tips
														const std::vector<long>		&anchor_nodes,
														const std::vector<long>		&anchor_max_ages){
	// determine parent clade for each clade
	std::vector<long> clade2parent;
	get_parent_per_clade(Ntips, Nnodes, Nedges, tree_edge, clade2parent);

	std::vector<double> max_node_ages(Nnodes,INFTY_D);
	for(long a=0; a<anchor_nodes.size(); ++a){
		max_node_ages[anchor_nodes[a]] = anchor_max_ages[a];
	}
	
	// propagate max_ages downstream
	for(long q=0, clade, node, pnode; q<traversal_queue.size(); ++q){
		clade	= traversal_queue[q];
		node	= clade - Ntips;
		pnode 	= clade2parent[clade] - Ntips;
		max_node_ages[pnode] = min(max_node_ages[pnode], max_node_ages[node]);
	}
	
	return max_node_ages;
}



#pragma mark -
#pragma mark Tree statistics
#pragma mark 





// Count the number of tips descending from each node
// The tree must be rooted; the root should be the unique node with no parent
void get_total_tip_count_per_node(	const long			Ntips,
									const long 			Nnodes,
									const long			Nedges,
									const IntegerVector &tree_edge, 			// (INPUT) 2D array (in row-major format) of size Nedges x 2
									std::vector<long>	&node2total_tip_count){	// (OUTPUT) array of size Nnodes, with each entry being the total number of tips descending from the node
	long clade;

	// determine parent clade for each clade
	std::vector<long> clade2parent;
	get_parent_per_clade(Ntips, Nnodes, Nedges, tree_edge, clade2parent);

	// find root using the mapping clade2parent
	const long root = get_root_from_clade2parent(Ntips, clade2parent);

	// get tree traversal route (tips --> root)
	// traversal_queue[] will be of size Nclades, and will have entries in 0:(Nclades-1)
	std::vector<long> traversal_queue, traversal_node2first_edge, traversal_node2last_edge, traversal_edges;
	get_tree_traversal_root_to_tips(	Ntips,
										Nnodes,
										Nedges,
										root,
										tree_edge,
										true,
										false,
										traversal_queue,
										traversal_node2first_edge,
										traversal_node2last_edge,
										traversal_edges,
										false,
										"");

	// calculate number of descending tips per node, traversing tips-->root (excluding the root)
	node2total_tip_count.assign(Nnodes,0);
	for(long q=traversal_queue.size()-1; q>=1; --q){
		clade = traversal_queue[q];
		node2total_tip_count[clade2parent[clade]-Ntips] += (clade<Ntips ? 1 : node2total_tip_count[clade-Ntips]);
	}
}


// Count the number of tips descending from each node
// This is an Rcpp wrapper function for get_total_tip_count_per_node(..)
// Requirements:
//   The tree must be rooted; the root should be the unique node with no parent
//   The tree can include multifurcations as well as monofurcations
// [[Rcpp::export]]
IntegerVector get_total_tip_count_per_node_CPP(	const long			Ntips,
												const long 			Nnodes,
												const long			Nedges,
												const IntegerVector &tree_edge){	// (INPUT) 2D array (in row-major format) of size Nedges x 2
	std::vector<long> node2total_tip_count;
	get_total_tip_count_per_node(	Ntips,
									Nnodes,
									Nedges,
									tree_edge, 
									node2total_tip_count);
	return Rcpp::wrap(node2total_tip_count);
}


// Pick random subsets of tips from a tree, by traversing from root-->tips and at each node choosing randomly between children
// The size of each random subset is Nrandoms, the number of independent subsets is Nsubsets
// [[Rcpp::export]]
std::vector<long> pick_random_tips_CPP(	const long			Ntips,
										const long 			Nnodes,
										const long			Nedges,
										const IntegerVector &tree_edge, 			// (INPUT) 2D array (in row-major format) of size Nedges x 2, or an empty std::vector (no tree available).
										const long			Nrandoms,				// (INPUT) number of random tips to pick at each experiment (i.e. in each independent subset)
										const long			Nsubsets,				// (INPUT) number of times the experiment should be repeated, i.e. each time drawing Nrandoms tips anew.
										const bool			with_replacement){		// (INPUT) pick tips with replacement. If false, then children with no descending tips left to pick from, are excluded from traversal; all other children of a node remain equally probable.
	
	const long Nclades = Ntips+Nnodes;
	long tip, clade;
	if((!with_replacement) && (Nrandoms>Ntips)) return std::vector<long>(); // this should not happen
	
	// determine parent clade for each clade
	std::vector<long> clade2parent;
	get_parent_per_clade(Ntips, Nnodes, Nedges, tree_edge, clade2parent);

	// find root using the mapping clade2parent
	const long root = get_root_from_clade2parent(Ntips, clade2parent);
		
	// prepare access structures for random tip selection
	std::vector<long> node2first_child, node2last_child, children, node2total_tip_count;
	get_children_per_node(	Ntips,
							Nnodes,
							Nedges,
							root,
							tree_edge,
							node2first_child,
							node2last_child,
							children);
	if(!with_replacement){
		get_total_tip_count_per_node(	Ntips,
										Nnodes,
										Nedges,
										tree_edge,
										node2total_tip_count);
	}
	std::vector<long> 	tips_remaining_per_node;
	std::vector<double> clade2weight;
	
	// pick Nrepeats random tip subsets of size Nrandoms
	std::vector<long> random_tips(Nsubsets*Nrandoms);
	for(long s=0; s<Nsubsets; ++s){
		if(!with_replacement){
			// re-initialize counters for this repeat
			tips_remaining_per_node = node2total_tip_count;
			clade2weight.resize(Nclades);
			for(clade=0; clade<Nclades; ++clade){
				clade2weight[clade] = (clade<Ntips ? 1.0 : (tips_remaining_per_node[clade-Ntips]>0 ? 1 : 0));
			}
		}
		for(long t=0; t<Nrandoms; ++t){
			if(with_replacement){
				tip = get_tip_by_random_uniform_traversal(	Ntips,
															root,
															node2first_child,
															node2last_child,
															children);
			}else{
				tip = get_tip_by_random_traversal(	Ntips,
													root,
													node2first_child,
													node2last_child,
													children,
													clade2weight);
				clade2weight[tip] = 0; // prevent re-inclusion of this tip in the future (i.e. don't replace)
				// propagate information upwards
				clade = tip;
				while(clade!=root){
					clade = clade2parent[clade];
					tips_remaining_per_node[clade-Ntips] -= 1;
					if(tips_remaining_per_node[clade-Ntips]<=0){
						// no more tips to draw from this clade, so set weight to zero
						clade2weight[clade] = 0.0;
					}
				}					
			}
			random_tips[s*Nrandoms + t] = tip;
		}
		// abort if the user has interrupted the calling R program
		Rcpp::checkUserInterrupt();
	}
	return random_tips;
}




// Calculate sum of all branch lengths, for each subtree in a tree
// This is equivalent to the 'phylogenetic diversity' measure introduced by Faith (1992).
// References: 
//    Faith (1992). Conservation evaluation and phylogenetic diversity. Biological Conservation 61:1-10.
//    Mark Vellend et al. Measuring phylogenetic biodiversity. Table 14.2 (presence/absence based phylogenetic diversity)
void get_cumulative_edge_lengths_per_node(	const long			Ntips,
											const long 			Nnodes,
											const long			Nedges,
											const long 			root, 				// (INPUT) index of root node, i.e. an integer in 0:(Ntips+Nnodes-1)
											const IntegerVector &tree_edge, 		// (INPUT) 2D array (in row-major format) of size Nedges x 2
											const NumericVector &edge_length, 		// (INPUT) 1D array of size Nedges, or an empty std::vector (all branches have length 1)
											std::vector<double>	&node2CBL){			// (OUTPUT) array of size Nnodes, with each entry being the cumulative branch length (phylogenetic diversity) of the subtree rooted in that node.
	const long Nclades = Ntips+Nnodes;
	long clade;

	// determine parent clade for each clade
	std::vector<long> clade2parent;
	get_parent_per_clade(Ntips, Nnodes, Nedges, tree_edge, clade2parent);
	
	// determine incoming edge per clade
	std::vector<long> incoming_edge_per_clade(Nclades,-1);
	for(long edge=0; edge<Nedges; ++edge){
		incoming_edge_per_clade[tree_edge[edge*2+1]] = edge;
	}

	// get tree traversal route (tips --> root)
	// traversal_queue[] will be of size Nclades, and will have entries in 0:(Nclades-1)
	std::vector<long> traversal_queue, traversal_node2first_edge, traversal_node2last_edge, traversal_edges;
	get_tree_traversal_root_to_tips(	Ntips,
										Nnodes,
										Nedges,
										root,
										tree_edge,
										true,
										false,
										traversal_queue,
										traversal_node2first_edge,
										traversal_node2last_edge,
										traversal_edges,
										false,
										"");
	reverse_array(traversal_queue); // make tips-->roots

	// calculate phylogenetic diversities (cumulative branch lengths), traversing tips-->root
	node2CBL.assign(Nnodes,0);
	for(long q=0; q<traversal_queue.size(); ++q){
		clade  = traversal_queue[q];
		if(clade==root) continue;
		node2CBL[clade2parent[clade]-Ntips] += (clade<Ntips ? 0.0 : node2CBL[clade-Ntips]) + (edge_length.size()==0 ? 1.0 : edge_length[incoming_edge_per_clade[clade]]);
	}
}



// Calculate the maximum & minimum distance of any tip to the root
// [[Rcpp::export]]
Rcpp::List get_min_max_tip_distance_from_root_CPP(	const long 			Ntips,
													const long 			Nnodes,
													const long 			Nedges,
													const IntegerVector &tree_edge,		// (INPUT) 2D array of size Nedges x 2 in row-major format
													const NumericVector &edge_length){	// (INPUT) 1D array of size Nedges, or an empty std::vector (all branches have length 1)
	const long Nclades = Ntips + Nnodes;
	long parent, clade;
										
	// determine root
	const long root = get_root_clade(Ntips, Nnodes, Nedges, tree_edge);
	
	// get tree traversal route (root --> tips)											
	std::vector<long> traversal_queue, node2first_edge, node2last_edge, edge_mapping;
	get_tree_traversal_root_to_tips(	Ntips,
										Nnodes,
										Nedges,
										root,
										tree_edge,
										true,
										false,
										traversal_queue,
										node2first_edge,
										node2last_edge,
										edge_mapping,
										false,
										"");
	
	// determine incoming edge per clade
	std::vector<long> incoming_edge_per_clade(Nclades,-1);
	for(long edge=0; edge<Nedges; ++edge){
		incoming_edge_per_clade[tree_edge[edge*2+1]] = edge;
	}
										
	// calculate distance of each node to its nearest descending tip
	// (traverse tips --> root)
	std::vector<double> min_tip_distance_per_node(Nnodes,INFTY_D), max_tip_distance_per_node(Nnodes,0);
	double min_distance, max_distance;
	for(long q=traversal_queue.size()-1; q>=0; --q){
		clade = traversal_queue[q];
		if(clade==root) continue;
		parent = tree_edge[incoming_edge_per_clade[clade]*2 + 0];
		min_distance = (clade<Ntips ? 0.0 : min_tip_distance_per_node[clade-Ntips]) + (edge_length.size()==0 ? 1.0 : edge_length[incoming_edge_per_clade[clade]]);
		max_distance = (clade<Ntips ? 0.0 : max_tip_distance_per_node[clade-Ntips]) + (edge_length.size()==0 ? 1.0 : edge_length[incoming_edge_per_clade[clade]]);
		min_tip_distance_per_node[parent-Ntips]	= min(min_tip_distance_per_node[parent-Ntips], min_distance);
		max_tip_distance_per_node[parent-Ntips]	= max(max_tip_distance_per_node[parent-Ntips], max_distance);
	}

	return Rcpp::List::create(	Rcpp::Named("min_distance")	= Rcpp::wrap(min_tip_distance_per_node[root-Ntips]),
								Rcpp::Named("max_distance")	= Rcpp::wrap(max_tip_distance_per_node[root-Ntips]));
}





// calculate distance from root, for each clade (tips+nodes)
// distance from root = cumulative branch length from root to the clade
void get_distances_from_root(	const long 			Ntips,
								const long 			Nnodes,
								const long 			Nedges,
								const IntegerVector &tree_edge,		// (INPUT) 2D array of size Nedges x 2 in row-major format
								const NumericVector &edge_length, 	// (INPUT) 1D array of size Nedges, or an empty std::vector (all branches have length 1)
								std::vector<double>	&distances){	// (OUTPUT) 1D array of size Nclades, listing the phylogenetic distance of each clade from the root
	const long Nclades = Ntips + Nnodes;
	long parent, clade;
										
	// determine parent clade for each clade
	std::vector<long> clade2parent;
	get_parent_per_clade(Ntips, Nnodes, Nedges, tree_edge, clade2parent);

	// find root using the mapping clade2parent
	const long root = get_root_from_clade2parent(Ntips, clade2parent);
	
	// get tree traversal route (root --> tips)											
	std::vector<long> traversal_queue, traversal_node2first_edge, traversal_node2last_edge, traversal_edges;
	get_tree_traversal_root_to_tips(	Ntips,
										Nnodes,
										Nedges,
										root,
										tree_edge,
										true,
										false,
										traversal_queue,
										traversal_node2first_edge,
										traversal_node2last_edge,
										traversal_edges,
										false,
										"");
	
	// determine incoming edge per clade
	std::vector<long> incoming_edge_per_clade(Nclades,-1);
	for(long edge=0; edge<Nedges; ++edge){
		incoming_edge_per_clade[tree_edge[edge*2+1]] = edge;
	}
										
	// calculate distance from root for each clade
	// (traverse root --> tips)
	distances.resize(Nclades);
	distances[root] = 0;
	for(long q=0; q<traversal_queue.size(); ++q){
		clade = traversal_queue[q];
		if(clade==root) continue;
		parent = clade2parent[clade];
		distances[clade] = (edge_length.size()==0 ? 1.0 : edge_length[incoming_edge_per_clade[clade]]) + distances[parent];
	}
}



// calculate distance from root, for each clade (tips+nodes)
// distance from root = cumulative branch length from root to the clade
// This is an Rcpp wrapper for the function get_distances_from_root(..)
// [[Rcpp::export]]
NumericVector get_distances_from_root_CPP(	const long 			Ntips,
											const long 			Nnodes,
											const long 			Nedges,
											const IntegerVector &tree_edge,			// (INPUT) 2D array of size Nedges x 2 in row-major format
											const NumericVector &edge_length){ 		// (INPUT) 1D array of size Nedges, or an empty std::vector (all branches have length 1)
	std::vector<double> distances;
	get_distances_from_root(Ntips,
							Nnodes,
							Nedges,
							tree_edge,
							edge_length,
							distances);
	return Rcpp::wrap(distances);
}









// For each clade (tip & node) in a tree, find the closest tip (in terms of cumulative branch length).
// Optionally, the search can be restricted to descending tips.
// Optionally, the search can also be restricted to a subset of target tips.
// If you want distances in terms of branch counts (instead of cumulative branch lengths), simply provide an empty edge_length[].
// Requirements:
//   The input tree must be rooted (root will be determined automatically, as the node that has no incoming edge)
//   The input tree can be multifurcating and/or monofurcating
// [[Rcpp::export]]
Rcpp::List get_closest_tip_per_clade_CPP(	const long 			Ntips,
											const long 			Nnodes,
											const long 			Nedges,
											const IntegerVector &tree_edge,				// 2D array of size Nedges x 2 in row-major format
											const NumericVector &edge_length, 			// 1D array of size Nedges, or an empty std::vector (all branches have length 1)
											const IntegerVector	&onlyToTips,			// 1D array listing target tips to restrict search to, or an empty std::vector (consider all tips as targets)
											bool				only_descending_tips,	// if true, then for each clade only descending tips are considered for nearest-distance. If false, some clades may have non-descending tips assigned as nearest tips.
											bool 				verbose,
											const std::string	&verbose_prefix){
	const long Nclades = Ntips + Nnodes;
	long parent, clade, tip, incoming_edge;
	double candidate_distance;

	// determine parent clade for each clade
	std::vector<long> clade2parent;
	get_parent_per_clade(Ntips, Nnodes, Nedges, tree_edge, clade2parent);
	
	// determine incoming edge per clade
	std::vector<long> incoming_edge_per_clade(Nclades,-1);
	for(long edge=0; edge<Nedges; ++edge){
		incoming_edge_per_clade[tree_edge[edge*2+1]] = edge;
	}
	
	// find root using the mapping clade2parent
	const long root = get_root_from_clade2parent(Ntips, clade2parent);
	
	// get tree traversal route (root --> tips)											
	std::vector<long> traversal_queue_root2tips, traversal_node2first_edge, traversal_node2last_edge, traversal_edges;
	get_tree_traversal_root_to_tips(	Ntips,
										Nnodes,
										Nedges,
										root,
										tree_edge,
										true,
										false,
										traversal_queue_root2tips,
										traversal_node2first_edge,
										traversal_node2last_edge,
										traversal_edges,
										verbose,
										verbose_prefix);

	// Step 1: calculate nearest descending tip per clade (traverse tips --> root)
	std::vector<long> nearest_descending_tip_per_clade(Nclades,-1);
	std::vector<double> distance_to_nearest_descending_tip_per_clade(Nclades,INFTY_D);
	if(onlyToTips.size()==0){
		// consider all tips as potential targets
		for(long tip=0; tip<Ntips; ++tip){
			nearest_descending_tip_per_clade[tip] = tip;
			distance_to_nearest_descending_tip_per_clade[tip] = 0;
		}
	}else{
		// only consider provided tips as targets
		for(long t=0; t<onlyToTips.size(); ++t){
			tip = onlyToTips[t];
			nearest_descending_tip_per_clade[tip] = tip;
			distance_to_nearest_descending_tip_per_clade[tip] = 0;
		}
	}
	for(long q=traversal_queue_root2tips.size()-1; q>=0; --q){
		clade = traversal_queue_root2tips[q];
		if(clade==root) continue;
		if(nearest_descending_tip_per_clade[clade]<0) continue; // no descending tip available from this clade
		parent			= clade2parent[clade];
		incoming_edge	= incoming_edge_per_clade[clade];
		// propagate information about nearest descending tip, to parent (if better than already saved for the parent)
		candidate_distance = (edge_length.size()==0 ? 1.0 : edge_length[incoming_edge]) + distance_to_nearest_descending_tip_per_clade[clade];
		if(candidate_distance<distance_to_nearest_descending_tip_per_clade[parent]){
			distance_to_nearest_descending_tip_per_clade[parent] = candidate_distance;
			nearest_descending_tip_per_clade[parent] = nearest_descending_tip_per_clade[clade];
		}
	}
	
	if(only_descending_tips){
		// only descending tips allowed, so we're finished
		return Rcpp::List::create(	Rcpp::Named("nearest_tips") 		= Rcpp::wrap(nearest_descending_tip_per_clade),
									Rcpp::Named("nearest_distances") 	= Rcpp::wrap(distance_to_nearest_descending_tip_per_clade));
	}
	
	// Step 2: calculate nearest tip per clade, regardless of whether descending or not (traverse root --> tips)
	std::vector<long> nearest_tip_per_clade(Nclades);
	std::vector<double> distance_to_nearest_tip_per_clade(Nclades);
	nearest_tip_per_clade[root] = nearest_descending_tip_per_clade[root];
	distance_to_nearest_tip_per_clade[root] = distance_to_nearest_descending_tip_per_clade[root];
	for(long q=0; q<traversal_queue_root2tips.size(); ++q){
		clade = traversal_queue_root2tips[q];
		if(clade==root) continue;
		parent				= clade2parent[clade];
		incoming_edge 		= incoming_edge_per_clade[clade];
		candidate_distance 	= (edge_length.size()==0 ? 1.0 : edge_length[incoming_edge]) + distance_to_nearest_tip_per_clade[parent];
		if(candidate_distance<distance_to_nearest_descending_tip_per_clade[clade]){
			// it's shorter to go upwards, than downwards
			distance_to_nearest_tip_per_clade[clade] = candidate_distance;
			nearest_tip_per_clade[clade] = nearest_tip_per_clade[parent];
		}else{
			// nearest descending tip is also nearest tip
			distance_to_nearest_tip_per_clade[clade] = distance_to_nearest_descending_tip_per_clade[clade];
			nearest_tip_per_clade[clade] = nearest_descending_tip_per_clade[clade];
		}
	}
	
	return Rcpp::List::create(	Rcpp::Named("nearest_tips") 		= Rcpp::wrap(nearest_tip_per_clade),
								Rcpp::Named("nearest_distances") 	= Rcpp::wrap(distance_to_nearest_tip_per_clade));
}





// Calculate phylogenetic distance matrix between all pairs of focal_clades
// Distance = cumulative branch length of both clades back to their most recent common ancestor (aka "patristic distance")
// This function is slightly different from get_distances_between_clades_CPP(), in that here the distances between all possible clade pairs are returned.
// The time complexity is O(Nfocals*Nfocals*Nanc + Ntips), where Nanc is the average number of ancestors per tip.
// Requirements:
//   The input tree must be rooted (root will be determined automatically, as the node that has no incoming edge)
//   The input tree can include multifurcations and monofurcations
// Attention: 0-based indexing is used for input and output variables, so make sure to shift indices in R before calling this function

// [[Rcpp::export]]
NumericMatrix get_distance_matrix_between_clades_CPP(	const long 			Ntips,
														const long 			Nnodes,
														const long 			Nedges,
														const IntegerVector &tree_edge,			// 2D array of size Nedges x 2 in row-major format
														const NumericVector &edge_length, 		// 1D array of size Nedges, or an empty std::vector (all branches have length 1)
														const IntegerVector &focal_clades,		// 1D array of size Nfocals, containing values in 0:(Nclades-1). These will correspond to the rows & columns of the returned distance matrix.
														bool 				verbose,
														const std::string 	&verbose_prefix){
	const long Nclades = Ntips + Nnodes;
	const long Nfocals = focal_clades.size();
	long parent, clade;

	// determine parent clade for each clade
	std::vector<long> clade2parent;
	get_parent_per_clade(Ntips, Nnodes, Nedges, tree_edge, clade2parent);
	
	// determine incoming edge per clade
	std::vector<long> incoming_edge_per_clade(Nclades,-1);
	for(long edge=0; edge<Nedges; ++edge){
		incoming_edge_per_clade[tree_edge[edge*2+1]] = edge;
	}
	
	// find root using the mapping clade2parent
	const long root = get_root_from_clade2parent(Ntips, clade2parent);
	
	// get tree traversal route (root --> tips)											
	std::vector<long> traversal_queue, traversal_node2first_edge, traversal_node2last_edge, traversal_edges;
	get_tree_traversal_root_to_tips(	Ntips,
										Nnodes,
										Nedges,
										root,
										tree_edge,
										true,
										false,
										traversal_queue,
										traversal_node2first_edge,
										traversal_node2last_edge,
										traversal_edges,
										verbose,
										verbose_prefix);
															
	// calculate number of ancestors and distance from root for each clade
	// (traverse root --> tips)
	std::vector<long> ancestor_count_per_clade(Nclades);
	std::vector<double> distance_from_root_per_clade(Nclades);
	ancestor_count_per_clade[root] 		= 0;
	distance_from_root_per_clade[root] 	= 0;
	for(long q=0; q<traversal_queue.size(); ++q){
		clade = traversal_queue[q];
		if(clade==root) continue;
		parent = clade2parent[clade];
		ancestor_count_per_clade[clade] 	= 1 + ancestor_count_per_clade[parent];
		distance_from_root_per_clade[clade]	=  (edge_length.size()==0 ? 1.0 : edge_length[incoming_edge_per_clade[clade]]) + distance_from_root_per_clade[parent];
	}
	const long total_ancestor_count = vector_sum(ancestor_count_per_clade);
	
	// calculate ancestry for each clade in a long list ancestors[]
	// (traverse root --> tips)
	std::vector<long> clade2first_ancestor(Nclades); // for each clade c, ancestors[clade2first_ancestor[c]..clade2last_ancestor[c]] will be the list of ancestor clades leading to the clade c
	std::vector<long> clade2last_ancestor(Nclades);
	clade2first_ancestor[0] = 0;
	clade2last_ancestor[0] = clade2first_ancestor[0] + ancestor_count_per_clade[0] - 1;
	for(long c=1; c<Nclades; ++c){
		clade2first_ancestor[c] = clade2last_ancestor[c-1] + 1;
		clade2last_ancestor[c]  = clade2first_ancestor[c] + ancestor_count_per_clade[c] - 1;
	}
	std::vector<long> ancestors(total_ancestor_count);
	for(long q=0; q<traversal_queue.size(); ++q){
		clade = traversal_queue[q];
		if(clade==root) continue;
		parent = clade2parent[clade];
		// step 1: copy the parent's ancestry to the child's ancestry
		for(long a=clade2first_ancestor[parent]; a<=clade2last_ancestor[parent]; ++a){
			ancestors[clade2first_ancestor[clade]+(a-clade2first_ancestor[parent])] = ancestors[a];
		}
		// step 2: append the parent to the clade's ancestry
		ancestors[clade2last_ancestor[clade]] = parent;
	}
	
	// calculate most-recent-common-ancestor and phylogenetic distance for each focal clade pair
	long cladeA, cladeB, mrca;
	NumericMatrix distances(Nfocals,Nfocals);
	for(long i=0; i<Nfocals; ++i){
		cladeA = focal_clades[i];
		distances(i,i) = 0;
		for(long j=i+1; j<Nfocals; ++j){
			cladeB = focal_clades[j];
			// step 1: determine most recent common ancestor
			// check for trivial case
			if(cladeA==cladeB){
				mrca = cladeA;
			}else{
				// follow ancestry of both clades in synchrony, until they diverge
				// note that the first ancestor of every clade will be the root
				long a,b;
				for(a=clade2first_ancestor[cladeA], b=clade2first_ancestor[cladeB]; (a<=clade2last_ancestor[cladeA]) && (b<=clade2last_ancestor[cladeB]); ++a, ++b){
					if(ancestors[a]!=ancestors[b]) break;
					else mrca = ancestors[a];
				}
				// check special case where one clade is descendant of the other (this would lead to a "premature" stop of the above loop)
				if((a<=clade2last_ancestor[cladeA]) && (ancestors[a]==cladeB)){
					mrca = cladeB;
				}else if((b<=clade2last_ancestor[cladeB]) && (ancestors[b]==cladeA)){
					mrca = cladeA;
				}
			}
			// step 2: calculate distance
			distances(i,j) = distance_from_root_per_clade[cladeA] + distance_from_root_per_clade[cladeB] - 2*distance_from_root_per_clade[mrca];
			distances(j,i) = distances(i,j);
		}
	}
	
	return(distances);
}
												


// Calculate phylogenetic distance for pairs of clades (cladesA[] vs cladesB[])
// Distance = cumulative branch length of both clades back to their most recent common ancestor (aka "patristic distance")
// There's some initial overhead involved with this function, but for large number of clade pairs this becomes more efficient
// Time complexity is O(Ntips + Npairs*log(Ntips)).
// Returns a NumericVector of size Npairs, with each entry being the distance between the two clades
// This function is slightly different from get_distance_matrix_between_clades_CPP(), in that here only distances between specific clade pairs are returned.
// Requirements:
//   The input tree must be rooted (root will be determined automatically, as the node that has no incoming edge)
//   The input tree can be multifurcating and/or monofurcating
// Attention: 0-based indexing is used for input and output variables, so make sure to shift indices in R before and after calling this function

// [[Rcpp::export]]
NumericVector get_distances_between_clades_CPP(	const long 			Ntips,
												const long 			Nnodes,
												const long 			Nedges,
												const IntegerVector &tree_edge,			// 2D array of size Nedges x 2 in row-major format
												const NumericVector &edge_length, 		// 1D array of size Nedges, or an empty std::vector (all branches have length 1)
												const IntegerVector &cladesA,			// 1D array of size Npairs, containing values in 0:(Nclades-1)
												const IntegerVector	&cladesB,			// 1D array of size Npairs, containing values in 0:(Nclades-1)
												bool 				verbose,
												const std::string 	&verbose_prefix){
	const long Npairs = cladesA.size();
	const long Nclades = Ntips + Nnodes;
	long parent, clade;

	// determine parent clade for each clade
	std::vector<long> clade2parent;
	get_parent_per_clade(Ntips, Nnodes, Nedges, tree_edge, clade2parent);
	
	// determine incoming edge per clade
	std::vector<long> incoming_edge_per_clade(Nclades,-1);
	for(long edge=0; edge<Nedges; ++edge){
		incoming_edge_per_clade[tree_edge[edge*2+1]] = edge;
	}
	
	// find root using the mapping clade2parent
	const long root = get_root_from_clade2parent(Ntips, clade2parent);
	
	// get tree traversal route (root --> tips)											
	std::vector<long> traversal_queue, traversal_node2first_edge, traversal_node2last_edge, traversal_edges;
	get_tree_traversal_root_to_tips(	Ntips,
										Nnodes,
										Nedges,
										root,
										tree_edge,
										true,
										false,
										traversal_queue,
										traversal_node2first_edge,
										traversal_node2last_edge,
										traversal_edges,
										verbose,
										verbose_prefix);
															
	// calculate number of ancestors and distance from root for each clade
	// (traverse root --> tips)
	std::vector<long> ancestor_count_per_clade(Nclades);
	std::vector<double> distance_from_root_per_clade(Nclades);
	ancestor_count_per_clade[root] 		= 0;
	distance_from_root_per_clade[root] 	= 0;
	for(long q=0; q<traversal_queue.size(); ++q){
		clade = traversal_queue[q];
		if(clade==root) continue;
		parent = clade2parent[clade];
		ancestor_count_per_clade[clade] 	= 1 + ancestor_count_per_clade[parent];
		distance_from_root_per_clade[clade]	=  (edge_length.size()==0 ? 1.0 : edge_length[incoming_edge_per_clade[clade]]) + distance_from_root_per_clade[parent];
	}
	const long total_ancestor_count = vector_sum(ancestor_count_per_clade);
	
	// calculate ancestry for each clade in a long list ancestors[]
	// (traverse root --> tips)
	std::vector<long> clade2first_ancestor(Nclades); // for each clade c, ancestors[clade2first_ancestor[c]..clade2last_ancestor[c]] will be the list of ancestor clades leading to the clade c
	std::vector<long> clade2last_ancestor(Nclades);
	clade2first_ancestor[0] = 0;
	clade2last_ancestor[0] = clade2first_ancestor[0] + ancestor_count_per_clade[0] - 1;
	for(long c=1; c<Nclades; ++c){
		clade2first_ancestor[c] = clade2last_ancestor[c-1] + 1;
		clade2last_ancestor[c]  = clade2first_ancestor[c] + ancestor_count_per_clade[c] - 1;
	}
	std::vector<long> ancestors(total_ancestor_count);
	for(long q=0; q<traversal_queue.size(); ++q){
		clade = traversal_queue[q];
		if(clade==root) continue;
		parent = clade2parent[clade];
		// step 1: copy the parent's ancestry to the child's ancestry
		for(long a=clade2first_ancestor[parent]; a<=clade2last_ancestor[parent]; ++a){
			ancestors[clade2first_ancestor[clade]+(a-clade2first_ancestor[parent])] = ancestors[a];
		}
		// step 2: append the parent to the clade's ancestry
		ancestors[clade2last_ancestor[clade]] = parent;
	}
	
	// calculate most-recent-common-ancestor for each clade pair
	std::vector<long> mrca_per_pair(Npairs);
	long cladeA, cladeB;
	for(long p=0; p<Npairs; ++p){
		cladeA = cladesA[p];
		cladeB = cladesB[p];
		// check for trivial case
		if(cladeA==cladeB){
			mrca_per_pair[p] = cladeA;
			continue;
		}
		// follow ancestry of both clades in synchrony, until they diverge
		// note that the first ancestor of every clade will be the root
		long a, b, mrca=-1;
		for(a=clade2first_ancestor[cladeA], b=clade2first_ancestor[cladeB]; (a<=clade2last_ancestor[cladeA]) && (b<=clade2last_ancestor[cladeB]); ++a, ++b){
			if(ancestors[a]!=ancestors[b]) break;
			else mrca = ancestors[a];
		}
		// check special case where one clade is descendant of the other (this would lead to a "premature" stop of the above loop)
		if((a<=clade2last_ancestor[cladeA]) && (ancestors[a]==cladeB)){
			mrca = cladeB;
		}else if((b<=clade2last_ancestor[cladeB]) && (ancestors[b]==cladeA)){
			mrca = cladeA;
		}
		mrca_per_pair[p] = mrca;
	}
	
	// calculate pairwise distances
	NumericVector distances(Npairs);
	for(long p=0; p<Npairs; ++p){
		distances[p] = distance_from_root_per_clade[cladesA[p]] + distance_from_root_per_clade[cladesB[p]] - 2*distance_from_root_per_clade[mrca_per_pair[p]];
	}
	return(distances);
}




// Count the number of extant clades per time point (depth point), where time points are taken on a regular grid
// The tree need not be ultrametric, although in general this function only makes sense for ultrametric trees (e.g. where edge lengths are time intervals)
// [[Rcpp::export]]
Rcpp::List count_clades_at_regular_times_CPP(	const long 			Ntips,
												const long 			Nnodes,
												const long 			Nedges,
												const IntegerVector	&tree_edge,		// (INPUT) 2D array of size Nedges x 2, flattened in row-major format
												const NumericVector	&edge_length, 	// (INPUT) 1D array of size Nedges, or an empty std::vector (all branches have length 1)
												const long			Ntimes,			// (INPUT) number of time points
												const bool			include_slopes){	// (INPUT) if true, slopes of the clades_per_time_point curve are also returned	
	// calculate clade distances from root
	const NumericVector clade_times = get_distances_from_root_CPP(	Ntips,
																	Nnodes,
																	Nedges,
																	tree_edge,
																	edge_length);
	const double max_time = get_array_max(clade_times);
	const double min_time = 0;
	
	// determine distance bins
	const double time_step = (1.0-1e-7)*(max_time-min_time)/(Ntimes-1);
	std::vector<double> time_points(Ntimes);
	for(long t=0; t<Ntimes; ++t){
		time_points[t] = min_time + time_step*t;
	}
	
	// calculate number of clades within each time point
	std::vector<long> diversities(Ntimes,0);
	for(long edge=0, child, parent; edge<Nedges; ++edge){
		parent = tree_edge[edge*2+0];
		child  = tree_edge[edge*2+1];
		const long last_time_point 	= max(0L,min(Ntimes-1,long(floor((clade_times[child]-min_time)/time_step))));
		const long first_time_point = (parent<0 ? last_time_point : max(0L,min(Ntimes-1,long(ceil((clade_times[parent]-min_time)/time_step)))));
		if(first_time_point==last_time_point){ ++diversities[first_time_point]; }
		else{ for(long t=first_time_point; t<=last_time_point; ++t) ++diversities[t]; }
	}
		
	// calculate slopes (symmetric difference coefficient)
	// use one-sided differences at the edges
	std::vector<double> slopes(include_slopes ? Ntimes : 0);
	std::vector<double> relative_slopes(include_slopes ? Ntimes : 0);
	if(include_slopes && (Ntimes>=2)){
		for(long t=0, left, right; t<Ntimes; ++t){
			left  				= max(t-1,0L);
			right 				= min(t+1,Ntimes-1);
			const double dt		= (time_points[right]-time_points[left]);
			const double CC		= ((left<t && t<right) ? (diversities[left]+diversities[t]+diversities[right])/3.0 : (diversities[left]+diversities[right])/2.0);
			slopes[t] 			= (diversities[right]-diversities[left])/dt;
			relative_slopes[t] 	= (CC==0 ? NAN_D : slopes[t]/CC);
		}
	}
	
	return Rcpp::List::create(	Rcpp::Named("time_points") 		= Rcpp::wrap(time_points),
								Rcpp::Named("diversities") 		= Rcpp::wrap(diversities),
								Rcpp::Named("slopes") 			= Rcpp::wrap(slopes),
								Rcpp::Named("relative_slopes") 	= Rcpp::wrap(relative_slopes));
}



// Count number of extant clades at arbitrary time points
// [[Rcpp::export]]
IntegerVector count_clades_at_times_CPP(const long			Ntips,
										const long 			Nnodes,
										const long			Nedges,
										IntegerVector 		tree_edge,			// (INPUT) 2D array (in row-major format) of size Nedges x 2
										const NumericVector	&edge_length,		// (INPUT) 1D array of size Nedges, or empty
										const NumericVector	&times){			// (INPUT) 1D array of size Ntimes
	const long Nclades = Ntips + Nnodes;

	// determine parent clade for each clade
	std::vector<long> clade2parent;
	get_parent_per_clade(Ntips, Nnodes, Nedges, tree_edge, clade2parent);

	// find root using the mapping clade2parent
	const long root = get_root_from_clade2parent(Ntips, clade2parent);	
	
	// calculate distances from root
	std::vector<double> distances_from_root(Nclades);
	get_distances_from_root(Ntips,
							Nnodes,
							Nedges,
							tree_edge,
							edge_length,
							distances_from_root);
	
	const long Ntimes = times.size();
	std::vector<long> diversities(Ntimes,0);
	for(long t=0; t<Ntimes; ++t){
		for(long clade=0; clade<Nclades; ++clade){
			if(clade==root) continue;
			if((distances_from_root[clade]>=times[t]) && (distances_from_root[clade2parent[clade]]<=times[t])){
				diversities[t] += 1;
			}
		}
	}
	return(Rcpp::wrap(diversities));
}



// returns true if the tree includes nodes that have more than 2 children
template<class ARRAY_TYPE>
bool tree_has_multifurcations(	const long			Ntips,
								const long 			Nnodes,
								const long			Nedges,
								const ARRAY_TYPE	&tree_edge){	// (INPUT) 2D array (in row-major format) of size Nedges x 2
	std::vector<long> child_count_per_node(Nnodes,0);
	for(long edge=0; edge<Nedges; ++edge) ++child_count_per_node[tree_edge[2*edge+0]-Ntips];
	for(long node=0; node<Nnodes; ++node){
		if(child_count_per_node[node]>2) return true;
	}
	return false;
}


template<class ARRAY_TYPE>
void count_monofurcations_and_multifurcations(	const long			Ntips,
												const long 			Nnodes,
												const long			Nedges,
												const ARRAY_TYPE	&tree_edge,				// (INPUT) 2D array (in row-major format) of size Nedges x 2
												long				&Nmonofurcations,		// (OUTPUT) number of monofurcating nodes
												long				&Nbifurcations,			// (OUTPUT) number of bifurcating nodes
												long				&Nmultifurcations){		// (OUTPUT) number of multifurcating nodes
	std::vector<long> child_count_per_node(Nnodes,0);
	for(long edge=0; edge<Nedges; ++edge) ++child_count_per_node[tree_edge[2*edge+0]-Ntips];
	Nmonofurcations = Nbifurcations = Nmultifurcations = 0;
	for(long node=0; node<Nnodes; ++node){
		if(child_count_per_node[node]==1) ++Nmonofurcations;
		else if(child_count_per_node[node]==2) ++Nbifurcations;
		else ++Nmultifurcations;
	}
}




// Extract speciation events and extinction events (with associated waiting times) from tree
//   Speciation event = non-root node
//   Extinction event = non-crown tip, i.e. tip that is not at maximum distance from the root
// The returned info may be used for fitting birth-death tree-generation models
// [[Rcpp::export]]
Rcpp::List get_speciation_extinction_events_CPP(const long				Ntips,
												const long 				Nnodes,
												const long				Nedges,
												IntegerVector 			tree_edge,			// (INPUT) 2D array (in row-major format) of size Nedges x 2
												const NumericVector		&edge_length,		// (INPUT) 1D array of size Nedges, or empty
												const double			min_age,			// (INPUT) min phylogenetic distance from the tree crown, to be considered. If <=0, this constraint is ignored.
												const double			max_age,			// (INPUT) max phylogenetic distance from the tree crown, to be considered. If <=0, this constraint is ignored.
												const IntegerVector		&only_clades,		// (INPUT) optional list of clade indices to consider. Can also be empty, in which case no filtering is done.
												const IntegerVector		&omit_clades){		// (INPUT) optional list of clade indices to omit. Can also be empty.
	const long Nclades = Ntips + Nnodes;
	long clade, parent;
										
	// determine parent clade for each clade
	std::vector<long> clade2parent;
	get_parent_per_clade(Ntips, Nnodes, Nedges, tree_edge, clade2parent);

	// find root using the mapping clade2parent
	const long root = get_root_from_clade2parent(Ntips, clade2parent);

	// get incoming edge for each clade
	std::vector<long> incoming_edge_per_clade;
	get_incoming_edge_per_clade(Ntips, Nnodes, Nedges, tree_edge, incoming_edge_per_clade);

	// get tree traversal route (tips --> root)											
	std::vector<long> queue_root2tips, node2first_edge, node2last_edge, edge_mapping;
	get_tree_traversal_root_to_tips(Ntips,
									Nnodes,
									Nedges,
									root,
									tree_edge,
									true,
									false,
									queue_root2tips,
									node2first_edge,
									node2last_edge,
									edge_mapping,
									false,
									"");
								
	// calculate distance from root for each clade (traverse root --> tips)
	std::vector<double> distances_from_root(Nclades);
	distances_from_root[root] = 0;
	for(long q=0; q<queue_root2tips.size(); ++q){
		clade = queue_root2tips[q];
		if(clade==root) continue;
		parent = clade2parent[clade];
		distances_from_root[clade] = (edge_length.size()==0 ? 1.0 : edge_length[incoming_edge_per_clade[clade]]) + distances_from_root[parent];
	}
	
	// sort clades in chronological order (increasing distance from root)
	std::vector<long> chronological_clade_order;
	qsortIndices(distances_from_root, chronological_clade_order);
	const double max_distance_from_root = distances_from_root[chronological_clade_order.back()];
	
	// count number of extant species at the time of each clade
	std::vector<long> diversity_at_clade(Nclades);
	long current_diversity = 1;
	for(long c=0; c<Nclades; ++c){
		clade = chronological_clade_order[c];
		diversity_at_clade[clade] = current_diversity;
		if(clade<Ntips) --current_diversity; // a tip marks an extinction event
		else current_diversity += node2last_edge[clade-Ntips] - node2first_edge[clade-Ntips]; // a node marks a speciation event, generating N-1 new species (where N is the number of children)
	}
	
	// figure out which clades to consider
	std::vector<bool> include_clade(Nclades,(only_clades.size()==0));
	for(long c=0; c<only_clades.size(); ++c) include_clade[only_clades[c]] = true;
	for(long c=0; c<omit_clades.size(); ++c) include_clade[omit_clades[c]] = false;
	if(max_age>0){
		// only consider clades above a certain distance from the root
		for(long c=0; c<Nclades; ++c){
			clade = chronological_clade_order[c];
			include_clade[clade] = include_clade[clade] && ((max_distance_from_root-distances_from_root[clade])<=max_age);
		}
	}
	if(min_age>0){
		// only consider clades below a certain distance from the root
		for(long c=0; c<Nclades; ++c){
			clade = chronological_clade_order[c];
			include_clade[clade] = include_clade[clade] && ((max_distance_from_root-distances_from_root[clade])>=min_age);
		}
	}
	
	// figure out number of speciation events
	long Nspeciations = 0;
	long Npoints = 0;
	for(long c=0; c<Nclades; ++c){
		clade = chronological_clade_order[c];
		if(include_clade[clade]) ++Npoints;
		if((clade>=Ntips) && include_clade[clade]) ++Nspeciations;
	}
	if(Nspeciations>0) --Nspeciations;  // the first event should be discarded, because it has unknown waiting time
	
	// extract speciation events
	// speciation_waiting_times[event] = waiting time to speciation event (counted from previous considered speciation event)
	// speciation_diversities[event] = number of extant clades during the speciation event
	std::vector<double> speciation_waiting_times(Nspeciations);
	std::vector<long> 	speciation_diversities(Nspeciations);
	std::vector<long> 	speciation_clades(Nspeciations);
	std::vector<double> speciation_times(Nspeciations);
	double previous_time = -1; // negative means undefined
	for(long c=0, event=0; c<Nclades; ++c){
		clade = chronological_clade_order[c];
		if(!include_clade[clade]) continue; // omit this clade
		if(clade<Ntips) continue; // not a speciation event
		if(previous_time<0){
			// this is the first speciation event encountered, so just record its time but don't include in returned events
			previous_time = distances_from_root[clade];
			continue;
		}
		speciation_diversities[event]	= diversity_at_clade[clade];
		speciation_clades[event]		= clade;
		speciation_times[event]			= distances_from_root[clade];
		speciation_waiting_times[event]	= speciation_times[event] - previous_time;
		previous_time 					= speciation_times[event];
		++event;
	}
	
	// extract extinction events (non-crown tips mark extinctions)
	std::vector<long> extinction_tips;
	previous_time = -1;
	for(long c=0; c<Nclades; ++c){
		clade = chronological_clade_order[c];
		if((clade>=Ntips) || (distances_from_root[clade]>=(1.0-RELATIVE_EPSILON)*max_distance_from_root)) continue; // non-crown tips correspond to extinction events
		if(previous_time<0){
			// this is the first extinction event encountered, so just record its time but don't include in returned events
			previous_time = distances_from_root[clade];
			continue;
		}
		extinction_tips.push_back(clade);
	}
	const long Nextinctions = extinction_tips.size();
	std::vector<double> extinction_waiting_times(Nextinctions);
	std::vector<long> 	extinction_diversities(Nextinctions);
	std::vector<double> extinction_times(Nextinctions);
	for(long event=0; event<Nextinctions; ++event){
		clade 							= extinction_tips[event];
		extinction_diversities[event]	= diversity_at_clade[clade];
		extinction_times[event]			= distances_from_root[clade];
		extinction_waiting_times[event]	= extinction_times[event] - previous_time;
		previous_time 					= extinction_times[event];
	}
	
	// all considered points (tips & nodes) in chronological order
	std::vector<double> times(Npoints);
	std::vector<long>   diversities(Npoints), clades(Npoints);
	for(long c=0, point=0; c<Nclades; ++c){
		clade = chronological_clade_order[c];
		if(!include_clade[clade]) continue;
		times[point] 		= distances_from_root[clade];
		diversities[point]	= diversity_at_clade[clade];
		clades[point]		= clade;
		++point;
	}
		
	return Rcpp::List::create(	Rcpp::Named("Nspeciations")				= Nspeciations,							// number of speciation events included
								Rcpp::Named("speciation_waiting_times")	= Rcpp::wrap(speciation_waiting_times),	// waiting time until each speciation event (distance from the previous speciation event)
								Rcpp::Named("speciation_diversities") 	= Rcpp::wrap(speciation_diversities), 	// number of clades just prior to the speciation
								Rcpp::Named("speciation_times")			= Rcpp::wrap(speciation_times),			// time at the speciation event
								Rcpp::Named("speciation_clades")		= Rcpp::wrap(speciation_clades),		// clade marking each speciation event
								Rcpp::Named("Nextinctions")				= Nextinctions,							// number of extinction events included
								Rcpp::Named("extinction_waiting_times")	= Rcpp::wrap(extinction_waiting_times),	// waiting time until each extinction event (distance from the previous extinction event)
								Rcpp::Named("extinction_diversities") 	= Rcpp::wrap(extinction_diversities), 	// number of clades just prior to the extinction
								Rcpp::Named("extinction_times")			= Rcpp::wrap(extinction_times),			// time at the extinction event
								Rcpp::Named("extinction_tips")			= Rcpp::wrap(extinction_tips),			// tip marking each extinction event
								Rcpp::Named("max_distance_from_root")	= max_distance_from_root,				// maximum distance of any tip to the root. Note that this may be outside of the span of times[]
								Rcpp::Named("times")					= Rcpp::wrap(times),					// time at each clade, in chronological order (root is first)
								Rcpp::Named("clades")					= Rcpp::wrap(clades),					// clade associated with each time point
								Rcpp::Named("diversities")				= Rcpp::wrap(diversities));				// number of clades just prior to the splitting or extinction of each clade, in chronological order (root is first)
}






#pragma mark -
#pragma mark Writing trees to file
#pragma mark 


// convert a tree to a string in Newick (parenthetic) format
// If the tree is not rooted, it is first rooted 
// [[Rcpp::export]]
std::string tree_to_Newick_string_CPP(	const long			Ntips,
										const long 			Nnodes,
										const long			Nedges,
										IntegerVector 		tree_edge,			// (INPUT) 2D array (in row-major format) of size Nedges x 2
										const NumericVector	&edge_length,		// (INPUT) 1D array of size Nedges, or empty
										const StringVector	&tip_labels,		// (INPUT) 1D array of size Ntips, or empty
										const StringVector	&node_labels,		// (INPUT) 1D array of size Nnodes, or empty
										const long			digits,				// (INPUT) number of digits used for printing edge lengths
										const double		root_edge_length){	// (INPUT) optional edge length leading into the root. Not really an edge of the tree. Ignored if negative.
	const bool has_tip_labels  	= (tip_labels.size()>0);
	const bool has_node_labels  = (node_labels.size()>0);
	const bool has_edge_lengths	= (edge_length.size()>0);
	const long Nclades = Ntips + Nnodes;
	long child,node,clade;
	ostringstream output;
	output << std::setprecision(digits);

	// get incoming edge for each clade
	std::vector<long> incoming_edge_per_clade;
	get_incoming_edge_per_clade(Ntips, Nnodes, Nedges, tree_edge, incoming_edge_per_clade);

	// find root based on incoming_edge_per_clade
	// make sure there is exactly one root
	long root = -1;
	for(long clade=Ntips; clade<Nclades; ++clade){
		if(incoming_edge_per_clade[clade]<0){
			if(root>=0){
				// already encountered a root, so abort and re-root properly
				root = -1;
				break;
			}else{
				root = clade;
			}
		}
	}
	
	if(root<0){
		// tree is not rooted, so root at some node
		root = Ntips;
		root_tree_at_node(	Ntips,
							Nnodes,
							Nedges,
							tree_edge,	// will be modified in-situ
							root);
							
		// re-calculate incoming edge for each clade
		get_incoming_edge_per_clade(Ntips, Nnodes, Nedges, tree_edge, incoming_edge_per_clade);
	}

	// get edge mappings (for efficient tree traversal)
	std::vector<long> node2first_edge, node2last_edge, edge_mapping;
	get_node2edge_mappings(	Ntips,
							Nnodes,
							Nedges,
							tree_edge,
							node2first_edge,
							node2last_edge,
							edge_mapping);
	
	// create output queue in depth-first-search direction
	// also count how many brackets to close at each clade
	// use a scratch_stack for traversing nodes
	// note that the returned tree will actually be reverse version (i.e. tips-->root, in reversed debth-first-search)
	std::vector<long> scratch_stack, queue;
	std::vector<long> Nbrackets_to_close(Nclades,0);
	std::vector<bool> is_first_child(Nclades,false);
	scratch_stack.reserve(floor(2*log(Ntips)/log(2.0))); // rough estimate of typical tree depth x 2. scratch_stack may be resized along the way if needed.
	scratch_stack.push_back(root);
	queue.reserve(Nclades);
	is_first_child[root] = true;
	while(scratch_stack.size()>0){
		clade = scratch_stack.back();
		scratch_stack.pop_back();
		queue.push_back(clade);
		if(clade>=Ntips){
			node = clade - Ntips;
			for(long e=node2last_edge[node]; e>=node2first_edge[node]; --e){
				child = tree_edge[edge_mapping[e]*2+1];
				scratch_stack.push_back(child); // add child to stack for further exploration in the next iteration
				if(e==node2last_edge[node]){
					// delegate braket closing for this clade to last descending tip
					Nbrackets_to_close[child] = 1 + Nbrackets_to_close[clade];
					Nbrackets_to_close[clade] = 0;
				}
				is_first_child[child] = (e==node2first_edge[node]);
			}
		}
	}
		
	// traverse output queue in reverse direction
	for(long q=queue.size()-1; q>=0; --q){
		clade = queue[q];
		for(long b=0; b<Nbrackets_to_close[clade]; ++b) output << "(";
		if(clade<Ntips){
			if(has_tip_labels) output << tip_labels[clade];
		}else{
			node = clade - Ntips;
			output << ")";
			if(has_node_labels) output << node_labels[node];
		}
		if(has_edge_lengths && (clade!=root)) output << ":" << edge_length[incoming_edge_per_clade[clade]];
		if(!is_first_child[clade]) output << ",";
	}
	if(has_edge_lengths && (root_edge_length>=0)) output << ":" << root_edge_length;
	output << ";";

	return(output.str());
}



#pragma mark -
#pragma mark Statistics of trait distribution
#pragma mark 




// Auxiliary function to get_trait_depth_consenTRAIT_CPP()
// Assumes a pre-calculated traversal root (root-->tips)
void aux_get_trait_depth_consenTRAIT(	const long 					Ntips,
										const long 					Nnodes,
										const long 					Nedges,
										const long					root,				// integer in Ntips:(Nclades-1)
										const IntegerVector 		&tree_edge,			// 2D array of size Nedges x 2, flattened in row-major format
										const NumericVector 		&edge_length, 		// 1D array of size Nedges, or an empty std::vector (all branches have length 1)
										const std::vector<long>		&state_per_tip,		// 1D std::vector of integer states of size Ntips. <=0 means absence, >0 means presence.
										const double				threshold_fraction,	// minimum fraction of tips in a clade that must share the trait, in order for the clade to be counted towards tau_D. In the original paper by Martiny et al (2013), this was 0.9.
										const bool					count_singletons,	// if true, then singleton tips (i.e. having the trait but not counted towards any positive clade) are included by assuming a non-discovered sister tip having the trait with probability 0.5, and diverging just after the parent node [Martiny et al 2013].
										const bool					weighted,			// if true, then positive clades (i.e. counted towards tauD) are weighted according to the number of positive tips
										const std::vector<long>		&traversal_queue,				// (INPUT) 1D array of size Nclades, with values in 0:(Nclades-1). Traversal queue root-->tips. Generated using the function get_tree_traversal_root_to_tips(include_tips=true).
										const std::vector<long>		&traversal_node2first_edge,		// (INPUT) 1D array of size Nnodes, with values in 0:(Nedges-1). Generated using the function get_tree_traversal_root_to_tips().
										const std::vector<long>		&traversal_node2last_edge,		// (INPUT) 1D array of size Nnodes, with values in 0:(Nedges-1). Generated using the function get_tree_traversal_root_to_tips().
										const std::vector<long>		&traversal_edges,				// (INPUT) 1D array of size Nedges, with values in 0:(Nedges-1). Generated using the function get_tree_traversal_root_to_tips().
										const std::vector<long> 	&clade2parent,					// (INPUT) 1D array of size Nclades, with values in 0:(Nclades-1).
										const std::vector<long>		&incoming_edge_per_clade,		// (INPUT) 1D array of size Nclades, with values in 0:(Nedges-1).
										const std::vector<double>	&node2max_tip_distance,			// (INPUT) 1D array of size Nnodes, specifying the max phylogenetic distance for each node to its tips.
										const double				singleton_threshold,			// (INPUT) phylogenetic distance threshold for counting a node as a single tip if its max distance to its tips is equal to or below this threshold. For example, if this is 0, then nodes whose descendants are all identical, will be considered as singletons.
										const bool					include_list_of_positives,		// (INPUT) if true, a list of the clades found to be positive is returned.
										std::vector<long>			&tips_per_clade,				// (OUTPUT) 1D array of size Nclades, listing the number of tips descending from each clade
										std::vector<long>			&positives_per_clade,			// (OUTPUT) 1D array of size Nclades, listing the number of positive tips descending from each clade.
										std::vector<double> 		&mean_depth_per_clade,			// (OUTPUT) 1D array of size Nclades, listing the mean phylogenetic depth of each clade
										double						&mean_depth,					// (OUTPUT) mean clade depth at which trait is conserved. This is the original tau_D introduced by Martiny et al. (2013).
										double						&var_depth,						// (OUTPUT) variance of clade depth at which trait is conserved.
										double						&min_depth,						// (OUTPUT) minimum clade depth at which trait is conserved.
										double						&max_depth,						// (OUTPUT) maximum clade depth at which trait is conserved.
										long						&Npositives,					// (OUTPUT) number of positive clades counted towards the tauD statistic
										std::vector<long>			&positive_clades,				// (OUTPUT) clade indices that were found to be positive in the trait. Optional output (see option include_list_of_positives). Will include only those clades that are also included in the statistics; in particular, singletons are only included if count_singletons==true.
										bool 						verbose,
										const std::string			&verbose_prefix){
	const long Nclades = Ntips + Nnodes;
	long clade, parent, child, node, incoming_edge;
	positive_clades.clear();

	// count number of tips with the trait ("positives"), for each clade
	// also count cumulative phylogenetic distances to tips, for each clade
	tips_per_clade.assign(Nclades, 0);
	positives_per_clade.assign(Nclades, 0);
	mean_depth_per_clade.assign(Nclades, 0); // initially this will be the cumulative depth_per_clade, normalization is done at the end
	for(long tip=0; tip<Ntips; ++tip){
		tips_per_clade[tip] = 1;
		positives_per_clade[tip] = (state_per_tip[tip]<=0 ? 0 : 1);
	}
	// traverse tips-->root
	for(long q=traversal_queue.size()-1; q>=0; --q){
		clade = traversal_queue[q];
		if(clade==root) continue;
		parent 			= clade2parent[clade];
		incoming_edge 	= incoming_edge_per_clade[clade];
		positives_per_clade[parent] 	+= positives_per_clade[clade];
		tips_per_clade[parent] 			+= tips_per_clade[clade];
		mean_depth_per_clade[parent] 	+= mean_depth_per_clade[clade] + tips_per_clade[clade] * (edge_length.size()==0 ? 1 : edge_length[incoming_edge]);
	}
	for(clade=0; clade<Nclades; ++clade){
		mean_depth_per_clade[clade] /= tips_per_clade[clade];
	}
	
	// traverse through all nodes for which "almost all" tips share the trait (i.e. "positive" nodes), and calculate their mean depth (--> tau_D statistic)
	// traverse root-->tips, and whenever a node is counted as positive, also mark its sub-clades as counted (without actually counting them)
	double sum_depths 		= 0;
	double sum_sqr_depths 	= 0;
	double total_weight		= 0;
	min_depth 				= NAN_D;
	max_depth 				= NAN_D;
	Npositives		 		= 0;
	std::vector<long> clade_counted(Nclades, false);
	for(long q=0; q<traversal_queue.size(); ++q){
		clade = traversal_queue[q];
		const double fraction_positives = positives_per_clade[clade]/double(tips_per_clade[clade]);
		if((!clade_counted[clade]) && (fraction_positives>=threshold_fraction)){
			if((clade<Ntips) || (node2max_tip_distance[clade-Ntips]<=singleton_threshold)){
				// clade is a singleton, so treat in a special way (or omit)
				// singleton tip sensu Martiny et al (2013): tip having the trait, but not counted towards any "positive" node (because too few of it's neighbors have the trait)
				// according to Martiny, such cases may occur due to undersampling of the tree.
				// this is a modified version of Martiny's: Some nodes may also be counted as singletons, if their descendants are very similar to each other
				if(count_singletons){
					++Npositives;
					const double temp_depth = 0.5*(edge_length.size()==0 ? 1 : edge_length[incoming_edge_per_clade[clade]]);
					total_weight	+= 1;
					sum_depths 		+= temp_depth;
					sum_sqr_depths 	+= SQR(temp_depth);
					if(isnan(min_depth) || (min_depth>temp_depth)) min_depth = temp_depth;
					if(isnan(max_depth) || (max_depth<temp_depth)) max_depth = temp_depth;
					clade_counted[clade] = true;
					positive_clades.push_back(clade);
				}
			}else{
				// clade is a node with phylogenetic diversity above the singleton threshold
				++Npositives;
				const double temp_depth = mean_depth_per_clade[clade];
				const double weight = (weighted ? positives_per_clade[clade] : 1);
				total_weight	+= weight;
				sum_depths 		+= weight * temp_depth;
				sum_sqr_depths 	+= weight * SQR(temp_depth);
				if(isnan(min_depth) || (min_depth>temp_depth)) min_depth = temp_depth;
				if(isnan(max_depth) || (max_depth<temp_depth)) max_depth = temp_depth;
				clade_counted[clade] = true;
				positive_clades.push_back(clade);
			}
		}
		if(clade_counted[clade] && (clade>=Ntips)){
			// clade was counted (either because it's positive, or due to an ancestral positive node)
			// so mark its children as counted as well (the state of being counted will propagate all the way to the tips)
			node = clade - Ntips;
			for(long e=traversal_node2first_edge[node]; e<=traversal_node2last_edge[node]; ++e){
				child = tree_edge[traversal_edges[e]*2+1];
				clade_counted[child] = true;
			}
		}
	}
	
	mean_depth = (total_weight==0 ? NAN_D : sum_depths/total_weight);
	var_depth  = (total_weight==0 ? NAN_D : (Npositives==1 ? 0.0 : (sum_sqr_depths/total_weight - SQR(mean_depth))));
}




// Calculate phylogenetic depth of a binary trait (presence/absence) on a tree
// Reference: Martiny et al (2013). Phylogenetic conservatism of functional traits in microorganisms. ISME Journal. 7:830-838
// consenTRAIT: Consensus Analysis of Phylogenetic Trait Distribution
//
// Input: A phylogenetic tree, and the states of a binary trait on all tips of the tree.
// Output: Mean depth at which the trait varies ("trait depth").
// P-value is probability that random tauD (with randomly re-assigned states) would lead to an equal or greater tauD than observed. Traits are reassigned based on the empirical distribution of presence/absences.
// The time complexity of this routine is O(Nedges).
//
// Requirements:
//   Tree must be rooted (root will be determined automatically based on the edges)
//   Tree may include monofurcations or multifurcations.
// [[Rcpp::export]]
Rcpp::List get_trait_depth_consenTRAIT_CPP(	const long 			Ntips,
											const long 			Nnodes,
											const long 			Nedges,
											const IntegerVector &tree_edge,				// 2D array of size Nedges x 2, flattened in row-major format
											const NumericVector &edge_length, 			// 1D array of size Nedges, or an empty std::vector (all branches have length 1)
											const IntegerVector	&state_per_tip,			// 1D std::vector of integer states of size Ntips. <=0 means absence, >0 means presence.
											const double		threshold_fraction,		// minimum fraction of tips in a clade that must share the trait, in order for the clade to be counted towards tau_D. In the original paper by Martiny et al (2013), this was 0.9.
											const bool			count_singletons,		// if true, then singleton tips (i.e. having the trait but not counted towards any positive clade) are included by assuming a non-discovered sister tip having the trait with probability 0.5, and diverging just after the parent node [Martiny et al 2013]. If false, leftover singletons will be ignored. This may be useful if you suspect many of them to be false positives.
											const bool			weighted,				// if true, then positive clades (i.e. counted towards tauD) are weighted according to the number of positive tips
											const double		singleton_threshold,	// phylogenetic distance threshold (max distance from node to its tips) for counting a node as a single tip if its diversity is equal to or below this threshold. For example, if this is 0, then nodes whose descendants are all identical, will be considered as singletons.
											const long			Npermutations,			// number of random permutations (re-assignments of states) for estimating P-values
											bool 				verbose,
											const std::string 	&verbose_prefix){
	const long Nclades = Ntips + Nnodes;
	long clade, parent;
	
	// determine parent clade for each clade
	std::vector<long> clade2parent;
	get_parent_per_clade(Ntips, Nnodes, Nedges, tree_edge, clade2parent);

	// determine incoming edge per clade
	std::vector<long> incoming_edge_per_clade(Nclades,-1);
	for(long edge=0; edge<Nedges; ++edge){
		incoming_edge_per_clade[tree_edge[edge*2+1]] = edge;
	}	
	
	// find root using the mapping clade2parent
	const long root = get_root_from_clade2parent(Ntips, clade2parent);
	
	// get tree traversal route (root --> tips)											
	std::vector<long> traversal_queue, traversal_node2first_edge, traversal_node2last_edge, traversal_edges;
	get_tree_traversal_root_to_tips(	Ntips,
										Nnodes,
										Nedges,
										root,
										tree_edge,
										true,
										false,
										traversal_queue,
										traversal_node2first_edge,
										traversal_node2last_edge,
										traversal_edges,
										verbose,
										verbose_prefix);
										
	// calculate max phylogenetic distance for each node to any of its tips, traversing tips-->root
	// needed to identify nodes that are essentially tips (i.e. all of its descendants are closely related)
	std::vector<double> node2max_tip_distance(Nnodes,0);
	for(long q=traversal_queue.size()-1; q>=0; --q){
		clade  = traversal_queue[q];
		if(clade==root) continue;
		parent = clade2parent[clade];
		node2max_tip_distance[parent-Ntips] = max(node2max_tip_distance[parent-Ntips], (clade<Ntips ? 0.0 : node2max_tip_distance[clade-Ntips]) + (edge_length.size()==0 ? 1.0 : edge_length[incoming_edge_per_clade[clade]]));
	}
				
	// get observed tauD
	std::vector<long> current_state_per_tip = Rcpp::as< std::vector<long> >(state_per_tip);
	std::vector<long> positive_clades, positives_per_clade, tips_per_clade;
	std::vector<double> mean_depth_per_clade;
	double tauD, varD, minD, maxD;
	long Npositives;
	aux_get_trait_depth_consenTRAIT(Ntips,
									Nnodes,
									Nedges,
									root,
									tree_edge,
									edge_length,
									current_state_per_tip,
									threshold_fraction,
									count_singletons,
									weighted,
									traversal_queue,
									traversal_node2first_edge,
									traversal_node2last_edge,
									traversal_edges,
									clade2parent,
									incoming_edge_per_clade,
									node2max_tip_distance,
									singleton_threshold,
									true,
									tips_per_clade,
									positives_per_clade,
									mean_depth_per_clade,
									tauD,
									varD,
									minD,
									maxD,
									Npositives,
									positive_clades,
									verbose,
									verbose_prefix);
														
	// estimate P-value of observed tauD, by randomly re-assigning states to tips (based on the empirical distribution)
	double Pvalue = NAN_D;
	double mean_random_tauD = NAN_D;
	if(Npermutations>0){
		long count_deeper 				= 0;
		double sum_random_tauD 			= 0;
		long count_valid_permutations 	= 0;
		const double fraction_zeros 	= vector_count_zeros(current_state_per_tip)/double(Ntips);
		double random_tauD, random_varD, random_minD, random_maxD;
		long random_Npositives;
		std::vector<long> positive_clades, dummy_positive_clades, dummy_positives_per_clade, dummy_tips_per_clade;
		std::vector<double> dummy_mean_depth_per_clade;
		for(long p=0; p<Npermutations; ++p){
			for(long tip=0; tip<Ntips; ++tip) current_state_per_tip[tip] = (random_bernoulli(fraction_zeros) ? 0 : 1);
			aux_get_trait_depth_consenTRAIT(Ntips,
											Nnodes,
											Nedges,
											root,
											tree_edge,
											edge_length,
											current_state_per_tip,
											threshold_fraction,
											count_singletons,
											weighted,
											traversal_queue,
											traversal_node2first_edge,
											traversal_node2last_edge,
											traversal_edges,
											clade2parent,
											incoming_edge_per_clade,
											node2max_tip_distance,
											singleton_threshold,
											false,
											dummy_tips_per_clade,
											dummy_positives_per_clade,
											dummy_mean_depth_per_clade,
											random_tauD,
											random_varD,
											random_minD,
											random_maxD,
											random_Npositives,
											dummy_positive_clades,
											false,
											verbose_prefix);
			if(!isnan(random_tauD)){
				++count_valid_permutations;
				if(random_tauD>=tauD) ++count_deeper;
				sum_random_tauD += random_tauD;
			}
		}
		Pvalue = (count_valid_permutations>0 ? count_deeper/count_valid_permutations : NAN_D);
		mean_random_tauD = (count_valid_permutations>0 ? sum_random_tauD/count_valid_permutations : NAN_D);
	}

	return Rcpp::List::create(	Rcpp::Named("tauD") 					= tauD,
								Rcpp::Named("varD") 					= varD,
								Rcpp::Named("minD") 					= minD,
								Rcpp::Named("maxD") 					= maxD,
								Rcpp::Named("Npositives") 				= Npositives,
								Rcpp::Named("tips_per_clade")			= Rcpp::wrap(tips_per_clade),
								Rcpp::Named("positive_clades")			= Rcpp::wrap(positive_clades),
								Rcpp::Named("positives_per_clade")		= Rcpp::wrap(positives_per_clade),
								Rcpp::Named("mean_depth_per_clade")		= Rcpp::wrap(mean_depth_per_clade),
								Rcpp::Named("Pvalue") 					= Pvalue,
								Rcpp::Named("mean_random_tauD")			= mean_random_tauD);
}





// Calculate phylogenetic ("spatial") autocorrelation (AC) function of a continuous trait on tips of a tree
// AC(x) = correlation between the trait values of two random tips at distance x from each other
// Distance x = cumulative branch length between the two tips (aka "patristic distance")
// The autocorrelation function is estimated on a discrete distance grid (with x ranging from 0 to max)
//
// Input: 
//   A phylogenetic tree, and the states of a continuous trait on all tips of the tree.
// Output: 
//   The distance grid on which AC was estimated.
//	 Estimated AC(x) for each x on the distance grid
//   The number of tip pairs used to estimate AC(x) for each x on the distance grid
//	The mean absolute deviation between trait values of tips at distance x, for each x on the distance grid
//
// The time complexity of this routine is O(Ntips + Npairs*Nanc), where Nanc is the typical number of ancestors per tip.
// The memory complexity is O(Nedges+Npairs)
//
// Requirements:
//   Tree must be rooted (root will be determined automatically based on the edges)
//   Tree may include monofurcations or multifurcations.

// [[Rcpp::export]]
Rcpp::List autocorrelation_function_of_continuous_trait_CPP(const long 			Ntips,
															const long 			Nnodes,
															const long 			Nedges,
															const IntegerVector &tree_edge,			// 2D array of size Nedges x 2, flattened in row-major format
															const NumericVector &edge_length, 		// 1D array of size Nedges, or an empty std::vector (all branches have length 1)
															const NumericVector	&state_per_tip,		// 1D std::vector of numeric states of size Ntips.
															const long			Npairs,				// number of random tip pairs for estimating the correlation function
															const long			Nbins,				// number of different distances x at which to estimate the correlation function. This is the number of bins spanning all possible phylogenetic distances.
															bool 				verbose,
															const std::string 	&verbose_prefix){	
	// pick random tip pairs
	IntegerVector tipsA(Npairs);
	IntegerVector tipsB(Npairs);
	for(long p=0; p<Npairs; ++p){
		tipsA[p] = uniformIntWithin(0,Ntips-1);
		tipsB[p] = uniformIntWithin(0,Ntips-1);
	}

	// calculate distance for each tip pair
	const NumericVector distances = get_distances_between_clades_CPP(	Ntips,
																		Nnodes,
																		Nedges,
																		tree_edge,
																		edge_length,
																		tipsA,
																		tipsB,
																		verbose,
																		verbose_prefix);
																		
	// determine distance bins
	const double min_distance = get_array_min(distances);
	const double max_distance = get_array_max(distances);
	const double dx = (max_distance-min_distance)/Nbins;
	std::vector<double> distance_grid(Nbins);
	for(long b=0; b<Nbins; ++b){
		distance_grid[b] = min_distance + dx*(b+0.5);
	}
																		
	// calculate correlation within each distance-bin
	// use empirical means & variances for tipsA and tipsB separately, and separately for each bin, to ensure that empirical autocorrelations are always within [-1,1].
	std::vector<double> meanA_per_bin(Nbins,0);
	std::vector<double> meanB_per_bin(Nbins,0);
	std::vector<long>   pairs_per_bin(Nbins,0);
	for(long p=0, bin; p<Npairs; ++p){
		bin = max(0L,min(Nbins-1,long(floor((distances[p]-min_distance)/dx))));
		pairs_per_bin[bin] += 1;
		meanA_per_bin[bin] += state_per_tip[tipsA[p]];
		meanB_per_bin[bin] += state_per_tip[tipsB[p]];
	}
	for(long bin=0; bin<Nbins; ++bin){
		if(pairs_per_bin[bin]==0) continue;
		meanA_per_bin[bin] /= pairs_per_bin[bin];
		meanB_per_bin[bin] /= pairs_per_bin[bin];
	}
	double stateA, stateB;
	std::vector<double> varA_per_bin(Nbins,0);
	std::vector<double> varB_per_bin(Nbins,0);
	std::vector<double> autocorrelations(Nbins,0);
	std::vector<double> mean_abs_deviations(Nbins,0);
	std::vector<double> mean_rel_deviations(Nbins,0);
	for(long p=0, bin; p<Npairs; ++p){
		stateA 	= state_per_tip[tipsA[p]];
		stateB 	= state_per_tip[tipsB[p]];
		bin 	= max(0L,min(Nbins-1,long(floor((distances[p]-min_distance)/dx))));
		varA_per_bin[bin] += SQR(stateA-meanA_per_bin[bin]);
		varB_per_bin[bin] += SQR(stateB-meanB_per_bin[bin]);
		autocorrelations[bin] 	 += (stateA-meanA_per_bin[bin])*(stateB-meanB_per_bin[bin]);
		mean_abs_deviations[bin] += abs(stateA-stateB);
		mean_rel_deviations[bin] += (stateA==stateB ? 0.0 : abs(stateA-stateB)/(0.5*(abs(stateA)+abs(stateB))));
	}
	for(long bin=0; bin<Nbins; ++bin){
		if(pairs_per_bin[bin]==0){
			varA_per_bin[bin] = NAN_D;
			varB_per_bin[bin] = NAN_D;
			autocorrelations[bin] 	 = NAN_D;
			mean_abs_deviations[bin] = NAN_D;
			mean_rel_deviations[bin] = NAN_D;
		}else{
			varA_per_bin[bin] /= pairs_per_bin[bin];
			varB_per_bin[bin] /= pairs_per_bin[bin];
			autocorrelations[bin] 	 /= pairs_per_bin[bin] * sqrt(varA_per_bin[bin]*varB_per_bin[bin]);
			mean_abs_deviations[bin] /= pairs_per_bin[bin];
			mean_rel_deviations[bin] /= pairs_per_bin[bin];
		}
	}
	
	return Rcpp::List::create(	Rcpp::Named("distance_grid") 		= Rcpp::wrap(distance_grid),
								Rcpp::Named("N_per_grid_point")		= Rcpp::wrap(pairs_per_bin),
								Rcpp::Named("autocorrelations") 	= Rcpp::wrap(autocorrelations),
								Rcpp::Named("mean_rel_deviations") 	= Rcpp::wrap(mean_rel_deviations),
								Rcpp::Named("mean_abs_deviations") 	= Rcpp::wrap(mean_abs_deviations));
}




// For each node in a tree, calculate the empirical frequencies of states of a discrete trait, based on the states of descending tips
// This may be a very crude reconstruction of ancestral state probabilities (when normalized)
// Returns a 2D integer array of size Nnodes x Nstates, in row-major format
// [[Rcpp::export]]
IntegerVector get_empirical_state_frequencies_per_node(	const long			Ntips,
														const long			Nnodes,
														const long			Nedges,
														const long			Nstates,		// (INPUT) number of discrete states for the trait
														const IntegerVector &tree_edge,		// (INPUT) 2D array (in row-major format) of size Nedges x 2, or an empty std::vector (no tree available). A tree is needed if the tip_distribution relies on a tree structure.
														const IntegerVector	&tip_states){	// (INPUT) 1D array of size Ntips, listing the discrete state for each tip
	
	
	// determine parent clade for each clade
	std::vector<long> clade2parent;
	get_parent_per_clade(Ntips, Nnodes, Nedges, tree_edge, clade2parent);

	// find root using the mapping clade2parent
	const long root = get_root_from_clade2parent(Ntips, clade2parent);

	// get tree traversal route (tips --> root)
	// traversal_queue[] will be of size Nclades, and will have entries in 0:(Nclades-1)
	std::vector<long> traversal_queue, traversal_node2first_edge, traversal_node2last_edge, traversal_edges;
	get_tree_traversal_root_to_tips(	Ntips,
										Nnodes,
										Nedges,
										root,
										tree_edge,
										true,	// include tips
										false,	// edge mapping is not pre-computed
										traversal_queue,
										traversal_node2first_edge,
										traversal_node2last_edge,
										traversal_edges,
										false,
										"");
	
	// calculate empirical frequencies per node, traversing tips-->root (excluding the root)
	std::vector<long> frequencies_per_node(Nnodes*Nstates,0); // 2D array in row-major format
	long clade, parent;
	for(long q=traversal_queue.size()-1; q>=1; --q){
		clade	= traversal_queue[q];
		parent 	= clade2parent[clade];
		if(clade<Ntips){
			frequencies_per_node[(parent-Ntips)*Nstates+tip_states[clade]] += 1;
		}else{
			for(long s=0; s<Nstates; ++s) frequencies_per_node[(parent-Ntips)*Nstates + s] += frequencies_per_node[(clade-Ntips)*Nstates + s];
		}
	}
	
	return Rcpp::wrap(frequencies_per_node);
}





// [[Rcpp::export]]
Rcpp::List get_trait_richness_collectors_curve_CPP(	const long			Ntips,
													const long 			Nnodes,
													const long			Nedges,
													const long			Ntraits,
													const long 			root, 					// (INPUT) index of root node, i.e. an integer in 0:(Ntips+Nnodes-1)
													const IntegerVector &tree_edge, 			// (INPUT) 2D array (in row-major format) of size Nedges x 2, or an empty std::vector (no tree available). A tree is needed if the tip_distribution relies on a tree structure.
													const NumericVector &edge_length, 			// (INPUT) 1D array of size Nedges, or an empty std::vector (all branches have length 1). Only relevant if tip_distribution=="phylogenetic_diversity".
													const IntegerVector	&tip2first_trait,		// (INPUT) 1D array of size Nnodes, with values in 0:(NCtraits-1).
													const IntegerVector	&tip2last_trait,		// (INPUT) 1D array of size Nnodes, with values in 0:(NCtraits-1).
													const IntegerVector	&traits,				// (INPUT) 1D array of size NCtraits, with values in 0:(Ntraits-1)
													const IntegerVector	&rarefaction_depths,	// (INPUT) 1D array of size Ndepths, defining the different rarefaction depths at which to estimate the collector's curve.
													const long			Nrepeats,				// (INPUT) number of random repeats for calculating the average collector's curve at a given rarefaction depth
													const std::string	&tip_distribution,		// (INPUT) probability distribution for randomly choosing tips. Options are "uniform_tips" (each tip is equally likely), "uniform_children" (traversing root-->tips, with each child of a node being equally likely), “uniform_tips_without_replacement" (choose random tip without replacement), "uniform_children_without_replacement"
													const bool			use_realized_depths){	// (INPUT) if true, the rarefaction_depths are interpreted as centroids of realized rarefaction depth intervals, and the collector's curve is calculated as a function of realized rarefaction depth (rather than imposed rarefaction depth). Only relevant if tip_distribution is with replacement.
	
	const long Nclades = Ntips+Nnodes;
	const long Ndepths = rarefaction_depths.size();
	long tip, clade, count_tips_remaining;
	
	const bool distribution_uniform_tips 				= (tip_distribution=="uniform_tips");
	const bool distribution_uniform_tips_wr 			= (tip_distribution=="uniform_tips_without_replacement");
	const bool distribution_uniform_children 			= (tip_distribution=="uniform_children");
	const bool distribution_uniform_children_wr			= (tip_distribution=="uniform_children_without_replacement");
	const bool distribution_phylogenetic_diversity		= (tip_distribution=="phylogenetic_diversity");
	const bool need_tree 								= (distribution_uniform_children || distribution_uniform_children_wr || distribution_phylogenetic_diversity);

	// determine parent clade for each clade
	std::vector<long> clade2parent;
	if(need_tree) get_parent_per_clade(Ntips, Nnodes, Nedges, tree_edge, clade2parent);
		
	// prepare access structures for random tip selection
	std::vector<long> node2first_child, node2last_child, children;
	std::vector<long> node2total_tip_count;
	std::vector<double> clade2weight, node2diversity;
	if(distribution_uniform_children || distribution_uniform_children_wr || distribution_phylogenetic_diversity){
		get_children_per_node(	Ntips,
								Nnodes,
								Nedges,
								root,
								tree_edge,
								node2first_child,
								node2last_child,
								children);
	}
	if(distribution_uniform_children_wr){
		get_total_tip_count_per_node(	Ntips,
										Nnodes,
										Nedges,
										tree_edge,
										node2total_tip_count);
	}
	if(distribution_phylogenetic_diversity){
		// calculate phylogenetic diversities (cumulative branch lengths) for nodes
		get_cumulative_edge_lengths_per_node(	Ntips,
												Nnodes,
												Nedges,
												root,
												tree_edge,
												edge_length,
												node2diversity);
		// assign weights to children (for random traversal root-->tips) proportional to their phylogenetic diversity (PD)
		// set tip-children (and children with zero PD) to the average weight as their non-tip sister clades.
		clade2weight.resize(Nclades);
		for(long node=0; node<Nnodes; ++node){
			double cumulative_PD 	= 0;
			long tip_count 			= 0;
			const long child_count 	= node2last_child[node] - node2first_child[node] + 1;
			for(long c=node2first_child[node], child; c<=node2last_child[node]; ++c){
				child = children[c];
				if((child<Ntips) || (node2diversity[child-Ntips]==0)){
					tip_count += 1;
				}else{
					cumulative_PD += node2diversity[child-Ntips];
				}
			}
			const double PD_per_tip = ((tip_count<child_count) ? cumulative_PD/(child_count-tip_count) : 1.0);
			for(long c=node2first_child[node], child; c<=node2last_child[node]; ++c){
				child = children[c];
				clade2weight[child] = (((child<Ntips) || (node2diversity[child-Ntips]==0)) ? PD_per_tip : node2diversity[child-Ntips]);
			}			
		}
		clade2weight[root] = node2diversity[root-Ntips];
	}
	std::vector<long> tips_remaining_per_node;
								
	// calculate trait_richness collector's curve
	std::vector<char> 	tip_included;
	std::vector<char> 	trait_included(Ntraits);
	std::vector<double> 	trait_richness_means(Ndepths,0);
	std::vector<double> 	trait_richness_stds(Ndepths,0);
	std::vector<long> 	trait_richness_mins(Ndepths,Ntraits+1);
	std::vector<long> 	trait_richness_maxs(Ndepths,-1);
	std::vector<long>	Nrepeats_per_depth(Ndepths,0);
	for(long d=0; d<Ndepths; ++d){
		long rarefaction_depth = rarefaction_depths[d];
		if(distribution_uniform_tips_wr || distribution_uniform_children_wr) rarefaction_depth = min(Ntips,rarefaction_depth);
		for(long r=0; r<Nrepeats; ++r){
			// re-initialize some counters for this iteration
			trait_included.assign(Ntraits,false);
			tip_included.assign(Ntips,false);
			if(distribution_uniform_tips_wr){
				count_tips_remaining = Ntips;
			}else if(distribution_uniform_children_wr){
				tips_remaining_per_node = node2total_tip_count;
				clade2weight.resize(Nclades);
				for(clade=0; clade<Nclades; ++clade) clade2weight[clade] = (clade<Ntips ? 1.0 : (tips_remaining_per_node[clade-Ntips]>0 ? 1 : 0));
			}
			// choose multiple random tips and count trait richness
			for(long t=0; t<rarefaction_depth; ++t){
				if(distribution_uniform_tips){
					tip = uniformIntWithin(0,Ntips-1);

				}else if(distribution_uniform_tips_wr){
					// choose random tip from remaining tips (i.e. without replacement)
					long tip_counter = uniformIntWithin(0,count_tips_remaining-1);
					for(tip=0; tip<Ntips; ++tip){
						if(!tip_included[tip]) --tip_counter;
						if(tip_counter<0) break;
					}
					--count_tips_remaining;

				}else if(distribution_uniform_children){
					tip = get_tip_by_random_uniform_traversal(	Ntips,
																root,
																node2first_child,
																node2last_child,
																children);
				}else if(distribution_uniform_children_wr){
					tip = get_tip_by_random_traversal(	Ntips,
														root,
														node2first_child,
														node2last_child,
														children,
														clade2weight);
					clade2weight[tip] = 0; // prevent re-inclusion of this tip in the future (i.e. don't replace)
					// propagate information upwards
					clade = tip;
					while(clade!=root){
						clade = clade2parent[clade];
						tips_remaining_per_node[clade-Ntips] -= 1;
						if(tips_remaining_per_node[clade-Ntips]<=0){
							// no more tips to draw from this clade, so set weight to zero
							clade2weight[clade] = 0.0;
						}
					}
					
				}else if(distribution_phylogenetic_diversity){
					tip = get_tip_by_random_traversal(	Ntips,
														root,
														node2first_child,
														node2last_child,
														children,
														clade2weight);
				}
				tip_included[tip] = true;
				for(long i=tip2first_trait[tip]; i<=tip2last_trait[tip]; ++i){
					trait_included[traits[i]] = true;
				}
			}
			
			// count towards statistics of collector's curve
			const long trait_richness			= vector_sum(trait_included);
			const long effective_d				= (use_realized_depths ? get_nearest_index(rarefaction_depths, vector_sum(tip_included)) : d);
			trait_richness_means[effective_d] 	+= trait_richness;
			trait_richness_stds[effective_d]  	+= SQR(trait_richness);
			trait_richness_mins[effective_d] 	= min(trait_richness_mins[effective_d], trait_richness);
			trait_richness_maxs[effective_d] 	= max(trait_richness_maxs[effective_d], trait_richness);
			Nrepeats_per_depth[effective_d]  	+= 1;			
		}
	}

	for(long d=0; d<Ndepths; ++d){
		if(Nrepeats_per_depth[d]==0){
			trait_richness_means[d] = NAN_D;
			trait_richness_stds[d]  = NAN_D;
			trait_richness_mins[d]  = NAN_D;
			trait_richness_maxs[d]  = NAN_D;
		}else{
			trait_richness_means[d] /= Nrepeats_per_depth[d];
			trait_richness_stds[d]  = sqrt(trait_richness_stds[d]/Nrepeats_per_depth[d] - SQR(trait_richness_means[d]));
		}	
	}
	
	return Rcpp::List::create(	Rcpp::Named("trait_richness_means") = Rcpp::wrap(trait_richness_means),
								Rcpp::Named("trait_richness_stds") 	= Rcpp::wrap(trait_richness_stds),
								Rcpp::Named("trait_richness_mins") 	= Rcpp::wrap(trait_richness_mins),
								Rcpp::Named("trait_richness_maxs") 	= Rcpp::wrap(trait_richness_maxs),
								Rcpp::Named("Nrepeats_per_depth") 	= Rcpp::wrap(Nrepeats_per_depth));
}










// Calculate Phylogenetic Independent Contrasts (PIC) for multiple continuous traits on a tree [Felsenstein 1985, page 10]
// One PIC is returned for each non-monofurcating node (or each bifurcating node, if only_bifurcations=true).
// If only_bifurcations==false, then:
//    If multifurcations are present and, these are internally expanded to bifurcations and an additional PIC is returned for each such bifurcation.
//    Hence, the number of returned PICs is the number of bifurcations in the tree, after multifurcations have been expanded to bifurcations.
// Literature:
//    Felsenstein (1985). Phylogenies and the Comparative Method. The American Naturalist. 125:1-15.
// Requirements:
//   Tree can include multi- and mono-furcations.
//   Tree must be rooted. Root will be determined automatically as the node with no parent.
void get_phylogenetic_independent_contrasts(const long			Ntips,
											const long 			Nnodes,
											const long			Nedges,
											const long			Ntraits,
											const IntegerVector &tree_edge,			// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
											const NumericVector &edge_length,		// (INPUT) 1D array of size Nedges, or an empty std::vector (all edges have length 1)
											const NumericVector	&tip_states,		// (INPUT) 2D array of size Ntips x Ntraits, in row-major format, listing numeric states for each trait at each tip
											const bool			only_bifurcations,	// (INPUT) if true, PICs are only calculated for bifurcating nodes in the input tree, and multifurcations are not expanded.
											const bool			scaled,				// (INPUT) if true, then PICs are rescaled by the square-root of their expected variances (=the edge lengths connecting the compared nodes/tips)
											std::vector<double>	&PICs,				// (OUTPUT) 2D array of size Npics x Ntraits, in row-major format, listing PICs for each trait and for each considered node
											std::vector<double>	&distances,			// (OUTPUT) 1D array of size Npics, listing phylogenetic distances corresponding to the PICs. Under a Brownian motion mode, these are proportional to the variance of each PIC
											std::vector<long>	&PIC_nodes){		// (OUTPUT) 1D array of size Npics, listing node indices for which PICs were calculated. Negative values indicate nodes that were not actually in the original tree, but were created temporarily during expansion of multifurcations
	// check if tree has monofurcations & multifurcations
	long Nmonofurcations, Nbifurcations, Nmultifurcations;
	count_monofurcations_and_multifurcations(	Ntips,
												Nnodes,
												Nedges,
												tree_edge,
												Nmonofurcations,
												Nbifurcations,
												Nmultifurcations);
	
	std::vector<long> local_tree_edge;
	std::vector<double> local_edge_length;
	long Nlocal_edges, Nlocal_nodes;
	if((!only_bifurcations) && (Nmultifurcations>0)){
		// Tree has multifurcations, so expand them first to bifurcations
		// Note that the number of monofurcations will remain unchanged, but the number of bifurcations/nodes/edges will increase
		std::vector<long> dummy;
		multifurcations_to_bifurcations(Ntips,
										Nnodes,
										Nedges,
										tree_edge,
										edge_length,
										0,
										Nlocal_nodes,
										Nlocal_edges,
										local_tree_edge,
										local_edge_length,
										dummy);
	}else{
		local_tree_edge 	=  Rcpp::as<vector<long> >(tree_edge);
		local_edge_length 	=  Rcpp::as<vector<double> >(edge_length);
		Nlocal_nodes 		= Nnodes;
		Nlocal_edges 		= Nedges;
	}
	const long Nlocal_clades = Ntips + Nlocal_nodes;
	const long Npics = (only_bifurcations ? Nbifurcations : (Nlocal_nodes - Nmonofurcations));
	
	// get incoming edge for each clade
	std::vector<long> incoming_edge_per_clade;
	get_incoming_edge_per_clade(Ntips, Nlocal_nodes, Nlocal_edges, local_tree_edge, incoming_edge_per_clade);
	
	// get root
	const long root = get_root_from_incoming_edge_per_clade(Ntips, local_tree_edge, incoming_edge_per_clade);

	// prepare tree traversal route (root-->tips) and edge mappings
	std::vector<long> traversal_queue, node2first_edge, node2last_edge, edge_mapping;
	get_tree_traversal_root_to_tips(Ntips,
									Nlocal_nodes,
									Nlocal_edges,
									root,
									local_tree_edge,
									false,	// don't include tips
									false,	// edge mappings are not pre-calculated
									traversal_queue,
									node2first_edge,
									node2last_edge,	
									edge_mapping,
									false,
									"");
									
	// prepare incoming edge length per clade (will be modified later on as part of Felsenstein's algorithm)
	std::vector<double> incoming_length_per_clade(Nlocal_clades);
	if(local_edge_length.size()>0){
		for(long clade=0; clade<Nlocal_clades; ++clade){
			if(clade!=root) incoming_length_per_clade[clade] = local_edge_length[incoming_edge_per_clade[clade]];
		}
	}else{
		incoming_length_per_clade.assign(incoming_length_per_clade.size(),1);
		incoming_length_per_clade[root] = 0;
	}
	
									
	// calculate Felsenstein's PICs in a postorder traversal (tips-->root)
	const double edge_length_epsilon = RELATIVE_EPSILON * get_array_nonzero_min(local_edge_length); // substitute to be used for zero edge lengths
	std::vector<double> node_states(Nlocal_nodes*Ntraits,0);	// 2D numeric array of size Nlocal_nodes x Ntraits
	PICs.clear(); PICs.reserve(Npics*Ntraits);
	distances.clear(); distances.reserve(Npics);
	PIC_nodes.clear(); PIC_nodes.reserve(Npics);
	long edge1, edge2, child, child1, child2, node, clade, trait;
	double length, total_weight, X1, X2;
	for(long q=traversal_queue.size()-1; q>=0; --q){
		clade	= traversal_queue[q];
		node	= clade - Ntips;
		// calculate Felsenstein's X_k (node_states) and nu_k (incoming_length_per_clade)
		total_weight = 0;
		for(long e=node2first_edge[node]; e<=node2last_edge[node]; ++e){
			child 	= local_tree_edge[2*edge_mapping[e]+1];
			length 	= incoming_length_per_clade[child];
			if(length==0) length = edge_length_epsilon;
			for(trait=0; trait<Ntraits; ++trait){
				node_states[node*Ntraits+trait] += (1.0/length) * (child<Ntips ? tip_states[child*Ntraits+trait] : node_states[(child-Ntips)*Ntraits+trait]);
			}
			total_weight += (1.0/length);
		}
		for(trait=0; trait<Ntraits; ++trait) node_states[node*Ntraits+trait] /= total_weight;
		incoming_length_per_clade[clade] += 1.0/total_weight;
		
		// calculate PICs using Felsenstein's X_i & nu_i (skip over monofurcating nodes)
		if(1+node2last_edge[node]-node2first_edge[node] != 2) continue; // node is not bifurcating
		edge1		= edge_mapping[node2first_edge[node]];
		edge2		= edge_mapping[node2first_edge[node]+1];
		child1 		= local_tree_edge[2*edge1+1];
		child2 		= local_tree_edge[2*edge2+1];
		for(trait=0; trait<Ntraits; ++trait){
			X1 = (child1<Ntips ? tip_states[child1*Ntraits+trait] : node_states[(child1-Ntips)*Ntraits+trait]);
			X2 = (child2<Ntips ? tip_states[child2*Ntraits+trait] : node_states[(child2-Ntips)*Ntraits+trait]);
			PICs.push_back(X2 - X1);
		}
		distances.push_back(incoming_length_per_clade[child1] + incoming_length_per_clade[child2]);
		PIC_nodes.push_back(node<Nnodes ? node : -1); // keep track which node this PIC corresponds to. -1 means this temporary node did not exist in the original tree
	}
	
	// rescale if needed
	if(scaled){
		for(long p=0; p<Npics; ++p){
			for(trait=0; trait<Ntraits; ++trait){
				PICs[p*Ntraits+trait] /= sqrt(distances[p]);
			}
		}
	}
}




// Calculate Phylogenetic Independent Contrasts (PIC) for multiple continuous traits on a tree [Felsenstein 1985, page 10]
// This is an Rcpp wrapper for the function get_phylogenetic_independent_contrasts(..)
// [[Rcpp::export]]
Rcpp::List get_phylogenetic_independent_contrasts_CPP(	const long			Ntips,
														const long 			Nnodes,
														const long			Nedges,
														const long			Ntraits,
														const IntegerVector &tree_edge,			// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
														const NumericVector &edge_length,		// (INPUT) 1D array of size Nedges, or an empty std::vector (all edges have length 1)
														const NumericVector	&tip_states,		// (INPUT) 2D array of size Ntips x Ntraits, in row-major format, listing numeric states for each trait at each tip
														const bool			only_bifurcations,	// (INPUT) if true, PICs are only calculated for bifurcating nodes in the input tree, and multifurcations are not expanded.
														const bool			scaled){			// (INPUT)if true, then PICs are rescaled by the square-root of their expected variances (=the edge lengths connecting the compared nodes/tips)
		std::vector<double> PICs, distances;
		std::vector<long> PIC_nodes;
		get_phylogenetic_independent_contrasts(	Ntips,
												Nnodes,
												Nedges,
												Ntraits,
												tree_edge,
												edge_length,
												tip_states,
												only_bifurcations,
												scaled,
												PICs,
												distances,
												PIC_nodes);
		return Rcpp::List::create(	Rcpp::Named("PICs")  		= Rcpp::wrap(PICs),
									Rcpp::Named("distances") 	= Rcpp::wrap(distances),
									Rcpp::Named("nodes") 		= Rcpp::wrap(PIC_nodes));
}





// Fit a multivariate Brownian motion model for multiple correlated continuous traits evolving under Brownian motion
// Estimates the diffusivity matrix D[i,j], so that exp(-X^T*D^{-1}*X/(4*L))/sqrt(det(2*pi*D)) is the probability density for the multidimensional trait vector X after phylogenetic distance L, if initially located at the origin.
void fit_Brownian_motion_model(	const long			Ntips,
								const long 			Nnodes,
								const long			Nedges,
								const long			Ntraits,
								const IntegerVector &tree_edge,			// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
								const NumericVector &edge_length,		// (INPUT) 1D array of size Nedges, or an empty std::vector (all edges have length 1)
								const NumericVector	&tip_states,		// (INPUT) 2D array of size Ntips x Ntraits, in row-major format, listing numeric states for each trait at each tip
								std::vector<double>	&diffusivity){		// (OUTPUT) 2D array of size Ntraits x Ntraits, in row-major format, listing the fitted diffusion matrix D
											
	// calculate phylogenetic independent contrasts (PIC)
	// PICs correspond to independent increments of a multivariate Brownian motion process
	std::vector<double> PICs, distances;
	std::vector<long> PIC_nodes;
	get_phylogenetic_independent_contrasts(	Ntips,
											Nnodes,
											Nedges,
											Ntraits,
											tree_edge,
											edge_length,
											tip_states,
											false,
											true,	// rescale PICs by phylogenetic distances
											PICs,
											distances,
											PIC_nodes);
											
	// estimate diffusivity matrix based on independent contrasts
	// maximum-likelihood estimator on the intrinsic geometry of positive-definite matrices
	const long Npics = distances.size();
	diffusivity.assign(Ntraits*Ntraits,0); // 2D matrix of size Ntraits x Ntraits, in row-major format
	for(long t1=0; t1<Ntraits; ++t1){
		for(long t2=t1; t2<Ntraits; ++t2){
			for(long p=0; p<Npics; ++p){
				diffusivity[t1*Ntraits+t2] += PICs[p*Ntraits+t1]*PICs[p*Ntraits+t2];
			}
			diffusivity[t1*Ntraits+t2] /= (2*Npics);
			diffusivity[t2*Ntraits+t1]  = diffusivity[t1*Ntraits+t2];
		}
	}
}




// Calculate mean & standard deviation of a numeric trait across all extant clades over time
// This function requires that the trait is known for all tips and nodes of the tree
// [[Rcpp::export]]
Rcpp::List get_trait_stats_at_times_CPP(const long			Ntips,
										const long 			Nnodes,
										const long			Nedges,
										IntegerVector 		tree_edge,			// (INPUT) 2D array (in row-major format) of size Nedges x 2
										const NumericVector	&edge_length,		// (INPUT) 1D array of size Nedges, or empty
										const NumericVector	&times,				// (INPUT) 1D array of size Ntimes
										const NumericVector	&states){			// (INPUT) 1D array of size Nclades, listing the trait's value on each tip & node
	const long Nclades = Ntips + Nnodes;

	// determine parent clade for each clade
	std::vector<long> clade2parent;
	get_parent_per_clade(Ntips, Nnodes, Nedges, tree_edge, clade2parent);

	// find root using the mapping clade2parent
	const long root = get_root_from_clade2parent(Ntips, clade2parent);	
	
	// calculate distances from root
	std::vector<double> distances_from_root(Nclades);
	get_distances_from_root(Ntips,
							Nnodes,
							Nedges,
							tree_edge,
							edge_length,
							distances_from_root);
	
	const long Ntimes = times.size();
	std::vector<double> state_means(Ntimes,0), state_stds(Ntimes,0);
	std::vector<long> clade_counts(Ntimes,0); // number of clades over which the trait statistics were calculated, at each time point
	for(long t=0; t<Ntimes; ++t){
		for(long clade=0; clade<Nclades; ++clade){
			if(clade==root) continue;
			if((distances_from_root[clade]>=times[t]) && (distances_from_root[clade2parent[clade]]<=times[t]) && (!std::isnan(states[clade]))){
				clade_counts[t] += 1;
				state_means[t]  += states[clade];
				state_stds[t]   += SQ(states[clade]);
			}
		}
		state_means[t] /= clade_counts[t];
		state_stds[t] = sqrt(state_stds[t]/clade_counts[t] - SQ(state_means[t]));
	}
	return Rcpp::List::create(	Rcpp::Named("clade_counts")	= Rcpp::wrap(clade_counts),
								Rcpp::Named("means")		= Rcpp::wrap(state_means),
								Rcpp::Named("stds")			= Rcpp::wrap(state_stds));
}



#pragma mark -
#pragma mark Most recent common ancestors
#pragma mark 



// for each clade in a list of clades (nodes or tips), find a minimal set of tips so that the clade is the most recent common ancestor of those tips.
// This routine has an overhead at the beginning, but becomes more efficient when the number of given MRCAs is high.
// Time complexity = O(Ntips+Nnodes) + O(mrcas)
// Requirements:
//   The input tree must be rooted
//   The input tree can be multifurcating and/or monofurcating (except for the given MRCAs, each of which must have at least 2 children)
// Returns:
//   mrca2first_tip[]: 1D array of size Nmrcas, with values being indices in mrca_tips[]
//   mrca2last_tip[]: 1D array of size Nmrcas, with values being indices in mrca_tips[]
//   mrca_tips[]: 1D array with values in 0:(Ntips-1), so that mrca_tips[mrca2first_tip[m]],..,mrca_tips[mrca2last_tip[m]] are tips whose MRCA is clade mrcas[m]

// [[Rcpp::export]]
Rcpp::List get_mrca_defining_tips_CPP(	const long 			Ntips,
										const long 			Nnodes,
										const long 			Nedges,
										const IntegerVector &tree_edge,		// (INPUT) 2D array of size Nedges x 2 in row-major format, with values in 0:(Nclades-1)
										const IntegerVector	&mrcas,			// (INPUT) 1D array of size Nmrcas, with values in 0,..,(Nclades-1). For each mrcas[i], a set of tips is returned such that mrcas[i] is the MRCA of those tips
										bool				verbose,
										const std::string	&verbose_prefix){
	/* indexing conventions:
		parent, child: always run within 0,..,Nclades-1
		node: always runs within 0,..,Nnodes-1
		tip: always runs within 0,..,Ntips
	*/
	const long Nmrcas = mrcas.size();
	long clade, node, child;

	// determine parent clade for each clade
	std::vector<long> clade2parent;
	get_parent_per_clade(Ntips, Nnodes, Nedges, tree_edge, clade2parent);

	// find root using the mapping clade2parent
	const long root = get_root_from_clade2parent(Ntips, clade2parent);
	
	// get tree traversal route (tips --> root)
	std::vector<long> traversal_queue, traversal_node2first_edge, traversal_node2last_edge, traversal_edges;
	get_tree_traversal_root_to_tips(	Ntips,
										Nnodes,
										Nedges,
										root,
										tree_edge,
										true,
										false,
										traversal_queue,
										traversal_node2first_edge,
										traversal_node2last_edge,
										traversal_edges,
										verbose,
										verbose_prefix);
	reverse_array(traversal_queue); // make tips-->roots
										
	// to each node, assign one tip descending from it (traverse tips --> root)
	std::vector<long> node2tip(Nnodes, -1);
	for(long q=0, parent; q<traversal_queue.size(); ++q){
		clade = traversal_queue[q];
		if(clade==root) continue;
		parent = clade2parent[clade];
		if(clade<Ntips){
			// clade is a tip
			node2tip[parent-Ntips] = clade;
		}else{
			// clade is a node, so propagate its tip upwards
			node2tip[parent-Ntips] = node2tip[clade-Ntips];
		}
	}
	
	// for each MRCA, collect one descending tip per child
	std::vector<long> mrca2first_tip(Nmrcas);
	std::vector<long> mrca2last_tip(Nmrcas);
	std::vector<long> mrca_tips;
	mrca_tips.reserve(2*Nmrcas); // will need 2 tips per mrca
	for(long m=0, mrca; m<Nmrcas; ++m){
		mrca = mrcas[m];
		if(mrca<Ntips){
			// mrca is a tip
			mrca2first_tip[m] = mrca_tips.size();
			mrca2last_tip[m]  = mrca_tips.size();
			mrca_tips.push_back(mrca);
		}else{
			// mrca is a node, so collect one tip per child (for up to 2 children, no more needed)
			node = mrca-Ntips;
			mrca2first_tip[m] = mrca_tips.size();
			for(long ei=traversal_node2first_edge[node]; ei<=min(traversal_node2last_edge[node],traversal_node2first_edge[node]+1); ++ei){
				child = tree_edge[traversal_edges[ei]*2+1];
				if(child<Ntips) mrca_tips.push_back(child);
				else mrca_tips.push_back(node2tip[child-Ntips]);
			}
			mrca2last_tip[m] = mrca_tips.size()-1;
		}
	}
	
	return Rcpp::List::create(	Rcpp::Named("mrca2first_tip") 	= Rcpp::wrap(mrca2first_tip),
								Rcpp::Named("mrca2last_tip") 	= Rcpp::wrap(mrca2last_tip),
								Rcpp::Named("mrca_tips")		= Rcpp::wrap(mrca_tips));
}



// determine pairwise ancestry relationships between a subset of focal clades
// i.e. determine all pairs c,a for which focal_clades[a] is an ancestor of focal_clades[c]
// Requirements:
//   The input tree must be rooted
//   The input tree can be multifurcating and/or monofurcating
// Returns a 2D array in row-major format, listing pairs (c,a)

// [[Rcpp::export]]
IntegerVector get_pairwise_ancestries_CPP(	const long 			Ntips,
											const long 			Nnodes,
											const long 			Nedges,
											const long			root,				// (INPUT) integer within Ntips:(Ntips+Nnodes-1)
											const IntegerVector &tree_edge,			// (INPUT) 2D array of size Nedges x 2 in row-major format, with values in 0:(Nclades-1)
											const IntegerVector	&focal_clades){		// (INPUT) 1D array of size Nfocals, with values in 0,..,(Nclades-1)

	const long Nclades = Ntips + Nnodes;
	long ancestor, clade;

	// determine parent clade for each clade
	std::vector<long> clade2parent(Ntips+Nnodes, -1);
	for(long edge=0; edge<Nedges; ++edge){
		clade2parent[tree_edge[edge*2+1]] = tree_edge[edge*2+0];
	}
	
	// create mapping all_clades --> focal_clades (use -1 if not mapped)
	std::vector<long> all2focal_clade_index(Nclades,-1);
	for(long fc=0; fc<focal_clades.size(); ++fc){
		all2focal_clade_index[focal_clades[fc]] = fc;
	}

	// traverse upwards from each of the focal clades and keep track of encountered focal ancestors
	std::vector<long> descendant_ancestor_pairs; // list of descendant-ancestor pairs found, in row-major format. I.e. descendant_ancestor_pairs[2*p+0] is the descendant of descendant_ancestor_pairs[2*p+1].
	for(long fc=0; fc<focal_clades.size(); ++fc){
		clade = focal_clades[fc];
		if(clade==root) continue;
		ancestor = clade;
		do{
			ancestor = clade2parent[ancestor];
			if(all2focal_clade_index[ancestor]>=0){
				descendant_ancestor_pairs.push_back(fc);
				descendant_ancestor_pairs.push_back(all2focal_clade_index[ancestor]);
			}
		}while(ancestor!=root);
	}
	
	return  Rcpp::wrap(descendant_ancestor_pairs);
}




// Calculate most-recent-common-ancestors (MRCA) for pairs of clades (cladesA[] vs cladesB[])
// there's some initial overhead involved with this function, but for large number of clade pairs this becomes more efficient
// Time complexity is O(Ntips+Nnodes).
// Returns an IntegerVector of size Npairs, with each entry being the clade index of the MRCA of the pair
// If one clade is descendant of the other clade, the latter will be returned as MRCA
// Requirements:
//   The input tree must be rooted (root will be determined automatically, as the node that has no incoming edge)
//   The input tree can be multifurcating and/or monofurcating
// Attention: 0-based indexing is used for input and output variables, so make sure to shift indices in R before and after calling this function
// [[Rcpp::export]]
IntegerVector get_most_recent_common_ancestors_CPP(	const long 			Ntips,
													const long 			Nnodes,
													const long 			Nedges,
													const IntegerVector &tree_edge,			// 2D array of size Nedges x 2 in row-major format
													const IntegerVector &cladesA,			// 1D array of size Npairs, containing values in 0:(Nclades-1)
													const IntegerVector	&cladesB,			// 1D array of size Npairs, containing values in 0:(Nclades-1)
													bool 				verbose,
													const std::string 	&verbose_prefix){
	const long Npairs = cladesA.size();
	const long Nclades = Ntips + Nnodes;
	long parent, clade;

	// determine parent clade for each clade
	std::vector<long> clade2parent;
	get_parent_per_clade(Ntips, Nnodes, Nedges, tree_edge, clade2parent);
	
	// determine incoming edge per clade
	std::vector<long> incoming_edge_per_clade(Nclades,-1);
	for(long edge=0; edge<Nedges; ++edge){
		incoming_edge_per_clade[tree_edge[edge*2+1]] = edge;
	}
	
	// find root using the mapping clade2parent
	const long root = get_root_from_clade2parent(Ntips, clade2parent);
	
	// get tree traversal route (root --> tips)											
	std::vector<long> traversal_queue, traversal_node2first_edge, traversal_node2last_edge, traversal_edges;
	get_tree_traversal_root_to_tips(	Ntips,
										Nnodes,
										Nedges,
										root,
										tree_edge,
										true,
										false,
										traversal_queue,
										traversal_node2first_edge,
										traversal_node2last_edge,
										traversal_edges,
										verbose,
										verbose_prefix);
															
	// calculate number of ancestors and distance from root for each clade
	// (traverse root --> tips)
	std::vector<long> ancestor_count_per_clade(Nclades);
	ancestor_count_per_clade[root] = 0;
	for(long q=0; q<traversal_queue.size(); ++q){
		clade = traversal_queue[q];
		if(clade==root) continue;
		parent = clade2parent[clade];
		ancestor_count_per_clade[clade] = 1 + ancestor_count_per_clade[parent];
	}
	const long total_ancestor_count = vector_sum(ancestor_count_per_clade);
	
	// calculate ancestry for each clade in a long list ancestors[]
	// (traverse root --> tips)
	std::vector<long> clade2first_ancestor(Nclades); // for each clade c, ancestors[clade2first_ancestor[c]..clade2last_ancestor[c]] will be the list of ancestor clades leading to the clade c
	std::vector<long> clade2last_ancestor(Nclades);
	clade2first_ancestor[0] = 0;
	clade2last_ancestor[0] = clade2first_ancestor[0] + ancestor_count_per_clade[0] - 1;
	for(long c=1; c<Nclades; ++c){
		clade2first_ancestor[c] = clade2last_ancestor[c-1] + 1;
		clade2last_ancestor[c]  = clade2first_ancestor[c] + ancestor_count_per_clade[c] - 1;
	}
	std::vector<long> ancestors(total_ancestor_count);
	for(long q=0; q<traversal_queue.size(); ++q){
		clade = traversal_queue[q];
		if(clade==root) continue;
		parent = clade2parent[clade];
		// step 1: copy the parent's ancestry to the child's ancestry
		for(long a=clade2first_ancestor[parent]; a<=clade2last_ancestor[parent]; ++a){
			ancestors[clade2first_ancestor[clade]+(a-clade2first_ancestor[parent])] = ancestors[a];
		}
		// step 2: append the parent to the clade's ancestry
		ancestors[clade2last_ancestor[clade]] = parent;
	}
	
	// calculate most-recent-common-ancestor for each clade pair
	IntegerVector mrca_per_pair(Npairs);
	long cladeA, cladeB;
	for(long p=0; p<Npairs; ++p){
		cladeA = cladesA[p];
		cladeB = cladesB[p];
		// check for trivial case
		if(cladeA==cladeB){
			mrca_per_pair[p] = cladeA;
			continue;
		}
		// follow ancestry of both clades in synchrony, until they diverge
		// note that the first ancestor of every clade will be the root
		long a,b, mrca=-1;
		for(a=clade2first_ancestor[cladeA], b=clade2first_ancestor[cladeB]; (a<=clade2last_ancestor[cladeA]) && (b<=clade2last_ancestor[cladeB]); ++a, ++b){
			if(ancestors[a]!=ancestors[b]) break;
			else mrca = ancestors[a];
		}
		// check special case where one clade is descendant of the other (this would lead to a "premature" stop of the above loop)
		if((a<=clade2last_ancestor[cladeA]) && (ancestors[a]==cladeB)){
			mrca = cladeB;
		}else if((b<=clade2last_ancestor[cladeB]) && (ancestors[b]==cladeA)){
			mrca = cladeA;
		}
		mrca_per_pair[p] = mrca;
	}
	return(mrca_per_pair);
}




// Given a set of clades (tips & nodes, "descendants"), find their most recent common ancestor
// [[Rcpp::export]]
long get_most_recent_common_ancestor_CPP(	const long 			Ntips,
											const long 			Nnodes,
											const long 			Nedges,
											const IntegerVector &tree_edge,			// 2D array of size Nedges x 2 in row-major format
											const IntegerVector &descendants){		// 1D array of size ND, containing values in 0:(Nclades-1)
	const long Nclades = Ntips + Nnodes;
	const long ND = descendants.size();
	if(ND==0) return 0;
	if(ND==1) return descendants[0];
	
	// determine parent clade for each clade
	std::vector<long> clade2parent;
	get_parent_per_clade(Ntips, Nnodes, Nedges, tree_edge, clade2parent);

	// find root using the mapping clade2parent
	const long root = get_root_from_clade2parent(Ntips, clade2parent);

	// traverse upwards from each descendant and count the number of visits to each clade in the tree
	long clade;
	std::vector<long> visits_per_clade(Nclades, 0);
	for(long d=0; d<ND; ++d){
		clade = descendants[d];
		++visits_per_clade[clade];
		while(clade2parent[clade]>=0){
			clade = clade2parent[clade];
			++visits_per_clade[clade];
		}
	}
	
	// traverse from one of the descendants back to the root once more, to find the first clade from which all descendants passed
	clade = descendants[0];
	while(clade2parent[clade]>=0){
		if(visits_per_clade[clade]==ND) return clade; // all descendants passed by this clade, so it's the MRCA
		clade = clade2parent[clade];
	}
	return root;
}



// Given a rooted tree, determine whether a set of tips is monophyletic
// Requirements:
//   The tree must be rooted
//   The tree can include monofurcations and multifurcations
// [[Rcpp::export]]
bool is_monophyletic_tip_set_CPP(	const long 			Ntips,
									const long 			Nnodes,
									const long 			Nedges,
									const IntegerVector &tree_edge,			// (INPUT) 2D array of size Nedges x 2 in row-major format
									const IntegerVector &focal_tips){		// (INPUT) 1D array of size Nfocals, listing focal tip indices to test for monophyly
	if(focal_tips.size()<=1) return true;
	// find MRCA of focal tips
	const long mrca = get_most_recent_common_ancestor_CPP(	Ntips,
															Nnodes,
															Nedges,
															tree_edge,
															focal_tips);
	if(mrca<Ntips) return true; // MRCA is a tip, so must be monophyletic
	
	// get edge mappings for tree traversal
	std::vector<long> node2first_edge, node2last_edge, edges;
	get_node2edge_mappings(	Ntips,
							Nnodes,
							Nedges,
							tree_edge,
							node2first_edge,
							node2last_edge,
							edges);

	// traverse from MRCA to tips, to check if any of the reached tips are not in the focal set
	// use a scratch_stack for depth-first-search traversal
	std::vector<bool> tip_is_focal(Ntips, false);
	for(long t=0; t<focal_tips.size(); ++t) tip_is_focal[focal_tips[t]] = true;
	long child,node;
	std::vector<long> scratch_stack;
	scratch_stack.reserve(floor(2*log(Ntips)/log(2.0))); // rough estimate of typical tree depth x 2. scratch_stack may be resized along the way if needed.
	scratch_stack.push_back(mrca);
	while(scratch_stack.size()>0){
		node = scratch_stack.back() - Ntips;
		scratch_stack.pop_back();
		for(long e=node2first_edge[node]; e<=node2last_edge[node]; ++e){
			child = tree_edge[edges[e]*2+1];
			if((child<Ntips) && (!tip_is_focal[child])) return false; // reached a non-focal child
			if(child>=Ntips) scratch_stack.push_back(child); // add child node to stack for further exploration in the next iteration
		}
	}

	return true;
}





#pragma mark -
#pragma mark Ancestral state reconstruction
#pragma mark





// calculate the best (lowest) cost of any transition parent-->child, given a particular parent state and a particular child cost table (extracted from master_cost_table)
// this function assumes that the cost table for the child has already been calculated (hence, you should move tips-->root)
double aux_get_cost_of_parent_state_transitioning_to_one_child(	const long					Nstates,
																const long					parent_state, 
																const long 					edge,
																const double				edge_weight,
																const long 					child,
																const NumericVector			transition_costs, 		// (INPUT) 2D array of size Nstates x Nstates (in row-major format)
																const std::vector<double>	&master_cost_table,		// (INPUT) 2D array of size (Ntips+Nnodes) x Nstates (in row-major format)
																std::vector<double>			&scratch_space,			// temporary space for intermediate operations
																std::vector<long>			&master_transitions,				// (INPUT/OUTPUT) 1D array (preferably reserved up to size Nnodes*Nstates*Nstates)
																std::vector<long>			&edge_and_state2first_transition,	// (INPUT/OUTPUT) 1D array of size Nedges*Nstates.
																std::vector<long>			&edge_and_state2last_transition){	// (INPUT/OUTPUT) 1D array of size Nedges*Nstates.
	std::vector<double> &choice_costs = scratch_space;
	choice_costs.resize(Nstates);
	for(long state=0; state<Nstates; ++state){
		choice_costs[state] = transition_costs[parent_state*Nstates + state]*edge_weight + master_cost_table[child*Nstates + state];
	}
	const double best_cost = get_array_min(choice_costs);
	edge_and_state2first_transition[edge*Nstates + parent_state] = master_transitions.size();
	for(long transition=0; transition<Nstates; ++transition){
		if(abs(choice_costs[transition]-best_cost)<=RELATIVE_EPSILON*best_cost){
			master_transitions.push_back(transition);
		}
	}
	edge_and_state2last_transition[edge*Nstates + parent_state]	= master_transitions.size()-1;	
	return best_cost;
}




// calculate the cost of a particular state in a particular node, best on the best (lowest) costs of transitions to each of the children
// this function assumes that the cost table of each child has already been calculated (hence, you should move tips-->root)
double aux_get_cost_of_parent_state_transitioning_to_all_children(	const long					Nstates,
																	const long 					node, 							// (INPUT) integer in 0:(Nnodes-1)
																	const long 					parent_state,					// (INPUT) integer in 0:(Nstates-1)
																	const double				branch_length_exponent,			// (INPUT) non-negative number
																	const NumericVector			transition_costs,	 			// (INPUT) 2D array of size Nstates x Nstates (in row-major format)
																	const std::vector<double>	&master_cost_table,				// (INPUT) 2D array of size (Ntips+Nnodes) x Nstates (in row-major format)
																	const IntegerVector			&tree_edge,						// (INPUT) 2D array of size Nedges x 2 (in row-major format), in similar format as tree$edge in R "phylo" trees.
																	const NumericVector			&edge_length,				// (INPUT) 1D array of size Nedges
																	const std::vector<long>		&traversal_node2first_edge,		// (INPUT) 1D array of size Nnodes, with values in 0:(Nedges-1). Generated using the function get_tree_traversal_root_to_tips().
																	const std::vector<long>		&traversal_node2last_edge,		// (INPUT) 1D array of size Nnodes, with values in 0:(Nedges-1). Generated using the function get_tree_traversal_root_to_tips().
																	const std::vector<long>		&traversal_edges,				// (INPUT) 1D array of size Nedges, with values in 0:(Nedges-1). Generated using the function get_tree_traversal_root_to_tips().
																	std::vector<double>			&scratch_space,						// temporary space for intermediate operations
																	std::vector<long>			&master_transitions,				// (INPUT/OUTPUT) 1D array (preferably reserved up to size Nnodes*Nstates*Nstates)
																	std::vector<long>			&edge_and_state2first_transition,	// (INPUT/OUTPUT) 1D array of size Nedges*Nstates.
																	std::vector<long>			&edge_and_state2last_transition){	// (INPUT/OUTPUT) 1D array of size Nedges*Nstates.
	double S = 0;
	double edge_weight;
	long edge, child;
	for(long ei=traversal_node2first_edge[node]; ei<=traversal_node2last_edge[node]; ++ei){
		edge  = traversal_edges[ei];
		child = tree_edge[edge*2+1];
		edge_weight = (branch_length_exponent==0 ? 1.0 : 1/pow((edge_length.size()==0 ? 1.0 : edge_length[edge]),branch_length_exponent));
		S += aux_get_cost_of_parent_state_transitioning_to_one_child(	Nstates,
																		parent_state, 
																		edge, 
																		edge_weight, 
																		child,
																		transition_costs,
																		master_cost_table,
																		scratch_space,
																		master_transitions,
																		edge_and_state2first_transition,
																		edge_and_state2last_transition);
	}
	return S;
}
	
	

// Weighted maximum parsimony ansestral state reconstruction for discrete traits.
// Modification of Sankoff algorithm for reconstructing discrete ancestral states (Weighted Small Parsimony Problem)
// Sankoff's algorithm allows the inclusion of a cost matrix: transition_costs[i,j] is the cost of transitioning i-->j (ignoring edge length)
// The modification of this function is that optionally, edge lengths can be used to weight the transition costs:
// 	Longer edges imply smaller transition costs between states
// 	Specifically, the cost of transitioning is transition_cost[i,j]/(edge_length^branch_length_exponent)
//	where branch_length_exponent can be e.g. 0 (Sankoff's original algorithm), 1 (linear weighting) or 0.5 (square-root weighting, corresponding to a Brownian motion)
// Requirements:
//	Tree can be multifurcating, and can also include nodes with a single child
//	If (branch_length_exponent!=0) then: All branches must have length > 0
// For a description of the original Sankoff algorithm, see: 
//	http://telliott99.blogspot.ca/2010/03/fitch-and-sankoff-algorithms-for.html
//	(page 11) https://cs.brown.edu/courses/csci1950-z/slides/CSCI1950ZFall09_Lecture2.pdf
// The function returns a (non-flattened) NumericMatrix of size Nnodes x Nstates.
//
// Attention: Be carefull to use the C++ style indexing (0-based) when passing index-variables or index arrays to this function.
// For example, root must be a 0-based index, and tree_edge[] must have values in 0:(Ntips+Nnodes-1) instead of 1:(Ntips+Nnodes)

// [[Rcpp::export]]
NumericMatrix WMPR_ASR_CPP(	const long			Ntips,
							const long 			Nnodes,
							const long			Nedges,
							const long			Nstates, 				// (INPUT) number of possible states
							const IntegerVector	&tree_edge, 			// (INPUT) 2D array of size Nedges x 2 (in row-major format), in similar format as tree$edge in R "phylo" trees. This array holds the topology of the tree (apart from branch lengths).
							const NumericVector	&edge_length,			// (INPUT) 1D array of size Nedges, synchronized with the rows of tree_edge[,], i.e. with edge_length[e] being the length of edge e. Can also be an empty vector (all edges have length 1.0).
							const IntegerVector	&tip_states, 			// (INPUT) 1D array of size Ntips, with values being in 0:(Nstates-1)
							const NumericVector	&transition_costs,	 	// (INPUT) 2D array of size Nstates x Nstates (in row-major format), with transition_costs[i,j] being the cost of transition i-->j. Normally transition_cost[i,i]=0 for all i. Some transitions may be vorbitten, in which case the transition cost should be set to infinity (INFTY_D).
							const double 		branch_length_exponent, // (INPUT) exponent for weighting transition costs by branch length. To ignore branch lengths (i.e. to obtain the non-weighted MPR algorithm), set this to 0.
							bool				weight_posteriors_by_scenario_counts,	// (INPUT) if true, then the posterior_probability of a state (in a specific node) is proportional to the number of scenarios in which the trait is at that state
							bool				verbose,
							const std::string	&verbose_prefix){
	// Terminology in this function:
	// 	'node' runs from 0:(Nnodes-1)
	// 	'tip' runs from 0:(Ntips-1)
	// 	'parent' and 'child' runs from 0:(Ntips+Nnodes-1)
	// 	'edge' runs from 0:(Nedges-1)
	// 	'state' runs from 0:(Nstates-1)
	const long Nclades = Ntips+Nnodes;
	long node, state, parent, transition, child, edge;
	std::vector<double> scratch_space;

	// determine root
	const long root = get_root_clade(Ntips, Nnodes, Nedges, tree_edge);
	
	// create tree-access structures and determine order in which to traverse tree
	std::vector<long> traversal_queue_root2tips, traversal_node2first_edge, traversal_node2last_edge, traversal_edges;
	get_tree_traversal_root_to_tips(	Ntips,
										Nnodes,
										Nedges,
										root,
										tree_edge,
										false,
										false,
										traversal_queue_root2tips,
										traversal_node2first_edge,
										traversal_node2last_edge,
										traversal_edges,
										verbose,
										verbose_prefix);
	
	// get traversal route tips --> root
	std::vector<long> traversal_queue_tips2root = traversal_queue_root2tips;
	reverse_array(traversal_queue_tips2root); 
	
	
	// master_cost_table[,] should be a 2D numeric array of size (Ntips+Nnodes) x Nstates (in row-major-format)
	// the row master_cost_table[r,] is the cost table of tip or node r
	std::vector<double> master_cost_table(Nclades*Nstates, 0.0);
	
	// fill costs for tips (this is easy, since states are known)
	for(long tip=0; tip<Ntips; ++tip){
		for(state=0; state<Nstates; ++state){
			master_cost_table[tip*Nstates + state] = INFTY_D;
		}
		master_cost_table[tip*Nstates + tip_states[tip]] = 0.0;
	}

	// edge_and_state2first_transition[,] and edge_and_state2last_transition[,] contain indices mapping to master_transitions[]
	// for any edge e connecting parent node n to some child, and assuming n is at state s, the integers
	//	master_transitions[edge_and_state2first_transition[e,s]:edge_and_state2last_transition[e,s]]
	// are within in 1:Nstates, and are the optimal states to which node n switched during edge e.
	// Typically there will be only one "best transition" (i.e. edge_and_state2first_transition[e,s]==edge_and_state2last_transition[e,s]), 
	// but in case of multiple MPR solutions some edges may have multiple alternative best transitions (i.e. edge_and_state2first_transition[e,s] > edge_and_state2last_transition[e,s]

	// pre-allocate space at upper bound of possible need
	std::vector<long> master_transitions;
	master_transitions.reserve(Nnodes*Nstates*Nstates);
	std::vector<long> edge_and_state2first_transition(Nedges*Nstates);
	std::vector<long> edge_and_state2last_transition(Nedges*Nstates);
	
	
	// traverse tree (tips-->root) and build master_cost_table
	for(long parent_i=0; parent_i<Nnodes; ++parent_i){
		parent = traversal_queue_tips2root[parent_i];
		// calculate the cost associated with any state in this particular node
		for(state=0; state<Nstates; ++state){
			master_cost_table[parent*Nstates + state] = aux_get_cost_of_parent_state_transitioning_to_all_children(	Nstates,
																													(parent - Ntips),
																													state,
																													branch_length_exponent,
																													transition_costs,
																													master_cost_table,
																													tree_edge,
																													edge_length,
																													traversal_node2first_edge,
																													traversal_node2last_edge,
																													traversal_edges,
																													scratch_space,
																													master_transitions,
																													edge_and_state2first_transition,
																													edge_and_state2last_transition);
		}
	}
	
	// count number of scenarios (MPR solutions) implying each state in each node (based on lowest cost in the root, and then the associated transitions to the children)
	// scenario_count_per_node_and_state[n,s] will be the number of MPR solutions ("scenarios") in which node n is at state s
	// scenario_count_per_node_and_state[,] will be filled in the order root-->tips
	// See pages 18-19 in: https://cs.brown.edu/courses/csci1950-z/slides/CSCI1950ZFall09_Lecture2.pdf
	// Note: This should be floating point, not int, because in the latter case you risk integer overflow and thus the spontaneous introduction of negative values! This cost me 1 day of bug-hunting!
	std::vector<double> scenario_count_per_node_and_state(Nnodes*Nstates, 0.0); // scenario_count_per_node_and_state[Nnodes x Nstates] in row-major format
	
	const double best_root_cost = get_array_min(master_cost_table, root*Nstates, (root*Nstates+Nstates-1));
	for(state=0; state<Nstates; ++state){
		if(abs(master_cost_table[root*Nstates+state]-best_root_cost)<=RELATIVE_EPSILON*best_root_cost){
			scenario_count_per_node_and_state[(root-Ntips)*Nstates + state] = 1;
		}
	}
	

	for(long q=0; q<traversal_queue_root2tips.size(); ++q){
		parent 	= traversal_queue_root2tips[q];
		node	= parent-Ntips;
		for(long ei=traversal_node2first_edge[node]; ei<=traversal_node2last_edge[node]; ++ei){
			edge  = traversal_edges[ei];
			child = tree_edge[edge*2+1];
			if(child<Ntips) continue;
			for(state=0; state<Nstates; ++state){
				if(scenario_count_per_node_and_state[node*Nstates+state]>0){
					// examine all optimal transitions parent --> child, when parent is at this particular state
					for(long transition_i=edge_and_state2first_transition[edge*Nstates+state]; transition_i<=edge_and_state2last_transition[edge*Nstates+state]; ++transition_i){
						transition = master_transitions[transition_i];
						// increment scenario_count for the corresponding state in this particular child
						scenario_count_per_node_and_state[(child-Ntips)*Nstates + transition] += scenario_count_per_node_and_state[node*Nstates + state];
					}
				}
			}
		}
	}

		
	// For a given tree, there may be multiple alternative scenarios (MPR solutions) for the ancestral states
	// based on the scenario count per node and per state, define posterior_probabilities for nodes
	double mass;
	NumericMatrix posterior_probabilities(Nnodes,Nstates); // this will be a 2D array of size Nnodes x Nstates. Note that in this case we're not flattening, for convenience, because we're returning this to R and there we like to have a non-flattened 2D matrix. 
	for(node=0; node<Nnodes; ++node){
		mass = 0;
		if(weight_posteriors_by_scenario_counts){
			// weight states proportional to the number of scenarios
			for(state=0, mass=0; state<Nstates; ++state) mass += scenario_count_per_node_and_state[node*Nstates + state];
			if(mass==0){
				//if(verbose) Rcout << verbose_prefix << "WARNING: Node " << node << " (clade " << (node+Ntips) << ") has max-parsimony mass = 0 (i.e. no predicted state). It's posterior probabilities will all be set to NaN\n";
				for(state=0; state<Nstates; ++state) posterior_probabilities(node,state) = NAN_D;
			}else{
				for(state=0; state<Nstates; ++state) posterior_probabilities(node,state) = scenario_count_per_node_and_state[node*Nstates + state]/mass;
			}
		}else{
			// all states with non-zero scenario count are weighted equally
			for(state=0, mass=0; state<Nstates; ++state) mass += (scenario_count_per_node_and_state[node*Nstates + state]>0 ? 1.0 : 0.0);
			for(state=0; state<Nstates; ++state) posterior_probabilities(node,state) = (scenario_count_per_node_and_state[node*Nstates + state]>0 ? 1.0 : 0.0)/mass;			
		}
	}
	
	return(posterior_probabilities);
}



// calculate log-likelihood and posterior probability at the root, for a fixed-rates continuous-time Markov model for discrete character evolution
// A major computational bottleneck is the exponentiation of the transition matrix along each edge, i.e. exp(edge_length*transition_matrix)
// Exponentiation of the transition matrix can happen in one of 3 ways:
//    1. By providing all pre-calculated exponentials, via expQ_per_edge. This uses more RAM, but is much faster than the other methods below.
//    2. By providing an eigendecomposition of the transition_matrix, via eigenvalues, EVmatrix and inverse_EVmatrix. This is slower than 1 and similarly fast as 3. Uses lower RAM than the other options, but an accurate eigendecomposition may not alwasy be available.
//	  3. By providing polynomials of the transition matrix, via transition_polynomials, transition_polynomial_norms, exponentiation_balances & exponentiation_scaling_power. Slowest, but works well and can always be used as a last resort.
// The above options are checked and utilized in the listed order. Whenever the associated arrays of an option are empty, the next option is checked.
void aux_ASR_with_fixed_rates_Markov_model(	const long					Ntips,
											const long 					Nnodes,
											const long					Nedges,
											const long					Nstates,
											const long					root,								// (INPUT) root of the tree, an integer in Ntips,..,Nnodes+Ntips-1
											const IntegerVector			&tree_edge,							// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
											const NumericVector 		&edge_length, 						// (INPUT) 1D array of size Nedges, or an empty std::vector (all branches have length 1)
											const NumericVector			&prior_probabilities_per_tip, 		// (INPUT) 2D array of size Ntips x Nstates, in row-major format, listing prior probability distributions for each tip
											const NumericVector			&prior_probabilities_for_root,		// (INPUT) 1D array of size Nstates, listing prior probability distribution for root. Which prior you use for the root is often a matter of choice.
											const long					Npolynomials,						// (INPUT) Number of transition polynomials provided (polynomials of the rescaled transition matrix) for exponentiating the transition matrix
											const std::vector<double>	&transition_polynomials,			// (INPUT) Optional array of size Npolynomials x Nstates x Nstates, containing the pre-computed polynomials of the rescaled transition matrix (in layer-row-major format). Only relevant if use_precalculated_expQ==false. Can be empty if eigendecomposition is available instead.
											const std::vector<double>	&transition_polynomial_norms,		// (INPUT) Optional array of size Npolynomials, containing the Hilbert-Schmidt L2 norm of each transition_polynomial. Can be empty if eigendecomposition is available instead.
											const double				exponentiation_accuracy,			// (INPUT) desired accuracy when exponentiating using polynomials
											const long					min_polynomials,					// (INPUT) minimum number of polynomials to include during exponentiation, regardless of the pursued exponentiation_accuracy.
											const std::vector<double>	&exponentiation_balances,			// (INPUT) array of size Nstates, storing the diagonal elements of the balancing transformation that was applied to the polynomials (and should be reversed for the exponential)
											const long					exponentiation_scaling_power,		// (INPUT) base-2 scaling power that was applied prior to calculation of the polynomials (and should be reversed for the exponential)
											const std::vector<cdouble>	&eigenvalues,						// (INPUT) Optional 1D vector of size Nstates, listing the eigenvalues of the transition_matrix (corresponding to some diagonalization). Can also be an empty vector if eigendecomposition not available.
											const std::vector<cdouble>	&EVmatrix,							// (INPUT) Optional 2D array of size Nstates x Nstates, in row-major format, whose columns are the eigenvectors of the transition_matrix. Can also be an empty vector if eigendecomposition not available.
											const std::vector<cdouble>	&inverse_EVmatrix,					// (INPUT) Optional 2D array of size Nstates x Nstates, in row-major format, the inverse of EVmatrix. Can also be an empty vector if eigendecomposition not available.
											const std::vector<double>	&expQ_per_edge,						// (INPUT) 3D array of size Nedges x Nstates x Nstates, in layer-row-major format, listing the exponentiated transition matrix along each edge. Only relevant if use_precalculated_expQ==true.
											const std::vector<long>		&traversal_queue,					// (INPUT) 1D array of size Nnodes, with values in Ntips:(Nclades-1). Traversal queue root-->tips (not including tips). Generated using the function get_tree_traversal_root_to_tips(include_tips=true).
											const std::vector<long>		&traversal_node2first_edge,			// (INPUT) 1D array of size Nnodes, with values in 0:(Nedges-1). Generated using the function get_tree_traversal_root_to_tips().
											const std::vector<long>		&traversal_node2last_edge,			// (INPUT) 1D array of size Nnodes, with values in 0:(Nedges-1). Generated using the function get_tree_traversal_root_to_tips().
											const std::vector<long>		&traversal_edges,					// (INPUT) 1D array of size Nedges, with values in 0:(Nedges-1). Generated using the function get_tree_traversal_root_to_tips().
											std::vector<double>			&posteriors,						// (OUTPUT) 1D array of size Nnodes x Nstates, listing the posterior probabilities at each node. This is used both as internal scratch space as well as to return the final posteriors.
											double						&loglikelihood){					// (OUTPUT) log-likelihood for the full tree
	long clade, edge, child, node;
	const bool has_eigendecomposition = (eigenvalues.size()>0 && EVmatrix.size()>0 && inverse_EVmatrix.size()>0);
	const bool use_precalculated_expQ = (expQ_per_edge.size()>0);
	const double max_edge_length 	  = (edge_length.size()==0 ? 1.0 : get_array_max(edge_length));
							
							
	// calculate probability distribution on each node (traverse tips-->root)
	posteriors.assign(Nnodes*Nstates,1.0);
	std::vector<double> Y, expQ;
	std::vector<cdouble> exponentiation_scratch;
	double const *expQ_pointer;
	loglikelihood = 0;
	for(long q=traversal_queue.size()-1; q>=0; --q){
		clade  = traversal_queue[q];
		node   = clade-Ntips;
		// set probability distribution of clade to the element-wise product of its children's probability distributions
		for(long e=traversal_node2first_edge[node]; e<=traversal_node2last_edge[node]; ++e){
			edge  = traversal_edges[e];
			child = tree_edge[2*edge+1];
			if(use_precalculated_expQ){
				expQ_pointer = &expQ_per_edge[edge*Nstates*Nstates];
			}else if(has_eigendecomposition){
				get_matrix_exponential_using_eigendecomposition(Nstates,
																eigenvalues,
																EVmatrix,
																inverse_EVmatrix,
																(edge_length.size()==0 ? 1.0 : edge_length[edge]),
																exponentiation_scratch,
																expQ);
				expQ_pointer = &expQ[0];
			}else{
				// calculate exponential of transition matrix along edge
				get_matrix_exponential_using_balanced_polynomials(	Nstates,
																	Npolynomials,
																	transition_polynomials,
																	transition_polynomial_norms,
																	(edge_length.size()==0 ? 1.0 : edge_length[edge])/max_edge_length,
																	exponentiation_accuracy,
																	min_polynomials,
																	exponentiation_balances,
																	exponentiation_scaling_power,
																	expQ);
				expQ_pointer = &expQ[0];
			}
						
			// use exponentiated transition matrix to propagate probabilities from children to parent
			// probabilities[parent] = product_{child in children} exp(edge_length * Q^T) * probabilities[child]
			if(child<Ntips) multiply_matrix_with_vector(Nstates, Nstates, expQ_pointer, &prior_probabilities_per_tip[child*Nstates], Y);
			else multiply_matrix_with_vector(Nstates, Nstates, expQ_pointer, &posteriors[(child-Ntips)*Nstates], Y);
			for(long s=0; s<Nstates; ++s) posteriors[node*Nstates+s] *= max(0.0,Y[s]); // factor Y into the posterior of this node. Avoid negative values from rounding errors			
		}
	
		// multiply root's probability distribution with its prior
		if(clade==root){
			for(long s=0; s<Nstates; ++s) posteriors[node*Nstates+s] *= prior_probabilities_for_root[s];
		}
	
		// normalize clade's probability distribution
		double S = 0;
		for(long s=0; s<Nstates; ++s) S += posteriors[node*Nstates+s];
		for(long s=0; s<Nstates; ++s) posteriors[node*Nstates+s] /= S;
		loglikelihood += log(S);
	}
}







// re-root the tree and update the posterior probabilities at each node (since each node's posterior is determined by the posteriors of its children)
// since rerooting at a new node does not change the tip/node/edge indices, nor the total number of inout edges per clade, we just need to update a few tree-access data structures and update the posteriors for the affected nodes.
// Note: This algorithm cannot easily be generalized to rerooting at tips, because this would change the total number of tips & nodes and hence mess up the indexing of tips {0,..,Ntips-1} and nodes {Ntips,...,Ntips+Nnodes-1}.
// Exponentiation of the transition matrix can happen in one of 3 ways:
//    1. By providing all pre-calculated exponentials, via expQ_per_edge. This uses more RAM, but is much faster than the other methods below.
//    2. By providing an eigendecomposition of the transition_matrix, via eigenvalues, EVmatrix and inverse_EVmatrix. This is slower than 1 but faster than 3. Uses lower RAM than the other options, but an accurate eigendecomposition may not alwasy be available.
//	  3. By providing polynomials of the transition matrix, via transition_polynomials, transition_polynomial_norms, exponentiation_balances & exponentiation_scaling_power. Slowest, but works well and can always be used as a last resort.
// The above options are checked and utilized in the listed order. Whenever the associated arrays of an option are empty, the next option is checked.
void aux_reroot_and_update_ASR_with_fixed_rates_Markov_model(	const long					Ntips,
																const long 					Nnodes,
																const long					Nedges,
																const long					Nstates,
																const long					old_root,							// (INPUT) old (current) root of the tree, an integer in Ntips,..,Nnodes+Ntips-1
																const long					new_root,							// (INPUT) new root for the tree, an integer in Ntips,..,Nnodes+Ntips-1
																const NumericVector 		&edge_length, 						// (INPUT) 1D array of size Nedges, or an empty std::vector (all branches have length 1)
																const double				max_edge_length,					// (INPUT) 
																const std::vector<long>		&clade2first_inout_edge,			// (INPUT) 1D array of size Nnodes, with values in 0:(Nedges-1), mapping clades to their first incoming or outgoing edge.
																const std::vector<long>		&clade2last_inout_edge,				// (INPUT) 1D array of size Nnodes, with values in 0:(Nedges-1), mapping clades to their last incoming or outgoing edge.
																const std::vector<long>		&inout_edges,						// (INPUT) 1D array of size 2*Nedges, with values in 0:(Nedges-1). Maps internal edge indices (i.e. as listed in clade2first_inout_edge[] and clade2last_inout_edge[]) to original edge indices.
																const NumericVector			&prior_probabilities_per_tip, 		// (INPUT) 2D array of size Ntips x Nstates, in row-major format, listing prior probability distributions for each tip
																const NumericVector			&prior_probabilities_for_root,		// (INPUT) 1D array of size Nstates, listing prior probability distribution for root. Which prior you use for the root is often a matter of choice.
																const long					Npolynomials,						// (INPUT) Number of transition polynomials provided (polynomials of the rescaled transition matrix) for exponentiating the transition matrix
																const std::vector<double>	&transition_polynomials,			// (INPUT) Optional array of size Npolynomials x Nstates x Nstates, containing the pre-computed polynomials of the rescaled transition matrix (in layer-row-major format). Only relevant if use_precalculated_expQ==false. Can be empty if eigendecomposition is available instead.
																const std::vector<double>	&transition_polynomial_norms,		// (INPUT) Optional array of size Npolynomials, containing the Hilbert-Schmidt L2 norm of each transition_polynomial. Can be empty if eigendecomposition is available instead.
																const double				exponentiation_accuracy,			// (INPUT) Desired accuracy when exponentiating transition matrix using polynomials.
																const long					min_polynomials,					// (INPUT) minimum number of polynomials to include during exponentiation, regardless of the pursued exponentiation_accuracy.
																const std::vector<double>	&exponentiation_balances,			// (INPUT) array of size Nstates, storing the diagonal elements of the balancing transformation that was applied prior to calculation of the polynomials (and should be reversed for the exponential)
																const long					exponentiation_scaling_power,		// (INPUT) base-2 scaling power that was applied prior to calculation of the polynomials (and should be reversed for the exponential)
																const std::vector<cdouble>	&eigenvalues,						// (INPUT) Optional 1D vector of size Nstates, listing the eigenvalues of the transition_matrix (corresponding to some diagonalization). Can also be an empty vector if eigendecomposition not available.
																const std::vector<cdouble>	&EVmatrix,							// (INPUT) Optional 2D array of size Nstates x Nstates, in row-major format, whose columns are the eigenvectors of the transition_matrix. Can also be an empty vector if eigendecomposition not available.
																const std::vector<cdouble>	&inverse_EVmatrix,					// (INPUT) Optional 2D array of size Nstates x Nstates, in row-major format, the inverse of EVmatrix. Can also be an empty vector if eigendecomposition not available.
																const std::vector<double>	&expQ_per_edge,						// (INPUT) Optional 3D array of size Nedges x Nstates x Nstates, in layer-row-major format, listing the exponentiated transition matrix along each edge. Only relevant if use_precalculated_expQ==true. Can be empty, if exp(Q) is to be calculated using eigendecomposition or polynomials.
																std::vector<long> 			&tree_edge,							// (INPUT/OUTPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1). Will be updated after the rerooting.		
																std::vector<long>			&incoming_edge_per_clade,			// (INPUT/OUTPUT) 1D array of size Nclades, with elements in 0,..,Nedges-1. Will be updated after the rerooting.
																std::vector<double>			&posteriors){						// (INPUT/OUTPUT) 1D array of size Nnodes x Nstates, listing the posterior probabilities at each node. Should be pre-computed for the current tree, and will be updated after the rerooting.
	const bool has_eigendecomposition = (eigenvalues.size()>0 && EVmatrix.size()>0 && inverse_EVmatrix.size()>0);
	const bool use_precalculated_expQ = (expQ_per_edge.size()>0);
	if(new_root==old_root) return; // nothing to do
	
	// reroot (this changes edge directions, but keeps tip/node/edge indices the same)
	reroot_tree_at_node(Ntips, Nnodes, Nedges, old_root, new_root, tree_edge, incoming_edge_per_clade);
	
	// update posteriors of nodes that have been traversed by the rerooting
	std::vector<double> Y, expQ;
	std::vector<cdouble> exponentiation_scratch;
	double const *expQ_pointer; // will point to the location of the exponentiated transition matrix (which may be different for each edge)
	long clade = old_root;
	while(true){
		// re-calculate posterior for this clade (based on the posteriors of its children)
		const long node = clade-Ntips;
		for(long s=0; s<Nstates; ++s) posteriors[node*Nstates+s] = 1.0; // initialize posteriors for this node, repopulate below based on children
		for(long e=clade2first_inout_edge[clade], edge, child; e<=clade2last_inout_edge[clade]; ++e){
			edge  = inout_edges[e];
			if(tree_edge[2*edge+0]!=clade) continue; // this edge is not outgoing from this clade
			child = tree_edge[2*edge+1];
			if(use_precalculated_expQ){
				expQ_pointer = &expQ_per_edge[edge*Nstates*Nstates];
			}else if(has_eigendecomposition){
				get_matrix_exponential_using_eigendecomposition(Nstates,
																eigenvalues,
																EVmatrix,
																inverse_EVmatrix,
																(edge_length.size()==0 ? 1.0 : edge_length[edge]),
																exponentiation_scratch,
																expQ);
				expQ_pointer = &expQ[0];
			}else{
				// calculate exponential of transition matrix along edge
				get_matrix_exponential_using_balanced_polynomials(	Nstates,
																	Npolynomials,
																	transition_polynomials,
																	transition_polynomial_norms,
																	(edge_length.size()==0 ? 1.0 : edge_length[edge])/max_edge_length,
																	exponentiation_accuracy,
																	min_polynomials,
																	exponentiation_balances,
																	exponentiation_scaling_power,
																	expQ);
				expQ_pointer = &expQ[0];
			}
			// use exponentiated transition matrix to propagate probabilities from children to parent
			// probabilities[parent] = product_{child in children} exp(edge_length * Q^T) * probabilities[child]
			if(child<Ntips) multiply_matrix_with_vector(Nstates, Nstates, expQ_pointer, &prior_probabilities_per_tip[child*Nstates], Y);
			else multiply_matrix_with_vector(Nstates, Nstates, expQ_pointer, &posteriors[(child-Ntips)*Nstates], Y);
			for(long s=0; s<Nstates; ++s) posteriors[node*Nstates+s] *= max(0.0,Y[s]); // factor Y into the posterior of this node. Avoid negative values from rounding errors
		}
	
		// multiply root's probability distribution with its prior
		if(clade==new_root){
			for(long s=0; s<Nstates; ++s) posteriors[node*Nstates+s] *= prior_probabilities_for_root[s];
		}
			
		// normalize clade's probability distribution
		double S = 0;
		for(long s=0; s<Nstates; ++s) S += posteriors[node*Nstates+s];
		for(long s=0; s<Nstates; ++s) posteriors[node*Nstates+s] /= S;
		
		// move on to parent
		if(clade!=new_root) clade = tree_edge[incoming_edge_per_clade[clade]*2+0];
		else break;
	}
}




// calculate the loglikelihood of a fixed-rates Markov model for discrete character evolution on a phylogenetic tree, provided a fixed transition matrix
// Optionally, the marginal ancestral states (likelihoods) can be computed for all nodes, using the rerooting algorithm by [Yang et al. (1995). Genetics 141:1641-1650]
// Reconstructing marginal ancestral states substantially increases computation time, so don't request this if you only care about the loglikelihood (e.g. for fitting purposes)
// Optionally, the eigendecomposition of the transition matrix (eigenvalues & eigenvectors) can be provided to speed up exponentiations. If provided, this is used blindly for all exponentiations.
// If an eigendecomposition for the transition_matrix is not provided, then a Taylor series (polynomials) approximation is used instead. This includes preconditioning steps and seems to work quite well, and is similarly fast as using an eigendecomposition.
// Requirements:
//   Tree can include multi- and mono-furcations.
//   Tree must be rooted. Root will be determined automatically as the node with no parent.
// [[Rcpp::export]]
Rcpp::List ASR_with_fixed_rates_Markov_model_CPP(	const long			Ntips,
													const long 			Nnodes,
													const long			Nedges,
													const long			Nstates,
													const IntegerVector &tree_edge,						// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
													const NumericVector &edge_length, 					// (INPUT) 1D array of size Nedges, or an empty std::vector (all edges have length 1)
													const NumericVector &transition_matrix,				// (INPUT) 2D array of size Nstates x Nstates, in row-major format. Transition-rate matrix Q in row-major format, i.e. Q[r*Nstates + c] is (r,c)-th-element of Q and equal to the transition rate r-->c.
													const ComplexVector	&eigenvalues,					// (INPUT) Optional 1D vector of size Nstates, listing the eigenvalues of the transition_matrix (corresponding to some diagonalization). Can also be an empty vector if eigendecomposition not available.
													const ComplexVector	&EVmatrix,						// (INPUT) Optional 2D array of size Nstates x Nstates, in row-major format, whose columns are the eigenvectors of the transition_matrix. Can also be an empty vector if eigendecomposition not available.
													const ComplexVector	&inverse_EVmatrix,				// (INPUT) Optional 2D array of size Nstates x Nstates, in row-major format, the inverse of EVmatrix. Can also be an empty vector if eigendecomposition not available.
													const NumericVector	&prior_probabilities_per_tip, 	// (INPUT) 2D array of size Ntips x Nstates, in row-major format, listing prior probability distributions for each tip
													const NumericVector	&prior_probabilities_for_root,	// (INPUT) 1D array of size Nstates, listing prior probability distribution for root. Which prior you use for the root is often a matter of choice.
													bool				include_ancestral_likelihoods,	// (INPUT) if true, then the marginal ancestral states estimates (conditional scaled likelihoods as if each node was a root) of all nodes are also returned as an array of size Nnodes x Nstates. Only use if needed, since it's computationally expensive.
													const double		exponentiation_accuracy,		// (INPUT) maximum allowed error when exponentiating the transition matrix via polynomials, in terms of the Hilbert-Schmidt L2 norm. Only relevant if exponentiation is done using the polynomials.
													const long			max_polynomials,				// (INPUT) maximum possible number of polynomials to use for exponentiating the transition_matrix via polynomials, regardless of the pursued accuracy epsilon. Used as safety vault, but may break the guaranteed accuracy. A value ~100 is usually enough.
													const bool			store_exponentials){			// (INPUT) if True, then exponentials are pre-calculated and stored for the calculation of ancestral_likelihoods. This may save time because each exponential is only calculated once, but will use up more memory since all exponentials are stored. Only relevant if include_ancestral_likelihoods==TRUE, otherwise exponentials are never stored.
	const long Nclades 					= Ntips + Nnodes;
	const bool use_precalculated_expQ 	= (include_ancestral_likelihoods && store_exponentials);
	const bool has_eigendecomposition	= (eigenvalues.size()>0 && EVmatrix.size()>0 && inverse_EVmatrix.size()>0);
	const double max_edge_length 		= (edge_length.size()==0 ? 1.0 : get_array_max(edge_length));
	const long root 					= get_root_clade(Ntips, Nnodes, Nedges, tree_edge);


	// transform some of the R vectors to C++ vectors
	std::vector<cdouble> eigenvaluesCPP, EVmatrixCPP, inverse_EVmatrixCPP;
	cast_ComplexVector_to_CPP(eigenvalues, eigenvaluesCPP);
	cast_ComplexVector_to_CPP(EVmatrix, EVmatrixCPP);
	cast_ComplexVector_to_CPP(inverse_EVmatrix, inverse_EVmatrixCPP);

	// prepare tree traversal route (root-->tips)
	std::vector<long> traversal_queue, traversal_node2first_edge, traversal_node2last_edge, traversal_edges;
	get_tree_traversal_root_to_tips(Ntips,
									Nnodes,
									Nedges,
									root,
									tree_edge,
									false,	// don't include tips
									false,	// edge mappings are not pre-calculated
									traversal_queue,
									traversal_node2first_edge,
									traversal_node2last_edge,	
									traversal_edges,
									false,
									"");
									

	// prepare data structures for exponentiations of transition matrix
	std::vector<double> transition_polynomials, transition_polynomial_norms, exponentiation_balances;
	long Npolynomials;
	const long min_polynomials = min_polynomials_for_positive_exponential_of_irreducible_matrix(Nstates, transition_matrix);
	long exponentiation_scaling_power;
	if(!has_eigendecomposition){
		calculate_balanced_matrix_polynomials(	Nstates,
												std::vector<double>(transition_matrix.begin(), transition_matrix.end()),
												max_edge_length,
												exponentiation_accuracy,
												min_polynomials,
												max_polynomials,
												transition_polynomials,
												transition_polynomial_norms,
												Npolynomials,
												exponentiation_balances,
												exponentiation_scaling_power);
	}
											
	// pre-calculate exponentials of transition_matrix along edges if needed
	// This is faster because it avoids repeated calculations of the same exponential, but needs more RAM if Nedges is large
	std::vector<double> expQ_per_edge;
	if(use_precalculated_expQ){
		std::vector<cdouble> exponentiation_scratch;
		std::vector<double> scratch_expQ;
		expQ_per_edge.resize(Nedges*Nstates*Nstates);
		for(long edge=0; edge<Nedges; ++edge){
			if(has_eigendecomposition){
				get_matrix_exponential_using_eigendecomposition(Nstates,
																eigenvaluesCPP,
																EVmatrixCPP,
																inverse_EVmatrixCPP,
																(edge_length.size()==0 ? 1.0 : edge_length[edge]),
																exponentiation_scratch,
																scratch_expQ);
			}else{
				get_matrix_exponential_using_balanced_polynomials(	Nstates,
																	Npolynomials,
																	transition_polynomials,
																	transition_polynomial_norms,
																	(edge_length.size()==0 ? 1.0 : edge_length[edge])/max_edge_length,
																	exponentiation_accuracy,
																	min_polynomials,
																	exponentiation_balances,
																	exponentiation_scaling_power,
																	scratch_expQ);
			}
			for(long r=0; r<Nstates; ++r){
				for(long c=0; c<Nstates; ++c){
					expQ_per_edge[edge*Nstates*Nstates + r*Nstates + c] = scratch_expQ[r*Nstates + c];
				}
			}
		}
	}

								
	// calculate loglikelihood & posteriors for all nodes (this is not the same as the ancestral likelihoods)
	std::vector<double> posteriors;
	double loglikelihood;
	aux_ASR_with_fixed_rates_Markov_model(	Ntips,
											Nnodes,
											Nedges,
											Nstates,
											root,
											tree_edge,
											edge_length,
											prior_probabilities_per_tip,
											prior_probabilities_for_root,
											Npolynomials,
											transition_polynomials,
											transition_polynomial_norms,
											exponentiation_accuracy,
											min_polynomials,
											exponentiation_balances,
											exponentiation_scaling_power,
											eigenvaluesCPP,
											EVmatrixCPP,
											inverse_EVmatrixCPP,
											expQ_per_edge,
											traversal_queue,
											traversal_node2first_edge,
											traversal_node2last_edge,
											traversal_edges,
											posteriors,
											loglikelihood);


	// calculate marginal ancestral states (posterior probabilities) at each node, as if that node was the root [Yang et al. 1995]
	// note that the original edge mappings (e.g. traversal_node2first_edge[]) will no longer be consistent with current_tree_edge after rerooting
	// Notation: current_(..) refers to tree-access structures that are updated at each rerooting. They will be consistent with each other, but not necessarily with the original tree structure.
	std::vector<double> ancestral_likelihoods(Nnodes*Nstates);
	if(include_ancestral_likelihoods){
		// prepare some data structures, which will be updated everytime we reroot
		long current_root = root;
		std::vector<long> current_tree_edge(tree_edge.begin(), tree_edge.end());
		std::vector<long> clade2first_inout_edge, clade2last_inout_edge, inout_edges; // will be invariant to rerootings
		get_inout_edges_per_clade(Ntips, Nnodes, Nedges, tree_edge, clade2first_inout_edge, clade2last_inout_edge, inout_edges);
		std::vector<long> current_incoming_edge_per_clade(Nclades,-1); // will change with each rerooting
		for(long edge=0; edge<Nedges; ++edge) current_incoming_edge_per_clade[current_tree_edge[edge*2+1]] = edge;
		// calculate depth-first-search traversal route
		// this minimizes the distance between successive old_root-new_root pairs and thus the computation required to update the posteriors on each rerooting
		// similarly, depth_first_search_queue[] will still list successive nodes, but not necessarily in downstream order
		std::vector<long> depth_first_search_queue;
		get_tree_traversal_depth_first_search(	Ntips,
												Nnodes,
												Nedges,
												root,
												tree_edge,
												false,	// don't include tips
												true, 	// edge mappings are pre-calculated
												depth_first_search_queue,
												traversal_node2first_edge,
												traversal_node2last_edge,
												traversal_edges);
		
		// reroot at each node, updating the tree-access structures and updating the posteriors at each node
		// edge directions will change, but tip/node/edge indices remain the same
		for(long q=0, new_root; q<depth_first_search_queue.size(); ++q){
			new_root = depth_first_search_queue[q];
			aux_reroot_and_update_ASR_with_fixed_rates_Markov_model(Ntips,
																	Nnodes,
																	Nedges,
																	Nstates,
																	current_root,
																	new_root,
																	edge_length,
																	max_edge_length,
																	clade2first_inout_edge,
																	clade2last_inout_edge,
																	inout_edges,
																	prior_probabilities_per_tip,
																	prior_probabilities_for_root,
																	Npolynomials,
																	transition_polynomials,
																	transition_polynomial_norms,
																	exponentiation_accuracy,
																	min_polynomials,
																	exponentiation_balances,
																	exponentiation_scaling_power,
																	eigenvaluesCPP,
																	EVmatrixCPP,
																	inverse_EVmatrixCPP,
																	expQ_per_edge,
																	current_tree_edge,
																	current_incoming_edge_per_clade,
																	posteriors);
			// the posteriors of the root are equal to its ancestral_likelihoods, so extract those
			for(long s=0; s<Nstates; ++s) ancestral_likelihoods[(new_root-Ntips)*Nstates+s] = posteriors[(new_root-Ntips)*Nstates + s];
			current_root = new_root;
			// abort if the user has interrupted the calling R program
			Rcpp::checkUserInterrupt();
		}
	}
	
	if(include_ancestral_likelihoods){
		return Rcpp::List::create(	Rcpp::Named("loglikelihood") = loglikelihood,
									Rcpp::Named("ancestral_likelihoods") = Rcpp::wrap(ancestral_likelihoods));
	}else{
		return Rcpp::List::create(	Rcpp::Named("loglikelihood") = loglikelihood);
	}
}



/* OLD CODE. WORKS, BUT LESS EFFICIENTLY.

// calculate the loglikelihood of a fixed-rates Markov model for discrete character evolution on a phylogenetic tree, provided a fixed transition matrix
// Optionally, the marginal ancestral states (likelihoods) can be computed for all nodes, using the rerooting algorithm by [Yang et al. (1995). Genetics 141:1641-1650]
// Reconstructing marginal ancestral states substantially increases computation time, so don't request this if you only care about the loglikelihood (e.g. for fitting purposes)
// Requirements:
//   Tree must be rooted. Root will be determined automatically as the node with no parent.
// [[Rcpp::export]]
Rcpp::List ASR_with_fixed_rates_Markov_model_CPP(	const long			Ntips,
													const long 			Nnodes,
													const long			Nedges,
													const long			Nstates,
													const IntegerVector &tree_edge,						// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
													const NumericVector &edge_length, 					// (INPUT) 1D array of size Nedges, or an empty std::vector (all branches have length 1)
													const NumericVector &transition_matrix,				// (INPUT) 2D array of size Nstates x Nstates, in row-major format. Transition-rate matrix Q in row-major format, i.e. Q[r*Nstates + c] is (r,c)-th-element of Q
													const NumericVector	&prior_probabilities_per_tip, 	// (INPUT) 2D array of size Ntips x Nstates, in row-major format, listing prior probability distributions for each tip
													const NumericVector	&prior_probabilities_for_root,	// (INPUT) 1D array of size Nstates, listing prior probability distribution for root. Which prior you use for the root is often a matter of choice.
													bool				include_ancestral_likelihoods){	// (INPUT) if true, then the marginal ancestral states (posterior probabilities as if each node was a root) of all nodes are also returned as an array of size Nnodes x Nstates. Only use if needed, since it's computationally expensive.
	const long Nclades = Ntips+Nnodes;
	const bool use_precalculated_expQ = include_ancestral_likelihoods;
	std::vector<long> current_tree_edge(tree_edge.begin(), tree_edge.end());
	const double max_edge_length = (edge_length.size()==0 ? 1.0 : get_array_max(edge_length));
	long current_root = get_root_clade(Ntips, Nnodes, Nedges, current_tree_edge);
	
	// prepare data structures for exponentiations of transition matrix
	std::vector<double> transition_polynomials, transition_polynomial_norms;
	long Npolynomials;
	const long NPmax = 1000; // max number of matrix polynomials for calculating the exponential
	const double exponentiation_accuracy = 1e-3;
	calculate_matrix_polynomials(Nstates,
								std::vector<double>(transition_matrix.begin(), transition_matrix.end()),
								max_edge_length,
								exponentiation_accuracy,
								NPmin,
								NPmax,
								transition_polynomials,
								transition_polynomial_norms,
								Npolynomials);
								
	// pre-calculate exponentials of transition_matrix along edges if needed
	// This is faster but needs more RAM if Nedges is large
	std::vector<double> expQ_per_edge;
	if(use_precalculated_expQ){
		expQ_per_edge.resize(Nedges*Nstates*Nstates);
		std::vector<double> scratch_expQ;
		for(long edge=0; edge<Nedges; ++edge){
			get_matrix_exponential_using_polynomials(	Nstates,
													Npolynomials,
													transition_polynomials,
													transition_polynomial_norms,
													(edge_length.size()==0 ? 1.0 : edge_length[edge])/max_edge_length,
													exponentiation_accuracy,
													scratch_expQ);
			for(long r=0; r<Nstates; ++r){
				for(long c=0; c<Nstates; ++c){
					expQ_per_edge[edge*Nstates*Nstates + r*Nstates + c] = scratch_expQ[r*Nstates + c];
				}
			}
		}
	}
								
	// calculate loglikelihood
	std::vector<double> posteriors;
	std::vector<long> scratch_traversal_queue, scratch_traversal_node2first_edge, scratch_traversal_node2last_edge, scratch_traversal_edges;
	double loglikelihood;
	aux_ASR_with_fixed_rates_Markov_model(	Ntips,
											Nnodes,
											Nedges,
											Nstates,
											current_root,
											current_tree_edge,
											edge_length,
											prior_probabilities_per_tip,
											prior_probabilities_for_root,
											Npolynomials,
											transition_polynomials,
											transition_polynomial_norms,
											exponentiation_accuracy,
											NPmin,
											expQ_per_edge,
											use_precalculated_expQ,
											scratch_traversal_queue,
											scratch_traversal_node2first_edge,
											scratch_traversal_node2last_edge,
											scratch_traversal_edges,
											posteriors,
											loglikelihood);

	// calculate marginal ancestral states (posterior probabilities) at each node, as if that node was the root [Yang et al. 1995]
	std::vector<double> ancestral_likelihoods(Nnodes*Nstates);
	if(include_ancestral_likelihoods){
		std::vector<long> clade2first_inout_edge, clade2last_inout_edge, inout_edges;
		get_inout_edges_per_clade(Ntips, Nnodes, Nedges, current_tree_edge, clade2first_inout_edge, clade2last_inout_edge, inout_edges);
		double dummy_loglikelihood;
		for(long new_root=Ntips; new_root<Nclades; ++new_root){
			//cout << "  debug: rerooting at clade # " << new_root << endl; // debug
			// reroot tree at this node (edge directions will change, but tip/node/edge indices remain the same)
			reroot_tree_at_node(Ntips, Nnodes, Nedges, current_root, new_root, current_tree_edge);
			current_root = new_root;
			// calculate posterior as if this node was the root
			aux_ASR_with_fixed_rates_Markov_model(	Ntips,
													Nnodes,
													Nedges,
													Nstates,
													new_root,
													current_tree_edge,
													edge_length,
													prior_probabilities_per_tip,
													prior_probabilities_for_root,
													Npolynomials,
													transition_polynomials,
													transition_polynomial_norms,
													exponentiation_accuracy,
													NPmin,
													expQ_per_edge,
													use_precalculated_expQ,
													scratch_traversal_queue,
													scratch_traversal_node2first_edge,
													scratch_traversal_node2last_edge,
													scratch_traversal_edges,
													posteriors,
													dummy_loglikelihood);
			for(long s=0; s<Nstates; ++s) ancestral_likelihoods[(new_root-Ntips)*Nstates+s] = posteriors[(new_root-Ntips)*Nstates + s];
		}
	}
	
	if(include_ancestral_likelihoods){
		return Rcpp::List::create(	Rcpp::Named("loglikelihood") = loglikelihood,
									Rcpp::Named("ancestral_likelihoods") = Rcpp::wrap(ancestral_likelihoods));
	}else{
		return Rcpp::List::create(	Rcpp::Named("loglikelihood") = loglikelihood);
	}
}
*/


// forward-project marginal likelihoods of a subset of nodes to their descending tips by applying the exponentiated transition rate matrix
// [[Rcpp::export]]
NumericVector apply_fixed_rate_Markov_model_to_missing_clades_CPP(	const long			Ntips,
																	const long 			Nnodes,
																	const long			Nedges,
																	const long			Nstates,
																	const IntegerVector &tree_edge,				// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
																	const NumericVector &edge_length, 			// (INPUT) 1D array of size Nedges, or an empty std::vector (all branches have length 1)
																	const NumericVector	&transition_matrix,		// (INPUT) 2D matrix of size Nstates x Nstates, in row-major format. Transition rate matrix of the Markov model. Sum-per-row should be zero.
																	const double		exponentiation_accuracy,// (INPUT) maximum allowed error when exponentiating the transition matrix, in terms of the Hilbert-Schmidt L2 norm.
																	const long			max_polynomials,		// (INPUT) maximum possible number of polynomials to use for exponentiating the transition_matrix, regardless of the pursued accuracy epsilon. Used as safety vault, but may break the guaranteed accuracy. A value ~100 is usually enough.
																	LogicalVector		likelihoods_known,		// (INPUT) 1D array of size Nclades, indicating whether the likelihoods for a particular clade are known (1) or unknown/to be determined (0).
																	NumericVector 		likelihoods){			// (INPUT) 2D matrix of size Nclades x Nstates, in row-major format. Likelihoods of each state in each clade (tip & node) of the tree.

	// determine root
	const long root = get_root_clade(Ntips, Nnodes, Nedges, tree_edge);

	// get tree traversal route (tips --> root)
	std::vector<long> traversal_queue, traversal_node2first_edge, traversal_node2last_edge, traversal_edges;
	get_tree_traversal_root_to_tips(	Ntips,
										Nnodes,
										Nedges,
										root,
										tree_edge,
										false,
										false,
										traversal_queue,
										traversal_node2first_edge,
										traversal_node2last_edge,
										traversal_edges,
										false,
										"");
										
	// prepare data structures for exponentiations of transition matrix
	const double max_edge_length = (edge_length.size()==0 ? 1.0 : get_array_max(edge_length));
	std::vector<double> transition_polynomials, transition_polynomial_norms, exponentiation_balances;
	const long min_polynomials = min_polynomials_for_positive_exponential_of_irreducible_matrix(Nstates, transition_matrix);
	long Npolynomials;
	long exponentiation_scaling_power;
	calculate_balanced_matrix_polynomials(	Nstates,
											std::vector<double>(transition_matrix.begin(), transition_matrix.end()),
											max_edge_length,
											exponentiation_accuracy,
											min_polynomials,
											max_polynomials,
											transition_polynomials,
											transition_polynomial_norms,
											Npolynomials,
											exponentiation_balances,
											exponentiation_scaling_power);

	// calculate unknown likelihoods based on parents with known likelihoods (traverse root --> tips)
	std::vector<double> expQ, Y;
	for(long q=0, node, clade; q<traversal_queue.size(); ++q){
		clade = traversal_queue[q];
		node = clade - Ntips;
		for(long e=traversal_node2first_edge[node], edge, child; e<=traversal_node2last_edge[node]; ++e){
			edge  = traversal_edges[e];
			child = tree_edge[edge*2+1];
			if(likelihoods_known[child]) continue; // skip children with known likelihoods
			// exponentiate transition matrix along this edge
			get_matrix_exponential_using_balanced_polynomials(	Nstates,
																Npolynomials,
																transition_polynomials,
																transition_polynomial_norms,
																(edge_length.size()==0 ? 1.0 : edge_length[edge])/max_edge_length,
																exponentiation_accuracy,
																min_polynomials,
																exponentiation_balances,
																exponentiation_scaling_power,
																expQ);
			// propagate clade's likelihoods to child by multiplying from the right with exponentiated transition matrix
			multiply_vector_with_matrix(Nstates, Nstates, &expQ[0], &likelihoods[clade*Nstates], Y);
			for(long s=0; s<Nstates; ++s) likelihoods[child*Nstates+s] = Y[s];
			likelihoods_known[child] = true;
		}
	}
										
	return likelihoods;														
}






// [[Rcpp::export]]
NumericVector apply_MPR_to_missing_clades_CPP(	const long			Ntips,
												const long 			Nnodes,
												const long			Nedges,
												const long			Nstates,
												const IntegerVector &tree_edge,				// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
												LogicalVector		likelihoods_known,		// (INPUT) 1D array of size Nclades, indicating whether the likelihoods for a particular clade are known (1) or unknown/to be determined (0).
												NumericVector 		likelihoods){			// (INPUT) 2D matrix of size Nclades x Nstates, in row-major format. Likelihoods of each state in each clade (tip & node) of the tree.

	// determine root
	const long root = get_root_clade(Ntips, Nnodes, Nedges, tree_edge);

	// get tree traversal route (tips --> root)
	std::vector<long> traversal_queue, traversal_node2first_edge, traversal_node2last_edge, traversal_edges;
	get_tree_traversal_root_to_tips(	Ntips,
										Nnodes,
										Nedges,
										root,
										tree_edge,
										false,
										false,
										traversal_queue,
										traversal_node2first_edge,
										traversal_node2last_edge,
										traversal_edges,
										false,
										"");

	// set unknown likelihoods based on parents with known likelihoods (traverse root --> tips)
	for(long q=0, node, clade; q<traversal_queue.size(); ++q){
		clade = traversal_queue[q];
		node = clade - Ntips;
		for(long e=traversal_node2first_edge[node], edge, child; e<=traversal_node2last_edge[node]; ++e){
			edge  = traversal_edges[e];
			child = tree_edge[edge*2+1];
			if(likelihoods_known[child]) continue; // skip children with known likelihoods
			// propagate clade's likelihoods to child
			for(long s=0; s<Nstates; ++s) likelihoods[child*Nstates+s] = likelihoods[clade*Nstates+s];
			likelihoods_known[child] = true;
		}
	}
										
	return likelihoods;														
}





// calculate (or update) quadratic parameters for squared-change parsimony [Maddison 1991], for a single focal node (based on its children)
template<class ARRAY_TYPE>
void aux_get_quadratic_parameters_for_squared_change_parsimony(	const long				Ntips,
																const long 				Nnodes,
																const ARRAY_TYPE 		&tree_edge,							// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
																const NumericVector 	&edge_length,						// (INPUT) 1D array of size Nedges, or an empty std::vector (all edges have length 1)
																const NumericVector		&tip_states,						// (INPUT) 1D array of size Ntips
																const std::vector<long>	&node2first_edge,					// (INPUT) 1D array of size Nnodes, with values in 0:(Nedges-1), mapping nodes to their first outgoing edge (or incoming/outgoing edge if edge_mapping_inout==true).
																const std::vector<long>	&node2last_edge,					// (INPUT) 1D array of size Nnodes, with values in 0:(Nedges-1), mapping nodes to their last outgoing edge (or incoming/outgoing edge if edge_mapping_inout==true).
																const std::vector<long>	&edge_mapping,						// (INPUT) 1D array of size Nedges (or 2*Nedges if edge_mapping_inout==true), with values in 0:(Nedges-1). Maps internal edge indices (i.e. as listed in node2first_edge[] and node2last_edge[]) to original edge indices.
																const bool				edge_mapping_inout,					// (INPUT) whether the provided edge mapping lists inout edges. If false, then edge mapping must only list outgoing edges.
																const double			edge_length_epsilon,				// (INPUT) replacement length to use for edges with length 0, so that children are weighted much more strongly compared to children with edge_length>0
																const long				focal_node,							// (INPUT) focal node for which to calculate/update the quadratic parameters
																std::vector<double>		&quadratic_parameters_per_node){	// (INPUT/OUTPUT) 2D array of size Nnodes x 3, in row-major format, storing the squared-change parsimony quadratic parametsrs [Maddison 1991]. For focal_node, these parameters will be updated based on its children.
		const long focal_clade = focal_node + Ntips;
		double p1, p2, p3, length;
		long edge, child;
		
		// initialize
		quadratic_parameters_per_node[focal_node*3+0] = 0;
		quadratic_parameters_per_node[focal_node*3+1] = 0;
		quadratic_parameters_per_node[focal_node*3+2] = 0;
		
		// traverse through the focal node's children
		for(long e=node2first_edge[focal_node]; e<=node2last_edge[focal_node]; ++e){
			edge  	= edge_mapping[e];
			if(edge_mapping_inout && (tree_edge[2*edge+0]!=focal_clade)) continue; // this is not an outgoing edge
			child 	= tree_edge[2*edge+1];
			length 	= (edge_length.size()==0 ? 1.0 : edge_length[edge]);
			if(length==0) length = edge_length_epsilon;
			if(child>=Ntips){
				p1 = quadratic_parameters_per_node[(child-Ntips)*3+0];
				p2 = quadratic_parameters_per_node[(child-Ntips)*3+1];
				p3 = quadratic_parameters_per_node[(child-Ntips)*3+2];
			}
			quadratic_parameters_per_node[focal_node*3+0] += (child<Ntips ? 1.0/length						: p1/(length*p1+1));
			quadratic_parameters_per_node[focal_node*3+1] += (child<Ntips ? -2*tip_states[child]/length 	: p2/(length*p1+1));
			quadratic_parameters_per_node[focal_node*3+2] += (child<Ntips ? SQR(tip_states[child])/length 	: p3 - length*SQR(p2)/(4*(length*p1+1)));
		}
}




// Reconstruct of continuous ancestral states via squared change parsimony
// This is related to Phylogenetically Independent Contrasts ASR), with the difference that PIC reconstruction for a node only takes into account the subtree descending from a node.
// Note that the "mean node value" X_i of a node used to calculate the phylogenetic independent contrasts (PICs) by Felsenstein, 
//		is in fact the value that would have been reconstructed at P by (weighted) squared-change parsimony were P's clade is the whole tree; 
//		that is, it minimizes locally (over P's clade) the sum of (weighted) squared changes [as explained by Maddison 1991].
//		For the root, this is also the global optimum.
//		Hence, to obtain ancestral states for non-root nodes, you need to reroot at each node.
//		Maddison (1991) provides an alternative postorder-traversal algorithm to Felsenstein (whose original intention was not ASR), for arriving at the same local estimates for each node (=global estimate for the root) 
//		by keeping track of "quadratic parameters" at each tip/node. The quadratic parameters of each node can be calculated purely based on the  quadratic parameters of its children, and branch lengths can be included for weighting.
//		It turns out that the calculation of quadratic parameters sensu Maddison is more easily generalizable to multifurcations.
// To obtain the classical PIC estimates, set global=false. To obtain Maddison's global estimates (via rerooting) set global=true.
// Literature:
//    Felsenstein (1985). Phylogenies and the Comparative Method. The American Naturalist. 125:1-15.
//    Maddison (1991). Squared-change parsimony reconstructions of ancestral states for continuous-valued characters on a phylogenetic tree. Systematic Zoology. 40:304-314.
// Requirements:
//   Tree can include multi- and mono-furcations.
//   Tree must be rooted. Root will be determined automatically as the node with no parent.
// [[Rcpp::export]]
Rcpp::List ASR_via_squared_change_parsimony_CPP(const long			Ntips,
												const long 			Nnodes,
												const long			Nedges,
												const IntegerVector &tree_edge,			// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
												const NumericVector &edge_length,		// (INPUT) 1D array of size Nedges, or an empty std::vector (all edges have length 1)
												const NumericVector	&tip_states,		// (INPUT) numeric states of the trait at each tip
												bool				global){			// (INPUT) if true, then global squared-change parsimony state estimates are returned (as recommended by Maddison). This requires rerooting at each node and updating the quadratic parameters. If false, then the local estimate of each node (i.e. only considering its descending subtree) are returned. This corresponds to classical PIC.
	long node, clade;
	const double edge_length_epsilon = RELATIVE_EPSILON * get_array_nonzero_min(edge_length); // substitute to be used for zero edge lengths
	std::vector<double> ancestral_states(Nnodes); // will be populated later on
	
	// determine root
	const long root = get_root_clade(Ntips, Nnodes, Nedges, tree_edge);

	// prepare tree traversal route (root-->tips)
	std::vector<long> traversal_queue, traversal_node2first_edge, traversal_node2last_edge, traversal_edges;
	get_tree_traversal_root_to_tips(Ntips,
									Nnodes,
									Nedges,
									root,
									tree_edge,
									false,	// don't include tips
									false,	// edge mappings are not pre-calculated
									traversal_queue,
									traversal_node2first_edge,
									traversal_node2last_edge,	
									traversal_edges,
									false,
									"");
								
	// calculate quadratic parameters [Maddison 1991] in a postorder traversal (tips-->root)
	std::vector<double> quadratic_parameters_per_node(Nnodes*3); // 2D array of size Nnodes x 3, in row-major format
	for(long q=traversal_queue.size()-1; q>=0; --q){
		aux_get_quadratic_parameters_for_squared_change_parsimony(	Ntips,
																	Nnodes,
																	tree_edge,
																	edge_length,
																	tip_states,
																	traversal_node2first_edge,
																	traversal_node2last_edge,
																	traversal_edges,
																	false, // egde mappings are not inout
																	edge_length_epsilon,
																	traversal_queue[q]-Ntips, // focal node for which to calculate quadratic parameters (based on its children)
																	quadratic_parameters_per_node);
	}
	const double TSS = quadratic_parameters_per_node[(root-Ntips)*3+2] - SQR(quadratic_parameters_per_node[(root-Ntips)*3+1])/(4*quadratic_parameters_per_node[(root-Ntips)*3+0]); // minimized total sum of squared changes over the tree [Maddison 1991, Formula 7]
	
	if(!global){
		// return the local state estimate for each node, i.e. only taking into account its descending subtree
		// this is equivalent to phylogenetic independent contrasts
		// however we used Maddison's quadratic parameters because it's convenient and has been generalized to multifurcations [Maddison 1991]
		for(node=0; node<Nnodes; ++node){
			ancestral_states[node] = -quadratic_parameters_per_node[node*3+1]/(2*quadratic_parameters_per_node[node*3+0]);
		}
		
	}else{
		// calculate global estimates: This requires rerooting at each node and updating affected quadratic parameters.
		
		// prepare some data structures for tree traversal, which will be updated everytime we reroot
		// Note that the original edge mappings (e.g. traversal_node2first_edge[]) will no longer be consistent with current_tree_edge after rerooting
		// Notation: current_(..) refers to tree-access structures that are updated at each rerooting. They will be consistent with each other, but not necessarily with the original tree structure.
		long current_root = root;
		std::vector<long> current_tree_edge(tree_edge.begin(), tree_edge.end());
		std::vector<long> current_incoming_edge_per_clade; // will change with each rerooting
		get_incoming_edge_per_clade(Ntips, Nnodes, Nedges, current_tree_edge, current_incoming_edge_per_clade);
		
		// get inout edges for each node
		// these will be invariant to rerootings
		std::vector<long> node2first_inout_edge, node2last_inout_edge, inout_edges;
		get_inout_edges_per_node(Ntips, Nnodes, Nedges, tree_edge, node2first_inout_edge, node2last_inout_edge, inout_edges);

		// calculate depth-first-search traversal route
		// this minimizes the distance between successive old_root-new_root pairs and thus the computation required to update the posteriors on each rerooting
		// similarly, depth_first_search_queue[] will still list successive nodes, but not necessarily in downstream order
		std::vector<long> depth_first_search_queue;
		get_tree_traversal_depth_first_search(	Ntips,
												Nnodes,
												Nedges,
												root,
												tree_edge,
												false,	// don't include tips
												true, 	// edge mappings are pre-calculated
												depth_first_search_queue,
												traversal_node2first_edge,
												traversal_node2last_edge,
												traversal_edges);
	
		// estimate ancestral states at each node from its quadratic parameters, as if that node was the root
		// While rerooting at each node, update the tree-access structures and update the quadratic parameters of each node
		// Edge directions will change, but tip/node/edge indices remain the same
		ancestral_states[root-Ntips] = -quadratic_parameters_per_node[(root-Ntips)*3+1]/(2*quadratic_parameters_per_node[(root-Ntips)*3+0]); // for the root, the local estimate is also the global estimate
		for(long q=0, new_root; q<depth_first_search_queue.size(); ++q){
			new_root = depth_first_search_queue[q];
			if(new_root==current_root) continue; // nothing to do
			reroot_tree_at_node(Ntips, Nnodes, Nedges, current_root, new_root, current_tree_edge, current_incoming_edge_per_clade); // reroot (this changes edge directions, but keeps tip/node/edge indices the same)
	
			// update quadratic parameters for nodes that have been traversed by the rerooting
			clade = current_root;
			while(true){
				// re-calculate the quadratic parameters for this clade (based on the quadratic parameters of its children)
				aux_get_quadratic_parameters_for_squared_change_parsimony(	Ntips,
																			Nnodes,
																			current_tree_edge,
																			edge_length,
																			tip_states,
																			node2first_inout_edge,
																			node2last_inout_edge,
																			inout_edges,
																			true, // edge mappings are inout
																			edge_length_epsilon,
																			clade-Ntips, // focal node for which to calculate quadratic parameters (based on its children)
																			quadratic_parameters_per_node);
				if(clade!=new_root) clade = current_tree_edge[current_incoming_edge_per_clade[clade]*2+0]; // move on to parent
				else break;
			}

			// calculate the global (=local) estimate of the root from its quadratic parameters
			ancestral_states[new_root-Ntips] = -quadratic_parameters_per_node[(new_root-Ntips)*3+1]/(2*quadratic_parameters_per_node[(new_root-Ntips)*3+0]);
			current_root = new_root;
			// abort if the user has interrupted the calling R program
			Rcpp::checkUserInterrupt();
		}
	}
	
	return Rcpp::List::create(	Rcpp::Named("TSS")				= TSS,	// total sum of squared changes
								Rcpp::Named("ancestral_states") = Rcpp::wrap(ancestral_states));
}




// based on ML estimates of a continuous trait for a subset of tips, extrapolate those estimates to the remaining tips under a Brownian Motion model
// This essentially sets the state of an unknown tip to the state of its most recent known ancestor.
// [[Rcpp::export]]
NumericVector apply_BM_parsimony_to_missing_clades_CPP(	const long			Ntips,
														const long 			Nnodes,
														const long			Nedges,
														const IntegerVector &tree_edge,			// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
														LogicalVector		states_known,		// (INPUT) 1D array of size Nclades, indicating whether the states for a particular clade are known/already estimated (1) or unknown/to be determined (0).
														NumericVector 		states){			// (INPUT) 1D array of size Nclades, listing the state of each clade in the tree. May contain NA/NaN in those cases where states_known[i]=false.

	// determine root
	const long root = get_root_clade(Ntips, Nnodes, Nedges, tree_edge);

	// get tree traversal route (tips --> root)
	std::vector<long> traversal_queue, traversal_node2first_edge, traversal_node2last_edge, traversal_edges;
	get_tree_traversal_root_to_tips(	Ntips,
										Nnodes,
										Nedges,
										root,
										tree_edge,
										false,
										false,
										traversal_queue,
										traversal_node2first_edge,
										traversal_node2last_edge,
										traversal_edges,
										false,
										"");

	// set unknown statets to the estimates state of their parent (traverse root --> tips)
	for(long q=0, node, clade; q<traversal_queue.size(); ++q){
		clade = traversal_queue[q];
		node = clade - Ntips;
		for(long e=traversal_node2first_edge[node], edge, child; e<=traversal_node2last_edge[node]; ++e){
			edge  = traversal_edges[e];
			child = tree_edge[edge*2+1];
			if(states_known[child]) continue; // skip children with known state
			// propagate clade's state to child
			states[child] 		= states[clade];
			states_known[child] = true;
		}
	}
										
	return states;														
}




#pragma mark -
#pragma mark Simulate models of trait evolution
#pragma mark








// Perform random simulations of a fixed-rates continuous-time Markov model of discrete character evolution
// Starting with a specified vector of root_probabilities, and moving from root to tips, each node is assigned a random state according to its parent's state and according to the markov transition matrix.
// Optionally, multiple independent simulations can be performed using the same model (e.g. as part of some Monte Carlo integration)
// [[Rcpp::export]]
Rcpp::List simulate_fixed_rates_Markov_model_CPP(	const long			Ntips,
													const long 			Nnodes,
													const long			Nedges,
													const long			Nstates,
													const IntegerVector &tree_edge,				// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
													const NumericVector &edge_length, 			// (INPUT) 1D array of size Nedges, or an empty std::vector (all branches have length 1)
													const NumericVector &transition_matrix,		// (INPUT) 2D array of size Nstates x Nstates, in row-major format. Transition-rate matrix Q in row-major format, i.e. Q[r*Nstates + c] is (r,c)-th-element of Q, which is the transition rate r-->c.
													const NumericVector	&root_probabilities,	// (INPUT) probability distribution of states at the root. sum(root_probabilities) must be 1.0.
													const bool			include_tips,			// (INPUT) include states for tips in the output
													const bool			include_nodes,			// (INPUT) include states for nodes in the output
													const long			Nsimulations){			// (INPUT) number of random simulations (draws) of the model on the tree. If 1, then a single simulation is performed, yielding a single random state for each node and/or tip.
	if((Nsimulations<=0) || ((!include_tips) && (!include_nodes))){
		return	Rcpp::List::create(	Rcpp::Named("tip_states")  = IntegerVector(),
									Rcpp::Named("node_states") = IntegerVector());
	}

	// get incoming edge for each clade
	std::vector<long> incoming_edge_per_clade;
	get_incoming_edge_per_clade(Ntips, Nnodes, Nedges, tree_edge, incoming_edge_per_clade);

	// find root based on mapping clade-->incoming_edge
	const long root = get_root_from_incoming_edge_per_clade(Ntips, tree_edge, incoming_edge_per_clade);

	// get tree traversal route (tips --> root)
	std::vector<long> traversal_queue, traversal_node2first_edge, traversal_node2last_edge, traversal_edges;
	get_tree_traversal_root_to_tips(	Ntips,
										Nnodes,
										Nedges,
										root,
										tree_edge,
										include_tips,
										false,
										traversal_queue,
										traversal_node2first_edge,
										traversal_node2last_edge,
										traversal_edges,
										false,
										"");
										
	// prepare data structures for exponentiations of transition matrix
	const double max_edge_length = (edge_length.size()==0 ? 1.0 : get_array_max(edge_length));
	std::vector<double> transition_polynomials, transition_polynomial_norms, exponentiation_balances;
	long Npolynomials;
	const double exponentiation_accuracy = 1e-4;
	const long NPmin = min_polynomials_for_positive_exponential_of_irreducible_matrix(Nstates, transition_matrix);
	const long NPmax = 1000;
	long exponentiation_scaling_power;
	calculate_balanced_matrix_polynomials(	Nstates,
											std::vector<double>(transition_matrix.begin(), transition_matrix.end()),
											max_edge_length,
											exponentiation_accuracy,
											NPmin,
											NPmax,
											transition_polynomials,
											transition_polynomial_norms,
											Npolynomials,
											exponentiation_balances,
											exponentiation_scaling_power);
	
	// traverse root-->tips and draw random states, conditional upon their parent's state
	vector<double> expQ;
	vector<long> tip_states, node_states;
	if(include_tips) tip_states.resize(Nsimulations*Ntips);
	node_states.resize(Nsimulations*Nnodes); // always store node states, since needed for moving root-->tips
	long clade, edge, parent, parent_state, state;
	for(long q=0; q<traversal_queue.size(); ++q){
		clade = traversal_queue[q];
		if(clade!=root){
			edge 	= incoming_edge_per_clade[clade];
			parent 	= tree_edge[edge*2+0];
			// exponentiate transition matrix along incoming edge
			get_matrix_exponential_using_balanced_polynomials(	Nstates,
																Npolynomials,
																transition_polynomials,
																transition_polynomial_norms,
																(edge_length.size()==0 ? 1.0 : edge_length[edge])/max_edge_length,
																exponentiation_accuracy,
																NPmin,
																exponentiation_balances,
																exponentiation_scaling_power,
																expQ);
		}
		for(long r=0; r<Nsimulations; ++r){
			if(clade==root){
				state = random_int_from_distribution<double>(&root_probabilities[0], Nstates);
			}else{
				// use row of exponentiated transition matrix corresponding to the parent's state, as probability vector for the child's state
				parent_state = node_states[r*Nnodes + (parent-Ntips)];
				state = random_int_from_distribution<double>(&expQ[parent_state*Nstates+0], Nstates);
			}
			if((clade<Ntips) && include_tips) tip_states[r*Ntips + clade] = state;
			else if(clade>=Ntips) node_states[r*Nnodes + (clade-Ntips)] = state;
		}
	}
	if(!include_nodes) node_states.clear(); // clear memory if content is not to be returned
			
	return Rcpp::List::create(	Rcpp::Named("tip_states")  = Rcpp::wrap(tip_states),
								Rcpp::Named("node_states") = Rcpp::wrap(node_states));
}





// simulate a specific Ornstein-Uhlenbeck model of continuous trait evolution on a tree, starting from the root and moving towards the tips
// [[Rcpp::export]]
Rcpp::List simulate_Ornstein_Uhlenbeck_model_CPP(	const long			Ntips,
													const long 			Nnodes,
													const long			Nedges,
													const IntegerVector &tree_edge,				// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
													const NumericVector &edge_length, 			// (INPUT) 1D array of size Nedges, or an empty std::vector (all branches have length 1)
													const double		stationary_mean,		// (INPUT) mean of the stationary distribution
													const double		stationary_std,			// (INPUT) standard deviation of the stationary distribution
													const double 		decay_rate,				// (INPUT) exponential decay rate (or equilibration rate), in units 1/edge_length
													const bool			include_tips,			// (INPUT) include states for tips in the output
													const bool			include_nodes,			// (INPUT) include states for nodes in the output
													const long			Nsimulations){				// (INPUT) number of random simulations (draws) of the model on the tree. If 1, then a single simulation is performed, yielding a single random state for each node and/or tip.
	if((Nsimulations<=0) || ((!include_tips) && (!include_nodes))){
		// nothing to do
		return	Rcpp::List::create(	Rcpp::Named("tip_states")  = NumericVector(),
									Rcpp::Named("node_states") = NumericVector());
	}

	// get incoming edge for each clade
	std::vector<long> incoming_edge_per_clade;
	get_incoming_edge_per_clade(Ntips, Nnodes, Nedges, tree_edge, incoming_edge_per_clade);

	// find root based on mapping clade-->incoming_edge
	const long root = get_root_from_incoming_edge_per_clade(Ntips, tree_edge, incoming_edge_per_clade);

	// get tree traversal route (tips --> root)
	std::vector<long> traversal_queue, traversal_node2first_edge, traversal_node2last_edge, traversal_edges;
	get_tree_traversal_root_to_tips(	Ntips,
										Nnodes,
										Nedges,
										root,
										tree_edge,
										include_tips,
										false,
										traversal_queue,
										traversal_node2first_edge,
										traversal_node2last_edge,
										traversal_edges,
										false,
										"");
										
	
	// traverse root-->tips and draw random states, conditional upon their parent's state
	std::vector<double> tip_states, node_states;
	if(include_tips) tip_states.resize(Nsimulations*Ntips);
	node_states.resize(Nsimulations*Nnodes); // always store node states, since needed for moving root-->tips
	long clade, edge, parent;
	double parent_state, state;
	for(long r=0; r<Nsimulations; ++r){
		for(long q=0; q<traversal_queue.size(); ++q){
			clade = traversal_queue[q];
			if(clade==root){
				// set root to random state drawn from the stationary distribution
				state = stationary_mean + stationary_std*random_standard_normal();
			}else{
				edge 			= incoming_edge_per_clade[clade];
				parent 			= tree_edge[edge*2+0];
				parent_state 	= node_states[r*Nnodes + (parent-Ntips)];
				state 			= get_next_OU_sample(stationary_mean, decay_rate, stationary_std, (edge_length.size()==0 ? 1.0 : edge_length[edge]), parent_state);
			}
			if((clade<Ntips) && include_tips) tip_states[r*Ntips + clade] = state;
			else if(clade>=Ntips) node_states[r*Nnodes + (clade-Ntips)] = state;
		}
	}
	if(!include_nodes) node_states.clear(); // clear memory if content is not to be returned
		
	return Rcpp::List::create(	Rcpp::Named("tip_states")  = Rcpp::wrap(tip_states),
								Rcpp::Named("node_states") = Rcpp::wrap(node_states));
}






// simulate a Reflected Ornstein-Uhlenbeck (ROU) model of continuous trait evolution on a tree, starting from the root and moving towards the tips
// The ROU process is reflected at some minimum, which also happens to be its deterministic equilibrium.
// Note that in general the reflection point and deterministic equilibrium of an ROU may differ, but simulating such a more general process is much harder due to a reduced symmetry
// For details on the ROU process see:
//    http://www-bcf.usc.edu/~amyward/ROUproperties.pdf
//    http://hkumath.hku.hk/~jiansong/HLLS-9-10-12.pdf
// The root's state is drawn from the ROU's stationary distribution.
// [[Rcpp::export]]
Rcpp::List simulate_reflected_Ornstein_Uhlenbeck_model_CPP(	const long			Ntips,
															const long 			Nnodes,
															const long			Nedges,
															const IntegerVector &tree_edge,				// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
															const NumericVector &edge_length, 			// (INPUT) 1D array of size Nedges, or an empty std::vector (all branches have length 1)
															const double		reflection_point,		// (INPUT) the point of reflection, i.e. the ROU's minimum. This is also the deterministic equilibrium (see explanation above).
															const double		spread,					// (INPUT) standard deviation of the stationary distribution of the corresponding unreflected OU process
															const double 		decay_rate,				// (INPUT) exponential decay rate (or equilibration rate), in units 1/edge_length
															const bool			include_tips,			// (INPUT) include states for tips in the output
															const bool			include_nodes,			// (INPUT) include states for nodes in the output
															const long			Nsimulations){			// (INPUT) number of random simulations (draws) of the model on the tree. If 1, then a single simulation is performed, yielding a single random state for each node and/or tip.
	if((Nsimulations<=0) || ((!include_tips) && (!include_nodes))){
		// nothing to do
		return	Rcpp::List::create(	Rcpp::Named("tip_states")  = NumericVector(),
									Rcpp::Named("node_states") = NumericVector());
	}

	// get incoming edge for each clade
	std::vector<long> incoming_edge_per_clade;
	get_incoming_edge_per_clade(Ntips, Nnodes, Nedges, tree_edge, incoming_edge_per_clade);

	// find root based on mapping clade-->incoming_edge
	const long root = get_root_from_incoming_edge_per_clade(Ntips, tree_edge, incoming_edge_per_clade);

	// get tree traversal route (tips --> root)
	std::vector<long> traversal_queue, traversal_node2first_edge, traversal_node2last_edge, traversal_edges;
	get_tree_traversal_root_to_tips(	Ntips,
										Nnodes,
										Nedges,
										root,
										tree_edge,
										include_tips,
										false,
										traversal_queue,
										traversal_node2first_edge,
										traversal_node2last_edge,
										traversal_edges,
										false,
										"");
										
	
	// traverse root-->tips and draw random states, conditional upon their parent's state
	std::vector<double> tip_states, node_states;
	if(include_tips) tip_states.resize(Nsimulations*Ntips);
	node_states.resize(Nsimulations*Nnodes); // always store node states, since needed for moving root-->tips
	long clade, edge, parent;
	double parent_state, state;
	for(long r=0; r<Nsimulations; ++r){
		for(long q=0; q<traversal_queue.size(); ++q){
			clade = traversal_queue[q];
			if(clade==root){
				// set root to random state drawn from the stationary distribution
				state = reflection_point + abs(spread*random_standard_normal());
			}else{
				edge 			= incoming_edge_per_clade[clade];
				parent 			= tree_edge[edge*2+0];
				parent_state 	= node_states[r*Nnodes + (parent-Ntips)];
				state 			= reflection_point + abs(get_next_OU_sample(0, decay_rate, spread, (edge_length.size()==0 ? 1.0 : edge_length[edge]), parent_state-reflection_point));
			}
			if((clade<Ntips) && include_tips) tip_states[r*Ntips + clade] = state;
			else if(clade>=Ntips) node_states[r*Nnodes + (clade-Ntips)] = state;
		}
	}
	if(!include_nodes) node_states.clear(); // clear memory if content is not to be returned
		
	return Rcpp::List::create(	Rcpp::Named("tip_states")  = Rcpp::wrap(tip_states),
								Rcpp::Named("node_states") = Rcpp::wrap(node_states));
}


// simulate a Brownian motion model of continuous trait evolution for a scalar continuous trait, starting from the root and moving towards the tips
// The BM model is specified by means of its diffusivity D
// [[Rcpp::export]]
Rcpp::List simulate_scalar_Brownian_motion_model_CPP(	const long			Ntips,
														const long 			Nnodes,
														const long			Nedges,
														const IntegerVector &tree_edge,				// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
														const NumericVector &edge_length, 			// (INPUT) 1D array of size Nedges, or an empty std::vector (all branches have length 1)
														const NumericVector	&root_states,			// (INPUT) 1D array of arbitrary size, specifying root states for each simulation. If smaller than Nsimulations, values are recycled in rotation. If empty, zero is used as root state.
														const double		diffusivity,			// (INPUT) diffusion coefficient of the BM model, dX = sqrt(2D) * dB
														const bool			include_tips,			// (INPUT) include states for tips in the output
														const bool			include_nodes,			// (INPUT) include states for nodes in the output
														const long			Nsimulations){			// (INPUT) number of random simulations (draws) of the model on the tree. If 1, then a single simulation is performed, yielding a single random state for each node and/or tip.
	if((Nsimulations<=0) || ((!include_tips) && (!include_nodes))){
		return	Rcpp::List::create(	Rcpp::Named("tip_states")  = NumericVector(),
									Rcpp::Named("node_states") = NumericVector());
	}

	// get incoming edge for each clade
	std::vector<long> incoming_edge_per_clade;
	get_incoming_edge_per_clade(Ntips, Nnodes, Nedges, tree_edge, incoming_edge_per_clade);

	// find root based on mapping clade-->incoming_edge
	const long root = get_root_from_incoming_edge_per_clade(Ntips, tree_edge, incoming_edge_per_clade);

	// get tree traversal route (tips --> root)
	std::vector<long> traversal_queue, traversal_node2first_edge, traversal_node2last_edge, traversal_edges;
	get_tree_traversal_root_to_tips(	Ntips,
										Nnodes,
										Nedges,
										root,
										tree_edge,
										include_tips,
										false,
										traversal_queue,
										traversal_node2first_edge,
										traversal_node2last_edge,
										traversal_edges,
										false,
										"");
										
	
	// traverse root-->tips and draw random states, conditional upon their parent's state
	std::vector<double> tip_states, node_states;
	if(include_tips) tip_states.resize(Nsimulations*Ntips);
	node_states.resize(Nsimulations*Nnodes); // always store node states, since needed for moving root-->tips
	long clade, edge, parent;
	double parent_state, state;
	for(long r=0; r<Nsimulations; ++r){
		for(long q=0; q<traversal_queue.size(); ++q){
			clade = traversal_queue[q];
			if(clade==root){
				state = (root_states.size()==0 ? 0.0 : root_states[r % root_states.size()]);
			}else{
				edge 			= incoming_edge_per_clade[clade];
				parent 			= tree_edge[edge*2+0];
				parent_state 	= node_states[r*Nnodes + (parent-Ntips)];
				state 			= parent_state + sqrt(2 * diffusivity * (edge_length.size()==0 ? 1.0 : edge_length[edge])) * random_standard_normal();
			}
			if((clade<Ntips) && include_tips) tip_states[r*Ntips + clade] = state;
			else if(clade>=Ntips) node_states[r*Nnodes + (clade-Ntips)] = state;
		}
	}
	if(!include_nodes) node_states.clear(); // clear memory if content is not to be returned
		
	return Rcpp::List::create(	Rcpp::Named("tip_states")  = Rcpp::wrap(tip_states),
								Rcpp::Named("node_states") = Rcpp::wrap(node_states));
}




// simulate a Brownian motion model of co-evolution of multiple continuous traits on a tree, starting from the root and moving towards the tips.
// Model:
//   dX = sigma * dW
// where sigma is a noise amplitude matrix sigma of size Ntraits x Ndegrees (number of traits x degrees of freedom)
// Note that sigma * sigma^T = 2*diffusivity, i.e. sigma is related to the Cholesky decomposition of the diffusivity matrix
// [[Rcpp::export]]
Rcpp::List simulate_multivariate_Brownian_motion_model_CPP(	const long					Ntips,
															const long 					Nnodes,
															const long					Nedges,
															const long					Ntraits,				// (INPUT) number of traits, i.e. dimensionality of the space in which the Brownian motion random walk travels
															const long 					Ndegrees,				// (INPUT) degrees of freedom, i.e. number of independent Brownian motions driving the random walk
															const IntegerVector 		&tree_edge,				// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
															const NumericVector			&edge_length, 			// (INPUT) 1D array of size Nedges, or an empty std::vector (all branches have length 1)
															const NumericVector			&root_states,			// (INPUT) 2D array of size NR x Ntraits in row-major format, where NR can be arbitrary, specifying root states for each simulation. If NR is smaller than Nsimulations, values are recycled in rotation. If empty, zero is used as root state for all traits.
															const std::vector<double>	&sigma,					// (INPUT) 2D array of size Ntraits x Ndegrees, in row-major format. Noise amplitude matrix. If Ntraits = Ndegrees, then sigma is related to the cholesky decomposition of the diffusivity matrix, i.e. sigma satisfies: 2 * diffusivity = sigma * sigma^T.
															const bool					include_tips,			// (INPUT) include states for tips in the output
															const bool					include_nodes,			// (INPUT) include states for nodes in the output
															const long					Nsimulations){			// (INPUT) number of random simulations (draws) of the model on the tree. If 1, then a single simulation is performed, yielding a single random state for each node and/or tip.
	if((Nsimulations<=0) || ((!include_tips) && (!include_nodes))){
		return	Rcpp::List::create(	Rcpp::Named("tip_states")  = NumericVector(),
									Rcpp::Named("node_states") = NumericVector());
	}
	const long Nroot_states = root_states.size()/Ntraits;

	// get incoming edge for each clade
	std::vector<long> incoming_edge_per_clade;
	get_incoming_edge_per_clade(Ntips, Nnodes, Nedges, tree_edge, incoming_edge_per_clade);

	// find root based on mapping clade-->incoming_edge
	const long root = get_root_from_incoming_edge_per_clade(Ntips, tree_edge, incoming_edge_per_clade);

	// get tree traversal route (tips --> root)
	std::vector<long> traversal_queue, traversal_node2first_edge, traversal_node2last_edge, traversal_edges;
	get_tree_traversal_root_to_tips(	Ntips,
										Nnodes,
										Nedges,
										root,
										tree_edge,
										include_tips,
										false,
										traversal_queue,
										traversal_node2first_edge,
										traversal_node2last_edge,
										traversal_edges,
										false,
										"");
										
	
	// traverse root-->tips and draw random states, conditional upon their parent's state
	// node_states[] will be a 3D array in layer-row-major format, i.e. indexed as simulation*Nnodes*Ntraits+node*Ntraits+trait
	// similarly for tip_states
	std::vector<double> tip_states, node_states, state(Ntraits), standard_normals(Ndegrees);
	if(include_tips) tip_states.resize(Nsimulations*Ntips*Ntraits);
	node_states.resize(Nsimulations*Nnodes*Ntraits); // always store node states, since needed for moving root-->tips
	long clade, edge, parent, trait, length;
	for(long r=0; r<Nsimulations; ++r){
		for(long q=0; q<traversal_queue.size(); ++q){
			clade = traversal_queue[q];
			if(clade==root){
				for(trait=0; trait<Ntraits; ++trait) state[trait] = (Nroot_states==0 ? 0.0 : root_states[(r % Nroot_states)*Ntraits+trait]);
			}else{
				edge	= incoming_edge_per_clade[clade];
				parent	= tree_edge[edge*2+0];
				length	= (edge_length.size()==0 ? 1.0 : edge_length[edge]);
				// generate multivariate (potentially correlated) normally distributed numbers as:
				//    mean + sigma * standard_normals
				for(long d=0; d<Ndegrees; ++d){
					standard_normals[d]	= random_standard_normal();
				}
				multiply_matrix_with_vector(Ntraits, Ndegrees, sigma, standard_normals, state);
				for(trait=0; trait<Ntraits; ++trait){
					state[trait] = node_states[r*Nnodes*Ntraits + (parent-Ntips)*Ntraits+trait] + sqrt(length) * state[trait];
				}
			}
			for(trait=0; trait<Ntraits; ++trait){
				if((clade<Ntips) && include_tips) tip_states[r*Ntips*Ntraits + clade*Ntraits + trait] = state[trait];
				else if(clade>=Ntips) node_states[r*Nnodes*Ntraits + (clade-Ntips)*Ntraits + trait]   = state[trait];
			}
		}
	}
	if(!include_nodes) node_states.clear(); // clear memory if content is not to be returned
		
	return Rcpp::List::create(	Rcpp::Named("tip_states")  = Rcpp::wrap(tip_states),	// 3D array of size Nsimulations x Ntips x Ntraits, in layer-row-major format, i.e. indexed as simulation*Ntips*Ntraits+tip*Ntraits+trait
								Rcpp::Named("node_states") = Rcpp::wrap(node_states));	// 3D array of size Nsimulations x Nnodes x Ntraits, in layer-row-major format, i.e. indexed as simulation*Nnodes*Ntraits+node*Ntraits+trait
}





// WARNING: THIS FUNCTION IS NOT DEBUGGED YET. IT LIKELY WORKS FINE.
// Simulate a neutrally evolving gene along a phylogenetic tree.
// Each site of the gene is subject to a fixed mutation rate (mutations per site per edge_length_unit).
// Mutations are assumed to be independent for each site, and happen at an exponential rate, at equal probability to each state (all-to-all).
// Optionally, distances between alleles can be used to calculate new edge lengths (edge length = number of site differences between parent & child allele)
// [[Rcpp::export]]
Rcpp::List simulate_neutral_gene_evolution_CPP(	const long				Ntips,
												const long 				Nnodes,
												const long				Nedges,
												const long				Nsites,					// (INPUT) number of sites (e.g. nucleotides) at which the gene can vary neutrally
												const long 				Nstates,				// (INPUT) number of states each site can take (e.g. 4 for nucleotides). States are indexed 0,..,Nstates-1
												const IntegerVector 	&tree_edge,				// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
												const NumericVector		&edge_length, 			// (INPUT) 1D array of size Nedges, or an empty std::vector (all branches have length 1)
												const IntegerVector		&root_states,			// (INPUT) 2D array of size NR x Nsites in row-major format, with values in 0,..,Nstates-1, where NR can be arbitrary, specifying root states for each simulation. If NR is smaller than Nsimulations, values are recycled in rotation. If empty, zero is used as root state for all traits.
												const double			mutation_rate,			// (INPUT) mutation probability rate (mutations per site per edge_length_unit). 
												const bool				include_tips,			// (INPUT) include states for tips in the output
												const bool				include_nodes,			// (INPUT) include states for nodes in the output
												const bool				include_gene_distances,	// (INPUT) if true, then for each edge the distance between the parent & child allele is returned (for each simulation). This may be useful if you want to replace the original edge lengths with new lengths corresponding to the evolved gene.
												const long				Nsimulations){			// (INPUT) number of random simulations (draws) of the model on the tree. If 1, then a single simulation is performed, yielding a single random state for each node and/or tip.
	if((Nsimulations<=0) || ((!include_tips) && (!include_nodes))){
		return	Rcpp::List::create(	Rcpp::Named("tip_states")  = NumericVector(),
									Rcpp::Named("node_states") = NumericVector(),
									Rcpp::Named("gene_distances") = NumericVector());
	}
	const long Nroot_states = root_states.size()/Nsites;
	
	// get incoming edge for each clade
	std::vector<long> incoming_edge_per_clade;
	get_incoming_edge_per_clade(Ntips, Nnodes, Nedges, tree_edge, incoming_edge_per_clade);

	// find root based on mapping clade-->incoming_edge
	const long root = get_root_from_incoming_edge_per_clade(Ntips, tree_edge, incoming_edge_per_clade);

	// get tree traversal route (tips --> root)
	std::vector<long> traversal_queue, traversal_node2first_edge, traversal_node2last_edge, traversal_edges;
	get_tree_traversal_root_to_tips(	Ntips,
										Nnodes,
										Nedges,
										root,
										tree_edge,
										include_tips,
										false,
										traversal_queue,
										traversal_node2first_edge,
										traversal_node2last_edge,
										traversal_edges,
										false,
										"");	
	
	// traverse root-->tips and draw random states, conditional upon their parent's state
	// node_states[] will be a 3D array in layer-row-major format, i.e. indexed as simulation*Nnodes*Nsites+node*Nsites+site
	// similarly for tip_states
	std::vector<long> tip_states, node_states, gene(Nsites);
	if(include_tips) tip_states.resize(Nsimulations*Ntips*Nsites);
	node_states.resize(Nsimulations*Nnodes*Nsites); // always store node states, since needed for moving root-->tips
	std::vector<long> gene_distances(include_gene_distances ? Nsimulations*Nedges : 0, 0);
	long clade, edge, parent, site, length;
	for(long r=0; r<Nsimulations; ++r){
		for(long q=0; q<traversal_queue.size(); ++q){
			clade = traversal_queue[q];
			if(clade==root){
				for(site=0; site<Nsites; ++site){
					gene[site] = (Nroot_states==0 ? 0l : root_states[(r % Nroot_states)*Nsites+site]);
				}
			}else{
				edge	= incoming_edge_per_clade[clade];
				parent	= tree_edge[edge*2+0];
				length	= (edge_length.size()==0 ? 1.0 : edge_length[edge]);
				// mutate each site at Poisson probability rate
				for(site=0; site<Nsites; ++site){
					const double q = uniformWithinInclusiveLeft(0,1);
					const double lambda = mutation_rate*length;
					const double parent_state = node_states[r*Nnodes*Nsites + (parent-Ntips)*Nsites+site];
					if(q<=exp(-lambda)){
						// no mutation (same state as parent)
						gene[site] = parent_state;
					}else if(q<=exp(-lambda)*(1+lambda)){
						// exactly one mutation for this site (all states other than the original are equally probable)
						gene[site] = uniformIntWithin(0,Nstates-2);
						gene[site] = (gene[site]<parent_state ? gene[site] : 1+gene[site]);
					}else{
						// two or more mutations for this site (all states equally probable)
						gene[site] = uniformIntWithin(0,Nstates-1); 
					}
					if(include_gene_distances) gene_distances[r*Nedges+edge] += (gene[site]==parent_state ? 0 : 1);
				}
			}
			for(site=0; site<Nsites; ++site){
				if((clade<Ntips) && include_tips) tip_states[r*Ntips*Nsites + clade*Nsites + site] = gene[site];
				else if(clade>=Ntips) node_states[r*Nnodes*Nsites + (clade-Ntips)*Nsites + site]   = gene[site];
			}
		}
	}
	if(!include_nodes) node_states.clear(); // clear memory if content is not to be returned
	
	return Rcpp::List::create(	Rcpp::Named("tip_states")  		= Rcpp::wrap(tip_states),		// 3D array of size Nsimulations x Ntips x Nsites, in layer-row-major format, i.e. indexed as [simulation*Ntips*Nsites+tip*Nsites+site]
								Rcpp::Named("node_states") 		= Rcpp::wrap(node_states),		// 3D array of size Nsimulations x Nnodes x Nsites, in layer-row-major format, i.e. indexed as [simulation*Nnodes*Nsites+node*Nsites+site]
								Rcpp::Named("gene_distances")	= Rcpp::wrap(gene_distances));	// 2D array of size Nsimulations x Nedges, in row-major format, i.e. indexed as [simulation*Nedges+edge]
}




