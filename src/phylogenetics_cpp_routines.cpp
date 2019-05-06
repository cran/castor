/*
This C++ code library includes routines for efficiently scratching with huge phylogenetic trees in R (using the Rcpp interface).

Most code supports multifurcating trees, as well as trees containing monofurcations (i.e. some nodes having only one child).
In most cases, the input tree must be rooted.
The library is meant for large (>100,000 tips) trees, structured in the conventional "phylo" format in R.
The computational complexity of most routines is O(Nedges).

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
*/

#include <new>
#include <limits>
#include <cmath>
#include <string>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <complex>
#include <algorithm>
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
#define RELATIVE_EPSILON 1e-10
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
typedef std::vector<double> dvector;
typedef std::vector<long> lvector;

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



// returns the number of seconds since an arbitrary (but consistent) point in the past
// The timer is monotonic (i.e. not affected by manual changes in the system time), and is specific to the calling thread
// Note that for Windows this function is not thread-specific 
double get_thread_monotonic_walltime_seconds(){
	#if __MACH__ 
		mach_timebase_info_data_t info;
		int error_code = mach_timebase_info(&info);
		if (error_code != KERN_SUCCESS) return 0.0;
		return 1e-9 * mach_absolute_time() * double(info.numer)/double(info.denom);
	#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__linux__)
		// POSIX code
		// For details on clock_gettime() see: http://www.tin.org/bin/man.cgi?section=3&topic=clock_gettime
		struct timespec T;
		clock_gettime(CLOCK_MONOTONIC, &T); // Note that CLOCK_MONOTONIC_RAW is not available on all Linux distros
		return double(T.tv_sec) + 1e-9*T.tv_nsec;
	#elif defined(IS_WINDOWS)
		//return GetTickCount()*1e-6; // note that this is not thread-specific, but it's the best you can get on Windows. Update: It requires the windows.h library, which causes problems with Rcpp on CRAN.
		return clock()/double(CLOCKS_PER_SEC);
	#else
		return 0; // not implemented for other systems
	#endif
}



inline double string2Double(const string &number){
	return strtod(number.c_str(), NULL);
}

template<class TYPE> 
string makeString(const TYPE &data){
	ostringstream stream;
	stream << data;
	return stream.str();
}


// Formated string creation
string vstringprintf(const char *format, va_list args){
	va_list temp;
	va_copy(temp, args);
	char *buffer = new char[vsnprintf(NULL, 0, format, temp) + 1];
	va_end(temp);
	
	vsprintf(buffer, format, args);
	string s(buffer);
	delete [] buffer;
	return s;
}

// Formated string creation
string stringprintf(const char *format, ...){
	string s;
	va_list args;
	va_start(args, format);
	s = vstringprintf(format, args);
	va_end(args);
	return s;
}



string trim_whitespace(const std::string &haystack){
	long right = haystack.length()-1;
	long left = 0;
	while(((haystack[right]==' ') || (haystack[right]=='\t') || (haystack[right]=='\n')) && (right>=0)){
		--right;
	}
	while(((haystack[left]==' ') || (haystack[left]=='\t') || (haystack[left]=='\n')) && (left<right)){
		++left;
	}
	return haystack.substr(left,right-left+1);
}

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
inline int sgn(const TYPE value){
	return (value<0 ? -1 : 1);
}


// calculate result = a*X + b*Y, for vectors X & Y and scalars a & b
template<class TYPE>
inline void linear_combination(const double a, const std::vector<TYPE> &X, const double b, const std::vector<TYPE> &Y, std::vector<TYPE> &result){
	result.resize(X.size());
	for(long i=0; i<X.size(); ++i) result[i] = a*X[i] + b*Y[i];
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


// searches for needle in a list of ascending values
template<class TYPE1, class TYPE2>
long find_in_ascending_list(const std::vector<TYPE1> &haystack, const TYPE2 needle, const long start){
	for(long n=start; n<haystack.size(); ++n){
		if(haystack[n]>needle) return -1;
		if(haystack[n]==needle) return n;
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
inline double get_array_nonzero_min(const ARRAY_TYPE &X){
	const long N = X.size();
	double minX = NAN_D;
	for(long n=0; n<N; ++n){
		if((X[n]!=0) && (isnan(minX) || (X[n]<minX))) minX = X[n];
	}
	return minX;
}


template<class ARRAY_TYPE>
inline double get_array_max(const ARRAY_TYPE &X){
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


template<class ARRAY_TYPE>
inline bool arrays_are_equal(const ARRAY_TYPE &A, const ARRAY_TYPE &B){
	if(A.size()!=B.size()) return false;
	for(long i=0; i<A.size(); ++i){
		if(A[i]!=B[i]) return false;
	}
	return true;
}


template<class TYPE>
inline TYPE vector_sum(const std::vector<TYPE> &values){
	TYPE S = 0;
	for(long i=0; i<values.size(); ++i) S += values[i];
	return S;
}


inline long vector_sum(const std::vector<char> &values){
	long S = 0;
	for(long i=0; i<values.size(); ++i) S += values[i];
	return S;
}


double smallest_nonzero_step(const std::vector<double> &times){
	double S = INFTY_D;
	for(long i=0; i<times.size(); ++i){
		if(times[i+1]>times[i]){
			S = min(S, times[i+1]-times[i]);
		}
	}
	return S;
}


inline double vector_mean(const std::vector<double> &values){
	double S = 0;
	for(long i=0; i<values.size(); ++i) S += values[i];
	return (S/values.size());
}


inline double vector_abs_mean(const std::vector<double> &values){
	double S = 0;
	for(long i=0; i<values.size(); ++i) S += abs(values[i]);
	return (S/values.size());
}



inline double vector_mean(const std::vector<double> &values, const long first, const long last){
	double S = 0;
	for(long i=first; i<=last; ++i) S += values[i];
	return (S/(last-first+1.0));
}


// calculate sum of a single row in a 2D matrix of size NR x NC
// matrix must be provided in row-major format, i.e. matrix[r*NC+c] is the entry in row r & column c
inline double row_sum(const std::vector<double> &matrix, const long NC, const long row){
	double S = 0;
	for(long c=0; c<NC; ++c){
		S += matrix[row*NC + c];
	}
	return S;
}


// make sure no entry is negative
void make_vector_positive(std::vector<double> &values){
	for(long i=0; i<values.size(); ++i) values[i] = max(0.0, values[i]);
}

// replace negative entries
void replace_negatives(std::vector<double> &values, const double replacement){
	for(long i=0; i<values.size(); ++i){
		if(values[i]<0) values[i] = replacement;
	}
}

// replace non-strictly positive entries
void replace_non_positives(std::vector<double> &values, const double replacement){
	for(long i=0; i<values.size(); ++i){
		if(values[i]<=0) values[i] = replacement;
	}
}


// make sure entries in a vector are within the specified limits [min_value:max_value]
void cap_values(const double 		min_value,
				const double 		max_value,
				std::vector<double> &values){ // (INPUT/OUTPUT) the vector to be modified in-situ
	for(long i=0; i<values.size(); ++i){
		values[i] = max(min_value, min(max_value, values[i]));
	}
}


// extract a specific row from a 2D matrix in row-major format
template<class TYPE>
void extract_row(const std::vector<TYPE> &matrix, const long NC, const long row, std::vector<TYPE> &extracted_row){
	extracted_row.resize(NC);
	for(long c=0; c<NC; ++c){
		extracted_row[c] = matrix[row*NC+c];
	}
}


template<class TYPE>
void flatten_matrix(const std::vector<std::vector<TYPE> > 	&matrix,		// (INPUT) 2D matrix, where matrix[r][c] is the element in row r and column c. The number of rows is assumed to be matrix.size(). The number of columns is determined assumed to be matrix[0].size().
					std::vector<TYPE> 						&flattened){	// (OUTPUT) Flattened 2D matrix in row-major format, where flattened[r*NC + c] is the element in row r and column c.
	const long NR = matrix.size();
	if(NR==0){
		flattened.resize(0);
		return;
	}
	const long NC = matrix[0].size();
	flattened.resize(NR*NC);
	for(long r=0; r<NR; ++r){
		for(long c=0; c<NC; ++c){
			flattened[r*NC + c] = matrix[r][c];
		}
	}		
}


template<class TYPE>
inline bool contains_nan(const std::vector<TYPE> &values){
	for(long i=0; i<values.size(); ++i){
		if(std::isnan(values[i])) return true;
	}
	return false;
}


template<class TYPE>
inline bool contains_inf(const std::vector<TYPE> &values){
	for(long i=0; i<values.size(); ++i){
		if(std::isinf(values[i])) return true;
	}
	return false;
}

inline long vector_count_zeros(const std::vector<long> &values){
	long S = 0;
	for(long i=0; i<values.size(); ++i) S += (values[i]==0 ? 1 : 0);
	return S;
}


template<class ARRAY_TYPE>
inline long count_values_below_threshold(const ARRAY_TYPE &values, double threshold){
	long N = 0;
	for(long i=0; i<values.size(); ++i){
		if(values[i]<=threshold) ++N;
	}
	return N;
}


// remove an item from a vector, by replacing it with the last item in the list and then removing the last item
// if the order of items in the vector do not matter, then this is more efficient than removing an item from within a vector
// index is assumed to be a valid location in the vector
template<class TYPE>
inline void remove_item_from_vector(std::vector<TYPE> &list, long index){
	if(index==list.size()-1){
		// the item is the last one in the list, or nothing else left in the list, so just remove item
		list.pop_back();
	}else{
		list[index] = list.back();
		list.pop_back();
	}
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
inline void aux_qsortIndices_swap(std::vector<TYPE> &X, long A, long B, TYPE &temp){
	temp = X[A];
	X[A] = X[B];
	X[B] = temp;
}


template<typename TYPE>
long aux_qsortIndices_partition(const std::vector<TYPE> &values, std::vector<long> &indices, long start, long end){
	long temp;
	
	// choose pivot at middle
	const long pi = start+(end-start)/2;
	
	// Alternative: choose pivot via median-of-three rule
	/*
	long pi = start+(end-start)/2;
	if(values[indices[end]] < values[indices[start]]){
		aux_qsortIndices_swap(indices, start, end, temp);
	}
	if(values[indices[pi]] < values[indices[start]]){
		aux_qsortIndices_swap(indices, start, pi, temp);
	}
	if(values[indices[end]] < values[indices[pi]]){
		aux_qsortIndices_swap(indices, end, pi, temp);
	}
	*/
	const TYPE pv = values[indices[pi]];
	
	// swap (put pivot temporarily at end)
	aux_qsortIndices_swap(indices, pi, end, temp);
	
	long si = start;
	for(long i=start; i<end; ++i){
		if((values[indices[i]] < pv) || ((values[indices[i]]==pv) && (i%2==0))){ // modified condition of incrementing partition point, to avoid worst-case-scenarios when all (or many) values are equal
			// swap
			aux_qsortIndices_swap(indices, i, si, temp);
			++si;
		}
	}
	// swap (place pivot onto si)
	aux_qsortIndices_swap(indices, si, end, temp);
	
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




// inverse cumulative distribution function of the Student's t distribution with n degrees of freedom
// Returns t such that P(x<=t) = p
// Approximation according to:
//    Abramowitz and Stegun (1970). Handbook of mathematical functions. Page 949
//    Voutier (2010). A New Approximation to the Normal Distribution Quantile Function. arXiv:1002.0567
double quantile_Students_t(double p, long n){
	if(p<0.5) return -quantile_Students_t(1-p,n); // the approximation formula below is only valid for p>0.5 (i.e. q<0.5)
	double q = 1-p;
	const double A = sqrt(-2*log(q));
	//const double xq = A - (2.515517 + 0.802853*A + 0.010328*SQ(A))/(1 + 1.432788*A + 0.189269*SQ(A) + 0.001308*Qube(A));// [Abramowitz and Stegun (1970). Page 933, formula 26.2.23]
	const double xq = A - (2.653962002601684482 + 1.561533700212080345*A + 0.061146735765196993*SQ(A))/(1 + 1.904875182836498708*A + 0.454055536444233510*SQ(A) + 0.009547745327068945*Qube(A));// [Voutier (2010). A New Approximation to the Normal Distribution Quantile Function. arXiv:1002.0567. Pages 5-6]
	const double g1 = 0.25 * (Qube(xq) + xq);
	const double g2 = (1.0/96) * (5*pow(xq,5) + 16*Qube(xq) + 3*xq);
	const double g3 = (1.0/384) * (3*pow(xq,7) + 19*pow(xq,5) + 17*pow(xq,3) - 15*xq);
	const double g4 = (1.0/92160) * (79*pow(xq,9) + 776*pow(xq,7) + 1482*pow(xq,5) - 1920*Qube(xq) - 945*xq);
	double tq = xq + g1/n + g2/SQ(n) + g3/Qube(n) + g4/QTR(n); // [Abramowitz and Stegun (1970). Page 949]
	return tq;
}





#pragma mark -
#pragma mark Vectorized basic arithmetics
#pragma mark


template<class TYPE>
inline vector<TYPE> operator*(vector<TYPE> x, const vector<TYPE> &y){
	for(long i=0; i<x.size(); ++i){
		x[i] *= y[i];
	}
	return x;
}


template<class TYPE>
inline vector<TYPE>& operator*=(vector<TYPE> &x, const vector<TYPE> &y){
	for(long i=0; i<x.size(); ++i){
		x[i] *= y[i];
	}
	return x;
}


template<class TYPE>
inline vector<TYPE> operator*(vector<TYPE> x, double scalar){
	for(long i=0; i<x.size(); ++i){
		x[i] *= scalar;
	}
	return x;
}


template<class TYPE>
inline vector<TYPE> operator*(double scalar, vector<TYPE> x){
	for(long i=0; i<x.size(); ++i){
		x[i] *= scalar;
	}
	return x;
}


template<class TYPE>
vector<TYPE>& operator*=(vector<TYPE> &x, double scalar) {
	for(long i=0; i<x.size(); ++i){
		x[i] *= scalar;
	}
	return x;
}


template<class TYPE>
inline vector<TYPE> operator+(vector<TYPE> x, const vector<TYPE> &y){
	for(long i=0; i<x.size(); ++i){
		x[i] += y[i];
	}
	return x;
}


template<class TYPE>
inline vector<TYPE> &operator+=(vector<TYPE> &x, const vector<TYPE> &y){
	for(long i=0; i<x.size(); ++i){
		x[i] += y[i];
	}
	return x;
}


template<class TYPE>
inline vector<TYPE> operator-(vector<TYPE> x, const vector<TYPE> &y){
	for(long i=0; i<x.size(); ++i){
		x[i] -= y[i];
	}
	return x;
}


template<class TYPE>
inline vector<TYPE> &operator-=(vector<TYPE> &x, const vector<TYPE> &y){
	for(long i=0; i<x.size(); ++i){
		x[i] -= y[i];
	}
	return x;
}




template<class TYPE>
inline vector<TYPE> operator/(vector<TYPE> x, const vector<TYPE> &y){
	for(long i=0; i<x.size(); ++i){
		x[i] /= y[i];
	}
	return x;
}



template<class TYPE>
inline vector<TYPE> operator/(vector<TYPE> x, double scalar){
	for(long i=0; i<x.size(); ++i){
		x[i] /= scalar;
	}
	return x;
}


template<class TYPE>
inline vector<TYPE> &operator/=(vector<TYPE> &x, double scalar){
	for(long i=0; i<x.size(); ++i){
		x[i] /= scalar;
	}
	return x;
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
	return sqrt(-2.0*log(uniformWithinInclusiveRight(0, 1)))*cos(2.0*M_PI*uniformWithinInclusiveRight(0,1));
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


// pick an index within 0:(N-1) at probability proportional to weights[i]
// for efficiency, the caller guarantees that total_weight = sum_i weights[i]
long random_int_from_distribution(const std::vector<double> &weights, const double total_weight){
	const long N = weights.size();
	double p = R::runif(0.0,1.0);
	for(long i=0; i<N; ++i){
		if(p<=weights[i]/total_weight) return i;
		p -= weights[i]/total_weight;
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


// pick an index within 0:(N-1) at probability proportional to weights[weight_indices[index_pool[i]]]
// for efficiency, the caller guarantees that total_weight = sum_i weights[index_pool[i]]
long random_int_from_distribution(const std::vector<long> &index_pool, const std::vector<double> &weights, const std::vector<long> &weight_indices, const double total_weight){
	const long N = index_pool.size();
	double probability = R::runif(0.0,1.0);
	double dp;
	for(long i=0; i<N; ++i){
		dp = weights[weight_indices[index_pool[i]]]/total_weight;
		if(probability<=dp) return i;
		probability -= dp;
	}
	return N-1;
}



// given a list of index pools pool_1, pool_2, ..., pool_NP, and a probability weight associated with each pool (e.g. each item in pool_i has weight weight_i), pick a random item among all pools
// for efficiency, the caller guarantees that total_weight = 1, where total_weight := sum_p weights[p] * pools[p].size()
// The time complexity of this function is O(NP)
// if the number items in each pool is much larger than the number of distinct pools, this function is much more efficient than considering all items within a single pool
// this function guarantees that only non-empty pools are picked, or returns false if all pools are actually empty
// returns true upon success (guaranteed, provided that at least one pooll is non-empty)
bool random_int_from_pools(	const std::vector<lvector> 	&pools, 		// (INPUT) array of size NP, each element of which is a pool of integers (indices)
							const std::vector<double> 	&weights, 		// (INPUT) array of size NP, listing probability weights associated with each pool. Hence, weights[p] is the weight assigned to each element in index_pools[p], and hence weights[p]*pools[p].size() is the probability weight of landing in pool p.
							const double 				total_weight,	// (INPUT) the total weight across all items in all pools. Provided by the caller, for efficiency
							long						&p,				// (OUTPUT) the random pool chosen, i.e. a random integer between 0,..,NP
							long						&i){			// (OUTPUT) the random item in the chosen pool, i.e. a random integer between 0,...,pools[p].size()-1
	const long NP = pools.size();
	// step 1: choose a random pool based on their total probabilities
	double probability = R::runif(0.0,1.0);
	double dp;
	p = 0;
	long last_valid_p = -1; // keep track of the last non-empty pool
	while(p<NP){
		if(!pools[p].empty()) last_valid_p = p;
		dp = weights[p]*pools[p].size()/total_weight;
		if((probability<=dp) && (!pools[p].empty())) break; // pick this pool
		probability -= dp;
		++p; // try the next pool
	}
	if(last_valid_p<0){
		// all pools were empty, so return failure code
		p = -1; i = -1;
		return false;
	}
	if(p>=NP) p = last_valid_p; // probability overflow, likely due to numerical rounding errors. So pick the last non-empty pool

	// step 2: choose a random item in the chosen pool
	// this step is efficient, because all items in this pool have the same probability weight
	i = uniformIntWithin(0,pools[p].size()-1);
	return true;
}





// generate exponentially distributed random variable, with PDF f(x) = lambda*exp(-lambda*x)
double random_exponential_distribution(double lambda){
	return -log(R::runif(0.0,1.0))/lambda;
}


inline bool random_bernoulli(double p){
	//return ((double(rand())/RAND_MAX) <= p); // rand() is discouraged by R package builders
	return (R::runif(0.0,1.0)<=p);
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


template<class ARRAY_TYPE>
void print_as_matrix_column_major(const long NR, const long NC, const ARRAY_TYPE &A){
	for(long r=0; r<NR; ++r){
		for(long c=0; c<NC; ++c){
			Rcout << (c>0 ? ", " : "") << A[r + c*NR];
		}
		Rcout << "\n";
	}
	Rcout << "\n";
}


// ##########################################################
// MATRIX ALGEBRA

#pragma mark -
#pragma mark Matrix algebra
#pragma mark 


void get_identity_matrix(	const long 			NR,
							std::vector<double>	&A){	// (INPUT/OUTPUT) will be the identity matrix of size NR x NR, in row-major format
	A.assign(NR*NR,0);
	for(long r=0; r<NR; ++r){
		A[r*NR + r] = 1;
	}
}


// calculate the dot-product between two vectors, as needed for QR decomposition
double QR_dot_product(	const long 		N, 			 // (INPUT) the number of entries in the vectors
						const double	X[], 		 // (INPUT) vector of size N
						const long 		xincrement,  // (INPUT) the increment between successive entries in X
						const double 	Y[], 		 // (INPUT) vector of size N
						const long 		yincrement){ // (INPUT) the increment between successive entries in Y
	long i, ix, iy, m;
	if (N<=0) return 0;

	double S = 0.0;
	if(xincrement != 1 || yincrement != 1){
		if(0 <= xincrement){
			ix = 0;
		}else{
			ix = (-N + 1) * xincrement;
		}

		if(0 <= yincrement){
			iy = 0;
		}else{
			iy = (-N + 1) * yincrement;
		}

		for(i = 0; i<N; i++){
			S = S + X[ix] * Y[iy];
			ix = ix + xincrement;
			iy = iy + yincrement;
		}
	}else{
		m = N % 5;
		for (i = 0; i<m; i++){
			S = S + X[i] * Y[i];
		}
		for(i = m; i<N; i = i + 5){
			S = S + X[i] * Y[i] + X[i+1] * Y[i+1] + X[i+2] * Y[i+2] + X[i+3] * Y[i+3] + X[i+4] * Y[i+4];
		}
	}
	return S;
}


// swap two vectors in-situ, as needed for QR decomposition
void QR_swap_vectors(	const long 	N, 				// (INPUT) the length of the vectors
						double 		x[],			// (INPUT/OUTPUT) one of the vectors to swap
						const long 	xincrement, 	// (INPUT) increment between successive elements in Y
						double 		y[], 			// (INPUT/OUTPUT) one of the vectors to swap
						const long 	yincrement){	// (INPUT) increment between successive elements in Y
	long i, ix, iy, m;
	double temp;
	if(N<=0) return;
	
	if((xincrement == 1) && (yincrement == 1)){
		m = N % 3;
		for(i = 0; i<m; i++){
			temp = x[i];
			x[i] = y[i];
			y[i] = temp;
		}

		for(i = m; i<N; i = i + 3){
			temp 	= x[i];
			x[i] 	= y[i];
			y[i] 	= temp;

			temp 	= x[i+1];
			x[i+1] 	= y[i+1];
			y[i+1] 	= temp;

			temp 	= x[i+2];
			x[i+2] 	= y[i+2];
			y[i+2] 	= temp;
		}
	}else{
		if(0 <= xincrement){
			ix = 0;
		}else{
			ix = (-N + 1) * xincrement;
		}
		if(0 <= yincrement){
			iy = 0;
		}else{
			iy = (-N + 1) * yincrement;
		}
		for(i = 0; i<N; i++){
			temp 	= x[ix];
			x[ix] 	= y[iy];
			y[iy] 	= temp;
			ix 		= ix + xincrement;
			iy 		= iy + yincrement;
		}
	}
}


// calculate euclidean norm of vector X
double euclidean_norm(	const long 		N,				// (INPUT) the number of elements in the vector
						const double	X[], 			// (INPUT) 1D vector of size N
						const long 		xincrement){	// (INPUT) the increment between successive entries of X
	long i, ix;
	double absxi, norm, scale, sum_squares;

	if((N<1) || (xincrement<1)){
		norm = 0;
	}else if(N== 1){
		norm = abs(X[0]);
	}else{
		scale 	= 0.0;
		sum_squares 	= 1.0;
		ix 		= 0;
		for(i = 0; i < N; i++){
			if(X[ix] != 0.0){
				absxi = abs(X[ix]);
				if(scale < absxi){
					sum_squares = 1.0 + sum_squares * (scale/absxi) * (scale/absxi);
					scale = absxi;
				}else{
					sum_squares = sum_squares + (absxi/scale) * (absxi/scale);
				}
			}
			ix = ix + xincrement;
		}
		norm  = scale * sqrt(sum_squares);
	}
	return norm;
}



void QR_scale_vector(	const long 	N, 				// (INPUT) the number of elements in the vector
						double 		scaling_factor, // (INPUT) scaling factor
						double 		X[], 			// (INPUT/OUTPUT) 1D vector of size N
						const long	xincrement){			// (INPUT) the increment between successive entries of X.
	long i, ix, m;
	if(N<=0) return;

	if(xincrement == 1){
		m = N % 5;
		for(i = 0; i < m; i++){
			X[i] = scaling_factor * X[i];
		}
		for(i = m; i<N; i = i + 5){
			X[i]   = scaling_factor * X[i];
			X[i+1] = scaling_factor * X[i+1];
			X[i+2] = scaling_factor * X[i+2];
			X[i+3] = scaling_factor * X[i+3];
			X[i+4] = scaling_factor * X[i+4];
		}
	}else{
		if(0 <= xincrement){
			ix = 0;
		}else{
			ix = (-N + 1) * xincrement;
		}
		for(i = 0; i < N; i++){
			X[ix] = scaling_factor * X[ix];
			ix = ix + xincrement;
		}
	}
	return;
}



// add two vectors, as used for QR decomposition
void QR_add_vectors(const long 		N, 				// (INPUT) the number of elements in the vectors
					const double 	alpha, 			// (INPUT) an optional multiplicative factor for X
					const double	X[], 			// (INPUT) 1D vector of size N
					const long 		xincrement, 	// (INPUT) the increment between successive entries of X
					double 			Y[], 			// (INPUT/OUTPUT) 1D vector of size N. Upon return, this will store the result alpha*X+Y
					const long		yincrement){	// (INPUT) the increment between successive entries of Y
  long i, ix, iy, m;
  if(N<=0) return;
  if(alpha==0) return;

	if((xincrement != 1) || (yincrement != 1)){
		if(0 <= xincrement){
			ix = 0;
		}else{
			ix = (-N + 1) * xincrement;
		}

		if(0 <= yincrement){
			iy = 0;
		}else{
			iy = (-N + 1) * yincrement;
		}

		for(i = 0; i<N; i++){
			Y[iy] = Y[iy] + alpha * X[ix];
			ix = ix + xincrement;
			iy = iy + yincrement;
		}
	}else{
		m = N % 4;

		for(i = 0; i<m; i++){
			Y[i] = Y[i] + alpha * X[i];
		}

		for(i = m; i<N; i = i + 4){
			Y[i]   = Y[i]   + alpha * X[i];
			Y[i+1] = Y[i+1] + alpha * X[i+1];
			Y[i+2] = Y[i+2] + alpha * X[i+2];
			Y[i+3] = Y[i+3] + alpha * X[i+3];
		}
	}
	return;
}


// QR decomposition of a real 2D matrix of size NR x NC
// code adjusted from: https://people.sc.fsu.edu/~jburkardt/c_src/qr_solve/qr_solve.html
void QR_decomposition (	double 			A[], 			// (INPUT/OUTPUT) on input, the matrix to be decomposed, in column-major format. On output, A contains in its upper triangle the upper triangular matrix R of the QR factorization.  Below its diagonal A contains information from which the orthogonal part of the decomposition can be recovered.
						const long 		LDA, 			// (INPUT) the leading dimension of A, i.e. linear index periodicity between subsequent columns. Typically this will be equal to NR, but may also be larger (e.g. when processing a sub-matrix)
						const long 		NR,				// (INPUT) the number of rows of A
						const long 		NC,				// (INPUT) the number of columns of A
						const bool 		pivoting,		// (INPUT) whether to perform pivoting
						double 			scratch[], 		// (SCRATCH) pre-allocated scratch space of size >=NC. Only needed if pivoting=true.
						double 			QRaux[], 		// (OUTPUT) 1D array of size NC, must be pre-allocated. Will contain information required to recover the orthogonal part of the decomposition.
						long 			pivots[],		// (OUTPUT) 1D array of size NC, must be pre-allocated. pivots[k] specifies the index of the column of A that has been interchanged into the K-th column, if pivoting was requested.
						long			&rank){			// (OUTPUT) the estimated rank of the matrix A
	long j, jp, l, lup, maxj;
	double maxnrm, nrmxl;
	long pl, pu, swapj;
	double t, tt;

	// initialize pivots
	if(pivoting){
		for(long j = 0; j<NC; j++) pivots[j] = 0;
	}

	pl = 1;
	pu = 0;
	if(pivoting){	
		for(j = 1; j <= NC; j++){
			swapj = (0 < pivots[j-1]);
			if(pivots[j-1] < 0){
				pivots[j-1] = -j;
			}else{
				pivots[j-1] = j;
			}

			if(swapj){
				if(j != pl) QR_swap_vectors(NR, A+0+(pl-1)*LDA, 1, A+0+(j-1), 1);
				pivots[j-1] = pivots[pl-1];
				pivots[pl-1] = j;
				pl = pl + 1;
			}
		}
		pu = NC;
		
		for(j = NC; 1<=j; j--){
			if(pivots[j-1] < 0){
				pivots[j-1] = -pivots[j-1];
				if(j != pu){
					QR_swap_vectors(NR, A+0+(pu-1)*LDA, 1, A+0+(j-1)*LDA, 1);
					jp = pivots[pu-1];
					pivots[pu-1] = pivots[j-1];
					pivots[j-1] = jp;
				}
				pu = pu - 1;
			}
		}
	}

	// Compute column norms
	for(j = pl; j<=pu; j++){
		QRaux[j-1] = euclidean_norm(NR, A+0+(j-1)*LDA, 1);
	}

	for(j = pl; j<=pu; j++){
		scratch[j-1] = QRaux[j-1];
	}

	// Householder reduction of A
	lup = min(NR, NC);
	for(l = 1; l<=lup; l++){
		if(pl<=l && l<pu){
			maxnrm = 0.0;
			maxj = l;
			for(j = l; j<=pu; j++){
				if(maxnrm < QRaux[j-1]){
					maxnrm 	= QRaux[j-1];
					maxj 	= j;
				}
			}

			if(maxj != l){
				QR_swap_vectors(NR, A+0+(l-1)*LDA, 1, A+0+(maxj-1)*LDA, 1);
				QRaux[maxj-1] 	= QRaux[l-1];
				scratch[maxj-1] = scratch[l-1];
				jp 				= pivots[maxj-1];
				pivots[maxj-1] 	= pivots[l-1];
				pivots[l-1] 	= jp;
			}
		}

		// Householder transformation for column L.
		QRaux[l-1] = 0.0;
		if(l != NR){
			nrmxl = euclidean_norm(NR-l+1, A+l-1+(l-1)*LDA, 1);
			if(nrmxl != 0.0){
				if(A[l-1+(l-1)*LDA] != 0.0){
					nrmxl = nrmxl * (A[l-1+(l-1)*LDA]<0 ? -1 : +1);
				}
				QR_scale_vector(NR-l+1, 1.0 / nrmxl, A+l-1+(l-1)*LDA, 1);
				A[l-1+(l-1)*LDA] = 1.0 + A[l-1+(l-1)*LDA];
				for(j = l + 1; j <= NC; j++){
					t = -QR_dot_product(NR-l+1, A+l-1+(l-1)*LDA, 1, A+l-1+(j-1)*LDA, 1)/A[l-1+(l-1)*LDA];
					QR_add_vectors(NR-l+1, t, A+l-1+(l-1)*LDA, 1, A+l-1+(j-1)*LDA, 1);

					if(pl <= j && j <= pu){
						if(QRaux[j-1] != 0.0){
							tt = 1.0 - pow(abs(A[l-1+(j-1)*LDA])/QRaux[j-1], 2);
							tt = max(tt, 0.0);
							t = tt;
							tt = 1.0 + 0.05 * tt * pow(QRaux[j-1]/scratch[j-1], 2);

							if(tt != 1.0){
								QRaux[j-1] = QRaux[j-1] * sqrt(t);
							}else{
								QRaux[j-1] = euclidean_norm(NR-l, A+l+(j-1)*LDA, 1);
								scratch[j-1] = QRaux[j-1];
							}
						}
					}
				}
				// save transformation
				QRaux[l-1] = A[l-1+(l-1)*LDA];
				A[l-1+(l-1)*LDA] = -nrmxl;
			}
		}
	}
	
	// compute the rank
	const double tol = 1e-6;
	rank = 0;
  	long k = min(NR,NC);
	for(j = 0; j < k; j++){
		if(abs(A[j+j*LDA]) <= tol * abs(A[0+0*LDA])){
			break;
		}
		rank = j + 1;
	}	
}



// perform various operations using the QR-decomposition of some real 2D matrix A
// this function does not actually use the matrix A, but instead its pre-computed QR-decomposition
// code adjusted from: https://people.sc.fsu.edu/~jburkardt/c_src/qr_solve/qr_solve.html
long QR_operation(	double			QRA[], 		// (INPUT/OUTPUT) matrix of size LDA*NC, storing the QR-decomposition of A. Not actually modified permanently by this function, but used temporarily internally. So its output is the same as input.
					const long		LDA, 		// (INPUT) the leading dimension of the matrix A
					const long		NR, 		// (INPUT) number of rows in A
					const long 		NK, 		// (INPUT) number of columns in the AK matrix, formed during QR decomposition. NK will always be <=min(NR,NC)
					const double	QRaux[], 	// (INPUT) 1D array of size NC, containing information regarding the QR decomposition. Must be as calculated using QR_decomposition(..).
					const double	Y[], 		// (INPUT) 1D vector of size NR, to be operated upon by the QR-decomposed A
					double 			QY[], 		// (OUTPUT) the product Q*Y, if requested. In that case, must be preallocated of size NR.
					double 			QTY[], 		// (OUTPUT) the product Q^T*Y, if requested. In that case, must be preallocated of size NR.
					double 			X[], 		// (OUTPUT) the solution to the least squares problem: minimize ||AK*X - Y||_2, if requested. In that case, it must be preallocated of size NK
					double 			residuals[], 	// (OUTPUT) the least-squares residuals Y - AK*X, if requested. In that case, must be preallocated of size NR.
					double 			AX[], 			// (OUTPUT) the least-squares approximation AX:=AK*X, if requested. In that case, it must be preallocated of size NR.
					const string 	&job){			// (INPUT) string of size 5 and of format ABCDE, where each character is either 1 or 0, specifying the various jobs to perform. For example, "00110" performs jobs C and D.
	long i, info, j, jj, ju;
	double t, temp;
	info = 0;

	// determine jobs
	bool cQY 	= (job[0]!='0'); // compute QY
	bool cQTY	= (job[1]!='0'); // compute QTY
	bool cx 	= (job[2]!='0'); // compute QTY and X
	bool cr 	= (job[3]!='0'); // compute QTY and RSD
	bool cax 	= (job[4]!='0'); // compute QTY and AX
	
	cQTY = cQTY || cx || cr || cax; // always compute QTY if one of cx, cr, cax is requested

	ju = min(NK, NR-1);
	if(ju == 0){
		if(cQY) QY[0] = Y[0];
		if(cQTY) QTY[0] = Y[0];
		if(cax) AX[0] = Y[0];
		if(cx){
			if(QRA[0+0*LDA] == 0.0){
				info = 1;
			}else{
				X[0] = Y[0] / QRA[0+0*LDA];
			}
		}
		if(cr) residuals[0] = 0.0;
		return info;
	}

	// prepare computation of QY or QTY
	if(cQY){
		for(i = 1; i <= NR; i++){
			QY[i-1] = Y[i-1];
		}
	}
	if(cQTY){
		for(i = 1; i <= NR; i++){
			QTY[i-1] = Y[i-1];
		}
	}
	
	// compute QY
	if(cQY){
		for(jj = 1; jj <= ju; jj++){
			j = ju - jj + 1;
			if(QRaux[j-1] != 0.0){
				temp = QRA[j-1+(j-1)*LDA];
				QRA[j-1+(j-1)*LDA] = QRaux[j-1];
				t = -QR_dot_product(NR-j+1, QRA+j-1+(j-1)*LDA, 1, QY+j-1, 1 ) / QRA[j-1+(j-1)*LDA];
				QR_add_vectors(NR-j+1, t, QRA+j-1+(j-1)*LDA, 1, QY+j-1, 1 );
				QRA[j-1+(j-1)*LDA] = temp;
			}
		}
	}
	
	// compute Q'*Y.
	if(cQTY){
		for(j = 1; j <= ju; j++){
			if(QRaux[j-1] != 0.0){
				temp = QRA[j-1+(j-1)*LDA];
				QRA[j-1+(j-1)*LDA] = QRaux[j-1];
				t = -QR_dot_product(NR-j+1, QRA+j-1+(j-1)*LDA, 1, QTY+j-1, 1 ) / QRA[j-1+(j-1)*LDA];
				QR_add_vectors(NR-j+1, t, QRA+j-1+(j-1)*LDA, 1, QTY+j-1, 1 );
				QRA[j-1+(j-1)*LDA] = temp;
			}
		}
	}

	// prepare computation of X, RSD, or AX.
	if(cx){
		for(i = 1; i <= NK; i++){
			X[i-1] = QTY[i-1];
		}
	}
	if(cax){
		for(i = 1; i <= NK; i++){
			AX[i-1] = QTY[i-1];
		}
	}
	if(cr && NK < NR){
		for(i = NK+1; i <= NR; i++){
			residuals[i-1] = QTY[i-1];
		}
	}
	if(cax && NK+1 <= NR){
		for(i = NK+1; i <= NR; i++){
			AX[i-1] = 0.0;
		}
	}
	if(cr){
		for(i = 1; i <= NK; i++){
			residuals[i-1] = 0.0;
		}
	}

	// compute X
	if(cx){
		for(jj = 1; jj <= NK; jj++){
			j = NK - jj + 1;
			if(QRA[j-1+(j-1)*LDA] == 0.0){
				info = j;
				break;
			}
			X[j-1] = X[j-1] / QRA[j-1+(j-1)*LDA];
			if(j != 1){
				t = -X[j-1];
				QR_add_vectors(j-1, t, QRA+0+(j-1)*LDA, 1, X, 1);
			}
		}
	}

	// compute residuals and AX, if needed
	if(cr || cax){
		for(jj = 1; jj <= ju; jj++){
			j = ju - jj + 1;
			if(QRaux[j-1] != 0.0){
				temp = QRA[j-1+(j-1)*LDA];
				QRA[j-1+(j-1)*LDA] = QRaux[j-1];
				if(cr){
					t = -QR_dot_product(NR-j+1, QRA+j-1+(j-1)*LDA, 1, residuals+j-1, 1) / QRA[j-1+(j-1)*LDA];
					QR_add_vectors(NR-j+1, t, QRA+j-1+(j-1)*LDA, 1, residuals+j-1, 1);
				}
				if(cax){
					t = -QR_dot_product(NR-j+1, QRA+j-1+(j-1)*LDA, 1, AX+j-1, 1) / QRA[j-1+(j-1)*LDA];
					QR_add_vectors(NR-j+1, t, QRA+j-1+(j-1)*LDA, 1, AX+j-1, 1 );
				}
				QRA[j-1+(j-1)*LDA] = temp;
			}
		}
	}

	return info;
}



// solve a linear system of equations of the format:
//	  A*X = B
// in a least squares sence, i.e. minimize ||A*X - B||_2
// Here, A is a real matrix of size NR x NC, and X is a 1D vector of size NC.
// Uses QR-decomposition, which must be performed beforehand using the function QR_decomposition(..)
// Hence, this function does not actually use the matrix A, but its QR decomposition QRA (and some other auxiliary variables).
// Code adjusted from: https://people.sc.fsu.edu/~jburkardt/c_src/qr_solve/qr_solve.html
void QR_linear_least_squares( 	double			QRA[], 			// (INPUT/OUTPUT) array of size LDA*NC, containing the QR factorization, as computed using QR_decomposition(..).  Not actually modified permanently by this function, but used temporarily internally. So its output is the same as input.
								const long 		LDA,			// (INPUT) the leading dimension of A, i.e. linear index periodicity between subsequent columns. Typically this will be equal to NR, but may also be larger (e.g. when processing a sub-matrix)
								const long 		NR, 			// (INPUT) the number of rows of A
								const long 		NC,				// (INPUT) the number of columns of A
								const long		rank, 			// (INPUT) the rank of the matrix, for example estimated via QR_decomposition
								const double	B[], 			// (INPUT) 1D vector of size NR
								long 			pivots[],		// (INPUT) 1D array of size NC, specifying the index of the column of A that has been interchanged into the K-th column, if pivoting was requested during QR. Must be as calculated using QR_decomposition(..).
								double 			QRaux[],		// (INPUT) 1D array of size NC, containing information regarding the QR decomposition. Must be as calculated using QR_decomposition(..).
								double 			X[], 			// (OUTPUT) 1D vector of size NC, containing a least-squares solution. Must be preallocated.
								double 			residuals[]){	// (OUTPUT) residuals, B - A*X
	long i, info, j, k;
	double t;
	double *dummyAB=NULL, *dummyQY=NULL;

	if(rank != 0){
		info = QR_operation(QRA, LDA, NR, rank, QRaux, B, dummyQY, residuals, X, residuals, dummyAB, "00110");
	}
	for(i = 0; i<NC; i++){
		pivots[i] = - pivots[i];
	}
	for(i = rank; i<NC; i++){
		X[i] = 0.0;
	}
	for(j = 1; j<=NC; j++){
		if(pivots[j-1] <= 0){
			k = - pivots[j-1];
			pivots[j-1] = k;

			while(k != j){
				t = X[j-1];
				X[j-1] = X[k-1];
				X[k-1] = t;
				pivots[k-1] = -pivots[k-1];
				k = pivots[k-1];
			}
		}
	}
	return;
}



// Solver the linear least squares problem:
//	minimize ||A*X-B||_2
// for some matrix A and some vector or matrix X, using QR-decomposition of A
// If B (and thus X) has NCb columns, this corresponds to NCb independent least-squares problems, which are solved separately (but profiting from the same QR decomposition)
// Most of the computation time goes into the QR-decomposition of the matrix A
void QR_linear_least_squares(	const long		NRa,		// (INPUT) number of rows in A
								const long 		NCa,		// (INPUT) number of columns in A, and number of rows in X
								const long		NCb,		// (INPUT) number of columns in B. This is the number of independent problems to solve.
								const std::vector<double> 	&A,			// (INPUT) 2D matrix of size NRa x NCa, in row-major or column-major format
								const std::vector<double> 	&B,			// (INPUT) 1D column vector of length NRa = NRb, or a 2D matrix of size NRa x NCb, in row-major or column-major format
								const bool					row_major,	// (INPUT) indicating whether A and B are stored in row-major format (instead of column-major). The same applies to the output X.
								std::vector<double>			&QRA,		// (SCRATCH) scratch space for internal computations and for storing the QR-decomposition of A. Will be resized as needed up to size NRa x NCa
								std::vector<double>			&scratch,	// (SCRATCH) scratch space for internal computations. Will be resized as needed up to size max(NCa,NRa)
								std::vector<double>			&X,			// (OUTPUT) vector or matrix of size NCa x NCb, storing the solutitions to the least-squares problems. Will be in row-major format if row_major=true, otherwise it will be in column-major format.
								long						&rank){		// (OUTPUT) the estimated rank of the matrix A
	const long NRx = NCa;
	const long NCx = NCb;
	long r,c;

	// prepare and work on a copy of A, in column-major format
	const long LDA = NRa; // the leading dimension of A, i.e. linear index periodicity between subsequent columns. Here we assume that A is stored compactly in memory, i.e. LDA=NRa (after A is turned into column-major format, if needed).
	std::vector<double> BCM(NRa*NCb); // store a copy of B in column-major format
	if(row_major){
		// transform A into column-major format and store result in QRA
		QRA.resize(NRa*NCa);
		for(r=0; r<NRa; ++r){
			for(c=0; c<NCa; ++c){
				QRA[c*NRa + r] = A[r*NCa + c];
			}
		}
		// transform B into column-major format and store result in BCM
		for(r=0; r<NRa; ++r){
			for(c=0; c<NCb; ++c){
				BCM[c*NRa + r] = B[r*NCb + c];
			}
		}
	}else{
		// A and B are already in column-major format
		QRA = A;
		BCM = B;
	}
						
	// compute QR decomposition of A, and store results directly in QRA
	scratch.resize(NCa);
	std::vector<double> QRaux(NCa);
	std::vector<long> pivots(NCa);
	QR_decomposition (	&QRA[0],
						LDA,
						NRa,
						NCa,
						true, // pivoting
						&scratch[0],
						&QRaux[0],
						&pivots[0],
						rank);
	
	// use QR decomposition to minimize ||A*X - B||_2, separately for each column in X
	// QRA, BCM, XCM are all treated in column-major format at this point
	std::vector<double> XCM(NCa*NCb);
	scratch.resize(NRa); // use to store residuals
	for(long j=0; j<NCb; ++j){
		QR_linear_least_squares(&QRA[0], LDA, NRa, NCa, rank, &BCM[j*NRa], &pivots[0], &QRaux[0], &XCM[j*NCa], &scratch[0]);
	}
		
	// transform XCM to row-major format and store in X, if needed
	if(row_major){
		X.resize(NRx*NCx);
		for(r=0; r<NRx; ++r){
			for(c=0; c<NCx; ++c){
				X[r*NCx + c] = XCM[c*NRx + r];
			}
		}
	}else{
		X = XCM;
	}
}



// calculate the inverse of a square matrix A, using QR decomposition
// A can be in row-major or column-major format; the internal computation is the same, because confusing majority is equivalernt to temporary transposing
void QR_matrix_inverse(	const long					N,			// (INPUT) number of rows & columns in A
						const std::vector<double> 	&A,			// (INPUT) 2D matrix of size N x N, in row-major or column-major format (which, does not matter)
						std::vector<double>			&QRA,		// (SCRATCH) scratch space for internal computations and for storing the QR-decomposition of A. Will be resized as needed up to size N^2
						std::vector<double>			&Ainv,		// (OUTPUT) 2D matrix of size N x N, storing the inverse of A (or an approximation thereof, if A is non-invertible). Ainv will be in row-major format iff A was in row-major format.
						long						&rank){		// (OUTPUT) an estimate of the rank of A. If A is invertible, this will be equal to N.
	dvector identity, scratch;
	get_identity_matrix(N,identity);
	// solve linear system A*Ainv = identity, in the least squares sense
	// note that we pretent as if A was in column-major format for efficiency (QR would otherwise transform everything temporarily)
	//   If A is actually in row-major format, then pretending as if it's column-major is equivalent to transposing A, 
	//   hence the obtained Ainv will just be the inverse of A^T in column-major format, or equivalently, the inverse of A in row-major format.
	QR_linear_least_squares(N,N,N,A,identity,false,QRA,scratch,Ainv,rank);
}




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



// Solve linear system Ax=b for vector x, using LU decomposition of the matrix A
// scratchSpace[] is only needed for the internal calculations, and should be allocated prior to calling (size at least N*N).
// maxError defines the maximum acceptable L2 norm E(X):=|AX-b|, where X is the approximate solution.
// As long as E(X)>maxError, the solution X is improved iteratively.
// If maxError<=0 or maxImprovements==0, no iterative improvement is performed, i.e. X is set to the solution obtained from regular LU decomposition.
// x[] will contain the solution vector upon return.
// x[] should be allocated prior to calling (size at least N).
// Returns false on error (e.g. when matrix is singular).
template<class TYPE>
bool LUsolveLinearSystem(	const TYPE 			matrix[], 			// (INPUT) in row major format
							TYPE				scratchSpace[],		// (SCRATCH) pre-allocated of size N*N
							unsigned long 		N,
							const TYPE 			b[],
							double 				maxError, 			//max allowed L2 error |Ax-b|
							unsigned int 		maxImprovements, 	//max allowed number of iterative improvements of solution. Use this to cap computing time if your matrices are unpredictably *pathological*
							TYPE 				x[]){				// (OUTPUT) the solution vector x. Must be pre-allocated of size N.
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


// Calculate the inverse of a rectangular non-singular matrix.
// Matrix needs to be in LU decomposition format (LUmatrix), as returned by the routine LUDecomposition(..) defined above.
// IPIV[] should store the pivoting indices from the LU decomposition, as returned by LUDecomposition(..)
// inverse[] (output) will contain the inverse matrix in row-major format, and should already be allocated before calling (size at least N*N).
template<class TYPE>
void LUInverse(	const TYPE			LUmatrix[], 
				unsigned long 		N, 
				const unsigned long IPIV[], 
				TYPE				inverse[]){
	TYPE *col = new TYPE[N];
	long i,j;
	for(i=0; i<N; ++i){
		for(j=0; j<N; ++j){ col[j] = 0; }
		col[i] = 1;
		LUSolveLinearSystem(LUmatrix, N, IPIV, col);
		for(j=0; j<N; ++j){ inverse[j*N+i] = col[j]; }
	}
	delete[] col;
}



// Calculate the inverse of a rectangular, non-singular matrix.
// Uses LU decomposition, forward substitution and back substitution.
// Merely a wrapper for the routines LUDecomposition() and LUInverse() defined above.
// Returns false on error (e.g. if matrix is singular).
// matrix[] (input) will be modified after this call (so back it up if you still need it).
// inverse[] (output) should already be allocated before calling (at least size N*N).
template<class TYPE>
bool inverseMatrix(	TYPE 			matrix[], 
					unsigned long 	N, 
					TYPE			inverse[]){
	if(N==0) return false;
	unsigned long *IPIV = new unsigned long[N];
	int Psign;
	if(!LUDecomposition(matrix, N, IPIV, Psign)){
		delete[] IPIV;
		return false; // error on LU decomposition
	}
	LUInverse(matrix, N, IPIV, inverse);
	delete[] IPIV;
	return true;
}



// same as above, but LU decomposition method is supplemented by iterative improvement of solution
template<class TYPE>
bool inverseMatrix(	const TYPE			matrix[], // (input)
					unsigned long 		N, 
					TYPE		 		inverse[], // (output) should be allocated to size at least N*N
					TYPE				scratchSpace[], // should be allocated to size at least N*N	
					double 				maxError, //max allowed L2 error |Ax - b|, where x is any column of the inverse and b is the corresponding unit vector (0,..,1,..,0)
					unsigned int 		maxImprovements){ //max allowed number of iterative improvements of inverse. Use this to cap computing time if your matrices are unpredictably *pathological*. Set to 0 for unlimited number of improvements.
	if(N==0) return false;
	unsigned long *IPIV = new unsigned long[N];
	int Psign;
	long i,j;
	int impCount;	
	
	//scratch with a copy of the matrix
	for(i=0; i<N*N; ++i){
		scratchSpace[i] = matrix[i];
	}
	
	if(!LUDecomposition(scratchSpace, N, IPIV, Psign)){
		delete[] IPIV;
		return false; // error on LU decomposition
	}
	
	TYPE *b = new TYPE[N];
	TYPE *x = new TYPE[N];
	for(j=0; j<N; ++j){
		for(i=0; i<N; ++i){ x[i] = b[i] = 0; }
		x[j] = b[j] = 1;
		LUSolveLinearSystem(scratchSpace, N, IPIV, x);
		if((maxError>0) && (maxImprovements>0)){
			impCount = 0;
			while((errorInLinearSolution(matrix, N, b, x) > maxError) && ((impCount<maxImprovements) || (maxImprovements==0))){
				LUImproveSolutionToLinearSystem(matrix, scratchSpace, N, IPIV, b, x);
				++impCount;
			}
		}
		for(i=0; i<N; ++i){ inverse[i*N+j] = x[i]; }
	}
	
	
	delete[] b;
	delete[] x;
	delete[] IPIV;
	return true;
}



double get_matrix_trace(const long NR, const std::vector<double> &matrix){
	double T = 0;
	for(long r=0; r<NR; ++r){
		T +=matrix[r*NR+r];
	}
	return T;
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


double sum_of_column(	const long 					NR,			// (INPUT) number of rows in the matrix A
						const long 					NC,			// (INPUT) number of columns in the matrix A
						const std::vector<double>	&A,			// (INPUT) matrix of size NR x NC, in row-major format
						const long 					column){	// (INPUT) focal column, the entries of which are to be summed
	double sum = 0;
	for(long r=0; r<NR; ++r){
		sum += A[r*NC+column];
	}
	return sum;
}

double sum_of_row(	const long 					NR,			// (INPUT) number of rows in the matrix A
					const long 					NC,			// (INPUT) number of columns in the matrix A
					const std::vector<double>	&A,			// (INPUT) matrix of size NR x NC, in row-major format
					const long 					row){	// (INPUT) focal row, the entries of which are to be summed
	double sum = 0;
	for(long c=0; c<NC; ++c){
		sum += A[row*NC+c];
	}
	return sum;
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
	if((NR==2) && (NC==2)){
		// 2 x 2 matrix, treat as special case for computational efficiency
		// this is useful for example in BiSSE models, where a matrix is multiplied with a vector numerous times
		Y.resize(NR);
		Y[0] = A[0]*X[0] + A[1]*X[1];
		Y[1] = A[2]*X[0] + A[3]*X[1];
	}else{
		Y.assign(NR,0);
		for(long r=0; r<NR; ++r){
			for(long c=0; c<NC; ++c){
				Y[r] += A[r*NC+c] * X[c];
			}
		}
	}
}


// multiply matrix A with column-vector X
// matrix is assumed in row-major format
template<class TYPE1, class TYPE2, class TYPE3>
inline void multiply_matrix_with_vector(const long					NR,
										const long					NC,
										const std::vector<TYPE1>	&A,		// (INPUT) matrix of size NR*NC, in row-major format
										const std::vector<TYPE2>	&X,		// (INPUT) std::vector of size NC
										std::vector<TYPE3>			&Y){	// (OUTPUT) product A*X, of size NR
	if((NR==2) && (NC==2)){
		// 2 x 2 matrix, treat as special case for computational efficiency
		// this is useful for example in BiSSE models, where a matrix is multiplied with a vector numerous times
		Y.resize(NR);
		Y[0] = A[0]*X[0] + A[1]*X[1];
		Y[1] = A[2]*X[0] + A[3]*X[1];
	}else{
		Y.assign(NR,0);
		long r, c;
		for(r=0; r<NR; ++r){
			for(c=0; c<NC; ++c){
				Y[r] += A[r*NC+c] * X[c];
			}
		}
	}
}


// Calculate the product Y = X^T * A
template<class TYPE1,class TYPE2,class TYPE3>
void multiply_vector_with_matrix(	const long			NR,
									const long			NC,
									TYPE2				X[],	// (INPUT) pre-allocated array of size NR or greater
									TYPE1				A[],	// (INPUT) array of size NR*NC, in row-major format
									std::vector<TYPE3>	&Y){	// (OUTPUT) product X^T*A, of size NC
	Y.assign(NC,0);
	for(long r=0; r<NR; ++r){
		for(long c=0; c<NC; ++c){
			Y[c] += X[r] * A[r*NC+c];
		}
	}
}




// Calculate the product Y = X^T * A
// matrix is assumed in row-major format
template<class TYPE1, class TYPE2, class TYPE3>
void multiply_vector_with_matrix(	const long					NR,
									const long					NC,
									const std::vector<TYPE2>	&X,		// (INPUT) std::vector of size NR
									const std::vector<TYPE1>	&A,		// (INPUT) matrix of size NR*NC, in row-major format
									std::vector<TYPE3>			&Y){	// (OUTPUT) product X^T*A, of size NC
	Y.assign(NC,0);
	long r,c;
	for(r=0; r<NR; ++r){
		for(c=0; c<NC; ++c){
			Y[c] += X[r] * A[r*NC+c];
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
	AB.assign(NR*NC,0);
	long r,c,k;
	for(r=0; r<NR; ++r){
		for(c=0; c<NC; ++c){
			for(k=0; k<NCa; ++k){
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



// wrapper class for preparing exponentation of a matrix A and exponentiating A using various scalar factors
class matrix_exponentiator{
private:
	// variables for exponentiation using matrix polynomials
	long 				NP;					// number of polynomials available for calculating the exponential of the input matrix
	std::vector<double> polynomials;		// array of size NP x NR x NR, containing the pre-computed polynomials of matrix: Cp:=A^p/p! for p=0,..,NP-1. Polynomials are stored in layer-row-major format, polynomials[p*NR*NR + r*NR + c] is (r,c)-th-element of A^p/p!
	std::vector<double>	polynomial_norms;	// array of size NP, containing the Hilbert-Schmidt L2 norm of each polynomial Cp, ||Cp||_2. Used to reduce the number of incorporated polynomials to the only the necessary number (for optimization reasons).
	long				NPmin;
	double				epsilon;
	bool				balanced;			// whether the input matrix was balanced prior to calculating the polynomials
	std::vector<double>	balances;			// 1D array of size NR, storing the diagonal elements of a diagonal matrix D that was applied for balancing A prior to polynomial calculation. This transformation will be reversed after exponentiation. Only relevant if balanced==true.
	long				scaling_power;		// base-2 scaling power that was applied to the matrix prior to calculating the polynomials. This scaling will be reversed after exponentiation, via repeated squaring. Only relevant if balanced==true.
	
	// variables for exponentiation using eigendecomposition
	bool					use_eigendecomposition;	// whether exponentiation should be performed using eigendecomposition
	mutable std::vector<cdouble>	exponentiation_scratch;
	std::vector<cdouble> 	eigenvalues;
	std::vector<cdouble> 	EVmatrix;
	std::vector<cdouble> 	inverse_EVmatrix;
public:
	bool				initialized;
	long 				NR;					// number of rows & columns in the input matrix
	
	// Constructors
	matrix_exponentiator(){ initialized = false; }
	matrix_exponentiator(	const long 			_NR,
							std::vector<double>	A,					// (INPUT) 2D array of size NR x NR, in row-major format
							const double		rescaling,			// (INPUT) optional scalar scaling factor for input matrix (i.e. use rescaling*A instead of A in all calculations). Set to 1.0 for no rescaling.
							double				_epsilon,			// (INPUT) norm threshold for calculated polynomials C_p=A^p/p!, i.e. stop calculating polynomials as soon as ||exp(A)-sum_{p=0}^{NP-1}C_p||<epsilon. Norm refers to the Hilbert-Schmidt L2 norm.
							const long			_NPmin,				// (INPUT) minimum number of polynomials to calculate if possible (including A^0), regardless of the pursued accuracy epsilon. For sparse Markov transition matrix it is recommended to set this to NR+1, so that the matrix exponential does not contain zeros that it shouldn't contain (assuming A is irreducible). The presence of false zeros in exp(A) can mess up ancestral state reconstruction algorithms.
							const long			NPmax,				// (INPUT) maximum possible number of polynomials to calculate, regardless of the pursued accuracy epsilon. Used as safety vault, but may break the guaranteed accuracy.
							bool				_balanced){
		initialize(_NR, A, rescaling, _epsilon, _NPmin, NPmax, _balanced);
	}
	matrix_exponentiator(	const long 					_NR,
							const std::vector<cdouble>	&_eigenvalues,
							const std::vector<cdouble>	&_EVmatrix,
							const std::vector<cdouble>	&_inverse_EVmatrix,
							const double				rescaling){			// (INPUT) optional scalar scaling factor for input matrix (i.e. use rescaling*A instead of A in all calculations). Set to 1.0 for no rescaling.
		initialize(_NR, _eigenvalues, _EVmatrix, _inverse_EVmatrix, rescaling);
	}
	
	
	// prepare exponentiation of matrix A, by pre-calculating polynomials of A
	// if balanced==true then preps include balancing the matrix. This may be needed for some weird matrixes
	void initialize(const long 			_NR,
					std::vector<double>	A,					// (INPUT) 2D array of size NR x NR, in row-major format
					const double		rescaling,			// (INPUT) optional scalar scaling factor for input matrix (i.e. use rescaling*A instead of A in all calculations). Set to 1.0 for no rescaling.
					double				_epsilon,			// (INPUT) norm threshold for calculated polynomials C_p=A^p/p!, i.e. stop calculating polynomials as soon as ||exp(A)-sum_{p=0}^{NP-1}C_p||<epsilon. Norm refers to the Hilbert-Schmidt L2 norm.
					const long			_NPmin,				// (INPUT) minimum number of polynomials to calculate if possible (including A^0), regardless of the pursued accuracy epsilon. For sparse Markov transition matrix it is recommended to set this to NR+1, so that the matrix exponential does not contain zeros that it shouldn't contain (assuming A is irreducible). The presence of false zeros in exp(A) can mess up ancestral state reconstruction algorithms.
					const long			NPmax,				// (INPUT) maximum possible number of polynomials to calculate, regardless of the pursued accuracy epsilon. Used as safety vault, but may break the guaranteed accuracy.
					bool				_balanced){
		balanced 	= _balanced;
		NR 			= _NR;
		NPmin 		= _NPmin;
		epsilon 	= _epsilon;
		initialized	= true;
		use_eigendecomposition = false;
		if(balanced){
			calculate_balanced_matrix_polynomials(	NR,
													A,
													rescaling,
													epsilon,
													NPmin,
													NPmax,
													polynomials,
													polynomial_norms,
													NP,
													balances,
													scaling_power);
		}else{
			calculate_matrix_polynomials(	NR,
											A,
											rescaling,
											epsilon,
											NPmin,
											NPmax,
											polynomials,
											polynomial_norms,
											NP);
		}
	}
	
	// prepare exponentiation of matrix A based on eigendecomposition
	void initialize(const long 					_NR,
					const std::vector<cdouble>	&_eigenvalues,
					const std::vector<cdouble>	&_EVmatrix,
					const std::vector<cdouble>	&_inverse_EVmatrix,
					const double				rescaling){			// (INPUT) optional scalar scaling factor for input matrix (i.e. use rescaling*A instead of A in all calculations). Set to 1.0 for no rescaling.
		NR 						= _NR;
		initialized				= true;
		use_eigendecomposition 	= true;
		eigenvalues 			= _eigenvalues;
		EVmatrix 				= _EVmatrix;
		inverse_EVmatrix 		= _inverse_EVmatrix;
		if(rescaling!=1.0){
			for(long r=0; r<eigenvalues.size(); ++r) eigenvalues[r] *= rescaling;
		}
	}
	
	
	// calculate exp(tau*A)
	void get_exponential(	double				tau,					// (INPUT) scaling factor in exponent
							std::vector<double>	&exponential) const{	// (OUTPUT) array of size NR x NR, containing the exponentiated matrix exp(tau*A), in row-major format.
		if(use_eigendecomposition){
			get_matrix_exponential_using_eigendecomposition(NR, eigenvalues, EVmatrix, inverse_EVmatrix, tau, exponentiation_scratch, exponential);
		}else if(balanced){
			get_matrix_exponential_using_balanced_polynomials(NR, NP, polynomials, polynomial_norms, tau, epsilon, NPmin, balances, scaling_power, exponential);
		}else{
			get_matrix_exponential_using_polynomials(NR, NP, polynomials, polynomial_norms, tau, epsilon, NPmin, exponential);
		}
	}
};



// get an approximation for the expression exp(scaling*A)*X, where A is a quadratic matrix of size NR x NR, and X is a vector or matrix of size NR x NC
// the approximation is obtained using the expression sum_n=0^order (scaling*A)^n/n! * X
void apply_approximate_matrix_exponential(	const long					NR, 		// (INPUT) number of rows and columns in A, also equal to the number of rows in X
											const long					NC, 		// (INPUT) number of columns in X
											const std::vector<double> 	&A, 		// (INPUT) 2D matrix of size NR x NR, in row-major format. The matrix to be exponentially applied.
											const double				&scaling, 	// (INPUT) an optional multiplicative factor for A, i.e. compute exp(scaling*A) instead of exp(A)*X. Set this to 1 for no scaling.
											const std::vector<double>	&X, 		// (INPUT) the vector or matrix to which exp(A) should be applied
											const long					&order,		// (INPUT) how many polynomials to include in the expansion of exp(A), i.e. 1 + A + A^2/2! + A^3/3! ... A^order/order!
											std::vector<double>			&scratch1,	// (SCRATCH) scratch space needed for internal computation
											std::vector<double>			&scratch2,	// (SCRATCH) scratch space needed for internal computation
											std::vector<double>			&Y){		// (OUTPUT) the resulting vector or matrix of size NR x NC, an approximation to Y=exp(scaling*A)*X
	long n, r, c, k;
	scratch1 = X;
	scratch2.resize(NR*NC);
	Y = X;
	std::vector<double> *source, *target;
	for(n=1; n<=order; ++n){
		source = (n%2==1 ? &scratch1 : &scratch2);
		target = (n%2==0 ? &scratch1 : &scratch2);
		// compute the polynomial term of order n, i.e. ((scaling*A)^n/n!)*X
		// store results in target matrix
		target->assign(NR*NC,0);
		for(r=0; r<NR; ++r){
			for(c=0; c<NC; ++c){
				// compute target := (scaling*A)/n * source
				for(k=0; k<NR; ++k){
					(*target)[r*NC+c] += (scaling/n)*A[r*NR+k]*(*source)[k*NC+c];
				}
			}
		}
		// add target to Y
		for(r=0; r<NR; ++r){
			for(c=0; c<NC; ++c){
				Y[r*NC+c] += (*target)[r*NC+c];
			}
		}
	}
}




// get an approximation for the expression [exp(scaling*A) - Id]*A^{-1}*X, where A is a quadratic matrix of size NR x NR, and X is a vector or matrix of size NR x NC
// the approximation is obtained using the expression sum_n=0^order (scaling/(n+1)) * (scaling*A)^n/n! * X
// This routine is mainly used for the Rosenbrock-Euler ODE solver
void apply_approximate_RosenbrockEuler_exponential(	const long					NR, 		// (INPUT) number of rows and columns in A, also equal to the number of rows in X
													const long					NC, 		// (INPUT) number of columns in X
													const std::vector<double> 	&A, 		// (INPUT) 2D matrix of size NR x NR, in row-major format. The matrix to be exponentially applied.
													const double				&scaling, 	// (INPUT) an optional multiplicative factor for A, i.e. compute exp(scaling*A) instead of exp(A)*X. Set this to 1 for no scaling.
													const std::vector<double>	&X, 		// (INPUT) the vector or matrix to which exp(A) should be applied
													const long					&order,		// (INPUT) how many polynomials to include in the expansion of exp(A), i.e. 1 + A + A^2/2! + A^3/3! ... A^order/order!
													std::vector<double>			&scratch1,	// (SCRATCH) scratch space needed for internal computation
													std::vector<double>			&scratch2,	// (SCRATCH) scratch space needed for internal computation
													std::vector<double>			&Y){		// (OUTPUT) the resulting vector or matrix of size NR x NC, an approximation to Y=exp(scaling*A)*X
	long n, r, c, k;
	scratch1 = X;
	scratch2.resize(NR*NC);
	
	// zeroth-order term
	Y.resize(X.size());
	for(k=0; k<Y.size(); ++k) Y[k] = scaling*X[k];
	
	double REfactor;
	std::vector<double> *source, *target;
	for(n=1; n<=order; ++n){
		source = (n%2==1 ? &scratch1 : &scratch2);
		target = (n%2==0 ? &scratch1 : &scratch2);
		// compute the polynomial term of order n, i.e. ((scaling*A)^n/n!)*X
		// store results in target matrix
		target->assign(NR*NC,0);
		for(r=0; r<NR; ++r){
			for(c=0; c<NC; ++c){
				// compute target := (scaling*A)/n * source
				for(k=0; k<NR; ++k){
					(*target)[r*NC+c] += (scaling/n)*A[r*NR+k]*(*source)[k*NC+c];
				}
			}
		}
		// add target to Y
		REfactor =  (scaling/(n+1.0)); // modifying factor to get Rosenbrock-Euler form, instead of classical exponential
		for(r=0; r<NR; ++r){
			for(c=0; c<NC; ++c){
				Y[r*NC+c] += REfactor * (*target)[r*NC+c];
			}
		}
	}
}



inline double dot_product(const dvector &X, const dvector &Y){
	double S = 0;
	for(long i=0; i<X.size(); ++i) S += X[i]*Y[i];
	return S;
}



// estimate the dominant eigenvalue and corresponding eigenvector of a square matrix, using power iteration
// returns true if convergence was successful after two trials
bool get_dominant_eigenvalue(	const long		N,				// (INPUT) the number of rows & columns in A
								const dvector 	&A,				// (INPUT) square 2D matrix of size N x N, in row-major format
								const long		max_iterations, // (INPUT) maximum number of iterations (matrix multiplications) to perform
								const double	tolerance,		// (INPUT) relative error tolerance for detecting convergence of X (in terms of direction). This is the sine of the angle between two successive iterations of X; if the angle is small, it means the direction of X changes little. A tolerance = 1e-3 is usually sufficient.
								dvector			&X,				// (OUTPUT) an estimate for the dominant eigenvector, normalized to norm 1. This vector is only uniquely defined in terms of direction and norm, but not sign.
								double			&lambda){		// (OUTPUT) an estimate for the dominant eigenvalue
	double XAX, error, best_lambda, best_error;
	dvector AX, Xnew(N), best_X;
	X.resize(N);
	for(int trial=0; trial<min(3L,N); ++trial){ // perform multiple trials, in case the first happens to start with a "bad" direction X
		if(trial==0){
			// generate a random start vector X
			for(long i=0; i<N; ++i) X[i] = R::runif(0.0,1.0);
		}else{
			// pick a new sart vector X that is perpendicular to the previous X
			// step 1: find largest (in magnitude) element in X
			long largest = 0;
			for(long i=0; i<N; ++i){
				if(abs(X[i])>abs(X[largest])) largest = i;
			}
			// step 2: generate N-1 random entries for Xnew, except at the position [largest]
			for(long i=0; i<N; ++i) Xnew[i] = (i==largest ? 0.0 : R::runif(0.0,1.0));
			// step 3: choose Xnew[largest] such that the dot-product X^T*Xnew is zero
			Xnew[largest] = (-1/X[largest]) * dot_product(X,Xnew);
			// step 4: adopt Xnew, and repeat power iteration below
			X = Xnew;
		}

		// normalize X
		double S = sqrt(dot_product(X,X));
		for(long i=0; i<N; ++i) X[i] /= S;
	
		// iterate over X(i+1) = A*X(i)/lambda(i), where lambda(i):=X(i)^T*A*X(i)/||X(i)||_2
		for(long n=0; n<max_iterations; ++n){
			// compute A*X
			AX.assign(N,0);
			for(long r=0; r<N; ++r){
				for(long c=0; c<N; ++c){
					AX[r] += A[r*N+c]*X[c];	
				}
			}
			// compute X^TX
			S = dot_product(X,X);
			// compute X^T*(A*X)
			XAX = 0;
			for(long i=0; i<N; ++i) XAX += X[i]*AX[i];
			if(XAX==0){
				// landed on the zero-eigenvalue, so we're stuck (converged) on this eigenspace
				lambda = 0;
				error  = 0;
				break;
			}
			// compute lambda = X^T*A*X/(X^T*X), the current estimate of the dominant eigenvalue
			lambda = XAX/S;
			// compute next iteration, Xnext = A*X/lambda
			for(long i=0; i<N; ++i) Xnew[i] = AX[i]/lambda;
			// normalize Xnew, just in case
			S = sqrt(dot_product(Xnew,Xnew));
			for(long i=0; i<N; ++i) Xnew[i] /= S;
			// compute error := sin(angle between X and Xnew)
			// an error << 1 implies convergence
			error = sqrt(1 - min(1.0,SQ(dot_product(X,Xnew)))); // this formula relies on the fact that X and Xnew are normalized
			// adopt new X
			X = Xnew;
			if(error<tolerance) break; // reached convergence within relative tolerance
		}
		if((trial==0) || (abs(lambda)>abs(best_lambda))){
			// record the outcome of this trial, regardless of convergence, if lambda is better than the previous ones
			best_lambda = lambda;
			best_error  = error;
			best_X		= X;
		}
	}
	lambda  = best_lambda;
	X		= best_X;
	return (best_error<tolerance);
}




// find the smallest eigenvalue (by magnitude) of a square matrix A
// this method first inverts A using QR decomposition and then uses the power method to compute the dominant eigenvalue of A
// if the matrix is non-invertible, its weakest eigenvalue is estimated to be 0
bool get_weakest_eigenvalue(const long		N,				// (INPUT) the number of rows & columns in A
							const dvector 	&A,				// (INPUT) square 2D matrix of size N x N, in row-major format
							const long		max_iterations, // (INPUT) maximum number of iterations (matrix multiplications) to perform during the power method
							const double	tolerance,		// (INPUT) relative error tolerance for detecting convergence during the power method
							double			&lambda){
	long rank;
	dvector QRA,Ainv,X;
	QR_matrix_inverse(N,A,QRA,Ainv,rank);
	if(rank<N){
		lambda = 0;
		return true;
	}
	// get dominant eigenvalue of Ainv; this will be the inverse of the weakest eigenvalue of A
	const bool converged = get_dominant_eigenvalue(N,Ainv,max_iterations, tolerance, X, lambda);
	lambda = 1/lambda;
	return converged;
}





// balance a square real matrix prior to eigendecomposition
void EIG_balance_matrix(const long 	N, 			// (INPUT) the number of rows & columns in A
						double 		A[], 		// (INPUT/OUTPUT) 2D matrix of size N x N, in column-major format. On input, the matrix to be balanced. On output, the matrix is modified (balanced) in situ.
						long 		&low, 		// (OUTPUT) information on the location of zeros in A
						long 		&igh,		// (OUTPUT) information on the location of zeros in A
						double 		scale[]){	// (OUTPUT) technical information on the permutations and scalings used. Must be preallocated to size N
	double b2, c, f, g, r, radix, s, t;
	bool done, noconv, swap;
	long i, j, k, l, m;

	radix 	= 16.0;
	b2 		= radix * radix;
	j 		= -1;
	m 		= -1;
	k 		= 0;
	l 		= N - 1;

	done = false;
	while(!done){
		for(j = l; 0 <= j; j--){
			swap = true;
			for(i = 0; i <= l; i++){
				if(i != j){
					if(A[j+i*N] != 0.0){
						swap = false;
						break;
					}
				}
			}
			if(swap){
				m = l;
				scale[m] = (double) j;

				if(j != m){
					for(i=0; i <= l; i++){
						t        = A[i+j*N];
						A[i+j*N] = A[i+m*N];
						A[i+m*N] = t;
					}
					for(i=k; i < N; i++){
						t        = A[j+i*N];
						A[j+i*N] = A[m+i*N];
						A[m+i*N] = t;
					}
				}
				if(l == 0){
					low = k;
					igh = l;
					return;
				}

				l = l - 1;
				if(l < 0){
					done = true;
				}
				break;
			}else if(j == 0){
				done = true;
				break;
			}
		}
	}

	done = false;
	while(!done){
		for(j = k; j <= l; j++){
			swap = true;
			for(i = k; i <= l; i++){
				if(i != j){
					if(A[i+j*N] != 0.0){
						swap = false;
						break;
					}
				}
			}
			if(swap){
				m = k;
				scale[m] = (double) j;
				if(j != m){
					for(i=0; i <= l; i++){
						t        = A[i+j*N];
						A[i+j*N] = A[i+m*N];
						A[i+m*N] = t;
					}
					for(i=k; i < N; i++){
						t        = A[j+i*N];
						A[j+i*N] = A[m+i*N];
						A[m+i*N] = t;
					}
				}
				k = k + 1;
				if(l < k){
					done = true;
				}
				break;
			}else{
				if(j == l){
					done = true;
					break;
				}
			}
		}
	}

	// balance submatrix in rows K to L.
	for(i = k; i <= l; i++){
		scale[i] = 1.0;
	}

	// norm reduction.
	noconv = true;
	while(noconv){
		noconv = false;
		for(i = k; i <= l; i++){
			c = 0.0;
			r = 0.0;
			for(j = k; j <= l; j++){
				if(j != i){
					c = c + abs(A[j+i*N]);
					r = r + abs(A[i+j*N]);
				}
			}
	
			// deal with zero C or R due to underflow.
			if(c != 0.0 && r != 0.0){
				g = r / radix;
				f = 1.0;
				s = c + r;
				while(c < g){
					f = f * radix;
					c = c * b2;
				}
				g = r * radix;
				while(g <= c){
					f = f / radix;
					c = c / b2;
				}
				if(( c + r ) / f < 0.95 * s){
					g = 1.0 / f;
					scale[i] = scale[i] * f;
					noconv = true;

					for(j = k; j < N; j++){
						A[i+j*N] = A[i+j*N] * g;
					}
					for(j = 0; j <= l; j++){
						A[j+i*N] = A[j+i*N] * f;
					}
				}
			}
		}
	}

	low = k;
	igh = l;
}




void EIG_reverse_balancing(	const long 		N, 			// (INPUT) the number of rows & columns in the original matrix A
							const long 		low, 		// (INPUT) the low index, as returned by EIG_balance_matrix
							const long 		igh, 		// (INPUT) the igh index, as returned by EIG_balance_matrix
							const double 	scale[], 	// (INPUT) the scale vector of size N, as returned by EIG_balance_matrix
							const long 		M,			// (INPUT) the number of columns of Z to be back-transformed
							double 			Z[]){		// (INPUT/OUTPUT) 2D matrix of size N x M, in column-major format, containing the real & imaginary parts of eigenvectors. Upon return, these will have been back-transformed
	long i, ii, j, k;
	double t;
	if (M <= 0) return;

	if(igh != low){
		for(i = low; i <= igh; i++){
			for(j = 0; j < M; j++){
				Z[i+j*N] = Z[i+j*N] * scale[i];
			}
		}
	}

	for(ii = 0; ii < N; ii++){
		i = ii;
		if(i<low || igh<i){
			if(i < low){
				i = low - ii;
			}
			k = long(scale[i]);
			if(k != i){
				for(j=0; j < M; j++){
					t        = Z[i+j*N];
					Z[i+j*N] = Z[k+j*N];
					Z[k+j*N] = t;
				}
			}
		}
	}
}




// auxiliary function used for EIG_eigendecomposition
void EIG_accumulate_similarity(	const long 	N, 		// (INPUT) the number of rows & column in the original matrix
								const long 	low, 	// (INPUT) The low index, as returned by EIG_balance_matrix. If balancing was not performed, set low = 0, igh = N - 1.
								const long 	igh, 	// (INPUT) The igh index, as returned by EIG_balance_matrix. If balancing was not performed, set low = 0, igh = N - 1.
								double 		H[], 	// (INPUT) 2D matrix of size N x N, in column-major format, containing HELMES multiplies
								const long	ind[], 	// (INPUT) technical information on the ELMHES reduction
								double 		Z[]){	// (OUTPUT) 2D matrix of size N x N, in column-major format, containing the transformation matrix produced by ELMHES. Must be preallocated to size N*N.
	long i,j, mp;
	
	// initialize to the identity matrix
	for(i=0; i<N; ++i){
		for(long j=0; j<N; ++j){
			Z[i + j*N] = 0;	
		}
		Z[i + i*N] = 1;
	}

	if(igh < low + 2){
		return;
	}
	for(mp = igh - 1; low + 1 <= mp; mp--){
		for(i = mp + 1; i <= igh; i++){
			Z[i+mp*N] = H[i+(mp-1)*N];
		}
		i = ind[mp];

		if(i != mp){
			for(j = mp; j <= igh; j++){
				Z[mp+j*N] = Z[i+j*N];
			}
			Z[i+mp*N] = 1.0;
			for(j = mp + 1; j <= igh; j++){
				Z[i+j*N] = 0.0;
			}
		}
	}
}




// transform a real square matrix to upper Hessenberg form, e.g. for subsequent eigendecomposition
// Code adjusted from: https://people.sc.fsu.edu/~jburkardt/c_src/eispack/eispack.html
void EIG_ELMHES(const long 	N, 		// (INPUT) the number of rows & column in the original matrix
				const long 	low, 	// (INPUT) The low index, as returned by EIG_balance_matrix. If balancing was not performed, set low = 1, igh = N.
				const long 	igh, 	// (INPUT) The igh index, as returned by EIG_balance_matrix. If balancing was not performed, set low = 1, igh = N.
				double 		A[], 	// (INPUT/OUTPUT) 2D matrix of size N x N, in column-major format. On input, the matrix to be reduced. On output, the Hessenberg matrix.
				long 		ind[]){	// (OUTPUT) 1D array of size N, containing informationon the rows & columns interchanged. Must be preallocated to size N.
	long i, j, m;
	double t, x, y;

	for(m = low + 1; m <= igh - 1; m++){
		x = 0.0;
		i = m;
		for(j = m; j <= igh; j++){
			if(abs(x) < abs(A[j+(m-1)*N] )){
				x = A[j+(m-1)*N];
				i = j;
			}
		}

		ind[m] = i;
		if(i != m){
			for(j = m - 1; j < N; j++){
				t        = A[i+j*N];
				A[i+j*N] = A[m+j*N];
				A[m+j*N] = t;
			}
			for(j = 0; j <= igh; j++){
				t        = A[j+i*N];
				A[j+i*N] = A[j+m*N];
				A[j+m*N] = t;
			}
		}

		if(x != 0.0){
			for(i = m + 1; i <= igh; i++){
				y = A[i+(m-1)*N];

				if(y != 0.0){
					y = y / x;
					A[i+(m-1)*N] = y;
					for(j = m; j < N; j++){
						A[i+j*N] = A[i+j*N] - y * A[m+j*N];
					}
					for(j = 0; j <= igh; j++){
						A[j+m*N] = A[j+m*N] + y * A[j+i*N];
					}
				}
			}
		}
	}
}




// compute the eigenvalues of a a real square upper Hessenberg matrix, using the QR method
// Returns 0 upon success, or J>0 if convergence was not reached at the J-th eigenvalue
// Code adjusted from: https://people.sc.fsu.edu/~jburkardt/c_src/eispack/eispack.html
long EIG_eigenvalues_RUH(	const long 	N, 			// (INPUT) the number of rows & columns in the matrix
							const long 	low, 		// (INPUT) The low index, as returned by EIG_balance_matrix. If balancing was not performed, set low = 0, igh = N-1
							const long 	igh, 		// (INPUT) The igh index, as returned by EIG_balance_matrix. If balancing was not performed, set low = 0, igh = N-1
							double 		H[],		// (INPUT/OUTPUT) 2D matrix of size N x N, in column-major format, storing the upper Hessenberg matrix. On output, this matrix will be modified, to store information about internal transformations
							double 		lambdaR[], 	// (OUTPUT) 1D array of size N, storing the real part of the eigenvalues. Conjugate pairs are listed consequitively, with the eigenvalue having positive imaginary part listed first. Must be preallocated to size N.
							double		lambdaI[]){	// (OUTPUT) 1D array of size N, storing the imaginary part of the eigenvalues. Must be preallocated to size N.
	long en, enm2, i, ierr, itn, its, j, k, l, m, na;
	double norm, p, q, r, s, t, tst1, tst2, w, x, y, zz;
	bool notlas;

	ierr = 0;
	norm = 0.0;
	k 	 = 0;
	for(i = 0; i < N; i++){
		for(j = k; j < N; j++){
			norm = norm + abs(H[i+j*N]);
		}
		k = i;
		if((i<low) || (igh<i)){
			lambdaR[i] = H[i+i*N];
			lambdaI[i] = 0.0;
		}
	}

	en 	= igh;
	t 	= 0.0;
	itn = 30 * N;

	if(igh < low) return ierr;

	its 	= 0;
	na 		= igh - 1;
	enm2 	= igh - 2;
	while(true){
		for(l = en; low <= l; l--){
			if(l == low) break;
			s = abs(H[l-1+(l-1)*N]) + abs(H[l+l*N]);
			if(s == 0.0) s = norm;
			tst1 = s;
			tst2 = tst1 + abs(H[l+(l-1)*N]);
			if(tst2 == tst1) break;
		}
		x = H[en+en*N];
		if(l == en){
			lambdaR[en] = x + t;
			lambdaI[en] = 0.0;
			en = na;
			if(en < low){
				return ierr;
			}
			its  = 0;
			na   = en - 1;
			enm2 = na - 1;
			continue;
		}
		y = H[na+na*N];
		w = H[en+na*N] * H[na+en*N];
		if(l == na){
			p = (y - x) / 2.0;
			q = p * p + w;
			zz = sqrt(abs(q));
			x = x + t;
			if(0.0 <= q){
				zz = p + abs(zz ) * (p<0 ? -1 : 1);
				lambdaR[na] = x + zz;
				if(zz == 0.0){
					lambdaR[en] = lambdaR[na];
				}else{
					lambdaR[en] = x - w / zz;
				}
				lambdaI[na] = 0.0;
				lambdaI[en] = 0.0;
			}else{
				lambdaR[na] = x + p;
				lambdaR[en] = x + p;
				lambdaI[na] = zz;
				lambdaI[en] = - zz;
			}
			en = enm2;
			if(en < low) return ierr;
			its = 0;
			na = en - 1;
			enm2 = na - 1;
			continue;
		}

		if(itn == 0){
			ierr = en;
			return ierr;
		}
		// exceptional shift.
		if(its == 10 || its == 20){
			t = t + x;
			for(i = low; i <= en; i++){
				H[i+i*N] = H[i+i*N] - x;
			}
			s = abs(H[en+na*N] ) + abs(H[na+enm2*N] );
			x = 0.75 * s;
			y = x;
			w = - 0.4375 * s * s;
		}

		its = its + 1;
		itn = itn - 1;
		for(m = enm2; l <= m; m--){
			zz = H[m+m*N];
			r = x - zz;
			s = y - zz;
			p = ( r * s - w ) / H[m+1+m*N] + H[m+(m+1)*N];
			q = H[m+1+(m+1)*N] - zz - r - s;
			r = H[m+2+(m+1)*N];
			s = abs(p) + abs(q) + abs(r);
			p = p / s;
			q = q / s;
			r = r / s;
			if(m == l) break;
			tst1 = abs(p) * (abs(H[m-1+(m-1)*N]) + abs(zz) + abs(H[m+1+(m+1)*N]));
			tst2 = tst1 + abs(H[m+(m-1)*N]) * (abs(q) + abs(r));
			if(tst2 == tst1) break;
		}

		for(i = m + 2; i <= en; i++){
			H[i+(i-2)*N] = 0.0;
			if(i != m + 2){
				H[i+(i-3)*N] = 0.0;
			}
		}

		// double QR step
		for(k = m; k <= na; k++){
			notlas = (k != na);
			if(k != m){
				p = H[k+(k-1)*N];
				q = H[k+1+(k-1)*N];
				if(notlas){
					r = H[k+2+(k-1)*N];
				}else{
					r = 0.0;
				}
				x = abs(p) + abs(q) + abs(r);
				if(x == 0.0) continue;
				p = p / x;
				q = q / x;
				r = r / x;
			}
			s = sqrt( p * p + q * q + r * r ) * (p<0 ? -1 : 1);
			if(k != m){
				H[k+(k-1)*N] = - s * x;
			}else if(l != m){
				H[k+(k-1)*N] = - H[k+(k-1)*N];
			}
			p  = p + s;
			x  = p / s;
			y  = q / s;
			zz = r / s;
			q  = q / p;
			r  = r / p;
		
			// row modification
			if(! notlas){
				for(j = k; j < N; j++){
					p = H[k+j*N] + q * H[k+1+j*N];
					H[k+j*N] = H[k+j*N] - p * x;
					H[k+1+j*N] = H[k+1+j*N] - p * y;
				}
				j = min(en, k + 3);
				// column modification
				for(i = 0; i <= j; i++){
					p = x * H[i+k*N] + y * H[i+(k+1)*N];
					H[i+k*N] = H[i+k*N] - p;
					H[i+(k+1)*N] = H[i+(k+1)*N] - p * q;
				}
			}else{
				for(j = k; j < N; j++){
					p = H[k+j*N] + q * H[k+1+j*N] + r * H[k+2+j*N];
					H[k+j*N]   = H[k+j*N] - p * x;
					H[k+1+j*N] = H[k+1+j*N] - p * y;
					H[k+2+j*N] = H[k+2+j*N] - p * zz;
				}
				j = min(en, k + 3);
				// column modification
				for(i = 0; i <= j; i++){
					p = x * H[i+k*N] + y * H[i+(k+1)*N] + zz * H[i+(k+2)*N];
					H[i+k*N] 	 = H[i+k*N] - p;
					H[i+(k+1)*N] = H[i+(k+1)*N] - p * q;
					H[i+(k+2)*N] = H[i+(k+2)*N] - p * r;
				}
			}
		}
	}
	return ierr;
}



// divide two complex numbers a and b, each represented as a real tuples (real & imaginary part)
// the result, c:=a/b is also returned as a real tuple
void divide_complex(const double ar, const double ai, const double br, const double bi, double &cr, double &ci){
	double ais, ars, bis, brs, s;
	s = abs(br) + abs(bi);

	ars = ar / s;
	ais = ai / s;
	brs = br / s;
	bis = bi / s;

	s = brs * brs + bis * bis;
	cr = (ars * brs + ais * bis) / s;
	ci = (ais * brs - ars * bis) / s;
}



// calculate eigenvalues & eigenvectors of a square real upper Hessenberg matrix, using the QR method
// returns 0 upon success, or an integer J>0 if convergence failed for the J-th eigenvalue.
// Code adjusted from: https://people.sc.fsu.edu/~jburkardt/c_src/eispack/eispack.html
long EIG_eigenvalues_RUH2(	const long 	N, 			// (INPUT) the number of rows & columns in the matrix
							const long 	low, 		// (INPUT) The low index, as returned by EIG_balance_matrix. If balancing was not performed, set low = 0, igh = N-1
							const long 	igh, 		// (INPUT) The igh index, as returned by EIG_balance_matrix. If balancing was not performed, set low = 0, igh = N-1
							double 		H[],		// (INPUT/OUTPUT) 2D matrix of size N x N, in column-major format, storing the upper Hessenberg matrix. On output, this matrix will be modified, to store information about internal transformations
							double 		lambdaR[], 	// (OUTPUT) 1D array of size N, storing the real part of the eigenvalues. Conjugate pairs are listed consequitively, with the eigenvalue having positive imaginary part listed first. Must be preallocated to size N. If an error occurred (ierr>0), the eigenvalues will be correct for indices ierr+1, .., N.
							double		lambdaI[],	// (OUTPUT) 1D array of size N, storing the imaginary part of the eigenvalues. Must be preallocated to size N.
							double		Z[]){		// (INPUT/OUTPUT) 2D matrix of size N x N, in column-major format. On input, contains the ELTRAN or ORTRAN trasformation, if performed. If the eigenvectors of H are requested, this must be the identity matrix. On output, Z contains the real & imaginary part of the eigenvectors. If the function returns with an error code, none of the returned eigenvectors are valid.
	long en, enm2, i, ierr, itn, its, j, k, l, m, na;
	double norm, p, q, r, ra, s, sa, t, ti, tr, tst1, tst2, vi, vr, w, x, y, zz;
	bool notlas;

	ierr = 0;
	norm = 0.0;
	k = 0;

	for(i = 0; i < N; i++){
		for(j = k; j < N; j++){
			norm = norm + abs(H[i+j*N]);
		}
		k = i;
		if(i < low || igh < i){
			lambdaR[i] = H[i+i*N];
			lambdaI[i] = 0.0;
		}
	}

	en = igh;
	t = 0.0;
	itn = 30 * N;
	while(low <= en){
		its = 0;
		na = en - 1;
		enm2 = na - 1;
		while(true){
			for(l = en; low <= l; l--){
				if(l == low) break;
				s = abs(H[l-1+(l-1)*N]) + abs(H[l+l*N]);
				if(s == 0.0) s = norm;
				tst1 = s;
				tst2 = tst1 + abs(H[l+(l-1)*N] );
				if(tst2 == tst1) break;
			}

			x = H[en+en*N];
			// found one root
			if(l == en){
				H[en+en*N] = x + t;
				lambdaR[en] = H[en+en*N];
				lambdaI[en] = 0.0;
				en = na;
				break;
			}
			y = H[na+na*N];
			w = H[en+na*N] * H[na+en*N];

			// found 2 roots
			if(l == na){
				p = (y - x) / 2.0;
				q = p * p + w;
				zz = sqrt(abs(q ) );
				H[en+en*N] = x + t;
				x = H[en+en*N];
				H[na+na*N] = y + t;

				if(q < 0.0){
					lambdaR[na] = x + p;
					lambdaR[en] = x + p;
					lambdaI[na] = zz;
					lambdaI[en] = - zz;
				}else{
					zz = p + abs(zz ) * (p<0 ? -1 : 1);
					lambdaR[na] = x + zz;
					lambdaR[en] = lambdaR[na];
					if(zz != 0.0 ){
						lambdaR[en] = x - w / zz;
					}
					lambdaI[na] = 0.0;
					lambdaI[en] = 0.0;
					x = H[en+na*N];
					s = abs(x ) + abs(zz );
					p = x / s;
					q = zz / s;
					r = sqrt(p * p + q * q );
					p = p / r;
					q = q / r;
					// row modification.
					for(j = na; j < N; j++){
						zz = H[na+j*N];
						H[na+j*N] = q * zz + p * H[en+j*N];
						H[en+j*N] = q * H[en+j*N] - p * zz;
					}
					// column modification.
					for(i = 0; i <= en; i++){
						zz = H[i+na*N];
						H[i+na*N] = q * zz + p * H[i+en*N];
						H[i+en*N] = q * H[i+en*N] - p * zz;
					}
					// accumulate transformations.
					for(i = low; i <= igh; i++){
						zz = Z[i+na*N];
						Z[i+na*N] = q * zz + p * Z[i+en*N];
						Z[i+en*N] = q * Z[i+en*N] - p * zz;
					}
				}
				en = enm2;
				break;
			}

			if(itn == 0 ){
				ierr = en;
				return ierr;
			}

			// exceptional shift.
			if(its == 10 || its == 20){
				t = t + x;
				for(i = low; i <= en; i++){
					H[i+i*N] = H[i+i*N] - x;
				}
				s = abs(H[en+na*N] ) + abs(H[na+enm2*N] );
				x = 0.75 * s;
				y = x;
				w = - 0.4375 * s * s;
			}

			its = its + 1;
			itn = itn - 1;
			for(m = enm2; l <= m; m--){
				zz = H[m+m*N];
				r = x - zz;
				s = y - zz;
				p = ( r * s - w ) / H[m+1+m*N] + H[m+(m+1)*N];
				q = H[m+1+(m+1)*N] - zz - r - s;
				r = H[m+2+(m+1)*N];
				s = abs(p ) + abs(q ) + abs(r );
				p = p / s;
				q = q / s;
				r = r / s;
				if(m == l) break;
				tst1 = abs(p ) * ( abs(H[m-1+(m-1)*N] ) + abs(zz ) 
				+ abs(H[m+1+(m+1)*N] ) );
				tst2 = tst1 + abs(H[m+(m-1)*N] ) * ( abs(q ) + abs(r ) );
				if(tst2 == tst1) break;
			}

			for(i = m + 2; i <= en; i++){
				H[i+(i-2)*N] = 0.0;
				if(i != m + 2 ){
					H[i+(i-3)*N] = 0.0;
				}
			}

			// double QR step
			for(k = m; k <= na; k++){
				notlas = (k != na);
				if(k != m){
					p = H[k+(k-1)*N];
					q = H[k+1+(k-1)*N];
					r = 0.0;
					if(notlas) r = H[k+2+(k-1)*N];
					x = abs(p ) + abs(q ) + abs(r );
					if(x == 0.0) continue;
					p = p / x;
					q = q / x;
					r = r / x;
				}
				s = sqrt(p * p + q * q + r * r) * (p<0 ? -1 : 1);
				if(k != m){
					H[k+(k-1)*N] = - s * x;
				}else if(l != m){
					H[k+(k-1)*N] = - H[k+(k-1)*N];
				}

				p = p + s;
				x = p / s;
				y = q / s;
				zz = r / s;
				q = q / p;
				r = r / p;

				if(!notlas){
					// row modification
					for(j = k; j < N; j++){
						p = H[k+j*N] + q * H[k+1+j*N];
						H[k+j*N] = H[k+j*N] - p * x;
						H[k+1+j*N] = H[k+1+j*N] - p * y;
					}
					j = min(en, k + 3);
					// column modification
					for(i = 0; i <= j; i++){
						p = x * H[i+k*N] + y * H[i+(k+1)*N];
						H[i+k*N] = H[i+k*N] - p;
						H[i+(k+1)*N] = H[i+(k+1)*N] - p * q;
					}
					// accumulate transformations
					for(i = low; i <= igh; i++){
						p = x * Z[i+k*N] + y * Z[i+(k+1)*N];
						Z[i+k*N] = Z[i+k*N] - p;
						Z[i+(k+1)*N] = Z[i+(k+1)*N] - p * q;
					}
				}else{
					// row modification
					for(j = k; j < N; j++){
						p = H[k+j*N] + q * H[k+1+j*N] + r * H[k+2+j*N];
						H[k+j*N] = H[k+j*N] - p * x;
						H[k+1+j*N] = H[k+1+j*N] - p * y;
						H[k+2+j*N] = H[k+2+j*N] - p * zz;
					}
					j = min(en, k + 3);
					// column modification
					for(i = 0; i <= j; i++){
						p = x * H[i+k*N] + y * H[i+(k+1)*N] + zz * H[i+(k+2)*N];
						H[i+k*N] = H[i+k*N] - p;
						H[i+(k+1)*N] = H[i+(k+1)*N] - p * q;
						H[i+(k+2)*N] = H[i+(k+2)*N] - p * r;
					}
					// accumulate transformations
					for(i = low; i <= igh; i++){
						p = x * Z[i+k*N] + y * Z[i+(k+1)*N] + zz * Z[i+(k+2)*N];
						Z[i+k*N] = Z[i+k*N] - p;
						Z[i+(k+1)*N] = Z[i+(k+1)*N] - p * q;
						Z[i+(k+2)*N] = Z[i+(k+2)*N] - p * r;
					}
				}
			}
		}
	}

	// all roots found, now backsubstitute
	if(norm == 0.0) return ierr;
	for(en = N - 1; 0 <= en; en--){
		p = lambdaR[en];
		q = lambdaI[en];
		na = en - 1;

		if(0.0 < q ){
			continue;
		}else if(q == 0.0){
			m = en;
			H[en+en*N] = 1.0;

			for(i = en - 1; en - na - 1 <= i; i--){
				w = H[i+i*N] - p;
				r = 0.0;
				for(j = m; j <= en; j++){
					r = r + H[i+j*N] * H[j+en*N];
				}
				if(lambdaI[i] < 0.0 ){
					zz = w;
					s = r;
					continue;
				}
				m = i;
				if(lambdaI[i] == 0.0){
					t = w;
					if(t == 0.0 ){
						tst1 = norm;
						t = tst1;
						while(true){
							t = 0.01 * t;
							tst2 = norm + t;
							if(tst2 <= tst1) break;
						}
					}
					H[i+en*N] = - r / t;
				}else{
					x = H[i+(i+1)*N];
					y = H[i+1+i*N];
					q = (lambdaR[i] - p) * (lambdaR[i] - p) + lambdaI[i] * lambdaI[i];
					t = (x * s - zz * r) / q;
					H[i+en*N] = t;
					if(abs(zz) < abs(x)){
						H[i+1+en*N] = ( - r - w * t ) / x;
					}else{
						H[i+1+en*N] = ( - s - y * t ) / zz;
					}
				}
				// overflow control
				t = abs(H[i+en*N]);
				if(t != 0.0 ){
					tst1 = t;
					tst2 = tst1 + 1.0 / tst1;
					if(tst2 <= tst1){
						for(j = i; j <= en; j++){
							H[j+en*N] = H[j+en*N] / t;
						}
					}
				}
			}
		}else if(q < 0.0){
			// complex vector
			m = na;
			if(abs(H[na+en*N]) < abs(H[en+na*N])){
				H[na+na*N] = q / H[en+na*N];
				H[na+en*N] = - ( H[en+en*N] - p ) / H[en+na*N];
			}else{
				divide_complex(0.0, -H[na+en*N], H[na+na*N] - p, q, tr, ti);
				H[na+na*N] = tr;
				H[na+en*N] = ti;
			}
			H[en+na*N] = 0.0;
			H[en+en*N] = 1.0;
			enm2 = na - 1;
			for(i = na - 1; na - enm2 <= i; i--){
				w = H[i+i*N] - p;
				ra = 0.0;
				sa = 0.0;
				for(j = m; j <= en; j++){
					ra = ra + H[i+j*N] * H[j+na*N];
					sa = sa + H[i+j*N] * H[j+en*N];
				}
				if(lambdaI[i] < 0.0 ){
					zz = w;
					r = ra;
					s = sa;
				}
				m = i;
				if(lambdaI[i] == 0.0 ){
					divide_complex(-ra, -sa, w, q, tr, ti);
					H[i+na*N] = tr;
					H[i+en*N] = ti;
				}else{
					x = H[i+(i+1)*N];
					y = H[i+1+i*N];
					vr = (lambdaR[i] - p) * (lambdaR[i] - p) + lambdaI[i] * lambdaI[i] - q * q;
					vi = (lambdaR[i] - p) * 2.0 * q;

					if(vr == 0.0 && vi == 0.0){
						tst1 = norm * (abs(w) + abs(q) + abs(x) + abs(y) + abs(zz));
						vr = tst1;
						while(true ){
							vr = 0.01 * vr;
							tst2 = tst1 + vr;
							if(tst2 <= tst1) break;
						}
					}

					divide_complex(x * r - zz * ra + q * sa, x * s - zz * sa - q * ra, vr, vi, tr, tr);
					H[i+na*N] = tr;
					H[i+en*N] = ti;
					if(abs(zz) + abs(q) < abs(x)){
						H[i+1+na*N] = (- ra - w * H[i+na*N] + q * H[i+en*N]) / x;
						H[i+1+en*N] = (- sa - w * H[i+en*N] - q * H[i+na*N]) / x;
					}else{
						divide_complex(- r - y * H[i+na*N], - s - y * H[i+en*N], zz, q, tr, ti);
						H[i+1+na*N] = tr;
						H[i+1+en*N] = ti;
					}
				}
				// overflow control.
				t = max(abs(H[i+na*N]), abs(H[i+en*N]));

				if(t != 0.0){
					tst1 = t;
					tst2 = tst1 + 1.0 / tst1;
					if(tst2 <= tst1){
						for(j = i; j <= en; j++){
							H[j+na*N] = H[j+na*N] / t;
							H[j+en*N] = H[j+en*N] / t;
						}
					}
				}
			}
		}
	}
	
	// end back substitution
	for(i = 0; i < N; i++){
		if(i < low || igh < i){
			for(j = i; j < N; j++){
				Z[i+j*N] = H[i+j*N];
			}
		}
	}
	
	// multiply by transformation matrix to obtain vectors of original matrix
	for(j = N - 1; low <= j; j--){
		m = min(j, igh);
		for(i = low; i <= igh; i++){
			zz = 0.0;
			for(k = low; k <= m; k++){
				zz = zz + Z[i+k*N] * H[k+j*N];
			}
			Z[i+j*N] = zz;
		}
	}

	return ierr;
}






// calculate eigenvalues and (optionally) eigenvectors of a real square matrix A
// Internally, A is first transformed into upper Hessenberg form, and then eigenvalues & eigenvectors are computed using the QR method
// returns 0 upon success, otherwise the error code is as returned by EIG_eigenvalues_RUH() or EIG_eigenvalues_RUH2()
// Code adjusted from: https://people.sc.fsu.edu/~jburkardt/c_src/eispack/eispack.html
long EIG_eigendecomposition(const long				N,						// (INPUT) the number of rows & columns in A
							const dvector 			&A,						// (INPUT) 2D matrix of size N x N, in row-major or column-major format
							const bool				row_major,				// (INPUT) whether the matrix A is in row-major or column-major format. The same formatting will also apply to the output. Note that internally this function transforms everything into column-major format if needed.
							const bool				include_eigenvectors,	// (INPUT) whether to also compute eigenvectors
							dvector					&scratchA,				// (SCRATCH) scratch space, will be resized to N x N
							dvector					&scratchZ,				// (SCRATCH) scratch space for eigenvector computation, will be resized to N x N. Only relevant if include_eigenvectors==true.
							dvector					&eigenvaluesR,			// (OUTPUT) the real part of the eigenvalues
							dvector					&eigenvaluesI,			// (OUTPUT) the imaginary part of the eigenvalues
							std::vector<cdouble>	&eigenvectors){			// (OUTPUT) 2D matrix of size N x N, containing the eigenvectors (one eigenvector per column). Only relevant if include_eigenvectors==true. Will be in row-major format iff row_major=true, otherwise it will be in column-major format.
	if(row_major){
		// transform A into column-major format and store result in scratchA
		scratchA.resize(N*N);
		for(long r=0; r<N; ++r){
			for(long c=0; c<N; ++c){
				scratchA[r + c*N] = A[r*N + c];
			}
		}
	}else{
		// A is already in column-major format
		scratchA = A;
	}
	
	int ierr;
	long is1, is2;
	dvector fv1(N);
	lvector iv1(N);
	eigenvaluesR.resize(N);
	eigenvaluesI.resize(N);
	EIG_balance_matrix(N, &scratchA[0], is1, is2, &fv1[0]);
	EIG_ELMHES(N, is1, is2, &scratchA[0], &iv1[0]);

	if(!include_eigenvectors){
		ierr = EIG_eigenvalues_RUH(N, is1, is2, &scratchA[0], &eigenvaluesR[0], &eigenvaluesI[0]);
		if(ierr!=0) return ierr;  // an error occurred during HQR
	}else{
		scratchZ.resize(N*N); // extra scratch space for eigenvectors
		EIG_accumulate_similarity(N, is1, is2, &scratchA[0], &iv1[0], &scratchZ[0]);
		ierr = EIG_eigenvalues_RUH2(N, is1, is2, &scratchA[0], &eigenvaluesR[0], &eigenvaluesI[0], &scratchZ[0]);
		if(ierr!=0) return ierr; // an error occurred during HQR2
		EIG_reverse_balancing(N, is1, is2, &fv1[0], N, &scratchZ[0]);
		
		// compactify eigenvectors to complex numbers
		// note that scratchZ is structured in column-major format, and each column (or pair of columns) corresponds to an eigenvector
		eigenvectors.resize(N*N);
		for(long j=0; j<N; ++j){
			if(eigenvaluesI[j]==0){
				// this eigenvalue is real, so the eigenvector is stored in the j-th column of Z
				for(long r=0; r<N; ++r){
					eigenvectors[r+j*N] = complex<double>(scratchZ[r+j*N],0);
				}
			}else{
				// this eigenvalue is complex, so the eigenvector is stored in the j-th and (j+1)-th columns of Z (real & imaginary part, respectively). 
				// the conjugate of this eigenvector is the eigenvector for the conjugate eigenvalue
				for(long r=0; r<N; ++r){
					eigenvectors[r+j*N]		= complex<double>(scratchZ[r+j*N],+scratchZ[r+(j+1)*N]);
					eigenvectors[r+(j+1)*N]	= complex<double>(scratchZ[r+j*N],-scratchZ[r+(j+1)*N]);
				}				
				++j; // skip next eigenvalue, since it is the complex-conjugate of the j-th eigenvalue
			}
		}
	}
	
	// transform eigenvectors from column-major to row-major format if needed
	// this basically means transposing the eigenvector matrix
	if(include_eigenvectors && row_major){
		cdouble temp;
		for(long r=1; r<N; ++r){
			for(long c=0; c<r; ++c){
				temp = eigenvectors[r*N + c];
				eigenvectors[r*N + c] = eigenvectors[r + c*N];
				eigenvectors[r + c*N] = temp;
			}
		}
	}
	return ierr;
}



// get spectral range of a real square general matrix, i.e. the difference between the largest (most positive) real part of an eigenvalue and the smallest (most negative) real part of an eigenvalue
// That is, if (R1,I1),...,(RN,IN) are the eigenvalues of the matrix (where Ri and Ii is the real & imaginary part of the i-th eigenvalue), calculate max_{i,j} (Ri-Rj).
// The input matrix A can be in row-major or column-major format; the result is the same regardless (since confusing majority simply corresponds to transposing A)
// returns NAN_D upon failure
double get_spectral_range(	const long 		N,		// (INPUT) the number of rows & columns in A
							const dvector	&A){	// (INPUT) 2D matrix of size N x N, in row-major or column-major format (it does not matter which)
	// get all eigenvalues of A
	std::vector<cdouble> eigenvectors;
	dvector scratchA, scratchZ, eigenvaluesR, eigenvaluesI;
	const long error = EIG_eigendecomposition(N, A, false, false, scratchA, scratchZ, eigenvaluesR, eigenvaluesI, eigenvectors);
	if(error!=0) return NAN_D;
	
	// find most positive and most negative part of any eigenvalue
	double maxR = eigenvaluesR[0], minR = eigenvaluesR[0];
	for(long i=0; i<N; ++i){
		maxR = max(maxR, eigenvaluesR[i]);
		minR = min(minR, eigenvaluesR[i]);
	}
	
	return (maxR-minR);
}



#pragma mark -
#pragma mark Stochastic processes
#pragma mark 



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




long get_next_Mk_state(	const matrix_exponentiator 	&transition_matrix_exponentiator,
						std::vector<double>			&scratch_exp,					// (SCRATCH) scratch space used to store the temporary exponential
						const double 				dt,
						const long 					previous_state){
	const long Nstates = transition_matrix_exponentiator.NR;
	transition_matrix_exponentiator.get_exponential(dt, scratch_exp);
	// use row of exponentiated transition matrix corresponding to the previous state, as probability vector for the next state
	const long next_state = random_int_from_distribution<double>(&scratch_exp[previous_state*Nstates+0], Nstates);
	return next_state;
}



// get next discrete state of a Markov chain, based on the transition rate matrix, and conditionally upon a single transition having occurred
// for efficiency, the caller guarantees that total_transition_rate equals the sum of transition rates from old_state
long get_next_Mk_state(	const long 					Nstates,					// (INPUT) number of discrete states
						const std::vector<double>	&transition_matrix,			// (INPUT) transition rate matrix. 2D array of size Nstates x Nstates, in row-major format, such that transition_matrix[r*Nstates+c] is the transition rate r-->c
						const double				total_transition_rate,		// (INPUT) equal to sum_{c != old_state} transition_matrix[old_state*Nstates+c]
						const long					old_state){					// (INPUT) old state, from which a single transition is to be simulated
	double p = R::runif(0.0,1.0);
	for(long c=0; c<Nstates; ++c){
		if(p<=transition_matrix[old_state*Nstates+c]/total_transition_rate) return c;
		p -= transition_matrix[old_state*Nstates+c]/total_transition_rate;
	}
	return Nstates-1;
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


double get_average(double values[], long count){
	double S = 0;
	for(long i=0; i<count; ++i) S += values[i];
	return S/count;
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
TYPE_Y interpolateBetween2PointsLinear(const TYPE_X &x1, const TYPE_Y &y1, const TYPE_X &x2, const TYPE_Y &y2, const TYPE_X &x){
	if(x1 == x2){
		return TYPE_Y(NAN_D);
	}else{
		return y1 + (x-x1)*(y2-y1)/(x2-x1);
	}
}




// interpolates an 'old' time series at new time points
// only values from includedNewTimesStart to includedNewTimesEnd (inclusive) will be well-defined; outside values are unpredictable
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
	
	valuesAtNewTimes.resize(newTimes.size());
	
	long tpOld = oldStart;
	double time;
	for(; tpNew<=tpLastNew; ++tpNew){
		time = newTimes[tpNew];
		
		// at this point it is guaranteed that   oldTimes[tpOld]<=time   and that oldTimes[oldEnd]>=time
		
		//find left-closest oldTime to time
		while((tpOld+1<=oldEnd) && (oldTimes[tpOld+1]<=time)){ ++tpOld; }
		
		// at this point it is guaranteed that either (tpOld==oldEnd and oldTimes[tpOld]==time), or (tpOld<oldEnd and oldTimes[tpOld]<=time and oldTimes[tpOld+1]>time)
				
		// interpolate old series at focal time
		valuesAtNewTimes[tpNew] = (tpOld==oldEnd ? valuesAtOldTimes[oldEnd] : interpolateBetween2PointsLinear(oldTimes[tpOld], valuesAtOldTimes[tpOld], oldTimes[tpOld+1], valuesAtOldTimes[tpOld+1], time));
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
	
	if(!LUsolveLinearSystem(	&matrix[0],
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
	
	if(!LUsolveLinearSystem(	&matrix[0],
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
	
	if(!LUsolveLinearSystem(	&matrix[0],
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
// for details see: http://www.mathscratchs.com/help/curvefit/smoothing-data.html
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
										const bool			allow_less_data_at_edges, // only relevant if windowIndexSpan is used
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
			if(!allow_less_data_at_edges){
				if(nl<0) nr += (0-nl); // shift to the right since we're at the left edge
				else if(nr>=N) nl -= (nr-N+1); // shift to the left since we're at the right edge
			}
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
			case 0:
				smoothenedData[n] = get_average(&localData[0], M);
				break;
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
	const bool success = smoothenTimeSeriesSavitzkyGolay(times, data, windowTimeSpan, windowIndexSpan, order, true, smoothened_data);
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


template<class TYPE>
inline void convex_combination(const TYPE &y1, const TYPE &y2, const double lambda, TYPE &y){
	y = (1-lambda)*y1 + lambda*y2;
}


// specialize for vector-valued y
template<>
inline void convex_combination(const std::vector<double> &y1, const std::vector<double> &y2, const double lambda, std::vector<double> &y){
	const long N = y1.size();
	y.resize(N);
	for(long i=0; i<N; ++i){
		y[i] = (1-lambda)*y1[i] + lambda*y2[i];
	}
}




//Class for defining piecewiese linear function using a regular or irregular reference grid
//In the case of a regular grid, value retrieval is of time complexity O(1)
//In the case of an irregular grid, value retrieval is fastest when subsequent requested points are close to each other, because the search for the enclosing grid cell starts from the previous requested grid cell
template<class VALUE_TYPE>
class LinearInterpolationFunctor{
private:
	std::vector<double> 	referencePoints; // store reference points, in the case of an irregular grid. Will be empty iff the grid was regular during setup.
	std::vector<VALUE_TYPE>	referenceValues;
	
	double		domain_min, domain_max;
	double		domainStep;
	double		lengthScale;
	
	bool		periodic;
	VALUE_TYPE	outlier_value_left, outlier_value_right;

	// keep track of the referencePoint matched to the last requested time point
	// This helps find the next reference point faster in the case of non-regular grids, 
	//   assuming that adjacent requests are for times in close proximity (as is the case in ODE solvers)
	// last_requested_reference will correspond to the referencePoint at or immediately to the left (i.e. below) the last requested time
	mutable long last_requested_reference;
	
	void 	set_to_regular_grid_values(	long 				referenceCount,
										const double		domain_min, 
										const double		domain_max, 
										VALUE_TYPE 			referenceValues[],
										bool				periodic,
										const VALUE_TYPE	&outlier_value_left,	// value to use for extending time series to the left if needed, if periodic==false. Irrelevant if periodic==true.
										const VALUE_TYPE	&outlier_value_right);	// value to use for extending time series to the right if needed, if periodic==false. Irrelevant if periodic==true.
public:

	LinearInterpolationFunctor(): domain_min(0), domain_max(0), domainStep(0), lengthScale(1), periodic(false), last_requested_reference(-1) {}
	
	// define piecewise linear function, on a regular or irregular gird
	// domain is assumed to span from lowest to highest referencePoints
	// If coordinate is periodic, the last reference point is identified with the first one (its reference value is taken to be (referenceValues[0]+referenceValues[last])/2)
	// Optionally, the function can be internally re-mapped onto a regular grid, for more efficient later value retrieval
	// In the case of an irregular grid (preInterpolateOnRegularGrid=false), later value retrieval is fastest if sequentially requested points are close to each other, because the search starts from the last requested point
	//   Hence, for example, if an ODE solver requests values at non-decreasing points, then value retrieval will be of time-complexity O(N/R) where N is the number of stored reference points and R is the total number of requested points
	LinearInterpolationFunctor(	const std::vector<double> 		&referencePoints, // domain grid points, in ascending order
								const std::vector<VALUE_TYPE> 	&referenceValues,
								bool							periodic,	// if true, time series are extended periodically if needed (e.g. if the reference domain is too small). Otherwise they are extended with zeros.
								const VALUE_TYPE				&outlier_value_left,
								const VALUE_TYPE				&outlier_value_right,
								bool							preInterpolateOnRegularGrid, // if true, then the function is internally pre-interpolated on a regular domain grid. In that case, value retrieval later on will be more efficient.
								double							regularGridStep);	// regular reference-grid step to use, if preInterpolateOnRegularGrid==true. If <=0, then the regular grid step is chosen as the mean step in the (irregular) input grid.
	
	//define piecewise linear function on a regular grid
	//referenceValues[i] should correspond to domain_min + i * intervalWidth(domain)/(referenceCount-1)
	//in particular, referenceValue[referenceCount-1] should correspond to domain_max
	//if coordinate is periodic, the last reference point is identified with the first one (its reference value is taken to be (referenceValues[0]+referenceValues[last])/2)
	//since the domain grid is assumed to be regular, interpolations defined this way are much for efficient
	LinearInterpolationFunctor(	long 				referenceCount, 
								const double		domain_min, 
								const double		domain_max, 
								VALUE_TYPE 			referenceValues[],
								bool				periodic,
								const VALUE_TYPE	&outlier_value_left,
								const VALUE_TYPE	&outlier_value_right);
							
	VALUE_TYPE operator()(double x) const;
	void getValue(double x, VALUE_TYPE &y) const; // equivalent to operator(), but avoiding copy operators
	
	bool isNULL() const{ return referenceValues.empty(); }
	bool isPeriodic() const{ return periodic; }
	void getDomain(double &_domain_min, double &_domain_max) const{ _domain_min=domain_min; _domain_max=domain_max; }
	long getReferenceCount() const{ return referenceValues.size(); }
	
	void setTypicalLengthScale(double _lengthScale){ lengthScale = _lengthScale; } 
	double typicalLengthScale() const{ return lengthScale; }
	
	const std::vector<double> &getReferencePoints() const{ return referencePoints; }
	const std::vector<double> &getReferenceValues() const{ return referenceValues; }
	
	inline VALUE_TYPE getFirstReferenceValue() const{ return referenceValues[0]; }
	inline VALUE_TYPE getLastReferenceValue() const{ return referenceValues.back(); }
	inline const VALUE_TYPE* getFirstReferenceValuePointer() const{ return &referenceValues[0]; }
	inline const VALUE_TYPE* getLastReferenceValuePointer() const{ return &referenceValues[referenceValues.size()-1]; }
};





template<class VALUE_TYPE>
LinearInterpolationFunctor<VALUE_TYPE>::LinearInterpolationFunctor(	const std::vector<double> 		&_referencePoints, 
																	const std::vector<VALUE_TYPE> 	&_referenceValues,
																	bool							_periodic,
																	const VALUE_TYPE				&_outlier_value_left,
																	const VALUE_TYPE				&_outlier_value_right,
																	bool							preInterpolateOnRegularGrid,
																	double							regularGridStep){
	periodic 			= _periodic;
	outlier_value_left 	= _outlier_value_left;
	outlier_value_right	= _outlier_value_right;
	referencePoints.clear(); 
	referenceValues.clear();
	last_requested_reference = -1;
	if(_referencePoints.empty()) return;

	if(preInterpolateOnRegularGrid){
		// pre-interpolate on regular domain grid, then setup interpolator functor
		regularGridStep 			= (regularGridStep<=0 ? vector_mean(_referencePoints) : regularGridStep);
		const long NR 				= 1 + (_referencePoints.back() - _referencePoints.front())/regularGridStep;
		const double _domain_min 	= _referencePoints.front();
		std::vector<double> regular_reference_points(NR);
		std::vector<VALUE_TYPE> regular_reference_values;
		for(long i=0; i<NR; ++i) regular_reference_points[i] = _domain_min + i*regularGridStep;
		long includedNewTimesStart, includedNewTimesEnd;
		interpolateTimeSeriesAtTimes<VALUE_TYPE,vector<double> >(	_referencePoints,
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
									_outlier_value_left,
									_outlier_value_right);
		
	}else{
		// use provided irregular grid
		referencePoints = _referencePoints;
		referenceValues = _referenceValues;
		const long referenceCount = referencePoints.size();
		
		// figure out domain
		domain_min	= referencePoints.front();
		domain_max	= referencePoints.back();
		lengthScale	= domain_max - domain_min;
		
		// force periodic boundaries if necessary
		if(periodic){
			referenceValues[0] = referenceValues[referenceCount-1] = 0.5*(referenceValues[0] + referenceValues[referenceCount-1]);
		}
	}
}


template<class VALUE_TYPE>
void LinearInterpolationFunctor<VALUE_TYPE>::set_to_regular_grid_values(long 				referenceCount, 
																		const double 		_domain_min, 
																		const double 		_domain_max, 
																		VALUE_TYPE 			*_referenceValues,
																		bool				_periodic,
																		const VALUE_TYPE	&_outlier_value_left,
																		const VALUE_TYPE	&_outlier_value_right){

	periodic			= _periodic;
	domain_min			= _domain_min;
	domain_max			= _domain_max;
	lengthScale			= (domain_max-domain_min);
	domainStep			= (domain_max-domain_min)/max(1.0, referenceCount-1.0);
	outlier_value_left 	= _outlier_value_left;
	outlier_value_right	= _outlier_value_right;
	referencePoints.clear(); 
	referenceValues.clear();
	last_requested_reference = -1;
	if(referenceCount == 0) return;
		
	referenceValues.resize(referenceCount);
	for(long i=0; i<referenceCount; ++i){
		referenceValues[i] = _referenceValues[i];
	}
		
	// enforce periodic boundary values if necessary
	if(periodic){
		referenceValues[0] = referenceValues[referenceCount-1] = (referenceValues[0] + referenceValues[referenceCount-1])/2.0;
	}
}


template<class VALUE_TYPE>
LinearInterpolationFunctor<VALUE_TYPE>::LinearInterpolationFunctor(	long 				referenceCount, 
																	const double 		_domain_min, 
																	const double 		_domain_max, 
																	VALUE_TYPE 			_referenceValues[],
																	bool				_periodic,
																	const VALUE_TYPE	&_outlier_value_left,
																	const VALUE_TYPE	&_outlier_value_right){
	set_to_regular_grid_values(referenceCount, _domain_min, _domain_max, _referenceValues, _periodic, _outlier_value_left, _outlier_value_right);
}



template<class VALUE_TYPE>
void LinearInterpolationFunctor<VALUE_TYPE>::getValue(double x, VALUE_TYPE &y) const{
	if(referenceValues.empty()){
		//nothing available
		y = outlier_value_left;
		return;
	}
	const long referenceCount = referenceValues.size();
	if(periodic){
		x = moduloInterval(x, domain_min, domain_max);
	}else if(x < domain_min){
		// requested x is outside (left) of the reference domain
		y = outlier_value_left;
		last_requested_reference = 0;
		return;
	}else if(x > domain_max){
		// requested x is outside (right) of the reference domain
		y = outlier_value_right;
		last_requested_reference = referenceCount-1;
		return;
	}
	if(referenceCount == 1){
		y = referenceValues[0];
		last_requested_reference = 0;
		return;
	}
	
	if(referencePoints.empty()){
		long j = floor((x - domain_min)/domainStep); //left border of relevant domain-interval
		j = min(j, referenceCount - 1); //avoid out-of-bounds errors caused by numerical inaccuracies
		last_requested_reference = j;
		if(j == referenceCount - 1){
			y = referenceValues[referenceCount - 1];
		}else{
			//linear interpolation between referenceValues
			const double lambda = (x - (domain_min + j*domainStep))/domainStep;
			convex_combination(referenceValues[j], referenceValues[j+1], lambda, y);
		}
	}else{
		long j;
		if(last_requested_reference<0) last_requested_reference = 0;
		if(referencePoints[last_requested_reference]<=x){
			// search for reference point, towards the right
			for(j=last_requested_reference; j<referenceCount-1; ++j){
				if(referencePoints[j+1] > x) break;
			}
		}else{
			// search for reference point, towards the left
			for(j=last_requested_reference; j>=0; --j){
				if(referencePoints[j] <= x) break;
			}
		}
		if(j>=referenceCount-1){
			y = referenceValues[referenceCount - 1];
			last_requested_reference = referenceCount - 1;
		}else if(j<=0){
			y = referenceValues[0];
			last_requested_reference = 0;			
		}else{
			const double lambda = (x - referencePoints[j])/(referencePoints[j+1] - referencePoints[j]);
			convex_combination(referenceValues[j], referenceValues[j+1], lambda, y);
			last_requested_reference = j;
		}
	}
}


template<class VALUE_TYPE>
VALUE_TYPE LinearInterpolationFunctor<VALUE_TYPE>::operator()(double x) const{
	VALUE_TYPE V;
	getValue(x,V);
	return V;
}





#pragma mark -
#pragma mark Numerical solvers
#pragma mark



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
					long						maxRecordedPoints, 				// (INPUT) if small, some intermediate trajectory points are skipped (i.e. not recorded). This might be useful if accurate integration requires a small time step, but would produce a needlesly large time series. If maxRecordedPoints==2, then only the start and end points are recorded. If maxRecordedPoints==1, then only the end point is recorded.
					long						maxTimeStepRefinements,			// (INPUT) max allowed number of refinements of local time step when encountering invalid states. Only relevant if abortOnInvalidState==false.
					const double				refinement_factor,				// (INPUT) factor by which to refine time steps. Typical values are 2-10.
					const PROGRESS_REPORTER		&reporter,						// (INPUT) callback functor to handle progress reporting during simulation
					const double				runtime_out_seconds,			// (INPUT) max allowed runtime in seconds. If <=0, this is ignored.
					string						&warningMessage){				// (OUTPUT) will be non-empty in case of error, or if non-fatal problems occurred
	COORDINATE currentPoint, candidatePoint, point2, jumpPoint;
	COORDINATE k1, k2, kConsensus;
	double t=startTime, t2, current_dt1, current_dt2;
	long recorded, iterations;
	RequestedDynamics dynamics;
	CrossedBoundary crossedBoundary;
	const double simulationTime 	= endTime-startTime; // target duration of simulation
	const double start_runtime  	= get_thread_monotonic_walltime_seconds();
	const double min_dt 			= dt/pow(refinement_factor,maxTimeStepRefinements);
	
	// keep track of warnings
	bool crossedBoundaryButFixedByReducingTimeStep = false;
	bool crossedBoundaryYesButFixedBruteForce = false;
	bool crossedBoundaryYesButFixedUsingEuler = false;

	//preliminary error checking
	warningMessage = "";
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
	
	const double recordingTimeStep = simulationTime/max(1L,maxRecordedPoints-1);
	model.reserveSpaceForTimeSeries(maxRecordedPoints);
	bool k1AlreadyCalculated = false;
	
	//initialization
	if(!model.getInitialState(t, currentPoint)){
		warningMessage = "Failed to get initial state";
		return false;
	}
	double lastRecordedTime = -INFTY_D;
	if(maxRecordedPoints>1){
		model.registerState(t, currentPoint);
		lastRecordedTime = t;
	}
	
	//run simulation
	//default time step is dt, but temporarily dt might be reduced to a smaller value
	for(recorded=1, current_dt1=current_dt2=dt, iterations=0; t<endTime; /* increment of t handled in loop */ ){
		// at this stage currentPoint is guaranteed to be a valid state
		// t should always correspond to currentPoint
		// t2 should always correspond to point2
		// current_dt1 is used for moving from currentPoint-->point2
		// current_dt2 is used for moving from currentPoint-->candidatePoint (= potential next currentPoint)
		++iterations;
		
		// check runtime-out
		if((runtime_out_seconds>0) && (iterations%100==0) && (get_thread_monotonic_walltime_seconds()-start_runtime>=runtime_out_seconds)){
			warningMessage += string(warningMessage=="" ? "" : "\n") + "Simulation was aborted prematurely after "+makeString(iterations)+" iterations because it reached the maximum allowed processing time.";
			return true;
		}
		
		// reduce time step to not exceed end time, if needed
		current_dt1 = min(min(current_dt1, current_dt2), endTime-t);
				
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
		crossedBoundary = model.checkCrossedDomainBoundaryAndFix(t, currentPoint, current_dt1, point2, true);
		if(crossedBoundary == CrossedBoundaryYesButFixedBruteForce){ crossedBoundaryYesButFixedBruteForce = true; }
		else if(crossedBoundary == CrossedBoundaryYesButFixedByReducingTimeStep){ crossedBoundaryButFixedByReducingTimeStep = true; }
		
		// check if the time step should be reduced (as suggested by the model)
		while((current_dt1>min_dt) && model.checkShouldRefineTimeStep(t, currentPoint, current_dt1, point2)){
			current_dt1 /= refinement_factor;
			point2 = currentPoint + k1*current_dt1;
			t2 = t + current_dt1; // t2 should always correspond to point2
		}
		
		// stage-2 time-step should not be larger than the stage-1 time-step, but also not below min_dt
		current_dt2 = max(min(min_dt,endTime-t), min(current_dt1,current_dt2));
		
		// get dynamics (k2) at point2
		dynamics = model.getRateOfChangeAtState(t2, point2, k2, jumpPoint);
		if(dynamics==RequestedDynamicsForceJumpToState){
			currentPoint = jumpPoint;
			t = t2;
			goto RK2_REGISTER_STATE; 
		}

		// use k1 & k2 to move the currentPoint forward, potentially after refining the time step multiple times
		kConsensus 		= (k1+k2)*0.5;
		candidatePoint 	= currentPoint + kConsensus*current_dt2;
		crossedBoundary = model.checkCrossedDomainBoundaryAndFix(t, currentPoint, current_dt2, candidatePoint, false);
		if(crossedBoundary==CrossedBoundaryYesButFixedBruteForce){
			if(current_dt2>min_dt){
				current_dt2 /= refinement_factor;
				crossedBoundaryButFixedByReducingTimeStep = true;
				continue; // repeat the whole iteration with a smaller time step
			}else{
				// use simple Euler scheme for this time step, since kConsensus throws us out of the boundary at even arbitrarily small time steps
				currentPoint = point2;
				t 			 = t2;
				crossedBoundaryYesButFixedUsingEuler = true;
				goto RK2_REGISTER_STATE;
			}
		}
				
		// check if the time step should be reduced (as suggested by the model)
		if((current_dt2>min_dt) && model.checkShouldRefineTimeStep(t, currentPoint, current_dt2, candidatePoint)){
			current_dt2 /= refinement_factor; 
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
			warningMessage += string(warningMessage=="" ? "" : "\n") + "Simulation reached NaN.";
			return (recorded>1);
		}
		k1AlreadyCalculated = false;
		if((t-lastRecordedTime > recordingTimeStep) && (recorded<maxRecordedPoints-1)){ // don't record if maxRecordedPoints has been almost exhausted (i.e. 1 remaining), since the final state will be recorded outside of the loop
			model.registerState(t, currentPoint);
			++recorded;
			lastRecordedTime = t;
			reporter(recorded, maxRecordedPoints, 1-(endTime-t)/simulationTime);
		}
		
		// initialize counters for next iteration (note: this step may be skipped by a 'continue' statement)
		// only gradually increase the time step, to avoid repeatedly wasting time in refinements
		current_dt1  = min(refinement_factor*current_dt1,dt);
		current_dt2  = min(refinement_factor*current_dt2,dt);
	}
	
	// register final state if needed
	if(t>lastRecordedTime){
		model.registerState(t, currentPoint);
	}
	
	// construct warnings message if needed
	if(crossedBoundaryYesButFixedUsingEuler) warningMessage += string(warningMessage=="" ? "" : "\n") + "Simulation crossed domain boundary and was fixed by locally using forward Euler.";
	if(crossedBoundaryYesButFixedBruteForce) warningMessage += string(warningMessage=="" ? "" : "\n") + "Simulation crossed domain boundary and was fixed by brute-force adjusting the trajectory.";
	if(crossedBoundaryButFixedByReducingTimeStep) warningMessage += string(warningMessage=="" ? "" : "\n") + "Simulation crossed domain boundary and was fixed by locally reducing the time step.";

	return true;
}







// RungeKutta2 integrator, handling:
//    1. Suggestions by the model for temporarily refined time steps
//    2. Situations where the domain boundary is crossed during an iteration (in which case the time step is temporarily decreased as needed in order to not overshoot)
//    3. Requests by the model to jump discontinuously to another state
// requested times might reverse (i.e. times at which model dynamics are requested might not be strictly increasing), but recorded time series will be strictly forward in time
// This is similar to RungeKutta2 above, the difference being that here points are recorded when deemed appropriate based on the difference from the previous time point (or equivalently, the local rate of change)
// Hence, instead of specifying the max number of recorded time points, the caller specifies the minRecordingTimeStep and the recordingRelValueStep
template<class COORDINATE, class MODEL, class PROGRESS_REPORTER>
bool RungeKutta2(	double						startTime,						// (INPUT) simulation start time
					double 						endTime, 						// (INPUT) simulation end time
					double 						dt, 							// (INPUT) default integration time step
					MODEL 						&model,							// (INPUT/OUTPUT) object defining model dynamics, also handling time series storage
					double						minRecordingTimeStep, 			// (INPUT) the minimum time difference between subsequent recordings
					double						recordingRelValueStep, 			// (INPUT) the minimum relative value difference to the previous recording, before triggering a new recording. This is sub-ordinate to minRecordingTimeStep, i.e. minRecordingTimeStep will never be violated (except perhaps at the end, because the last time point is always recorded). The "relative difference" between two states is defined by the model, not the solver. Typically this will be the maximum relative difference between components.
					long						maxTimeStepRefinements,			// (INPUT) max allowed number of refinements of local time step when encountering invalid states. Only relevant if abortOnInvalidState==false.
					const double				refinement_factor,				// (INPUT) factor by which to refine time steps. Typical values are 2-10.
					const PROGRESS_REPORTER		&reporter,						// (INPUT) callback functor to handle progress reporting during simulation
					const double				runtime_out_seconds,			// (INPUT) max allowed runtime in seconds. If <=0, this is ignored.
					string						&warningMessage){				// (OUTPUT) will be non-empty in case of error, or if non-fatal problems occurred
	COORDINATE currentPoint, candidatePoint, point2, jumpPoint;
	COORDINATE k1, k2, kConsensus;
	double t, t2, current_dt1, current_dt2;
	long recorded, iterations;
	RequestedDynamics dynamics;
	CrossedBoundary crossedBoundary;
	const double simulationTime 	= endTime-startTime; // target duration of simulation
	const double start_runtime  	= get_thread_monotonic_walltime_seconds();
	const double min_dt 			= dt/pow(refinement_factor,maxTimeStepRefinements);
	
	// keep track of warnings
	bool crossedBoundaryButFixedByReducingTimeStep = false;
	bool crossedBoundaryYesButFixedBruteForce = false;
	bool crossedBoundaryYesButFixedUsingEuler = false;

	//preliminary error checking
	warningMessage = "";
	if(dt<simulationTime*RELATIVE_EPSILON){
		warningMessage = "Time step too small";
		return false;
	}
	if(simulationTime < dt){
		warningMessage = "Time step exceeds simulation time";
		return false;
	}
	if(recordingRelValueStep < 0){
		warningMessage = "recordingRelValueStep is negative";
		return false;
	}
	
	model.reserveSpaceForTimeSeries(simulationTime,minRecordingTimeStep,recordingRelValueStep);
	bool k1AlreadyCalculated = false;
	
	//initialization
	t = startTime;
	if(!model.getInitialState(t, currentPoint)){
		warningMessage = "Failed to get initial state";
		return false;
	}
	model.registerState(t, currentPoint);
	COORDINATE lastRecordedPoint = currentPoint;
	double lastRecordedTime = t;
	
	//run simulation
	//default time step is dt, but temporarily dt might be reduced to a smaller value
	for(recorded=1, current_dt1=current_dt2=dt, iterations=0; t<endTime; /* increment of t handled in loop */ ){
		// at this stage currentPoint is guaranteed to be a valid state
		// t should always correspond to currentPoint
		// t2 should always correspond to point2
		// lastRecordedTime should always correspond to lastRecordedPoint
		// current_dt1 is used for moving from currentPoint-->point2
		// current_dt2 is used for moving from currentPoint-->candidatePoint (= potential next currentPoint)
		++iterations;
		
		// check runtime-out
		if((runtime_out_seconds>0) && (iterations%100==0) && (get_thread_monotonic_walltime_seconds()-start_runtime>=runtime_out_seconds)){
			warningMessage += string(warningMessage=="" ? "" : "\n") + "Simulation was aborted prematurely after "+makeString(iterations)+" iterations because it reached the maximum allowed processing time.";
			return true;
		}
		
		// reduce time step to not exceed end time, if needed
		current_dt1 = min(min(current_dt1, current_dt2), endTime-t);
				
		// get dynamics (k1) at currentPoint if needed
		if(!k1AlreadyCalculated){
			dynamics = model.getRateOfChangeAtState(t, currentPoint, k1, jumpPoint);
			if(dynamics==RequestedDynamicsForceJumpToState){
				currentPoint = jumpPoint;
				t 			 = t + current_dt1;
				goto ARK2_REGISTER_STATE; 
			}
			k1AlreadyCalculated = true;
		}
		
		// get point2 (in forward k1-direction)
		linear_combination(1.0,currentPoint,current_dt1,k1,point2); // point2 = currentPoint + k1*current_dt1;
		t2 = t + current_dt1; // t2 should always correspond to point2
		
		// check if point2 crossed the boundary and move point2 backward as much as needed
		// point2 and current_dt1 will be adjusted if this routine returns true
		crossedBoundary = model.checkCrossedDomainBoundaryAndFix(t, currentPoint, current_dt1, point2, true);
		if(crossedBoundary == CrossedBoundaryYesButFixedBruteForce){ crossedBoundaryYesButFixedBruteForce = true; }
		else if(crossedBoundary == CrossedBoundaryYesButFixedByReducingTimeStep){ crossedBoundaryButFixedByReducingTimeStep = true; }
		
		// check if the time step should be reduced (as suggested by the model)
		while((current_dt1>min_dt) && model.checkShouldRefineTimeStep(t, currentPoint, current_dt1, point2)){
			current_dt1 /= refinement_factor;
			linear_combination(1.0,currentPoint,current_dt1,k1,point2); // point2 = currentPoint + k1*current_dt1;
			t2 = t + current_dt1; // t2 should always correspond to point2
		}
		
		// stage-2 time-step should not be larger than the stage-1 time-step, but also not below min_dt
		current_dt2 = max(min(min_dt,endTime-t), min(current_dt1,current_dt2));
		
		// get dynamics (k2) at point2
		dynamics = model.getRateOfChangeAtState(t2, point2, k2, jumpPoint);
		if(dynamics==RequestedDynamicsForceJumpToState){
			currentPoint = jumpPoint;
			t			 = t2;
			goto ARK2_REGISTER_STATE; 
		}

		// use k1 & k2 to move the currentPoint forward, potentially after refining the time step multiple times
		linear_combination(0.5, k1, 0.5, k2, kConsensus); // kConsensus = 0.5*(k1+k2)
		linear_combination(1.0, currentPoint, current_dt2, kConsensus, candidatePoint); // candidatePoint 	= currentPoint + kConsensus*current_dt2
		crossedBoundary = model.checkCrossedDomainBoundaryAndFix(t, currentPoint, current_dt2, candidatePoint, false);
		if(crossedBoundary==CrossedBoundaryYesButFixedBruteForce){
			if(current_dt2>min_dt){
				current_dt2 /= refinement_factor;
				crossedBoundaryButFixedByReducingTimeStep = true;
				continue; // repeat the whole iteration with a smaller time step
			}else{
				// use simple Euler scheme for this time step, since kConsensus throws us out of the boundary at even arbitrarily small time steps
				currentPoint = point2;
				t	 		 = t2;
				crossedBoundaryYesButFixedUsingEuler = true;
				goto ARK2_REGISTER_STATE;
			}
		}
				
		// check if the time step should be reduced (as suggested by the model)
		if((current_dt2>min_dt) && model.checkShouldRefineTimeStep(t, currentPoint, current_dt2, candidatePoint)){
			current_dt2 /= refinement_factor; 
			continue; // repeat the last part with a smaller time step
		}

		// seems everything worked out as normal, so set the new currentPoint
		currentPoint = candidatePoint;
		t 			 = t + current_dt2; // t should always correspond to currentPoint
				
		// register new state if needed	
		ARK2_REGISTER_STATE:
				
		// check sanity (some ODEs may explode)		
		if(std::isnan(t) || model.stateIsNaN(currentPoint)){
			warningMessage += string(warningMessage=="" ? "" : "\n") + "Simulation reached NaN.";
			return (recorded>1);
		}
		k1AlreadyCalculated = false;
		if((t-lastRecordedTime >= minRecordingTimeStep) && (model.getRelativeChange(lastRecordedPoint,currentPoint)>recordingRelValueStep)){ // don't record if maxRecordedPoints has been almost exhausted (i.e. 1 remaining), since the final state will be recorded outside of the loop			
			model.registerState(t, currentPoint);
			++recorded;
			lastRecordedTime  = t;
			lastRecordedPoint = currentPoint;
			reporter((t-startTime), simulationTime);
		}
		
		// initialize counters for next iteration (note: this step may be skipped by a 'continue' statement)
		// only gradually increase the time step, to avoid repeatedly wasting time in refinements
		current_dt1  = min(refinement_factor*current_dt1,dt);
		current_dt2  = min(refinement_factor*current_dt2,dt);
	}
	
	// register final state if needed
	if(t>lastRecordedTime){
		model.registerState(t, currentPoint);
	}
	
	// construct warnings message if needed
	if(crossedBoundaryYesButFixedUsingEuler) warningMessage += string(warningMessage=="" ? "" : "\n") + "Simulation crossed domain boundary and was fixed by locally using forward Euler.";
	if(crossedBoundaryYesButFixedBruteForce) warningMessage += string(warningMessage=="" ? "" : "\n") + "Simulation crossed domain boundary and was fixed by brute-force adjusting the trajectory.";
	if(crossedBoundaryButFixedByReducingTimeStep) warningMessage += string(warningMessage=="" ? "" : "\n") + "Simulation crossed domain boundary and was fixed by locally reducing the time step.";

	return true;
}






// Two-step Runge-Kutta solver for ODEs where the "scale" and the "shape" (rescaled-version) of the trajectory should be simulated instead.
// Consider the following ODE:
//    dX/dt = f(t,X)
// Here, instead of directly solving the ODE for X(t), we solve two auxiliary ODEs for S(t) and Y(t), where S is the "scale" (in some abstract sense) of X, and Y is a rescaled version of X (i.e. the "shape" of X)
// For example, if X is a classical vector, then S may be defined as log(mean(X)) and Y = X/S.
// Hence, any current state X is represented by a tuple (S,Y), where S is scalar and Y is of similar nature as X (the exact format of Y is handled transparently by the model, not the solver).
// A separation between "shape" Y and "scale" S may be especially useful if X(t) decays to zero exponentially, resulting in numerical underflow, or explodes to infinity exponentially, resulting in numerical overflow.
//
// The solver can handle:
//    1. Suggestions by the model for temporarily refined time steps
//    2. Situations where the domain boundary is crossed during an iteration (in which case the time step is temporarily decreased as needed in order to not overshoot)
//    3. Requests by the model to jump discontinuously to another state
//
// The shape Y must be of type COORDINATE, and the scale S must be of type double.
// The model must be able to provide initial conditions and the dynamics for S and Y at each time point, as well as storing of the computed trajectory (defined by S and Y).
// Of course, the model may internally store the trajectory X instead of S and Y, but the integrator will always compute and provide S and Y.
//
// Requested times might reverse (i.e. times at which model dynamics are requested might not be strictly increasing), but recorded time series will be strictly forward in time
template<class COORDINATE, class MODEL, class PROGRESS_REPORTER>
bool ScaledRungeKutta(	const double				startTime,						// (INPUT) simulation start time
						const double				endTime, 						// (INPUT) simulation end time
						double 						dt, 							// (INPUT) default integration time step
						MODEL 						&model,							// (INPUT/OUTPUT) object defining model dynamics (rates of change of S & Y), and handling time series storage
						long						maxRecordedPoints, 				// (INPUT) if small, some intermediate trajectory points are skipped (i.e. not recorded). This might be useful if accurate integration requires a small time step, but would produce a needlesly large time series. If maxRecordedPoints==2, then only the start and end points are recorded. If maxRecordedPoints==1, then only the end point is recorded.
						long						maxTimeStepRefinements,			// (INPUT) max allowed number of refinements of local time step when encountering invalid states. Only relevant if abortOnInvalidState==false.
						const PROGRESS_REPORTER		&reporter,						// (INPUT) callback functor to handle progress reporting during simulation
						const double				runtime_out_seconds,			// (INPUT) max allowed runtime in seconds. If <=0, this is ignored.
						string						&warningMessage){				// (OUTPUT) will be non-empty in case of error, or if non-fatal problems occurred
	COORDINATE currentY, candidateY, Y2, jumpY;
	double currentS, candidateS, S2, jumpS;
	COORDINATE kY1, kY2, kYConsensus;
	double kS1, kS2, kSConsensus;
	double t=startTime, t2, current_dt1, current_dt2;
	long recorded, refinements, iterations;
	RequestedDynamics dynamics;
	CrossedBoundary crossedBoundary;
	const double simulationTime = endTime-startTime; // target duration of simulation
	const double start_runtime  = get_thread_monotonic_walltime_seconds();
	
	// keep track of warnings
	bool crossedBoundaryButFixedByReducingTimeStep = false;
	bool crossedBoundaryYesButFixedBruteForce = false;
	bool crossedBoundaryYesButFixedUsingEuler = false;
	
	//preliminary error checking
	warningMessage = "";
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
	model.reserveSpaceForScaledTimeSeries(maxRecordedPoints);
	bool k1AlreadyCalculated = false;
	
	// get initial state (shape and scale)
	if(!model.getScaledInitialState(t, currentY, currentS)){
		warningMessage = "Failed to get initial state";
		return false;
	}
	
	// record initial state
	double lastRecordedTime = -INFTY_D;
	if(maxRecordedPoints>1){
		model.registerScaledState(t, currentY, currentS);
		lastRecordedTime = t;
	}
	
	//run simulation
	//default time step is dt, but temporarily dt might be reduced to a smaller value
	for(recorded=1, current_dt1=current_dt2=dt, refinements=0, iterations=0; t<endTime; /* increment of t handled in loop */ ){
		// at this stage currentY is guaranteed to be a valid state
		// t should always correspond to currentY
		// t2 should always correspond to Y2
		// current_dt1 is used for moving from currentY-->Y2
		// current_dt2 is used for moving from currentY-->candidateY (= potential next currentY)
		++iterations;
		
		// check runtime-out
		if((runtime_out_seconds>0) && (iterations%100==0) && (get_thread_monotonic_walltime_seconds()-start_runtime>=runtime_out_seconds)){
			warningMessage += string(warningMessage=="" ? "" : "\n") + "Simulation was aborted prematurely after "+makeString(iterations)+"  because it reached the maximum allowed processing time.";
			return true;
		}
		
		// reduce time step to not exceed end time, if needed
		current_dt1 = min(min(current_dt1, current_dt2), endTime-t);
				
		// get dynamics (k1) at the current state if needed
		// the dynamics will be the rate of change for the "shape" Y and for the "scale" S, i.e. symbolically k1 = (dY/dt, dS/dt)
		if(!k1AlreadyCalculated){
			dynamics = model.getRateOfChangeAtScaledState(t, currentY, currentS, kY1, kS1, jumpY, jumpS);
			if(dynamics==RequestedDynamicsForceJumpToState){
				currentY = jumpY;
				currentS = jumpS;
				t = t + current_dt1;
				goto SRK2_REGISTER_STATE; 
			}
			k1AlreadyCalculated = true;
		}
		
		// get Y2 & S2 (in forward k1-direction)
		Y2  = currentY + kY1*current_dt1;
		S2  = currentS + kS1*current_dt1;
		t2  = t + current_dt1; // t2 should always correspond to Y2
		
		// check if (Y2,S2) crossed the boundary and move (Y2,S2) backward as much as needed
		// (Y2,S2) and current_dt1 will be adjusted if this routine returns true
		crossedBoundary = model.checkCrossedDomainBoundaryAndFix(t, currentY, currentS, current_dt1, Y2, S2);
		if(crossedBoundary == CrossedBoundaryYesButFixedBruteForce){ crossedBoundaryYesButFixedBruteForce = true; }
		else if(crossedBoundary == CrossedBoundaryYesButFixedByReducingTimeStep){ crossedBoundaryButFixedByReducingTimeStep = true; }
				
		// check if the time step should be reduced (as suggested by the model)
		while((refinements<maxTimeStepRefinements) && model.checkShouldRefineTimeStep(t, currentY, currentS, current_dt1, Y2, S2)){
			current_dt1 /= 2;
			++refinements;
			Y2 = currentY + kY1*current_dt1;
			S2 = currentS + kS1*current_dt1;
			t2 = t + current_dt1; // t2 should always correspond to (Y2,S2)
		}		
		
		// get dynamics (k2) at (Y2,S2)
		dynamics = model.getRateOfChangeAtScaledState(t2, Y2, S2, kY2, kS2, jumpY, jumpS);
		if(dynamics==RequestedDynamicsForceJumpToState){
			currentY = jumpY;
			currentS = jumpS;
			t = t2;
			goto SRK2_REGISTER_STATE; 
		}

		// use k1 & k2 to move the currentY forward, potentially after refining the time step multiple times
		kYConsensus = (kY1+kY2)*0.5;
		kSConsensus = (kS1+kS2)*0.5;
		candidateY  = currentY + kYConsensus*current_dt2;
		candidateS  = currentS + kSConsensus*current_dt2;
		crossedBoundary = model.checkCrossedDomainBoundaryAndFix(t, currentY, currentS, current_dt2, candidateY, candidateS);
		if(crossedBoundary==CrossedBoundaryYesButFixedBruteForce){
			if(refinements<maxTimeStepRefinements){
				++refinements; 
				current_dt2 /= 2;
				crossedBoundaryButFixedByReducingTimeStep = true;
				continue; // repeat the whole iteration with a smaller time step
			}else{
				// use simple Euler scheme for this time step, since kConsensus throws us out of the boundary at even arbitrarily small time steps
				currentY = Y2;
				currentS = S2;
				t = t2;
				crossedBoundaryYesButFixedUsingEuler = true;
				goto SRK2_REGISTER_STATE;
			}
		}
				
		// check if the time step should be reduced (as suggested by the model)
		if((refinements<maxTimeStepRefinements) && model.checkShouldRefineTimeStep(t, currentY, currentS, current_dt2, candidateY, candidateS)){
			++refinements; 
			current_dt2 /= 2; 
			continue; // repeat the last part with a smaller time step
		}
		
		// seems everything worked out as normal, so set the new (currentY, currentS)
		currentY = candidateY;
		currentS = candidateS;
		t = t + current_dt2; // t should always correspond to (currentY, currentS)
		goto SRK2_REGISTER_STATE;
		
		// check and register new state if needed	
		SRK2_REGISTER_STATE:
		// check sanity (some ODEs may explode)		
		if(std::isnan(t) || model.scaledStateIsNaN(currentY,currentS)){
			warningMessage += string(warningMessage=="" ? "" : "\n") + "Simulation reached NaN.";
			return (recorded>1);
		}
		k1AlreadyCalculated = false;
		if((t-lastRecordedTime > recordingTimeStep) && (recorded<maxRecordedPoints-1)){ // don't record if maxRecordedPoints has been almost exhausted (i.e. 1 remaining), since the final state will be recorded outside of the loop
			model.registerScaledState(t, currentY, currentS);
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
		model.registerScaledState(t, currentY, currentS);
	}
	
	// construct warnings message if needed
	if(crossedBoundaryYesButFixedUsingEuler) warningMessage += string(warningMessage=="" ? "" : "\n") + "Simulation crossed domain boundary and was fixed by locally using forward Euler.";
	if(crossedBoundaryYesButFixedBruteForce) warningMessage += string(warningMessage=="" ? "" : "\n") + "Simulation crossed domain boundary and was fixed by brute-force adjusting the trajectory.";
	if(crossedBoundaryButFixedByReducingTimeStep) warningMessage += string(warningMessage=="" ? "" : "\n") + "Simulation crossed domain boundary and was fixed by locally reducing the time step.";

	return true;
}





// Two-step Runge-Kutta solver for linear ODEs where the log-scale and a scaled-version of the trajectory should be simulated instead.
// Consider the following linear ODE:
//    dX/dt = A(t)*X(t),
// where the dynamic variable X(t) is a vector or a matrix of size NR x NC, and where A(t) is an explicitly given matrix of size NR x NR.
// Here, the variable X(t) is split into two auxiliary variables, the log-scale S (S=log(mean(X))) and the "shape" Y=X/exp(S).
// Hence, any current state X is represented by a tuple (S,Y), where S is scalar and Y is of similar nature as X.
// A separation between "shape" Y and "scale" S may be especially useful if X(t) decays to zero exponentially, resulting in numerical underflow, or explodes to infinity exponentially, resulting in numerical overflow.
//
// The solver can handle:
//    1. Suggestions by the model for temporarily refined time steps
//    2. Situations where the domain boundary is crossed during an iteration (in which case the time step is temporarily decreased as needed in order to not overshoot)
//
// The shape Y must be of type COORDINATE (either a vector or matrix in row-major format), and the scale S will be of type double.
// The model must be able to provide initial conditions for X, the dynamics for X in the form of A(t) at each time point, as well as storing of the computed trajectory (defined by S and Y).
// Of course, the model may internally store the trajectory X instead of S and Y, but the integrator will always compute and provide S and Y.
//
// Requested times might reverse (i.e. times at which model dynamics are requested might not be strictly increasing), but recorded time series will be strictly forward in time
template<class COORDINATE, class MODEL, class PROGRESS_REPORTER>
bool LinearScaledRungeKutta(const long					NR,								// (INPUT) number of rows in A and X
							const long					NC,								// (INPUT) number of columns in X
							const double				startTime,						// (INPUT) simulation start time
							const double				endTime, 						// (INPUT) simulation end time
							const double 				dt, 							// (INPUT) default integration time step
							MODEL 						&model,							// (INPUT/OUTPUT) object defining model dynamics (via the matrix A), and handling time series storage
							const long					maxRecordedPoints, 				// (INPUT) if small, some intermediate trajectory points are skipped (i.e. not recorded). This might be useful if accurate integration requires a small time step, but would produce a needlesly large time series. If maxRecordedPoints==2, then only the start and end points are recorded. If maxRecordedPoints==1, then only the end point is recorded.
							const long					maxTimeStepRefinements,			// (INPUT) max allowed number of refinements of local time step when encountering invalid states. Only relevant if abortOnInvalidState==false.
							const double				refinement_factor,				// (INPUT) factor by which to refine time steps. Typical values are 2-10.
							const long 					max_exp_order,					// (INPUT) maximum polynomial order for approximating exponentials (short-term propagators). Typical values are 2-5
							const PROGRESS_REPORTER		&reporter,						// (INPUT) callback functor to handle progress reporting during simulation
							const double				runtime_out_seconds,			// (INPUT) max allowed runtime in seconds. If <=0, this is ignored.
							string						&warningMessage){				// (OUTPUT) will be non-empty in case of error, or if non-fatal problems occurred
	COORDINATE currentY, candidateY, Y2, initialX;
	double currentS, candidateS, S2, rescaling;
	std::vector<double> A1, A2, Aconsensus, scratch1, scratch2;
	double t=0, t2, current_dt1, current_dt2, original_dt2;
	long recorded, iterations;
	CrossedBoundary crossedBoundary;
	const double DeltaTime 			= endTime-startTime; // target duration of simulation
	const double start_runtime  	= get_thread_monotonic_walltime_seconds();
	const double min_dt 			= dt/pow(refinement_factor,maxTimeStepRefinements);
	
	// note that we internally measure time in relative time "t", i.e. starting at 0 until DeltaTime
	// this is needed to avoid floating point rounding errors when DeltaT is much smaller than startTime
	// whenever we "speak" to the model (e.g. requesting dynamics or recording a trajectory point), we transform it to real time
	
	
	// keep track of warnings
	bool crossedBoundaryButFixedByReducingTimeStep = false;
	bool crossedBoundaryYesButFixedBruteForce = false;
	bool crossedBoundaryYesButFixedUsingEuler = false;
		
	//preliminary error checking
	warningMessage = "";
	if(dt<DeltaTime*RELATIVE_EPSILON){
		warningMessage = "Time step too small";
		return false;
	}
	if(DeltaTime < dt){
		warningMessage = "Time step exceeds simulation time";
		return false;
	}
	if(maxRecordedPoints < 1){
		warningMessage = "Requested zero recorded points";
		return false;
	}
	
	const double recordingTimeStep = DeltaTime/max(1L,maxRecordedPoints-1);
	model.reserveSpaceForScaledTimeSeries(maxRecordedPoints);
	bool k1AlreadyCalculated = false;
	
	// get initial state (shape and scale)
	if(!model.getInitialState(startTime+t, initialX)){
		warningMessage = "Failed to get initial state";
		return false;
	}
	
	// create scaled version of initial state
	const double initial_mean = vector_mean(initialX);
	if(initial_mean<=0){
		warningMessage = "Initial state as negative or zero scale";
		return false;
	}
	currentY = initialX/initial_mean;
	currentS = log(initial_mean);
	
	// record initial state
	double lastRecordedTime = -INFTY_D;
	if(maxRecordedPoints>1){
		model.registerScaledState(startTime+t, currentY, currentS);
		lastRecordedTime = t;
	}
		
	//run simulation
	//default time step is dt, but temporarily dt might be reduced to a smaller value
	for(recorded=1, current_dt1=current_dt2=dt, iterations=0; t<DeltaTime; /* increment of t handled in loop */ ){
		// at this stage currentY is guaranteed to be a valid state
		// t should always correspond to currentY
		// t2 should always correspond to Y2
		// current_dt1 is used for moving from currentY-->Y2
		// current_dt2 is used for moving from currentY-->candidateY (= potential next currentY)
		++iterations;
		
		// check runtime-out
		if((runtime_out_seconds>0) && (iterations%100==0) && (get_thread_monotonic_walltime_seconds()-start_runtime>=runtime_out_seconds)){
			warningMessage += string(warningMessage=="" ? "" : "; ") + "Simulation was aborted prematurely after "+makeString(iterations)+" iterations because it reached the maximum allowed processing time.";
			return true;
		}
		
		// reduce time step dt1 to not exceed end time nor dt2, if needed
		current_dt1 = min(min(current_dt1,current_dt2), DeltaTime-t);
				
		// get dynamics at the current time t, in the form of the matrix A(t)
		if(!k1AlreadyCalculated){
			model.getLinearDynamics(startTime+t, A1);
			k1AlreadyCalculated = true;
		}
		
		// get Y2 by applying the dynamics A: Y2 = exp(current_dt1*A)*currentY
		// don't update S yet, do so at the end of this iteration
		apply_approximate_matrix_exponential(NR, NC, A1, current_dt1, currentY, max_exp_order, scratch1, scratch2, Y2);
		S2  = currentS; // don't change scale until later
		t2  = t + current_dt1; // t2 should always correspond to Y2

		// check if Y2 crossed the boundary and move Y2 backward as much as needed
		// (Y2,S2) and current_dt1 will be adjusted if this routine returns true
		crossedBoundary = model.checkCrossedDomainBoundaryAndFix(startTime+t, currentY, currentS, current_dt1, Y2, S2, true);
		if(crossedBoundary == CrossedBoundaryYesButFixedBruteForce){ crossedBoundaryYesButFixedBruteForce = true; }
		else if(crossedBoundary == CrossedBoundaryYesButFixedByReducingTimeStep){ crossedBoundaryButFixedByReducingTimeStep = true; }
						
		// check if the time step should be reduced (as suggested by the model)
		while((current_dt1>min_dt) && model.checkShouldRefineTimeStep(startTime+t, currentY, currentS, current_dt1, Y2, S2)){
			// recalculate currentY --> Y2 with a smaller time step
			current_dt1 /= refinement_factor;
			apply_approximate_matrix_exponential(NR, NC, A1, current_dt1, currentY, max_exp_order, scratch1, scratch2, Y2);
			S2 = currentS;
			t2 = t + current_dt1; // t2 should always correspond to (Y2,S2)
		}
		
		// stage-2 time-step should not be larger than the stage-1 time-step, but also not below min_dt
		current_dt2 = max(min(min_dt,endTime-t), min(current_dt1,current_dt2));
		
		// get dynamics (matrix A) at t2
		model.getLinearDynamics(startTime+t2, A2);

		// use A1 & A2 to move the currentY forward (t --> t+current_dt2), and check if time step was sufficiently small
		linear_combination(0.5, A1, 0.5, A2, Aconsensus);
		apply_approximate_matrix_exponential(NR, NC, Aconsensus, current_dt2, currentY, max_exp_order, scratch1, scratch2, candidateY);
		candidateS = currentS; // increment scale later on

		// check if we crossed the domain boundary, and correct if needed		
		original_dt2 = current_dt2;
		crossedBoundary = model.checkCrossedDomainBoundaryAndFix(startTime+t, currentY, currentS, current_dt2, candidateY, candidateS, true);
		if(crossedBoundary==CrossedBoundaryYesButFixedBruteForce){
			if(current_dt2>min_dt){
				current_dt2 /= refinement_factor;
				crossedBoundaryButFixedByReducingTimeStep = true;
				continue; // repeat the whole iteration with a smaller time step
			}else{
				// use simple Euler scheme for this time step, since kConsensus throws us out of the boundary at even arbitrarily small time steps
				rescaling 	= vector_mean(Y2);
				currentY 	= Y2/rescaling;
				currentS 	= S2+log(rescaling); // since we are rescaling Y2, correct for this rescaling in currentS		
				t 			= t2;
				crossedBoundaryYesButFixedUsingEuler = true;
				goto SRK2_REGISTER_STATE;
			}
		}else if(crossedBoundary == CrossedBoundaryYesButFixedByReducingTimeStep){
			crossedBoundaryButFixedByReducingTimeStep = true;
		}

		// check if the time step should be reduced (as suggested by the model)
		if((current_dt2>min_dt) && model.checkShouldRefineTimeStep(startTime+t, currentY, currentS, current_dt2, candidateY, candidateS)){
			current_dt2 /= refinement_factor;
			continue; // repeat this iteration with a smaller time step, forget about the current candidateY
		}
		
		// seems everything worked out as normal, so set the new (currentY, currentS)
		rescaling = vector_mean(candidateY);
		currentY = candidateY/rescaling;
		currentS = candidateS + log(rescaling); // since we are rescaling currentY, correct the scale for this rescaling
		t = t + current_dt2; // t should always correspond to (currentY, currentS)
		goto SRK2_REGISTER_STATE;
		
		// check and register new state if needed	
		SRK2_REGISTER_STATE:
		// check sanity (some ODEs may explode)		
		if(std::isnan(t) || model.scaledStateIsNaN(currentY,currentS)){
			warningMessage += string(warningMessage=="" ? "" : "; ") + "Simulation reached NaN.";
			return (recorded>1);
		}
		k1AlreadyCalculated = false;
		if((t-lastRecordedTime > recordingTimeStep) && (recorded<maxRecordedPoints-1)){ // don't record if maxRecordedPoints has been almost exhausted (i.e. 1 remaining), since the final state will be recorded outside of the loop
			model.registerScaledState(startTime+t, currentY, currentS);
			++recorded;
			lastRecordedTime = t;
			reporter(recorded, maxRecordedPoints, 1-(DeltaTime-t)/DeltaTime);
		}
		
		// initialize counters for next iteration (note: this step may be skipped by a 'continue' statement)
		// only gradually increase the time step, to avoid repeatedly wasting time in refinements
		current_dt1  = min(refinement_factor*current_dt1,dt);
		current_dt2  = min(refinement_factor*current_dt2,dt);
	}
	
	// register final state if needed
	if(t>lastRecordedTime){
		model.registerScaledState(startTime+t, currentY, currentS);
	}

	// construct warnings message if needed
	if(crossedBoundaryYesButFixedUsingEuler) warningMessage += string(warningMessage=="" ? "" : "; ") + "Simulation crossed domain boundary and was fixed by locally using forward Euler.";
	if(crossedBoundaryYesButFixedBruteForce) warningMessage += string(warningMessage=="" ? "" : "; ") + "Simulation crossed domain boundary and was fixed by brute-force adjusting the trajectory.";
	if(crossedBoundaryButFixedByReducingTimeStep) warningMessage += string(warningMessage=="" ? "" : "; ") + "Simulation crossed domain boundary and was fixed by locally reducing the time step.";
	
	return true;
}






// Rosenbrock-Euler order-2 ODE solver for equations of the form:
//    dX/dt = f(t,X) = A(t)*X(t) + N(t,X(t)),
// where the dynamic variable X(t) is a vector or a matrix of size NR x NC, and where A(t) ("linearity") is an explicitly given matrix of size NR x NR and N(t,X(t)) is the non-linear residual.
// This is of the simplest exponential integrators. Exponential integrators are particularly suited for stiff problems, when the stiffness is caused by the linear dynamics.
// At each iteration, the integration proceeds with the following step:
//	 X(t+eps) = X(t) + [exp(eps*A(t)) - Id] * A(t)^{-1} * f(t,X(t))
// where:
//   f(t,X(t)) = [A(t)*X(t) + N(t,X(t))]
// The exponential * inverse_A term can be approximated as:
//   [exp(eps*A(t)) - Id] * A(t)^{-1} = sum_{k=0}^n (eps/(k+1)) * (eps*A(t))^k/k!
//  Here, we shall call this the Rosenbrock-Euler propagator.
//
// The solver can handle:
//    1. Suggestions by the model for temporarily refined time steps
//    2. Situations where the domain boundary is crossed during an iteration (in which case the time step is temporarily decreased as needed in order to not overshoot)
//
// The model must be able to provide initial conditions for X, the linearity A(t) and the nonlinearity N(t,X(t)) at each time point, as well as storing of the computed trajectory
// For a quick derivation of the Rosenbrock-Euler scheme see: 
//   Chen et al (2018). Exponential Rosenbrock-Euler integrators for elastodynamic simulation. IEEE transactions on visualization and computer graphics. 24:2702-2713
//
// Requested times might reverse (i.e. times at which model dynamics are requested might not be strictly increasing), but recorded time series will be strictly forward in time
template<class COORDINATE, class MODEL, class PROGRESS_REPORTER>
bool RosenbrockEuler(	const long					NR,								// (INPUT) number of rows in A and X
						const long					NC,								// (INPUT) number of columns in X
						const double				startTime,						// (INPUT) simulation start time
						const double				endTime, 						// (INPUT) simulation end time
						const double 				dt, 							// (INPUT) default integration time step
						MODEL 						&model,							// (INPUT/OUTPUT) object defining model dynamics (via the matrix A), and handling time series storage
						const long					maxRecordedPoints, 				// (INPUT) if small, some intermediate trajectory points are skipped (i.e. not recorded). This might be useful if accurate integration requires a small time step, but would produce a needlesly large time series. If maxRecordedPoints==2, then only the start and end points are recorded. If maxRecordedPoints==1, then only the end point is recorded.
						const long					maxTimeStepRefinements,			// (INPUT) max allowed number of refinements of local time step when encountering invalid states. Only relevant if abortOnInvalidState==false.
						const long 					max_exp_order,					// (INPUT) maximum polynomial order for approximating exponentials (short-term propagators). Typical values are 2-5
						const PROGRESS_REPORTER		&reporter,						// (INPUT) callback functor to handle progress reporting during simulation
						const double				runtime_out_seconds,			// (INPUT) max allowed runtime in seconds. If <=0, this is ignored.
						string						&warningMessage){				// (OUTPUT) will be non-empty in case of error, or if non-fatal problems occurred
	COORDINATE currentX, candidateX, X2, nonlinearity, rate;
	std::vector<double> linearity, scratch1, scratch2, default_RE_propagator;
	double t=startTime, t2, current_dt, candidate_dt;
	long recorded;
	CrossedBoundary crossedBoundary;
	const double simulationTime = endTime-startTime; // target duration of simulation
	const double start_runtime  = get_thread_monotonic_walltime_seconds();
	const double min_dt = dt/pow(2.0,maxTimeStepRefinements);
	const bool constant_linearity = model.linearDynamicsAreConstant();
	
	// keep track of warnings
	bool crossedBoundaryButFixedByReducingTimeStep = false;
	bool crossedBoundaryYesButFixedBruteForce = false;
	bool crossedBoundaryYesButFixedUsingEuler = false;
		
	//preliminary error checking
	warningMessage = "";
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
	
	const double recordingTimeStep = simulationTime/max(1L,maxRecordedPoints-1);
	model.reserveSpaceForTimeSeries(maxRecordedPoints);
	bool k1AlreadyCalculated = false;
	
	// get initial state (shape and scale)
	if(!model.getInitialState(t, currentX)){
		warningMessage = "Failed to get initial state";
		return false;
	}
		
	// record initial state
	double lastRecordedTime = -INFTY_D;
	if(maxRecordedPoints>1){
		model.registerState(t, currentX);
		lastRecordedTime = t;
	}
	
	if(constant_linearity){
		// since linearity is always the same, precompute the default Rosenbrock-Euler propagator (i.e. for the default time step)
		model.getLinearAndNonlinearDynamics(t, currentX, linearity, nonlinearity);
		std::vector<double> identity;
		get_identity_matrix(NR,identity);
		apply_approximate_RosenbrockEuler_exponential(NR, NR, linearity, dt, identity, max_exp_order, scratch1, scratch2, default_RE_propagator);
	}
	
	//run simulation
	//default time step is dt, but temporarily dt might be reduced to a smaller value
	for(recorded=1, current_dt=dt; t<endTime; /* increment of t handled in loop */ ){
		// at this stage currentY is guaranteed to be a valid state
		// t should always correspond to currentY
		
		// check runtime-out
		if((runtime_out_seconds>0) && (get_thread_monotonic_walltime_seconds()-start_runtime>=runtime_out_seconds)){
			warningMessage += string(warningMessage=="" ? "" : "\n") + "Simulation was aborted prematurely because it reached the maximum allowed processing time.";
			return true;
		}
		
		// reduce time step to not exceed end time, if needed
		current_dt = min(current_dt, endTime-t);
						
		// get dynamics at the current time t, in the form of the linearity (matrix) and nonlinearity
		if(!k1AlreadyCalculated){
			model.getLinearAndNonlinearDynamics(t, currentX, linearity, nonlinearity);
			k1AlreadyCalculated = true;
		}
		
		// get rate of change: rate = linearity*currentX + nonlinearity
		multiply_matrices(NR,NR,NC,linearity,currentX,rate);
		rate += nonlinearity;
		
		// get X2 by applying the Rosenbrock-Euler propagator
		// X2 = currentX + [exp(current_dt*linearity)-Id]*linearity^{-1}*rate
		// t2 should always correspond to X2
		if(constant_linearity && (abs(current_dt-dt)<RELATIVE_EPSILON*dt)){
			multiply_matrices(NR,NR,NC,default_RE_propagator,rate,X2);
		}else{
			apply_approximate_RosenbrockEuler_exponential(NR, NC, linearity, current_dt, rate, max_exp_order, scratch1, scratch2, X2);
		}
		X2 += currentX;
		t2  = t + current_dt; 
						
		// check if the time step should be reduced (as suggested by the model)
		while((current_dt>min_dt) && model.checkShouldRefineTimeStep(t, currentX, current_dt, X2)){
			// recalculate currentX --> X2 with a smaller time step
			current_dt /= 2;
			apply_approximate_RosenbrockEuler_exponential(NR, NC, linearity, current_dt, rate, max_exp_order, scratch1, scratch2, X2);
			X2 += currentX;
			t2  = t + current_dt; // t2 should always correspond to X2
		}
		
		// check if we crossed the domain boundary, and correct if needed
		// candidate_dt should always correspond to candidateX
		candidateX 		= X2; // candidateX and current_dt may be modified by the model below
		candidate_dt 	= current_dt;
		crossedBoundary = model.checkCrossedDomainBoundaryAndFix(t, currentX, candidate_dt, candidateX);
		if(crossedBoundary==CrossedBoundaryYesButFixedBruteForce){
			if(current_dt>min_dt){
				current_dt /= 2;
				crossedBoundaryButFixedByReducingTimeStep = true;
				continue; // repeat the whole iteration with a smaller time step
			}else{
				// no more time step refinements allowed, and crossing the boundary at even infinitesimal step, so just register candidateX
				currentX 	= candidateX;
				t 			= t2;
				crossedBoundaryYesButFixedUsingEuler = true;
				goto RBE_REGISTER_STATE;
			}
		}else if(crossedBoundary==CrossedBoundaryYesButFixedByReducingTimeStep){
		}
		
		// seems everything worked out as normal, so set the new currentX
		// t should always correspond to currentX
		currentX = candidateX;
		t = t + candidate_dt;
		
		// check and register new state if needed	
		RBE_REGISTER_STATE:
		// check sanity (some ODEs may explode)		
		if(std::isnan(t) || model.stateIsNaN(currentX)){
			warningMessage += string(warningMessage=="" ? "" : "\n") + "Simulation reached NaN.";
			return (recorded>1);
		}
		k1AlreadyCalculated = false;
		if((t-lastRecordedTime > recordingTimeStep) && (recorded<maxRecordedPoints-1)){ // don't record if maxRecordedPoints has been almost exhausted (i.e. 1 remaining), since the final state will be recorded outside of the loop
			model.registerState(t, currentX);
			++recorded;
			lastRecordedTime = t;
			reporter(recorded, maxRecordedPoints, 1-(endTime-t)/simulationTime);
		}
		
		// initialize counters for next iteration (note: this step may be skipped by a 'continue' statement)
		current_dt  = dt;
	}
	
	// register final state if needed
	if(t>lastRecordedTime){
		model.registerState(t, currentX);
	}

	// construct warnings message if needed
	if(crossedBoundaryYesButFixedUsingEuler) warningMessage += string(warningMessage=="" ? "" : "\n") + "Simulation crossed domain boundary and was fixed by locally using forward Euler.";
	if(crossedBoundaryYesButFixedBruteForce) warningMessage += string(warningMessage=="" ? "" : "\n") + "Simulation crossed domain boundary and was fixed by brute-force adjusting the trajectory.";
	if(crossedBoundaryButFixedByReducingTimeStep) warningMessage += string(warningMessage=="" ? "" : "\n") + "Simulation crossed domain boundary and was fixed by locally reducing the time step.";
	
	return true;
}






// Rosenbrock-Euler integrator nested into two-stage Runge-Kutta scheme for solving ODEs of the format:
//    dX/dt = f(t,X) = A(t)*X(t) + N(t,X(t)),
// where the dynamic variable X(t) is a vector or a matrix of size NR x NC, and where A(t) ("linearity") is an explicitly given matrix of size NR x NR and N(t,X(t)) is the non-linear residual.
// This solver uses an exponential integrator for each stage, and includes two stages per iteration. 
// Exponential integrators are particularly suited for stiff problems, when the stiffness is caused by the linear dynamics.
// At each stage, the integration proceeds with the following step:
//	 X(t+eps) = X(t) + [exp(eps*A(t)) - Id] * A(t)^{-1} * f(t,X(t))
// where:
//   f(t,X(t)) = [A(t)*X(t) + N(t,X(t))]
// The exponential * inverse_A term can be approximated as:
//   [exp(eps*A(t)) - Id] * A(t)^{-1} = sum_{k=0}^n (eps/(k+1)) * (eps*A(t))^k/k!
//  Here, we shall call this the Rosenbrock-Euler propagator.
//
// The solver can handle:
//    1. Suggestions by the model for temporarily refined time steps
//    2. Situations where the domain boundary is crossed during an iteration (in which case the time step is temporarily decreased as needed in order to not overshoot)
//
// The model must be able to provide initial conditions for X, the linearity A(t) and the nonlinearity N(t,X(t)) at each time point, as well as storing of the computed trajectory
// For a quick derivation of the Rosenbrock-Euler scheme see: 
//   Chen et al (2018). Exponential Rosenbrock-Euler integrators for elastodynamic simulation. IEEE transactions on visualization and computer graphics. 24:2702-2713
//
// Requested times might reverse (i.e. times at which model dynamics are requested might not be strictly increasing), but recorded time series will be strictly forward in time
template<class COORDINATE, class MODEL, class PROGRESS_REPORTER>
bool RosenbrockEulerRungeKutta2(	const long					NR,								// (INPUT) number of rows in A and X
									const long					NC,								// (INPUT) number of columns in X
									const double				startTime,						// (INPUT) simulation start time
									const double				endTime, 						// (INPUT) simulation end time
									const double 				dt, 							// (INPUT) default integration time step
									MODEL 						&model,							// (INPUT/OUTPUT) object defining model dynamics (via the matrix A), and handling time series storage
									const long					maxRecordedPoints, 				// (INPUT) if small, some intermediate trajectory points are skipped (i.e. not recorded). This might be useful if accurate integration requires a small time step, but would produce a needlesly large time series. If maxRecordedPoints==2, then only the start and end points are recorded. If maxRecordedPoints==1, then only the end point is recorded.
									const long					maxTimeStepRefinements,			// (INPUT) max allowed number of refinements of local time step when encountering invalid states. Only relevant if abortOnInvalidState==false.
									const double				refinement_factor,				// (INPUT) factr by which to refine time steps
									const long 					max_exp_order,					// (INPUT) maximum polynomial order for approximating exponentials (short-term propagators). Typical values are 2-5
									const PROGRESS_REPORTER		&reporter,						// (INPUT) callback functor to handle progress reporting during simulation
									const double				runtime_out_seconds,			// (INPUT) max allowed runtime in seconds. If <=0, this is ignored.
									string						&warningMessage){				// (OUTPUT) will be non-empty in case of error, or if non-fatal problems occurred
	COORDINATE currentX, candidateX, X2, NL1, NL2, rate1, rate2, Rconsensus;
	std::vector<double> A1, A2, Aconsensus, scratch1, scratch2, default_RE_propagator;
	double t=startTime, t2, current_dt1, current_dt2;
	long recorded;
	CrossedBoundary crossedBoundary;
	const double simulationTime 	= endTime-startTime; // target duration of simulation
	const double start_runtime  	= get_thread_monotonic_walltime_seconds();
	const double min_dt 			= dt/pow(refinement_factor,maxTimeStepRefinements);
	const bool constant_linearity 	= model.linearDynamicsAreConstant();
	
	// keep track of warnings
	bool crossedBoundaryButFixedByReducingTimeStep = false;
	bool crossedBoundaryYesButFixedBruteForce = false;
	bool crossedBoundaryYesButFixedUsingEuler = false;
		
	//preliminary error checking
	warningMessage = "";
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
	
	const double recordingTimeStep = simulationTime/max(1L,maxRecordedPoints-1);
	model.reserveSpaceForTimeSeries(maxRecordedPoints);
	bool k1AlreadyCalculated = false;
	
	// get initial state (shape and scale)
	if(!model.getInitialState(t, currentX)){
		warningMessage = "Failed to get initial state";
		return false;
	}
		
	// record initial state
	double lastRecordedTime = -INFTY_D;
	if(maxRecordedPoints>1){
		model.registerState(t, currentX);
		lastRecordedTime = t;
	}
	
	if(constant_linearity){
		// since linearity is always the same, precompute the default Rosenbrock-Euler propagator (i.e. for the default time step)
		model.getLinearAndNonlinearDynamics(t, currentX, A1, NL1);
		std::vector<double> identity;
		get_identity_matrix(NR,identity);
		apply_approximate_RosenbrockEuler_exponential(NR, NR, A1, dt, identity, max_exp_order, scratch1, scratch2, default_RE_propagator);
	}
	
	//run simulation
	//default time step is dt, but temporarily dt might be reduced to a smaller value
	for(recorded=1, current_dt1=current_dt2=dt; t<endTime; /* increment of t handled in loop */ ){
		// at this stage currentY is guaranteed to be a valid state
		// t should always correspond to currentY
		// t2 should always correspond to X2
		// current_dt1 is used for moving from currentX-->X2, which is an intermediate point for obtaining a second "take" on the local dynamics
		// current_dt2 is used for moving from currentX-->candidateX (= potential next currentX), based on the dynamics averaged between currentX and X2
		
		// check runtime-out
		if((runtime_out_seconds>0) && (get_thread_monotonic_walltime_seconds()-start_runtime>=runtime_out_seconds)){
			warningMessage += string(warningMessage=="" ? "" : "\n") + "Simulation was aborted prematurely because it reached the maximum allowed processing time.";
			return true;
		}
		
		// reduce time step dt1 to not exceed dt2, nor overshoot past endTime
		current_dt1 = min(min(current_dt1,current_dt2), endTime-t);
		
		// get dynamics at the current time t, in the form of the linearity (matrix) and nonlinearity
		// the rate of change will be: rate = linearity*currentX + nonlinearity
		if(!k1AlreadyCalculated){
			model.getLinearAndNonlinearDynamics(t, currentX, A1, NL1);
			k1AlreadyCalculated = true;
			multiply_matrices(NR,NR,NC,A1,currentX,rate1);
			rate1 += NL1;
		}
		
		// get X2 by applying the Rosenbrock-Euler propagator
		if(constant_linearity && (abs(current_dt1-dt)<RELATIVE_EPSILON*dt)){
			multiply_matrices(NR,NR,NC,default_RE_propagator,rate1,X2);
		}else{
			apply_approximate_RosenbrockEuler_exponential(NR, NC, A1, current_dt1, rate1, max_exp_order, scratch1, scratch2, X2);
		}
		X2 += currentX;
		t2  = t + current_dt1; // t2 should always correspond to X2

		// check if X2 crossed the boundary and move X2 backward as much as needed
		// X2 and current_dt1 may be adjusted if this routine returns something other than CrossedBoundaryNo
		crossedBoundary = model.checkCrossedDomainBoundaryAndFix(t, currentX, current_dt1, X2);
		if(crossedBoundary == CrossedBoundaryYesButFixedBruteForce){ crossedBoundaryYesButFixedBruteForce = true; }
		else if(crossedBoundary == CrossedBoundaryYesButFixedByReducingTimeStep){ crossedBoundaryButFixedByReducingTimeStep = true; }
						
		// check if the time step should be reduced (as suggested by the model)
		while((current_dt1>min_dt) && model.checkShouldRefineTimeStep(t, currentX, current_dt1, X2)){
			// recalculate currentX --> X2 with a smaller time step
			current_dt1 /= refinement_factor;
			apply_approximate_RosenbrockEuler_exponential(NR, NC, A1, current_dt1, rate1, max_exp_order, scratch1, scratch2, X2);
			X2 += currentX;
			t2  = t + current_dt1; // t2 should always correspond to X2
		}
		
		// stage-2 time-step should not be larger than the stage-1 time-step, but also not below min_dt
		current_dt2 = max(min(min_dt,endTime-t), min(current_dt1,current_dt2));
				
		// get dynamics at t2
		model.getLinearAndNonlinearDynamics(t2, X2, A2, NL2);
		multiply_matrices(NR,NR,NC,A2,X2,rate2);
		rate2 += NL2;

		// use rate1 & rate2 to move the currentX forward, and check if time step was sufficiently small
		linear_combination(0.5, A1, 0.5, A2, Aconsensus);
		linear_combination(0.5, rate1, 0.5, rate2, Rconsensus);
		if(constant_linearity && (abs(current_dt2-dt)<RELATIVE_EPSILON*dt)){
			multiply_matrices(NR,NR,NC,default_RE_propagator,Rconsensus,candidateX);
		}else{
			apply_approximate_RosenbrockEuler_exponential(NR, NC, Aconsensus, current_dt2, Rconsensus, max_exp_order, scratch1, scratch2, candidateX);
		}
		candidateX += currentX;

		// check if we crossed the domain boundary, and correct if needed
		crossedBoundary = model.checkCrossedDomainBoundaryAndFix(t, currentX, current_dt2, candidateX);
		if(crossedBoundary==CrossedBoundaryYesButFixedBruteForce){
			if(current_dt2>min_dt){
				current_dt2 /= refinement_factor;
				crossedBoundaryButFixedByReducingTimeStep = true;
				continue; // repeat the whole iteration with a smaller time step
			}else{
				// use prediction from 1-stage Rosenbrock-Euler scheme for this time step, since kConsensus throws us out of the boundary at even arbitrarily small time steps
				currentX	= X2;
				t 			= t2;
				crossedBoundaryYesButFixedUsingEuler = true;
				goto SRERK2_REGISTER_STATE;
			}
		}

		// check if the time step should be reduced (as suggested by the model)
		if((current_dt2>min_dt) && model.checkShouldRefineTimeStep(t, currentX, current_dt2, candidateX)){
			current_dt2 /= refinement_factor;
			continue; // repeat this iteration with a smaller time step, forget about the current candidateX
		}
		
		// seems everything worked out as normal, so set the new currentX
		currentX = candidateX;
		t = t + current_dt2; // t should always correspond to (currentX, currentS)
		goto SRERK2_REGISTER_STATE;
		
		// check and register new state if needed	
		SRERK2_REGISTER_STATE:
		// check sanity (some ODEs may explode)		
		if(std::isnan(t) || model.stateIsNaN(currentX)){
			warningMessage += string(warningMessage=="" ? "" : "\n") + "Simulation reached NaN.";
			return (recorded>1);
		}
		k1AlreadyCalculated = false;
		if((t-lastRecordedTime > recordingTimeStep) && (recorded<maxRecordedPoints-1)){ // don't record if maxRecordedPoints has been almost exhausted (i.e. 1 remaining), since the final state will be recorded outside of the loop
			model.registerState(t, currentX);
			++recorded;
			lastRecordedTime = t;
			reporter(recorded, maxRecordedPoints, 1-(endTime-t)/simulationTime);
		}
		
		// initialize counters for next iteration (note: this step may be skipped by a 'continue' statement)
		current_dt1  = dt;
		current_dt2  = dt;
	}
	
	// register final state if needed
	if(t>lastRecordedTime){
		model.registerState(t, currentX);
	}

	// construct warnings message if needed
	if(crossedBoundaryYesButFixedUsingEuler) warningMessage += string(warningMessage=="" ? "" : "\n") + "Simulation crossed domain boundary and was fixed by locally using 1-stage Rosenbrock-Euler.";
	if(crossedBoundaryYesButFixedBruteForce) warningMessage += string(warningMessage=="" ? "" : "\n") + "Simulation crossed domain boundary and was fixed by brute-force adjusting the trajectory.";
	if(crossedBoundaryButFixedByReducingTimeStep) warningMessage += string(warningMessage=="" ? "" : "\n") + "Simulation crossed domain boundary and was fixed by locally reducing the time step.";
	
	return true;
}





class ProgressReporter{
private:
	long reportCount;
	string prefix, suffix;
	bool asPercentage;
	mutable long lastReportedCase;
	mutable double lastReportedFraction;
	bool silent;
public:
	ProgressReporter(const bool _silent){ reportCount = 10; lastReportedCase=-1; lastReportedFraction = 0; silent = _silent; }
	ProgressReporter(long _reportCount){ reportCount = _reportCount; lastReportedCase=-1; lastReportedFraction=0; }
	ProgressReporter(long _reportCount, const string &_prefix, const string &_suffix, bool _asPercentage){ reportCount = _reportCount; prefix = _prefix; suffix = _suffix; asPercentage = _asPercentage; lastReportedCase = -1; lastReportedFraction = 0; }
	
	void setReportCount(long count){ reportCount = count; }
	void setPrefix(const string &_prefix){ prefix = _prefix; }
	void setSuffix(const string &_suffix){ suffix = _suffix; }
	void setAsPercentage(bool _asPercentage){ asPercentage = _asPercentage; }
	void reset(){ lastReportedCase = -1; lastReportedFraction = 0; }
	
	void operator()(long casesFinished, long totalCases, double exactFraction) const{
		if(reportCount<=0) return;
		if(silent) return;
		const double p = (1.0*casesFinished)/totalCases;
		const double last_p = (1.0*lastReportedCase)/totalCases;
		const double dp = 1.0/reportCount;
		if((lastReportedCase>=0) && (floor(p/dp)<=floor(last_p/dp))) return; //don't report this case
		if(floor(p/dp)==0) return;
		lastReportedCase 	 = casesFinished;
		lastReportedFraction = exactFraction;
		const long rounder = pow(10.0,1+log10(max(reportCount,1l)));
		Rcout << prefix;
		if(asPercentage){ Rcout << long(rounder*100.0*exactFraction)/rounder << " %"; }
		else{ Rcout << long(floor(exactFraction*totalCases)) << " out of " << totalCases; }
		Rcout << suffix;
	}
	
	
	void operator()(double rangeFinished, double totalRange) const{
		if(reportCount<=0) return;
		if(silent) return;
		const double fractionFinished = rangeFinished/totalRange;
		const double dfraction = 1.0/reportCount; // targetted fraction step between recordings
		if((floor(fractionFinished/dfraction)<=floor(lastReportedFraction/dfraction))) return; //don't report this case
		if(floor(fractionFinished/dfraction)==0) return;
		lastReportedFraction = fractionFinished;
		const long rounder = pow(10.0,1+log10(max(reportCount,1l))); // auxiliary factor for rounding reported fractions on printout
		Rcout << prefix;
		if(asPercentage){ Rcout << long(rounder*100.0*fractionFinished)/rounder << " %"; }
		else{ Rcout << rangeFinished << " out of " << totalRange; }
		Rcout << suffix;
	}
};




#pragma mark -
#pragma mark Deterministic diversity dynamics
#pragma mark



// Class for representing simulated diversity at any time point, including cumulative speciation & extinction events
class TreeStateHistory{
public:
	double diversity;				// current number of extant species
	double coalescent_diversity;	// number of species in coalescent final tree
	double Nbirths;					// cumulative number of speciation events since the beginning of a simulation
	double Ndeaths;					// cumulative number of extinction events since the beginning of a simulation
	double Pextinction;				// extinction probability of a size-1 clade until max time. Only makes sense for simulations performed in reverse.
	double Pmissing;				// probability of a size-1 clade missing at max time (either extinction or not discovered). Only makes sense for simulations performed in reverse.
	TreeStateHistory(){ diversity = coalescent_diversity = Nbirths = Ndeaths = Pextinction = 0; }
	TreeStateHistory(const double start_diversity){ diversity = coalescent_diversity = start_diversity; Nbirths = Ndeaths = 0; Pextinction = 0; Pmissing = 0; }
	bool isnan() const{ return (std::isnan(diversity) || std::isnan(coalescent_diversity) || std::isnan(Nbirths) || std::isnan(Ndeaths) || std::isnan(Pextinction) || std::isnan(Pmissing)); }
};


inline TreeStateHistory operator*(TreeStateHistory x, double scalar){
	x.diversity 			*= scalar;
	x.coalescent_diversity 	*= scalar;
	x.Nbirths 				*= scalar;
	x.Ndeaths 				*= scalar;
	x.Pextinction 			*= scalar;
	x.Pmissing	 			*= scalar;
	return x;
}


TreeStateHistory& operator*=(TreeStateHistory &x, double scalar) {
	x.diversity 			*= scalar;
	x.coalescent_diversity	*= scalar;
	x.Nbirths 				*= scalar;
	x.Ndeaths 				*= scalar;
	x.Pextinction			*= scalar;
	x.Pmissing				*= scalar;
	return x;
}


inline TreeStateHistory operator+(TreeStateHistory x, const TreeStateHistory &y){
	x.diversity 			+= y.diversity;
	x.coalescent_diversity	+= y.coalescent_diversity;
	x.Nbirths 				+= y.Nbirths;
	x.Ndeaths 				+= y.Ndeaths;
	x.Pextinction			+= y.Pextinction;
	x.Pmissing				+= y.Pmissing;
	return x;
}


inline TreeStateHistory &operator+=(TreeStateHistory &x, const TreeStateHistory &y){
	x.diversity 			+= y.diversity;
	x.coalescent_diversity	+= y.coalescent_diversity;
	x.Nbirths 				+= y.Nbirths;
	x.Ndeaths 				+= y.Ndeaths;
	x.Pextinction			+= y.Pextinction;
	x.Pmissing				+= y.Pmissing;
	return x;
}


inline TreeStateHistory operator-(TreeStateHistory x, const TreeStateHistory &y){
	x.diversity 			-= y.diversity;
	x.coalescent_diversity 	-= y.coalescent_diversity;
	x.Nbirths 				-= y.Nbirths;
	x.Ndeaths 				-= y.Ndeaths;
	x.Pextinction			-= y.Pextinction;
	x.Pmissing				-= y.Pmissing;
	return x;
}


inline TreeStateHistory &operator-=(TreeStateHistory &x, const TreeStateHistory &y){
	x.diversity 			-= y.diversity;
	x.coalescent_diversity	-= y.coalescent_diversity;
	x.Nbirths 				-= y.Nbirths;
	x.Ndeaths 				-= y.Ndeaths;
	x.Pextinction			-= y.Pextinction;
	x.Pmissing				-= y.Pmissing;
	return x;
}


inline TreeStateHistory operator/(TreeStateHistory x, const TreeStateHistory &y){
	x.diversity 			/= y.diversity;
	x.coalescent_diversity	/= y.coalescent_diversity;
	x.Nbirths 				/= y.Nbirths;
	x.Ndeaths 				/= y.Ndeaths;
	x.Pextinction			/= y.Pextinction;
	x.Pmissing				/= y.Pmissing;
	return x;
}


inline TreeStateHistory operator/(TreeStateHistory x, double scalar){
	x.diversity 			/= scalar;
	x.coalescent_diversity	/= scalar;
	x.Nbirths 				/= scalar;
	x.Ndeaths 				/= scalar;
	x.Pextinction			/= scalar;
	x.Pmissing				/= scalar;
	return x;
}


inline TreeStateHistory &operator/=(TreeStateHistory &x, double scalar){
	x.diversity 			/= scalar;
	x.coalescent_diversity	/= scalar;
	x.Nbirths 				/= scalar;
	x.Ndeaths 				/= scalar;
	x.Pextinction			/= scalar;
	x.Pmissing				/= scalar;
	return x;
}



// class for species speciation-extinction model for tree generation
// includes integration of cumulative birth & death events
// can be integrated in reverse (time counted backwards) as well as in forward direction
// When integrated in reverse, the extinction probability Pextinction(t,T) and the probability of missing Pmissing are also simulated
// In any case, the simulated diversity is the true (total) extant diversity at any time. To make the simulated time series coalescent, use calculate_coalescent_diversities() afterwards.
class TreeSpeciationExtinctionModel{
	double 	min_valid_diversity;
	bool 	reverse;
	double 	reflection_time;
	bool	has_probabilities;
public:
	std::vector<TreeStateHistory> trajectory;
	std::vector<double> times;
	double initial_diversity;
	double rarefaction;

	TreeSpeciationExtinctionModel(){ min_valid_diversity = 0; reverse = false; rarefaction = 1; has_probabilities = false; }
	double birth_rate_intercept, birth_rate_factor, birth_rate_exponent;
	double death_rate_intercept, death_rate_factor, death_rate_exponent;
	long Nsplits;
	
	// use parameters from another model instance
	void adopt_parameters(const TreeSpeciationExtinctionModel &model2){
		birth_rate_intercept 	= model2.birth_rate_intercept;
		birth_rate_factor 		= model2.birth_rate_factor;
		birth_rate_exponent 	= model2.birth_rate_exponent;
		death_rate_intercept 	= model2.death_rate_intercept;
		death_rate_factor 		= model2.death_rate_factor;
		death_rate_exponent 	= model2.death_rate_exponent;
		Nsplits 				= model2.Nsplits;
		min_valid_diversity 	= model2.min_valid_diversity;
		reflection_time 		= model2.reflection_time;
		reverse	 				= model2.reverse;
	}
		
	// model is warned that a time series of size count will need to be stored
	void reserveSpaceForTimeSeries(long Ntimes){ 
		trajectory.clear();
		trajectory.reserve(Ntimes); 
		times.clear();
		times.reserve(Ntimes); 
	}
	
	// reverse the dynamics of the model, i.e. return vector field -V(_reflection_time-t) whenever V(t) is requested.
	// this should be called prior to simulating
	void set_reverse(double _reflection_time){
		reflection_time = _reflection_time;
		reverse = true;
		has_probabilities = true;
	}
	
	// reverse the time course of a simulated trajectory
	// this should be called after simulating
	void reverse_trajectory(const double final_time){
		TreeStateHistory state;
		const long NMT = times.size();
		double time;
		for(long t=0, s; t<NMT/2; ++t){
			s = NMT-t-1;
			state = trajectory[t];
			trajectory[t] = trajectory[s];
			trajectory[s] = state;
			time = times[t];
			times[t] = times[s];
			times[s] = time;
		}
		// correct direction of accumulation of time, Nbirths & Ndeaths
		const double max_Nbirths = trajectory[0].Nbirths;
		const double max_Ndeaths = trajectory[0].Ndeaths;
		for(long t=0; t<NMT; ++t){
			times[t] = final_time - times[t];
			trajectory[t].Nbirths = max_Nbirths - trajectory[t].Nbirths;
			trajectory[t].Ndeaths = max_Ndeaths - trajectory[t].Ndeaths;
		}
	}
	
	// calculate probabilities fo extinction and missing, by integrating in reverse
	// this should be done after the end of a simulation, since the ODE's has initial condition and rate of change depends on the final trajectory of diversity
	// assumes that times[] are in ascending order
	void calculate_probabilities(){
		if(has_probabilities) return;
		const long NMT = times.size();
		trajectory[NMT-1].Pextinction = 0; 					// probability of extinction begins at 0
		trajectory[NMT-1].Pmissing = (1.0-rarefaction); 	// probability of missing begins at (1-rarefacton_fraction)
		for(long t=NMT-2; t>=0; --t){
			const double birth_rate_pc = get_speciation_rate_at_state(times[t+1],trajectory[t+1].diversity)/trajectory[t+1].diversity;
			const double death_rate_pc = get_extinction_rate_at_state(times[t+1],trajectory[t+1].diversity)/trajectory[t+1].diversity;
			const double dt = times[t+1]-times[t];
			trajectory[t].Pextinction = trajectory[t+1].Pextinction + dt*(death_rate_pc - trajectory[t+1].Pextinction*(birth_rate_pc+death_rate_pc) + pow(trajectory[t+1].Pextinction,Nsplits)*birth_rate_pc);
			trajectory[t].Pmissing	  = trajectory[t+1].Pmissing + dt*(death_rate_pc - trajectory[t+1].Pmissing*(birth_rate_pc+death_rate_pc) + pow(trajectory[t+1].Pmissing,Nsplits)*birth_rate_pc);
		}
		has_probabilities = true;
	}
	
	
	void get_coalescent_trajectory(	double 							resolution, 	// (INPUT)
									double 							rarefaction, 	// (INPUT)
									std::vector<TreeStateHistory> 	&coalescent) const{	// (OUTPUT)
		const long NMT = times.size();
		const double final_time = times.back();
		coalescent = trajectory;
		coalescent[NMT-1].Pextinction = 0; // probability of extinction begins at 0
		coalescent[NMT-1].Pmissing = (resolution<=0 ? (1.0-rarefaction) : 0.0);
		long resolution_jump = NMT-1; // will be updated (decreased) when we actually find the resolution jump
		double total_diversity_at_resolution_age, coalescent_diversity_at_resolution_age;
		for(long t=NMT-2; t>=0; --t){
			const double birth_rate_pc = get_speciation_rate_at_state(times[t+1],coalescent[t+1].diversity)/coalescent[t+1].diversity;
			const double death_rate_pc = get_extinction_rate_at_state(times[t+1],coalescent[t+1].diversity)/coalescent[t+1].diversity;
			double dt = times[t+1]-times[t];
			coalescent[t].Pextinction = coalescent[t+1].Pextinction + dt*(death_rate_pc - coalescent[t+1].Pextinction*(birth_rate_pc+death_rate_pc) + pow(coalescent[t+1].Pextinction,Nsplits)*birth_rate_pc);
			coalescent[t].Pmissing	  = coalescent[t+1].Pmissing + dt*(death_rate_pc - coalescent[t+1].Pmissing*(birth_rate_pc+death_rate_pc) + pow(coalescent[t+1].Pmissing,Nsplits)*birth_rate_pc);
			const double last_age = final_time-times[t+1];
			const double new_age  = final_time-times[t];
			if((new_age>=resolution) && (last_age<resolution)){
				// we just jumped over rarefaction_age, so apply rarefaction to Pmissing
				resolution_jump = t+1;
				total_diversity_at_resolution_age = trajectory[t+1].diversity + (trajectory[t].diversity-trajectory[t+1].diversity)*(resolution-last_age)/(new_age-last_age);
				double Pmissing_at_resolution_age = coalescent[t+1].Pmissing + (coalescent[t].Pmissing-coalescent[t+1].Pmissing)*(resolution-last_age)/(new_age-last_age); // linearly interpolate to get Pmissing at resolution age
				Pmissing_at_resolution_age = 1 - rarefaction*(1-Pmissing_at_resolution_age); // apply rarefaction at resolution age
				coalescent_diversity_at_resolution_age = total_diversity_at_resolution_age * (1-Pmissing_at_resolution_age);
				dt = new_age-resolution;
				coalescent[t].Pmissing = Pmissing_at_resolution_age + dt*(death_rate_pc - Pmissing_at_resolution_age*(birth_rate_pc+death_rate_pc) + pow(Pmissing_at_resolution_age,Nsplits)*birth_rate_pc);
			}
		}
		// flatten coalescent diversity for ages < resolution_age
		for(long t=NMT-2; t>=0; --t){
			if((final_time-times[t])<resolution){
				coalescent[t].coalescent_diversity = coalescent_diversity_at_resolution_age;
				coalescent[t].Pmissing = coalescent_diversity_at_resolution_age/coalescent[t].diversity;
			}else{
				coalescent[t].coalescent_diversity = coalescent[t].diversity * (1-coalescent[t].Pmissing);
			}
		}
	}
	
	
	// calculate coalescent diversities time series, based on the diversities & Pmissing
	void calculate_coalescent_diversities(){
		calculate_probabilities();  // make sure probabilities have been calculated
		for(long t=0; t<times.size(); ++t){
			trajectory[t].coalescent_diversity = trajectory[t].diversity * (1-trajectory[t].Pmissing);
		}
	}

	bool getInitialState(double time, TreeStateHistory &state) const{ 
		state = TreeStateHistory(initial_diversity); 
		if(reverse){
			state.Pextinction = 0;
			state.Pmissing = (1.0-rarefaction);
		}else{
			// not implemented for forward integration; Pextinction & Pmissing can be calculated afterwards (in reverse) using calculate_probabilities()
			state.Pextinction = 0;
			state.Pmissing = 0;
		}
		return true; 
	}

	// record time series point
	void registerState(double time, const TreeStateHistory &state){
		trajectory.push_back(state); 
		times.push_back(time); 
		if(reverse){
			trajectory.back().Pextinction = min(1.0,max(0.0,trajectory.back().Pextinction));
			trajectory.back().Pmissing = min(1.0,max(0.0,trajectory.back().Pmissing));
		}
	}

	double get_speciation_rate_at_state(double time, const double diversity) const{
		return birth_rate_intercept + birth_rate_factor*pow(diversity,birth_rate_exponent);
	}
	
	double get_extinction_rate_at_state(double time, const double diversity) const{
		return death_rate_intercept + death_rate_factor*pow(diversity,death_rate_exponent);
	}

	RequestedDynamics getRateOfChangeAtState(double time, const TreeStateHistory &current_state, TreeStateHistory &rate_of_change, TreeStateHistory &jump_state){
		if(reverse) time = reflection_time-time;
		const double N 	= current_state.diversity;
		const double Pe = current_state.Pextinction;
		const double Pm = current_state.Pmissing;
		const double B 	= get_speciation_rate_at_state(time,current_state.diversity);
		const double D 	= get_extinction_rate_at_state(time,current_state.diversity);
		rate_of_change.Nbirths 		= B;
		rate_of_change.Ndeaths 		= D;
		rate_of_change.diversity 	= (Nsplits-1)*B - D;
		if(reverse){
			rate_of_change.Pextinction	= D/N - Pe*(B+D)/N + pow(Pe,Nsplits)*B/N;
			rate_of_change.Pmissing		= D/N - Pm*(B+D)/N + pow(Pm,Nsplits)*B/N;
		}else{
			rate_of_change.Pextinction = 0; // not implemented for forward integration
			rate_of_change.Pmissing = 0; // not implemented for forward integration			
		}
		if(reverse) rate_of_change.diversity = -rate_of_change.diversity;
		return RequestedDynamicsRateOfChange; 
	}

	// returns true if for some reason the time step should be refined, e.g. if the domain boundary is crossed
	bool checkShouldRefineTimeStep(double time, const TreeStateHistory &current_state, double dt, const TreeStateHistory &candidate_state) const{ 
		if(reverse){
			return (candidate_state.diversity<min_valid_diversity) || (candidate_state.Pextinction<0) || (candidate_state.Pextinction>1) || (candidate_state.Pmissing<0) || (candidate_state.Pmissing>1); 
		}else{
			return (candidate_state.diversity<min_valid_diversity);
		}
	}
	
	// check if candidate_state is outside of the domain boundaries.
	// In that case, tries to correct the candidate_state to be the "last" valid state on the linear trajectory
	// If this is not possible (e.g. because previous_state was already on the boundary), the problematic components are brute-force adjusted to be within the domain. In this case, the routine returns CrossedBoundaryYesButFixedBruteForce.
	CrossedBoundary checkCrossedDomainBoundaryAndFix(	double					previous_time,
														const TreeStateHistory	&previous_state,		// previous state (assumed to be valid!)
														double 					&dt,					// (INPUT/OUTPUT) will be adjusted (reduced) if candidate_state crossed the boundary. The modified value is guaranteed to be within (0,dt]
														TreeStateHistory		&candidate_state,		// (INPUT/OUTPUT) if candidate_state is outside of domain boundaries, then this may become the "last" state (on the linear trajectory from previous_state to candidate_state) within the domain (if possible).
														const bool				intermediate) const{	// (INPUT) is the candidate point an intermediate point (i.e. as used in Runge-Kutta schemes), or a final point (i.e., the final prediction for the candidate time)
		if(candidate_state.diversity>=min_valid_diversity){
			return CrossedBoundaryNo;
		}else if(previous_state.diversity>min_valid_diversity){
			double lambda = 1;
			if(candidate_state.diversity<min_valid_diversity) lambda = min(lambda,(previous_state.diversity-min_valid_diversity)/(previous_state.diversity-candidate_state.diversity));
			candidate_state = previous_state*(1-lambda) + candidate_state*lambda;
			dt *= lambda;
			return CrossedBoundaryYesButFixedByReducingTimeStep;
		}else{
			candidate_state.diversity 	= max(min_valid_diversity,candidate_state.diversity);
			candidate_state.Pextinction 	= min(1.0,max(0.0,candidate_state.Pextinction));
			candidate_state.Pmissing 	= min(1.0,max(0.0,candidate_state.Pmissing));
			return CrossedBoundaryYesButFixedBruteForce;
		}
		/*
		if((candidate_state.diversity>=min_valid_diversity) && ((!reverse) || ((candidate_state.Pextinction>=0) && (candidate_state.Pextinction<=1)))){
			return CrossedBoundaryNo;
		}else if((previous_state.diversity>min_valid_diversity) && ((!reverse) || ((previous_state.Pextinction>=0) && (previous_state.Pextinction<=1)))){
			double lambda = 1;
			if(candidate_state.diversity<min_valid_diversity) lambda = min(lambda,(previous_state.diversity-min_valid_diversity)/(previous_state.diversity-candidate_state.diversity));
			if(reverse && (candidate_state.Pextinction<0)) lambda = min(lambda,(previous_state.Pextinction-0)/(previous_state.Pextinction-candidate_state.Pextinction));
			if(reverse && (candidate_state.Pextinction>1)) lambda = min(lambda,(1-previous_state.Pextinction)/(candidate_state.Pextinction-previous_state.Pextinction));
			candidate_state = previous_state*(1-lambda) + candidate_state*lambda;
			dt *= lambda;
			return CrossedBoundaryYesButFixedByReducingTimeStep;
		}else{
			candidate_state.diversity 	= max(min_valid_diversity,candidate_state.diversity);
			candidate_state.Pextinction 	= min(1.0,max(0.0,candidate_state.Pextinction));
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
														const double		resolution,				// (INPUT) optional resoluton at which to collapse final tree (prior to any rarefaction). Set this to 0 for no collapsing. Note that collapsing affects the last part of the coalescent diversity curve.
														const double		rarefaction,			// (INPUT) optional rarefaction fraction to apply, i.e. fraction of tips remaining after rarefaction. Rarefaction is assumed to occur after collapsing at the given resolution (i.e. the collapsed tree is rarefied). Set this to 1 for no rarefaction. Note that rarefaction affects the full coalescent diversity curve in a non-uniform manner.
														const long			Nsplits,				// (INPUT) number of children to create at each diversification event. Must be at least 2. For a bifurcating tree this should be set to 2. If >2, then the tree will be multifurcating.
														const NumericVector	&times,					// (INPUT) times at which to calculate deterministic diversities, in ascending order
														const double		start_time,				// (INPUT) Time at which to start simulation. This should be the beginning of the tree (<=times[0]).
														const double		final_time,				// (INPUT) Time at which to end simulation. This should be the ending time of the tree (>= times.back()).
														const double		start_diversity,		// (INPUT) diversity of extant clades at start_time. If reverse==true, this is the visible diversity after rarefaction. Only relevant if reverse=false.
														const double		final_diversity,		// (INPUT) total extant diversity at final_time, without any rarefaction or collapsing at some resolution. Only relevant if reverse==true.
														const bool			reverse,				// (OUTPUT) if true, then the tree model is integrated in backward time direction. In that case, start_diversity is interpreted as the true diversity at times.back()
														const bool			include_coalescent,			// (INPUT) if true, the coalescent diversities are also returned. These are the diversities that would remain in a coalescent tree (i.e. only including ancestors of extant tips)
														const bool			include_probabilities,		// (INPUT) if true, then Prepresentation (for each time point) will also be returned. This only makes sense for reverse integrations
														const bool			include_birth_rates,		// (INPUT) if true, the speciation rates corresponding to times[] will also be returned
														const bool			include_death_rates,		// (INPUT) if true, the extinction rates corresponding to times[] will also be returned
														const bool			include_Nevents,			// (INPUT) if true, then the total predicted birth events (starting from times[0] and up to each time point) will also be calculated and returned
														const double		runtime_out_seconds){		// (INPUT) max allowed runtime in seconds. If <=0, this option is ignored.
	const long NT = times.size();
	
	// prepare returned data structures, fill later
	std::vector<double> birth_rates, death_rates, Nbirths, Ndeaths;
	std::vector<double> coalescent_diversities, total_diversities, Psurvival, Prepresentation;
	
	if((birth_rate_intercept==0) && (birth_rate_exponent==1) && (death_rate_intercept==0) && (death_rate_exponent==1) && (Nsplits==2)){
		// special case: constant per-capita speciation & extinction rates
		// use known analytical solutions
		total_diversities.resize(NT);
		if(include_coalescent) coalescent_diversities.resize(NT);
		if(include_probabilities){
			Psurvival.resize(NT);
			Prepresentation.resize(NT);
		}
		if(include_birth_rates) birth_rates.resize(NT);
		if(include_death_rates) death_rates.resize(NT);
		if(include_Nevents){ 
			Nbirths.resize(NT); 
			Ndeaths.resize(NT); 
		}
		const double diversification_rate = birth_rate_factor-death_rate_factor;
		double Pmissing;
		double Pmissing_at_resolution_age, total_diversity_at_resolution_age, coalescent_diversity_at_resolution_age;
		// calculate total & coalescent diversity at resolution_age, as well as Pmissing at resolution_age.
		if(resolution>0){
			if(birth_rate_factor==death_rate_factor){
				Pmissing_at_resolution_age = 1 - 1/(1+resolution*birth_rate_factor);
			}else{
				Pmissing_at_resolution_age = 1 - diversification_rate/(birth_rate_factor - death_rate_factor*exp(-resolution*diversification_rate));
			}
			if(reverse){
				total_diversity_at_resolution_age = final_diversity * exp(-resolution*diversification_rate);
			}else{
				total_diversity_at_resolution_age = start_diversity * exp((final_time-resolution-start_time)*diversification_rate);
			}
			coalescent_diversity_at_resolution_age = total_diversity_at_resolution_age * rarefaction * (1-Pmissing_at_resolution_age);
			
		}
		// calculate total & coalescnet diversity and other variables for all time points
		for(long t=0; t<NT; ++t){
			const double age = final_time - times[t];
			if(reverse){
				total_diversities[t] = final_diversity * exp(-age*diversification_rate);
			}else{
				total_diversities[t] = start_diversity * exp((times[t]-start_time)*diversification_rate);
			}
			if(age>=resolution){
				const double effective_rarefaction = rarefaction * (1-Pmissing_at_resolution_age);
				const double effective_age = age-resolution;
				if(birth_rate_factor==death_rate_factor){
					Pmissing = 1 - effective_rarefaction/(1+effective_rarefaction*effective_age*birth_rate_factor);
				}else{
					Pmissing = 1 - effective_rarefaction*diversification_rate/(effective_rarefaction*birth_rate_factor + ((1-effective_rarefaction)*birth_rate_factor - death_rate_factor)*exp(-effective_age*diversification_rate));
				}
			}else{
				Pmissing = 1.0 - coalescent_diversity_at_resolution_age/total_diversities[t]; // effective Pmissing, based on coalescent and total diversity after collapsing and rarefaction
			}
			if(include_probabilities){
				if(birth_rate_factor==death_rate_factor){
					Psurvival[t] = 1/(1+age*birth_rate_factor);
				}else{
					Psurvival[t] = diversification_rate/(birth_rate_factor - death_rate_factor*exp(-age*diversification_rate));
				}
				Prepresentation[t] = 1-Pmissing;
			}
			if(include_birth_rates) birth_rates[t] = birth_rate_factor * total_diversities[t];
			if(include_death_rates) death_rates[t] = death_rate_factor * total_diversities[t];
			if(include_Nevents){
				Nbirths[t] = total_diversities[0]*(birth_rate_factor/diversification_rate) * (exp(diversification_rate*(times[t]-times[0])) - 1);
				Ndeaths[t] = total_diversities[0]*(death_rate_factor/diversification_rate) * (exp(diversification_rate*(times[t]-times[0])) - 1);
			}
			if(include_coalescent) coalescent_diversities[t] = total_diversities[t] * (1-Pmissing);
		}
		
	}else{
		// general case: integrate ODEs using Runge-Kutta
		
		// prepare data structures for storing raw simulation results (will be interpolated onto requested time points afterwards)
		std::vector<double> model_times;
		std::vector<TreeStateHistory> model_trajectory;
		long NMT;

		// initialize model for tree growth
		TreeSpeciationExtinctionModel model;
		model.birth_rate_intercept 	= birth_rate_intercept;
		model.birth_rate_factor 	= birth_rate_factor;
		model.birth_rate_exponent 	= birth_rate_exponent;
		model.death_rate_intercept 	= death_rate_intercept;
		model.death_rate_factor 	= death_rate_factor;
		model.death_rate_exponent 	= death_rate_exponent;
		model.Nsplits				= Nsplits;
		model.initial_diversity		= (reverse ? final_diversity : start_diversity);
		if(reverse) model.set_reverse(final_time); // time-reverse the dynamics (vector field) of the model, in order to the simulate backwards in time using the conventional Runge-Kitta integrator
		string warningMessage;
		const double default_dt = min((final_time-times[0])/NT, (reverse ? 1e-1*final_diversity/model.estimate_max_rate_of_change(final_time, final_diversity,false) : (final_time-times[0])/NT)); // guess appropriate simulation time step
		
		if((resolution>0) && reverse){
			// run simulation in 2 chunks, one for ages<resolution and then again for ages>resolution
			model.rarefaction = 1; // rarefaction occurs after collapsing tree at resolution, so simulate first part as if there was no rarefaction
			bool success = RungeKutta2<TreeStateHistory,TreeSpeciationExtinctionModel,ProgressReporter>
										(0,
										resolution,
										min(resolution/100, 1e-1*final_diversity/model.estimate_max_rate_of_change(final_time, final_diversity,false)),
										model,
										3, // number of points to record
										2,
										2, // refinement_factor
										ProgressReporter(true),
										runtime_out_seconds,
										warningMessage);
			if((!success) || (model.times.back()<resolution)) return Rcpp::List::create(Rcpp::Named("success") = false, Rcpp::Named("error") = "Could not integrate backwards until resolution age: "+warningMessage ); // simulation failed
			// extract some information about the resolution_age (last time point of simulation)
			const long NMTr = model.times.size();
			const double total_diversity_at_resolution_age 		= model.trajectory[NMTr-1].diversity;
			const double Pmissing_at_resolution_age 			= model.trajectory[NMTr-1].Pmissing;
			const double coalescent_diversity_at_resolution_age = total_diversity_at_resolution_age * rarefaction * (1-Pmissing_at_resolution_age);
			
			// integrate further back in time, starting at resolution_age
			TreeSpeciationExtinctionModel model2;
			model2.adopt_parameters(model);
			model2.initial_diversity = total_diversity_at_resolution_age;
			model2.rarefaction = rarefaction*(1-Pmissing_at_resolution_age); // effective rarefaction at resolution age (i.e. after collapsing and rarefaction)
			success = RungeKutta2<TreeStateHistory,TreeSpeciationExtinctionModel,ProgressReporter>
										(resolution,
										default_dt + (final_time-times[0]),
										default_dt,
										model2,
										10*times.size(),
										2,
										2, // refinement_factor
										ProgressReporter(true),
										runtime_out_seconds,
										warningMessage);
			if(!success) return Rcpp::List::create(Rcpp::Named("success") = false, Rcpp::Named("error") = "Could not integrate backwards past resolution age: "+warningMessage ); // simulation failed

			// reverse time course of simulated trajectories
			model.reverse_trajectory(final_time); 
			model2.reverse_trajectory(final_time);
			
			// until resolution_age, coalescent diversities and Pmissing are special, due to collapsing of the tree at the given resolution
			for(long t=0; t<NMTr; ++t){
				model.trajectory[t].coalescent_diversity = coalescent_diversity_at_resolution_age;
				model.trajectory[t].Pmissing = coalescent_diversity_at_resolution_age/model.trajectory[t].diversity;
			}

			if(include_probabilities) model2.calculate_probabilities(); // make sure probabilities have been calculated (may not be the case if simulation was not in reverse). This must come after reversing the trajectory.
			if(include_coalescent) model2.calculate_coalescent_diversities(); // make sure coalescent diversities have been calculated for the trajectory (may not be the case if simulation was not in reverse). This must come after reversing the trajectory.
			
			// adjust Nbirths & Ndeaths
			if(include_Nevents){
				for(long t=0; t<NMTr; ++t){
					model.trajectory[t].Nbirths += model2.trajectory.back().Nbirths;
					model.trajectory[t].Ndeaths += model2.trajectory.back().Ndeaths;
				}
			}

			// extract and concatenate trajectories from the two simulations
			NMT = NMTr+model2.times.size();
			model_times.reserve(NMT);
			model_times.insert(model_times.end(), model2.times.begin(), model2.times.end());
			model_times.insert(model_times.end(), model.times.begin(), model.times.end());
			model_trajectory.reserve(NMT);
			model_trajectory.insert(model_trajectory.end(), model2.trajectory.begin(), model2.trajectory.end());
			model_trajectory.insert(model_trajectory.end(), model.trajectory.begin(), model.trajectory.end());


		}else{
			// run simulation for all ages
			// note that simulation a priori calculates total extant diversities over time, we make the diversities coalescent afterwards
			model.rarefaction = rarefaction;
			const bool success = RungeKutta2<TreeStateHistory,TreeSpeciationExtinctionModel,ProgressReporter>
										((reverse ? 0 : start_time),
										default_dt+(reverse ? final_time-times[0] : final_time),
										default_dt,
										model,
										10*times.size(),
										2,
										2, // refinement_factor
										ProgressReporter(true),
										runtime_out_seconds,
										warningMessage);
			if(!success) return Rcpp::List::create(Rcpp::Named("success") = false, Rcpp::Named("error") = warningMessage ); // simulation failed

			if(reverse) model.reverse_trajectory(final_time); // reverse time course of simulated trajectory
			
			if((!reverse) && (resolution>0)){
				model.get_coalescent_trajectory(resolution, rarefaction, model_trajectory);
			}else{
				if(include_probabilities) model.calculate_probabilities(); // make sure probabilities have been calculated (may not be the case if simulation was not in reverse). This must come after reversing the trajectory.
				if(include_coalescent) model.calculate_coalescent_diversities(); // make sure coalescent diversities have been calculated for the trajectory (may not be the case if simulation was not in reverse). This must come after reversing the trajectory.
				model_trajectory = model.trajectory;
			}
			model_times = model.times;
			NMT = model_times.size();
		}
		
		// extract total diversities (extinct+extant) from simulated trajectory
		long includedNewTimesStart, includedNewTimesEnd;
		std::vector<double> model_total_diversities(NMT);
		for(long t=0; t<NMT; ++t) model_total_diversities[t] = model_trajectory[t].diversity;
		interpolateTimeSeriesAtTimes(model_times, model_total_diversities, 0, NMT-1, times, 0, NT-1, includedNewTimesStart, includedNewTimesEnd, total_diversities);

		if(include_coalescent){
			// extract coalescent diversities from simulated trajectory
			std::vector<double> model_coalescent_diversities(NMT);
			for(long t=0; t<NMT; ++t) model_coalescent_diversities[t] = model_trajectory[t].coalescent_diversity;
			interpolateTimeSeriesAtTimes(model_times, model_coalescent_diversities, 0, NMT-1, times, 0, NT-1, includedNewTimesStart, includedNewTimesEnd, coalescent_diversities);
		}
	
		if(include_probabilities){
			// extract probabilities of lineage survival & representation from simulated trajectory
			std::vector<double> model_Prepresentation(NMT);
			std::vector<double> model_Psurvival(NMT);
			for(long t=0; t<NMT; ++t){
				model_Psurvival[t] = 1-model_trajectory[t].Pextinction;
				model_Prepresentation[t] = 1-model_trajectory[t].Pmissing;
			}
			interpolateTimeSeriesAtTimes(model_times, model_Prepresentation, 0, NMT-1, times, 0, NT-1, includedNewTimesStart, includedNewTimesEnd, Prepresentation);
			interpolateTimeSeriesAtTimes(model_times, model_Psurvival, 0, NMT-1, times, 0, NT-1, includedNewTimesStart, includedNewTimesEnd, Psurvival);
		}

		// extract cumulative speciation & extinction events, if needed
		if(include_Nevents){
			// extract from model trajectory
			std::vector<double> model_Nbirths(NMT), model_Ndeaths(NMT);
			for(long t=0; t<NMT; ++t){
				model_Nbirths[t] = model_trajectory[t].Nbirths;
				model_Ndeaths[t] = model_trajectory[t].Ndeaths;
			}
			interpolateTimeSeriesAtTimes(model_times, model_Nbirths, 0, NMT-1, times, 0, NT-1, includedNewTimesStart, includedNewTimesEnd, Nbirths);
			interpolateTimeSeriesAtTimes(model_times, model_Ndeaths, 0, NMT-1, times, 0, NT-1, includedNewTimesStart, includedNewTimesEnd, Ndeaths);
		}

		// calculate speciation rates for requested times, if needed
		if(include_birth_rates){
			birth_rates.resize(NT);
			for(long t=0; t<NT; ++t) birth_rates[t] = model.get_speciation_rate_at_state(times[t],total_diversities[t]);
		}
		// calculate extinction rates for requested times, if needed
		if(include_death_rates){
			death_rates.resize(NT);
			for(long t=0; t<NT; ++t) death_rates[t] = model.get_extinction_rate_at_state(times[t],total_diversities[t]);
		}
	}
	

	 		
	return Rcpp::List::create(	Rcpp::Named("success")					= true,
								Rcpp::Named("coalescent_diversities") 	= Rcpp::wrap(coalescent_diversities),	// coalescent diversity at each time point (if requested)
								Rcpp::Named("total_diversities")		= Rcpp::wrap(total_diversities),	// total (extant+extinct) diversity at each time point
								Rcpp::Named("Psurvival")				= Rcpp::wrap(Psurvival),
								Rcpp::Named("Prepresentation")			= Rcpp::wrap(Prepresentation),
								Rcpp::Named("birth_rates")				= Rcpp::wrap(birth_rates),
								Rcpp::Named("death_rates")				= Rcpp::wrap(death_rates),
								Rcpp::Named("Nbirths")					= Nbirths,
								Rcpp::Named("Ndeaths")					= Ndeaths);
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
		if(!smoothenTimeSeriesSavitzkyGolay(times, raw_coalescent_diversities, 0.0, smoothing_span, min(smoothing_span-2,smoothing_order), true, smooth_coalescent_diversities)){
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
	std::vector<double> Prepresentation(NT), total_diversities(NT), nu(NT), pulled_total_diversities(NT);
	for(long t=NT-1, tl, tr; t>=0; --t){
		tl = (t==0 ? 0 : t-1);
		tr = (t==0 ? 1 : t); // use left-sided derivative, to avoid the discontinuity at age=0
		const double birth_rate_pc = (constant_birth_rate_pc ? last_birth_rate_pc : birth_rates_pc[t]);
		nu[t] = -log(coalescent_diversities[tr]/coalescent_diversities[tl])/(times[tr]-times[tl]); // (1/Nc)*dNc/dtau = dlog(Nc)/dtau, where Nc: coalescent_diversity
		Prepresentation[t] = -nu[t]/birth_rate_pc;
		if((Prepresentation[t]<=1e-8) && (coalescent_diversities[t]>0)){
			// diversity not defined (Prepresentation seems spuriously low)
			total_diversities[t] = NAN_D;
			Prepresentation[t] 	= max(0.0, Prepresentation[t]);
		}else if(Prepresentation[t]>1){
			// diversity not defined (Prepresentation > 1)
			Prepresentation[t]  = 1;
			total_diversities[t] = NAN_D;			
		}else{
			total_diversities[t] = coalescent_diversities[t]/Prepresentation[t];
		}
		pulled_total_diversities[t] = - last_birth_rate_pc*coalescent_diversities[t]/nu[t];
	}
	total_diversities[NT-1] = coalescent_diversities[NT-1]/rarefaction; // coalescent diversity last time point is assumed to be equal to true diversity multiplied by rarefaction fraction
	
	// calculate birth & death rates & pulled diversification rate
	std::vector<double> death_rates(NT), death_rates_pc(NT), birth_rates(NT), pulled_diversification_rates(NT), pulled_extinction_rates(NT);
	for(long t=NT-1, tl, tr; t>=0; --t){
		tl = (t==0 ? 0 : t-1);
		tr = (t==0 ? 1 : t); // use left-sided derivative, to avoid the discontinuity at age=0
		const double birth_rate_pc = (constant_birth_rate_pc ? last_birth_rate_pc : birth_rates_pc[t]);
		birth_rates[t] 		= birth_rate_pc*total_diversities[t];
		death_rates[t] 		= -(total_diversities[tr] - total_diversities[tl])/(times[tr] - times[tl]) + birth_rates[t];
		death_rates_pc[t] 	= death_rates[t]/(0.5*(total_diversities[tr] + total_diversities[tl]));
		pulled_diversification_rates[t] = (log(-coalescent_diversities[tr]/nu[tr]) - log(-coalescent_diversities[tl]/nu[tl]))/(times[tr] - times[tl]);
		pulled_extinction_rates[t]		= last_birth_rate_pc - pulled_diversification_rates[t];
	}
	
	// calculate Psurvival (based on estimate per-capita birth & death rates)
	std::vector<double> Psurvival(NT);
	Psurvival[NT-1] = 1;
	for(long t=NT-2; t>=0; --t){
		const double dt = times[t+1]-times[t];
		const double birth_rate_pc = (constant_birth_rate_pc ? last_birth_rate_pc : birth_rates_pc[t]);
		// 2-step Runge-Kutta
		double rate1 = - (death_rates_pc[t] - (1-Psurvival[t+1])*(birth_rate_pc+death_rates_pc[t]) + pow((1-Psurvival[t+1]),Nsplits)*birth_rate_pc);
		double temp_Psurvival = Psurvival[t+1] + dt*rate1;
		double rate2 = - (death_rates_pc[t] - (1-temp_Psurvival)*(birth_rate_pc+death_rates_pc[t]) + pow((1-temp_Psurvival),Nsplits)*birth_rate_pc);
		Psurvival[t] = max(0.0, min(1.0, Psurvival[t+1] + dt*0.5*(rate1+rate2)));
	}
	
	
	// calculate total number of births & deaths
	const double total_births = integrate1D(times,birth_rates,first_time_point,NT-1,true);
	const double total_deaths = integrate1D(times,death_rates,first_time_point,NT-1,true);
	
	return Rcpp::List::create(	Rcpp::Named("success")						= true,
								Rcpp::Named("total_diversities")			= Rcpp::wrap(total_diversities),
								Rcpp::Named("birth_rates") 					= Rcpp::wrap(birth_rates),
								Rcpp::Named("death_rates") 					= Rcpp::wrap(death_rates),
								Rcpp::Named("Prepresentation")				= Rcpp::wrap(Prepresentation),
								Rcpp::Named("Psurvival")					= Rcpp::wrap(Psurvival),
								Rcpp::Named("total_births")					= total_births,
								Rcpp::Named("total_deaths")					= total_deaths,
								Rcpp::Named("last_birth_rate_pc")			= last_birth_rate_pc,
								Rcpp::Named("last_death_rate_pc")			= death_rates[NT-1]/total_diversities[NT-1],
								Rcpp::Named("pulled_diversification_rates")	= Rcpp::wrap(pulled_diversification_rates),
								Rcpp::Named("pulled_extinction_rates")		= Rcpp::wrap(pulled_extinction_rates),
								Rcpp::Named("pulled_total_diversities")		= Rcpp::wrap(pulled_total_diversities));
}






// Estimate past diversity, birth & death rates, based on a time series of coalescent diversities
// This reconstruction is non-parametric, i.e. no particular model is assumed (except for knowledge or constancy of the per-capita birth rate)
// Similar to reconstruct_past_diversity_from_biased_coalescent_CPP, but this one allows specification of an arbitrary discovery_fraction for each time, as opposed to assuming random unbiased sampling at the end
// Note: This function is currently only implemented for bifurcating trees. 
// Input:
//	A time series of coalescent diversities, i.e. as would be visible in a coalescent phylogenetic tree.
// 	Corresponding assumed per-capita birth rates. Alternatively, these can be assumed to be constant and be estimated directly from the coalescent_diversity time series.
//  A time series of discovery_fractions, at each age specifying the fraction of lineages at that age, with extant discovered descendants
// Output:
//  The estimated true past diversities (N(t))
//  The estimated corresponding death (extinction) rates (delta(t))
//  The probability of a size-1 clade surviving from each time point to the present (P(t))
//  The estimated total number of speciation & extinction events, over the considered age interval
// [[Rcpp::export]]
Rcpp::List reconstruct_past_diversity_from_biased_coalescent_CPP(	const std::vector<double>	&times,						// (INPUT) 1D numeric array of size NT, listing times in ascending order. The last time point corresponds to the "present" (age=0)
																	const std::vector<double>	&raw_coalescent_diversities,// (INPUT) 1D numeric array of size NT, listing coalescent diversities as visible after rarefaction. Should be unsmoothened.
																	const std::vector<double> 	&birth_rates_pc,			// (INPUT) 1D numeric array of size NT, listing known or assumed per-capita birth rates. Can also be of size 1, in which case the same per-capita birth rate is assumed throughout. Can also be empty, in which case a constant per-capita birth rate is assumed and estimated from the last slope of the coalescent_diversity curve.
																	const std::vector<double> 	&discovery_fractions,		// (INPUT) 1D numeric array of size NT, listing discovery fractions for each age and synchronized with times[]. For example, discovery_fractions.back() corresponds to the fraction of discovered extant species.
																	const std::vector<double> 	&discovery_fraction_slopes,	// (INPUT) 1D numeric array of size NT, listing the 1st derivative of the discovery_fractions (w.r.t. time) at times[]
																	const double 				max_age,					// (INPUT) max age (distance from last time point) to consider for integrating total births & deaths. If <=0, all times are considered.
																	const long					smoothing_span,				// (INPUT) Integer. Optional number of time points for smoothening the diversities time series via Savitzky-Golay-filter. If <=2, no smoothing is done. Smoothening the coalescent diversity can reduce the noise in the non-parametric reconstruction. 
																	const long					smoothing_order){			// (INPUT) Integer. Optional polynomial order of the smoothing model.
	const long NT 			= times.size();
	const double max_time 	= times[NT-1];
	
	// smoothen time series if needed
	// using a Savitzky-Golay-filter of 2nd-order (1st order if smoothen==3)
	const bool smoothen = (smoothing_span>2);
	std::vector<double> smooth_coalescent_diversities(NT);
	if(smoothen){
		if(!smoothenTimeSeriesSavitzkyGolay(times, raw_coalescent_diversities, 0.0, smoothing_span, min(smoothing_span-2,smoothing_order), true, smooth_coalescent_diversities)){
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
		// correct for the effects of incomplete discovery
		last_birth_rate_pc = max(0.0, last_birth_rate_pc*discovery_fractions[NT-1] - discovery_fraction_slopes[NT-1]);
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
		survival_chances[t] = -(discovery_fractions[t]*nu[t] + discovery_fraction_slopes[t])/birth_rate_pc;
		if((survival_chances[t]<=1e-2/coalescent_diversities[NT-1]) && (coalescent_diversities[t]>0)){
			// diversity not defined (survival seems spuriously low)
			true_diversities[t] = NAN_D;
		}else if(survival_chances[t]>1){
			// diversity not defined (survival chance > 1)
			true_diversities[t] = NAN_D;			
		}else{
			true_diversities[t] = coalescent_diversities[t]/(survival_chances[t] * discovery_fractions[t]);
		}
	}
	true_diversities[NT-1] = coalescent_diversities[NT-1]/discovery_fractions[NT-1]; // coalescent diversity last time point is assumed to be equal to true diversity multiplied by discovery fraction (since survival_chance = 1)
	
	// calculate birth & death rates & per-capita (exponential) pulled growth rate
	std::vector<double> death_rates(NT), birth_rates(NT), pulled_diversification_rates(NT);
	for(long t=NT-1, tl, tr; t>=0; --t){
		tl = (t==0 ? 0 : t-1);
		tr = (t==0 ? 1 : t); // use left-sided derivative, to avoid the discontinuity at age=0
		const double birth_rate_pc = (constant_birth_rate_pc ? last_birth_rate_pc : birth_rates_pc[t]);
		birth_rates[t] = birth_rate_pc*true_diversities[t];
		death_rates[t] = max(0.0, -(true_diversities[tr] - true_diversities[tl])/(times[tr] - times[tl]) + birth_rates[t]); // death rate must be non-negative
		pulled_diversification_rates[t] = (log(-coalescent_diversities[tr]/(nu[tr]*SQ(discovery_fractions[tr]) - discovery_fractions[tr]*discovery_fraction_slopes[tr])) - log(-coalescent_diversities[tl]/(nu[tl]*SQ(discovery_fractions[tl]) - discovery_fractions[tl]*discovery_fraction_slopes[tl])))/(times[tr] - times[tl]);
	}
	
	// calculate total number of births & deaths
	const double total_births = integrate1D(times,birth_rates,first_time_point,NT-1,true);
	const double total_deaths = integrate1D(times,death_rates,first_time_point,NT-1,true);
	
	return Rcpp::List::create(	Rcpp::Named("success")						= true,
								Rcpp::Named("true_diversities")				= Rcpp::wrap(true_diversities),
								Rcpp::Named("birth_rates") 					= Rcpp::wrap(birth_rates),
								Rcpp::Named("death_rates") 					= Rcpp::wrap(death_rates),
								Rcpp::Named("Psurvival")					= Rcpp::wrap(survival_chances),
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
		if(!smoothenTimeSeriesSavitzkyGolay(times, raw_diversities, 0.0, smoothing_span, min(smoothing_span-2,smoothing_order), true, smooth_diversities)){
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
	
	// integrate backwards to calculate probability of survival & representation
	// probability of representation at age=0 begins at (1-rarefaction), probability of survival at age=0 begins at 1
	std::vector<double> Psurvival(NT), Prepresentation(NT);
	Prepresentation[NT-1] = rarefaction;
	Psurvival[NT-1] = 1;
	for(long t=NT-2; t>=0; --t){
		const double dt = times[t+1]-times[t];
		const double birth_rate_pc = (const_birth_rate_pc ? birth_rates_pc[0] : birth_rates_pc[t]);
		const double death_rate_pc = 0.5*(death_rates[t+1]/true_diversities[t+1] + death_rates[t]/true_diversities[t]);
		Psurvival[t] = Psurvival[t+1] - dt*(death_rate_pc - (1-Psurvival[t+1])*(birth_rate_pc+death_rate_pc) + pow((1-Psurvival[t+1]),Nsplits)*birth_rate_pc);
		Prepresentation[t] = Prepresentation[t+1] - dt*(death_rate_pc - (1-Prepresentation[t+1])*(birth_rate_pc+death_rate_pc) + pow((1-Prepresentation[t+1]),Nsplits)*birth_rate_pc);
	}
	
	// use Prepresentation & Psurvival to calculate coalescent diversity and Pdiscovery
	std::vector<double> coalescent_diversities(NT);
	std::vector<double> Pdiscovery(NT);
	for(long t=0; t<NT; ++t){
		coalescent_diversities[t] = true_diversities[t]*Prepresentation[t];	
		Pdiscovery[t] = Prepresentation[t]/Psurvival[t];
	}

	
	// calculate total number of births & deaths
	const double total_births = integrate1D(times,birth_rates,first_time_point,NT-1,true);
	const double total_deaths = integrate1D(times,death_rates,first_time_point,NT-1,true);
	
	return Rcpp::List::create(	Rcpp::Named("success")					= true,
								Rcpp::Named("birth_rates") 				= Rcpp::wrap(birth_rates),
								Rcpp::Named("death_rates") 				= Rcpp::wrap(death_rates),
								Rcpp::Named("Psurvival")				= Rcpp::wrap(Psurvival),		// probability of a lineage surviving until today
								Rcpp::Named("Pdiscovery")				= Rcpp::wrap(Pdiscovery),		// probability of an extant lineage being discovered
								Rcpp::Named("Prepresentation")			= Rcpp::wrap(Prepresentation),	// probability of a lineage surviving and being discovered (=Psurvival*Pdiscovery)
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


template<class ARRAY_TYPE>
void get_incoming_edge_per_tip(	const long			Ntips,
								const long			Nedges,
								const ARRAY_TYPE	&tree_edge,					// (INPUT) 2D array of size Nedges x 2, in row-major format
								std::vector<long>	&incoming_edge_per_tip){	// (OUTPUT) 1D array of size Ntips, with values in 0,..,Nedges-1. 
	incoming_edge_per_tip.assign(Ntips,-1);
	for(long edge=0, child; edge<Nedges; ++edge){
		child = tree_edge[edge*2+1];
		if(child<Ntips) incoming_edge_per_tip[child] = edge;
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





// determine root of a tree
// Assuming that the tree is connected (when edge directions are ignored), this function will return -1 if the tree is not properly rooted
// Hence, this function can also be used to check if the tree is properly rooted (provided that it is connected)
template<class ARRAY_TYPE>
long get_root_clade(const long			Ntips,
					const long 			Nnodes,
					const long			Nedges,
					const ARRAY_TYPE	&tree_edge){			// (INPUT) 2D array (in row-major format) of size Nedges x 2
	const long Nclades = Ntips+Nnodes;
	std::vector<long> Nparents_per_clade(Nclades,0);
	for(long edge=0; edge<Nedges; ++edge){
		Nparents_per_clade[tree_edge[edge*2+1]] += 1;
	}
	long root = -1;
	for(long c=0; c<Nclades; ++c){
		if(Nparents_per_clade[c]>1) return -1; // found a clade with multiple parents, which cannot be in a rooted tree
		if(Nparents_per_clade[c]==0){
			// found clade with no parents, so this may be root
			if(root>=0) return -1; // multiple roots found, which cannot be
			root = c;
		}
	}
	return root;
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



class tree_traversal{
public:
	bool includes_tips;
	long Ntips, Nnodes, Nedges;
	std::vector<long> queue;
	std::vector<long> node2first_edge, node2last_edge;
	std::vector<long> edge_mapping;
	
	// Constructor/initializer
	template<class ARRAY_TYPE>
	tree_traversal(	const long			_Ntips,
					const long 			_Nnodes,
					const long			_Nedges,
					const long 			root, 							// (INPUT) index of root node, i.e. an integer in 0:(Ntips+Nnodes-1)
					const ARRAY_TYPE	&tree_edge, 					// (INPUT) 2D array (in row-major format) of size Nedges x 2
					const bool			include_tips,					// (INPUT) if true, then tips are included in the returned queue[]. This does not affect the returned arrays node2first_edge[], node2last_edge[], edges[].
					const bool			precalculated_edge_mappings){	// (INPUT) if true, then the edge mapping tables node2first_edge[], node2last_edge[] and edges[] are taken as is. Otherwise, they are calculated from scratch.
		includes_tips = include_tips;
		Ntips  = _Ntips;
		Nnodes = _Nnodes;
		Nedges = _Nedges;
		get_tree_traversal_root_to_tips(Ntips,
										Nnodes,
										Nedges,
										root,
										tree_edge,
										include_tips,
										precalculated_edge_mappings,
										queue,
										node2first_edge,
										node2last_edge,
										edge_mapping,
										false,
										"");
	}
};






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

	// get tree traversal route (root --> tips)
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





// For each node, calculate the mean phylogenetic distance to its descending tips
// Requirements:
//   The tree must be rooted; the root should be the unique node with no parent
//   The tree can include multifurcations as well as monofurcations
// [[Rcpp::export]]
NumericVector get_mean_depth_per_node_CPP(	const long			Ntips,
											const long 			Nnodes,
											const long			Nedges,
											const IntegerVector &tree_edge, 	// (INPUT) 2D array (in row-major format) of size Nedges x 2
											const NumericVector &edge_length){ 	// (INPUT) 1D array of size Nedges, or an empty std::vector (all branches have length 1)

	// determine parent clade for each clade
	std::vector<long> clade2parent;
	get_parent_per_clade(Ntips, Nnodes, Nedges, tree_edge, clade2parent);

	// get incoming edge for each clade
	std::vector<long> incoming_edge_per_clade;
	get_incoming_edge_per_clade(Ntips, Nnodes, Nedges, tree_edge, incoming_edge_per_clade);

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
										true, 	// include tips
										false, 	// no precalculated edge mappings
										traversal_queue,
										traversal_node2first_edge,
										traversal_node2last_edge,
										traversal_edges,
										false,
										"");

	// calculate number of descending tips per node & mean distance to tips, traversing tips-->root (excluding the root)
	std::vector<long> node2tip_count(Nnodes,0);
	std::vector<double> node2tip_depth(Nnodes,0);
	for(long q=traversal_queue.size()-1, clade, parent, cnode; q>=1; --q){
		clade	= traversal_queue[q];
		cnode	= clade-Ntips;
		parent 	= clade2parent[clade];
		node2tip_count[parent-Ntips] += (clade<Ntips ? 1 : node2tip_count[cnode]);
		node2tip_depth[parent-Ntips] += (clade<Ntips ? 0 : node2tip_depth[cnode]) + (clade<Ntips ? 1 : node2tip_count[cnode]) * (edge_length.size()==0 ? 1 : edge_length[incoming_edge_per_clade[clade]]);
	}
	for(long node=0; node<Nnodes; ++node){
		node2tip_depth[node] /= node2tip_count[node];
	}
	return Rcpp::wrap(node2tip_depth);
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
template<class ARRAY_TYPE_INT, class ARRAY_TYPE_D>
void get_distances_from_root(	const long 				Ntips,
								const long 				Nnodes,
								const long 				Nedges,
								const ARRAY_TYPE_INT 	&tree_edge,		// (INPUT) 2D array of size Nedges x 2 in row-major format
								const ARRAY_TYPE_D		&edge_length, 	// (INPUT) 1D array of size Nedges, or an empty std::vector (all branches have length 1)
								std::vector<double>		&distances){	// (OUTPUT) 1D array of size Nclades, listing the phylogenetic distance of each clade from the root
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






// For each clade (tip & node) in a tree, find the most distant tip (in terms of cumulative branch length).
// Optionally, the search can be restricted to descending tips.
// Optionally, the search can also be restricted to a subset of target tips.
// If you want distances in terms of branch counts (instead of cumulative branch lengths), simply provide an empty edge_length[].
// Requirements:
//   The input tree must be rooted (root will be determined automatically, as the node that has no incoming edge)
//   The input tree can be multifurcating and/or monofurcating
// [[Rcpp::export]]
Rcpp::List get_farthest_tip_per_clade_CPP(	const long 			Ntips,
											const long 			Nnodes,
											const long 			Nedges,
											const IntegerVector &tree_edge,				// 2D array of size Nedges x 2 in row-major format
											const NumericVector &edge_length, 			// 1D array of size Nedges, or an empty std::vector (all branches have length 1)
											const IntegerVector	&onlyToTips,			// 1D array listing target tips to restrict search to, or an empty std::vector (consider all tips as targets)
											bool				only_descending_tips,	// if true, then for each clade only descending tips are considered for farthest-distance. If false, some clades may have non-descending tips assigned as farthest tips.
											bool 				verbose,
											const std::string	&verbose_prefix){
	const long Nclades = Ntips + Nnodes;
	long parent, clade, tip, incoming_edge;
	double candidate_distance;
	const bool unit_edge_lengths = (edge_length.size()==0);

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

	// Step 1: calculate farthest descending tip per clade (traverse tips --> root)
	std::vector<long> farthest_descending_tip_per_clade(Nclades,-1);
	std::vector<double> distance_to_farthest_descending_tip_per_clade(Nclades,0);
	if(onlyToTips.size()==0){
		// consider all tips as potential targets
		for(long tip=0; tip<Ntips; ++tip){
			farthest_descending_tip_per_clade[tip] = tip;
		}
	}else{
		// only consider provided tips as targets
		for(long t=0; t<onlyToTips.size(); ++t){
			tip = onlyToTips[t];
			farthest_descending_tip_per_clade[tip] = tip;
		}
	}
	for(long q=traversal_queue_root2tips.size()-1; q>=0; --q){
		clade = traversal_queue_root2tips[q];
		if(clade==root) continue;
		if(farthest_descending_tip_per_clade[clade]<0) continue; // no descending tip available from this clade
		parent			= clade2parent[clade];
		incoming_edge	= incoming_edge_per_clade[clade];
		// propagate information about farthest descending tip, to parent (if more distant than already saved for the parent)
		candidate_distance = (unit_edge_lengths ? 1.0 : edge_length[incoming_edge]) + distance_to_farthest_descending_tip_per_clade[clade];
		if((candidate_distance>distance_to_farthest_descending_tip_per_clade[parent]) || (farthest_descending_tip_per_clade[parent]<0)){
			distance_to_farthest_descending_tip_per_clade[parent] = candidate_distance;
			farthest_descending_tip_per_clade[parent] = farthest_descending_tip_per_clade[clade];
		}
	}
	
	if(only_descending_tips){
		// only descending tips allowed, so we're finished
		return Rcpp::List::create(	Rcpp::Named("farthest_tips") 		= Rcpp::wrap(farthest_descending_tip_per_clade),
									Rcpp::Named("farthest_distances") 	= Rcpp::wrap(distance_to_farthest_descending_tip_per_clade));
	}
	

	// Step 2: calculate farthest upstream tip per clade
	std::vector<long> farthest_upstream_tip_per_clade(Nclades,-1);
	std::vector<double> distance_to_farthest_upstream_tip_per_clade(Nclades,0);
	for(long q=1; q<traversal_queue_root2tips.size(); ++q){
		clade	= traversal_queue_root2tips[q];
		parent	= clade2parent[clade];
		incoming_edge = incoming_edge_per_clade[clade];
		for(long e=traversal_node2first_edge[parent-Ntips], edge, child; e<=traversal_node2last_edge[parent-Ntips]; ++e){
			edge = traversal_edges[e];
			if(edge==incoming_edge) continue;
			child = tree_edge[2*edge+1];
			if(farthest_descending_tip_per_clade[child]<0) continue;
			candidate_distance = (unit_edge_lengths ? 1.0+1.0 : edge_length[edge]+edge_length[incoming_edge]) + distance_to_farthest_descending_tip_per_clade[child];
			if((candidate_distance>distance_to_farthest_upstream_tip_per_clade[clade]) || (farthest_upstream_tip_per_clade[clade]<0)){
				farthest_upstream_tip_per_clade[clade] = farthest_descending_tip_per_clade[child];
				distance_to_farthest_upstream_tip_per_clade[clade] = candidate_distance;
			}
		}
		// check if going further up than the parrent leads to an even farther target tip
		if(farthest_upstream_tip_per_clade[parent]>=0){
			candidate_distance = (unit_edge_lengths ? 1.0 : edge_length[incoming_edge]) + distance_to_farthest_upstream_tip_per_clade[parent];
			if((candidate_distance>distance_to_farthest_upstream_tip_per_clade[clade]) || (farthest_upstream_tip_per_clade[clade]<0)){
				farthest_upstream_tip_per_clade[clade] = farthest_upstream_tip_per_clade[parent];
				distance_to_farthest_upstream_tip_per_clade[clade] = candidate_distance;
			}
		}
	}

	// Step 3: calculate farthest tip per clade, regardless of whether descending or not (traverse root --> tips)
	std::vector<long> farthest_tip_per_clade(Nclades,-1);
	std::vector<double> distance_to_farthest_tip_per_clade(Nclades,0);
	for(long q=0; q<traversal_queue_root2tips.size(); ++q){
		clade = traversal_queue_root2tips[q];
		if((farthest_upstream_tip_per_clade[clade]<0) || (distance_to_farthest_descending_tip_per_clade[clade]>distance_to_farthest_upstream_tip_per_clade[clade])){
			// farthest tip for this clade is downstream
			farthest_tip_per_clade[clade] = farthest_descending_tip_per_clade[clade];
			distance_to_farthest_tip_per_clade[clade] = distance_to_farthest_descending_tip_per_clade[clade];
		}else if((farthest_descending_tip_per_clade[clade]<0) || (distance_to_farthest_descending_tip_per_clade[clade]<distance_to_farthest_upstream_tip_per_clade[clade])){
			// farthest tip for this clade is upstream
			farthest_tip_per_clade[clade] = farthest_upstream_tip_per_clade[clade];
			distance_to_farthest_tip_per_clade[clade] = distance_to_farthest_upstream_tip_per_clade[clade];
		}
	}
	
	return Rcpp::List::create(	Rcpp::Named("farthest_tips") 		= Rcpp::wrap(farthest_tip_per_clade),
								Rcpp::Named("farthest_distances") 	= Rcpp::wrap(distance_to_farthest_tip_per_clade));
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
												const IntegerVector	&tree_edge,			// (INPUT) 2D array of size Nedges x 2, flattened in row-major format
												const NumericVector	&edge_length, 		// (INPUT) 1D array of size Nedges, or an empty std::vector (all branches have length 1)
												const long			Ntimes,				// (INPUT) number of time points
												double				min_time,			// (INPUT) minimum time (distance from root) to consider. If negative, will be set to the minimum possible.
												double				max_time,			// (INPUT) maximum time (distance from root) to consider. If Infinite, will be set to the maximum possible.
												const bool			include_slopes){	// (INPUT) if true, slopes of the clades_per_time_point curve are also returned	
	// calculate clade distances from root
	const NumericVector clade_times = get_distances_from_root_CPP(Ntips, Nnodes, Nedges, tree_edge, edge_length);
	max_time = min(max_time, get_array_max(clade_times));
	min_time = max(0.0, min_time);
	
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
		const long last_time_point 	= min(Ntimes-1,long(floor((clade_times[child]-min_time)/time_step)));
		if(last_time_point<0) continue; // edge is outside of considered time span
		const long first_time_point = (parent<0 ? last_time_point : max(0L,long(ceil((clade_times[parent]-min_time)/time_step))));
		if(first_time_point>Ntimes-1) continue; // edge is outside of considered time span
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
	get_distances_from_root(Ntips, Nnodes, Nedges, tree_edge, edge_length, distances_from_root);
	
	const long Ntimes = times.size();
	std::vector<long> diversities(Ntimes,0);
	for(long t=0; t<Ntimes; ++t){
		if(times[t]==0){
			diversities[t] = 1; // by convention, only one clade at root
			continue;
		}
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
												const IntegerVector		&tree_edge,			// (INPUT) 2D array (in row-major format) of size Nedges x 2
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

	// get tree traversal route (root --> tips)											
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
	std::vector<long> chronological_clade_order(Nclades);
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



// Calculate relative evolutionary divergences (RED) of nodes, similarly to the PhyloRank v0.0.27 package [Parks et al. 2018]
// The RED of a node is a measure of its relative placement between the root and its descending tips
// Hence, the RED is always between 0 and 1, with the root having an RED of 0 and all tips having an RED of 1. REDs for tips are not returned here, since they will always be 1.
// Requirements:
//   The tree must be rooted; the root should be the unique node with no parent
//   The tree can include multifurcations as well as monofurcations
void get_relative_evolutionary_divergences(	const long				Ntips,
											const long 				Nnodes,
											const long				Nedges,
											const IntegerVector		&tree_edge,			// (INPUT) 2D array (in row-major format) of size Nedges x 2
											const NumericVector		&edge_length,		// (INPUT) 1D array of size Nedges, or empty in which case each edge is interpreted as having length 1)
											std::vector<double>		&REDs){				// (OUTPUT) 1D array of size Nnodes, listing the RED for each node	
	const bool unit_edge_lengths = (edge_length.size()==0);

	// determine parent clade for each clade
	std::vector<long> clade2parent;
	get_parent_per_clade(Ntips, Nnodes, Nedges, tree_edge, clade2parent);

	// find root using the mapping clade2parent
	const long root = get_root_from_clade2parent(Ntips, clade2parent);

	// get incoming edge for each clade
	std::vector<long> incoming_edge_per_clade;
	get_incoming_edge_per_clade(Ntips, Nnodes, Nedges, tree_edge, incoming_edge_per_clade);

	// get tree traversal route (root --> tips), including tips										
	tree_traversal traversal(Ntips, Nnodes, Nedges, root, tree_edge, true, false);
									
	// get mean distance of each node to its descending tips (traversing tips-->root, exclude root)
	std::vector<long> node2tip_count(Nnodes,0);
	std::vector<double> node2tip_depth(Nnodes,0);
	for(long q=traversal.queue.size()-1, clade, pnode, cnode; q>=1; --q){
		clade	= traversal.queue[q];
		cnode	= clade-Ntips;
		pnode 	= clade2parent[clade] - Ntips;
		node2tip_count[pnode] += (clade<Ntips ? 1 : node2tip_count[cnode]);
		node2tip_depth[pnode] += (clade<Ntips ? 0 : node2tip_depth[cnode]) + (clade<Ntips ? 1 : node2tip_count[cnode]) * (unit_edge_lengths ? 1 : edge_length[incoming_edge_per_clade[clade]]);
	}
	for(long node=0; node<Nnodes; ++node){
		node2tip_depth[node] /= node2tip_count[node];
	}
		
	// calculate RED for each node (traverse root --> tips)
	REDs.resize(Nnodes);
	REDs[root-Ntips] = 0;
	for(long q=1, clade, pnode, cnode; q<traversal.queue.size(); ++q){
		clade = traversal.queue[q];
		if(clade<Ntips) continue; // skip tips
		pnode 	= clade2parent[clade] - Ntips;
		cnode	= clade - Ntips;
		const double mean_distance_to_tips  = node2tip_depth[cnode];
		const double incoming_edge_length 	= (unit_edge_lengths ? 1.0 : edge_length[incoming_edge_per_clade[clade]]);
		if((mean_distance_to_tips + incoming_edge_length)==0){
			REDs[cnode] = REDs[pnode];
		}else{
			REDs[cnode] = min(1.0, REDs[pnode] + (incoming_edge_length/(incoming_edge_length+mean_distance_to_tips)) * (1.0-REDs[pnode]));
		}
	}
}


// Rcpp wrapper for the homonymous base function
// Returns relative evolutionary divergence (RED) values for each node in the tree [Parks et al. 2018]
// The tree must be rooted, but may include monofurcations and multifurcations
// [[Rcpp::export]]
NumericVector get_relative_evolutionary_divergences_CPP(const long				Ntips,
														const long 				Nnodes,
														const long				Nedges,
														const IntegerVector		&tree_edge,		// (INPUT) 2D array (in row-major format) of size Nedges x 2
														const NumericVector		&edge_length){	// (INPUT) 1D array of size Nedges, or empty in which case each edge is interpreted as having length 1)
	std::vector<double> REDs;
	get_relative_evolutionary_divergences(Ntips, Nnodes, Nedges, tree_edge, edge_length, REDs);
	
	return Rcpp::wrap(REDs);
}


// Date (make ultrametric) a phylogenetic tree based on relative evolutionary divergences (RED)
// The RED of a node measures its relative placement between the root and its descending tips.
// For each edge, the RED difference between child & parent is used to set the new length of that edge (times some scaling factor to reproduce the anchor age).
// The tree must be rooted, but may include monofurcations and multifurcations
// [[Rcpp::export]]
Rcpp::List date_tree_via_RED_CPP(	const long				Ntips,
									const long 				Nnodes,
									const long				Nedges,
									const IntegerVector		&tree_edge,		// (INPUT) 2D array (in row-major format) of size Nedges x 2
									const NumericVector		&edge_length,	// (INPUT) 1D array of size Nedges, or empty in which case each edge is interpreted as having length 1)
									const long				anchor_node,	// (INPUT) Index of node to be used as anchor, an integer between 0,..,Nnodes-1. Can also be negative, in which case the root is used as anchor.
									const long				anchor_age){	// (INPUT) Age of the anchor node, a positive real number. If anchor_node<0, then this specifies the age of the root.
	// get node REDs
	std::vector<double> node_REDs;
	get_relative_evolutionary_divergences(Ntips, Nnodes, Nedges, tree_edge, edge_length, node_REDs);
	
	// find scaling factor (time_units/RED_units)
	const double anchor_RED = (anchor_node<0 ? 0 : node_REDs[anchor_node]);
	if(anchor_RED==1) return Rcpp::List::create(Rcpp::Named("success") = false, Rcpp::Named("error") = "Anchor is essentially a tip (its relative evolutionary divergence is 1).");

	const double scaling = anchor_age/(1-anchor_RED);
	std::vector<double> edge_times(Nedges);
	for(long edge=0, parent, child; edge<Nedges; ++edge){
		parent = tree_edge[2*edge+0];
		child  = tree_edge[2*edge+1];
		edge_times[edge] = scaling*max(0.0,((child<Ntips ? 1.0 : node_REDs[child-Ntips]) - node_REDs[parent-Ntips]));
	}	

	return Rcpp::List::create(	Rcpp::Named("edge_times") 	= edge_times, 	// the new edge lengths in time units, such that the 
								Rcpp::Named("node_REDs") 	= node_REDs,
								Rcpp::Named("success") 		= true);
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
											long					&new_root,			// (OUTPUT) index of the root in the collapsed tree. In newer implementations, this is actually guaranteed to be = Ntips (as is common convention).
											double					&root_shift){		// (OUTPUT) phylogenetic distance of new root from the old root (will be 0 if the old root is kept)
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





// Given a set of target tips in an unrooted tree, find the node (or node) in the tree such that when that node is made root, the target tips form a monophyletic group whose MRCA is a child of that node (if as_MRCA==false) or whose MRCA is that node (if as_MRCA==true).
// This node can also be defined as follows (if as_MRCA==false): It is the single node, for which exactly one connected edge satisfies "all tips on the other side are targets", and all other connected edges satisfy "all tips on the other side are non-targets".
// Returns -1 on failure, otherwise it will return a clade index
// [[Rcpp::export]]
long find_root_for_monophyletic_clade_CPP(	const long				Ntips,
											const long 				Nnodes,
											const long				Nedges,
											IntegerVector 			tree_edge,		// (INPUT) 2D array (in row-major format) of size Nedges x 2
											const bool				is_rooted,		// (INPUT) if true, the input tree is guaranteed to already be rooted. Otherwise, it will be temporarily rooted internally at some arbitrary node.
											const std::vector<long>	&target_tips,	// (INPUT) 1D array of tip indices, listing target tips to be made monophyletic
											const bool				as_MRCA){		// (INPUT) if true, the MRCA of the target tips is returned, otherwise the parent of the MRCA is returned
	const long Nclades = Ntips+Nnodes;
	long clade, node, parent;
	if(target_tips.empty()) return -1;
			
	// temporarily root tree if needed (for purposes of traversal)
	// all tip/node/edge indices remain the same
	if(!is_rooted){
		root_tree_at_node(	Ntips,
							Nnodes,
							Nedges,
							tree_edge,	// will be modified in-situ
							1);
	}
	
	// determine parent clade for each clade
	std::vector<long> clade2parent;
	get_parent_per_clade(Ntips, Nnodes, Nedges, tree_edge, clade2parent);
		
	// find root using the mapping clade2parent
	const long root = get_root_from_clade2parent(Ntips, clade2parent);

	// get tree traversal route (root --> tips)		
	std::vector<long> traversal_queue_root2tips, traversal_node2first_edge, traversal_node2last_edge, traversal_edges;
	get_tree_traversal_root_to_tips(	Ntips,
										Nnodes,
										Nedges,
										root,
										tree_edge,
										true, 	// include tips
										false, 	// no precalcuated edge mapping
										traversal_queue_root2tips,
										traversal_node2first_edge,
										traversal_node2last_edge,
										traversal_edges,
										false,
										"");
	
	// Step 1: for each clade, determine the number of children with descending (or being) targets and descending (or being) non-targets (traverse tips --> root)
	std::vector<long> Nchildren_with_descending_targets(Nnodes, 0), Nchildren_with_descending_nontargets(Nnodes, 0), Nchildren_per_node(Nnodes,0), Nnonmonofurcating_children_with_descending_targets(Nnodes,0);
	std::vector<bool> tip_is_target(Ntips,false);
	for(long t=0; t<target_tips.size(); ++t){
		tip_is_target[target_tips[t]] = true;
	}
	for(long q=traversal_queue_root2tips.size()-1, cnode, pnode; q>=1; --q){
		clade 	= traversal_queue_root2tips[q];
		parent 	= clade2parent[clade];
		cnode 	= clade-Ntips;
		pnode 	= parent-Ntips;
		Nchildren_per_node[pnode] += 1;
		// propagate information about descending targets & non-targets to parent
		if(clade<Ntips){
			Nchildren_with_descending_targets[pnode] 	+= (tip_is_target[clade] ? 1 : 0);
			Nchildren_with_descending_nontargets[pnode] += (tip_is_target[clade] ? 0 : 1);
			Nnonmonofurcating_children_with_descending_targets[pnode] += (tip_is_target[clade] ? 1 : 0);
		}else{
			Nchildren_with_descending_targets[pnode] 	+= (Nchildren_with_descending_targets[cnode]>0 ? 1 : 0);
			Nchildren_with_descending_nontargets[pnode] += (Nchildren_with_descending_nontargets[cnode]>0 ? 1 : 0);
			if(traversal_node2last_edge[cnode]-traversal_node2first_edge[cnode]>0) Nnonmonofurcating_children_with_descending_targets[pnode] += (Nchildren_with_descending_targets[cnode]>0 ? 1 : 0);
		}
	}	
	
	// Step 2: determine which clades have upstream targets (traverse root --> tips)
	std::vector<bool> clade_has_upstream_targets(Nclades, false), clade_has_upstream_nontargets(Nclades, false);
	for(long q=1, cnode, pnode; q<traversal_queue_root2tips.size(); ++q){
		clade  	= traversal_queue_root2tips[q];
		parent 	= clade2parent[clade];
		cnode 	= clade-Ntips;
		pnode 	= parent-Ntips;
		if(clade_has_upstream_targets[parent]) clade_has_upstream_targets[clade] = true;
		else if((clade<Ntips ? (tip_is_target[clade] ? 1 : 0) : (Nchildren_with_descending_targets[cnode]>0 ? 1 : 0)) < Nchildren_with_descending_targets[pnode]) clade_has_upstream_targets[clade] = true;
		if(clade_has_upstream_nontargets[parent]) clade_has_upstream_nontargets[clade] = true;
		else if((clade<Ntips ? (tip_is_target[clade] ? 0 : 1) : (Nchildren_with_descending_nontargets[cnode]>0 ? 1 : 0)) < Nchildren_with_descending_nontargets[pnode]) clade_has_upstream_nontargets[clade] = true;
	}
	
	if(as_MRCA){
		// Step 3: Find clade for which at most one inout edge has non-targets on the other side, and such that that edge has only non-targets
		// Monofurcations need to be accounted for in special ways, i.e. make sure the new root is the MRCA of the target tips and not further upstream connected via monofurcations
		for(clade=0; clade<Nclades; ++clade){
			node = clade-Ntips;
			const long Nedges_with_targets = (clade_has_upstream_targets[clade] ? 1 : 0) + (clade<Ntips ? 0 : Nchildren_with_descending_targets[node]);
			const long Nedges_with_nontargets = (clade_has_upstream_nontargets[clade] ? 1 : 0) + (clade<Ntips ? 0 : Nchildren_with_descending_nontargets[node]);
			const long Nedges_total = (clade<Ntips ? 0 : Nchildren_per_node[node])+(clade==root ? 0 : 1);
			if(Nedges_with_nontargets>1) continue; // clade has more than one inout edges with non-targets on the other side
			if(Nedges_with_nontargets==0) return clade;
			if((Nedges_with_nontargets==1) && (Nedges_with_targets==Nedges_total-1) && (Nedges_total>2)) return clade; // if clade only has 2 edges, then if it were made root it would be a monofurcation, in which case the proper MRCA would actually be further downstream
		}
	}else{
		// Step 3: Find clade that has exactly one inout edge, on the other side of which are all targets and only targets
		// Monofurcations need to be accounted for in special ways, i.e. make sure the MRCA of the target tips is a direct (not indirect) descendant of the new root
		for(clade=0; clade<Nclades; ++clade){
			node 	= clade-Ntips;
			parent 	= clade2parent[clade];
			if((clade<Ntips) && (tip_is_target[clade])) continue;
			if(clade_has_upstream_targets[clade] && clade_has_upstream_nontargets[clade]) continue; // there's targets as well as non-targets upstream
			if((clade>=Ntips) && ((Nchildren_with_descending_targets[node]+Nchildren_with_descending_nontargets[node])>Nchildren_per_node[node])) continue; // some children include targets as well as non-targets
			const bool monofurcating_parent = (parent<0 ? false : (parent==root ? (traversal_node2last_edge[parent-Ntips]-traversal_node2first_edge[parent-Ntips]==1) : (traversal_node2last_edge[parent-Ntips]-traversal_node2first_edge[parent-Ntips]==0))); // would the parent node become a monofurcation if clade was made root?
			if(clade<Ntips){
				if(clade_has_upstream_targets[clade] && (!monofurcating_parent)) return clade; // all targets are upstream, but make sure that immediately upstream node is not a monofurcation
			}else{
				if(clade_has_upstream_targets[clade] && (Nchildren_with_descending_targets[node]==0) && (!monofurcating_parent)) return clade; // all targets are upstream, but make sure that immediately upstream node is not a monofurcation
				if((!clade_has_upstream_targets[clade]) && (Nchildren_with_descending_targets[node]==1) && (Nnonmonofurcating_children_with_descending_targets[node]==1)) return clade; // targets are on the other side of exactly one descending edge that does not lead to a monofurcation
			}
		}
	}
	return -1;
}








// Given a set of target tips in an unrooted tree, find the "separator" edge in the tree such that all (or most) target tips are on the one side of the edge, and all (or most) non-target tips are on the other side.
// Specifically, for any edge e (with some arbitrary direction) let N_e^u & N_e^d be the number of targets upstream & downstream, respectively, and let M_e^u & M_e^d be the number of non-targets upstream and downstream, respectively.
// Define E(e) := (N_e^u>N_e^d ? M_e^u+N_e^d : M_e^d+N_e^u). Then the "separator" edge is the edge e that minimizes E(e).
// This function can be used to determine the root of an unrooted tree that would make a set of target tips monophyletic or at least "almost monophyletic" (i.e. with minimal deviation from monophyly).
// This function returns an edge index (in 0,..,Nedges-1), or -1 in the case of failure.
// [[Rcpp::export]]
Rcpp::List find_edge_splitting_tree_CPP(const long				Ntips,
										const long 				Nnodes,
										const long				Nedges,
										IntegerVector 			tree_edge,			// (INPUT) 2D array (in row-major format) of size Nedges x 2
										const bool				is_rooted,			// (INPUT) if true, the input tree is guaranteed to already be rooted. Otherwise, it will be temporarily rooted internally at some arbitrary node.
										const std::vector<long>	&target_tips,		// (INPUT) 1D array of tip indices, listing target tips by which to separate the tree. Can contain duplicates (will be ignored)
										const bool				include_misplaced){	// (INPUT) if true, then the misplaced tips (corresponding to the "optimal" splitting) are also returned as lists. This requires some extra computation.
	const long Nclades = Ntips+Nnodes;
	long clade, child, parent;
	if(target_tips.empty()) return Rcpp::List::create(Rcpp::Named("edge") = -1);
			
	// temporarily root tree if needed (for purposes of traversal)
	// all tip/node/edge indices remain the same
	if(!is_rooted){
		root_tree_at_node(	Ntips,
							Nnodes,
							Nedges,
							tree_edge,	// will be modified in-situ
							1);
	}
	
	// prepare auxiliary lookup tables
	std::vector<long> clade2parent, incoming_edge_per_clade;;
	get_parent_per_clade(Ntips, Nnodes, Nedges, tree_edge, clade2parent);
	get_incoming_edge_per_clade(Ntips, Nnodes, Nedges, tree_edge, incoming_edge_per_clade);
		
	// find root using the mapping clade2parent
	const long root = get_root_from_clade2parent(Ntips, clade2parent);

	// get tree traversal route (root --> tips)		
	std::vector<long> traversal_queue_root2tips, traversal_node2first_edge, traversal_node2last_edge, traversal_edges;
	get_tree_traversal_root_to_tips(	Ntips,
										Nnodes,
										Nedges,
										root,
										tree_edge,
										true, 	// include tips
										false, 	// no precalculated edge mapping
										traversal_queue_root2tips,
										traversal_node2first_edge,
										traversal_node2last_edge,
										traversal_edges,
										false,
										"");
	
	// Step 1: for each clade, determine the number of downstream targets and non-targets (traverse tips-->root)
	std::vector<long> Ntargets_downstream_per_clade(Nclades,0), Nnontargets_downstream_per_clade(Nclades,0);
	for(long tip=0; tip<Ntips; ++tip){
		Nnontargets_downstream_per_clade[tip] = 1;
	}
	for(long t=0, tip; t<target_tips.size(); ++t){
		tip = target_tips[t];
		Ntargets_downstream_per_clade[tip] 		= 1;
		Nnontargets_downstream_per_clade[tip] 	= 0;
	}
	for(long q=traversal_queue_root2tips.size()-1; q>=1; --q){
		clade 	= traversal_queue_root2tips[q];
		parent	= clade2parent[clade];
		Ntargets_downstream_per_clade[parent] 	 += Ntargets_downstream_per_clade[clade];
		Nnontargets_downstream_per_clade[parent] += Nnontargets_downstream_per_clade[clade];
	}
	
	// Step 2: for each clade, determine the number of upstream targets and non-targets (traverse root-->tips)
	std::vector<long> Ntargets_upstream_per_clade(Nclades,0), Nnontargets_upstream_per_clade(Nclades,0);
	for(long q=1; q<traversal_queue_root2tips.size(); ++q){
		clade  	= traversal_queue_root2tips[q];
		parent 	= clade2parent[clade];
		Ntargets_upstream_per_clade[clade] 		= Ntargets_upstream_per_clade[parent] + (Ntargets_downstream_per_clade[parent]-Ntargets_downstream_per_clade[clade]);
		Nnontargets_upstream_per_clade[clade] 	= Nnontargets_upstream_per_clade[parent] + (Nnontargets_downstream_per_clade[parent]-Nnontargets_downstream_per_clade[clade]);
	}
	
	// Step 3: find edge with minimal error in target monophyly (i.e. most targets on one side and most non-targets on the other side)
	long best_edge = -1, best_Ntargets_upstream, best_Ntargets_downstream, best_Nmisplaced_targets, best_Nmisplaced_nontargets;
	for(long edge=0; edge<Nedges; ++edge){
		child 	= tree_edge[edge*2+1];
		const long Ntargets_upstream 		= Ntargets_upstream_per_clade[child];
		const long Ntargets_downstream 		= Ntargets_downstream_per_clade[child];
		const long Nnontargets_upstream 	= Nnontargets_upstream_per_clade[child];
		const long Nnontargets_downstream 	= Nnontargets_downstream_per_clade[child];
		const long Nmisplaced_targets		= (Ntargets_upstream>Ntargets_downstream ? Ntargets_downstream : Ntargets_upstream);
		const long Nmisplaced_nontargets	= (Ntargets_upstream>Ntargets_downstream ? Nnontargets_upstream : Nnontargets_downstream);
		if((best_edge<0) || (Nmisplaced_targets+Nmisplaced_nontargets<best_Nmisplaced_targets+best_Nmisplaced_nontargets)){
			best_edge 	= edge;
			best_Nmisplaced_targets		= Nmisplaced_targets;
			best_Nmisplaced_nontargets	= Nmisplaced_nontargets;
			best_Ntargets_upstream 		= Ntargets_upstream;
			best_Ntargets_downstream 	= Ntargets_downstream;
		}
	}
	
	// determine misplaced tips if requested
	std::vector<long> misplaced_targets, misplaced_nontargets;
	if(include_misplaced){
		misplaced_targets.reserve(best_Nmisplaced_targets);
		misplaced_nontargets.reserve(best_Nmisplaced_nontargets);
		
		// determine which clades (especially tips) descend from best_edge (traverse root-->tips)
		std::vector<bool> descends_from_best_edge(Nclades,false);
		descends_from_best_edge[tree_edge[2*best_edge+1]] = true;
		for(long q=1; q<traversal_queue_root2tips.size(); ++q){
			clade 	= traversal_queue_root2tips[q];
			parent 	= clade2parent[clade];
			if(!descends_from_best_edge[clade]) descends_from_best_edge[clade] = descends_from_best_edge[parent];
		}
		
		// collect misplaced target & non-target tips
		const bool targets_should_be_upstream = (best_Ntargets_upstream>best_Ntargets_downstream);
		for(long tip=0; tip<Ntips; ++tip){
			if(descends_from_best_edge[tip] && (Ntargets_downstream_per_clade[tip]==1) && targets_should_be_upstream) misplaced_targets.push_back(tip);					// misplaced downstream target
			else if(descends_from_best_edge[tip] && (Nnontargets_downstream_per_clade[tip]==1) && (!targets_should_be_upstream)) misplaced_nontargets.push_back(tip);	// misplaced downstream non-target
			else if((!descends_from_best_edge[tip]) && (Ntargets_downstream_per_clade[tip]==1) && (!targets_should_be_upstream)) misplaced_targets.push_back(tip);		// misplaced upstream target
			else if((!descends_from_best_edge[tip]) && (Nnontargets_downstream_per_clade[tip]==1) && targets_should_be_upstream) misplaced_nontargets.push_back(tip); 	// misplaced upstream non-target
		}
	}
	
	return Rcpp::List::create(	Rcpp::Named("edge") 					= best_edge,
								Rcpp::Named("Nmisplaced_targets")		= best_Nmisplaced_targets,
								Rcpp::Named("Nmisplaced_nontargets")	= best_Nmisplaced_nontargets,
								Rcpp::Named("Ntargets_upstream") 		= best_Ntargets_upstream,
								Rcpp::Named("Ntargets_downstream") 		= best_Ntargets_downstream,
								Rcpp::Named("misplaced_targets") 		= misplaced_targets,
								Rcpp::Named("misplaced_nontargets") 	= misplaced_nontargets);
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
	
	// determine new root & root shift
	// traverse from new to old root and calculate cumulative distance
	// use old tree structure
	const long new_root = old2new_clade[root];
	clade = root;
	long root_shift = 0;
	while(incoming_edge_per_clade[clade]>=0){
		edge = incoming_edge_per_clade[clade];
		root_shift += (unit_edge_lengths ? 1.0 : edge_length[edge]);
		clade = tree_edge[edge*2 + 0];
	}

	return Rcpp::List::create(	Rcpp::Named("new_tree_edge") 		= Rcpp::wrap(new_tree_edge),
								Rcpp::Named("new_edge_length") 		= Rcpp::wrap(new_edge_length), // only relevant if (shorten==false)
								Rcpp::Named("new2old_clade") 		= Rcpp::wrap(new2old_clade),
								Rcpp::Named("new2old_edge") 		= Rcpp::wrap(new2old_edge),
								Rcpp::Named("old2new_clade") 		= Rcpp::wrap(old2new_clade),
								Rcpp::Named("collapsed_nodes") 		= Rcpp::wrap(collapsed_nodes),
								Rcpp::Named("representative_tips") 	= Rcpp::wrap(representative_tips),
								Rcpp::Named("new_root") 			= new_root, // in newer implementations this is actually guaranteed to be = Ntips_new
								Rcpp::Named("Ntips_new") 			= Ntips_new,
								Rcpp::Named("Nnodes_new") 			= Nnodes_new,
								Rcpp::Named("Nedges_new") 			= Nedges_new,
								Rcpp::Named("root_shift") 			= root_shift);
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
	std::vector<double> new_edge_length(Nedges_new);
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




// extend terminal edges (edges leading to tips) so that each tip has the same fixed distance (new_height) from the root
// if a tip already extends beyond the specified new_height, its incoming edge remains unchanged
// this is a quick-and-dirty way to make the tree ultrametric
// [[Rcpp::export]]
Rcpp::List extend_tree_to_height_CPP(	const long			Ntips,
										const long 			Nnodes,
										const long			Nedges,
										const IntegerVector	&tree_edge,		// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
										const NumericVector &edge_length, 	// (INPUT) 1D array of size Nedges, or an empty std::vector (all branches have length 1)
										double				new_height){	// (INPUT) phylogenetic distance from the root, to which all tips are to be extended. If negative, this is set to the max_distance_to_root of the input tree		
	// get incoming edge per tip (terminal edges)
	std::vector<long> incoming_edge_per_tip;
	get_incoming_edge_per_tip(Ntips, Nedges, tree_edge, incoming_edge_per_tip);
	
	// calculate all current distances from root
	std::vector<double> distances_from_root;
	get_distances_from_root(Ntips, Nnodes, Nedges, tree_edge, edge_length, distances_from_root);
	if(new_height<0) new_height = get_array_max(distances_from_root);
	
	// extend terminal edges to new_height
	double max_extension = 0;
	std::vector<double> new_edge_length(Nedges);
	if(edge_length.size()==0){
		new_edge_length.assign(Nedges,1);
	}else{
		for(long edge=0; edge<Nedges; ++edge) new_edge_length[edge] = edge_length[edge];
	}
	for(long tip=0; tip<Ntips; ++tip){
		const double extension = new_height - distances_from_root[tip];
		if(extension>0) new_edge_length[incoming_edge_per_tip[tip]] += extension;
		max_extension = max(max_extension, extension);
	}

	return Rcpp::List::create(	Rcpp::Named("new_edge_length")	= Rcpp::wrap(new_edge_length),
								Rcpp::Named("max_extension")	= Rcpp::wrap(max_extension)); // max length that was added to any edge
}


// Eliminate multifurcations in a tree by replacing them with multiple descending bifurcations
// Tips indices remain the same, but edge indices and the total number of nodes/edges may increase (if the tree includes multifurcations)
// Old nodes retain their index, and new nodes will have indices Nnodes,...,(Nnew_nodes-1)
// New nodes will always descend from the old multifurcating nodes, that is, for every multifurcating old node that is split into bifurcations, the newly added nodes will descend from the old node (that kept its original index)
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







// assign tips & nodes of a tree to groups, such that each group is monophyletic (a "taxon") represented by exactly one of given representative tips
// this is the "reverse" operation of picking one representative from each taxon, for a given partitioning of tips into taxa
// The tree must be rooted; the root should be the unique node with no parent
void assign_clades_to_taxa(	const long				Ntips,
							const long 				Nnodes,
							const long				Nedges,
							const IntegerVector 	&tree_edge,			// (INPUT) 2D array (in row-major format) of size Nedges x 2
							const std::vector<long>	&representatives,	// (INPUT) 1D array of size NR, each element of which is the index of a tip representing a distinct taxon
							std::vector<long>		&clade2taxon){		// (OUTPUT) 1D array of size Nclades, mapping each tip & node of the tree to one of NR taxa. In particular, tip2taxon[representatives[r]] = r for all r=1,..,NR. Nodes with more than 1 descending representatives (and thus not part of a specific taxon) will have value -1. If NR==0, then all clades will be assigned to value -1. Also clades with ambiguous taxon assignment (as can occur due to multifurcations) will have value -1
	const long Nclades = Ntips+Nnodes;
	const long NR = representatives.size();
	
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
										true,		// include tips in traversal
										false,		// edge mappings are not yet computed
										traversal_queue,
										traversal_node2first_edge,
										traversal_node2last_edge,
										traversal_edges,
										false,
										"");
	
	// traverse tips-->root and keep track of which nodes have a descending representative
	clade2taxon.assign(Nclades,-1);
	std::vector<long> clade2Nrep_tips(Nclades,0);
	std::vector<long> clade2Nrep_children(Nclades,0);
	for(long r=0; r<NR; ++r){
		clade2taxon[representatives[r]] = r;
		clade2Nrep_tips[representatives[r]] = 1;
		clade2Nrep_children[representatives[r]] = 1;
	}
	for(long q=traversal_queue.size()-1, parent, clade; q>=1; --q){
		clade  = traversal_queue[q];
		parent = clade2parent[clade];
		if((clade2taxon[clade]>=0) && (clade2taxon[parent]<0)) clade2taxon[parent] = clade2taxon[clade];
		clade2Nrep_tips[parent] += clade2Nrep_tips[clade];
		if(clade2Nrep_children[clade]>0) clade2Nrep_children[parent] += 1;
	}
	
	
	// traverse root-->tips and assign taxa & status to clades
	// status = -1 means multiple descending representatives
	// status = 0 means an unambiguous taxon assignment
	// status = 1 means that the taxon assignment - although performed - was ambiguous
	for(long q=0, clade; q<traversal_queue.size(); ++q){
		clade = traversal_queue[q];
		if(clade2Nrep_tips[clade]>1){
			// node contains multiple descending representatives, so no taxon can be assigned
			clade2taxon[clade] = -1;
		}else if(clade2Nrep_tips[clade]==1){
			// node contains exactly one descending representative, so keep whatever taxon was assigned
			continue;
		}else if(clade==root){
			// clade is root and contains no descending representative. This case is pathological, and only occurs if NR==0.
			continue;
		}else{
			// this non-root clade contains no descending representatives, so need to assign taxon based on parent
			clade2taxon[clade] = clade2taxon[clade2parent[clade]];
		}
	}
}




// Rcpp wrapper for the homonymous function above
// [[Rcpp::export]]
Rcpp::List assign_clades_to_taxa_CPP(	const long				Ntips,
										const long 				Nnodes,
										const long				Nedges,
										const IntegerVector 	&tree_edge,			// (INPUT) 2D array (in row-major format) of size Nedges x 2
										const std::vector<long>	&representatives){	// (INPUT) 1D array of size NR, each element of which is the index of a tip representing a distinct taxon
	std::vector<long> clade2taxon;
	assign_clades_to_taxa(	Ntips,
							Nnodes,
							Nedges,
							tree_edge,
							representatives,
							clade2taxon);
	return Rcpp::List::create(Rcpp::Named("clade2taxon") = clade2taxon);
}




#pragma mark -
#pragma mark Comparing trees
#pragma mark 


// Congruify trees (match nodes) as described by [Eastman et al (2013). Congruification: support for time scaling large phylogenetic trees. Methods in Ecology and Evolution. 4:688-691]
// This function essentially finds nodes in the "target" tree (T) that are equivalent ("concordant") branching points to nodes in the "reference" tree (R)
// In case multiple T-nodes are concordant to the same R-node, preference is given to the one closest to the tips. The same holds for R-nodes.
// This may be useful if R is dated (time-calibrated) and T is to be dated using information on branching times in R.
// [[Rcpp::export]]
Rcpp::List congruify_trees_CPP(	const long				RNtips,
								const long 				RNnodes,
								const long				RNedges,
								const IntegerVector 	&Rtree_edge,
								const long				TNtips,
								const long 				TNnodes,
								const long				TNedges,
								const IntegerVector 	&Ttree_edge,
								const IntegerVector		&mapping){		// 2D array of size NM x 2, in row-major format, mapping a subset of T-tips to a subset of R-tips (mapping[m,0]-->mapping[m,1]). The mapping need not be one-to-one, but T-->R must be a valid mapping. In particular, every T-tip can appear at most once in the mapping.
	const long NM = mapping.size()/2;
	const long RNclades = RNtips + RNnodes;
	const long TNclades = TNtips + TNnodes;
	
	// determine parent clade for each clade
	std::vector<long> Rclade2parent, Tclade2parent;
	get_parent_per_clade(RNtips, RNnodes, RNedges, Rtree_edge, Rclade2parent);
	get_parent_per_clade(TNtips, TNnodes, TNedges, Ttree_edge, Tclade2parent);

	// find root using the mapping clade2parent
	const long Rroot = get_root_from_clade2parent(RNtips, Rclade2parent);
	const long Troot = get_root_from_clade2parent(TNtips, Tclade2parent);

	// get tree traversal route (root --> tips)
	std::vector<long> Rtraversal_queue, Rtraversal_node2first_edge, Rtraversal_node2last_edge, Rtraversal_edges;
	get_tree_traversal_root_to_tips(	RNtips,
										RNnodes,
										RNedges,
										Rroot,
										Rtree_edge,
										true,		// include tips in traversal
										false,		// edge mappings are not yet computed
										Rtraversal_queue,
										Rtraversal_node2first_edge,
										Rtraversal_node2last_edge,
										Rtraversal_edges,
										false,
										"");
	std::vector<long> Ttraversal_queue, Ttraversal_node2first_edge, Ttraversal_node2last_edge, Ttraversal_edges;
	get_tree_traversal_root_to_tips(	TNtips,
										TNnodes,
										TNedges,
										Troot,
										Ttree_edge,
										true,		// include tips in traversal
										false,		// edge mappings are not yet computed
										Ttraversal_queue,
										Ttraversal_node2first_edge,
										Ttraversal_node2last_edge,
										Ttraversal_edges,
										false,
										"");
	
	// create mapping from T-tips to focal tips
	// focals = domain(mapping)
	std::vector<long> Ttip2focal(TNtips,-1); // Ttip2Focal[r] maps the T-tip r to the focal tip index. Can be -1, if T-tip r is not included in the focals (i.e. not in image(mapping))
	long next_focal=0;
	for(long m=0, Ttip; m<NM; ++m){
		Ttip = mapping[2*m+0];
		if(Ttip2focal[Ttip]<0){
			Ttip2focal[Ttip] = next_focal;
			++next_focal;
		}
	}
	const long Nfocals = next_focal;
	
	
	// create membership tables (2D array of size Nclades x Nfocals, in row-major format)
	// memberships[c,f] specifies whether clade c includes (has a descendant) the focal tip f
	std::vector<bool> Rmemberships(RNclades*Nfocals,false), Tmemberships(TNclades*Nfocals,false);
	// set membership of tips included in mapping
	for(long m=0, Rtip, Ttip, Ftip; m<NM; ++m){
		Ttip = mapping[2*m+0];
		Rtip = mapping[2*m+1];
		Ftip = Ttip2focal[Ttip];
		Rmemberships[Nfocals*Rtip+Ftip] = true;
		Tmemberships[Nfocals*Ttip+Ftip] = true;
	}
	// propagate memberships (inclusion of focals) upwards (tips-->root)
	for(long q=Rtraversal_queue.size()-1, parent, clade; q>=1; --q){
		clade  = Rtraversal_queue[q];
		parent = Rclade2parent[clade];
		for(long f=0; f<Nfocals; ++f){
			Rmemberships[Nfocals*parent+f] = (Rmemberships[Nfocals*parent+f] || Rmemberships[Nfocals*clade+f]);
		}
	}
	for(long q=Ttraversal_queue.size()-1, parent, clade; q>=1; --q){
		clade  = Ttraversal_queue[q];
		parent = Tclade2parent[clade];
		for(long f=0; f<Nfocals; ++f){
			Tmemberships[Nfocals*parent+f] = (Tmemberships[Nfocals*parent+f] || Tmemberships[Nfocals*clade+f]);
		}
	}
	
	
	// find equivalent R & T nodes based on membership tables
	// traverse T-tree tips-->roots, so that Tclades closer to the tips are found first (in case of multiple matches)
	std::vector<long> mapped_Tnodes, mapped_Rnodes;
	for(long q=Ttraversal_queue.size()-1, Tclade; q>=0; --q){
		Tclade = Ttraversal_queue[q];
		if(Tclade<TNtips) continue;
		
		// find equivalent node in R-tree, if possible
		// start searching at root, moving towards tips
		// at each branching point, at most one of the children will be a valid next step
		long Rclade = Rroot;
		while(Rclade>=RNtips){
			// requirement "T_in_R": at this point it is guaranteed that Rmemberships[Rclade,:] includes all focal tips that are included in Tmemberships[Tclade,:]
			// At most one of the children of Rclade will still satisfy "T_in_R", so go to that one
			// If a child violates "T_in_R", then all of its descendants also violate "T_in_R", so there is no point in continuing in children violating "T_in_R"
			long Rnode = Rclade-RNtips;
			long next_Rclade = -1;
			for(long e=Rtraversal_node2first_edge[Rnode], Rchild; e<=Rtraversal_node2last_edge[Rnode]; ++e){
				Rchild = Rtree_edge[2*Rtraversal_edges[e]+1];
				bool OK = true;
				for(long f=0; f<Nfocals; ++f){
					if(Tmemberships[Nfocals*Tclade+f] && (!Rmemberships[Nfocals*Rchild+f])){
						OK = false;
						break;
					}
				}
				if(OK){
					next_Rclade = Rchild;
					break;
				}
			}
			if(next_Rclade<0){
				// none of the children of Rclade satisfy requirement "T_in_R"
				// so if there exists an equivalent R-clade to Tclade, it must be Rclade, so check
				bool OK = true;
				for(long f=0; f<Nfocals; ++f){
					if(Tmemberships[Nfocals*Tclade+f] != Rmemberships[Nfocals*Rclade+f]){
						OK = false;
						break;
					}
				}
				if(OK){
					// found equivalent clade
					mapped_Tnodes.push_back(Tclade-TNtips);
					mapped_Rnodes.push_back(Rclade-RNtips);
				}
				break; // move to next Tclade, regardless of success
			}else{
				// move one step closer to tips
				Rclade = next_Rclade;
			}
		}
		if(q%100==0) Rcpp::checkUserInterrupt(); // abort if the user has interrupted the calling R program
	}

	// remove duplicate mapped Tclades
	std::vector<long> mapped_Tnodes_deduplicated, mapped_Rnodes_deduplicated;
	std::vector<bool> Rincluded(RNnodes);
	mapped_Tnodes_deduplicated.reserve(mapped_Tnodes.size());
	mapped_Rnodes_deduplicated.reserve(mapped_Rnodes.size());
	for(long m=0, Tnode, Rnode; m<mapped_Tnodes.size(); ++m){
		Tnode = mapped_Tnodes[m];
		Rnode = mapped_Rnodes[m];
		if(!Rincluded[Rnode]){
			mapped_Tnodes_deduplicated.push_back(Tnode);
			mapped_Rnodes_deduplicated.push_back(mapped_Rnodes[m]);
			Rincluded[Rnode] = true;
		}
	}

	return Rcpp::List::create(	Rcpp::Named("mapped_Tnodes") = mapped_Tnodes_deduplicated,
								Rcpp::Named("mapped_Rnodes") = mapped_Rnodes_deduplicated);
}




// Match nodes from one tree to another, assuming that the tree topologies are the same (but indexed differently) and that both have the same root
// This may be useful if nodes and/or tips were re-indexed, and the only way to match old to new nodes is based on topology (e.g. node names are missing)
// This function returns an error if the trees don't have equivalent topologies, so it can also be used as a simple equivalence test
// If you are dealing with different trees, consider using congruify_trees_CPP(..)
// [[Rcpp::export]]
Rcpp::List match_tree_nodes_CPP(const long				Ntips,
								const long 				Nnodes,
								const long				Nedges,
								const IntegerVector 	&tree_edgeA,
								const IntegerVector 	&tree_edgeB,
								const IntegerVector		&tipsA2B){		// 1D array of size Ntips, mapping A-tip indices to B-tip indices (tipsA2B[a] is the B-tip corresponding to a-th A-tip)
	const long Nclades = Ntips + Nnodes;
	
	// determine parent clade for each clade
	std::vector<long> clade2parentA, clade2parentB;
	get_parent_per_clade(Ntips, Nnodes, Nedges, tree_edgeA, clade2parentA);
	get_parent_per_clade(Ntips, Nnodes, Nedges, tree_edgeB, clade2parentB);

	// find root using the mapping clade2parent
	const long rootA = get_root_from_clade2parent(Ntips, clade2parentA);
	const long rootB = get_root_from_clade2parent(Ntips, clade2parentB);

	// get tree traversal route (root --> tips)
	std::vector<long> traversal_queueA, traversal_node2first_edgeA, traversal_node2last_edgeA, traversal_edgesA;
	get_tree_traversal_root_to_tips(	Ntips,
										Nnodes,
										Nedges,
										rootA,
										tree_edgeA,
										true,		// include tips in traversal
										false,		// edge mappings are not yet computed
										traversal_queueA,
										traversal_node2first_edgeA,
										traversal_node2last_edgeA,
										traversal_edgesA,
										false,
										"");

	// map clades A-->B
	// traverse tips-->root and propagate information from child to parent
	std::vector<long> cladesA2B(Nclades,-1);
	std::vector<bool> matchedB(Nnodes,false);
	long Nmatched = 0;
	for(long tip=0; tip<Ntips; ++tip) cladesA2B[tip] = tipsA2B[tip];
	for(long q=traversal_queueA.size()-1, cladeA, parentA, cladeB, parentB; q>=1; --q){
		cladeA 	= traversal_queueA[q];
		parentA = clade2parentA[cladeA];
		if(cladesA2B[parentA]>=0) continue; // parentA already mapped, so skip
		// assuming cladeA is already mapped to some B-clade, map parentA to the B-clade's parent.
		cladeB = cladesA2B[cladeA];
		if(cladeB==rootB){
			// something went wrong (non-rootA mapped to rootB. This can only happen if the two trees are rooted differently
			return Rcpp::List::create(Rcpp::Named("success") = false, Rcpp::Named("error") = "Tree roots don't match");
		}else{
			parentB = clade2parentB[cladeB];
			cladesA2B[parentA] = parentB;
			++Nmatched;
			if(matchedB[parentB-Ntips]) return Rcpp::List::create(Rcpp::Named("success") = false, Rcpp::Named("error") = "Some nodes in tree B were matched more than once");
			matchedB[parentB-Ntips] = true;
		}
	}
	if(Nmatched<Nnodes) return Rcpp::List::create(Rcpp::Named("success") = false, Rcpp::Named("error") = "Some nodes in tree A could not be matched");
		
	// extract mapped nodes
	std::vector<long> nodesA2B(Nnodes);
	for(long node=0; node<Nnodes; ++node) nodesA2B[node] = cladesA2B[Ntips+node]-Ntips;

	return Rcpp::List::create(	Rcpp::Named("success") 	= true,
								Rcpp::Named("nodesA2B") = nodesA2B,
								Rcpp::Named("rootA")	= rootA,
								Rcpp::Named("rootB")	= rootB);
}





// Calculate Robinson-Foulds distance between two rooted trees
// The trees must share the same tips, but may exhibit different topologies
// This metric quantifies the disagreement in topologies, but does not take into account edge lengths
// [William Day (1985). Optimal algorithms for comparing trees with labeled leaves]
// [[Rcpp::export]]
Rcpp::List get_Robinson_Foulds_distance_CPP(const long				Ntips,
											const long 				NnodesA,
											const long				NedgesA,
											const IntegerVector 	&tree_edgeA,
											const long 				NnodesB,
											const long				NedgesB,
											const IntegerVector 	&tree_edgeB,
											const IntegerVector		&tipsA2B){		// 1D array of size Ntips, mapping A-tip indices to B-tip indices (tipsA2B[a] is the B-tip corresponding to a-th A-tip)
	const long NcladesA = Ntips + NnodesA;
	const long NcladesB = Ntips + NnodesB;
	
	// determine parent clade for each clade
	std::vector<long> clade2parentA, clade2parentB;
	get_parent_per_clade(Ntips, NnodesA, NedgesA, tree_edgeA, clade2parentA);
	get_parent_per_clade(Ntips, NnodesB, NedgesB, tree_edgeB, clade2parentB);

	// find root using the mapping clade2parent
	const long rootA = get_root_from_clade2parent(Ntips, clade2parentA);
	const long rootB = get_root_from_clade2parent(Ntips, clade2parentB);

	// get tree traversal route (root --> tips) in depth-first-search mode (DFS is important, to ensure a certain traversal order of tips and nodes)
	std::vector<long> traversal_queueA, traversal_node2first_edgeA, traversal_node2last_edgeA, traversal_edgesA;
	get_tree_traversal_depth_first_search(	Ntips,
											NnodesA,
											NedgesA,
											rootA,
											tree_edgeA,
											true,		// include tips in traversal
											false,		// edge mappings are not yet computed
											traversal_queueA,
											traversal_node2first_edgeA,
											traversal_node2last_edgeA,
											traversal_edgesA);
	std::vector<long> traversal_queueB, traversal_node2first_edgeB, traversal_node2last_edgeB, traversal_edgesB;
	get_tree_traversal_depth_first_search(	Ntips,
											NnodesB,
											NedgesB,
											rootB,
											tree_edgeB,
											true,		// include tips in traversal
											false,		// edge mappings are not yet computed
											traversal_queueB,
											traversal_node2first_edgeB,
											traversal_node2last_edgeB,
											traversal_edgesB);
										
	
	// create mappings of tipsA-->focal_tips and tipsB-->focal_tips (where focal tips are just a re-indexing of tips, so that they are in the same order as reached by the DFS tree traversal of treeA)
	// requirement: A2F[a] = B2F[tipsA2B[a]]
	std::vector<long> tipsA2F(Ntips,-1), tipsB2F(Ntips,-1);
	long next_focal = 0;
	for(long q=traversal_queueA.size()-1, clade; q>=0; --q){
		clade = traversal_queueA[q];
		if(clade<Ntips) tipsA2F[clade] = (next_focal++);
	}
	for(long tipA=0; tipA<Ntips; ++tipA){
		tipsB2F[tipsA2B[tipA]] = tipsA2F[tipA];
	}
	

	//count the number of tips descending from each clade (traverse tips-->root)
	std::vector<long> cladeA2tip_counts(NcladesA,0), cladeB2tip_counts(NcladesB,0);
	for(long tip=0; tip<Ntips; ++tip){
		cladeA2tip_counts[tip] = 1;
		cladeB2tip_counts[tip] = 1;
	}
	for(long q=traversal_queueA.size()-1, cladeA; q>=1; --q){
		cladeA = traversal_queueA[q];
		cladeA2tip_counts[clade2parentA[cladeA]] += cladeA2tip_counts[cladeA];
	}
	for(long q=traversal_queueB.size()-1, cladeB; q>=1; --q){
		cladeB = traversal_queueB[q];
		cladeB2tip_counts[clade2parentB[cladeB]] += cladeB2tip_counts[cladeB];
	}
		
	// create membership tables (list of focal tips descending from each node)
	// A memberships will be sorted (i.e. each membershipsA[a][] will be a list of ascending focal tip indices)
	// this is achieved because we're traversing the tree in reverse-depth-first-search, and A tips are mapped to focal tips in ascending order, and node traversal is consistent with tip order
	std::vector<std::vector<long> > membershipsA(NnodesA), membershipsB(NnodesB);
	for(long node=0; node<NnodesA; ++node) membershipsA[node].reserve(cladeA2tip_counts[node+Ntips]); // preallocate space
	for(long node=0; node<NnodesB; ++node) membershipsB[node].reserve(cladeB2tip_counts[node+Ntips]); // preallocate space
	for(long q=traversal_queueA.size()-1, clade, cnode, pnode; q>=1; --q){
		clade = traversal_queueA[q];
		pnode = clade2parentA[clade] - Ntips;
		if(clade<Ntips){
			membershipsA[pnode].push_back(tipsA2F[clade]);
		}else{
			cnode = clade-Ntips;
			membershipsA[pnode].insert(membershipsA[pnode].end(), membershipsA[cnode].begin(), membershipsA[cnode].end());
		}
		if(q%100==0) Rcpp::checkUserInterrupt(); // abort if the user has interrupted the calling R program
	}
	for(long q=traversal_queueB.size()-1, clade, cnode, pnode; q>=1; --q){
		clade = traversal_queueB[q];
		pnode = clade2parentB[clade] - Ntips;
		if(clade<Ntips){
			membershipsB[pnode].push_back(tipsB2F[clade]);
		}else{
			cnode = clade-Ntips;
			membershipsB[pnode].insert(membershipsB[pnode].end(), membershipsB[cnode].begin(), membershipsB[cnode].end());
		}
		if(q%100==0) Rcpp::checkUserInterrupt(); // abort if the user has interrupted the calling R program
	}

	// also sort B memberships (A memberships are already sorted)
	for(long node=0; node<NnodesB; ++node){
		std::sort(membershipsB[node].begin(), membershipsB[node].end());
	}
	
	
	// find equivalent nodes between the two trees, by traversing tips-->roots
	long Nmatches=0;
	std::vector<long> cladeA2B(NcladesA); // cladeA2B[a] points to a clade in treeB that is fully contained (but not necessarily equal) to A-clade a
	std::vector<bool> matchedB(NcladesB,false); // keep track of B-clades that were matched before, to avoid duplicate matching
	for(long tipA=0; tipA<Ntips; ++tipA) cladeA2B[tipA] = tipsA2B[tipA];
	for(long q=traversal_queueA.size()-1, cladeA, nodeA; q>=0; --q){
		cladeA = traversal_queueA[q];
		if(cladeA<Ntips) continue;
		nodeA = cladeA - Ntips;
		// at this point, each child of cladeA is mapped to a cladeB which it fully contains (i.e. all of whose tips also descends the child)
		// check if any of the mapped cladeBs can be moved upstream and still be contained in cladeA
		for(long e=traversal_node2first_edgeA[nodeA], childA, cladeB, nodeB; e<=traversal_node2last_edgeA[nodeA]; ++e){
			childA = tree_edgeA[2*traversal_edgesA[e]+1];
			cladeB = cladeA2B[childA];
			cladeA2B[cladeA] = cladeB; // since 
			bool OK = true;
			while(OK && (cladeB!=rootB)){
				cladeB 	= clade2parentB[cladeB];
				nodeB	= cladeB-Ntips;
				// membershipsA[nodeA][] and membershipsB[nodeB][] and guaranteed to be sorted in ascending order
				for(long fb=0, fa=-1, f; fb<membershipsB[nodeB].size(); ++fb){
					f  = membershipsB[nodeB][fb];
					fa = find_in_ascending_list(membershipsA[nodeA],f,fa+1);
					if(fa<0){
						OK = false;
						break;
					}
				}
				if(OK) cladeA2B[cladeA] = cladeB; // still OK, so update matched clade
			}
			cladeB = cladeA2B[cladeA];
			if(cladeA2tip_counts[cladeA]==cladeB2tip_counts[cladeB]){
				// found fully matching B-clade (since cladeA2B[cladeA] is contained in cladeA, and actually has the same size)
				if(!matchedB[cladeB]){
					++Nmatches; // found equivalent B-clade that wasn't matched before
					matchedB[cladeB] = true;
				}
				break;
			}
		}
	}

	return Rcpp::List::create(Rcpp::Named("Nmatches") = Nmatches);
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
#pragma mark Plotting trees
#pragma mark 


// calculate the geometric placement of tips & nodes for plotting a tree as a phylogram
// The root is placed on the left end, tips are placed on the right end, edges extend horizontally left to right
// tips y-coordinates of all clades will be within 0 and Ntips
// [[Rcpp::export]]
Rcpp::List get_phylogram_geometry_CPP(	const long			Ntips,
										const long 			Nnodes,
										const long			Nedges,
										IntegerVector 		tree_edge,			// (INPUT) 2D array (in row-major format) of size Nedges x 2
										const NumericVector	&edge_length){		// (INPUT) 1D array of size Nedges, or empty
	const long Nclades = Ntips + Nnodes;
										
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
										true,	// include tips
										false,	// edge mappings are not yet calculated
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
	// (traverse root --> tips, excluding the root)
	std::vector<double> distances_from_root(Nclades);
	distances_from_root[root] = 0;
	for(long q=1, clade, parent; q<traversal_queue.size(); ++q){
		clade = traversal_queue[q];
		parent = clade2parent[clade];
		distances_from_root[clade] = (edge_length.size()==0 ? 1.0 : edge_length[incoming_edge_per_clade[clade]]) + distances_from_root[parent];
	}

	// calculate number of descending tips per node, traversing tips-->root (excluding the root)
	std::vector<long> node2total_tip_count(Nnodes,0);
	for(long q=traversal_queue.size()-1, clade; q>=1; --q){
		clade = traversal_queue[q];
		node2total_tip_count[clade2parent[clade]-Ntips] += (clade<Ntips ? 1 : node2total_tip_count[clade-Ntips]);
	}
	
	// calculate y-intervals of clades (traverse root-->tips)
	std::vector<double> clade_min_y(Nclades), clade_max_y(Nclades);
	clade_min_y[root] = 0;
	clade_max_y[root] = Ntips;
	for(long q=0, clade, node; q<traversal_queue.size(); ++q){
		clade = traversal_queue[q];
		if(clade<Ntips) continue;
		node = clade - Ntips;
		double cumulative_fraction = 0, fraction;
		for(long e=traversal_node2first_edge[node], child; e<=traversal_node2last_edge[node]; ++e){
			child 				= tree_edge[2*traversal_edges[e]+1];
			fraction			= (child<Ntips ? 1l : node2total_tip_count[child-Ntips])/double(node2total_tip_count[node]);
			clade_min_y[child] 	= clade_min_y[clade] + cumulative_fraction*(clade_max_y[clade]-clade_min_y[clade]);
			clade_max_y[child] 	= clade_min_y[clade] + (cumulative_fraction+fraction)*(clade_max_y[clade]-clade_min_y[clade]);
			cumulative_fraction += fraction;
		}
	}
	
	// calculate y-coordinates of clades (centers of y-intervals)
	std::vector<double> clade_y(Nclades);
	for(long clade=0; clade<Nclades; ++clade){
		clade_y[clade] = 0.5*(clade_min_y[clade] + clade_max_y[clade]);
	}
										
	return Rcpp::List::create(	Rcpp::Named("clade_x") 	= Rcpp::wrap(distances_from_root),
								Rcpp::Named("clade_y") 	= Rcpp::wrap(clade_y),
								Rcpp::Named("min_x")	= 0.0,
								Rcpp::Named("max_x")	= get_array_max(distances_from_root),
								Rcpp::Named("min_y")	= 0.5,
								Rcpp::Named("max_y")	= Ntips-0.5,
								Rcpp::Named("root_y")	= clade_y[root]);
}



#pragma mark -
#pragma mark Writing/reading trees
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


// auxiliary routine for parsing a single edge in a Newick string
// returns false on failure
bool aux_Newick_extract_next_edge(	const std::string 	&input,
									long 				&pointer,		// (INPUT/OUTPUT) will move towards the left
									string 				&name,			// (OUTPUT) child name. Will be empty ("") if not available
									double				&length,		// (OUTPUT) edge length. Will be NAN_D if not available
									string				&error){		// (OUTPUT) error description in case of failure
	long left = -1, split=-1;
	for(long i=pointer; i>=0; --i){
		if(input[i]==':') split = i;
		if((input[i]=='(') || (input[i]==')') || (input[i]==',')){
			left = i;
			break;
		}
	}
	if(left<0){
		error = "Missing terminal character '(', ')' or ','";
		return false;
	}
	if(left==pointer){
		// no name nor length available
		name = "";
		length = NAN_D;
		return true;
	}
	if(split<0){
		// no length available, interpret whole specifier as name
		name   = input.substr(left+1,pointer-left);
		length = NAN_D;
	}else{
		name   = input.substr(left+1,split-left-1);
		length = string2Double(input.substr(split+1,pointer-split));
	}
	pointer = left;
	return true;
}


// read a phylogenetic tree from a Newick-formatted string
// [[Rcpp::export]]
Rcpp::List read_Newick_string_CPP(	std::string	input,
									const bool	underscores_as_blanks){
	// remove any newline characters
	input.erase(std::remove(input.begin(), input.end(), '\n'), input.end());
	
	// trim any whitespace
	input = trim_whitespace(input);
	
	// replace underscores with blanks if needed
	if(underscores_as_blanks){
		std::replace(input.begin(), input.end(), '_', ' ');
	}
	
	// estimate number of tips, nodes & edges for pre-allocation purposes
	const long estimated_Nclades  = std::count(input.begin(), input.end(), ',') + std::count(input.begin(), input.end(), ')');
	const long estimated_Nedges = estimated_Nclades - 1;
	
	
	// pre-allocate space
	std::vector<std::string> clade_names;
	std::vector<double> edge_length;
	std::vector<long> tree_edge;
	clade_names.reserve(estimated_Nclades);
	edge_length.reserve(estimated_Nedges);
	tree_edge.reserve(2*estimated_Nedges);
	
	// prepare auxiliary data structures
	std::vector<long> clade_stack; // keep track of which node we are currently in. clade_stack[n+1] is a child of clade_stack[n]
	long pointer = input.length()-1;
	if(input[pointer]==';') --pointer;
	std::string error, name;
	double length, root_edge;
	
	// read input left<--right
	while(pointer>=0){
		if(clade_stack.empty() && (!clade_names.empty())){
			return Rcpp::List::create(Rcpp::Named("success") = false, Rcpp::Named("error") = "Tree appears to have multiple roots: Reached top tree level prior to reaching left-end of input string, at position "+makeString(pointer+1));
		}
		if(!aux_Newick_extract_next_edge(input, pointer, name, length, error)){
			return Rcpp::List::create(Rcpp::Named("success") = false, Rcpp::Named("error") = "Invalid child specifier to the left of position "+makeString(pointer+1)+": "+error);
		}
		clade_names.push_back(name);
		if(clade_stack.empty()){
			// clade is root
			root_edge = length;
		}else{
			edge_length.push_back(length);
			tree_edge.push_back(clade_stack.back());
			tree_edge.push_back(clade_names.size()-1);
		}
		if(input[pointer]==')'){
			// moving one level deeper, into a new child
			clade_stack.push_back(clade_names.size()-1);
			--pointer;
		}else if(input[pointer]=='('){
			// finished at this level, moving up to parents
			while((pointer>=0) && (input[pointer]=='(')){
				if(clade_stack.empty()){
					return Rcpp::List::create(Rcpp::Named("success") = false, Rcpp::Named("error") = "Unbalanced parentheses, found redundant closing '(' at position "+makeString(pointer+1));
				}
				clade_stack.pop_back();
				--pointer;
			}
			if((pointer>=0) && (input[pointer]==',')) --pointer;
			else if((pointer>=0) && (input[pointer]==')')) return Rcpp::List::create(Rcpp::Named("success") = false, Rcpp::Named("error") = "Unexpected opening paranthesis ')' at position "+makeString(pointer+1));
		}else{
			// more clades to be extracted at this level
			--pointer;
		}
		if(clade_names.size()%1000==0) Rcpp::checkUserInterrupt(); // abort if the user has interrupted the calling R program
	}
	
	// nothing left to parse, so check if we came back to level 0
	if(!clade_stack.empty()){
		return Rcpp::List::create(Rcpp::Named("success") = false, Rcpp::Named("error") = "Unbalanced parentheses, missing "+makeString(clade_stack.size())+" closing parentheses '(' on the left end");
	}
		
	// re-index clades (tips & nodes) consistent with the phylo format
	const long Nclades = clade_names.size();
	const long Nedges  = edge_length.size();
	std::vector<long> old2new_clade;
	long Ntips, Nnodes;
	reindex_clades(	Nclades,
					Nedges,
					tree_edge,
					true,
					Ntips,
					Nnodes,
					old2new_clade);
	const long root = Ntips;
	for(long edge=0; edge<Nedges; ++edge){
		tree_edge[2*edge+0] = old2new_clade[tree_edge[2*edge+0]];
		tree_edge[2*edge+1] = old2new_clade[tree_edge[2*edge+1]];
	}
	vector<string> tip_names(Ntips), node_names(Nnodes);
	for(long clade=0, new_clade; clade<Nclades; ++clade){
		new_clade = old2new_clade[clade];
		if(new_clade<Ntips)  tip_names[new_clade] = clade_names[clade];
		if(new_clade>=Ntips) node_names[new_clade-Ntips] = clade_names[clade];
	}
	
	return Rcpp::List::create(	Rcpp::Named("Ntips") 		= Ntips,
								Rcpp::Named("Nnodes") 		= Nnodes,
								Rcpp::Named("Nedges") 		= Nedges,
								Rcpp::Named("tip_names") 	= Rcpp::wrap(tip_names),
								Rcpp::Named("node_names") 	= Rcpp::wrap(node_names),
								Rcpp::Named("tree_edge")	= Rcpp::wrap(tree_edge),
								Rcpp::Named("edge_length")	= Rcpp::wrap(edge_length),
								Rcpp::Named("root")			= root,
								Rcpp::Named("root_edge")	= root_edge, // length of dummy "edge" (lacking a parent) leading into root
								Rcpp::Named("success")		= true);
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
		Pvalue = (count_valid_permutations>0 ? count_deeper/double(count_valid_permutations) : NAN_D);
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
Rcpp::List get_empirical_state_frequencies_per_node_CPP(	const long			Ntips,
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
	
	return Rcpp::List::create(Rcpp::Named("frequencies_per_node") = Rcpp::wrap(frequencies_per_node));
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
											const IntegerVector &tree_edge,				// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
											const NumericVector &edge_length,			// (INPUT) 1D array of size Nedges, or an empty std::vector (all edges have length 1)
											const NumericVector	&tip_states,			// (INPUT) 2D array of size Ntips x Ntraits, in row-major format, listing numeric states for each trait at each tip
											const bool			only_bifurcations,		// (INPUT) if true, PICs are only calculated for bifurcating nodes in the input tree, and multifurcations are not expanded.
											const bool			scaled,					// (INPUT) if true, then PICs are rescaled by the square-root of their expected variances (=the edge lengths connecting the compared nodes/tips)
											std::vector<double>	&PICs,					// (OUTPUT) 2D array of size Npics x Ntraits, in row-major format, listing PICs for each trait and for each considered node
											std::vector<double>	&distances,				// (OUTPUT) 1D array of size Npics, listing phylogenetic distances corresponding to the PICs. Under a Brownian motion mode, these are proportional to the variance of each PIC
											std::vector<long>	&PIC_nodes,				// (OUTPUT) 1D array of size Npics, listing node indices for which PICs were calculated. Negative values indicate nodes that were not actually in the original tree, but were created temporarily during expansion of multifurcations
											std::vector<double>	&root_state,			// (OUTPUT) 1D array of size Ntraits, holding the root's globally reconstructed state (X_k sensu Felsenstein)
											std::vector<double>	&root_standard_error){	// (OUTPUT) 1D array of size Ntraits, listing the standard errors for the root's state (under a Brownian motion model) [Garland et al. (1999). An introduction to phylogenetically based statistical methods, with a new method for confidence intervals on ancestral values]
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
	const long local_root = get_root_from_incoming_edge_per_clade(Ntips, local_tree_edge, incoming_edge_per_clade);

	// prepare tree traversal route (root-->tips) and edge mappings
	std::vector<long> traversal_queue, node2first_edge, node2last_edge, edge_mapping;
	get_tree_traversal_root_to_tips(Ntips,
									Nlocal_nodes,
									Nlocal_edges,
									local_root,
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
			if(clade!=local_root) incoming_length_per_clade[clade] = local_edge_length[incoming_edge_per_clade[clade]];
		}
	}else{
		incoming_length_per_clade.assign(incoming_length_per_clade.size(),1);
		incoming_length_per_clade[local_root] = 0;
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
		// note that monofurcating nodes acquire the same state as their child, and their modified incoming_length is the same as their child plus the length of their incoming edge
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
	
	// extract estimated root state & calculate standard error
	// this should come before the scaling further below
	// Standard error formula according to: [Garland et al. (1999). Page 377]
	root_state.resize(Ntraits);
	root_standard_error.assign(Ntraits,0);
	for(trait=0; trait<Ntraits; ++trait){
		root_state[trait] = node_states[(local_root - Ntips)*Ntraits+trait];
		for(long p=0; p<Npics; ++p){
			root_standard_error[trait] += SQ(PICs[p*Ntraits+trait])/distances[p];
		}
		root_standard_error[trait] *= incoming_length_per_clade[local_root] / Npics;
		root_standard_error[trait] = sqrt(root_standard_error[trait]);
	}
	
	// rescale returned PICs if needed
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
		std::vector<double> PICs, distances, root_state, root_standard_error;
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
												PIC_nodes,
												root_state,
												root_standard_error);
		return Rcpp::List::create(	Rcpp::Named("PICs")  				= Rcpp::wrap(PICs),
									Rcpp::Named("distances") 			= Rcpp::wrap(distances),
									Rcpp::Named("nodes") 				= Rcpp::wrap(PIC_nodes),
									Rcpp::Named("root_state") 			= Rcpp::wrap(root_state),
									Rcpp::Named("root_standard_error") 	= Rcpp::wrap(root_standard_error));
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
	std::vector<double> PICs, distances, root_state, root_standard_error;
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
											PIC_nodes,
											root_state,
											root_standard_error);
											
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
Rcpp::List WMPR_ASR_CPP(const long			Ntips,
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

	return Rcpp::List::create(	Rcpp::Named("posterior_probabilities") 	= Rcpp::wrap(posterior_probabilities),
								Rcpp::Named("best_root_cost")			= best_root_cost);
}



// calculate log-likelihood and posterior probability at the root, for a fixed-rates continuous-time Markov model for discrete character evolution
// A major computational bottleneck is the exponentiation of the transition matrix along each edge, i.e. exp(edge_length*transition_matrix)
// Exponentiation of the transition matrix can happen in one of 3 ways:
//    1. By providing all pre-calculated exponentials, via expQ_per_edge. This uses more RAM, but is much faster than the other methods below.
//	  2. By providing a transition_exponentiator object.
// The above options are checked and utilized in the listed order. Whenever the associated arrays of an option are empty, the next option is checked.
void aux_ASR_with_fixed_rates_Markov_model(	const long					Ntips,
											const long 					Nnodes,
											const long					Nedges,
											const long					Nstates,
											const long					root,								// (INPUT) root of the tree, an integer in Ntips,..,Nnodes+Ntips-1
											const IntegerVector			&tree_edge,							// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
											const NumericVector 		&edge_length, 						// (INPUT) 1D array of size Nedges, or an empty std::vector (all branches have length 1)
											const NumericVector			&prior_probabilities_per_tip, 		// (INPUT) 2D array of size Ntips x Nstates, in row-major format, listing prior probability distributions for each tip
											const NumericVector			&prior_probabilities_for_root,		// (INPUT) 1D array of size Nstates, listing prior probability distribution for root, to be used for combining the root's probabilities into a single likelihood. Which prior you use for the root is often a matter of choice.
											const matrix_exponentiator	&transition_exponentiator,			// (INPUT) initialized exponentiator object for transition matrix
											const std::vector<double>	&expQ_per_edge,						// (INPUT) 3D array of size Nedges x Nstates x Nstates, in layer-row-major format, listing the exponentiated transition matrix along each edge. Only relevant if use_precalculated_expQ==true.
											const std::vector<long>		&traversal_queue,					// (INPUT) 1D array of size Nnodes, with values in Ntips:(Nclades-1). Traversal queue root-->tips (not including tips). Generated using the function get_tree_traversal_root_to_tips(include_tips=false).
											const std::vector<long>		&traversal_node2first_edge,			// (INPUT) 1D array of size Nnodes, with values in 0:(Nedges-1). Generated using the function get_tree_traversal_root_to_tips().
											const std::vector<long>		&traversal_node2last_edge,			// (INPUT) 1D array of size Nnodes, with values in 0:(Nedges-1). Generated using the function get_tree_traversal_root_to_tips().
											const std::vector<long>		&traversal_edges,					// (INPUT) 1D array of size Nedges, with values in 0:(Nedges-1). Generated using the function get_tree_traversal_root_to_tips().
											std::vector<double>			&posteriors,						// (OUTPUT) 1D array of size Nnodes x Nstates, listing the posterior likelihoods at each node (rescaled to sum to 1 for each node). This is used both as internal scratch space as well as to return the final posteriors.
											double						&loglikelihood){					// (OUTPUT) log-likelihood for the full tree
	long clade, edge, child, node;
	const bool use_precalculated_expQ = (expQ_per_edge.size()>0);
	const double max_edge_length 	  = (edge_length.size()==0 ? 1.0 : get_array_max(edge_length));
							
	// calculate probability distribution on each node (traverse tips-->root)
	posteriors.assign(Nnodes*Nstates,1.0);
	std::vector<double> Y, expQ;
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
			}else{
				// calculate exponential of transition matrix along edge
				transition_exponentiator.get_exponential((edge_length.size()==0 ? 1.0 : edge_length[edge])/max_edge_length, expQ);
				expQ_pointer = &expQ[0];
			}
						
			// use exponentiated transition matrix to propagate likelihoods ("posteriors") from child to parent
			// posteriors[parent] = product_{child in children} exp(edge_length * Q) * posteriors[child]
			// this corresponds to solving the Kolmogorov backward equation along each child-edge in backward-time direction, with initial condition the likelihoods at the child
			if(child<Ntips) multiply_matrix_with_vector(Nstates, Nstates, expQ_pointer, &prior_probabilities_per_tip[child*Nstates], Y);
			else multiply_matrix_with_vector(Nstates, Nstates, expQ_pointer, &posteriors[(child-Ntips)*Nstates], Y);
			for(long s=0; s<Nstates; ++s) posteriors[node*Nstates+s] *= max(0.0,Y[s]); // factor Y into the posterior likelihoods of this node. Avoid negative values from rounding errors			
		}
						
		// rescale (normalize) clade's posterior likelihoods
		// this is necessary due to scaling issues for very large trees; the non-normalized posterior likelihoods tend to 0 for older nodes
		// note that since the Mk propagator (exp(t*Q)) is linear, rescaling just rescales the overall likelihood (this is corrected for below)
		double S = 0;
		for(long s=0; s<Nstates; ++s) S += posteriors[node*Nstates+s];
		for(long s=0; s<Nstates; ++s) posteriors[node*Nstates+s] /= S;
		
		// incorporate rescaling factor (used to normalize this node's posterior) into tree's loglikelihood
		// if we weren't rescaling each node's posterior, this would not be necessary
		loglikelihood += log(S);
	}
	
	// use root's posterior (combined with it's prior) to calculate the model's loglikelihood
	for(long s=0; s<Nstates; ++s) loglikelihood += log(posteriors[(root-Ntips)*Nstates+s]*prior_probabilities_for_root[s]);
}







// re-root the tree and update the posterior probabilities at each node (since each node's posterior is determined by the posteriors of its children)
// since rerooting at a new node does not change the tip/node/edge indices, nor the total number of inout edges per clade, we just need to update a few tree-access data structures and update the posteriors for the affected nodes.
// Note: This algorithm cannot easily be generalized to rerooting at tips, because this would change the total number of tips & nodes and hence mess up the indexing of tips {0,..,Ntips-1} and nodes {Ntips,...,Ntips+Nnodes-1}.
// Exponentiation of the transition matrix can happen in one of 3 ways:
//    1. By providing all pre-calculated exponentials, via expQ_per_edge. This uses more RAM, but is much faster than the other methods below.
//	  2. By providing a transition_exponentiator object.
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
																const matrix_exponentiator	&transition_exponentiator,			// (INPUT) initialized exponentiator object for transition matrix
																const std::vector<double>	&expQ_per_edge,						// (INPUT) Optional 3D array of size Nedges x Nstates x Nstates, in layer-row-major format, listing the exponentiated transition matrix along each edge. Only relevant if use_precalculated_expQ==true. Can be empty, if exp(Q) is to be calculated using eigendecomposition or polynomials.
																std::vector<long> 			&current_tree_edge,					// (INPUT/OUTPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1). Will be updated after the rerooting.		
																std::vector<long>			&current_incoming_edge_per_clade,	// (INPUT/OUTPUT) 1D array of size Nclades, with elements in 0,..,Nedges-1. Will be updated after the rerooting.
																std::vector<double>			&posteriors){						// (INPUT/OUTPUT) 1D array of size Nnodes x Nstates, listing the posterior probabilities at each node. Should be pre-computed for the current tree, and will be updated after the rerooting.
	const bool use_precalculated_expQ = (expQ_per_edge.size()>0);
	if(new_root==old_root) return; // nothing to do
	
	// reroot (this changes edge directions, but keeps tip/node/edge indices the same)
	reroot_tree_at_node(Ntips, Nnodes, Nedges, old_root, new_root, current_tree_edge, current_incoming_edge_per_clade);
	
	// update posteriors of nodes that have been traversed by the rerooting
	std::vector<double> Y, expQ;
	double const *expQ_pointer; // will point to the location of the exponentiated transition matrix (which may be different for each edge)
	long clade = old_root;
	while(true){
		// re-calculate posterior for the currently focal clade (based on the posteriors of its children)
		const long node = clade-Ntips;
		for(long s=0; s<Nstates; ++s) posteriors[node*Nstates+s] = 1.0; // initialize posteriors for this node, repopulate below based on children
		for(long e=clade2first_inout_edge[clade], edge, child; e<=clade2last_inout_edge[clade]; ++e){
			edge  = inout_edges[e];
			if(current_tree_edge[2*edge+0]!=clade) continue; // this edge is not outgoing from this clade (in the rerooted tree)
			child = current_tree_edge[2*edge+1];
			if(use_precalculated_expQ){
				expQ_pointer = &expQ_per_edge[edge*Nstates*Nstates];
			}else{
				// calculate exponential of transition matrix along edge
				transition_exponentiator.get_exponential((edge_length.size()==0 ? 1.0 : edge_length[edge])/max_edge_length, expQ);
				expQ_pointer = &expQ[0];
			}
			// use exponentiated transition matrix to propagate probabilities from children to parent
			// probabilities[parent] = product_{child in children} exp(edge_length * Q^T) * probabilities[child]
			if(child<Ntips) multiply_matrix_with_vector(Nstates, Nstates, expQ_pointer, &prior_probabilities_per_tip[child*Nstates], Y);
			else multiply_matrix_with_vector(Nstates, Nstates, expQ_pointer, &posteriors[(child-Ntips)*Nstates], Y);
			for(long s=0; s<Nstates; ++s) posteriors[node*Nstates+s] *= max(0.0,Y[s]); // factor Y into the posterior of this node. Avoid negative values from rounding errors
		}
			
		// rescale (normalize) clade's probability distribution
		double S = 0;
		for(long s=0; s<Nstates; ++s) S += posteriors[node*Nstates+s];
		for(long s=0; s<Nstates; ++s) posteriors[node*Nstates+s] /= S;
		
		// move on to parent
		if(clade!=new_root) clade = current_tree_edge[current_incoming_edge_per_clade[clade]*2+0];
		else break;
	}
}




// calculate the loglikelihood of a fixed-rates Markov model for discrete character evolution on a phylogenetic tree, provided a fixed transition matrix
// Optionally, the marginal ancestral likelihoods can be computed for all nodes, using the rerooting algorithm by [Yang et al. (1995). Genetics 141:1641-1650]
// Calculating marginal ancestral likelihoods substantially increases computation time, so don't request this if you only care about the loglikelihood (e.g. for fitting purposes)
// Optionally, the eigendecomposition of the transition matrix (eigenvalues & eigenvectors) can be provided to potentially speed up exponentiations. If provided, this is used blindly for all exponentiations.
// If an eigendecomposition for the transition_matrix is not provided, then a Taylor series (matrix polynomials) approximation is used instead. This includes preconditioning steps and seems to scratch quite well, and is similarly fast as using an eigendecomposition.
// Requirements:
//   Tree can include multi- and mono-furcations.
//   Tree must be rooted. Root will be determined automatically as the node with no parent.
// [[Rcpp::export]]
Rcpp::List ASR_with_fixed_rates_Markov_model_CPP(	const long					Ntips,
													const long 					Nnodes,
													const long					Nedges,
													const long					Nstates,
													const IntegerVector 		&tree_edge,						// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
													const NumericVector 		&edge_length, 					// (INPUT) 1D array of size Nedges, or an empty std::vector (all edges have length 1)
													const std::vector<double> 	&transition_matrix,				// (INPUT) 2D array of size Nstates x Nstates, in row-major format. Transition-rate matrix Q in row-major format, i.e. Q[r*Nstates + c] is (r,c)-th-element of Q and equal to the transition rate r-->c.
													const ComplexVector			&eigenvalues,					// (INPUT) Optional 1D vector of size Nstates, listing the eigenvalues of the transition_matrix (corresponding to some diagonalization). Can also be an empty vector if eigendecomposition not available.
													const ComplexVector			&EVmatrix,						// (INPUT) Optional 2D array of size Nstates x Nstates, in row-major format, whose columns are the eigenvectors of the transition_matrix. Can also be an empty vector if eigendecomposition not available.
													const ComplexVector			&inverse_EVmatrix,				// (INPUT) Optional 2D array of size Nstates x Nstates, in row-major format, the inverse of EVmatrix. Can also be an empty vector if eigendecomposition not available.
													const NumericVector			&prior_probabilities_per_tip, 	// (INPUT) 2D array of size Ntips x Nstates, in row-major format, listing prior probability distributions for each tip
													const NumericVector			&prior_probabilities_for_root,	// (INPUT) 1D array of size Nstates, listing prior probability distribution for root. Which prior you use for the root is often a matter of choice.
													bool						include_ancestral_likelihoods,	// (INPUT) whether to include ancestral likelihoods in return variables
													bool						reroot,							// (INPUT) if true, then the marginal ancestral states estimates (conditional scaled likelihoods as if each node was a root) of all nodes are also returned as an array of size Nnodes x Nstates. Only use if needed, since it's computationally expensive.
													const double				exponentiation_accuracy,		// (INPUT) maximum allowed error when exponentiating the transition matrix via polynomials, in terms of the Hilbert-Schmidt L2 norm. Only relevant if exponentiation is done using the polynomials.
													const long					max_polynomials,				// (INPUT) maximum possible number of polynomials to use for exponentiating the transition_matrix via polynomials, regardless of the pursued accuracy epsilon. Used as safety vault, but may break the guaranteed accuracy. A value ~100 is usually enough.
													const bool					store_exponentials){			// (INPUT) if True, then exponentials are pre-calculated and stored for the calculation of ancestral_likelihoods. This may save time because each exponential is only calculated once, but will use up more memory since all exponentials are stored. Only relevant if reroot==TRUE, otherwise exponentials are never stored.
	const long Nclades 					= Ntips + Nnodes;
	const bool use_precalculated_expQ 	= (reroot && store_exponentials);
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

	// prepare data structure for exponentiations of transition matrix, if needed
	matrix_exponentiator transition_exponentiator;
	if(has_eigendecomposition){
		// prepare exponentiator using eigen-decomposition
		transition_exponentiator.initialize(Nstates, eigenvaluesCPP, EVmatrixCPP, inverse_EVmatrixCPP, max_edge_length);
	}else{
		// prepare exponentiator using matrix polynomials
		const long min_polynomials = min_polynomials_for_positive_exponential_of_irreducible_matrix(Nstates, transition_matrix);
		transition_exponentiator.initialize(Nstates, transition_matrix, max_edge_length, 1e-4, min_polynomials, max_polynomials, true);
	}
											
	// pre-calculate exponentials of transition_matrix along edges if needed
	// This is faster because it avoids repeated calculations of the same exponential, but needs more RAM if Nedges is large
	std::vector<double> expQ_per_edge;
	if(use_precalculated_expQ){
		std::vector<double> scratch_expQ;
		expQ_per_edge.resize(Nedges*Nstates*Nstates);
		for(long edge=0; edge<Nedges; ++edge){
			transition_exponentiator.get_exponential((edge_length.size()==0 ? 1.0 : edge_length[edge])/max_edge_length, scratch_expQ);
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
											transition_exponentiator,
											expQ_per_edge,
											traversal_queue,
											traversal_node2first_edge,
											traversal_node2last_edge,
											traversal_edges,
											posteriors,
											loglikelihood);
	if(!include_ancestral_likelihoods){
		return Rcpp::List::create(	Rcpp::Named("loglikelihood") = loglikelihood);
	}else if(!reroot){
		return Rcpp::List::create(	Rcpp::Named("loglikelihood") = loglikelihood,
									Rcpp::Named("ancestral_likelihoods") = Rcpp::wrap(posteriors));
	}

	// calculate marginal ancestral states (posterior probabilities) at each node, as if that node was the root [Yang et al. 1995]
	// note that the original edge mappings (e.g. traversal_node2first_edge[]) will no longer be consistent with current_tree_edge after rerooting
	// Notation: current_(..) refers to tree-access structures that are updated at each rerooting. They will be consistent with each other, but not necessarily with the original tree structure.
	std::vector<double> ancestral_likelihoods(Nnodes*Nstates);
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
																transition_exponentiator,
																expQ_per_edge,
																current_tree_edge,
																current_incoming_edge_per_clade,
																posteriors);
		// the posteriors of the root are equal to its ancestral_likelihoods, so extract those
		for(long s=0; s<Nstates; ++s) ancestral_likelihoods[(new_root-Ntips)*Nstates+s] = posteriors[(new_root-Ntips)*Nstates+s];
		current_root = new_root;
		// abort if the user has interrupted the calling R program
		Rcpp::checkUserInterrupt();
	}
	
	return Rcpp::List::create(	Rcpp::Named("loglikelihood") = loglikelihood,
								Rcpp::Named("ancestral_likelihoods") = Rcpp::wrap(ancestral_likelihoods));
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
NumericVector apply_fixed_rate_Markov_model_to_missing_clades_CPP(	const long					Ntips,
																	const long 					Nnodes,
																	const long					Nedges,
																	const long					Nstates,
																	const IntegerVector 		&tree_edge,				// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
																	const NumericVector 		&edge_length, 			// (INPUT) 1D array of size Nedges, or an empty std::vector (all branches have length 1)
																	const std::vector<double> 	&transition_matrix,		// (INPUT) 2D array of size Nstates x Nstates, in row-major format. Transition-rate matrix Q in row-major format, i.e. Q[r*Nstates + c] is (r,c)-th-element of Q and equal to the transition rate r-->c. Sum-per-row should be zero.
																	const double				exponentiation_accuracy,// (INPUT) maximum allowed error when exponentiating the transition matrix, in terms of the Hilbert-Schmidt L2 norm.
																	const long					max_polynomials,		// (INPUT) maximum possible number of polynomials to use for exponentiating the transition_matrix, regardless of the pursued accuracy epsilon. Used as safety vault, but may break the guaranteed accuracy. A value ~100 is usually enough.
																	LogicalVector				likelihoods_known,		// (INPUT) 1D array of size Nclades, indicating whether the likelihoods for a particular clade are known (1) or unknown/to be determined (0).
																	NumericVector 				likelihoods,			// (INPUT) 2D matrix of size Nclades x Nstates, in row-major format. Likelihoods of each state in each clade (tip & node) of the tree.
																	const bool					unknown_likelihoods_as_priors){	// (INPUT) use unknown likelihoods (i.e. likelihoods[r,:] when likelihoods_known[r]==false) as priors. If false, unknown likelihoods are ignored, and effectively a flat prior is used

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
	matrix_exponentiator transition_exponentiator;
	const long min_polynomials = min_polynomials_for_positive_exponential_of_irreducible_matrix(Nstates, transition_matrix);
	transition_exponentiator.initialize(Nstates, transition_matrix, max_edge_length, exponentiation_accuracy, min_polynomials, max_polynomials, true);
	
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
			transition_exponentiator.get_exponential((edge_length.size()==0 ? 1.0 : edge_length[edge])/max_edge_length, expQ);
			// propagate clade's likelihoods to child by multiplying with exponentiated transition matrix
			multiply_matrix_with_vector(Nstates, Nstates, &expQ[0], &likelihoods[clade*Nstates], Y);
			if(unknown_likelihoods_as_priors){
				for(long s=0; s<Nstates; ++s) likelihoods[child*Nstates+s] = Y[s];
			}else{
				for(long s=0; s<Nstates; ++s) likelihoods[child*Nstates+s] *= Y[s];
			}
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
		// this is equivalent to the X_k from Felsenstein's phylogenetic independent contrasts
		// however we used Maddison's quadratic parameters because it's convenient and has been generalized to multifurcations [Maddison 1991]
		for(node=0; node<Nnodes; ++node){
			ancestral_states[node] = -quadratic_parameters_per_node[node*3+1]/(2*quadratic_parameters_per_node[node*3+0]); // [Maddison 1991, page 312]
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




// [[Rcpp::export]]
Rcpp::List get_mean_state_per_node_CPP(	const long					Ntips,
										const long 					Nnodes,
										const long					Nedges,
										const IntegerVector 		&tree_edge,					// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
										const NumericVector			&edge_length,				// (INPUT) 1D array of size Nedges, or an empty std::vector (all edges have length 1)
										const std::vector<double>	&tip_states){
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
	
	// calculate average & standard deviation per node, traversing tips-->root (excluding the root)
	std::vector<double> means(Nnodes,0), stds(Nnodes,0), counts(Nnodes,0);
	long clade, parent, pnode, cnode;
	for(long q=traversal_queue.size()-1; q>=1; --q){
		clade	= traversal_queue[q];
		parent 	= clade2parent[clade];
		cnode	= clade - Ntips;
		pnode	= parent - Ntips;
		if(clade<Ntips){
			means[pnode] 	+= tip_states[clade];
			stds[pnode]		+= SQ(tip_states[clade]);
			counts[pnode]	+= 1;
		}else{
			means[pnode] 	+= means[cnode];
			stds[pnode] 	+= stds[cnode];
			counts[pnode] 	+= counts[cnode];
		}
	}
	for(long node=0; node<Nnodes; ++node){
		means[node] /= counts[node];
		stds[node]  = sqrt(stds[node]/counts[node] - SQ(means[node]));
	}
	
	return	Rcpp::List::create(	Rcpp::Named("means")  	= Rcpp::wrap(means),
								Rcpp::Named("stds") 	= Rcpp::wrap(stds),
								Rcpp::Named("counts")	= Rcpp::wrap(counts));										
}



// ASR via Phylogenetic Independent Contrasts (PIC) for a scalar continuous trait on a tree [Felsenstein 1985, page 10]
// Confidence intervals can be calculated as described by [Garland et al. (1999), page 377]
// Literature:
//    Felsenstein (1985). Phylogenies and the Comparative Method. The American Naturalist. 125:1-15.
//    Garland et al. (1999). An introduction to phylogenetically based statistical methods, with a new method for confidence intervals on ancestral values. American Zoologist. 39:374-388.
// Requirements:
//   Tree can include multi- and mono-furcations (multifurcations are automatically expanded into descending bifurcations).
//   Tree must be rooted. Root will be determined automatically as the node with no parent.
// [[Rcpp::export]]
Rcpp::List ASR_via_independent_contrasts_CPP(	const long					Ntips,
												const long 					Nnodes,
												const long					Nedges,
												const IntegerVector 		&tree_edge,					// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
												const NumericVector			&edge_length,				// (INPUT) 1D array of size Nedges, or an empty std::vector (all edges have length 1)
												const std::vector<double>	&tip_states,				// (INPUT) 1D array of size Ntips, listing the numeric state at each tip
												const bool					include_standard_errors){	// (INPUT) if true, then standard errors of the local estimates are also calculated according to [Garland et al (1999)]
	long child, node, clade;
	double length, total_weight;

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
	if(Nmultifurcations>0){
		// Tree has multifurcations, so expand them first to bifurcations
		// Note that the number of monofurcations will remain unchanged, but the number of bifurcations/nodes/edges will increase
		// dummy nodes will always descend from the original nodes, and will be indexed Nnodes, Nnodes+1, ...
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
	
	// get incoming edge for each clade
	std::vector<long> incoming_edge_per_clade;
	get_incoming_edge_per_clade(Ntips, Nlocal_nodes, Nlocal_edges, local_tree_edge, incoming_edge_per_clade);
	
	// get root
	const long local_root = get_root_from_incoming_edge_per_clade(Ntips, local_tree_edge, incoming_edge_per_clade);

	// prepare tree traversal route (root-->tips) and edge mappings
	std::vector<long> traversal_queue, node2first_edge, node2last_edge, edge_mapping;
	get_tree_traversal_root_to_tips(Ntips,
									Nlocal_nodes,
									Nlocal_edges,
									local_root,
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
	const double edge_length_epsilon = RELATIVE_EPSILON * get_array_nonzero_min(local_edge_length); // substitute to be used for zero edge lengths
	std::vector<double> incoming_length_per_clade(Nlocal_clades);
	if(local_edge_length.size()>0){
		for(long clade=0; clade<Nlocal_clades; ++clade){
			if(clade!=local_root){
				length = local_edge_length[incoming_edge_per_clade[clade]];
				incoming_length_per_clade[clade] = (length==0 ? edge_length_epsilon : length);
			}
		}
	}else{
		incoming_length_per_clade.assign(incoming_length_per_clade.size(),1);
		incoming_length_per_clade[local_root] = 0;
	}
									
	// calculate Felsenstein's X_k and PICs in a postorder traversal (tips-->root)
	std::vector<double> node_states(Nlocal_nodes,0);
	for(long q=traversal_queue.size()-1; q>=0; --q){
		clade	= traversal_queue[q];
		node	= clade - Ntips;
		// calculate Felsenstein's X_k (node_states) and nu_k (incoming_length_per_clade)
		total_weight = 0;
		for(long e=node2first_edge[node]; e<=node2last_edge[node]; ++e){
			child 	= local_tree_edge[2*edge_mapping[e]+1];
			length 	= incoming_length_per_clade[child];
			node_states[node] 	+= (1.0/length) * (child<Ntips ? tip_states[child] : node_states[child-Ntips]);
			total_weight 		+= (1.0/length);
		}
		node_states[node] /= total_weight;
		incoming_length_per_clade[clade] += 1.0/total_weight;
	}
	
	// calculate standard errors as described by [Garland et al. (1999). Page 377]
	// traverse tips-->root in order to calculate cumulative sum of squared standardized PICs descending from each node
	std::vector<double> node_standard_errors;
	std::vector<double> sum_squared_standardized_PICs; // cumulative sum of squared standardized PICs descending from each node
	std::vector<long> NPICs_per_node;
	std::vector<double> node_CI95s;
	if(include_standard_errors){
		node_standard_errors.resize(Nnodes);
		sum_squared_standardized_PICs.assign(Nnodes,0);
		NPICs_per_node.assign(Nnodes,0);
		node_CI95s.resize(Nnodes);
		double X1, X2, distance;
		long child1, child2;
		for(long q=traversal_queue.size()-1; q>=0; --q){
			clade	= traversal_queue[q];
			node	= clade - Ntips;
			if(node2last_edge[node]==node2first_edge[node]){
				// Treat monofurcating nodes in a special way. 
				// Note that monofurcating nodes acquire the same state as their child,
				// while their modified incoming_length is the same as their child plus the length of their incoming edge
				child1 = local_tree_edge[2*edge_mapping[node2first_edge[node]]+1];
				if(child1>=Ntips){
					sum_squared_standardized_PICs[node] = sum_squared_standardized_PICs[child1-Ntips];
					NPICs_per_node[node] = NPICs_per_node[child1-Ntips];
				}
				node_standard_errors[node] = (NPICs_per_node[node]==0 ? NAN_D : sqrt((sum_squared_standardized_PICs[node]/NPICs_per_node[node]) * incoming_length_per_clade[child1]));
			}else{
				child1		= local_tree_edge[2*edge_mapping[node2first_edge[node]]+1];
				child2		= local_tree_edge[2*edge_mapping[node2first_edge[node]+1]+1];
				X1 			= (child1<Ntips ? tip_states[child1] : node_states[child1-Ntips]);
				X2 			= (child2<Ntips ? tip_states[child2] : node_states[child2-Ntips]);
				distance	= incoming_length_per_clade[child1] + incoming_length_per_clade[child2];
				sum_squared_standardized_PICs[node] = SQ(X2 - X1)/distance;
				NPICs_per_node[node] = 1;
				if(child1>=Ntips){
					sum_squared_standardized_PICs[node] += sum_squared_standardized_PICs[child1-Ntips];
					NPICs_per_node[node] += NPICs_per_node[child1-Ntips];
				}
				if(child2>=Ntips){
					sum_squared_standardized_PICs[node] += sum_squared_standardized_PICs[child2-Ntips];
					NPICs_per_node[node] += NPICs_per_node[child2-Ntips];
				}
				node_standard_errors[node] = (NPICs_per_node[node]==0 ? NAN_D : sqrt((sum_squared_standardized_PICs[node]/NPICs_per_node[node]) * (incoming_length_per_clade[child1]*incoming_length_per_clade[child2]/(incoming_length_per_clade[child1]+incoming_length_per_clade[child2]))));
			}
			node_CI95s[node] = quantile_Students_t(0.975, NPICs_per_node[node]) * node_standard_errors[node];
		}
	}
	
	// omit reconstructed states for dummy nodes, i.e. only keep original (multifurcating) nodes
	if(node_states.size()>Nnodes) node_states.resize(Nnodes);
	if(node_standard_errors.size()>Nnodes) node_standard_errors.resize(Nnodes);
	if(node_CI95s.size()>Nnodes) node_CI95s.resize(Nnodes);

	return	Rcpp::List::create(	Rcpp::Named("node_states")  		= Rcpp::wrap(node_states),
								Rcpp::Named("node_standard_errors") = Rcpp::wrap(node_standard_errors),
								Rcpp::Named("node_CI95s")			= Rcpp::wrap(node_CI95s));
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
#pragma mark MuSSE model fitting
#pragma mark


// MuSSE: Multiple State Speciation Extinction dynamics
// This is an extension of the BiSSE model by Maddison (2007)
// The model can account for tips with unknown state as well as for sampling biases (tips with known state are biased towards certain states)
// References:
//  [Maddison et al (2007). Estimating a binary character's effect on speciation and extinction. Systematic Biology. 56:701-710] 
//  [FitzJohn et al. (2009). Estimating trait-dependent speciation and extinction rates from incompletely resolved phylogenies. Systematic Biology. 58:595-611]
//  [FitzJohn (2012). Diversitree: comparative phylogenetic analyses of diversification in R. Methods in Ecology and Evolution. 3:1084-1092]

typedef std::vector<double> MuSSEstateD;
typedef std::vector<double> MuSSEstateE;


// class for integrating D-part of a MuSSE model (Multiple State Speciation Extinction) along a single edge
// this model should be integrated in reverse age (i.e. increasing age)
// this class defines the ODE for calculating D over increasing age, given some MuSSE model and some initial condition D(0)
// the class can also be used to calculate the dynamics of the mapping Phi(t):D(0)-->D(t), by treating the current_state as a matrix (with Dmap(0)=Id).
// the class can also be used to calculate the inverse mapping over increasing age, i.e. Theta(t):=Phi(t)^{-1}. Use inverse=true and matrix_form=true, to enable this.
class MuSSEmodelD{
	// scratch space
	mutable std::vector<double> mapping; // 2D matrix of size Nstates x Nstates, used to map a current state to its linear rate of change
	mutable std::vector<double> AY, YA; // 1D vector or 2D matrix, depending on whether matrix_form is true or false
	
public:
	std::vector<MuSSEstateD> trajectory;
	std::vector<double> ages; // integration normally proceeds in increasing age direction (i.e. backward time).
	std::vector<double> initial; // vector of size Nstates, or 2D matrix of size Nstates x Nstates (if matrix_form==true), specifying the initial D. Must be set by whoever creates the model instance

	// alternative storage space, for storing the rescaled trajectory instead of the actual trajectory
	// Hence, trajectory[t] = trajectory_shape[t] * exp(trajectory_scale[t])
	// Such a storage format may be needed to avoid overflow or underflow (i.e. if the trajectory explodes to infinity or implodes to zero)
	std::vector<MuSSEstateD> trajectory_shape;
	std::vector<double> trajectory_scale; // log-scale of trajectory mean
	
	bool matrix_form; // if true, the model is a model for 2D square matrices, not for vectors, i.e. dM/dt = A*M instead of dV/dt = A*V (internally, state variables are still stored as vectors, in row-major form)
	bool inverse; // if true, the model describes the evolution of the inverse mapping. That is, if the model normally describes the evolution of Phi(t), where dPhi(t)/dT = A(T)*Phi(t) and Phi(t=0)=Id, then the inverse model describes the evolution of Theta(t):=Phi^{-1}(t), i.e. dTheta(t)/ds = -Theta(t)*A(t) and Theta(0)=Id
		
	// Pre-computed extinction probabilities (E) over age, as a functor
	// Must be set by whoever creates the model instance
	LinearInterpolationFunctor<MuSSEstateE> E;
	
	
	void clear(){
		trajectory.clear();
		ages.clear();
		trajectory_shape.clear();
		trajectory_scale.clear();
	}
	
	// model parameters
	// should be manually set by whoever creates the MuSSEmodel instance
	std::vector<double> transition_rates; // 2D array of size Nstates x Nstates, in row-major format, listing Markov transition rates between states. transition_rates[r,c] is the transition rate r-->c. Non-diagonal entries must be positive, the sum of each row must be zero.
	std::vector<double> speciation_rates; // 1D array of size Nstates, listing extinction rates at each state
	std::vector<double> extinction_rates; // 1D array of size Nstates, listing extinction rates at each state
	long Nstates;
	
	// default constructor
	MuSSEmodelD(){
		matrix_form 	= false;
		inverse 		= false;
		Nstates 		= 0;
	}
	
	// use parameters from another model instance
	template<class MODEL_TYPE>
	void adopt_parameters(const MODEL_TYPE &model2){
		transition_rates 	= model2.transition_rates;
		speciation_rates 	= model2.speciation_rates;
		extinction_rates 	= model2.extinction_rates;
		Nstates 			= model2.Nstates;
	}
		
	// model is notified by the numerical solver that a time series of size count will need to be stored
	void reserveSpaceForTimeSeries(long Nages){ 
		trajectory.clear();
		trajectory.reserve(Nages); 
		ages.clear();
		ages.reserve(Nages); 
	}

	void reserveSpaceForScaledTimeSeries(long Nages){ 
		trajectory_shape.clear();
		trajectory_shape.reserve(Nages); 
		trajectory_scale.clear();
		trajectory_scale.reserve(Nages); 
		ages.clear();
		ages.reserve(Nages); 
	}
	
	// provide initial state to numerical solver
	bool getInitialState(double age, MuSSEstateD &state) const{ 
		state = initial;
		return true; 
	}
	
	// provide scaled initial state to numerical solver
	bool getScaledInitialState(double age, MuSSEstateD &shape, double &scale) const{ 
		const double initial_mean = vector_mean(initial);
		if(initial_mean<=0) return false;
		scale = log(initial_mean);
		shape = initial/initial_mean;
		return true; 
	}

	// record a new time series point, provided by the numerical solver
	void registerState(double age, const MuSSEstateD &state){
		trajectory.push_back(state); 
		ages.push_back(age); 

		// make sure entries are in [0,1]
		const long i = trajectory.size()-1;
		for(long s=0; s<trajectory[i].size(); ++s) trajectory[i][s] = max(0.0, min(1.0, trajectory[i][s]));
	}
	
	
	// record a new trajectory point, provided by the numerical solver
	// The state X to be recorded is provided in rescaled format, i.e. state = exp(scale) * shape
	// You can either record shape and scale separately, or combine them to obtain the actual state
	void registerScaledState(double age, const MuSSEstateD &shape, const double scale){
		trajectory_shape.push_back(shape);
		trajectory_scale.push_back(scale);
		ages.push_back(age); 
		
		// make sure entries are >0
		const long i = trajectory_shape.size()-1;
		for(long s=0; s<trajectory_shape[i].size(); ++s) trajectory_shape[i][s] = max(0.0, trajectory_shape[i][s]);
	}
	
	// provide matrix encoding linear rates of change of the current state X, i.e. return A(t), where:
	//   dX/dt = A(t)*X(t) (if inverse==false)
	// or:
	//   dX/dt = X(t)*A(t) (if inverse==true)
	// note that, in principle, A may also depend on the current state X, i.e. A=A(t,X(t))
	// The returned A must be in row-major format
	void getLinearDynamics(double age, std::vector<double> &A) const{
		const MuSSEstateE current_E = E(age);
		// The mapping A is essentially the transition_matrix, plus some additional terms on the diagonal
		A = transition_rates;
		for(long r=0; r<Nstates; ++r){
			A[r*Nstates+r] += - (speciation_rates[r]+extinction_rates[r]) + 2*speciation_rates[r]*current_E[r]; // add more terms to diagonal
		}
		if(inverse) A *= -1;
	}
	
	
	// return the rate of change of S & Y, where S is the "scale" of the trajectory and Y is the "shape" (scaled version) of the trajectory
	// I.e. instead of returning dX/dt, return dS/dt and dY/dt
	// The correspondence between X and (Y,S) should be handled transparently by the model class, not by the calling ODE solver
	// In this particular case, "scale" S is defined as S := log(mean(X)) and Y := X/S, so that X = exp(S)*Y
	RequestedDynamics getRateOfChangeAtScaledState(double age, const MuSSEstateD &currentY, const double currentS, MuSSEstateD &rateY, double &rateS, MuSSEstateD &jump_scaled_state, double &jump_scale) const{
		// first construct mapping matrix, then apply to currentY and currentS
		getLinearDynamics(age,mapping);
		// apply mapping to current_state
		// if dX/dt = A*X (inverse=false), then:
		//    dS/dt = mean(A*Y) and dY/dt = A*Y - Y*mean(A*Y)
		// if dX/dt = X^T*A (inverse=true), then:
		//    dS/dt = mean(Y^T*A) and dY/dt = Y^T*A - Y^T*mean(Y^T*A)
		const double corrective_rate = abs(get_matrix_trace(Nstates,mapping)) * (vector_mean(currentY) - 1.0);
		if(matrix_form){
			// current_state is a 2D matrix of size Nstates x Nstates, so the rate of change should also be a 2D matrix of size Nstates x Nstates
			if(inverse){
				multiply_matrices(Nstates,Nstates,Nstates,currentY,mapping,YA);
				rateS = vector_mean(YA) + corrective_rate;
				rateY = YA - currentY*rateS;
			}else{
				multiply_matrices(Nstates,Nstates,Nstates,mapping,currentY,AY);
				rateS = vector_mean(AY) + corrective_rate;
				rateY = AY - currentY*rateS;
			}
		}else{
			// current_state is a 1D vector of size Nstates
			if(inverse){
				multiply_vector_with_matrix(Nstates,Nstates,currentY,mapping,YA);
				rateS = vector_mean(YA) + corrective_rate;
				rateY = YA - currentY*rateS;
			}else{
				multiply_matrix_with_vector(Nstates,Nstates,mapping,currentY,AY);
				rateS = vector_mean(AY) + corrective_rate;
				rateY = AY - currentY*rateS;
			}
		}
		return RequestedDynamicsRateOfChange; 
	}


	// provide rates of change (dD/dt) to the numerical solver (where t=age)
	RequestedDynamics getRateOfChangeAtState(double age, const MuSSEstateD &current_state, MuSSEstateD &rate_of_change, MuSSEstateD &jump_state) const{
		// first construct mapping matrix, then apply to current state
		getLinearDynamics(age,mapping);
		// apply mapping to current_state
		// dX/dt = A*X (if inverse=false), or dX/dt = X^T*A (if inverse=true)
		if(matrix_form){
			// current_state is a 2D matrix of size Nstates x Nstates, so the rate of change should also be a 2D matrix of size Nstates x Nstates
			if(inverse){
				multiply_matrices(Nstates,Nstates,Nstates,current_state,mapping,rate_of_change);
			}else{
				multiply_matrices(Nstates,Nstates,Nstates,mapping,current_state,rate_of_change);
			}
		}else{
			// current_state is a 1D vector of size Nstates
			if(inverse){
				multiply_vector_with_matrix(Nstates,Nstates,current_state,mapping,rate_of_change);
			}else{
				multiply_matrix_with_vector(Nstates,Nstates,mapping,current_state,rate_of_change);
			}
		}
		return RequestedDynamicsRateOfChange; 
	}

	// returns true if for some reason the age step should be refined, e.g. if the domain boundary is crossed
	// this check may be requested by the numerical solver
	// note that current_state & candidate_state may be a vector of size Nstates (if matrix_form==false) or a 2D matrix of size Nstates x Nstates (if matrix_form==true).
	bool checkShouldRefineTimeStep(double age, const MuSSEstateD &current_state, double dt, const MuSSEstateD &candidate_state) const{ 
		if(inverse){
			// since the inverted mapping Theta(t) = Phi^{-1}(t) may validly map outside of [0,1]^Nstates, we cannot insist on candidate_state being within [0,1]^Nstates
			return false;
		}else{
			for(long i=0; i<candidate_state.size(); ++i){
				if(candidate_state[i]>1) return true;
				if(candidate_state[i]<0) return true;
			}
		}
		return false;
	}
	
	// similar to checkShouldRefineTimeStep above, but handling rescaled states, i.e. split into "shape" Y and "scale" S
	// If X is the actual simulated state (e.g. likelihoods D), then S:=log(mean(X)) and Y:=X/exp(S)
	bool checkShouldRefineTimeStep(double age, const MuSSEstateD &currentY, const double currentS, double dt, const MuSSEstateD &candidateY, const double candidateS) const{ 
		if(inverse){
			// since the inverted mapping Theta(t) = Phi^{-1}(t) may validly map outside of [0,1]^Nstates, we cannot insist on candidate_state being within [0,1]^Nstates
			for(long i=0; i<candidateY.size(); ++i){
				if((currentY[0]>0) && (candidateY[0]>0) && (abs(currentY[i]/currentY[1] - candidateY[i]/candidateY[0])>0.01)) return(true); // change in shape seems to drastic for one time step
			}
			return false;
		}else{
			for(long i=0; i<candidateY.size(); ++i){
				if((candidateY[i]>0) && (candidateS+log(candidateY[i])>0)) return true; // check if candidate_state>1
				if(candidateY[i]<0) return true; // checking if candidateY<0
				if((currentY[0]>0) && (candidateY[0]>0) && (abs(currentY[i]/currentY[1] - candidateY[i]/candidateY[0])>0.01)) return(true); // change in shape seems to drastic for one time step
			}
		}
		return false;
	}
	
	// check if candidate_state is outside of the domain boundaries.
	// In that case, tries to correct the candidate_state to be the "last" valid state on the linear trajectory
	// If this is not possible (e.g. because previous_state was already on the boundary), the problematic components are brute-force adjusted to be within the domain. In this case, the routine returns CrossedBoundaryYesButFixedBruteForce.
	CrossedBoundary checkCrossedDomainBoundaryAndFix(	double				previous_age,
														const MuSSEstateD	&previous_state,				// previous state (assumed to be valid!)
														double 				&dt,							// (INPUT/OUTPUT) will be adjusted (reduced) if candidate_state crossed the boundary. The modified value is guaranteed to be within (0,dt]
														MuSSEstateD			&candidate_state,				// (INPUT/OUTPUT) if candidate_state is outside of domain boundaries, then this may become the "last" state (on the linear trajectory from previous_state to candidate_state) within the domain (if possible).		
														const bool			intermediate) const{			// (INPUT) is the candidate point an intermediate point (i.e. as used in Runge-Kutta schemes), or a final point (i.e., the final prediction for the candidate time)
		if(inverse) return CrossedBoundaryNo; // since the inverted mapping Theta(t) = Phi^{-1}(t) may validly map outside of [0,1]^Nstates, we cannot insist on candidate_state being within [0,1]^Nstates
		double lambda = 1;
		const double min_lambda = 0.00001;
		for(long s=0; s<candidate_state.size(); ++s){
			if(candidate_state[s]<0){
				if(previous_state[s]<=0){
					// even refining the step would probably not help
					lambda = 0; break;
				}
				lambda = min(lambda,(previous_state[s]-0)/(previous_state[s]-candidate_state[s]));
			}
			if((!intermediate) && (candidate_state[s]>1)){
				if(previous_state[s]>=1){
					// even refining the step would probably not help
					lambda = 0; break;
				}
				lambda = min(lambda,(1-previous_state[s])/(candidate_state[s]-previous_state[s]));
			}
		}
		if((lambda<1) && (lambda>min_lambda)){
			candidate_state = previous_state*(1-lambda) + candidate_state*lambda;
			dt *= lambda;
			return CrossedBoundaryYesButFixedByReducingTimeStep;
		}else if((lambda<1) && (lambda<=min_lambda)){
			// at least one state variable has to be fixed brute-force, so do this for all
			for(long s=0; s<candidate_state.size(); ++s){
				candidate_state[s] = max(0.0, min(1.0, candidate_state[s]));
			}
			return CrossedBoundaryYesButFixedBruteForce;
		}else{
			return CrossedBoundaryNo;
		}
	}

	// similar to checkCrossedDomainBoundaryAndFix above, but handling rescaled states, i.e. (Y,S) instead of the actual simulated state
	// S is the "scale" of the trajectory and Y is the "shape" of the trajectory (a rescaled variant of the trajectory)
	// If X is the actual simulated state (e.g. likelihoods D), then S:=log(mean(X)) and Y:=X/exp(S)
	CrossedBoundary checkCrossedDomainBoundaryAndFix(	double				previous_age,
														const MuSSEstateD	&previousY,			// (INPUT) shape of the previous state (assumed to be valid!)
														const double		previousS,			// (INPUT) scale of the previous state
														double 				&dt,				// (INPUT/OUTPUT) will be adjusted (reduced) if candidate_state crossed the boundary. The modified value is guaranteed to be within (0,dt]
														MuSSEstateD			&candidateY,		// (INPUT/OUTPUT) if candidate_state is outside of domain boundaries, then this may become the "last" scaled state (on the linear trajectory from previousY to candidateY) within the domain (if possible).		
														double				&candidateS,		// (INPUT/OUTPUT) scale of the candidate state
														const bool			intermediate) const{// (INPUT) is the candidate point an intermediate point (i.e. as used in Runge-Kutta schemes), or a final point (i.e., the final prediction for the candidate time)
		if(inverse) return CrossedBoundaryNo; // since the inverted mapping Theta(t) = Phi^{-1}(t) may validly map outside of [0,1]^Nstates, we cannot insist on candidate_state being within [0,1]^Nstates
				
		double lambda = 1;
		const double min_lambda = 0.00001;
		for(long s=0; s<candidateY.size(); ++s){
			if(candidateY[s]<0){
				if(previousY[s]<=0){
					// even refining the step would probably not help
					lambda = 0; break;
				}
				lambda = min(lambda,(previousY[s]-0)/(previousY[s]-candidateY[s]));
			}
			if((!intermediate) && (candidateY[s]>0) && (candidateS+log(candidateY[s])>0)){ // check if candidate_state>1
				if((previousY[s]<0) || (previousS+log(previousY[s])>=0)){ // check if previous_state was outside of [0,1)
					// even refining the step would probably not help
					lambda = 0; break;
				}
				lambda = min(lambda,(0-(previousS+log(previousY[s])))/((candidateS+log(candidateY[s]))-(previousS+log(previousY[s])))); // calculate necessary step refinement in log-scale space
			}
			
		}
		if((lambda<1) && (lambda>min_lambda)){
			candidateY 	= previousY*(1-lambda) + candidateY*lambda;
			candidateS 	= previousS*(1-lambda) + candidateS*lambda;
			dt *= lambda;
			return CrossedBoundaryYesButFixedByReducingTimeStep;
		}else if((lambda<1) && (lambda<=min_lambda)){
			// at least one state variable has to be fixed brute-force, so do this for all
			for(long s=0; s<candidateY.size(); ++s){
				candidateY[s] = max(0.0, candidateY[s]);
				candidateY[s] = exp(min(-candidateS,log(candidateY[s]))); // force candidate_state[s] to be <= 1
			}
			return CrossedBoundaryYesButFixedBruteForce;
		}else{
			return CrossedBoundaryNo;
		}
	}

	
	bool stateIsNaN(const MuSSEstateD &state) const{
		return contains_nan(state);
	}
	
	bool scaledStateIsNaN(const MuSSEstateD &scaled_state, const double scale) const{
		return (std::isnan(scale) || contains_nan(scaled_state));
	}
	
	// estimate the maximum absolute rate of change or provide a reasonable upper bound for it
	// may be useful for choosing the age step for simulations
	double estimate_max_rate_of_change() const{
		double max_rate = 0;
		for(long s=0; s<Nstates; ++s){
			const double sumT = row_sum(transition_rates,Nstates,s)-transition_rates[s*Nstates + s]; // sum of non-diagonal transition rates from state s
			// get maximum absolute rate of change of D[s]
			// The sum of multiple functions f1+f2+.. is always contained within [min(f1)+min(f2)+..., max(f1)+max(f2)+...]
			//   so calculate minima and maxima of individual summands, then use these to bound the maximum absolute value
			const double minD1 = - (speciation_rates[s] + extinction_rates[s] + sumT);
			const double maxD1 = 0;
			const double minD2 = 0;
			const double maxD2 = 2*speciation_rates[s];
			const double minD3 = 0;
			const double maxD3 = sumT;
			const double maxD = max(abs(minD1+minD2+minD3),abs(maxD1+maxD2+maxD3));
			// combine maximum absolute rate of change of D[s]
			max_rate = max(max_rate, maxD);
		}
		return max_rate;
	}
};





// class for integrating E-part of a MuSSE model (Multiple State Speciation Extinction) along a single edge
// this model should be integrated in reverse age (i.e. increasing age)
class MuSSEmodelE{
protected:
	// model parameters
	std::vector<double> transition_rates; // 2D array of size Nstates x Nstates, in row-major format, listing Markov transition rates between states. transition_rates[r,c] is the transition rate r-->c. Non-diagonal entries must be positive, the sum of each row must be zero.
	std::vector<double> speciation_rates; // 1D array of size Nstates, listing extinction rates at each state
	std::vector<double> extinction_rates; // 1D array of size Nstates, listing extinction rates at each state
	long Nstates;

	std::vector<double> linear_dynamics; // store the linear part of the E-dynamics. This is a constant matrix, under BiSSE/MuSSE/HiSSE/SecSSE models
	friend class MuSSEmodelD;
public:
	std::vector<MuSSEstateE> trajectory;
	std::vector<double> ages;
	std::vector<double> initial; // vector of size Nstates, specifying the initial E (at the start of the edge). Must be manually set by whoever creates the model instance
	bool matrix_form;
	
	void clear(){
		trajectory.clear();
		ages.clear();
	}
	
	void setup(	const long		_Nstates,
				const dvector 	&_transition_rates, 
				const dvector 	&_speciation_rates,
				const dvector 	&_extinction_rates){
		Nstates 		 = _Nstates;
		transition_rates = _transition_rates;
		speciation_rates = _speciation_rates;	
		extinction_rates = _extinction_rates;
		
		// setup linear part of dynamics, as a square matrix
		// this is basically the transition matrix, plus some stuff on the diagonal
		linear_dynamics = transition_rates;
		for(long s=0; s<Nstates; ++s){
			linear_dynamics[s*Nstates+s] -= (speciation_rates[s]+extinction_rates[s]);
		}
	}
	
	// use parameters from another model instance
	template<class MODEL_TYPE>
	void adopt_parameters(const MODEL_TYPE &model2){
		setup(model2.Nstates, model2.transition_rates, model2.speciation_rates, model2.extinction_rates);
		matrix_form = model2.matrix_form;
	}
		
	// model is notified by the numerical solver that a time series of size count will need to be stored
	void reserveSpaceForTimeSeries(long Nages){ 
		trajectory.clear();
		trajectory.reserve(Nages); 
		ages.clear();
		ages.reserve(Nages); 
	}

	// model is notified by the numerical solver that a time series will need to be stored
	void reserveSpaceForTimeSeries(const double simulationTime, const double minRecordingTimeStep, const double recordingValueStep){ 
		const double max_rate = estimate_max_rate_of_change();
		const long Nrecordings = 2 + min(simulationTime/minRecordingTimeStep, simulationTime * max_rate/recordingValueStep); // estimate total number of recording that will be performed
		trajectory.clear();
		trajectory.reserve(Nrecordings); 
		ages.clear();
		ages.reserve(Nrecordings); 
	}
	
	// provide initial state to numerical solver
	bool getInitialState(double age, MuSSEstateE &state) const{ 
		state = initial;
		return true; 
	}
	
	bool linearDynamicsAreConstant() const{ return true; }

	// record a new time series point, provided by the numerical solver
	void registerState(double age, const MuSSEstateE &state){
		trajectory.push_back(state); 
		ages.push_back(age); 
	}

	// provide rates of change to the numerical solver
	RequestedDynamics getRateOfChangeAtState(double age, const MuSSEstateE &current_state, MuSSEstateE &rate_of_change, MuSSEstateE &jump_state) const{
		if(matrix_form){
			// current_state is a 2D matrix of size Nstates x Nstates, so the rate of change should also be a 2D matrix of size Nstates x Nstates
			// treat every column (c) as an independent trajectory
			rate_of_change.resize(Nstates*Nstates);
			for(long c=0; c<Nstates; ++c){
				for(long s=0, j; s<Nstates; ++s){
					j = s*Nstates + c;
					rate_of_change[j] = extinction_rates[s] - (speciation_rates[s]+extinction_rates[s])*current_state[j] + speciation_rates[s]*SQ(current_state[j]);
					for(long z=0; z<Nstates; ++z){
						rate_of_change[j] += transition_rates[s*Nstates+z] * current_state[z*Nstates+c];
					}
				}
			}
		}else{
			rate_of_change.resize(Nstates);
			for(long s=0; s<Nstates; ++s){
				rate_of_change[s] = extinction_rates[s] - (speciation_rates[s]+extinction_rates[s])*current_state[s] + speciation_rates[s]*SQ(current_state[s]);
				for(long j=0; j<Nstates; ++j){
					rate_of_change[s] += transition_rates[s*Nstates+j] * current_state[j];
				}
			}
		}
		return RequestedDynamicsRateOfChange; 
	}


	// provide rates of change to the numerical solver, split into a linear and non-linear part
	// the rate of change of X will be: linearity*X + nonlinearity
	RequestedDynamics getLinearAndNonlinearDynamics(double age, const MuSSEstateE &current_state, dvector &linearity, MuSSEstateE &nonlinearity) const{
		linearity = linear_dynamics;
		if(matrix_form){
			// current_state is a 2D matrix of size Nstates x Nstates, so the rate of change should also be a 2D matrix of size Nstates x Nstates
			// treat every column (c) as an independent trajectory
			nonlinearity.resize(Nstates*Nstates);
			for(long c=0; c<Nstates; ++c){
				for(long s=0; s<Nstates; ++s){
					nonlinearity[s*Nstates + c] = extinction_rates[s] + speciation_rates[s]*SQ(current_state[s*Nstates + c]);
				}
			}
		}else{
			nonlinearity.resize(Nstates);
			for(long s=0; s<Nstates; ++s){
				nonlinearity[s] = extinction_rates[s] + speciation_rates[s]*SQ(current_state[s]);
			}
		}
		return RequestedDynamicsRateOfChange; 
	}



	// returns true if for some reason the age step should be refined, e.g. if the domain boundary is crossed
	// this check may be requested by the numerical solver
	bool checkShouldRefineTimeStep(double age, const MuSSEstateE &current_state, double dt, const MuSSEstateE &candidate_state) const{ 
		for(long s=0; s<candidate_state.size(); ++s){
			if(candidate_state[s]>1) return true;
			if(candidate_state[s]<0) return true;
			if(abs(candidate_state[s]-current_state[s])>0.0005*0.5*(candidate_state[s]+current_state[s])) return true; // relative change of state seems too drastic, so recommend reducing the time step
		}
		return false;
	}
	
	
	// calculate a measure of relative difference between two E-states
	// this may be requested by the numerical ODE solver, to decide whether a new point should be recorded
	double getRelativeChange(const MuSSEstateE &stateA, const MuSSEstateE &stateB) const{ 
		double relativeChange = 0;
		for(long s=0; s<stateA.size(); ++s){
			if((stateA[s]!=0) || (stateB[s]!=0)){
				relativeChange = max(relativeChange, abs(stateA[s]-stateB[s])*2/(abs(stateA[s])+abs(stateB[s])));
			}
		}
		return relativeChange;
	}
	
	// check if candidate_state is outside of the domain boundaries.
	// In that case, tries to correct the candidate_state to be the "last" valid state on the linear trajectory
	// If this is not possible (e.g. because previous_state was already on the boundary), the problematic components are brute-force adjusted to be within the domain. In this case, the routine returns CrossedBoundaryYesButFixedBruteForce.
	CrossedBoundary checkCrossedDomainBoundaryAndFix(	double				previous_age,
														const MuSSEstateE	&previous_state,		// previous state (assumed to be valid!)
														double 				&dt,					// (INPUT/OUTPUT) will be adjusted (reduced) if candidate_state crossed the boundary. The modified value is guaranteed to be within (0,dt]
														MuSSEstateE			&candidate_state,		// (INPUT/OUTPUT) if candidate_state is outside of domain boundaries, then this may become the "last" state (on the linear trajectory from previous_state to candidate_state) within the domain (if possible).
														const bool			intermediate) const{	// (INPUT) is the candidate point an intermediate point (i.e. as used in Runge-Kutta schemes), or a final point (i.e., the final prediction for the candidate time)
		double lambda = 1;
		for(long s=0; s<candidate_state.size(); ++s){
			if(candidate_state[s]>1){
				if(previous_state[s]>=1){
					// even refining the step would probably not help
					lambda = 0; break;
				}
				lambda = min(lambda,(1-previous_state[s])/(candidate_state[s]-previous_state[s]));
			}
			if((!intermediate) && (candidate_state[s]<0)){
				if(previous_state[s]<=0){
					// even refining the step would probably not help
					lambda = 0; break;
				}
				lambda = min(lambda,(previous_state[s]-0)/(previous_state[s]-candidate_state[s]));
			}
		}
		if((lambda<1) && (lambda>0)){
			candidate_state = previous_state*(1-lambda) + candidate_state*lambda;
			dt *= lambda;
			return CrossedBoundaryYesButFixedByReducingTimeStep;
		}else if((lambda<1) && (lambda<=0)){
			// at least one state variable has to be fixed brute-force, so do this for all
			for(long s=0; s<candidate_state.size(); ++s){
				candidate_state[s] = max(0.0, min(1.0, candidate_state[s]));
			}
			return CrossedBoundaryYesButFixedBruteForce;
		}else{
			return CrossedBoundaryNo;
		}
	}
	
	bool stateIsNaN(const MuSSEstateE &state) const{
		return contains_nan(state);
	}
	
	// estimate the maximum absolute rate of change or provide a reasonable upper bound for it
	// may be useful for choosing the age step for simulations
	double estimate_max_rate_of_change() const{
		double max_rate = 0;
		for(long s=0; s<Nstates; ++s){
			const double sumT = row_sum(transition_rates,Nstates,s)-transition_rates[s*Nstates + s]; // sum of non-diagonal transition rates from state s
			// get maximum absolute rate of change of E[s]
			// The sum of multiple functions f1+f2+.. is always contained within [min(f1)+min(f2)+..., max(f1)+max(f2)+...]
			//   so calculate minima and maxima of individual summands, then use these to bound the maximum absolute value
			const double minE1 = 0;
			const double maxE1 = extinction_rates[s];
			const double minE2 = - (speciation_rates[s] + extinction_rates[s] + sumT);
			const double maxE2 = 0;
			const double minE3 = 0;
			const double maxE3 = speciation_rates[s];
			const double minE4 = 0;
			const double maxE4 = sumT;
			const double maxE = max(abs(minE1+minE2+minE3+minE4),abs(maxE1+maxE2+maxE3+maxE4));
			// combine maximum absolute rate of change of E[s]
			max_rate = max(max_rate, maxE);
		}
		return max_rate;
	}
};


// DEPRECATED: Original MuSSE algorithm, but somewhat slower than the newer version below
// calculate log-likelihood of MuSSE model (Multiple State Speciation Extinction) on a tree
// initial conditions for D and E should be provided by the caller
// Requirements:
//   Tree can include multi- and mono-furcations.
//   Tree must be rooted. Root will be determined automatically as the node with no parent.
//   Tree must be ultrametric (e.g. a timetree of extant species). In particular, all tips are assumed to have age 0.
// [[Rcpp::export]]
Rcpp::List get_MuSSE_loglikelihood_classic_CPP(	const long					Ntips,
												const long 					Nnodes,
												const long					Nedges,
												const long					Nstates,
												const IntegerVector 		&tree_edge,					// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
												const std::vector<double> 	&node_ages, 				// (INPUT) 1D array of size Nclades, specifying the age (time before present) of each tip & node. All tips are assumed to have age 0.
												const std::vector<double> 	&transition_rates,			// (INPUT) 2D array of size Nstates x Nstates, in row-major format. Transition-rate matrix Q in row-major format, i.e. Q[r*Nstates + c] is (r,c)-th-element of Q and equal to the transition rate r-->c.
												const std::vector<double> 	&speciation_rates,			// (INPUT) 1D array of size Nstate, specifying the speciation rate at each state.
												const std::vector<double> 	&extinction_rates,			// (INPUT) 1D array of size Nstate, specifying the extinction rate at each state.
												const std::vector<double>	&initial_D_per_tip, 		// (INPUT) 2D array of size Ntips x Nstates, in row-major format, listing initial conditions for D (clade likelihoods) for each tip
												const std::vector<double>	&initial_E_per_state, 		// (INPUT) 1D array of Nstates, listing initial conditions for E (extinction probabilities) conditional upon the state
												const std::vector<double>	&root_prior,				// (INPUT) 1D array of size Nstates, listing prior probability distribution for root. Used to combine the root's clade likelihoods (D) into an overall log-likelihood of the model										
												const double				runtime_out_seconds){		// (INPUT) max allowed MuSSE integration runtime in seconds, per edge. If <=0, this option is ignored.
	const long root 		= get_root_clade(Ntips, Nnodes, Nedges, tree_edge);
	const long root_node 	= root - Ntips;
	const double root_age 	= node_ages[root_node];

	// prepare tree traversal route (root-->tips)
	// Note: This seems to only have a minuscule contribution to the total runtime
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
									
	// prepare MuSSE model for E (extinction probabilities) and D (clade likelihoods)
	MuSSEmodelD modelD;
	MuSSEmodelE modelE;
	modelE.setup(Nstates, transition_rates, speciation_rates, extinction_rates);
	modelE.matrix_form 		= false;
	modelD.adopt_parameters(modelE);
	modelD.matrix_form 		= false;
	modelD.inverse 			= false;
	
	// integrate MuSSE model for E
	string warningMessage;
	modelE.initial = initial_E_per_state;
	const double dt = min(0.1*root_age,0.1/modelE.estimate_max_rate_of_change());
	const long Nrecords = min(100000.0,max(100.0,root_age*modelE.estimate_max_rate_of_change())); // number of points to record
	bool success = RungeKutta2<MuSSEstateE,MuSSEmodelE,ProgressReporter>
								(0,
								root_age,
								dt,
								modelE,
								Nrecords, 	// number of points to record
								2,	// maxTimeStepRefinements
								2, // refinement_factor
								ProgressReporter(true),
								runtime_out_seconds,
								warningMessage);
	if(!success) return Rcpp::List::create(Rcpp::Named("success") = false, Rcpp::Named("error") = "Could not integrate MuSSE model (extinction probabilities E): "+warningMessage ); // simulation failed
	
	// define interpolator functions for E and provide to D-model
	modelD.E = LinearInterpolationFunctor<MuSSEstateE>(modelE.ages,modelE.trajectory,false,modelE.trajectory[0],modelE.trajectory.back(),true,0);
		
	// traverse tips-->root, to calculate D at each node
	std::vector<MuSSEstateD> posteriors(Nnodes,speciation_rates); // 1D array storing MuSSE posterior clade likelihoods (D) for each node
	double loglikelihood = 0;
	for(long q=traversal_queue.size()-1, clade, node; q>=0; --q){
		clade  = traversal_queue[q];
		node   = clade-Ntips;
		// set MuSSE likelihoods of clade to the element-wise product of its children's MuSSE likelihoods
		for(long e=traversal_node2first_edge[node], edge, child; e<=traversal_node2last_edge[node]; ++e){
			edge  = traversal_edges[e];
			child = tree_edge[2*edge+1];
			const double child_age = (child<Ntips ? 0.0 : node_ages[child-Ntips]);

			// prepare MuSSE model ODE integration along this edge (specify initial condition)
			modelD.clear();
			if(child<Ntips){
				// use tip's prior as initial condition
				extract_row(initial_D_per_tip, Nstates, child, modelD.initial);
			}else{
				// use child's MuSSE posteriors as initial condition
				modelD.initial = posteriors[child-Ntips];
			}
			
			// solve MuSSE model ODE along edge and save final point as node's posterior
			if(node_ages[node]==child_age){
				// child has same age as parent node, so just take child posterior as-is (no need to integrate over zero time)
				posteriors[node] *= modelD.initial;
			}else{
				success = RungeKutta2<MuSSEstateD,MuSSEmodelD,ProgressReporter>
											(child_age,
											node_ages[node],
											min(0.5*(node_ages[node]-child_age),0.1/modelD.estimate_max_rate_of_change()), // dt
											modelD,
											1, 	// number of points to record
											2,	// maxTimeStepRefinements
											2, // refinement_factor
											ProgressReporter(true),
											runtime_out_seconds,
											warningMessage);
				if((!success) || modelD.ages.empty()) return Rcpp::List::create(Rcpp::Named("success") = false, Rcpp::Named("error") = "Could not integrate MuSSE model (clade likelihoods D) along edge: "+warningMessage ); // simulation failed
				posteriors[node] *= modelD.trajectory.back();
			}
		}
				
		// rescale (normalize) node's posterior D
		// this is necessary due to scaling issues for very large trees; the non-normalized posterior D tends to 0 for older nodes
		// note that since the MuSSE ODE is linear in D, rescaling just rescales the overall model likelihood (this is corrected for below)
		const double S = vector_sum(posteriors[node]);
		posteriors[node] /= S;
		
		// incorporate rescaling factor (used to normalize this node's posterior D) into tree's loglikelihood
		// if we weren't rescaling each node's posterior D, this would not be necessary
		loglikelihood += log(S);
	}
	
	// calculate model's log-likelihood from root's posterior
	loglikelihood += log(vector_sum(posteriors[root_node]*root_prior));
	return Rcpp::List::create(Rcpp::Named("success") = true, Rcpp::Named("loglikelihood") = loglikelihood);
}






// calculate log-likelihood of MuSSE model (Multiple State Speciation Extinction) on a tree
// initial conditions for D and E should be provided by the caller
// Requirements:
//   Tree can include multi- and mono-furcations.
//   Tree must be rooted. Root will be determined automatically as the node with no parent.
//   Tree must be ultrametric (e.g. a timetree of extant species). In particular, all tips are assumed to have age 0.
// [[Rcpp::export]]
Rcpp::List get_MuSSE_loglikelihood_CPP(	const long					Ntips,
										const long 					Nnodes,
										const long					Nedges,							// (INPUT) number of edges in the tree
										const long					Nstates,						// (INPUT) number of discrete states that the modeled trait can have
										const IntegerVector 		&tree_edge,						// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
										const std::vector<double> 	&node_ages, 					// (INPUT) 1D array of size Nclades, specifying the age (time before present) of each tip & node. All tips are assumed to have age 0.
										const std::vector<double> 	&transition_rates,				// (INPUT) 2D array of size Nstates x Nstates, in row-major format. Transition-rate matrix Q in row-major format, i.e. Q[r*Nstates + c] is (r,c)-th-element of Q and equal to the transition rate r-->c.
										const std::vector<double> 	&speciation_rates,				// (INPUT) 1D array of size Nstate, specifying the speciation rate at each state.
										const std::vector<double> 	&extinction_rates,				// (INPUT) 1D array of size Nstate, specifying the extinction rate at each state.
										const std::vector<double>	&initial_D_per_tip, 			// (INPUT) 2D array of size Ntips x Nstates, in row-major format, listing initial conditions for D (clade likelihoods) for each tip
										const std::vector<double>	&initial_E_per_state, 			// (INPUT) 1D array of Nstates, listing initial conditions for E (extinction probabilities) conditional upon the state
										std::vector<double>			&root_prior,					// (INPUT) 1D array of size Nstates, listing prior probability distribution for root. Used to combine the root's clade likelihoods (D) into an overall log-likelihood of the model. Can also be an empty vector, in which case the computed state-likelihoods at the root (D[s]) are used as prior probability distribution.
										const std::string			&root_conditioning,				// (INPUT) either "none", "madfitz" or "herr_als", specifying how to condition the root's state-likelihoods prior to averaging. "none" corresponds to the original BiSSE model by Maddison (2007), "madfitz" and "herr_als" are options introduced by the hisse R package.
										const bool					include_ancestral_likelihoods,	// (INPUT) whether to also return the state likelihoods (D) for each node. This may be used as "local" ancestral state reconstructions.
										const bool					include_warnings,				// (INPUT) whether to also return all warning messages that occurred
										const double				max_condition_number,			// (INPUT) unitless number, the maximum acceptable condition number for the Gmap (as estimated from the linearized dynamics), when choosing the integration interval size. A larger max_condition number leads to fewer age-splits, thus faster computation but also lower accuracy. Hence, this number controls the trade-off between speed and accuracy. Typical values are 1e4 (slower, more accurate) up to 1e8 (faster, less accurate).
										const double				relative_ODE_step,				// (INPUT) unitless number, default relative integration time step for the ODE solvers. Relative to the typical time scales of the dynamics, as estimated from the theoretically maximum possible rate of change of D or E. Typical values are 0.01 - 0.1.
										const double				E_value_step,					// (INPUT) unitless number, relative step for interpolating E over time. So a E_value_step of 0.001 means that E is recorded and interpolated between points between which E differs by roughy 0.001. Typical values are 0.01-0.0001. A smaller E_value_step increases interpolation accuracy, but also increases memory requirements and adds runtime (scales with the tree's age span, not Ntips).
										const double				D_temporal_resolution,			// (INPUT) unitless number, relative resolution for interpolating Gmap over time. This is relative to the "typical" time scales at which E and Gmap vary. So a resolution of 10 means for every typical time scale there will be 10 interpolation points. Typical values are 1-100. A greater resolution increases interpolation accuracy, but also increases memory requirements and adds runtime (scales with the tree's age span, not Ntips).
										const double				runtime_out_seconds){			// (INPUT) max allowed MuSSE integration runtime in seconds, per edge. If <=0, this option is ignored.
	const long root 			= get_root_clade(Ntips, Nnodes, Nedges, tree_edge);
	const long root_node 		= root - Ntips;
	const double root_age 		= node_ages[root_node];
	const double start_runtime 	= (runtime_out_seconds>0 ? get_thread_monotonic_walltime_seconds() : 0.0);
	std::vector<string> warnings;
	dvector scratchDV, dummyDV; // scratch space or dummy variables in the form of a vector of doubles
	
	// calculate rescaled birth rates (e.g. normalized relative to the mean birth rate)
	const double speciation_rate_log_scale 	= log(vector_mean(speciation_rates));
	const dvector scaled_speciation_rates	= speciation_rates/exp(speciation_rate_log_scale);
		
	// prepare tree traversal route (root-->tips)
	// Note: This seems to only have a minuscule contribution to the total runtime
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
									
	// prepare MuSSE model for E (extinction probabilities) and D (clade likelihoods)
	MuSSEmodelD modelD;
	MuSSEmodelE modelE;
	modelE.setup(Nstates, transition_rates, speciation_rates, extinction_rates);
	modelE.matrix_form = false;
	modelD.adopt_parameters(modelE);
	
	// integrate MuSSE model for E
	string warningMessage;
	modelE.initial = initial_E_per_state;
	bool success = RungeKutta2<MuSSEstateE,MuSSEmodelE,ProgressReporter>
								(0, // start_time
								root_age, // end_time
								max(0.000001*root_age,min(0.2*root_age,relative_ODE_step/modelE.estimate_max_rate_of_change())), // default time step
								modelE,
								1e-6/modelE.estimate_max_rate_of_change(), // minRecordingTimeStep
								E_value_step, 		// recordingRelValueStep
								5,			// maxTimeStepRefinements
								4, 			// refinement_factor
								ProgressReporter(true),
								(runtime_out_seconds>0 ? max(runtime_out_seconds*0.01, runtime_out_seconds+start_runtime-get_thread_monotonic_walltime_seconds()) : 0.0),
								warningMessage);	
	if(!success) return Rcpp::List::create(Rcpp::Named("success") = false, Rcpp::Named("error") = "Could not integrate MuSSE model (extinction probabilities E): "+warningMessage ); // simulation failed
	if(include_warnings && (warningMessage!="")) warnings.push_back("Numerical integration of MuSSE-ODE for E was problematic: "+warningMessage);
	// define interpolator functions for E and provide to D-model
	// interpolate on irregular recording grid. Since later value requests by the ODE solver will be for nearby (i.e. slightly varying) time points, value retrieval will be fast
	modelD.E = LinearInterpolationFunctor<MuSSEstateE>(modelE.ages,modelE.trajectory,false,modelE.trajectory[0],modelE.trajectory.back(),false,0);
	const MuSSEstateE root_E = modelE.trajectory.back();


	// integrate D model (mapping D(0)-->D(t)) from age 0 to root_age in incremental time intervals, and store as interpolators
	// We split the interval [0:root_age] into sub-intervals, because the integrated Gmaps converge to singular matrices over time, 
	//   so we need to "renew" them regularly (i.e. start again at identity matrix)
	// The extent of DeltaT is chosen adaptively depending on the "dissipation rate" of the linear dynamics, i.e. how fast the condition number of Gmap increases over time
	// The spectral range of the linear dynamics (i.e. the difference between the most positive and most negative eigenvalue in terms of their real part) dictates the dissipation rate of the D-process.
	// Specifically, the spectral range is the exponential rate at which the condition number of Gmap (k(Gmap)) will increase over time, since Gmap ~ exp(t * dynamics)
	// Hence, a larger spectral range necessitates a smaller DeltaT, i.e. splitting [0:root_age] into smaller sub-intervals
	modelD.inverse 		= false;
	modelD.matrix_form 	= true;
	// determine linear dynamics (matrix form) of D at various representative ages, and keep track of the worst spectral range encountered
	dvector dynamics, ages(2);
	ages[0] = 0; ages[1] = root_age;
	double max_spectral_range = 0;
	for(long a=0; a<ages.size(); ++a){
		modelD.getLinearDynamics(ages[a], dynamics); // get linear dynamics at this age
		// calculate spectral range (most positive minus most negative eigenvalue) of the dynamics at this age
		double spectral_range = get_spectral_range(Nstates, dynamics);
		if(std::isnan(spectral_range)){
			// failed to calculate spectral range for some reason
			if(include_warnings) warnings.push_back(stringprintf("Failed to get spectral range of D-dynamics at age %g; attempting to estimate upper bound using power-method",ages[a]));
			// try to get an upper bound, based on the dominant eigenvalue (i.e. with largest modulus)
			// If lambda is the dominant eigenvalue, then we know that the spectral range is not greater than 2*|lambda|
			double dominant_eigenvalue;
			bool eigenvalue_converged = get_dominant_eigenvalue(Nstates,dynamics,1000,1e-3,dummyDV,dominant_eigenvalue); // get dominant eigenvalue (by magnitude) of the linear dynamics
			if(include_warnings && (!eigenvalue_converged)) warnings.push_back(stringprintf("Power iteration of D-dynamics at age %g did not converge, so dominant eigenvalue could not be accurately estimated; using the provisional eigenvalue %g",ages[a],dominant_eigenvalue));
			spectral_range = 2*abs(dominant_eigenvalue);
		}
		max_spectral_range = max(max_spectral_range,spectral_range);
	}
	// choose DeltaT (integration interval for partial D-maps) and the corresponding number of age-intervals according to the spectral range of the dynamics
	// the condition number of Gmap after time DeltaT is ~ exp(SR*DeltaT), where SR is the typical spectral range of the dynamics
	// A greater number of intervals means that more matrix multiplications will be needed at each edge during the postorder traversal later on
	const long max_Nintervals 	= 100000; // hard limit on the number of age intervals allowed
	const double DeltaT 		= max(root_age/max_Nintervals, min(1.0000001*root_age,log(max_condition_number)/max_spectral_range));
	const long Nintervals 		= ceil(root_age/DeltaT);
	if(include_warnings && (Nintervals==max_Nintervals)) warnings.push_back("Number of age intervals for computing Gmap reached the upper limit; perhaps the rate values are too high compared to the tree's time scales?");
	// calculate the D-maps by integrating the ODEs for D across each age interval
	std::vector<LinearInterpolationFunctor<MuSSEstateD> > Gmap_shape_functors(Nintervals);
	std::vector<LinearInterpolationFunctor<double> > Gmap_scale_functors(Nintervals);
	const double Dmax_rate_of_change = modelD.estimate_max_rate_of_change(); // estimated maximum rate of change of D under the model's dynamics
	for(long n=0; n<Nintervals; ++n){
		get_identity_matrix(Nstates,modelD.initial);
		warningMessage = "";
		const double start_age = n*DeltaT;
		const double end_age = start_age+DeltaT;
		const double max_runtime_for_integration = (runtime_out_seconds>0 ? runtime_out_seconds+start_runtime-get_thread_monotonic_walltime_seconds() : 0.0);
		if((runtime_out_seconds>0) && (max_runtime_for_integration<=0)) return Rcpp::List::create(Rcpp::Named("success") = false, Rcpp::Named("error") = stringprintf("Aborted prematurely during pre-calculation of %d D-mappings, because we reached the maximum allowed processing time",Nintervals));
		const long NDpoints = max(10.0,min(1000000.0/Nintervals,D_temporal_resolution*DeltaT*Dmax_rate_of_change)); 	// number of points to record for this age interval		
		success = LinearScaledRungeKutta<MuSSEstateD,MuSSEmodelD,ProgressReporter>
									(Nstates,
									Nstates,
									start_age,
									end_age,
									max(0.00001*DeltaT,min(0.5*DeltaT,relative_ODE_step/Dmax_rate_of_change)), // default time step
									modelD,
									NDpoints, 	// number of points to record
									4,			// maxTimeStepRefinements
									2,			// refinement_factor
									4, 			// max_exp_order
									ProgressReporter(true),
									max_runtime_for_integration,
									warningMessage);
		if((!success) || (modelD.ages.back()<end_age - 0.001*DeltaT)) return Rcpp::List::create(Rcpp::Named("success") = false, Rcpp::Named("error") = stringprintf("Could not integrate MuSSE D-model along age interval %g --> %g: ",start_age,end_age)+warningMessage); // simulation failed
		if(include_warnings && (warningMessage!="")) warnings.push_back(stringprintf("Numerical integration of MuSSE-ODE for D between ages %g --> %g was problematic: ",start_age,end_age)+warningMessage);
		Gmap_shape_functors[n] = LinearInterpolationFunctor<MuSSEstateD>(modelD.ages,modelD.trajectory_shape,false,modelD.trajectory_shape[0],modelD.trajectory_shape.back(),true,DeltaT/(NDpoints-1));
		Gmap_scale_functors[n] = LinearInterpolationFunctor<double>(modelD.ages,modelD.trajectory_scale,false,modelD.trajectory_scale[0],modelD.trajectory_scale.back(),true,DeltaT/(NDpoints-1));
		Rcpp::checkUserInterrupt(); // abort if the user has interrupted the calling R program
	}
	

	// traverse tips-->root, to calculate D at each node
	std::vector<dvector> posteriors(Nnodes,dvector(Nstates,1.0)); // 1D array storing MuSSE posterior clade likelihoods (D) for each node
	double loglikelihood = 0;
	dvector child_D, clade_D;
	dvector Gmap_to_clade_shape, Gmap_to_child_shape;
	double Gmap_to_child_scale, scaling;
	dvector X(Nstates), Y(Nstates);
	dvector inversion_scratch(Nstates*Nstates);	
	long Gmap_to_child_rank, start_functor, end_functor, Nsplits;
	const long max_rank_warnings = 100; // maximum number of Gmap-rank associated warnings to include. Used to prevent returning huge lists of largely similar warnings.
	long Nrank_warnings = 0; // keep track of the number of Gmap-rank associated warnings
	for(long q=traversal_queue.size()-1, clade, node; q>=0; --q){
		clade = traversal_queue[q];
		node  = clade-Ntips;
		const double clade_age = node_ages[node];
		
		// account for splitting event, i.e. initialize the likelihoods at the node to lambda[]^(Nsplits-1)
		// Note that lambda^(Nsplits-1) is the leading-order term (in dt) for the probability of a Yule process (pure birth-process), with a per-capita birth-rate lambda, to have Nsplits splits after time interval dt (when starting with a single lineage)
		// Alternatively, a multiplication with lambda^(Nsplits-1) can be justified by first breaking a multifurcation up into nearby bifurcations, apply the original MuSSE formula, and then taking the limit where those bifurcations are infinitesimally close to each other
		// To prevent numerical under- or over-flow, we only multiply by the rescaled lambdas, and correct for the rescaling in the model's overall loglikelihood
		Nsplits = traversal_node2last_edge[node]-traversal_node2first_edge[node] + 1;
		for(long s=0; s<Nstates; ++s) posteriors[node][s] = pow(scaled_speciation_rates[s],Nsplits-1.0);
		loglikelihood += (Nsplits-1) * speciation_rate_log_scale;

		// set MuSSE likelihoods of clade to the element-wise product of its children's MuSSE likelihoods
		for(long e=traversal_node2first_edge[node], edge, child; e<=traversal_node2last_edge[node]; ++e){
			edge  = traversal_edges[e];
			child = tree_edge[2*edge+1];
			const double child_age = (child<Ntips ? 0.0 : node_ages[child-Ntips]);
			const double age_interval = clade_age - child_age;

			// specify initial condition (at child) for MuSSE model ODE integration along this edge
			if(child<Ntips){
				// child is a tip, so use child's prior as initial condition
				extract_row(initial_D_per_tip, Nstates, child, child_D);
			}else{
				// use child's MuSSE posteriors as initial condition
				child_D = posteriors[child-Ntips];
			}
			
			// map child_D --> clade_D
			// this is equivalent to solving the MuSSE model ODE along edge
			if(age_interval<=root_age*RELATIVE_EPSILON){
				// child has same age as parent node, so just take child posterior as-is (no need to integrate over zero time)
				posteriors[node] *= child_D;			
			}else{			
			
				start_functor 	= min(Nintervals-1,long(child_age/DeltaT)); // Gmap functor defined on the interval that includes child_age
				end_functor  	= min(Nintervals-1,long(clade_age/DeltaT)); // Gmap functor defined on the interval that includes clade_age. May be the same as start_functor.
				// map child_D --> X (defined at start_functor's start time)
				Gmap_shape_functors[start_functor].getValue(child_age,Gmap_to_child_shape);
				Gmap_to_child_scale = Gmap_scale_functors[start_functor](child_age);
				//LUsolveLinearSystem(&Gmap_to_child_shape[0],&inversion_scratch[0],Nstates,&child_D[0],1e-6*vector_abs_mean(child_D),10,&X[0]);
				QR_linear_least_squares(Nstates,Nstates,1,Gmap_to_child_shape,child_D,true,scratchDV,inversion_scratch,X,Gmap_to_child_rank);
				if(include_warnings && (Gmap_to_child_rank<Nstates)){
					if(Nrank_warnings<max_rank_warnings) warnings.push_back(stringprintf("G-map from age %g to %g is rank deficient (has estimated rank %d), and hence its inversionn is numerically unstable",start_functor*DeltaT,child_age,Gmap_to_child_rank));
					++Nrank_warnings;
				}
				// map X --> clade_D
				// multiply Gmaps for all age-intervals between child & clade
				clade_D = X;
				loglikelihood -= Gmap_to_child_scale; // correct LL for scaling of inverted Gmap_to_child
				for(long n=start_functor; n<end_functor; ++n){
					multiply_matrix_with_vector(Nstates,Nstates,Gmap_shape_functors[n].getLastReferenceValue(),clade_D,Y);
					make_vector_positive(Y);
					scaling = vector_mean(Y);
					clade_D  = Y/scaling;
					loglikelihood += log(scaling);
					loglikelihood += Gmap_scale_functors[n].getLastReferenceValue();							
				}
				Gmap_shape_functors[end_functor].getValue(clade_age,Gmap_to_clade_shape);
				multiply_matrix_with_vector(Nstates,Nstates,&Gmap_to_clade_shape[0],&clade_D[0],Y);
				clade_D  = Y;
				loglikelihood += Gmap_scale_functors[end_functor](clade_age);
				
			
				/* CODE V0. DENOVO INTEGRATION ALONG EDGE
				// map child-->clade
				modelD.clear();
				warningMessage = "";
				modelD.initial = child_D;
				const double dt = max(0.001*age_interval,min(0.5*age_interval,0.1/modelD.estimate_max_rate_of_change()));
				success = LinearScaledRungeKutta<MuSSEstateD,MuSSEmodelD,ProgressReporter>
								(Nstates,
								1,
								child_age,
								clade_age,
								dt, // default time step
								modelD,
								2, 	// number of points to record
								2,	// maxTimeStepRefinements
								2,	// refinement_factor
								2, 	// max_exp_order
								ProgressReporter(true),
								runtime_out_seconds,
								warningMessage);
				if((!success) || (modelD.ages.back()<(clade_age-age_interval*1e-4))) return Rcpp::List::create(	Rcpp::Named("success") = false, Rcpp::Named("error") = "Failed to integrate D-model along edge during postorder traversal (age "+makeString(child_age)+" --> "+makeString(clade_age)+"): "+warningMessage);
				clade_D = modelD.trajectory_shape.back();
				const double scaling = modelD.trajectory_scale.back();
				*/
				
				/* CODE V1. USING A SINGLE FLOW FROM AGE 0 --> CLADE-AGE
				// determine MuSSE mapping along this edge child-->clade (use pre-computed mappings D & D_inverse)
				// Gmap(child-->clade) = Gmap(0-->clade) * Gmap(0-->child)^{-1}
				Gmap_to_clade_shape = functor_Gmap_shape(clade_age); // scaled version ("shape") of Gmap_to_clade
				Gmap_to_clade_scale = functor_Gmap_scale(clade_age); // log-scale of Gmap_to_clade
				Gmap_to_child_shape = functor_Gmap_shape(child_age); // scaled version ("shape") of Gmap_to_child
				Gmap_to_child_scale = functor_Gmap_scale(child_age); // log-scale of Gmap_to_child

				// solve linear system: Find vector X such that: Gmap_to_child_shape * X = child_D
				LUsolveLinearSystem(&Gmap_to_child_shape[0],&inversion_scratch[0],Nstates,&child_D[0],1e-6*get_array_nonzero_min(child_D),100,&X[0]);
				
				// calculate: clade_D = Gmap_to_clade_shape * X
				multiply_matrix_with_vector(Nstates,Nstates,Gmap_to_clade_shape,X,clade_D);
				loglikelihood += Gmap_to_clade_scale - Gmap_to_child_scale; // since the scales of Gmap_to_child_shape & Gmap_to_clade_shape were not included in the mapping child-->clade (nor does it need to be, since posteriors[node] will be normalized anyway), they must be incorporated into the loglikelihood

				// at this point Gmap_shape is a scaled version of the actual Gmap:child-->clade that we would like
				// Specifically, Gmap_shape = Gmap/exp(Gmap_scale)
				*/

				replace_non_positives(clade_D, 1e-8*vector_abs_mean(clade_D)); // replace non-positive values with a relatively small value, to avoid NaNs in the loglikelihood
				posteriors[node] *= clade_D;
				
				// rescale (normalize) node's posterior D
				// this is necessary due to scaling issues for very large trees; the non-normalized posterior D tends to 0 for older nodes
				// note that since the MuSSE ODE is linear in D, rescaling just rescales the overall model likelihood (this is corrected for below)
				// normalization should be done at every sub-iteration (i.e. for every child), because for some very degenerate trees some nodes can have hundreds of children
				scaling = vector_sum(posteriors[node]);
				for(long s=0; s<Nstates; ++s) posteriors[node][s] /= scaling;
				loglikelihood += log(scaling); // correct model's loglikelihood for rescaling of this node's posterior

				// check validity of likelihoods so far
				if((scaling==0) || std::isnan(scaling) || std::isinf(scaling)) return Rcpp::List::create(	Rcpp::Named("success") = false, Rcpp::Named("error") = "Likelihood reached NaN or Inf during postorder traversal", Rcpp::Named("warnings") = Rcpp::wrap(warnings));
			}
		}
		
		// abort if the user has interrupted the calling R program, or if we ran out of time
		if(q%100==0){
			Rcpp::checkUserInterrupt();
			if((runtime_out_seconds>0) && (get_thread_monotonic_walltime_seconds()-start_runtime>=runtime_out_seconds)){
				return Rcpp::List::create(	Rcpp::Named("success") = false, Rcpp::Named("error") = "Aborted prematurely during postorder traversal, because we reached the maximum allowed processing time", Rcpp::Named("warnings") = Rcpp::wrap(warnings));
			}
		}
	}
	
	// include an additional summarizing warning if some rank-warnings were omitted
	if(include_warnings && (Nrank_warnings>max_rank_warnings)){
		warnings.push_back(stringprintf("An additional %d Gmap-rank-related warnings have been omitted",(Nrank_warnings-max_rank_warnings)));
	}

	// determine maximum-likelihood root stage, based on state-likelihoods
	const long ML_root_state = get_array_max(posteriors[root_node]);	
	
	// calculate root-prior, if needed
	if(root_prior.empty()){
		// use state-likelihoods at the root to define a probability distribution
		root_prior = posteriors[root_node]/vector_sum(posteriors[root_node]);
	}
		
	// condition root's state-likelihoods, if needed
	if(root_conditioning=="madfitz"){
		// this is the same as root.type="madfitz" in the hisse package, and condition.surv=TRUE in diversitree (function 'rootfunc.musse' in file 'model-musse.R')
		double scaling = 0;
		for(long s=0; s<Nstates; ++s) scaling += root_prior[s]*speciation_rates[s]*SQ(1 - root_E[s]);
		posteriors[root_node] /= scaling;
	}else if(root_conditioning=="herr_als"){
		double scaling;
		for(long s=0; s<Nstates; ++s){
			scaling = speciation_rates[s]*SQ(1 - root_E[s]);
			if(scaling>0) posteriors[root_node][s] /= scaling;
		}
	}
	
	// calculate model's log-likelihood from root's posterior
	loglikelihood += log(vector_sum(posteriors[root_node]*root_prior));
	
	// prepare return values
	Rcpp::List results = Rcpp::List::create(Rcpp::Named("success") 			= true, 
											Rcpp::Named("warnings") 		= Rcpp::wrap(warnings),
											Rcpp::Named("NErecordings") 	= modelE.ages.size(),
											Rcpp::Named("Nintervals") 		= Nintervals,
											Rcpp::Named("loglikelihood") 	= loglikelihood,
											Rcpp::Named("ML_root_state") 	= ML_root_state);
	if(include_ancestral_likelihoods){
		dvector ancestral_likelihoods;
		flatten_matrix(posteriors, ancestral_likelihoods); // flatten posteriors into row-major format
		results.push_back(Rcpp::wrap(ancestral_likelihoods), "ancestral_likelihoods");
	}
	return results;

}




#pragma mark -
#pragma mark Simulate models of trait evolution
#pragma mark



// Perform random simulations of a fixed-rates continuous-time Markov model of discrete character evolution
// Starting with a specified vector of root_probabilities, and moving from root to tips, each node is assigned a random state according to its parent's state and according to the markov transition matrix.
// Optionally, multiple independent simulations can be performed using the same model (e.g. as part of some Monte Carlo integration)
// [[Rcpp::export]]
Rcpp::List simulate_fixed_rates_Markov_model_CPP(	const long					Ntips,
													const long 					Nnodes,
													const long					Nedges,
													const long					Nstates,
													const IntegerVector 		&tree_edge,				// (INPUT) 2D array of size Nedges x 2, in row-major format, with elements in 0,..,(Nclades-1)				
													const NumericVector		 	&edge_length, 			// (INPUT) 1D array of size Nedges, or an empty std::vector (all branches have length 1)
													const std::vector<double> 	&transition_matrix,		// (INPUT) 2D array of size Nstates x Nstates, in row-major format. Transition-rate matrix Q in row-major format, i.e. Q[r*Nstates + c] is (r,c)-th-element of Q, which is the transition rate r-->c.
													const NumericVector			&root_probabilities,	// (INPUT) probability distribution of states at the root. sum(root_probabilities) must be 1.0.
													const bool					include_tips,			// (INPUT) include states for tips in the output
													const bool					include_nodes,			// (INPUT) include states for nodes in the output
													const long					Nsimulations){			// (INPUT) number of random simulations (draws) of the model on the tree. If 1, then a single simulation is performed, yielding a single random state for each node and/or tip.
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
	const long NPmin = min_polynomials_for_positive_exponential_of_irreducible_matrix(Nstates, transition_matrix);
	const matrix_exponentiator transition_exponentiator(Nstates, transition_matrix, max_edge_length, 1e-4, NPmin, 1000, true);
	
	// traverse root-->tips and draw random states, conditional upon their parent's state
	vector<double> expQ;
	vector<long> tip_states, node_states;
	if(include_tips) tip_states.resize(Nsimulations*Ntips);
	node_states.assign(Nsimulations*Nnodes,0); // always store node states, since needed for moving root-->tips. Assign default value so that valgrind memcheck does not complain about uninitialized values.
	long clade, edge, parent, parent_state, state=0;
	for(long q=0; q<traversal_queue.size(); ++q){
		clade = traversal_queue[q];
		if(clade!=root){
			edge 	= incoming_edge_per_clade[clade];
			parent 	= tree_edge[edge*2+0];
			// exponentiate transition matrix along incoming edge
			transition_exponentiator.get_exponential((edge_length.size()==0 ? 1.0 : edge_length[edge])/max_edge_length, expQ);
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
	double parent_state, state = 0;
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
	double parent_state, state = 0;
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
	double parent_state, state = 0;
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








#pragma mark -
#pragma mark Generating random trees via cladogenic models
#pragma mark 



// auxiliary function used to generate random phylogenetic trees based on some cladogenic model
void aux_finalize_generated_random_tree(const double				time,					// (INPUT)
										const bool					as_generations,			// (INPUT)
										const bool					coalescent,				// (INPUT)
										const bool					include_rates,			// (INPUT)
										const std::vector<double>	&clade2end_time,		// (INPUT)
										const std::vector<double>	&clade2birth_rate_pc,	// (INPUT) optional input, listing pc birth rate per clade. Only relevant if include_rates==true.
										const std::vector<double>	&clade2death_rate_pc,	// (INPUT) optional input, listing pc birth rate per clade. Only relevant if include_rates==true.
										const bool					tree_had_deaths,		// (INPUT) whether deaths (extinctions) were part of cladogenesis
										long						&Ntips,					// (INPUT/OUTPUT)
										long						&Nclades,				// (INPUT/OUTPUT)
										long 						&Nedges,				// (INPUT/OUTPUT)
										long						&root,					// (INPUT/OUTPUT)
										double						&root_time,				// (OUTPUT)
										std::vector<long>			&tree_edge,				// (INPUT/OUTPUT)
										std::vector<double>			&edge_length,			// (OUTPUT)
										std::vector<long>			&extant_tips,			// (INPUT/OUTPUT)
										std::vector<long>			&new2old_clade,			// (OUTPUT)
										std::vector<double> 		&birth_rates_pc,		// (OUTPUT) optional output of size Nclades, only returned if include_rates==true
										std::vector<double> 		&death_rates_pc){		// (OUTPUT) optional output of size Nclades, only returned if include_rates==true
	long Nnodes = Nclades - Ntips;
	edge_length.resize(Nedges);
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
	root_time = clade2end_time[root]; // root_time = time at which the root split
	
	// identify tips as the clades with no outgoing edge
	std::vector<bool> clade_is_tip(Nclades,true);
	for(long edge=0; edge<Nedges; ++edge){
		clade_is_tip[tree_edge[edge*2+0]] = false;
	}
	
	// re-number tip & node indices to conform with the phylo format, where tips are indexed first (0,..,Ntips-1) and nodes last (Ntips,..,Ntips+Nnodes-1).
	std::vector<long> old2new_clade(Nclades,-1);
	new2old_clade.resize(Nclades);
	long next_new_tip  = 0;
	long next_new_node = 0;
	for(long clade=0; clade<Nclades; ++clade){
		if(clade_is_tip[clade]) old2new_clade[clade] = (next_new_tip++);
		else old2new_clade[clade] = Ntips + (next_new_node++);
		new2old_clade[old2new_clade[clade]] = clade;
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
	if(coalescent && tree_had_deaths){
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
		
		// update new2old_clade & old2new_clade
		const std::vector<long> old2older_clade(new2old_clade);
		old2new_clade.assign(old2new_clade.size(),-1);
		new2old_clade.resize(Nclades);
		for(long new_clade=0; new_clade<pruning_new2old_clade.size(); ++new_clade){
			new2old_clade[new_clade] = old2older_clade[pruning_new2old_clade[new_clade]];
			old2new_clade[new2old_clade[new_clade]] = new_clade;
		}
	}
	
	// also archive pc birth & death rates (using new clade indices) if needed
	if(include_rates){
		birth_rates_pc.resize(Nclades);
		death_rates_pc.resize(Nclades);
		for(long c=0; c<Nclades; ++c){
			birth_rates_pc[c] = clade2birth_rate_pc[new2old_clade[c]];
			death_rates_pc[c] = clade2death_rate_pc[new2old_clade[c]];
		}
	}
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
	if(has_added_birth_rates) added_birth_rates_pc = LinearInterpolationFunctor<double>(additional_rates_times, additional_birth_rates_pc,additional_periodic,0.0,0.0,true,0);
	if(has_added_death_rates) added_death_rates_pc = LinearInterpolationFunctor<double>(additional_rates_times, additional_death_rates_pc,additional_periodic,0.0,0.0,true,0);
	
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
	long Nevents		= 0;
	double time 				= 0;
	double total_rate			= INFTY_D;
	double equilibrium_time 	= INFTY_D;
	double initial_growth_rate 	= NAN_D;
	while(((max_tips<=0) || ((coalescent ? extant_tips.size() : Ntips)<max_tips)) && ((max_time<=0) || (time+1/total_rate<max_time)) && ((max_time_since_equilibrium<0) || (time-equilibrium_time+1/total_rate<max_time_since_equilibrium))){
		// determine time of next speciation or extinction event
		// prevent deaths if only one tip is left
		const double NEtips 	= extant_tips.size();
		const double birth_rate = max(0.0, birth_rate_intercept + birth_rate_factor * pow(NEtips, birth_rate_exponent) + (has_added_birth_rates ? added_birth_rates_pc(time)*NEtips : 0.0));
		const double death_rate = (NEtips<=1 ? 0 : max(0.0, (death_rate_intercept + death_rate_factor * pow(NEtips, death_rate_exponent)) + (has_added_death_rates ? added_death_rates_pc(time)*NEtips : 0.0)));
		if(std::isnan(initial_growth_rate)) initial_growth_rate = birth_rate-death_rate;
		if(((birth_rate-death_rate)*initial_growth_rate<0) && (equilibrium_time>time)) equilibrium_time = time; // first crossing of equilibrium state, so keep record
		total_rate = birth_rate+death_rate;
		
		// determine next event
		const double dt = random_exponential_distribution(total_rate);
		time += dt;
		const bool birth = random_bernoulli(birth_rate/total_rate);
		++Nevents;
				
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
		if(Nevents%100==0) Rcpp::checkUserInterrupt();
	}
		
	if((coalescent ? extant_tips.size() : Ntips)<=1){
		// something went wrong (e.g. zero birth & death rates)
		return Rcpp::List::create(Rcpp::Named("success") = false, Rcpp::Named("error")="Generated tree is empty or has only one tip");
	}

	// add a small dt at the end to make all edges non-zero length
	time += random_exponential_distribution(total_rate);
	if(max_time>0) time = min(time, max_time); // prevent going past max_time
	if(max_time_since_equilibrium>=0) time = min(time, max_time_since_equilibrium+equilibrium_time); // prevent going past max_time_since_equilibrium
	
	// finalize tree (make valid phylo structure, make coalescent if needed)
	std::vector<double> edge_length, birth_rates_pc, death_rates_pc, vdummy1, vdummy2, vdummy3, vdummy4;
	std::vector<long> new2old_clade;
	double root_time;
	aux_finalize_generated_random_tree(	time,
										as_generations,
										coalescent,
										false, // don't return rates-per-clade
										clade2end_time,
										vdummy1,
										vdummy2,
										((death_rate_intercept!=0) || (death_rate_factor!=0) || has_added_death_rates), // whether deaths were included
										Ntips,
										Nclades,
										Nedges,
										root,
										root_time,
										tree_edge,
										edge_length,
										extant_tips,
										new2old_clade,
										vdummy3,
										vdummy4);
	long Nnodes = Nclades - Ntips;
	
	return Rcpp::List::create(	Rcpp::Named("success") 			= true,
								Rcpp::Named("tree_edge") 		= Rcpp::wrap(tree_edge),
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
												const double	root_birth_rate_pc,			// (INPUT) initial pc birth rate of root
												const double	root_death_rate_pc,			// (INPUT) initial pc death rate of root
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
	clade2birth_rate_pc.push_back(root_birth_rate_pc);
	clade2death_rate_pc.push_back(root_death_rate_pc);
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
			// pick per-capita birth & death rates for this child and update total birth & death rates
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
				// pick per-capita birth & death rates for this child and update total birth & death rates
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

	if(Ntips<=1){
		// something went wrong (e.g. zero birth & death rates)
		return Rcpp::List::create(Rcpp::Named("success") = false, Rcpp::Named("error")="Generated tree is empty or has only one tip");
	}
	
	// add a small dt at the end to make all edges non-zero length
	time += random_exponential_distribution(total_rate);
	if(max_time>0) time = min(time, max_time); // prevent going past max_time
	if(max_time_since_equilibrium>=0) time = min(time, max_time_since_equilibrium+equilibrium_time); // prevent going past max_time_since_equilibrium

	// finalize tree (make valid phylo structure, make coalescent if needed)
	std::vector<double> edge_length, birth_rates_pc, death_rates_pc;
	std::vector<long> new2old_clade;
	double root_time;
	aux_finalize_generated_random_tree(	time,
										as_generations,
										coalescent,
										include_rates,
										clade2end_time,
										clade2birth_rate_pc,
										clade2death_rate_pc,
										(max_death_rate_pc>0),
										Ntips,
										Nclades,
										Nedges,
										root,
										root_time,
										tree_edge,
										edge_length,
										extant_tips,
										new2old_clade,
										birth_rates_pc,
										death_rates_pc);
	long Nnodes = Nclades - Ntips;
	
	return Rcpp::List::create(	Rcpp::Named("success") 			= true,
								Rcpp::Named("tree_edge") 		= Rcpp::wrap(tree_edge),
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


// Generate a random phylogenetic tree under a speciation/extinction model, where species are born or go extinct as a Poissonian process
// New species are added by splitting one of the currently extant tips (chosen randomly) into Nsplits new tips
// This function is similar to generate_random_tree_CPP(..) above, with one important difference:
//    Per-capita speciation and extinction rates are modelled as discrete-state continuous-time Markov chains evolving along the tree edges, with transitions occuring along edges (but not at branching points)
//	  Hence per-capita speciation/extinction rates are scalar traits specific to each node & tip.
// [[Rcpp::export]]
Rcpp::List generate_random_tree_Mk_rates_CPP(	const long 	 				max_tips,					// (INPUT) max number of tips (extant tips, if coalescent==true). If <=0, this constraint is ignored.
												const double				max_time,					// (INPUT) max simulation time. If <=0, this constraint is ignored.
												const double				max_time_since_equilibrium,	// (INPUT) max simulation time, counting from the first time point where death_rate-birth_rate changed sign. This may be used as an alternative to (or in conjunction with) max_time to ensure the tree has reached speciation/extinction equilibrium. If <0, this constraint is ignored.
												const long					Nstates,					// (INPUT) number of possible states
												const std::vector<double>	&state_birth_rates,			// (INPUT) list of per-capita birth rates for each state
												const std::vector<double>	&state_death_rates,			// (INPUT) list of per-capita death rates for each state
												const long					root_state,					// (INPUT) root state (0,..,Nstates-1)
												const std::vector<double>	&transition_matrix,			// (INPUT) 2D array of size Nstates x Nstates, in row-major format. Transition-rate matrix Q for the possible states. In row-major format, i.e. Q[r*Nstates + c] is (r,c)-th-element of Q, which is the transition rate r-->c.
												const bool					coalescent,					// (INPUT) whether to return only the coalescent tree (i.e. including only extant tips)
												const long					Nsplits,					// (INPUT) number of children to create at each diversification event. Must be at least 2. For a bifurcating tree this should be set to 2. If >2, then the tree will be multifurcating.
												const bool					as_generations,				// (INPUT) if false, then edge lengths correspond to time. If true, then edge lengths correspond to generations (hence if coalescent==false, all edges will have unit length).
												const bool					all_transitions,			// (INPUT) if true, then all transitions between states are simulated along edges over time. This is the exact version of the model. If false, an approximation is used whereby transitions only occur at branching points (i.e. during speciation events).
												const bool					no_full_extinction,	// (INPUT) if true, then extinction of the entire tree is prevented. This is done by temporarily disabling extinctions when the number of extant tips is 1.
												const bool					include_birth_times,		// (INPUT) if true, then the times of speciations (in order of occurrence) will also be returned
												const bool					include_death_times,		// (INPUT) if true, then the times of extinctions (in order of occurrence) will also be returned
												const bool					include_rates){				// (INPUT) if true, then the per-capita birth & death rates for each clade will also be returned
	const double max_death_rate_pc = get_array_max(state_death_rates);
	const long expected_Nclades = (max_tips<0 ? 2l : max_tips);
	long next_Nsplits = max(2l, Nsplits);
	std::vector<long> tree_edge;
	std::vector<lvector> extant_tips(Nstates,std::vector<long>());
	std::vector<long> clade2parent;
	std::vector<double> clade2end_time;
	std::vector<long> clade2state; 		// state for each node/tip (i.e. determining the waiting times until speciation & extinction)
	std::vector<double> birth_times, death_times;

	// the following lookup tables are updated continuously as the tree grows:
	//   tree_edge: edge structure of the tree, with entries being clade indices
	//   extant_tips: keeps track of which clades are extant tips, grouped by state. For any state s and extant-tip i (in state s), the value of extant_tips[s][i] will be an index for clade2state[], clade2parent[] and clade2state[]
	//   clade2parent: keeps track of parent of each clade ever created (including nodes and dead tips)
	//   clade2end_time: keeps track of time at which each clade ever created split or went extinct (negative if clade is an extant tip)
	//   clade2state: keeps track of the state of each clade ever created
		
	// reserve memory based on rough expectations
	tree_edge.reserve(expected_Nclades*2);
	clade2parent.reserve(expected_Nclades);
	clade2end_time.reserve(expected_Nclades);
	clade2state.reserve(expected_Nclades);
	for(long state=0; state<Nstates; ++state) extant_tips[state].reserve(ceil(expected_Nclades/2/Nstates));
		
	// prepare exponentiation of birth & death rate matrix
	const double exponentiation_rescaling = max(1.0,max(max_time_since_equilibrium,max_time));
	const long NPmin = min_polynomials_for_positive_exponential_of_irreducible_matrix(Nstates, transition_matrix);
	const matrix_exponentiator transition_exponentiator(Nstates, transition_matrix, exponentiation_rescaling, 1e-4, NPmin, 1000, true);
	vector<double> scratch_exp;
	
	// pre-calculate auxiliary info about transition matrix
	std::vector<double> state_total_transition_rates(Nstates,0);
	for(long state=0; state<Nstates; ++state){
		state_total_transition_rates[state] = sum_of_row(Nstates, Nstates, transition_matrix, state) - transition_matrix[state*Nstates+state];
	}
		
	// create the first tip (which is also the root)
	long Nclades 	= 0;	// current number of clades
	long root 		= 0;
	extant_tips[root_state].push_back(Nclades++);
	clade2parent.push_back(-1); // root has no parent
	clade2end_time.push_back(-1);
	clade2state.push_back(root_state);
	long Ntips 		 = 1; // current number of extant + extinct tips
	long NextantTips = 1; // current number of extant tips (i.e. that can split or die). This must always be equal to sum_s extant_tips[s].size()
	
			
	// create additional tips, by splitting existing tips at each step (turning the split parent tip into a node)
	long Nedges 					= 0;
	double time 					= 0;
	long Nevents					= 0; // number of births/deaths/transitions so far
	double total_rate				= INFTY_D;
	double total_birth_rate			= state_birth_rates[clade2state[root]];
	double total_death_rate			= state_death_rates[clade2state[root]];
	double total_transition_rate 	= state_total_transition_rates[clade2state[root]];
	double equilibrium_time 		= INFTY_D;
	double initial_growth_rate 		= NAN_D; // keep track of the net growth rate (birth rate - death rate) at the beginning of the simulation
	std::vector<long> Ntransitions(Nstates*Nstates,0); // keep track of the number of transitions between each pair of states
	std::vector<long> Nbirths(Nstates,0);
	std::vector<long> Ndeaths(Nstates,0);
	long stip=0, clade=0, tip_state=0, old_state=0, child_state=0;
	while((NextantTips>0) && ((max_tips<=0) || ((coalescent ? NextantTips : Ntips)<max_tips)) && ((max_time<=0) || (time+1/total_rate<max_time)) && ((max_time_since_equilibrium<0) || (time-equilibrium_time+1/total_rate<max_time_since_equilibrium))){
		// determine time of next speciation/extinction/transition event
		// prevent deaths if only one tip is left
		const double restricted_death_rate = (no_full_extinction && (NextantTips<=1) ? 0 : total_death_rate);
		if(std::isnan(initial_growth_rate)) initial_growth_rate = total_birth_rate - restricted_death_rate;
		if(((total_birth_rate - restricted_death_rate)*initial_growth_rate<0) && (equilibrium_time>time)) equilibrium_time = time; // first crossing of equilibrium state, so keep record
		total_rate = total_birth_rate + restricted_death_rate + (all_transitions ? total_transition_rate : 0.0);
		
		// draw random (exponentially distributed) time lag to next event
		time += random_exponential_distribution(total_rate);
		++Nevents;
		
		const bool cladogenic = random_bernoulli((total_birth_rate+restricted_death_rate)/total_rate);
		if(cladogenic){
			// speciation or extinction event, decide which of the two
			const bool birth = ((restricted_death_rate==0) || random_bernoulli(total_birth_rate/(total_birth_rate+restricted_death_rate)));
						
			// randomly pick an existing tip to split or kill, proportionally to their state-dependent birth or death rates
			// returns the state of the picked tip, and its index in extant_tips[tip_state][]
			random_int_from_pools(extant_tips, (birth ? state_birth_rates : state_death_rates), (birth ? total_birth_rate : total_death_rate), tip_state, stip);
			clade = extant_tips[tip_state][stip];
			clade2end_time[clade] = time;
			const double edge_length = clade2end_time[clade]-(clade==root ? 0 : clade2end_time[clade2parent[clade]]);
		
			// update total birth & death rates for the removal of the chosen clade from the pool of tips
			total_birth_rate -= state_birth_rates[tip_state];
			total_death_rate -= state_death_rates[tip_state];
				
			// update total transition rate for the removal of the chosen clade from the pool of tips
			total_transition_rate -= state_total_transition_rates[tip_state];

			if(birth){
				// determine number of splits (typically Nsplits, but may be smaller to prevent complete extinction)
				if(max_tips>0) next_Nsplits = min(1+max_tips-long(coalescent ? NextantTips : Ntips), max(2l, Nsplits)); // cap Ntips if needed to stay within limits
				
				// remove tip from pool of extant tips (will remain in pool of clades, essentially as an internal node)
				remove_item_from_vector(extant_tips[tip_state], stip);
				
				// keep track of when this clade split
				if(include_birth_times) birth_times.push_back(time);

				// split chosen tip into next_Nsplits daughter-tips & create new edges leading into those tips
				// also choose the state of each child
				for(long c=0; c<next_Nsplits; ++c){
					// choose child state
					if(all_transitions){
						// child state is the same as parent, since transitions are modelled separately
						child_state = tip_state;
					}else{
						// pick random state for this child according to Markov transition rates
						child_state = get_next_Mk_state(transition_exponentiator,scratch_exp,edge_length/exponentiation_rescaling,tip_state);
					}
					// add child to lookup tables
					// its clade index is given by the current Nclades (since we append it to the end of the clades tables)
					tree_edge.push_back(clade);
					tree_edge.push_back(Nclades);
					extant_tips[child_state].push_back(Nclades);
					clade2parent.push_back(clade);
					clade2end_time.push_back(-1);
					clade2state.push_back(child_state);
					++Nclades;
					
					// update total birth & death rates
					total_birth_rate += state_birth_rates[child_state];
					total_death_rate += state_death_rates[child_state];
					
					// update total transition rate
					total_transition_rate += state_total_transition_rates[child_state];
				}
				
				// update some summary variables
				NextantTips += next_Nsplits-1;
				Nedges 		+= next_Nsplits;
				Ntips 		+= next_Nsplits-1;
				Nbirths[tip_state] += 1;					
				
			}else{
				// kill chosen tip (remove from pool of extant tips)
				// note that it still remains a tip (e.g. listed in clade2state[]), but it can't diversify or die anymore
				
				// remove from pool of extant tips
				remove_item_from_vector(extant_tips[tip_state], stip);

				// keep track of when this clade died
				if(include_death_times) death_times.push_back(time);

				// update some summary variables
				--NextantTips;
				Ndeaths[tip_state] += 1;
			}

		}else{
			// transition event
			// randomly pick an existing extant tip to transition
			random_int_from_pools(extant_tips, state_total_transition_rates, total_transition_rate, old_state, stip);
			clade = extant_tips[old_state][stip];
			
			// randomly pick new state
			const long new_state = get_next_Mk_state(Nstates,transition_matrix,state_total_transition_rates[old_state],old_state);
			clade2state[clade] 	 = new_state;

			// move tip into the appropriate pool of extant_tips
			remove_item_from_vector(extant_tips[old_state], stip);
			extant_tips[new_state].push_back(clade);

			// update total transition rate
			total_transition_rate += state_total_transition_rates[new_state] - state_total_transition_rates[old_state];
			
			// update total birth & death rate
			total_birth_rate += state_birth_rates[new_state] - state_birth_rates[old_state];
			total_death_rate += state_death_rates[new_state] - state_death_rates[old_state];
			
			// keep record of this transition event (old_state-->new_state)
			Ntransitions[old_state*Nstates + new_state] += 1;			
		}
		// abort if the user has interrupted the calling R program
		if(Nevents%100 == 0) Rcpp::checkUserInterrupt();
	}
	
	if(Ntips<1){
		// something went wrong (there should always be at least one tip
		return Rcpp::List::create(Rcpp::Named("success") = false, Rcpp::Named("error")="Something went wrong: Generated tree is empty");
	}else if((NextantTips<1) && coalescent){
		// the tree seems to have gone extinct, and only the coalescent tree was requested (i.e. we would need to return an empty tree)
		if((!no_full_extinction) && (max_death_rate_pc>0)){
			// full extinction is a plausible scenario since we did not try to prevent extinctions and death_rates were positive
			return Rcpp::List::create(Rcpp::Named("success") = false, Rcpp::Named("error")="Tree went fully extinct");
		}else{
			return Rcpp::List::create(Rcpp::Named("success") = false, Rcpp::Named("error")="Something went wrong: Generated tree has no extant tips");
		}
	}
		
	// add a small dt at the end to make all edges non-zero length
	time += random_exponential_distribution(total_rate);
	if(max_time>0) time = min(time, max_time); // prevent going past max_time
	if(max_time_since_equilibrium>=0) time = min(time, max_time_since_equilibrium+equilibrium_time); // prevent going past max_time_since_equilibrium

	// translate clade states to clade-specific birth & death rates
	std::vector<double> clade2birth_rate_pc(Nclades), clade2death_rate_pc(Nclades);
	for(long c=0; c<Nclades; ++c){
		clade2birth_rate_pc[c] = state_birth_rates[clade2state[c]];
		clade2death_rate_pc[c] = state_death_rates[clade2state[c]];
	}
	
	// flatten state-dependent pools of extant tips, into a single list
	std::vector<long> all_extant_tips(NextantTips);
	for(long state=0, tip=0; state<Nstates; ++state){
		for(long stip=0; stip<extant_tips[state].size(); ++stip){
			all_extant_tips[tip++] = extant_tips[state][stip];
		}
	}

	// finalize tree (make valid phylo structure, make coalescent if needed)
	std::vector<double> edge_length, birth_rates_pc, death_rates_pc;
	std::vector<long> new2old_clade;
	double root_time;
	aux_finalize_generated_random_tree(	time,
										as_generations,
										coalescent,
										include_rates,
										clade2end_time,
										clade2birth_rate_pc,
										clade2death_rate_pc,
										(max_death_rate_pc>0),
										Ntips,
										Nclades,
										Nedges,
										root,
										root_time,
										tree_edge,
										edge_length,
										all_extant_tips,
										new2old_clade,
										birth_rates_pc,
										death_rates_pc);
	long Nnodes = Nclades - Ntips;
	
	// update clade2state to new clade indices
	const std::vector<long> old_clade2state(clade2state);
	clade2state.resize(Nclades);
	for(long clade=0; clade<Nclades; ++clade) clade2state[clade] = old_clade2state[new2old_clade[clade]];
	
	return Rcpp::List::create(	Rcpp::Named("success") 			= true,
								Rcpp::Named("tree_edge") 		= Rcpp::wrap(tree_edge),
								Rcpp::Named("edge_length") 		= Rcpp::wrap(edge_length),
								Rcpp::Named("Nnodes") 			= Nnodes,
								Rcpp::Named("Ntips") 			= Ntips,
								Rcpp::Named("Nedges") 			= Nedges,
								Rcpp::Named("root")				= root, // this is actually guaranteed to be = Ntips
								Rcpp::Named("Nbirths")			= Nbirths,
								Rcpp::Named("Ndeaths")			= Ndeaths,
								Rcpp::Named("Ntransitions")		= Ntransitions,
								Rcpp::Named("root_time")		= root_time,
								Rcpp::Named("final_time")		= time,
								Rcpp::Named("equilibrium_time")	= equilibrium_time,
								Rcpp::Named("clade_states")		= clade2state,
								Rcpp::Named("birth_times")		= Rcpp::wrap(birth_times),
								Rcpp::Named("death_times")		= Rcpp::wrap(death_times),
								Rcpp::Named("birth_rates_pc")	= Rcpp::wrap(birth_rates_pc),	// birth_rates_pc[c] will be the per-capita birth rate of clade c (prior to its extinction or splitting)
								Rcpp::Named("death_rates_pc")	= Rcpp::wrap(death_rates_pc));	// death_rates_pc[c] will be the per-capita death rate of clade c (prior to its extinction or splitting)
}








