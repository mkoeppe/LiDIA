// -*- C++ -*-
//==============================================================================================
//
//	This file is part of LiDIA --- a library for computational number theory
//
//	Copyright (c) 1994--2001 the LiDIA Group.  All rights reserved.
//
//	See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//----------------------------------------------------------------------------------------------
//
//	$Id$
//
//	Author	: Michael Jacobson (MJ)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_QO_SIEVE_H_GUARD_
#define LIDIA_QO_SIEVE_H_GUARD_


#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_SORT_VECTOR_H_GUARD_
# include	"LiDIA/sort_vector.h"
#endif
#ifndef LIDIA_MATRIX_H_GUARD_
# include	"LiDIA/matrix.h"
#endif
#ifndef LIDIA_INDEXED_HASH_TABLE_H_GUARD_
# include	"LiDIA/indexed_hash_table.h"
#endif
#ifndef LIDIA_PARTIAL_RELATION_H_GUARD_
# include	"LiDIA/number_fields/partial_relation.h"
#endif
#ifndef LIDIA_QUADRATIC_FORM_H_GUARD_
# include	"LiDIA/quadratic_form.h"
#endif
#ifndef LIDIA_QUADRATIC_NUMBER_STANDARD_H_GUARD_
# include	"LiDIA/quadratic_number_standard.h"
#endif
#ifndef LIDIA_ARITH_INL_GUARD_
# include       "LiDIA/arith.inl"
#endif
#include	<fstream>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class qi_class;


//
// Class: qo_sieve
//
// This class sieves quadratic polynomials and produces relations for the
//    subexponential class group algorithms based on the MPQS factoring
//    algorithm.  Depending on the initialization, it either uses regular
//    MPQS or SIQS.
// The relations are stored in a matrix as they are generated, and the user
//    can access them one at a time.  If the user asks for a relation and
//    the matrix has some still available, no sieving is done.
//

#define CAND_NUMBER 10000
#define MAX_FB 9001
//#define MAX_SIEB 500010
#define MAX_SIEB 10000010  // MJJ

typedef unsigned char SIEBTYP;

class qo_sieve
{
private:

	bigint Delta; // discriminant of polynomials to sieve
	int BachBound; // Bach's bound for a generating system

	bool is_init; // true if sieve is initialized
	bool is_si; // type of polynomial generation

	// sieveing parameters
	int size_FB; // size of the factor base
	int M; // sieve radius (-M, M)
	int P_ONCE; // number of distinct primes to form leading coeff.
	int P_TOTAL; // number of primes from which to choose
	int smallstart; // first prime to sieve with
	int POLY; // number of B's for one A coefficient
	long lpbound; // bound for large prime relations
	double Tval; // tolerance value
	double d_wurz; // double approx of sqrt(Delta)
	double size_A; // size of leading coefficient
	int pBound; // max prime for non-self-init

	sort_vector< bigint > D_primes; // factor base primes which divide
	// Delta

	// relation array
	matrix< lidia_size_t > relation_array; // array of previously computed
	// relations
	base_vector< quadratic_number_standard > minima_array; // the corresponding minima


	// sieving arrays
	int *FB; // norms of prime ideals in factor base
	bigint *FB_b; // b values of the polynomials
	SIEBTYP *LOGP; // logarithms of the norms
	int *SQRTkN; // square root of Delta mod the norms
	int *START1; // sieving starting points for each norm
	int *START2; // sieving starting points (second ones)

	SIEBTYP *sieb; // the sieve array
	SIEBTYP *ende; // the last entry in the sieve array
	int *CANDIDATE; // list of candidate x values in the sieve array

	// second set of starting points (when mixing si and normal)
	int *START1_2; // sieving starting points for each norm
	int *START2_2; // sieving starting points (second ones)


	// self-initialization variables
	int start_fb; // first prime to take for leading coefficient
	int bin_index; // index of current small prime product
	int max_bin; // maximum possible with current paramters
	int index_i; // index of current B value

	// self-initialization arrays
	sort_vector< int > Q_prime_glob; // all primes which can be used to
	// construct the leading coeff

	int *a_inv; // inverse of leading coeffiecient mod each norm
	int **vorb;
	bigint *BG;


	// large prime relation files
	char faktor[LINE_LEN];
	std::ofstream lp_file;
	char PART[50];
	char TEMP[50];
	char ZUSA[50];

	// hash tables for avoiding duplicate relations and polynomials
	hash_table< bigint > htable;
	hash_table< bigint > used_polys;
	hash_table< int > lp_table;


	// polynomial info
	bool next_is_si; // true if next poly to sieve was generated with si
	quadratic_form F_si; // the polynomial to sieve from si
	math_vector< int > e_si; // vector over FB representing F_si;

	quadratic_form F_ra; // the polynomial to sieve from random generation
	math_vector< int > e_ra; // vector over FB representing F_ra;

	sort_vector< lidia_size_t > Q_prime_si; // indices of primes dividing a
	sort_vector< lidia_size_t > Q_prime_ra; // indices of primes dividing a




	//
	// PRIVATE FUNCTIONS
	//

	void qo_get_name(char *s);

	void new_polynomial(lidia_size_t, int);
	void new_polynomial_si();

	void compute_roots();

	void qo_sieve_interval();

	bool test_candidates(int, lidia_size_t, bool);

	bool test_single(math_vector< int > &, quadratic_number_standard &);




public:

	// constructor and destructor
	qo_sieve();
	~qo_sieve();



	//
	// initialization
	//

	void init(bigint &, base_vector< qi_class > &, double, int, int,
		  int, int, int, long);
	void init_self_initialization();
	void reset();
	void restart();



	//
	// state queries
	//

	bool is_initialized();
	bool is_si_initialized();
	void no_self_init();

	bigint discriminant_used();
	int FB_size();
	void dense_bound(int &);




	//
	// generate next polynomial
	//

	void next_polynomial();
	void next_polynomial(lidia_size_t);
	void next_polynomial(const qi_class &);


	//
	// sieve current polynomial (all solutions)
	//

	void sieve();



	//
	// sieve current polynomial (first solution only)
	//

	bool sieve_one(lidia_size_t, math_vector< lidia_size_t > &,
		       quadratic_number_standard &);
	bool sieve_one(const qi_class &, math_vector< lidia_size_t > &,
		       quadratic_number_standard &);
	bool sieve_one(const qi_class &);



	//
	// large prime handling
	//

	lidia_size_t count_lp_relations();
	lidia_size_t get_lp_relations();


	//
	// access to relations
	//

	bool get_relation(math_vector< lidia_size_t > &, quadratic_number_standard &);


};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifndef LIDIA_QI_CLASS_H_GUARD_
# include	"LiDIA/qi_class.h"
#endif



#endif	// LIDIA_QO_SIEVE_H_GUARD_
