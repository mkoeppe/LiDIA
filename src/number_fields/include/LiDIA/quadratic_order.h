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
//	Changes	:
//
//==============================================================================================


#ifndef LIDIA_QUADRATIC_ORDER_H_GUARD_
#define LIDIA_QUADRATIC_ORDER_H_GUARD_



#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_XBIGFLOAT_H_GUARD_
# include	"LiDIA/xbigfloat.h"
#endif
#ifndef LIDIA_BIGFLOAT_INT_H_GUARD_
# include	"LiDIA/bigfloat_int.h"
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
#ifndef LIDIA_RATIONAL_FACTORIZATION_H_GUARD_
# include	"LiDIA/rational_factorization.h"
#endif
#ifndef LIDIA_QUADRATIC_NUMBER_STANDARD_H_GUARD_
# include	"LiDIA/quadratic_number_standard.h"
#endif
#ifndef LIDIA_QUADRATIC_NUMBER_POWER_PRODUCT_H_GUARD_
# include	"LiDIA/quadratic_number_power_product.h"
#endif
#ifndef LIDIA_QO_SIEVE_H_GUARD_
# include	"LiDIA/number_fields/qo_sieve.h"
#endif
#ifndef LIDIA_QO_LIST_H_GUARD_
# include	"LiDIA/number_fields/qo_list.h"
#endif
#ifndef LIDIA_QO_UTIL_H_GUARD_
# include	"LiDIA/number_fields/qo_util.h"
#endif
#ifndef LIDIA_QI_CLASS_H_GUARD_
# include	"LiDIA/qi_class.h"
#endif
#ifndef LIDIA_INFO_H_GUARD_
# include	"LiDIA/info.h"
#endif
#ifndef LIDIA_PRIME_LIST_H_GUARD_
# include	"LiDIA/prime_list.h"
#endif

#include	<iomanip>
#include	<cstdlib>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



//
// Class: quadratic_order
//
// This class represents a quadratic order Z + (Delta + sqrt(Delta))/2 Z,
//    where Delta is the discriminant of the order.
// A quadratic order is represented internally by its discriminant.  The
//    following invariants are stored as they are computed:
//
//    R = regulator
//    h = class number
//    L = approximation of L(1, X)
//    Cfunc = approximation of C(Delta)
//    CL = vector of class group invariants
//    gens = set of generators of each cyclic subgroup
//
// The following information is stored for use in some computations:
//
//     hfact = integer factorization of h
//     disc_fact = integer factorization of Delta
//     fact_base = factor base of prime ideals
//     U = transformation matrix from generators to factor base
//     UINV = transformation matrix form factor base to generators
//     prin_list = list of principal ideals
//

#define RBJTB 7			// for regulator, use BJT if D < 10^7
#define RSHANKSB 13		// for regulator, use shanks if D < 10^13
#define CNBJT_RB 3 		// for class number, use BJT if D < 10^3
#define CNSHANKS_RB 11		// for class number, use Shanks if D < 10^11
#define CGBJT_RB 9 		// for class group, use BJT if D < 10^9
#define CGRAND_RB 12		// for class group, use randexp if D < 10^12

#define CNBJT_IB 11 		// for class number, use BJT if |D| < 10^11
#define CNSHANKS_IB 16		// for class number, use Shanks if |D+ < 10^16
#define CGBJT_IB 11 		// for class group, use BJT if |D| < 10^11
#define CGSHANKS_IB 14		// for class group, use Shanks if |D| < 10^14



class quadratic_order
{
private:

	bigint Delta; // discriminant
	bigfloat R; // floating point approximation of regulator
	xbigfloat R_xbig; // same with xbigfloat
	long R_xbig_prec; // R_xbig is an absolute R_xbig_prec approx.
	bigint h; // class number
	bigfloat L; // approximation of L(1, \chi)
	bigfloat Cfunc; // approximation of C(Delta)
	base_vector< bigint > CL; // invariants of the class groups
	base_vector< qi_class > gens; // system of generators of CL

	quadratic_number_power_product fund_unit; // fundamental unit

	rational_factorization hfact; // factorization of the class number
	rational_factorization disc_fact; // factorization of Delta
	base_vector< qi_class > fact_base; // factor base used to compute CL
	base_vector< lidia_size_t > contributors; // factor base elements which
	// contribute to CL

	math_vector< quadratic_number_standard > minima_qn; // and minima itselves

	matrix< long > AMAT; // matrix of all relations generated
	trans_matrix TR; // transformation matrix

	bigint DET, fu_mult;

	qo_sieve sieve; // used for sieve-bases relation generation
	prime_list PL; // list of prime numbers
	prime_list SPL; // list of small prime numbers

	static int info; // whether to print times

public:

	matrix< bigint > RHNF; // hnf of relation matrix
	matrix< bigint > U; // transformation from gens to FB
	matrix< bigint > UINV; // transformation from FB to gens


	bool subused; // true if the class group was computed
				// using the subexponential algorithm

	bigint nucomp_bound; // bound for nucomp and nudupl (if Delta < 0)

	bigfloat FI; // L(1) approximation using p < 2Q (Bach method)
	long Q;

	bigfloat min_log; // minimum regulator possible
	static int do_verify; // whether to verify computations


private:
	//
	// PRIVATE FUNCTIONS
	//

	void set_transformations();

	// regulator
	void rBJT();
	void rshanks();
	bigfloat find_hstar(const bigfloat & hR, const bigfloat & E,
			    const bigfloat & u, const bigfloat & Lval);

	// class number, Shanks method
	void cns_imag();
	void cns_real();

	// class group, BJT method
	void cgBJT_imag();
	void cgBJT_real();

	// class group, Shanks method
	void cgs_imag(const bigfloat & tFI, const bigfloat & tF, bigint & h1);
	void cgs_real(const bigfloat & tFI, bigint & h1);

	// class group, when h is known
	void cgh_imag();
	void cgh_real();

	// class group, subexponential routines (file qo_subexp.c)
	void rBJT_bound(const bigfloat &);
	void qo_get_fb(int & size_FB);
	void generate_free_relations(lidia_size_t &, sort_vector< lidia_size_t > &);
	void qo_hnf_and_det(matrix< bigint > & A1, sort_vector< lidia_size_t > & S,
			    xbigfloat & LL, bool input_reg = false);

	// Duellmann's algorithm
	void qo_read_parameters(int decimal_digits, int &size_FB);
	void qo_init_powers(int & k0, qi_class ***powers,
			    quadratic_number_standard ***mins,
			    base_vector< lidia_size_t > & power_index);
	bool factor_over_FB(qi_class & AA, math_vector< lidia_size_t > &v);
	void qo_relations_randexp(int size_FB, int n, qi_class **powers,
				  quadratic_number_standard **mins,
				  base_vector< lidia_size_t > & power_index,
				  sort_vector< lidia_size_t > & S,
				  matrix< bigint > & A1,
				  lidia_size_t & curr_A,
				  lidia_size_t & curr_A1);
	void cg_randexp();



	// MPQS algorithm
	void qo_read_parameters(int decimal_digits, double &T, int &M, int &size_FB,
				int &P_ONCE, int &P_TOTAL, int &smallstart, int &POLY);
	void qo_relations_mpqs(int size_FB,
			       int n,
			       sort_vector< lidia_size_t > & S,
			       matrix< bigint > & A1,
			       lidia_size_t & curr_A,
			       lidia_size_t & curr_A1,
			       bool self_init);
	void cg_mpqs(bool self_init, bool large_primes);

	// representing over factor base and generators
	math_vector< bigint > represent_over_generators(math_vector< bigint > & vec);





	// factoring the discriminant using the 2-Sylow subgroup of the class group
	void factor_imag();
	void factor_real();
	void find_ambiguous(qi_class_real & AMB, sort_vector< bigint > & facts,
			    int & numfacts);

	// computing smallest power of diagonal element for relation
	bigint exact_power(const bigint & omult, const qi_class & G,
			   indexed_hash_table< ideal_node > & RT,
			   indexed_hash_table< ideal_node > & QT,
			   lidia_size_t numRpr, long & r, long & q);
	bigint exact_power_real(const bigint & omult, const qi_class & G,
				indexed_hash_table< ideal_node > & RT,
				indexed_hash_table< ideal_node > & QT,
				lidia_size_t numRpr, long & r, long & q,
				qi_class_real & Gstep);



public:

	long prec; // floating point approximation precision
	long xprec; // xbigfloat precision

	hash_table< qi_class_real > prin_list; // principal ideal list

	static qo_list& qo_l(); // list of current quadratic orders
	static quadratic_order zero_QO; // order with discriminant zero

	static bool special; // FOR TABLE GENERATION
	static bool special2; // FOR TABLE GENERATION
	static bool noCL; // FOR TABLE GENERATION
	static char matfile[]; // FOR TABLE GENERATION


	//
	// constructors and destructor
	//

	quadratic_order();
	quadratic_order(const long D);
	quadratic_order(const bigint & D);
	quadratic_order(const quadratic_order & QO);

	~quadratic_order();



	//
	// initialization
	//

	static void verbose(int state);
	static void verification(int level);



	//
	// assignments
	//

	bool assign(const long D);
	bool assign(const bigint & D);
	void assign(const quadratic_order & QO2);

	quadratic_order & operator = (const quadratic_order & QO);



	//
	// access functions
	//

	bigint discriminant() const;
#ifndef HEADBANGER
	bigint nu_bound() const;
#endif


	//
	// comparisons
	//

	bool is_zero() const;
	bool is_equal(const quadratic_order & QO2) const;
	bool is_subset(quadratic_order & QO2);
	bool is_proper_subset(quadratic_order & QO2);

	friend bool operator == (const quadratic_order & QO1,
				 const quadratic_order & QO2);
	friend bool operator != (const quadratic_order & QO1,
				 const quadratic_order & QO2);

	friend bool operator <= (quadratic_order & QO1, quadratic_order & QO2);
	friend bool operator < (quadratic_order & QO1, quadratic_order & QO2);
	friend bool operator >= (quadratic_order & QO1, quadratic_order & QO2);
	friend bool operator > (quadratic_order & QO1, quadratic_order & QO2);

	friend bool operator ! (const quadratic_order & QO);


	//
	// basic functions
	//

	friend bool is_quadratic_discriminant(const long D);
	friend bool is_quadratic_discriminant(const bigint & D);
	bool is_imaginary() const;
	bool is_real() const;
#ifndef HEADBANGER
	bool is_R_computed() const;
	bool is_h_computed() const;
	bool is_L_computed() const;
	bool is_C_computed() const;
	bool is_CL_computed() const;
	bool is_subexp_computed() const;
	bool is_fundamental_unit_computed () const; //MM
#endif

	friend void swap(quadratic_order & Q1, quadratic_order & Q2);


	//
	// high level functions
	//

	int bach_bound();
	int bach_bound(const bigint & f);

	bigint conductor();
	bool is_maximal();
	void maximize();
	void maximize(quadratic_order & max_ord);



	//
	// L functions, Littlewood indices, C function
	//

	friend int kronecker(const bigint & D, const bigint & p);

	// optimal value of Q for h* < h < 2h* approximation
	long generate_optimal_Q();
	long get_optimal_Q();

	// optimal value of Q to prove that h = 3
	long generate_optimal_Q_cnum(const bigfloat & h2, const bigfloat & t);
	long get_optimal_Q_cnum();

	// optimal value of Q to compute C(Delta)
	long generate_optimal_Q_cfunc();
	long get_optimal_Q_cfunc();

	// truncated product estimate of C(Delta)
	bigfloat estimate_C(const long nQ);

	// truncated product estimate of L(s)
	bigfloat estimate_L(const int s, const long nQ);

	// Bach's method for estimating L(1)
	bigfloat estimate_L1(const long nQ);
	bigfloat estimate_L1_error(const long nQ) const;


	// L functions, C function, and Littlewood indices
	bigfloat Lfunction();
	bigfloat LDfunction();
	bigfloat LLI();
	bigfloat LLI_D();
	bigfloat ULI();
	bigfloat ULI_D();
	bigfloat Cfunction();
	bigfloat_int Cfunction_bigfloat(int pmult);

	bigfloat LLI(const bigint &);
	bigfloat LLI(const bigint &, const bigfloat &);
	bigfloat ULI(const bigint &);
	bigfloat ULI(const bigint &, const bigfloat &);
	bigfloat Cfunction(const bigint &, const rational_factorization &,
			   int pmult = 0);
	bigfloat Cfunction(const bigint &, const bigfloat &,
			   const rational_factorization &, int pmult = 0);



	//
	// regulator, class number, class group
	//

	bigfloat regulator();
	xbigfloat regulator(long k); // MM
	bigint class_number();
	base_vector< bigint > class_group();

#ifndef HEADBANGER
	// forcing subexp algorithm
	bigfloat regulator(bool usesub);
	base_vector< bigint > class_group(bool usesub);
#endif

	// O(D^1/4) regulator computation
	bigfloat regulator_BJT();

	// O(D^1/5) regulator computation (uses estimate of L(1))
	bigfloat regulator_shanks();


	// verification
	bool verify_regulator();
	// < MM >
	bool could_be_regulator_multiple (const xbigfloat & l);
	// < /MM >

	// O(D^1/5) class number computation
	bigint class_number_shanks();



	// from L(1) estimate and BJT algorithm, class number not known
	base_vector< bigint > class_group_BJT();

	// from L(1) estimate and my algorithm, class number not known
	base_vector< bigint > class_group_shanks();

#ifndef HEADBANGER
	// class number known
	base_vector< bigint > class_group_h();
#endif

	// subexponential class group algorithm
	base_vector< bigint > class_group_randexp();
	base_vector< bigint > class_group_mpqs();
	base_vector< bigint > class_group_siqs();
	base_vector< bigint > class_group_lp();

#ifndef HEADBANGER
	void class_group_siqs(int new_size_FB, int new_M = 0, lidia_size_t c = 0, bool in_reg = false);
#endif

	// verification
	int count_primes_to_verify(int BachBound);
	bool verify_class_group(int level);
	bool verify_class_group();
	bool verify_cg(bool use_sieve);

#ifndef HEADBANGER
	long cg_rel_test(lidia_size_t n,
			 int size_FB_in,
			 int M_in,
			 int & size_FB_used,
			 int & M_used,
			 int & PTOTAL_used);


	// principality test (subexponential)
	bool is_in_lattice(const qi_class_real & AA);
	bool is_in_lattice(math_vector< bigint > & vec);
	bool is_in_lattice(const qi_class_real & AA, bigfloat & dist);
#endif

	// representing over FB and generators
	math_vector< bigint > represent_over_FB(const qi_class & AA,
						quadratic_number_standard & q);
	math_vector< bigint > represent_over_generators(qi_class &AA);


	math_vector< bigint > represent_over_FB_sieve(const qi_class & AA,
						      quadratic_number_standard & q);
	bool represent_over_FB_sieve(const qi_class & AA);
	math_vector< bigint > represent_over_generators_sieve(qi_class &AA);


	// misc
	rational_factorization factor_h();
	bigint exponent();
	int p_rank(const bigint & p);
	base_vector< qi_class > generators();
	rational_factorization factor_discriminant();

	// fundamental unit
	const quadratic_number_power_product & fundamental_unit ();

	void compute_HNF_and_R(char *inname, char *outname);


	//
	// input/output
	//

	friend std::istream & operator >> (std::istream & in, quadratic_order & QO);
	friend std::ostream & operator << (std::ostream & out, const quadratic_order & QO);


	// < MM >
	//
	// The interface for the following functions is not yet stable
	// and may change in the future.
	//



	// L(1, chi_D) approximation
	//
private:

	void bachs_L1_error_triplet(const bigint & Q, xbigfloat & A, xbigfloat
				    & B) const;

	bigint number_of_terms (const bigint & Delta, long k, int info = 0) const;

	void ell (xbigfloat & l, unsigned long n, const bigint & Delta, long k, int info = 0) const;


public:

	void relative_L1chi_approximation (xbigfloat & L, long k, int info = 0) const;

	xbigfloat relative_L1chi_approximation (long k, int info = 0) const;


private:

	static void test_L1chi();


public:

	void approximate_exp_l_n_delta(xbigfloat & L, const bigint & Delta, int info = 0) const;


	//
	// class number verification
	//

	bool verify_h(const bigint & Delta, const bigint & h, int info = 0);

	bool verify_h(xbigfloat & L, const bigint & Delta, const bigint & h, int info = 0);

	bool verify_hR(const bigint & Delta, const bigint & h, const xbigfloat & r, int info = 0);

	bool verify_hR(xbigfloat & L, const bigint & Delta, const bigint & h, const xbigfloat & r, int info = 0);


	//
	// lower bound on the regulator
	//
	long lower_regulator_bound(const bigint & Delta);
	long lower_regulator_bound(const bigint & Delta, const bigint & h, int info = 0);
	long lower_regulator_bound(xbigfloat & L, const bigint & Delta, const bigint & h, int info = 0);


	// DO NOT USE THE FOLLOWING FUNCTIONS ANYMORE.
	//
	// The functions below will be removed in the future,
	// because they are member functions of quadratic_number_power_product
	// now. There are still here only for debugging purposes.
	//

	//
	// approximation of Ln of power products of quadratic numbers
	//

	void absolute_Ln_approximation (xbigfloat & l,
					base_vector< xbigfloat > & log_q,
					base_vector< long > & prec_log_q,
					const bigint_matrix & M,
					const base_vector< quadratic_number_standard > & q,
					lidia_size_t j,
					long k,
					int info = 0);

	void absolute_Ln_approximation_with_check (xbigfloat & l,
						   base_vector< xbigfloat > & log_q,
						   base_vector< long > & prec_log_q,
						   const bigint_matrix & M,
						   const base_vector< quadratic_number_standard > & q,
						   lidia_size_t j,
						   long k,
						   int info = 0);

	void absolute_log_approximation (bigfloat & l,
					 const matrix< bigint > & M,
					 const base_vector< quadratic_number_standard > & q,
					 lidia_size_t j,
					 int info = 0);

	void absolute_Ln_approximation (bigfloat & l,
					const matrix< bigint > & M,
					const base_vector< quadratic_number_standard > & q,
					lidia_size_t j,
					int info = 0);

	void absolute_Ln_approximation (xbigfloat & l,
					base_vector< xbigfloat > & log_q,
					base_vector< long > & prec_log_q,
					const base_vector< bigint > & exponents,
					const base_vector< quadratic_number_standard > & q,
					long k,
					int info = 0);

	void absolute_Ln_approximation (xbigfloat & l,
					const base_vector< quadratic_number_standard > & q,
					const base_vector< bigint > & e,
					long k,
					int info = 0);

	void relative_Ln_approximation (xbigfloat & l,
					const base_vector< quadratic_number_standard > & q,
					const base_vector< bigint > & e,
					long k,
					long m,
					int info = 0);

	void reduce_modulo_regulator (bigint & q,
				      const base_vector< quadratic_number_standard > & base_w,
				      const base_vector< bigint > & exp_w,
				      base_vector< xbigfloat > & log_base_w,
				      base_vector< long > & prec_log_base_w,
				      const base_vector< quadratic_number_standard > & base_r,
				      const base_vector< bigint > & exp_r,
				      base_vector< xbigfloat > & log_base_r,
				      base_vector< long > & prec_log_base_r,
				      long br);

	//
	// finding a generating unit from a set of units
	//
	void remove_one_from_unit_array (base_vector< xbigfloat > & log_q,
					 base_vector< long > & prec_log_q,
					 bigint_matrix & M,
					 const base_vector< quadratic_number_standard > & q,
					 long m,
					 int info = 0);

	bool is_one(base_vector< xbigfloat > & log_q,
		    base_vector< long > & prec_log_q,
		    const base_vector< quadratic_number_standard > & q,
		    const base_vector< bigint > & tmp_exponents,
		    lidia_size_t m,
		    int info = 0);

	int  is_generator (base_vector< xbigfloat > & log_q,
			   base_vector< long > & prec_log_q,
			   const matrix< bigint > & M,
			   const base_vector< quadratic_number_standard > & q,
			   lidia_size_t i,
			   lidia_size_t j,
			   const bigint & x,
			   const bigint & y,
			   const bigint & M1,
			   const bigint & M2,
			   long m,
			   int info = 0);

	void precompute_logarithm_bases (base_vector< xbigfloat > & log_q,
					 base_vector< long > & prec_log_q,
					 const bigint_matrix & M,
					 const base_vector< quadratic_number_standard > & q,
					 long k);

	long estimate_pp_accuracy (const matrix< bigint > & M,
				   long m,
				   long c);

	void cfrac_expansion_bounded_den (bigint & conv_num,
					  bigint & conv_den,
					  bigint num,
					  bigint den,
					  const bigint & S,
					  int info);

	void rgcd (bigint & x,
		   bigint & y,
		   bigint & M1,
		   bigint & M2,
		   base_vector< xbigfloat > & log_q,
		   base_vector< long > & prec_log_q,
		   const bigint_matrix & M,
		   const base_vector< quadratic_number_standard > & q,
		   lidia_size_t j1,
		   lidia_size_t j2,
		   long m,
		   int info = 0);

	void find_generating_unit (bigint_matrix & M,
				   const base_vector< quadratic_number_standard > & q,
				   long m,
				   int info = 0);


	bigfloat find_generating_unit2 (bigint_matrix & M,
					const base_vector< quadratic_number_standard > & q,
					long m,
					bigint & Delta);
};


// friend functions of quadratic order

bool operator == (const quadratic_order & QO1,
		  const quadratic_order & QO2);
bool operator != (const quadratic_order & QO1,
		  const quadratic_order & QO2);

bool operator <= (quadratic_order & QO1, quadratic_order & QO2);
bool operator < (quadratic_order & QO1, quadratic_order & QO2);
bool operator >= (quadratic_order & QO1, quadratic_order & QO2);
bool operator > (quadratic_order & QO1, quadratic_order & QO2);

bool operator ! (const quadratic_order & QO);


    //
    // basic functions
    //

bool is_quadratic_discriminant(const long D);
bool is_quadratic_discriminant(const bigint & D);

void swap(quadratic_order & Q1, quadratic_order & Q2);

    //
    // L functions, Littlewood indices, C function
    //

int kronecker(const bigint & D, const bigint & p);

    //
    // input/output
    //

std::istream & operator >> (std::istream & in, quadratic_order & QO);
std::ostream & operator << (std::ostream & out, const quadratic_order & QO);




#ifndef HEADBANGER

inline void
nugcd(bigint & u, bigint & d, bigint & v1, bigint & v3, const bigint & PEA_L)
{
	bigint t1, t3, q;
	bigint h, v1_x, v3_x, d_x, u_x;
	bool flag;

	remainder(v3, v3, d);
	if (v3.is_lt_zero())
		add(v3, v3, d);

	subtract(u, d, v3);
	if (abs(u) < abs(v3)) {
		v3.assign(u);
		flag = true;
	}
	else
		flag = false;

	v1.assign_one();
	u.assign_zero();

	// trivial case
	if (abs(v3) <= PEA_L) {
		if (flag)
			v3.negate();
		return;
	}

	// iteration
	u_x.assign(u);
	d_x.assign(d);
	v1_x.assign(v1);
	v3_x.assign(v3);

	do {
		div_rem(q, t3, d_x, v3_x);

		multiply(t1, q, v1_x);
		subtract(t1, u_x, t1);

		h.assign(u_x); u_x.assign(v1_x); v1_x.assign(t1); t1.assign(h);
		h.assign(d_x); d_x.assign(v3_x); v3_x.assign(t3); t3.assign(h);
	} while (abs(v3_x) > PEA_L);

	// store results
	if (u_x == v1) {
		u.assign(u_x); v1.assign(v1_x); d.assign(d_x); v3.assign(v3_x);
	}
	else if (u_x != u) {
		v1.assign(v1_x); u.assign(u_x); v3.assign(v3_x); d.assign(d_x);
	}

	// modify signs
	if (flag) {
		v1.negate();
		u.negate();
	}

	if (v1.is_lt_zero()) {
		v1.negate();
		v3.negate();
	}
}

#endif // HEADBANGER



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_QUADRATIC_ORDER_H_GUARD_
