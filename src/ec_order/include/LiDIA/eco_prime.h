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
//	Author	: Frank Lehmann, Markus Maurer, Volker Mueller
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_ECO_PRIME_H_GUARD_
#define LIDIA_ECO_PRIME_H_GUARD_



#ifndef LIDIA_BIGMOD_H_GUARD_
# include "LiDIA/bigmod.h"
#endif
#ifndef LIDIA_FF1_H_GUARD_
# include "LiDIA/ff1.h"
#endif
#ifndef LIDIA_FF2_H_GUARD_
# include "LiDIA/ff2.h"
#endif
#ifndef LIDIA_MEQ_PRIME_H_GUARD_
# include "LiDIA/meq_prime.h"
#endif
#ifndef LIDIA_DENSE_POWER_SERIES_H_GUARD_
# include "LiDIA/dense_power_series.h"
#endif
#ifndef LIDIA_FP_RATIONAL_FUNCTION_H_GUARD_
# include "LiDIA/Fp_rational_function.h"
#endif
#ifndef LIDIA_WEP_RAT_FUNCTION_H_GUARD_
# include "LiDIA/wep_rat_function.h"
#endif
#ifndef LIDIA_FP_POLY_MULTIPLIER_H_GUARD_
# include "LiDIA/Fp_poly_multiplier.h"
#endif
#ifndef LIDIA_BASE_VECTOR_H_GUARD_
# include "LiDIA/base_vector.h"
#endif
#ifndef LIDIA_POINT_H_GUARD_
# include "LiDIA/point.h"
#endif
#ifndef LIDIA_GF_ELEMENT_H_GUARD_
# include "LiDIA/gf_element.h"
#endif
#ifdef TIMING                    // to allow precise timings of all different
# include "LiDIA/timer.h"         // parts of the program
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class eco_prime : public meq_prime
{
private :

	// ***** types used in the class eco_prime *****

	typedef bigmod  	       ff_element;
	typedef Fp_polynomial        ff_pol;
	typedef Fp_rational_function ff_rat;
	typedef Fp_poly_modulus      ff_polmod;
	typedef Fp_poly_multiplier   ff_polmult;
	typedef ff1     	       ff1_element;
	typedef ff2  	               ff2_element;

        // the forward and friend declarations are required since 
        // PsiPowers needs access to the private type spcifier ff_pol.
        struct PsiPowers;
        friend struct PsiPowers;

	struct PsiPowers {
		ff_pol *pow1;
		ff_pol *pow2;
		ff_pol *pow3;
	};


	// ***** constants *****


	static const char UNKNOWN; 		      // values for sp_type
	static const char ALL_SUBGRPS_INV;
	static const char ONE_SUBGRP_INV;
	static const char TWO_SUBGRPS_INV;
	static const char NO_SUBGRP_INV;

public:

	static const char COMPUTE_SP_DEGREE; // values for degree_mode
	static const char DONT_COMPUTE_SP_DEGREE;
	static const char COMPUTE_SP_DEGREE_IF_ATKIN;
	static const char COMPUTE_SP_DEGREE_IF_ELKIES;
	static const char COMPUTE_SP_DEGREE_AND_ROOTS;

	static const char EV_RATIONAL_FUNCTION; // values for ev_strategy below
	static const char EV_DIVISION_POLYNOMIAL;
	static const char EV_RATIONAL_FUNCTION_TABLE;
	static const char EV_BABYSTEP_GIANTSTEP;
	static const char EV_FUNNY_BABYSTEP_GIANTSTEP;
	static const char EV_OPTIMAL;

	// in the highest bit, we store what coordinate is tested --> for 0, use
	// y-coordinate (default), for 1 use X-coordinate only.

	static const char TEST_X_COORDINATE; // test with X or Y-coordinate
	static const char TEST_Y_COORDINATE; // in Elkies part


private:
	
	// ***** class variables, statics *****

	bigint p; 		// characteristic of finite field
	bigint pn; 	 	  	// number of elements of finite field
	ff_element A, B; // elliptic curve over finite field
	ff_element jinv; // j - invariant of E=(A, B)

	elliptic_curve< gf_element > E_; // the curve, _ to mark member variables


	// ***** class variables, non statics *****

	ff_pol        meq_pol; // l - th modular polynomial
	ff_polmod     meq_pol_mod; // PolyModulus of it
	ff_element    ftau[2]; // roots of meq_pol
	char          sp_type; // splitting type of meq_pol
	char          degree_mode;
	udigit        sp_degree; // contains the splitting degree of meq_pol
	bool          sp_degree_known; // if sp_degree_known is true
	char          ev_strategy; // algorithm to find eigenvalues.
	unsigned int  schoof_bound; // for smaller Atkin primes, use Schoof
	unsigned int  atkin_bound; //  no splt.type computed for l >..
        unsigned int  lower_bound_for_BG;   // don't use BG, FBG algorithms for
        unsigned int  lower_bound_for_FBG;  // small l
	bool          use_corr_test; // use prob. correctness test iff true

	udigit              q; // pn mod l, l might be a prime power
	sort_vector< udigit > c; // contains c mod l, c is trace of Frobenius
	udigit              l;

	int info;


private:

	bool test_x_coordinate()
	{
		if (ev_strategy == EV_FUNNY_BABYSTEP_GIANTSTEP)
			return true;

		if (ev_strategy == EV_OPTIMAL) {
			if (l % 4 == 3 && l >= lower_bound_for_FBG)
				return true;
			else
				return false;
		}

		if (ev_strategy == EV_BABYSTEP_GIANTSTEP)
			return false;

		if (ev_strategy & 0x80)   // highest bit is 1 !!
			return true;
		else
			return false;
	}

	bool test_x_coordinate_uniquely()
	{
		return ((test_x_coordinate() && (l % 4 == 3)) || 
                    (ev_strategy == EV_FUNNY_BABYSTEP_GIANTSTEP 
                     && l >= lower_bound_for_FBG));
	}

	bool test_y_coordinate()
	{
		return !test_x_coordinate();
	}

	// ---------------------------------------------------------
	// now we have a lot of private functions used for Elkies
	// computation

	//******* File elkies/div_of_divpol.cc ******************************
	
	void compute_divisor_of_division_polynomial(ff_pol & div_of_div_pol,
						    ff_element  & E2a,
						    ff_element  & E2b);
	void guess_jltau (base_vector< ff_element > &);

	void coefficient_comparison(ff_pol &, const ff_element &,
				    const ff_element &, const ff_element &);

	void compute_coefficients_of_weierstrass_p (base_vector< ff_element > &,
						    lidia_size_t, const ff_element &,
						    const ff_element &);

	void compute_weierstrass_p (dense_power_series< ff_element > &,
				    const base_vector< ff_element > &, lidia_size_t);

	void compute_coeffseries(dense_power_series< ff_element > &,
				 const base_vector< ff_element > &,
				 const base_vector< ff_element > &,
				 const ff_element &, lidia_size_t);
	void precompute_table(dense_power_series< ff_element > * &,
			      const dense_power_series< ff_element > &, lidia_size_t);

	void set_odd_weierstrass(dense_power_series< ff_element >* const,
				 lidia_size_t);


	//****** File elkies/tildeE.cc ********************************

	int tildeEA (ff_element & tildeEa, ff_element & tildeEb,
		     ff_element & P1, const ff_element & A_tau,
		     const ff_element & j_ltau, const ff_element & Ea,
		     const ff_element & Eb);

	int tildeEf (ff_element & tildeEa, ff_element & tildeEb, ff_element & P1,
		     const ff_element & A_tau, const ff_element & Ea,
		     const ff_element & Eb);

	//******* File elkies/compute_psi.cc ***************************

public:

	void compute_psi (ff_pol &res, lidia_size_t k, const ff_polmod &f);
	void compute_psi (ff_pol &res, lidia_size_t k);
	void build_plan  (lidia_size_t ** &to_use, lidia_size_t k);
	void print_plan  (lidia_size_t **to_use, lidia_size_t k);


	//********* File elkies/find_eigenvalue.cc ***************************

	bool find_eigenvalue (lidia_size_t & alpha, const ff_pol & fC,
			      lidia_size_t prev_ev, lidia_size_t lpower,
			      bool test_all);

	bool eigenv_mod3 (lidia_size_t & alpha, const ff_pol  & fC);


	//****** File elkies/schoofpart.cc *****************************

	bool schoofpart (lidia_size_t & eigenval, const ff_pol &F,
			 lidia_size_t *klist, bool test_all);
	bool schoofpart_rat_function (lidia_size_t & eigenval, const ff_pol &F,
				      lidia_size_t *klist, bool test_all);

	//****** File elkies/schoof_algorithm.cc ***********************

	bool schoof_algorithm();

	//********** File eco/eco_prime.cc ****************************



	void Ytop_f (ff_pol &res, const ff_polmod &f);
	void CurveEqn (ff_pol & pol, const bigint & p, const ff_element & a,
		       const ff_element & b, const ff_polmod & f);

	int extract_sqrt(point< gf_element > * &, const point< gf_element > &);
	int max_two_power(const point< gf_element > &);

	//********* File elkies/compute_psi.cc ********************

	void init_psi_comp (PsiPowers* & psi_pow, lidia_size_t & top,
			    lidia_size_t nmax, const ff_polmod &  f);

	void next_psi (PsiPowers* & psi_pow, lidia_size_t & top,
		       const ff_polmod & f, const ff_pol & sqrYY,
		       const ff_element& inv2);
							
	void free_psi (PsiPowers * &psi_pow, lidia_size_t b, lidia_size_t e);

	//********* File elkies/compute_list.cc *************************

	void compute_d_list (lidia_size_t * & dlist);
	void compute_ev_list(lidia_size_t * dlist, lidia_size_t * & ev_list);
	void refine_ev_list(const ff_pol & fC, lidia_size_t * & klist);

	//********* File elkies/BG_algorithms.cc *************************

	void mult_matrix(base_vector< bigint > &, const ff_pol &,
			 const ff_polmod &);

	int search_in_table(const ff_rat & f, ff_rat * table, const ff_polmod & F,
			    int size);

	bool schoofpart_BG (lidia_size_t & ev, const ff_pol & fC,
			    lidia_size_t* klist);

	void double_point_x_only(ff_rat &, const ff_rat &, const ff_polmod &);

	void find_number_of_FBG_steps(lidia_size_t &, lidia_size_t &,
				      double &, lidia_size_t *);
	bool schoofpart_FBG (lidia_size_t & ev, const ff_pol & fC,
			     lidia_size_t* klist);

	//********* File elkies/compute_sign.cc *************************

	void sign_dewaghe(lidia_size_t &, const ff_pol &);


	//******************************************************


	// ---------------------------------------------------------
	// the following public functions are implemented in
	// eco/eco_prime.cc

public :

	//
	// constructors /destructor
	//

	eco_prime ();
	eco_prime (const gf_element & a4, const gf_element & a6);
	~eco_prime ();


	//
	// assignment and access
	//

	eco_prime & operator = (const eco_prime &);

	bool set_prime (udigit prime);

	void set_curve(const gf_element & a, const gf_element & b);

	void reset();
	void set_strategy (char m);
	void set_ev_search_strategy(char m);

	void set_schoof_bound(unsigned int  m)
	{
		schoof_bound = m;
	}

	void set_atkin_bound(unsigned int  m)
	{
		atkin_bound = m;
	}

        void set_BG_bound(unsigned int  m = 40)
        {  
          lower_bound_for_BG = m; 
        }

        void set_FBG_bound(unsigned int  m = 40)
        {  
          lower_bound_for_FBG = m; 
        }

	void use_correctness_test (bool y)
	{
		use_corr_test = y;
	}

	void set_info_mode(int i = 0)
	{
		info = i;
	}

	const base_vector< udigit > & get_relation () const
	{
		return c;
	}

	long get_prime() const
	{
		return static_cast<long>(l);
	}

	//
	// some utility routines
	//

	static void compute_jinv  (ff_element & jinv, const ff_element &a4,
				   const ff_element & a6);

	void compute_jinv(ff_element & jinv)
	{
		eco_prime::compute_jinv(jinv, A, B);
	}

	//
	// splitting type and trace of Frobenius
	//

	void compute_splitting_type();
	bool is_elkies();

private:
	// introduced this new function for unknown sp_type
	void compute_c_for_d(lidia_size_t);

public:
	void compute_trace_atkin();
	void compute_trace_elkies();

	void compute_mod_2_power();

	//
	// group order
	//

	bigint compute_group_order(int info = 0);


	//
	// special cases
	//
	bool is_supersingular(bigint &);
	bool check_j_0_1728(bigint &);

	//
	// verification
	//

	bool probabilistic_correctness_proof(const bigint & res,
					     lidia_size_t tests = 5);

};



bigint compute_group_order(const gf_element & a4, const gf_element & a6, int info = 0);



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_ECO_PRIME_H_GUARD_
