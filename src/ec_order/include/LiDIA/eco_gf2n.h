// -*- C++ -*-
//==============================================================================================
//
//      This file is part of LiDIA --- a library for computational number theory
//
//      Copyright (c) 1994--2001 the LiDIA Group.  All rights reserved.
//
//      See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//----------------------------------------------------------------------------------------------
//
//      $Id$
//
//      Author  :
//      Changes : See CVS log
//
//==============================================================================================


#ifndef LIDIA_ECO_GF2N_H_GUARD_
#define LIDIA_ECO_GF2N_H_GUARD_



#ifndef LIDIA_GF2N_POLY_MODULUS_H_GUARD_
# include	"LiDIA/gf2n_poly_modulus.h"
#endif
#ifndef LIDIA_GF2N_RATIONAL_FUNCTION_H_GUARD_
# include	"LiDIA/gf2n_rational_function.h"
#endif
#ifndef LIDIA_POWER_TABLE_H_GUARD_
# include	"LiDIA/power_table.h"
#endif
#ifndef LIDIA_MV_POLY_H_GUARD_
# include	"LiDIA/mv_poly.h"
#endif
#ifndef LIDIA_UDIGIT_H_GUARD_
# include	"LiDIA/udigit.h"
#endif
#ifndef LIDIA_FF1_H_GUARD_
# include	"LiDIA/ff1.h"
#endif
#ifndef LIDIA_FF2_H_GUARD_
# include	"LiDIA/ff2.h"
#endif
#ifndef LIDIA_UDIGIT_MOD_H_GUARD_
# include	"LiDIA/udigit_mod.h"
#endif
#ifndef LIDIA_BASE_VECTOR_H_GUARD_
# include	"LiDIA/base_vector.h"
#endif
#ifndef LIDIA_ELLIPTIC_CURVE_H_GUARD_
# include	"LiDIA/elliptic_curve.h"
#endif
#ifndef LIDIA_POINT_H_GUARD_
# include	"LiDIA/point.h"
#endif
#ifndef LIDIA_GF_ELEMENT_H_GUARD_
# include	"LiDIA/gf_element.h"
#endif

#ifdef TIMING			// to allow precise timings of all different
# include	"LiDIA/timer.h"	// parts of the program
#endif


#define MEQ2_PREFIX "meq_2_"



#ifdef LIDIA_NAMESPACE
namespace LiDIA
{
# define IN_NAMESPACE_LIDIA
#endif



  class eco_gf2n
  {
    private:
      // ***** types used in the class eco_gf2n *****
    typedef gf2n ff_element;
    typedef gf2n_polynomial ff_pol;
    typedef gf2n_rational_function ff_rat;
    typedef gf2n_poly_modulus ff_polmod;
    typedef ff1 ff1_element;
    typedef ff2 ff2_element;

      // the forward and friend declarations are required since 
      // PsiPowers needs access to the private type spcifier ff_pol.
    struct PsiPowers;
    friend struct PsiPowers;

    struct PsiPowers
    {
      ff_pol *pow1;
      ff_pol *pow2;
      ff_pol *pow3;
    };



    // ***** constants *****

    static const char UNKNOWN;	// values for sp_type
    static const char ALL_SUBGRPS_INV;
    static const char ONE_SUBGRP_INV;
    static const char TWO_SUBGRPS_INV;
    static const char NO_SUBGRP_INV;

  public:

    static const char COMPUTE_SP_DEGREE;	// values for degree_mode
    static const char DONT_COMPUTE_SP_DEGREE;
    static const char COMPUTE_SP_DEGREE_IF_ATKIN;
    static const char COMPUTE_SP_DEGREE_IF_ELKIES;
    static const char COMPUTE_SP_DEGREE_AND_ROOTS;

    static const char EV_RATIONAL_FUNCTION;	// values for ev_strategy
    static const char EV_DIVISION_POLYNOMIAL;
    static const char EV_RATIONAL_FUNCTION_TABLE;
    static const char EV_BABYSTEP_GIANTSTEP;
    static const char EV_FUNNY_BABYSTEP_GIANTSTEP;
    static const char EV_FUNNY_BABYSTEP_GIANTSTEP_SIGN;

    static const char PLAIN_LERCIER;	// use Lerciers idea or the Velu extension
    static const char WITH_VELU;

    // in the highest bit, we store what coordinate is tested --> for 0, use
    // y-coordinate (default), for 1 use X-coordinate only.

    static const char TEST_X_COORDINATE;	// test with X or Y-coordinate
    static const char TEST_Y_COORDINATE;	// in Elkies part

  private:

    // ***** class variables

      bigint pn;		// number of elements of finite field
    ff_element A6;		// elliptic curve Y^2 + X*Y = X^3 + A6
    ff_element jinv;		// j - invariant of E
    ff_pol meq_pol;		// l - th modular polynomial
    ff_element ftau[2];		// roots of meq_pol
    char sp_type;		// splitting type of meq_pol
    char degree_mode;
    udigit sp_degree;		// contains the splitting degree of meq_pol
    bool sp_degree_known;	// if sp_degree_known is true
    char ev_strategy;		// algorithm to find eigenvalues.
    unsigned int schoof_bound;	// for smaller Atkin primes, use Schoof
    unsigned int atkin_bound;	//  no splt.type computed for l >...
    unsigned int power_bound;	// degree bound for using powers of l
    unsigned int lower_bound_for_BG;	// don't use BG, FBG algorithms for
    unsigned int lower_bound_for_FBG;	// small l
    unsigned int degree_max_2_power;	// degree of polynomial for 2^i TP
    bool use_corr_test;		// use prob. correctness test iff true

    int info;
    udigit q;			// pn mod l, l might be a prime power
      sort_vector < udigit > c;	// contains c mod l, c is trace of Frobenius
    udigit l;


    bool test_x_coordinate ()
    {
      if (((ev_strategy == EV_FUNNY_BABYSTEP_GIANTSTEP ||
	    ev_strategy == EV_FUNNY_BABYSTEP_GIANTSTEP_SIGN) &&
	   l > lower_bound_for_FBG) || (ev_strategy & 0x80))
	return true;
      else
	return false;
    }

    bool test_x_coordinate_uniquely ()
    {
      if (is_prime (static_cast < udigit > (l)) &&
	  ev_strategy == EV_FUNNY_BABYSTEP_GIANTSTEP_SIGN &&
	  l >= lower_bound_for_FBG && sp_type != ONE_SUBGRP_INV
	  && sp_type != ALL_SUBGRPS_INV)
	return true;

      if (!is_prime (static_cast < udigit > (l)) && c.size () == 1)
	return true;

      return false;
    }

    bool test_y_coordinate ()
    {
      return !test_x_coordinate ();
    }


    //********* File lercier/lercier_functions.cc ****************

    void compute_sqrtP (ff_element * &p_k,
			const power_tab < ff_element > &alpha_pow,
			const ff_element & sqrt4_alpha,
			const ff_element & sqrt4_beta,
			const ff_element & beta, lidia_size_t d,
			char strat = eco_gf2n::WITH_VELU);

    void eqnsys2n_1 (mv_poly * &p_mv, lidia_size_t k,
		     const power_tab < ff_element > &ap, lidia_size_t d);

    void eqnsys2n_2 (mv_poly * &p_mv, lidia_size_t k,
		     const power_tab < ff_element > &ap, lidia_size_t d);

    int n_over_k (lidia_size_t n, lidia_size_t k);

    void build_eqn1 (mv_poly & EQ, mv_poly * &p_mv, lidia_size_t k,
		     const power_tab < ff_element > &ap, lidia_size_t d);

    void build_eqn2 (mv_poly & EQ, mv_poly * &p_mv, lidia_size_t k,
		     const power_tab < ff_element > &ap, lidia_size_t d);

    bool velu_improvement (mv_poly * &p_mv, lidia_size_t K, lidia_size_t d,
			   const power_tab < ff_element > &ap,
			   const ff_element & sqrt_alpha,
			   const ff_element & sqrt_beta);

    bigint solve_eqs (mv_poly * &EQS, lidia_size_t number_var,
		      lidia_size_t number_eq, lidia_size_t bit_mask,
		      lidia_size_t number_of_guessed_var);

    void compute_sqrtP_guess (int current_K,
			      int d, mv_poly * &p_mv,
			      gf2n * &p_k,
			      const power_tab < gf2n > &alpha_pow,
			      const gf2n & sqrt4_alpha,
			      const gf2n & sqrt4_beta, char strategy);


    //******* File elkies/compute_psi_eco2.cc *****************************

    void init_psi_comp (PsiPowers * &, lidia_size_t &, lidia_size_t,
			ff_polmod &);
    void compute_psi (ff_pol & res, lidia_size_t k, ff_polmod &);
    void compute_psi (ff_pol & res, lidia_size_t k);
    void build_plan (lidia_size_t ** &to_use, lidia_size_t k);
    void print_plan (lidia_size_t ** to_use, lidia_size_t k);
    void next_psi (PsiPowers * &, lidia_size_t &, ff_polmod &);
    void free_psi (PsiPowers * &, lidia_size_t, lidia_size_t);


    //********* File elkies/find_eigenvalue_eco2.cc ************************


    bool find_eigenvalue (lidia_size_t & alpha, const ff_pol & fC,
			  lidia_size_t prev_ev, lidia_size_t lpower,
			  bool test_all);

    bool eigenv_mod3 (lidia_size_t & alpha, const ff_pol & fC);

    bool eigenv_rat_oneposs (lidia_size_t & eigenval, const ff_pol & fC,
			     lidia_size_t ev);


    //****** File elkies/schoofpart_eco2.cc *******************************


    bool schoofpart (lidia_size_t & eigenval, const ff_pol & F,
		     lidia_size_t * klist, bool test_all);
    bool schoofpart_rat_function (lidia_size_t & eigenval, const ff_pol & F,
				  lidia_size_t * klist, bool test_all);


    //****** File elkies/schoof_algorithm_eco2.cc **********************

    bool schoof_algorithm (int split_degree = -1);

    //********* File elkies/compute_list_eco2.cc *********************


    void compute_d_list (lidia_size_t * &dlist);
    void compute_ev_list (lidia_size_t * dlist, lidia_size_t * &ev_list);

    //********* File elkies/Ytop_eco2.cc ******************************

    void square_y_term (ff_pol & res0, ff_pol & res1,
			const ff_pol & curve, ff_polmod & f);
    void power_block (ff_pol & g, ff_pol * table, int r);
    void Ytop_f (ff_pol & y0, ff_pol & y1, ff_polmod &);
    void Ytop_f_plain (ff_pol & y0, ff_pol & y1, ff_polmod &);

    //********** File elkies/two_power_eco2.cc ***********************

    void compose_2_torsion_pol (ff_pol &, const ff_pol &, const ff_pol &);

    int sqrt_of_point (point < gf_element > &res,
		       const point < gf_element > &);

    //*********** eco/eco_gf2n.cc *************************************

    void read_meq (ff_pol &, int l);

    void compute_c_for_d (lidia_size_t d);

    //********* File elkies/BG_algorithms_eco2.cc *************************

    void mult_matrix (base_vector < ff_element > &, const ff_pol &,
		      ff_polmod &);

    void inner_product (ff_element & x, const base_vector < ff_element > &a,
			const ff_pol & b);

    int search_in_table (const ff_rat & f, ff_rat * table, ff_polmod & F,
			 int size);

    bool schoofpart_BG (lidia_size_t & ev, const ff_pol & fC,
			lidia_size_t * klist);

    void double_point_x_only (ff_rat &, const ff_rat &, ff_polmod &);
    void add_different_point_x_only (ff_rat &, const ff_rat &,
				     const ff_rat &, const ff_rat &,
				     ff_polmod &);

    void find_number_of_FBG_steps (lidia_size_t &, lidia_size_t &,
				   double &, lidia_size_t *);

    bool schoofpart_FBG (lidia_size_t & ev, const ff_pol & fC,
			 lidia_size_t * klist, bool t = false);



  public:

    //
    // constructors /destructor
    //

    eco_gf2n ();
    eco_gf2n (const ff_element & a6);
    ~eco_gf2n ();


    //
    // assignment and access
    //

    eco_gf2n & operator = (const eco_gf2n &);

    void reset ();
    bool set_prime (udigit prime);
    void set_curve (const ff_element & a6);
    void set_strategy (char m);
    void set_ev_search_strategy (char m);

    void set_schoof_bound (unsigned int m = 15)
    {
      schoof_bound = m;
    }

    void set_atkin_bound (unsigned int m = 20)
    {
      atkin_bound = m;
    }

    void set_power_bound (unsigned int m = 100)
    {
      power_bound = m;
    }

    void set_BG_bound (unsigned int m = 40)
    {
      lower_bound_for_BG = m;
    }

    void set_FBG_bound (unsigned int m = 40)
    {
      lower_bound_for_FBG = m;
    }


    void set_info_mode (int i = 1)
    {
      info = i;
    }

    void use_correctness_test (bool y = true)
    {
      use_corr_test = y;
    }

    const base_vector < udigit > &get_relation () const
    {
      return c;
    }

    long get_prime () const
    {
      return static_cast < long >(l);
    }

    //
    // some utility routines
    //

    static void compute_jinv (ff_element & jinv, const ff_element & a6)
    {
      if (!a6.is_zero ())
	invert (jinv, a6);
      else
	lidia_error_handler ("compute_jinv", "elliptic curve is singular");
    }

    void compute_jinv (ff_element & jinv)
    {
      eco_gf2n::compute_jinv (jinv, A6);
    }


    //
    // splitting type and trace of Frobenius
    //

    void compute_splitting_type ();

    bool is_elkies ();

    void compute_trace_atkin ();

    void compute_trace_elkies ();

    void compute_powers_of_l (const ff_pol & fC, lidia_size_t ev,
			      lidia_size_t ev1_order,
			      lidia_size_t ev2_order,
			      lidia_size_t & ev_lpower,
			      lidia_size_t & lpower);

    void compute_sign (lidia_size_t &, const ff_pol &);

    void compute_mod_2_power ();

    bigint compute_group_order ();


    //
    // verification
    //

    bool probabilistic_correctness_proof (const bigint & res,
					  lidia_size_t tests = 5);

    //******* lercier/lercier_isogeny.cc ***************************

    void divisor_of_divpol_via_lercier (ff_pol & Q,
					char strategy = eco_gf2n::WITH_VELU);

    void isogeny_via_lercier (ff_pol & numerator, ff_pol & denominator,
			      const ff_element & A6_1,
			      const ff_element & A6_2, char strat =
			      eco_gf2n::WITH_VELU);

    void prime_square (ff_pol & divpol_divisor, const ff_pol & fC,
		       const ff_element & A6_2);


  public:static gf_element convertToGFElement (const gf2n & a);
  public:static gf2n convertToGF2n (const gf_element & a);
  public:static elliptic_curve < gf_element >
      generateCurve (const gf2n & a6);
  public:static gf2n j_invariant (const gf2n & a1, const gf2n & a2,
			     const gf2n & a3, const gf2n & a4,
			     const gf2n & a6);
  };



  bigint compute_group_order (const gf2n & a6, int info = 0);
  bigint compute_group_order (const gf2n & a2, const gf2n & a6, int info = 0);
  bigint compute_group_order (const gf2n & a1, const gf2n & a2,
			      const gf2n & a3, const gf2n & a4,
			      const gf2n & a6, int info = 0);



#ifdef LIDIA_NAMESPACE
}				// end of namespace LiDIA

# undef IN_NAMESPACE_LIDIA
#endif



#endif // LIDIA_ECO_GF2H_H_GUARD_
