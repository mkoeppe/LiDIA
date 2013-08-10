// -*- C++ -*-
//==============================================================================================
//
//      This file is part of LiDIA --- a library for computational number theory
//
//      Copyright (c) 1994--2003 the LiDIA Group.  All rights reserved.
//
//      See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//----------------------------------------------------------------------------------------------
//
//
//  File    : prime_proof.h
//  Author  : Jochen Hechler (JH)
//      Changes : See CVS log
//
//==============================================================================================



#ifndef LIDIA_PRIME_PROOF_H_GUARD_
#define LIDIA_PRIME_PROOF_H_GUARD_

#include        "LiDIA/bigint.h"
#include        "LiDIA/base_vector.h"
#include        "LiDIA/bigfloat.h"
#include        "LiDIA/LiDIA.h"
#include        "LiDIA/quadratic_form.h"
#include        "LiDIA/quadratic_order.h"
#include        "LiDIA/bigcomplex.h"
#include        "LiDIA/polynomial.h"
#include        "LiDIA/gf_polynomial.h"
#include	"LiDIA/certificate.h"
#include 	"LiDIA/elliptic_curve.h"
#include	"LiDIA/point.h"


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif


class prime_proof
{

//
// INTERNAL VARIABLES
//

private:


        bigint	                        n; // the number which will be tested
	bigint				divisior_of_n; // if a divisor is found it will be stored in 
        bool                            verbose; // if true it will display some useful information
	bool				top_level; // if true, this test is the first in the downstep
	certificate			down_cert; // holds the certificate
	int				ecpp_order_mode; // the mode for finding an suitable order in the ecpp algorithm 
	base_vector <bigint>            r_i; // stores the used exponent for the certificate for spp test
        bigint                          p; // a divisor of n-1 or n+1 depends on the test which is used
        bigint                          q; // n-1 = p*q or n+1 = p*q depends on the test which is used
        bigint                          a; // a base which is used in the n_minus_one part
        bigint                          npo_r; // special variables for n_plus_one test. stored for certificate
        bigint                          npo_s; // special variables for n_plus_one test. stored for certificate
	base_vector < bigint >          cert_vector; // vector with certification information
        bigint                          s; // lower bound of divisor for the order 
        bigint                          D; // discriminante used in the ecpp algorithm
        bigint                          start_list; // start prime for the primelist 
        base_vector <int>               prime_list; // the primelist for the ecpp algorithm
	int 				pl_length; // primelist length
        bigint                          order; // order of the elliptic curve
        bigint                          prime_factor; // prime factor of the order
        bigint                          gcd; // if n is not prime, gcd is a factor of n, found bz the ecpp algorihtm
        short                           delta_case;
        long                            class_number;
        base_vector<quadratic_form>     class_group;
        static const bigint             default_delta;
        static const bigint             default_delta_1;
        static const bigint             default_delta_2;
        static const bigint             default_delta_4;
        static const short              default_delta_case;
        static const long               bitlength_no_of_intervals;
        static                          const bigfloat alpha;
        unsigned int                    generation_mode;
        long                            complex_precision;
        polynomial < bigint >           class_polynomial;
        lidia_size_t                    word_length; // the word length of the processor
        static                          const lidia_size_t default_word_length; // its default value
        static                          const lidia_size_t default_degree_oef;
        bool                            mode_one_tried;
	int				cl_min; // minimal classnumber for discriminants, should be one
	int				cl_max; // maximal classnumber for discriminants, sould be 50
	int				pollard_roh_start; // bitlength for pollard rho factor search
	int				pollard_roh_rounds; // number of pollard rho rounds used to find a factor 
	protected: galois_field F_q; // the finite field F_q
//	protected: point<gf_element> G; // the base point of order r
	
	

	
 
 

	
	
	
//
// PUBLIC FUNCTIONS
//

public:
        prime_proof();
	prime_proof (const bigint & p);
        ~prime_proof();
//Changed by G.A - no qualificator prime_proof::
	bool		prove_prime(const bigint & new_p);
	bool		prove_prime();
	bool 		verify_certificate(certificate & cert);
	certificate 	get_certificate();
	void 		set_verbose(bool verbose);
	void		set_classnumbers(const int min, const int max);
	void		set_pollard_rho_parameter(const int start_length,const int rounds);
	void		set_prime(const bigint & p);
	void		set_top(bool is_top);
	void		set_ecpp_mode(int mode);
	bool		spp(const bigint & spp_n);
	bool		spp_verify_proof_prime(const base_vector <bigint> cert_vector);
	bool		n_minus_one( const bigint & new_q, const bigint & new_p);
	bool		nmo_verify_proof_prime(base_vector <bigint> cert_vector);
	bool		n_plus_one( const bigint & new_q, const bigint & new_p);
	bool		npo_verify_proof_prime(base_vector <bigint> cert_vector);
	bool		ecpp();
	bool		ecpp_verify_proof(base_vector <bigint> cert_vector);

	
	

private:
	bool		factorize_n_pm_one();
	bool		minus_one_test(bigint & biggest_prime_factor);
	bool 		plus_one_test(bigint & biggest_prime_factor);
	bool		ecpp_test();
	certificate	spp_get_certificate();
	certificate	nmo_get_certificate();
	certificate	npo_get_certificate();
	certificate	ecpp_get_certificate();
	bigmod		lucas_sequence(const bigint & r, const bigint & s, const bigint & d);
	bigint		found_divisor();
	bool		prove_n_with_curve ();
	void		calculate_prime_list();
	void		calculate_lower_bound_s();
	bool		prove_downstep(const bigint & r);
	bool		assign_class_invariant_to_curve( const gf_element & invariant );
	bool		assign_class_invariant_to_curve_case_minus_three();
	bool		assign_class_invariant_to_curve_case_minus_four();
	bigcomplex	class_invariant( const quadratic_form & Q );
	void		product_real( base_vector<bigfloat> & c,const base_vector<bigfloat> & a,const base_vector<bigfloat> & b );
	void		product_karatsuba( base_vector<bigfloat> & result, const base_vector<bigfloat> & f,const base_vector<bigfloat> & g);
	void		compute_class_polynomial();
	void		compute_class_polynomial_real();
	void		compute_class_polynomial_karatsuba_high_precision();
	void		compute_class_polynomial_karatsuba_low_precision();
	bool		find_next_discriminant(long start_h, long end_h);
	bool		is_good_order(bigint & t, bigint & y);
	bigint		pollard_rho(bigint & t);
	void		set_complex_precision();
	void		set_generation_mode();
	bool 		generate_curve();

	

	
	
};

#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif


#endif
