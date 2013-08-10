// -*- C++ -*-
//==============================================================================================
//
//	This file is part of LiDIA --- a library for computational number theory
//
//	Copyright (c) 1994--2002 the LiDIA Group.  All rights reserved.
//
//	See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//----------------------------------------------------------------------------------------------
//
//
//  File    : gec_complex_multiplication.h
//  Author  : Harald Baier (HB)
//	Changes	: See CVS log
//
//==============================================================================================

#ifndef LIDIA_GEC_COMPLEX_MULTIPLICATION_H_GUARD_
#define LIDIA_GEC_COMPLEX_MULTIPLICATION_H_GUARD_

#include <string>

#ifndef LIDIA_BIGCOMPLEX_H_GUARD_
# include "LiDIA/bigcomplex.h"
#endif
#ifndef LIDIA_GEC_H_GUARD_
# include "LiDIA/gec.h"
#endif
#ifndef LIDIA_POLYNOMIAL_H_GUARD_
# include "LiDIA/polynomial.h"
#endif
#ifndef LIDIA_GF_POLYNOMIAL_H_GUARD_
# include "LiDIA/gf_polynomial.h"
#endif


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif




class 
gec_complex_multiplication : 
	public gec
{
	/*************************************************
	 *
	 * Class attributes 
	 *
	 **************************************************/

	//
	// delta_case and class group of End(E)
	//
	// Currently we distinguish the cases 
	//       delta = 0 mod 4  <=>    delta_case = 0
	//       delta = 1 mod 4  <=>    delta_case = 1
protected: short delta_case; 
protected: base_vector<quadratic_form> class_group;
	
	//
	// Static variables: default discriminant, its case, ... (see *.cc file)
	// 
protected: static const bigint default_delta;
protected: static const bigint default_delta_1;
protected: static const bigint default_delta_2;
protected: static const bigint default_delta_4;
protected: static const short default_delta_case;
protected: static const long bitlength_no_of_intervals;
protected: static const bigfloat alpha;

	// In case of the Fixed Field Approach or the OEF,
	// default values of interval of h where discriminants are tested
protected: static const long default_lower_bound_class_number;
protected: static const long default_upper_bound_class_number;

	//
	// generation_mode stores the class invariant used for generation
	// Currently we support the following class invariants:
	// generation_mode = 1 <=> j
	// generation_mode = 3 <=> gammma_2
	// generation_mode = 4 <=> weber_f
	//
protected: unsigned int generation_mode;
	// complex_precision stores the floating point precision	
protected: long complex_precision;
	// class_polynomial stores the minimal polynomial of class invariant
protected: polynomial < bigint > class_polynomial;
	// stores whether class polynomial is already computed
protected: bool is_polynomial_set;
	
	// to set lower bound and upper bound of h for Fixed Field Approach
	// and OEF case: Only discriminants in this interval are considered.
protected: long lower_bound_class_number;
protected: long upper_bound_class_number;

	//
	// some attributes for the OEF case
	//
protected: lidia_size_t word_length; // the word length of the processor
protected: static const lidia_size_t default_word_length; // its default value
protected: static const lidia_size_t default_degree_oef;


	/*************************************************
	 *
	 * Class methods 
	 *
	 **************************************************/

	//
	// constructor and destructor
	//
public: gec_complex_multiplication();
public: gec_complex_multiplication( const bigint & );
public: gec_complex_multiplication( lidia_size_t d );
public: gec_complex_multiplication( std::ostream & out );

public: ~gec_complex_multiplication();


	//
	// accessors
	//
public: unsigned int get_generation_mode() const;
public: const bigint & get_delta() const;
public: long get_complex_precision() const;
public: polynomial < bigint > get_class_polynomial() const;

	//
	// mutators
	//
public: void set_generation_mode( unsigned int );
public: void set_delta( const bigint & );
public: void set_complex_precision( long );
public: void set_lower_bound_class_number( const long );
public: void set_upper_bound_class_number( const long );
public: void set_class_polynomial( const polynomial < bigint > & );
public: void set_class_polynomial(const bigint &, const lidia_size_t, char * );

	//
	// assignment
	//
public: void assign( const gec_complex_multiplication & );

 
	//
	// High level functions
	//
public: void generate();
public: void generate_oef( lidia_size_t, lidia_size_t, lidia_size_t );

public: virtual bool find_good_prime();
public: bool find_discriminant();
public: bool find_discriminant_inefficient();
protected: bool find_good_prime_0_mod_4();
protected: bool find_good_prime_1_mod_8();
protected: bool find_good_prime_5_mod_8();
protected: bool find_good_prime_power();

protected: bool find_oef();

protected: bool assign_class_invariant_to_curve( const gf_element & );
protected: void assign_efficient_curve_parameters( const Fp_polynomial & );

public: void compute_class_polynomial();
protected: void compute_class_polynomial_real();
protected: void compute_class_polynomial_karatsuba_low_precision();
protected: void compute_class_polynomial_karatsuba_high_precision();
protected: void product_karatsuba( base_vector<bigfloat> & result,
								   const base_vector<bigfloat> & f,
								   const base_vector<bigfloat> & g);
protected: void product_real( base_vector<bigfloat> & result,
								   const base_vector<bigfloat> & f,
								   const base_vector<bigfloat> & g);


protected: bigcomplex class_invariant( const quadratic_form & );

public: static std::string path_to_discriminant_files();
};


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif

#endif	// LIDIA_GEC_COMPLEX_MULTIPLICATION_H_GUARD_
