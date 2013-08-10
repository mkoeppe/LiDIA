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
//  File    : gec_E_and_twist_prime.h
//  Author  : Harald Baier (HB)
//	Changes	: See CVS log
//
//==============================================================================================

#ifndef LIDIA_GEC_E_AND_TWIST_PRIME_H_GUARD_
#define LIDIA_GEC_E_AND_TWIST_PRIME_H_GUARD_


#ifndef LIDIA_GEC_COMPLEX_MULTIPLICATION_H_GUARD_
# include "LiDIA/gec_complex_multiplication.h"
#endif


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif


class 
gec_E_and_twist_prime:
public gec_complex_multiplication
{
	/*************************************************
	 *
	 * Class attributes 
	 *
	 **************************************************/

	//
	// A default discriminant respecting the BSI requirements
	// 
private: static bigint default_delta_E_and_twist_prime;
private: static short default_delta_case_E_and_twist_prime;

private: bigint r_tw; // the cryptographic prime factor of twist
private: bigint k_tw; // the cofactor: #E_tw = r_tw * k_tw
private: point< gf_element > G_tw; // the base point of order r_tw 
private: gf_element a2_tw, a4_tw, a6_tw;

	/*************************************************
	 *
	 * Class methods 
	 *
	 **************************************************/

	//
	// constructor and destructor
	//
public: gec_E_and_twist_prime();
public:	gec_E_and_twist_prime( std::ostream & out );
public: ~gec_E_and_twist_prime();

public: const bigint & get_r_tw() const;
public: const bigint & get_k_tw() const;
public: const point<gf_element> & get_G_tw () const;

	//
	// Accessors for the twist parameters
	//
public: const gf_element & get_a2_tw() const;
public: const gf_element & get_a4_tw() const;
public: const gf_element & get_a6_tw() const;

	//
	// High level functions
	//
private: bool assign_class_invariant_to_curve_and_twist( const gf_element & );
private: void assign_efficient_curve_parameters_and_twist( const Fp_polynomial & );
public: void generate();
public: bool find_good_prime();
private: bool is_cryptographically_strong();
};

#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif

#endif	// LIDIA_GEC_E_AND_TWIST_PRIME_H_GUARD_
