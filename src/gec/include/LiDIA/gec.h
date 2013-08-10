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
//  File    : gec.h
//  Author  : Harald Baier (HB)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_GEC_H_GUARD_
#define LIDIA_GEC_H_GUARD_

#ifndef LIDIA_BASE_VECTOR_H_GUARD_
# include "LiDIA/base_vector.h"
#endif
#ifndef LIDIA_BIGFLOAT_H_GUARD_
# include "LiDIA/bigfloat.h"
#endif
#ifndef LIDIA_BIGINT_H_GUARD_
# include "LiDIA/bigint.h"
#endif
#ifndef LIDIA_ELLIPTIC_CURVE_H_GUARD_
# include "LiDIA/elliptic_curve.h"
#endif
#ifndef LIDIA_GF_ELEMENT_H_GUARD_
# include "LiDIA/gf_element.h"
#endif
#ifndef LIDIA_POINT_H_GUARD_
# include "LiDIA/point.h"
#endif
#ifndef LIDIA_QUADRATIC_FORM_H_GUARD_
# include "LiDIA/quadratic_form.h"
#endif
#ifndef LIDIA_QUADRATIC_ORDER_H_GUARD_
# include "LiDIA/quadratic_order.h"
#endif
#ifndef LIDIA_RATIONAL_FACTORIZATION_H_GUARD_
# include "LiDIA/rational_factorization.h"
#endif
#ifdef TIMING                    // to allow precise timings of all different
# include "LiDIA/timer.h"         // parts of the program
#endif

#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif


/*************************************************

gec is a base class for the
classes which provide the generation methods.
This class declares the common attributes and
methods.

*************************************************/
   
class gec
{
	/*************************************************
	 *
	 * The attributes of this class are named according
	 * to the P1363 standard of IEEE.
	 *
	 **************************************************/

protected: bigint q; // the cardinality of the finite field F_q
protected: bigint p; // the characteristic of the finite field F_q
protected: lidia_size_t degree; // the degree of F_q over F_p
protected: galois_field F_q; // the finite field F_q
protected: bigint r; // the cryptographic prime factor
protected: bigint k; // the cofactor: #E = r * k
protected: point<gf_element> G; // the base point of order r 

	//
	// Upper bound of k and lower bound of bitlength of r 
	// and their static default values.
	// The default values are set in gec.cc and may be changed there.
	// 
protected: lidia_size_t lower_bound_bitlength_r;
protected: static lidia_size_t default_lower_bound_bitlength_r;
protected: bigint upper_bound_k;
protected: static bigint default_upper_bound_k;

	//
	// The coefficients of the elliptic curve: We have
	// 
	// y^2      = x^3 +          a4*x + a6  if char(F_q)!=2,3
	// y^2 + xy = x^3 + a2*x^2 +        a6  if char(F_q) = 2,
	// y^2      = x^3 + a2*x^2 + a4*x + a6  if char(F_q) = 3.
	//
protected: gf_element a2, a4, a6; 
	
	//
	// the boolean is_initialized is equal to true if and only if
	// all attributes are initialized. This can be achieved
	// by either invoking the method generate() or 
	// by using the appropriate constructor.
	//
protected: bool is_initialized;

	//
	// Static variable to store the upper bound of n 
	// in case of GF(2^n): For larger values point counting is not possible.
	//
public: static lidia_size_t upper_bound_of_degree;

	//
	// Static variable to define number of Miller-Rabin-Tests
	//
public: static int nr_of_prob_prime_tests;

	//
	// the boolean VERBOSE gives a lot information on the status
	// of the program if it is set to true.
	//
protected: bool VERBOSE;

	//
	// the boolean according_to_BSI is equal to true if and only if
	// the generated curve respects the BSI requirements
	// 
protected: bool according_to_BSI;
	// lower bound of the class number of maximal order containing End(E)
protected: static lidia_size_t BSI_lower_bound_h_field;
	// lower bound of extension degree to avoid MOV/FR-attacks
protected: static lidia_size_t BSI_lower_bound_extension_degree;

	//
	// Default value of the bitlength of fields in which
	// MOV/FR-attacks are considered to be intractable.
	// This value differs from the BSI requirements and is only
	// of relevance if according_to_BSI is equal to false.
	//
protected: static lidia_size_t default_lower_bound_extension_bitlength;

	//
	// the discriminant and the lower bound of the class number of the field 
	// containing the endomorphism ring of the curve.
	// * delta_field is only computed iff according_to_BSI == true.
	// * h_field is only computed iff according_to_BSI == true and the CM-Flag
	//   is set
	// * lower_bound_h_field is only computed iff according_to_BSI == true 
	//   and the PC-Flag is set
protected: bigint delta_field;
protected: lidia_size_t h_field;
protected: lidia_size_t lower_bound_h_field;

	//
	// discriminant delta, its absolute value, and class number h of End(E)
	//
protected: bigint delta;
protected: bigint abs_delta;
protected: lidia_size_t h;

	//
	// flags to handle different cases in the inheriting classes
	//
protected: bool CM_FLAG; // gec is used as base class for CM-method
protected: bool PC_FLAG; // gec is used as base class for random approach

	// if it is set to true curves with a = -3 are generated
protected: bool efficient_curve_parameters;

	//
	// file handling: output files for timings and parameters
	//
public: char output_timings[ 100 ];
public: char output_parameters[ 100 ];
public: bool timings;
public: timer measure_time;
public: std::ostream & out_;

	/*************************************************
	 *
	 * Class methods 
	 *
	 **************************************************/
	
	//
	// constructor and destructor
	//
public: gec();
public: gec( std::ostream & out ); // constructor to write data to out
public: virtual ~gec();

	//
	// Accessors for primes, cofactor, base point:
	// The methods return 0 if the corresponding attribute is not set.
	//
public: const bigint & get_p() const;
public: const bigint & get_q() const;
public: lidia_size_t get_degree() const;
public: const bigint & get_r() const;
public: const bigint & get_k() const;
public: const point<gf_element> & get_G () const;
	
	//
	// Accessors for the curve parameters:
	// The lidia_error_handler is invoked if curve is not initialized.
	//
public: const gf_element & get_a2() const;
public: const gf_element & get_a4() const;
public: const gf_element & get_a6() const;

	//
	// Accessors for default values of bounds of bitlength_r and upper_bound_k
	//
public: lidia_size_t get_default_lower_bound_bitlength_r() const;
public: const bigint & get_default_upper_bound_k() const;

	//
	// Accessors for bounds of bitlength_r and upper_bound_k:
	// If they are not set, the method returns 0.
	//
public: lidia_size_t get_lower_bound_bitlength_r() const;
public: const bigint & get_upper_bound_k() const;

	//
	// Accessors for BSI settings
	//
public: bool get_according_to_BSI() const;
public: lidia_size_t get_BSI_lower_bound_h_field() const;
public: lidia_size_t get_BSI_lower_bound_extension_degree() const;

	//
	// Accessor for default_lower_bound_extension_bitlength
	//
public: lidia_size_t get_default_lower_bound_extension_bitlength() const;

	//
	// Accessors for discriminant, class number and lower bound of 
	// class number of Quot(End(E))
	// They return 0 if the corresponding attribute is not set.
	//
public: const bigint & get_delta_field() const;
public: lidia_size_t get_h_field() const;
public: lidia_size_t get_lower_bound_h_field() const;

	//
	// Accessors for discriminant and class number of End(E).
	// If they are not computed (e.g. because of PC-method), 0 is returned.
	//
public: const bigint & get_delta() const;
public: lidia_size_t get_h() const;

	//
	// Accessors of verbosity and efficient curve parameters
	//
public: bool get_verbose_level();
public: bool get_efficient_curve_parameters();



	//
	// Mutators of field and degree
	//
public: void set_field( const bigint & ); // to initialize with prime field
public: void set_degree( lidia_size_t ); // to initialize degree of field

	//
	// Mutators of security level
	//
public: void set_lower_bound_bitlength_r( lidia_size_t );
public: void set_upper_bound_k( const bigint & );

	//
	// Mutators of default security bounds
	//
public: void set_default_lower_bound_bitlength_r( lidia_size_t );
public: void set_default_upper_bound_k( const bigint & );
public: void set_default_lower_bound_extension_bitlength( lidia_size_t );

	//
	// Mutators of BSI-attributes
	//
public: void set_according_to_BSI( bool );
public: void set_BSI_lower_bound_h_field( lidia_size_t );
public: void set_BSI_lower_bound_extension_degree( lidia_size_t );

	//
	// Mutators of verbosity and efficient curve parameters
	//
public: void set_verbose_level( bool );
public: void set_efficient_curve_parameters( bool );

	//
	// High level functions: 
    // * is_cryptographically_strong: gec, gec_point_counting, 
	//   gec_complex_multiplication:
	//      returns true if the curve defined by a2, a4, a6
	//      is cryptographically strong.
	// * is_cryptographically_strong: gec_E_and_twist_prime:
	//      returns true if the curve defined by a2, a4, a6
	//      and its twist are both cryptographically strong.
	//
protected: bool is_cryptographically_strong( const bigint & order );

	//
	// High level functions: 
    // * compute_delta_field: computes the discriminant
	//   of the maximal order containing End(E)
	// * compute_lower_bound_h_field: computes a lower
	//   bound of the class number of the maximal order containing End(E)
	//
protected: bool compute_delta_field();
protected: lidia_size_t compute_lower_bound_h_field();

	//
	// High level functions: 
	// * pure virtual function generate(), which has to 
	//   be implemented by the inheriting class
	// * check_the_settings(): invoked by generate() to control if the
	//   set values are proper
	//
public: virtual void generate() = 0;
protected: virtual void check_the_settings();
};

#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif

#endif	// LIDIA_GEC_H_GUARD_
