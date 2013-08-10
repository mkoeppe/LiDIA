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
//  File    : gec_complex_multiplication.cc
//  Author  : Harald Baier (HB)
//	Changes	: See CVS log
//
//==============================================================================================

#ifdef HAVE_CONFIG_H
# include       "config.h"
#endif
# include        "LiDIA/gec_complex_multiplication.h"
#ifdef VERBOSE
# include        <cassert>
#endif


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif


//
// Initialize static variables
//

// default discriminant (fundamental of class number 200), if no delta
// is set by the user 
const bigint gec_complex_multiplication::default_delta = -21311;
const bigint gec_complex_multiplication::default_delta_4 = -21311;
const bigint gec_complex_multiplication::default_delta_1 = -125579;
const bigint gec_complex_multiplication::default_delta_2 = -53444;
// as default_delta is odd, the corresponding case is set to 1
const short gec_complex_multiplication::default_delta_case = 1;
// the default number of subintervals is 2^40
const long gec_complex_multiplication::bitlength_no_of_intervals = 40;
// the precision for computation of Weber polynomials is alpha * L (see PhD)
const bigfloat gec_complex_multiplication::alpha = 0.015;
// the default lower bound of the class number is 200:
// in case of Fixed Field Approach or OEF, this is the default value
// for the lower bound of the class number
const long gec_complex_multiplication::default_lower_bound_class_number = 200;
// the default upper bound of the class number is 1000:
// in case of Fixed Field Approach or OEF, this is the default value
// for the upper bound of the class number
const long gec_complex_multiplication::default_upper_bound_class_number = 1000;

// some constant static variables for the OEF case
const lidia_size_t gec_complex_multiplication::default_word_length = 32;
const lidia_size_t gec_complex_multiplication::default_degree_oef = 5;


//
// Constructor / Destructor
//
gec_complex_multiplication::
gec_complex_multiplication() :
	gec()
{
	CM_FLAG = true;

	class_polynomial.assign_zero();
	complex_precision = 0;
	degree = 1;
	generation_mode = 0;

	lower_bound_class_number = upper_bound_class_number = 0;
	word_length = 0;

	is_initialized = false;
	is_polynomial_set = false;
	efficient_curve_parameters = false;
}

gec_complex_multiplication::
gec_complex_multiplication( const bigint & Delta ) :
	gec()
{
	CM_FLAG = true;

	delta.assign( Delta );
	class_polynomial.assign_zero();
	complex_precision = 0;
	degree = 1;
	generation_mode = 0;

	lower_bound_class_number = upper_bound_class_number = 0;
	word_length = 0;

	is_initialized = false;
	is_polynomial_set = false;
	efficient_curve_parameters = false;
}


gec_complex_multiplication::
gec_complex_multiplication( lidia_size_t d ) :
	gec()
{
	CM_FLAG = true;
  
	class_polynomial.assign_zero();
	complex_precision = 0;
	degree = d;
	generation_mode = 0;

	lower_bound_class_number = upper_bound_class_number = 0;
	word_length = 0;

	is_initialized = false;
	is_polynomial_set = false;
	efficient_curve_parameters = false;

	lower_bound_class_number = default_lower_bound_class_number;
}

gec_complex_multiplication::
gec_complex_multiplication( std::ostream & out ) : gec( out )
{
	CM_FLAG = true;

	class_polynomial.assign_zero();
	complex_precision = 0;
	degree = 1;
	generation_mode = 0;

	lower_bound_class_number = upper_bound_class_number = 0;
	word_length = 0;

	is_initialized = false;
	is_polynomial_set = false;
	efficient_curve_parameters = false;
}


gec_complex_multiplication::
~gec_complex_multiplication()
{
}

//
// Accessors
//
unsigned int
gec_complex_multiplication::get_generation_mode() const
{
	return( generation_mode );
}

const bigint & 
gec_complex_multiplication::get_delta() const
{
	return( delta );
}

long
gec_complex_multiplication::get_complex_precision() const
{
	return( complex_precision );
}

polynomial < bigint >
gec_complex_multiplication::get_class_polynomial() const
{
	return( class_polynomial );
}

//
// Mutators
//

// Mutator for generation_mode: see file set_generation_mode.cc

// Mutator for discriminant of End(E)
void
gec_complex_multiplication::set_delta( const bigint & d )
{
	if( d.is_ge_zero() )
		lidia_error_handler("gec_complex_multiplication", "set_delta( const bigint & delta): delta is non negative");

	abs_delta.assign( d );
	abs_delta.negate();
	
	if( ( abs_delta.bit(0) == 0 && abs_delta.bit(1) == 0 ) )
	{
		delta_case = 0;
		delta = d;
	}
	else if( ( abs_delta.bit(0) == 1 && abs_delta.bit(1) == 1 ) )
	{
		delta_case = 1;
		delta = d;
	}
	else
		lidia_error_handler("gec_complex_multiplication", "set_delta( const bigint & delta): delta is not = 0,1 mod 4");
}

// Mutator for complex_precision: see file set_complex_precision.cc


// Mutator for lower_bound_class_number
void
gec_complex_multiplication::set_lower_bound_class_number( 
	const long lower_bound )
{
	if( lower_bound < default_lower_bound_class_number ||
		lower_bound > default_upper_bound_class_number )
		lidia_error_handler("gec_complex_multiplication", "set_lower_bound_class_number( const long h_): h_ is not in the correct interval. Check the default values.");
	else if( upper_bound_class_number != 0 &&
			 upper_bound_class_number <= lower_bound )
		lidia_error_handler("gec_complex_multiplication", "set_lower_bound_class_number( const long h_): h_ has to be smaller than upper_bound_class_number.");

	lower_bound_class_number = lower_bound;
}

// Mutator for upper_bound_class_number
void
gec_complex_multiplication::set_upper_bound_class_number( 
	const long upper_bound )
{
	if( upper_bound < default_lower_bound_class_number ||
		upper_bound > default_upper_bound_class_number )
		lidia_error_handler("gec_complex_multiplication", "set_upper_bound_class_number( const long h_): h_ is not in the correct interval. Check the default values.");
	else if( lower_bound_class_number >= upper_bound )
		lidia_error_handler("gec_complex_multiplication", "set_upper_bound_class_number( const long h_): h_ has to be larger than lower_bound_class_number.");

	upper_bound_class_number = upper_bound;
}


// Mutator for class_polynomial: The user responsible for the right choice.
void
gec_complex_multiplication::set_class_polynomial(
	const polynomial < bigint > &  class_pol )
{
	class_polynomial.assign( class_pol );
	is_polynomial_set = true;
}

// Mutator for class_polynomial: The user responsible for the right choice.
// In addition, the function requires the coefficients to be stored
// in the file <file_directory>/h_/delta_ starting with the
// coefficient of X^h_ = 1 in descending order.
void
gec_complex_multiplication::set_class_polynomial(
	const bigint & delta_, const lidia_size_t h_, char * file_directory )
{
	delta.assign( delta_ );
	class_polynomial.set_degree( h_ );

	char input_file[ 200 ], delta_as_string[ 100 ];
	strcpy( input_file, file_directory );
	strcat( input_file, "/" );
	bigint_to_string( delta_, delta_as_string );
	strcat( input_file, delta_as_string );

	std::ifstream in( input_file );
	
	for( lidia_size_t kk = h_; kk >= 0; kk-- )
		in >> class_polynomial[ kk ];

	is_polynomial_set = true;
}

//
// Assignment by other instance of this class
//
void
gec_complex_multiplication::assign( const gec_complex_multiplication & I )
{
	if( this != &I )
	{
		is_initialized = I.is_initialized;
      
		if( is_initialized )
		{
			q   = I.q;
			F_q = I.F_q;
			r   = I.r;
			k   = I.k;
			G   = I.G;
	  
			lower_bound_bitlength_r = I.lower_bound_bitlength_r;
			upper_bound_k           = I.upper_bound_k;

			a2 = I.a2;
			a4 = I.a4;
			a6 = I.a6;

			according_to_BSI     = I.according_to_BSI;
			delta_field          = I.delta_field;
			lower_bound_h_field  = I.lower_bound_h_field;
		}
	}
}


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
