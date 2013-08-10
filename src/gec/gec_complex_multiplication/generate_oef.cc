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
//  File    : generate_oef.cc
//  Author  : Harald Baier (HB)
//	Changes	: See CVS log
//
//==============================================================================================

#ifdef HAVE_CONFIG_H
# include       "config.h"
#endif
# include        "LiDIA/gec_complex_multiplication.h"
# include        "LiDIA/Fp_poly_modulus.h"

#ifdef VERBOSE
# include        <cassert>
#endif


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif


//
// generate_oef(lidia_size_t word_size_, lidia_size_t degree_, lidia_size_t )
// implements algorithm oefCurve of PhD (11.6)
// the user has to ensure word_size_ * degree_ >= lower_bound_bitlength_r
//
void
gec_complex_multiplication::
generate_oef( lidia_size_t word_size_, lidia_size_t degree_, 
			  lidia_size_t lower_bound_h_ )
{
	if( VERBOSE )
		std::cout << "Generating a curve over an OEF...." << std::endl;

	if( word_size_ > 1 )
		word_length = word_size_;
	else
		lidia_error_handler( "gec_complex_multiplication",
							 "generate_oef(): word_size_ has to be > 1");

	if( degree_ > 1 && is_prime( degree_ ) )
		degree = degree_;
	else
		lidia_error_handler( "gec_complex_multiplication",
							 "generate_oef(): degree has to be a prime");

   	// Discriminants are searched in the class number interval specified by
	// the following bounds; if the bound is too low,
	// the default value is used.
	if( lower_bound_h_ < default_lower_bound_class_number )
	{
		std::cout << "We set lower_bound_h_ = "
				  << default_lower_bound_class_number << std::endl;

		lower_bound_class_number = default_lower_bound_class_number;
	}
	else if( lower_bound_h_ > default_upper_bound_class_number )
	{
		std::cout << "We only support discriminants up to class number "
				  << default_upper_bound_class_number << std::endl;
		std::cout << "Aborting." << std::endl;
		
		exit( 1 );
	}
	else
		lower_bound_class_number = lower_bound_h_;

	if( upper_bound_class_number < lower_bound_class_number ||
		upper_bound_class_number > default_upper_bound_class_number )
		upper_bound_class_number = default_upper_bound_class_number;

	//
	// Assign default values if no user defined values are given
	//
	if( delta.is_zero() && q.is_zero() )
	{
		delta = default_delta;
		delta_case = default_delta_case;
	}

	// we want to generate a curve of prime order
	upper_bound_k = 1;

	if( lower_bound_bitlength_r == 0 )
		lower_bound_bitlength_r = default_lower_bound_bitlength_r;

	if( degree * word_length < lower_bound_bitlength_r )
	{
		std::cout << "You have to ensure degree * word_length ";
		std::cout << " >= lower_bound_bitlength_r." << std::endl;
		std::cout << "Aborting." << std::endl;
		
		exit( 1 );		
	}

	if( ! find_oef() )
	{
		std::cout << "find_oef() returns false. Exiting." << std::endl;
		
		exit( 1 );
	}
	
	//
	// find an irreducible polynomial X^degree - omega:
	// it implements algorithm findBinomial (11.4) of PhD
	//
	bigint omega( 1 );
	bigint degree_as_bigint( degree ), rem( 0 ), R, degree_to_exp_min_one;
	lidia_size_t exponent_of_degree = 0;

	R = p - 1;

	while( rem.is_zero() )
	{
		divide( R, R, degree_as_bigint );
		remainder( rem, R, degree_as_bigint );

		exponent_of_degree++;
	}

	if( exponent_of_degree > 1 )
		power( degree_to_exp_min_one, degree_as_bigint, exponent_of_degree - 1 );

	while( true )
	{
		omega = next_prime( omega );

		power_mod( rem, omega, R, p );
		if( rem.is_one() )
			continue;

		if( exponent_of_degree > 1 )
		{
			power_mod( rem, rem, degree_to_exp_min_one, p );

			if( rem.is_one() )
				continue;
		}

		// an irreducible binomial is found
		break;
	}

	if( VERBOSE )
		std::cout << "Computing the class group....  ";

	class_group.set_mode( EXPAND );
	class_group = compute_class_group( delta );
	h = class_group.size();

	if( VERBOSE )
		std::cout << "h = " << h << std::endl;

	if( generation_mode == 0 )
		set_generation_mode( generation_mode );

	if( complex_precision == 0 )
		set_complex_precision( 0 ); 

	if( VERBOSE )
		std::cout << "Computing a class polynomial...." << std::endl;
   
	compute_class_polynomial();

	if( VERBOSE )
	{
		std::cout << "Computing a factor of degree " << degree
				  << " of minimal polynomial mod p...." << std::endl;
	}
	
	// initialize the class polynomial mod p
	Fp_polynomial class_polynomial_mod_p( class_polynomial, p );

	// case degree > 1
	if( degree > 1 )
	{
		Fp_polynomial field_polynomial;
		field_polynomial.set_modulus( p );
		field_polynomial.set_coefficient( degree );
		field_polynomial[ 0 ] = - omega;
		
		Fp_polynomial factor_of_degree_d =
			find_factor( class_polynomial_mod_p, degree );

		std::cout << "field_polynomial = " << field_polynomial << std::endl;
		std::cout << "factor_of_degree_d = "
				  << factor_of_degree_d << std::endl;

		if( det_irred_test( field_polynomial ) )
			std::cout << "field_polynomial is irreducible." << std::endl;
		else
		{
			std::cout << "field_polynomial is not irreducible." << std::endl;
			exit( 1 );
		}

		F_q.assign( galois_field( field_polynomial ) );

		gf_polynomial factor_of_degree_d_gf( F_q );
		factor_of_degree_d_gf.set_degree( degree );

		gf_element tmp_gf( F_q );

		tmp_gf.assign( bigint( 1 ) );
		factor_of_degree_d_gf[ degree ].assign( tmp_gf );

		for( lidia_size_t kk = 0; kk < degree; kk++ )
		{
			tmp_gf.assign( factor_of_degree_d[ kk ] );
			factor_of_degree_d_gf[ kk ].assign( tmp_gf );
		}
	   
		std::cout << "factor_of_degree_d = " << factor_of_degree_d_gf << std::endl;

		gf_element class_invariant_Fq( F_q );
		class_invariant_Fq = find_root( factor_of_degree_d_gf );

		gf_element::set_output_format( 1 );
		std::cout << "class_invariant_Fq = "
				  << class_invariant_Fq << std::endl;

		if( ! assign_class_invariant_to_curve( class_invariant_Fq ) )
		{
			std::cout << "Could not assign class_invariant_mod_p to curve" << std::endl;
			exit( 1 );
		}
	}
	else
	{
		std::cout << "Degree is " << degree << ". However,";
		std::cout << "it has to be a prime. Exiting." << std::endl;

		exit( 1 );
	}

	if( VERBOSE )
	{
		std::cout << std::endl << "*******************"
				  << std::endl << std::endl;

		std::cout << "q = " << q << std::endl;
		std::cout << "q.bit_length = " << q.bit_length() << std::endl;
		std::cout << "p = " << p << std::endl;
		std::cout << "p.bit_length = " << p.bit_length() << std::endl;
		std::cout << "degree = " << degree << std::endl;
		std::cout << "k = " << k << std::endl;
		std::cout << "r = " << r << std::endl;
		std::cout << "r.bit_length = " << r.bit_length() << std::endl;
		std::cout << "delta = " << delta << std::endl;
		std::cout << "h = " << h << std::endl;
		std::cout << "delta_field = " << delta_field << std::endl;
		std::cout << "lower_bound_h_field = " << lower_bound_h_field << std::endl;
		std::cout << "a4 = " << a4 << std::endl;
		std::cout << "a6 = " << a6 << std::endl;
		std::cout << "G = " << G << std::endl;
	}

	is_initialized = true;

}

#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
