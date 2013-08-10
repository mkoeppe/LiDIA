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
//  File    : generate_twist.cc
//  Author  : Harald Baier (HB)
//	Changes	: See CVS log
//
//==============================================================================================

#ifdef HAVE_CONFIG_H
# include       "config.h"
#endif
# include        "LiDIA/gec_E_and_twist_prime.h"
#ifdef VERBOSE
# include        <cassert>
#endif


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// generate() implements the virtual function of the base class
// in the scope of gec_E_and_twist_prime
//
void
gec_E_and_twist_prime::generate()
{
	if( VERBOSE )
	{
		std::cout << "Generating curve and twist of prime order ";
		std::cout << "via complex multiplication...." << std::endl;
	}

	check_the_settings();

	//
	// Assign default values if no user defined values are given
	//
	if( delta.is_zero() )
	{
		delta = default_delta_E_and_twist_prime;
		delta_case = default_delta_case_E_and_twist_prime;
	}
	else
	{
		bigint remainder_mod_8, remainder_mod_3;

		remainder( remainder_mod_8, delta, bigint( 8 ) );
		remainder( remainder_mod_3, delta, bigint( 3 ) );

		if( remainder_mod_8 != -3 || remainder_mod_3 != -1 )
		{
			std::cout << "delta is not = 5 mod 8 and = 2 mod 3. Exiting." << std::endl;
			exit( 1 );
		}
	}

	//
	// We can already check if the set discriminant meets the 
	// BSI requirements
	//
	if( according_to_BSI )
	{
	   	compute_delta_field();

		quadratic_order QO;
		QO.assign( delta_field );

		bigint h_field_as_bigint;
		h_field_as_bigint = QO.class_number();
		if( h_field_as_bigint.intify( h_field ) )
			lidia_error_handler( "gec_E_and_twist_prime",
								 "generate: cannot assign h_field_as_bigint to h_field" );

		lower_bound_h_field = h_field;
	   
		if( VERBOSE )
		{
			std::cout << "delta = " << delta << std::endl;
			std::cout << "h = " << h << std::endl;
			std::cout << "delta_field = " << delta_field << std::endl;
			std::cout << "h_field = " << h_field << std::endl;
			std::cout << "lower_bound_h_field = " << lower_bound_h_field << std::endl;
		}

		if( h_field < BSI_lower_bound_h_field )
		{
			std::cout << "h_field = " << h_field << " for your discriminant.";
			std::cout << std::endl << "As the BSI requires delta_field >= " <<
				BSI_lower_bound_h_field << " " << std::endl;
			std::cout << "we assign delta = " << default_delta_E_and_twist_prime << std::endl;

			h_field = lower_bound_h_field = BSI_lower_bound_h_field; 
			delta = delta_field = default_delta_E_and_twist_prime;
			delta_case = default_delta_case_E_and_twist_prime;
		}
	}

	if( VERBOSE )
	{
		std::cout << std::endl << "*************************" << std::endl;
		std::cout << "We use the following parameters: " << std::endl;
		std::cout << "r.bit_length() >= " << lower_bound_bitlength_r << std::endl;
		std::cout << "k <= " << upper_bound_k <<std::endl;
		std::cout << "delta = " << delta << std::endl;
		std::cout << "delta_case = " << delta_case << std::endl;
		std::cout << "degree = " << degree << std::endl;
		std::cout << std::endl;
	}

	if( degree <= 1 )
	{
		if( ! find_good_prime() )
			std::cout << "Cannot find prime number." << std::endl;
	}

	if( VERBOSE )
	{
		std::cout << "find_good_prime() finished successfully: " << std::endl;
		std::cout << "q = " << q << std::endl;
		std::cout << "p = " << p << std::endl;
		std::cout << "degree = " << degree << std::endl;
		std::cout << "r = " << r << std::endl;
		std::cout << "r.bit_length() = " << r.bit_length() << std::endl;       
		std::cout << "k = " << k << std::endl;
		std::cout << "k.bit_length() = " << k.bit_length() << std::endl;
	}

	if( VERBOSE )
		std::cout << "Computing the class group...." << std::endl;

	class_group.set_mode( EXPAND );
	class_group = compute_class_group( delta );
	h = class_group.size();

	if( VERBOSE )
		std::cout << "h = " << h << std::endl;

	if( generation_mode == 0 )
		set_generation_mode( generation_mode );

	if( complex_precision == 0 && is_polynomial_set == false )
		set_complex_precision( 0 );

	if( is_polynomial_set == false )
	{
		if( VERBOSE )
			std::cout << "Computing a class polynomial...." << std::endl;
		
		measure_time.start_timer();	
		compute_class_polynomial();
		measure_time.stop_timer();
		
		if( timings == true )
			out_ << ( float ) measure_time.real_time() / 100.0 << " & ";
	}

	if( VERBOSE )
		std::cout << "Computing a factor of degree " 
				  << degree << " of minimal polynomial mod p...." << std::endl;

	// initialize the class polynomial mod p
	Fp_polynomial class_polynomial_mod_p( class_polynomial, p );

	if( degree <= 1 && efficient_curve_parameters == false )
	{
		measure_time.start_timer();
		bigint root_of_class_polynomial_mod_p = 
			find_root( class_polynomial_mod_p );
		measure_time.stop_timer();

		if( timings == true )
		{
			out_ << ( float ) measure_time.real_time() / 100.0 << " & ";
//			out_ << "\\\\ \\hline" << std::endl;
		}
		
		F_q.assign( galois_field( q ) );

		//
		// cast bigint root_of_class_polynomial_mod_p to 
		// gf_element class_invariant_mod_p
		//
		gf_element class_invariant_mod_p( F_q );
		class_invariant_mod_p.assign( root_of_class_polynomial_mod_p );
	   
		if( VERBOSE )
			std::cout << "class_invariant_mod_p = " << class_invariant_mod_p << std::endl;
	   
		//
		// Determine the curve parameters and the base points:
		//
		// This is accomplished by the function 
		// assign_class_invariant_to_curve.
		//
		if( ! assign_class_invariant_to_curve_and_twist( 
			class_invariant_mod_p ) )
		{
			std::cout << "Could not assign class_invariant_mod_p to curve" << std::endl;
			exit( 1 );
		}
	}
	//
	// case of q = p: We search for a curve of the form (-3,b)
	//
	else if( degree <= 1 && efficient_curve_parameters == true )
		assign_efficient_curve_parameters_and_twist( class_polynomial_mod_p );

	//
	// case degree > 1
	//
	if( degree > 1 )
	{
		Fp_polynomial field_polynomial = 
			find_factor( class_polynomial_mod_p, degree );

		std::cout << "field_polynomial = " << field_polynomial << std::endl;

		if( det_irred_test( field_polynomial ) )
			std::cout << "field_polynomial is irreducible." << std::endl;
		else
		{
			std::cout << "field_polynomial is not irreducible." << std::endl;
			exit( 1 );
		}

		F_q.assign( galois_field( field_polynomial ) );
	   
		Fp_polynomial class_invariant_pol;
		class_invariant_pol.set_modulus( p );
		class_invariant_pol.assign_x();

		gf_element class_invariant_Fq( F_q );
		class_invariant_Fq.set_polynomial_rep( class_invariant_pol );

		gf_element::set_output_format( 1 );
		std::cout << "class_invariant_Fq = " << class_invariant_Fq << std::endl;

		if( ! assign_class_invariant_to_curve( class_invariant_Fq ) )
			std::cout << "Could not assign class_invariant_mod_p to curve" << std::endl;

	}

	if( VERBOSE )
	{
		std::cout << "*******************" << std::endl;

		std::cout << "q = " << q << std::endl;
		std::cout << "-----------------" << std::endl;
		std::cout << "delta = " << delta << std::endl;
		std::cout << "h = " << h << std::endl;
		std::cout << "delta_field = " << delta_field << std::endl;
		std::cout << "lower_bound_h_field = "
				  << lower_bound_h_field << std::endl;
		std::cout << "generation_mode = " << generation_mode << std::endl;
		std::cout << "-----------------" << std::endl;
		std::cout << "k = " << k << std::endl;
		std::cout << "r = " << r << std::endl;
		std::cout << "a4 = " << a4 << std::endl;
		std::cout << "a6 = " << a6 << std::endl;
		std::cout << "G = " << G << std::endl;
		std::cout << "-----------------" << std::endl;
		std::cout << "k_tw = " << k_tw << std::endl;
		std::cout << "r_tw = " << r_tw << std::endl;
		std::cout << "a4_tw = " << a4_tw << std::endl;
		std::cout << "a6_tw = " << a6_tw << std::endl;
		std::cout << "G_tw = " << G_tw << std::endl;
	}

	is_initialized = true;

}


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
