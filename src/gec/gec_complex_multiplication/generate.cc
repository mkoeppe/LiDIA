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
//  File    : generate.cc
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
// generate() implements the virtual function of the base class
// in the scope of the complex multiplication approach
//
void
gec_complex_multiplication::generate()
{
	if( VERBOSE )
	{
		std::cout << "Generating a curve via complex multiplication...." 
				  << std::endl;
	}

	//
	// * check whether delta = 0,1 mod 4
	// * check whether delta and upper_bound_k fit
	//
	check_the_settings();

	//
	// Assign default values if no user defined values are given
	//
	if( delta.is_zero() && q.is_zero() )
	{
		delta = default_delta;
		delta_case = default_delta_case;
	}

	// if the prime field is set, find an appropriate discriminant
	// using the Fixed Field Approach: It is implemented in
	// the function find_discriminant().
	if( ! q.is_zero() )
		if( ! find_discriminant() )
		{
			std::cout << "Cannot find discriminant. Exiting" << std::endl;
			
			exit( 1 );
		}

	//
	// We can already check if the discriminant meets the 
	// BSI/GISA requirements: If not, use the default values
	//
	if( according_to_BSI )
	{
	   	compute_delta_field();

		quadratic_order QO;
		QO.assign( delta_field );

		bigint h_field_as_bigint;
		h_field_as_bigint = QO.class_number();
		if( h_field_as_bigint.intify( h_field ) )
			lidia_error_handler( "gec_complex_multiplication",
								 "generate: cannot assign h_field_as_bigint to h_field" );

		lower_bound_h_field = h_field;
	   
		if( h_field < BSI_lower_bound_h_field )
		{
			std::cout << "h_field = " << h_field << " for your discriminant.";
			std::cout << std::endl << "As the BSI requires delta_field >= " <<
				BSI_lower_bound_h_field << " " << std::endl;
			std::cout << "we assign field discriminant of class number" << BSI_lower_bound_h_field << "." << std::endl;

			h_field = lower_bound_h_field = BSI_lower_bound_h_field;
			if( upper_bound_k == 1 )
			  {
			    delta = delta_field = default_delta_1;
			    delta_case = default_delta_case;
			  }
			else if( upper_bound_k < 4 )
			  {
			    delta = delta_field = default_delta_2;
			    delta_case = 0;
			  }
			else
			  {
			    delta = delta_field = default_delta_4;
			    delta_case = default_delta_case;
			  }
		}
	}

	// If no field is set, use the Fixed Discriminant Approach
	if( q.is_zero() )
	{
		// find an appropriate prime field
		if( degree == 1 )
		{
			if( ! find_good_prime() )
			{
				std::cout << "Cannot find prime. Exiting" << std::endl;
				
				exit( 1 );
			}
		}
		// find an appropriate prime power
		else
		{
			if( ! find_good_prime_power() )
			{
				std::cout << "Cannot find prime power. Exiting" << std::endl;

				exit( 1 );
			}
		}
	}

	// assign the field
	F_q.assign( galois_field( q ) );
	
	// if VERBOSE = true, print some information to standard output
	if( VERBOSE )
	{
		std::cout << "q = " << q << std::endl;
		std::cout << "p = " << p << std::endl;
		std::cout << "degree = " << degree << std::endl;
		std::cout << "r = " << r << std::endl;
		std::cout << "r.bit_length() = " << r.bit_length() << std::endl;       
		std::cout << "k = " << k << std::endl;
		std::cout << "k.bit_length() = " << k.bit_length() << std::endl;
		std::cout << "delta = " << delta << std::endl;
	}

	//
	// Compute the reduced representatives of discriminant delta
	//
	if( VERBOSE )
		std::cout << "Computing the class group....  ";

	// the reduced representatives are stored in the base_vector class_group
	class_group.set_mode( EXPAND );
	class_group = compute_class_group( delta ); // see chapter 6 of PhD
	h = class_group.size();

	if( VERBOSE )
		std::cout << "h = " << h << std::endl;

	//
	// Set the generation mode
	//
	if( generation_mode == 0 )
		set_generation_mode( generation_mode );

	//
	// Set the complex precision (only relevant, if no class polynomial
	// is set before)
	//
	if( complex_precision == 0 && is_polynomial_set == false )
		set_complex_precision( 0 );

	//
	// If no class polynomial is set before, compute it according
	// to generation_mode within the complex precision set before.
	//
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

	//
	// Compute some irreducible factor of degree degree of class polynomial
	//
	if( VERBOSE )
		std::cout << "Computing a factor of degree " << degree 
				  << " of minimal polynomial mod p...." << std::endl;

	// initialize the class polynomial mod p
	Fp_polynomial class_polynomial_mod_p( class_polynomial, p );
	
	// case of q = p: We do not search for a curve of the form (-3,b)
	if( degree == 1 && efficient_curve_parameters == false )
	{
		measure_time.start_timer();
		bigint root_of_class_polynomial_mod_p = 
			find_root( class_polynomial_mod_p );
		measure_time.stop_timer();

		if( timings == true )
		{
			out_ << ( float ) measure_time.real_time() / 100.0 << " & ";
		}
		
		// cast bigint root_of_class_polynomial_mod_p to 
		// gf_element class_invariant_mod_p
		gf_element class_invariant_mod_p( F_q );
		class_invariant_mod_p.assign( root_of_class_polynomial_mod_p );
	   
		if( VERBOSE )
			std::cout << "class_invariant_mod_p = " 
					  << class_invariant_mod_p << std::endl;
	   
		// Determine the curve parameters and a base point:
		// This is accomplished by the function 
		// assign_class_invariant_to_curve.
		if( ! assign_class_invariant_to_curve( class_invariant_mod_p ) )
		{
			std::cout << "Could not assign class_invariant_mod_p to curve" 
					  << std::endl;
			exit( 1 );
		}
	}
	// case of q = p: We search for a curve of the form (-3,b)
	else if( degree == 1 && efficient_curve_parameters == true )
		assign_efficient_curve_parameters( class_polynomial_mod_p );
	// case degree > 1
	else if( degree > 1 )
	{
		Fp_polynomial field_polynomial = 
			find_factor( class_polynomial_mod_p, degree );

		if( VERBOSE )
			std::cout << "field_polynomial = "
					  << field_polynomial << std::endl;

		if( ! det_irred_test( field_polynomial ) )
		{
			std::cout << "field_polynomial is not irreducible." << std::endl;
			exit( 1 );
		}

		if( VERBOSE )
			std::cout << "field_polynomial is irreducible." << std::endl;

		F_q.assign( galois_field( field_polynomial ) );

		Fp_polynomial class_invariant_pol;
		class_invariant_pol.set_modulus( p );
		class_invariant_pol.assign_x();

		gf_element class_invariant_Fq( F_q );
		class_invariant_Fq.set_polynomial_rep( class_invariant_pol );

		gf_element::set_output_format( 1 );
		std::cout << "class_invariant_Fq = "
				  << class_invariant_Fq << std::endl;

		if( ! assign_class_invariant_to_curve( class_invariant_Fq ) )
		{
			std::cout << "Could not assign class_invariant_mod_p to curve"
					  << std::endl;
			exit( 1 );
		}
	}

	// if VERBOSE = true, print the results to standard output
	if( VERBOSE )
	{
		std::cout << std::endl << std::endl;

		std::cout << "q = " << q << std::endl;
		std::cout << "k = " << k << std::endl;
		std::cout << "r = " << r << std::endl;
		std::cout << "r*k = " << r*k << std::endl;
		std::cout << "delta = " << delta << std::endl;
		std::cout << "h = " << h << std::endl;
		std::cout << "delta_field = " << delta_field << std::endl;
		std::cout << "lower_bound_h_field = "
				  << lower_bound_h_field << std::endl;
		std::cout << "generation_mode = " << generation_mode << std::endl;
		std::cout << "a4 = " << a4 << std::endl;
		std::cout << "a6 = " << a6 << std::endl;
		std::cout << "G = " << G << std::endl;
	}
	
	// set is_initialized to true
	is_initialized = true;
}


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
