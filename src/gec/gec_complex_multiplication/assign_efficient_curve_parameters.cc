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
//  File    : assign_efficient_curve_parameters.cc
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
// assign_efficient_curve_parameters() tries to
// assign an elliptic curve of the form (-3,b).
// This function assumes F_q to be a prime field.
//
void
gec_complex_multiplication::assign_efficient_curve_parameters
( const Fp_polynomial & class_polynomial_mod_p )
{
	//
	// declare and assign some variables
	//

	// a boolean to decide whether a good curve is found
	bool good_curve = false;
	// a counting variable
	int counter = 0;
	// an upper bound of the degree of the factor of the class polynomial
	const lidia_size_t degree_of_factor = 30;
	// the actual degree of the current factor of the class polynomial
	lidia_size_t degree_of_current_factor;
	// some temporary variables
	gf_element tmp_gf_1( F_q ), tmp_gf_2( F_q ), tmp_gf_3( F_q );
	// the j-invariant
	gf_element j( F_q );
	// for testing if the current curve is isomorphic to a curve (-3,b)
	gf_element beta( F_q );
	// a root of the class polynomial, that is a class invariant
	gf_element invariant( F_q );
	// kappa is used to assign curve parameters for given j-invariant
	gf_element kappa( F_q );

	// the constants -3 and 1728 in F_q
	gf_element minus_3( F_q ), gf_element_1728( F_q );
	minus_3.assign( bigint( -3 ) );
	gf_element_1728.assign( bigint( 1728 ) );

	// the curve candidates are (a4_, a6_) and a twisted curve (a4_tw, a6_tw)
	gf_element a4_( F_q ), a6_( F_q ), a4_tw( F_q ), a6_tw( F_q );

	// a root of the class polynomial modulo p
	bigint root_of_class_polynomial_mod_p;

	//
	// some polynomials needed during the computation
	//

	// the current polynomial: class polynomial / current_whole_factor_mod_p
	Fp_polynomial current_polynomial_mod_p( class_polynomial_mod_p );
	// the current linear factor
	Fp_polynomial linear_polynomial_mod_p;
	linear_polynomial_mod_p.assign_x( p );
	// the current factor of the class polynomial
	Fp_polynomial current_whole_factor_mod_p;
	// current_factor_mod_p stores current_whole_factor_mod_p divided
	// by all linear factors
	Fp_polynomial current_factor_mod_p;

	// two elliptic curves and two points on each curve
	elliptic_curve<gf_element> e_0, e_1;
	point<gf_element> P_0, P_1, Q_0, Q_1;
	
	// Search a quadratic nonresidue in F_q
	gf_element q_non_res( F_q );
	bigint exponent;
	do
	{
		q_non_res.randomize();
		divide( exponent, q - bigint(1), bigint(2));
		power( tmp_gf_1, q_non_res, exponent );
	}
	while( tmp_gf_1.is_one() || tmp_gf_1.is_zero() );

	// to measure the time for finding the parameters
	measure_time.start_timer();
	
	// counter stores the number of tested invariants
	while( counter < h )
	{
		// compte a factor of degree at most degree_of_factor of
		// current_polynomial_mod_p
		current_whole_factor_mod_p = find_factor_degree_d( 
			current_polynomial_mod_p, degree_of_factor );

		current_factor_mod_p.assign( current_whole_factor_mod_p );
		degree_of_current_factor = current_factor_mod_p.degree();

		lidia_size_t counter_2 = 0;
		while( counter_2 < degree_of_current_factor )
		{
			root_of_class_polynomial_mod_p = 
				find_root( current_factor_mod_p );

			//
			// cast bigint root_of_class_polynomial_mod_p to 
			// gf_element invariant
			//
			invariant.assign( root_of_class_polynomial_mod_p );
		
			if( VERBOSE )
				std::cout << "current class invariant mod p = "
					 << invariant << std::endl;

			// assign root_of_class_polynomial_mod_p to j-invariant
			if( generation_mode == 4 )
			{
				tmp_gf_3.assign( bigint( 16 ) );
			
				power(tmp_gf_1, invariant, bigint(24));
				subtract(tmp_gf_3, tmp_gf_1, tmp_gf_3);
				power(tmp_gf_2, tmp_gf_3, bigint(3));
				divide(tmp_gf_2, tmp_gf_2, tmp_gf_1);
			
				j = tmp_gf_2;
			}     
			else if ( generation_mode == 3 )
			{
				power(tmp_gf_1, invariant, bigint(3));
			
				j = tmp_gf_1;
			}
			else
				j = invariant;
			
			//
			// exclude the cases j == 0 and j == 1728
			//
			if( j.is_zero() )
			{
				std::cout << "j-invariant is zero. Aborting." << std::endl;
				exit( 1 );
			}
			if( j == gf_element_1728 )
			{
				std::cout << "j-invariant is 1728. Aborting." << std::endl;
				exit( 1 );
			}

			//
			// Compute two twisted curves with j-invariant j.
			//
			
			// assign the coefficients of the first curve
			subtract(tmp_gf_1, gf_element_1728, j);
			divide(kappa, j, tmp_gf_1); // kappa = j / (1728 - j)
			multiply(a4_, bigint( 3 ), kappa); // a4_ = 3*kappa
			multiply(a6_, bigint( 2 ), kappa); // a6_ = 2*kappa

			// assign coefficients of twisted curve
			square(tmp_gf_1, q_non_res);
			multiply(a4_tw, a4_, tmp_gf_1);
			multiply(tmp_gf_1, tmp_gf_1, q_non_res);
			multiply(a6_tw, a6_, tmp_gf_1);

			// assign curve parameters
			e_0.set_coefficients( a4_, a6_ );
			e_1.set_coefficients( a4_tw, a6_tw );

			if( VERBOSE )
			{
				std::cout << std::endl << "Curve candidates: " << std::endl;
				std::cout << "Curve 1: y^2 = x^3 + " <<  a4_ << "x + " 
					 << a6_ << std::endl;
				std::cout << "Curve 2: y^2 = x^3 + " <<  a4_tw << "x + " 
					 << a6_tw << std::endl;
				std::cout << std::endl;
			}

			// we allow at most 5 tests to find curve of order r*k
			for( int i = 1; i <= 5; i++ )
			{
				int tester_0 = 0, tester_1 = 0; // to assign the proper curve
			
				// generate random points P_0 on e_0 and P_1 on e_1
				// with k*P_0 != 0 and k*P_1 != 0
				do
				{
					P_0 = e_0.random_point();
					P_1 = e_1.random_point();
				
					multiply(Q_0, k, P_0); // Q_0 = k * P_0
					multiply(Q_1, k, P_1); // Q_1 = k * P_1
				} while( Q_0.is_zero() || Q_1.is_zero() );
			
				multiply(P_0, r, Q_0); // P_0 = r * Q_0
				multiply(P_1, r, Q_1); // P_1 = r * Q_1      
			
				if( P_0.is_zero() )
					tester_0 = 1;
				if( P_1.is_zero() )
					tester_1 = 1;
			
				if( ( tester_0 == 1 ) && ( tester_1 == 0 ) )
				{
					a4 = a4_;
					a6 = a6_;
				
					break;
				}
			
				if( ( tester_0 == 0 ) && ( tester_1 == 1 ) )
				{
					a4 = a4_tw;
					a6 = a6_tw;
				
					break;
				}
			}

			//
			// test whether curve (a4, a6) \sim (-3,b),
			// that is -3 = a4 * beta^2 with a square beta in F_q^*
			// hence -3 / a4 = beta^2
			//

			// compute -3 / a4
			divide( tmp_gf_1, minus_3, a4 );

			// decide if -3 / a4 is a square in F_q^*
			if( tmp_gf_1.is_square() )
			{
				// compute square root of -3 / a4 and decide whether
				// it is a square in F_q^*
				beta = sqrt( tmp_gf_1 );
				if( beta.is_square() )
				{
					a4.assign( minus_3 );
					power( tmp_gf_1, beta, 3 );
					multiply( a6, a6, tmp_gf_1 );

					// we are done and leave the while-loop
					good_curve = true;
					break;
				}
				// if p = 3 mod 4, -1 is not a square; hence -beta is a square
				else if( p.bit( 1 ) == 1 )
				{
					beta.negate();

					if( beta.is_square() )
					{
						a4.assign( minus_3 );
						power( tmp_gf_1, beta, 3 );
						multiply( a6, a6, tmp_gf_1 );

						// we are done and leave the while-loop
						good_curve = true;
						break;
					}
				}
			}

			// if we are not successful, we divide the current
			// polynomial by the tested linear factor
			linear_polynomial_mod_p[ 0 ] = -root_of_class_polynomial_mod_p;
			divide( current_factor_mod_p, current_factor_mod_p,
					linear_polynomial_mod_p );

			// increase both counters
			counter++;
			counter_2++;
		}
		
		// we are done and leave the while-loop
		if( good_curve == true )
			break;
		
		// otherwise we divide the class polynomial by the current factor
		divide( current_polynomial_mod_p, current_polynomial_mod_p,
				current_whole_factor_mod_p );
	}

	// stop timing
	measure_time.stop_timer();

	// set curve to (-3,b)
	e_0.set_coefficients( a4, a6 );


	// find a base point of order r
	do
	{
		P_0 = e_0.random_point();
		multiply(Q_0, k, P_0); // Q_0 = k * P_0
	}
	while( Q_0.is_zero() );

	multiply(P_0, r, Q_0); // P_0 = r * Q_0
			
	if( P_0.is_zero() && counter < h )
	{
		G = Q_0;

		if( VERBOSE )
			std::cout << "We tested " << counter << " roots." << std::endl;
	}
	else if( counter == h )
	{
		std::cout << "There is no curve with a4 = -3." << std::endl;
		exit( 1 );
	}
	else
	{
		std::cout << "Not the right order." << std::endl;
		exit( 1 );
	}	

	// if timings are relevant, write them to output stream out_
	if( timings == true )
		out_ << ( float ) measure_time.real_time() / 100.0 << " & ";
}


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
