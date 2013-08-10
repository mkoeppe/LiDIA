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
//  File    : assign_efficient_curve_parameters_and_twist.cc
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
// assign_efficient_curve_parameters_and_twist() tries to
// assign an elliptic curves of the form (-3,b) and (-3,b')
//

void
gec_E_and_twist_prime::
assign_efficient_curve_parameters_and_twist
( const Fp_polynomial & class_polynomial_mod_p )
{
	bool good_twisted_pair = false;
	bool first_curve_good;
	int counter = 0;
	const lidia_size_t degree_of_factor = 30;
	lidia_size_t degree_of_current_factor;

	F_q.assign( galois_field( p ) );
	gf_element minus_3( F_q );
	minus_3.assign( bigint( -3 ) );

	gf_element dummy_1( F_q ), dummy_2( F_q ), dummy_3( F_q );
	gf_element j( F_q );
	gf_element beta( F_q );
	gf_element invariant( F_q );

	gf_element gf_element_1728( F_q );
	gf_element_1728.assign( bigint( 1728 ) );

	base_vector<gf_element> curve_parameters( 4, FIXED );

	bigint root_of_class_polynomial_mod_p;

	Fp_polynomial current_polynomial_mod_p( class_polynomial_mod_p );
	Fp_polynomial linear_polynomial_mod_p;
	Fp_polynomial current_whole_factor_mod_p, current_factor_mod_p;
	linear_polynomial_mod_p.assign_x( p );

	point<gf_element> P_0, P_1, Q_0, Q_1;
	elliptic_curve<gf_element> e_0, e_1;
	
	// Search a quadratic nonresidue in F_q
	gf_element q_non_res( F_q );
	do
	{
		q_non_res.randomize();
		bigint exponent;
		divide( exponent, q - bigint(1), bigint(2));
		power( dummy_1, q_non_res, exponent );
	}
	while( dummy_1.is_one() || dummy_1.is_zero() );

	measure_time.start_timer();

	while( counter < h )
	{
		current_whole_factor_mod_p = find_factor_degree_d( 
			current_polynomial_mod_p, degree_of_factor );

		current_factor_mod_p.assign( current_whole_factor_mod_p );
		degree_of_current_factor = current_factor_mod_p.degree();

		lidia_size_t counter_2 = 0;

		while( counter_2 < degree_of_current_factor )
		{
			first_curve_good = false;

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

			if( generation_mode == 4 )
			{
				dummy_3.assign( bigint( 16 ) );
			
				power(dummy_1, invariant, bigint(24));
				subtract(dummy_3, dummy_1, dummy_3);
				power(dummy_2, dummy_3, bigint(3));
				divide(dummy_2, dummy_2, dummy_1);
			
				j = dummy_2;
			}     
			else if ( generation_mode == 3 )
			{
				power(dummy_1, invariant, bigint(3));
			
				j = dummy_1;
			}
			else
				j = invariant;
		
//
// exclude the cases j == 0 and j == 1728
//
			if( j.is_zero() )
			{
				std::cout << "j-invariant is zero. Aborting." << std::endl;
				exit(1);
			}
  
// define 1728 as gf_element
			if( j == gf_element_1728 )
			{
				std::cout << "j-invariant is 1728. Aborting." << std::endl;
				exit(1);
			}

//
// Compute the representatives of the two
// possible equivalence classes.
//
		
			// assign the coefficients of the first curve
			gf_element kappa( F_q );

			subtract(dummy_1, gf_element_1728, j);
			divide(kappa, j, dummy_1);
		
			multiply(curve_parameters[0], bigint(3), kappa); // a4=3*kappa
			multiply(curve_parameters[1], bigint(2), kappa); // a6=2*kappa
		
			square(dummy_1, q_non_res);
			multiply(curve_parameters[2], curve_parameters[0], dummy_1);
		
			multiply(dummy_1, dummy_1, q_non_res);
			multiply(curve_parameters[3], curve_parameters[1], dummy_1);
		
			e_0.set_coefficients( curve_parameters[0], curve_parameters[1] );
			e_1.set_coefficients( curve_parameters[2], curve_parameters[3] );

			if( VERBOSE )
			{
				std::cout << std::endl << "Curve candidates: " << std::endl;
				std::cout << "Curve 1: y^2 = x^3 + " <<  curve_parameters[0] << "x + " 
					 << curve_parameters[1] << std::endl;
				std::cout << "Curve 2: y^2 = x^3 + " <<  curve_parameters[2] << "x + " 
					 << curve_parameters[3] << std::endl;
				std::cout << std::endl;
			}
		
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
			
				if(P_0.is_zero())
					tester_0 = 1;
				if(P_1.is_zero())
					tester_1 = 1;
				
				if( (tester_0 == 1) && (tester_1 == 0) )
				{
					G = Q_0;
					a4 = curve_parameters[0];
					a6 = curve_parameters[1];
				
					G_tw = Q_1;
					a4_tw = curve_parameters[2];
					a6_tw = curve_parameters[3];
					
					break;
				}
			
				if((tester_0 == 0)&&(tester_1 == 1))
				{
					G = Q_1;
					a4 = curve_parameters[2];
					a6 = curve_parameters[3];
				
					G_tw = Q_0;
					a4_tw = curve_parameters[0];
					a6_tw = curve_parameters[1];
					
					break;
				}
			}
			
			// test whether curve is isomorph to one of the form (-3,b)
			divide( dummy_1, minus_3, a4 );
		
			if( dummy_1.is_square() )
			{
				beta = sqrt( dummy_1 );
				if( beta.is_square() )
				{
					a4.assign( minus_3 );
					power( dummy_1, beta, 3 );
					multiply( a6, a6, dummy_1 );

					first_curve_good = true;
				}
				else
				{
					beta.negate();

					if( beta.is_square() )
					{
						a4.assign( minus_3 );
						power( dummy_1, beta, 3 );
						multiply( a6, a6, dummy_1 );
					
						first_curve_good = true;
					}
				}
			}

			if( first_curve_good == true )
			{
				divide( dummy_1, minus_3, a4_tw );
		
				if( dummy_1.is_square() )
				{
					beta = sqrt( dummy_1 );
					if( beta.is_square() )
					{
						a4_tw.assign( minus_3 );
						power( dummy_1, beta, 3 );
						multiply( a6_tw, a6_tw, dummy_1 );

						good_twisted_pair = true;
						break;
					}
					else
					{
						beta.negate();
						
						if( beta.is_square() )
						{
							a4_tw.assign( minus_3 );
							power( dummy_1, beta, 3 );
							multiply( a6_tw, a6_tw, dummy_1 );
					
							good_twisted_pair = true;
							break;
						}
					}
				}
			}

			linear_polynomial_mod_p[ 0 ] = -root_of_class_polynomial_mod_p;

			divide( current_factor_mod_p, current_factor_mod_p,
					linear_polynomial_mod_p );

			counter++;
			counter_2++;

			std::cout << "counter = " << counter << std::endl;
		}

		if( good_twisted_pair == true )
			break;
		
		divide( current_polynomial_mod_p, current_polynomial_mod_p,
				current_whole_factor_mod_p );
	}
	
	measure_time.stop_timer();
	
	e_0.set_coefficients( a4, a6 );
	e_1.set_coefficients( a4_tw, a6_tw );

	do{
		P_0 = e_0.random_point();
		P_1 = e_1.random_point();
	}
	while( P_0.is_zero() || P_1.is_zero() );

	multiply(Q_0, k, P_0); // Q_0 = k * P_0
	multiply(P_0, r, Q_0); // P_0 = r * Q_0

	multiply(Q_1, k_tw, P_1); // Q_1 = k_tw * P_1
	multiply(P_1, r_tw, Q_1); // P_1 = r_tw * Q_1
			
	if(P_0.is_zero() && P_1.is_zero() && counter < h )
	{
		G = Q_0;
		G_tw = Q_1;

		if( VERBOSE )
			std::cout << "counter = " << counter << std::endl;
	}
	else if( counter == h )
	{
		std::cout << "There is no curve with a4 = a4_tw = -3." << std::endl;
		exit( 1 );
	}
	else
	{
		std::cout << "Not the right order." << std::endl;
		exit( 1 );
	}	
	
	if( timings == true )
		out_ << ( float ) measure_time.real_time() / 100.0 << " & ";
}


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
