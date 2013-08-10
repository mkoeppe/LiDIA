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
//  File    : assign_class_invariant_to_curve_and_twist.cc
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


bool
gec_E_and_twist_prime::
assign_class_invariant_to_curve_and_twist( const gf_element & invariant )
{ 
	//
	// Compute a j-invariant j mod q of a curve of order r*k
	//
	gf_element dummy_1( F_q ), dummy_2( F_q ), dummy_3( F_q );
	gf_element j( F_q );

	if( generation_mode == 4 )
	{
		dummy_3.assign( bigint(16) );

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
	gf_element gf_element_1728( F_q );
	gf_element_1728.assign( bigint( 1728 ) );

	if( j == gf_element_1728 )
	{
		std::cout << "j-invariant is 1728. Aborting." << std::endl;
		exit(1);
	}

	//
	// Compute the representatives of the two
	// possible equivalence classes.
	//
	base_vector<gf_element> curve_parameters( 4, FIXED );

	// assign the coefficients of the first curve
	gf_element kappa( F_q );

	subtract(dummy_1, gf_element_1728, j);
	divide(kappa, j, dummy_1);

	multiply(curve_parameters[0], bigint(3), kappa); // a4=3*kappa
	multiply(curve_parameters[1], bigint(2), kappa); // a6=2*kappa

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

	square(dummy_1, q_non_res);
	multiply(curve_parameters[2], curve_parameters[0], dummy_1);

	multiply(dummy_1, dummy_1, q_non_res);
	multiply(curve_parameters[3], curve_parameters[1], dummy_1);

	point<gf_element> P_0, P_1, Q_0, Q_1;
	elliptic_curve<gf_element> e_0( curve_parameters[0], curve_parameters[1] );
	elliptic_curve<gf_element> e_1( curve_parameters[2], curve_parameters[3] );

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

		if( P_0.is_zero() )
			tester_0 = 1;
		if( P_1.is_zero() )
			tester_1 = 1;
        
		if( (tester_0 == 1) && (tester_1 == 0) )
		{
			multiply( P_1, r_tw, Q_1 );

			if( P_1.is_zero() )
			{
				G = Q_0;
				a4 = curve_parameters[0];
				a6 = curve_parameters[1];
				
				G_tw = Q_1;
				a4_tw = curve_parameters[2];
				a6_tw = curve_parameters[3];
				
				return( 1 );
			}
		}

		if((tester_0 == 0)&&(tester_1 == 1))
		{
			multiply( P_0, r_tw, Q_0 );

			if( P_0.is_zero() )
			{
				G = Q_1;
				a4 = curve_parameters[2];
				a6 = curve_parameters[3];

				G_tw = Q_0;
				a4_tw = curve_parameters[0];
				a6_tw = curve_parameters[1];

				return( 1 );
			}
		}
	}
	
	return( 0 );
}


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
