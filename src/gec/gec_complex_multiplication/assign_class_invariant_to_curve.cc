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
//  File    : assign_class_invariant_to_curve.cc
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
// the function assign_class_invariant_to_curve computes
// two twisted curves of j-invariant invariant. It tests
// if one of the curves is of order r*k. If yes, the parameters
// are assigned to a4 and a6, and true is returned.
//
// Otherwise, the method returns false.
//
bool
gec_complex_multiplication::
assign_class_invariant_to_curve( const gf_element & invariant )
{ 
	// three temporary variables of type gf_element
	gf_element tmp_gf_1( F_q ), tmp_gf_2( F_q ), tmp_gf_3( F_q );
	// the j-invariant
	gf_element j( F_q );

	//
	// Compute a j-invariant j in F_q of a curve of order r*k
	//

	// if f(tau) is a class invariant, compute corresponding j
	if( generation_mode == 4 )
	{
		tmp_gf_3.assign( bigint(16) );

		power(tmp_gf_1, invariant, bigint(24));
		subtract(tmp_gf_3, tmp_gf_1, tmp_gf_3);
		power(tmp_gf_2, tmp_gf_3, bigint(3));
		divide(tmp_gf_2, tmp_gf_2, tmp_gf_1);
      
		j = tmp_gf_2;
	}

	// if gamma_2(tau) is a class invariant, compute corresponding j
	else if( generation_mode == 3 )
		power( j, invariant, bigint( 3 ) );

	// otherwise we have j = invariant
	else
		j = invariant;

	//
	// We have to exclude the cases j == 0 and j == 1728
	//
	
	// ensure j != 0 in F_q
	if( j.is_zero() )
	{
		std::cout << "j-invariant is zero. Aborting." << std::endl;
		exit( 1 );
	}
	// ensure j != 1728 in F_q
	tmp_gf_1.assign( bigint( 1728 ) );
	if( j == tmp_gf_1 )
	{
		std::cout << "j-invariant is 1728. Aborting." << std::endl;
		exit(1);
	}

	//
	// Compute two twisted elliptic curves; one of them is of order r*k
	//
	
	// the curves are (a4_, a6_) and a twisted curve (a4_tw, a6_tw)
	gf_element a4_( F_q ), a6_( F_q ), a4_tw( F_q ), a6_tw( F_q );
	gf_element kappa( F_q );

	// compute kappa = j / ( 1728 - j )
	subtract( tmp_gf_2, tmp_gf_1, j );
	divide( kappa, j, tmp_gf_2 );

	// assign the coefficients of the first curve
	multiply( a4_, bigint( 3 ), kappa ); // a4_ = 3*kappa
	multiply( a6_, bigint(2), kappa ); // a6_ = 2*kappa

	// Compute a quadratic nonresidue in F_q by trial-and-error
	gf_element q_non_res( F_q );
	bigint e_; // exponent e_
	do
	{
		q_non_res.randomize(); // select a random element
		divide( e_, q - bigint( 1 ), bigint( 2 ) ); // e_ = (q-1)/2
		power( tmp_gf_1, q_non_res, e_ );
	}
	while( tmp_gf_1.is_one() || tmp_gf_1.is_zero() );

	// compute the parameters of the twisted curve
	square(tmp_gf_1, q_non_res);
	multiply(a4_tw, a4_, tmp_gf_1); // a4_tw = a4_ * q_non_res^2
	multiply(tmp_gf_1, tmp_gf_1, q_non_res);
	multiply(a6_tw, a6_, tmp_gf_1); // a6_tw = a6_ * q_non_res^3

	// assign the elliptic curves
	elliptic_curve<gf_element> e_0( a4_, a6_ );
	elliptic_curve<gf_element> e_1( a4_tw, a6_tw );

	// if VERBOSE == true, write the parameters to standard output
	if( VERBOSE )
	{
		std::cout << std::endl << "Curve candidates: " << std::endl;
		std::cout << "Curve 1: y^2 = x^3 + " <<  a4_ << "x + " 
				  << a6_ << std::endl;
		std::cout << "Curve 2: y^2 = x^3 + " <<  a4_tw << "x + " 
				  << a6_tw << std::endl;
		std::cout << std::endl;
	}

	// four points for finding the curve of order r*k
	point<gf_element> P_0, P_1, Q_0, Q_1;

	// we allow at most 5 tests; otherwise none of the curves is of order r*k
	for( int i = 1; i <= 5; i++ )
	{
		int tester_0 = 0, tester_1 = 0; // to assign the proper curve

		// generate a random point P_0 on e_0 with k * P_0 != O
		do
		{
			P_0 = e_0.random_point();
			multiply(Q_0, k, P_0); // Q_0 = k * P_0
		}
		while( Q_0.is_zero() );

		// generate a random point P_1 on e_1 with k * P_1 != O
		do
		{
			P_1 = e_1.random_point();
			multiply(Q_1, k, P_1); // Q_1 = k * P_1
		}
		while( Q_1.is_zero() );

		// compute r*Q_0 and r*Q_1
		multiply(P_0, r, Q_0); // P_0 = r * Q_0
		multiply(P_1, r, Q_1); // P_1 = r * Q_1      

		// if order of Q_0 is r set tester_0 = 1
		if(P_0.is_zero())
			tester_0 = 1;
		// if order of Q_1 is r set tester_1 to 1
		if(P_1.is_zero())
			tester_1 = 1;

		// assign proper curve of order r*k, if possible
		if( ( tester_0 == 1 ) && ( tester_1 == 0 ) )
		{
			G = Q_0; // base point of order r is Q_0
			a4 = a4_;
			a6 = a6_;
	   
			return( true );
		}

		if( (tester_0 == 0 ) && ( tester_1 == 1 ) )
		{
			G = Q_1; // base point of order r is Q_1
			a4 = a4_tw;
			a6 = a6_tw;
	   
			return( true );
		}
	}

	// after 5 failures we are convinced that none of the two
	// curves is of order r*k
	return( false );
}


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
