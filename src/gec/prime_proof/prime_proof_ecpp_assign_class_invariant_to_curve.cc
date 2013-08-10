// -*- C++ -*-
//==============================================================================================
//
//      This file is part of LiDIA --- a library for computational number theory
//
//      Copyright (c) 1994--2003 the LiDIA Group.  All rights reserved.
//
//      See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//----------------------------------------------------------------------------------------------
//
//
//  File    : prime_proof_ecpp_assign_class_invariant_to_curve.cc
//  Author  : Jochen Hechler (JH)
//      Changes : See CVS log
//
//==============================================================================================


#include       	"LiDIA/prime_proof.h"

#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



    //
    // the function assign_class_invariant_to_curve computes
    // two twisted curves of j-invariant invariant. 
    //

    bool prime_proof::assign_class_invariant_to_curve( const gf_element & invariant )
    { 
	gf_element tmp_gf_1( F_q ), tmp_gf_2( F_q ), tmp_gf_3( F_q );
	gf_element j( F_q );

	//
	// Compute a j-invariant j 
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

        else if( generation_mode == 3 )
	    power( j, invariant, bigint( 3 ) );
        else
	    j = invariant;


        tmp_gf_1.assign( bigint( 1728 ) );


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

	// Compute a quadratic nonresidue by trial-and-error
	bigint e_; // exponent e_



        gf_element q_non_res( F_q );
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



	// we allow at most 5 tests; otherwise none of the curves is of order
        point<gf_element> P_0, P_1, Q_0, Q_1;

	bigint div_order;
	divide(div_order,order,prime_factor);
	bool b1,b2,b3;
	gcd = 1;
	for( int i = 1; i <= 5; i++ )
	{
	    int tester_0 = 0, tester_1 = 0; // to assign the proper curve
	    // generate a random point P_0 on e_0 with P_0 != O
	    do
	    {
		P_0 = e_0.random_point();
		multiply(Q_0, order, P_0); // Q_0 = k * P_0
	    }
	    while( Q_0.is_zero() );


	    // generate a random point P_1 on e_1 with P_1 != O
	    do
	    {
		P_1 = e_1.random_point();
		multiply(Q_1, order, P_1); // Q_1 = k * P_1
	    }
	    while( Q_1.is_zero() );


	    // compute order*P_0 and order*P_1
	    multiply(P_0, order, Q_0); // P_0 = r * Q_0
	    multiply(P_1, order, Q_1); // P_1 = r * Q_1      

	    // if order of Q_0 is r set tester_0 = 1
	    if(P_0.is_zero())
		tester_0 = 1;
	    // if order of Q_1 is r set tester_1 to 1
	    if(P_1.is_zero())
		tester_1 = 1;


	    if( ( tester_0 == 1 ) && ( tester_1 == 0 ) )
	    {
		do
		{
		    P_0 = e_0.random_point();	
		    multiply(Q_0, order, P_0);
		    b1= Q_0.is_zero();
		    multiply(Q_0, div_order, P_0);
		    b2= Q_0.is_zero();
		    b3= Q_0.on_curve();
		}
		while( (!b1 || b2) && b3 );
		if(verbose && b1 && !b2 && b3)std::cout<<"ecpp:Test success"<<std::endl;
		if(verbose && (!b1 || b2) )std::cout<<"ecpp: Test was not successfull"<<std::endl;
		if(!b3) gcd = 2;
		if(!b3 || !b1)return(false);
		cert_vector[0]=D; // type of test an discriminant
		cert_vector[1]=n; // prime to test
		cert_vector[2]=order; // order of elliptic curve
		cert_vector[3]=prime_factor; //primefactor of order
		cert_vector[4]=a4_.lift_to_Z(); // y^2=x^3+x*a4_bm+a6_bm
		cert_vector[5]=a6_.lift_to_Z();
		cert_vector[6]=P_0.get_x().lift_to_Z();
		cert_vector[7]=P_0.get_y().lift_to_Z();
		return( true );
	    }

	    if( (tester_0 == 0 ) && ( tester_1 == 1 ) )
	    {
		do
		{
		    P_1 = e_1.random_point();
		    multiply(Q_1, order, P_1);
		    b1= Q_1.is_zero();
		    multiply(Q_1, div_order, P_1);
		    b2= Q_1.is_zero();
		    b3=Q_1.on_curve();
		}
		while( (!b1 || b2) && b3 );
		if(verbose && b1 && !b2 && b3)std::cout<<"ecpp:Test success"<<std::endl;
		if(verbose && (!b1 || b2))std::cout<<"ecpp: Test was not successfull"<<std::endl;
		if(!b3) gcd = 2;
		if(!b3 || !b1)return(false);
		cert_vector[0]=D; // type of test an discriminant
		cert_vector[1]=n; // prime to test
		cert_vector[2]=order; // order of elliptic curve
		cert_vector[3]=prime_factor; //primefactor of order
		cert_vector[4]=a4_tw.lift_to_Z(); // y^2=x^3+x*a4_bm+a6_bm
		cert_vector[5]=a6_tw.lift_to_Z();
		cert_vector[6]=P_1.get_x().lift_to_Z();
		cert_vector[7]=P_1.get_y().lift_to_Z();
		return( true );
	    }
	}
	if(verbose)std::cout<<"ecpp: was not successfull 1"<<std::endl;
	return( false );
    }


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
