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
//  File    : prime_proof_ecpp_set_complex_precision.cc
//  Author  : Harald Baier (HB),
//	      small changes by Jochen Hechler (JH)
//      Changes : See CVS log
//
//==============================================================================================

#include       "LiDIA/prime_proof.h"


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif

// This Code is taken from the gec package.

//
// Determine the decimal precision "complex_precision" 
// for computation of class invariant.
// Currently we use the precision proposed in [Lay, Zimmer]
// for generation_mode != 4.
//



    void prime_proof::set_complex_precision()
    {

	// declare some bigfloats
	bigfloat precision_as_bigfloat; 
	bigfloat L( 0.0 ); 
	bigfloat tmp_f; 
	
  
	// compute sum of inverses of a of all quadratic forms (a,b,c)
	// of discriminant delta
	for( int i = 0; i < class_number; i++ )
	    add( L, L, inverse( bigfloat( class_group[ i ].get_a() ) ) );

	//
	// Compute L = Pi() * sqrt( delta ) / log(10) * sum ( 1/a )
	//
	multiply( L, L, Pi() );
	multiply( L, L, sqrt( bigfloat( abs( D ) ) ) );
	divide( L, L, log( bigfloat( 10.0 ) ) );

	//
	// Set precision according to the used class invariant
	// If generation_mode != 4, we make use of the precision
	// proposed by Lay and Zimmer.
	//
	// Otherwise, we set precision = alpha * L .
	// alpha is a static variable set in ecpp.cc
	//
	
	if( generation_mode != 4 ) 
	{  
	    add( precision_as_bigfloat, L, bigfloat( 5.0 ) );
	    divide( tmp_f, bigfloat( class_number ), bigfloat( 4.0 ) );
	    add( precision_as_bigfloat, precision_as_bigfloat, tmp_f );

	    if( generation_mode == 3 )
	    {
		divide( precision_as_bigfloat, 
			precision_as_bigfloat, bigfloat( 3.0 ) ); 
	    }
	}
	else
	    multiply( precision_as_bigfloat, L, alpha );

	// if precision is too low, we set it to 20
	if(precision_as_bigfloat < 5)precision_as_bigfloat = 20.0;

	//
	// assign precision_as_bigfloat to complex_precision
	//
	truncate( precision_as_bigfloat, precision_as_bigfloat );
	inc( precision_as_bigfloat );
	if( precision_as_bigfloat.longify( complex_precision ) )
	    lidia_error_handler( "gec_complex_multiplication",
				 "set_complex_precision: "
				 "cannot assign bigfloat to attribute "
				 "complex_precision");

    }

#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
