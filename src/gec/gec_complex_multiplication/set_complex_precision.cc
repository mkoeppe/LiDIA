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
//  File    : set_complex_precision.cc
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
// Determine the decimal precision "complex_precision" 
// for computation of class invariant.
// Currently we use the precision proposed in [Lay, Zimmer]
// for generation_mode != 4.
//
// If comp_prec > 0, the precision is set to comp_prec.
// The user is responsible that the chosen precision yields a correct result.
//
void
gec_complex_multiplication::set_complex_precision( long comp_prec )
{
	// if comp_prec > 0, use this precision
	if( comp_prec > 0 )
	{
		complex_precision = comp_prec;
		return;
	}

	// if either delta = 0 or h = 0, invoke the lidia_error_handler
	if( delta == 0 || h == 0 )
		lidia_error_handler("gec_complex_multiplication",
							"set_complex_precision(): either delta or h is not initialized");
	
	// declare some bigfloats
	bigfloat precision_as_bigfloat; // the precision as a bigfloat
	bigfloat L( 0.0 ); // as described in PhD
	bigfloat tmp_f; // temporary variable of type bigfloat
  
	// compute sum of inverses of a of all quadratic forms (a,b,c)
	// of discriminant delta
	for( int i = 0; i < h; i++ )
		add( L, L, inverse( bigfloat( class_group[ i ].get_a() ) ) );

	//
	// Compute L = Pi() * sqrt( delta ) / log(10) * sum ( 1/a )
	//
	multiply( L, L, Pi() );
	multiply( L, L, sqrt( bigfloat( abs( delta ) ) ) );
	divide( L, L, log( bigfloat( 10.0 ) ) );

	//
	// Set precision according to the used class invariant
	// If generation_mode != 4, we make use of the precision
	// proposed by Lay and Zimmer.
	//
	// Otherwise, we set precision = alpha * L .
	// alpha is a static variable set in gec_complex_multiplication.cc
	//
	if( generation_mode != 4 ) 
	{  
		add( precision_as_bigfloat, L, bigfloat( 5.0 ) );
		divide( tmp_f, bigfloat( h ), bigfloat( 4.0 ) );
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
	if( precision_as_bigfloat < 5 )
		precision_as_bigfloat = 20.0;

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

	if( VERBOSE )
		std::cout << std::endl << "The floating point precision is " 
				  << complex_precision << std::endl;
}


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
