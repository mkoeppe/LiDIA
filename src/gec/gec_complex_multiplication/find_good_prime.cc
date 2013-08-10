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
//  File    : find_good_prime.cc
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


bool
gec_complex_multiplication::find_good_prime()
{
	if( VERBOSE )
	{
		std::cout << "Fixed Discriminant Approach: Searching the prime..."
				  << std::endl;
	}

	// compute |delta| for distinguishing three cases:
	abs_delta.assign( delta );
	abs_delta.negate();

	// we only use efficient approach, if lower_bound_bitlength_r >=100
	// otherwise, we use simple approach
	bool use_efficient_find_prime = true;
	if( lower_bound_bitlength_r < 100)
		use_efficient_find_prime = false;

	if( use_efficient_find_prime )
	{
		// delta = 1 mod 8  <=>  |delta| = 7 mod 8
		// In addition, the efficient implementation assumes upper_bound_k = 4.
		if( (abs_delta.bit(0) == 1) && (abs_delta.bit(1) == 1) && 
			(abs_delta.bit(2) == 1) && upper_bound_k == 4 )
		{
			if( find_good_prime_1_mod_8() )
				return( true ) ;
			else
				lidia_error_handler( "gec_complex_multiplication", 
									 "find_good_prime: unable to find prime. " );
		}

		// delta = 5 mod 8  <=>  |delta| = 3 mod 8
		// In addition, the efficient implementation assumes upper_bound_k = 4.
		if( (abs_delta.bit(0) == 1) && (abs_delta.bit(1) == 1) && 
			(abs_delta.bit(2) == 0) && upper_bound_k == 1 )
		{
			if( find_good_prime_5_mod_8() )
				return( true ) ;
			else
				lidia_error_handler( "gec_complex_multiplication", 
									 "find_good_prime: unable to find prime. " );
		}

		// delta = 0 mod 4  <=>  |delta| = 0 mod 4
		// In addition, the efficient implementation assumes upper_bound_k <= 4.
		if( abs_delta.bit( 0 ) == 0 && abs_delta.bit( 1 ) == 0 &&
			upper_bound_k <= 4 )
		{
			if( find_good_prime_0_mod_4() )
				return( true ) ;
			else
				lidia_error_handler( "gec_complex_multiplication", 
									 "find_good_prime: unable to find prime. " );
		}
	}

/*****************************************************

Otherwise, we use some non-optimized and slow approach:

see algorithm generatePrimeRandomTrace (algorithm 4.3) from PhD

*****************************************************/
   
	//
	// lower_bound_q = 2^(lower_bound_bitlength_r - 1) * upper_bound_k
	// upper_bound_q = 2 * lower_bound_q
	//
	bigint lower_bound_q(1), upper_bound_q(1);
	shift_left( lower_bound_q, lower_bound_q, lower_bound_bitlength_r - 1 );
	multiply( lower_bound_q, lower_bound_q, upper_bound_k );
	shift_left( upper_bound_q, lower_bound_q, 1 );

	//
	// trace t and y for representation 4q=t^2-delta*y^2
	//
	bigint t;
	bigint y;
	bigint order; // order of curve: order = q + 1 - t

	//
	// temporary variables of type bigint and bigfloat
	//
	bigint tmp_i;
	bigfloat tmp_f;

	//
	// upper bound of t: |t| <= \lceil 2*sqrt(q) \rceil
	//
	bigint upper_bound_t;
	sqrt( tmp_f, bigfloat( upper_bound_q ) );
	shift_left( tmp_f, tmp_f, 1 );
	truncate( tmp_f, tmp_f );
	tmp_f.bigintify( upper_bound_t );

	//
	// lower and upper bound of y for the chosen t
	//
	bigint lower_bound_y, upper_bound_y;

	//
	// we set lower_bound_q <- 4*lower_bound_q and
	//        upper_bound_q <- 4*upper_bound_q
	// because of the factor 4 in 4q=t^2-delta*y^2
	shift_left( lower_bound_q, lower_bound_q, 2);
	shift_left( upper_bound_q, upper_bound_q, 2);
	       
	while( true )
	{
		t = randomize( upper_bound_t );

		//
		// compute the lower bound of y
		//
		tmp_i = lower_bound_q - t*t;
		if( tmp_i <= 0 )
			lower_bound_y = 0;
		else
		{
			divide( tmp_f, bigfloat( tmp_i ), bigfloat( -delta ) );
			sqrt( tmp_f, tmp_f );
			truncate( tmp_f, tmp_f );
			tmp_f += 1;
			tmp_f.bigintify( lower_bound_y );
		}
	
		//
		// compute the upper bound of y
		//
		tmp_i = upper_bound_q - t*t;
		divide( tmp_f, bigfloat( tmp_i ), bigfloat( -delta ) );
		sqrt( tmp_f, tmp_f );
		truncate( tmp_f, tmp_f );
		tmp_f.bigintify( upper_bound_y );

		if( upper_bound_y != lower_bound_y )
			y = randomize( upper_bound_y - lower_bound_y );
		else
			y.assign_zero();

		y += lower_bound_y;

		if( delta_case == 0 )
		{
			if( t.bit(0) )
				t--;
		}      

		if( delta_case == 1 )
		{
			if( t.bit(0) != y.bit(0) )
				t--;
		}

		q = t*t - delta*y*y;

		shift_right( q, q, 2 );

		if( is_prime( q, nr_of_prob_prime_tests ) )
		{
			if( VERBOSE )
				std::cout << "Testing for the prime q = " << q << std::endl;
	  
			p = q;	
    
			order = q + 1 - t;
			if( is_cryptographically_strong( order ) )
				return( true );    

			order = q + 1 + t;
			if( is_cryptographically_strong( order ) )
				return( true );
		}
	}
}


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif









