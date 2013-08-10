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
//  File    : find_good_prime_twist.cc
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
// the function find_good_prime for finding a twisted pair
//
bool
gec_E_and_twist_prime::find_good_prime()
{
	const int T = 2500; // definition: see PhD
  
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
	bigint order; // order of curve: order = q + 1 +- t

	//
	// temporary variables of type bigint and bigfloat
	//
	bigint tmp_i;
	bigfloat tmp_f;

	//
	// some literal constants for computations
	const bigint three( 3 );
	const bigint three_eight_five( 385 );
	const bigint seven_seven_zero( 770 );

	//
	// we ensure q \in [start_of_search, end_of_search]
	// currently we use 2^40 intervals (see bitlength_no_of_intervals,
	// defined in class gec_complex_multiplication)
	bigint start_of_search, end_of_search;

	bigint K( 1 ); // number of intervals
	shift_left(K, K, bitlength_no_of_intervals ); // K = 2^40
	
	bigint length_of_interval;
	divide( length_of_interval, lower_bound_q, K );

	bigint random_interval;
	random_interval.randomize( K );

	multiply( tmp_i, random_interval, length_of_interval );
	add( start_of_search, lower_bound_q, tmp_i );
	add( end_of_search, start_of_search, length_of_interval );

   //
   // upper bound of t: |t| <= 2\lceil sqrt(q) \rceil
   //
	bigint upper_bound_t;
	sqrt( tmp_f, bigfloat( start_of_search ) );
	shift_left( tmp_f, tmp_f, 1 ); // because 4*q = t^2 - Dy^2 
	truncate( tmp_f, tmp_f );
	tmp_f.bigintify( upper_bound_t );

	//
	// lower and upper bound of y for the chosen t
	//
	bigint lower_bound_y, upper_bound_y;

	bool in_interval;
	int addend_1, addend_2, factor, c;
	int j;

   // variables for statistics
	int no_of_initialization = 0;
	int c_r = 0; // to store cycles in which q left unchanged
   
   //
   // Variables to check if p or cardinality of curve 
   // is divisible by small primes.
   // We choose type long because of LiDIA function
   // long remainder(bigint &, long &).
   //
   
   // q mod 13, q mod 17, q mod 19, ...
	long pmod13, pmod17, pmod19,
		pmod23, pmod29, pmod31, pmod37, pmod41, 
		pmod43, pmod47, pmod53, pmod59,
		pmod61, pmod67, pmod71, pmod73, pmod79, pmod83,
		pmod89, pmod97;
   
	// t mod 5, t mod 7, ...
	long tmod5, tmod7, tmod11, tmod13, tmod17, tmod19,
		tmod23, tmod29, tmod31, tmod37, tmod41, 
		tmod43, tmod47, tmod53, tmod59,
		tmod61, tmod67, tmod71, tmod73, tmod79, tmod83,
		tmod89, tmod97;
   
	while( true )
	{
		in_interval = 1;
	   
		c_r = 0;
		no_of_initialization++;

		j = 0;

		// t random number in the given range
		t.randomize( upper_bound_t );
	   
		// y \approx sqrt( (t^2  - 4 * start_of_search ) / delta )
		square( tmp_i, t ); // tmp_i = t^2
		subtract( tmp_i, tmp_i, 4 * start_of_search );
		divide( tmp_f, bigfloat( tmp_i ), delta );
		sqrt( tmp_f, tmp_f );
		tmp_f.bigintify( y ); 

		// initialize t = 1 mod 770
		remainder( tmp_i, t, seven_seven_zero );
		subtract( tmp_i, seven_seven_zero, tmp_i );
		add( t, t, tmp_i );
		inc( t );
		
		// initialize y = 385 mod 770, y != 0 mod 3
		remainder( tmp_i, y, seven_seven_zero );
		subtract( tmp_i, seven_seven_zero, tmp_i );
		add( y, y, tmp_i );
		add( y, y, three_eight_five );

		remainder( tmp_i, y, three );
		if( tmp_i.is_zero() )
			add( y, y, seven_seven_zero );

		//
		// Set   q = ( t^2 - D * y^2 ) / 4
		//
		square( q, t );
		square( tmp_i, y );
		multiply( tmp_i, tmp_i, delta );
		subtract( q, q, tmp_i );

		shift_right( q, q, 2 );

		//
		// Ensure that q is in the interval
		//
		while( q < start_of_search)
	    {
			shift_left( tmp_i, t, 2 ); // tmp_i = 4 * t
			add( tmp_i, tmp_i, bigint( 4 ) ); // tmp_i = 4*t + 4
			add( q, q, tmp_i ); 
			add( t, t, bigint( 2 ) );
	    }

		//
		// Initialize variables to check if q 
		// or cardinality of curve is divisible by small primes 
		//
		pmod13 = remainder(q, 13);
		pmod17 = remainder(q, 17);
		pmod19 = remainder(q, 19);
		pmod23 = remainder(q, 23);
		pmod29 = remainder(q, 29);
		pmod31 = remainder(q, 31);
		pmod37 = remainder(q, 37);  
		pmod41 = remainder(q, 41);
		pmod43 = remainder(q, 43);
		pmod47 = remainder(q, 47);
		pmod53 = remainder(q, 53);
		pmod59 = remainder(q, 59);
		pmod61 = remainder(q, 61);
		pmod67 = remainder(q, 67);
		pmod71 = remainder(q, 71);
		pmod73 = remainder(q, 73);  
		pmod79 = remainder(q, 79);
		pmod83 = remainder(q, 83);
		pmod89 = remainder(q, 89);
		pmod97 = remainder(q, 97);

		tmod5 = remainder(t, 5);  
		tmod7 = remainder(t, 7);
		tmod11 = remainder(t, 11);
		tmod13 = remainder(t, 13);
		tmod17 = remainder(t, 17);
		tmod19 = remainder(t, 19);
		tmod23 = remainder(t, 23);
		tmod29 = remainder(t, 29);
		tmod31 = remainder(t, 31);
		tmod37 = remainder(t, 37);  
		tmod41 = remainder(t, 41);
		tmod43 = remainder(t, 43);
		tmod47 = remainder(t, 47);
		tmod53 = remainder(t, 53);
		tmod59 = remainder(t, 59);
		tmod61 = remainder(t, 61);
		tmod67 = remainder(t, 67);
		tmod71 = remainder(t, 71);
		tmod73 = remainder(t, 73);  
		tmod79 = remainder(t, 79);
		tmod83 = remainder(t, 83);
		tmod89 = remainder(t, 89);
		tmod97 = remainder(t, 97);

		do
	    {
			if( pmod13 && pmod17 && 
				pmod19 && pmod23 && pmod29 &&
				pmod31 && pmod37 && pmod41 && 
				pmod43 && pmod47 && pmod53 && pmod59 &&
				pmod61 && pmod67 && pmod71 && pmod73 && 
				pmod79 && pmod83 &&
				pmod89 && pmod97 &&
				((pmod13 + 1 - tmod13)%13) &&
				((pmod17 + 1 - tmod17)%17) &&
				((pmod19 + 1 - tmod19)%19) &&
				((pmod23 + 1 - tmod23)%23) &&
				((pmod29 + 1 - tmod29)%29) &&
				((pmod31 + 1 - tmod31)%31) &&
				((pmod37 + 1 - tmod37)%37) &&
				((pmod41 + 1 - tmod41)%41) &&
				((pmod43 + 1 - tmod43)%43) &&
				((pmod47 + 1 - tmod47)%47) &&
				((pmod53 + 1 - tmod53)%53) &&
				((pmod59 + 1 - tmod59)%59) &&
				((pmod61 + 1 - tmod61)%61) &&
				((pmod67 + 1 - tmod67)%67) &&
				((pmod71 + 1 - tmod71)%71) &&
				((pmod73 + 1 - tmod73)%73) &&
				((pmod79 + 1 - tmod79)%79) &&
				((pmod83 + 1 - tmod83)%83) &&
				((pmod89 + 1 - tmod89)%89) &&
				((pmod97 + 1 - tmod97)%97) &&
				((pmod13 + 1 + tmod13)%13) &&
				((pmod17 + 1 + tmod17)%17) &&
				((pmod19 + 1 + tmod19)%19) &&
				((pmod23 + 1 + tmod23)%23) &&
				((pmod29 + 1 + tmod29)%29) &&
				((pmod31 + 1 + tmod31)%31) &&
				((pmod37 + 1 + tmod37)%37) &&
				((pmod41 + 1 + tmod41)%41) &&
				((pmod43 + 1 + tmod43)%43) &&
				((pmod47 + 1 + tmod47)%47) &&
				((pmod53 + 1 + tmod53)%53) &&
				((pmod59 + 1 + tmod59)%59) &&
				((pmod61 + 1 + tmod61)%61) &&
				((pmod67 + 1 + tmod67)%67) &&
				((pmod71 + 1 + tmod71)%71) &&
				((pmod73 + 1 + tmod73)%73) &&
				((pmod79 + 1 + tmod79)%79) &&
				((pmod83 + 1 + tmod83)%83) &&
				((pmod89 + 1 + tmod89)%89) &&
				((pmod97 + 1 + tmod97)%97) )
			{
				
				// We only compute p and t if p is not 
				// divisible by small primes.
				// We have 
				// t(c_r) = t(0) + 2*c_r,
				// p(c_r) = p(0) + c_r*(t(0) + c_r).
				
				if( c_r )
				{
					add( tmp_i, t, bigint( c_r ) );
					multiply( tmp_i, tmp_i, bigint( c_r ) );
					add( q, q, tmp_i ); // adapt q
					add( t, t, bigint( 2 * c_r ) ); // adapt t
					
					c_r = 0;
				}
				
				if( q > end_of_search )
					in_interval = 0;
				
				if( in_interval && ( efficient_curve_parameters == false || 
									 q.bit( 1 ) == 1 ) &&
					is_prime( q, nr_of_prob_prime_tests ) )
				{
					// compute r = q + 1 - t

					subtract( r, q, t ); // r = q - t
					inc( r );

						// compute r_tw = q + 1 + t
					add( r_tw, q, t ); // r_tw = q + t
					inc( r_tw );

					p.assign( q );

					if( is_prime( r, nr_of_prob_prime_tests ) && 
						is_prime( r_tw, nr_of_prob_prime_tests ) && 
						is_cryptographically_strong() )
					{
						return( true );
					}
				}
			}

			if( j >= T )
				break;
			
			// Adapt the remainder of q and t modulo small primes.
			// We have to change both pmod and tmod because of the 
			// inequation
			// p(c_r) = p(0) + 4*c_r*(t(0) + c_r) !=
            // p(0) + 4*c_r*(t(0) + 1)
			// in case of c_r > 1.
			// 
			// According to the current state of tmod5, tmod7 and tmod11
			// we can skip several rounds.

			if( tmod5 == 1 && tmod11 == 3 && ( tmod7 == 1 || tmod7 == 3 ) )
				c = 19;

			else if( tmod5 == 1 && 
					 ( ( tmod7 == 1 && ( tmod11==1 || tmod11==10 ) ) ||
					   ( tmod7 == 6 && ( tmod11==4 || tmod11==6 ) ) ||
					   ( tmod7 == 3 && tmod11 == 1 ) ||
					   ( tmod7 == 4 && tmod11 == 4 ) ) )
				c = 14;

			else if( tmod5 == 4 && tmod7 == 3 && 
					 ( tmod11 == 1 || tmod11 == 10 ) )
				c = 11;

			else if( tmod7 == 4 && 
					 ( ( tmod5 == 1 && tmod11 == 6 ) || 
					   ( tmod5 == 4 && tmod11 == 7 ) ) )
				c = 10;

			else if( tmod5 == 1 && ( tmod7 == 4 || tmod7 == 6 ) )
				c = 9;

			else if( tmod5 == 4 && 
					 ( ( tmod7 == 3 && tmod11 == 3 ) || 
					   ( tmod7 == 6 && tmod11 == 7 ) ) )
				c = 6;

			else if( ( tmod5 == 1 && 
					   ( tmod7 == 1 || ( tmod7 == 3 && tmod11 == 5 ) ) ) || 
					 ( tmod5 == 4 && 
					   ( tmod7 == 3 || ( tmod7 == 1 && tmod11 == 7 ) ) ) ) 
				c = 5;

			else if( tmod5 == 1 && tmod7 == 3 )
				c = 4;

			else
				c = 1;

			addend_2 = 2 * c;
			addend_1 = c * c;
			factor = c;
			c_r += c;

			pmod13 = (pmod13 + factor * tmod13 + addend_1)%13;
			pmod17 = (pmod17 + factor * tmod17 + addend_1)%17;
			pmod19 = (pmod19 + factor * tmod19 + addend_1)%19;
			pmod23 = (pmod23 + factor * tmod23 + addend_1)%23;
			pmod29 = (pmod29 + factor * tmod29 + addend_1)%29;
			pmod31 = (pmod31 + factor * tmod31 + addend_1)%31;
			pmod37 = (pmod37 + factor * tmod37 + addend_1)%37;
			pmod41 = (pmod41 + factor * tmod41 + addend_1)%41;
			pmod43 = (pmod43 + factor * tmod43 + addend_1)%43;
			pmod47 = (pmod47 + factor * tmod47 + addend_1)%47;
			pmod53 = (pmod53 + factor * tmod53 + addend_1)%53;
			pmod59 = (pmod59 + factor * tmod59 + addend_1)%59;
			pmod61 = (pmod61 + factor * tmod61 + addend_1)%61;
			pmod67 = (pmod67 + factor * tmod67 + addend_1)%67;
			pmod71 = (pmod71 + factor * tmod71 + addend_1)%71;
			pmod73 = (pmod73 + factor * tmod73 + addend_1)%73;
			pmod79 = (pmod79 + factor * tmod79 + addend_1)%79;
			pmod83 = (pmod83 + factor * tmod83 + addend_1)%83;
			pmod89 = (pmod89 + factor * tmod89 + addend_1)%89;
			pmod97 = (pmod97 + factor * tmod97 + addend_1)%97;
				
			tmod5 = (tmod5 + addend_2)%5;
			tmod7 = (tmod7 + addend_2)%7;
			tmod11 = (tmod11 + addend_2)%11;
			tmod13 = (tmod13 + addend_2)%13;
			tmod17 = (tmod17 + addend_2)%17;
			tmod19 = (tmod19 + addend_2)%19;
			tmod23 = (tmod23 + addend_2)%23;
			tmod29 = (tmod29 + addend_2)%29;
			tmod31 = (tmod31 + addend_2)%31;
			tmod37 = (tmod37 + addend_2)%37;
			tmod41 = (tmod41 + addend_2)%41;
			tmod43 = (tmod43 + addend_2)%43;
			tmod47 = (tmod47 + addend_2)%47;
			tmod53 = (tmod53 + addend_2)%53;	     
			tmod59 = (tmod59 + addend_2)%59;
			tmod61 = (tmod61 + addend_2)%61;
			tmod67 = (tmod67 + addend_2)%67;
			tmod71 = (tmod71 + addend_2)%71;
			tmod73 = (tmod73 + addend_2)%73;
			tmod79 = (tmod79 + addend_2)%79;
			tmod83 = (tmod83 + addend_2)%83;
			tmod89 = (tmod89 + addend_2)%89;
			tmod97 = (tmod97 + addend_2)%97;

			j++;
			
		} while( in_interval );
	}
}


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
