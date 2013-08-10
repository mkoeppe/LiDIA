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
//  File    : find_good_prime_0_mod_4.cc
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
// This function implements algorithm findPrime0Mod4 from PhD.
// For details of the implementation, we refer to the corresponding
// chapter of PhD.
//
bool
gec_complex_multiplication::find_good_prime_0_mod_4()
{
	// number of subsequent tested pairs (t,y);
	// if not successful after T tries, we initialize a new pair.
	const int T = 2000;

	int t_bit_1; // bit of 2 in binary expansion of t

	//
	// handle the 2 different cases:
	//
	//  ( t mod 2, t_bit_1 )   |   delta mod 16
	// ----------------------------------------
	//                         |
	//       (  0, 1 )         |      0, 8
	//                         |
	//       (  0, 0 )         |      4, 12

	
	abs_delta.assign( delta );
	abs_delta.negate();

	if( abs_delta.bit( 2 ) == 0 ) // delta = 0, 8 mod 16
		t_bit_1 = 1;
	else // delta = 4, 12 mod 16
		t_bit_1 = 0;

	//
	// temporary variables of type bigint and bigfloat
	//
	bigint tmp_i;
	bigfloat tmp_f;
	bigint one_zero_five( 105 );

	//
	// if delta mod 3 = 2 we have more chances
	//
	bool delta_mod_3_2 = false;
	remainder( tmp_i, delta, bigint( 3 ) );	
	if( tmp_i == -1 )
		delta_mod_3_2 = true;
	//
	// lower_bound_q = 2^(lower_bound_bitlength_r - 1) * upper_bound_k
	// upper_bound_q = 2 * lower_bound_q
	//
	bigint lower_bound_q(1), upper_bound_q(1);
	shift_left( lower_bound_q, lower_bound_q, lower_bound_bitlength_r - 1 );
	multiply( lower_bound_q, lower_bound_q, upper_bound_k );
	shift_left( upper_bound_q, lower_bound_q, 1 );

	//
	// trace t and y for representation q=t^2-delta*y^2
	// 
	bigint t;
	bigint y;
	bigint order; // order of curve: order = q + 1 +- t

	//
	// we ensure q \in [start_of_search, end_of_search]
	// currently we use 2^40 intervals
	bigint start_of_search, end_of_search;

	bigint K( 1 ); // number of intervals
	shift_left(K, K, bitlength_no_of_intervals ); // =2^40
	
	bigint length_of_interval( 1 );
	divide( length_of_interval, lower_bound_q, K );

	bigint random_interval;
	random_interval.randomize( K );

	multiply( tmp_i, random_interval, length_of_interval );
	add( start_of_search, lower_bound_q, tmp_i );
	add( end_of_search, start_of_search, length_of_interval );

	//
	// upper bound of t: |t| <= 2\lceil sqrt(start_of_search) \rceil
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

	//
	// some more variables
	//
	bool in_interval;
	int addend_1, addend_2, factor, c; // for updating p mod m, t mod m
	bigint cycles_times_2; // 2 * cycles
	int j; // counts how often t is changed
	// to keep track whether we pass through p+1-t and p is prime
	bool minus_good, minus_prime_true, test_q;

	// variables for statistics
	int no_of_initialization = 0;
	int c_r = 0; // to store cycles in which q left unchanged

	//
	// Variables to check if p or cardinality of curve 
	// is divisible by small primes.
	// We choose type long because of LiDIA function
	// long remainder(bigint &, long &).
	//
	
	// q mod 11, q mod 13, ...
	// we do not need q mod 3, q mod 5, and q mod 7 as
	// all these residues are nonzero
	long pmod11, pmod13, pmod17, pmod19,
		pmod23, pmod29, pmod31, pmod37, pmod41, 
		pmod43, pmod47, pmod53, pmod59,
		pmod61, pmod67, pmod71, pmod73, pmod79, pmod83,
		pmod89, pmod97;
   
	// t mod 3, t mod 5, ...
	long tmod3, tmod5, tmod7, tmod11, tmod13, tmod17, tmod19,
		tmod23, tmod29, tmod31, tmod37, tmod41, 
		tmod43, tmod47, tmod53, tmod59,
		tmod61, tmod67, tmod71, tmod73, tmod79, tmod83,
		tmod89, tmod97;
   
	while( true )
	{
		// initialize some variables
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

		//
		// initialize t as described in PhD
		//

		// initialize t = 1 mod 105
		remainder( tmp_i, t, one_zero_five );
		subtract( tmp_i, one_zero_five, tmp_i );
		add( t, t, tmp_i );
		inc( t );

		// ensure t mod 4 as in PhD
		while( t.bit( 0 ) == 1 || t.bit( 1 ) != t_bit_1 )
			add( t, t, one_zero_five ); // 105 = 1 mod 4

		//
		// initialize y as specified in PhD
		//
		
		if( ! delta_mod_3_2 )
		{
			// we ensure y mod 210 = 105
			remainder( tmp_i, y, bigint( 210 ) );
			subtract( tmp_i, bigint( 210 ), tmp_i );
			add( y, y, tmp_i );
			add( y, y, one_zero_five );
		}
		else
		{
			// we ensure y mod 70 = 35, y mod 3 != 0
			remainder( tmp_i, y, bigint( 70 ) );
			subtract( tmp_i, bigint( 70 ), tmp_i );
			add( y, y, tmp_i );
			add( y, y, bigint( 35 ) ); // y mod 70 = 35

			remainder( tmp_i, y, bigint( 3 ) );

			if( tmp_i == 0 )
				add( y, y, bigint( 70 ) ); // y mod 3 = 1
		}


		//
		// Set   q = ( t^2 - D * y^2 ) / 4
		//
		square( q, t );
		square( tmp_i, y );
		multiply( tmp_i, tmp_i, delta );
		subtract( q, q, tmp_i );

		shift_right( q, q, 2 );

		//
		// Initialize variables to check if q 
		// or cardinality of curve is divisible by small primes 
		//
		pmod11 = remainder(q, 11);
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

		tmod3 = remainder(t, 3);
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
			if( pmod11 && 
				pmod13 && pmod17 && 
				pmod19 && pmod23 && pmod29 &&
				pmod31 && pmod37 && pmod41 && 
				pmod43 && pmod47 && pmod53 && pmod59 &&
				pmod61 && pmod67 && pmod71 && pmod73 && 
				pmod79 && pmod83 &&
				pmod89 && pmod97)
			{
				minus_good = minus_prime_true = false;
				
				if( ( tmod3 == 1 || delta_mod_3_2 ) &&
					((pmod11 + 1 - tmod11)%11) &&
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
					((pmod97 + 1 - tmod97)%97) )
				{
					minus_good = true;

					// We only compute p and t if p is not 
					// divisible by small primes.
					// We have 
					// t(c_r) = t(0) + 4*c_r,
					// p(c_r) = p(0) + 2*c_r*t(0) + 4*c_r*c_r
					//           = p(0) + 2*c_r*( t(0) + 2*c_r)
					
					cycles_times_2.assign( bigint( 2 * c_r ) );
				
					// adapt q
					add( tmp_i, t, cycles_times_2 );
					multiply( tmp_i, tmp_i, cycles_times_2 );
					add( q, q, tmp_i );
					
					// adapt t
					shift_left( cycles_times_2, cycles_times_2, 1 );
					add( t, t, cycles_times_2 );
					
					c_r = 0;
					
					if( q > end_of_search )
						in_interval = 0;
					
					if( in_interval && is_prime( q ) )
					{
						minus_prime_true = true;

						// compute order = q + 1 - t
						subtract( order, q, t ); // order = q - t
						inc( order );

						p.assign( q );

						if( is_cryptographically_strong( order ) )
							return( true );
					}
				}

				if( ( tmod3 == 2 || delta_mod_3_2 ) &&
					((pmod11 + 1 + tmod11)%11) &&
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
					if( c_r != 0 )
					{
						cycles_times_2.assign( bigint( 2 * c_r ) );
						
						// adapt q
						add( tmp_i, t, cycles_times_2 );
						multiply( tmp_i, tmp_i, cycles_times_2 );
						add( q, q, tmp_i );
						
						// adapt t
						shift_left( cycles_times_2, cycles_times_2, 1 );
						add( t, t, cycles_times_2 );
						
						c_r = 0;
					}

					if( q > end_of_search )
						in_interval = 0;

					if( minus_good == true && minus_prime_true == false )
						test_q = false;
					else
						test_q = true;

					if( in_interval && test_q && is_prime( q ) )
					{
						// compute order = q + 1 + t
						add( order, q, t ); // order = q + t
						inc( order );
						
						p.assign( q );
						
						if( is_cryptographically_strong( order ) )
							return( true );
					}
				}
			}

			if( j > T )
				break;

			// Adapt the remainder of q and t modulo small primes.
			// We have to change both pmod and tmod because of the 
			// inequation
			// p(c_r) = p(0) + 4*c_r*(t(0) + c_r) !=
            // p(0) + 4*c_r*(t(0) + 1)
			// in case of c_r > 1.
			// 
			// p(c_r) = p(0) + 2*c_r*t + 4*c_r*c_r
			// t(c_r) = t(0) + 4*c_r
			//
			// factor = 2*c_r; addend_1 = 4*c_r*c_r
			// addend_2 = 4*c_r

			if( tmod3 == 2 && 
				( ( tmod5 == 1 && tmod7 == 1 ) ||
				  ( tmod5 == 2 && tmod7 == 2 ) ) )
				
				c = 8;

			else if( tmod3 == 1 && 
					 ( ( tmod5 == 3 && tmod7 == 3 ) ||
					   ( tmod5 == 4 && tmod7 == 4 ) ) )
				
				c = 7;

			else if( tmod3 == 2 && 
					 ( ( tmod5 == 1 && tmod7 == 6 ) ||
					   ( tmod5 == 4 && tmod7 == 2 ) ) )
				
				c = 5;

			else if( tmod3 == 1 && 
					 ( ( tmod5 == 1 && tmod7 == 4 ) ||
					   ( tmod5 == 3 && tmod7 == 1 ) ) )
				
				c = 4;

			else if( ( tmod3 == 1 && ( tmod5 == 1 || tmod5 == 4 ) ) || 
					 ( tmod3 == 2 && ( tmod5 == 2 || tmod5 == 4 ) ) )
				
				c = 3;

			else if( tmod3 == 2 && tmod5 == 1 ) 

				c = 2;

			else
				
				c = 1;
				

			addend_2 = 4 * c;
			addend_1 = addend_2 * c;
			factor = 2 * c;
			c_r += c;

			pmod11 = (pmod11 + factor * tmod11 + addend_1)%11;
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
				
			tmod3 = (tmod3 + addend_2)%3;
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
