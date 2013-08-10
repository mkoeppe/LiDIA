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
//  File    : find_good_prime_1_mod_8.cc
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
// This function implements algorithm findPrime1Mod8 from PhD.
// For details of the implementation, we refer to the corresponding
// chapter of PhD.
//
bool
gec_complex_multiplication::find_good_prime_1_mod_8()
{
	// number of subsequent tested pairs (t,y);
	// if not successful after T tries, we initialize a new pair.
	const int T = 3000;

	//
	// lower_bound_q = 2^(lower_bound_bitlength_r - 1) * upper_bound_k
	// upper_bound_q = 2 * lower_bound_q
	//
	bigint lower_bound_q(1), upper_bound_q(1);
	shift_left( lower_bound_q, lower_bound_q, lower_bound_bitlength_r - 1 );
	multiply( lower_bound_q, lower_bound_q, upper_bound_k );
	shift_left( upper_bound_q, lower_bound_q, 1 );

	//
	// trace t_1_2 and y_1_2 for representation q=t_1_2^2-delta*y_1_2^2
	// Remark that we only deal with half the trace and y!!!!
	//
	bigint t_1_2;
	bigint y_1_2;
	bigint order; // order of curve: order = q + 1 +- 2 * t_1_2

	//
	// temporary variables of type bigint and bigfloat
	//
	bigint tmp_i;
	bigfloat tmp_f;
	bigint one_zero_five( 105 );
	
	//
	// we ensure q \in [start_of_search, end_of_search]
	// currently we use 2^40 intervals
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
	// upper bound of t_1_2: |t_1_2| <= \lceil sqrt(q) \rceil
	//
	bigint upper_bound_t;
	sqrt( tmp_f, bigfloat( start_of_search ) );
	truncate( tmp_f, tmp_f );
	tmp_f.bigintify( upper_bound_t );

	//
	// lower and upper bound of y_1_2 for the chosen t_1_2
	//
	bigint lower_bound_y, upper_bound_y;

	//
	// variables for sieving modifier
	//
	int addend_1, addend_2, factor, c;
	bool in_interval;
	int j;
	// to keep track whether we pass through p+1-t and p is prime
	bool b_m, b_p;

	// variables for statistics
	int no_of_initialization = 0;
	int c_r = 0; // to store cycles in which q left unchanged
   
	//
	// Variables to check if p or cardinality of curve 
	// is divisible by small primes.
	// We choose type long because of LiDIA function
	// long remainder(bigint &, long &).
	//
	
	// q mod 3, q mod 5, ...
	long pmod3, pmod5, pmod7, pmod11, pmod13, pmod17, pmod19,
		pmod23, pmod29, pmod31, pmod37, pmod41, 
		pmod43, pmod47, pmod53, pmod59,
		pmod61, pmod67, pmod71, pmod73, pmod79, pmod83,
		pmod89, pmod97;
	
	// t mod 3, t mod 5, ...
	long tmod3, tmod5, tmod7, tmod11, tmod13, tmod17, tmod19,
		tmod23, tmod29, tmod31, tmod37, tmod41, 
		tmod43, tmod47, tmod53, tmod59,
		tmod61, tmod67, tmod71, tmod73, tmod79, tmod83,
		tmod89, tmod97;//t mod 3,...
   
	while( true )
	{
		in_interval = 1;

		j = 0;
		c_r = 0;
		no_of_initialization++;

		// t_1_2 random number in the given range
		t_1_2.randomize( upper_bound_t );
	   
		// y_1_2 \approx sqrt( (t_1_2^2  - start_of_search ) / delta )
		square( tmp_i, t_1_2 ); // tmp_i = t_1_2^2
		subtract( tmp_i, tmp_i, start_of_search );
		divide( tmp_f, bigfloat( tmp_i ), delta );
		sqrt( tmp_f, tmp_f );
		tmp_f.bigintify( y_1_2 ); 

		//
		// we have to ensure t_1_2 = 1 mod 3*5*7 
		//
		remainder( tmp_i, t_1_2, one_zero_five );
		subtract( tmp_i, one_zero_five, tmp_i );
		add( t_1_2, t_1_2, tmp_i );
		inc( t_1_2 );

		// initialize t_1_2 = 1 mod 2 and 1 mod 105
		if( t_1_2.bit( 0 ) == 0 )
			add( t_1_2, t_1_2, one_zero_five );

		//
		// we have to ensure that y_1_2 is divisible by 3*5*7=105
		//
		remainder( tmp_i, y_1_2, one_zero_five );
		subtract( tmp_i, one_zero_five, tmp_i );
		add( y_1_2, y_1_2, tmp_i );
		
		// initialize y_1_2 = 0 mod 2
		if( y_1_2.bit( 0 ) == 1 )
			add( y_1_2, y_1_2, one_zero_five );
		
		// initialize (t_1_2 mod 4, y_1_2 mod 4) = (1,0), (3,2)
		if( ( t_1_2.bit( 1 ) == 0 && y_1_2.bit( 1 ) == 1 ) ||
			( t_1_2.bit( 1 ) == 1 && y_1_2.bit( 1 ) == 0 ) )
			add( y_1_2, y_1_2, bigint( 210 ) );
		
		//
		// Set   q = t_1_2^2 - D * y_1_2^2
		//
		square( q, t_1_2 );
		square( tmp_i, y_1_2 );
		multiply( tmp_i, tmp_i, delta );
		subtract( q, q, tmp_i );
	  

		//
		// Initialize variables to check if q 
		// or cardinality of curve is divisible by small primes 
		//
		pmod3 = remainder(q, 3);
		pmod5 = remainder(q, 5);  
		pmod7 = remainder(q, 7);
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

		tmod3 = remainder(t_1_2, 3);
		tmod5 = remainder(t_1_2, 5);  
		tmod7 = remainder(t_1_2, 7);
		tmod11 = remainder(t_1_2, 11);
		tmod13 = remainder(t_1_2, 13);
		tmod17 = remainder(t_1_2, 17);
		tmod19 = remainder(t_1_2, 19);
		tmod23 = remainder(t_1_2, 23);
		tmod29 = remainder(t_1_2, 29);
		tmod31 = remainder(t_1_2, 31);
		tmod37 = remainder(t_1_2, 37);  
		tmod41 = remainder(t_1_2, 41);
		tmod43 = remainder(t_1_2, 43);
		tmod47 = remainder(t_1_2, 47);
		tmod53 = remainder(t_1_2, 53);
		tmod59 = remainder(t_1_2, 59);
		tmod61 = remainder(t_1_2, 61);
		tmod67 = remainder(t_1_2, 67);
		tmod71 = remainder(t_1_2, 71);
		tmod73 = remainder(t_1_2, 73);  
		tmod79 = remainder(t_1_2, 79);
		tmod83 = remainder(t_1_2, 83);
		tmod89 = remainder(t_1_2, 89);
		tmod97 = remainder(t_1_2, 97);

		do
	    {
			if( pmod3 && pmod5 && 
				pmod7 && pmod11 && 
				pmod13 && pmod17 && 
				pmod19 && pmod23 && pmod29 &&
				pmod31 && pmod37 && pmod41 && 
				pmod43 && pmod47 && pmod53 && pmod59 &&
				pmod61 && pmod67 && pmod71 && pmod73 && 
				pmod79 && pmod83 &&
				pmod89 && pmod97)
			{
				b_m = b_p = false;

				if( (tmod3 == 2) &&
					((pmod11 + 1 - 2*tmod11)%11) &&
					((pmod13 + 1 - 2*tmod13)%13) &&
					((pmod17 + 1 - 2*tmod17)%17) &&
					((pmod19 + 1 - 2*tmod19)%19) &&
					((pmod23 + 1 - 2*tmod23)%23) &&
					((pmod29 + 1 - 2*tmod29)%29) &&
					((pmod31 + 1 - 2*tmod31)%31) &&
					((pmod37 + 1 - 2*tmod37)%37) &&
					((pmod41 + 1 - 2*tmod41)%41) &&
					((pmod43 + 1 - 2*tmod43)%43) &&
					((pmod47 + 1 - 2*tmod47)%47) &&
					((pmod53 + 1 - 2*tmod53)%53) &&
					((pmod59 + 1 - 2*tmod59)%59) &&
					((pmod61 + 1 - 2*tmod61)%61) &&
					((pmod67 + 1 - 2*tmod67)%67) &&
					((pmod71 + 1 - 2*tmod71)%71) &&
					((pmod73 + 1 - 2*tmod73)%73) &&
					((pmod79 + 1 - 2*tmod79)%79) &&
					((pmod83 + 1 - 2*tmod83)%83) &&
					((pmod89 + 1 - 2*tmod89)%89) &&
					((pmod97 + 1 - 2*tmod97)%97) )
				{
					b_m = true;

					// We only compute p and t if p is not 
					// divisible by small primes.
					// We have 
					// t_1_2(c_r) = t_1_2(0) + 2*c_r,
					// p(c_r) = p(0) + 4*c_r*(t(0) + c_r).
					
					add( tmp_i, t_1_2, bigint( c_r ) );
					multiply( tmp_i, tmp_i, bigint( 4 * c_r ) );
					add( q, q, tmp_i ); // adapt q
					add( t_1_2, t_1_2, bigint( 2 * c_r ) ); // adapt t_1_2
					
					c_r = 0;

					if( q > end_of_search )
						in_interval = 0;

					if( in_interval && is_prime( q ) )
					{
						b_p = true;

						// compute order = q + 1 - 2 * t_1_2
						shift_left( order, t_1_2, 1 ); // order = 2*t_1_2
						subtract( order, q, order ); // order = q - 2*t_1_2
						inc( order );

						p.assign( q );
						
						if( is_cryptographically_strong( order ) )
						{
							return( true );
						}
					}
				}

				if( (tmod3 == 1) &&
					((pmod11 + 1 + 2*tmod11)%11) &&
					((pmod13 + 1 + 2*tmod13)%13) &&
					((pmod17 + 1 + 2*tmod17)%17) &&
					((pmod19 + 1 + 2*tmod19)%19) &&
					((pmod23 + 1 + 2*tmod23)%23) &&
					((pmod29 + 1 + 2*tmod29)%29) &&
					((pmod31 + 1 + 2*tmod31)%31) &&
					((pmod37 + 1 + 2*tmod37)%37) &&
					((pmod41 + 1 + 2*tmod41)%41) &&
					((pmod43 + 1 + 2*tmod43)%43) &&
					((pmod47 + 1 + 2*tmod47)%47) &&
					((pmod53 + 1 + 2*tmod53)%53) &&
					((pmod59 + 1 + 2*tmod59)%59) &&
					((pmod61 + 1 + 2*tmod61)%61) &&
					((pmod67 + 1 + 2*tmod67)%67) &&
					((pmod71 + 1 + 2*tmod71)%71) &&
					((pmod73 + 1 + 2*tmod73)%73) &&
					((pmod79 + 1 + 2*tmod79)%79) &&
					((pmod83 + 1 + 2*tmod83)%83) &&
					((pmod89 + 1 + 2*tmod89)%89) &&
					((pmod97 + 1 + 2*tmod97)%97) )
				{
					if( b_m == false )
					{
						add( tmp_i, t_1_2, bigint( c_r ) );
						multiply( tmp_i, tmp_i, bigint( 4 * c_r ) );
						add( q, q, tmp_i ); // adapt q
						add( t_1_2, t_1_2, bigint( 2 * c_r ) ); // adapt t_1_2
						c_r = 0;
					}

					if( q > end_of_search )
						in_interval = 0;

					if( in_interval && ( b_p == true ||
										 ( b_m == false && is_prime( q ) ) ) )
					{
						// compute order = q + 1 + 2 * t_1_2
						shift_left( order, t_1_2, 1 ); // order = 2*t_1_2
						add( order, q, order ); // order = q + 2*t_1_2
						inc( order );
						
						p.assign( q );

						if( is_cryptographically_strong( order ) )
						{
							return( true );
						}
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

			if( tmod3 == 1 && tmod5 == 3 && tmod7 == 4 )
				
				c = 17;

			else if( tmod5 == 3 && 
					 ( ( tmod3 == 2 && tmod7 == 6 ) ||
					   ( tmod3 == 1 && tmod7 == 5 ) ) )
				
				c = 12;

			else if( tmod3 == 1 && 
					 ( ( tmod5 == 1 && ( tmod7 == 1 || tmod7 == 2 ) ) ||
					   ( tmod5 == 2 && ( tmod7 == 4 || tmod7 == 5 ) ) ) )

				c = 11;
				
			else if( tmod3 == 2 && 
					 ( ( tmod5 == 2 && ( tmod7 == 2 || tmod7 == 3 ) ) ||
					   tmod5 == 3 || 
					   ( tmod5 == 4 && ( tmod7 == 4 || tmod7 == 5 ) ) ) )

				c = 7;

			else if( ( tmod3 == 2 && tmod5 == 2 ) ||
					 ( tmod3 == 1 && tmod5 == 1 ) )

				c = 6;

			else if( tmod3 == 1 && ( tmod5 == 2 || tmod5 == 3 ) )

				c = 5;

			else
			
				c = 1;
			

			addend_2 = 2 * c;
			addend_1 = addend_2 * addend_2;
			factor = 4 * c;
			c_r += c;

			pmod3  = (pmod3 + factor * tmod3 + addend_1) % 3;
			pmod5  = (pmod5 + factor * tmod5 + addend_1)%5;
			pmod7  = (pmod7 + factor * tmod7 + addend_1)%7;
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
