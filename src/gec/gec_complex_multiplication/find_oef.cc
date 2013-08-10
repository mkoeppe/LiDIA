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
//  File    : find_oef.cc
//  Author  : Harald Baier (HB)
//	Changes	: See CVS log
//
//==============================================================================================

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include "LiDIA/gec_complex_multiplication.h"
#include "LiDIA/osstream.h"
#include "LiDIA/error.h"

#include <string>
#include <map>
#include <fstream>
#ifdef VERBOSE
# include        <cassert>
#endif


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif


//
// this function implements algorithm findOEFField (11.3) from
// PhD
//
    bool
    gec_complex_multiplication::
    find_oef()
    {
	std::cout << upper_bound_k << std::endl;
	std::cout << lower_bound_bitlength_r << std::endl;
	std::cout << degree << std::endl;
	std::cout << word_length << std::endl;
	std::cout << lower_bound_class_number << std::endl;
	std::cout << upper_bound_class_number << std::endl;

	if( VERBOSE )
	    std::cout << "Searching the OEF....." << std::endl;

	long decreaser = 2 * degree;
   
	//
	// Some constants for reading discriminants
	//
	long remain;

	// ensure that degree | lower_bound_class_number
	remain = lower_bound_class_number % degree;
	if( remain != 0 )
	    lower_bound_class_number += ( degree - remain );
	
	// ensure that degree | upper_bound_class_number   
	remain = upper_bound_class_number % degree;
	upper_bound_class_number -=  remain;

	const int no_of_rounds = 6; // how often read files
	const bigint delta_interval( -1000000 );

	bigint rem;

	//
	// Some other variables
	//
	bigint delta_lower_bound, delta_upper_bound;
	long class_number;
	char class_number_as_string[ 10 ];
   
	bool h_is_even;

	//
	// small primes for sieving
	const lidia_size_t no_of_small_primes = 7;
	int small_primes[] = { 3, 5, 7, 11, 13, 17, 19 };
	long remainder_mod_small_primes[ no_of_small_primes ];
	lidia_size_t kk;
	bool bad_prime, bad_mod;
   
	//
	// define the interval of p
	// 
	bigint lower_bound_p, upper_bound_p;
	upper_bound_p.assign_one();
	shift_left( upper_bound_p, upper_bound_p, word_length );

	sqrt( lower_bound_p, upper_bound_p );
	subtract( lower_bound_p, upper_bound_p, lower_bound_p );

	// prepare upper_bound_p for efficient sieving
	remainder( rem, upper_bound_p, degree );
	subtract( upper_bound_p, upper_bound_p, rem - 1 );

	if( upper_bound_p.bit( 0 ) == 0 )
	    subtract( upper_bound_p, upper_bound_p, degree );

	// set the bounds appropriate
	add( upper_bound_p, upper_bound_p, decreaser );
	add( lower_bound_p, lower_bound_p, decreaser );

	//
	// trace t and y for representation 4q=t^2-delta*y^2
	//
	bigint t;
	bigint y;
	bigint order; // order of curve: order = q + 1 - t

	//
	// file handling
	//
	std::string path = path_to_discriminant_files();
	std::ifstream in;
   
	// to store positions in files
	std::map<long, std::streampos> position;
   
	for( int count_inits = 0; count_inits < no_of_rounds; count_inits++ )
	{
	    delta_upper_bound = count_inits * delta_interval;
	    delta_lower_bound = ( count_inits + 1 ) * delta_interval;

	    class_number = lower_bound_class_number;

	    while( class_number < upper_bound_class_number )
	    {
		h_is_even = (class_number % 2 == 0);

		osstream oss;
		oss << "h_" << class_number;
		std::string in_filename = path + extractString(oss);
		
		in.open( in_filename.c_str() );
		if(!in.good()) {
		    std::string msg("Can't open file ");
		    msg += in_filename;
		    lidia_error_handler("gec_complex_multiplication::"
					"find_oef()",
					msg.c_str());
		}
		
		// set the read position if count_inits > 0
		if( count_inits > 0 ) {
		    in.seekg( position[ class_number ] );
		    if(in.fail()) {
			std::string msg("Can't seek pos within file ");
			msg += in_filename;
			lidia_error_handler("gec_complex_multiplication::"
					    "find_oef()",
					    msg.c_str());
		    }
		}
		do
		{
		    bool bad_delta = true;

		    do
		    {
			in >> delta;
			if(in.fail()) {
			    std::string msg("Error while reading file ");
			    msg += in_filename;
			    lidia_error_handler("gec_complex_multiplication::"
						"find_oef()",
						msg.c_str());
			}

			remainder( rem, delta, 8 );

			if( rem == -3 )
			{
			    remainder( rem, delta, 3 );
					
			    std::cout << rem;

			    if( ! rem.is_zero() )
			    {
				bad_delta = false;
			    }
			}
					
			if( h_is_even )
			    in.ignore( INT_MAX,'\n' );

		    } while( bad_delta );

		    //
		    // assign p and initialize remainders
		    p.assign( upper_bound_p );

		    for( kk = 0; kk < no_of_small_primes; kk++ )
			remainder_mod_small_primes[ kk ]
			    = remainder( p, small_primes[ kk ] );

		    while( p > lower_bound_p + degree )
		    {
			bad_prime = true;
					
			while( bad_prime )
			{
			    bad_mod = true;

			    while( bad_mod )
			    {
				bad_mod = false;
				p = p - decreaser;

				for( kk = 0; kk < no_of_small_primes; kk++ )
				{
				    remainder_mod_small_primes[ kk ] = 
					(remainder_mod_small_primes[kk] 
					 - decreaser ) % small_primes[ kk ];
						
				    if( remainder_mod_small_primes[kk] == 0 )
					bad_mod = true;
				}
			    }
						
			    if( is_prime( p ) )
				bad_prime = false;
			}					

			//
			// check that p does not split to two principal ideals
			// HB 030515: first condition introduced
			// as proposed by Per Christoffersen
			if( 4 * p > -delta )
			    if( cornacchia( t, y, delta, p ) )
				continue;
					
			//
			// check that p^degree splits to two principal ideals
			//
			if( ! cornacchia_prime_power( t, y, delta, p, degree ) )
			    continue;

			power( q, p, degree );
	   
			if( VERBOSE )
			{
			    std::cout << "t = " << t << std::endl;
			    std::cout << "y = " << y << std::endl;
			    std::cout << "delta = " << delta << std::endl;
			    std::cout << "degree = " << degree << std::endl;
			}

			order = q + 1 - t;
			if( is_cryptographically_strong( order ) )  {
			    return( true );
			}
					
			order = q + 1 + t;
			if( is_cryptographically_strong( order ) ) {
			    return( true );
			}
		    }
				
		    if( h_is_even )
			in.ignore( INT_MAX,'\n' );

		} while( delta > delta_lower_bound );

		position[ class_number ] = in.tellg();
		in.clear();
		in.close();
		if(in.fail()) {
		    std::string msg("Error while closing file ");
		    msg += in_filename;
		    lidia_error_handler("gec_complex_multiplication::"
					"find_oef()",
					msg.c_str());
		}
			
		// ensure : degree | class_number 
		class_number += degree;
	    }
	}
	
	std::cout << "Returning false." << std::endl;

	return( false );
    }


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
