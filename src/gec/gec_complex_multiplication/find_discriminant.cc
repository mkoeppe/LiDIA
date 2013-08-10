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
//  File    : find_discriminant.cc
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
#include <climits>
#include <cstdlib>
#ifdef VERBOSE
# include        <cassert>
#endif


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif

    std::string
    gec_complex_multiplication::path_to_discriminant_files() {
	std::string path;

	// did the user provide the path in the environment?
	if( getenv( "LIDIA_DISCRIMINANTS" ) != NULL ) {
	    path = getenv( "LIDIA_DISCRIMINANTS" );
	}

	// alternatively, take the value defined in path.h
	if(path.length() == 0) {
	    path = LIDIA_DISCRIMINANTS_DB; 
	    if(path.length() == 0) {
		lidia_error_handler("gec_complex_multiplication::"
				    "path_to_discriminant_files()",
				    "LIDIA_DISCRIMINANTS_DB (defined "
				    "in LiDIA/path.h) is empty");
	    }
	}

	// enforce a trailing '/'
	char const dirsep = '/';
	if(path[path.length()-1] != dirsep) {
	    path += dirsep;
	}

	return path;
    }
	

//
// This function implements algorithm findDiscriminant (algorithm 5.1) of PhD.
// For details, we refer to PhD.
//
    bool
    gec_complex_multiplication::find_discriminant() {
	//
	// ensure that q is prime
	//
	if( degree == 1 )
	    p.assign( q );
	else
	    lidia_error_handler( "gec_complex_multiplication::find_discriminant",
				 "degree has to equal 1");

	//
	// Some constants which have to be adapted if source changes
	//
	const int no_of_small_primes = 25; // #{l prime: 3 <= l <= 101}

	// Discriminants are searched in the interval specified by
	// the following bounds; if no bound was set before,
	// the default values are used.
	if( lower_bound_class_number == 0 )
	    lower_bound_class_number = default_lower_bound_class_number;

	if( upper_bound_class_number == 0 )
	    upper_bound_class_number = default_upper_bound_class_number;

	// to decide whether curve of prime order is searched
	bool curve_of_prime_order = false;
	if( upper_bound_k == 1 )
	    curve_of_prime_order = true;

	// remainder mod 8 (of p or delta)
        long rem_mod_8;

	//
	// other variables
	//
	bigint t, y, order;
	bigint B_delta;
	int delta_no_of_prime_factors;
	bigint tmp_i, local_delta;
	bigint delta_lower_bound, delta_upper_bound;

	long class_number = lower_bound_class_number;
	lidia_size_t count;
	char class_number_as_string[ 10 ];

	bool good_delta, h_is_even;

	// 
	// compute B_p
	//
	bigint small_prime, B_p( 0 ), addend( 1 );

	int small_primes[] = { 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 
			       41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 
			       83, 89, 97, 101 };

	for( int ii = 0; ii < no_of_small_primes; ii++ ) 
	{	
	    small_prime.assign( bigint( small_primes[ ii ] ) ); 
	    if( jacobi( p, small_prime ) == -1 )
		add( B_p, B_p, addend );

	    shift_left( addend, addend, 1 );
	}

	// to consider congruence condition due to residue of p mod 8:
	remainder( rem_mod_8, p, 8 );

	addend.assign_one();
	if( rem_mod_8 == 3 )
	{
	    shift_left( addend, addend, no_of_small_primes );
	    add( B_p, B_p, addend );
	}
	else if( rem_mod_8 == 5 )
	{
	    shift_left( addend, addend, no_of_small_primes + 1 );
	    add( B_p, B_p, addend );
	}
	else if( rem_mod_8 == 7 )
	{
	    shift_left( addend, addend, no_of_small_primes + 2 );
	    add( B_p, B_p, addend );
	}

	//
	// file handling
	//

	std::string path = path_to_discriminant_files();
	std::ifstream in;
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
				    "find_discriminant()",
				    msg.c_str());
	    }
		
	    while( true )
	    {
		in >> delta;
		if(in.eof() && in.fail() ) {
		    break;
		}
		if(in.fail()) {
		    std::string msg("Error while reading file ");
		    msg += in_filename;
		    lidia_error_handler("gec_complex_multiplication::"
					"find_discriminant()",
					msg.c_str());
		}

		// if upper_bound_k = 1 ensure delta = 5 mod 8
		if( curve_of_prime_order == true)
		{
		    remainder( rem_mod_8, delta, 8 );
		    if( rem_mod_8 != -3 )
		    {
			// discard current line
			in.ignore( INT_MAX,'\n' );
			continue;
		    }
		}

		abs_delta.assign( delta );
		abs_delta.negate();

		good_delta = true;

		if( h_is_even )
		{
		    in >> B_delta;
		    if(in.fail()) {
			std::string msg("Error while reading file ");
			msg += in_filename;
			if(in.eof()) {
			    msg += " (unexpected EOF)";
			}
			lidia_error_handler("gec_complex_multiplication::"
					    "find_discriminant()",
					    msg.c_str());
		    }

		    bitwise_and( tmp_i, B_p, B_delta );
		    if( ! tmp_i.is_zero() )
		    {
			// discard line
			in.ignore( INT_MAX,'\n' );
			continue;
		    }

		    in >> delta_no_of_prime_factors;
		    if(in.fail()) {
			std::string msg("Error while reading file ");
			msg += in_filename;
			if(in.eof()) {
			    msg += " (unexpected EOF)";
			}
			lidia_error_handler("gec_complex_multiplication::"
					    "find_discriminant()",
					    msg.c_str());
		    }
		    for( count = 0; 
			 count < delta_no_of_prime_factors; count++ ) {
			in >> tmp_i;
			if(in.fail()) {
			    std::string msg("Error while reading file ");
			    msg += in_filename;
			    if(in.eof()) {
				msg += " (unexpected EOF)";
			    }
			    lidia_error_handler("gec_complex_multiplication::"
						"find_discriminant()",
						msg.c_str());
			}
					
			if( jacobi( p, tmp_i ) != 1 ) {
			    good_delta = false;
			    break;
			}
		    }
		}
				
		if( good_delta == false || jacobi( delta, p) != 1 ) {
		    in.ignore( INT_MAX ,'\n' );
		    continue;
		}
			
		if( ! cornacchia( t, tmp_i, delta, p ) ) {
		    in.ignore( INT_MAX ,'\n' );
		    continue;
		}

		add( order, p, t );
		inc( order );
			
		if( is_cryptographically_strong( order ) ) {
		    h = class_number;

		    return( true );
		}
			
		shift_left( t, t, 1 );
		subtract( order, order, t );

		if( is_cryptographically_strong( order ) ) {
		    h = class_number;
					
		    return( true );
		}
	    }
	    in.clear();
	    in.close();
	    if(in.fail()) {
		std::string msg("Error while closing file ");
		msg += in_filename;
		lidia_error_handler("gec_complex_multiplication::"
				    "find_discriminant()",
				    msg.c_str());
	    }
		
	    class_number++;
	}

	return( false );
    }


    bool
    gec_complex_multiplication::find_discriminant_inefficient()
    {
	//
	// ensure that q is prime
	//
	if( degree == 1 )
	    p.assign( q );
	else
	    lidia_error_handler( "gec_complex_multiplication::find_discriminant",
				 "degree has to equal 1");

	// Discriminants are searched in the interval specified by
	// the following bounds; if no bound was set before,
	// the default values are used.
	if( lower_bound_class_number == 0 )
	    lower_bound_class_number = default_lower_bound_class_number;

	if( upper_bound_class_number == 0 )
	    upper_bound_class_number = default_upper_bound_class_number;

	// to read discriminant from file
	int discriminant;

	// to decide whether curve of prime order is searched
	bool curve_of_prime_order = false;
	if( upper_bound_k == 1 )
	    curve_of_prime_order = true;

	// remainder mod 8 (of p or delta)
	long rem_mod_8;

	//
	// other variables
	//
	bigint t, y, order;
	bigint B_delta;
	bigint tmp_i, local_delta;
	bigint delta_lower_bound, delta_upper_bound;

	long class_number = lower_bound_class_number;
	char class_number_as_string[ 10 ];
	char in_file[ 150 ], path[ 150 ];

	//
	// file handling
	//
	if( getenv( "LIDIA_DISCRIMINANTS" ) != NULL )
	    strcpy( path, getenv( "LIDIA_DISCRIMINANTS" ) );
	else
	    lidia_error_handler( "gec_complex_multiplication::find_discriminant",
				 "env LIDIA_DISCRIMINANTS not set.");

	strcat( path, "/h_" );
	std::ifstream in;

	while( class_number < upper_bound_class_number )
	{
	    sprintf( class_number_as_string, "%ld", class_number );
	    strcpy( in_file, path );
	    strcat( in_file, class_number_as_string );
			
	    in.open( in_file );

	    while( true )
	    {
		in >> discriminant;

		if( in.eof() )
		    break;

		delta.assign( bigint( discriminant ) );

		// if upper_bound_k = 1 ensure delta = 5 mod 8
		if( curve_of_prime_order == true)
		{
		    remainder( rem_mod_8, delta, 8 );
		    if( rem_mod_8 != -3 )
		    {
			in.ignore( INT_MAX,'\n' );
			continue;
		    }
		}

		abs_delta.assign( delta );
		abs_delta.negate();

		if( ! cornacchia( t, tmp_i, delta, p ) )
		{
		    in.ignore( INT_MAX ,'\n' );
		    continue;
		}

		add( order, p, t );
		inc( order );
			
		if( is_cryptographically_strong( order ) )
		{
		    h = class_number;

		    return( true );
		}
			
		shift_left( t, t, 1 );
		subtract( order, order, t );

		if( is_cryptographically_strong( order ) )
		{
		    h = class_number;
					
		    return( true );
		}
	    }

	    in.close();
		
	    class_number++;
	}

	return( false );
    }


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
