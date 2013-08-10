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
//  File    : prime_proof_ecpp_discriminant.cc
//  Author  : Harald Baier (HB),
//	      small changes by Jochen Hechler (JH)
//      Changes : See CVS log
//
//==============================================================================================



#include	<iostream>
#include	"LiDIA/LiDIA.h"
#include	"LiDIA/error.h"
#include        "LiDIA/gec_complex_multiplication.h"
#include        "LiDIA/prime_proof.h"
#include	"LiDIA/bigint.h"
#include        "LiDIA/osstream.h"

#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif

    bool prime_proof::find_next_discriminant(long start_h, long end_h)
    {

	const int no_of_small_primes = 25; // #{l prime: 3 <= l <= 101}

	// remainder mod 8 (of n or delta)
	long rem_mod_8;

	// other variables
	bigint abs_D;
	bigint t,y;
	bigint B_delta;
	int delta_no_of_prime_factors;
	bigint tmp_i, local_delta;
	class_number = cl_min;
	long upper_bound_class_number = cl_max;
	lidia_size_t count;
	char class_number_as_string[ 10 ];
	
	bool h_is_even;


	// compute B_p
	bigint small_prime, B_n( 0 ), addend( 1 );
	int small_primes[] = { 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37,
			       41, 43, 47, 53, 59, 61, 67, 71, 73, 79,
			       83, 89, 97, 101 };
	for( int ii = 0; ii < no_of_small_primes; ii++ )
	{
	    small_prime.assign( bigint( small_primes[ ii ] ) );
	    if( jacobi( n, small_prime ) == -1 )
		add( B_n, B_n, addend );
	    shift_left( addend, addend, 1 );
	}

	// to consider congruence condition due to residue of n mod 8:
	remainder( rem_mod_8, n, 8 );

	// file handling
	std::string path(gec_complex_multiplication::
			 path_to_discriminant_files());
	std::ifstream in;
	int discriminant;
	while( class_number < upper_bound_class_number )
	{
	    h_is_even = (class_number % 2 == 0);

	    osstream oss;
	    oss << "h_" << class_number;
	    std::string in_filename = path + extractString(oss);

	    in.open(in_filename.c_str());
	    if(!in.good()) {
		std::string msg("Can't open file ");
		msg += in_filename;
		lidia_error_handler("prime_proof::find_next_discriminant",
				    msg.c_str());
	    }
	    while( true )
	    {
		in >> discriminant;
		if(in.eof() && in.fail()) {
		    break;
		}
		if(in.fail()) {
		    std::string msg("Error while reading file ");
		    msg += in_filename;
		    lidia_error_handler("prime_proof::find_next_discriminant",
					msg.c_str());
		}
                 
		abs_D.assign( bigint( discriminant ) );
		remainder( rem_mod_8, abs_D, 8 );
		if(rem_mod_8 == 3)
		{
		    // skip this case, polynomials might be two high
		    in.ignore( INT_MAX,'\n' );
		    continue;
		}

		D.assign(abs_D);
		abs_D.negate();
		if( h_is_even )
		{
		    in >> B_delta; 
		    if(in.fail()) {
			std::string msg("Error while reading file ");
			msg += in_filename;
			if(in.eof()) {
			    msg += " (unexpected EOF)";
			}
			lidia_error_handler("prime_proof::"
					    "find_next_discriminant",
					    msg.c_str());
		    }
		    
		    bitwise_and( tmp_i, B_n, B_delta );
		    if( ! tmp_i.is_zero() )
		    {
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
			lidia_error_handler("prime_proof::"
					    "find_next_discriminant",
					    msg.c_str());
		    }
		    for( count = 0;
			 count < delta_no_of_prime_factors; count++ )
		    {
			in >> tmp_i;
			if(in.fail()) {
			    std::string msg("Error while reading file ");
			    msg += in_filename;
			    if(in.eof()) {
				msg += " (unexpected EOF)";
			    }
			    lidia_error_handler("prime_proof::"
						"find_next_discriminant",
						msg.c_str());
			}
			if( jacobi( n, tmp_i ) != 1 )
			{
			    break;
			}
		    }
		}
		if( abs_D < 5) {
		    continue;
		}
		if(  jacobi( D, n) != 1 )
		{
		    in.ignore( INT_MAX ,'\n' );
		    continue;
		}
		if( 4*n < abs_D) {
		    continue;
		}
		if( ! cornacchia( t, y, D, n ) )
		{
		    in.ignore( INT_MAX ,'\n' );
		    continue;
		}
		if( is_good_order(t,y))
		{
		    return( true );
		}	
		if(gcd>1) {
		    return( false );
		}
	    }

	    in.clear();
	    in.close();
	    if(in.fail()) {
		std::string msg("Error while closing file ");
		msg += in_filename;
		lidia_error_handler("prime_proof::"
				    "find_next_discriminant",
				    msg.c_str());
	    }

	    class_number++;
	    if(verbose) {
		std::cout<<"ecpp:next class: "<<class_number<<std::endl;
	    }
	    if(class_number==50 && ecpp_order_mode==1){
		class_number = 1;
		mode_one_tried = true;
	    }
	}
	
	return( false );
    }

#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
