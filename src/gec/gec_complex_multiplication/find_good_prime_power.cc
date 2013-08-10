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
//  File    : find_good_prime_power.cc
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
// the algorithm tries to find a prime power suitable for
// the CM-method of discriminant Delta. However, the implementation
// at hand is naive and not efficient.
//
bool
gec_complex_multiplication::
find_good_prime_power()
{
   if( VERBOSE )
   {
       std::cout << std::endl << "************************" << std::endl;
       std::cout << "Generating the prime power....." << std::endl;
   }

   //
   // find an appropriate begin for prime search
   //
   lidia_size_t bitlength_p = lower_bound_bitlength_r + 
	   upper_bound_k.bit_length();
   
   bitlength_p /= degree;
   bitlength_p++;
   if( bitlength_p % 2 )
	   bitlength_p++;

   p.assign_one();
   shift_left( p, p, bitlength_p );

   //
   // trace t and y for representation 4q=t^2-delta*y^2
   //
   bigint t;
   bigint y;
   bigint order; // order of curve: order = q + 1 - t

   bool good_prime; // to test the order of the prime ideal in quadratic order
   int upper_bound_of_d = degree/2;
   int d;

   while( true )
   {
	   good_prime = false;

	   p = next_prime( p );

	   //
	   // Check that the order of the prime ideal in O_D containing p 
	   // has order degree.
	   // We do this by checking that for all divisors d of degree 
	   // p^d does not factor as principal ideals in O_D
	   //
	   for( d = 1; d <= upper_bound_of_d; d++ )
	   {
		   if( degree % d == 0 && cornacchia_prime_power( t, y, delta, p, d ) )
			   break;
	   }

	   if( d == upper_bound_of_d + 1 )
		   good_prime = true;

	   if( VERBOSE )
		   std::cout << "p = " << p << std::endl;

	   if( good_prime && cornacchia_prime_power( t, y, delta, p, degree ) )
	   {
		   power( q, p, degree );

		   if( VERBOSE )
		   {
			   std::cout << "t = " << t << std::endl;
			   std::cout << "y = " << y << std::endl;
			   std::cout << "delta = " << delta << std::endl;
			   std::cout << "degree = " << degree << std::endl;
		   }

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
