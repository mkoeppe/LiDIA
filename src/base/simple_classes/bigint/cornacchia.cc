//==============================================================================================
//
//	This file is part of LiDIA --- a library for computational number theory
//
//	Copyright (c) 1994--2001 the LiDIA Group.  All rights reserved.
//
//	See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//----------------------------------------------------------------------------------------------
//
//	$Id$
//
//	Author	: Volker Mueller (VM), Harald Baier (HB)
//	Changes	: See CVS log
//
//==============================================================================================

// functions: * cornacchia
//            * cornacchia_prime_power
//



// Description: determines (if possible) solutions for
//              x^2 + |D| * y^2 = 4 * p, where D = 0,1 mod 4, D < 0,
//              |D| < 4*p, p odd.
// Output: sets  x,y und returns true, if solution found
//	   otherwise return false


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigint.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



bool cornacchia (bigint & x, bigint & y, const bigint & DD, const bigint & p)
{
	bigint x0, a, b, l, r;
	bigint D, tmp2, D_abs, p_four;
	bigint rr;
	bool r_is_sqr;
	long s;

	if (!DD.is_negative())
		lidia_error_handler("cornacchia", "D not negative");

	s = (4 - (DD.least_significant_digit() & 3)) & 3; // DD mod 4

	shift_left(p_four, p, 2);
	D.assign(DD);

	if (!is_prime(p) || !p.is_positive() || p == 2)
		lidia_error_handler("cornacchia", "no odd prime number");

	if ((-D) >= p_four)
		lidia_error_handler("cornacchia", "|D| >= 4*p");

	if (s != 1 && s != 0)
		lidia_error_handler("cornacchia", "D != 0 or 1 mod 4");

	if (jacobi(D, p) == -1) {
		return false; // (D/p) = -1 -->no solution
	}
	else {
		ressol(x0, p+D, p);
		if (x0.is_even() != D.is_even())
			subtract(x0, p, x0);

		shift_left(a, p, 1);
		b.assign(x0);

		shift_left(l, p, 2);
		sqrt(l, l); // l = floor(2*sqrt(p)) = floor(sqrt(4p))

		while (b > l) {
			remainder(r, a, b);
			a.assign(b);
			b.assign(r);
		}

		square(l, b);
		subtract(a, p_four, l);
		// a = 4p - b^2 now

		D.negate();
		remainder(r, a, D);

		if (!r.is_zero()) {
			return false;
		}
		else {
			divide(r, a, D);
			r_is_sqr = is_square (rr, r);
			if (!r_is_sqr) {
				return false;
			}
			else {
				x.assign(b);
				y.assign(rr);
				return true;
			}
		}
	}
}


// Function    : cornacchia_prime_power 
// Author      : Harald Baier
// Remark      : extends cornacchia to prime powers with odd exponent
// Algorithm   : bases on Book of Cohen, Algorithm 1.5.3 
// Description: determines (if possible) solutions for 
//              x^2 + |D| * y^2 = 4 * p^exp, where D = 0,1 mod 4, D < 0,
//              |D| < 4*p^exp, p and exp odd 
// Output: sets  x,y und returns true, if solution found
//	   otherwise return false
//

bool cornacchia_prime_power( bigint & x, bigint & y, 
							 const bigint & DD, 
							 const bigint & p, const int exp )
{
  bigint x0, a, b, l, r ;
  bigint D, tmp2, D_abs, modulus_four;
  bigint rr;
  bigint modulus;
  bool r_is_sqr;
  long s;

  if( ! DD.is_negative() )
     lidia_error_handler("cornacchia_prime_power","D not negative");

  s = (4 - (DD.least_significant_digit() & 3)) & 3;    // DD mod 4

//  if( ! is_prime(p)  || ! p.is_positive() || p == 2)
//    lidia_error_handler("cornacchia_prime_power","p no odd prime number");
  
  if( exp <= 0 || ( exp & 1 == 0 ) )
	  lidia_error_handler("cornacchia_prime_power","exp no odd integer");

  power( modulus, p, exp );
  shift_left(modulus_four, modulus, 2);
  D.assign(DD);

  if ((-D) >= modulus_four)
	  lidia_error_handler("cornacchia","|D| >= 4*modulus");

  if (s != 1 && s != 0)
	  lidia_error_handler("cornacchia","D != 0 or 1 mod 4");

  // as exp is odd, we have in case of gcd(D,p) = 1: 
  // D is a square mod p^exp  <=>  D is a square mod p
  if( jacobi( D, p) != 1 )  
  {
      return false ;   // (D/p) = -1: no solution; (D/p) = 0: (D,p) = p
  }
  else 
  {             
      ressol_prime_power( x0, modulus + D, p, exp );   
      if( x0.is_even() != D.is_even() )
         subtract(x0, modulus, x0);

     shift_left(a, modulus, 1);
     b.assign( x0 );
      
     shift_left( l, modulus, 2 );
     sqrt( l, l );        // l = floor(2*sqrt(p^exp)) = floor(sqrt(4p^exp))
 
     while( b > l )
	 {
		 remainder(r, a, b);
		 a.assign(b);
		 b.assign(r);
	 }

     square( l, b );
     subtract( a, modulus_four, l );// a = 4*p^exp - b^2
    
     D.negate();
     remainder( r, a, D );

     if( !r.is_zero() )
	 {
		 return false;
	 }
     else
	 {
		 divide(r, a, D);
		 r_is_sqr = is_square( rr, r ); 
		 if ( ! r_is_sqr )  
		 {
			 return false;
		 }
		 else 
		 {
			 x.assign(b);
			 y.assign(rr);
			 return true;
		 }
	 }
  } 
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
