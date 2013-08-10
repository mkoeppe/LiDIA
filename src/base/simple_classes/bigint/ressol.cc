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
//	Author	: Thomas Denny (TD), Andreas M"uller (AM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigint.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// Input:
//    prime number p
//    integer aa
//
// Output:
//    square root r of aa mod p
//    0 <= r < p
//
// Condition: p = 2^s * (2k+1) + 1
//	      where s fits in a long
//
// Algorithm: Shanks-Ressol
//


void ressol(bigint & x, const bigint & aa, const bigint & p)
{
	bigint v, one, a(aa);

	if (aa.is_negative()) {
		add(a, p, a % p);
		remainder(a, a, p);
	}

	if (a.is_zero()) {
		x.assign_zero();
		return;
	}

	if (p == bigint(2)) {
		x.assign(a);
		return;
	}

	//MM: used to avoid repeated casts bigint(1)
	one.assign_one();

	// (1) p = 3 mod 4
	if ((p.least_significant_digit() & 0x3) == 3) {
		if (jacobi(a, p) == 1) {
			add (v, p, one);
			shift_right(v, v, 2);
			power_mod(x, a, v, p);
			return;
		}
		else
			lidia_error_handler("ressol", "a is not a quadratic residue mod p");
	}

	//
	// MM: To be honest declare s and t as longs,
	// otherwise you have to use intify in the
	// "while (n > one) loop"
	//

	bigint r, n, c, z, k;
	long   s, t;

	// (2) initialisation:
	// compute k, s : p = 2^s (2k+1) + 1
	// (MM: shift_right as part of the loop)
	subtract(k, p, one);
	s = 0;
	while (k.is_even()) {
		s++;
		k.divide_by_2();
	}

	dec(k);
	k.divide_by_2();

	// (3) initial values
	power_mod(r, a, k, p);

	square(n, r);
	remainder(n, n, p); // n = (r * r) % p;
	multiply(n, n, a);
	remainder(n, n, p); // n = (n * a) % p;
	multiply(r, r, a);
	remainder(r, r, p); // r = (r * a) % p;

	if (n.is_one()) {
		x.assign(r);
		return;
	}

	// (4) non-quadratic residue
	z.assign(bigint(2));
	while (jacobi(z, p) == 1) inc(z);

	v.assign(k);
	v.multiply_by_2();
	inc(v); // v = (k << 1) + 1;
	power_mod(c, z, v, p);

	// (5) iteration
	while (n > one) {
		k.assign(n);
		t = s;
		s = 0;

		while (!k.is_one()) {
			square(k, k);
			remainder(k, k, p); // k = (k * k) % p;
			s++;
		}

		t = t - s;
		if (t == 0)
			lidia_error_handler("ressol", "a is not a quadratic residue mod p");

		shift_left(v, one, t-1); // v = 1 << (t-1);
		power_mod(c, c, v, p);
		multiply(r, r, c);
		remainder(r, r, p); // r = (r * c) % p;
		square(c, c);
		remainder(c, c, p); // c = (c * c) % p;
		multiply(n, n, c);
		remainder(n, n, p); // n = (n * c) % p;
	}

	x.assign(r);
	return;
}



void ressol_p(bigint & x, const bigint & a, const bigint & p)
{
	if (!p.is_prime())
		lidia_error_handler("ressol", "p is not prime!");
	else
		ressol(x, a, p);
}

// function: ressol_prime_power (HB, November 2000)
//
// Input:
//    integer aa
//    prime number p
//    exponent exp (it has to be an odd positive integer)
//
// Output:
//    square root r of aa mod p^exp
//    0 <= r < p^exp
// 
// Conditions: 
//   * p^exp = 2^s * (2k+1) + 1 , where s fits in a long
//	 * aa and p are coprime     
//
// Algorithm: Extension of Shanks-RESSOL to prime powers
//


void ressol_prime_power(bigint & x, const bigint & aa, 
						const bigint & p, const int exp)
{
  bool DEBUG = false;

  // check some conditions  
  if( exp <= 0 )
	  lidia_error_handler( "ressol_prime_power", "exponent is not a positive integer" );

  if( exp % 2 == 0 )
	  lidia_error_handler( "ressol_prime_power", "exponent is not odd" );

  bigint v, one(1), a(aa), modulus, euler_phi, tmp;

  // initialize:
  //   * the modulus as p^exp
  //   * the euler_phi value: euler_phi = p^(exp - 1) * ( p - 1 )
  power( tmp, p, exp - 1 );
  multiply( modulus, tmp, p );
  subtract( euler_phi, modulus, tmp );

/*  if( DEBUG )
/  {
	  cout << "p = " << p << endl;
	  cout << "exp = " << exp << endl;
	  cout << "p^exp = " << modulus << endl;
	  cout << "euler_phi = " << euler_phi << endl;
  }
*/

  if (aa.is_negative())
  {
	  add( a, p, a % modulus );
	  remainder( a, a, modulus );
  }


  /*
  MM: To be honest declare s and t as longs,
      otherwise you have to use intify in the
      "while (n > one) loop"
  */

  bigint r, n, c, z, k;
  long   s, t;

  // (2) initialisation:
  // compute k, s : euler_phi = 2^s (2k+1) 
  // (MM: shift_right as part of the loop)
  k.assign( euler_phi );
  s = 0;

  while( k.is_even() )
  {
	  s++;
	  k.divide_by_2(); 
  }
  
  dec( k );
  k.divide_by_2();
  
  // (3) initial values: 
  // n = a^(2k+1) mod p^exp,
  // r = a^(k+1) mod p^exp   (initial value of root)
  power_mod( r, a, k, modulus );
  
  square(n, r);
  remainder(n, n, modulus); // n = (r * r) % modulus;
  multiply(n, n, a);
  remainder(n, n, modulus); // n = (n * a) % modulus;
  multiply(r, r, a);
  remainder(r, r, modulus); // r = (r * a) % modulus;
  
  if( n.is_one() )
  {
	  x.assign( r );
	  return;
  }
  
  // (4) find a quadratic non-residue mod p^exp
  // as exp is odd it suffices to compute the legendre symbol modulo p
  z.assign( ( bigint ) 2 );
  while ( jacobi( z, p ) == 1 ) 
	  inc(z);
  
  v.assign(k);
  v.multiply_by_2();
  inc(v); // v = (k << 1) + 1;
  power_mod(c, z, v, modulus); // c generates the 2-Sylow subgroup
  
  // (5) iteration
  while( n > one )
  {
      k.assign(n);
      t  = s;
      s  = 0;
      
      while ( ! k.is_one() )
	  {
		  square(k, k);
		  remainder(k, k, modulus); // k = (k * k) % modulus; 
		  s++; 
	  }
      
      t = t - s; 
      if (t == 0)
		  lidia_error_handler("ressol_prime_power", "a is not a quadratic residue mod p^k");
       
      shift_left(v, one, t-1); // v = 1 << (t-1);
      power_mod(c, c, v, modulus);
      multiply(r, r, c);
      remainder(r, r, modulus); // r = (r * c) % modulus;
      square(c, c);
      remainder(c, c, modulus); // c = (c * c) % modulus;
      multiply(n, n, c);
      remainder(n, n, modulus); // n = (n * c) % modulus;
    } 
  
  x.assign(r);
  return;
} 


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
