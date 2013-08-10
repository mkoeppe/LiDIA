// -*- C++ -*-
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
//	Author	: Thomas Pfahler (TPf), Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_UDIGIT_H_GUARD_
#define LIDIA_UDIGIT_H_GUARD_



//
// Exactly the following functions has to be
// defined and implemented in
// udigit_def.h and udigit_def.c:
//

// udigit max_udigit()
        //RV = maximal udigit
        //RADIX = max_udigit()+1

// unsigned int bits_per_udigit()
         // RV = number of bits a udigit consists of

// udigit max_udigit_modulus()
        //RV = maximal udigit which can be used in
        //     the _mod - functions below

// unsigned int bits_per_udigit_modulus()
         // RV = number of bits a udigit modulus consists of

// udigit add(udigit & sum, udigit a, udigit b, udigit carry = 0)
        //RV * RADIX + sum = a + b + carry

// udigit subtract(udigit & sum, udigit a, udigit b, udigit carry = 0)
        //-RV * RADIX + diff = a - b - carry

// udigit multiply(udigit & prod, udigit a, udigit b)
        //RV * RADIX + prod = a * b

// udigit divide(udigit & quot, udigit a1, udigit a0, udigit b)
        //a1 * RADIX + a0 = quot * b + RV
	//a1 must be less than b


//
// The input to the modular functions must
// be reduced modulo p.
//

// udigit add_mod(udigit a, udigit b, udigit p)
	//RV = (a + b) mod p
	//assumes a, b < p, p > 1

// udigit subtract_mod(udigit a, udigit b, udigit p)
	//RV = (a - b) mod p
	//assumes a, b < p, p > 1

// udigit multiply_mod(udigit a, udigit b, udigit p)
        // RV = (a * b) mod p
        // assumes a, b < p, p > 1


#ifndef LIDIA_LIDIA_H_GUARD_
# include	"LiDIA/LiDIA.h"
#endif



#include	"LiDIA/kernel/udigit_def.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



#ifndef HEADBANGER



udigit	max_udigit ();
udigit	max_udigit_modulus ();
unsigned int	bits_per_udigit ();
unsigned int	bits_per_udigit_modulus ();
udigit	add (udigit & sum, udigit a, udigit b, udigit carry = 0);
udigit	subtract (udigit & diff, udigit a, udigit b, udigit carry = 0);
udigit	multiply (udigit & prod, udigit a, udigit b);
udigit	divide (udigit & quot, udigit a1, udigit a0, udigit b);
udigit	add_mod (udigit a, udigit b, udigit p);
udigit	subtract_mod (udigit a, udigit b, udigit p);
udigit	multiply_mod (udigit a, udigit b, udigit p);



#if INLINE_INTERFACE
# include	"LiDIA/kernel/udigit_interface.h"
#endif



//
// The following functions are implemented in udigit.c
// and need not be changed for use with a different kernel:
//

udigit	power (udigit a, unsigned int e);
// compute a^e, without taking care of carrys

udigit	sqrt_mod (udigit a, udigit p);
// p is assumed to be prime

udigit	power_mod (udigit a, udigit e, udigit p);
// compute (a^e) mod p

bool	is_prime (udigit n, unsigned int trials = 9);
// returns true if n passes the Miller-Rabin test,
// otherwise it returns false

int	jacobi (udigit a, udigit b);
// returns 1 if a is a square mod b, -1 if not
// and zero if b divides a.

udigit	next_prime (udigit a);
// returns next prime which is larger than a.
// UEBERLAUF ?????

udigit	previous_prime (udigit a);
// returns the next prime smaller than a.

udigit	gcd (udigit a, udigit b);
// binary gcd algorithm

inline
udigit	lcm (udigit a, udigit b)
{
	udigit c;
	divide(c, 0, a, gcd(a, b));
	multiply(c, c, b);
	return c;
}

udigit	xgcd (udigit &u, udigit &v, udigit a, udigit b);
// extended binary gcd algorithm, we take care of overflows
// RV = u * a - v * b, |RV| = gcd(a, b)

udigit	negate_mod (udigit a, udigit p);
// RV = -a mod p
// assumes a < p, p > 1

udigit	invert_mod (udigit a, udigit p);
// invert a mod p, binary gcd algorithm,
// only coefficient for a is computed

inline
udigit	divide_mod (udigit a, udigit b, udigit p)
{
	udigit binv = invert_mod(b, p);
	return multiply_mod(a, binv, p);
}



#endif	// HEADBANGER



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_UDIGIT_H_GUARD_
