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



//
// Arithmetical operations for type udigit.
// Use of libI kernel functions.
// RV = return value.
//



inline udigit
max_udigit ()
{
	return ~(static_cast<udigit>(0));
}



inline udigit
max_udigit_modulus ()
{
	return ~(static_cast<udigit>(0));
}



inline unsigned int
bits_per_udigit ()
{
	return UDIGIT_NBITS;
}



inline unsigned int
bits_per_udigit_modulus ()
{
	return UDIGIT_NBITS;
}



inline udigit
add (udigit & sum,
	    udigit a, udigit b, udigit carry)
        //RV * 2^base + sum = a + b + carry
{
	return DigitAdd(&sum, a, b, carry);
}



inline udigit
subtract (udigit & diff,
		 udigit a, udigit b,
		 udigit carry)
        //-RV * 2^base + diff = a - b - carry
{
	return DigitSub(&diff, a, b, carry);
}



inline udigit
multiply (udigit & prod,
		 udigit a, udigit b)
        //RV * 2^base + prod = a * b
{
	return DigitMult(& prod, a, b);
}



inline udigit
divide (udigit & quot,
	       udigit a1, udigit a0, udigit b)
	//a1 * 2^base + a0 = quot * b + RV
	//a1 must be less than b
{
	return DigitDiv(& quot, a1, a0, b);
}



//
// modular functions, we assume that input is reduced modulo modul
//

inline udigit
add_mod(udigit a, udigit b, udigit p)
	//RV = (a + b) % p
	//assumes a,b < p
{
	udigit carry, rem;

	carry = DigitAdd(&rem, a, b, 0);

	if (carry) {
		rem = (max_udigit() - p) + 1 + rem;
	}
	else {
		if (rem >= p) {
			rem = rem - p;
		}
	}

	return rem;
}



inline udigit
subtract_mod (udigit a, udigit b, udigit p)
	//RV = (a - b) % p
	//assumes a,b < p
{
	udigit rem;

	if (a >= b) {
		rem = a - b;
	}
	else {
		rem = (p - b) + a;
	}
	return rem;
}



inline udigit
multiply_mod (udigit a, udigit b, udigit p)
	// a,b must be less than p
{
	udigit high, low;

	high = DigitMult(&low, a, b);
	return DigitDiv(&low, high, low, p);
}
