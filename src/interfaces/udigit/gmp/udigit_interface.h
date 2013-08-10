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
// Use of gmp kernel functions.
// RV = return value.
//



inline udigit
max_udigit ()
	// RV = maximal udigit
{
	return ~(static_cast<udigit>(0));
}



inline udigit
max_udigit_modulus ()
	// RV = maximal udigit that can be used as a modulus
{
	return ~(static_cast<udigit>(0));
}



inline unsigned int
bits_per_udigit ()
	// RV = the number of bits in a udigit
{
	return UDIGIT_NBITS;
}



inline unsigned int
bits_per_udigit_modulus ()
	// RV = the number of bits in a udigit modulus
{
	return UDIGIT_NBITS;
}



inline udigit
add (udigit & sum, udigit a, udigit b, udigit carry)
        //RV * 2^base + sum = a + b + carry
{
	udigit res[2] = { a, 0UL };

	mpn_add_1(res, res, 2, b);
	if (carry != 0) {
		mpn_add_1(res, res, 2, carry);
	}
	sum = res[0];
	return res[1];
}



inline udigit
subtract (udigit & diff, udigit a, udigit b, udigit carry)
        //-RV * 2^base + diff = a - b - carry
{
	udigit res[2] = { a, 0UL };

	mpn_sub_1(res, res, 2, b);
	if (carry != 0) {
		mpn_sub_1(res, res, 2, carry);
	}
	diff = res[0];
	return -res[1];
}



inline udigit
multiply (udigit & prod, udigit a, udigit b)
        //RV * 2^base + prod = a * b
{
	udigit res[2] = { a, 0UL };

	mpn_mul_1(res, res, 2, b);
	prod = res[0];
	return res[1];
}



inline udigit
divide (udigit & quot, udigit a1, udigit a0, udigit b)
	//a1 * 2^base + a0 = quot * b + RV
	//a1 must be less than b
{
	udigit a[2] = { a0, a1 };
        udigit tmp[2];
        udigit rv = mpn_divrem_1(tmp, 0, a, 2, b);
        quot = tmp[0];
        return rv;
}



//
// modular functions, we assume that input is reduced modulo modul
//

inline udigit
add_mod (udigit a, udigit b, udigit p)
	//RV = (a + b) % p
	//assumes a,b < p
{
	udigit res[2] = { a, 0UL };

	mpn_add_1(res, res, 2, b);
	if (res[1] != 0) {
		mpn_sub_1(res, res, 2, p);
	}
	else {
		if (res[0] >= p) {
			res[0] -= p;
		}
	}

	return res[0];
}



inline udigit
subtract_mod (udigit a, udigit b, udigit p)
	//RV = (a - b) % p
	//assumes a,b < p
{
	udigit res;

	if (a >= b) {
	  	res = a - b;
	}
	else {
	  	res = (p - b) + a;
	}
	return res;
}



inline udigit
multiply_mod (udigit a, udigit b, udigit p)
	// a,b must be less than p
{
	udigit res[2] = { a, 0UL };

	mpn_mul_1(res, res, 2, b);
	return mpn_mod_1(res, 2, p);
}
