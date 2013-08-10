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
//	$Id: udigit_interface.h,v 2.2 2001/05/29 09:07:09 hamdy Exp $
//
//	Author	:
//	Changes	: See CVS log
//
//==============================================================================================


//
// Arithmetical operations for type udigit.
// Use of piologie kernel functions.
// RV = return value.
//


class uNumberBase : public NumberBase
{
};

static uNumberBase udigit_NumberBase;



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
	    udigit a, udigit b, udigit carry = 0)
        //RV * 2^base + sum = a + b + carry
{
	a += carry;
	carry = (a < carry);
	a += b;
	carry += (a < b);
	sum = a;
	return carry;
}



inline udigit
subtract (udigit & diff,
		 udigit a, udigit b,
		 udigit carry = 0)
        //-RV * 2^base + diff = a - b - carry
{
	udigit tmp = a - carry;


	carry = (tmp > a);
	a = tmp - b;
	carry += (a > tmp);
	diff = a;
	return carry;
}



inline udigit
multiply (udigit & prod,
		 udigit a, udigit b)
        //RV * 2^base + prod = a * b
{
	udigit c;


	udigit_NumberBase.digitmul(a, b, c, prod);
	return c;
}



inline udigit
divide (udigit & quot,
	       udigit a1, udigit a0, udigit b)
	//a1 * 2^base + a0 = quot * b + RV
	//a1 must be less than b
{
	udigit r;


	udigit_NumberBase.digitdiv(a1, a0, b, quot, r);
	return r;
}



//
// modular functions, we assume that input is reduced modulo modul
//

inline udigit
add_mod (udigit a, udigit b, udigit p)
	//RV = (a + b) % p
	//assumes a,b < p
{
	udigit carry, rem;


	carry = add(rem, a, b);

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
	udigit high, low, q, r;


	udigit_NumberBase.digitmul(a, b, high, low);
	udigit_NumberBase.digitdiv(high, low, p, q, r);
	return r;
}
