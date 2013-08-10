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
    sum = a + b;  
    udigit new_carry = (sum < a) ? 1 : 0;

    sum += carry;
    if(sum < carry) {
	++ new_carry;
    }

    return new_carry;
}



inline udigit
subtract (udigit & diff, udigit a, udigit b, udigit carry)
        //-RV * 2^base + diff = a - b - carry
{
    udigit new_carry = (b > a) ? 1 : 0;
    diff = a - b;

    if(carry > diff) {
	++new_carry;
    }
    diff -= carry;

    return new_carry;
}



inline udigit
multiply (udigit & prod, udigit a, udigit b)
        //RV * 2^base + prod = a * b
{
    unsigned long shifts = intDsize / 2;
    udigit const threshold = udigit(1) << shifts;
    udigit const lowBits = threshold - 1;

    udigit a0 = a & lowBits;
    udigit a1 = a >> shifts;
    udigit b0 = b & lowBits;
    udigit b1 = b >> shifts;

    // (a0 + a1*2^{N/2})(b0+b1*2^{N/2}) =
    //   a0*b0 + (a0*b1 + a1*b0)*2^{N/2} + a1*b1*2^N
    prod = a0 * b0;
    udigit tmp;
    udigit highDigit = add(tmp, a0 * b1, a1 * b0, 0);

    // highDigit won't overflow!
    highDigit <<= shifts;
    highDigit += (tmp >> shifts);
    highDigit += add(prod, prod, tmp << shifts, 0);
    highDigit += a1*b1;

    return highDigit;
}



inline udigit
divide (udigit & quot, udigit a1, udigit a0, udigit b)
	//a1 * 2^base + a0 = quot * b + RV
	//a1 must be less than b
{
    cln::cl_I a(a1);
    a = (a << intDsize) + a0;

    cln::cl_I_div_t q_and_rem = truncate2(a, b);
    quot = cl_I_to_ulong(q_and_rem.quotient);
    return cl_I_to_ulong(q_and_rem.remainder);
}


//
// modular functions, we assume that input is reduced modulo modul
//

inline udigit
add_mod (udigit a, udigit b, udigit p)
	//RV = (a + b) % p
	//assumes a,b < p
{
    udigit sum;
    udigit carry = add(sum, a, b, udigit(0));

    if(carry) { // ==> sum > p
	subtract(sum, sum, p, udigit(0));
    }
    else if(p <= sum) {
	sum -= p;
    }

    return sum;
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
    cln::cl_I prod = cln::cl_I(a) * cln::cl_I(b);
    return cl_I_to_ulong(mod(prod, cln::cl_I(p)));
}
