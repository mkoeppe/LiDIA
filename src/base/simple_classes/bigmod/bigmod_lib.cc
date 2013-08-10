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
//	Author	: Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigmod.h"
#include	"LiDIA/rational_factorization.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//**** get generator of finite field *********************************

bigmod get_generator(const bigmod &x)
{
	bigmod a, apow, ahelp;
	bigint *v, h1, gord, h;
	bool found = false;
	unsigned int i, length;
	rational_factorization rf;

	if (!is_prime(x.modulus()))
		lidia_error_handler("file bigmod_lib.c::get_generator",
				    "no prime modulus");

	gord = x.modulus(); // order of mult. group
	dec(gord);

	rf = factor_completely(gord);
	length = rf.no_of_comp();
	v = new bigint[length];

	// compute in v the difference between ord/p_i for prime factors p_i
	// of ord

	divide(v[0], gord, rf.base(length-1));
	h.assign(v[0]);
	i = 1;

	while (i < length) {
		divide(h1, gord, rf.base(length-1-i));
		subtract(v[i++], h1, h);
		h.assign(h1);
	}

	while (found == false)	{
		// check for generator
		do {
			a.randomize();
		} while (a.is_zero() || a.is_one());

		power(apow, a, v[0]);
		if (apow.is_one())
			continue;

		i = 1;
		while (i < length) {
			power(ahelp, a, v[i++]);
			multiply(apow, apow, ahelp);
			if (apow.is_one())
				break;
		}
		if ((i == length) && (apow.is_one() == false))
			found = true;
	}
	delete[] v;

	return a;
}



//**** hash functions ************************************************

udigit hash(const bigmod & x)
{
	return udigit(x.mantissa().least_significant_digit());
}



//**** is_square *****************************************************

bool is_square(const bigmod & a)
{
	if (!is_prime(a.modulus()))
		lidia_error_handler("file bigmod_lib.c::is_square",
				    "no prime modulus");

	return (jacobi(a.mantissa(), a.modulus()) != -1);
}



//**** sqrt **********************************************************

bigmod sqrt(const bigmod & a)
{
	bigmod t;
	bigmod res1, res0, h, h2;
	bigint exp = a.modulus();
	unsigned int i;

	if (a.is_zero())
		return a;

	do {
		t.randomize();
	} while (is_square(t*t-4*a));

	inc(exp); shift_right(exp, exp, 1);

	res1.assign_one(); res0.assign_zero();

	for (i = exp.bit_length()-2; i > 0; i--) {
		// fast exponentiation left-right
		h.assign(res0);
		square(h2, h);
		square(res0, res1);
		multiply(res0, res0, a);
		subtract(res0, h2, res0); // res0 = h*h - a*res1*res1;

		h.multiply_by_2();
		multiply(h2, res1, t);
		add(h, h, h2);
		multiply(res1, res1, h); // res1 = res1*(res1*t + 2*h);

		if (exp.bit(i)) {
			h.assign(res0);
			multiply(res0, res1, a);
			negate(res0, res0); // res0 = -res1*a;

			multiply(res1, res1, t);
			add(res1, res1, h); // res1 = res1*t + h;
		}
	}

	h.assign(res0);
	square(h2, h);
	square(res0, res1);
	multiply(res0, res0, a);
	subtract(res0, h2, res0); // res0 = h*h - a*res1*res1;

	h.multiply_by_2();
	multiply(h2, res1, t);
	add(h, h, h2);
	multiply(res1, res1, h); // res1 = res1*(res1*t + 2*h);

	if (exp.is_odd()) {
		multiply(res0, res1, a);
		negate(res0, res0); // res0 = -res1*a;
	}

	return (res0);
}



//***** solve_quadratic **********************************************
//  Input: polynomial X^2 + a1*X + a0
//  Output: true and set root, if polynomial splits mod p
//          false otherwise
//********************************************************************

bool solve_quadratic(bigmod & root, const bigmod & a1, const bigmod & a0)
{
	bigmod delta;

	if (a0.is_zero()) {
		root.assign_zero();
		return true;
	}

	square(delta, a1);
	subtract(delta, delta, 4*a0);

	if (!is_square(delta))
		return false;

	root = (-a1 + sqrt(delta))/2;
	return true;
}



//**** characteristic ************************************************

bigint characteristic(const bigmod & a)
{
	if (is_prime(a.modulus()))
		return a.modulus();
	else {
		lidia_error_handler("file bigmod_lib.c::characteristic",
				    "no prime modulus");
		return 0;
	}
}



//**** number_of_elements ********************************************

bigint number_of_elements(const bigmod & a)
{
	if (is_prime(a.modulus()))
		return a.modulus();
	else {
		lidia_error_handler("file bigmod_lib.c::number_of_elements",
				    "no prime modulus");
		return 0;
	}
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
