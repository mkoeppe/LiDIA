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
//	Author	: Markus Maurer (MM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/base/base_bigmod.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



residue_class_list< bigint > * base_bigmod::L = NULL;



//
// normalize
//

void
base_bigmod::normalize (const bigint & m)
{
	if (I.is_zero() || I.is_one()) {
		// nothing to do, regardless of m (even if m is zero)
		return;
	}

	if (m.is_zero())
		lidia_error_handler ("base_bigmod::normalize",
				     "The modulus is set to zero!");

	bigint q, r;

	div_rem(q, r, I, m);
	I.assign(r);
	if (I.is_negative())
		add(I, I, m);
}



//
// modifiers
//

bigint
base_bigmod::invert(const bigint & m, int verbose)
{
	bigint u;
	bigint d(xgcd_left(u, I, m));


	if (d.is_one()) {
		I.assign(u);
		normalize(m);
	}
	else {
		if (verbose) {
			warning_handler("base_bigmod", "invert::inverse does not exist.");
		}
		else {
			lidia_error_handler("base_bigmod", "invert::inverse does not exist.");
		}
	}
	return d;
}



void
base_bigmod::multiply_by_2 (const bigint & m)
{
	I.multiply_by_2();
	if (I.compare(m) >= 0)
		subtract(I, I, m);
}



void
base_bigmod::divide_by_2 (const bigint & m)
{
	if (m.is_even() && !I.is_even())
		lidia_error_handler("base_bigmod", "divide_by_2(bigint&)::even modulus");
	if (!I.is_even())
		add(I, I, m);
	I.divide_by_2();
}



//
// arithmetic procedures
//

void
divide (base_bigmod & c, const base_bigmod & a, const base_bigmod & b, const bigint & m)
{
	bigint u;
	bigint d(xgcd_left(u, b.I, m));


	if (!d.is_one())
		lidia_error_handler("base_bigmod", "divide::inverse undefined");
	if (u.is_negative())
		add(u, u, m);
	multiply(c.I, a.I, u);
	remainder(c.I, c.I, m);
}



void
divide (base_bigmod & c, const base_bigmod & a, const bigint & b, const bigint & m)
{
	bigint u;
	bigint d(xgcd_left(u, b, m));


	if (!d.is_one())
		lidia_error_handler("base_bigmod", "operator/::inverse undefined");
	if (u.is_negative())
		add(u, u, m);
	multiply(c.I, a.I, u);
	remainder(c.I, c.I, m);
}



void
power (base_bigmod & c, const base_bigmod & a, const bigint & b, const bigint & m)
{
	base_bigmod multiplier = a;
	bigint exponent;


	if (b.is_zero() || a.is_one())
		c.assign_one();
	else {
		exponent.assign(b);
		if (exponent.is_negative()) {
			bigint d, u;
			d = xgcd_left(u, a.I, m);
			if (!d.is_one())
				lidia_error_handler("base_bigmod", "power::inverse undefined");
			if (u.is_negative())
				add(u, u, m);
			multiplier.I.assign(u);
			exponent.negate();
		}
		else
			multiplier.I.assign(a.I);
		c.assign_one();
		while (exponent.is_gt_zero()) {
			if (!exponent.is_even())
				multiply(c, c, multiplier, m);
			square(multiplier, multiplier, m);
			exponent.divide_by_2();
		}
	}
}



void
power (base_bigmod & c, const base_bigmod & a, long b, const bigint & m)
{
	base_bigmod multiplier;
	long exponent;


	if ((b == 0) || a.is_one())
		c.assign_one();
	else {
		exponent = b;
		if (exponent < 0) {
			bigint d, u;
			d = xgcd_left(u, a.I, m);
			if (!d.is_one())
				lidia_error_handler("base_bigmod", "operator^::inverse undefined");
			if (u.is_negative())
				add(u, u, m);
			multiplier.I.assign(u);
			exponent = -exponent;
		}
		else
			multiplier.I = a.I;
		c.assign_one();
		while (exponent > 0) {
			if (exponent & 1)
				multiply(c, c, multiplier, m);
			square(multiplier, multiplier, m);
			exponent >>= 1;
		}
	}
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
