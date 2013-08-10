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
//	Author	: Thomas Papanikolaou (TP)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigrational.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// assigners
//

void
bigrational::assign (double d)
{
	int exponent, sign = (d < 0.0);
	double fraction;
	base_digit m;
	long ndigits = 0;

	num.assign_zero();
	den.assign_one();

	if (d) {
		if (sign)
			d = -d;

		fraction = std::frexp(d, &exponent);

		while (fraction != 0) {
			shift_left(num, num, bigint::bits_per_digit());
			fraction = std::ldexp(fraction, bigint::bits_per_digit());
			m = static_cast<base_digit>(fraction);
			fraction -= m;
			ndigits++;
			add(num, num, m);
		}
		// the product on the right side is never > MAXINT
		exponent -= static_cast<int>(ndigits * bigint::bits_per_digit());
		if (exponent < 0)
			shift_left(den, den, -exponent);
		else
			shift_left(num, num, exponent);
		if (sign)
			num.negate();
	}
	normalize();
}



//
// converters
//

double
dbl (const bigrational & a)
{
	if (a.numerator().is_zero()) {
		return 0.0;
	}

	long ln = a.numerator().bit_length();
	long ld = a.denominator().bit_length();
	long en = 0, ed = 0;

	if (ln > 1023) {
		en = ln - 1023;
	}
	if (ld > 1023) {
		ed = ld - 1023;
	}


	bigint an(a.numerator() >> en);
	bigint ad(a.denominator() >> ed);

	// an = (a/2^en)*2^en
	// ad = (a/2^ed)*2^ed

	double d = dbl(an) / dbl(ad);
	d = std::ldexp(d, static_cast<int>(en - ed));

	return d;
}



//
// modifiers
//

void
bigrational::normalize ()
{
	int s = den.sign();
	if (s == 0) {
		lidia_error_handler("bigrational", "normalize()::division by zero.");
		assign_zero();
		return;
	}
	if (s < 0) {
		num.negate();
		den.negate();
	}
	bigint g = bgcd(num, den);
	num /= g;
	den /= g;
}



void
bigrational::invert ()
{
	if (num.is_zero()) {
		lidia_error_handler("bigrational", "invert()::division by zero.");
		assign_zero();
		return;
	}
	num.swap(den);
	if (den.is_negative()) {
		num.negate();
		den.negate();
	}
}



//
// arithmetics
//

void
add (bigrational & c, const bigrational & a, const bigrational & b)
{
	if (a.is_zero()) {
		c.num.assign(b.num);
		c.den.assign(b.den);
	}
	else if (b.is_zero()) {
		c.num.assign(a.num);
		c.den.assign(a.den);
	}
	else {
		bigint g(gcd(a.den, b.den)), h;
		if (g.is_one()) {
			// c.num   a.num * b.den + a.den * b.num;
			// ----- = -----------------------------
			// c.den         a.den * b.den;

			multiply(h, a.num, b.den);
			multiply(c.num, a.den, b.num);
			add(c.num, c.num, h);
			multiply(c.den, a.den, b.den);
		}
		else {
			// bigint t = a.num * (b.den / g) + b.num * (a.den / g);
			// bigint h = gcd(t, g);
			// c.num = t / h;
			// c.den = (a.den / g) * (b.den / h);

			bigint t, s, ss;
			divide(s, b.den, g);
			multiply(t, a.num, s);
			divide(s, a.den, g);
			multiply(ss, s, b.num);
			add(t, t, ss);
			h = gcd(t, g);
			divide(c.num, t, h);
			divide(t, b.den, h);
			multiply(c.den, s, t);
		}
		if (c.den.is_negative()) {
			c.num.negate();
			c.den.negate();
		}
	}
}



void
subtract (bigrational & c, const bigrational & a, const bigrational & b)
{
	if (a.is_zero()) {
		c.num.assign(b.num);
		c.num.negate();
		c.den.assign(b.den);
	}
	else if (b.is_zero()) {
		c.num.assign(a.num);
		c.den.assign(a.den);
	}
	else {
		bigint g(gcd(a.den, b.den)), h;
		if (g.is_one()) {
			// c.num   a.num * b.den - a.den * b.num;
			// ----- = -----------------------------
			// c.den         a.den * b.den;

			multiply(h, a.num, b.den);
			multiply(g, a.den, b.num);
			subtract(c.num, h, g);
			multiply(c.den, a.den, b.den);
		}
		else {
			// bigint t = a.num * (b.den / g) - b.num * (a.den / g);
			// bigint h = gcd(t, g);
			// c.num = t / h;
			// c.den = (a.den / g) * (b.den / h);

			bigint t, s, ss;
			divide(s, b.den, g);
			multiply(t, a.num, s);
			divide(s, a.den, g);
			multiply(ss, s, b.num);
			subtract(t, t, ss);
			h = gcd(t, g);
			divide(c.num, t, h);
			divide(t, b.den, h);
			multiply(c.den, s, t);
		}
		if (c.den.is_negative()) {
			c.num.negate();
			c.den.negate();
		}
	}
}



void
multiply (bigrational & c, const bigrational & a, const bigrational & b)
{
	bigint g(gcd(a.num, b.den));
	bigint h(gcd(a.den, b.num));
	bigint s, t;

	switch (g.is_one() * 2 + h.is_one()) {
	case 0:
		// c.num   (a.num / g) * (b.num / h)
		// ----- = -------------------------
		// c.den   (a.den / h) * (b.den / g)
		divide(s, a.num, g);
		divide(t, b.num, h);
		multiply(c.num, s, t);

		divide(s, a.den, h);
		divide(t, b.den, g);
		multiply(c.den, s, t);
		break;
	case 1:
		// c.num   (a.num / g) * b.num
		// ----- = -------------------
		// c.den   a.den * (b.den / g)
		divide(s, a.num, g);
		multiply(c.num, s, b.num);
		divide(t, b.den, g);
		multiply(c.den, a.den, t);
		break;
	case 2:
		// c.num   a.num * (b.num / h)
		// ----- = -------------------
		// c.den   (a.den / h) * b.den
		divide(t, b.num, h);
		multiply(c.num, a.num, t);
		divide(s, a.den, h);
		multiply(c.den, s, b.den);
		break;
	case 3:
		// c.num   a.num * b.num
		// ----- = -------------
		// c.den   a.den * b.den
		multiply(c.num, a.num, b.num);
		multiply(c.den, a.den, b.den);
		break;
	}
	if (c.den.is_negative()) {
		c.num.negate();
		c.den.negate();
	}
}



void
divide (bigrational & c, const bigrational & a, const bigrational & b)
{
	bigint g(gcd(a.num, b.num));
	bigint h(gcd(a.den, b.den));
	bigint s, t;

	switch (g.is_one() * 2 + h.is_one()) {
	case 0:
		// c.num   (a.num / g) * (b.den / h)
		// ----- = -------------------------
		// c.den   (a.den / h) * (b.num / g)
		divide(s, a.num, g);
		divide(t, b.den, h);
		multiply(c.num, s, t);

		divide(s, a.den, h);
		divide(t, b.num, g);
		multiply(c.den, s, t);
		break;
	case 1:
		// c.num   (a.num / g) * b.den
		// ----- = -------------------
		// c.den   a.den * (b.num / g)
		divide(s, a.num, g);
		multiply(c.num, s, b.den);
		divide(t, b.num, g);
		multiply(c.den, a.den, t);
		break;
	case 2:
		// c.num   a.num * (b.den / h)
		// ----- = -------------------
		// c.den   (a.den / h) * b.num
		divide(t, b.den, h);
		multiply(c.num, a.num, t);
		divide(s, a.den, h);
		multiply(c.den, s, b.num);
		break;
	case 3:
		// c.num   a.num * b.den
		// ----- = -------------
		// c.den   a.den * b.num
		multiply(c.num, a.num, b.den);
		multiply(c.den, a.den, b.num);
		break;
	}
	if (c.den.is_negative()) {
		c.num.negate();
		c.den.negate();
	}
}



void
shift_left (bigrational & c, const bigrational & a, long ui)
{
	if (ui < 0) {
		lidia_error_handler("bigrational", "shift_left()::index is negative.");
		c.num.assign(a.num);
		c.den.assign(a.den);
		return;
	}
	c.den.assign(a.den);
	shift_left(c.num, a.num, ui);
	c.normalize();
}



void
shift_right (bigrational & c, const bigrational & a, long ui)
{
	if (ui < 0) {
		lidia_error_handler("bigrational", "shift_right()::index is negative.");
		c.num.assign(a.num);
		c.den.assign(a.den);
		return;
	}
	c.num.assign(a.num);
	shift_left(c.den, a.den, ui);
	c.normalize();
}



void
power (bigrational & c, const bigrational & a, const bigint & b)
{
	bigint n(1UL), d(1UL);

	if (!b.is_zero()) {
		if (b.is_gt_zero()) {
			power(n, a.numerator(), b);
			power(d, a.denominator(), b);
		}
		else {
			bigint abs_b(-b);

			power(n, a.denominator(), abs_b);
			power(d, a.numerator(), abs_b);
			if (d.is_negative()) {
				n.negate();
				d.negate();
			}
		}
	}
	c.assign_numerator(n);
	c.assign_denominator(d);
}



void
power (bigrational & c, const bigrational & a, long i)
{
	bigint n(1UL), d(1UL);

	if (i) {
		if (i > 0) {
			power(n, a.numerator(), i);
			power(d, a.denominator(), i);
		}
		else {
			i = -i;
			power(n, a.denominator(), i);
			power(d, a.numerator(), i);
			if (d.is_negative()) {
				n.negate();
				d.negate();
			}
		}
	}
	c.assign_numerator(n);
	c.assign_denominator(d);
}



bigint
round (const bigrational & a)
{
	bigrational c;
	bigint rn(a.numerator()), rd(a.denominator());

	rn.multiply_by_2();
	add(rn, rn, rd);
	rd.multiply_by_2();
	c.assign_numerator(rn);
	c.assign_denominator(rd);
	return floor(c);
}



bigint
floor (const bigrational & a)
{
	bigint q, r;

	div_rem(q, r, abs(a.numerator()), abs(a.denominator()));
	if (a.sign() < 0 && r.sign() != 0) {
		if (a.sign() < 0)
			q.negate();
		dec(q);
	}
	else {
		if (a.sign() < 0)
			q.negate();
	}
	return q;
}



bigint
ceiling (const bigrational & a)
{
	bigint q, r;

	div_rem(q, r, abs(a.numerator()), abs(a.denominator()));
	if (a.sign() >= 0 && r.sign() != 0)
		inc(q);
	if (a.sign() < 0)
		q.negate();
	return q;
}



bigint
truncate (const bigrational & a)
{
	bigint q, r;

	div_rem(q, r, abs(a.numerator()), abs(a.denominator()));
	if (a.sign() < 0)
		q.negate();
	return q;
}



//
// non friend functions that use bigrationals
//

bool
square_root (bigrational & root, const bigrational & s)
{
	if (s.is_lt_zero()) {
		return false;
	}

	bigint ns, ds;

	sqrt(ns, s.numerator());
	sqrt(ds, s.denominator());
	root = bigrational(ns, ds);
	return (root*root == s);
}



bool
cube_root (bigrational & root, const bigrational & s)
{
	bigint ns(s.numerator()), ds;

	if (ns.is_lt_zero()) {
		ns = -ns;
	}
	newton_root(ns, ns, 3);
	newton_root(ds, s.denominator(), 3);
	root = bigrational(ns, ds);
	if (s.is_lt_zero()) {
		root.negate();
	}
	return ((root*root*root) == s);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
