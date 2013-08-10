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
//	Author	:
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigfloat.h"
#include	"LiDIA/bigfloat_config.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



void
acos (bigfloat & y, const bigfloat & x)
{
	bigfloat p1(x), p2(1), p3;

	if (p1.compare(p2) == 0) {
		y.assign_zero();
		return;
	}

	p2.negate();
	if (p1.compare(p2) == 0) {
		constant_Pi(y);
		return;
	}

	multiply(p2, x, p1);
	dec(p2);
	p2.negate();
	if (p2.is_lt_zero())
		lidia_error_handler("bigfloat", "acos::cannot handle arguments > 1");
	sqrt(p2, p2);
	divide(p1, x, p2);
	atan(p3, p1);
	constant_Pi(p2);
	p2.e--;
	subtract(y, p2, p3);
}



void
acot (bigfloat & y, const bigfloat & x)
{
	bigfloat p1, p2;
	atan(p1, x);
	constant_Pi(p2);
	p2.e--;
	subtract(y, p2, p1);
}



void
asin (bigfloat & y, const bigfloat & x)
{
	bigfloat p1(x), p2(1);

	if (p1.compare(p2) == 0) {
		constant_Pi(y);
		y.e--;
		return;
	}

	p2.negate();
	if (p1.compare(p2) == 0) {
		constant_Pi(y);
		y.e--;
		y.negate();
		return;
	}

	multiply(p2, x, p1);
	dec(p2);
	p2.negate();
	if (p2.is_lt_zero())
		lidia_error_handler("bigfloat", "asin::cannot handle arguments > 1");
	sqrt(p2, p2);
	divide(p1, x, p2);
	atan(y, p1);
}



void
atan (bigfloat & y, const bigfloat & x)
{
	long i, j, ex, t = bigfloat::binary_precision;
	long m, u, f;

	if (x.is_zero()) {
		y.assign_zero();
		return;
	}

	ex = x.e + x.prec;
	if (ex > bigfloat::binary_precision)
		bigfloat::binary_precision = ex + bits_per_base_digit
			- ex % bits_per_base_digit;
	bigfloat a(1);
	bigfloat tmp(x);

	m = 0;
	if (tmp.is_lt_zero()) {
		m = 1;
		tmp.negate();
	}

	if (tmp.compare(a) == 0) {
		constant_Pi(y);
		y.e -= 2;
		if (m)
			y.negate();
		return;
	}

	tmp.normalize();

	ex = tmp.e + bigfloat::binary_precision;
	u = 0;
	if (ex > 0) {
		divide(tmp, 1, tmp);
		u = 1;
	}

	ex = tmp.e + bigfloat::binary_precision;
	f = 0;
	bigfloat q(tmp);
	if (ex > -10)
		while (tmp.e + bigfloat::binary_precision > -10) {
			multiply(q, q, tmp);
			add(q, q, a);
			sqrt(q, q);
			add(q, q, a);
			divide(tmp, tmp, q);
			q.assign(tmp);
			f++;
		}
	square(a, tmp);

	ex = tmp.e + bigfloat::binary_precision;
	if (ex < 0)
		ex = -ex;
	ex <<= 1;
	j = bigfloat::binary_precision / ex;
	if (j & 1)
		j++;
	divide(y, 1, 2 * j + 1);
	bigfloat::binary_precision = 4 * ex;

	for (i = j; i >= 1; i--) {
		multiply(y, y, a);
		divide(q, 1, 2 * i - 1);
		bigfloat::binary_precision += 2 * ex;
		if (bigfloat::binary_precision > t)
			bigfloat::binary_precision = t;
		y.negate();
		add(y, y, q);
	}
	bigfloat::binary_precision = t;
	multiply(y, y, tmp);
	y.e += f;

	if (u) {
		constant_Pi(a);
		a.e--;
		subtract(y, y, a);
		y.negate();
	}
	if (m)
		y.negate();

}



void
atan2 (bigfloat & z, const bigfloat & y, const bigfloat & x)
{
	bigfloat w;
	int ys = y.sign(), xs = x.sign();
	int yss = (ys < 0), xss = (xs < 0);
	char code = yss + (xss << 1);

	if (xs == 0) {
		if (ys == 0)
			z.assign_zero();
		else {
			constant_Pi(z);
			z.divide_by_2();
			if (ys < 0)
				z.negate();
		}
		return;
	}

	if (ys == 0) {
		if (xs < 0)
			constant_Pi(z);
		else
			z.assign_zero();
		return;
	}

	switch(code) {
	case 0:
	case 1:
		w.assign_zero();
		break;
	case 2:
		constant_Pi(w);
		break;
	case 3:
		constant_Pi(w);
		w.negate();
		break;
	}

	z.assign(y);
	divide(z, z, x);
	atan(z, z);
	add(z, z, w);
}



void
cos (bigfloat & y, const bigfloat & x)
{
	bigfloat tmp;
	double d, df;
	long ex, f = 0, s = 0, t = bigfloat::digit_precision;
	long i, k, n, s1, t0, t1, t2;

	//
	// test for zero
	//

	if (x.is_zero()) {
		y.assign_one();
		return;
	}

	//
	// p0 = x; increase the precision to bigfloat::binary_precision + log2(p0)
	//

	bigfloat p0(x);
	ex = p0.e + p0.prec;
	if (ex > 0) {
		bigfloat::binary_precision += ex + bits_per_base_digit
			- (ex % bits_per_base_digit);
		bigfloat::digit_precision = bigfloat::binary_precision / bits_per_base_digit;
	}

	//
	// normalize p0 and set p0 = abs(p0)
	//

	p0.normalize();
	p0.absolute_value();


	//
	// calculate pi and reduce p0 according to the equation
	//     cos(p0) = cos(p0' + k pi) = (-1)^k cos(p0')
	//

	constant_Pi(y);
	if (p0.compare(bigfloat::C_Pi) >= 0) {
		divide(y, p0, bigfloat::C_Pi);
		truncate(y.m, y);
		y.e = 0;
		s = y.m.is_odd();
		y.normalize();
		multiply(y, y, bigfloat::C_Pi);
		subtract(p0, p0, y);
	}

	//
	// calculate ex = log2(p0) exactly (as double); the calculate the
	// reduction range k and find the number f of reduction steps needed
	// so that
	//     p0 * 2 ^ f < p0 * 2 ^ k = p0 * 2 ^ (-sqrt(t/4)).
	// Then reduce p0.
	//

	ex = static_cast<long>(std::log(std::fabs(p0.m.most_significant_digit() / max_base_digit_2)) *
			       INVLOG2 + (p0.e + p0.length() * bits_per_base_digit));
	k = static_cast<long>(-std::sqrt(bigfloat::binary_precision / 4.0));
	if (ex > k) {
		n = static_cast<long>(1 + std::sqrt(4.0 * bigfloat::binary_precision));
		f = 1 + ex - k;
		p0.e -= f;
	}
	else {
		n = static_cast<long>(1 + (bigfloat::binary_precision * 1.0) / ((ex > 0) ? ex : -ex));
		f = 0;
	}

	//
	// reduce further the number of steps using a fix-point iteration
	//

	d = 0.0;
	for (i = 2; i <= n; i++)
		d += std::log(static_cast<double>(i)) * INVLOG2;
	i = 0;
	df = i * k + d;
	while (df > 0) {
		i++;
		d -= std::log(static_cast<double>(n)) * INVLOG2;
		df = i * k + d;
		n--;
	}
	n++;

	//
	// make n even
	//

	if (n & 1)
		n++;

	//
	// Set p0 = -p0^2, y = 1; then evaluate the taylor series for
	//           cos(x)-1 = -x^2/2! + x^4/4! -+ ...
	// using increasing precision at each step
	//

	square(p0, p0);
	p0.negate();
	y.assign(1);

	//
	// starting precision
	//

	t1 = ((-k) / bits_per_base_digit) << 1;
	t2 = bigfloat::digit_precision - 2;
	s1 = 0;
	bigfloat::binary_precision = (t1 + 2) * bits_per_base_digit;

	//
	// Use increasing precision at each step. l1 is the number
	// of correct Digits at the current point. l2 ist the precision
	// to be set before the next iteration occurs. The code for precision
	// increasing is adapted from the PARI system.
	//

	for (i = n; i > 2; i -= 2) {
		divide(y, y, i * (i - 1));
		if (p0.length() > t1) {
			tmp.cut(p0, t1 * bits_per_base_digit);
			multiply(y, y, tmp);
		}
		else
			multiply(y, y, p0);
		ex = y.e + bigfloat::binary_precision;

		// GP/PARI : Start
		s1 = s1 - ex;
		t0 = s1 / bits_per_base_digit;
		t1 += t0;
		if (t1 > t2)
			t1 = t2;
		s1 %= bits_per_base_digit;
		// GP/PARI : End

		bigfloat::binary_precision = (t1 + 2) * bits_per_base_digit;
		inc(y);
	}
	divide(y, y, i * (i - 1));
	multiply(y, y, p0);

	//
	// we have calculated y = cos(x/2^f) - 1. We will calculate now
	// cos(x/2^(f-1)) - 1 = 2 * y * (y + 2). y is exact, (y+2) is exact
	// so the final result (cos(x)-1) is exact. Error can only occur if
	// x was very near to pi/2.
	//

	if (f) {
		add(p0, y, 2);
		for (i = 0; i < f; i++) {
			multiply(y, y, p0);
			y.e++;
			add(p0, y, 2);
		}
	}

	//
	// Before adding 1 we check the result for very small y and increase
	// the precision accordingly so that the final result after the addition
	// is correct at bigfloat::binary_precision + abs(log2(y)) digits
	//

	ex = y.e;
	if (ex < 0)
		ex = -ex;
	if (ex > bigfloat::binary_precision) {
		bigfloat::binary_precision = ex + bits_per_base_digit
			- ex % bits_per_base_digit;
		bigfloat::digit_precision = bigfloat::binary_precision / bits_per_base_digit;
	}
	inc(y);

	//
	// correct the sign and normalize the answer
	//

	if (s)
		y.negate();

	//
	// set the old precision
	//

	bigfloat::digit_precision = t;
	bigfloat::binary_precision = t * bits_per_base_digit;
}



void
cot (bigfloat & y, const bigfloat & x)
{
	long s = 0, ex, t = bigfloat::digit_precision;

	ex = x.e + x.prec;
	if (ex < 0) {
		ex = -ex;
		bigfloat::binary_precision += 4 * (ex + bits_per_base_digit
						   - ex % bits_per_base_digit);
		bigfloat::digit_precision = bigfloat::binary_precision / bits_per_base_digit;
	}
	bigfloat p0(x);
	if (p0.is_lt_zero()) {
		s = 1;
		p0.negate();
	}
	constant_Pi(y);
	if (p0.compare(bigfloat::C_Pi) >= 0) {
		divide(y, p0, bigfloat::C_Pi);
		truncate(y, y);
		multiply(y, y, bigfloat::C_Pi);
		subtract(p0, p0, y);
	}
	cos(y, p0);
	square(p0, y);
	p0.negate();
	inc(p0);
	sqrt(p0, p0);
	divide(y, y, p0);
	if (s)
		y.negate();
	bigfloat::digit_precision = t;
	bigfloat::binary_precision = t * bits_per_base_digit;
}



void
sin (bigfloat & y, const bigfloat & x)
{
	bigfloat p0;
	long ex, t = bigfloat::digit_precision;

	ex = x.e + x.prec;
	if (ex <= 0) {
		int xs = x.sign();
		cos(y, x);
		p0.assign(y);
		inc(p0);

		if (bigfloat::binary_precision < (y.length() * bits_per_base_digit))
			bigfloat::binary_precision = y.length() * bits_per_base_digit;
		ex = -ex;
		bigfloat::binary_precision += 2 * (ex + bits_per_base_digit
						   - (ex % bits_per_base_digit));
		bigfloat::digit_precision = bigfloat::binary_precision / bits_per_base_digit;
		dec(y);
		y.negate();
		bigfloat::digit_precision = t;
		bigfloat::binary_precision = t * bits_per_base_digit;

		multiply(y, y, p0);
		sqrt(y, y);
		if (xs < 0)
			y.negate();
		return;
	}

	constant_Pi(p0);
	p0.e--;
	add(y, x, p0);
	cos(y, y);
	y.negate();

	bigfloat::digit_precision = t;
	bigfloat::binary_precision = t * bits_per_base_digit;
}



void
tan (bigfloat & y, const bigfloat & x)
{
	long s = 0, ex, t = bigfloat::digit_precision;

	ex = x.e + x.prec;
	if (ex < 0) {
		ex = -ex;
		bigfloat::binary_precision += 4 * (ex + bits_per_base_digit
						   - ex % bits_per_base_digit);
		bigfloat::digit_precision = bigfloat::binary_precision / bits_per_base_digit;
	}
	bigfloat p0(x);
	if (p0.is_lt_zero()) {
		s = 1;
		p0.negate();
	}
	constant_Pi(y);
	if (p0.compare(bigfloat::C_Pi) >= 0) {
		divide(y, p0, bigfloat::C_Pi);
		truncate(y, y);
		multiply(y, y, bigfloat::C_Pi);
		subtract(p0, p0, y);
	}
	cos(p0, p0);
	square(y, p0);
	y.negate();
	inc(y);
	sqrt(y, y);
	divide(y, y, p0);
	if (s)
		y.negate();
	bigfloat::digit_precision = t;
	bigfloat::binary_precision = t * bits_per_base_digit;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
