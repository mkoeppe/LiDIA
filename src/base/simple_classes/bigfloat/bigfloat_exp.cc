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
//==============================================================================================



#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigfloat.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



void
exp (bigfloat & y, const bigfloat & x)
{
	long s, f = 0, i, ex, n, t = bigfloat::binary_precision - rounding_bits;
	double a, b;
	double reg;
	long l0, l1, l2, s1, ex1;

	//
	// if abs(x) < 2^(-bigfloat::binary_precision) then x = 1 to full precision
	//

	if (x.is_zero()) {
		y.assign_one();
		return;
	}

	//
	// save the sign of x and set y = abs(x)
	//

	s = x.sign();
	y.assign(x);
	if (s < 0)
		y.negate();

	//
	// remove leading zeros and calculate the length of the integer
	// part of y
	//

	ex1 = y.e + y.length() * bits_per_base_digit;

	a = static_cast<double>(y.m.most_significant_digit());
	a /= max_base_digit_2;
	a = std::log(a) * INVLOG2 + ex1;

	//
	// if x was a small number choose a better convergence domain
	//

	if (y.length() > 2)
		reg = 1.246294;
	else
		reg = 5.0;

	//
	// increase the precision to avoid cancelation; calculate
	// the number f of reduction steps needed so that
	// y / (2 ^ f) < 2 ^ (-sqrt(t / reg)). Reduce y.
	//

	t = t + static_cast<long>(std::sqrt(static_cast<double>(t))) - (ex1 < 0) * ex1;
	b = -std::sqrt(t / reg);
	if (a > b) {
		n = 1 - static_cast<long>(reg * b);
		f = 1 + static_cast<long>(a - b);
		y.e -= f;
	}
	else
		n = 1 - static_cast<long>(t / a);

	//
	// begin of the main calculation loop
	//

	bigfloat tmp(1), tmp1;

	//
	// Use increasing precision at each step. l1 is the number
	// of correct Digits at the current point. l2 ist the precision
	// to be set before the next iteration occurs. The code for precision
	// increasing is adapted from the PARI system.
	//

	l1 = 2;
	l2 = 1 + (t / bits_per_base_digit);
	s1 = 0;
	bigfloat::binary_precision = (l1 + 2) * bits_per_base_digit;

	for (i = n; i; i--) {
		divide(tmp, tmp, i);
		if (y.length() > l1) {
			tmp1.cut(y, l1 * bits_per_base_digit);
			multiply(tmp, tmp, tmp1);
		}
		else
			multiply(tmp, tmp, y);

		// GP/PARI : start

		ex = tmp.e + bigfloat::binary_precision;
		s1 = s1 - ex;
		l0 = s1 / bits_per_base_digit;
		l1 += l0;
		if (l1 > l2)
			l1 = l2;
		s1 %= bits_per_base_digit;

		// GP/PARI : end

		bigfloat::binary_precision = (l1 + 2) * bits_per_base_digit;
		inc(tmp);
	}

	//
	// reconstruct the result by repeated squaring
	//

	if (f) {
		while (f--)
			square(tmp, tmp);
	}

	//
	// if x was negativ we must invert the result
	//

	if (s < 0) {
		divide(y, 1, tmp);
	}
	else {
		y.assign(tmp);
		y.normalize();
	}

	//
	// set the old precision
	//

	bigfloat::binary_precision = bigfloat::digit_precision * bits_per_base_digit;
}



void
log (bigfloat & y, const bigfloat & x)
{
	long i;
	long dt = bigfloat::digit_precision;
	long t = bigfloat::binary_precision;

	if (x.is_le_zero())
		lidia_error_handler("bigfloat", "log()::argument is <= 0");

	//
	// if (x - 1) == 0 then the result is 0
	//

	bigfloat pr(x);
	dec(pr);

	long ex1 = pr.e + pr.prec; // , ex2;

	if (pr.is_zero()) {
		y.assign_zero();
		return;
	}

	//
	// set y = x; if expo(x-1) < 0 increase the precision to
	// bigfloat::binary_precision + abs(log2(x-1)) to avoid cancellation when
	// summing the series for arctanh()
	//

	y.assign(x);
	if (ex1 < 0) {
		ex1 = -ex1;
		ex1 = ex1 + bits_per_base_digit - ex1 % bits_per_base_digit;
		t = (bigfloat::binary_precision += ex1);
		bigfloat::digit_precision += (ex1 / bits_per_base_digit);
	}
#if 0
	else if ((ex2 = y.e + y.prec) > bigfloat::binary_precision) {
		ex2 = ex2 + bits_per_base_digit - ex2 % bits_per_base_digit;
		t = (bigfloat::binary_precision += ex2);
		bigfloat::digit_precision += (ex2 / bits_per_base_digit);
	}
#endif

	//
	// use the identity log(1/y) = -log(y) to make y > 0
	//

	long s = 0;
	if ((y.e + y.prec) < 1) {
		s = 1;
		divide(y, 1, x);
	}

	//
	// first reduction using the identity log(x) = 2 * log(sqrt(x))
	//

	long f1 = 0;
	while ((y.e + y.prec) > 1) {
		f1++;
		sqrt(y, y);
	}

	//
	// since we are calculating the arctanh() we calculate the number of
	// reduction/calculation steps (f/n) (using the above identity) through (y-1)
	//

	long f = 0, n;

	add(pr, y, -1);

	double xc = pr.m.most_significant_digit() / max_base_digit_2;
	xc = std::log(xc) * INVLOG2 + (pr.e + t);

	long l = bigfloat::digit_precision;
	double aa = -xc, bb = 2.75 * std::sqrt(l / 3.0);

	if (aa <= bb) {
		n = 1 + static_cast<long>((bits_per_base_digit >> 1) / 2.75 * std::sqrt(3.0 * l));
		f = 1 + static_cast<long>(bb - aa);
	}
	else {
		double a = aa * LOG2, b = (bits_per_base_digit >> 1) * l * LOG2;
		n = 2 + static_cast<long>(b / a);
	}

	//
	// second reduction using the identity log(x) = 2 * log(sqrt(x))
	//

	for (i = 0; i < f; i++)
		sqrt(y, y);

	//
	// now calculate log((y - 1) / (y + 1))
	//

	dec(y);
	add(pr, y, 2);
	divide(y, y, pr);

	bigfloat sum(y);
	square(y, y);

	long ex = y.e + y.prec;

	//
	// Use increasing precision at each step. l1 is the number
	// of correct Digits at the current point. l2 ist the precision
	// to be set before the next iteration occurs. The code for precision
	// increasing is adapted from the PARI system.
	//

	bigfloat term, tmp;

	n = (n << 1) + 1;
	long l0, l1 = (f + f1) / (bits_per_base_digit >> 1) + 1;
	long hlp = 0;
	long l2 = bigfloat::digit_precision - 3 + (f / bits_per_base_digit);
	bigfloat::binary_precision = (l1 + 2) * bits_per_base_digit;
	divide(pr, pr, n);
	n -= 2;

	for (i = n; i >= 1; i -= 2) {
		if (y.length() > l1) {
			tmp.cut(y, l1 * bits_per_base_digit);
			multiply(pr, pr, tmp);
		}
		else
			multiply(pr, pr, y);
		divide(term, 1, i);
		// GP/PARI : start
		hlp -= ex;
		l0 = hlp / bits_per_base_digit;
		l1 += l0;
		if (l1 > l2)
			l1 = l2;
		hlp %= bits_per_base_digit;
		// GP/PARI : end
		bigfloat::binary_precision = (l1 + 2) * bits_per_base_digit;
		add(pr, pr, term);
	}

	//
	// restore precision and reconstruct (shift); set the sign;
	// delete local variables
	//

	bigfloat::digit_precision = dt;
	bigfloat::binary_precision = dt * bits_per_base_digit;
	multiply(y, sum, pr);
	y.e += 1 + f + f1;
	if (s)
		y.negate();
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
