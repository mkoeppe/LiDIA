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



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



bigfloat bigfloat::C_E;
bigfloat bigfloat::C_Catalan;
bigfloat bigfloat::C_Euler;
bigfloat bigfloat::C_Pi;



void
constant_Catalan (bigfloat & y)
{

	if (bigfloat::C_Catalan.length() >= bigfloat::digit_precision) {
		y.assign(bigfloat::C_Catalan);
		y.normalize();
	}
	else {
		bigfloat temp;
		long i = 1, j = 3;

		y.assign(0.5);
		bigfloat term(0.5);
		bigfloat p(0.5);

		while (!term.is_approx_zero()) {
			multiply(p.m, p.m, i);
			p.normalize();
			divide(p, p, j);
			multiply(term.m, term.m, i);
			term.normalize();
			divide(term, term, j);
			divide(temp, p, j);
			add(term, term, temp);
			add(y, y, term);
			i++;
			j += 2;
		}
		bigfloat::C_Catalan.assign(y);
	}
}



#ifndef USE_OLD_E

void
bigint_exp (long a, long b, bigint & first, bigint & second)
{
	bigint tmp1_first, tmp1_second;
	bigint tmp2_first, tmp2_second;
	bigint tmp;

	if (b == a + 1) {
		first.assign(b);
		second.assign(b);
		return;
	}

	bigint_exp(a, (a + b) >> 1, tmp1_first, tmp1_second);
	bigint_exp((a + b) >> 1, b, tmp2_first, tmp2_second);

	multiply(tmp, tmp1_first, tmp2_second);
	add(first, tmp, tmp2_first);
	multiply(second, tmp1_second, tmp2_second);
}



void
constant_E (bigfloat & x)
{
	static long lgamma_c[] = {40, 112, 289, 709, 1675, 3866, 8759, 19569,
				  43238, 94672, 205733, 444239, 954020,
				  2039119, 4340390, 9205076, 19458735,
				  41014632, 86223576, 180835770, 378448767,
				  790451976};

	if (bigfloat::C_E.length() >= bigfloat::digit_precision) {
		x.assign(bigfloat::C_E);
		x.normalize();
	}
	else {
		bigint ans_first;
		bigint ans_second;
		double y = 1.0 + static_cast<double>(bigfloat::binary_precision) * LOG2;
		long l, h;

		h = 16;
		l = 0;
		do {
			h <<= 1;
			l++;
		} while (lgamma_c[l] < y);
		bigint_exp(0, h, ans_first, ans_second);
		x.assign(ans_first, ans_second);
		bigfloat::C_E.assign(x);
	}
}



#else

void
constant_E (bigfloat & x)
{
	double f, t;
	long i, n;

	if (bigfloat::C_E.length() >= bigfloat::digit_precision) {
		x.assign(bigfloat::C_E);
		x.normalize();
	}
	else {
		n = 2;
		f = 0.0;
		t = bigfloat::binary_precision;
		while (f < t) {
			f += log(static_cast<double>(n)) * INVLOG2;
			n++;
		}
		n -= 2;
		x.assign(1 + n);
		bigfloat fac(n);

		for (i = 1; i < n; i++) {
			multiply(fac.m, fac.m, n - i);
			add(x.m, x.m, fac.m);
		}

		divide(x, x, fac);
		bigfloat::C_E.assign(x);
	}
}

#endif	// USE_OLD_E



void
constant_Euler (bigfloat & y)
{
	long l, n, k, x;
	bigfloat u, v, a, b, c;

	if (bigfloat::C_Euler.length() >= bigfloat::digit_precision) {
		y.assign(bigfloat::C_Euler);
		y.normalize();
	}
	else {
		l = bigfloat::digit_precision;

		x = 1 + static_cast<long>((0.25 * (l - 3)) * (bits_per_base_digit * LOG2));
		n = 1 + static_cast<long>(3.591 * x);

		a.assign(x);
		log(u, a);
		if (u.sign() > 0)
			u.negate();
		a.assign(u);
		b.assign_one();
		v.assign_one();

		for (k = 1; k <= n; k++) {
			multiply(b, b, x);
			multiply(b, b, x);
			divide(b, b, (k * k));
			multiply(a, a, x);
			multiply(a, a, x);
			divide(a, a, k);
			add(c, a, b);
			divide(a, c, k);
			add(u, u, a);
			add(v, v, b);
		}
		divide(y, u, v);
		bigfloat::C_Euler.assign(y);
	}
}



void
constant_Pi (bigfloat & x)
{
	long l, n, n1, hprec = bigfloat::binary_precision, prec;
	double alpha;

	const long k1 = 545140134;
	const long k2 = 13591409;
	const long k3 = 640320;
	const double alpha2 = (47.110413 / bits_per_base_digit);

	prec = bigfloat::digit_precision;
	if (bigfloat::C_Pi.length() >= prec) {
		x.assign(bigfloat::C_Pi);
		x.normalize();
	}
	else {
		bigfloat p1(k1), p2, p3, h;

		prec -= 2;
		n = static_cast<long>(1 + prec * (1 / alpha2));
		n1 = 6 * n - 1;

		multiply(h, p1, n);
		add(p2, h, k2);
		p1.assign(p2);

		if (prec >= 4)
			l = 4;
		else
			l = prec;
		alpha = l;
		p1.cut(l * bits_per_base_digit);

		while (n) {
			if (n > 1290) {
				if (n1 > 46340) {
					multiply(h, p1, n1 - 2);
					multiply(h, h, n1);
					multiply(h, h, n1 - 4);
					divide(h, h, (n * n));
					divide(p3, h, n);
				}
				else {
					multiply(h, p1, (n1 * (n1 - 2)));
					multiply(h, h, (n1 - 4));
					divide(h, h, (n * n));
					divide(p3, h, n);
				}
			}
			else {
				multiply(h, p1, (n1 * (n1 - 2)));
				multiply(h, h, (n1 - 4));
				divide(p3, h, (n * n * n));
			}
			divide(p3, p3, 100100025);
			divide(p3, p3, 327843840);
			subtract(p2, p2, k1);
			subtract(p1, p2, p3);
			alpha += alpha2;
			l = static_cast<long>(1 + alpha);
			if (l > prec)
				l = prec;
			p1.cut(l * bits_per_base_digit);
			n--;
			n1 -= 6;
		}
		divide(x, 53360, p1);
		h.assign(k3);
		sqrt(h, h);
		bigfloat::binary_precision = hprec;
		multiply(x, x, h);
		bigfloat::C_Pi.assign(x);
	}
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
