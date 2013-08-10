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
#include	"LiDIA/bigfloat.h"
#include	"LiDIA/base/b_value.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//#define BIGFLOAT_EXP_DEBUG


int  bigfloat::old_rnd_mode;
long bigfloat::old_bin_prec;



void
sqrt (bigfloat & y, const bigfloat & x, long k)
{
	if (x.is_lt_zero()) {
		lidia_error_handler("bigfloat", "sqrt()::argument must be >= 0");
		y.assign_zero();
		return;
	}
	if (x.is_zero())
		y.assign_zero();
	else {
		if (k < 1) {

			//lidia_error_handler("bigfloat::sqrt(bigfloat, bigfloat, long)",
			//"k must be >= 1");
			//return;
			k = 1;
		}

		int old_rnd_mode = bigfloat::get_mode();

		long eps, ex, i, n, q, t;
		bigfloat p1, p2, p3;
		bigint   tmp1;

		// Normalization is important to make sure
		// that the most significant digit of the
		// mantissa only contains significant bits
		// and to guarantee that length() * bits_per_digit
		// equals bit_length() below.

		p1.assign(x);
		p1.base_digit_normalize();

#ifdef BIGFLOAT_SQRT_DEBUG
		std::cout << "x = " << x << std::endl;
		std::cout << "x.exponent() = " << x.exponent() << std::endl;
		std::cout << "x.mantissa() = " << x.mantissa() << std::endl;
		std::cout << "p1 = " << p1 << std::endl;
		std::cout << "p1.exponent() = " << p1.exponent() << std::endl;
		std::cout << "p1.mantissa() = " << p1.mantissa() << std::endl;
#endif

		if (k < 3)
			t = 3;
		else
			t = k+1;


		//
		// calculate the exponent ex of the result. If ex is odd
		// multiply the first approximation by 2
		//
		// x = 0.m * 2^eps * 2^(2ex), eps = 0, 1
		//

		ex = p1.e + p1.prec;
		eps = ex % 2; // eps is b
		ex /= 2;

		if (eps < 0) {
			eps += 2;
			ex--;
		}

#ifdef BIGFLOAT_SQRT_DEBUG
		std::cout << "exponent/2 = " << ex << std::endl;
		std::cout << "exponent%2 = " << eps << std::endl;
#endif

		//
		// Normalize before calculating the start value for
		// the newton iteration. Set p1 = mantissa(x). Then
		// take the square root of the leading word and store
		// it in p2. Our approximation is at (bits_per_double - 1)
		// bits correct, so 1 + log2(digit_precision) steps
		// are sufficient
		//

		p1.e = -(p1.length() * bits_per_base_digit) + eps;

		double beta = static_cast<double>(p1.m.most_significant_digit());
		beta /= max_base_digit_2;
		beta = std::sqrt((eps + 1) * beta);

		bigfloat::set_mode(MP_EXACT);
		y.assign(beta);

#ifdef BIGFLOAT_SQRT_DEBUG
		std::cout << "start approximation y = " << y << std::endl;
		std::cout << "y.exponent() = " << y.exponent() << std::endl;
		std::cout << "y.mantissa() = " << y.mantissa() << std::endl;
#endif

		// set the starting precision
		q = bits_per_base_digit - eps -2;

		// but we always need q-2, so compute it
		q -= 2;

		// Determine the number of iterations
		n = 1 + static_cast<long>(LiDIA::log2(static_cast<double>(1 + (t-2)/q)));

#ifdef BIGFLOAT_SQRT_DEBUG
		std::cout << "q = " << q << std::endl;
		std::cout << "number of iterations n = " << n << std::endl;
#endif

		//
		// check maximal precision:
		// 2^n q + 8 + 2 <= max_long
		//

		tmp1.assign(q);
		shift_left (tmp1, tmp1, n);
		add (tmp1, tmp1, bigint(10));

		if (tmp1 > bigint(max_long)) {
			lidia_error_handler ("bigfloat2.c::bigfloat::sqrt(bigfloat, bigfloat, long)",
					     "Precision overflow.");
			return;
		}

		//
		// start Newton iteration
		//

		for (i = 1; i <= n; i++) {
			// eps = 2^i q + 8
			q <<= 1;
			eps = q + 8;

#ifdef BIGFLOAT_SQRT_DEBUG
			std::cout << "eps in step " << i << " : " << eps << std::endl;
			std::cout << "y.bit_length() " << y.bit_length() << std::endl;
			std::cout << "y vor cut : " << y << std::endl;
			std::cout << "y.exponent() = " << y.exponent() << std::endl;
			std::cout << "y.mantissa() = " << y.mantissa() << std::endl;
#endif

			// truncating to eps+1 yields relative eps approx,
			if (y.bit_length() > eps+1)
				y.cut(eps+1);

#ifdef BIGFLOAT_SQRT_DEBUG
			std::cout << "p1.bit_length() " << p1.bit_length() << std::endl;
			std::cout << "y after cut : " << y << std::endl;
			std::cout << "y.exponent() = " << y.exponent() << std::endl;
			std::cout << "y.mantissa() = " << y.mantissa() << std::endl;
#endif

			// truncating to eps+1 yields relative eps approx,
			if (p1.bit_length() > eps+1) {
				p2.cut (p1, eps+1);
				divide (p3, p2, y, eps);
			}
			else
				divide (p3, p1, y, eps);

#ifdef BIGFLOAT_SQRT_DEBUG
			std::cout << "y after division : " << y << std::endl;
			std::cout << "p3 after division : " << p3 << std::endl;
#endif
			add (y, y, p3, eps);
			y.e--;

#ifdef BIGFLOAT_SQRT_DEBUG
			std::cout << "y after iteration " << i << " : " << y << std::endl;
#endif
		}

		//
		// restore the precision and the exponent and
		// normalize the result
		//

		y.e += ex;

		// y is a relative k+1 approx. to sqrt(x);
		// truncating to k+3 yields a relative k approx.
		if (y.bit_length() > k+3)
			y.cut(k+3);

		bigfloat::set_mode(old_rnd_mode);
	}
}



void
log (bigfloat & y, const bigfloat & xx, long t)
{

	int old_rnd_mode = bigfloat::get_mode();

	int  s = 0;
	long d;
	long f1, k, c, n, m, i;

	bigfloat x;
	bigfloat pr;

	//
	// Error, if x <= 0
	//

	if (xx.is_le_zero()) {
		lidia_error_handler("void log(bigfloat&, const bigfloat&, long)",
				    "argument is <= 0");
		return;
	}


	//
	// if xx == 1 then the result is 0
	//

	if (xx.m.is_one() && xx.e == 0) {
		y.assign_zero();
		return;
	}


	bigfloat::set_mode(MP_EXACT);
	x = xx;

	//
	// We need t >= 1
	//

	if (t < 1)
		t = 1;

#ifdef BIGFLOAT_LOG_DEBUG
	std::cout << "Computing a relative " << t;
	std::cout << " approximation to " << x << std::endl;
	std::cout << "x.exponent() = " << x.exponent() << std::endl;
	std::cout << "x.mantissa() = " << x.mantissa() << std::endl;
#endif

	//
	// if 0 < x < 1 then use ln(x) = - ln(1/x)
	//

	if (x.e + x.prec <= 0) {
		//
		// For 0 < x < 1
		// find minimal d with 0 < x <= 1 - 2^{-d}
		// x-1 = - 0.|m| * 2^{e+prec} <= - (1/2) * 2^{e+prec}
		//

		bigfloat::set_mode (MP_EXACT);
		pr = x;
		dec(pr);
		d = -(pr.e + pr.prec - 1);

		//
		// invert x with machine precision
		// eps = 2^{-(t+d+4)}
		//

		s = 1;
		divide(x, 1, x, t+d+4);

		t += 2;

#ifdef BIGFLOAT_LOG_DEBUG
		std::cout << "x is less than 1." << std::endl;
		std::cout << "Found 0 < x <= 1 - 2^{-d}, d = " << d << std::endl;
		std::cout << "Inverting x with machine precision " << t+d+2 << std::endl;
		std::cout << "Approximation of x = 1/x = " << x << std::endl;
		std::cout << "x.exponent() = " << x.exponent() << std::endl;
		std::cout << "x.mantissa() = " << x.mantissa() << std::endl;
		std::cout << "Computing relative " << t << " approx to new x" << std::endl;
#endif
	}


	//
	// if x >= 2 then use ln(x) = 1/2 ln(sqrt(x))
	//

	if (x.e + x.prec > 1) {
		//
		// k >= ceil(log2(ceil(log2(x))))
		//
		f1 = x.e + x.prec;
		k = 0;
		while (f1 != 0) {
			k++;
			f1 >>= 1;
		}

#ifdef BIGFLOAT_LOG_DEBUG
		std::cout << "x >= 2, doing SQUARE ROOT REDUCTION" << std::endl;
		std::cout << "Found k >= ceil(log2(ceil(log2(x)))) = " << k << std::endl;
#endif

		//
		// square root reduction with machine precision
		// eps 2^{-max{t+k+6,8}}
		//
		k = t+k+6;
		if (k < 8)
			k = 8;

		f1 = 0;
		while (x.e + x.prec > 1) {
			f1++;
			sqrt(x, x, k);

#ifdef BIGFLOAT_LOG_DEBUG
			std::cout << "SQRT: new x = " << x << std::endl;
			std::cout << "x.exponent() = " << x.exponent() << std::endl;
			std::cout << "x.mantissa() = " << x.mantissa() << std::endl;
#endif

		}

		t += f1 + 2;

#ifdef BIGFLOAT_LOG_DEBUG
		std::cout << "Taking " << f1 << " square roots with machine ";
		std::cout << " precision k = max{t+k+6, 8} = " << k << std::endl;
		std::cout << " new x = " << x << std::endl;
		std::cout << "x.exponent() = " << x.exponent() << std::endl;
		std::cout << "x.mantissa() = " << x.mantissa() << std::endl;
		std::cout << "Computing relative " << t << " approx to new x" << std::endl;
#endif
	}
	else
		f1 = 0;


	//
	// For 1 < x < 2, we compute a relative
	// t-approximation to ln(x) now.
	//

	//
	// For 1 < x < 2
	// find minimal d >= 1 with < 1 + 2^{-d} <= x
	// x-1 = 0,m * 2^{e+prec} >= 1/2 * 2^{e+prec}
	//

	bigfloat::set_mode (MP_EXACT);
	pr = x;
	dec(pr);
	d = -(pr.e + pr.prec) + 1;

	//
	// Compute an absolute t+d+1-approximation
	// to obtain the relative t-approximation
	// to ln(x)
	//

	t += d+1;

#ifdef BIGFLOAT_LOG_DEBUG
	std::cout << "1 < x < 2" << std::endl;
	std::cout << "Found 1 + 2^{-d} <= x, d = " << d << std::endl;
	std::cout << "Computing a relative " << t << " approximation to x" << std::endl;
#endif

	//
	// Find c >= 1 with 1 < x < 2^{-c+1}
	//

	c = d;

	//
	// n = ceil((t+1)/2c)
	//

	t += 1;
	n = (t/(2*c)) + 1;

#ifdef BIGFLOAT_LOG_DEBUG
	std::cout << "HORNER SCHEME" << std::endl;
	std::cout << "1 < x < 2^{-c+1}, d = c = " << c << std::endl;
	std::cout << "Finding an absolute " << t << " approximation." << std::endl;
	std::cout << "Taking " << n << " terms of the series." << std::endl;
#endif

	//
	// relative eps0 = 2^{-max{t-c+5,8}} approximation z0
	// to z = ((x-1)/(x+1))^2
	//

	bigfloat h1, h2;
	bigfloat z0;

	long eps0 = (t+5)-c;

	if (eps0 < 8)
		eps0 = 8;

	subtract (h1, x, 1, eps0 + 3);
	add      (h2, x, 1, eps0 + 4);
	divide   (pr, h1, h2, eps0 + 3);
	square   (z0, pr, eps0 + 3);

#ifdef BIGFLOAT_LOG_DEBUG
	std::cout << "Computing a relative " << eps0 << " approximation ";
	std::cout << "to z0 = ((x-1)/(x+1))^2" << std::endl;
	std::cout << "z0 = " << z0 << std::endl;
	std::cout << "z0.exponent() = " << z0.exponent() << std::endl;
	std::cout << "z0.mantissa() = " << z0.mantissa() << std::endl;
#endif

	//
	// Horner Scheme to evaluate sum(i=0,n) z / (2i+1)
	// with precisions epsn,...,eps0
	//

	bigfloat zr;

	long c2m1 = (c*2-1);
	long epsn = eps0 - n * c2m1;
	long epsr;
	long b = 7 + static_cast<long>(std::ceil(LiDIA::log2(3.0*n+2.0)));

#ifdef BIGFLOAT_LOG_DEBUG
	std::cout << "lower bound 7+log2(3n+2) for the machine ";
	std::cout << "precision is " << b << std::endl;
	std::cout << "epsn before comparison with b, epsn = " << epsn << std::endl;
#endif

	if (epsn < b)
		epsr = b;
	else
		epsr = epsn;

#ifdef BIGFLOAT_LOG_DEBUG
	std::cout << "Choosing machine precision " << epsr;
	std::cout << " in step " << n << std::endl;
#endif

	m = 2*n + 1;
	invert (y, m, epsr);

#ifdef BIGFLOAT_LOG_DEBUG
	std::cout << "inverted coefficient " << m << " is y = " << y << std::endl;
	std::cout << "y.exponent() = " << y.exponent() << std::endl;
	std::cout << "y.mantissa() = " << y.mantissa() << std::endl;
#endif

	for (i = n-1, m -= 2; i >= 0; i--, m -= 2) {
		// multiplication of an relative
		// eps(i+1)-approximation to z
		// with machine precision eps(i+1)
		//

		if (z0.bit_length() > (epsr+2)) {
			zr.cut   (z0, epsr+2);
			multiply (y, y, zr, epsr);
		}
		else
			multiply (y, y, z0, epsr);

#ifdef BIGFLOAT_LOG_DEBUG
		std::cout << "Multiplication yields y = " << y << std::endl;
		std::cout << "y.exponent() = " << y.exponent() << std::endl;
		std::cout << "y.mantissa() = " << y.mantissa() << std::endl;
#endif

		// determine the larger machine precision eps(i)

		epsn = epsr + c2m1;

		if (epsn < b)
			epsr = b;
		else
			epsr = epsn;

#ifdef BIGFLOAT_LOG_DEBUG
		std::cout << "Changing to new machine precision " << epsr;
		std::cout << " in step " << i << std::endl;
#endif

		// add a relative eps(i)-approx. of 1/(2i+1)
		// to y with machine precision eps(i)

		invert (h1, m, epsr);

#ifdef BIGFLOAT_LOG_DEBUG
		std::cout << "inverted coefficient " << m << " is y = " << h1 << std::endl;
		std::cout << "y.exponent() = " << y.exponent() << std::endl;
		std::cout << "y.mantissa() = " << y.mantissa() << std::endl;
#endif

		add    (y, y, h1, epsr);

#ifdef BIGFLOAT_LOG_DEBUG
		std::cout << "Addition yields y = " << y << std::endl;
		std::cout << "y.exponent() = " << y.exponent() << std::endl;
		std::cout << "y.mantissa() = " << y.mantissa() << std::endl;
#endif
	}

	//
	// 2pr is a relative eps0-approx. to 2(x-1)/(x+1);
	// multiply y and pr with machine precision eps0
	//

	pr.multiply_by_2 ();
	multiply (y, y, pr, eps0);

#ifdef BIGFLOAT_LOG_DEBUG
	std::cout << "Multiplication with 2pr yields y = " << y << std::endl;
	std::cout << "y.exponent() = " << y.exponent() << std::endl;
	std::cout << "y.mantissa() = " << y.mantissa() << std::endl;
#endif

	//
	// multiply the logarithm in case
	// of inversion or square root
	// reduction
	//

	y.e += f1;

	if (s)
		y.negate ();

	//
	// restore old rounding mode
	//

	bigfloat::set_mode (old_rnd_mode);

	if (y.bit_length() > t+3)
		y.cut(t+3);
}



//
// Returns an absolute k-approximation to ln(x)
//

void
log_absolute (bigfloat & l, const bigfloat & x, long k)
{
	bigfloat y;
	long     blnx; // b(|b(x)|)+1
	int      old_rnd_mode;

	old_rnd_mode = bigfloat::get_mode();
	bigfloat::set_mode(MP_EXACT);

	// guarantee 1/2 <= mantissa(x) < 1
	y = x;
	y.base_digit_normalize();

	// b(x) = exponent, if 1/2 <= mantissa(x) < 1
	blnx = y.e + y.prec;

	if (blnx < 0)
		blnx = -blnx;

	blnx = b_value(blnx);
	++blnx;

	// l relative k+1+blnx-approximation to ln(x)
	log (l, y, k+1+blnx);

	// truncate l after k+1 bits right to the comma;
	// hence l is a absolute k-approximation ln(x)
	l.cut_absolute (k+1);

	bigfloat::set_mode(old_rnd_mode);
}



void
exp (bigfloat & y, const bigfloat & xx, long k)
{
	bigfloat x;
	int  old_rnd_mode = bigfloat::get_mode();
	int  sign;
	long l, t, n, f, i, c;

	//
	// If xx == 0, return 1.
	//

	if (xx.m.is_zero()) {
		y.m.assign_one();
		y.e = 0;
		return;
	}


	bigfloat::set_mode (MP_EXACT);
	x = xx;


	//
	// Determine a relative k approximation
	// to exp(x)
	//


	// If x is negative, determine
	// relative k+4 approximation to exp(-x)
	// and invert that approximation.

	if (x.is_lt_zero()) {
		sign = -1;
		x.negate ();
		l = k + 4;

#ifdef BIGFLOAT_EXP_DEBUG
		std::cout << "Negating x, now x = " << x << std::endl;
		std::cout << "Determing relative " << l << " approximation now." << std::endl;
#endif
	}
	else {
		sign = 1;
		l = k;
	}



	//
	// Determine a relative l
	// approximation to exp(x), x > 0
	//

	// x.prec = b(m);
	// f >= 0  iff  x >= 1/2
	// If x >= 1/2 set x = x / 2^{e+b(m)+1}

	f = x.e + x.prec + 1;

	if (f >= 0) {
		t = l + f + b_value(f+1) + 3;
		x.e = -x.prec - 1;

#ifdef BIGFLOAT_EXP_DEBUG
		std::cout << "Transforming x to 0 < x < 1/2, x = " << x << std::endl;
		std::cout << "Determing relative " << t << " approximation now." << std::endl;
#endif
	}
	else
		t = l;



	//
	// Determine a relative t-approximation
	// to exp(x) for 0 < x < 1/2
	//

	t++;

	// For 0 < x < 1/2
	// find maximal c >= 1 with 2^{-c} > x
	// x = 0,m * 2^{e+prec}

	bigfloat::set_mode (MP_EXACT);
	c = -(x.e + x.prec);

	if (c < 1)
		lidia_error_handler ("bigfloat::exp(bigfloat&, const bigfloat&, long)",
				     "c < 1.");


	//
	// n = max{2, ceil((t+1)/(c+1))
	//

	n = (t+1)/(c+1) + 1;
	if (n < 2)
		n = 2;

#ifdef BIGFLOAT_EXP_DEBUG
	std::cout << "HORNER SCHEME" << std::endl;
	std::cout << "0 < x < 2^{-c}, c = " << c << std::endl;
	std::cout << "Finding an absolute " << t << " approximation." << std::endl;
	std::cout << "Taking " << n << " terms of the series." << std::endl;
#endif

	//
	// Horner Scheme to evaluate sum(i=0,n) x^i / i!
	// with precisions epsn,...,eps0
	//

	bigfloat z;
	long eps0 = t+4+b_value(n+1);
	long epsn = eps0 - n * c;
	long epsr;
	long b = 1 + static_cast<long>(std::ceil(LiDIA::log2(3.0*n+2.0)));

#ifdef BIGFLOAT_EXP_DEBUG
	std::cout << "lower bound 1+log2(3n+2) for the machine ";
	std::cout << "precision is " << b << std::endl;
	std::cout << "epsn before comparison with b, epsn = " << epsn << std::endl;
#endif

	if (epsn < b)
		epsr = b;
	else
		epsr = epsn;

#ifdef BIGFLOAT_EXP_DEBUG
	std::cout << "Choosing machine precision " << epsr;
	std::cout << " in step " << n << std::endl;
#endif


	bigfloat::set_mode (MP_EXACT);
	y.assign_one();


#ifdef BIGFLOAT_EXP_DEBUG
	std::cout << "y = " << y << std::endl;
	std::cout << "y.exponent() = " << y.exponent() << std::endl;
	std::cout << "y.mantissa() = " << y.mantissa() << std::endl;
#endif

	for (i = n-1; i >= 0; i--) {
		// multiplication of an relative
		// epsr - approximation to x/(i+1)
		// with machine precision epsr
		//

		if (x.prec > (epsr+3)) {
			z.cut  (x, epsr+3);
			divide (z, z, bigfloat(i+1), epsr+2);
		}
		else
			divide (z, x, bigfloat(i+1), epsr+2);

		multiply (y, y, z, epsr);


#ifdef BIGFLOAT_EXP_DEBUG
		std::cout << "Multiplication yields y = " << y << std::endl;
		std::cout << "y.exponent() = " << y.exponent() << std::endl;
		std::cout << "y.mantissa() = " << y.mantissa() << std::endl;
#endif

		// determine the larger machine precision eps(i)

		epsn = epsr + c;

		if (epsn < b)
			epsr = b;
		else
			epsr = epsn;

#ifdef BIGFLOAT_EXP_DEBUG
		std::cout << "Changing to new machine precision " << epsr;
		std::cout << " in step " << i << std::endl;
#endif

		// add 1 to y exact (we could use machine precision eps(i))
		inc (y);

#ifdef BIGFLOAT_EXP_DEBUG
		std::cout << "Addition of 1 yields y = " << y << std::endl;
		std::cout << "y.exponent() = " << y.exponent() << std::endl;
		std::cout << "y.mantissa() = " << y.mantissa() << std::endl;
#endif
	}

	if (y.prec > t+3)
		y.cut (y, t+3);


	//
	// Now, y is a relative t approximation.
	// Raise it to power f in case of f > 0
	// and obtain a relative l approximation.
	//

	if (f > 0) {
		t = l + f + b_value(f+1) + 3;

		for (i = 1; i <= f; i++) {
			square (y, y, t-i);
		}

		if (y.prec > l+3)
			y.cut (y, l+3);

#ifdef BIGFLOAT_EXP_DEBUG
		std::cout << "after " << f << " squarings" << std::endl;
		std::cout << "y = " << y << std::endl;
#endif
	}


	//
	// Now, y is a relative l approximation.
	// Invert it in case of negative sign
	// and obtain a relative k approximation.
	//

	if (sign == -1) {
		invert (y, y, k+3);

		if (y.prec > k+3)
			y.cut (y, k+3);

#ifdef BIGFLOAT_EXP_DEBUG
		std::cout << "after inversion, y = " << y << std::endl;
#endif
	}

	bigfloat::set_mode(old_rnd_mode);
}



void
add (bigfloat & c, const bigfloat & a,
     const bigfloat & b,
     long             t)
{
	bigfloat::old_rnd_mode = bigfloat::rounding_mode;
	bigfloat::old_bin_prec = bigfloat::binary_precision;

	bigfloat::rounding_mode = MP_RND;
	bigfloat::binary_precision = t;

	add (c, a, b);

	bigfloat::rounding_mode = bigfloat::old_rnd_mode;
	bigfloat::binary_precision = bigfloat::old_bin_prec;
}



void
subtract (bigfloat & c, const bigfloat & a,
	  const bigfloat & b,
	  long             t)
{
	bigfloat::old_rnd_mode = bigfloat::rounding_mode;
	bigfloat::old_bin_prec = bigfloat::binary_precision;

	bigfloat::rounding_mode = MP_RND;
	bigfloat::binary_precision = t;

	subtract (c, a, b);

	bigfloat::rounding_mode = bigfloat::old_rnd_mode;
	bigfloat::binary_precision = bigfloat::old_bin_prec;
}



void
multiply (bigfloat & c, const bigfloat & a, const bigfloat & b, long t)
{
	bigfloat::old_rnd_mode = bigfloat::rounding_mode;
	bigfloat::old_bin_prec = bigfloat::binary_precision;

	bigfloat::rounding_mode = MP_RND;
	bigfloat::binary_precision = t;

	multiply(c, a, b);

	bigfloat::rounding_mode = bigfloat::old_rnd_mode;
	bigfloat::binary_precision = bigfloat::old_bin_prec;
}



void
square (bigfloat & c, const bigfloat & a, long t)
{
	bigfloat::old_rnd_mode = bigfloat::rounding_mode;
	bigfloat::old_bin_prec = bigfloat::binary_precision;

	bigfloat::rounding_mode = MP_RND;
	bigfloat::binary_precision = t;

	square (c, a);

	bigfloat::rounding_mode = bigfloat::old_rnd_mode;
	bigfloat::binary_precision = bigfloat::old_bin_prec;
}



void
divide (bigfloat & c, const bigfloat & a, const bigfloat & b, long t)
{
	bigfloat::old_rnd_mode = bigfloat::rounding_mode;
	bigfloat::old_bin_prec = bigfloat::binary_precision;

	bigfloat::rounding_mode = MP_RND;
	bigfloat::binary_precision = t+2;	// +2, because rounding gives
						//     a further error
	divide (c, a, b);

	bigfloat::rounding_mode = bigfloat::old_rnd_mode;
	bigfloat::binary_precision = bigfloat::old_bin_prec;
}



void
invert (bigfloat & c, const bigfloat & a, long t)
{
	bigfloat::old_rnd_mode = bigfloat::rounding_mode;
	bigfloat::old_bin_prec = bigfloat::binary_precision;

	bigfloat::rounding_mode = MP_RND;
	bigfloat::binary_precision = t+2;

	invert (c, a);

	bigfloat::rounding_mode = bigfloat::old_rnd_mode;
	bigfloat::binary_precision = bigfloat::old_bin_prec;
}



void
bigfloat::print_binary (std::ostream & out) const
{
	lidia_size_t i;

	for (i = bit_length()-1; i >= 0; i--)
		out << m.bit(i);

	out << " * 2^" << e;
}



// truncates the mantissa, such that they are
// at most k bits right to the comma

void
bigfloat::cut_absolute (long k)
{
	// If there are no bits right to
	// the comma, then do nothing.

	if (e >= 0)
		return;

	// Otherwise, there are -e bits
	// right to the comma; save all
	// bits of the mantissa left to
	// the comma and k bits right to
	// the comma.

	// There are e+prec bits left
	// to the comma.

	if (e + prec >= 0)
		k += e + prec;

	// There are -(e+prec) leading
	// zeros right to the comma.
	// Keep k-(e+prec) bits of the
	// mantissa.

	else
		k -= (e+prec);

	this->cut(k);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
