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
#include	"LiDIA/bigfloat.h"
#include	"LiDIA/bigfloat_config.h"
#include	<cstdlib>
#include	<cstring>
#include	<cctype>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// assigners
//

void
bigfloat::assign (const bigint & y)
{
	m.assign(y);
	e = 0;
	prec = m.bit_length();
	base_digit_normalize();
	if (prec == 0)
		flag() = PlusZero;
	else
		flag() = Exact;
}



void
bigfloat::assign (const bigint & x, const bigint & y)
{
	bigint q, r;

	div_rem(q, r, x, y);
	if (r.is_zero()) {
		m.assign(q);
		e = 0;
		prec = m.bit_length();
		base_digit_normalize();
		if (prec == 0)
			flag() = PlusZero;
		else
			flag() = Exact;
	}
	else {
		bigfloat tmp(y);
		assign(x);
		divide(*this, *this, tmp);
	}
}



void
bigfloat::assign (const bigfloat & y)
{
	if (&y != this) {
		m.assign(y.m);
		e = y.e;
		prec = y.prec;
		flag() = y.flag();
	}
}



void
bigfloat::assign (double d)
{
	base_digit mant;
	int expo;
	int s = 0;
	int ndigits = 0;

	if (d == 0.0) {
		assign_zero();
	}
	else {
		if (d < 0.0) {
			s = 1;
			d = -d;
		}

		d = std::frexp(d, &expo);
		m.assign_zero();
		while (d != 0) {
			d = std::ldexp(d, bits_per_base_digit);
			mant = static_cast<base_digit>(d);
			d -= mant;
			ndigits++;
			shift_left(m, m, bits_per_base_digit);
			add(m, m, bigint(mant));
		}
		if (s)
			m.negate();
		e = expo - (ndigits * bits_per_base_digit);
		prec = m.bit_length();
		base_digit_normalize();
		flag() = NotExact;
	}
}



void
bigfloat::assign (const xdouble & xd)
{
	bigfloat temp(xd.l());
	assign(xd.h());
	add(*this, *this, temp);
//   char s[200];
//   xdouble_to_string(s, xd);
//   string_to_bigfloat(s, *this);
}



void
bigfloat::assign (int i)
{
	m.assign(i);
	e = 0;
	prec = m.bit_length();
	base_digit_normalize();
	if (prec == 0)
		flag() = PlusZero;
	else
		flag() = Exact;
}



void
bigfloat::assign (long i)
{
	m.assign(i);
	e = 0;
	prec = m.bit_length();
	base_digit_normalize();
	if (prec == 0)
		flag() = PlusZero;
	else
		flag() = Exact;
}



void
bigfloat::assign (unsigned long ui)
{
	m.assign(ui);
	e = 0;
	prec = m.bit_length();
	base_digit_normalize();
	if (prec == 0)
		flag() = PlusZero;
	else
		flag() = Exact;
}



void
bigfloat::assign_exact_zero ()
{
	m.assign_zero();
	e = 0;
	prec = 0;
	flag() = PlusZero;
}



void
bigfloat::assign_zero ()
{
	m.assign_zero();
	e = -bigfloat::binary_precision;
	prec = 0;
	flag() = NotExact;
}



void
bigfloat::assign_one ()
{
	m.assign_one();
	e = 0;
	prec = 1;
	base_digit_normalize();
	flag() = Exact;
}



//
// comparators
//

int
bigfloat::abs_compare (const bigfloat & y) const
{
	long sx, sy, ex, ey;

	if (&y == this) {
		return 0;
	}

	sx = m.sign();
	sy = y.m.sign();
	if (sx == 0) {
		return ((sy == 0) ? 0 : -1);
	}
	if (sy == 0) {
		return 1;
	}
	// sx and sy are non zero

	ex = e + m.bit_length();
	ey = y.e + y.m.bit_length();

	if (ex > ey) {
		// |x| > |y|
		return 1;
	}
	if (ex < ey) {
		// |x| < |y|
		return -1;
	}

	// ex = ey, i.e. x and y have the same magnitude
	bigfloat tmp;

	if (sx == sy) {
		subtract(tmp, *this, y);
	}
	else {
		add(tmp, *this, y);
	}
	return ((sx > 0) ? tmp.sign() : -tmp.sign());
}



int
bigfloat::compare (const bigfloat & y) const
{
	long sx, sy, ex, ey;

	if (&y == this) {
		return 0;
	}

	sx = m.sign();
	sy = y.m.sign();
	if (sx > sy) {
		return 1;
	}
	if (sx < sy) {
		return -1;
	}
	// sx = sy, i.e. x and y have the same sign

	if (sx == 0) {
		// since sx = sy = 0, x = y = 0
		return 0;
	}

	ex = e + m.bit_length();
	ey = y.e + y.m.bit_length();
	if (ex > ey) {
		// |x| > |y|
		return ((sx > 0) ? 1 : -1);
	}
	if (ex < ey) {
		// |x| < |y|
		return ((sx < 0) ? 1 : -1);
	}

	// ex = ey, i.e. x and y have the same magnitude
	bigfloat tmp;

	subtract(tmp, *this, y);
	return tmp.sign();
}



//
// predicates
//

bool
bigfloat::is_char () const
{
	if (is_zero())
		return true;

	bigfloat tmp;
	truncate(tmp, *this);
	bool b = (compare(tmp) == 0 && (tmp.precision() <= bits_per_(char)));

	return b;
}



bool
bigfloat::is_uchar () const
{
	if (is_zero())
		return true;

	bigfloat tmp;
	truncate(tmp, *this);
	bool b = (compare(tmp) == 0 && (tmp.precision() <= bits_per_(unsigned char)));

	return b;
}



bool
bigfloat::is_short () const
{
	if (is_zero())
		return true;

	bigfloat tmp;
	truncate(tmp, *this);
	bool b = (compare(tmp) == 0 && (tmp.precision() <= bits_per_(short)));

	return b;
}



bool
bigfloat::is_ushort () const
{

	if (is_zero())
		return true;

	bigfloat tmp;
	truncate(tmp, *this);
	bool b = (compare(tmp) == 0 && (tmp.precision() <= bits_per_(unsigned short)));

	return b;
}



bool
bigfloat::is_int () const
{

	if (is_zero())
		return true;

	bigfloat tmp;
	truncate(tmp, *this);
	bool b = (compare(tmp) == 0 && (tmp.precision() <= bits_per_(int)));

	return b;
}



bool
bigfloat::is_uint () const
{

	if (is_zero())
		return true;

	bigfloat tmp;
	truncate(tmp, *this);
	bool b = (compare(tmp) == 0 && (tmp.precision() <= bits_per_(unsigned int)));

	return b;
}



bool
bigfloat::is_long () const
{

	if (is_zero())
		return true;

	bigfloat tmp;
	truncate(tmp, *this);
	bool b = (compare(tmp) == 0 && (tmp.precision() <= bits_per_(long)));

	return b;
}



bool
bigfloat::is_ulong () const
{

	if (is_zero())
		return true;

	bigfloat tmp;
	truncate(tmp, *this);
	bool b = (compare(tmp) == 0 && (tmp.precision() <= bits_per_(unsigned long)));

	return b;
}



bool
bigfloat::is_double () const
{
	long ex = precision();

	if (ex == 0)
		return true;
	ex += exponent();
	return ((ex <= 0x3ff) && (ex > -1023));
}



//
// converters
//

void
bigfloat::bigintify (bigint & y) const
{
	bigfloat r = round(*this);
	y.assign(r.m);
	if (r.e > 0)
		shift_left(y, y, r.e);
	else
		shift_right(y, y, -r.e);
}



bool
bigfloat::doublify (double &d) const
{
	long l = m.bit_length(), expo = e + l;

	if (m.is_zero() || (expo < -1023)) {
		d = 0.0;
		return false;
	}

	if (expo >= 0x3ff) {
		warning_handler("bigfloat", "impossible assignment double = bigfloat");
		return true;
	}

	if (l > (SIZEOF_DOUBLE * 8)) {
		long ed = l - (SIZEOF_DOUBLE * 8);
		bigint bd = m >> ed; // bd = (m/2^ed)*2^ed
		d = dbl(bd) * std::ldexp(2.0, static_cast<int>(e + ed - 1));
	}
	else
		d = dbl(m) * std::ldexp(2.0, static_cast<int>(e - 1));
	return false;
}



bool
bigfloat::xdoublify (xdouble &xd) const
{
	long l = m.bit_length(), expo = e + l;

	if (m.is_zero() || (expo < -996)) {
		xd = 0.0;
		return false;
	}

	if (expo >= 996) {
		warning_handler("bigfloat", "impossible assignment xdouble = bigfloat"
			);
		return true;
	}

	if (l > (2*SIZEOF_DOUBLE * 8)) {
		long ed = l - (2*SIZEOF_DOUBLE * 8);
		bigint bd = m >> ed; // bd = (m/2^ed)*2^ed
		xd = xdbl(bd) * exp(log(static_cast<xdouble>(2.0)) * static_cast<int>(e + ed));
	}
	else
		xd = xdbl(m) * exp(log(static_cast<xdouble>(2.0)) * static_cast<int>(e));
	return false;
}



bool
bigfloat::longify (long &i) const
{
	double d;
	int flag;
	flag = doublify(d);
	if (flag == 1 || (d > max_long)) {
		warning_handler("bigfloat", "impossible assignment long = bigfloat");
		return true;
	}
	i = static_cast<long>(d);
	return false;
}



bool
bigfloat::intify (int &i) const
{
	double d;
	int flag;

	flag = doublify(d);
	if (flag == 1 || (d > max_int)) {
		warning_handler("bigfloat", "impossible assignment int = bigfloat");
		return true;
	}
	i = static_cast<int>(d);
	return false;
}



//
// modifiers
//

void
bigfloat::invert ()
{
	bigint one(1UL), mant(m);
	long t = bigfloat::binary_precision;

	t += ((mant.length() + 1) * bits_per_base_digit);
	LiDIA::shift_left(one, one, t);
	div_rem(m, one, one, mant);
	e = -t - e;
	prec = m.bit_length();
	flag() = NotExact;
	normalize();
}



void
bigfloat::swap (bigfloat & x)
{
	long t = x.e;
	x.e = e;
	e = t;

	t = x.prec;
	x.prec = prec;
	prec = t;

	bigfloat_flag f;
	f.flag() = x.flag();
	x.flag() = flag();
	flag() = f.flag();

	m.swap(x.m);
}



//
// procedural arithmetic
//

void
add (bigfloat & sum, const bigfloat & x, const bigfloat & y)
{
	const bigfloat *px, *py;
	long log_x = x.e + x.prec, log_y = y.e + y.prec;
	long log_dif = log_x - log_y;
	long ed, e1, e2;

	if (log_dif < 0) {
		px = &y;
		py = &x;
		log_dif = -log_dif;
	}
	else {
		px = &x;
		py = &y;
	}


	e1 = px->e;
	e2 = py->e;

	int exp_error = check_overflow(ed, e1, -e2);
	// ed = |e1 - e2|

	if (exp_error) {
		if (exp_error > 0)
			lidia_error_handler("add", "exponent overflow");
		else
			lidia_error_handler("add", "exponent undeflow.");
	}

	if (bigfloat::rounding_mode != MP_EXACT)
		if (log_dif > (bigfloat::binary_precision + rounding_bits)) {
			sum.assign(*px);
			return;
		}

	if (ed > 0) {
		bigint tmp(px->m);
		shift_left(tmp, tmp, ed);
		add(sum.m, tmp, py->m);
		sum.e = e2;
	}
	else {
		bigint tmp(py->m);
		shift_left(tmp, tmp, -ed);
		add(sum.m, tmp, px->m);
		sum.e = e1;
	}
	sum.prec = sum.m.bit_length();
	sum.flag() = NotExact;
	sum.normalize();
}



void
subtract (bigfloat & dif, const bigfloat & x, const bigfloat & y)
{
	const bigfloat *px, *py;
	long log_x = x.e + x.prec, log_y = y.e + y.prec;
	long log_dif = log_x - log_y;
	long ed, e1, e2, i;

	i = log_dif;
	if (log_dif < 0) {
		px = &y;
		py = &x;
		log_dif = -log_dif;
	}
	else {
		px = &x;
		py = &y;
	}


	e1 = px->e;
	e2 = py->e;

	int exp_error = check_overflow(ed, e1, -e2);
	// ed = |e1 - e2|

	if (exp_error) {
		if (exp_error > 0)
			lidia_error_handler("subtract", "exponent overflow");
		else
			lidia_error_handler("subtract", "exponent undeflow.");
	}

	if (bigfloat::rounding_mode != MP_EXACT)
		if (log_dif > (bigfloat::binary_precision + rounding_bits)) {
			dif.assign(*px);
			if (i < 0)
				dif.m.negate();
			return;
		}

	if (ed > 0) {
		bigint tmp(px->m);
		shift_left(tmp, tmp, ed);
		subtract(dif.m, tmp, py->m);
		dif.e = e2;
	}
	else {
		bigint tmp(py->m);
		shift_left(tmp, tmp, -ed);
		subtract(dif.m, px->m, tmp);
		dif.e = e1;
	}
	if (i < 0)
		dif.m.negate();
	dif.prec = dif.m.bit_length();
	dif.flag() = NotExact;
	dif.normalize();
}



void
multiply (bigfloat & prod, const bigfloat & x, const bigfloat & y)
{
	long ex = x.e, ey = y.e;

	multiply(prod.m, x.m, y.m);
	prod.e = ex + ey;
	int exp_error = check_overflow(prod.e, ex, ey);
	if (exp_error) {
		if (exp_error > 0)
			lidia_error_handler("multiply", "exponent overflow");
		else
			lidia_error_handler("multiply", "exponent undeflow.");
	}
	prod.prec = prod.m.bit_length();
	prod.flag() = NotExact;
	prod.normalize();
}



void
multiply (bigfloat & prod, const bigfloat & x, long i)
{
	multiply(prod.m, x.m, i);
	prod.e = x.e;
	prod.prec = prod.m.bit_length();
	prod.flag() = NotExact;
	prod.normalize();
}



void
multiply (bigfloat & prod, long i, const bigfloat & x)
{
	multiply(prod.m, x.m, i);
	prod.e = x.e;
	prod.prec = prod.m.bit_length();
	prod.flag() = NotExact;
	prod.normalize();
}



void
square (bigfloat & c, const bigfloat & a)
{
	long ea = a.e;
	square(c.m, a.m);
	c.e = 2 * ea;
	int exp_error = check_overflow(c.e, ea, ea);
	if (exp_error) {
		if (exp_error > 0)
			lidia_error_handler("square", "exponent overflow");
		else
			lidia_error_handler("square", "exponent undeflow.");
	}
	c.prec = c.m.bit_length();
	c.flag() = NotExact;
	c.normalize();
}



void
divide (bigfloat & quot, const bigfloat & x, const bigfloat & y)
{
	long ex = x.e, ey = y.e;
	bigint tmp1(1), tmp2;
	long t = bigfloat::binary_precision;
	t += ((y.m.length() + 1) * bits_per_base_digit);
	shift_left(tmp1, tmp1, t);
	div_rem(tmp1, tmp2, tmp1, y.m);
	multiply(quot.m, tmp1, x.m);
	quot.e = -t + ex - ey;
	quot.prec = quot.m.bit_length();
	quot.flag() = NotExact;
	quot.normalize();
}



void
divide (bigfloat & quot, const bigfloat & x, long i)
{
	quot.assign(x);
	quot.flag() = NotExact;
	quot.normalize();
	shift_left(quot.m, quot.m, bits_per_base_digit);
	quot.e -= bits_per_base_digit;
	divide(quot.m, quot.m, i);
	quot.prec = quot.m.bit_length();
	quot.normalize();
}



void
divide (bigfloat & quot, long i, const bigfloat & x)
{
	long t = bigfloat::binary_precision;
	bigint tmp(i);
	t += ((x.m.length() + 1) * bits_per_base_digit);
	shift_left(tmp, tmp, t);
	div_rem(quot.m, tmp, tmp, x.m);
	quot.e = -t - x.e;
	quot.prec = quot.m.bit_length();
	quot.flag() = NotExact;
	quot.normalize();
}



void
divide (bigfloat & quot, long i, long j)
{
	long t = bigfloat::binary_precision + bits_per_base_digit;
	quot.m = i;
	shift_left(quot.m, quot.m, t);
	divide(quot.m, quot.m, j);
	quot.e = -t;
	quot.prec = quot.m.bit_length();
	quot.flag() = NotExact;
	quot.normalize();
}



void
ceil (bigfloat & y, const bigfloat & x)
{
	bigfloat fraction;
	if (&y == &x) {
		bigfloat tmp;
		truncate(tmp, x);
		subtract(fraction, x, tmp);
		y.assign(tmp);
	}
	else {
		truncate(y, x);
		subtract(fraction, x, y);
	}
	if (fraction.is_gt_zero())
		inc(y);
}



void
ceil (bigint & y, const bigfloat & x)
{
	truncate(y, x);

	bigfloat fraction(y);
	subtract(fraction, fraction, x);
	fraction.negate();

	if (fraction.is_gt_zero())
		inc(y);
}



void
floor (bigint & y, const bigfloat & x)
{
	truncate(y, x);

	bigfloat fraction(y);
	subtract(fraction, fraction, x);
	fraction.negate();

	if (fraction.is_lt_zero())
		dec(y);
}



void
floor (bigfloat & y, const bigfloat & x)
{
	bigfloat fraction;
	if (&y == &x) {
		bigfloat tmp;
		truncate(tmp, x);
		subtract(fraction, x, tmp);
		y.assign(tmp);
	}
	else {
		truncate(y, x);
		subtract(fraction, x, y);
	}
	if (fraction.is_lt_zero())
		dec(y);
}



void
round (bigint & y, const bigfloat & x)
{
	if (x.is_zero()) {
		y.assign_zero();
		return;
	}

	bigfloat half_plus_x(0.5);
	add(half_plus_x, half_plus_x, x);
	floor(y, half_plus_x);
}



void
round (bigfloat & y, const bigfloat & x)
{
	if (x.is_zero()) {
		y.assign_zero();
		return;
	}

	bigfloat half_plus_x(0.5);
	add(half_plus_x, half_plus_x, x);
	floor(y, half_plus_x);
}



void
truncate (bigint & y, const bigfloat & x)
{
	if (x.e >= 0)
		shift_left(y, x.m, x.e);
	else
		shift_right(y, x.m, -x.e);
}



void
truncate (bigfloat & y, const bigfloat & x)
{
	if (x.e >= 0)
		shift_left(y.m, x.m, x.e);
	else
		shift_right(y.m, x.m, -x.e);
	y.e = 0;
	y.prec = y.m.bit_length();
	y.normalize();
}



int
bigfloat::leading_zeros () const
{
	int l = m.bit_length() % bits_per_base_digit;

	if (l == 0)
		return l;
	else
		return bits_per_base_digit - l;
}



void
bigfloat::base_digit_normalize ()
{
	if (!m.is_zero()) {
		int i = leading_zeros();
		shift_left(m, m, i);
		e -= i;
		prec += i;
	}
}



void
bigfloat::normalize ()
{
	// exact numbers are not normalized;
	// MP_EXACT mode implies no rounding
	if (is_exact() || bigfloat::rounding_mode == MP_EXACT)
		return;

	long dif = prec - bigfloat::binary_precision;
	bigint tmp;

	// length of the number is greater than current precision
	if (dif > 0) {
		// find and store the last bit before cutting
		int addt;
		dif--;
		shift_right(tmp, m, dif);
		e += dif;
		addt = tmp.is_odd();

		// cut
		tmp.divide_by_2();
		e++;

		// if the last bit is equal to 1
		if (addt) {
			int lx = tmp.bit_length();

			// round
			switch (bigfloat::rounding_mode) {
			case MP_RND:
				if (tmp.is_odd()) {
					if (tmp.is_lt_zero())
						tmp.dec();
					else
						tmp.inc();
				}
				break;
			case MP_RND_UP:
				if (tmp.is_gt_zero())
					tmp.inc();
				break;
			case MP_RND_DOWN:
				if (tmp.is_lt_zero())
					tmp.dec();
				break;
			}
			if (tmp.bit_length() > lx) {
				bigint tmp2;
				shift_right(tmp2, tmp, 1);
				tmp2.swap(tmp);
				e++;
			}
		}
		m.swap(tmp);
	}
	// length of the number is less than current precision; shift left
	else {
		dif = -dif;
		shift_left(m, m, dif);
		e -= dif;
	}
	prec = bigfloat::binary_precision;
}



//
// cuts the precision to l bits
//

void
bigfloat::cut (long l)
{
	if (prec <= l) return;
	long dif = prec - l;
	shift_right(m, m, dif);
	e += dif;
	prec -= dif;
}



void
bigfloat::cut (const bigfloat & x, long l)
{
	flag() = x.flag();
	if (x.prec <= l) {
		m.assign(x.m);
		e = x.e;
		prec = x.prec;
		return;
	}
	long dif = x.prec - l;
	shift_right(m, x.m, dif);
	e = x.e + dif;
	prec = x.prec - dif;
}



void
power(bigfloat & z, const bigfloat & x, long y)
{
	if (x.is_zero()) {
		z.assign_zero();
		return;
	}
	if (y == 0) {
		z.assign_one();
		return;
	}

	// tmp stores the current 2^n-power of basis x
	bigfloat tmp(x);
	z.assign_one();

	// use fast exponentiation for computing the power
	long yy = y;

	// if y is negative use absolute value of y
	if (yy < 0)
		yy *= -1;

	if (yy & 1)
		multiply(z, z, tmp);
	yy >>= 1;

	while (yy > 0) {
		square(tmp, tmp);
		if (yy & 1)
			multiply(z, z, tmp);
		yy >>= 1;
	}
	
	// if y is negative we have to invert the result
	if (y < 0)
		invert(z, z);
}



void
power(bigfloat & z, const bigfloat & x, const bigfloat & y)
{
	if (x.is_zero()) {
		z.assign_zero();
		return;
	}

	if (y.is_zero()) {
		z.assign_one();
		return;
	}
	
	long yy = 0;
	y.longify(yy);

	if (y == yy)
		power(z, x, yy);

	else if (x.sign() > 0) {
		bigfloat tmp;
		
		log(tmp, x);
		multiply(z, y, tmp);
		exp(z, z);
	}

	else // i.e. if (x.sign() < 0) {
		lidia_error_handler("bigfloat",
				    "power(bigfloat & z, const bigfloat & x, const bigfloat & y): x < 0 and y not an integer");
}



void
sqrt (bigfloat & y, const bigfloat & x)
{

	if (x.is_lt_zero())
		lidia_error_handler("bigfloat", "sqrt()::argument must be >= 0");
	if (x.is_zero())
		y.assign_zero();
	else {
		long t = bigfloat::binary_precision, l = bigfloat::digit_precision;
		long eps, ex, i, n, l1, l0, l2, l02;
		bigfloat p1, p2, p3;

		p1.assign(x);

		//
		// calculate the exponent ex of the result. If ex is odd
		// multiply the first approximation by 2
		//

		ex = x.e + x.prec;
		eps = ex % 2;
		ex /= 2;

		if (eps < 0) {
			eps += 2;
			ex--;
		}

		//
		// normalize before calculating the start value for
		// the newton iteration. Set p1 = mantissa(x). Then
		// take the square root of the leading words and store
		// it in p2. Our approximation is at (bits_per_double - 1)
		// bits correct, so 1 + log2(digit_precision) steps are sufficient
		//

		p1.e = -(p1.length() * bits_per_base_digit) + eps;

		double beta = static_cast<double>(p1.m.most_significant_digit());
		beta /= max_base_digit_2;
		beta = std::sqrt((eps + 1) * beta);
		y.assign(beta);

		l -= 2;
		n = 1 + static_cast<long>(std::log(static_cast<double>(l)) * INVLOG2);

		//
		// Use increasing precision at each step. l1 is the number
		// of correct Digits at the current point. l2 ist the precision
		// to be set before the next iteration occurs. The code for
		// precision increasing is adapted from the PARI system.
		//

		l1 = 1;
		l2 = 3;
		for (i = 1; i <= n; i++) {
			// from GP/PARI : start

			l0 = l1 << 1;
			if (l0 <= l) {
				l2 += l1;

				//
				// Use only the first l0 + 2 Digits of p3. Remember
				// that the newton iteration is self correcting
				//

				l02 = l0 + 2;
				p3.cut(l02 * bits_per_base_digit);

				//
				// Use only the first l0 + 2 Digits of y. Remember
				// that the newton iteration is self correcting
				//

				y.cut(l02 * bits_per_base_digit);
				l1 = l0;
			}
			else {
				l2 += -l1 + l + 1;
				l1 = l + 1;
			}
			bigfloat::binary_precision = l2 * bits_per_base_digit;

			// from GP/PARI : end

			//
			// Use only the first l2 Digits of p1. Remember
			// that the newton iteration is self correcting
			//

			if (p1.length() > l2) {
				p2.cut(p1, l2 * bits_per_base_digit);
				divide(p3, p2, y);
			}
			else
				divide(p3, p1, y);
			add(y, y, p3);
			y.e--;
		}

		//
		// restore the precision and the exponent and
		// normalize the result
		//

		bigfloat::binary_precision = t;
		y.e += ex;
		y.normalize();
	}
}



void
bigfloat::randomize (const bigint & n, long f)
{
	random_generator rg;
	long ff;

	assign(LiDIA::randomize(n));

	if (f == 0)
		e = 0;
	else {
		rg >> ff;
		e = ff % f;
		if (ff & 1)
			e = -e;
	}
}



long bigfloat::decimal_precision = static_cast<long>(std::ceil(2 * bits_per_base_digit * L2B10));
long bigfloat::binary_precision = 5 * bits_per_base_digit;
long bigfloat::digit_precision = 5;



void
bigfloat::set_precision (long t)
{
	if (t > 0) {
		bigfloat::decimal_precision = t;
		bigfloat::digit_precision = 3 + static_cast<long>(std::ceil(t / (L2B10 * bits_per_base_digit)));
		bigfloat::binary_precision = bits_per_base_digit * bigfloat::digit_precision;
	}
	else {
		bigfloat::decimal_precision = static_cast<long>(std::ceil(2*bits_per_base_digit * L2B10));
		bigfloat::binary_precision = 5 * bits_per_base_digit;
		bigfloat::digit_precision = 5;
	}
}



int bigfloat::rounding_mode = MP_RND;



void
bigfloat::set_mode (int m)
{
	switch (m) {
	case MP_TRUNC:
		bigfloat::rounding_mode = MP_TRUNC;
		break;
	case MP_RND:
		bigfloat::rounding_mode = MP_RND;
		break;
	case MP_RND_UP:
		bigfloat::rounding_mode = MP_RND_UP;
		break;
	case MP_RND_DOWN:
		bigfloat::rounding_mode = MP_RND_DOWN;
		break;
	case MP_EXACT:
		bigfloat::rounding_mode = MP_EXACT;
		break;
	default:
		lidia_error_handler("bigfloat", "mode()::invalid rounding mode");
		break;
	}
}



int
string_to_bigfloat (char *s, bigfloat & n)
{
	int mi, ei, eii, count;
	int j, l = strlen(s), expo = 0;
	char *mant, *exs, c;

	//
	// allocate temporary space for the mantissa
	// and the short exponent; clear both arrays
	//

	mant = new char[l + 2];
	exs = new char[10];
	for (mi = 0; mi <= l + 1; mi++)
		mant[mi] = '\0';
	for (ei = 0; ei < 10; ei++)
		exs[ei] = '\0';

	//
	// mi    -->counts the digits of the mantissa (including sign)
	// ei    -->counts the digits of the exponent (including sign)
	// eii   -->counts the digits after the decimal point
	// count -->is the counter on the array s
	//

	mi = 0;
	ei = 0;
	eii = 0;
	count = 0;

	// remove leading spaces
	c = s[count];
	while (isspace(c) || iscntrl(c)) {
		count++;
		c = s[count];
	}

	// mantissa mign
	if (c == '-') {
		mant[mi] = c;
		count++;
		mi++;
	}
	else if (c == '+')
		count++;

	// mart x
	c = s[count];
	while (isdigit(c)) {
		mant[mi] = c;
		count++;
		mi++;
		c = s[count];
	}

	// radixpoint
	if (c == '.') {
		count++;
		c = s[count];
	}

	// part y
	while (isdigit(c)) {
		mant[mi] = c;
		count++;
		mi++;
		eii--;
		c = s[count];
	}

	// exponent
	if (c == 'E' || c == 'e') {
		count++;
		c = s[count];

		// sign
		if (c == '-') {
			exs[ei] = c;
			count++;
			ei++;
		}
		else if (c == '+')
			count++;

		// digits
		c = s[count];
		while (isdigit(c)) {
			exs[ei] = c;
			count++;
			ei++;
			c = s[count];
		}
	}

        // take care of leading zeros
        size_t leading_zeros = strspn(mant, "0");
        if(mant[leading_zeros] == '\0') {
          if(leading_zeros > 0) {
            --leading_zeros;
          }
          else {
            // this implies mi == 0
            mant[0] = '0';
            mi = 1;
          }
        }
        mi -= leading_zeros;
            
	// store mantissa and exponent in basis 10
	string_to_bigint(mant + leading_zeros, n.m);
	expo = static_cast<int>(atol(exs) + eii);
	n.e = expo;

	// if exponent >= 0 we have an integer
	if (expo >= 0) {
		// calculate n * 10^expo
		long q = expo, r = q % log_max_pow_10;
		for (j = 0; j < r; j++)
			multiply(n.m, n.m, 10L);
		q -= r;
		while (q > 0) {
			multiply(n.m, n.m, max_pow_10);
			q -= log_max_pow_10;
		}
		n.e = 0;
		n.prec = n.m.bit_length();
		n.base_digit_normalize();
		n.flag() = Exact;
	}
	else {
		// we must convert a real
		bigint tmp(1);
		long h, q, r;

		h = static_cast<long>((strlen(mant) + (-expo)) / L2B10);
		q = -expo;
		r = q % log_max_pow_10;
		for (j = 0; j < r; j++)
			multiply(tmp, tmp, 10L);
		q -= r;
		while (q > 0) {
			multiply(tmp, tmp, max_pow_10);
			q -= log_max_pow_10;
		}
		if (h < bigfloat::binary_precision)
			h = bigfloat::binary_precision;
		h <<= 1;
		shift_left(n.m, n.m, h);
		divide(n.m, n.m, tmp);
		n.e = -h;
		n.prec = n.m.bit_length();
		n.flag() = NotExact;
		n.normalize();
	}
	delete[] mant;
	delete[] exs;
	return l;
}



int
bigfloat_to_string (const bigfloat & n, char *s)
{
	if (n.is_zero()) {
		s[0] = '0';
		s[1] = '\0';
		return 1;
	}

	//
	// Normalize the result, so that we have a uniform output
	//

	bigfloat h(n);
	h.flag() = NotExact;
	h.normalize();

	//
	// get the exponent of h and its sign
	//

	long e1 = h.e, sign_exp = 0, e2;
	if (e1 == 0)
		return static_cast<int>(bigint_to_string(h.m, s));
	if (e1 > 0)
		sign_exp = 1;
	else
		sign_exp = -1;

	// set e1 = abs(e1), e2 = e1 * log10(2)
	// case 1:
	//
	// e1 > 0 <=  =  > sign_exp = 1
	//
	// binary_mantissa(h) * 2^e1 = decimal_mantissa(h) * 10^e2, i.e
	// decimal_mantissa(h) = binary_mantissa(h) * 2^e1 / 10^e2
	//
	// e1< 0 <=  =  > sign_exp = -1
	//
	// binary_mantissa(h) * 2^e1 = decimal_mantissa(h) * 10^e2. Now
	// e1< 0 ==  > binary_mantissa(h) * 2^e1 = binary_mantissa(h) / 2^abs(e1)
	// e2< 0 ==  > decimal_mantissa(h) * 10^e2 = decimal_mantissa(h) / 10^abs(e2)
	//
	// ==  > binary_mantissa(h) / 2^abs(e1) = decimal_mantissa(h) / 10^abs(e2)
	// ==  > decimal_mantissa(h) = binary_mantissa(h) * 10^abs(e2) / 2^abs(e1)
	//
	//

	e1 = sign_exp * e1;
	e2 = static_cast<long>(e1 * L2B10);

	long q = e2, r = q % log_max_pow_10, i;
	if (sign_exp == 1) {
		shift_left(h.m, h.m, e1);
		for (i = 0; i < r; i++)
			divide(h.m, h.m, 10L);
		q -= r;
		while (q > 0) {
			divide(h.m, h.m, max_pow_10);
			q -= log_max_pow_10;
		}
	}
	else {
		for (i = 0; i < r; i++)
			multiply(h.m, h.m, 10L);
		q -= r;
		while (q > 0) {
			multiply(h.m, h.m, max_pow_10);
			q -= log_max_pow_10;
		}
		shift_right(h.m, h.m, e1);
	}
	h.e = sign_exp * e2;

	// get the sign ms and convert the decimal_mantissa(h) to a string
	// Adjust the decimal exponent appropriately
	//
	long ms = 0;
	if (h.m.is_lt_zero()) {
		ms = 1;
		h.m.negate();
	}
	bigint_to_string(h.m, s);
	int sl = strlen(s);
	e1 = sl + sign_exp * e2;

	//
	// Now we are going to cut the output to have at most
	// bigfloat::decimal_precision places. We do this fi
	// and only if the string has a length sl greater than
	// the current decimal precision.
	//

	long cut;
	char ch, E[10];

	i = 0;
	if (sl > bigfloat::decimal_precision) {
		//
		// round by taking the last two places in account. This
		// probably causes a carry (i = 1).
		//
		cut = bigfloat::decimal_precision;
		ch = s[cut + 1];
		if (ch > '5')
			s[cut] += (s[cut] > '5');
		ch = s[cut];
		if (ch >= '5')
			i = 1;
		else
			i = 0;
		s[cut] = '\0';
		if (i) {
			//
			// We have a carry which add stringbolically
			//
			cut--;
			while (cut >= 0 && i) {
				if (s[cut] == '9') {
					s[cut] = '0';
					cut--;
					i = 1;
				}
				else {
					s[cut] += 1;
					i = 0;
				}
			}
			if (cut == -1 && i) {
				int m, l = strlen(s) + 1;
				for (m = l; m > 0; m--)
					s[m] = s[m - 1];
				s[0] = '1';
			}
		}
		//
		// Adjust the exponent since there might have been a carry
		//
		sl = static_cast<int>(strlen(s) + i);
		e1 += (sl > bigfloat::decimal_precision);
	}

	//
	// remove trailing zeroes to make ouput shorter
	//

	i = sl - 1;
	while (s[i] == '0' && i > (e1 - 1)) {
		s[i] = '\0';
		i--;
	}

	//
	// check if the decimal point lies at the end of the right side
	// If true, do not print it.
	//

	if ((i + 1) == e1) {
		s[i + 1] = '\0';
		sl = static_cast<int>(i + 1);
		if (ms) {
			for (i = sl; i >= 0; i--)
				s[i + 1] = s[i];
			s[0] = '-';
			sl++;
		}
		return static_cast<int>(strlen(s));
	}

	//
	// if the decimal point lies in the middle of the number
	// shift the fractional part of the number so that we
	// get a character to put the decimal point in it.
	// Else print 0. and the mantissa of the number. If we have
	// an exponent != 0, print it in schientific notation.
	//
	sl = strlen(s);
	if (e1 > 0 && e1 < sl) {
		for (i = sl; i >= e1; i--) {
			s[i + 1] = s[i];
		}
		s[e1] = '.';
		sl++;
		if (ms) {
			for (i = sl; i >= 0; i--)
				s[i + 1] = s[i];
			s[0] = '-';
			sl++;
		}
		return static_cast<int>(strlen(s));
	}
	else {
		for (i = sl + 3 - (!ms); i >= 3 - (!ms); i--)
			s[i] = s[i - 3 + (!ms)];
		if (ms) {
			s[0] = '-';
			s[1] = '0';
			s[2] = '.';
			sl += 3;
		}
		else {
			s[0] = '0';
			s[1] = '.';
			sl += 2;
		}
		if (e1 != 0) {
			s[sl] = 'e';
			s[++sl] = '\0';
			sprintf(E, "%ld", e1);
			strcat(s, E);
		}
	}
	return static_cast<int>(strlen(s));
}



#ifndef HEADBANGER

int
check_overflow(long &z, long x, long y)
{
	if (x == 0) {
		z = y;
		return 0;
	}
	if (y == 0) {
		z = x;
		return 0;
	}

	unsigned long uz = 0;
	int over_under = 0;
	int condition = (x > 0) * 2 + (y > 0);
	int sign = 0;

	switch (condition) {
	case 3:
		// x > 0 and y > 0
		uz = x + y;
		break;
	case 2:
		// x > 0 and y <= 0
		if (x >= -y)
			uz = x - (-y);
		else {
			uz = (-y) - x;
			sign = 1;
		}
		// uz = |x + y|
		break;
	case 1:
		// x <= 0 and y > 0
		if (-x >= y) {
			uz = (-x) - y;
			sign = 1;
		}
		else
			uz = y - (-x);
		// uz = |x + y|
		break;
	case 0:
		// x <= 0 and y <= 0
		uz = (-x) + (-y);
		sign = 1;
		break;
	}
	z = static_cast<long>(uz);
	over_under = static_cast<int>(uz >> bits_per_base_digit_minus_1);
	// over_under will be 1 if the highest bit of uz is set,
	// i.e. if uz doesn't fit into a signed long.  `sign' determines
	// whether this is an over- or under-flow.
	if (sign) {
		z = -z;
		over_under = -over_under;
	}
	return over_under;
}

#endif	// HEADBANGER



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
