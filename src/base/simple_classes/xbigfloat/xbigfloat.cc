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
#include	"LiDIA/xbigfloat.h"
#include	"LiDIA/xdouble.h"
#include	"LiDIA/bigint.h"
//#include	"LiDIA/bigrational.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// assigners
//

void
xbigfloat::assign_zero ()
{
	int old_rnd_mode;

	debug_handler("xbigfloat", "assign(int)");
	old_rnd_mode = bigfloat::get_mode();
	bigfloat::set_mode_no_check(MP_EXACT);
	x.assign_zero();
	bigfloat::set_mode_no_check(old_rnd_mode);
}



void
xbigfloat::assign_one ()
{
	int old_rnd_mode;

	debug_handler("xbigfloat", "assign(int)");
	old_rnd_mode = bigfloat::get_mode();
	bigfloat::set_mode_no_check(MP_EXACT);
	x.assign_one();
	bigfloat::set_mode_no_check(old_rnd_mode);
}



void
xbigfloat::assign (int i)
{
	int old_rnd_mode;

	debug_handler("xbigfloat", "assign(int)");
	old_rnd_mode = bigfloat::get_mode();
	bigfloat::set_mode_no_check(MP_EXACT);
	x.assign(i);
	bigfloat::set_mode_no_check(old_rnd_mode);
}



void
xbigfloat::assign (long i)
{
	int old_rnd_mode;

	debug_handler("xbigfloat", "assign(long)");
	old_rnd_mode = bigfloat::get_mode();
	bigfloat::set_mode_no_check(MP_EXACT);
	x.assign(i);
	bigfloat::set_mode_no_check(old_rnd_mode);
}



void
xbigfloat::assign (unsigned long i)
{
	int old_rnd_mode;

	debug_handler("xbigfloat", "assign(unsigned long)");
	old_rnd_mode = bigfloat::get_mode();
	bigfloat::set_mode_no_check(MP_EXACT);
	x.assign(i);
	bigfloat::set_mode_no_check(old_rnd_mode);
}



void
xbigfloat::assign (double d)
{
	int old_rnd_mode;

	debug_handler("xbigfloat", "assign(double)");
	old_rnd_mode = bigfloat::get_mode();
	bigfloat::set_mode_no_check(MP_EXACT);
	x.assign(d);
	bigfloat::set_mode_no_check(old_rnd_mode);
}



void
xbigfloat::assign (const xdouble & xd)
{
	int old_rnd_mode;

	debug_handler("xbigfloat", "assign(xdouble)");
	old_rnd_mode = bigfloat::get_mode();
	bigfloat::set_mode_no_check(MP_EXACT);
	x.assign(xd);
	bigfloat::set_mode_no_check(old_rnd_mode);
}



void
xbigfloat::assign (const bigint & bi)
{
	int old_rnd_mode;

	debug_handler("xbigfloat", "assign(const bigint &)");
	old_rnd_mode = bigfloat::get_mode();
	bigfloat::set_mode_no_check(MP_EXACT);
	x.assign(bi);
	bigfloat::set_mode_no_check(old_rnd_mode);
}



#if 0
void xbigfloat::assign(const bigint & n, const bigint & d)
{
	int old_rnd_mode;

	debug_handler("xbigfloat", "assign(const bigint&, const bigint&)");
	old_rnd_mode = bigfloat::get_mode();
	bigfloat::set_mode_no_check(MP_EXACT);
	x.assign(n, d);
	bigfloat::set_mode_no_check(old_rnd_mode);
}
#endif



void
xbigfloat::assign (const bigfloat & bf)
{
	int old_rnd_mode;

	debug_handler("xbigfloat", "assign(const bigfloat&)");
	old_rnd_mode = bigfloat::get_mode();
	bigfloat::set_mode_no_check(MP_EXACT);
	x.assign(bf);
	bigfloat::set_mode_no_check(old_rnd_mode);
}



void
xbigfloat::assign (const xbigfloat & xbf)
{
	int old_rnd_mode;

	debug_handler("xbigfloat", "assign(const xbigfloat&)");
	old_rnd_mode = bigfloat::get_mode();
	bigfloat::set_mode_no_check(MP_EXACT);
	x.assign(xbf.x);
	bigfloat::set_mode_no_check(old_rnd_mode);
}



//
// comparators
//

bool
xbigfloat::is_one () const
{
	int old_rnd_mode;

	debug_handler("xbigfloat", "is_one()const");
	old_rnd_mode = bigfloat::get_mode();
	bigfloat::set_mode_no_check(MP_EXACT);

	bigfloat tmp(1L);

	subtract(tmp, tmp, x);
	bigfloat::set_mode_no_check(old_rnd_mode);

	if (tmp.mantissa().is_zero())
		return true;
	return false;
}



int
xbigfloat::compare (const xbigfloat & a) const
{
	int t;
	int old_rnd_mode;


	debug_handler("xbigfloat", "xbigfloat::compare (const xbigfloat&) const");
	if (&a == this) {
		t = 0;
	}
	else {
		old_rnd_mode = bigfloat::get_mode();
		bigfloat::set_mode_no_check(MP_EXACT);
		t = x.compare(a.x);
		bigfloat::set_mode_no_check(old_rnd_mode);
	}

	return t;
}



//
// modifiers
//

void
xbigfloat::inc ()
{
	int old_rnd_mode;

	debug_handler("xbigfloat", "inc()");
	old_rnd_mode = bigfloat::get_mode();
	bigfloat::set_mode_no_check(MP_EXACT);
	x.inc();
	bigfloat::set_mode_no_check(old_rnd_mode);
}



void
xbigfloat::dec ()
{
	int old_rnd_mode;

	debug_handler("xbigfloat", "dec()");
	old_rnd_mode = bigfloat::get_mode();
	bigfloat::set_mode_no_check(MP_EXACT);
	x.dec();
	bigfloat::set_mode_no_check(old_rnd_mode);
}



//
// Procedural versions
//

void
add (xbigfloat & c, const xbigfloat & a, const xbigfloat & b)
{
	int old_rnd_mode;

	debug_handler("xbigfloat", "add(xbigfloat&, const xbigfloat&, const xbigfloat&)");
	old_rnd_mode = bigfloat::get_mode();
	bigfloat::set_mode_no_check(MP_EXACT);
	add(c.x, a.x, b.x);
	bigfloat::set_mode_no_check(old_rnd_mode);
}



void
subtract (xbigfloat & c, const xbigfloat & a, const xbigfloat & b)
{
	int old_rnd_mode;

	debug_handler("xbigfloat", "subtract(xbigfloat&, const xbigfloat&, const xbigfloat&)");
	old_rnd_mode = bigfloat::get_mode();
	bigfloat::set_mode_no_check(MP_EXACT);
	subtract(c.x, a.x, b.x);
	bigfloat::set_mode_no_check(old_rnd_mode);
}



void
multiply (xbigfloat & c, const xbigfloat & a, const xbigfloat & b)
{
	int old_rnd_mode;

	debug_handler("xbigfloat", "multiply(xbigfloat&, const xbigfloat&, const xbigfloat&)");
	old_rnd_mode = bigfloat::get_mode();
	bigfloat::set_mode_no_check(MP_EXACT);
	multiply(c.x, a.x, b.x);
	bigfloat::set_mode_no_check(old_rnd_mode);
}



void
square (xbigfloat & c, const xbigfloat & a)
{
	int old_rnd_mode;

	debug_handler("xbigfloat", "square(xbigfloat&, const xbigfloat&, const xbigfloat&)");
	old_rnd_mode = bigfloat::get_mode();
	bigfloat::set_mode_no_check(MP_EXACT);
	square(c.x, a.x);
	bigfloat::set_mode_no_check(old_rnd_mode);
}



void
divide (xbigfloat & c, const xbigfloat & a,
	const xbigfloat & b, long k)
{
	debug_handler("xbigfloat", "divide(xbigfloat&, const xbigfloat&, const xbigfloat&, long)");
	divide(c.x, a.x, b.x, k);
}



void
divide (bigint & q, const xbigfloat & a, const xbigfloat & b)

	//
	// determines q = floor (a/b).
	//
{
	debug_handler ("xbigfloat", "divide(bigint, xbigfloat, xbigfloat");

	long k;

	//
	// a = m * 2^e, b = n * 2^f
	//
	// a/b = m/n * 2^(e-f),
	//

	if (check_overflow(k, a.x.exponent(), -b.x.exponent())) {
		lidia_error_handler ("xbigfloat::divide",
				     "Exponent overflow.");
		return;
	}

	bigint num(a.x.mantissa());
	bigint den(b.x.mantissa());
	bigint r;

	if (k > 0)
		shift_left (num, num, k);
	else
		shift_left (den, den, -k);

	div_rem (q, r, num, den);
}



//
//  procedures
//

void
shift_left (xbigfloat & c, const xbigfloat & a, long i)
{
	int old_rnd_mode;

	debug_handler("xbigfloat", "shift_left(long)");
	old_rnd_mode = bigfloat::get_mode();
	bigfloat::set_mode_no_check(MP_EXACT);
	LiDIA::shift_left(c.x, a.x, i);
	bigfloat::set_mode_no_check(old_rnd_mode);
}



void
shift_right (xbigfloat & c, const xbigfloat & a, long i)
{
	int old_rnd_mode;

	debug_handler("xbigfloat", "shift_right(long)");
	old_rnd_mode = bigfloat::get_mode();
	bigfloat::set_mode_no_check(MP_EXACT);
	LiDIA::shift_right(c.x, a.x, i);
	bigfloat::set_mode_no_check(old_rnd_mode);
}



void
ceil (bigint & c, const xbigfloat & a)
{
	int old_rnd_mode;

	debug_handler("xbigfloat", "ceil(bigint&, const xbigfloat&)");
	old_rnd_mode = bigfloat::get_mode();
	bigfloat::set_mode_no_check(MP_EXACT);
	LiDIA::ceil(c, a.x);
	bigfloat::set_mode_no_check(old_rnd_mode);
}



void
floor (bigint & c, const xbigfloat & a)
{
	int old_rnd_mode;

	debug_handler("xbigfloat", "floor(bigint&, const xbigfloat&)");
	old_rnd_mode = bigfloat::get_mode();
	bigfloat::set_mode_no_check(MP_EXACT);
	LiDIA::floor(c, a.x);
	bigfloat::set_mode_no_check(old_rnd_mode);
}



//
//   input / output
//

//
// operator<<(std::ostream out)
//
// Prints the xbigfloat to stream out.
//
// Format:   (MANTISSA,EXPONENT)
//
// Internally we use the bigfloat format m * 2^e. Hence,
// we must output (m,e+b(m)).
//

std::ostream &
operator << (std::ostream & out, const xbigfloat & a)
{
	debug_handler("xbigfloat", "operator << (std::ostream&, const xbigfloat&)");
	out << "(" << a.get_mantissa() << ", " << a.get_exponent() << ")";
	return out;
}



//
// operator>>(std::istream in)
//
// Reads the xbigfloat from stream in.
//
// Format:   (MANTISSA,EXPONENT)
//
// Internally we use the bigfloat format m * 2^e. Hence,
// we must read MANTISSA and EXPONENT and assign
//
// m = MANTISSA, e = EXPONENT - b(MANTISSA).
//

std::istream &
operator >> (std::istream & in, xbigfloat & a)
{
	debug_handler("xbigfloat", "operator >> (std::istream&, xbigfloat&)");
	bigint m;
	long e;
	char c;
	int old_rnd_mode;

	// read white spaces
	//
	in >> c;
	while (c == ' ') in >> c;

	// read "("
	if (c != '(') {
		lidia_error_handler("xbigfloat::operator >> (std::istream, a)",
				    "(expected.");
		return in;
	}

	// read mantissa
	in >> m;

	// read white spaces
	in >> c;
	while (c == ' ') in >> c;

	// read ','
	if (c != ',') {
		lidia_error_handler ("xbigfloat::operator(std::istream & in, ...)",
				     "mantissa read, ',' expected.");
		return in;
	}

	// read exponent
	in >> e;

	// read white spaces
	//
	in >> c;
	while (c == ' ') in >> c;

	// read ")"
	if (c != ')') {
		lidia_error_handler("xbigfloat::operator >> (std::istream, a)",
				    "mantissa, exponent read, ')' expected.");
		return in;
	}


	// assign
	old_rnd_mode = bigfloat::get_mode();
	bigfloat::set_mode_no_check(MP_EXACT);

	if (m == 0)
		a.assign(m);
	else {
		a.x.assign(m);
		e -= m.bit_length();
		if (e >= 0)
			shift_left(a.x, a.x, e);
		else
			shift_right(a.x, a.x, -e);
	}
	bigfloat::set_mode_no_check(old_rnd_mode);

	// return
	return in;
}



void
xbigfloat::print_as_bigfloat () const
{
	debug_handler("xbigfloat", "print_as_bigfloat()const");
	this->print_as_bigfloat(std::cout);
}



void
xbigfloat::print_as_bigfloat (std::ostream & out) const
{
	debug_handler("xbigfloat", "print_as_bigfloat(std::ostream&)const");
	out << x;
}



//
// Given x,y, k >= 1, and c >= 0, this function returns true,
// iff |x-y| < v, where
//
//     v = 2^{-k+1} * (2^{b(x)} + 2^{b(y)-c}).
//
// This can be used to check the correctness of relative
// approximations, i.e., if
//
//     k >= 1, c >= 0,
//     |x/z -1| < 2^{-k},
//     |y/z -1| < 2^{-k-c},
//
// then |x-y| < v for z != 0.
//

bool
xbigfloat::check_relative_error (const xbigfloat & x,
				 const xbigfloat & y,
				 long k,
				 long c)
{
	xbigfloat u, v;

	v.assign_one();
	if (x.b_value() < 0)
		v >>= -x.b_value();
	else
		v <<= x.b_value();

	u.assign_one();
	if (y.b_value() - c < 0)
		u >>= -(y.b_value()-c);
	else
		u <<= y.b_value()-c;

	v += u;
	v >>= (k-1);

	subtract(u, x, y);
	u.absolute_value();
	return (u < v);
}



//
// Given x, y, k, and c, this function returns true,
// iff |x-y| < v, where
//
//     v = 2^{-k} + 2^{-k-c}.
//
// This can be used to check the correctness of absolute
// approximations, i.e., if
//
// |x - z| < 2^{-k}
// |y - z| < 2^{-k-c}
//
// then |x-y| < v.
//

bool
xbigfloat::check_absolute_error (const xbigfloat & x,
				 const xbigfloat & y,
				 long k,
				 long c)
{
	xbigfloat u, v;

	v.assign_one();
	if (k < 0)
		v >>= -k;
	else
		v <<= k;

	u.assign_one();
	if (-k-c < 0)
		u >>= -(-k-c);
	else
		u <<= -k-c;

	v += u;

	subtract(u, x, y);
	u.absolute_value();
	return (u < v);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
