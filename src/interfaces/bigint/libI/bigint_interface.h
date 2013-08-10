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
//	Author	:Thomas Papanikolaou (TP)
//	Changes	: See CVS log
//
//==============================================================================================


//
// c'tors and d'tor
//

inline
bigint::bigint ()
{
	cI(&I);
}



inline
bigint::bigint (int i)
{
	cIasint(&I, i);
}



inline
bigint::bigint (long l)
{
	cIaslong(&I, l);
}



inline
bigint::bigint (unsigned long ul)
{
	cIasulong(&I, ul);
}



inline
bigint::bigint (double d)
{
	cI(&I);
	Iasdbl(&I, d);
}



inline
bigint::bigint (const bigint_rep_t & a)
{
	cIasI(&I, &a);
}



inline
bigint::bigint (const bigint & a)
{
	cIasI(&I, &a.I);
}



inline
bigint::~bigint()
{
	dI(&I);
}



//
// accessors
//

inline int
bigint::bit (unsigned long i) const
{
	if (i >= static_cast<unsigned long>(Ilog(&I) + 1)) {
		return 0;
	}
	return ((I.vec[i/BitsPerDigit] >> (i % BitsPerDigit)) & 1);
}



inline unsigned long
bigint::most_significant_digit () const
{
	long l = I.length;

	if (l == 0) {
		return 0;
	}

	return I.vec[l - 1];
}



inline unsigned long
bigint::least_significant_digit () const
{
	if (I.length == 0) {
		return 0;
	}

	return I.vec[0];
}



//
// assigners
//

inline void
bigint::assign_zero ()
{
	I.sign = PLUS;
	I.length = 0;
}



inline void
bigint::assign_one ()
{
	I.sign = PLUS;
	I.length = 1;
	*(I.vec) = 1;
}



inline void
bigint::assign (int i)
{
	Iasint(&I, i);
}



inline void
bigint::assign (long ui)
{
	Iaslong(&I, ui);
}



inline void
bigint::assign (unsigned long ui)
{
	Iasulong(&I, ui);
}



inline void
bigint::assign (double d)
{
	Iasdbl(&I, d);
}



inline void
bigint::assign (const bigint_rep_t & a)
{
	if (&a != &I) {
		IasI(&I, &a);
	}
}



inline void
bigint::assign (const bigint & a)
{
	if (&a != this) {
		IasI(&I, &a.I);
	}
}



//
// comparators
//

inline int
bigint::abs_compare (const bigint & a) const
{
	if (&a == this) {
		return 0;
	}

	int lg_I = I.length, lg_a = a.I.length;

	if ((lg_a == lg_I) && DigitVecEq(a.I.vec, I.vec, lg_I)) {
		return 0;
	}
	if ((lg_a > lg_I) || ((lg_a == lg_I) && DigitVecGt(a.I.vec, I.vec, lg_I))) {
		return -1;
	}
	else {
		return 1;
	}
}



inline int
bigint::abs_compare (unsigned long a) const
{
	int lg_I = I.length;

	if ((lg_I == 1) && I.vec[0] == a) {
		return 0;
	}
	if ((lg_I == 0) && a == 0) {
		return 0;
	}
	if ((lg_I == 0) || ((lg_I == 1) && I.vec[0] < a)) {
		return -1;
	}
	else {
		return 1;
	}
}



inline int
bigint::compare (const bigint & a) const
{
	if (&a == this) {
		return 0;
	}

	if (IgtI(&a.I, &I)) {
		return -1;
	}
	else {
		return !IeqI(&a.I, &I);
	}
}



inline int
bigint::compare (const long a) const
{
	bigint _a(a);

	if (IgtI(&_a.I, &I)) {
		return -1;
	}
	else {
		return !IeqI(&_a.I, &I);
	}
}



inline int
bigint::compare (const unsigned long a) const
{
	bigint	_a(a);

	if (IgtI(&_a.I, &I)) {
		return -1;
	}
	else {
		return !IeqI(&_a.I, &I);
	}
}



inline int
bigint::sign () const
{
	return ((I.sign) ? -1 : (I.length != 0));
}



//
// properties
//

inline lidia_size_t
bigint::length () const
{
	return I.length;
}



inline lidia_size_t
bigint::bit_length () const
{
	return Ilog(&I) + 1;
}



//
// predicates
//

inline bool
bigint::is_even () const
{
	return ((I.length == 0) || (!(*(I.vec) & 1)));
}



inline bool
bigint::is_odd () const
{
	return ((I.length != 0) && (*(I.vec) & 1));
}



//
// converters
//

inline bool
bigint::intify (int & i) const
{
	if (Iisint(&I)) {
		i = intasI(&I);
		return false;
	}
	else
		return true;
}



inline bool
bigint::longify (long & i) const
{
	if (Iislong(&I)) {
		i = longasI(&I);
		return false;
	}
	else
		return true;
}



//
// misc. functions
//

inline double
bigint::radix ()
{
	return std::ldexp(1.0, BitsPerDigit);
}



inline int
bigint::bits_per_digit ()
{
	return (BitsPerDigit);
}



//
// modifiers
//

inline void
bigint::absolute_value ()
{
	I.sign = PLUS;
}



inline void
bigint::abs ()
{
	I.sign = PLUS;
}



inline void
bigint::negate ()
{
	if (I.length) {
		I.sign ^= MINUS;
	}
}



inline void
bigint::inc ()
{
	Iinc(&I);
}



inline void
bigint::dec ()
{
	Idec(&I);
}



inline void
bigint::multiply_by_2 ()
{
	IslasD(&I, 1);
}



inline void
bigint::divide_by_2 ()
{
	Isr1(&I);
	if (I.length == 0) {
		I.sign = PLUS;
	}
}



//
// type predicates
//

inline bool
bigint::is_char () const
{
	return ((sign() >= 0) ?
		(compare(static_cast<long>(SCHAR_MAX)) <= 0) :
		(compare(static_cast<long>(SCHAR_MIN)) >= 0));
}



inline bool
bigint::is_uchar () const
{
	return ((sign() >= 0) && (compare(static_cast<unsigned long>(UCHAR_MAX)) <= 0));
}



inline bool
bigint::is_short () const
{
	return ((sign() >= 0) ?
		(compare(static_cast<long>(SHRT_MAX)) <= 0) :
		(compare(static_cast<long>(SHRT_MIN)) >= 0));

}



inline bool
bigint::is_ushort () const
{
	return ((sign() >= 0) && (compare(static_cast<unsigned long>(USHRT_MAX)) <= 0));
}



inline bool
bigint::is_int () const
{
	return Iisint(&I);
}



inline bool
bigint::is_uint () const
{
	return Iisuint(&I);
}



inline bool
bigint::is_long () const
{
	return Iislong(&I);
}



inline bool
bigint::is_ulong () const
{
	return Iisulong(&I);
}



//
// arithmetic procedures
//

inline void
negate (bigint & a, const bigint & b)
{
	IasI(&a.I, &b.I);
	if (a.I.length) {
		a.I.sign ^= MINUS;
	}
}



inline void
add (bigint & c, const bigint & a, const bigint & b)
{
	IasIplI(&c.I, &a.I, &b.I);
}



inline void
add (bigint & c, const bigint & a, long b)
{
	bigint _b(b);

	IasIplI(&c.I, &a.I, &_b.I);
}



inline void
add (bigint & c, const bigint & a, unsigned long b)
{
	bigint _b(b);

	IasIplI(&c.I, &a.I, &_b.I);
}



inline void
subtract (bigint & c, const bigint & a, const bigint & b)
{
	IasImiI(&c.I, &a.I, &b.I);
}



inline void
subtract (bigint & c, const bigint & a, long b)
{
	bigint _b(b);

	IasImiI(&c.I, &a.I, &_b.I);
}



inline void
subtract (bigint & c, const bigint & a, unsigned long b)
{
	bigint _b(b);

	IasImiI(&c.I, &a.I, &_b.I);
}



inline void
multiply (bigint & c, const bigint & a, const bigint & b)
{
	IasImuI(&c.I, &a.I, &b.I);
}



inline void
multiply (bigint & c, const bigint & a, long b)
{
	int s = 0;

	if (b < 0) {
		s = 1;
		b = -b;
	}
	IasImuD(&c.I, &a.I, b);
	if (s && c.I.length) {
		c.I.sign ^= MINUS;
	}
}



inline void
multiply (bigint & c, const bigint & a, unsigned long b)
{
	IasImuD(&c.I, &a.I, b);
}



inline void
square (bigint & a, const bigint & b)
{
	IasImuI(&a.I, &b.I, &b.I);
}



inline void
divide (bigint & c, const bigint & a, const bigint & b)
{
	bigint r;
	int sa = a.I.sign;
	int sb = b.I.sign;

	uIdiv(&c.I, &r.I, &a.I, &b.I);
	if (c.I.length) {
		c.I.sign = sa ^ sb;
	}
}



inline void
divide (bigint & c, const bigint & a, long b)
{
	int s = 0;
	int sa = a.I.sign;

	if (b < 0) {
		s = 1;
		b = -b;
	}
	uIasIdiD(&c.I, &a.I, b);
	if (c.I.length) {
		c.I.sign = sa ^ s;
	}
}



inline void
divide (bigint & c, const bigint & a, unsigned long b)
{
	int sa = a.I.sign;

	uIasIdiD(&c.I, &a.I, b);
	if (c.I.length) {
		c.I.sign = sa;
	}
}



inline void
remainder (bigint & c, const bigint & a, const bigint & b)
{
	bigint q;
	int sa = a.I.sign;

	uIdiv(&q.I, &c.I, &a.I, &b.I);
	if (c.I.length) {
		c.I.sign = sa;
	}
}



inline void
remainder (bigint & c, const bigint & a, long b)
{
	bigint q;
	bigint _b(b);
	int sa = a.I.sign;

	uIdiv(&q.I, &c.I, &a.I, &_b.I);
	if (c.I.length) {
		c.I.sign = sa;
	}
}



inline void
remainder (bigint & c, const bigint & a, unsigned long b)
{
	bigint q;
	bigint _b(b);
	int sa = a.I.sign;

	uIdiv(&q.I, &c.I, &a.I, &_b.I);
	if (c.I.length) {
		c.I.sign = sa;
	}
}



inline void
remainder (long & r, const bigint & a, long b)
{
	bigint q;

	if (b < 0) {
		b = -b;
	}
	r = uIasIdiD(&q.I, &a.I, b);
	if (r != 0) {
		if (a.I.sign)
			r = -r;
	}
}



inline void
remainder (unsigned long & r, const bigint & a, unsigned long b)
{
	bigint q;

	r = uIasIdiD(&q.I, &a.I, b);
	if (r != 0) {
		if (a.I.sign)
			r = -r;
	}
}



inline long
remainder (const bigint & a, long b)
{
	long r;

	remainder(r, a, b);
	return r;
}



inline void
div_rem (bigint & q, bigint & r, const bigint & a, const bigint & b)
{
	int sa = a.I.sign;
	int sb = b.I.sign;

	uIdiv(&q.I, &r.I, &a.I, &b.I);
	if (q.I.length) {
		q.I.sign = sa ^ sb;
	}
	if (r.I.length) {
		r.I.sign = sa;
	}
}



inline void
div_rem (bigint & q, bigint & r, const bigint & a, long b)
{
	bigint _b(b);
	int sa = a.I.sign;
	int sb = _b.I.sign;

	uIdiv(&q.I, &r.I, &a.I, &_b.I);
	if (q.I.length) {
		q.I.sign = sa ^ sb;
	}
	if (r.I.length) {
		r.I.sign = sa;
	}
}



inline void
div_rem (bigint & q, bigint & r, const bigint & a, unsigned long b)
{
	bigint _b(b);
	int sa = a.I.sign;
	int sb = _b.I.sign;

	uIdiv(&q.I, &r.I, &a.I, &_b.I);
	if (q.I.length) {
		q.I.sign = sa ^ sb;
	}
	if (r.I.length) {
		r.I.sign = sa;
	}
}



inline void
div_rem (bigint & q, long & r, const bigint & a, long b)
{
	int s = 0;
	int sa = a.I.sign;

	if (b < 0) {
		s = 1;
		b = -b;
	}
	r = uIasIdiD(&q.I, &a.I, b);
	if (q.I.length)
		q.I.sign = sa ^ s;
	if (r) {
		if (sa)
			r = -r;
	}
}



inline void
invert (bigint & a, const bigint & b)
{
	if ((b.I.length == 1) && (*(b.I.vec) == 1)) {
		IasI(&a.I, &b.I);
	}
	else {
		lidia_error_handler("bigint", "invert::inverting a non-unit.");
	}
}



inline void
sqrt (bigint & a, const bigint & b)
{
	if (b.is_lt_zero()) {
		lidia_error_handler("bigint", "sqrt(bigint&a, const bigint&b):: b < 0");
	}
	else {
		IassqrtI(&a.I, &b.I);
	}
}



inline void
shift_left (bigint & c, const bigint & a, long ui)
{
	if (ui < 0) {
		lidia_error_handler("bigint", "shift_left()::index is negative.");
	}
	IasIslD(&c.I, &a.I, static_cast<unsigned long>(ui));
}



inline void
shift_right (bigint & c, const bigint & a, long ui)
{
	if (ui < 0) {
		lidia_error_handler("bigint", "shift_right()::index is negative.");
	}
	IasIsrD(&c.I, &a.I, static_cast<unsigned long>(ui));
}



inline void
bitwise_and (bigint & c, const bigint & a, const bigint & b)
{
	IasIandI(&c.I, &a.I, &b.I);
}



inline void
bitwise_or (bigint & c, const bigint & a, const bigint & b)
{
	IasIorI(&c.I, &a.I, &b.I);
}



inline void
bitwise_xor (bigint & c, const bigint & a, const bigint & b)
{
	IasIxorI(&c.I, &a.I, &b.I);
}



inline void
bitwise_not (bigint & b, const bigint & a)
{
	IasnotI(&b.I, &a.I);
}



//
// gcd's
//


inline bigint
gcd (const bigint & a, const bigint & b)
{
	bigint c;

	Igcd(&c.I, &a.I, &b.I);
	return c;
}



inline bigint
bgcd (const bigint & a, const bigint & b)
{
	bigint c;

	Ibgcd(&c.I, &a.I, &b.I);
	return c;
}



inline bigint
dgcd (const bigint & a, const bigint & b)
{
	bigint c;

	Idgcd(&c.I, &a.I, &b.I);
	return c;
}



inline bigint
xgcd (bigint & u, bigint & v, const bigint & a, const bigint & b)
{
	bigint c;

	Ixgcd(&c.I, &u.I, &v.I, &a.I, &b.I);
	if (a.is_zero() && b.is_zero()) {
		v.assign_zero();
	}
	return c;
}



inline bigint
xgcd_left (bigint & u, const bigint & a, const bigint & b)
{
	bigint c;

	Idxgcd_left(&c.I, &u.I, &a.I, &b.I);
	return c;
}



inline bigint
xgcd_right (bigint & v, const bigint & a, const bigint & b)
{
	bigint c;

	Idxgcd_left(&c.I, &v.I, &b.I, &a.I);
	return c;
}



inline bigint
lcm (const bigint & a, const bigint & b)
{
	bigint c;
	bigint d;

	Igcd(&d.I, &a.I, &b.I);
	IasIdiI(&c.I, &a.I, &d.I);
	ImuasI(&c.I, &b.I);
	c.I.sign = PLUS;
	return c;
}



//
// I/O
//


inline int
bigint_to_string (const bigint & a, char *s)
{
	return Itoa(&a.I, s);
}
