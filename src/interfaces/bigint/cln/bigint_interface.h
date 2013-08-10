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
//	Author	: Thomas Papanikolaou (TP), Bruno Haible (HB)
//	Changes	: See CVS log
//
//==============================================================================================


//
// The bit operations assume a sign/magnitude representation. (The doc
// seems to indicate this, and this is supported by the fact that libI,
// lip, gmp all work on sign/magnitude representation. Thomas Papanikolaou
// says the contrary, but I think he's not right about this detail.) HB
//
#define BITOPS_ASSUME_ABS

//
// c'tors and d'tor
//

inline
bigint::bigint ()
	: I ()
{
	// nothing to do
}



inline
bigint::bigint (int i)
	: I (static_cast<long>(i))
{
	// nothing to do
}



inline
bigint::bigint (long l)
	: I (l)
{
	// nothing to do
}



inline
bigint::bigint (unsigned long ul)
	: I (ul)
{
	// nothing to do
}



inline
bigint::bigint (double d)
	: I (round1(cln::cl_DF(d)))
{
	// nothing to do
}



inline
bigint::bigint (const bigint_rep_t & a)
	: I (a)
{
	// nothing to do
}



inline
bigint::bigint (const bigint & a)
	: I (a.I)
{
	// nothing to do
}



inline
bigint::~bigint()
{
	// nothing to do
}



//
// accessors
//

inline int
bigint::bit (unsigned long idx) const
{
#ifdef BITOPS_ASSUME_ABS
	if (minusp(I))
		return static_cast<int>(logbitp(idx, -I));
	else
#endif
		return static_cast<int>(logbitp(idx, I));
}



inline unsigned long
bigint::most_significant_digit () const
{
	if (zerop(I))
		return 0;
	cln::cl_I aI = cln::abs(I);
	return cl_I_to_UL(ldb(aI, cln::cl_byte(intDsize,
					       ((integer_length(aI) - 1) /
						intDsize)*intDsize)));
}



inline unsigned long
bigint::least_significant_digit () const
{
	cln::cl_I aI = cln::abs(I);
	return cl_I_to_UL(ldb(aI, cln::cl_byte(intDsize, 0)));
}



//
// assigners
//

inline void
bigint::assign_zero ()
{
	I = 0;
}



inline void
bigint::assign_one ()
{
	I = 1;
}



inline void
bigint::assign (int i)
{
	I = static_cast<long>(i);
}



inline void
bigint::assign (long i)
{
	I = i;
}



inline void
bigint::assign (unsigned long ui)
{
	I = ui;
}



inline void
bigint::assign (double d)
{
	I = round1(cln::cl_DF(d));
}



inline void
bigint::assign (const bigint_rep_t & a)
{
	if (&a != &I) {
		I = a;
	}
}



inline void
bigint::assign (const bigint & a)
{
	if (&a != this) {
		I = a.I;
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

	return static_cast<int>(cln::compare(cln::abs(I), cln::abs(a.I)));
}



inline int
bigint::abs_compare (unsigned long a) const
{
	return static_cast<int>(cln::compare(cln::abs(I), cln::cl_I(a)));
}



inline int
bigint::compare (const bigint & a) const
{
	if (&a == this) {
		return 0;
	}
	return static_cast<int>(cln::compare(I, a.I));
}



inline int
bigint::compare (unsigned long a) const
{
	return static_cast<int>(cln::compare(I, cln::cl_I(a)));
}



inline int
bigint::compare (long a) const
{
	return static_cast<int>(cln::compare(I, cln::cl_I(a)));
}



inline int
bigint::sign () const
{
	if (minusp(I)) {
		return -1;
	}
	if (zerop(I)) {
		return 0;
	}
	else {
		return 1;
	}
}



//
// properties
//

inline lidia_size_t
bigint::length () const
{
#ifdef BITOPS_ASSUME_ABS
	if (minusp(I))
		return (integer_length(-I)+intDsize-1)/intDsize;
	else
#endif
		return (integer_length(I)+intDsize-1)/intDsize;
}



inline lidia_size_t
bigint::bit_length () const
{
#ifdef BITOPS_ASSUME_ABS
	if (minusp(I))
		return integer_length(-I);
	else
#endif
		return integer_length(I);
}



//
// predicates
//

inline bool
bigint::is_even () const
{
	return evenp(I);
}



inline bool
bigint::is_odd () const
{
	return oddp(I);
}



//
// converters
//

inline bool
bigint::intify (int & i) const
{
	if (integer_length(I) >= int_bitsize)
		return true;
	i = cln::cl_I_to_L(I);
	return false;
}



inline bool
bigint::longify (long & i) const
{
	if (integer_length(I) >= long_bitsize)
		return true;
	i = cl_I_to_long(I);
	return false;
}



inline double
bigint::dbl () const
{
	return cln::double_approx(I);
}



inline xdouble
bigint::xdbl () const
{
	double d1 = dbl();
	bigint a1;

	a1.assign(d1);
	double d2 = cln::double_approx(I - a1.I);
	return xdouble(d1) + xdouble(d2);
}



//
// misc. functions
//

inline double
bigint::radix ()
{
	return std::ldexp(1.0, intDsize);
}



inline int
bigint::bits_per_digit ()
{
	return intDsize;
}



//
// modifiers
//

inline void
bigint::absolute_value ()
{
	if (minusp(I)) {
		I = -I;
	}
}



inline void
bigint::abs ()
{
	if (minusp(I)) {
		I = -I;
	}
}



inline void
bigint::negate ()
{
	I = -I;
}



inline void
bigint::inc ()
{
	I = plus1(I);
}



inline void
bigint::dec ()
{
	I = minus1(I);
}



inline void
bigint::multiply_by_2 ()
{
	I = ash(I, 1);
}



inline void
bigint::divide_by_2 ()
{
	if (minusp(I)) {
		I = ash(plus1(I), -1);
	}
	else {
		I = ash(I, -1);
	}
}



inline void
bigint::swap (bigint & a)
{
	void* tmp = a.I.pointer;

	a.I.pointer = I.pointer;
	I.pointer = tmp;
}



//
// predicates
//

inline bool
bigint::is_char () const
{
	return (!minusp(I) ? (I <= SCHAR_MAX) : (I >= SCHAR_MIN));
}



inline bool
bigint::is_uchar () const
{
	return (!minusp(I) && (I <= UCHAR_MAX));
}



inline bool
bigint::is_short () const
{
	return (!minusp(I) ? (I <= SHRT_MAX) : (I >= SHRT_MIN));
}



inline bool
bigint::is_ushort () const
{
	return (!minusp(I) && (I <= USHRT_MAX));
}



inline bool
bigint::is_int () const
{
	return (!minusp(I) ? (I <= INT_MAX) : (I >= INT_MIN));
}



inline bool
bigint::is_uint () const
{
	return (!minusp(I) && (I <= UINT_MAX));
}



inline bool
bigint::is_long () const
{
	return (!minusp(I) ? (I <= LONG_MAX) : (I >= LONG_MIN));
}



inline bool
bigint::is_ulong () const
{
	return (!minusp(I) && (I <= ULONG_MAX));
}



//
// arithmetic procedures
//

inline void
negate (bigint & a, const bigint & b)
{
	a.I = -b.I;
}



inline void
add (bigint & c, const bigint & a, const bigint & b)
{
	c.I = a.I + b.I;
}



inline void
add (bigint & c, const bigint & a, unsigned long b)
{
	c.I = a.I + b;
}



inline void
add (bigint & c, const bigint & a, long b)
{
	c.I = a.I + b;
}



inline void
subtract (bigint & c, const bigint & a, const bigint & b)
{
	c.I = a.I - b.I;
}



inline void
subtract (bigint & c, const bigint & a, unsigned long b)
{
	c.I = a.I - b;
}



inline void
subtract (bigint & c, const bigint & a, long b)
{
	c.I = a.I - b;
}



inline void
multiply (bigint & c, const bigint & a, const bigint & b)
{
	c.I = a.I * b.I;
}



inline void
multiply (bigint & c, const bigint & a, unsigned long b)
{
	c.I = a.I * b;
}



inline void
multiply (bigint & c, const bigint & a, long b)
{
	c.I = a.I * b;
}



inline void
square (bigint & a, const bigint & b)
{
	a.I = square(b.I);
}



inline void
divide (bigint & c, const bigint & a, const bigint & b)
{
	c.I = truncate1(a.I, b.I);
}



inline void
divide (bigint & c, const bigint & a, unsigned long b)
{
	c.I = truncate1(a.I, b);
}



inline void
divide (bigint & c, const bigint & a, long b)
{
	c.I = truncate1(a.I, b);
}



inline void
remainder (bigint & c, const bigint & a, const bigint & b)
{
	c.I = rem(a.I, b.I);
}



inline void
remainder (bigint & c, const bigint & a, unsigned long b)
{
	c.I = rem(a.I, b);
}



inline void
remainder (bigint & c, const bigint & a, long b)
{
	c.I = rem(a.I, b);
}



inline void
remainder (long & r, const bigint & a, long b)
{
	r = cl_I_to_long(rem(a.I, b));
}



inline void
remainder (unsigned long & r, const bigint & a, unsigned long b)
{
	cln::cl_I t = rem(a.I, b);
	if (minusp(t))
		t = t + b;
	r = cl_I_to_ulong(t);
}



inline long
remainder (const bigint & a, long b)
{
	return cl_I_to_long(rem(a.I, b));
}



inline void
div_rem (bigint & q, bigint & r, const bigint & a, const bigint & b)
{
	cln::cl_I_div_t q_r = truncate2(a.I, b.I);
	q.I = q_r.quotient;
	r.I = q_r.remainder;
}



inline void
div_rem (bigint & q, bigint & r, const bigint & a, unsigned long b)
{
	cln::cl_I_div_t q_r = truncate2(a.I, b);
	q.I = q_r.quotient;
	r.I = q_r.remainder;
}



inline void
div_rem (bigint & q, bigint & r, const bigint & a, long b)
{
	cln::cl_I_div_t q_r = truncate2(a.I, b);
	q.I = q_r.quotient;
	r.I = q_r.remainder;
}



inline void
div_rem (bigint & q, long & r, const bigint & a, long b)
{
	cln::cl_I_div_t q_r = truncate2(a.I, b);
	q.I = q_r.quotient;
	r = cl_I_to_long(q_r.remainder);
}



inline void
invert (bigint & a, const bigint & b)
{
	if ((b.I == 1) || (b.I == -1))
		a.I = b.I;
	else
		lidia_error_handler("bigint", "invert::inverting of a non-unit.");
}



inline void
sqrt (bigint & a, const bigint & b)
{
	a.I = isqrt(b.I);
}



inline void
shift_left (bigint & c, const bigint & a, long ui)
{
	if (ui < 0)
		lidia_error_handler("bigint", "shift_left()::index is negative.");
	c.I = ash(a.I, cln::cl_I(ui));
}



inline void
shift_right (bigint & c, const bigint & a, long ui)
{
	if (ui < 0)
		lidia_error_handler("bigint", "shift_right()::index is negative.");
	c.I = ash(a.I, -cln::cl_I(ui));
}



inline void
bitwise_and (bigint & c, const bigint & a, const bigint & b)
{
	c.I = logand(a.I, b.I);
}



inline void
bitwise_or (bigint & c, const bigint & a, const bigint & b)
{
	c.I = logior(a.I, b.I);
}



inline void
bitwise_xor (bigint & c, const bigint & a, const bigint & b)
{
	c.I = logxor(a.I, b.I);
}



inline void
bitwise_not (bigint & b, const bigint & a)
{
	b.I = lognot(a.I);
}



//
// gcd's
//

inline bigint
gcd (const bigint & a, const bigint & b)
{
	return bigint(gcd(a.I, b.I));
}



inline bigint
bgcd (const bigint & a, const bigint & b)
{
	return bigint(gcd(a.I, b.I));
}



inline bigint
dgcd (const bigint & a, const bigint & b)
{
	return bigint(gcd(a.I, b.I));
}



inline bigint
xgcd_left (bigint & u, const bigint & a, const bigint & b)
{
	cln::cl_I v;
	return bigint(xgcd(a.I, b.I, &u.I, &v));
}



inline bigint
xgcd_right (bigint & v, const bigint & a, const bigint & b)
{
	cln::cl_I u;
	return bigint(xgcd(a.I, b.I, &u, &v.I));
}



inline bigint
lcm (const bigint & a, const bigint & b)
{
	return bigint(lcm(a.I, b.I));
}



// end of bigint_interface.h
