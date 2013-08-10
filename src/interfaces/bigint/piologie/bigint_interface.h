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
//	$Id: bigint_interface.h,v 2.6 2001/07/03 13:10:17 hamdy Exp $
//
//	Author	: Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


//
// constructors and destructor; we could leave out some of these
//

inline
bigint::bigint ()
	: I ()
{
	// nothing to do
}



inline
bigint::bigint (int i)
	: I (static_cast<SignDigit>(i))
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
	: I (Natural(static_cast<Digit>(ul)))
{
	// nothing to do
}



inline
bigint::bigint (double d)
{
	assign(d);
}



inline
bigint::bigint (const Integer & a)
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
bigint::~bigint ()
{
	// nothing to do
}



//
// accessors
//

inline int
bigint::bit (unsigned long i) const
{
	return I.testbit(i);
}



inline unsigned long
bigint::most_significant_digit () const
{
	return I.highest();
}



inline unsigned long
bigint::least_significant_digit () const
{
	return I.lowest();
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
	I = static_cast<SignDigit>(i);
}



inline void
bigint::assign (long i)
{
	I = i;
}



inline void
bigint::assign (unsigned long ui)
{
	I = Natural(static_cast<Digit>(ui));
}



inline void
bigint::assign (const bigint_rep_t & a)
{
	I = a;
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

	int lg_I = I.length();
	int lg_a = a.I.length();

	if ((lg_a == lg_I) && (::abs(a.I) == ::abs(I))) {
		return 0;
	}
	if (lg_a > lg_I || ((lg_a == lg_I) && (::abs(a.I) >::abs(I)))) {
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

	if (a.I > I) {
		return -1;
	}
	else {
		return (a.I != I);
	}
}



inline int
bigint::sign () const
{
	return ::sign(I);
}



//
// properties
//

inline lidia_size_t
bigint::length () const
{
	return I.length();
}



inline lidia_size_t
bigint::bit_length () const
{
	return log2(I)+1;
}



//
// predicates
//

inline bool
bigint::is_even () const
{
	return I.even();
}



inline bool
bigint::is_odd () const
{
	return I.odd();
}



//
// converters
//

inline bool
bigint::intify (int & i) const
{
	if (log2(I) > (SIZEOF_INT * CHAR_BIT-1)) {
		return true;
	}
	i = static_cast<int>(I.lowest());
	if (::sign(I) < 0) {
		i = -i;
	}
	return false;
}



inline bool
bigint::longify (long & i) const
{
	if (log2(I) > (SIZEOF_LONG * CHAR_BIT-1)) {
		return true;
	}
	i = static_cast<long>(I.lowest());
	if (::sign(I) < 0) {
		i = -i;
	}
	return false;
}



//
// misc. functions
//

inline double
bigint::radix ()
{
	return std::ldexp(1.0, BETA);
}



inline int
bigint::bits_per_digit ()
{
	return BETA;
}



//
// modifiers
//

inline void
bigint::absolute_value ()
{
	if (::sign(I) < 0) {
		I = -I;
	}
}



inline void
bigint::abs ()
{
	if (::sign(I) < 0) {
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
	++I;
}



inline void
bigint::dec ()
{
	--I;
}



inline void
bigint::multiply_by_2 ()
{
	I <<= 1;
}



inline void
bigint::divide_by_2 ()
{
	I >>= 1;
}



inline void
bigint::swap (bigint & a)
{
	swap(I, a.I);
}



//
// type predicates
//

inline bool
bigint::is_char () const
{
	return ((::sign(I) >= 0) ? (I <= SCHAR_MAX) : (I >= SCHAR_MIN));
}



inline bool
bigint::is_uchar () const
{
	return ((::sign(I) >= 0) && (I <= UCHAR_MAX));
}



inline bool
bigint::is_short () const
{
	return ((::sign(I) >= 0) ? (I <= SHRT_MAX) : (I >= SHRT_MIN));
}



inline bool
bigint::is_ushort () const
{
	return ((::sign(I) >= 0) && (I <= USHRT_MAX));
}



inline bool
bigint::is_int () const
{
	return ((::sign(I) >= 0) ? (I <= INT_MAX) : (I >= INT_MIN));
}



inline bool
bigint::is_uint () const
{
	return ((::sign(I) >= 0) && (I <= USHRT_MAX));
}


inline bool
bigint::is_long () const
{
	return ((::sign(I) >= 0) ? (I <= LONG_MAX) : (I >= LONG_MIN));
}



inline bool
bigint::is_ulong () const
{
	return ((::sign(I) >= 0) && (I <= USHRT_MAX));
}



//
// Procedural versions
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
add (bigint & c, const bigint & a, long b)
{
	c.I = a.I + b;
}



inline void
add (bigint & c, const bigint & a, unsigned long b)
{
	c.I = a.I + b;
}



inline void
subtract (bigint & c, const bigint & a, const bigint & b)
{
	c.I = a.I - b.I;
}



inline void
subtract (bigint & c, const bigint & a, long b)
{
	c.I = a.I - b;
}



inline void
subtract (bigint & c, const bigint & a, unsigned long b)
{
	c.I = a.I - b;
}



inline void
multiply (bigint & c, const bigint & a, const bigint & b)
{
	c.I = a.I * b.I;
}



inline void
multiply (bigint & c, const bigint & a, long b)
{
	c.I = a.I * b;
}



inline void
multiply (bigint & c, const bigint & a, unsigned long b)
{
	c.I = a.I * b;
}



inline void
square (bigint & a, const bigint & b)
{
	a.I = b.I * b.I;
}



inline void
divide (bigint & c, const bigint & a, const bigint & b)
{
	c.I = a.I / b.I;
}



inline void
divide (bigint & c, const bigint & a, long b)
{
	c.I = a.I / b;
}



inline void
divide (bigint & c, const bigint & a, unsigned long b)
{
	c.I = a.I / b;
}



inline void
remainder (bigint & c, const bigint & a, const bigint & b)
{
	c.I = a.I % b.I;
}



inline void
remainder (bigint & c, const bigint & a, long b)
{
	c.I = a.I % b;
}



inline void
remainder (bigint & c, const bigint & a, unsigned long b)
{
	c.I = a.I % b;
}



inline void
remainder (long & r, const bigint & a, long b)
{
	r = a.I % b;
}



inline void
remainder (unsigned long & r, const bigint & a, unsigned long b)
{
	r = ::abs(a.I) % static_cast<Digit>(b);
}



inline long
remainder (const bigint & a, long b)
{
	return a.I % b;
}



inline void
div_rem (bigint & q, bigint & r, const bigint & a, const bigint & b)
{
	div(a.I, b.I, q.I, r.I);
}



inline void
div_rem (bigint & q, bigint & r, const bigint & a, long & b)
{
	div(a.I, b, q.I, r.I);
}



inline void
div_rem (bigint & q, bigint & r, const bigint & a, unsigned long b)
{
	div(a.I, b, q.I, r.I);
}



inline void
div_rem (bigint & q, long & r, const bigint & a, long b)
{
	div(a.I, b, q.I, r);
}



inline void
invert (bigint & a, const bigint & b)
{
	if ((b.I == 1) || (b.I == -1)) {
		a.I = b.I;
	}
	else {
		lidia_error_handler("bigint", "invert::inverting of a non-unit.");
	}
}



inline void
sqrt (bigint & a, const bigint & b)
{
	a.I = sqrt(b.I);
}



inline void
shift_left (bigint & c, const bigint & a, long ui)
{
	if (ui < 0) {
		lidia_error_handler("bigint", "shift_left()::index is negative.");
	}
	c.I = a.I << static_cast<unsigned long>(ui);
}



inline void
shift_right (bigint & c, const bigint & a, long ui)
{
	if (ui < 0) {
		lidia_error_handler("bigint", "shift_right()::index is negative.");
	}
	c.I = a.I >> static_cast<unsigned long>(ui);
}



inline void
bitwise_and (bigint & c, const bigint & a, const bigint & b)
{
	c.I = a.I & b.I;
}



inline void
bitwise_or (bigint & c, const bigint & a, const bigint & b)
{
	c.I = a.I | b.I;
}



inline void
bitwise_xor (bigint & c, const bigint & a, const bigint & b)
{
	c.I = a.I ^ b.I;
}



inline void
bitwise_not (bigint & b, const bigint & a)
{
	b.I = ~a.I;
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
	Integer v, g;

	gcd(a.I, b.I, u.I, v, g);
	return bigint(g);
}



inline bigint
xgcd_right (bigint & v, const bigint & a, const bigint & b)
{
	Integer u, g;

	gcd(a.I, b.I, u, v.I, g);
	return bigint(g);
}



//
// functions
//

inline bigint
abs (const bigint & a)
{
	return bigint(abs(a.I));
}



//
// I/O
//

inline int
string_to_bigint (const char *s, bigint & a)
{
	const char *e = a.I.atoI(s, 10);

	return (e-s);
}



inline int
bigint_to_string (const bigint & a, char *s)
{
	s = Itoa(a.I, s, 10);
	return strlen(s);
}
