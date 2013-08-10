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
//	Author	: Thomas Papanikolaou (TP)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_BIGCOMPLEX_H_GUARD_
#define LIDIA_BIGCOMPLEX_H_GUARD_



#ifndef LIDIA_BIGFLOAT_H_GUARD_
# include	"LiDIA/bigfloat.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class bigcomplex
{
	//
	// the C++ type we use to represent a bigcomplex
	//

	bigfloat re;
	bigfloat im;

public:

	//
	// c'tors and d'tor
	//

	bigcomplex();
	bigcomplex(const bigfloat &);
	bigcomplex(const bigfloat &, const bigfloat &);
	bigcomplex(const bigcomplex &);
	~bigcomplex();

#ifndef HEADBANGER

	//
	// accessors
	//

	const bigfloat & real() const;
	const bigfloat & imag() const;



	//
	// assigners
	//

	void assign_zero();
	void assign_one();

	void assign(const bigfloat & x);
	void assign(const bigfloat & x, const bigfloat & y);
	void assign(const bigcomplex & x);

	void assign_real(const bigfloat & x);
	void assign_imag(const bigfloat & y);

	bigcomplex & operator = (const bigfloat & y);
	bigcomplex & operator = (const bigcomplex & y);



	//
	// comparators
	//

	bool is_equal(const bigcomplex & a) const;
	bool is_equal(const bigfloat & a) const;
	//bool is_approx_equal (const bigcomplex & a) const;
	//bool is_approx_equal (const bigfloat & a) const;

	bool is_zero() const;
	bool is_approx_zero() const;
	bool is_one() const;
	//bool is_approx_one() const;

	bool is_real() const;
	bool is_approx_real() const;
	bool is_imaginary() const;
	bool is_approx_imaginary() const;

	void negate();
	void invert();



	//
	// modifiers
	//

	void swap(bigcomplex & x);


	static bigint characteristic();

	//
	// Precision and rounding mode setting
	//

	static void set_precision(long l);
	static void set_mode(int l);



	//
	// Procedural versions
	//

	friend void negate(bigcomplex & x, const bigcomplex & y);

	friend void add(bigcomplex & x, const bigcomplex & y, const bigcomplex & z);
	friend void add(bigcomplex & x, const bigcomplex & y, const bigfloat & z);
	friend void add(bigcomplex & x, const bigfloat & y, const bigcomplex & z);

	friend void subtract(bigcomplex & x, const bigcomplex & y, const bigcomplex & z);
	friend void subtract(bigcomplex & x, const bigcomplex & y, const bigfloat & z);
	friend void subtract(bigcomplex & x, const bigfloat & y, const bigcomplex & z);

	friend void multiply(bigcomplex & x, const bigcomplex & y, const bigcomplex & z);
	friend void multiply(bigcomplex & x, const bigcomplex & y, const bigfloat & z);
	friend void multiply(bigcomplex & x, const bigfloat & y, const bigcomplex & z);
	friend void square(bigcomplex & x, const bigcomplex & z);

	friend void divide(bigcomplex & x, const bigcomplex & y, const bigcomplex & z);
	friend void divide(bigcomplex & x, const bigcomplex & y, const bigfloat & z);
	friend void divide(bigcomplex & x, const bigfloat & y, const bigcomplex & z);

	//
	// Functions
	//

	friend void conj(bigcomplex & y, const bigcomplex & x);
	friend void polar(bigcomplex & y, const bigfloat & r, const bigfloat & t);

	friend void sin(bigcomplex & y, const bigcomplex & x);
	friend void cos(bigcomplex & y, const bigcomplex & x);

	friend void sinh(bigcomplex & y, const bigcomplex & x);
	friend void cosh(bigcomplex & y, const bigcomplex & x);

	friend void sqrt(bigcomplex & y, const bigcomplex & x);

	friend void exp(bigcomplex & y, const bigcomplex & x);
	friend void log(bigcomplex & y, const bigcomplex & x);

	friend void power(bigcomplex & y, const bigcomplex & x, const bigcomplex & p);
	friend void power(bigcomplex & y, const bigcomplex & x, const bigfloat & p);
	friend void invert(bigcomplex & x, const bigcomplex & z);

#endif	// HEADBANGER


	//
	// Input/Output
	//

	friend std::istream & operator >> (std::istream & s, bigcomplex & x);
	friend std::ostream & operator << (std::ostream & s, const bigcomplex & x);

};



//
// Constructors and destructor
//

inline
bigcomplex::bigcomplex ()
{
	// nothing to do
}



inline
bigcomplex::bigcomplex (const bigcomplex & y)
	: re(y.re),
	  im(y.im)
{
	// nothing to do
}



inline
bigcomplex::bigcomplex (const bigfloat & y)
	: re(y)
{
	im.assign_exact_zero();
}



inline
bigcomplex::bigcomplex (const bigfloat & r, const bigfloat & i)
	: re(r),
	  im(i)
{
	// nothing to do
}



inline
bigcomplex::~bigcomplex ()
{
	// nothing to do
}



//
// accessors
//

inline const bigfloat &
bigcomplex::real () const
{
	return re;
}



inline const bigfloat &
bigcomplex::imag () const
{
	return im;
}



inline const bigfloat &
real (const bigcomplex & x)
{
	return x.real();
}



inline const bigfloat &
imag (const bigcomplex & x)
{
	return x.imag();
}



//
// assigners
//

inline void
bigcomplex::assign_zero ()
{
	re.assign_zero();
	im.assign_zero();
}



inline void
bigcomplex::assign_one ()
{
	re.assign_one();
	im.assign_zero();
}



inline void
bigcomplex::assign (const bigfloat & x)
{
	re.assign(x);
	im.assign_zero();
}



inline void
bigcomplex::assign (const bigfloat & x, const bigfloat & y)
{
	re.assign(x);
	im.assign(y);
}



inline void
bigcomplex::assign (const bigcomplex & x)
{
	if (&x != this) {
		re.assign(x.re);
		im.assign(x.im);
	}
}



inline void
bigcomplex::assign_real (const bigfloat & x)
{
	re.assign(x);
}



inline void
bigcomplex::assign_imag (const bigfloat & y)
{
	im.assign(y);
}



inline bigcomplex &
bigcomplex::operator = (const bigfloat & y)
{
	assign(y);
	return *this;
}



inline bigcomplex &
bigcomplex::operator = (const bigcomplex & y)
{
	assign(y);
	return *this;
}



//
// comparators
//

inline bool
bigcomplex::is_equal (const bigcomplex & a) const
{
	return (&a == this || (re.compare(a.re) == 0 && im.compare(a.im) == 0));
}



inline bool
bigcomplex::is_equal (const bigfloat & a) const
{
	return (re.compare(a) == 0 && im.is_zero());
}



inline bool
is_equal (const bigcomplex & a, const bigcomplex & b)
{
	return a.is_equal(b);
}



inline bool
is_equal (const bigcomplex & a, const bigfloat & b)
{
	return a.is_equal(b);
}



inline bool
is_equal (const bigfloat & a, const bigcomplex & b)
{
	return b.is_equal(a);
}



inline bool
bigcomplex::is_zero () const
{
	return (re.is_zero() && im.is_zero());
}



inline bool
bigcomplex::is_approx_zero () const
{
	return (re.is_approx_zero() && im.is_approx_zero());
}



inline bool
bigcomplex::is_one () const
{
	return (re.is_one() && im.is_zero());
}



inline bool
operator == (const bigcomplex & x, const bigcomplex & y)
{
	return x.is_equal(y);
}



inline bool
operator == (const bigcomplex & x, const bigfloat & y)
{
	return x.is_equal(y);
}



inline bool
operator == (const bigfloat & x, const bigcomplex & y)
{
	return y.is_equal(x);
}



inline bool
operator != (const bigcomplex & x, const bigcomplex & y)
{
	return !x.is_equal(y);
}



inline bool
operator != (const bigcomplex & x, const bigfloat & y)
{
	return !x.is_equal(y);
}



inline bool
operator != (const bigfloat & x, const bigcomplex & y)
{
	return !y.is_equal(x);
}



//
// predicates
//

inline bool
bigcomplex::is_real () const
{
	return (im.is_zero());
}



inline bool
bigcomplex::is_approx_real () const
{
	return (im.is_approx_zero());
}



inline bool
bigcomplex::is_imaginary () const
{
	return (re.is_zero());
}


inline bool
bigcomplex::is_approx_imaginary () const
{
	return (re.is_approx_zero());
}



inline bool
is_bigfloat (const bigcomplex & a)
{
	return (a.imag().is_zero());
}



//
// modifiers
//

inline void
bigcomplex::negate ()
{
	re.negate();
	im.negate();
}



inline void
bigcomplex::swap (bigcomplex & x)
{
	re.swap(x.re);
	im.swap(x.im);
}



inline void
swap (bigcomplex & x, bigcomplex & y)
{
	x.swap(y);
}



//
// class modifiers
//

inline void
bigcomplex::set_precision (long l)
{
	bigfloat::set_precision(l);
}



inline void
bigcomplex::set_mode (int l)
{
	bigfloat::set_mode(l);
}



//
// procedural arithmetic
//

// (cr, ci) = (xr, xi) + (yr, yi) = (xr+yr, xi+yi)
inline void
add (bigcomplex & c, const bigcomplex & x, const bigcomplex & y)
{
	add(c.re, x.re, y.re);
	add(c.im, x.im, y.im);
}



inline void
add (bigcomplex & c, const bigfloat & x, const bigcomplex & y)
{
	add(c.re, x, y.re);
	c.im.assign(y.im);
}



// (cr, ci) = (xr, xi) + y = (xr+y, xi)
inline void
add (bigcomplex & c, const bigcomplex & x, const bigfloat & y)
{
	add(c.re, x.re, y);
	c.im.assign(x.im);
}



// (cr, ci) = (xr, xi) - (yr, yi) = (xr-yr, xi-yi)
inline void
subtract (bigcomplex & c, const bigcomplex & x, const bigcomplex & y)
{
	subtract(c.re, x.re, y.re);
	subtract(c.im, x.im, y.im);
}



// (cr, ci) = (xr, xi) - y = (xr-y, xi)
inline void
subtract (bigcomplex & c, const bigcomplex & x, const bigfloat & y)
{
	subtract(c.re, x.re, y);
	c.im.assign(x.im);
}



inline void
subtract (bigcomplex & c, const bigfloat & x, const bigcomplex & y)
{
	subtract(c.re, y.re, x);

	c.re.negate();
	c.im.assign(y.im);
	c.im.negate();
}



void multiply(bigcomplex & x, const bigcomplex & y, const bigcomplex & z);



// (cr, ci) = (xr, xi) * y = (xr*y, xi*y)
inline void
multiply (bigcomplex & c, const bigcomplex & x, const bigfloat & y)
{
	multiply(c.re, x.re, y);
	multiply(c.im, x.im, y);
}



inline void
multiply (bigcomplex & c, const bigfloat & x, const bigcomplex & y)
{
	multiply(c.re, y.re, x);
	multiply(c.im, y.im, x);
}



void square(bigcomplex & x, const bigcomplex & z);



void divide(bigcomplex & x, const bigcomplex & y, const bigcomplex & z);



// (cr, ci) = (xr, xi) / y = (xr/y, xi/y)
inline void
divide (bigcomplex & c, const bigcomplex & x, const bigfloat & y)
{
	if (y.is_zero()) {
		lidia_error_handler("bigcomplex", "operator/3::division by zero.");
		return;
	}
	divide(c.re, x.re, y);
	divide(c.im, x.im, y);
}



//
// arithmetic operators
//

inline bigcomplex
operator - (const bigcomplex & x)
{
	bigcomplex c(x);

	c.negate();
	return c;
}



inline bigcomplex
operator + (const bigcomplex & x, const bigcomplex & y)
{
	bigcomplex c;

	add(c, x, y);
	return c;
}



inline bigcomplex
operator + (const bigcomplex & x, const bigfloat & y)
{
	bigcomplex c;

	add(c, x, y);
	return c;
}



inline bigcomplex
operator + (const bigfloat & x, const bigcomplex & y)
{
	bigcomplex c;

	add(c, x, y);
	return c;
}



inline bigcomplex
operator - (const bigcomplex & x, const bigcomplex & y)
{
	bigcomplex c;

	subtract(c, x, y);
	return c;
}



inline bigcomplex
operator - (const bigcomplex & x, const bigfloat & y)
{
	bigcomplex c;

	subtract(c, x, y);
	return c;
}



inline bigcomplex
operator - (const bigfloat & x, const bigcomplex & y)
{
	bigcomplex c;

	subtract(c, x, y);
	return c;
}



inline bigcomplex
operator * (const bigcomplex & x, const bigcomplex & y)
{
	bigcomplex c;

	multiply(c, x, y);
	return c;
}



inline bigcomplex
operator * (const bigcomplex & x, const bigfloat & y)
{
	bigcomplex c;

	multiply(c, x, y);
	return c;
}



inline bigcomplex
operator * (const bigfloat & x, const bigcomplex & y)
{
	bigcomplex c;

	multiply(c, x, y);
	return c;
}



inline bigcomplex
operator / (const bigcomplex & x, const bigcomplex & y)
{
	bigcomplex c;

	divide(c, x, y);
	return c;
}



inline bigcomplex
operator / (const bigcomplex & x, const bigfloat & y)
{
	bigcomplex c;

	divide(c, x, y);
	return c;
}



inline bigcomplex
operator / (const bigfloat & x, const bigcomplex & y)
{
	bigcomplex c;

	divide(c, x, y);
	return c;
}



inline bigcomplex &
operator += (bigcomplex & x, const bigcomplex & y)
{
	add(x, x, y);
	return x;
}



inline bigcomplex &
operator += (bigcomplex & x, const bigfloat & y)
{
	add(x, x, y);
	return x;
}



inline bigcomplex &
operator -= (bigcomplex & x, const bigcomplex & y)
{
	subtract(x, x, y);
	return x;
}



inline bigcomplex &
operator -= (bigcomplex & x, const bigfloat & y)
{
	subtract(x, x, y);
	return x;
}



inline bigcomplex &
operator *= (bigcomplex & x, const bigcomplex & y)
{
	multiply(x, x, y);
	return x;
}



inline bigcomplex &
operator *= (bigcomplex & x, const bigfloat & y)
{
	multiply(x, x, y);
	return x;
}



inline bigcomplex &
operator /= (bigcomplex & x, const bigcomplex & y)
{
	divide(x, x, y);
	return x;
}



inline bigcomplex &
operator /= (bigcomplex & x, const bigfloat & y)
{
	divide(x, x, y);
	return x;
}



//
// misc. functions and procedures
//

// (cr, ci) = (-xr, -xi)
inline void
negate (bigcomplex & c, const bigcomplex & x)
{
	LiDIA::negate(c.re, x.re);
	LiDIA::negate(c.im, x.im);
}



inline void
invert (bigcomplex & c, const bigcomplex & y)
{
	c.assign(y);
	c.invert();
}



inline bigfloat
arg (const bigcomplex & x)
{
	bigfloat c;

	atan2(c, x.imag(), x.real());
	return c;
}



// (cr, ci) = ~(xr, xi) = (xr, -xi)
inline void
conj (bigcomplex & c, const bigcomplex & x)
{
	c.assign(x);
	c.im.negate();
}



inline bigcomplex
conj (const bigcomplex & x)
{
	bigcomplex c;

	conj(c, x);
	return c;
}



// (cr, ci) = (r * cos(t), r * sin(t))
inline void
polar (bigcomplex & c, const bigfloat & r, const bigfloat & t)
{
	bigfloat tmp;

	cos(tmp, t);
	multiply(c.re, r, tmp);
	sin(tmp, t);
	multiply(c.im, r, tmp);
}



inline bigcomplex
polar (const bigfloat & r, const bigfloat & t)
{
	bigcomplex c;

	polar(c, r, t);
	return c;
}



void sqrt(bigcomplex & y, const bigcomplex & x);



inline bigcomplex
sqrt (const bigcomplex & x)
{
	bigcomplex c;

	sqrt(c, x);
	return c;
}



void exp(bigcomplex & y, const bigcomplex & x);
void log(bigcomplex & y, const bigcomplex & x);



inline bigcomplex
exp (const bigcomplex & x)
{
	bigcomplex c;

	exp(c, x);
	return c;
}



inline bigcomplex
log (const bigcomplex & x)
{
	bigcomplex c;

	log(c, x);
	return c;
}



void power(bigcomplex & y, const bigcomplex & x, const bigcomplex & p);
void power(bigcomplex & y, const bigcomplex & x, const bigfloat & p);



inline bigcomplex
power (const bigcomplex & x, const bigcomplex & p)
{
	bigcomplex c;

	power(c, x, p);
	return c;
}



inline bigcomplex
power (const bigcomplex & x, const bigfloat & p)
{
	bigcomplex c;

	power(c, x, p);
	return c;
}



inline bigcomplex
inverse (const bigcomplex & x)
{
	bigcomplex c;

	invert(c, x);
	return c;
}



inline bigcomplex
square (const bigcomplex & x)
{
	bigcomplex c;

	square(c, x);
	return c;
}



inline bigfloat
abs (const bigcomplex & x)
{
	bigfloat c, tmp;

	square(c, x.real());
	square(tmp, x.imag());
	add(c, c, tmp);
	sqrt(c, c);
	return c;
}



inline bigfloat
norm (const bigcomplex & x)
{
	bigfloat c, tmp;

	square(c, x.real());
	square(tmp, x.imag());
	add(c, c, tmp);
	return c;
}



inline bigfloat
hypot (const bigcomplex & x)
{
	bigfloat c, tmp;

	square(c, x.real());
	square(tmp, x.imag());
	add(c, c, tmp);
	sqrt(c, c);
	return c;
}



//
// trigonometric functions
//

void sin(bigcomplex & y, const bigcomplex & x);
void cos(bigcomplex & y, const bigcomplex & x);



inline bigcomplex
sin (const bigcomplex& x)
{
	bigcomplex c;

	sin(c, x);
	return c;
}



inline bigcomplex
cos (const bigcomplex& x)
{
	bigcomplex c;

	cos(c, x);
	return c;
}



//
// hyperbolic functions
//

void sinh(bigcomplex & y, const bigcomplex & x);
void cosh(bigcomplex & y, const bigcomplex & x);



inline bigcomplex
sinh (const bigcomplex& x)
{
	bigcomplex c;

	sinh(c, x);
	return c;
}



inline bigcomplex
cosh (const bigcomplex& x)
{
	bigcomplex c;

	cosh(c, x);
	return c;
}



//
// Input/Output
//

std::istream & operator >> (std::istream & s, bigcomplex & x);
std::ostream & operator << (std::ostream & s, const bigcomplex & x);


//
// some properties of C
//

inline bigint
bigcomplex::characteristic ()
{
	return 0UL;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#define LIDIA_CLASS_BIGCOMPLEX

#include	"LiDIA/specialization/bigcomplex.special"



#endif	// LIDIA_BIGCOMPLEX_H_GUARD_
