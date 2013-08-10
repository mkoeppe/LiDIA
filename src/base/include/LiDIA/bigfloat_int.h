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
//	Author	: Nigel Smart (NiSm), Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_BIGFLOAT_INT_H_GUARD_
#define LIDIA_BIGFLOAT_INT_H_GUARD_



#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_BIGRATIONAL_H_GUARD_
# include	"LiDIA/bigrational.h"
#endif
#ifndef LIDIA_BIGFLOAT_H_GUARD_
# include	"LiDIA/bigfloat.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class bigfloat_int
{
private:

	//
	//  An interval is represented by two bigfloats
	//     (Upp, Low)
	//

	bigfloat Upp, Low;

	//
	// Private functions used for normalisation
	//

	void normalize(bigfloat &, bigfloat &, bigfloat &, bigfloat &);

public:

	//
	// c'tors and d'tor
	//

	bigfloat_int();
	bigfloat_int(int);
	bigfloat_int(long);
	bigfloat_int(unsigned long);
	bigfloat_int(double);
	bigfloat_int(const bigint &);
	bigfloat_int(const bigint &, const bigint &);
	bigfloat_int(const bigrational &);
	bigfloat_int(const bigfloat_int &);
	bigfloat_int(const bigfloat &);
	bigfloat_int(const bigfloat &, const bigfloat &);
	~bigfloat_int();



	//
	// accessors
	//

	const bigfloat & upper() const;
	const bigfloat & lower() const;



	//
	// assigners
	//

	void assign_zero();
	void assign_one();

	void assign(int);
	void assign(long);
	void assign(unsigned long);
	void assign(double);
	void assign(const bigint &);
	void assign(const bigint &, const bigint &);
	void assign(const bigrational &);
	void assign(const bigfloat &);
	void assign(const bigfloat_int &);
	void assign(const bigfloat &, const bigfloat &);

	bigfloat_int & operator = (int);
	bigfloat_int & operator = (long);
	bigfloat_int & operator = (unsigned long);
	bigfloat_int & operator = (double);
	bigfloat_int & operator = (const bigint &);
	bigfloat_int & operator = (const bigrational &);
	bigfloat_int & operator = (const bigfloat &);
	bigfloat_int & operator = (const bigfloat_int &);

	//
	// comparators
	//

	bool contains_zero() const;
	bool is_gt_zero() const;
	bool is_lt_zero() const;

	int sign() const;



	//
	// properties
	//
	void int_length(bigfloat &);



	//
	// converters
	//

	// methods for approximation and conversion
	void dble(double &);
	void Approx();
	void approx(bigfloat &);



	//
	// modifiers
	//

	void absolute_value();

	void inc();
	void dec();

	void negate();
	void invert();

	void swap(bigfloat_int &);



	//
	// arithmetic procedures
	//

	friend void add(bigfloat_int &, const bigfloat_int &, const bigfloat_int &);
	friend void add(bigfloat_int &, const bigfloat_int &, long);
	friend void add(bigfloat_int &, long, const bigfloat_int &);

	friend void subtract(bigfloat_int &, const bigfloat_int &, const bigfloat_int &);
	friend void subtract(bigfloat_int &, const bigfloat_int &, long);
	friend void subtract(bigfloat_int &, long, const bigfloat_int &);

	friend void multiply(bigfloat_int &, const bigfloat_int &, const bigfloat_int &);
	friend void multiply(bigfloat_int &, const bigfloat_int &, long);
	friend void multiply(bigfloat_int &, long, const bigfloat_int &);
	friend void square(bigfloat_int &, const bigfloat_int &);

	friend void divide(bigfloat_int &, const bigfloat_int &, const bigfloat_int &);
	friend void divide(bigfloat_int &, const bigfloat_int &, long);
	friend void divide(bigfloat_int &, long, const bigfloat_int &);

	friend void power(bigfloat_int &, const bigfloat_int &, const bigfloat_int &);
	friend void power(bigfloat_int &, const bigfloat_int &, long);
	friend void sqrt(bigfloat_int &, const bigfloat_int &);



	//
	// Transcedental procedures
	//

	friend void exp(bigfloat_int &, const bigfloat_int &);
	friend void log(bigfloat_int &, const bigfloat_int &);
	friend void besselj(bigfloat_int &, int, const bigfloat_int &);

	//
	// (Inverse) trigonometrical procedures
	//

	friend void sin(bigfloat_int &, const bigfloat_int &);
	friend void cos(bigfloat_int &, const bigfloat_int &);
	friend void tan(bigfloat_int &, const bigfloat_int &);
	friend void cot(bigfloat_int &, const bigfloat_int &);
	friend void asin(bigfloat_int &, const bigfloat_int &);
	friend void acos(bigfloat_int &, const bigfloat_int &);
	friend void atan(bigfloat_int &, const bigfloat_int &);
	friend void atan2(bigfloat_int &, const bigfloat_int &, const bigfloat_int &);
	friend void acot(bigfloat_int &, const bigfloat_int &);

	//
	// (Inverse) hyperbolic trigonometrical procedures
	//

	friend void sinh(bigfloat_int &, const bigfloat_int &);
	friend void cosh(bigfloat_int &, const bigfloat_int &);
	friend void tanh(bigfloat_int &, const bigfloat_int &);
	friend void coth(bigfloat_int &, const bigfloat_int &);
	friend void asinh(bigfloat_int &, const bigfloat_int &);
	friend void acosh(bigfloat_int &, const bigfloat_int &);
	friend void atanh(bigfloat_int &, const bigfloat_int &);
	friend void acoth(bigfloat_int &, const bigfloat_int &);

	//
	// I/O
	//

	friend std::istream & operator >> (std::istream & in, bigfloat_int & a);
	friend std::ostream & operator << (std::ostream & out, const bigfloat_int & a);

};



//
// c'tors and d'tor
//

inline
bigfloat_int::~bigfloat_int ()
{
	// nothing to do
}



inline
bigfloat_int::bigfloat_int ()
	: Upp(0UL),
	  Low(0UL)
{
	// nothing to do
}



inline
bigfloat_int::bigfloat_int (int i)
	: Upp(i),
	  Low(i)
{
	// nothing to do
}



inline
bigfloat_int::bigfloat_int (long l)
	: Upp(l),
	  Low(l)
{
	// nothing to do
}



inline
bigfloat_int::bigfloat_int (unsigned long ul)
	: Upp(ul),
	  Low(ul)
{
	// nothing to do
}



inline
bigfloat_int::bigfloat_int (double d)
	: Upp(d),
	  Low(d)
{
	// nothing to do
}



inline
bigfloat_int::bigfloat_int (const bigint & bi)
	: Upp(bi),
	  Low(bi)
{
	// nothing to do
}



inline
bigfloat_int::bigfloat_int (const bigfloat_int & bf)
	: Upp(bf.Upp),
	  Low(bf.Low)
{
	// nothing to do
}



inline
bigfloat_int::bigfloat_int (const bigfloat & bf1, const bigfloat & bf2)
	: Upp(bf2),
	  Low(bf1)
{
	// nothing to do
}



inline
bigfloat_int::bigfloat_int (const bigfloat & bf)
	: Upp(bf),
	  Low(bf)
{
	// nothing to do
}



inline
bigfloat_int::bigfloat_int (const bigint & n, const bigint & d)
	: Upp(n),
	  Low(n)
{
	divide(*this, *this, d);
}



inline
bigfloat_int::bigfloat_int (const bigrational & r)
	: Upp(r.numerator()),
	  Low(r.numerator())
{
	divide(*this, *this, r.denominator());
}



//
// accessors
//

inline const bigfloat &
bigfloat_int::upper () const
{
	return Upp;
}



inline const bigfloat &
bigfloat_int::lower () const
{
	return Low;
}



//
// assigners
//

inline void
bigfloat_int::assign (const bigint & y)
{
	Upp.assign(y);
	Low.assign(y);
}



inline void
bigfloat_int::assign (const bigint & x, const bigint & y)
{
	assign(x);
	divide(*this, *this, y);
}



inline void
bigfloat_int::assign (const bigrational & r)
{
	assign(r.numerator());
	divide(*this, *this, r.denominator());
}



inline void
bigfloat_int::assign (const bigfloat_int & y)
{
	Upp.assign(y.Upp);
	Low.assign(y.Low);
}



inline void
bigfloat_int::assign (const bigfloat & y)
{
	Upp.assign(y);
	Low.assign(y);
}



inline void
bigfloat_int::assign (const bigfloat & bf1, const bigfloat & bf2)
{
	Upp.assign(bf2);
	Low.assign(bf1);
}



inline void
bigfloat_int::assign (double d)
{
	Upp.assign(d);
	Low.assign(d);
}



inline void
bigfloat_int::assign (int i)
{
	Upp.assign(i);
	Low.assign(i);
}



inline void
bigfloat_int::assign (long i)
{
	Upp.assign(i);
	Low.assign(i);
}



inline void
bigfloat_int::assign (unsigned long ui)
{
	Upp.assign(ui);
	Low.assign(ui);
}



inline void
bigfloat_int::assign_zero ()
{
	Upp.assign_zero();
	Low.assign_zero();
}



inline void
bigfloat_int::assign_one ()
{
	Upp.assign_one();
	Low.assign_one();
}



inline bigfloat_int &
bigfloat_int::operator = (int i)
{
	assign(i);
	return *this;
}



inline bigfloat_int &
bigfloat_int::operator = (long l)
{
	assign(l);
	return *this;
}



inline bigfloat_int &
bigfloat_int::operator = (unsigned long ul)
{
	assign(ul);
	return *this;
}



inline bigfloat_int &
bigfloat_int::operator = (double d)
{
	assign(d);
	return *this;
}



inline bigfloat_int &
bigfloat_int::operator = (const bigint & bi)
{
	assign(bi);
	return *this;
}



inline bigfloat_int &
bigfloat_int::operator = (const bigrational & br)
{
	assign(br);
	return *this;
}



inline bigfloat_int &
bigfloat_int::operator = (const bigfloat & bf)
{
	assign(bf);
	return *this;
}



inline bigfloat_int &
bigfloat_int::operator = (const bigfloat_int & bf)
{
	assign(bf);
	return *this;
}



//
// comparators
//

inline bool
bigfloat_int::contains_zero () const
{
	return ((Upp.is_ge_zero()) && (Low.is_le_zero()));
}



inline bool
bigfloat_int::is_gt_zero () const
{
	return Low.is_gt_zero();
}



inline bool
bigfloat_int::is_lt_zero () const
{
	return Upp.is_lt_zero();
}



// Note: not `==' is not the same as `!=' !
inline bool
operator == (const bigfloat_int & x, const bigfloat_int & y)
{
	return (&x == &y);
}



inline bool
operator != (const bigfloat_int & x, const bigfloat_int & y)
{
	return ((x.lower() > y.upper()) || (x.upper() < y.lower()));
}



inline bool
operator > (const bigfloat_int & x, const bigfloat_int & y)
{
	return (x.lower() > y.upper());
}



inline bool
operator < (const bigfloat_int & x, const bigfloat_int & y)
{
	return (x.upper() < y.lower());
}



//
// modifiers
//

inline void
bigfloat_int::inc ()
{
	add(*this, *this, 1UL);
}



inline void
bigfloat_int::dec ()
{
	subtract(*this, *this, 1UL);
}



inline void
bigfloat_int::swap (bigfloat_int & a)
{
	Upp.swap(a.Upp);
	Low.swap(a.Low);
}



//
// arithmetic procedures
//

void add(bigfloat_int &, const bigfloat_int &, const bigfloat_int &);
void add(bigfloat_int &, const bigfloat_int &, long);
void add(bigfloat_int &, long, const bigfloat_int &);

void subtract(bigfloat_int &, const bigfloat_int &, const bigfloat_int &);
void subtract(bigfloat_int &, const bigfloat_int &, long);
void subtract(bigfloat_int &, long, const bigfloat_int &);

void multiply(bigfloat_int &, const bigfloat_int &, const bigfloat_int &);
void multiply(bigfloat_int &, const bigfloat_int &, long);
void multiply(bigfloat_int &, long, const bigfloat_int &);

void divide(bigfloat_int &, const bigfloat_int &, const bigfloat_int &);
void divide(bigfloat_int &, const bigfloat_int &, long);
void divide(bigfloat_int &, long, const bigfloat_int &);



inline void
inc (bigfloat_int & x)
{
	x.inc();
}



inline void
dec (bigfloat_int & x)
{
	x.dec();
}



inline void
negate (bigfloat_int & x, const bigfloat_int & y)
{
	x.assign(y);
	x.negate();
}



inline void
invert (bigfloat_int & y, const bigfloat_int & x)
{
	y.assign(x);
	y.invert();
}



inline bigfloat_int
inverse (const bigfloat_int & x)
{
	bigfloat_int y;

	y.invert();
	return y;
}



void square(bigfloat_int &, const bigfloat_int &);



inline void
square (bigfloat_int & c, const bigfloat_int & a)
{
	multiply(c, a, a);
}



inline bigfloat_int
square (const bigfloat_int & a)
{
	bigfloat_int c;

	square(c, a);
	return c;
}



//
// arithmetic operators
//

inline bigfloat_int
operator - (const bigfloat_int & a)
{
	bigfloat_int c(a);

	c.negate();
	return c;
}



inline bigfloat_int
operator + (const bigfloat_int & a, const bigfloat_int & b)
{
	bigfloat_int c;

	add(c, a, b);
	return c;
}



inline bigfloat_int
operator - (const bigfloat_int & a, const bigfloat_int & b)
{
	bigfloat_int c;

	subtract(c, a, b);
	return c;
}



inline bigfloat_int
operator * (const bigfloat_int & a, const bigfloat_int & b)
{
	bigfloat_int c;

	multiply(c, a, b);
	return c;
}



inline bigfloat_int
operator / (const bigfloat_int & a, const bigfloat_int & b)
{
	bigfloat_int c;

	divide(c, a, b);
	return c;
}



inline bigfloat_int &
operator += (bigfloat_int & a, const bigfloat_int & b)
{
	add(a, a, b);
	return a;
}



inline bigfloat_int &
operator -= (bigfloat_int & a, const bigfloat_int & b)
{
	subtract(a, a, b);
	return a;
}



inline bigfloat_int &
operator *= (bigfloat_int & a, const bigfloat_int & b)
{
	multiply(a, a, b);
	return a;
}



inline bigfloat_int &
operator /= (bigfloat_int & a, const bigfloat_int & b)
{
	divide(a, a, b);
	return a;
}



inline bigfloat_int &
operator ++ (bigfloat_int & a)
{
	a.inc();
	return a;
}



inline bigfloat_int
operator ++ (bigfloat_int & a, int)
{
	bigfloat_int c(a);

	a.inc();
	return c;
}



inline bigfloat_int &
operator -- (bigfloat_int & a)
{
	a.dec();
	return a;
}



inline bigfloat_int
operator -- (bigfloat_int & a, int)
{
	bigfloat_int c(a);

	a.dec();
	return c;
}



inline void
floor (bigint & y, const bigfloat_int & x)
{
	floor(x.lower()).bigintify(y);
}



inline void
ceil (bigint & y, const bigfloat_int & x)
{
	ceil(x.upper()).bigintify(y);
}



inline int
sign (const bigfloat_int & x)
{
	return x.sign();
}



inline bigfloat_int
abs (const bigfloat_int & x)
{
	bigfloat_int y(x);

	y.absolute_value();
	return y;
}



void truncate(bigint &, const bigfloat_int &);
void round(bigint &, const bigfloat_int &);



inline bigint
truncate (const bigfloat_int & x)
{
	bigint y;

	truncate(y, x);
	return y;
}



inline bigint
round (const bigfloat_int & x)
{
	bigint y;

	round(y, x);
	return y;
}



inline bigint
floor (const bigfloat_int & x)
{
	bigint y;

	floor(y, x);
	return y;
}



inline bigint
ceil (const bigfloat_int & x)
{
	bigint y;

	ceil(y, x);
	return y;
}



void sqrt(bigfloat_int &, const bigfloat_int &);



inline bigfloat_int
sqrt (const bigfloat_int & x)
{
	bigfloat_int y;

	sqrt(y, x);
	return y;
}



void power(bigfloat_int &, const bigfloat_int &, const bigfloat_int &);
void power(bigfloat_int &, const bigfloat_int &, long);



inline bigfloat_int
power (const bigfloat_int & x, const bigfloat_int & y)
{
	bigfloat_int z;

	power(z, x, y);
	return z;
}



inline bigfloat_int
power (const bigfloat_int & x, long n)
{
	bigfloat_int z;

	power(z, x, n);
	return z;
}



//
// Transcedental procedures and functions
//

void exp(bigfloat_int &, const bigfloat_int &);
void log(bigfloat_int &, const bigfloat_int &);
void besselj(bigfloat_int &, int, const bigfloat_int &);



inline bigfloat_int
exp (const bigfloat_int & x)
{
	bigfloat_int y;

	exp(y, x);
	return y;
}



inline bigfloat_int
log (const bigfloat_int & x)
{
	bigfloat_int y;

	log(y, x);
	return y;
}



inline bigfloat_int
besselj (int n, const bigfloat_int & x)
{
	bigfloat_int y;

	besselj(y, n, x);
	return y;
}



//
// (Inverse) trigonometrical procedures and functions
//

void sin(bigfloat_int &, const bigfloat_int &);
void cos(bigfloat_int &, const bigfloat_int &);
void tan(bigfloat_int &, const bigfloat_int &);
void cot(bigfloat_int &, const bigfloat_int &);
void asin(bigfloat_int &, const bigfloat_int &);
void acos(bigfloat_int &, const bigfloat_int &);
void atan(bigfloat_int &, const bigfloat_int &);
void atan2(bigfloat_int &, const bigfloat_int &, const bigfloat_int &);
void acot(bigfloat_int &, const bigfloat_int &);



inline bigfloat_int
sin (const bigfloat_int & x)
{
	bigfloat_int y;

	sin(y, x);
	return y;
}



inline bigfloat_int
cos (const bigfloat_int & x)
{
	bigfloat_int y;

	cos(y, x);
	return y;
}



inline bigfloat_int
tan (const bigfloat_int & x)
{
	bigfloat_int y;

	tan(y, x);
	return y;
}



inline bigfloat_int
cot (const bigfloat_int & x)
{
	bigfloat_int y;

	cot(y, x);
	return y;
}



inline bigfloat_int
asin (const bigfloat_int & x)
{
	bigfloat_int y;

	asin(y, x);
	return y;
}



inline bigfloat_int
acos (const bigfloat_int & x)
{
	bigfloat_int y;

	acos(y, x);
	return y;
}



inline bigfloat_int
atan (const bigfloat_int & x)
{
	bigfloat_int y;

	atan(y, x);
	return y;
}



inline bigfloat_int
atan2 (const bigfloat_int & y, const bigfloat_int & x)
{
	bigfloat_int z;

	atan2(z, y, x);
	return z;
}



inline bigfloat_int
acot (const bigfloat_int & x)
{
	bigfloat_int y;

	acot(y, x);
	return y;
}



//
// (Inverse) hyperbolic trigonometrical procedures and functions
//

void sinh(bigfloat_int &, const bigfloat_int &);
void cosh(bigfloat_int &, const bigfloat_int &);
void tanh(bigfloat_int &, const bigfloat_int &);
void coth(bigfloat_int &, const bigfloat_int &);
void asinh(bigfloat_int &, const bigfloat_int &);
void acosh(bigfloat_int &, const bigfloat_int &);
void atanh(bigfloat_int &, const bigfloat_int &);
void acoth(bigfloat_int &, const bigfloat_int &);



inline bigfloat_int
sinh (const bigfloat_int & x)
{
	bigfloat_int y;

	sinh(y, x);
	return y;
}



inline bigfloat_int
cosh (const bigfloat_int & x)
{
	bigfloat_int y;

	cosh(y, x);
	return y;
}



inline bigfloat_int
tanh (const bigfloat_int & x)
{
	bigfloat_int y;

	tanh(y, x);
	return y;
}



inline bigfloat_int
coth (const bigfloat_int & x)
{
	bigfloat_int y;

	coth(y, x);
	return y;
}



inline bigfloat_int
asinh (const bigfloat_int & x)
{
	bigfloat_int y;

	asinh(y, x);
	return y;
}



inline bigfloat_int
acosh (const bigfloat_int & x)
{
	bigfloat_int y;

	acosh(y, x);
	return y;
}



inline bigfloat_int
atanh (const bigfloat_int & x)
{
	bigfloat_int y;

	atanh(y, x);
	return y;
}



inline bigfloat_int
acoth (const bigfloat_int & x)
{
	bigfloat_int y;

	acoth(y, x);
	return y;
}



//
// misc. procedures
//

void constant_E (bigfloat_int & x);
void constant_Pi (bigfloat_int & x);
void constant_Euler (bigfloat_int & x);
void constant_Catalan (bigfloat_int & x);



inline void
swap (bigfloat_int& a, bigfloat_int& b)
{
	a.swap(b);
}



//
// I/O
//

std::istream & operator >> (std::istream & in, bigfloat_int & a);
std::ostream & operator << (std::ostream & out, const bigfloat_int & a);



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_BIGFLOAT_INT_H_GUARD_
