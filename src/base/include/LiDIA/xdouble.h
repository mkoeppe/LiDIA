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
//	Author	: Keith Briggs (KB)
//	Changes	: See CVS log
//
//==============================================================================================


/*
(Original Header)
C++ functions for doubledouble (i.e. double+double) precision.

These functions use techniques due to Dekker, Linnainmaa, Kahan, Knuth
and Priest.  I credit Kahan with the addition algorithm (the simplification
which permits the elimination of the tests and branches is due to Knuth);
Dekker and Linnainmaa with the multiply, divide, and square root routines,
and Priest for the initial transcription into C++.

A doubledouble x is represented as a pair of doubles, x.hi and x.lo,
such that the number represented by x is x.hi + x.lo, where

   |x.lo| <= 0.5*ulp(x.hi), (*)

and ulp(y) means `unit in the last place of y'.  For the software to
work correctly, IEEE Standard Arithmetic is sufficient.  That includes
just about every modern workstation.  Also sufficient is any platform
that implements arithmetic with correct rounding, i.e., given double
floating point numbers a and b, a op b is computed exactly and then
rounded to the nearest double.  The tie-breaking rule is not
important.

See:

T. J. Dekker
   Point Technique for Extending the Available Precision,
   Numer. Math. 18 (1971), pp. 224-242.
S. Linnainmaa
   Software for doubled-precision floating point computations
   ACM TOMS 7, 172-283 (1081).
D. Priest
  On properties of floating point arithmetics: numerical stability
  and the cost of accurate computations.   Ph.D. Diss, Berkeley 1992.

and more references in http://www.cs.wisc.edu/~shoup/ntl/quad_float.txt.
*/

#ifndef LIDIA_XDOUBLE_H_GUARD_
#define LIDIA_XDOUBLE_H_GUARD_


#ifndef LIDIA_LIDIA_H_GUARD_
# include	"LiDIA/LiDIA.h"
#endif
#include	<cmath>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



#ifndef __isinf
# define __isinf(x) ((x) != (x))
#endif

#ifdef x86
# define x86_FIX \
  unsigned short __old_cw, __new_cw; \
  asm volatile ("fnstcw %0":" = m" (__old_cw)); \
  __new_cw = (__old_cw & ~0x300) | 0x200; \
  asm volatile ("fldcw %0": :"m" (__new_cw));
# define END_x86_FIX  asm volatile ("fldcw %0": :"m" (__old_cw));
#else
# define x86_FIX
# define END_x86_FIX
#endif

//
// Some system use a copysign define
//
#ifdef copysign
# undef copysign
# undef HAVE_COPYSIGN
#endif

//
// 2^27 + 1, appropriate for IEEE double
//
//    static const double Split = 134217729.0L;

class xdouble
{
protected:
	double hi, lo;

	static const double Split;

	xdouble normalize();

public:

	//
	// Useful constants
	//
	static const xdouble Log2;
	static const xdouble Log10;
	static const xdouble Pi;
	static const xdouble TwoPi;
	static const xdouble Pion2;
	static const xdouble Pion4;
	static const xdouble _Pi;

	//
	// Public access to hi and lo
	//
	double h() const;
	double l() const;

	//
	// Constructors
	//
	inline xdouble();
	inline xdouble(const xdouble&);
	inline xdouble(double x, double y = 0.0);
	xdouble(char *);

	//
	// Conversions
	//
	friend double double_xdbl(const xdouble&);
	friend unsigned long ulong_xdbl(const xdouble&);
	friend long long_xdbl(const xdouble&);
	friend unsigned int uint_xdbl(const xdouble&);
	friend int int_xdbl(const xdouble&);
	friend unsigned short ushort_xdbl(const xdouble&);
	friend short short_xdbl(const xdouble&);
	friend unsigned char uchar_xdbl(const xdouble&);
	friend char char_xdbl(const xdouble&);

	//
	// Operators
	//
	xdouble& operator = (const char&);
	xdouble& operator = (const unsigned char&);
	xdouble& operator = (const short&);
	xdouble& operator = (const unsigned short&);
	xdouble& operator = (const int&);
	xdouble& operator = (const unsigned int&);
	xdouble& operator = (const long&);
	xdouble& operator = (const unsigned long&); // Get funny errors without this
	xdouble& operator = (const double&);
	xdouble& operator = (const xdouble&);
	xdouble& operator = (char*);

	//
	// Unary minus
	//
	friend xdouble operator - (const xdouble&);

	//
	// Addition
	//
	friend xdouble operator + (const xdouble&, const xdouble&);
	friend xdouble operator + (const double&, const xdouble&);
	friend xdouble operator + (const xdouble& x, const double& d);
	friend xdouble operator + (const unsigned long& ul, const xdouble& x);
	friend xdouble operator + (const xdouble& x, const unsigned long& ul);
	friend xdouble operator + (const long& l, const xdouble& x);
	friend xdouble operator + (const xdouble& x, const long& l);
	friend xdouble operator + (const unsigned int& ui, const xdouble& x);
	friend xdouble operator + (const xdouble& x, const unsigned int& ui);
	friend xdouble operator + (const int& i, const xdouble& x);
	friend xdouble operator + (const xdouble& x, const int& i);
	friend xdouble operator + (const unsigned short& us, const xdouble& x);
	friend xdouble operator + (const xdouble& x, const unsigned short& us);
	friend xdouble operator + (const short& s, const xdouble& x);
	friend xdouble operator + (const xdouble& x, const short& s);
	friend xdouble operator + (const unsigned char& uc, const xdouble& x);
	friend xdouble operator + (const xdouble& x, const unsigned char& uc);
	friend xdouble operator + (const char& c, const xdouble& x);
	friend xdouble operator + (const xdouble& x, const char& c);

	//
	// Assignment & Addition
	//
	xdouble& operator += (const xdouble& x);
	xdouble& operator += (const double& d);
	xdouble& operator += (const unsigned long& ul);
	xdouble& operator += (const long& l);
	xdouble& operator += (const unsigned int& ui);
	xdouble& operator += (const int& i);
	xdouble& operator += (const unsigned short& us);
	xdouble& operator += (const short& s);
	xdouble& operator += (const unsigned char& uc);
	xdouble& operator += (const char& c);

	//
	// Subtraction
	//
	friend xdouble operator - (const xdouble&, const xdouble&);
	friend xdouble operator - (const double&, const xdouble&);
	friend xdouble operator - (const xdouble& x, const double& d);
	friend xdouble operator - (const unsigned long& ul, const xdouble& x);
	friend xdouble operator - (const xdouble& x, const unsigned long& ul);
	friend xdouble operator - (const long& l, const xdouble& x);
	friend xdouble operator - (const xdouble& x, const long& l);
	friend xdouble operator - (const unsigned int& ui, const xdouble& x);
	friend xdouble operator - (const xdouble& x, const unsigned int& ui);
	friend xdouble operator - (const int& i, const xdouble& x);
	friend xdouble operator - (const xdouble& x, const int& i);
	friend xdouble operator - (const unsigned short& us, const xdouble& x);
	friend xdouble operator - (const xdouble& x, const unsigned short& us);
	friend xdouble operator - (const short& s, const xdouble& x);
	friend xdouble operator - (const xdouble& x, const short& s);
	friend xdouble operator - (const unsigned char& uc, const xdouble& x);
	friend xdouble operator - (const xdouble& x, const unsigned char& uc);
	friend xdouble operator - (const char& c, const xdouble& x);
	friend xdouble operator - (const xdouble& x, const char& c);

	//
	// Assignment & Subtraction
	//
	xdouble& operator -= (const xdouble& x);
	xdouble& operator -= (const double& d);
	xdouble& operator -= (const unsigned long& ul);
	xdouble& operator -= (const long& l);
	xdouble& operator -= (const unsigned int& ui);
	xdouble& operator -= (const int& i);
	xdouble& operator -= (const unsigned short& us);
	xdouble& operator -= (const short& s);
	xdouble& operator -= (const unsigned char& uc);
	xdouble& operator -= (const char& c);

	//
	// Multiplication
	//
	friend xdouble operator * (const xdouble&, const xdouble&);
	friend xdouble operator * (const double&, const xdouble&);
	friend xdouble operator * (const unsigned long&, const xdouble&);
	friend xdouble operator * (const long&, const xdouble&);
	friend xdouble operator * (const xdouble& x, const double& d);
	friend xdouble operator * (const xdouble& x, const unsigned long& ul);
	friend xdouble operator * (const xdouble& x, const long& l);
	friend xdouble operator * (const unsigned int& ui, const xdouble& x);
	friend xdouble operator * (const xdouble& x, const unsigned int& ui);
	friend xdouble operator * (const int& i, const xdouble& x);
	friend xdouble operator * (const xdouble& x, const int& i);
	friend xdouble operator * (const unsigned short& us, const xdouble& x);
	friend xdouble operator * (const xdouble& x, const unsigned short& us);
	friend xdouble operator * (const short& s, const xdouble& x);
	friend xdouble operator * (const xdouble& x, const short& s);
	friend xdouble operator * (const unsigned char& uc, const xdouble& x);
	friend xdouble operator * (const xdouble& x, const unsigned char& uc);
	friend xdouble operator * (const char& c, const xdouble& x);
	friend xdouble operator * (const xdouble& x, const char& c);

	//
	// Assignment & Multiplication
	//
	xdouble& operator *= (const xdouble& x);
	xdouble& operator *= (const double& d);
	xdouble& operator *= (const unsigned long& ul);
	xdouble& operator *= (const long& l);
	xdouble& operator *= (const unsigned int& ui);
	xdouble& operator *= (const int& i);
	xdouble& operator *= (const unsigned short& us);
	xdouble& operator *= (const short& s);
	xdouble& operator *= (const unsigned char& uc);
	xdouble& operator *= (const char& c);

	//
	// Division
	//
	friend xdouble operator / (const xdouble &, const xdouble &);
	friend xdouble operator / (const xdouble &, const double &);
	friend xdouble operator / (const double &, const xdouble &);
	friend xdouble operator / (const xdouble &, const long&);
	friend xdouble operator / (const xdouble& x, const unsigned long& ul);
	friend xdouble operator / (const unsigned long& ul, const xdouble& x);
	friend xdouble operator / (const long& l, const xdouble& x);
	friend xdouble operator / (const unsigned int& ui, const xdouble& x);
	friend xdouble operator / (const xdouble& x, const unsigned int& ui);
	friend xdouble operator / (const int& i, const xdouble& x);
	friend xdouble operator / (const xdouble& x, const int& i);
	friend xdouble operator / (const unsigned short& us, const xdouble& x);
	friend xdouble operator / (const xdouble& x, const unsigned short& us);
	friend xdouble operator / (const short& s, const xdouble& x);
	friend xdouble operator / (const xdouble& x, const short& s);
	friend xdouble operator / (const unsigned char& uc, const xdouble& x);
	friend xdouble operator / (const xdouble& x, const unsigned char& uc);
	friend xdouble operator / (const char& c, const xdouble& x);
	friend xdouble operator / (const xdouble& x, const char& c);

	//
	// Assignment & Division
	//
	xdouble& operator /= (const xdouble& x);
	xdouble& operator /= (const double& d);
	xdouble& operator /= (const unsigned long& ul);
	xdouble& operator /= (const long& l);
	xdouble& operator /= (const unsigned int& ui);
	xdouble& operator /= (const int& i);
	xdouble& operator /= (const unsigned short& us);
	xdouble& operator /= (const short& s);
	xdouble& operator /= (const unsigned char& uc);
	xdouble& operator /= (const char& c);

	//friend xdouble operator %(const xdouble&, const xdouble&);
	//friend xdouble operator | (const xdouble &, const xdouble &);

	//
	// Comparison
	//
	friend bool operator > (const xdouble &, const xdouble &);
	friend bool operator >= (const xdouble &, const xdouble &);
	friend bool operator < (const xdouble &, const xdouble &);
	friend bool operator <= (const xdouble &, const xdouble &);
	friend bool operator == (const xdouble &, const xdouble &);
	friend bool operator != (const xdouble &, const xdouble &);

	friend double dnorm(const xdouble &);

	//
	// member functions
	//
	void dump(char *); // debugging use only

	//
	// non inline functions (defined in xdouble.c)
	//
	// Input/Output
	//
	friend std::istream & operator >> (std::istream &, xdouble &);
	friend std::ostream & operator << (std::ostream &, const xdouble &);

	//
	// String conversion
	//
	friend void string_to_xdouble(const char*, xdouble&);
	friend void xdouble_to_string(const xdouble&, char*);

	//
	// for compatibility to double
	//
	friend xdouble hypot(const xdouble&, const xdouble&);
	friend xdouble recip(const xdouble&);
	friend xdouble sqrt(const xdouble&);
	friend xdouble sqr(const xdouble&);
	friend xdouble cub(const xdouble&);
	friend xdouble sqr_double(const double&);
	friend xdouble rint(const xdouble&);
	friend xdouble ceil(const xdouble&);
	friend xdouble floor(const xdouble&);
	friend xdouble trunc(const xdouble&);
	friend xdouble fmod(const xdouble&, const int);
	friend xdouble modf(const xdouble&, xdouble *ip);
	friend xdouble fabs(const xdouble&);
	friend xdouble exp(const xdouble&);
	friend xdouble log(const xdouble&);
	friend xdouble log10(const xdouble&);
	friend xdouble powint(const xdouble&, const long);
	friend xdouble pow(const xdouble&, const xdouble&);
	friend xdouble sin(const xdouble&);
	friend void sincos(const xdouble& x, xdouble& sinx, xdouble& cosx);
	friend xdouble cos(const xdouble&);
	friend xdouble atan(const xdouble&);
	friend xdouble atan2(const xdouble&, const xdouble&);
	friend xdouble asin(const xdouble&);
	friend xdouble sinh(const xdouble&);
	friend xdouble cosh(const xdouble&);
	friend xdouble tanh(const xdouble&);
	friend xdouble erf(const xdouble&);
	friend xdouble erfc(const xdouble);
	friend xdouble gamma(const xdouble&);
	friend int digits(const xdouble&, const xdouble&);
	friend xdouble modr(const xdouble& a, const xdouble& b, int& n, xdouble& rem);
	friend xdouble ldexp(const xdouble&, const int);
	friend int sign(const xdouble&);

#ifndef HAVE_COPYSIGN
	friend double copysign(double x, double y);
#endif
	friend xdouble copysign(const xdouble& x, const double y);
}; // end class xdouble


double double_xdbl(const xdouble&);
unsigned long ulong_xdbl(const xdouble&);
long long_xdbl(const xdouble&);
unsigned int uint_xdbl(const xdouble&);
int int_xdbl(const xdouble&);
unsigned short ushort_xdbl(const xdouble&);
short short_xdbl(const xdouble&);
unsigned char uchar_xdbl(const xdouble&);
char char_xdbl(const xdouble&);

//
// Unary minus
//
xdouble operator - (const xdouble&);

//
// Addition
//
xdouble operator + (const xdouble&, const xdouble&);
xdouble operator + (const double&, const xdouble&);
xdouble operator + (const xdouble& x, const double& d);
xdouble operator + (const unsigned long& ul, const xdouble& x);
xdouble operator + (const xdouble& x, const unsigned long& ul);
xdouble operator + (const long& l, const xdouble& x);
xdouble operator + (const xdouble& x, const long& l);
xdouble operator + (const unsigned int& ui, const xdouble& x);
xdouble operator + (const xdouble& x, const unsigned int& ui);
xdouble operator + (const int& i, const xdouble& x);
xdouble operator + (const xdouble& x, const int& i);
xdouble operator + (const unsigned short& us, const xdouble& x);
xdouble operator + (const xdouble& x, const unsigned short& us);
xdouble operator + (const short& s, const xdouble& x);
xdouble operator + (const xdouble& x, const short& s);
xdouble operator + (const unsigned char& uc, const xdouble& x);
xdouble operator + (const xdouble& x, const unsigned char& uc);
xdouble operator + (const char& c, const xdouble& x);
xdouble operator + (const xdouble& x, const char& c);

//
// Subtraction
//
xdouble operator - (const xdouble&, const xdouble&);
xdouble operator - (const double&, const xdouble&);
xdouble operator - (const xdouble& x, const double& d);
xdouble operator - (const unsigned long& ul, const xdouble& x);
xdouble operator - (const xdouble& x, const unsigned long& ul);
xdouble operator - (const long& l, const xdouble& x);
xdouble operator - (const xdouble& x, const long& l);
xdouble operator - (const unsigned int& ui, const xdouble& x);
xdouble operator - (const xdouble& x, const unsigned int& ui);
xdouble operator - (const int& i, const xdouble& x);
xdouble operator - (const xdouble& x, const int& i);
xdouble operator - (const unsigned short& us, const xdouble& x);
xdouble operator - (const xdouble& x, const unsigned short& us);
xdouble operator - (const short& s, const xdouble& x);
xdouble operator - (const xdouble& x, const short& s);
xdouble operator - (const unsigned char& uc, const xdouble& x);
xdouble operator - (const xdouble& x, const unsigned char& uc);
xdouble operator - (const char& c, const xdouble& x);
xdouble operator - (const xdouble& x, const char& c);

//
// Multiplication
//
xdouble operator * (const xdouble&, const xdouble&);
xdouble operator * (const double&, const xdouble&);
xdouble operator * (const unsigned long&, const xdouble&);
xdouble operator * (const long&, const xdouble&);
xdouble operator * (const xdouble& x, const double& d);
xdouble operator * (const xdouble& x, const unsigned long& ul);
xdouble operator * (const xdouble& x, const long& l);
xdouble operator * (const unsigned int& ui, const xdouble& x);
xdouble operator * (const xdouble& x, const unsigned int& ui);
xdouble operator * (const int& i, const xdouble& x);
xdouble operator * (const xdouble& x, const int& i);
xdouble operator * (const unsigned short& us, const xdouble& x);
xdouble operator * (const xdouble& x, const unsigned short& us);
xdouble operator * (const short& s, const xdouble& x);
xdouble operator * (const xdouble& x, const short& s);
xdouble operator * (const unsigned char& uc, const xdouble& x);
xdouble operator * (const xdouble& x, const unsigned char& uc);
xdouble operator * (const char& c, const xdouble& x);
xdouble operator * (const xdouble& x, const char& c);

//
// Division
//
xdouble operator / (const xdouble &, const xdouble &);
xdouble operator / (const xdouble &, const double &);
xdouble operator / (const double &, const xdouble &);
xdouble operator / (const xdouble &, const long&);
xdouble operator / (const xdouble& x, const unsigned long& ul);
xdouble operator / (const unsigned long& ul, const xdouble& x);
xdouble operator / (const long& l, const xdouble& x);
xdouble operator / (const unsigned int& ui, const xdouble& x);
xdouble operator / (const xdouble& x, const unsigned int& ui);
xdouble operator / (const int& i, const xdouble& x);
xdouble operator / (const xdouble& x, const int& i);
xdouble operator / (const unsigned short& us, const xdouble& x);
xdouble operator / (const xdouble& x, const unsigned short& us);
xdouble operator / (const short& s, const xdouble& x);
xdouble operator / (const xdouble& x, const short& s);
xdouble operator / (const unsigned char& uc, const xdouble& x);
xdouble operator / (const xdouble& x, const unsigned char& uc);
xdouble operator / (const char& c, const xdouble& x);
xdouble operator / (const xdouble& x, const char& c);

//
// Comparison
//
bool operator > (const xdouble &, const xdouble &);
bool operator >= (const xdouble &, const xdouble &);
bool operator < (const xdouble &, const xdouble &);
bool operator <= (const xdouble &, const xdouble &);
bool operator == (const xdouble &, const xdouble &);
bool operator != (const xdouble &, const xdouble &);

double dnorm(const xdouble &);

//
// Input/Output
//
std::istream & operator >> (std::istream &, xdouble &);
std::ostream & operator << (std::ostream &, const xdouble &);

//
// String conversion
//
void string_to_xdouble(const char*, xdouble&);
void xdouble_to_string(const xdouble&, char*);

//
// for compatibility to double
//
xdouble hypot(const xdouble&, const xdouble&);
xdouble recip(const xdouble&);
xdouble sqrt(const xdouble&);
xdouble sqr(const xdouble&);
xdouble cub(const xdouble&);
xdouble sqr_double(const double&);
xdouble rint(const xdouble&);
xdouble ceil(const xdouble&);
xdouble floor(const xdouble&);
xdouble trunc(const xdouble&);
xdouble fmod(const xdouble&, const int);
xdouble modf(const xdouble&, xdouble *ip);
xdouble fabs(const xdouble&);
xdouble exp(const xdouble&);
xdouble log(const xdouble&);
xdouble log10(const xdouble&);
xdouble powint(const xdouble&, const long);
xdouble pow(const xdouble&, const xdouble&);
xdouble sin(const xdouble&);
void sincos(const xdouble& x, xdouble& sinx, xdouble& cosx);
xdouble cos(const xdouble&);
xdouble atan(const xdouble&);
xdouble atan2(const xdouble&, const xdouble&);
xdouble asin(const xdouble&);
xdouble sinh(const xdouble&);
xdouble cosh(const xdouble&);
xdouble tanh(const xdouble&);
xdouble erf(const xdouble&);
xdouble erfc(const xdouble);
xdouble gamma(const xdouble&);
int digits(const xdouble&, const xdouble&);
xdouble modr(const xdouble& a, const xdouble& b, int& n, xdouble& rem);
xdouble ldexp(const xdouble&, const int);
int sign(const xdouble&);

#ifndef HAVE_COPYSIGN
double copysign(double x, double y);
#endif
xdouble copysign(const xdouble& x, const double y);

//
// inline members
//
inline double xdouble::h() const
{
	return (hi);
}

inline double xdouble::l() const
{
	return (lo);
}

//
// Constructors
//
inline xdouble::xdouble()
{
	hi = 0.0;
	lo = 0.0;
}

inline xdouble::xdouble(double x, double y)
{
	x86_FIX
	hi = x+y; lo = y+(x-hi); // normalize
	END_x86_FIX
}

inline xdouble::xdouble(const xdouble& x)
{
	hi = x.hi;
	lo = x.lo;
}

//
// Assignments
//
inline xdouble& xdouble::operator = (const char& c)
{
	hi = static_cast<double>(c);
	lo = 0.0;
	return (*this);
}

inline xdouble& xdouble::operator = (const unsigned char& uc)
{
	hi = static_cast<double>(uc);
	lo = 0.0;
	return (*this);
}

inline xdouble& xdouble::operator = (const short& s)
{
	hi = static_cast<double>(s);
	lo = 0.0;
	return (*this);
}

inline xdouble& xdouble::operator = (const unsigned short& us)
{
	hi = static_cast<double>(us);
	lo = 0.0;
	return (*this);
}

inline xdouble& xdouble::operator = (const int& i)
{
	hi = static_cast<double>(i);
	lo = 0.0;
	return (*this);
}

inline xdouble& xdouble::operator = (const unsigned int& ui)
{
	hi = static_cast<double>(ui);
	lo = 0.0;
	return (*this);
}

inline xdouble& xdouble::operator = (const long& l)
{
	hi = static_cast<double>(l);
	lo = 0.0;
	return (*this);
}

inline xdouble& xdouble::operator = (const unsigned long& ul)
{
	hi = static_cast<double>(ul);
	lo = 0.0;
	return (*this);
}

inline xdouble& xdouble::operator = (const double &d)
{
	hi = d;
	lo = 0.0;
	return (*this);
}

inline xdouble& xdouble::operator = (const xdouble& x)
{
	hi = x.hi;
	lo = x.lo;
	return (*this);
}

inline xdouble xdouble::normalize()
{
	double h = hi+lo;

	lo = lo+(hi-h);
	hi = h;
	return (*this);
}


#ifdef XDBL_NO_INLINE
# define XDBL_INLINE
#else
# define XDBL_INLINE inline
#endif

#ifndef XDBL_NO_INLINE
# include	"LiDIA/xdouble.inl"
#endif

//
// inline functions
//
inline double dnorm(const xdouble & x)
{
	return std::fabs(x.h());
}

//
// Conversions
//
inline double double_xdbl(const xdouble& xd)
{
	return (xd.hi);
}

inline unsigned long ulong_xdbl(const xdouble& xd)
{
	return (static_cast<unsigned long>(xd.hi));
}

inline long long_xdbl(const xdouble& xd)
{
	return (static_cast<long>(xd.hi));
}

inline unsigned int uint_xdbl(const xdouble& xd)
{
	return (static_cast<unsigned int>(xd.hi));
}

inline int int_xdbl(const xdouble& xd)
{
	return (static_cast<int>(xd.hi));
}

inline unsigned short ushort_xdbl(const xdouble& xd)
{
	return (static_cast<unsigned short>(xd.hi));
}

inline short short_xdbl(const xdouble& xd)
{
	return (static_cast<short>(xd.hi));
}

inline unsigned char uchar_xdbl(const xdouble& xd)
{
	return (static_cast<unsigned char>(xd.hi));
}

inline char char_xdbl(const xdouble& xd)
{
	return (static_cast<char>(xd.hi));
}

#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_XDOUBLE_H_GUARD_
