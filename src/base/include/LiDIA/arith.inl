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
//	Author	: Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_ARITH_INL_GUARD_
#define LIDIA_ARITH_INL_GUARD_



#ifndef LIDIA_XDOUBLE_H_GUARD_
# include	<LiDIA/xdouble.h>
#endif


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



//
// swap
//

inline void swap(int &a, int &b)
{
	int TMP;

	TMP = a;
	a = b;
	b = TMP;
}



inline void swap(unsigned int &a, unsigned int &b)
{
	int TMP;

	TMP = a;
	a = b;
	b = TMP;
}



inline void swap(long &a, long &b)
{
	long TMP;

	TMP = a;
	a = b;
	b = TMP;
}



inline void swap(unsigned long &a, unsigned long &b)
{
	long TMP;

	TMP = a;
	a = b;
	b = TMP;
}



inline void swap(float &a, float &b)
{
	float TMP;

	TMP = a;
	a = b;
	b = TMP;
}



inline void swap(double &a, double &b)
{
	double TMP;

	TMP = a;
	a = b;
	b = TMP;
}



inline void swap(xdouble &a, xdouble &b)
{
	xdouble TMP;

	TMP = a;
	a = b;
	b = TMP;
}



inline void swap(char &a, char &b)
{
	char TMP;

	TMP = a;
	a = b;
	b = TMP;
}



//
// addition
//

inline void add(int &c, int a, int b)
{
	c = a + b;
}



inline void add(long &c, long a, long b)
{
	c = a + b;
}



inline void add(float &c, float a, float b)
{
	c = a + b;
}



inline void add(double &c, double a, double b)
{
	c = a + b;
}



inline void add(xdouble &c, xdouble a, xdouble b)
{
	c = a + b;
}



//
// subtraction
//

inline void subtract(int &c, int a, int b)
{
	c = a - b;
}



inline void subtract(long &c, long a, long b)
{
	c = a - b;
}



inline void subtract(float &c, float a, float b)
{
	c = a - b;
}



inline void subtract(double &c, double a, double b)
{
	c = a - b;
}



inline void subtract(xdouble &c, xdouble a, xdouble b)
{
	c = a - b;
}



//
// multiplication
//

inline void multiply(int &c, int a, int b)
{
	c = a * b;
}



inline void multiply(long &c, long a, long b)
{
	c = a * b;
}



inline void multiply(float &c, float a, float b)
{
	c = a * b;
}



inline void multiply(double &c, double a, double b)
{
	c = a * b;
}



inline void multiply(xdouble &c, xdouble a, xdouble b)
{
	c = a * b;
}



//
// negation
//

inline void negate(int &c, int a)
{
	c = -a;
}



inline void negate(long &c, long a)
{
	c = -a;
}



inline void negate(float &c, float a)
{
	c = -a;
}



inline void negate(double &c, double a)
{
	c = -a;
}



inline void negate(xdouble &c, xdouble a)
{
	c = -a;
}



//
// divide
//

inline void divide(int &c, int a, int b)
{
	c = a / b;
}



inline void divide(long &c, long a, long b)
{
	c = a / b;
}



inline void divide(float &c, float a, float b)
{
	c = a / b;
}



inline void divide(double &c, double a, double b)
{
	c = a / b;
}



inline void divide(xdouble &c, xdouble a, xdouble b)
{
	c = a / b;
}


//
// binary logarithm
//

// some platforms define log2 as a macro that breaks LiDIA.
#ifdef log2
# undef log2
#endif

inline double log2(double a) {
  return std::log(a) / std::log(2.0);
}


//
//  is_zero
//

inline bool is_zero (int i)
{
	return (i == 0);
}



inline bool is_zero (long i)
{
	return (i == 0);
}



//
//  is_one
//

inline bool is_one (int i)
{
	return (i == 1);
}



inline bool is_one (long i)
{
	return (i == 1);
}



//
//  is_even
//

inline bool is_even (int i)
{
	return (!(i&1));
}



inline bool is_even (long i)
{
	return (!(i&1));
}



//
//  is_odd
//

inline bool is_odd (int i)
{
	return (i&1);
}



inline bool is_odd (long i)
{
	return (i&1);
}



//
//  is_negative
//

inline bool is_negative (int i)
{
	return (i < 0);
}



inline bool is_negative (long i)
{
	return (i < 0);
}



//
//  abs_compare
//

template< class T >
inline
int abs_compare(const T & a, const T & b)
{
	return a.abs_compare(b);
}



inline int abs_compare(long a, long b)
{
	if (a < 0) a = -a;
	if (b < 0) b = -b;
	if (a < b)
		return -1;
	else {
		if (a == b)
			return 0;
		else
			return 1;
	}
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_ARITH_INL_GUARD_
