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
//	Author	: Stefan Neis (SN)
//	Changes	: See CVS log
//
//==============================================================================================


//
// This  include  file  supports basic polynomial operations over
// all types
//
// Most   of   this   code (+, -, *, /, gcd, xgcd)  was
// originally  written by Victor Shoup (finite fields); Thomas Papanikolaou
// changed his implementation to use templates and the names of the
// functions to our standard. Then I extended this code to be able to
// treat polynomials over other types, especially over bigints
//
//                                             Stefan Neis
//


#ifndef LIDIA_FIELD_POLYNOMIAL_H_GUARD_
#define LIDIA_FIELD_POLYNOMIAL_H_GUARD_



#ifndef LIDIA_POLY_INTERN_H_GUARD_
# include	"LiDIA/base/poly_intern.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class bigint;



//
// The class in this file (field_polynomial) is for internal use only !!!
//
// The user is supposed to always use "polynomial < T >".
//

//
// debug level
//
//   0 : div_rem
//   1 : power_mod, integral
//   2 : gcd
//

template< class T >
class field_polynomial : public base_polynomial< T >
// Implements polynomials over fields.
// We require all functions needed for base_polynomial and
// additionally the functions divide(c, a, b) and invert(c, a).
{
public:

	//
	// constructors and destructor
	//

	field_polynomial(): base_polynomial< T > ()
	{ }

	field_polynomial(T x): base_polynomial< T > (x)
	{ }

	field_polynomial(const T * v, lidia_size_t d): base_polynomial< T > (v, d)
	{ }

	field_polynomial(const base_vector< T > v): base_polynomial< T > (v)
	{ }

	field_polynomial(const base_polynomial< T > &p): base_polynomial< T > (p)
	{ }

	~field_polynomial()
	{ }

	field_polynomial< T > &operator = (const base_polynomial< T > &a)
	{
		base_polynomial< T >::assign(a);
		return *this;
	}

	//
	// Division and related stuff
	//

	void div_rem(field_polynomial< T > &, const base_polynomial< T > &,
		     const base_polynomial< T > &);
	void divide(const base_polynomial< T > &, const T &);
	void power_mod(const base_polynomial< T > &, const bigint &,
		       const base_polynomial< T > &);
	void integral(const base_polynomial< T > &);
	void gcd(const base_polynomial< T > &, const base_polynomial< T > &);
	void xgcd(field_polynomial< T > &, field_polynomial< T > &,
		  const base_polynomial< T > &, const base_polynomial< T > &);

};



template< class T >
inline void
div_rem(field_polynomial< T > & q,
	field_polynomial< T > & r,
	const base_polynomial< T > & a,
	const base_polynomial< T > & b)
{
	q.div_rem(r, a, b);
}



template< class T >
inline void
divide(field_polynomial< T > & c,
       const base_polynomial< T > & a,
       const T & b)
{
	c.divide(a, b);
}



template< class T >
inline void
divide(field_polynomial< T > & q,
       const base_polynomial< T > & a,
       const base_polynomial< T > & b)
{
	field_polynomial< T > r;

	q.div_rem(r, a, b);
}



template< class T >
inline void
remainder(field_polynomial< T > & r,
	  const base_polynomial< T > & a,
	  const base_polynomial< T > & b)
{
	field_polynomial< T > q;

	q.div_rem(r, a, b);
}



template< class T >
inline void
power_mod(field_polynomial< T > & c,
	  const base_polynomial< T > & a,
	  const bigint & b,
	  const base_polynomial< T > & f)
{
	c.power_mod(a, b, f);
}



// template < class T >
// field_polynomial < T > operator / (const base_polynomial < T > &a,
//				    const base_polynomial < T > &b);

// template < class T >
// field_polynomial < T > operator / (const base_polynomial < T > &a,
//			    const T & b);

// template < class T >
// field_polynomial < T > operator % (const base_polynomial < T > &a,
//				    const base_polynomial < T > &b);

template< class T >
inline void
integral(field_polynomial< T > & c,
	 const base_polynomial< T > & a)
{
	c.integral(a);
}



//  template < class T >
//  field_polynomial < T > integral(const base_polynomial < T > & a);

template< class T >
inline void
gcd(field_polynomial< T > & d,
    const base_polynomial< T > & aa,
    const base_polynomial< T > & bb)
{
	d.gcd(aa, bb);
}



//  template < class T >
//  field_polynomial < T > gcd(const base_polynomial < T > &aa,
//			    const base_polynomial < T > &bb);

template< class T >
inline void
xgcd (field_polynomial< T > & d,
      field_polynomial< T > & x, field_polynomial< T > & y,
      const base_polynomial< T > & aa, const base_polynomial< T > & bb)
{
	d.xgcd(x, y, aa, bb);
}



//  template < class T >
//  field_polynomial < T > xgcd(field_polynomial < T > &x,
//			     field_polynomial < T > &y,
//			     const base_polynomial < T > &aa,
//			     const base_polynomial < T > &bb);




//
// operators
//

template< class T >
inline field_polynomial< T > &
operator /= (field_polynomial< T > & a, const base_polynomial< T > & b)
{
	field_polynomial< T > r;
	a.div_rem(r, a, b);
	return a;
}



template< class T >
inline field_polynomial< T > &
operator /= (field_polynomial< T > & a, const T & b)
{
	a.divide(a, b);
	return a;
}



template< class T >
field_polynomial< T > &
operator %= (field_polynomial< T > & a, const base_polynomial< T > & b)
{
	field_polynomial< T > q;

	q.div_rem(a, a, b);
	return a;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/field_polynomial.cc"
#endif



#endif	// LIDIA_FIELD_POLYNOMIAL_H_GUARD_
