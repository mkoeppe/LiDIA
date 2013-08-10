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
//	Author	: Frank Lehmann (FL), Markus Maurer (MM)
//                Thorsten Rottschaefer (TR)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_DENSE_POWER_SERIES_H_GUARD_
#define LIDIA_DENSE_POWER_SERIES_H_GUARD_


#ifndef LIDIA_BASE_DENSE_POWER_SERIES_H_GUARD_
# include	"LiDIA/finite_fields/base_dense_power_series.h"
#endif
#ifndef LIDIA_SPARSE_POWER_SERIES_H_GUARD_
# include	"LiDIA/sparse_power_series.h"
#endif
#ifndef LIDIA_FFT_PRIME_H_GUARD_
# include	"LiDIA/fft_prime.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



// fix the old incorrect spelling:
#define dense_power_serie dense_power_series

template< class T > class dense_power_series : public base_dense_power_series< T >
{
public:

	void square   (const dense_power_series< T > & a);
	void multiply (const dense_power_series< T > & a, const dense_power_series< T > & b);
	void invert   (const dense_power_series< T > & a);
	void power    (const dense_power_series< T > & a, long n);
	void divide   (const dense_power_series< T > & a, const dense_power_series< T > & b);
	void divide   (const T & b, const dense_power_series< T > & a);


public :

	//
	// c'tors and d'tor
	//

	dense_power_series ();
	dense_power_series (const T & a, lidia_size_t l);
	dense_power_series (const base_vector< T > & a, lidia_size_t f);
	dense_power_series (const base_dense_power_series< T > & a);

	~dense_power_series ();



	//
	// assigners
	//

	void assign (const dense_power_series< T > & a);
	void assign (const sparse_power_series< T > & a);
	dense_power_series< T > & operator = (const dense_power_series< T > & a);
	dense_power_series< T > & operator = (const sparse_power_series< T > & a);

};



//
// c'tors and d'tor
//

template< class T >
inline
dense_power_series< T >::dense_power_series ()
	: base_dense_power_series< T > ()
{
	debug_handler ("dense_power_series< T >",
		       "dense_power_series()");
}



template< class T >
inline
dense_power_series< T >::dense_power_series (const T & a, lidia_size_t l)
	: base_dense_power_series< T > (a, l)
{
	debug_handler ("dense_power_series< T >",
		       "dense_power_series(const T&, lidia_size_t)");
}



template< class T >
inline
dense_power_series< T >::dense_power_series (const base_vector< T > & a, lidia_size_t f)
	: base_dense_power_series< T > (a, f)
{
	debug_handler ("dense_power_series< T >",
		       "dense_power_series(const base_vector< T > &, lidia_size_t)");
}



template< class T >
inline
dense_power_series< T >::dense_power_series (const base_dense_power_series< T > & a)
	: base_dense_power_series< T > (a)
{
	debug_handler ("dense_power_series< T >",
		       "dense_power_series(const base_dense_power_series< T > &)");
}



template< class T >
inline
dense_power_series< T >::~dense_power_series ()
{
	debug_handler ("dense_power_series< T >",
		       "~dense_power_series()");
}



//
// assigners
//

template< class T >
inline void
dense_power_series< T >::assign (const dense_power_series< T > & a)
{
	debug_handler ("dense_power_series< T >",
		       "operator = (const dense_power_series< T > &)");

	if (&a != this) {
		base_dense_power_series< T >::assign(a);
	}
}



template< class T >
inline void
dense_power_series< T >::assign (const sparse_power_series< T > & a)
{
	debug_handler ("dense_power_series< T >",
		       "operator = (const sparse_power_series< T > &)");
	base_dense_power_series< T >::assign(a);
}



template< class T >
inline dense_power_series< T > &
dense_power_series< T >::operator = (const dense_power_series< T > & a)
{
	this->assign(a);
	return *this;
}



template< class T >
inline dense_power_series< T > &
dense_power_series< T >::operator = (const sparse_power_series< T > & a)
{
	this->assign(a);
	return *this;
}



//
// ***** arithmetic via functions *****
//

template< class T >
inline void
square (dense_power_series< T > & c, const dense_power_series< T > & a)
{
	c.square(a);
}



template< class T >
inline void
multiply (dense_power_series< T > & c, const dense_power_series< T > & a, const dense_power_series< T > & b)
{
	c.multiply(a, b);
}



template< class T >
inline void
invert (dense_power_series< T > & c, const dense_power_series< T > & a)
{
	c.invert(a);
}



template< class T >
inline void
power (dense_power_series< T > & c, const dense_power_series< T > & a, long n)
{
	c.power(a, n);
}



template< class T >
inline void
divide (dense_power_series< T > & c, dense_power_series< T > & a, dense_power_series< T > & b)
{
	c.divide(a, b);
}



template< class T >
inline void
divide (dense_power_series< T > & c, const T & b, const dense_power_series< T > & a)
{
	c.divide(b, a);
}



//
// ***** arithmetic via operators *****
//

template< class T >
inline dense_power_series< T >
operator * (const dense_power_series< T > & a,
	    const dense_power_series< T > & b)
{
	dense_power_series< T > c;

	c.multiply(a, b);
	return c;
}



template< class T >
inline dense_power_series< T > &
operator *= (dense_power_series< T > & a,
	     const dense_power_series< T > & b)
{
	a.multiply(a, b);
	return a;
}



template< class T >
inline dense_power_series< T >
operator / (const dense_power_series< T > & a,
	    const dense_power_series< T > & b)
{
	dense_power_series< T > c;

	c.divide(a, b);
	return c;
}



template< class T >
inline dense_power_series< T >
operator / (const T & b,
	    const dense_power_series< T > & a)
{
	dense_power_series< T > c;

	c.divide(b, a);
	return c;
}



template< class T >
inline dense_power_series< T > &
operator /= (dense_power_series< T > & a,
	     const dense_power_series< T > & b)
{
	a.divide(a, b);
	return a;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#define LIDIA_CLASS_DENSE_POWER_SERIES

#include	"LiDIA/specialization/dense_power_series.special"



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/dense_power_series.cc"
#endif



#endif	// LIDIA_DENSE_POWER_SERIES_H_GUARD_
