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
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_SPARSE_POWER_SERIES_H_GUARD_
#define LIDIA_SPARSE_POWER_SERIES_H_GUARD_


#ifndef LIDIA_BASE_SPARSE_POWER_SERIES_H_GUARD_
# include	"LiDIA/finite_fields/base_sparse_power_series.h"
#endif
#ifndef LIDIA_BIGMOD_H_GUARD_
# include	"LiDIA/bigmod.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



// **************************************************************************
// *
// *    class name :  sparse_power_series < T >
// *
// **************************************************************************

// fix the old incorrect spelling:
#define sparse_power_serie sparse_power_series


template< class T >
class sparse_power_series : public base_sparse_power_series< T >
{

//protected :
public:

	//
	// ***** protected member functions *****
	//

	void multiply (const sparse_power_series< T > & a, const sparse_power_series< T > & b);
	void invert   (const sparse_power_series< T > & a);
	void square   (const sparse_power_series< T > & a);
	void power    (const sparse_power_series< T > & a, long n);
	void divide   (const sparse_power_series< T > & a, const sparse_power_series< T > & b);
	void divide   (const T & b, const sparse_power_series< T > & a);


public :

	//
	// ****  constructor/destructor functions    ******
	//

	sparse_power_series ();
	sparse_power_series (const T & a, lidia_size_t l);
	sparse_power_series (const base_vector< T > & a, lidia_size_t f);
	sparse_power_series (const base_sparse_power_series< T > & x);
	~sparse_power_series ();


	//
	//  *****  assignment operator  *****
	//

	void assign (const base_sparse_power_series< T > & x);
	sparse_power_series< T > & operator = (const base_sparse_power_series< T > & x);


};



//
// ***** arithmetical procedures *****
//

template< class T >
inline 	void
multiply (sparse_power_series< T > & res,
	  const sparse_power_series< T > & a,
	  const sparse_power_series< T > & b)
{
	res.multiply(a, b);
}



template< class T >
inline 	void
invert (sparse_power_series< T > & res,
	const sparse_power_series< T > & a)
{
	res.invert(a);
}



template< class T >
inline 	void
square (sparse_power_series< T > & res,
	const sparse_power_series< T > & a)
{
	res.square(a);
}



template< class T >
inline 	void
power (sparse_power_series< T > & res,
       const sparse_power_series< T > & a,
       long n)
{
	res.power(a, n);
}



template< class T >
inline 	void
divide (sparse_power_series< T > & res,
	const sparse_power_series< T > & a,
	const sparse_power_series< T > & b)
{
	res.divide(a, b);
}



template< class T >
inline 	void
divide (sparse_power_series< T > & res,
	const T & b,
	const sparse_power_series< T > & a)
{
	res.divide(b, a);
}



//
//  *****  arithmetical operators  *****
//

template< class T >
inline 	sparse_power_series< T >
operator * (const sparse_power_series< T > & a,
	    const sparse_power_series< T > & b)
{
	sparse_power_series< T > c;

	c.multiply(a, b);
	return c;
}



template< class T >
inline 	sparse_power_series< T > &
operator *= (sparse_power_series< T > & a,
	     const sparse_power_series< T > & b)
{
	a.multiply(a, b);
	return a;
}



template< class T >
inline 	sparse_power_series< T >
operator / (const sparse_power_series< T > & a,
	    const sparse_power_series< T > & b)
{
	sparse_power_series< T > c;

	c.divide(a, b);
	return c;
}



template< class T >
inline 	sparse_power_series< T >
operator / (const T & a,
	    const sparse_power_series< T > & b)
{
	sparse_power_series< T > c;

	c.divide(a, b);
	return c;
}



template< class T >
inline 	sparse_power_series< T > &
operator /= (sparse_power_series< T > & a,
	     const sparse_power_series< T > & b)
{
	a.divide(a, b);
	return a;
}



#ifndef NO_PSR_BIGMOD

//******************************************************************************
//*************** Specialization : sparse_power_series < bigmod > **************
//******************************************************************************

template<>
class sparse_power_series< bigmod > : public base_sparse_power_series< bigmod >
{

protected :

	//
	// ***** protected member functions *****
	//

	void multiply (const sparse_power_series< bigmod > & a, const sparse_power_series< bigmod > & b);
	void invert   (const sparse_power_series< bigmod > & a);
	void square   (const sparse_power_series< bigmod > & a);
	void power    (const sparse_power_series< bigmod > & a, long n);
	void divide   (const sparse_power_series< bigmod > & a, const sparse_power_series< bigmod > & b);
	void divide   (const bigmod & b, const sparse_power_series< bigmod > & a);


public :

	//
	// ****  constructor/destructor functions    ******
	//

	sparse_power_series ();
	sparse_power_series (const bigmod & a, lidia_size_t l);
	sparse_power_series (const base_vector< bigmod > & a, lidia_size_t f);
	sparse_power_series (const base_sparse_power_series< bigmod > & x);
	~sparse_power_series () {};


	//
	//  *****  assignment operator  *****
	//

	sparse_power_series< bigmod > & operator = (const base_sparse_power_series< bigmod > & x);


	//
	// ***** arithmetical procedures *****
	//

	friend 	void
	multiply (sparse_power_series< bigmod > & res,
		  const sparse_power_series< bigmod > & a,
		  const sparse_power_series< bigmod > & b);

	friend 	void
	invert (sparse_power_series< bigmod > & res,
		const sparse_power_series< bigmod > & a);

	friend 	void
	square (sparse_power_series< bigmod > & res,
		const sparse_power_series< bigmod > & a);

	friend 	void
	power (sparse_power_series< bigmod > & res,
	       const sparse_power_series< bigmod > & a,
	       long n);

	friend 	void
	divide (sparse_power_series< bigmod > & res,
		const sparse_power_series< bigmod > & a,
		const sparse_power_series< bigmod > & b);

	friend 	void
	divide (sparse_power_series< bigmod > & res,
		const bigmod & b,
		const sparse_power_series< bigmod > & a);


	//
	//  *****  arithmetical operators  *****
	//

	friend 	sparse_power_series< bigmod >
	operator * (const sparse_power_series< bigmod > & a,
		    const sparse_power_series< bigmod > & b);

	friend 	sparse_power_series< bigmod > &
	operator *= (sparse_power_series< bigmod > & a,
		     const sparse_power_series< bigmod > & b);

	friend 	sparse_power_series< bigmod >
	operator / (const sparse_power_series< bigmod > & a,
		    const sparse_power_series< bigmod > & b);

	friend 	sparse_power_series< bigmod >
	operator / (const bigmod & a,
		    const sparse_power_series< bigmod > & b);

	friend 	sparse_power_series< bigmod > &
	operator /= (sparse_power_series< bigmod > & a,
		     const sparse_power_series< bigmod > & b);

};



//
// c'tors and d'tor
//

template< class T >
inline
sparse_power_series< T >::sparse_power_series ()
	: base_sparse_power_series< T > ()
{
	debug_handler ("sparse_power_series< T >",
		       "sparse_power_series< T > ()");
}



template< class T >
inline
sparse_power_series< T >::sparse_power_series (const T & a, lidia_size_t l)
	: base_sparse_power_series< T > (a, l)
{
	debug_handler ("sparse_power_series< T >",
		       "sparse_power_series< T > (const T&, lidia_size_t)");
}



template< class T >
inline
sparse_power_series< T >::sparse_power_series (const base_vector< T > & a, lidia_size_t f)
	: base_sparse_power_series< T > (a, f)
{
	debug_handler ("sparse_power_series< T >",
		       "sparse_power_series< T > (const base_vector< T > &, lidia_size_t)");
}



template< class T >
inline
sparse_power_series< T >::sparse_power_series (const base_sparse_power_series< T > & x)
	: base_sparse_power_series< T > (x)
{
	debug_handler ("sparse_power_series< T >",
		       "sparse_power_series< T > (const base_sparse_power_series< T > &)");
}



template< class T >
inline
sparse_power_series< T >::~sparse_power_series ()
{
	debug_handler ("sparse_power_series< T >",
		       "~sparse_power_series< T > ()");
}



//
// assigners
//

template< class T >
inline void
sparse_power_series< T >::assign (const base_sparse_power_series< T > & x)
{
	debug_handler ("sparse_power_series< T >",
		       "assign (const base_sparse_power_series< T > &");

	if (&x != this)
		base_sparse_power_series< T >::assign(x);
}



template< class T >
inline sparse_power_series< T > &
sparse_power_series< T >::operator = (const base_sparse_power_series< T > & x)
{
	assign(x);
	return *this;
}



inline 	void
multiply (sparse_power_series< bigmod > & res,
	  const sparse_power_series< bigmod > & a,
	  const sparse_power_series< bigmod > & b)
{
	res.multiply(a, b);
}



inline 	void
invert (sparse_power_series< bigmod > & res,
	const sparse_power_series< bigmod > & a)
{
	res.invert(a);
}



inline 	void
square (sparse_power_series< bigmod > & res,
	const sparse_power_series< bigmod > & a)
{
	res.square(a);
}



inline 	void
power (sparse_power_series< bigmod > & res,
       const sparse_power_series< bigmod > & a,
       long n)
{
	res.power(a, n);
}



inline 	void
divide (sparse_power_series< bigmod > & res,
	const sparse_power_series< bigmod > & a,
	const sparse_power_series< bigmod > & b)
{
	res.divide(a, b);
}



inline 	void
divide (sparse_power_series< bigmod > & res,
	const bigmod & b,
	const sparse_power_series< bigmod > & a)
{
	res.divide(b, a);
}



//
//  *****  arithmetical operators  *****
//

inline 	sparse_power_series< bigmod >
operator * (const sparse_power_series< bigmod > & a,
	    const sparse_power_series< bigmod > & b)
{
	sparse_power_series< bigmod > c;

	c.multiply(a, b);
	return c;
}



inline 	sparse_power_series< bigmod > &
operator *= (sparse_power_series< bigmod > & a,
	     const sparse_power_series< bigmod > & b)
{
	a.multiply(a, b);
	return a;
}



inline 	sparse_power_series< bigmod >
operator / (const sparse_power_series< bigmod > & a,
	    const sparse_power_series< bigmod > & b)
{
	sparse_power_series< bigmod > c;

	c.divide(a, b);
	return c;
}



inline 	sparse_power_series< bigmod >
operator / (const bigmod & a,
	    const sparse_power_series< bigmod > & b)
{
	sparse_power_series< bigmod > c;

	c.divide(a, b);
	return c;
}



inline 	sparse_power_series< bigmod > &
operator /= (sparse_power_series< bigmod > & a,
	     const sparse_power_series< bigmod > & b)
{
	a.divide(a, b);
	return a;
}



#endif	// NO_PSR_BIGMOD



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/sparse_power_series.cc"
#endif



#endif	// LIDIA_SPARSE_POWER_SERIES_H_GUARD_
