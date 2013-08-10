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
//                Thorsten Rottschaefer (TR), Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_DENSE_POWER_SERIES_BIGMOD_H_GUARD_
#define LIDIA_DENSE_POWER_SERIES_BIGMOD_H_GUARD_


#ifndef LIDIA_BIGMOD_H_GUARD_
# include	"LiDIA/bigmod.h"
#endif
#ifndef LIDIA_BASE_DENSE_POWER_SERIES_H_GUARD_
# include	"LiDIA/finite_fields/base_dense_power_series.h"
#endif
#ifndef LIDIA_DENSE_POWER_SERIES_H_GUARD_
# include	"LiDIA/dense_power_series.h"
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



#ifndef NO_PSR_BIGMOD


//******************************************************************************
//*************** Specialization : dense_power_series < bigmod > ***************
//******************************************************************************

template<>
class dense_power_series< bigmod > : public base_dense_power_series< bigmod >
{

protected :

	//
	// ***** protected member functions *****
	//

	void square   (const dense_power_series< bigmod > & a);
	void multiply (const dense_power_series< bigmod > & a, const dense_power_series< bigmod > & b);
	void invert   (const dense_power_series< bigmod > & a);
	void power    (const dense_power_series< bigmod > & a, long n);
	void divide   (const dense_power_series< bigmod > & a, const dense_power_series< bigmod > & b);
	void divide   (const bigmod & b, const dense_power_series< bigmod > & a);


public :

	//
	// ***** constructors / destructor *****
	//

	dense_power_series ();
	dense_power_series (const bigmod & a, lidia_size_t l);
	dense_power_series (const base_vector< bigmod > & a, lidia_size_t f);
	dense_power_series (const base_dense_power_series< bigmod > & a);
	~dense_power_series () {}


	//
	// ***** assignment - operator *****
	//

	dense_power_series< bigmod > & operator = (const dense_power_series< bigmod > & a);
	dense_power_series< bigmod > & operator = (const sparse_power_series< bigmod > & a);


	//
	// ***** arithmetic via functions *****
	//

	friend void square (dense_power_series< bigmod > & c,
			    const dense_power_series< bigmod > & a);

	friend void multiply (dense_power_series< bigmod > & c,
			      const dense_power_series< bigmod > & a,
			      const dense_power_series< bigmod > & b);

	friend void invert (dense_power_series< bigmod > & c,
			    const dense_power_series< bigmod > & a);

	friend void power (dense_power_series< bigmod > & c,
			   const dense_power_series< bigmod > & a,
			   long n);

	friend void divide (dense_power_series< bigmod > & c,
			    dense_power_series< bigmod > & a,
			    dense_power_series< bigmod > & b);

	friend void divide (dense_power_series< bigmod > & c,
			    const bigmod & b,
			    const dense_power_series< bigmod > & a);


	//
	// ***** arithmetic via operators *****
	//

	friend dense_power_series< bigmod >
	operator * (const dense_power_series< bigmod > & a,
		    const dense_power_series< bigmod > & b);

	friend dense_power_series< bigmod > &
	operator *= (dense_power_series< bigmod > & a,
		     const dense_power_series< bigmod > & b);

	friend dense_power_series< bigmod >
	operator / (const dense_power_series< bigmod > & a,
		    const dense_power_series< bigmod > & b);

	friend dense_power_series< bigmod >
	operator / (const bigmod & b,
		    const dense_power_series< bigmod > & a);

	friend dense_power_series< bigmod > &
	operator /= (dense_power_series< bigmod > & a,
		     const dense_power_series< bigmod > & b);

};



inline void
square (dense_power_series< bigmod > & c,
	const dense_power_series< bigmod > & a)
{
	c.square(a);
}



inline void
multiply (dense_power_series< bigmod > & c,
	  const dense_power_series< bigmod > & a,
	  const dense_power_series< bigmod > & b)
{
	c.multiply(a, b);
}



inline void
invert (dense_power_series< bigmod > & c,
	const dense_power_series< bigmod > & a)
{
	c.invert(a);
}



inline void
power (dense_power_series< bigmod > & c,
       const dense_power_series< bigmod > & a,
       long n)
{
	c.power(a, n);
}



inline void
divide (dense_power_series< bigmod > & c,
	dense_power_series< bigmod > & a,
	dense_power_series< bigmod > & b)
{
	c.divide(a, b);
}



inline void
divide (dense_power_series< bigmod > & c,
	const bigmod & b,
	const dense_power_series< bigmod > & a)
{
	c.divide(b, a);
}



//
// ***** arithmetic via operators *****
//

inline dense_power_series< bigmod >
operator * (const dense_power_series< bigmod > & a,
	   const dense_power_series< bigmod > & b)
{
	dense_power_series< bigmod > c;

	c.multiply(a, b);
	return c;
}



inline dense_power_series< bigmod > &
operator *= (dense_power_series< bigmod > & a,
	     const dense_power_series< bigmod > & b)
{
	a.multiply(a, b);
	return a;
}



inline dense_power_series< bigmod >
operator / (const dense_power_series< bigmod > & a,
	    const dense_power_series< bigmod > & b)
{
	dense_power_series< bigmod > c;

	c.divide(a, b);
	return c;
}



inline dense_power_series< bigmod >
operator / (const bigmod & b,
	    const dense_power_series< bigmod > & a)
{
	dense_power_series< bigmod > c;

	c.divide(b, a);
	return c;
}



inline dense_power_series< bigmod > &
operator /= (dense_power_series< bigmod > & a,
	     const dense_power_series< bigmod > & b)
{
	a.divide(a, b);
	return a;
}



#endif	// NO_PSR_BIGMOD



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_DENSE_POWER_SERIES_BIGMOD_H_GUARD_
