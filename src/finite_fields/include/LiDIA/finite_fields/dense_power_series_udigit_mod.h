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


#ifndef LIDIA_DENSE_POWER_SERIES_UDIGIT_MOD_H_GUARD_
#define LIDIA_DENSE_POWER_SERIES_UDIGIT_MOD_H_GUARD_


#ifndef LIDIA_UDIGIT_MOD_H_GUARD_
# include	"LiDIA/udigit_mod.h"
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



//#ifndef NO_PSR_BIGMOD

//*****************************************************************************
//************** Specialization : dense_power_series < udigit_mod > ***********
//*****************************************************************************

template<>
class dense_power_series< udigit_mod > : public base_dense_power_series< udigit_mod >
{
	static fft_prime q;
	static lidia_size_t fft_cross_over_point;

protected :

	//
	// ***** protected member functions *****
	//
	void square   (const dense_power_series< udigit_mod > & a);
	void multiply (const dense_power_series< udigit_mod > & a, const dense_power_series< udigit_mod > & b);
	void invert   (const dense_power_series< udigit_mod > & a);
	void power    (const dense_power_series< udigit_mod > & a, long n);
	void divide   (const dense_power_series< udigit_mod > & a, const dense_power_series< udigit_mod > & b);
	void divide   (const udigit_mod & b, const dense_power_series< udigit_mod > & a);

	void multiply_fft (const dense_power_series< udigit_mod > & a,
			   const dense_power_series< udigit_mod > & b);
	void multiply_plain (const dense_power_series< udigit_mod > & a,
			     const dense_power_series< udigit_mod > & b);


public :

	friend void square (dense_power_series< udigit_mod > & c,
			    const dense_power_series< udigit_mod > & a);
	friend void multiply (dense_power_series< udigit_mod > & c,
			      const dense_power_series< udigit_mod > & a,
			      const dense_power_series< udigit_mod > & b);
	friend void multiply_fft  (dense_power_series< udigit_mod > & c,
				   const dense_power_series< udigit_mod > & a,
				   const dense_power_series< udigit_mod > & b);
	friend void multiply_plain(dense_power_series< udigit_mod > & c,
				   const dense_power_series< udigit_mod > & a,
				   const dense_power_series< udigit_mod > & b);
	friend void invert (dense_power_series< udigit_mod > & c,
			    const dense_power_series< udigit_mod > & a);
	friend void power (dense_power_series< udigit_mod > & c,
			   const dense_power_series< udigit_mod > & a,
			   long n);
	friend void divide (dense_power_series< udigit_mod > & c,
			    dense_power_series< udigit_mod > & a,
			    dense_power_series< udigit_mod > & b);
	friend void divide (dense_power_series< udigit_mod > & c,
			    const udigit_mod & b,
			    const dense_power_series< udigit_mod > & a);


	friend dense_power_series< udigit_mod >
	operator * (const dense_power_series< udigit_mod > & a,
		    const dense_power_series< udigit_mod > & b);

	friend dense_power_series< udigit_mod > &
	operator *= (dense_power_series< udigit_mod > & a,
		     const dense_power_series< udigit_mod > & b);

	friend dense_power_series< udigit_mod >
	operator / (const dense_power_series< udigit_mod > & a,
		    const dense_power_series< udigit_mod > & b);

	friend dense_power_series< udigit_mod >
	operator/ (const udigit_mod & b,
		   const dense_power_series< udigit_mod > & a);

	friend dense_power_series< udigit_mod > &
	operator /= (dense_power_series< udigit_mod > & a,
		     const dense_power_series< udigit_mod > & b);


	//
	// ***** constructors / destructor *****
	//

	dense_power_series ();
	dense_power_series (const udigit_mod & a, lidia_size_t l);
	dense_power_series (const base_vector< udigit_mod > & a, lidia_size_t f);
	dense_power_series (const base_dense_power_series< udigit_mod > & a);
	~dense_power_series ();


	//
	// ***** assignment - operator *****
	//

	dense_power_series< udigit_mod > & operator = (const dense_power_series< udigit_mod > & a);
	dense_power_series< udigit_mod > & operator = (const sparse_power_series< udigit_mod > & a);


	static void set_fft_cross_over_point(lidia_size_t k)
	{
		dense_power_series< udigit_mod >::fft_cross_over_point = k;
	}

	static lidia_size_t get_fft_cross_over_point()
	{
		return dense_power_series< udigit_mod >::fft_cross_over_point;
	}
};



//
// ***** arithmetic via functions *****
//

inline
void
square (dense_power_series< udigit_mod > & c,
	const dense_power_series< udigit_mod > & a)
{
	c.square(a);
}



inline
void
multiply (dense_power_series< udigit_mod > & c,
	  const dense_power_series< udigit_mod > & a,
	  const dense_power_series< udigit_mod > & b)
{
	// call the intern fft function
	c.multiply(a, b);
}



inline
void
multiply_fft  (dense_power_series< udigit_mod > & c,
	       const dense_power_series< udigit_mod > & a,
	       const dense_power_series< udigit_mod > & b)
{
	// call the intern fft function
	c.multiply_fft(a, b);
}



inline
void
multiply_plain(dense_power_series< udigit_mod > & c,
	       const dense_power_series< udigit_mod > & a,
	       const dense_power_series< udigit_mod > & b)
{
	c.multiply_plain(a, b);
}



inline
void
invert (dense_power_series< udigit_mod > & c,
	const dense_power_series< udigit_mod > & a)
{
	c.invert(a);
}



inline
void
power (dense_power_series< udigit_mod > & c,
       const dense_power_series< udigit_mod > & a,
       long n)
{
	c.power(a, n);
}



inline
void
divide (dense_power_series< udigit_mod > & c,
	dense_power_series< udigit_mod > & a,
	dense_power_series< udigit_mod > & b)
{
	c.divide(a, b);
}



inline
void
divide (dense_power_series< udigit_mod > & c,
	const udigit_mod & b,
	const dense_power_series< udigit_mod > & a)
{
	c.divide(b, a);
}



//
// ***** arithmetic via operators *****
//

inline
dense_power_series< udigit_mod >
operator * (const dense_power_series< udigit_mod > & a,
	    const dense_power_series< udigit_mod > & b)
{
	dense_power_series< udigit_mod > c;

	c.multiply(a, b);
	return c;
}



inline
dense_power_series< udigit_mod > &
operator *= (dense_power_series< udigit_mod > & a,
	     const dense_power_series< udigit_mod > & b)
{
	a.multiply(a, b);
	return a;
}



inline
dense_power_series< udigit_mod >
operator / (const dense_power_series< udigit_mod > & a,
	    const dense_power_series< udigit_mod > & b)
{
	dense_power_series< udigit_mod > c;

	c.divide(a, b);
	return c;
}



inline
dense_power_series< udigit_mod >
operator/ (const udigit_mod & b,
	   const dense_power_series< udigit_mod > & a)
{
	dense_power_series< udigit_mod > c;

	c.divide(b, a);
	return c;
}



inline
dense_power_series< udigit_mod > &
operator /= (dense_power_series< udigit_mod > & a,
	     const dense_power_series< udigit_mod > & b)
{
	a.divide(a, b);
	return a;
}



//#endif



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_DENSE_POWER_SERIES_UDIGIT_MOD_H_GUARD_
