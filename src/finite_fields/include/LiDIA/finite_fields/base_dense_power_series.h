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


#ifndef LIDIA_BASE_DENSE_POWER_SERIES_H_GUARD_
#define LIDIA_BASE_DENSE_POWER_SERIES_H_GUARD_


#ifndef LIDIA_MATH_VECTOR_H_GUARD_
# include	"LiDIA/math_vector.h"
#endif
#ifndef LIDIA_BASE_SPARSE_POWER_SERIES_H_GUARD_
# include	"LiDIA/finite_fields/base_sparse_power_series.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



// fix the old incorrect spelling:
#define base_dense_power_serie base_dense_power_series


template< class T > class base_dense_power_series
{

protected :

	lidia_size_t    first;
	lidia_size_t    last;
	math_vector< T > *coeff;

	lidia_size_t    max_num_coeff;


	//
	// ***** protected member functions *****
	//

public:

	void add (const base_dense_power_series< T > & a, const base_dense_power_series< T > & b);
	void add (const base_dense_power_series< T > & a, const T & b);
	void add (const T & b, const base_dense_power_series< T > & a);

	void subtract (const base_dense_power_series< T > & a, const base_dense_power_series< T > & b);
	void subtract (const base_dense_power_series< T > & a, const T & b);
	void subtract (const T & b, const base_dense_power_series< T > & a);

	void multiply (const base_dense_power_series< T > & a, const T & b);
	void multiply (const T & b, const base_dense_power_series< T > & a);

	void divide (const base_dense_power_series< T > & a, const T & b);
	void negate (const base_dense_power_series< T > & a);

	void swap (base_dense_power_series< T > & a);


public :

	//
	// c'tors and d'tor
	//

	base_dense_power_series ();
	base_dense_power_series (const T & elem, lidia_size_t l);
	base_dense_power_series (const base_vector< T > & c, lidia_size_t f);
	base_dense_power_series (const base_dense_power_series< T > & a);

	~base_dense_power_series ();


	//
	// ***** utility routines *****
	//

protected :

	bool is_zero (lidia_size_t & non_zero_index) const;


public :

	bool is_zero () const;
	bool is_one () const;

	bool is_equal (const base_dense_power_series< T > & a) const;


protected :

	void allocate (lidia_size_t cp, math_vector< T > * c);


public :

	lidia_size_t get_first () const;
	lidia_size_t get_last  () const;
	void  reduce_last (lidia_size_t l);



	//
	// ***** coefficient handling *****
	//

	void get_coeff (T & elem, lidia_size_t e);
	void get (base_vector< T > & c);
	void get (T * & c, lidia_size_t & sz);
	void set_coeff (const T & elem, lidia_size_t e);
	void set_coeff (const T & elem, lidia_size_t e, const base_dense_power_series< T > & a);
	void set (base_vector< T > & c, lidia_size_t f);
	void set (const T *c, lidia_size_t prec, lidia_size_t f);
	void set (const T & elem, lidia_size_t l);


	//
	// ***** more utility functions *****
	//

	void clear ();
	void normalize ();
	void set_max_num_of_coeff (lidia_size_t sz);
	lidia_size_t get_max_num_of_coeff ();
	void  clear_max_num_of_coeff ();


	//
	// ***** subscripting *****
	//

	const T operator[] (lidia_size_t e) const;
	T &     operator() (lidia_size_t e);


	//
	// ***** assignment *****
	//

	void assign_zero (lidia_size_t f);
	void assign_one (lidia_size_t l);

	void assign (const base_dense_power_series< T > & a);
	void assign (const base_sparse_power_series< T > & a);

	base_dense_power_series< T > & operator = (const base_dense_power_series< T > & a);
	base_dense_power_series< T > & operator = (const base_sparse_power_series< T > & a);


	//
	// ***** cast - operator *****
	//

	operator base_sparse_power_series< T > ();


	//
	// ***** input / output *****
	//

	int read (std::istream & in);
	void write (std::ostream & out) const;



	//
	// ***** miscellaneous *****
	//

	void multiply_by_xn (lidia_size_t n);
	void compose (lidia_size_t n);

	//friend void randomize (base_dense_power_series < T > & a,
	//			    lidia_size_t f, lidia_size_t l,
	//			    const T & coeff_bound);


};



template< class T >
inline base_dense_power_series< T > &
base_dense_power_series< T >::operator = (const base_dense_power_series< T > & a)
{
	assign(a);
	return *this;
}



template< class T >
inline base_dense_power_series< T > &
base_dense_power_series< T >::operator = (const base_sparse_power_series< T > & a)
{
	assign(a);
	return *this;
}



template< class T >
inline bool
operator == (const base_dense_power_series< T > & a,
	     const base_dense_power_series< T > & b)
{
	return a.is_equal(b);
}



template< class T >
inline bool
operator != (const base_dense_power_series< T > & a,
	     const base_dense_power_series< T > & b)
{
	return !a.is_equal(b);
}



template< class T >
inline base_dense_power_series< T >
operator+ (const base_dense_power_series< T > & a,
	   const base_dense_power_series< T > & b)
{
	base_dense_power_series< T > c;

	c.add(a, b);
	return(c);
}



template< class T >
inline base_dense_power_series< T >
operator+ (const base_dense_power_series< T > & a,
	   const T & b)
{
	base_dense_power_series< T > c;
	c.add(a, b);
	return(c);
}



template< class T >
inline base_dense_power_series< T >
operator+ (const T & b,
	   const base_dense_power_series< T > & a)
{
	base_dense_power_series< T > c;

	c.add(b, a);
	return(c);
}



template< class T >
inline base_dense_power_series< T > &
operator += (base_dense_power_series< T > & a,
	     const base_dense_power_series< T > & b)
{
	a.add(a, b);
	return (a);
}



template< class T >
inline base_dense_power_series< T > &
operator += (base_dense_power_series< T > & a,
	     const T & b)
{
	a.add(a, b);
	return (a);
}



template< class T >
inline base_dense_power_series< T >
operator- (const base_dense_power_series< T > & a,
	   const base_dense_power_series< T > & b)
{
	base_dense_power_series< T > c;

	c.subtract (a, b);
	return(c);
}



template< class T >
inline base_dense_power_series< T >
operator- (const base_dense_power_series< T > & a,
	   const T & b)
{
	base_dense_power_series< T > c;

	c.subtract (a, b);
	return(c);
}



template< class T >
inline base_dense_power_series< T >
operator- (const T & b,
	   const base_dense_power_series< T > & a)
{
	base_dense_power_series< T > c;

	c.subtract(b, a);
	return(c);
}



template< class T >
inline base_dense_power_series< T > &
operator -= (base_dense_power_series< T > & a,
	     const base_dense_power_series< T > & b)
{
	a.subtract(a, b);
	return(a);
}



template< class T >
inline base_dense_power_series< T > &
operator -= (base_dense_power_series< T > & a,
	     const T & b)
{
	a.subtract(a, b);
	return(a);
}



template< class T >
inline base_dense_power_series< T >
operator* (const base_dense_power_series< T > & a,
	   const T & b)
{
	base_dense_power_series< T > c;

	c.multiply(a, b);
	return(c);
}



template< class T >
inline base_dense_power_series< T >
operator* (const T & b,
	   const base_dense_power_series< T > & a)
{
	base_dense_power_series< T > c;

	c.multiply(b, a);
	return(c);
}



template< class T >
inline base_dense_power_series< T > &
operator *= (base_dense_power_series< T > & a,
	     const T & b)
{
	a.multiply(a, b);
	return(a);
}



template< class T >
inline base_dense_power_series< T >
operator/ (const base_dense_power_series< T > & a,
	   const T & b)
{
	base_dense_power_series< T > c;
	T x;

	invert(x, b);
	c.multiply(a, x);
	return(c);
}



template< class T >
inline base_dense_power_series< T > &
operator /= (base_dense_power_series< T > & a,
	     const T & b)
{
	T x;

	invert(x, b);
	a.multiply(a, x);
	return(a);
}



//
// ***** input / output *****
//

template< class T >
inline std::istream &
operator >> (std::istream & in, base_dense_power_series< T > & a)
{
	a.read(in);
	return(in);
}



template< class T >
inline std::ostream &
operator << (std::ostream & out, const base_dense_power_series< T > & a)
{
	a.write(out);
	return(out);
}



template< class T >
inline void
add (base_dense_power_series< T > & c,
     const base_dense_power_series< T > & a,
     const base_dense_power_series< T > & b)
{
	c.add(a, b);
}



template< class T >
inline void
add (base_dense_power_series< T > & c,
     const base_dense_power_series< T > & a,
     const 	                 T  & b)
{
	c.add(a, b);
}



template< class T >
inline void
add (base_dense_power_series< T > & c ,
     const                   T  & b ,
     const base_dense_power_series< T > & a)
{
	c.add(b, a);
}



template< class T >
inline void
subtract (base_dense_power_series< T > & c,
	  const base_dense_power_series< T > & a,
	  const base_dense_power_series< T > & b)
{
	c.subtract(a, b);
}



template< class T >
inline void
subtract (base_dense_power_series< T > & c,
	  const base_dense_power_series< T > & a,
	  const              T  & b)
{
	c.subtract(a, b);
}



template< class T >
inline void
subtract (base_dense_power_series< T > & c,
	  const              T  & b,
	  const base_dense_power_series< T > & a)
{
	c.subtract(b, a);
}



template< class T >
inline void
multiply (base_dense_power_series< T > & c,
	  const base_dense_power_series< T > & a,
	  const  	             T  & b)
{
	c.multiply(a, b);
}



template< class T >
inline void
multiply (base_dense_power_series< T > & c,
	  const                   T  & b,
	  const base_dense_power_series< T > & a)
{
	c.multiply(b, a);
}



template< class T >
inline void
divide (base_dense_power_series< T > & c,
	const base_dense_power_series< T > & a,
	const T  & b)
{
	c.divide(a, b);
}



template< class T >
inline void
negate (base_dense_power_series< T > & c, const base_dense_power_series< T > & a)
{
	c.negate(a);
}



//
// ***** miscellaneous *****
//

template< class T >
inline void
swap (base_dense_power_series< T > & c, base_dense_power_series< T > & a)
{
	c.swap(a);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/finite_fields/base_dense_power_series.cc"
#endif



#endif	// LIDIA_BASE_DENSE_POWER_SERIES_H_GUARD_
