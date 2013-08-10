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


#ifndef LIDIA_BASE_SPARSE_POWER_SERIES_H_GUARD_
#define LIDIA_BASE_SPARSE_POWER_SERIES_H_GUARD_


#ifndef LIDIA_SORT_VECTOR_H_GUARD_
# include	"LiDIA/sort_vector.h"
#endif
#ifndef LIDIA_COEFF_SPARSE_POWER_SERIES_H_GUARD_
# include	"LiDIA/finite_fields/coeff_sparse_power_series.h"
#endif



// only used in randomize()
//# include	"LiDIA/bigint.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



// **************************************************************************
// *
// *    class name :  base_sparse_power_series < T >
// *
// **************************************************************************

// fix the old incorrect spelling:
#define base_sparse_power_serie base_sparse_power_series


template< class T >
class base_sparse_power_series
{

protected :

	lidia_size_t first;
	lidia_size_t last;
	sort_vector< spc< T > > *coeff;
	bool sorted;

	bool input_verification;

	//
	//  *****  some protected member functions  *****
	//

	lidia_size_t alloc_size (lidia_size_t sz)
	{
		return sz;
	}

	void sort ()
	{
		coeff->sort();
		sorted = true;
	}

	void sort_test () const
	{
		if (sorted == false)
			lidia_error_handler ("base_sparse_power_series< T >",
					     "invalid structure: please, use sort()");
	}

	void init_test (char * mess) const;
	void rebuild ();


	//
	// ***** protected arithmetical member functions *****
	//

public:

	void add (const base_sparse_power_series< T > & a, const base_sparse_power_series< T > & b);
	void add (const base_sparse_power_series< T > & a, const T & b);
	void add (const T & b, const base_sparse_power_series< T > & a);

	void subtract (const base_sparse_power_series< T > & a, const base_sparse_power_series< T > & b);
	void subtract (const base_sparse_power_series< T > & a, const T & b);
	void subtract (const T & b, const base_sparse_power_series< T > & a);

	void multiply (const base_sparse_power_series< T > & a, const T & b);
	void multiply (const T & b, const base_sparse_power_series< T > & a);

	void divide (const base_sparse_power_series< T > & a, const T & b);
	void negate (const base_sparse_power_series< T > & a);

	void swap (base_sparse_power_series< T > & b);



public :

	//
	// c'tors and d'tor
	//

	base_sparse_power_series ();
	base_sparse_power_series (const T & elem , lidia_size_t l);
	base_sparse_power_series (const base_vector< T > & c, lidia_size_t f);
	base_sparse_power_series (const base_sparse_power_series< T > & a);

	~base_sparse_power_series ();



	//
	// assigners
	//

	void assign_zero (lidia_size_t l);
	void assign_one (lidia_size_t l);

	void assign (const base_sparse_power_series< T > & a);

	base_sparse_power_series< T > & operator = (const base_sparse_power_series< T > & a);



	//
	// comparators
	//

	bool is_one  () const;
	bool is_zero () const;
	bool is_equal(const base_sparse_power_series< T > & a) const;



	//
	//  *****  initialization functions  *****
	//

	void set (const T & elem , lidia_size_t l);
	void set (const T * c, const lidia_size_t * e, lidia_size_t sz, lidia_size_t l);
	// void set (const base_vector < T > & c, const base_vector < lidia_size_t > & e , lidia_size_t l);
	void set_coeff (const T & elem , lidia_size_t exp);
	void reduce_last (lidia_size_t l);
	void clear ();
	void normalize ();


	//
	//  *****  information functions  *****
	//

	lidia_size_t get_first () const;
	lidia_size_t get_last () const;

	void get (T * & c, lidia_size_t * & e, lidia_size_t & sz) const;
	// void get (base_vector < T > & c , base_vector < lidia_size_t > & e) const;
	void get_coeff (T & elem , lidia_size_t exp) const;



	//
	//  *****  subscription operator  *****
	//

	const T operator [] (lidia_size_t exp) const;



	//
	//  *****  I/O functions and operators  *****
	//

	void set_input_verification(bool v);

	int read (std::istream & in);
	int write (std::ostream & out) const;
	void info (std::ostream & out) const;



	//
	//  *****  generating random series  *****
	//

	// friend
	// void
	// randomize (base_sparse_power_series < T > & a , lidia_size_t f, lidia_size_t l , lidia_size_t all , const T & max);


	//
	//  *****  arithmetical operations  *****
	//

	void multiply_by_xn (lidia_size_t n);
	void compose (lidia_size_t n);

};



template< class T >
inline base_sparse_power_series< T > &
base_sparse_power_series< T >::operator = (const base_sparse_power_series< T > & a)
{
	this->assign(a);
	return *this;
}



template< class T >
inline bool
operator == (const base_sparse_power_series< T > & a,
	     const base_sparse_power_series< T > & b)
{
	return a.is_equal(b);
}



template< class T >
inline bool
operator != (const base_sparse_power_series< T > & a,
	     const base_sparse_power_series< T > & b)
{
	return !a.is_equal(b);
}



template< class T >
inline void
add (base_sparse_power_series< T > & res,
     const base_sparse_power_series< T > & a,
     const base_sparse_power_series< T > & b)
{
	res.add(a, b);
}



template< class T >
inline void
subtract (base_sparse_power_series< T > & res,
	  const base_sparse_power_series< T > & a,
	  const base_sparse_power_series< T > & b)
{
	res.subtract(a, b);
}



template< class T >
inline void
negate   (base_sparse_power_series< T > & res,
	  const base_sparse_power_series< T > & a)
{
	res.negate(a);
}



//
//  *****  scalar - operations  *****
//

template< class T >
inline void
add (base_sparse_power_series< T > & res,
     const base_sparse_power_series< T > & a,
     const T & b)
{
	res.add(a, b);
}



template< class T >
inline void
add (base_sparse_power_series< T > & res,
     const T & b,
     const base_sparse_power_series< T > & a)
{
	res.add(a, b);
}



template< class T >
inline void
subtract (base_sparse_power_series< T > & res,
	  const base_sparse_power_series< T > & a,
	  const T & b)
{
	res.subtract(a, b);
}



template< class T >
inline void
subtract (base_sparse_power_series< T > & res,
	  const T & b,
	  const base_sparse_power_series< T > & a)
{
	res.subtract(b, a);
}



template< class T >
inline void
multiply (base_sparse_power_series< T > & res,
	  const base_sparse_power_series< T > & a,
	  const T & b)
{
	res.multiply(a, b);
}



template< class T >
inline void
multiply (base_sparse_power_series< T > & res,
	  const T & b,
	  const base_sparse_power_series< T > & a)
{
	res.multiply(b, a);
}



template< class T >
inline void
divide (base_sparse_power_series< T > & res,
	const base_sparse_power_series< T > & a,
	const T & b)
{
	res.divide(a, b);
}


//
//  *****  arithmetical operators  *****
//

//  -----  addition  -----

template< class T >
inline base_sparse_power_series< T >
operator + (const base_sparse_power_series< T > & a,
	    const base_sparse_power_series< T > & b)
{
	base_sparse_power_series< T > res;

	res.add(a, b);
	return res;
}



template< class T >
inline base_sparse_power_series< T >
operator + (const base_sparse_power_series< T > & a,
	    const T & b)
{
	base_sparse_power_series< T > res;

	res.add(a, b);
	return res;
}



template< class T >
inline base_sparse_power_series< T >
operator + (const T & a,
	    const base_sparse_power_series< T > & b)
{
	base_sparse_power_series< T > res;

	res.add(a, b);
	return res;
}



template< class T >
inline base_sparse_power_series< T > &
operator += (base_sparse_power_series< T > & a,
	     const base_sparse_power_series< T > & b)
{
	a.add(a, b);
	return a;
}



template< class T >
inline base_sparse_power_series< T > &
operator += (base_sparse_power_series< T > & a,
	     const T & b)
{
	a.add(a, b);
	return a;
}



//  -----  subtraction  -----

template< class T >
inline base_sparse_power_series< T >
operator - (const base_sparse_power_series< T > & a)
{
	base_sparse_power_series< T > res;

	res.negate (a);
	return res;
}



template< class T >
inline base_sparse_power_series< T >
operator - (const base_sparse_power_series< T > & a,
	    const base_sparse_power_series< T > & b)
{
	base_sparse_power_series< T > res;

	res.subtract (a, b);
	return res;
}



template< class T >
inline base_sparse_power_series< T >
operator - (const base_sparse_power_series< T > & a,
	    const T & b)
{
	base_sparse_power_series< T > res;

	res.subtract(a, b);
	return res;
}



template< class T >
inline base_sparse_power_series< T >
operator - (const T & a,
	    const base_sparse_power_series< T > & b)
{
	base_sparse_power_series< T > res;

	res.subtract(a, b);
	return res;
}



template< class T >
inline base_sparse_power_series< T > &
operator -= (base_sparse_power_series< T > & a,
	     const base_sparse_power_series< T > & b)
{
	a.subtract(a, b);
	return a;
}



template< class T >
inline base_sparse_power_series< T > &
operator -= (base_sparse_power_series< T > & a,
	     const T & b)
{
	a.subtract(a, b);
	return a;
}



//  -----  multiplication  -----

template< class T >
inline base_sparse_power_series< T >
operator * (const base_sparse_power_series< T > & a,
	    const T & b)
{
	base_sparse_power_series< T > res;

	res.multiply(a, b);
	return res;
}



template< class T >
inline base_sparse_power_series< T >
operator * (const T & a,
	    const base_sparse_power_series< T > & b)
{
	base_sparse_power_series< T > res;

	res.multiply(a, b);
	return res;
}



template< class T >
inline base_sparse_power_series< T > &
operator *= (base_sparse_power_series< T > & a,
	     const T & b)
{
	a.multiply(a, b);
	return a;
}



//  -----  division  -----

template< class T >
inline base_sparse_power_series< T >
operator / (const base_sparse_power_series< T > & a,
	    const T & b)
{
	base_sparse_power_series< T > res;

	res.divide(a, b);
	return res;
}



template< class T >
inline base_sparse_power_series< T > &
operator /= (base_sparse_power_series< T > & a,
	     const T & b)
{
	a.divide(a, b);
	return a;
}



//
//  *****  I/O functions and operators  *****
//

template< class T >
inline std::istream &
operator >> (std::istream & in , base_sparse_power_series< T > & a)
{
	a.read(in);
	return (in);
}



template< class T >
inline std::ostream &
operator << (std::ostream & out , const base_sparse_power_series< T > & a)
{
	a.write(out);
	return (out);
}



//
//  ***** swap - function *****
//

template< class T >
inline void
swap (base_sparse_power_series< T > & a , base_sparse_power_series< T > & b)
{
	a.swap(b);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/base_sparse_power_series.cc"
#endif



#endif	// LIDIA_SPARSE_POWER_SERIES_H_GUARD_
