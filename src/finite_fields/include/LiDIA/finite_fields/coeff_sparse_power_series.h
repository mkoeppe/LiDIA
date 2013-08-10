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


#ifndef LIDIA_COEFF_SPARSE_POWER_SERIES_H_GUARD_
#define LIDIA_COEFF_SPARSE_POWER_SERIES_H_GUARD_



#ifndef LIDIA_ARITH_INL_GUARD_
# include	"LiDIA/arith.inl"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



// **************************************************************************
// *
// *    class name :  spc < T >
// *
// **************************************************************************

template< class T >
struct spc
{
public:

        T     	     coeff;
        lidia_size_t exp;



public:

	//
	// c'tors and d'tor
	//

	spc ();
	spc (const spc< T > & coeff);
	spc (const T & c, lidia_size_t exp);
	~spc ();



	//
	// accessors
	//

	const T & get_coeff() const;
	lidia_size_t get_exp() const;



	//
	// assigners
	//

	void assign(const spc< T > & c);
	void assign(const T & c, lidia_size_t e);

	spc< T > & operator = (const spc< T > & c);



	//
	// comparators
	//

	bool operator < (const spc < T > & c) const;
	bool operator <= (const spc< T > & c) const;
	bool operator == (const spc< T > & c) const;

	int coeff_cmp_zero (const spc< T > & b) const;


	//
	// I/O
	//

	void read  (std::istream & in);
	void write (std::ostream & out) const;


	void swap (spc< T > & b);
};



//
// c'tors and d'tor
//

template< class T >
inline
spc< T >::spc ()
{
	// nothing to do
}


template< class T >
inline
spc< T >::spc (const spc< T > & c)
	: coeff(c.coeff),
	  exp(c.exp)
{
	// nothing to do
}



template< class T >
inline
spc< T >::spc (const T & c , lidia_size_t e)
	: coeff(c),
	  exp(e)
{
	// nothing to do
}



template< class T >
inline
spc< T >::~spc ()
{
	// nothing to do
}



//
// accessors
//

template< class T >
inline const T &
spc< T >::get_coeff () const
{
	return this->coeff;
}



template< class T >
inline lidia_size_t
spc< T >::get_exp () const
{
	return this->exp;
}



//
// assigners
//

template< class T >
inline void
spc< T >::assign (const spc< T > & c)
{
	if (&c != this) {
		this->coeff = c.coeff;
		this->exp = c.exp;
	}
}



template< class T >
inline void
spc< T >::assign (const T & c, lidia_size_t e)
{
	this->coeff = c;
	this->exp = e;
}



template< class T >
inline spc< T > &
spc< T >::operator = (const spc< T > & c)
{
	debug_handler ("sparse_powser::spc" , "operator = ");

	assign(c);
	return *this;
}



//
// comparators
//

template< class T >
inline bool
spc< T >::operator < (const spc < T > & c) const
{
	debug_handler ("sparse_powser::spc" , "operator < ");

	return (this->exp < c.exp);
}



template< class T >
inline bool
spc< T >::operator <= (const spc< T > & c) const
{
	debug_handler ("sparse_powser::spc< T >" , "operator <= ");

	return (this->exp <= c.exp);
}



template< class T >
inline bool
spc< T >::operator == (const spc< T > & c) const
{
	debug_handler ("sparse_powser::spc" , "operator == ");

	return (this->exp == c.exp);
}



template< class T >
inline void
spc< T >::swap (spc< T > & b)
{
	LiDIA::swap(this->coeff, b.coeff);
	LiDIA::swap(this->exp, b.exp);
}



template< class T >
inline std::istream &
operator >> (std::istream & in, spc< T > & c)
{
	c.read(in);
	return in;
}



template< class T >
inline std::ostream &
operator << (std::ostream & out, const spc< T > & c)
{
	c.write(out);
	return out;
}



template< class T >
inline void
swap (spc< T > & a, spc< T > & b)
{
	a.swap(b);
}



template< class T >
inline int
coeff_cmp_zero (const spc< T > & a, const spc< T > & b)
{
	return (a.coeff_cmp_zero(b));
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/coeff_sparse_power_series.cc"
#endif



#endif	// LIDIA_COEFF_SPARSE_POWER_SERIES_H_GUARD_
