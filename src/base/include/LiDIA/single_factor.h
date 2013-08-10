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
//	Author	:
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_SINGLE_FACTOR_H_GUARD_
#define LIDIA_SINGLE_FACTOR_H_GUARD_



#ifndef LIDIA_LIDIA_H_GUARD_
# include	"LiDIA/LiDIA.h"
#endif
#ifndef LIDIA_ARITH_INL_GUARD_
# include	"LiDIA/arith.inl"
#endif
#ifndef LIDIA_BASE_FACTOR_H_GUARD_
# include	"LiDIA/base/base_factor.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class T > class single_factor;
template< class T > class base_factor;
template< class T > class factorization;



//****************************************************************
//                 class single_factor< T >
//****************************************************************

// forward declarations required for template friends
template< class T >
lidia_size_t ord_divide(const single_factor< T > &a,
			const single_factor< T > &b);
template< class T >
void gcd(single_factor< T > & c, const single_factor< T > & a,
	 const single_factor< T > & b);


template< class T >
class single_factor : public base_factor< T >
{
	//****************************************************************************
	//
	// This is the default implementation for the class single_factor.
	// This code is used if there is no special instantiation for type T.
	//
	// There exist special instantiations for type T = {bigint, Fp_polynomial, ideal}
	// single_factor< bigint >
	// single_factor< Fp_polynomial >
	// single_factor< ideal >
	//
	// If you want to add another special instantiation (not the default one),
	// compare example below (Fp_polynomial, gf_polynomial, bigint)
	//
	//****************************************************************************


	friend class factorization< T >;

private:
	//additional information should be defined here



public:
	//
	// constructors, destructor
	//
	single_factor();
	//default value must be '1' (the neutral element of multiplication)
	single_factor(const single_factor< T > &);
	//copy constructor
	single_factor(const T&);
	//conversion for type T
	~single_factor();



	//
	// swap
	//

public :

	void swap(single_factor< T > & a);


	//
	// assignment
	//
	single_factor< T > & operator = (const single_factor< T > & x);
	single_factor< T > & operator = (const T & x);
#ifndef HEADBANGER
	void assign(const single_factor< T > & x);
	void assign(const T & x);
#endif	// HEADBANGER

	//
	// predicates
	//

	bool is_one() const;
	bool is_prime_factor() const;
	bool is_prime_factor(int test);
	//test == 0 ->no explicit primality test (only a flag is checked)
	//test != 0 ->explicit prime-test if prime_state() == unknown
	//see manual

	T extract_unit();
	friend lidia_size_t ord_divide< T >(const single_factor< T > &a,
					    const single_factor< T > &b);


	friend void gcd< T >(single_factor< T > & c, const single_factor< T > & a,
			const single_factor< T > & b);


	//
	// factorization algorithms
	//

	factorization< T > factor() const;
	// standard factorization algorithm, used in function
	// "factorization < T >::factor_all_components"

};


//
// c'tors and d'tor
//

template< class T >
inline
single_factor< T >::single_factor()
	: base_factor< T > ()
{
	//DEFAULT VALUE MUST BE '1', i.e. the neutral element of multiplication !!!
	this->rep = 1; //rep.assign_one();
}



template< class T >
inline
single_factor< T >::single_factor(const single_factor< T > & x)
	: base_factor< T > (x)
{
	// nothing to do
}



template< class T >
inline
single_factor< T >::single_factor(const T & x)
	: base_factor< T > (x)
{
	// nothing to do
}



template< class T >
inline
single_factor< T >::~single_factor()
{
	// nothing to do
}



//
// assigners
//

template< class T >
inline void
single_factor< T >::assign (const single_factor< T > & x)
{
	if (&x != this) {
		this->rep = x.rep;
		this->set_prime_flag(x.prime_flag());
	}
}



template< class T >
inline void
single_factor< T >::assign (const T & x)
{
	this->rep = x;
	this->set_prime_flag(decomposable_object::unknown);
}



template< class T >
inline single_factor< T > &
single_factor< T >::operator = (const single_factor< T > & x)
{
	assign(x);
	return *this;
}



template< class T >
inline single_factor< T > &
single_factor< T >::operator = (const T & x)
{
	assign(x);
	return *this;
}



template< class T >
inline T
single_factor< T >::extract_unit ()
{
	//DEFAULT VARIANT: always returns '1' and leaves rep unchanged
	//single_factor< T > one;
	//return one.rep;
	return T(1);
}



//
// predicated
//

template< class T >
inline bool
single_factor< T >::is_one () const
{
	//single_factor< T > ONE; // default value is '1'
	//return (*this == ONE);
	return (this->rep == T(1));
}



template< class T >
inline bool
single_factor< T >::is_prime_factor() const
{
	return (this->prime_flag() == decomposable_object::prime);
}



template< class T >
inline bool
single_factor< T >::is_prime_factor (int test)
{
	if (this->is_prime_factor())
		return true;

	if (this->test == 0)		// => no explicit primality test
		return false;


	// NO DEFAULT PRIMALITY TEST

	lidia_error_handler("factorization< T >", "is_prime_factor(int)::explicit test not implemented");
	return false;
}



template< class T >
inline factorization< T >
factor (const single_factor< T > & f)
{
	return f.factor();
}



template< class T >
inline void
single_factor< T >::swap (single_factor< T > & a)
{
	base_factor< T >::swap(a);
}



template< class T >
inline void
swap (single_factor< T > & a, single_factor< T > & b)
{
	a.swap(b);
}



template< class T >
inline void
gcd (single_factor< T > & c,
     const single_factor< T > & a,
     const single_factor< T > & b)
	//c = gcd(a, b)
{
	c.rep = gcd(a.rep, b.rep);
	c.set_prime_flag(single_factor< T >::unknown);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#define LIDIA_CLASS_SINGLE_FACTOR

#include	"LiDIA/specialization/single_factor.special"



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/single_factor.cc"
#endif



#endif	// LIDIA_SINGLE_FACTOR_H_GUARD_
