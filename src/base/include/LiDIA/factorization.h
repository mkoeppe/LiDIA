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


#ifndef LIDIA_FACTORIZATION_H_GUARD_
#define LIDIA_FACTORIZATION_H_GUARD_



#ifndef LIDIA_UDIGIT_H_GUARD_
# include	"LiDIA/udigit.h"
#endif
#ifndef LIDIA_SORT_VECTOR_H_GUARD_
# include	"LiDIA/sort_vector.h"
#endif
#ifndef LIDIA_PPAIR_H_GUARD_
# include	"LiDIA/base/ppair.h"
#endif
#ifndef LIDIA_SINGLE_FACTOR_H_GUARD_
# include	"LiDIA/single_factor.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



//template < class T > class single_factor;
//class decomposable_object;

class factorization_flags
{
protected:
	enum factorization_state {
		sorted = 1,
		normalized = 2,
		refined = 4
	};
};



template< class T >
class factorization: public factorization_flags, private decomposable_object
{

private:

	// data elements

	single_factor< T > *epsilon; // unit
	sort_vector< ppair< single_factor < T >, lidia_size_t > >
	prime_component, composite_component;


	udigit attributes; // information not concerning primality


	// access functions for 'attributes'
	void clear_attributes()
	{
		attributes &= ~ (sorted | normalized | refined);
	}

	void check_for_trivial_case();
	// sets information for simple factorizations


	static void
	refine2(sort_vector< ppair< single_factor< T >, lidia_size_t > > & v,
		ppair< single_factor< T >, lidia_size_t > & SF,
		const ppair< single_factor< T >, lidia_size_t > & a,
		const ppair< single_factor< T >, lidia_size_t > & b);

	static void
	normalize(sort_vector< ppair< single_factor< T >, lidia_size_t > > &v);

public:
	//
	// constructors, destructor
	//

	factorization();
	factorization(const factorization< T > & F);
	factorization(const single_factor< T > & f);
	~factorization();

	void reset(); // reset the factorization to one.
	void kill(); // reset the factorization and deallocate memory.

	//
	// assignment
	//

	factorization< T > & operator = (const factorization< T > & x);
	factorization< T > & operator = (const single_factor< T > & x);

#ifndef HEADBANGER
	void assign (const factorization< T > & F);
	void assign (const single_factor< T > & f);
#endif		//HEADBANGER

	//
	// comparisons
	//

	bool operator == (const factorization< T > & b) const;
	bool operator != (const factorization< T > & b)
	{
		return !(*this == b);
	}


	//
	// access functions
	//

	const single_factor< T > & prime_base(lidia_size_t) const;
	lidia_size_t prime_exponent(lidia_size_t) const;
	const single_factor< T > & composite_base(lidia_size_t) const;
	lidia_size_t composite_exponent(lidia_size_t) const;

	const T & unit() const;

	T value() const;

	lidia_size_t no_of_prime_components() const
	{
		return prime_component.size();
	}

	lidia_size_t no_of_composite_components() const
	{
		return composite_component.size();
	}

	lidia_size_t no_of_components() const
	{
		return prime_component.size() + composite_component.size();
	}


	//
	// predicates
	//

	bool is_prime_factorization();

	bool is_prime_factorization(int test);
	//  test == 0 ->no explicit primality test (only flags are checked)
	//  returns true iff (*this) is known to be a prime factorization
	//  test != 0 ->explicit primality test for unknown elements
	//  returns false iff (*this) has an element which is not_prime
	//  (see manual)

	bool is_sorted() const; // 'sorted' means: both lists are sorted

	bool is_normalized() const;

	bool is_refined() const;

	//
	// input / output
	//

	//
	// I/O - format :
	// [ [unit] ]
	// or
	// [ [unit], [base_1, exp_1], [base_2, exp_2], ..., [base_i, exp_i] ]
	// or
	// [ [base_1, exp_1], [base_2, exp_2], ..., [base_i, exp_i] ]
	//

	void read(std::istream &is);
	void write(std::ostream &os) const;


	//#####
	// NUR ZU TESTZWECKEN !!!
	void output() const;
	// **********************


	//
	//  functions for computing with factorizations
	//

#ifndef HEADBANGER
	void concat(const factorization< T > &a, const factorization< T > &b);
	// *this = a concat b
#endif		//HEADBANGER

	void invert(); //== power(-1)

	void square(); //== power(2)

	void power(lidia_size_t p);



	//
	//  high-level functions
	//

	void sort();

	void normalize(); //was: compose();

	void refine(); //factor refinement algorithm
	bool refine(const single_factor< T > & x);
	//returns true iff *this can be refined by computing gcd's
	//of all (composite) elements with 'x'

	void replace(lidia_size_t ix, const factorization< T > & F);
	// deletes composite_component[i] and appends the components of F

	void append(const single_factor< T > & a)
	{
		append(a, 1);
	}
	void append(const single_factor< T > & a, lidia_size_t exp);
	// appends 'a^exp' to the factorization

	void swapping_append(single_factor< T > & a)
	{
		swapping_append(a, 1);
	}
	void swapping_append(single_factor< T > & a, lidia_size_t exp);


	void append(const T& a, lidia_size_t exp);

	//
	// factoring
	//

	void factor_all_components();
	// calls 'single_factor < T >::factor()' for all composite components
	// and replaces the components by the resulting factorizations
};



template< class T >
inline
factorization< T >::~factorization ()
{
	delete epsilon;
}



//
// assigners
//

template< class T >
inline factorization< T > &
factorization< T >::operator = (const factorization< T > & f)
{
	assign(f);
	return *this;
}



template< class T >
inline factorization< T > &
factorization< T >::operator = (const single_factor< T > & f)
{
	assign(f);
	return *this;
}



//
// accessors
//

template< class T >
inline const single_factor< T > &
factorization< T >::prime_base (lidia_size_t index) const
{
	if ((index< 0) || (index >= no_of_prime_components()))
		lidia_error_handler("factorization< T >",
				    "prime_base::index out of range");

	return prime_component[index].left();
}



template< class T >
inline lidia_size_t
factorization< T >::prime_exponent (lidia_size_t index) const
{
	if ((index< 0) || (index >= no_of_prime_components()))
		lidia_error_handler("factorization< T >",
				    "prime_exponent::index out of range");

	return prime_component[index].right();
}



template< class T >
inline const single_factor< T > &
factorization< T >::composite_base (lidia_size_t index) const
{
	if ((index< 0) || (index >= no_of_composite_components()))
		lidia_error_handler("factorization< T >",
				    "composite_base::index out of range");

	return composite_component[index].left();
}



template< class T >
inline lidia_size_t
factorization< T >::composite_exponent (lidia_size_t index) const
{
	if ((index< 0) || (index >= no_of_composite_components()))
		lidia_error_handler("factorization< T >",
				    "composite_exponent::index out of range");

	return composite_component[index].right();
}



template< class T >
inline const T &
factorization< T >::unit () const
{
	return epsilon->base();
}



//
// predicates
//

template< class T >
inline bool
factorization< T >::is_prime_factorization ()
{
	return this->is_prime_factorization(0);
}



template< class T >
inline bool
factorization< T >::is_sorted () const
{
	return static_cast<bool>(attributes & sorted);
}



template< class T >
inline bool
factorization< T >::is_normalized () const
{
	return static_cast<bool>(attributes & normalized);
}



template< class T >
inline bool
factorization< T >::is_refined () const
{
	return static_cast<bool>(attributes & refined);
}



//
//
//

template< class T >
inline std::istream &
operator >> (std::istream &is, factorization< T > &F)
{
	F.read(is);
	return is;
}



template< class T >
inline std::ostream &
operator << (std::ostream &os, const factorization< T > & F)
{
	F.write(os);
	return os;
}



template< class T >
inline void
invert (factorization< T > & f,
	const factorization< T > & a)
{
	f.assign(a);
	f.invert();
}



template< class T >
inline factorization< T >
inverse (const factorization< T > &f)
{
	factorization< T > tmp(f);

	tmp.invert();
	return tmp;
}



template< class T >
inline void
square (factorization< T > & f,
	const factorization< T > & a)
{
	f.assign(a);
	f.square();
}



template< class T >
inline void
power (factorization< T > & f,
       const factorization< T > & a,
       lidia_size_t p)
{
	f.assign(a);
	f.power(p);
}



#ifndef HEADBANGER
template< class T >
inline void
multiply (factorization< T > &c,
	  const factorization< T > &a,
	  const factorization< T > &b)
{
	c.concat(a, b);
}



template< class T >
inline void
divide (factorization< T > &c,
	const factorization< T > &a,
	const factorization< T > &b)
{
	factorization< T > tmp(b);

	tmp.invert();
	multiply(c, a, tmp);
}
#endif	// HEADBANGER



template< class T >
inline factorization< T >
operator * (const factorization< T > & a,
	    const factorization< T > & b)
{
	factorization< T > c;

	multiply(c, a, b);
	return c;
}



template< class T >
inline factorization< T >
operator / (const factorization< T > & a,
	    const factorization< T > & b)
{
	factorization< T > c;

	divide(c, a, b);
	return c;
}



template< class T >
inline factorization< T >
factor_all_components (const factorization< T > &F)
{
	factorization< T > FF(F);

	FF.factor_all_components();
	return FF;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/factorization.cc"
#endif



#endif
