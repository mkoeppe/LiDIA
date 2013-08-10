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
// This  include  file  specializes the
// template polynomial class for bigints
//


#ifndef LIDIA_BIGINT_POLYNOMIAL_H_GUARD_
#define LIDIA_BIGINT_POLYNOMIAL_H_GUARD_



#ifndef LIDIA_POLY_INTERN_H_GUARD_
# include	"LiDIA/base/poly_intern.h"
#endif
#ifndef LIDIA_BIGRATIONAL_H_GUARD_
# include	"LiDIA/bigrational.h"
#endif
#ifndef LIDIA_BIGFLOAT_H_GUARD_
# include	"LiDIA/bigfloat.h"
#endif
#ifndef LIDIA_BIGCOMPLEX_H_GUARD_
# include	"LiDIA/bigcomplex.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class bigrational;
class bigfloat;
class bigcomplex;

//
// debug level
//
//   0 : div_rem
//   1 : power_mod
//   2 : gcd
//   3 : div_rem mod p
//   4 : power(_mod) mod p
//   5 : gcd mod p
//   6 : factoring mod p -- internals
//   7 : factoring mod p -- global
//   8 : number of real roots
//

template<>
class polynomial< bigint > : public base_polynomial< bigint >
{
	// For some routines, I would like to use Thomas Pfahler's
	// fp_polynomial, which is not yet finished...
public:

	//
	// constructors and destructor
	//

	polynomial(): base_polynomial< bigint > ()
	{ }

	polynomial(bigint x): base_polynomial< bigint > (x)
	{ }

	polynomial(const bigint * v, lidia_size_t d): base_polynomial< bigint > (v, d)
	{ }

	polynomial(const long * v, lidia_size_t d): base_polynomial< bigint > ()
	{
		if (d < 0 || v == NULL)
			lidia_error_handler_para(d, "d", "d >= 0",
						 PRT, "v", "v != NULL",
						 "base_polynomial< T >::"
						 "base_polynomial(const T * v, lidia_size_t d)",
						 "TEST", "JA AJ");
		deg = d;
		coeff = new bigint[d + 1];
		memory_handler(coeff, DM_BP, "base_polynomial(const T *, lidia_size_t)"
			       " :: Error in memory allocation (coeff)");

		for (lidia_size_t i = 0; i <= d; i++)
			coeff[i] = v[i];
		remove_leading_zeros();
	}

	polynomial(const base_vector< bigint > v): base_polynomial< bigint > (v)
	{ }


	polynomial(const base_polynomial< bigint > &p): base_polynomial< bigint > (p)
	{
		//std::cout << "p = " << p << std::endl;
		//std::cout << "IN: polynomial(const base_polynomial < bigint > &p)" << std::endl;
		//std::cout << "*this = " << *this << std::endl;
	}

	polynomial(const polynomial< bigint > &p): base_polynomial< bigint > (p)
	{
		//std::cout << "p = " << p << std::endl;
		//std::cout << "IN: polynomial(const polynomial < bigint > &p)" << std::endl;
		//std::cout << "*this = " << *this << std::endl;
	}

	~polynomial()
	{ }

	polynomial< bigint > &operator = (const base_polynomial< bigint > &a)
	{
		base_polynomial< bigint >::assign(a);
		return *this;
	}

	polynomial< bigint > &operator = (const polynomial< bigint > &a)
	{
		base_polynomial< bigint >::assign(a);
		return *this;
	}

	//
	// Cast operators:
	//

	operator base_polynomial< bigrational > () const;
	operator base_polynomial< bigfloat > () const;
	operator base_polynomial< bigcomplex > () const;


	//
	// Pseudo - division and related stuff
	//

	friend void div_rem(polynomial< bigint > &q, polynomial< bigint > &r,
			    const base_polynomial< bigint > &a,
			    const base_polynomial< bigint > &b);

	friend void divide(polynomial< bigint > & c,
			   const base_polynomial< bigint > & a, const bigint & b);

	friend void divide(polynomial< bigint > & q,
			   const base_polynomial< bigint > & a,
			   const base_polynomial< bigint > & b);

	friend void remainder(polynomial< bigint > & c,
			      const base_polynomial< bigint > & a, const bigint & b);

	friend void remainder(polynomial< bigint > & r,
			      const base_polynomial< bigint > & a,
			      const base_polynomial< bigint > & b);

	friend void power_mod(polynomial< bigint > & c,
			      const base_polynomial< bigint > & a, const bigint & b,
			      const base_polynomial< bigint > & f);

	friend polynomial< bigint > operator / (const base_polynomial< bigint > &a,
						const base_polynomial< bigint > &b);

	friend polynomial< bigint > operator / (const base_polynomial< bigint > &a,
						const bigint & b);

	friend polynomial< bigint > operator % (const base_polynomial< bigint > &a,
						const base_polynomial< bigint > &b);

	friend polynomial< bigint > operator % (const base_polynomial< bigint > &a,
						const bigint &b);

	polynomial< bigint > & operator /= (const base_polynomial< bigint > &a)
	{
		polynomial< bigint > r;

		div_rem(*this, r, *this, a);
		return *this;
	}

	polynomial< bigint > & operator /= (const bigint &a)
	{
		divide(*this, *this, a);
		return *this;
	}

	polynomial< bigint > & operator %= (const base_polynomial< bigint > &a)
	{
		polynomial< bigint > q;

		div_rem(q, *this, *this, a);
		return *this;
	}

	polynomial< bigint > & operator %= (const bigint &b)
	{
		remainder(*this, *this, b);
		return *this;
	}

	//
	// Gcd's and related stuff, i.e. content and primitive part.
	//

	friend bigint cont(const base_polynomial< bigint > &a);

	friend polynomial< bigint > pp(const base_polynomial< bigint > &a);

	friend polynomial< bigint > gcd(const base_polynomial< bigint > &aa,
					const base_polynomial< bigint > &bb);

	friend polynomial< bigint > xgcd(polynomial< bigint > &x,
					 polynomial< bigint > &y,
					 const base_polynomial< bigint > &aa,
					 const base_polynomial< bigint > &bb);

	//
	// Number of real roots
	//

	friend lidia_size_t no_of_real_roots(const base_polynomial< bigint > & poly_T);
};



inline void
divide (polynomial< bigint > & q,
	const base_polynomial< bigint > & a,
	const base_polynomial< bigint > & b)
{
	polynomial< bigint > r;

	div_rem(q, r, a, b);
}



inline void
remainder(polynomial< bigint > & r,
	  const base_polynomial< bigint > & a,
	  const base_polynomial< bigint > & b)
{
	polynomial< bigint > q;

	div_rem(q, r, a, b);
}



inline polynomial< bigint >
operator / (const base_polynomial< bigint > &a, const base_polynomial< bigint > &b)
{
	polynomial< bigint > q, r;

	div_rem(q, r, a, b);
	return q;
}



inline polynomial< bigint >
operator / (const base_polynomial< bigint > &a, const bigint & b)
{
	polynomial< bigint > q;

	divide(q, a, b);
	return q;
}



inline polynomial< bigint >
operator % (const base_polynomial< bigint > &a, const base_polynomial< bigint > &b)
{
	polynomial< bigint > q, r;

	div_rem(q, r, a, b);
	return r;
}



inline polynomial< bigint >
operator % (const base_polynomial< bigint > &a, const bigint &b)
{
	polynomial< bigint > r;

	remainder(r, a, b);
	return r;
}



//
// Resultant and Discriminant
//

bigint resultant(const polynomial< bigint > &aa,
		 const polynomial< bigint > &bb);
bigint discriminant(const polynomial< bigint > &a);



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_BIGINT_POLYNOMIAL_H_GUARD_
