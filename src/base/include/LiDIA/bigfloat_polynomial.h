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
// template polynomial class for bigfloats
//


#ifndef LIDIA_BIGFLOAT_POLYNOMIAL_H_GUARD_
#define LIDIA_BIGFLOAT_POLYNOMIAL_H_GUARD_



#ifndef LIDIA_POLY_INTERN_H_GUARD_
# include	"LiDIA/base/poly_intern.h"
#endif
#ifndef LIDIA_FIELD_POLYNOMIAL_H_GUARD_
# include	"LiDIA/field_polynomial.h"
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



template<>
class polynomial< bigfloat > : public field_polynomial< bigfloat >
{
public:

	polynomial(): field_polynomial< bigfloat > ()
	{ }

	polynomial(bigfloat x): field_polynomial< bigfloat > (x)
	{ }

	polynomial(const bigfloat * v, lidia_size_t d):
		field_polynomial< bigfloat > (v, d)
	{ }

	polynomial(const base_vector< bigfloat > v): field_polynomial< bigfloat > (v)
	{ }

	polynomial(const base_polynomial< bigfloat > &p):
		field_polynomial< bigfloat > (p)
	{ }

	~polynomial()
	{ }

	polynomial< bigfloat > & operator = (const base_polynomial< bigfloat > &a)
	{
		base_polynomial< bigfloat >::assign(a);
		return *this;
	}

	//
	// Cast operator:
	//

	operator base_polynomial< bigcomplex > () const;


	friend polynomial< bigfloat > operator / (const base_polynomial< bigfloat > &a,
						  const base_polynomial< bigfloat > &b);

	friend polynomial< bigfloat > operator / (const base_polynomial< bigfloat > &a,
						  const bigfloat & b);

	friend polynomial< bigfloat > operator % (const base_polynomial< bigfloat > &a,
						  const base_polynomial< bigfloat > &b);

	friend polynomial< bigfloat > integral(const base_polynomial< bigfloat > & a);

	friend polynomial< bigfloat > gcd(const base_polynomial< bigfloat > &aa,
					  const base_polynomial< bigfloat > &bb);

	friend polynomial< bigfloat > xgcd(polynomial< bigfloat > &x,
					   polynomial< bigfloat > &y,
					   const base_polynomial< bigfloat > &aa,
					   const base_polynomial< bigfloat > &bb);

};



inline polynomial< bigfloat >
operator / (const base_polynomial< bigfloat > &a,
	    const base_polynomial< bigfloat > &b)
{
	polynomial< bigfloat > q, r;

	div_rem (q, r, a, b);
	return q;
}



inline polynomial< bigfloat >
operator / (const base_polynomial< bigfloat > &a,
	    const bigfloat & b)
{
	polynomial< bigfloat > c;

	divide (c, a, b);
	return c;
}



inline polynomial< bigfloat >
operator % (const base_polynomial< bigfloat > &a,
	    const base_polynomial< bigfloat > &b)
{
	polynomial< bigfloat > q, r;

	div_rem (q, r, a, b);
	return r;
}



inline polynomial< bigfloat >
integral(const base_polynomial< bigfloat > & a)
{
	polynomial< bigfloat > c;

	integral(c, a);
	return c;
}



inline polynomial< bigfloat >
gcd(const base_polynomial< bigfloat > &aa,
    const base_polynomial< bigfloat > &bb)
{
	polynomial< bigfloat > g;

	gcd(g, aa, bb);
	return g;
}



inline polynomial< bigfloat >
xgcd(polynomial< bigfloat > &x,
     polynomial< bigfloat > &y,
     const base_polynomial< bigfloat > &aa,
     const base_polynomial< bigfloat > &bb)
{
	polynomial< bigfloat > g;

	xgcd(g, x, y, aa, bb);
	return g;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_BIGFLOAT_POLYNOMIAL_H_GUARD_
