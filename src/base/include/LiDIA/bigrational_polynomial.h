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
// template polynomial class for bigrationals
//


#ifndef LIDIA_BIGRATIONAL_POLYNOMIAL_H_GUARD_
#define LIDIA_BIGRATIONAL_POLYNOMIAL_H_GUARD_


#ifndef LIDIA_POLY_INTERN_H_GUARD_
# include	"LiDIA/base/poly_intern.h"
#endif
#ifndef LIDIA_FIELD_POLYNOMIAL_H_GUARD_
# include	"LiDIA/field_polynomial.h"
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



template<>
class polynomial< bigrational > : public field_polynomial< bigrational >
{
public:

	polynomial(): field_polynomial< bigrational > ()
	{ }

	polynomial(bigrational x): field_polynomial< bigrational > (x)
	{ }

	polynomial(const bigrational * v, lidia_size_t d):
		field_polynomial< bigrational > (v, d)
	{ }

	polynomial(const base_vector< bigrational > v):
		field_polynomial< bigrational > (v)
	{ }

	polynomial(const base_polynomial< bigrational > &p):
		field_polynomial< bigrational > (p)
	{ }

	~polynomial()
	{ }

	polynomial< bigrational > & operator = (const base_polynomial< bigrational >
						&a)
	{
		base_polynomial< bigrational >::assign(a);
		return *this;
	}


	//
	// Cast operators:
	//


	//
	// Cast operators:
	//

	operator base_polynomial< bigfloat > () const;
	operator base_polynomial< bigcomplex > () const;



	friend polynomial< bigrational > operator / (const base_polynomial< bigrational > &a,
						     const base_polynomial< bigrational > &b);

	friend polynomial< bigrational > operator / (const base_polynomial< bigrational > &a,
						     const bigrational & b);

	friend polynomial< bigrational > operator % (const base_polynomial< bigrational > &a,
						     const base_polynomial< bigrational > &b);

	friend polynomial< bigrational > integral(const base_polynomial< bigrational > & a);

	friend polynomial< bigrational > gcd(const base_polynomial< bigrational > &aa,
					     const base_polynomial< bigrational > &bb);

	friend polynomial< bigrational > xgcd(polynomial< bigrational > &x,
					      polynomial< bigrational > &y,
					      const base_polynomial< bigrational > &aa,
					      const base_polynomial< bigrational > &bb);

};



inline polynomial< bigrational >
operator / (const base_polynomial< bigrational > &a,
	    const base_polynomial< bigrational > &b)
{
	polynomial< bigrational > q, r;

	div_rem (q, r, a, b);
	return q;
}



inline polynomial< bigrational >
operator / (const base_polynomial< bigrational > &a,
	    const bigrational & b)
{
	polynomial< bigrational > c;

	divide (c, a, b);
	return c;
}



inline polynomial< bigrational >
operator % (const base_polynomial< bigrational > &a,
	    const base_polynomial< bigrational > &b)
{
	polynomial< bigrational > q, r;

	div_rem (q, r, a, b);
	return r;
}



inline polynomial< bigrational >
integral(const base_polynomial< bigrational > & a)
{
	polynomial< bigrational > c;

	integral(c, a);
	return c;
}



inline polynomial< bigrational >
gcd(const base_polynomial< bigrational > &aa,
    const base_polynomial< bigrational > &bb)
{
	polynomial< bigrational > g;

	gcd(g, aa, bb);
	return g;
}



inline polynomial< bigrational >
xgcd(polynomial< bigrational > &x,
     polynomial< bigrational > &y,
     const base_polynomial< bigrational > &aa,
     const base_polynomial< bigrational > &bb)
{
	polynomial< bigrational > g;

	xgcd(g, x, y, aa, bb);
	return g;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_BIGRATIONAL_POLYNOMIAL_H_GUARD_
