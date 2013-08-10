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
// template polynomial class for bigcomplexes
//


#ifndef LIDIA_BIGCOMPLEX_POLYNOMIAL_H_GUARD_
#define LIDIA_BIGCOMPLEX_POLYNOMIAL_H_GUARD_



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



//
// debug level
//
//   0 : remove_leading_zeros -- redefine.
//   1 : integral
//   2 : root
//   3 : cohen
//   4 : roots
//


template<>
class polynomial< bigcomplex > : public field_polynomial< bigcomplex >
{
public:

	polynomial(): field_polynomial< bigcomplex > ()
	{ }

	polynomial(bigcomplex x): field_polynomial< bigcomplex > (x)
	{ }

	polynomial(const bigcomplex * v, lidia_size_t d):
		field_polynomial< bigcomplex > (v, d)
	{ }

	polynomial(const base_vector< bigcomplex > v):field_polynomial< bigcomplex > (v)
	{ }

	polynomial(const base_polynomial< bigcomplex > &p):
		field_polynomial< bigcomplex > (p)
	{ }

	~polynomial()
	{ }

	polynomial< bigcomplex > & operator = (const base_polynomial< bigcomplex > &a)
	{
		base_polynomial< bigcomplex >::assign(a);
		return *this;
	}

	friend bigcomplex root(const base_polynomial< bigcomplex > &,
			       const bigcomplex &);
	friend void cohen(const base_polynomial< bigcomplex > &,
			  bigcomplex *, int, int &);
	friend void roots(const base_polynomial< bigcomplex > &,
			  bigcomplex *);
	friend bigcomplex integral(const bigfloat &, const bigfloat &,
				   const base_polynomial< bigcomplex > &);


	friend polynomial< bigcomplex > operator / (const base_polynomial< bigcomplex > &a,
						    const base_polynomial< bigcomplex > &b);

	friend polynomial< bigcomplex > operator / (const base_polynomial< bigcomplex > &a,
						    const bigcomplex & b);

	friend polynomial< bigcomplex > operator % (const base_polynomial< bigcomplex > &a,
						    const base_polynomial< bigcomplex > &b);

	friend polynomial< bigcomplex > integral(const base_polynomial< bigcomplex > & a);

	friend polynomial< bigcomplex > gcd(const base_polynomial< bigcomplex > &aa,
					    const base_polynomial< bigcomplex > &bb);

	friend polynomial< bigcomplex > xgcd(polynomial< bigcomplex > &x,
					     polynomial< bigcomplex > &y,
					     const base_polynomial< bigcomplex > &aa,
					     const base_polynomial< bigcomplex > &bb);

};

bigcomplex root(const base_polynomial< bigcomplex > &,
		const bigcomplex &);
void cohen(const base_polynomial< bigcomplex > &,
	   bigcomplex *, int, int &);
void roots(const base_polynomial< bigcomplex > &,
	   bigcomplex *);
bigcomplex integral(const bigfloat &, const bigfloat &,
		    const base_polynomial< bigcomplex > &);
polynomial< bigcomplex > operator / (const base_polynomial< bigcomplex > &a,
				     const base_polynomial< bigcomplex > &b);
polynomial< bigcomplex > operator / (const base_polynomial< bigcomplex > &a,
				     const bigcomplex & b);
polynomial< bigcomplex > operator % (const base_polynomial< bigcomplex > &a,
				     const base_polynomial< bigcomplex > &b);
polynomial< bigcomplex > integral(const base_polynomial< bigcomplex > & a);
polynomial< bigcomplex > gcd(const base_polynomial< bigcomplex > &aa,
			     const base_polynomial< bigcomplex > &bb);
polynomial< bigcomplex > xgcd(polynomial< bigcomplex > &x,
			      polynomial< bigcomplex > &y,
			      const base_polynomial< bigcomplex > &aa,
			      const base_polynomial< bigcomplex > &bb);



inline polynomial< bigcomplex >
operator / (const base_polynomial< bigcomplex > &a,
	    const base_polynomial< bigcomplex > &b)
{
	polynomial< bigcomplex > q, r;

	div_rem (q, r, a, b);
	return q;
}



inline polynomial< bigcomplex >
operator / (const base_polynomial< bigcomplex > &a,
	    const bigcomplex & b)
{
	polynomial< bigcomplex > c;

	divide (c, a, b);
	return c;
}



inline polynomial< bigcomplex >
operator % (const base_polynomial< bigcomplex > &a,
	    const base_polynomial< bigcomplex > &b)
{
	polynomial< bigcomplex > q, r;

	div_rem (q, r, a, b);
	return r;
}



inline polynomial< bigcomplex >
integral(const base_polynomial< bigcomplex > & a)
{
	polynomial< bigcomplex > c;

	integral(c, a);
	return c;
}



inline polynomial< bigcomplex >
gcd(const base_polynomial< bigcomplex > &aa,
    const base_polynomial< bigcomplex > &bb)
{
	polynomial< bigcomplex > g;

	gcd(g, aa, bb);
	return g;
}



inline polynomial< bigcomplex >
xgcd(polynomial< bigcomplex > &x,
     polynomial< bigcomplex > &y,
     const base_polynomial< bigcomplex > &aa,
     const base_polynomial< bigcomplex > &bb)
{
	polynomial< bigcomplex > g;

	xgcd(g, x, y, aa, bb);
	return g;
}



// Jenkins-Traub root finding
base_vector< bigcomplex > roots(const base_polynomial< bigfloat > &);
base_vector< bigcomplex > roots(const base_polynomial< bigcomplex > &);



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_BIGCOMPLEX_POLYNOMIAL_H_GUARD_
