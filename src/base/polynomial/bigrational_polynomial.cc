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


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigrational.h"
#include	"LiDIA/polynomial.h"

#ifndef LIDIA_INCLUDE_CC
# include	"LiDIA/base/poly_intern.cc"
# include	"LiDIA/field_polynomial.cc"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



// instantiate class base_polynomial
template class base_polynomial< bigrational >;

// instantiate friend functions -- arithmetical functions
template void negate(base_polynomial< bigrational > & c,
		     const base_polynomial< bigrational > &a);

template void add(base_polynomial< bigrational > & c,
		  const base_polynomial< bigrational > & a,
		  const base_polynomial< bigrational > & b);
template void add(base_polynomial< bigrational > & c,
		  const base_polynomial< bigrational > & a,
		  const bigrational & b);
template void add(base_polynomial< bigrational > & c,
		  const bigrational & b,
		  const base_polynomial< bigrational > & a);

template void subtract(base_polynomial< bigrational > & c,
		       const base_polynomial< bigrational > & a,
		       const base_polynomial< bigrational > & b);
template void subtract(base_polynomial< bigrational > & c,
		       const base_polynomial< bigrational > & a,
		       const bigrational & b);
template void subtract(base_polynomial< bigrational > & c,
		       const bigrational & b,
		       const base_polynomial< bigrational > & a);

template void multiply(base_polynomial< bigrational > & c,
		       const base_polynomial< bigrational > & a,
		       const base_polynomial< bigrational > & b);
template void multiply(base_polynomial< bigrational > & c,
		       const base_polynomial< bigrational > & a,
		       const bigrational & b);
template void multiply(base_polynomial< bigrational > & c,
		       const bigrational & b,
		       const base_polynomial< bigrational > & a);

template void power(base_polynomial< bigrational > & c,
		    const base_polynomial< bigrational > & a, const bigint & b);

template void derivative(base_polynomial< bigrational > &c,
			 const base_polynomial< bigrational > &a);
template base_polynomial< bigrational >
derivative(const base_polynomial< bigrational > &a);

// instantiate friend functions -- operators
template bool operator == (const base_polynomial< bigrational > &a,
			   const base_polynomial< bigrational > &b);
template bool operator != (const base_polynomial< bigrational > &a,
			   const base_polynomial< bigrational > &b);
template base_polynomial< bigrational >
operator -(const base_polynomial< bigrational > &a);

template base_polynomial< bigrational >
operator +(const base_polynomial< bigrational > &a,
	   const base_polynomial< bigrational > &b);

template base_polynomial< bigrational >
operator +(const base_polynomial< bigrational > &a,
	   const bigrational & b);

template base_polynomial< bigrational >
operator +(const bigrational & b,
	   const base_polynomial< bigrational > &a);

template base_polynomial< bigrational >
operator -(const base_polynomial< bigrational > &a,
	   const base_polynomial< bigrational > &b);

template base_polynomial< bigrational >
operator -(const base_polynomial< bigrational > &a,
	   const bigrational &b);

template base_polynomial< bigrational >
operator -(const bigrational &a,
	   const base_polynomial< bigrational > &b);

template base_polynomial< bigrational >
operator *(const base_polynomial< bigrational > &a,
	   const base_polynomial< bigrational > &b);

template base_polynomial< bigrational >
operator *(const base_polynomial< bigrational > &a,
	   const bigrational &b);

template base_polynomial< bigrational >
operator *(const bigrational &b,
	   const base_polynomial< bigrational > &a);

template std::istream & operator >> (std::istream &,
				     base_polynomial< bigrational > &);
template std::ostream & operator << (std::ostream &,
				     const base_polynomial< bigrational > &);

// instantiate friend functions -- access functions
template bigrational lead_coeff(const base_polynomial< bigrational > &a);
template bigrational const_term(const base_polynomial< bigrational > &a);

// instantiate class field_polynomial
template class field_polynomial< bigrational >;

// instantiate friend functions -- arithmetical functions
template void div_rem(field_polynomial< bigrational > &q,
		      field_polynomial< bigrational > &r,
		      const base_polynomial< bigrational > &a,
		      const base_polynomial< bigrational > &b);
template void divide(field_polynomial< bigrational > & c,
		     const base_polynomial< bigrational > & a,
		     const bigrational & b);
template void divide(field_polynomial< bigrational > & q,
		     const base_polynomial< bigrational > & a,
		     const base_polynomial< bigrational > & b);
template void remainder(field_polynomial< bigrational > & r,
			const base_polynomial< bigrational > & a,
			const base_polynomial< bigrational > & b);
template void power_mod(field_polynomial< bigrational > & c,
			const base_polynomial< bigrational > & a,
			const bigint & b,
			const base_polynomial< bigrational > & f);

template void integral(field_polynomial< bigrational > &c,
		       const base_polynomial< bigrational > &a);

template void gcd(field_polynomial< bigrational > &d,
		  const base_polynomial< bigrational > &aa,
		  const base_polynomial< bigrational > &bb);

template void xgcd(field_polynomial< bigrational > &d,
		   field_polynomial< bigrational > &x,
		   field_polynomial< bigrational > &y,
		   const base_polynomial< bigrational > &aa,
		   const base_polynomial< bigrational > &bb);



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
