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
#include	"LiDIA/bigfloat.h"
#include	"LiDIA/polynomial.h"

#ifndef LIDIA_INCLUDE_CC
# include	"LiDIA/base/poly_intern.cc"
# include	"LiDIA/field_polynomial.cc"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



#define DV_FLP LDBL_UNIPOL+20
#define DM_FLP "polynomial< bigfloat >"
#define LP_ERROR poly_error_msg


template <>
void base_polynomial< bigfloat >::remove_leading_zeros()
{
	debug_handler_c(DM_FLP, "in member - function remove_leading_zeros ()",
			DV_FLP, std::cout << "original degree is " << deg << std::endl << std::flush);
	bigfloat c, *np;
	lidia_size_t d, i;

	d = deg;
	np = coeff + d;
	while (d >= 0 && np->is_approx_zero())
		d--, np--;

	if (d < 0) {
		deg = d;
		delete[] coeff;
		coeff = NULL;
	}
	else if (d != deg) {
		debug_handler_c(DM_FLP, "in member - function remove_leading_zeros()",
				DV_FLP + 1, std::cout << "new degree is " << d << std::endl << std::flush);
		deg = d;
		np = new bigfloat[d + 1];
		memory_handler(np, DM_CP, "remove_leading_zeros() :: "
			       "Error in memory allocation (np)");
		for (i = 0; i <= d; i++)
			np[i] = coeff[i];
		delete[] coeff;
		coeff = np;
	}
}



// instantiate class base_polynomial
template class base_polynomial< bigfloat >;

// instantiate friend functions -- arithmetical functions
template void negate(base_polynomial< bigfloat > & c,
                     const base_polynomial< bigfloat > &a);

template void add(base_polynomial< bigfloat > & c,
                  const base_polynomial< bigfloat > & a,
                  const base_polynomial< bigfloat > & b);
template void add(base_polynomial< bigfloat > & c,
                  const base_polynomial< bigfloat > & a,
                  const bigfloat & b);
template void add(base_polynomial< bigfloat > & c,
                  const bigfloat & b,
                  const base_polynomial< bigfloat > & a);

template void subtract(base_polynomial< bigfloat > & c,
                       const base_polynomial< bigfloat > & a,
                       const base_polynomial< bigfloat > & b);
template void subtract(base_polynomial< bigfloat > & c,
                       const base_polynomial< bigfloat > & a,
                       const bigfloat & b);
template void subtract(base_polynomial< bigfloat > & c,
                       const bigfloat & b,
                       const base_polynomial< bigfloat > & a);

template void multiply(base_polynomial< bigfloat > & c,
                       const base_polynomial< bigfloat > & a,
                       const base_polynomial< bigfloat > & b);
template void multiply(base_polynomial< bigfloat > & c,
                       const base_polynomial< bigfloat > & a,
                       const bigfloat & b);
template void multiply(base_polynomial< bigfloat > & c,
                       const bigfloat & b,
                       const base_polynomial< bigfloat > & a);

template void power(base_polynomial< bigfloat > & c,
                    const base_polynomial< bigfloat > & a, const bigint & b);

template void derivative(base_polynomial< bigfloat > &c,
                         const base_polynomial< bigfloat > &a);

template base_polynomial< bigfloat >
derivative(const base_polynomial< bigfloat > &a);

// instantiate friend functions -- operators
template bool operator == (const base_polynomial< bigfloat > &a,
			   const base_polynomial< bigfloat > &b);
template bool operator != (const base_polynomial< bigfloat > &a,
			   const base_polynomial< bigfloat > &b);
template base_polynomial< bigfloat >
operator -(const base_polynomial< bigfloat > &a);

template base_polynomial< bigfloat >
operator +(const base_polynomial< bigfloat > &a,
	   const base_polynomial< bigfloat > &b);

template base_polynomial< bigfloat >
operator +(const base_polynomial< bigfloat > &a,
	   const bigfloat & b);

template base_polynomial< bigfloat >
operator +(const bigfloat & b,
	   const base_polynomial< bigfloat > &a);

template base_polynomial< bigfloat >
operator -(const base_polynomial< bigfloat > &a,
	   const base_polynomial< bigfloat > &b);

template base_polynomial< bigfloat >
operator -(const base_polynomial< bigfloat > &a,
	   const bigfloat &b);

template base_polynomial< bigfloat >
operator -(const bigfloat &a,
	   const base_polynomial< bigfloat > &b);

template base_polynomial< bigfloat >
operator *(const base_polynomial< bigfloat > &a,
	   const base_polynomial< bigfloat > &b);

template base_polynomial< bigfloat >
operator *(const base_polynomial< bigfloat > &a,
	   const bigfloat &b);

template base_polynomial< bigfloat >
operator *(const bigfloat &b,
	   const base_polynomial< bigfloat > &a);

template std::istream & operator >> (std::istream &,
				     base_polynomial< bigfloat > &);
template std::ostream & operator << (std::ostream &,
				     const base_polynomial< bigfloat > &);

// instantiate friend functions -- access functions
template bigfloat lead_coeff(const base_polynomial< bigfloat > &a);
template bigfloat const_term(const base_polynomial< bigfloat > &a);

// instantiate class field_polynomial
template class field_polynomial< bigfloat >;

// instantiate friend functions -- arithmetical functions
template void div_rem(field_polynomial< bigfloat > &q,
                      field_polynomial< bigfloat > &r,
                      const base_polynomial< bigfloat > &a,
                      const base_polynomial< bigfloat > &b);
template void divide(field_polynomial< bigfloat > & c,
                     const base_polynomial< bigfloat > & a,
                     const bigfloat & b);
template void divide(field_polynomial< bigfloat > & q,
		     const base_polynomial< bigfloat > & a,
		     const base_polynomial< bigfloat > & b);
template void remainder(field_polynomial< bigfloat > & r,
			const base_polynomial< bigfloat > & a,
			const base_polynomial< bigfloat > & b);
template void power_mod(field_polynomial< bigfloat > & c,
                        const base_polynomial< bigfloat > & a,
                        const bigint & b,
                        const base_polynomial< bigfloat > & f);

template void integral(field_polynomial< bigfloat > &c,
                       const base_polynomial< bigfloat > &a);

template void gcd(field_polynomial< bigfloat > &d, const base_polynomial< bigfloat > &aa,
		  const base_polynomial< bigfloat > &bb);

template void xgcd(field_polynomial< bigfloat > &d,
		   field_polynomial< bigfloat > &x,
		   field_polynomial< bigfloat > &y,
		   const base_polynomial< bigfloat > &aa,
		   const base_polynomial< bigfloat > &bb);



#undef DV_FLP
#undef DM_FLP
#undef LP_ERROR



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
