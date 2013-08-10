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
//	Author	: Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/gf_polynomial.h"
#include	"LiDIA/base/poly_intern.cc"
#include	"LiDIA/field_polynomial.cc"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



static bool inline is_undefined(const galois_field &K)
{
	return K.degree() == 0;
}


template<>
void base_polynomial< gf_element >::set_degree(lidia_size_t d)
{
	debug_handler_l("base_polynomial< gf_element >", "in member - function "
			"set_degree(lidia_size_t)", LDBL_UNIPOL + 2);

	if (d < -1)
		lidia_error_handler_para(d, "d", "d >= -1",
					 "void base_polynomial< gf_element >::"
					 "set_degree(lidia_size_t d)",
					 "base_polynomial< gf_element >", poly_error_msg[2]);

	if (d == deg)
		return;

	if (d < 0) {
		delete[] coeff; // Note: coeff != NULL, since otherwise d == deg !
		coeff = NULL;
		deg = d;
		return;
	}

	gf_element *tmp = coeff;
	coeff = new gf_element [d + 1];
	memory_handler(coeff, "base_polynomial< gf_element >",
		       "set_degree(lidia_size_t d) :: Error in memory allocation (coeff)");

	register lidia_size_t minimum = (d < deg)? d : deg;

	gf_element null(*gf_polynomial::FIELD);
	register lidia_size_t i;
	for (i = 0; i <= minimum; i++) {
		if (is_undefined(tmp[i].get_field()))
			LiDIA::add(coeff[i], tmp[i], null); // set field
		else
			coeff[i] = tmp[i];
	}
	for (i = minimum+1; i <= d; i++)
		coeff[i] = null;

	if (tmp != NULL)
		delete[] tmp;
	deg = d;
}



template<>
gf_element base_polynomial< gf_element >::lead_coeff() const
{
	if (deg < 0) return gf_element();
	else         return coeff[deg];
}



template<>
gf_element base_polynomial< gf_element >::const_term() const
{
	if (deg < 0) return gf_element();
	else         return coeff[0];
}



void derivative(base_polynomial< gf_element > &c,
		const base_polynomial< gf_element > & a)
{
	debug_handler_l("base_polynomial< gf_element >", "in member - function "
			"derivative (const base_polynomial< T > &)", LDBL_UNIPOL+5);

	lidia_size_t d = a.degree();

	if (d <= 0) {
		c.set_degree(-1);
		return;
	}

	c.set_degree(d - 1);
	//const gf_element *ap = a.coeff + 1;
	//gf_element *cp = c.coeff;
	for (lidia_size_t i = 1; i <= d; i++) {
		LiDIA::multiply(c[i], a[i], i);
	}
	c.remove_leading_zeros(); // necessary if characteristic != 0
}



void integral(field_polynomial< gf_element > &c,
	      const base_polynomial< gf_element > & a)
{
	debug_handler_l("field_polynomial< gf_element >", "in member - function "
			"integral (const base_polynomial< T > &)", LDBL_UNIPOL+8 +1);

	lidia_size_t d = a.degree();
	if (d < 0) {
		c.set_degree(-1);
		return;
	}

	c.set_degree(d + 1);

	//const gf_element *ap = ((field_polynomial< gf_element > *)(&a))->coeff;
	//gf_element *cp = c.coeff;
	c[0].assign_zero(*gf_polynomial::FIELD);
	//cp++;

	for (register lidia_size_t i = 0; i <= d; i++) {
		LiDIA::divide(c[i + 1], a[i], i + 1);
	}
}



template class base_polynomial< gf_element >;
template class field_polynomial< gf_element >;

// instantiate friend functions -- arithmetical functions
template void negate(base_polynomial< gf_element > & c,
		     const base_polynomial< gf_element > &a);

template void add(base_polynomial< gf_element > & c,
		  const base_polynomial< gf_element > & a,
		  const base_polynomial< gf_element > & b);
template void add(base_polynomial< gf_element > & c,
		  const base_polynomial< gf_element > & a,
		  const gf_element & b);
template void add(base_polynomial< gf_element > & c,
		  const gf_element & b,
		  const base_polynomial< gf_element > & a);

template void subtract(base_polynomial< gf_element > & c,
		       const base_polynomial< gf_element > & a,
		       const base_polynomial< gf_element > & b);
template void subtract(base_polynomial< gf_element > & c,
		       const base_polynomial< gf_element > & a,
		       const gf_element & b);
template void subtract(base_polynomial< gf_element > & c,
		       const gf_element & b,
		       const base_polynomial< gf_element > & a);

template void multiply(base_polynomial< gf_element > & c,
		       const base_polynomial< gf_element > & a,
		       const base_polynomial< gf_element > & b);
template void multiply(base_polynomial< gf_element > & c,
		       const base_polynomial< gf_element > & a,
		       const gf_element & b);
template void multiply(base_polynomial< gf_element > & c,
		       const gf_element & b,
		       const base_polynomial< gf_element > & a);

template void power(base_polynomial< gf_element > & c,
		    const base_polynomial< gf_element > & a, const bigint & b);

//template void derivative(base_polynomial <gf_element> &c,
//			 const base_polynomial <gf_element> &a);

// instantiate friend functions -- operators
template bool operator == (const base_polynomial< gf_element > &a,
			   const base_polynomial< gf_element > &b);
template bool operator != (const base_polynomial< gf_element > &a,
			   const base_polynomial< gf_element > &b);

// instantiate friend functions -- arithmetical functions
template void div_rem(field_polynomial< gf_element > &q,
		      field_polynomial< gf_element > &r,
		      const base_polynomial< gf_element > &a,
		      const base_polynomial< gf_element > &b);
template void divide(field_polynomial< gf_element > &q,
		     const base_polynomial< gf_element > &a,
		     const base_polynomial< gf_element > &b);
template void remainder(field_polynomial< gf_element > &r,
			const base_polynomial< gf_element > &a,
			const base_polynomial< gf_element > &b);
template void divide(field_polynomial< gf_element > & c,
		     const base_polynomial< gf_element > & a,
		     const gf_element & b);
template void power_mod(field_polynomial< gf_element > & c,
			const base_polynomial< gf_element > & a,
			const bigint & b,
			const base_polynomial< gf_element > & f);

//template void integral(field_polynomial <gf_element> &c,
//		       const base_polynomial <gf_element> &a);

template void gcd(field_polynomial< gf_element > &d,
		  const base_polynomial< gf_element > &aa,
		  const base_polynomial< gf_element > &bb);

template void xgcd(field_polynomial< gf_element > &d,
		   field_polynomial< gf_element > &x,
		   field_polynomial< gf_element > &y,
		   const base_polynomial< gf_element > &aa,
		   const base_polynomial< gf_element > &bb);

template std::istream & operator >> (std::istream &,
				     base_polynomial< gf_element > &);
template std::ostream & operator << (std::ostream &,
				     const base_polynomial< gf_element > &);

template void swap(base_polynomial< gf_element > &,
		   base_polynomial< gf_element > &);



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
