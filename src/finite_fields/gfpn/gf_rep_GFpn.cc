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
#include	"LiDIA/galois_field.h"
#include	"LiDIA/finite_fields/galois_field_rep.h"
#include	"LiDIA/gf_element.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//***********************************************************************
//*	functions for finite prime fields of odd characteristic		*
//***********************************************************************

#define ASSERT_FIELDS_MATCH(a, text) if (a.get_ff_rep() != this) lidia_error_handler("galois_field_rep", text": different fields");

#define ASSERT_NOT_NULLPTR(a, text) if (a.rep == 0) lidia_error_handler("galois_field_rep", text": null pointer");


void galois_field_rep::construct_GFpn(gf_element& a) const
{
	ASSERT_FIELDS_MATCH(a, "construct");
	if (a.rep != 0)
		lidia_error_handler("galois_field_rep",
				    "construct_GFpn: not null pointer");
	a.rep = new Fp_polynomial();
	memory_handler(a.rep, "gf_element", "construct_GFpn: out of memory");
	(static_cast<Fp_polynomial*>(a.rep))->set_modulus(p);
}



void galois_field_rep::destruct_GFpn(gf_element& a) const
{
	ASSERT_FIELDS_MATCH(a, "destruct_GFpn");
	delete (static_cast<Fp_polynomial*>(a.rep));
	a.rep = 0;
}



void galois_field_rep::copy_GFpn(gf_element& a, const gf_element& b) const
{
	if (a.get_ff_rep() != this || a.ff != b.ff)
		lidia_error_handler("galois_field_rep", "copy_GFpn: different fields");
	if (a.rep == 0) {
		a.rep = new Fp_polynomial(*(static_cast<Fp_polynomial*>(b.rep)));
		memory_handler(a.rep, "gf_element", "copy_GFpn: out of memory");
	}
	else
		(static_cast<Fp_polynomial*>(a.rep))->assign(*(static_cast<Fp_polynomial*>(b.rep)));
}



void galois_field_rep::as0_GFpn(gf_element& a) const
{
	ASSERT_FIELDS_MATCH(a, "as0_GFpn");
	ASSERT_NOT_NULLPTR(a, "as0_GFpn");

	(static_cast<Fp_polynomial*>(a.rep))->assign_zero();
}



void galois_field_rep::as1_GFpn(gf_element& a) const
{
	ASSERT_FIELDS_MATCH(a, "as1_GFpn");
	ASSERT_NOT_NULLPTR(a, "as1_GFpn");

	(static_cast<Fp_polynomial*>(a.rep))->assign_one();
}



bool galois_field_rep::iseq_GFpn(const gf_element& a, const gf_element& b) const
{
	if (a.get_ff_rep() != this || a.ff != b.ff)
		lidia_error_handler("galois_field_rep", "iseq_GFpn: different fields");
	if (a.rep == 0 || b.rep == 0)
		lidia_error_handler("galois_field_rep", "iseq_GFpn: null pointer");
	return (*(static_cast<Fp_polynomial*>(a.rep))) == (*(static_cast<Fp_polynomial*>(b.rep)));
}



bool galois_field_rep::is0_GFpn(const gf_element& a) const
{
	ASSERT_FIELDS_MATCH(a, "is0_GFpn");
	ASSERT_NOT_NULLPTR(a, "is0_GFpn");

	return (static_cast<Fp_polynomial*>(a.rep))->is_zero();
}



bool galois_field_rep::is1_GFpn(const gf_element& a) const
{
	ASSERT_FIELDS_MATCH(a, "is1_GFpn");
	ASSERT_NOT_NULLPTR(a, "is1_GFpn");

	return (static_cast<Fp_polynomial*>(a.rep))->is_one();
}



void galois_field_rep::neg_GFpn(gf_element& a) const
{
	ASSERT_FIELDS_MATCH(a, "neg_GFpn");
	ASSERT_NOT_NULLPTR(a, "neg_GFpn");

	(static_cast<Fp_polynomial*>(a.rep))->negate();
}



void galois_field_rep::add_GFpn(gf_element& c, const gf_element& a,
				const gf_element&b) const
{
	if (a.get_ff_rep() != this || a.ff != b.ff || a.ff != c.ff)
		lidia_error_handler("galois_field_rep", "add_GFpn: different fields");
	if (c.rep == 0 || a.rep == 0 || b.rep == 0)
		lidia_error_handler("galois_field_rep", "add_GFpn: null pointer");

	const Fp_polynomial& A = *(static_cast<Fp_polynomial*>(a.rep));
	const Fp_polynomial& B = *(static_cast<Fp_polynomial*>(b.rep));
	Fp_polynomial& C = *(static_cast<Fp_polynomial*>(c.rep));
	LiDIA::add(C, A, B);
}



void galois_field_rep::sub_GFpn(gf_element& c, const gf_element& a,
				const gf_element&b) const
{
	if (a.get_ff_rep() != this || a.ff != b.ff || a.ff != c.ff)
		lidia_error_handler("galois_field_rep", "sub_GFpn: different fields");
	if (c.rep == 0 || a.rep == 0 || b.rep == 0)
		lidia_error_handler("galois_field_rep", "sub_GFpn: null pointer");
	const Fp_polynomial& A = *(static_cast<Fp_polynomial*>(a.rep));
	const Fp_polynomial& B = *(static_cast<Fp_polynomial*>(b.rep));
	Fp_polynomial& C = *(static_cast<Fp_polynomial*>(c.rep));
	LiDIA::subtract(C, A, B);
}



void galois_field_rep::mul_GFpn(gf_element& c, const gf_element& a,
				const gf_element&b) const
{
	if (a.get_ff_rep() != this || a.ff != b.ff || a.ff != c.ff)
		lidia_error_handler("galois_field_rep", "mul_GFpn: different fields");
	if (c.rep == 0 || a.rep == 0 || b.rep == 0)
		lidia_error_handler("galois_field_rep", "mul_GFpn: null pointer");
	const Fp_polynomial& A = *(static_cast<Fp_polynomial*>(a.rep));
	const Fp_polynomial& B = *(static_cast<Fp_polynomial*>(b.rep));
	Fp_polynomial& C = *(static_cast<Fp_polynomial*>(c.rep));

	multiply(C, A, B, poly_mod);
}



void galois_field_rep::inv_GFpn(gf_element& a) const
{
	ASSERT_FIELDS_MATCH(a, "inv_GFpn");
	ASSERT_NOT_NULLPTR(a, "inv_GFpn");

	Fp_polynomial &A = *(static_cast<Fp_polynomial*>(a.rep));
	Fp_polynomial s, t, d;
	xgcd(d, s, t, A, poly_mod.modulus());
	// a s + b t = d = gcd(A, irred.poly)

	//this should never happen :
	if (!d.is_one()) lidia_error_handler("galois_field_rep",
					     "inv_GFpn()::cannot compute inverse");

	A.assign(s);
}



void galois_field_rep::sqr_GFpn(gf_element& c, const gf_element& a) const
{
	if (a.get_ff_rep() != this || a.ff != c.ff)
		lidia_error_handler("galois_field_rep", "sqr_GFpn: different fields");
	if (a.rep == 0 || c.rep == 0)
		lidia_error_handler("galois_field_rep", "sqr_GFpn: null pointer");

	square(*(static_cast<Fp_polynomial*>(c.rep)), *(static_cast<Fp_polynomial*>(a.rep)), poly_mod);
}



udigit galois_field_rep::hash_GFpn(const gf_element& a) const
{
	ASSERT_FIELDS_MATCH(a, "hash_GFpn");
	ASSERT_NOT_NULLPTR(a, "hash_GFpn");

	const Fp_polynomial &A = *(static_cast<Fp_polynomial*>(a.rep));
	udigit s = 0;
	lidia_size_t i;
	for (i = 0; i < A.degree(); i++)
		s ^= static_cast<udigit>(A[i].least_significant_digit());
	return s % max_udigit();
}



void galois_field_rep::set_pol_rep_GFpn(gf_element& a, const Fp_polynomial& b) const
{
	ASSERT_FIELDS_MATCH(a, "set_pol_rep_GFpn");
	ASSERT_NOT_NULLPTR(a, "set_pol_rep_GFpn");

	Fp_polynomial &A = *(static_cast<Fp_polynomial*>(a.rep));
	if (b.degree() < static_cast<lidia_size_t>(degree()))
		A.assign(b);
	else
		remainder(A, b, poly_mod);
}



const Fp_polynomial& galois_field_rep::get_pol_rep_GFpn(const gf_element& a)
	const
{
	ASSERT_FIELDS_MATCH(a, "get_pol_rep_GFpn");
	ASSERT_NOT_NULLPTR(a, "get_pol_rep_GFpn");

	return *(static_cast<Fp_polynomial*>(a.rep));
}



#if 0
void galois_field_rep::promote_GFpn(gf_element&a, const bigint& b) const
{
	ASSERT_FIELDS_MATCH(a, "promote_GFpn");
	ASSERT_NOT_NULLPTR(a, "promote_GFpn");

	(static_cast<Fp_polynomial*>(a.rep))->assign(b);
}
#endif



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
