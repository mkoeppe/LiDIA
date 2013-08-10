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
#include	"LiDIA/multi_bigmod.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//***********************************************************************
//*	functions for finite prime fields of odd characteristic		*
//***********************************************************************

#define ASSERT_FIELDS_MATCH(a, text) if (a.get_ff_rep() != this) lidia_error_handler("galois_field_rep", text": different fields");

#define ASSERT_NOT_NULLPTR(a, text) if (a.rep == 0) lidia_error_handler("galois_field_rep", text": null pointer");


void galois_field_rep::construct_GFp(gf_element& a) const
{
	ASSERT_FIELDS_MATCH(a, "construct");
	if (a.rep != 0)
		lidia_error_handler("galois_field_rep",
				    "construct_GFp: not null pointer");
	a.rep = new bigint();
	memory_handler(a.rep, "gf_element", "construct_GFp: out of memory");
}



void galois_field_rep::destruct_GFp(gf_element& a) const
{
	ASSERT_FIELDS_MATCH(a, "destruct_GFp");
	delete static_cast<bigint*>(a.rep);
	a.rep = 0;
}



void galois_field_rep::copy_GFp(gf_element& a, const gf_element& b) const
{
	if (a.get_ff_rep() != this || a.ff != b.ff)
		lidia_error_handler("galois_field_rep", "copy_GFp: different fields");
	if (a.rep == 0) {
		a.rep = new bigint(*(static_cast<bigint*>(b.rep)));
		memory_handler(a.rep, "gf_element", "copy_GFp: out of memory");
	}
	else
		(static_cast<bigint*>(a.rep))->assign(*(static_cast<bigint*>(b.rep)));
}



void galois_field_rep::as0_GFp(gf_element& a) const
{
	ASSERT_FIELDS_MATCH(a, "as0_GFp");
	ASSERT_NOT_NULLPTR(a, "as0_GFp");

	(static_cast<bigint*>(a.rep))->assign_zero();
}



void galois_field_rep::as1_GFp(gf_element& a) const
{
	ASSERT_FIELDS_MATCH(a, "as1_GFp");
	ASSERT_NOT_NULLPTR(a, "as1_GFp");

	(static_cast<bigint*>(a.rep))->assign_one();
}



bool galois_field_rep::iseq_GFp(const gf_element& a, const gf_element& b) const
{
	if (a.get_ff_rep() != this || a.ff != b.ff)
		lidia_error_handler("galois_field_rep", "iseq_GFp: different fields");
	if (a.rep == 0 || b.rep == 0)
		lidia_error_handler("galois_field_rep", "iseq_GFp: null pointer");
	return (*(static_cast<bigint*>(a.rep))) == (*(static_cast<bigint*>(b.rep)));
}



bool galois_field_rep::is0_GFp(const gf_element& a) const
{
	ASSERT_FIELDS_MATCH(a, "is0_GFp");
	ASSERT_NOT_NULLPTR(a, "is0_GFp");

	return (static_cast<bigint*>(a.rep))->is_zero();
}



bool galois_field_rep::is1_GFp(const gf_element& a) const
{
	ASSERT_FIELDS_MATCH(a, "is1_GFp");
	ASSERT_NOT_NULLPTR(a, "is1_GFp");

	return (static_cast<bigint*>(a.rep))->is_one();
}



void galois_field_rep::neg_GFp(gf_element& a) const
{
	ASSERT_FIELDS_MATCH(a, "neg_GFp");
	ASSERT_NOT_NULLPTR(a, "neg_GFp");

	bigint &A = *(static_cast<bigint*>(a.rep));
	if (!A.is_zero())
		LiDIA::subtract(A, p, A);
}



void galois_field_rep::add_GFp(gf_element& c, const gf_element& a,
			       const gf_element&b) const
{
	if (a.get_ff_rep() != this || a.ff != b.ff || a.ff != c.ff)
		lidia_error_handler("galois_field_rep", "add_GFp: different fields");
	if (c.rep == 0 || a.rep == 0 || b.rep == 0)
		lidia_error_handler("galois_field_rep", "add_GFp: null pointer");
	bigint &C = *(static_cast<bigint*>(c.rep));
	LiDIA::add(C, *(static_cast<bigint*>(a.rep)), *(static_cast<bigint*>(b.rep)));
	if (C.compare(p) >= 0)
		LiDIA::subtract(C, C, p);
}



void galois_field_rep::sub_GFp(gf_element& c, const gf_element& a,
			       const gf_element&b) const
{
	if (a.get_ff_rep() != this || a.ff != b.ff || a.ff != c.ff)
		lidia_error_handler("galois_field_rep", "sub_GFp: different fields");
	if (c.rep == 0 || a.rep == 0 || b.rep == 0)
		lidia_error_handler("galois_field_rep", "sub_GFp: null pointer");
	bigint &C = *(static_cast<bigint*>(c.rep));
	LiDIA::subtract(C, *(static_cast<bigint*>(a.rep)), *(static_cast<bigint*>(b.rep)));
	if (C.is_negative())
		LiDIA::add(C, C, p);
}



void galois_field_rep::mul_GFp(gf_element& c, const gf_element& a,
			       const gf_element&b) const
{
	if (a.get_ff_rep() != this || a.ff != b.ff || a.ff != c.ff)
		lidia_error_handler("galois_field_rep", "mul_GFp: different fields");
	if (c.rep == 0 || a.rep == 0 || b.rep == 0)
		lidia_error_handler("galois_field_rep", "mul_GFp: null pointer");
	bigint &C = *(static_cast<bigint*>(c.rep));
	LiDIA::multiply(C, *(static_cast<bigint*>(a.rep)), *(static_cast<bigint*>(b.rep)));
//	if (!p.is_zero())			// for uninitialized fields, p is zero
	LiDIA::remainder(C, C, p);
}



void galois_field_rep::inv_GFp(gf_element& a) const
{
	ASSERT_FIELDS_MATCH(a, "inv_GFp");
	ASSERT_NOT_NULLPTR(a, "inv_GFp");
	bigint &A = *(static_cast<bigint*>(a.rep));

//    if (p.is_zero())			// for uninitialized fields, p is zero
//    {
//	if (!abs(A).is_one())
//	    lidia_error_handler("galois_field_rep", "inv_GFp: cannot compute "
//		    "inverse (did you forget to initialize the field?)");
//	return;				// 1^{-1} = 1 and -1^{-1} = -1
//    }

	bigint d, u;
	d = xgcd_left(u, A, p);
	if (d.is_one()) {
		LiDIA::remainder(A, u, p);
		if (A.is_negative())
			LiDIA::add(A, A, p);
	}
	else
		lidia_error_handler("galois_field_rep",
				    "inv_GFp: inverse does not exist");
}



void galois_field_rep::sqr_GFp(gf_element& c, const gf_element& a) const
{
	if (a.get_ff_rep() != this || a.ff != c.ff)
		lidia_error_handler("galois_field_rep", "sqr_GFp: different fields");
	if (a.rep == 0 || c.rep == 0)
		lidia_error_handler("galois_field_rep", "sqr_GFp: null pointer");

	bigint &C = *(static_cast<bigint*>(c.rep));
	LiDIA::square(C, *(static_cast<bigint*>(a.rep)));
//    if (!p.is_zero())			// for uninitialized fields, p is zero
	LiDIA::remainder(C, C, p);
}



udigit galois_field_rep::hash_GFp(const gf_element& a) const
{
	ASSERT_FIELDS_MATCH(a, "hash_GFp");
	ASSERT_NOT_NULLPTR(a, "hash_GFp");

	return (static_cast<bigint*>(a.rep))->least_significant_digit();
}



const Fp_polynomial& galois_field_rep::get_pol_rep_GFp(const gf_element& a)
	const
{
	ASSERT_FIELDS_MATCH(a, "get_pol_rep_GFp");
	ASSERT_NOT_NULLPTR(a, "get_pol_rep_GFp");

//	if (p.is_zero())			// for uninitialized fields, p is zero
//		lidia_error_handler("galois_field_rep", "get_pol_rep_GFp: "
//				    "characteristic = 0 (did you forget to initialize the field?)");

	static Fp_polynomial pol;
	pol.set_modulus(p);
	pol.set_coefficient(*(static_cast<bigint*>(a.rep)), 0);
	return pol;
}



void galois_field_rep::set_pol_rep_GFp(gf_element&a, const Fp_polynomial& b) const
{
	ASSERT_FIELDS_MATCH(a, "set_pol_rep_GFp");
	ASSERT_NOT_NULLPTR(a, "set_pol_rep_GFp");

	(static_cast<bigint*>(a.rep))->assign(b.const_term());
}



#if 0
void galois_field_rep::promote_GFp(gf_element&a, const bigint& b) const
{
	ASSERT_FIELDS_MATCH(a, "promote_GFp");
	ASSERT_NOT_NULLPTR(a, "promote_GFp");

	(static_cast<bigint*>(a.rep))->assign(b);
}
#endif



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
