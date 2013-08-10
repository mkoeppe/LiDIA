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
#include	"LiDIA/random_generator.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//***********************************************************************
//*	functions for finite fields of even characteristic		*
//***********************************************************************

#define ASSERT_FIELDS_MATCH(a, text) if (a.get_ff_rep() != this) lidia_error_handler("galois_field_rep", text": different fields");

#define ASSERT_NOT_NULLPTR(a, text) if (a.rep == 0) lidia_error_handler("galois_field_rep", text": null pointer");


void galois_field_rep::construct_GF2n(gf_element& a) const
{
	ASSERT_FIELDS_MATCH(a, "construct");
	if (a.rep != 0)
		lidia_error_handler("galois_field_rep",
				    "construct_GF2n: not null pointer");

	a.rep = new udigit[I2.anzBI];
	memory_handler(a.rep, "gf_element", "construct_GF2n: out of memory");
	register unsigned int i;
	for (i = 0; i < I2.anzBI; i++)
		(static_cast<udigit*>(a.rep))[i] = static_cast<udigit>(0);
}



void galois_field_rep::destruct_GF2n(gf_element& a) const
{
	ASSERT_FIELDS_MATCH(a, "destruct_GF2n");
	delete[] static_cast<udigit*>(a.rep);
	a.rep = 0;
}



void galois_field_rep::copy_GF2n(gf_element& a, const gf_element& b) const
{
	if (a.get_ff_rep() != this || a.ff != b.ff)
		lidia_error_handler("galois_field_rep", "copy_GF2n: different fields");
	if (a.rep == 0) {
		a.rep = new udigit[I2.anzBI];
		memory_handler(a.rep, "gf_element", "construct_GF2n: out of memory");
	}
	register unsigned int i;
	for (i = 0; i < I2.anzBI; i++)
		(static_cast<udigit*>(a.rep))[i] = (static_cast<udigit*>(b.rep))[i];
}



void galois_field_rep::as0_GF2n(gf_element& a) const
{
	ASSERT_FIELDS_MATCH(a, "as0_GF2n");
	ASSERT_NOT_NULLPTR(a, "as0_GF2n");

	register unsigned int i;
	for (i = 0; i < I2.anzBI; i++)
		(static_cast<udigit*>(a.rep))[i] = static_cast<udigit>(0);
}



void galois_field_rep::as1_GF2n(gf_element& a) const
{
	ASSERT_FIELDS_MATCH(a, "as1_GF2n");
	ASSERT_NOT_NULLPTR(a, "as1_GF2n");

	register unsigned int i;
	udigit *ap = static_cast<udigit*>(a.rep);
	for (i = 0; i < I2.anzBI; i++)
		ap[i] = static_cast<udigit>(0);
	ap[0] = 1;
}



bool galois_field_rep::iseq_GF2n(const gf_element& a, const gf_element& b) const
{
	if (a.get_ff_rep() != this || a.ff != b.ff)
		lidia_error_handler("galois_field_rep", "iseq_GF2n: different fields");
	if (a.rep == 0 || b.rep == 0)
		lidia_error_handler("galois_field_rep", "iseq_GF2n: null pointer");

	udigit *ap = static_cast<udigit*>(a.rep),
		*bp = static_cast<udigit*>(b.rep);
	if (!I2.is_reduced(ap))
		I2.partial_reduce2(ap);
	if (!I2.is_reduced(bp))
		I2.partial_reduce2(bp);

	register int i = I2.anzBI-1;
	while (i >= 0 && (ap[i] == bp[i]))
		i--;

	if (i < 0) return true;
	else       return false;
}



bool galois_field_rep::is0_GF2n(const gf_element& a) const
{
	ASSERT_FIELDS_MATCH(a, "is0_GF2n");
	ASSERT_NOT_NULLPTR(a, "is0_GF2n");

	udigit *ap = static_cast<udigit*>(a.rep);

	if (!I2.is_reduced(ap))
		I2.partial_reduce2(ap);

	register int i = I2.anzBI-1;
	while ((i >= 0) && (ap[i] == static_cast<udigit>(0)))
		i--;

	if (i < 0) return true;
	else       return false;
}



bool galois_field_rep::is1_GF2n(const gf_element& a) const
{
	ASSERT_FIELDS_MATCH(a, "is1_GF2n");
	ASSERT_NOT_NULLPTR(a, "is1_GF2n");

	udigit *ap = static_cast<udigit*>(a.rep);
	if (!I2.is_reduced(ap))
		I2.partial_reduce2(ap);

	register int i = I2.anzBI-1;
	while (i >= 0 && ap[i] == static_cast<udigit>(0))
		i--;

	if (i != 0) return false;
	else        return (ap[0] == static_cast<udigit>(1));
}



void galois_field_rep::neg_GF2n(gf_element& a) const
{

	ASSERT_FIELDS_MATCH(a, "neg_GF2n");
	ASSERT_NOT_NULLPTR(a, "neg_GF2n");

	// nothing to be done since a = -a in GF(2^n)
}



void galois_field_rep::add_GF2n(gf_element& c, const gf_element& a,
				const gf_element&b) const
{
	if (a.get_ff_rep() != this || a.ff != b.ff || a.ff != c.ff)
		lidia_error_handler("galois_field_rep", "add_GF2n: different fields");
	if (c.rep == 0 || a.rep == 0 || b.rep == 0)
		lidia_error_handler("galois_field_rep", "add_GF2n: null pointer");

	register unsigned int i;
	register udigit *cp = static_cast<udigit*>(c.rep),
		*ap = static_cast<udigit*>(a.rep),
		*bp = static_cast<udigit*>(b.rep);

	for (i = 0; i < I2.anzBI; i++, cp++, ap++, bp++)
		*cp = (*ap) ^ (*bp);
}



void galois_field_rep::sub_GF2n(gf_element& c, const gf_element& a,
				const gf_element&b) const
{
	if (a.get_ff_rep() != this || a.ff != b.ff || a.ff != c.ff)
		lidia_error_handler("galois_field_rep", "sub_GF2n: different fields");
	if (c.rep == 0 || a.rep == 0 || b.rep == 0)
		lidia_error_handler("galois_field_rep", "sub_GF2n: null pointer");

	// identical to addition
	register unsigned int i;
	register udigit *cp = static_cast<udigit*>(c.rep),
		*ap = static_cast<udigit*>(a.rep),
		*bp = static_cast<udigit*>(b.rep);

	for (i = 0; i < I2.anzBI; i++, cp++, ap++, bp++)
		*cp = (*ap) ^ (*bp);
}



void galois_field_rep::mul_GF2n(gf_element& c, const gf_element& a,
				const gf_element&b) const
{
	if (a.get_ff_rep() != this || a.ff != b.ff || a.ff != c.ff)
		lidia_error_handler("galois_field_rep", "mul_GF2n: different fields");
	if (c.rep == 0 || a.rep == 0 || b.rep == 0)
		lidia_error_handler("galois_field_rep", "mul_GF2n: null pointer");

//	if (a.rep == b.rep) {
//		sqr_GF2n(c, a);
//		return;
//	}

	register unsigned int i;
	for (i = 0; i < 2*I2.anzBI; i++)
		I2.B[i] = static_cast<udigit>(0);

	I2.mul(I2.B, static_cast<udigit*>(a.rep), static_cast<udigit*>(b.rep));
	I2.partial_reduce1(I2.B);

	for (i = 0; i < I2.anzBI; i++)
		(static_cast<udigit*>(c.rep))[i] = I2.B[i];
}



void galois_field_rep::inv_GF2n(gf_element& a) const
{
	ASSERT_FIELDS_MATCH(a, "inv_GF2n");
	ASSERT_NOT_NULLPTR(a, "inv_GF2n");

	udigit *ap = static_cast<udigit*>(a.rep);
	register unsigned int i;
	for (i = 0; i < I2.anzBI; i++)
		I2.B[i] = ap[i];

	I2.inv(ap, I2.B);
}



void galois_field_rep::sqr_GF2n(gf_element& c, const gf_element& a) const
{
	if (a.get_ff_rep() != this || a.ff != c.ff)
		lidia_error_handler("galois_field_rep", "sqr_GF2n: different fields");
	if (a.rep == 0 || c.rep == 0)
		lidia_error_handler("galois_field_rep", "sqr_GF2n: null pointer");

	register unsigned int i;
	udigit *ap = static_cast<udigit*>(a.rep);

	I2.square(I2.B, ap);
	I2.partial_reduce1(I2.B);

	for (i = 0; i < I2.anzBI; i++)
		(static_cast<udigit*>(c.rep))[i] = I2.B[i];
}



#if 0
void galois_field_rep::rand_GF2n(gf_element& a) const
{
	ASSERT_FIELDS_MATCH(a, "rand_GF2n");
	ASSERT_NOT_NULLPTR(a, "rand_GF2n");

	static random_generator rg;
	long r;
	register unsigned int i;
	udigit *ap = static_cast<udigit*>(a.rep);

	for (i = 0; i < I2.anzBI; i++)
	{
		rg >> r;
		ap[i] = static_cast<udigit>(r);
	}
	ap[I2.anzBI-1] &= ~((~0) << (degree() % BITS_PER_LONG));
}



bigint galois_field_rep::lift_GF2n(const gf_element& a) const
{
	ASSERT_FIELDS_MATCH(a, "lift_GF2n");
	ASSERT_NOT_NULLPTR(a, "lift_GF2n");

	if (a.is_zero())
		return bigint(0);
	else if (a.is_one())
		return bigint(1);
	else {
// UIUIUIUIUIUIUIUIUIUIUIUIUIUIUI..............
//	lidia_error_handler("gf_element",
//		"lift_GF2n: element not in prime field");
//
		bigint h;
		int i;
		udigit *ap = static_cast<udigit*>(a.rep);

		if (!I2.is_reduced(ap))
			I2.partial_reduce2(ap);

		for (i = static_cast<int>(I2.anzBI)-1; i >= 0; i--) {
			// convert gf2n element into bigint
			shift_left(h, h, BITS_PER_LONG);
			::add(h, h, bigint(ap[i]));
		}
		return h;
	}
	return bigint(0);
}
#endif



udigit galois_field_rep::hash_GF2n(const gf_element& a) const
{
	ASSERT_FIELDS_MATCH(a, "hash_GF2n");
	ASSERT_NOT_NULLPTR(a, "hash_GF2n");

	udigit *ap = static_cast<udigit*>(a.rep);
	if (!I2.is_reduced(ap))
		I2.partial_reduce2(ap);

	return static_cast<udigit>(ap[0]);
}



#if 0
void galois_field_rep::print_GF2n(const gf_element& a, std::ostream& out) const
{
	ASSERT_FIELDS_MATCH(a, "print_GF2n");
	ASSERT_NOT_NULLPTR(a, "print_GF2n");

	udigit *ap = static_cast<udigit*>(a.rep);
	if (!I2.is_reduced(ap))
		I2.partial_reduce2(ap);

	out << "[";
	int i, j;
	udigit tmp;
	tmp = ap[I2.anzBI-1];
	int d = (degree()) % BITS_PER_LONG;
	tmp <<= BITS_PER_LONG-d;
	for (j = d-1; j >= 0; j--) {
		if ((tmp & (1 << (BITS_PER_LONG-1))) == 0) out << "0";
		else                                     out << "1";
		tmp <<= 1;
	}

	for (i = I2.anzBI-2; i >= 0; i--) {
		tmp = ap[i];
		for (j = BITS_PER_LONG-1; j >= 0; j--) {
			if ((tmp & (1 << (BITS_PER_LONG-1))) == 0) out << "0";
			else                                     out << "1";
			tmp <<= 1;
		}
	}
	out << "]_2";

	bigint h;

	if (!I2.is_reduced(ap))
		I2.partial_reduce2(ap);

	for (i = static_cast<int>(I2.anzBI)-1; i >= 0; i--) {
		// convert gf2n element into bigint
		shift_left(h, h, BITS_PER_LONG);
		::add(h, h, bigint(ap[i]));
	}

	out << ", Dec:" << h;
}



void galois_field_rep::promote_GF2n(gf_element&a, const bigint& b) const
{
	ASSERT_FIELDS_MATCH(a, "promote_GF2n");
	ASSERT_NOT_NULLPTR(a, "promote_GF2n");

	if (b.is_even())
		a.assign_zero();
	else
		a.assign_one();
}
#endif



void galois_field_rep::set_pol_rep_GF2n(gf_element& a, const Fp_polynomial& b)
	const
{
	ASSERT_FIELDS_MATCH(a, "set_pol_rep_GF2n");
	ASSERT_NOT_NULLPTR(a, "set_pol_rep_GF2n");

	const Fp_polynomial *pol;
	Fp_polynomial bb;
	if (b.degree() >= static_cast<int>(degree())) {
		remainder(bb, b, irred_polynomial());
		pol = &bb;
	}
	else
		pol = &b;

	int i, j;
	udigit tmp;
	udigit *ap = static_cast<udigit*>(a.rep);
	for (i = 0; i < static_cast<int>(I2.anzBI); i++) {
		tmp = 0;
		for (j = 0; j < BITS_PER_LONG; j++)
			if (!(*pol)[i*BITS_PER_LONG+j].is_zero())
				tmp |= (1 << j);
		ap[i] = tmp;
	}
}



const Fp_polynomial& galois_field_rep::get_pol_rep_GF2n(const gf_element& a)
	const
{
	ASSERT_FIELDS_MATCH(a, "get_pol_rep_GF2n");
	ASSERT_NOT_NULLPTR(a, "get_pol_rep_GF2n");

	static Fp_polynomial pol;
	pol.set_modulus(2);
	pol.set_max_degree(degree()-1);

	udigit *ap = static_cast<udigit*>(a.rep);
	if (!I2.is_reduced(ap))
		I2.partial_reduce2(ap);

	int i, j, deg = degree()-1;
	udigit tmp;
	tmp = ap[I2.anzBI-1];
	int d = (degree()) % BITS_PER_LONG;
	tmp <<= BITS_PER_LONG-d;
	for (j = d-1; j >= 0; j--) {
		if ((tmp & (1UL << (BITS_PER_LONG-1))) != 0)
			pol.set_coefficient(deg);
		deg--;
		tmp <<= 1;
	}

	for (i = I2.anzBI-2; i >= 0; i--) {
		tmp = ap[i];
		for (j = BITS_PER_LONG-1; j >= 0; j--) {
			if ((tmp & (1UL << (BITS_PER_LONG-1))) != 0)
				if ((tmp & (1UL << (BITS_PER_LONG-1))) != 0)
					pol.set_coefficient(deg);
			deg--;
			tmp <<= 1;
		}
	}

	return pol;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
