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


#ifndef LIDIA_FIELD_POLYNOMIAL_CC_GUARD_
#define LIDIA_FIELD_POLYNOMIAL_CC_GUARD_



#ifndef LIDIA_FIELD_POLYNOMIAL_H_GUARD_
# include	"LiDIA/field_polynomial.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



#define DV_FP LDBL_UNIPOL+8
#define DM_FP "field_polynomial"
#define LP_ERROR poly_error_msg

// FIELD_POLYNOMIAL FUNCTIONS

//
// Division
//

template< class T >
void
field_polynomial< T >::divide (const base_polynomial< T > & a, const T &b)
{
	debug_handler_l(DM_FP, "in member - function "
			"divide (const base_polynomial< T > &, "
			"const T &)", DV_FP);
	const T * ap;
	T * cp;
	lidia_size_t deg_a = a.degree();

	this->set_degree(deg_a);
	T b_inv;
	invert(b_inv, b);

	register lidia_size_t i = deg_a + 1;
	for (ap = ((field_polynomial< T > *)(&a))->coeff, cp = this->coeff;
	     i; i--, ap++, cp++)
		LiDIA::multiply(*cp, *ap, b_inv);
	// No remove_leading_zeros, since zero_divisors do not exist!!
}



template< class T >
void
field_polynomial< T >::div_rem (field_polynomial< T > &r,
				const base_polynomial< T > &aa, const base_polynomial< T > &bb)
{
	debug_handler_l(DM_FP, "in member - function "
			"div_rem (field_polynomial< T > &, "
			"const base_polynomial< T > &, "
			"const base_polynomial< T > &)", DV_FP);

	lidia_size_t deg_a = aa.degree(), deg_b = bb.degree();
	lidia_size_t deg_ab = deg_a - deg_b;

	if (deg_b < 0)
		lidia_error_handler("field_polynomial< T >", "div_rem::division by zero");
	if (deg_ab < 0) {
		r.assign(aa);
		this->assign_zero();
	}
	else {
		const T *ap, *bp;
		T *qp, *rp, x, y, z;
		lidia_size_t i, j;
		field_polynomial< T > a(aa), b(bb);
		this->set_degree(deg_ab);
		r.set_degree(deg_a);

		invert(x, b.coeff[deg_b]);

		for (i = deg_a + 1, rp = r.coeff, ap = a.coeff; i; i--, ap++, rp++)
			*rp = *ap;

		for (i = 0, qp = this->coeff+deg_ab; i <= deg_ab; i++, qp--) {
			rp = r.coeff + deg_a - i;
			LiDIA::multiply(y, x, *rp);
			*qp = y;
			for (j = deg_b + 1, bp = b.coeff + deg_b; j; j--, rp--, bp--) {
				LiDIA::multiply(z, y, *bp);
				LiDIA::subtract(*rp, *rp, z);
			}
		}
		this->remove_leading_zeros();
		r.remove_leading_zeros();
	}
}



template< class T >
void
field_polynomial< T >::power_mod (const base_polynomial< T > & a,
				  const bigint & b, const base_polynomial< T > & f)
{
	debug_handler_l(DM_FP, "in member - function "
			"power_mod (const base_polynomial< T > &, const bigint &, "
			"const base_polynomial< T > &)", DV_FP + 1);
	bigint exponent;
	field_polynomial< T > multiplier, garbage;
	if (b.is_negative())
		this->assign_zero();
	else if (b.is_zero() || a.is_one())
		this->assign_one();
	else {
		exponent.assign(b);
		multiplier = a;
		this->assign_one();
		while (exponent.is_gt_zero()) {
			if (!exponent.is_even()) {
				LiDIA::multiply(*this, *this, multiplier);
				garbage.div_rem(*this, *this, f);
			}
			LiDIA::multiply(multiplier, multiplier, multiplier);
			garbage.div_rem(multiplier, multiplier, f);
			exponent.divide_by_2();
		}
	}
}



//
// functions
//

template< class T >
void
field_polynomial< T >::integral (const base_polynomial< T > & a)
{
	debug_handler_l(DM_FP, "in member - function "
			"integral (const base_polynomial< T > &)", DV_FP + 1);

	lidia_size_t d = a.degree();
	if (d < 0) {
		this->set_degree(-1);
		return;
	}

	this->set_degree(d + 1);

	const T *ap = (static_cast<const field_polynomial< T > *>(&a))->coeff;
	T *cp = this->coeff;
	(*cp) = 0;
	cp++;

	T temp;

	for (register lidia_size_t i = 0; i <= d; i++, cp++, ap++) {
		temp = i + 1; // necessary, since bigcomplex does not
                                // support automatic cast !!
		LiDIA::divide(*cp, *ap, temp);
	}
}



//
// Gcd's
//

template< class T >
void
field_polynomial< T >::gcd (const base_polynomial< T > &aa,
			    const base_polynomial< T > &bb)
{
	debug_handler_l(DM_FP, "in member - function "
			"gcd (const base_polynomial< T > &, "
			"const base_polynomial< T > &)", DV_FP + 2);
	T lc_inv;

	if (bb.is_zero()) {
		LiDIA::invert(lc_inv, aa[aa.degree()]);
		LiDIA::multiply(*this, aa, lc_inv);
	}
	else {
		debug_handler_l("field_polynomial< T >", "gcd(...)::initialize; ", DV_FP-1);
		field_polynomial< T > q, r, a(aa), b(bb);
		debug_handler_l("field_polynomial< T >", "gcd(...)::initialized; ", DV_FP-1);
		do {
			q.div_rem(r, a, b);
			debug_handler_c("field_polynomial< T >", "gcd(...)", DV_FP - 1,
					std::cout << a << " = " << q << b << " + " << r << std::flush);
			a.assign(b);
			b.assign(r);
		} while (b.deg > 0);
		if (b.deg < 0) {
			invert(lc_inv, a.coeff[a.deg]);
			LiDIA::multiply(*this, a, lc_inv);
		}
		else {
			lc_inv = 1;
			this->assign(lc_inv);
		}
	}
}




template< class T >
void
field_polynomial< T >::xgcd (field_polynomial< T > &x,
			     field_polynomial< T > &y, const base_polynomial< T > &aa,
			     const base_polynomial< T > &bb)
{
	debug_handler_l(DM_FP, "in member - function "
			"xgcd (field_polynomial< T > &, "
			"field_polynomial< T > &, "
			"const base_polynomial< T >, "
			"const base_polynomial< T > &)", DV_FP + 2);

	field_polynomial< T > u0, v0, u2, v2, q, r, a = aa, b = bb;

	if (b.deg < 0) {
		*this = a;
		x.assign_one();
		y.assign_zero();
	}
	else {
		x.assign_one();
		y.assign_zero();
		u2.assign_zero();
		v2.assign_one();

		do {
			q.div_rem(r, a, b);
			a.assign(b);
			b.assign(r);
			u0.assign(u2);
			v0.assign(v2);
			LiDIA::multiply(r, q, u2);
			LiDIA::subtract(u2, x, r);
			LiDIA::multiply(r, q, v2);
			LiDIA::subtract(v2, y, r);
			x.assign(u0);
			y.assign(v0);
		} while (b.deg >= 0);
		this->assign(a);
	}

	if (this->deg < 0)
		return;
	T normalizer;
	invert(normalizer, this->lead_coeff());
	LiDIA::multiply(*this, *this, normalizer);
	LiDIA::multiply(x, x, normalizer);
	LiDIA::multiply(y, y, normalizer);
}



#undef DV_FP
#undef DM_FP
#undef LP_ERROR



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_FIELD_POLYNOMIAL_CC_GUARD_
