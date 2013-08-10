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
#include	"LiDIA/bigint.h"
#include	"LiDIA/polynomial.h"

#ifndef LIDIA_INCLUDE_CC
# include	"LiDIA/base/poly_intern.cc"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



#define DV_IP LDBL_UNIPOL+8
#define DM_IP "polynomial< bigint >"
#define LP_ERROR poly_error_msg



// instantiate class
template class base_polynomial< bigint >;

// instantiate friend functions -- arithmetical functions
template void negate(base_polynomial< bigint > & c,
                     const base_polynomial< bigint > &a);

template void add(base_polynomial< bigint > & c,
                  const base_polynomial< bigint > & a,
                  const base_polynomial< bigint > & b);
template void add(base_polynomial< bigint > & c,
                  const base_polynomial< bigint > & a, const bigint & b);
template void add(base_polynomial< bigint > & c,
                  const bigint & b, const base_polynomial< bigint > & a);

template void subtract(base_polynomial< bigint > & c,
                       const base_polynomial< bigint > & a,
                       const base_polynomial< bigint > & b);
template void subtract(base_polynomial< bigint > & c,
                       const base_polynomial< bigint > & a, const bigint & b);
template void subtract(base_polynomial< bigint > & c,
                       const bigint & b, const base_polynomial< bigint > & a);

template void multiply(base_polynomial< bigint > & c,
                       const base_polynomial< bigint > & a,
                       const base_polynomial< bigint > & b);
template void multiply(base_polynomial< bigint > & c,
                       const base_polynomial< bigint > & a, const bigint & b);
template void multiply(base_polynomial< bigint > & c,
                       const bigint & b, const base_polynomial< bigint > & a);

template void power(base_polynomial< bigint > & c,
                    const base_polynomial< bigint > & a, const bigint & b);

template void derivative(base_polynomial< bigint > &c,
                         const base_polynomial< bigint > &a);
template base_polynomial< bigint > derivative(const base_polynomial< bigint > &a);

// instantiate friend functions -- operators
template bool operator == (const base_polynomial< bigint > &a,
			   const base_polynomial< bigint > &b);
template bool operator != (const base_polynomial< bigint > &a,
			   const base_polynomial< bigint > &b);
template base_polynomial< bigint > operator -(const base_polynomial< bigint > &a);

template base_polynomial< bigint > operator +(const base_polynomial< bigint > &a,
					      const base_polynomial< bigint > &b);

template base_polynomial< bigint > operator +(const base_polynomial< bigint > &a,
					      const bigint & b);

template base_polynomial< bigint > operator +(const bigint & b,
					      const base_polynomial< bigint > &a);

template base_polynomial< bigint > operator -(const base_polynomial< bigint > &a,
					      const base_polynomial< bigint > &b);

template base_polynomial< bigint > operator -(const base_polynomial< bigint > &a,
					      const bigint &b);

template base_polynomial< bigint > operator -(const bigint &a,
					      const base_polynomial< bigint > &b);

template base_polynomial< bigint > operator *(const base_polynomial< bigint > &a,
					      const base_polynomial< bigint > &b);

template base_polynomial< bigint > operator *(const base_polynomial< bigint > &a,
					      const bigint &b);

template base_polynomial< bigint > operator *(const bigint &b,
					      const base_polynomial< bigint > &a);

template std::istream & operator >> (std::istream &, base_polynomial< bigint > &);
template std::ostream & operator << (std::ostream &, const base_polynomial< bigint > &);

// instantiate friend functions -- access functions
template bigint lead_coeff(const base_polynomial< bigint > &a);
template bigint const_term(const base_polynomial< bigint > &a);



//
// Division
//

bigint cont(const base_polynomial< bigint > &a)
{
	lidia_size_t d = a.degree();
	if (d < 0) return bigint(0);
	bigint g = a[d];
	for (lidia_size_t i = d - 1; i >= 0; g = gcd(a[i--], g));
	return g;
}



polynomial< bigint > pp(const base_polynomial< bigint > &a)
{
	if (a.is_zero())
		return a;
	else
		return (a / cont(a));
}



void div_rem(polynomial< bigint > &q, polynomial< bigint > &r,
	     const base_polynomial< bigint > &aa, const base_polynomial< bigint > &bb)
{
	debug_handler_l(DM_IP, "in friend - function "
			"div_rem (polynomial< bigint > &, "
			"polynomial< bigint > &, "
			"const base_polynomial< bigint > &, "
			"const base_polynomial< bigint > &)", DV_IP);

	lidia_size_t deg_a = aa.degree(), deg_b = bb.degree();
	lidia_size_t deg_ab = deg_a - deg_b;
	if (deg_b < 0) {
		lidia_error_handler("polynomial< bigint >", "div_rem::division by zero");
		return;
	}

	bigint pow(1);
	if (deg_ab < 0) {
		r.assign(aa);
		q.assign_zero();
	}
	else {
		const bigint *bp;
		bigint *qp, *rp, *tmp, x, y, z;
		lidia_size_t i, j, e = deg_ab + 1;

		polynomial< bigint > a = aa, b = bb;
		q.set_degree(deg_ab);
		r.assign(a);

		x = b.coeff[deg_b];

		for (i = 0, qp = q.coeff + deg_ab; i <= deg_ab; i++, qp--, e--) {
			rp = r.coeff + deg_a - i;
			y = (*rp);
			for (j = deg_ab+1, tmp = q.coeff; j; j--, tmp++)
				multiply(*tmp, x, *tmp);
			for (j = deg_a+1, tmp = r.coeff; j; j--, tmp++)
				multiply(*tmp, x, *tmp);
			*qp = y;
			for (j = deg_b + 1, bp = b.coeff + deg_b; j; j--, rp--, bp--) {
				multiply(z, y, *bp);
				subtract(*rp, *rp, z);
			}
		}
		q.remove_leading_zeros();
		power(pow, x, e);
		multiply(q, q, pow);
	}
	r.remove_leading_zeros();
	multiply(r, r, pow);
}



void divide(polynomial< bigint > & c,
            const base_polynomial< bigint > & a, const bigint & b)
{
	debug_handler_l(DM_IP, "in friend - function "
			"divide (polynomial< bigint > &, "
			"const base_polynomial< bigint > &, "
			"const bigint &)", DV_IP);
	const bigint *ap = ((polynomial< bigint > *)(&a))->coeff;
	lidia_size_t deg_a = a.degree();

	c.set_degree(deg_a);
	bigint r, *cp = c.coeff;

	for (register lidia_size_t i = deg_a + 1; i; i--, ap++, cp++) {
		div_rem ((*cp), r, (*ap), b);
		if (!(r.is_zero()))
			lidia_error_handler("polynomial< bigint >",
					    "divide(polynomial< bigint > &, "
					    "const base_polynomial< bigint > &, const bigint &)"
					    " :: division error");
	}
}



void power_mod(polynomial< bigint > & c,
               const base_polynomial< bigint > & a, const bigint & b,
               const base_polynomial< bigint > & f)
{
	debug_handler_l(DM_IP, "in friend - function "
			"power_mod (polynomial< bigint > &, "
			"const base_polynomial< bigint > &, const bigint &, "
			"const base_polynomial< bigint > &)", DV_IP + 1);
	bigint exponent;
	polynomial< bigint > multiplier;
	if (b.is_negative())
		c.assign_zero();
	else if (b.is_zero() || a.is_one())
		c.assign_one();
	else {
		exponent.assign(b);
		multiplier.assign(a);
		c.assign_one();
		while (exponent.is_gt_zero()) {
			if (!exponent.is_even()) {
				multiply(c, c, multiplier);
				remainder(c, c, f);
			}
			multiply(multiplier, multiplier, multiplier);
			remainder(multiplier, multiplier, f);
			exponent.divide_by_2();
		}
	}
}



//
// Gcd's
//

polynomial< bigint > gcd(const base_polynomial< bigint > &aa,
			 const base_polynomial< bigint > &bb)
// OOOh, sooooo slooowwww.... (Most naive version, should be improved
//                             to use either the sub-resultant-algorithm
//                             or modular techniques)
{
	debug_handler_l(DM_IP, "in friend - function "
			"gcd (const base_polynomial< bigint > &, "
			"const base_polynomial< bigint > &)", DV_IP + 2);
	if (bb.is_zero())
		return aa;
	polynomial< bigint > q, r, a = pp(aa), b = pp(bb);
	bigint cd = gcd(cont(aa), cont(bb));
	do {
		div_rem(q, r, a, b);
		a.assign(b);
		b.assign(pp(r));
	} while (b.deg >= 0);
	return cd * a;
}



polynomial< bigint > xgcd(polynomial< bigint > &x,
			  polynomial< bigint > &y,
			  const base_polynomial< bigint > &aa,
			  const base_polynomial< bigint > &bb)
// OOOh, sooooo slooowwww.... (Most naive version, should be improved
//                             to use sub-resultant-algorithm........)
{
	debug_handler_l(DM_IP, "in friend - function "
			"xgcd (polynomial< bigint > &, polynomial< bigint > &, "
			"const base_polynomial< bigint > &, "
			"const base_polynomial< bigint > &)", DV_IP + 2);
	polynomial< bigint > u0, v0, u2, v2, q, r, a = aa, b = bb;

	if (b.deg < 0) {
		x.assign_one();
		y.assign_zero();
		return a;
	}
	x.assign_one();
	y.assign_zero();
	u2.assign_zero();
	v2.assign_one();

	do {
		div_rem(q, r, a, b);
		a.assign(b);
		b.assign(pp(r));
		u0.assign(u2);
		v0.assign(v2);
		multiply(r, q, u2);
		subtract(u2, x, r);
		u2.assign(pp(u2));
		multiply(r, q, v2);
		subtract(v2, y, r);
		v2.assign(pp(v2));
		x.assign(u0);
		y.assign(v0);
	} while (b.deg >= 0);
	return a;
}



void remainder(polynomial< bigint > & c,
               const base_polynomial< bigint > & a, const bigint & b)
{
	debug_handler_l(DM_IP, "in friend - function "
			"remainder (polynomial< bigint > &, "
			"const base_polynomial< bigint > &, "
			"const bigint &)", DV_IP + 3);
	const bigint *ap;
	bigint *cp;
	bigint r;
	register lidia_size_t deg_a = a.degree();
	register lidia_size_t i;

	c.set_degree(deg_a);

	for (i = deg_a + 1, ap = ((polynomial< bigint > *)(&a))->coeff, cp = c.coeff;
	     i; i--, ap++, cp++) {
		remainder (*cp, *ap, b);
		if (cp->is_negative()) add(*cp, *cp, b);
	}
	c.remove_leading_zeros();
}



// The following function was contributed by Nigel Smart:

lidia_size_t no_of_real_roots(const base_polynomial< bigint > & poly_T)
{
	debug_handler_l(DM_IP, "in friend - function "
			"no_of_real_roots (const base_polynomial< bigint > &)",
			DV_IP + 8);

	if (poly_T.degree() <= 0)
		return 0;

	if (poly_T.degree() == 1)
		return 1;

	polynomial< bigint > A = pp(poly_T);
	polynomial< bigint > B = pp(derivative(poly_T));

	lidia_size_t n = A.deg, s = lead_coeff(A).sign(), t,
		r1 = 1, del, te, ll, degr = 1;

	bigint g = 1, h = 1;

	if (n%2)
		t = s;
	else
		t = -s;

	polynomial< bigint > R, Q;
	while (degr) {
		del = A.deg - B.deg;
		div_rem(Q, R, A, B);
		if (R.deg < 0) {
			lidia_error_handler("sturm", "no_of_real_roots::argument was not square free");
			return 0;
		}
		if (lead_coeff(B) > 0 || del%2 == 1) {
			R = -R;
		}
		degr = R.deg;
		te = lead_coeff(R).sign();
		if (te != s) {
			s = -s;
			r1 = r1-1;
		}
		if (degr%2)
			ll = -t;
		else
			ll = t;
		if (te != ll) {
			t = -t;
			r1++;
		}
		if (degr != 0) {
			A = B;
			bigint bi, bi2;
			power(bi, h, del);
			multiply(bi, g, bi);
			divide(B, R, bi);
			g = abs(lead_coeff(A));
			power(bi2, g, del);
			power(bi, h, del-1);
			divide(h, bi2, bi);
		}
	}
	return r1;
}



// The following two functions were contributed
//  by Roland Dreier (dreier@math.berkeley.edu)

// Calculate resultant of aa and bb via subresultant algorithm.
// Algorithm cribbed verbatim from Algorithm 3.3.7 of H. Cohen's "A
// course in computational algebraic number theory."  (Even the
// variables are named pretty much the same as in his book!)
//
// Author      : Roland Dreier (dreier@math.berkeley.edu)
//
bigint
resultant(const polynomial< bigint > &aa,
	  const polynomial< bigint > &bb)
{
	// Return zero if either polynomial is zero.
	if ((aa.is_zero()) || (bb.is_zero())) {
		return(bigint(0));
	}
	// otherwise...

	// Initialization steps:

	bigint acont, bcont;
	polynomial< bigint > a, b;
	bool neg = false;

	// Maybe skip one reduction by making sure deg(a) >= deg(b).
	if (aa.degree() >= bb.degree()) {
		a.assign(aa);
		b.assign(bb);
	} else {
		a.assign(bb);
		b.assign(aa);
		// Possibly change sign!!
		neg = ((a.degree() % 2) && (b.degree() % 2));
	}

	acont.assign(cont(a));
	bcont.assign(cont(b));

	divide(a, a, acont);
	divide(b, b, bcont);

	bigint g = 1;
	bigint h = 1;

	bigint t;
	bigint pow_temp;

	power(t, acont, int(b.degree()));
	power(pow_temp, bcont, int(a.degree()));
	multiply(t, t, pow_temp);

	// Loop over pseudo-division and reduction steps:

	lidia_size_t delta;
	polynomial< bigint > r;

	do {
		delta = a.degree() - b.degree();

		if ((a.degree() % 2) && (b.degree() % 2)) {
			neg = !neg;
		}
		remainder(r, a, b);

		a.assign(b);

		power(pow_temp, h, delta);
		divide(b, r, pow_temp);
		divide(b, b, g);

		g.assign(a.lead_coeff());

		power(pow_temp, g, delta--);
		if (delta <= 0) {
			power(h, h, -delta);
			multiply(h, h, pow_temp);
		}
		else {
			power(h, h, delta);
			divide(h, pow_temp, h);
		}
	} while (b.degree() > 0);

	// Finish up:

	power(pow_temp, b.lead_coeff(), a.degree());
	delta = a.degree()-1;
	if (delta <= 0) {
		power(h, h, -delta);
		multiply(h, h, pow_temp);
	}
	else {
		power(h, h, delta);
		divide(h, pow_temp, h);
	}
	if (neg)
		t.negate();
	return t*h;
}



//
// Calculate discriminant of a polynomial using the formula:
//   disc(a) = (-1)^(d*(d-1)/2) * resultant(a, a') / lead_coeff(a)
// where d = deg(f).
//
// Rather than raising -1 to a power, just look at d mod 4.
// Remember, d*(d-1)/2 is even iff d = 0 or 1 mod 4.
//
// Author      : Roland Dreier (dreier@math.berkeley.edu)
//
bigint
discriminant(const polynomial< bigint > &a)
{
	if (a.degree() <= 0) {
		return(bigint(0));
	}

	if ((int(a.degree()) % 4) <= 1) {
		return(resultant(a, derivative(a)) / a.lead_coeff());
	} else {
		return(-(resultant(a, derivative(a)) / a.lead_coeff()));
	}
}



#undef DV_IP
#undef DM_IP
#undef LP_ERROR



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
