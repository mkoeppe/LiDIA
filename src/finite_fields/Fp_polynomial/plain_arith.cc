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
//	Author	: Victor Shoup (VS), Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/Fp_polynomial.h"
#include	"LiDIA/udigit.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



// define USE_SINGLE_PREC, if you want to use single prec. subroutines for
// single prec. moduli (see negate, add, subtract)
#define USE_SINGLE_PREC
//#undef USE_SINGLE_PREC

// size of statically allocated array in div_rem routines
#define DIVREM_STATIC_ALLOC_SIZE 50


//***************************************************************
//	This File contains the implementation of the functions
//	-	plain_mul, plain_sqr
//	-	plain_div_rem, plain_div, plain_rem
//	-	plain_inv
//	These always use "classical" algorithms.
//****************************************************************

// compute ap[0..deg_a]*bp[0..deg_b] (polynomial multiplication) mod p
// store result in cp[0..deg_a+deg_b]
static inline void
polmul(bigint *cp, const bigint *ap, lidia_size_t deg_a,
       const bigint *bp, lidia_size_t deg_b, const bigint &p)
{
	lidia_size_t i, j, jmin, jmax, deg_ab = deg_a + deg_b;
	bigint temp, accum;

	for (i = 0; i <= deg_ab; i++) {
		accum.assign_zero();
		jmin = comparator< lidia_size_t >::max(0, i-deg_b);
		jmax = comparator< lidia_size_t >::min(deg_a, i);
		for (j = jmin; j <= jmax; j++) {
			multiply(temp, ap[j], bp[i-j]);
			add(accum, accum, temp);
		}
		Remainder(cp[i], accum, p);
	}
}



#ifdef USE_SINGLE_PREC
// same as above, but uses single prec.
static inline void
polmul_sgl(bigint *cp, const bigint *ap, lidia_size_t deg_a,
	   const bigint *bp, lidia_size_t deg_b, udigit pp)
{
	lidia_size_t i, j, jmin, jmax, deg_ab = deg_a + deg_b;

	udigit accum_l, accum_h, temp_h, temp_l, carry;
	for (i = 0; i <= deg_ab; i++) {
		accum_l = accum_h = 0;
		jmin = comparator< lidia_size_t >::max(0, i-deg_b);
		jmax = comparator< lidia_size_t >::min(deg_a, i);
		for (j = jmin; j <= jmax; j++) {
			// temp_h*B + temp_l = ap[j]*bp[i-j]
			temp_h = multiply(temp_l,
						 ap[j].least_significant_digit(),
						 bp[i-j].least_significant_digit());
			// accum += temp
			carry = add(accum_l, accum_l, temp_l, 0);
			carry = add(accum_h, accum_h, temp_h, carry);
			// if carry, reduce mod pp
			if (carry != 0) {
				accum_h = divide(carry, carry, accum_h, pp);
				accum_l = divide(carry, accum_h, accum_l, pp);
				accum_h = 0;
			}
		}
		accum_h = divide(carry, 0, accum_h, pp);
		cp[i].assign(divide(carry, accum_h, accum_l, pp));
	}
}
#endif	// USE_SINGLE_PREC



void plain_mul(Fp_polynomial &c, const Fp_polynomial &a, const Fp_polynomial &b)
{
	debug_handler_l("Fp_polynomial",
			"plain_mul (Fp_polynomial&, Fp_polynomial&, Fp_polynomial&)", 5);

	lidia_size_t deg_a = a.degree(), deg_b = b.degree();
	if (deg_a < 0 && deg_b < 0) {
		c.set_modulus(a); //assigns zero
		return;
	}
	const bigint *ap, *bp;
	bigint *cp;

	Fp_polynomial tmp_a, tmp_b;
	if (&c == &a) {
		tmp_a.assign(a);
		ap = tmp_a.coeff;
	}
	else
		ap = a.coeff;
	if (&c == &b) {
		tmp_b.assign(b);
		bp = tmp_b.coeff;
	}
	else
		bp = b.coeff;

	c.set_modulus(a);
	c.set_degree(deg_a+deg_b);
	cp = c.coeff;

	const bigint& p = a.modulus();
#ifdef USE_SINGLE_PREC
	register udigit pp = p.least_significant_digit();

	if (p.compare(pp) == 0 && pp < max_udigit_modulus())
		polmul_sgl(cp, ap, deg_a, bp, deg_b, pp);
	else
#endif // USE_SINGLE_PREC
		polmul(cp, ap, deg_a, bp, deg_b, p);

	c.remove_leading_zeros();
}



// ---------------------------------------------------------------------------

// compute ap[0..deg_a]^2 mod p (polynomial multiplication)
// store result in xp[0..2*deg_a]
static inline void
polsqr(bigint *xp, const bigint *ap, lidia_size_t deg_a, const bigint &p)
{
	lidia_size_t i, j, jmin, jmax;
	lidia_size_t m, m2;
	lidia_size_t d = 2*deg_a;
	bigint accum, temp;
	for (i = 0; i <= d; i++) {
		accum.assign_zero();
		jmin = comparator< lidia_size_t >::max(0, i-deg_a);
		jmax = comparator< lidia_size_t >::min(deg_a, i);
		m = jmax - jmin + 1;
		//i>deg_a: m = 2*deg_a - i + 1
		//i<deg_a: m = i + 1
		//i=deg_a: m = i + 1 = deg_a + 1
		//m(i) ist Anzahl der zu add. Produkte fuer Koeff. i
		m2 = m >> 1;
		//m2 = m/2
		//(wegen Symmetrie nur halb so viele; falls m ungerade, wird
		// nach der Schleitfe direkt quadriert)
		jmax = jmin + m2 - 1;
		for (j = jmin; j <= jmax; j++) {
			multiply(temp, ap[j], ap[i-j]);
			add(accum, accum, temp);
		}
		accum.multiply_by_2();
		if (m & 1) {
			square(temp, ap[jmax + 1]);
			add(accum, accum, temp);
		}
		Remainder(xp[i], accum, p);
	}
}



#ifdef USE_SINGLE_PREC
// same as above, but uses single prec.
static inline void
polsqr_sgl(bigint *xp, const bigint *ap, lidia_size_t deg_a, udigit pp)
{
	lidia_size_t i, j, jmin, jmax;
	lidia_size_t m, m2;
	lidia_size_t d = 2*deg_a;

	udigit accum_l, accum_h, temp_h, temp_l, carry, high;
	for (i = 0; i <= d; i++) {
		accum_h = accum_l = 0;
		jmin = comparator< lidia_size_t >::max(0, i-deg_a);
		jmax = comparator< lidia_size_t >::min(deg_a, i);
		m = jmax - jmin + 1; // for explanation of algorithm see above
		m2 = m >> 1;
		jmax = jmin + m2 - 1;
		for (j = jmin; j <= jmax; j++) {
			// temp_h*B + temp_l = ap[j]*ap[i-j]
			temp_h = multiply(temp_l,
						 ap[j].least_significant_digit(),
						 ap[i-j].least_significant_digit());
			// accum += temp
			carry = add(accum_l, accum_l, temp_l, 0);
			carry = add(accum_h, accum_h, temp_h, carry);
			// if carry, reduce mod pp
			if (carry != 0) {
				accum_h = divide(carry, carry, accum_h, pp);
				accum_l = divide(carry, accum_h, accum_l, pp);
				accum_h = 0;
			}
		}
		// accum *= 2
		carry = add(accum_l, accum_l, accum_l, 0);
		high = add(accum_h, accum_h, accum_h, carry);
		if (m & 1) {
			// temp = ap[jmax+1]^2
			temp_l = ap[jmax + 1].least_significant_digit();
			temp_h = multiply(temp_l, temp_l, temp_l);
			// accum += temp
			carry = add(accum_l, accum_l, temp_l, 0);
			carry = add(accum_h, accum_h, temp_h, carry);
			high += carry;
		}
		// reduce (high*B^2 + accum_h*B + accum_l) mod pp
		accum_h = divide(carry, high, accum_h, pp);
		accum_l = divide(carry, accum_h, accum_l, pp);

		xp[i].assign(accum_l);
	}
}
#endif // USE_SINGLE_PREC



void plain_sqr(Fp_polynomial &x, const Fp_polynomial &a)
{
	debug_handler_l ("Fp_polynomial",
			 "plain_sqr(Fp_polynomial&, Fp_polynomial&)", 5);

	if (&x != &a)
		x.set_modulus(a); //assigns zero

	lidia_size_t deg_a = a.degree();
	if (deg_a < 0) {
		x.assign_zero();
		return;
	}

	lidia_size_t d = 2*deg_a;
	const bigint *ap;

	Fp_polynomial temp_poly;
	if (&x == &a) {
		temp_poly.assign(a);
		ap = temp_poly.coeff;
	}
	else
		ap = a.coeff;

	x.set_degree(d);
	bigint *xp = x.coeff;

	const bigint& p = a.modulus();
#ifdef USE_SINGLE_PREC
	register udigit pp = p.least_significant_digit();

	if (p.compare(pp) == 0 && pp < max_udigit_modulus())
		polsqr_sgl(xp, ap, deg_a, pp);
	else
#endif // USE_SINGLE_PREC
		polsqr(xp, ap, deg_a, p);

	x.remove_leading_zeros();
}



// ---------------------------------------------------------------------------

// compute a = qb + r mod p (polynomial division with remainder)
// store result in qp[0..deg_a-deg_b], r[0..deg_b-1]
// (input must not alias output)
static inline void
poldivrem(bigint *qp, bigint *rp, const bigint *ap, lidia_size_t da,
	  const bigint *bp, lidia_size_t db, const bigint &p)
{
	//
	// compute inverse of b.lead_coeff(), if necessary
	//
	bigint LCInv, t, s;
	bool lc_is_one = bp[db].is_one();
	if (!lc_is_one)
		InvMod(LCInv, bp[db], p);


	lidia_size_t i, j;

	//
	// allocate temporary array if necessary
	//
	static bigint *xp = new bigint[DIVREM_STATIC_ALLOC_SIZE];
	bigint *xxp = 0;
	if (da + 1 > DIVREM_STATIC_ALLOC_SIZE) {
		xxp = xp;
		xp = new bigint[da + 1];
		memory_handler(xp, "Fp_polynomial",
			       "plain_div_rem :: Error in memory allocation");
	}

	for (i = 0; i <= da; i++)
		xp[i].assign(ap[i]);

	lidia_size_t dq = da - db;
	for (i = dq; i >= 0; i--) {
		Remainder(t, xp[i+db], p);
		if (!lc_is_one)
			MulMod(t, t, LCInv, p);
		qp[i].assign(t);

		negate(t, t);

		for (j = db-1; j >= 0; j--) {
			multiply(s, t, bp[j]);
			add(xp[i+j], xp[i+j], s);
		}
	}

	//
	// store remainder in r
	//
	for (i = 0; i < db; i++)
		Remainder(rp[i], xp[i], p);

	if (da + 1 > DIVREM_STATIC_ALLOC_SIZE) {
		delete[] xp;
		xp = xxp;
	}
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// same as above, but uses single prec.
#ifdef USE_SINGLE_PREC
static inline void
poldivrem_sgl(bigint *qp, bigint *rp, const bigint *ap, lidia_size_t da,
	      const bigint *bp, lidia_size_t db, udigit p)
{
	//
	// compute inverse of b.lead_coeff(), if necessary
	//
	udigit LCInv = 1;
	bool lc_is_one = bp[db].is_one();
	if (!lc_is_one)
		LCInv = invert_mod(bp[db].least_significant_digit(), p);


	lidia_size_t i, j;

	//
	// allocate temporary array
	//
	static udigit *xp = new udigit[2*DIVREM_STATIC_ALLOC_SIZE];
	udigit *xxp = 0;
	if (da + 1 > DIVREM_STATIC_ALLOC_SIZE) {
		xxp = xp;
		xp = new udigit[(da + 1) << 1]; // allocate twice the amount...
		memory_handler(xp, "Fp_polynomial",
			       "plain_div_rem :: Error in memory allocation");
	}

	udigit *xp_l = xp;
	udigit *xp_h = xp + (da + 1); // ...instead of allocating twice
	for (i = 0; i <= da; i++) {
		xp_l[i] = ap[i].least_significant_digit();
		xp_h[i] = static_cast<udigit>(0);
	}

	lidia_size_t dq = da - db;
	udigit t, s_h, s_l;
	udigit carry;
	for (i = dq; i >= 0; i--) {
		// we know that xp_h[i+db] < p, see (*)
		t = divide(carry, xp_h[i+db], xp_l[i+db], p);
		if (!lc_is_one)
			t = multiply_mod(t, LCInv, p);
		qp[i].assign(t);

		// t = -t
		t = subtract_mod(static_cast<udigit>(0), t, p);

		for (j = db-1; j >= 0; j--) {
			// s_h*B + s_l = t * bp[j]
			s_h = multiply(s_l, t, bp[j].least_significant_digit());
			// xp[i+j] += s
			carry = add(xp_l[i+j], xp_l[i+j], s_l, 0);
			carry = add(xp_h[i+j], xp_h[i+j], s_h, carry);
			// if carry or if xp_h[i+j]>=p, reduce mod p     (*)
			if (carry != 0 || xp_h[i+j] >= p) {
				s_h = divide(carry, carry, xp_h[i+j], p);
				xp_l[i+j] = divide(carry, s_h, xp_l[i+j], p);
				xp_h[i+j] = 0;
			}
		}
	}


	//
	// store remainder in r
	//
	udigit q_tmp;
	for (i = 0; i < db; i++) {
		// we know that xp_h[i] < p, see (*)
		rp[i].assign(divide(q_tmp, xp_h[i], xp_l[i], p));
	}

	if (da + 1 > DIVREM_STATIC_ALLOC_SIZE) {
		delete[] xp;
		xp = xxp;
	}
}
#endif // USE_SINGLE_PREC



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void plain_div_rem(Fp_polynomial& q, Fp_polynomial& r, const Fp_polynomial& a, const Fp_polynomial& b)
{
	debug_handler_l("Fp_polynomial", "plain_div_rem(Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_polynomial&)", 5);

	lidia_size_t da, db, dq;

	da = a.degree();
	db = b.degree();

	if (db < 0) {
		lidia_error_handler("Fp_polynomial", "plain_div_rem(Fp_polynomial&, "
				    "Fp_polynomial&, Fp_polynomial&, Fp_polynomial&)::"
				    "division by zero");
		return;
	}

	if (da < db) {
		r.assign(a);
		q.set_modulus(a); //assigns zero
		return;
	}

	dq = da - db;

	bigint *qp, *rp;
	const bigint *ap, *bp;
	ap = a.coeff;
	bp = b.coeff;

	Fp_polynomial tmp;


	// precautions to take if output aliases input:
	// 1) if q aliases input
	bool use_tmp = false;
	if (&q != &a && &q != &b) {
		q.set_modulus(a);
		q.set_degree(dq);
		qp = q.coeff;
	}
	else {
		tmp.set_modulus(a);
		tmp.set_degree(dq);
		qp = tmp.coeff;
		use_tmp = true;
	}
	// 2) if r aliases input
	bool truncate_r = false;
	if (&r != &a && &r != &b) {
		r.set_modulus(a);
		r.set_degree(db - 1);
	}
	else {
		truncate_r = true; // we have deg(a) >= deg(b) > deg(r)
		// relies on the fact that r is set after all other computations are
		// done !!!
	}
	rp = r.coeff;


	// now start doing the real work...
	const bigint& p = a.modulus();
#ifdef USE_SINGLE_PREC
	register udigit pp = p.least_significant_digit();
	bool sgl_prec = (p.compare(pp) == 0 && pp < max_udigit_modulus());
	if (sgl_prec)
		poldivrem_sgl(qp, rp, ap, da, bp, db, pp);
	else
#endif
		poldivrem(qp, rp, ap, da, bp, db, p);


	// tidy up if input aliases output...
	if (truncate_r)
		r.set_degree(db - 1);
	r.remove_leading_zeros();

	if (use_tmp)
		q.assign(tmp);
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// should re-write the following
void plain_rem(Fp_polynomial& r, const Fp_polynomial& a, const Fp_polynomial& b)
{
	debug_handler_l("Fp_polynomial", "plain_rem (Fp_polynomial&, Fp_polynomial&, Fp_polynomial&)", 5);

	Fp_polynomial tmp;
	plain_div_rem(tmp, r, a, b);
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void plain_div(Fp_polynomial &q, const Fp_polynomial &a, const Fp_polynomial &b)
{
	lidia_size_t deg_b = b.degree(), deg_a = a.degree();
	lidia_size_t deg_q = deg_a - deg_b;
	const bigint &p = a.modulus();
#ifdef USE_SINGLE_PREC
	register udigit pp = p.least_significant_digit();
	bool sgl_prec = (p.compare(pp) == 0 && pp < max_udigit_modulus());
	if (sgl_prec) {
		Fp_polynomial r;
		plain_div_rem(q, r, a, b);
		return;
	}
#endif

	if (deg_b < 0) {
		lidia_error_handler("Fp_polynomial", "plain_div (Fp_polynomial&, "
				    "Fp_polynomial&, Fp_polynomial&)::division by zero");
		return;
	}

	if (deg_q < 0) {
		q.set_modulus(a); //assigns zero
		return;
	}

	if (deg_q <= deg_b/2) {
		// variant 1 (multiplying with 'inverse'); faster if (deg_a/deg_b < 1.5)
		Fp_polynomial P1, P2, P3;
		copy_reverse(P3, b, 0, deg_b);
		invert(P2, P3, deg_q + 1);
		copy_reverse(P1, P2, 0, deg_q);

		shift_right(P3, a, deg_b);
		multiply(P2, P3, P1);
		trunc(P2, P2, 2*deg_q+1);
		shift_right(q, P2, deg_q);
	}// end variant 1
	else {
		// variant 2 (plain); cmp. plain_div_rem(...)
		lidia_size_t i, j;
		bool lc_is_one;
		bigint LCInv, t, s;
		const bigint *bp;

		if (b.lead_coeff().is_one())
			lc_is_one = true;
		else {
			lc_is_one = false;
			InvMod(LCInv, b.lead_coeff(), p);
		}

		Fp_polynomial lb;
		if (&q == &b) {
			lb = b;
			bp = lb.coeff;
		}
		else
			bp = b.coeff;

		bigint *xp = new bigint[deg_a+1];
		memory_handler(xp, "Fp_polynomial",
			       "plain_div :: Error in memory allocation");

		for (i = 0; i <= deg_a; i++)
			xp[i].assign(a.coeff[i]);

		q.set_modulus(a);
		q.set_degree(deg_q);
		bigint *qp = q.coeff;

		for (i = deg_q; i >= 0; i--) {
			Remainder(t, xp[i+deg_b], p);
			if (!lc_is_one)
				MulMod(t, t, LCInv, p);
			Remainder(qp[i], t, p);

			if (i > 0) {
				negate(t, t);
				for (j = deg_b - 1; j >= 0; j--) {
					multiply(s, t, bp[j]);
					add(xp[i+j], xp[i+j], s);
				}
			}
		}// end for
		delete[] xp;
	}// end variant 2
}



// ---------------------------------------------------------------------------

// x = (1/a) % X^m, input not output, constant term a is nonzero
void plain_inv(Fp_polynomial& x, const Fp_polynomial& a, lidia_size_t m)
{
	debug_handler("Fp_polynomial", "plain_inv(Fp_polynomial&, Fp_polynomial&, lidia_size_t)");

	lidia_size_t i, k, n, lb;
	bigint v, t, s;
	const bigint* ap;
	bigint* xp;

	n = a.degree();
	if (n < 0) {
		lidia_error_handler("Fp_polynomial", "plain_inv(Fp_polynomial&, "
				    "Fp_polynomial&, lidia_size_t)::division by zero");
		return;
	}

	x.set_modulus(a);

	const bigint &p = a.modulus();
	InvMod(s, a.const_term(), p);

	if (n == 0) {
		x.set_degree(0);
		x.coeff[0].assign(s);
		return;
	}

	ap = a.coeff;
	x.set_degree(m-1);
	xp = x.coeff;

	xp[0] = s;
	bool ct_is_one = s.is_one();

	for (k = 1; k < m; k++) {
		v.assign_zero();
		lb = comparator< lidia_size_t >::max(k-n, 0);
		for (i = lb; i <= k-1; i++) {
			multiply(t, xp[i], ap[k-i]);
			add(v, v, t);
		}
		Remainder(xp[k], v, p);
		NegateMod(xp[k], xp[k], p);
		if (!ct_is_one)
			MulMod(xp[k], xp[k], s, p);
	}
	x.remove_leading_zeros();
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
