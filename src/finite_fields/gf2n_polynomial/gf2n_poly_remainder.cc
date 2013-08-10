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
//	Author	: Volker Mueller
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/gf2n_polynomial.h"
#include        "LiDIA/arith.inl"
#include	"crossover.tbl"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//=====================================================================
// first the standard function for computing the remainder

void plain_rem(gf2n_polynomial & c, const gf2n_polynomial & a,
	       const gf2n_polynomial & bb)
{
	if (bb.is_zero())
		lidia_error_handler("gf2n_polynomial", "remainder::divisor is zero");

	gf2n_polynomial b (bb); // if c aliases bb -->problem without copying
	c.assign(a);

	if (a.deg < b.deg)
		return;

	gf2n *cp, *bp, t, h;
	const int bd = b.degree();
	int i, j;
	gf2n blc = inverse(b.lead_coeff());

	for (i = a.deg; i >= bd; i--) {
		cp = c.coeff + i;

		if (cp->is_zero())
			continue;

		multiply(t, *cp, blc);
		cp->assign_zero();

		for (j = bd-1, bp = b.coeff + bd - 1, cp --;
		     j >= 0; j--, cp--, bp--) {
			multiply(h, *bp, t);
			add(*cp, *cp, h);
		}
	}

	c.deg = bd-1;
	cp = c.coeff + c.deg;

	while (cp->is_zero() && c.deg > 0) {
		cp --;
		c.deg --;
	}

	if (c.deg == 0 && c.coeff[0].is_zero())
		c.deg = -1;
}


//==============================================================
// several functions needed for kara_remainder
//
// copy_reverse :
// x[0..hi-lo+1] = reverse(a[lo..hi]), with zero fill
// input may not alias output

void copy_reverse(gf2n_polynomial & x, const gf2n_polynomial & a,
		  lidia_size_t lo, lidia_size_t hi)
{
	debug_handler("gf2n_poly_remainder.c", "copy_reverse(gf2n_polynomial &, "
		      "gf2n_polynomial&, lidia_size_t, lidia_size_t) ");

	lidia_size_t i, j, n, m;
	n = hi - lo + 1;
	m = a.degree() + 1; // = a.c_length

	if (x.size < n)
		x.set_size(n);

	x.deg = n - 1;

	gf2n* xp = x.coeff;
	const gf2n* ap = a.coeff;

	for (i = 0; i < n; i++) {
		j = hi - i;
		if (j< 0 || j >= m)
			xp[i].assign_zero();
		else
			xp[i].assign(ap[j]);
	}
}




//=======================================================================
// the invert function just chooses the better of two implmented algorithms
//

void invert(gf2n_polynomial & x, const gf2n_polynomial & a, lidia_size_t m)
{
	debug_handler_l("gf2n_polynomial", "inv(gf2n_polynomial&, "
			"gf2n_polynomial&, lidia_size_t)", 5);



	if (&x == &a) {
		gf2n_polynomial la(a);

		if (a.degree() < CROSSOVER_PLAININV_NEWTON)
			plain_inv(x, la, m);
		else
			newton_inv(x, la, m);
	}
	else
		if (a.degree() < CROSSOVER_PLAININV_NEWTON)
			plain_inv(x, a, m);
		else
			newton_inv(x, a, m);
}

//====================================================================
// the standard inversion function
// x = (1/a) % X^m, input not output, constant term a is nonzero

void plain_inv(gf2n_polynomial& x, const gf2n_polynomial& a, lidia_size_t m)
{
	debug_handler("gf2n_polynomial", "plain_inv(gf2n_polynomial&, "
		      "gf2n_polynomial&, lidia_size_t)");

	lidia_size_t i, k, n, lb;
	gf2n v, t, s;
	const gf2n* ap;
	gf2n* xp;

	n = a.degree();
	if (n < 0) {
		lidia_error_handler("gf2n_polynomial", "plain_inv(gf2n_polynomial&, "
				    "gf2n_polynomial&, lidia_size_t)::division by "
				    "zero");
		return;
	}

	if (a.const_term().is_zero()) {
		lidia_error_handler("gf2n_polynomial", "plain_inv(gf2n_polynomial&, "
				    "gf2n_polynomial&, lidia_size_t)::const term "
				    "is zero");
		return;
	}

	if (x.size < m)
		x.set_size(m);

	invert(s, a.const_term());

	if (n == 0) {
		x.deg = 0;
		x.coeff[0].assign(s);
		return;
	}

	ap = a.coeff;
	x.deg = m-1;
	xp = x.coeff;

	xp[0].assign(s);
	bool ct_is_one = s.is_one();

	for (k = 1; k < m; k++) {
		v.assign_zero();
		lb = comparator< lidia_size_t >::max(k-n, 0);
		for (i = lb; i <= k-1; i++) {
			multiply(t, xp[i], ap[k-i]);
			add(v, v, t);
		}
		if (!ct_is_one)
			multiply(xp[k], v, s);
		else
			xp[k].assign(v);
	}
}

//===============================================================
// x = (1/a) % X^m, input not output, constant term a is nonzero

void newton_inv(gf2n_polynomial & x, const gf2n_polynomial& a, lidia_size_t m)
{
	debug_handler("gf2n_poly_remainder.c", "newton_inv(gf2n_polynomial&, "
		      "gf2n_polynomial&, lidia_size_t)");

	lidia_size_t l = 1, s = 1 << (static_cast<int>(LiDIA::log2(static_cast<double>(m)) + 1));
	gf2n_polynomial t, h;

	h.set_size(s);
	x.set_size(s);
	t.set_size(s);

	if (a.const_term().is_zero()) {
		lidia_error_handler("gf2n_poly_remainder.c", "newton_invert::constant"
				    "term is zero");
		return;
	}
	else {
		invert(x.coeff[0], a.const_term());
		x.deg = 0;
	}

	while (l < m) {
		square(t, x);
		trunc(h, a, l << 1);
		trunc(t, t, l << 1);
		multiply(t, t, h);
		trunc(t, t, l << 1);

		for (int i = l; i < l << 1; i++)
			x.coeff[i].assign(t.coeff[i]);

		x.deg = (l << 1) - 1;

		while (x.coeff[x.deg].is_zero() && x.deg >= 0)
			x.deg --;

		l <<= 1;
	}
	x.deg = m - 1;
	while (x.coeff[x.deg].is_zero() && x.deg >= 0)
		x.deg --;
}


//=================================================================
// shift right and delete lower coefficients

void floor (gf2n_polynomial & g, const gf2n_polynomial & f,
	    unsigned int d)

{
	register int i, fd = f.deg;
	gf2n *ap, *bp;

	if (d == 0 || f.is_zero()) {
		g.assign(f);
		return;
	}

	if (g.size <= static_cast<int>(fd-d))
		g.set_size (fd-d+1);

	for (i = d, ap = g.coeff, bp = f.coeff+d; i <= fd; i++, ap++, bp++)
		(*ap).assign(*bp);

	g.deg = fd-d;
}




//==================================================================
// the main reduction function

void kara_div(gf2n_polynomial & q, const gf2n_polynomial & a,
	      const gf2n_polynomial &b)
{
	lidia_size_t deg_b = b.degree(), deg_a = a.degree();
	lidia_size_t deg_q = deg_a - deg_b;

	if (deg_b < 0) {
		lidia_error_handler("gf2n_polynomial", "plain_div (gf2n_polynomial &, "
				    "gf2n_polynomial&, gf2n_polynomial&)::division "
				    "by zero");
		return;
	}

	if (deg_q < 0) {
		q.assign_zero();
		return;
	}

	gf2n_polynomial P1(deg_q+1), P2(deg_q+1), P3(deg_b+1);

	copy_reverse(P3, b, 0, deg_b);
	invert(P2, P3, deg_q + 1);
	copy_reverse(P1, P2, 0, deg_q);

	floor(P3, a, deg_b);
	multiply(P2, P3, P1);
	trunc(P2, P2, 2*deg_q+1);
	floor(q, P2, deg_q);
}

void kara_rem(gf2n_polynomial & r, const gf2n_polynomial & a,
	      const gf2n_polynomial & b)
{
	if (a.deg < b.deg) {
		r.assign(a);
		return;
	}

	gf2n_polynomial q(a.deg-b.deg);

	kara_div(q, a, b);
	multiply(q, q, b);
	add(r, a, q);
}


void remainder (gf2n_polynomial &r, const gf2n_polynomial &a,
		const gf2n_polynomial &b)
{
	if (b.is_zero()) {
		lidia_error_handler("gf2n_polynomial", "remainder::modulus is zero");
		return;
	}

	if (b.degree() == 0) {
		r.assign_zero();
		return;
	}


	if (b.degree() >= CROSSOVER_REDUCE)
		kara_rem(r, a, b);
	else
		plain_rem(r, a, b);
}


void divide (gf2n_polynomial &r, const gf2n_polynomial &a,
	     const gf2n_polynomial &b)
{
	if (b.is_zero()) {
		lidia_error_handler("gf2n_polynomial", "divide::modulus is zero");
		return;
	}

	if (b.degree() == 0) {
		r.assign_zero();
		return;
	}


	if (b.degree() >= CROSSOVER_DIVIDE)
		kara_div(r, a, b);
	else
		plain_div(r, a, b);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
