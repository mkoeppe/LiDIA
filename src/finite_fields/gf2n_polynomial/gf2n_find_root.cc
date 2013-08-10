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
#include	"LiDIA/base_vector.h"
#include	"LiDIA/gf2n_polynomial.h"
#include	"LiDIA/gf2n_poly_modulus.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//*****************************************************************
//
// determine all roots of f in GF(2^d), append these roots to vector v.
// It is assumed that f splits completely and all the roots are in GF(2^d)
// Algorithm: see Lidl-Niederreiter book.
//
//****************************************************************


void find_all_roots(base_vector< gf2n > & v, const gf2n_polynomial & f,
		    unsigned int d)
{
	if (f.degree() == 0)
		return;

	gf2n_polynomial h(f);
	h.make_monic();

	if (h.degree() == 1) {
		lidia_size_t k = v.size();

		if (v.capacity() < k+1)
			v.set_capacity(k+1);
		v.set_size(k+1);

		v[k] = h.const_term();
		if (v[k].relative_degree() > d)
			lidia_error_handler("gf2n_polynomial", "find_all_roots::root not in"
					    " field");

#ifdef DEBUG
		assert(f(v[k]).is_zero());
#endif
		return;
	}

	if (h.degree() == 2) {
		gf2n x1, x2;

		do {
			h.get_coefficient(x2, 1);
		}
		while (!solve_quadratic(x1, x2, h.const_term()));

		divide(x2, h.const_term(), x1);

		if (x1.relative_degree() > d || x2.relative_degree() > d)
			lidia_error_handler("gf2n_polynomial", "find_all_roots::root not in"
					    " field");

		lidia_size_t k = v.size();

		if (v.capacity() < k+2)
			v.set_capacity(k+2);
		v.set_size(k+2);

		v[k] = x1;
		v[k+1] = x2;

#ifdef DEBUG
		assert(f(x1).is_zero());
		assert(f(x2).is_zero());
#endif

		return;
	}

	gf2n_polynomial h2;
	bool end = false;

	gf2n_poly_modulus F(h);
	gf2n_polynomial a, b, gamma;
	gf2n t;
	unsigned int i;

	gamma.assign_x();
	t.randomize(d);
	gamma.set_coefficient(t, 0);

	do {
		a.assign(gamma);
		b.assign(a);
		shift_left(gamma, gamma, F);

		for (i = 1; i < d; i++) {
			// trace computation
			square(b, b, F);
			add(a, a, b);
		}
		gcd(h2, a, h);

		if (!h2.is_one() && h2 != h && !h2.is_zero()) {
			end = true;
			find_all_roots(v, h2, d);
			find_all_roots(v, h / h2, d);
		}
	} while (!end);
}



//*********************************************************************
// function computes one root of f and returns this root.
// Algorithm: see Lidl-Niederreiter book.
//*********************************************************************

gf2n find_root(const gf2n_polynomial & f, unsigned int d)
{
	if (degree(f) == 0)
		lidia_error_handler("find_root", "constant polynomial");

	gf2n_polynomial h(f), h2;
	h.make_monic();

	if (h.degree() == 1)
		if (h.const_term().relative_degree() > static_cast<unsigned int>(d))
			lidia_error_handler("gf2n_polynomial", "find_root::root not in field");
		else
			return h.const_term();

	if (h.degree() == 2) {
		gf2n x1, x2;
		do {
			h.get_coefficient(x2, 1);
		} while (!solve_quadratic(x1, x2, h.const_term()));
		if (x1.relative_degree() > static_cast<unsigned int>(d))
			lidia_error_handler("gf2n_polynomial", "find_root::root not in field");
		return x1;
	}

	gf2n_poly_modulus F(h);
	gf2n_polynomial a, b, gamma;
	gf2n t;
	unsigned int i;

	gamma.assign_x();
	t.randomize(d);
	gamma.set_coefficient(t, 0);

	do {
		a.assign(gamma);
		b.assign(a);
		shift_left(gamma, gamma, F);

		for (i = 1; i < d; i++) {
			square(b, b, F);
			add(a, a, b);
		}
		gcd(h2, a , h);

		if (!h2.is_one() && h2 != h && !h2.is_zero()) {
			// splitting found
			if (h2.degree() < h.degree() / 2)
				return find_root(h2, d);
			else
				return find_root(h / h2, d);
		}
	} while (true);
	return gf2n();
}



//
// Task:        computes w = a+a^q+...+^{q^{d-1}} mod f    (f = F.modulus())
//
// Conditions:  d >= 0,
//              b = x^q mod f,
//              q a power of p
//
// Algorithm:   trace map
//
//      The main idea is the following:
//      Let     w_m = sum_{0<=i<m} a^{q^i}
//      and     z_m = x^{q^m}.
//
//      Then    w_{2m} = w_m + w_m(z_m)
//      and     z_{2m} = z_m(z_m)
//

void trace_map(gf2n_polynomial & w, const gf2n_polynomial& a,
	       lidia_size_t d, gf2n_poly_modulus & F,
	       const gf2n_polynomial & b)
{
	gf2n_polynomial y, z, t;

	z.assign(b);
	y.assign(a);
	w.assign_zero();

	// Let d = sum_{0<=i<=l} d_i*2^i.
	// Set i=0.

	while (d) {
		if (d == 1) {
			if (w.is_zero())
				w.assign(y); // w = w(z) + y = y
			else {
				compose(w, w, z, F);
				add(w, w, y); // w = w(z) + y
			}
			// we don't need to compute new values for y and z
			// because we're done
		}
		else {
			if ((d & 1) == 0) {
				compose(t, y, z, F);
				compose(z, z, z, F); // z = z(z)
				add(y, t, y); // y = y(z) + y
			}
			else {
				if (w.is_zero()) {
					w.assign(y); // w = w(z) + y = y
					compose(t, y, z, F); // y = y(z) + y
					compose(z, z, z, F); // z = z(z)
					add(y, t, y);
				}
				else {
					compose(w, w, z, F); // w = w(z) + y = y
					compose(t, y, z, F); // y = y(z) + y
					compose (z, z, z, F); // z = z(z)
					add(w, w, y);
					add(y, t, y);
				}
			}
		}

		// Now, we have checked bit i of (the original) d.
		//      z = x^{q^{2^i}}
		//      y = sum_{1<=j<=2^i} a^{q^j}
		//      w = sum_{0<=j<=D_i} a^{q^j}, where D_i = sum_{0<=j<=i} d_j*2^j
		// Set i = i+1

		d = d >> 1;
	}
}



// the polynomial stored in F is a equal degree polynomial of degree d.
// The procedure tries to find one irreducible factor of degree d in F.

void EDF_one_factor(gf2n_polynomial & w, const gf2n_polynomial & Xq,
		    gf2n_poly_modulus & Fpm, int d)
{
	int Fd = Fpm.modulus().degree();

	if (Fd < d)
		lidia_error_handler("gf2n_polynomial", "EDF_one_factor::degree of input"
				    " smaller than degree of factor");

	if (Fd == d) {
		w.assign(Fpm.modulus());
		return;
	}

	if (d == 1) {
		gf2n t = find_root(Fpm.modulus());
		w = gf2n_polynomial(1);
		w.set_coefficient(t, 0);
		return;
	}

	gf2n_polynomial a, xq, b;
	gf2n_poly_modulus F(Fpm);
	int degree_field = gf2n::get_absolute_degree();

	if (Xq.degree() >= Fd)
		remainder(xq, Xq, F.modulus());
	else
		xq.assign(Xq);

	do {
		a.randomize(d-1);
		trace_map(b, a, d, F, xq);
		w.assign(b);
		for (int i = 1; i < degree_field; i++) {
			square(b, b, F);
			add(w, w, b);
		}

		gcd(w, w, F.modulus());
		if (!w.is_one() && w.degree() < Fd) {
			if (Fd - w.degree() < w.degree())
				divide(w, F.modulus(), w);

			if (w.degree() > d) {
				F.build(w);
				remainder(xq, xq, F.modulus());
				Fd = w.degree();
			}
		}
	} while (w.degree() != d);

#ifdef DEBUG
	assert((Fpm.modulus() % w).is_zero());
#endif
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
