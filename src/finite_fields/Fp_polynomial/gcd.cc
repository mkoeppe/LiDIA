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
#include	"LiDIA/finite_fields/Fp_polynomial_util.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//***************************************************************
//	This File contains the implementation of the functions
//	- gcd, xgcd
//	- plain_gcd, plain_xgcd (both always use "classical" algorithms)
//	- and some auxiliary functions
//	- resultant
//***************************************************************


//***************************************************************
//				    auxiliary functions
//	all assume that all input parameters have correct modulus
//***************************************************************

#if 0
// TPf: speedups in plain_rem make plain_rem_vec obsolete...

void plain_rem_vec(Fp_polynomial& r, const Fp_polynomial& a, const Fp_polynomial& b, bigint *x)
	//special version of plain_rem, uses memory offered by vector x
	//about 10% faster
	//assumes x.size() > a.degree()+1
	//cmp. plain_div_rem_vec in file poly_matrix.c
{
	debug_handler("Fp_polynomial", "plain_rem_vec (Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, bigint*)");

	lidia_size_t da, db, dq, i, j;
	bigint LCInv, t, s;

	da = a.degree();
	db = b.degree();
	dq = da - db;

	if (db < 0) {
		lidia_error_handler("Fp_polynomial", "plain_rem_vec(...)::"
				    "division by zero");
		return;
	}

	if (dq < 0) {
		r.assign(a);
		return;
	}

	bool lc_is_one;
	const bigint &p = a.modulus();
	const bigint *bp = b.coeff;
	if (bp[db].is_one())
		lc_is_one = true;
	else {
		lc_is_one = false;
		InvMod(LCInv, bp[db], p);
	}

	for (i = 0; i <= da; i++)
		x[i].assign(a.coeff[i]);

	for (i = dq; i >= 0; i--) {
		Remainder(t, x[i+db], p);
		if (!lc_is_one)
			MulMod(t, t, LCInv, p);
		negate(t, t);

		for (j = db-1; j >= 0; j--) {
			multiply(s, t, bp[j]);
			add(x[i+j], x[i+j], s);
		}
	}

	r.set_modulus(a);
	r.set_degree(db-1);
	for (i = 0; i < db; i++)
		Remainder(r.coeff[i], x[i], p);
	r.remove_leading_zeros();
}
#endif



static
void half_gcd(Fp_polynomial& U, Fp_polynomial& V)
{
	debug_handler("gcd.c", "half_gcd(Fp_polynomial&, Fp_polynomial&)");

	lidia_size_t du = U.degree();
	lidia_size_t d_red = (du+1)/2;
	if (V.is_zero() || V.degree() <= du - d_red)
		return;

	lidia_size_t d1 = (d_red + 1)/2;
	if (d1 < 1)
		d1 = 1;
	if (d1 >= d_red)
		d1 = d_red - 1;

	poly_matrix M1(V.modulus());
	M1.half_gcd(U, V, d1);
	M1.multiply(U, V);

	lidia_size_t d2 = V.degree() - du + d_red;

	if (V.is_zero() || d2 <= 0)
		return;

	M1.kill();

	Fp_polynomial Q;
	div_rem(Q, U, U, V);
	swap(U, V);

	poly_matrix M2(V.modulus());
	M2.half_gcd(U, V, d2);
	M2.multiply(U, V);
}



//***************************************************************
//
//			    gcd
//
//***************************************************************

void gcd(Fp_polynomial& d, const Fp_polynomial& u, const Fp_polynomial& v)
{
	debug_handler("Fp_polynomial", "gcd(Fp_polynomial&, Fp_polynomial&, Fp_polynomial&)");

	u.comp_modulus(v, "gcd");

	Fp_polynomial u1(u);
	Fp_polynomial v1(v);


	if (u1.degree() == v1.degree()) {
		if (u1.is_zero()) {
			d.set_modulus(v); //assigns zero
			return;
		}
		remainder(v1, v1, u1);
	}
	else
		if (u1.degree() < v1.degree())
			swap(u1, v1);


	// deg(u1) > deg(v1)
	lidia_size_t crov =
		Fp_polynomial::crossovers.gcd_crossover(u.modulus());
	while (u1.degree() >= crov && !v1.is_zero()) {
		half_gcd(u1, v1);
		if (!v1.is_zero()) {
			remainder(u1, u1, v1);
			swap(u1, v1);
		}
	}

	plain_gcd(d, u1, v1);
}



// should re-write the following
void plain_gcd(Fp_polynomial& x, const Fp_polynomial& a, const Fp_polynomial& b)
{
	debug_handler("Fp_polynomial", "plain_gcd(Fp_polynomial&, Fp_polynomial&, Fp_polynomial&");

	a.comp_modulus(b, "plain_gcd");

	if (b.is_zero())
		x.assign(a);
	else
		if (a.is_zero())
			x.assign(b);
		else {
			Fp_polynomial u(a);
			Fp_polynomial v(b);

			do
			{
				plain_rem(u, u, v);
				// TPf: original code:     plain_rem_vec(u, u, v, TEMP);
				//      (speedups in plain_rem make plain_rem_vec obsolete)
				swap(u, v);
			} while (!v.is_zero());

			x.assign(u);
		}

	if (x.is_zero())
		return;

	// make gcd monic
	x.make_monic();
}



//***************************************************************
//
//			    xgcd
//
//***************************************************************

void xgcd(Fp_polynomial& d, Fp_polynomial& s, Fp_polynomial& t, const Fp_polynomial& a, const Fp_polynomial& b)
{
	debug_handler("Fp_polynomial", "xgcd(Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_polynomial&)");

	a.comp_modulus(b, "xgcd");

	if (a.is_zero() && b.is_zero()) {
		d.set_modulus(a); //assigns zero
		s.set_modulus(a); //assigns zero
		t.set_modulus(a); //assigns zero
		return;
	}

	if (comparator< lidia_size_t >::max(a.degree(), b.degree()) <
	    Fp_polynomial::crossovers.xgcd_crossover(a.modulus())) {
		plain_xgcd(d, s, t, a, b);
		return;
	}

	Fp_polynomial U(a), V(b), Q;
	lidia_size_t flag = 0;

	if (U.degree() == V.degree()) {
		div_rem(Q, U, U, V);
		swap(U, V);
		flag = 1;
	}
	else
		if (U.degree() < V.degree()) {
			swap(U, V);
			flag = 2;
		}

	const bigint &p = a.modulus();
	poly_matrix M(p);

	M.xhalf_gcd(U, V, U.degree()+1);

	d = U;

	if (flag == 0) {
		s = M(0, 0);
		t = M(0, 1);
	}
	else
		if (flag == 1) {
			s = M(0, 1);
			multiply(t, Q, M(0, 1));
			subtract(t, M(0, 0), t);
		}
		else
		{  // flag == 2
			s = M(0, 1);
			t = M(0, 0);
		}

	if (d.is_monic())
		return;

	// normalize

	bigint w; //static
	InvMod(w, d.lead_coeff(), p);
	multiply_by_scalar(d, d, w);
	multiply_by_scalar(s, s, w);
	multiply_by_scalar(t, t, w);
}



void xgcd_left(Fp_polynomial& d, Fp_polynomial& s,
               const Fp_polynomial& a, const Fp_polynomial& b)
{
	debug_handler("Fp_polynomial", "xgcd_left(Fp_polynomial&, "
		      "Fp_polynomial&, Fp_polynomial&, Fp_polynomial&)");

	a.comp_modulus(b, "xgcd_left");

	if (a.is_zero() && b.is_zero()) {
		d.set_modulus(a); //assigns zero
		s.set_modulus(a); //assigns zero
		return;
	}

	if (comparator< lidia_size_t >::max(a.degree(), b.degree()) <
	    Fp_polynomial::crossovers.xgcd_crossover(a.modulus())) {
		plain_xgcd_left(d, s, a, b);
		return;
	}

	Fp_polynomial U(a), V(b), Q;
	lidia_size_t flag = 0;

	if (U.degree() == V.degree()) {
		div_rem(Q, U, U, V);
		swap(U, V);
		flag = 1;
	}
	else
		if (U.degree() < V.degree()) {
			swap(U, V);
			flag = 2;
		}

	const bigint &p = a.modulus();
	poly_matrix M(p);

	M.xhalf_gcd(U, V, U.degree()+1);

	d = U;

	if (flag == 0)
		s = M(0, 0);
	else
		if (flag == 1)
			s = M(0, 1);
		else        // flag == 2
			s = M(0, 1);

	if (d.is_monic())
		return;

	// normalize

	bigint w; //static
	InvMod(w, d.lead_coeff(), p);
	multiply_by_scalar(d, d, w);
	multiply_by_scalar(s, s, w);
}



// should re-write the following
void plain_xgcd(Fp_polynomial& d, Fp_polynomial& s, Fp_polynomial& t, const Fp_polynomial& a, const Fp_polynomial& b)
{
	debug_handler("Fp_polynomial", "plain_xgcd(Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_polynomial&");

	a.comp_modulus(b, "plain_xgcd");

	const bigint &p = a.modulus();
	if (b.is_zero()) {
		t.set_modulus(a); //assigns zero
		s.set_modulus(a);
		s.assign_one();
		d.assign(a);
	}
	else if (a.is_zero()) {
		s.set_modulus(a); //assigns zero
		t.set_modulus(a);
		t.assign_one();
		d.assign(b);
	}
	else {
		lidia_size_t e = comparator< lidia_size_t >::max(a.degree(), b.degree());

		Fp_polynomial u(a);
		Fp_polynomial v(b);
		Fp_polynomial temp, u0, v0, u1, v1, u2, v2, q, r;

		temp.set_max_degree(e);
		u0.set_max_degree(e);
		v0.set_max_degree(e);
		u1.set_max_degree(e);
		v1.set_max_degree(e);
		u2.set_max_degree(e);
		v2.set_max_degree(e);
		q.set_max_degree(e);
		r.set_max_degree(e);

		u1.set_modulus(a);
		u2.set_modulus(a);
		v1.set_modulus(a);
		v2.set_modulus(a);

		u1.assign_one(); v1.assign_zero();
		u2.assign_zero(); v2.assign_one();

		do
		{
			div_rem(q, r, u, v);
			u.assign(v); //u = v;
			v.assign(r); //v = r;
			u0.assign(u2); //u0 = u2;
			v0.assign(v2); //v0 = v2;
			multiply(temp, q, u2);
			subtract(u2, u1, temp);
			multiply(temp, q, v2);
			subtract(v2, v1, temp);
			u1.assign(u0); //u1 = u0;
			v1.assign(v0); //v1 = v0;
		} while (!v.is_zero());

		d.assign(u); //d = u;
		s.assign(u1); //s = u1;
		t.assign(v1); //t = v1;
	}

	if (d.is_zero() || d.is_monic())
		return;

	// make gcd monic
	bigint z; //static

	InvMod(z, d.lead_coeff(), p);
	multiply_by_scalar(d, d, z);
	multiply_by_scalar(s, s, z);
	multiply_by_scalar(t, t, z);
}



void plain_xgcd_left(Fp_polynomial& d, Fp_polynomial& s,
                     const Fp_polynomial& a, const Fp_polynomial& b)
{
	debug_handler("Fp_polynomial", "plain_xgcd_left(Fp_polynomial&, "
		      "Fp_polynomial&, Fp_polynomial&, "
		      "Fp_polynomial&");

	a.comp_modulus(b, "plain_xgcd_left");

	const bigint &p = a.modulus();
	if (b.is_zero()) {
		s.set_modulus(a);
		s.assign_one();
		d.assign(a);
	}
	else
		if (a.is_zero()) {
			s.set_modulus(a); //assigns zero
			d.assign(b);
		}
		else {
			lidia_size_t e = comparator< lidia_size_t >::max(a.degree(),
									  b.degree());

			Fp_polynomial u(a);
			Fp_polynomial v(b);
			Fp_polynomial temp, u0, v0, u1, v1, u2, v2, q, r;

			temp.set_max_degree(e);
			u0.set_max_degree(e);
			v0.set_max_degree(e);
			u1.set_max_degree(e);
			v1.set_max_degree(e);
			u2.set_max_degree(e);
			v2.set_max_degree(e);
			q.set_max_degree(e);
			r.set_max_degree(e);

			u1.set_modulus(a);
			u2.set_modulus(a);
			v1.set_modulus(a);
			v2.set_modulus(a);

			u1.assign_one(); v1.assign_zero();
			u2.assign_zero(); v2.assign_one();

			do
			{
				plain_div_rem(q, r, u, v);
				u.assign(v); //u = v;
				v.assign(r); //v = r;
				u0.assign(u2); //u0 = u2;
				multiply(temp, q, u2);
				subtract(u2, u1, temp);
				u1.assign(u0); //u1 = u0;
			}
			while (!v.is_zero());

			d.assign(u); //d = u;
			s.assign(u1); //s = u1;
		}

	if (d.is_zero() || d.is_monic())
		return;

	// make gcd monic
	bigint z; //static

	InvMod(z, d.lead_coeff(), p);
	multiply_by_scalar(d, d, z);
	multiply_by_scalar(s, s, z);
}



void resultant (bigint &r, const Fp_polynomial & ff, const Fp_polynomial & gg)
{
	Fp_polynomial f(ff), g(gg);
	bigint l, lead;
	const bigint & p = f.modulus();
	int m, d;

	r.assign_one();
	m = g.degree();

	do {
		remainder(g, g, f);
		d = g.degree();

		if (d < m) {
			l.assign(f.lead_coeff());
			power_mod(l, l, m-d, p);
			MulMod(r, r, l, p);
		}

		if (g.degree() <= 0)
			break;

		swap(f, g);

		m = g.degree();
		if ((f.degree() & 1 == 1) && (m & 1) == 1)
			subtract(r, p, r);
	} while (true);

	if (g.is_zero())
		r.assign_zero();
	else
		MulMod(r, r, g.const_term(), p);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
