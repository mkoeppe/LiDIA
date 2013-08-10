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



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



int gf2n_polynomial::default_size = 4;
int gf2n_polynomial::refc = 0;



void
gf2n_polynomial::set_size(unsigned int d)
{
	if (size > -1)
		delete[] coeff;
	if (d >= static_cast<unsigned int>(default_size))
		size = d+1;
	else
		size = default_size;

	coeff = new gf2n[size];
}



gf2n_polynomial::gf2n_polynomial()                 // = 0
{
	deg = -1;
	refc++;
	size = -1;
	coeff = NULL;
}



gf2n_polynomial::gf2n_polynomial(unsigned int d)   // = 1 * X^d
{
	register int i;
	refc++;
	deg = d;
	if (d >= static_cast<unsigned int>(default_size))
		size = d+1;
	else
		size = default_size;

	coeff = new gf2n[size];
	for (i = 0; i < static_cast<int>(d); i++)
		coeff[i].assign_zero();
	coeff[d].assign_one();
}



gf2n_polynomial::gf2n_polynomial(const gf2n & f)   // = f * X^0
{
	deg = 0;
	refc++;
	size = default_size;
	coeff = new gf2n[size];
	coeff[0].assign(f);
}



gf2n_polynomial::gf2n_polynomial(const gf2n_polynomial & p)  // = p
{
	gf2n *ap, *cp;
	register int i;

	refc++;
	deg = p.deg;
	if (p.size > default_size)
		size = p.size;
	else
		size = default_size;
	coeff = new gf2n[size];

	for (i = deg+1, ap = coeff, cp = p.coeff; i; i--, ap++, cp++)
		*ap = *cp;
}



//****  comparisons of polynomials ******************************


bool operator == (const gf2n_polynomial &a,
		  const gf2n_polynomial &b)
{
	gf2n *ap, *bp;
	register int i, da = a.deg, db = b.deg;

	if (da != db)
		return false;
	for (i = da + 1, ap = a.coeff, bp = b.coeff; i; i--, ap++, bp++)
		if (*ap != *bp)
			return false;
	return true;
}



void gf2n_polynomial::make_monic()
{
	if (deg < 0)
		lidia_error_handler("gf2n_polynomial",
				    "make_modic: zero gf2n_polynomial");

	if (coeff[deg].is_one())
		return;

	gf2n s = inverse(coeff[deg]);
	for (register int i = 0; i < deg; i++)
		multiply(coeff[i], s, coeff[i]);
	coeff[deg].assign_one();
}



unsigned int gf2n_polynomial::degree_of_definition() const
{
	if (is_zero())
		return 1;

	unsigned int d = 1;

	for (int i = 0; i <= deg; i++) {
		d = lcm(d, static_cast<udigit>(coeff[i].relative_degree()));
		if (d == gf2n::extension_degree())
			return d;
	}
	return d;
}



bool gf2n_polynomial::is_square() const
{
	if (deg & 1)
		return false;
	for (int i = 1; i < deg; i += 2)
		if (!(coeff[i]).is_zero())
			return false;
	return true;
}



void gf2n_polynomial::assign(const gf2n_polynomial & f)
{
	register int i;
	gf2n * fp, *gp;

	if (static_cast<int>(size) < f.deg+1)
		set_size(f.deg + 1);
	deg = f.deg;
	for (i = 0, gp = coeff, fp = f.coeff; i <= deg; i++, fp++, gp++)
		(*gp).assign(*fp);
}



void gf2n_polynomial::assign (const base_vector< gf2n > & a)
{
	lidia_size_t j;
	deg = a.size() - 1;
	set_size(deg);

	for (j = 0; j <= deg; j++)
		coeff[j].assign(a[j]);

	while (coeff[deg].is_zero() && deg > 0)
		deg --;

	if (deg == 0 && coeff[0].is_zero())
		deg = -1;
}



void gf2n_polynomial::randomize(unsigned int degree)
{
	int i;

	set_size(degree+1);
	deg = degree;

	for (i = 0; i <= deg; i++)
		coeff[i].randomize();

	while (coeff[deg].is_zero())
		coeff[deg].randomize();
}



void gf2n_polynomial::set_coefficient(const gf2n & x, unsigned int i)
{
	if (static_cast<int>(i) > deg) {
		if (static_cast<int>(i) > size)
			set_size(i+1);
		deg = i;
	}

	coeff[i].assign(x);

	if ((static_cast<int>(i) == deg) && (x.is_zero())) {
		while (coeff[deg].is_zero() && deg > 0)
			deg --;
		if (deg == 0 && coeff[0].is_zero())
			deg = -1;
	}
}



void gf2n_polynomial::set_coefficient(unsigned int i)
{
	if (static_cast<int>(i) > deg) {
		if (static_cast<int>(i) > size)
			set_size(i+1);
		deg = i;
	}

	coeff[i].assign_one();
}



void gf2n_polynomial::get_coefficient(gf2n & x, unsigned int i) const
{
	if (static_cast<int>(i) > deg) {
		lidia_error_handler("gf2n_polynomial",
				    "get_coefficient::index out of range");
		x.assign_zero();
	}
	else
		x.assign(coeff[i]);
}



void multiply_by_scalar(gf2n_polynomial & g, const gf2n & b,
			const gf2n_polynomial & a)
{
	gf2n *ap, *cp;
	register int i;

	if (b.is_zero()) {
		g.deg = -1;
		return;
	}

	if (b.is_one()) {
		g.assign(a);
		return;
	}

	register int deg_a = a.deg;

	if (deg_a >= static_cast<int>(g.size))
		g.set_size (deg_a + 1 > static_cast<int>(gf2n_polynomial::default_size) ?
			    deg_a+ 1 : gf2n_polynomial::default_size);
	g.deg = deg_a;

	for (i = 0, ap = a.coeff, cp = g.coeff; i <= deg_a;
	     i++, ap++, cp++)
		multiply(*cp, b, *ap);
}



//----------------------------------------------------------------------
// compute and return the value f(x)

gf2n eval(const gf2n_polynomial & f, const gf2n & x)
{
	int i, fd = f.deg;
	gf2n t, *ap;

	if (f.is_zero()) {
		t.assign_zero();
		return (t);
	}

	t.assign(f.coeff[fd]);

	if (fd == 0)
		return t;

	for (i = fd-1, ap = f.coeff+fd-1; i >= 0; i--, ap--) {
		multiply(t, t, x);
		add(t, t, *ap);
	}
	return (t);
}



void swap (gf2n_polynomial & f, gf2n_polynomial & g)
{
	gf2n_polynomial h(g);
	g.assign(f);
	f.assign(h);
}



//*****************************************************************
// compute f = sqrt(g), if g is a square and return TRUE, otherwise
// return false and f is undefined.

bool sqrt(gf2n_polynomial & f, const gf2n_polynomial & g)
{
	int i, d = g.deg >> 1;

	if (!g.is_square())
		return false;

	if (f.size < d+1)
		f.set_size(d+1);

	f.deg = d;
	for (i = 0; i <= d; i++)
		sqrt(f.coeff[i], g.coeff[i << 1]);

	return true;
}



//*****************************************************************
// set g = f * X^d,   g, f may be aliases

void shift_left (gf2n_polynomial & g, const gf2n_polynomial & f,
		 unsigned int d)

{
	gf2n *cp, *ap;

	if (d == 0 || f.is_zero()) {
		g.assign(f);
		return;
	}

	register int dg = f.deg + d, i, fd = f.deg;

	if (dg >= g.size) {
		g.size = dg+1;
		cp = new gf2n[g.size];
		cp = cp + dg; // f, g may be aliases in this case -->no free
	}
	else
		cp = g.coeff + dg;

	g.deg = dg;

	for (i = fd+1, ap = f.coeff+fd; i; i--, cp--, ap--)
		(*cp).assign(*ap);

	for (i = d; i; i--, cp--)
		(*cp).assign_zero();

	cp ++; // cp should point to the beginning of array now
	if (cp != g.coeff) {
		delete[] g.coeff;
		g.coeff = cp;
	}
}



//*********************************************************************
//   g = f / X^d,   g, f may be the same gf2n_polynomial, if f != 0 mod X^d,
//                  the error handler is invoked

void shift_right (gf2n_polynomial & g, const gf2n_polynomial & f,
		  unsigned int d)

{
	register int i, fd = f.deg;
	gf2n *ap, *bp;

	if (d == 0 || f.is_zero()) {
		g.assign(f);
		return;
	}

	for (i = 0, ap = f.coeff; i < static_cast<int>(d); i++, ap++)
		if (!(*ap).is_zero())
			lidia_error_handler("gf2n_polynomial", "divide_by_x: X^d no divisor");

	if (g.size <= static_cast<int>(fd-d))
		g.set_size (fd-d+1);

	for (i = d, ap = g.coeff, bp = f.coeff+d; i <= fd; i++, ap++, bp++)
		(*ap).assign(*bp);

	g.deg = fd-d;
}



//*********************************************************************
// g = f * (X + a)  ,  f,g may be aliases

void multiply_by_linear(gf2n_polynomial & g, const gf2n_polynomial & f,
			const gf2n & a)
{
	gf2n *ap, *cp, t;

	if (f.is_zero()) {
		g.assign(f);
		return;
	}

	if (a.is_zero()) {
		shift_left(g, f, 1);
		return;
	}

	register int dg = f.deg+1, i, fd = f.deg;

	if (dg >= static_cast<int>(g.size)) {
		g.size = dg+1;
		cp = new gf2n[g.size];
		cp = cp + dg; // now free, since f, g may be aliases
	}
	else
		cp = g.coeff + dg;

	g.deg = dg;

	ap = f.coeff+fd; // highest coefficient
	*cp = *ap;
	cp --;

	for (i = fd; i; i--, cp--) {
		multiply(t, a, *ap);
		ap --;
		add(*cp, t, *ap);
	}

	multiply(*cp, a, *ap);

	if (cp != g.coeff) {
		delete[] g.coeff;
		g.coeff = cp;
	}
}



//*********** operations:  add *****************************************


void add(gf2n_polynomial &c, const gf2n_polynomial & a,
	 const gf2n_polynomial & b)
{
	register int deg_a = a.deg, deg_b = b.deg;
	register int i, min_deg_ab, max_deg_ab;
	gf2n *ap, *bp, *cp, *cpp, t;

	if (deg_a < 0) {
		c.assign(b);
		return;
	}

	if (deg_b < 0) {
		c.assign(a);
		return;
	}

	if (deg_a > deg_b) {
		max_deg_ab = deg_a;
		min_deg_ab = deg_b;
	}
	else {
		max_deg_ab = deg_b;
		min_deg_ab = deg_a;
	}

	if (max_deg_ab >= static_cast<int>(c.size)) {
		if (max_deg_ab >= static_cast<int>(gf2n_polynomial::default_size))
			c.size = max_deg_ab + 1;
		else
			c.size = gf2n_polynomial::default_size;

		cpp = new gf2n[c.size];
	}
	else cpp = c.coeff;

	c.deg = max_deg_ab;
	ap = a.coeff;
	bp = b.coeff;
	cp = cpp;

	for (i = min_deg_ab + 1; i; i--, ap++, bp++, cp++)
		add(*cp, *ap, *bp);

	if (deg_a > min_deg_ab)
		for (i = deg_a - min_deg_ab; i; i--, cp++, ap++)
			*cp = *ap;
	else
		if (deg_b > min_deg_ab)
			for (i = deg_b - min_deg_ab; i; i--, cp++, bp++)
				*cp = *bp;
		else {
			cp --;
			while ((*cp).is_zero()) {
				max_deg_ab --;
				cp --;
				if (max_deg_ab < 0)
					break;
			}
			c.deg = max_deg_ab;
		}

	if (cpp != c.coeff) {
		delete[] c.coeff;
		c.coeff = cpp;
	}
}



//******** operations: multiply **************************************
// plain version, karatzuba version (default) is in gf2n_poly_karatzuba.cc
//

void plain_mul(gf2n_polynomial & c, const gf2n_polynomial & a,
	       const gf2n_polynomial & b)
{
	register int deg_a = a.deg, deg_b = b.deg;
	register int i, j, deg_ab = deg_a + deg_b;
	gf2n *ap, *bp, t, *cp, *cpp;

	if (deg_a < 0 || deg_b < 0) {
		c.deg = -1;
		return;
	}

	if (&a == &b) {
		square(c, a);
		return;
	}

	if (deg_ab >= static_cast<int>(c.size) || c.coeff == a.coeff || c.coeff == b.coeff) {
		if (deg_ab >= static_cast<int>(gf2n_polynomial::default_size))
			c.size = deg_ab + 1;
		else
			c.size = gf2n_polynomial::default_size;

		cpp = new gf2n[c.size];
	}
	else cpp = c.coeff;

	cp = cpp;
	c.deg = deg_ab;

	for (i = deg_ab + 1; i; i--, cp++)
		(*cp).assign_zero();

	for (i = 0, ap = a.coeff; i <= deg_a; i++, ap++)
		for (j = deg_b + 1, bp = b.coeff, cp = cpp + i; j; j--, bp++, cp++) {
			multiply(t, *ap, *bp);
			add(*cp, *cp, t);
		}

	if (cpp != c.coeff) {
		delete[] c.coeff;
		c.coeff = cpp;
	}
}



//************* operations: square *********************************

void square (gf2n_polynomial & c, const gf2n_polynomial & a)
{
	register int i, dg = a.deg, dgg = (a.deg << 1);
	gf2n * cp, *cpp, *ap;

	if (dg < 0) {
		c.deg = -1;
		return;
	}

	if (dgg >= static_cast<int>(c.size) || c.coeff == a.coeff) {
		if (dgg >= static_cast<int>(gf2n_polynomial::default_size))
			c.size = dgg + 1;
		else
			c.size = gf2n_polynomial::default_size;
		cpp = new gf2n[c.size];
	}
	else cpp = c.coeff;

	c.deg = dgg;

	for (i = 0, ap = a.coeff, cp = cpp; i < dg; i++, ap++, cp++) {
		square(*cp, *ap);
		cp ++;
		(*cp).assign_zero();
	}
	square(*cp, *ap); // leading coefficient

	if (cpp != c.coeff) {
		delete[] c.coeff;
		c.coeff = cpp;
	}
}



//************** output ***********************************************

std::ostream & operator << (std::ostream & s, const gf2n_polynomial  & a)
{
	int d = a.deg;

	if (d < 0)
		s << "0";
	else
		if (d == 0)
			s << a.coeff[0];
		else
			if (d == 1) {
				if (a.coeff[1].is_one())
					s << "X";
				else
					s << a.coeff[1] << " * X";
				if (!a.coeff[0].is_zero())
					s << " + " << a.coeff[0];
			}
			else {
				if (a.coeff[d].is_one())
					s << "X^" << d;
				else
					s << a.coeff[d] << " * X^" << d;
				for (register int i = d - 1; i > 1; i--)
					if (a.coeff[i].is_one())
						s << " + X^" << i;
					else if (!a.coeff[i].is_zero())
						s << " + " << a.coeff[i] << " * X^" << i;
				if (a.coeff[1].is_one())
					s << " + X";
				else if (!a.coeff[1].is_zero())
					s << " + " << a.coeff[1] << " * X";
				if (!a.coeff[0].is_zero())
					s << " + " << a.coeff[0];
			}
	return s;
}



//*******************************************************************

void plain_div_rem(gf2n_polynomial & div, gf2n_polynomial & rem,
	           const gf2n_polynomial & a, const gf2n_polynomial & b)
{
	if (b.deg < 0)
		lidia_error_handler("gf2n_polynomial", "div::divisor is zero");

	if (a.deg < b.deg) {
		div.assign_zero();
		rem.assign(a);
		return;
	}

	gf2n_polynomial res (a.deg - b.deg), tt;
	gf2n *cp, *bp, t;

	tt.assign(a);
	gf2n blc = inverse(b.lead_coeff()), h;
	int i;

	do {
		multiply(t, tt.lead_coeff(), blc);
		res.coeff[tt.deg-b.deg].assign(t);

		for (i = 0, bp = b.coeff, cp = tt.coeff+tt.deg-b.deg;
		     i < b.deg; i++, bp++, cp++) {
			multiply(h, t, *bp);
			add(*cp, *cp, h);
		}

		do {
			tt.deg --;
			if (tt.deg < 0)
				break;
			cp --;
		} while ((*cp).is_zero());
	} while (tt.deg >= b.deg);

	div.assign(res);
	rem.assign(tt);
}



//*******************************************************************

void plain_div(gf2n_polynomial & div,
	       const gf2n_polynomial & a, const gf2n_polynomial & b)
{
	if (b.deg < 0)
		lidia_error_handler("gf2n_polynomial", "plain_div::divisor is zero");

	if (a.deg < b.deg) {
		div.assign_zero();
		return;
	}

	gf2n_polynomial res (a.deg - b.deg), tt;
	gf2n *cp, *bp, t;

	tt.assign(a);
	gf2n blc = inverse(b.lead_coeff()), h;
	int i;

	do {
		multiply(t, tt.lead_coeff(), blc);
		res.coeff[tt.deg-b.deg].assign(t);

		for (i = 0, bp = b.coeff, cp = tt.coeff+tt.deg-b.deg;
		     i < b.deg; i++, bp++, cp++) {
			multiply(h, t, *bp);
			add(*cp, *cp, h);
		}

		do {
			tt.deg --;
			if (tt.deg < 0)
				break;
			cp --;
		} while ((*cp).is_zero());
	} while (tt.deg >= b.deg);

	div.assign(res);
}



//******* GCDs ************************************************************

void gcd (gf2n_polynomial & r, const gf2n_polynomial & a,
	  const gf2n_polynomial & b)
{
	gf2n_polynomial u, v;

	if (b.is_zero()) {
		r.assign(a);
		r.make_monic();
		return;
	}

	if (a.is_zero()) {
		r.assign(b);
		r.make_monic();
		return;
	}

	if (a.deg > b.deg) {
		u.assign(a); v.assign(b);
	}
	else {
		u.assign(b); v.assign(a);
	}

	gf2n_polynomial tt;

	do {
		remainder(tt, u, v);
		u.assign(v);
		v.assign(tt);
	} while (!v.is_zero());

	u.make_monic();
	r.assign(u);
}



gf2n_polynomial xgcd_left (gf2n_polynomial & s, const gf2n_polynomial & a,
			   const gf2n_polynomial & b)

{
	// return gcd of a, b, gcd = s*a + % * b  
	gf2n_polynomial h, u, v, u1, u2;
	int ud;
	gf2n t;

	if (b.is_zero()) {
		t.assign(a.lead_coeff());
		multiply(u, a, inverse(t));
		s = gf2n_polynomial(0);
		s.set_coefficient(t, 0);
		return (u);
	}

	if (a.is_zero()) {
		s = gf2n_polynomial();
		v.assign(b);
		v.make_monic();
		return (v);
	}

	//***** Invariant:  u = u1 * a + ?? * b,   v = u2 * a + ?? * b *

	u.assign(a);
	u1 = gf2n_polynomial(0);
	v.assign(b);
	u2 = gf2n_polynomial();

	do {
		if (u.deg > v.deg) {
			multiply(t, lead_coeff(u), inverse(lead_coeff(v)));
			ud = u.deg - v.deg;
			shift_left(h, v, ud);
			multiply(h, t, h);
			add(u, h, u);

			shift_left(h, u2, ud);
			multiply(h, t, h);
			add(u1, h, u1);
		}
		else {
			multiply(t, lead_coeff(v), inverse(lead_coeff(u)));
			ud = v.deg - u.deg;
			shift_left(h, u, ud);
			multiply(h, t, h);
			add(v, h, v);

			shift_left(h, u1, ud);
			multiply(h, t, h);
			add(u2, h, u2);
		}
	} while ((!v.is_zero()) && (!u.is_zero()));

	if (v.is_zero())     // result in v and s
	{
		v.assign(u);
		s.assign(u1);
	}
	else s.assign(u2);

	if (!v.is_monic()) {
		multiply(s, inverse(lead_coeff(v)), s);
		v.make_monic();
	}
	return(v);
}



void xgcd(gf2n_polynomial& d, gf2n_polynomial& s, gf2n_polynomial& t,
	  const gf2n_polynomial& a, const gf2n_polynomial& b)
{
	if (b.is_zero()) {
		t.assign_zero();
		s.assign_one();
		d.assign(a);
	}
	else
		if (a.is_zero()) {
			s.assign_zero();
			t.assign_one();
			d.assign(b);
		}
		else {
			gf2n_polynomial u(a);
			gf2n_polynomial v(b);
			gf2n_polynomial temp, u0, v0, u1, v1, u2, v2, q, r;

			u1.assign_one(); v1.assign_zero();
			u2.assign_zero(); v2.assign_one();

			do {
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
	gf2n z;

	invert(z, d.lead_coeff());
	multiply_by_scalar(d, d, z);
	multiply_by_scalar(s, s, z);
	multiply_by_scalar(t, t, z);
}



void resultant (gf2n &r, const gf2n_polynomial & ff,
		const gf2n_polynomial & gg)
{
	gf2n_polynomial f(ff), g(gg);
	gf2n l, lead;
	int m, d;

	r.assign_one();
	m = g.degree();

	do {
		remainder(g, g, f);
		d = g.degree();

		if (d < m) {
			l.assign(f.lead_coeff());
			power(l, l, m-d);
			multiply(r, r, l);
		}

		if (g.degree() <= 0)
			break;

		swap(f, g);

		m = g.degree();
	} while (true);

	if (g.is_zero())
		r.assign_zero();
	else
		multiply(r, r, g.const_term());
}



//*********************************************************************
//
//                          inner_product
//
//*********************************************************************

void inner_product(gf2n & x, const base_vector< gf2n > & a,
		   const gf2n_polynomial &b, lidia_size_t offset)
{
	debug_handler("gf2n_polynomial", "inner_product(gf2n &, "
		      "base_vector< gf2n > &, gf2n_polynomial &, lidia_size_t)");

	lidia_size_t n = comparator< lidia_size_t >::min(a.size(),
							  b.degree() + 1 + offset);
	lidia_size_t i;
	gf2n t;

	x.assign_zero();

	for (i = offset; i < n; i++) {
		multiply(t, a[i], b[i-offset]);
		add(x, x, t);
	}
}



//------------------------------------------------------------------
//************ modular operations ********************************

void multiply_mod (gf2n_polynomial & c, const gf2n_polynomial & a,
                   const gf2n_polynomial & b, const gf2n_polynomial & f)
{
	if (&c != &f) {
		multiply(c, a, b);
		remainder(c, c, f);
	}
	else {
		gf2n_polynomial d;
		multiply(d, a, b);
		remainder(c, d, f);
	}
}



void square_mod (gf2n_polynomial & c, const gf2n_polynomial & a,
                 const gf2n_polynomial & f)
{
	if (&c != &f) {
		square(c, a);
		remainder(c, c, f);
	}
	else {
		gf2n_polynomial d;
		square(d, a);
		remainder(c, d, f);
	}
}



void invert_mod(gf2n_polynomial& x, const gf2n_polynomial& a,
		const gf2n_polynomial& f)
{
	gf2n_polynomial d;
	d = xgcd_left(x, a, f);
	if (!d.is_one())
		lidia_error_handler("gf2n_polynomial", "invert_mod(gf2n_polynomial&, "
				    "gf2n_polynomial&, gf2n_polynomial&)::can't compute"
				    "multiplicative inverse");
}



bool invert_mod_status(gf2n_polynomial& x, const gf2n_polynomial& a,
		       const gf2n_polynomial& f)
{
	gf2n_polynomial d;
	d = xgcd_left(x, a, f);
	if (!d.is_one())
		x.assign(d);
	return (d.is_one());
}



void derivative(gf2n_polynomial &x, const gf2n_polynomial &a)
{
	lidia_size_t d = a.degree();
	lidia_size_t i;

	if (d <= 0) {
		x.assign_zero();
		return;
	}

	x.assign(gf2n_polynomial(d));

	for (i = 0; i < d; i += 2) {
		x.coeff[i].assign(a.coeff[i+1]);
		x.coeff[i+1].assign_zero();
	}
	x.deg = d-1;
}



//*********************************************************************
//  uses squaring trick to compute x^i for i = 0, ..., d -1
//  no deallocation of memory !!

gf2n_polynomial * power_table(const gf2n_polynomial & x, int d)
{
	gf2n_polynomial * table;
	int i, j;

	if (d < 2)
		lidia_error_handler("gf2n_polynomial", "power_table::index d < 2");

	table = new gf2n_polynomial [d];
	table[0].assign_one();
	table[1].assign(x);

	if (d == 3) {
		square(table[2], table[1]);
		return table;
	}

	i = 2; j = 3;

	while (j < d) {
		while (i < d) {
			square(table[i], table[i/2]);
			i <<= 1;
		}
		i = j; j += 2;
		multiply(table[i], table[i-1], x);
		i <<= 1;
	}
	return table;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
