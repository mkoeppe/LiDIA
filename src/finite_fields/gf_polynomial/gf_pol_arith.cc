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
#include	"LiDIA/gf_polynomial.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//******************************************************************
//		classical arithmetic
//******************************************************************

void plain_power(polynomial< gf_element > & c,
		 const polynomial< gf_element > & a, const bigint & b)
{
	//"classical" version
	c.ffield = a.ffield;
	polynomial< gf_element >::build_frame(a.ffield);
	power(c.pol, a.pol, b);
	c.delete_frame();
}



void plain_gcd(polynomial< gf_element > &d,
	       const polynomial< gf_element > &aa, const polynomial< gf_element > &bb)
{
	//"classical" version
	d.ffield = gf_polynomial::common_field(aa.ffield, bb.ffield);
	polynomial< gf_element >::build_frame(d.ffield);
	gcd(d.pol, aa.pol, bb.pol);
	polynomial< gf_element >::delete_frame();
	if (d.pol.degree() == 0)			// if gcd = 1, pol[0] is not
		d.pol[0].assign_one(d.ffield); // initialized
}



void plain_multiply(gf_polynomial &c, const gf_polynomial &a,
		    const gf_polynomial &b)
{
	//"classical" version
	c.ffield = gf_polynomial::common_field(a.ffield, b.ffield);
	polynomial< gf_element >::build_frame(c.ffield);
	multiply(c.pol, a.pol, b.pol);
	polynomial< gf_element >::delete_frame();
}



void plain_square(polynomial< gf_element > & c,
		  const polynomial< gf_element > & a)
{
	//"classical" version
	c.ffield = a.ffield;
	polynomial< gf_element >::build_frame(a.ffield);
	multiply(c.pol, a.pol, a.pol);
	polynomial< gf_element >::delete_frame();
}



void plain_div_rem(gf_polynomial &q, gf_polynomial &r,
		   const gf_polynomial &a, const gf_polynomial &b)
	//r = a - q*b, "classical" version
{
	q.ffield = gf_polynomial::common_field(a.ffield, b.ffield);
	r.ffield = q.ffield;
	polynomial< gf_element >::build_frame(q.ffield);
	div_rem(q.pol, r.pol, a.pol, b.pol);
	polynomial< gf_element >::delete_frame();
}



void plain_divide(gf_polynomial &q,
		  const gf_polynomial &a, const gf_polynomial &b)
	//r = a - q*b, "classical" version
{
	gf_polynomial r;
	q.ffield = gf_polynomial::common_field(a.ffield, b.ffield);
	r.ffield = q.ffield;
	polynomial< gf_element >::build_frame(q.ffield);
	div_rem(q.pol, r.pol, a.pol, b.pol);
	polynomial< gf_element >::delete_frame();
}



void plain_remainder(gf_polynomial &r,
		     const gf_polynomial &a, const gf_polynomial &b)
	//r = a - q*b, "classical" version
{
	gf_polynomial q;
	q.ffield = gf_polynomial::common_field(a.ffield, b.ffield);
	r.ffield = q.ffield;
	polynomial< gf_element >::build_frame(q.ffield);
	div_rem(q.pol, r.pol, (base_polynomial< gf_element > )a.pol,
		(base_polynomial< gf_element > )b.pol);
	polynomial< gf_element >::delete_frame();
}



//******************************************************************
//  the following 2 functions are the key to our fast multiplication
//  and division routines :
//******************************************************************

void to_Kronecker(Fp_polynomial &g, const gf_polynomial &f,
		  lidia_size_t lo, lidia_size_t hi)
	//computes Kronecker substitution:
	//y -> x^(2n-1)  for  g in GF(p^n)[y], f in Fp[x]
	//only g[lo..hi] is converted
{
	if (lo < 0)
		lidia_error_handler("gf_polynomial", "to_Kronecker(...)::bad indices");

	hi = comparator< lidia_size_t >::min(hi, f.degree());

	lidia_size_t m = comparator< lidia_size_t >::max(hi - lo + 1, 0);
	lidia_size_t n = f.get_field().degree();
	g.set_modulus(f.get_field().characteristic());
	g.set_max_degree((2*n-1)*m);
	//Kronecker substitution:
	//y -> x^(2n-1)
	lidia_size_t i, j;
	//###	std::cout<<"to, n = "<<n<<"  m = "<<m<<std::endl;
	for (i = lo; i <= hi; i++) {
		for (j = 0; j < n; j++) {
	//###	std::cout<<"g["<<(i-lo)*(2*n-1) + j<<"] = f["<<i<<"]["<<j<<"]"<<std::endl;
			g[(i-lo)*(2*n-1) + j].assign((f[i].polynomial_rep2())[j]);
		}
	}
}



void from_Kronecker(gf_polynomial &f, const Fp_polynomial &g,
		    lidia_size_t lo, lidia_size_t hi)
	//computes "backwards" Kronecker substitution:
	//X^(2n-1) -> Y  for  g in GF(p^n)[Y], f in Fp[X]
	//only coefficients [lo..hi] (in Y) are converted
{
	if (lo < 0)
		lidia_error_handler("gf_polynomial", "from_Kronecker(...)::bad indices");

	lidia_size_t n = f.get_field().degree();
	lidia_size_t m = comparator< lidia_size_t >::max(hi - lo + 1, 0);
	if (m == 0) {
		f.assign_zero();
		return;
	}
	const bigint &p = f.get_field().characteristic();
	Fp_polynomial tmp;
	tmp.set_modulus(p);
	tmp.set_max_degree(n-1);
	f.set_degree(m-1);
	lidia_size_t i, j;
	//### std::cout<<"from, n = "<<n<<"  m = "<<m<<std::endl;

	//Kronecker substitution:
	//x^(2n-1) -> y
	for (i = lo; i <= hi; i++) {
		for (j = 0; j < (2*n-1); j++) {
			tmp[j] = g[i*(2*n-1) + j];
		}
	//### std::cout<<"i-lo = "<<i-lo<<"   tmp = "<<tmp<<std::endl;
		f[i-lo].set_polynomial_rep(tmp);
	}
	f.remove_leading_zeros();
}



//******************************************************************
//		"Kronecker" based arithmetic
//******************************************************************


void copy_reverse(gf_polynomial &x, const gf_polynomial &a,
		  lidia_size_t lo, lidia_size_t hi)
	// x[0..hi-lo+1] = reverse(a[lo..hi]), with zero fill
	// input may not alias output
{
	debug_handler("gf_polynomial", "copy_reverse(gf_polynomial&, gf_polynomial&, lidia_size_t, lidia_size_t) ");
	lidia_size_t i, j, n, m;

	n = hi-lo+1;
	m = a.degree()+1;

	x.ffield = a.ffield;
	polynomial< gf_element >::build_frame(a.ffield);
	x.set_degree(n-1);

	for (i = 0; i < n; i++) {
		j = hi-i;
		if (j< 0 || j >= m)  x[i].assign_zero();
		else                  x[i].assign(a[j]);
	}

	x.remove_leading_zeros();
	polynomial< gf_element >::delete_frame();
}



// x = (1/a) % X^m, input not output, constant term a is nonzero
void invert(gf_polynomial& x, const gf_polynomial& a, lidia_size_t m)
{
	debug_handler("gf_polynomial", "invert(gf_polynomial&, gf_polynomial&, lidia_size_t)");

	const galois_field K = a.ffield;

	lidia_size_t i, k, n, lb;
	gf_element s(K);

	n = a.degree();
	if (n < 0) lidia_error_handler("gf_polynomial", "invert(gf_polynomial&, gf_polynomial&, lidia_size_t)::division by zero");

	invert(s, a.const_term());

	if (n == 0) {
		x.assign(s);
		return;
	}

	x.ffield = K;
	polynomial< gf_element >::build_frame(K);
	x.set_degree(m-1);
	x[0].assign(s);
	bool ct_is_one = s.is_one();

	if (K.characteristic() == 2) {
		gf_element v(K), t(K);
		for (k = 1; k < m; k++) {
			v.assign_zero();
			lb = comparator< lidia_size_t >::max(k-n, 0);
			for (i = lb; i <= k-1; i++) {
				multiply(t, x[i], a[k-i]);
				add(v, v, t);
			}
			v.negate();
			if (!ct_is_one)  multiply(v, v, s);
			x[k].assign(v);
		}
	}
	else {
		Fp_polynomial v, t;
		v.set_modulus(K.characteristic());
		for (k = 1; k < m; k++) {
			v.assign_zero();
			lb = comparator< lidia_size_t >::max(k-n, 0);
			for (i = lb; i <= k-1; i++) {
				multiply(t, x[i].polynomial_rep(), a[k-i].polynomial_rep());
				add(v, v, t);
			}
			v.negate();
			if (!ct_is_one)  multiply(v, v, s.polynomial_rep());
			x[k].set_polynomial_rep(v);
		}
	}
	x.remove_leading_zeros();
	polynomial< gf_element >::delete_frame();

#if 0
	// KONTROLLE
	gf_polynomial w1, w2, w3, w4;
	w4.build_frame(a.ffield);
	w4.set_degree(m);
	w4[m].assign_one();
	xgcd(w1, w2, w3, w4, a);
	if (w3 != x)  std::cout << "invert: error" << std::endl;
	else          std::cout << "invert: ok" << std::endl;
#endif
}



void div_rem(gf_polynomial &q, gf_polynomial &r,
	     const gf_polynomial &a, const gf_polynomial &b)
	//r = a - q*b
{
	lidia_size_t deg_b = b.degree(), deg_a = a.degree();
	if (deg_b < 0)
		lidia_error_handler("polynomial< gf_element >", "div_rem::division by zero");

	const galois_field &K = gf_polynomial::common_field(a.ffield, b.ffield);
	if (deg_a < deg_b) {
		q.assign_zero(K);
		r.assign(a);
		return;
	}

	if (K.characteristic() == 2
	    && (K.degree() > 1 || K.degree() == 1 && deg_a < 2*deg_b))
	{
		//"classical" version
		q.ffield = K;
		r.ffield = K;
		gf_polynomial::build_frame(K);
		div_rem(q.pol, r.pol, a.pol, b.pol);
		gf_polynomial::delete_frame();
		return;
	}

	gf_polynomial p1, p2, p3;
	Fp_polynomial r1, r2;
	copy_reverse(p3, b, 0, deg_b);
	invert(p2, p3, deg_a-deg_b+1);
	copy_reverse(p1, p2, 0, deg_a-deg_b);

	to_Kronecker(r1, p1, 0, p1.degree());
	to_Kronecker(r2, a, deg_b, deg_a);
//###	std::cout<<"p1:\n"<<p1<<std::endl<<r1<<std::endl;
//###	std::cout<<"a:  db="<<deg_b<<"  da="<<deg_a<<std::endl<<a<<std::endl<<r2<<std::endl;
	multiply(r1, r1, r2);
	from_Kronecker(p1, r1, deg_a-deg_b, 2*(deg_a-deg_b));

#if 0
	// KONTROLLE
	gf_polynomial tmp, tmp2;
	tmp.build_frame(a.ffield);
	tmp.set_degree(deg_a-deg_b);
	for (lidia_size_t i = deg_b; i <= deg_a; i++)
		tmp[i-deg_b] = a[i];
	multiply(tmp2, tmp, p1);
	std::cout << "tmp:\n" << tmp << std::endl << "tmp2:\n" << tmp2 << std::endl;
#endif

	multiply(p2, b, p1);
	subtract(r, a, p2);
	q.assign(p1);
}



void divide(gf_polynomial &q, const gf_polynomial &a, const gf_polynomial &b)
	//q = a / b
{
	lidia_size_t deg_b = b.degree(), deg_a = a.degree();
	if (deg_b < 0)
		lidia_error_handler("polynomial< gf_element >", "divide::division by zero");
	q.ffield = gf_polynomial::common_field(a.ffield, b.ffield);
	gf_polynomial::build_frame(q.ffield);
	if (deg_a < deg_b) {
		q.assign_zero();
		gf_polynomial::delete_frame();
		return;
	}

	if (q.ffield.characteristic() == 2 && q.ffield.degree() > 1) {
		//"classical" version
		divide(q.pol, a.pol, b.pol);
		gf_polynomial::delete_frame();
		return;
	}

	gf_polynomial r;

	gf_polynomial p1, p2, p3;
	Fp_polynomial r1, r2;
	copy_reverse(p3, b, 0, deg_b);
	invert(p2, p3, deg_a-deg_b+1);
	copy_reverse(p1, p2, 0, deg_a-deg_b);

	to_Kronecker(r1, p1, 0, p1.degree());
	to_Kronecker(r2, a, deg_b, deg_a);
	multiply(r1, r1, r2);
	from_Kronecker(q, r1, deg_a-deg_b, 2*(deg_a-deg_b));
	polynomial< gf_element >::delete_frame();
}



void remainder(gf_polynomial &r, const gf_polynomial &a, const gf_polynomial &b)
	//r = a - (a/b)*b
{
	lidia_size_t deg_b = b.degree(), deg_a = a.degree();
	if (deg_b < 0)
		lidia_error_handler("polynomial< gf_element >", "remainder::division by zero");
	gf_polynomial q;
	q.ffield = gf_polynomial::common_field(a.ffield, b.ffield);
	gf_polynomial::build_frame(q.ffield);
	if (deg_a < deg_b) {
		q.assign_zero();
		r.assign(a);
		gf_polynomial::delete_frame();
		return;
	}

	if (q.ffield.characteristic() == 2 && q.ffield.degree() > 1) {
		//"classical" version
		r.ffield = q.ffield;
		remainder(r.pol, a.pol, b.pol);
		gf_polynomial::delete_frame();
		return;
	}

	gf_polynomial p1, p2, p3;
	Fp_polynomial r1, r2;
	copy_reverse(p3, b, 0, deg_b);
	invert(p2, p3, deg_a-deg_b+1);
	copy_reverse(p1, p2, 0, deg_a-deg_b);

	to_Kronecker(r1, p1, 0, p1.degree());
	to_Kronecker(r2, a, deg_b, deg_a);
	multiply(r1, r1, r2);
	from_Kronecker(q, r1, deg_a-deg_b, 2*(deg_a-deg_b));
	polynomial< gf_element >::delete_frame();

	multiply(p1, b, q);
	subtract(r, a, p1);
}



void multiply(gf_polynomial &c, const gf_polynomial &a,
	      const gf_polynomial &b)
{
	c.ffield = gf_polynomial::common_field(a.ffield, b.ffield);
	polynomial< gf_element >::build_frame(c.ffield);

	if (c.ffield.characteristic() == 2 && c.ffield.degree() > 2) {
		//"classical" version
		multiply(c.pol, a.pol, b.pol);
	}
	else {
		Fp_polynomial r1, r2, r3;

		to_Kronecker(r1, a, 0, a.degree());
		to_Kronecker(r2, b, 0, b.degree());
		multiply(r3, r1, r2);
		from_Kronecker(c, r3, 0, a.degree()+b.degree());
	}

	polynomial< gf_element >::delete_frame();
}



void square(polynomial< gf_element > & c,
	    const polynomial< gf_element > & a)
{
	c.ffield = a.ffield;
	polynomial< gf_element >::build_frame(a.ffield);

	if (c.ffield.characteristic() == 2 && c.ffield.degree() > 2) {
		//"classical" version
		multiply(c.pol, a.pol, a.pol);
	}
	else {
		Fp_polynomial r1, r2;

		to_Kronecker(r1, a, 0, a.degree());
		square(r2, r1);
		from_Kronecker(c, r2, 0, 2*a.degree());
	}

	polynomial< gf_element >::delete_frame();
}



void gcd(polynomial< gf_element > &d,
	 const polynomial< gf_element > &aa, const polynomial< gf_element > &bb)
{
	debug_handler("gf_polynomial", "gcd(...)");

	if (bb.is_zero())
		d.assign(aa);
	else if (aa.is_zero())
		d.assign(bb);
	else {
		gf_polynomial r, a(aa), b(bb);
		do
		{
			remainder(a, a, b);
			swap(a, b);
		} while (!b.is_zero());
		d.assign(a);
	}

	if (d.is_zero())
		return;

	if (d.degree() == 0) {
		d.assign_one();
		return;
	}

	if (!d.lead_coeff().is_one()) {
		gf_element lc_inv;
		invert(lc_inv, d.lead_coeff());
		multiply(d, d, lc_inv);
	}
}



void power(polynomial< gf_element > & c,
	   const polynomial< gf_element > & a, const bigint & b)
{
	bigint exponent;
	gf_polynomial multiplier;
	if (b.is_negative())
		lidia_error_handler("gf_polynomial", "power(...)::negative exponent");

	c.assign_one(a.get_field());
	if (b.is_zero() || a.is_one())
		return;

	exponent.assign(b);
	multiplier.assign(a);
	while (exponent.is_gt_zero()) {
		if (!exponent.is_even())
			multiply(c, c, multiplier);
		square(multiplier, multiplier);
		exponent.divide_by_2();
	}
}



//***************************************************************
//
//		Modular Arithmetic without pre-conditioning
//
//***************************************************************

// arithmetic mod f.
// all inputs and outputs are polynomials of degree less than deg(f).
// ASSUMPTION: f is assumed monic, and deg(f) > 0.
// ALIAS RESTRICTIONS: f (and exponent e) can not alias an output.


void
multiply_mod(gf_polynomial & x, const gf_polynomial & a, const gf_polynomial & b, const gf_polynomial & f)
{
	debug_handler_l("gf_polynomial", "mul_mod(gf_polynomial&, gf_polynomial&, gf_polynomial&, gf_polynomial&)", 6);

	gf_polynomial t;
	multiply(t, a, b);
	remainder(x, t, f);
}



void
square_mod(gf_polynomial & x, const gf_polynomial & a, const gf_polynomial & f)
{
	debug_handler_l("gf_polynomial", "sqr_mod(gf_polynomial&, gf_polynomial&, gf_polynomial&)", 6);

	gf_polynomial t;
	square(t, a);
	remainder(x, t, f);
}



void
multiply_by_x_mod(gf_polynomial & h, const gf_polynomial & a, const gf_polynomial & f)
{
	debug_handler_l("gf_polynomial", "mul_by_x_mod(gf_polynomial&, gf_polynomial&, gf_polynomial&)", 6);

	h.ffield = gf_polynomial::common_field(a.ffield, f.ffield);
	polynomial< gf_element >::build_frame(h.ffield);

	lidia_size_t i, n, m;
	gf_element t, z;

	n = f.degree();
	m = a.degree();

	if (m < n - 1) {
		h.set_degree(m + 1);
		for (i = m + 1; i >= 1; i--)
			h[i].assign(a[i - 1]);
		h[0].assign_zero(f.get_field());
	}
	else {
		if (&f == &h) {
			//allows f to alias output

			h.assign_zero();
			h.delete_frame();
			return;
		}
		if (m >= n || !f.is_monic()) {
			gf_polynomial tmp;
			tmp.assign_x(f.get_field());
			multiply(h, a, tmp);
			remainder(h, h, f);
			h.delete_frame();
			return;
		}
		//now, m = n-1

		h.set_degree(n - 1);
		negate(z, a[n - 1]);
		for (i = n - 1; i >= 1; i--) {
			multiply(t, z, f[i]);
			add(h[i], a[i - 1], t);
		}
		multiply(h[0], z, f[0]);
		h.remove_leading_zeros();
	}
	polynomial< gf_element >::delete_frame();
}



void
invert_mod(gf_polynomial & x, const gf_polynomial & a, const gf_polynomial & f)
{
	debug_handler_l("gf_polynomial", "inv_mod(gf_polynomial&, gf_polynomial&, gf_polynomial&)", 6);

	gf_polynomial d, t;
	xgcd(d, x, t, a, f);
	if (!d.is_one())
		lidia_error_handler_c("gf_polynomial", "inv_mod(...)::"
				      "can't compute multiplicative inverse",
				      std::cout << "d = " << d << "a = " << a << "f = " << f << std::endl;
			);
}



bool
invert_mod_status(gf_polynomial & x, const gf_polynomial & a, const gf_polynomial & f)
{
	debug_handler_l("gf_polynomial", "inv_mod_status(gf_polynomial&, gf_polynomial&, gf_polynomial&)", 6);

	gf_polynomial d, t;
	xgcd(d, x, t, a, f);
	if (!d.is_one())
		x.assign(d);
	return (d.is_one());
}



void
power_mod(gf_polynomial & h, const gf_polynomial & g, const bigint & e, const gf_polynomial & f)
{
	debug_handler_l("gf_polynomial", "power_mod(gf_polynomial&, gf_polynomial&, bigint&, gf_polynomial&)", 6);

	gf_polynomial lg(g);
	lidia_size_t n = e.bit_length();

	h.ffield = gf_polynomial::common_field(g.ffield, f.ffield);
	polynomial< gf_element >::build_frame(h.ffield);
	h.assign_one();

	for (lidia_size_t i = n - 1; i >= 0; i--) {
		square_mod(h, h, f);
		if (e.bit(i))
			multiply_mod(h, h, lg, f);
	}
	polynomial< gf_element >::delete_frame();
}



void
power_x_mod(gf_polynomial & h, const bigint & e, const gf_polynomial & f)
{
	debug_handler_l("gf_polynomial", "power_x_mod(gf_polynomial&, bigint&, gf_polynomial&)", 6);

	if (&h == &f)
		lidia_error_handler("gf_polynomial", "power_x_mod(gf_polynomial&, bigint&, gf_polynomial&)::no alias allowed");

	h.ffield = f.ffield;
	polynomial< gf_element >::build_frame(f.ffield);
	h.set_degree(-1);

	lidia_size_t n;
	if (e < f.degree()) {
		e.sizetify(n);
		h.set_degree(n);
		h[n] = 1;
		//###	std::cout<<"1)  X^"<<e<<" mod "<<f<<" = "<<h<<std::endl;
	}
	else {
		n = e.bit_length();
		h.assign_one();

		for (lidia_size_t i = n - 1; i >= 0; i--) {
			square_mod(h, h, f);
			if (e.bit(i))
				multiply_by_x_mod(h, h, f);
		}
		//###	std::cout<<"2)  X^"<<e<<" mod "<<f<<" = "<<h<<std::endl;
	}

	polynomial< gf_element >::delete_frame();
}



// WARNING: obsolete.  Use power_mod with Fp_poly_modulus (see below).
void
power_x_plus_a_mod(gf_polynomial & h, const gf_element & a, const bigint & e, const gf_polynomial & f)
{
	debug_handler_l("gf_polynomial", "power_x_plus_a_mod(gf_polynomial&, bigint&, bigint&, gf_polynomial&)", 6);

	h.ffield = f.ffield;
	polynomial< gf_element >::build_frame(f.ffield);

	lidia_size_t i;
	gf_polynomial t1, t2;
	gf_element la(a);
	lidia_size_t n = e.bit_length();

	h.assign_one();

	for (i = n - 1; i >= 0; i--) {
		square_mod(h, h, f);
		if (e.bit(i)) {
			multiply_by_x_mod(t1, h, f);
			multiply(t2, h, la);
			add(h, t1, t2);
		}
	}

	polynomial< gf_element >::delete_frame();
}



// computes x = a mod X^m-1
void
cyclic_reduce(gf_polynomial & x, const gf_polynomial & a, lidia_size_t m)
{
	debug_handler_l("gf_polynomial", "cyclic_reduce(gf_polynomial&, gf_polynomial&, lidia_size_t)", 6);


	lidia_size_t n = a.degree();
	lidia_size_t i, j;
	gf_element accum; //static

	if (n < m) {
		x.assign(a);
		return;
	}

	if (&x == &a) {
		x.set_degree(m - 1);
		return;
	}

	x.ffield = a.ffield;
	polynomial< gf_element >::build_frame(a.ffield);
	x.set_degree(m - 1);

	for (i = 0; i < m; i++) {
		accum.assign(a[i]);
		for (j = i + m; j <= n; j += m)
			add(accum, accum, a[j]);
		x[i].assign(accum);
	}

	x.remove_leading_zeros();
	polynomial< gf_element >::delete_frame();
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
