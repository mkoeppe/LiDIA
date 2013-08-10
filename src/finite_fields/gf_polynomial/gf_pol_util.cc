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
//	Author	: Thomas Pfahler (TPf), Harald Baier (HB)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/gf_polynomial.h"
#include	"LiDIA/finite_fields/Fp_polynomial_fft.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif

//
// HB: computes a root of f over GF(p^m)
//
gf_element
find_root( const gf_polynomial & f )
{
// f is monic, and has deg(f) distinct roots.
// returns  roots

    if (f.degree() == 0)
		lidia_error_handler( "gf_pol_util::find_root(const gf_polynomial)", "polynomial is constant." );
	
    if (f.degree() == 1)
		return( - f[ 0 ] );

    gf_polynomial h, local_f( f );

    const galois_field & K = f.get_field();
    const bigint &p = K.characteristic();

    if( p != 2 )
    {
		gf_element one(K), r(K);
		one.assign_one();
		const bigint & q = K.number_of_elements();
		bigint q1( q );
		q1.divide_by_2();
		
		while( local_f.degree() > 1 )
		{
			do
			{
				r.randomize();
				power_x_plus_a_mod(h, r, q1, local_f);
				subtract(h, h, one);
				gcd(h, h, local_f);
			}
			while( h.degree() <= 0 || h.degree() == local_f.degree() );

			if( h.degree() <= local_f.degree() / 2 )
				local_f.assign( h );
			else
				divide( local_f, local_f, h );
		}
    }
    else
		lidia_error_handler( "gf_pol_util::find_root(const gf_polynomial)", "case p = 2 not implemented." );

	return( - local_f[ 0 ] );
}


void
find_roots(base_vector< gf_element > &x, const gf_polynomial & f)
{
// f is monic, and has deg(f) distinct roots.
// returns the list of roots
	debug_handler("factoring.c", "find_roots(base_vector< gf_element > &, gf_polynomial&)");

	if (f.degree() == 0)
		return;

	if (f.degree() == 1) {
		lidia_size_t k = x.size();

		if (x.capacity() < k + 1)
			x.set_capacity(k + 1);
		x.set_size(k + 1);

		negate(x[k], f.const_term());
		return;
	}

	gf_polynomial h;

	const galois_field &K = f.get_field();
	const bigint &p = K.characteristic();

	if (p != 2) {
		gf_element one(K), r(K);
		one.assign_one();
		const bigint &q = K.number_of_elements();
		bigint q1(q);
		q1.divide_by_2();

		{
			do {
				r.randomize();
				power_x_plus_a_mod(h, r, q1, f);
				subtract(h, h, one);
				gcd(h, h, f);
			} while (h.degree() <= 0 || h.degree() == f.degree());
		}
	}
	else {
		// p == 2
		do {
			lidia_size_t i, k = K.degree(), n = f.degree();
			gf_polynomial g;
			g = randomize(K, n-1);
			h = g;
			for (i = 1; i < k; i++) {
				square(g, g, f);
				add(h, h, g);
			}
			gcd(h, h, f);
		} while (h.degree() <= 0 || h.degree() == f.degree());
	}

	find_roots(x, h);
	divide(h, f, h);
	find_roots(x, h);
}



base_vector< gf_element > find_roots(const gf_polynomial &f)
{
	base_vector< gf_element > v;
	find_roots(v, f);
	lidia_size_t i;
	gf_polynomial h, g;
	h.assign_x(f.get_field());
	g.assign_one(f.get_field());
	for (i = 0; i < v.size(); i++) {
		h[0] = -v[i];
		multiply(g, g, h);
	}
	if (g != f)
		std::cout << "Mist." << v << std::endl << f << std::endl << g << std::endl;
	return v;
}



void compose(gf_polynomial& x, const gf_polynomial& g,
	     const gf_polynomial& h, const gf_poly_modulus& F)
// x = g(h) mod f
{
	lidia_size_t m = square_root(g.degree() + 1);
	if (m == 0) {
		x.assign_zero(h.get_field());
		return;
	}

	gf_poly_argument A;
	A.build(h, F, m);
	A.compose(x, g, F);
}



void compose2(gf_polynomial& x1, gf_polynomial& x2,
	      const gf_polynomial& g1, const gf_polynomial& g2,
	      const gf_polynomial& h, const gf_poly_modulus& F)
	// xi = gi(h) mod f (i=1,2)
	// ALIAS RESTRICTION:  xi may not alias gj, for i != j
{
	lidia_size_t m = square_root(g1.degree() + g2.degree() + 2);
	if (m == 0) {
		x1.assign_zero(h.get_field());
		x2.assign_zero(h.get_field());
		return;
	}

	gf_poly_argument A;
	A.build(h, F, m);
	A.compose(x1, g1, F);
	A.compose(x2, g2, F);
}



void compose3(gf_polynomial& x1, gf_polynomial& x2, gf_polynomial& x3,
	      const gf_polynomial& g1, const gf_polynomial& g2, const gf_polynomial& g3,
	      const gf_polynomial& h, const gf_poly_modulus& F)
	// xi = gi(h) mod f (i=1..3)
	// ALIAS RESTRICTION:  xi may not alias gj, for i != j
{
	lidia_size_t m = square_root(g1.degree() + g2.degree() + g3.degree() + 3);
	if (m == 0) {
		x1.assign_zero(h.get_field());
		x2.assign_zero(h.get_field());
		x3.assign_zero(h.get_field());
		return;
	}

	gf_poly_argument A;
	A.build(h, F, m);
	A.compose(x1, g1, F);
	A.compose(x2, g2, F);
	A.compose(x3, g2, F);
}



void
trace_map(gf_polynomial & w, const gf_polynomial & a, lidia_size_t d,
	  const gf_poly_modulus & F, const gf_polynomial & b)
	// w = a+a^q+...+^{q^{d-1}} mod f;
	// it is assumed that d >= 0, and b = X^q mod f, q a power of p
{
	debug_handler("gf_polynomial", "trace_map(gf_polynomial&, gf_polynomial&, lidia_size_t, gf_poly_modulus&, gf_polynomial&)");

	w.ffield = (gf_polynomial::common_field(a.ffield, b.ffield));
	gf_polynomial::build_frame(w.ffield);


	gf_polynomial z(b), y(a), t;
	w.assign_zero();

	while (d) {
		if (d == 1) {
			if (w.is_zero())
				w.assign(y);
			else {
				compose(w, w, z, F);
				add(w, w, y);
			}
		}
		else {
			if ((d & 1) == 0) {
				compose2(z, t, z, y, z, F);
				add(y, t, y);
			}
			else {
				if (w.is_zero()) {
					w.assign(y);
					compose2(z, t, z, y, z, F);
					add(y, t, y);
				}
				else {
					compose3(z, t, w, z, y, w, z, F);
					add(w, w, y);
					add(y, t, y);
				}
			}
		}
		d = d >> 1;
	}

	gf_polynomial::delete_frame();
}



#if 0

void
power_compose(gf_polynomial & y, const gf_polynomial & h, lidia_size_t q, const gf_polynomial & f)
	// w = X^{q^d} mod f;
	// it is assumed that d >= 0, and b = X^q mod f, q a power of p
{
	debug_handler("gf_polynomial", "power_compose(gf_polynomial&, gf_polynomial&, lidia_size_t, gf_polynomial&)");

	y.ffield = (gf_polynomial::common_field(h.ffield, f.ffield));
	gf_polynomial::build_frame(y.ffield);

	gf_polynomial z(h);
	lidia_size_t sw;

	y.assign_x();

	while (q) {
		sw = 0;

		if (q > 1)
			sw = 2;
		if (q & 1) {
			if (y.is_x())
				y.assign(z);
			else
				sw = sw | 1;
		}

		switch (sw) {
		case 0:
			break;

		case 1:
			compose(y, y, z, F);
			break;

		case 2:
			compose(z, z, z, F);
			break;

		case 3:
			compose2(y, z, y, z, z, F);
			break;
		}

		q = q >> 1;
	}

	gf_polynomial::delete_frame();
}

#endif



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
