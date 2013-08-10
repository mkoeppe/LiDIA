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
#include	"LiDIA/Fp_poly_modulus.h"
#include	"LiDIA/Fp_poly_multiplier.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//**********************************************************
//
//		   Modular Composition
//	    algorithms for computing g(h) mod f
//
//**********************************************************

//
// Task:        Computes x = g(h) mod f.
//

void compose(Fp_polynomial& x, const Fp_polynomial& g,
             const Fp_polynomial& h, const Fp_poly_modulus& F)
{
	debug_handler("compose.c", "compose(Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_poly_modulus&)");

	g.comp_modulus(h, "compose");

	lidia_size_t m = square_root(g.degree() + 1);

	if (m == 0) {
		x.set_modulus(h); //assigns zero
		return;
	}

	poly_argument A;
	A.build(h, F, m);
	A.compose(x, g, F);
}



//
// Task:        Computes xi = gi(h) mod f (i=1,2).
//
// Conditions:  ALIAS RESTRICTION:  xi may not alias gj, for i != j
//

void compose2(Fp_polynomial& x1, Fp_polynomial& x2,
              const Fp_polynomial& g1, const Fp_polynomial& g2,
              const Fp_polynomial& h, const Fp_poly_modulus& F)
{
	debug_handler("compose.c", "compose2(Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_poly_modulus&)");

	g1.comp_modulus(h, "compose2");
	g2.comp_modulus(h, "compose2");

	lidia_size_t m = square_root(g1.degree() + g2.degree() + 2);

	if (m == 0) {
		x1.set_modulus(h); //assigns zero
		x2.set_modulus(h); //assigns zero
		return;
	}

	if (&x1 == &g2) {
		lidia_error_handler("compose.c",
				    "compose2 (...)::alias restriction - x1 may not alias g2");
		return;
	}

	poly_argument A;
	A.build(h, F, m);

	A.compose(x1, g1, F);
	A.compose(x2, g2, F);
}



//
// Task:        Computes xi = gi(h) mod f (i=1,2,3).
//
// Conditions:  ALIAS RESTRICTION:  xi may not alias gj, for i != j
//

void compose3(Fp_polynomial& x1, Fp_polynomial& x2, Fp_polynomial& x3,
	      const Fp_polynomial& g1, const Fp_polynomial& g2, const Fp_polynomial& g3,
	      const Fp_polynomial& h, const Fp_poly_modulus& F)
{
	debug_handler("compose.c", "compose3(Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_poly_modulus&)");

	g1.comp_modulus(h, "compose3");
	g2.comp_modulus(h, "compose3");
	g3.comp_modulus(h, "compose3");

	lidia_size_t m = square_root(g1.degree() + g2.degree() + g3.degree() + 3);

	if (m == 0) {
		x1.set_modulus(h); //assigns zero
		x2.set_modulus(h); //assigns zero
		x3.set_modulus(h); //assigns zero
		return;
	}

	if (&x1 == &g2 || &x1 == &g3 || &x2 == &g3) {
		lidia_error_handler("compose.c",
				    "compose3 (...)::alias restriction - xi may not alias gj for i < j");
		return;
	}

	poly_argument A;
	A.build(h, F, m);

	A.compose(x1, g1, F);
	A.compose(x2, g2, F);
	A.compose(x3, g3, F);
}



//*********************************************************************
//
//			project_powers
//
//*********************************************************************

//
// Task:        Computes (a, 1), (a, h), ..., (a, h^{k-1} % f).
//              This operation is really a "transposed" of the
//              modular composition operation.
//

void project_powers(base_vector< bigint > & x, const base_vector< bigint > & a,
		    lidia_size_t k, const Fp_polynomial& h, const Fp_poly_modulus& F)
{
	debug_handler("Fp_polynomial", "project_powers(base_vector< bigint > &, base_vector< bigint > &, lidia_size_t, Fp_polynomial&, Fp_poly_modulus&)");

	if (k < 0) {
		lidia_error_handler("Fp_polynomial", "project_powers("
				    "base_vector< bigint > &, base_vector< bigint > &, lidia_size_t, "
				    "Fp_polynomial&, Fp_poly_modulus&)::bad args");
		return;
	}

	if (k == 0) {
		x.set_size(0);
		return;
	}

	lidia_size_t m = square_root(k);

	poly_argument H;
	H.build(h, F, m);

	project_powers(x, a, k, H, F);
}



//
// Task:        Same as above, but uses a pre-computed poly_argument
//

void project_powers(base_vector< bigint > & x, const base_vector< bigint > & a,
		    lidia_size_t k, const poly_argument& H, const Fp_poly_modulus& F)
{
	debug_handler("Fp_polynomial", "project_powers(base_vector< bigint > &, base_vector< bigint > &, lidia_size_t, poly_argument&, Fp_poly_modulus&)");

	lidia_size_t n = F.modulus().degree();
	lidia_size_t m = H.H.size()-1;
	lidia_size_t i, j, l = (k+m-1)/m - 1;

	Fp_poly_multiplier M(H.H[m], F);

	base_vector< bigint > s(n, n);
	s.assign(a);

	if (x.capacity() < k)
		x.set_capacity(k);
	x.set_size(k);

	for (i = 0; i <= l; i++) {
		lidia_size_t m1 = comparator< lidia_size_t >::min(m, k-i*m);
		bigint* w = &x[i*m];
		for (j = 0; j < m1; j++)
			inner_product(w[j], s, H.H[j]);
		if (i < l)
			update_map(s, s, M, F);
	}
}



//*********************************************************************
//
//	    		    inner_product
//
//*********************************************************************

void inner_product(bigint& x, const base_vector< bigint > & a,
		   const Fp_polynomial &b, lidia_size_t offset)
{
	debug_handler("Fp_polynomial", "inner_product(bigint&, base_vector< bigint > &, Fp_polynomial&, lidia_size_t)");

	lidia_size_t n = comparator< lidia_size_t >::min(a.size(),
							  b.degree()+1+offset);
	lidia_size_t i;
	bigint accum, t;

	for (i = offset; i < n; i++) {
		multiply(t, a[i], b[i-offset]);
		add(accum, accum, t);
	}
	Remainder(x, accum, b.modulus());
}



//*********************************************************************
//
//			    update_map
//
//*********************************************************************


//
// Task:        Computes (a, b), (a, (b*X)%f), ..., (a, (b*X^{n-1})%f),
//              where ( , ) denotes the vector inner product.
//              This is really a "transposed" MulMod by B.

void update_map(base_vector< bigint > & x, const base_vector< bigint > & a,
		const Fp_poly_multiplier& B, const Fp_poly_modulus& F)
{
	debug_handler("Fp_polynomial", "update_map(base_vector< bigint > &, base_vector< bigint > & Fp_poly_multiplier&, Fp_poly_modulus&)");

	if (B.poly_mod_ptr != &F) {
		lidia_error_handler("Fp_polynomial",
				    "update_map(base_vector< bigint > &, base_vector< bigint > & "
				    "Fp_poly_multiplier&, Fp_poly_modulus&)::wrong moduli");
		return;
	}

	lidia_size_t n = F.n;
	lidia_size_t i;

	if (!B.use_FFT) {
		update_map(x, a, B.b, F.f);
		return;
	}

	fft_rep R1, R2;
	R1.init(F.fd);
	R2.init(F.fd);

	base_vector< bigint > V1(n, n);

	R1.set_size(F.k);
	R2.set_size(F.k);
	R1.rev_to_fft_rep(a, 0, a.size()-1, 0);
	multiply(R2, R1, F.FRep);
	R2.rev_from_fft_rep(V1, 0, n-2);

	for (i = 0; i <= n-2; i++)
		NegateMod(V1[i], V1[i], F.f.modulus());

	R2.set_size(F.l);
	R2.rev_to_fft_rep(V1, 0, n-2, n-1);
	multiply(R2, B.B1, R2);
	multiply(R1, B.B2, R1);

	R2.add_expand(R1);
	R2.rev_from_fft_rep(x, 0, n-1);
}



//
// Task: Same as above, but uses only classical arithmetic
//

void update_map(base_vector< bigint > & xx, const base_vector< bigint > & a,
		const Fp_polynomial& b, const Fp_polynomial& f)
{
	debug_handler("Fp_polynomial", "update_map(base_vector< bigint > &, base_vector< bigint > &, Fp_polynomial&, Fp_polynomial&)");

	b.comp_modulus(f, "update_map");

	lidia_size_t n = f.degree();
	lidia_size_t i, m;

	if (b.is_zero()) {
		if (xx.capacity() < n)
			xx.set_capacity(n);
		xx.set_size(n);
		for (i = 0; i < n; i++)
			xx[i].assign_zero();
		return;
	}

	m = n-1 - b.degree();

	base_vector< bigint > x(n, n);

	for (i = 0; i <= m; i++)
		inner_product(x[i], a, b, i);

	if (b.degree() != 0) {
		Fp_polynomial c;
		c.set_max_degree(n-1);
		shift_left(c, b, m);

		for (i = m+1; i < n; i++) {
			multiply_by_x_mod(c, c, f);
			inner_product(x[i], a, c);
		}
	}

	xx = x;
}



//*********************************************************************
//
//		prob_min_poly, min_poly, irred_poly
//
//*********************************************************************

static
void min_poly_work(Fp_polynomial &h, const base_vector< bigint > &R,
		   lidia_size_t m, const Fp_polynomial &g, const Fp_poly_modulus &F)
{
	debug_handler("compose.c", "min_poly_work(Fp_polynomial&, base_vector< bigint > &, lidia_size_t, Fp_polynomial&, Fp_poly_modulus&)");

	base_vector< bigint > x;
	project_powers(x, R, 2*m, g, F);

	Fp_polynomial a, b;
	a.set_modulus(g);
	lidia_size_t i;
	for (i = 2*m-1; i >= 0; i--)
		a.set_coefficient(x[i], 2*m-i-1);

	b.set_modulus(g);
	b.set_coefficient(2*m);

	const bigint &p = g.modulus();
	poly_matrix M(p);
	M.half_gcd(b, a, m+1);

	// make monic
	bigint t;
	InvMod(t, (M(1, 1)).lead_coeff(), p);
	multiply_by_scalar(h, M(1, 1), t);
}



//
// Literature:  Shoup, J. Symbolic Comp. 17:371-391, 1994
//              Shoup, J. Symbolic Comp. 20:363-397, 1995
//
// Task:        Computes the monic minimal polynomial of (g mod f)
//              (f = F.modulus()).
//              The algorithm is probabilistic, always returns a divisor of
//              the minimal polynomial, and returns a proper divisor with
//              probability at most m/p. (p = f.modulus())
//
// Conditions:  m = a bound on the degree of the minimal polynomial
//

void prob_min_poly(Fp_polynomial& h, const Fp_polynomial& g, lidia_size_t m,
		   const Fp_poly_modulus& F)
{
	debug_handler("Fp_polynomial", "prob_min_poly(Fp_polynomial&, Fp_polynomial&, lidia_size_t, Fp_poly_modulus&)");

	const Fp_polynomial &f = F.modulus();
	g.comp_modulus(f, "prob_min_poly");

	lidia_size_t n = f.degree();
	const bigint &p = g.modulus();

	base_vector< bigint > R(n, n);
	lidia_size_t i;
	for (i = 0; i < n; i++)
		R[i].assign(randomize(p));

	min_poly_work(h, R, m, g, F);
}



//
// Task:        Same as above, but guarantees that result is correct
//

void min_poly(Fp_polynomial& h, const Fp_polynomial& g, lidia_size_t m,
	      const Fp_poly_modulus& F)
{
	debug_handler("Fp_polynomial", "min_poly(Fp_polynomial&, Fp_polynomial&, lidia_size_t, Fp_poly_modulus&)");

	g.comp_modulus(F.modulus(), "min_poly");

	Fp_polynomial h1;

	for (;;) {
		prob_min_poly(h, g, m, F);

		if (h.degree() == m)
			return;
		compose(h1, h, g, F);

		if (h1.is_zero())
			return;
	}
}



//
// Task:        Same as above, but assumes that f is irreducible
//
// Conditions:  Assumes that f is irreducible, or at least that the minimal
//              polynomial of g is itself irreducible.
//
// Algorithm:   The algorithm is deterministic (and is always correct).
//

void irred_poly(Fp_polynomial& h, const Fp_polynomial& g, lidia_size_t m,
		const Fp_poly_modulus& F)
{
	debug_handler("Fp_polynomial", "irred_poly(Fp_polynomial&, Fp_polynomial&, lidia_size_t, Fp_poly_modulus&)");

	g.comp_modulus(F.modulus(), "irred_poly");

	base_vector< bigint > R(lidia_size_t(1), lidia_size_t(1));

	R[0].assign_one();
	min_poly_work(h, R, m, g, F);
}



#if 0

//********************************************************
//
//		 Multiplication of several polynomials
//
//********************************************************


void multiply(Fp_polynomial& x, base_vector< Fp_polynomial > & a)
	// x = product of all a[i]'s, contents of a[i] are destroyed
	//***********************************
	// The algorithm works by successively taking the smallest two
	// elements from a, and replacing them by their product.
	// This calls for a priority-queue, of course, but the current
	// implementation uses a very dumb priority queue.
	// For long lists, a better priority-queue should perhaps be implemented.
	//************************************
{
	debug_handler("Fp_polynomial", "multiply(Fp_polynomial&, base_vector< Fp_polynomial > &)");

	lidia_size_t n = a.size();

	// first, deal with some trivial cases
	if (n == 0) {
		//WARNING: no modulus set for x !!!!!!!!!!!!!!
		lidia_error_handler("compose.c", "multiply(Fp_polynomial&, "
				    "base_vector< Fp_polynomial > &)::vector is empty - "
				    "no modulus could be set for result");
		return;
	}

	if (n == 1) {
		x.assign(a[0]);
		return;
	}

	lidia_size_t i, j;

	for (i = 0; i < n; i++) {
		if ((a[i]).is_zero()) {
			x.assign(a[i]); // sets modulus, too
			return;
		}
	}

	// assume n > 1 and all a[i]'s are nonzero

	// sort into non-increasing degrees

	for (i = 1; i <= n - 1; i++)
		for (j = 0; j <= n - i - 1; j++)
			if (a[j].degree() < a[j+1].degree())
				swap(a[j], a[j+1]);

	Fp_polynomial g;

	while (n > 1) {
		// replace smallest two poly's by their product
		multiply(g, a[n-2], a[n-1]);
		a[n-2].kill();
		a[n-1].kill();
		swap(g, a[n-2]);
		n--;

		// re-establish order
		i = n-1;
		while (i > 0 && a[i-1].degree() < a[i].degree()) {
			swap(a[i-1], a[i]);
			i--;
		}
	}

	x = a[0];

	a[0].kill();
	a.set_size(0);
}
#endif



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
