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

#include	"LiDIA/single_factor.h"
#include	"LiDIA/factorization.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//***********************************************************************
//
//			root_edf
//
//***********************************************************************


//
// Task:	compute EDF for d == 1
//
// Algorithm:	uses find_roots
//

void
root_edf(factorization< Fp_polynomial > &factors, const Fp_polynomial & f)
{
	debug_handler("Fp_polynomial", "root_edf(factorization< Fp_polynomial > &, Fp_polynomial&)");

	base_vector< bigint > roots;
	my_timer t;
	bool verbose = single_factor< Fp_polynomial >::verbose();

	if (verbose) t.start("finding roots...");
	rec_find_roots(roots, f);
	if (verbose) t.stop();

	lidia_size_t r = roots.size();

	Fp_polynomial tmp;
	tmp.set_modulus(f);
	tmp.assign_x();

	factors.kill();
	lidia_size_t j;
	for (j = 0; j < r; j++) {
		tmp.set_coefficient(-roots[j], 0); //tmp = x - roots[j]
		append_irred_factor(factors, tmp);
	}
}



//***********************************************************************
//
//		several methods for performing the edf
//
//***********************************************************************

static
void edf_minpoly_method(factorization< Fp_polynomial > &factors,
			const Fp_poly_modulus & F, const Fp_polynomial & b, lidia_size_t d)
{
	debug_handler("Fp_polynomial", "edf_minpoly_method(factorization< Fp_polynomial > &, Fp_poly_modulus&, Fp_polynomial&, lidia_size_t)");

	const Fp_polynomial & f = F.modulus();
	lidia_size_t n = f.degree();
	lidia_size_t r = n / d;

	Fp_polynomial a, g, h;
	const bigint &p = b.modulus();

	bool verbose = single_factor< Fp_polynomial >::verbose();
	my_timer t, tt;
	if (verbose) tt.start("computing edf...");
	do {
		randomize(a, p, n - 1);

		if (verbose) t.start("computing trace map...");
		trace_map(g, a, d, F, b);
		if (verbose) t.stop();

		if (verbose) t.start("computing splitting poly...");
		min_poly(h, g, r, F);
		if (verbose) t.stop();
	} while (h.degree() < r);

	base_vector< bigint > roots;

	if (verbose) t.start("finding roots...");
	rec_find_roots(roots, h);
	if (verbose) t.stop();

	if (verbose) t.start("finding factors...");
	find_irred_factors(factors, f, g, roots);
	if (verbose) t.stop();
	tt.stop();
}



static
void edf_dirpower_method(factorization< Fp_polynomial > &factors,
			 const Fp_poly_modulus & F, const Fp_polynomial & b, lidia_size_t d)
{
	debug_handler("Fp_polynomial", "edf_dirpower_method(factorization< Fp_polynomial > &, Fp_poly_modulus&, Fp_polynomial&, lidia_size_t)");

	const Fp_polynomial & f = F.modulus();
	lidia_size_t n = f.degree();
	lidia_size_t r = n / d;

	const bigint &p = b.modulus();
	bigint p2(p);
	p2.divide_by_2(); //p2 = (p-1)/2

	factorization< Fp_polynomial > Fact;
	Fact.append(f);
	Fp_polynomial a, g;

	bool verbose = single_factor< Fp_polynomial >::verbose();
	my_timer t;
	if (verbose) t.start("computing edf...");

	do {
		randomize(a, p, n-1);
		trace_map(g, a, d, F, b);
		if (Fact.refine(g) ? verbose : 0) std::cerr << "+";

		if (p != 2 && Fact.no_of_composite_components() < r) {
			power(g, g, p2, F);
			add(g, g, 1);
			if (Fact.refine(g) ? verbose : 0) std::cerr << "+";
		}
	} while (Fact.no_of_composite_components() < r);

	lidia_size_t i;
	for (i = 0; i < r; i++)
		append_irred_factor(factors, Fact.composite_base(i).base());

	if (verbose) t.stop();
}


//
// Literature:	J. von zur Gathen and V. Shoup
//		Computing Frobenius maps and factoring polynomials
//		Computational Complexity, 2, 1992, 187-224
//
// Task:	Performs equal-degree factorization.
//
// Conditions:	b = x^p mod f.
//		d = degree of irreducible factors of f
//
// Algorithm:	(Space for the trace-map computation can be controlled via
// 		ComposeBound.)
//

void
edf(factorization< Fp_polynomial > &factors, const Fp_poly_modulus & F,
    const Fp_polynomial & b, lidia_size_t d)
{
	debug_handler("Fp_polynomial", "edf(base_vector< Fp_polynomial > &, Fp_poly_modulus&, Fp_polynomial&, lidia_size_t)");

	const Fp_polynomial & f = F.modulus();
	const bigint & p = f.modulus();
	bool verbose = single_factor< Fp_polynomial >::verbose();

	f.comp_modulus(b, "edf");

	lidia_size_t n = f.degree();
	lidia_size_t r = n / d;

	if (verbose) std::cerr << "edf for d = " << d << ", r = " << r << " :\n";

	factors.kill();
	if (r == 0)
		return;

	if (r == 1) {
		append_irred_factor(factors, f);
		return;
	}

	if (d == 1) {
		root_edf(factors, f);
		return;
	}


	if ((p > 2*r) && (p > 20)) {
		edf_minpoly_method(factors, F, b, d);
	}
	else {
		// p <= 2*r || p <= 20
		edf_dirpower_method(factors, F, b, d);
	}
}



//***********************************************************************
//
//		special variants for edf
//
//***********************************************************************

//
// Task:	like edf, but returns just a single factor
//
// SHOULD REWRITE THIS FUNCTION
//

void
edf1(Fp_polynomial & factor, const Fp_poly_modulus & F, const Fp_polynomial & b,
     lidia_size_t d)
{
	debug_handler("Fp_polynomial", "edf1(Fp_polynomial&, Fp_poly_modulus&, Fp_polynomial&, lidia_size_t)");

	const Fp_polynomial & f = F.modulus();
	const bigint & p = f.modulus();
	bool verbose = single_factor< Fp_polynomial >::verbose();

	f.comp_modulus(b, "edf1");

	lidia_size_t n = f.degree();
	lidia_size_t r = n / d;
	my_timer t;

	if (r == 1) {
		factor = f;
		return;
	}

	bigint root;

	if (d == 1) {
		if (verbose) t.start("finding one root...");
		root.assign(find_root(f));
		if (verbose) t.stop();

		factor.set_modulus(f);
		factor.assign_x();
		subtract(factor, factor, root);
		return;
	}

	Fp_polynomial a;
	Fp_polynomial g, h;

	if (p > 2*r) {
		do {
			randomize(a, p, n - 1);

			if (verbose) t.start("computing trace map...");
			trace_map(g, a, d, F, b);
			if (verbose) t.stop();

			if (verbose) t.start("computing splitting poly...");
			min_poly(h, g, r, F);
			if (verbose) t.stop();

		} while (h.degree() < r);

		if (verbose) t.start("finding one root...");
		root.assign(find_root(h));
		if (verbose) t.stop();

		if (verbose) t.start("finding one factor...");
		subtract(g, g, root);
		gcd(factor, f, g);
		if (verbose) t.stop();
	}
	else {
		// p <= 2*r  
		if (verbose) t.start("computing a single factor");

		bigint p2(p);
		p2.divide_by_2(); //p2 = (p-1)/2

		for (;;) {
			randomize(a, p, n-1);
			trace_map(g, a, d, F, b);

			gcd(h, f, g);
			if (h.degree() == d) {
				factor = h;
				break;
			}
			if (h.degree() == n-d) {
				divide(factor, f, h);
				break;
			}

			if (p != 2) {
				power(g, g, p2, F);
				add(g, g, 1);
				gcd(h, f, g);
				if (h.degree() == d) {
					factor = h;
					break;
				}
				if (h.degree() == n-d) {
					divide(factor, f, h);
					break;
				}
			}
		}
		if (verbose) t.stop();
	}

}



//
// Conditions:	b = x^p mod f and p > n
//
// Algorithm:	this a probabilistic algorithm that computes a monic factor of
// 		f = F.modulus()
//		It never returns 1, and with very high probability (at least
//		for large p) it actually returns an irreducible factor.
//

void
mystery_edf(Fp_polynomial & factor, const Fp_poly_modulus & F,
	    const Fp_polynomial & b)
{
	debug_handler("Fp_polynomial", "mystery_edf(Fp_polynomial&, Fp_poly_modulus&, Fp_polynomial&)");

	my_timer t;
	bool verbose = single_factor< Fp_polynomial >::verbose();
	const Fp_polynomial & f = F.modulus();

	f.comp_modulus(b, "mystery_edf");

	lidia_size_t n = f.degree();

	bigint root;
	const bigint & p = f.modulus();

	if (b.is_x()) {
		if (verbose) t.start("finding one root...");
		root.assign(find_root(f));
		if (verbose) t.stop();

		factor.assign(b); // b == x !!!

		subtract(factor, factor, root);
		return;
	}

	Fp_polynomial a;
	Fp_polynomial g, h;

	randomize(a, p, n);

	if (verbose) t.start("computing trace map...");
	trace_map(g, a, n, F, b);
	if (verbose) t.stop();


	if (verbose) t.start("computing splitting poly...");
	min_poly(h, g, n / 2, F);
	if (verbose) t.stop();

	if (h.degree() == 0) {
		factor = f;
		return;
	}

	if (verbose) t.start("finding one root...");
	root.assign(find_root(h));
	if (verbose) t.stop();

	if (verbose) t.start("finding one factor...");
	subtract(g, g, root);
	gcd(factor, f, g);
	if (verbose) t.stop();

}



//***********************************************************************
//
//			    "interface"
//
//***********************************************************************

factorization< Fp_polynomial >
edf(const Fp_polynomial &f, lidia_size_t d)
{
	debug_handler("Fp_polynomial", "edf(Fp_polynomial&, lidia_size_t)");

	bool verbose = single_factor< Fp_polynomial >::verbose();
	my_timer t;

	factorization< Fp_polynomial > factors;
	Fp_polynomial b;
	Fp_poly_modulus F(f);

	if (verbose) t.start("computing x^p...");
	power_x(b, f.modulus(), F);
	if (verbose) t.stop();

	edf(factors, F, b, d);
	return factors;
}



factorization< Fp_polynomial >
single_factor< Fp_polynomial >::edf(lidia_size_t d) const
{
	debug_handler("single_factor< Fp_polynomial >", "edf(lidia_size_t)");

	return LiDIA::edf(rep, d);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
