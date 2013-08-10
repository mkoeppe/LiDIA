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
#include	"LiDIA/factorization.h"
#include	"LiDIA/finite_fields/Fp_polynomial_fft.h"
#include	"LiDIA/finite_fields/Fp_polynomial_util.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



void add_factor(factorization< gf_polynomial > & factors, const gf_polynomial &f,
		lidia_size_t d)
{
	debug_handler("gf_polynomial", "add_factor(...)");

	if (single_factor< gf_polynomial >::verbose())
		std::cerr << "degree = " << d << ", number = " << f.degree() / d << "\n";
	factors.append(f, d);
}



void tbl_gcd_refine(gf_polynomial &f, factorization< gf_polynomial > &factors,
		    const gf_poly_modulus &F, int limit, const gf_polynomial *tbl,
		    lidia_size_t d)
{
	debug_handler("gf_polynomial", "tbl_gcd_refine(...)");

	if (limit == 0) return;

	bool verbose = single_factor< gf_polynomial >::verbose();
	if (verbose) std::cerr << "+";

	gf_polynomial t, s;
	lidia_size_t i;

	if (limit == 1) {
		gcd(t, f, tbl[0]);
		if (t.degree() > 0) {
			add_factor(factors, t, d);
			divide(f, f, t);
		}
		return;
	}

	t.assign(tbl[0]);

	for (i = 1; i < limit; i++)
		multiply(t, t, tbl[i], F);

	gcd(t, f, t);
	if (t.degree() == 0) return;

	divide(f, f, t);

	i = 0;
	d = d - limit + 1;
	while (2 * d <= t.degree()) {
		gcd(s, tbl[i], t);
		if (s.degree() > 0) {
			add_factor(factors, s, d);
			divide(t, t, s);
		}
		i++;
		d++;
	}
	if (t.degree() > 0)
		add_factor(factors, t, t.degree());
}



void ddf(factorization< gf_polynomial > &factors, const gf_polynomial &ff,
	 const gf_polynomial &X_pow_q)
	//X_pow_q = X^q mod ff
	// the exponents of 'factors' show the degrees of the irred. polynomials
{
	debug_handler("gf_polynomial", "ddf(factorization< gf_polynomial > &, gf_polynomial&, gf_polynomial&)");

	factors.kill();

	if (ff.degree() == 1) {
		add_factor(factors, ff, 1);
		return;
	}

	lidia_size_t i;
	gf_poly_modulus F(ff);


	// build table of X_pow_q, X_pow_q^2, X_pow_q^3, ... X_pow_q^(l-1)  
	lidia_size_t L, l = 0;
	gf_poly_argument H;
	bool use_gf_poly_argument;
	const bigint &q = ff.get_field().number_of_elements();
	if (q > ff.degree()) {
		use_gf_poly_argument = true;
		L = 2 * square_root(ff.degree());
		l = comparator< lidia_size_t >::min(L, ff.degree());
		H.build(X_pow_q, F, l);
	}
	else
		use_gf_poly_argument = false;

	lidia_size_t d = 1; // degree of irred. factors currently tested  

	const int DDF_GCD_TABLE_SIZE = 4; // blocking size  
	gf_polynomial *tbl = new gf_polynomial[DDF_GCD_TABLE_SIZE];
	if (!tbl) lidia_error_handler("gf_polynomial", "ddf(...)::out of memory");
	i = 0; // index for tbl  

	gf_polynomial h(X_pow_q), g(X_pow_q), f(ff);
	gf_polynomial X(f.get_field());
	X.assign_x();
	lidia_size_t old_n;


	while (2*d <= f.degree()) {
		// while not all factors found  
		old_n = f.degree();
		subtract(tbl[i], g, X);
		i++;
		if (i == DDF_GCD_TABLE_SIZE) {
			// compute gcd and refine  
			tbl_gcd_refine(f, factors, F, DDF_GCD_TABLE_SIZE, tbl, d);
			i = 0;
		}

		d++;
		if (2*d <= f.degree()) {
			// compute X_pow_q^(q^d)  
			if (f.degree() < old_n) {
				F.build(f);
				remainder(h, h, F);
				remainder(g, g, F);
				if (use_gf_poly_argument) {
					l = comparator< lidia_size_t >::min(l, f.degree());
					H.build(h, F, l);
				}
			}
			if (use_gf_poly_argument)
				H.compose(g, g, F);
			else
				power(g, g, q, F);
		}
	}
	tbl_gcd_refine(f, factors, F, i, tbl, d - 1);

	delete[] tbl;
	if (f.degree() > 0)
		add_factor(factors, f, f.degree());
}



void
root_edf(factorization< gf_polynomial > &factors, const gf_polynomial & f)
	// edf for d==1
{
	debug_handler("gf_polynomial", "root_edf(...)");

	base_vector< gf_element > roots;
	my_timer t;
	bool verbose = single_factor< gf_polynomial >::verbose();

	if (verbose) t.start("finding roots...");
	roots.assign(find_roots(f));
	if (verbose) t.stop();

	lidia_size_t r = roots.size();

	gf_polynomial tmp(f.get_field());
	tmp.assign_x();

	factors.kill();
	for (lidia_size_t j = 0; j < r; j++) {
		negate(tmp[0], roots[j]); // tmp = X - roots[j]  
		append_irred_factor(factors, tmp);
	}
}



void edf(factorization< gf_polynomial > &factors, const gf_poly_modulus &F,
	 const gf_polynomial &X_pow_q, lidia_size_t d)
{
	debug_handler("gf_polynomial", "edf(...)");

	factors.kill();
	const gf_polynomial &f = F.modulus();
	lidia_size_t i, n = f.degree();
	lidia_size_t r = n / d;
	bool verbose = single_factor< gf_polynomial >::verbose();
	if (verbose) std::cerr << "edf for d = " << d << ", r = " << r << " :\n";

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

	my_timer t;
	if (verbose) t.start("computing edf...");

	const galois_field &K = f.get_field();
	bigint q2(K.number_of_elements());

	gf_polynomial a, g, h;
	gf_element ONE(K);
	ONE.assign_one();

	if ((q2 > 20) && (q2 > 2*r)) {
		do {
			a = randomize(K, n-1);
			trace_map(g, a, d, F, X_pow_q);
			if (verbose) std::cerr << "." << std::flush;
		} while (!checked_min_poly(h, g, r, F));
		base_vector< gf_element > roots;
		roots.assign(find_roots(h));
		find_irred_factors(factors, f, g, roots);
	}
	else {
		q2.divide_by_2(); // q2 = (q-1)/2  
		factorization< gf_polynomial > Fact;
		Fact.assign(f);
		do {
			a = randomize(K, n-1);
			trace_map(g, a, d, F, X_pow_q);
			if (Fact.refine(g) ? verbose : 0) std::cerr << "+";

			if (Fact.no_of_composite_components() < r) {
				if (K.characteristic() != 2) {
					power(g, g, q2, F);
					add(g, g, ONE);
				}
				else	// p == 2
					trace_map(g, a, K.degree(), F, X_pow_q);
				if (Fact.refine(g) ? verbose : 0) std::cerr << "+";
			}
		} while (Fact.no_of_composite_components() < r);

		for (i = 0; i < r; i++)
			append_irred_factor(factors, Fact.composite_base(i).base());
	}
}



void
sf_can_zass(factorization< gf_polynomial > &factors, const gf_polynomial & f)
	// Assumes f is square-free and monic
	// returns list of factors of f.
{
	debug_handler("gf_polynomial", "sf_can_zass(factorization< gf_polynomial > &, gf_polynomial& f)");

	factors.kill();

	my_timer t;
	bool verbose = single_factor< gf_polynomial >::verbose();

	const bigint & q = f.get_field().number_of_elements();

	gf_poly_modulus F(f);
	gf_polynomial h;

	if (verbose) t.start("computing X^q...");
	power_x(h, q, F);
	if (verbose) t.stop();

	factorization< gf_polynomial > ddf_fact, edf_fact;
	if (verbose) {
		std::cerr << "computing ddf...\n";
		t.start("ddf time: ");
	}

	ddf(ddf_fact, f, h);
	//the exponents of 'ddf_fact' are the degrees of the irred. polynomials

	if (verbose) t.stop();

	gf_polynomial hh;
	gf_poly_modulus FF;

	lidia_size_t i;
	for (i = 0; i < ddf_fact.no_of_composite_components(); i++) {
		const gf_polynomial & g = ddf_fact.composite_base(i).base();
		lidia_size_t d = ddf_fact.composite_exponent(i);
		lidia_size_t r = g.degree() / d;

		if (r == 1) {
			// g is already irreducible
			append_irred_factor(factors, g);
		}
		else {
			// must perform edf

			if (verbose)
				std::cerr << "edf, degree = " << d << ", number = " << r << std::endl;

			if (d == 1) {
				// root finding

				root_edf(edf_fact, g);
			}
			else {
				// general case

				FF.build(g);
				remainder(hh, h, g);
				edf(edf_fact, FF, hh, d);
			}
			multiply(factors, factors, edf_fact);
		}
	}
}



void
can_zass(factorization< gf_polynomial > &factors, const gf_polynomial & f)
	// returns a list of factors, with multiplicities.
{
	debug_handler("gf_polynomial", "can_zass(factorization< gf_polynomial > &, gf_polynomial&)");

	if (f.is_zero())
		lidia_error_handler("gf_polynomial", "can_zass(...)::input is zero polynomial");

	my_timer t;
	bool verbose = single_factor< gf_polynomial >::verbose();
	factorization< gf_polynomial > sfd, x;
	gf_polynomial ff;
	gf_element lc = f.lead_coeff();
	gf_element lc_inv = inverse(lc);
	multiply(ff, f, lc_inv);

	if (verbose) t.start("square-free decomposition...");
	square_free_decomp(sfd, ff);
	if (verbose) t.stop();

	lidia_size_t i;
	factors.kill();

	for (i = 0; i < sfd.no_of_composite_components(); i++) {
		if (verbose) {
			std::cerr << "factoring multiplicity " << sfd.composite_exponent(i);
			std::cerr << ", deg = " << (sfd.composite_base(i).base()).degree() << "\n";
		}

		sf_can_zass(x, sfd.composite_base(i).base());

		x.power(sfd.composite_exponent(i));
		multiply(factors, factors, x);
	}
	ff.set_degree(0);
	ff[0].assign(lc);
	factors.append(ff);
}



factorization< gf_polynomial > sf_can_zass(const gf_polynomial &f)
{
	factorization< gf_polynomial > factors;
	sf_can_zass(factors, f);
	return factors;
}



factorization< gf_polynomial > single_factor< gf_polynomial >::sf_can_zass() const
{
	return LiDIA::sf_can_zass(rep);
}



factorization< gf_polynomial > can_zass(const gf_polynomial &f)
{
	factorization< gf_polynomial > factors;
	can_zass(factors, f);
	return factors;
}



factorization< gf_polynomial > single_factor< gf_polynomial >::can_zass() const
{
	return LiDIA::can_zass(rep);
}



factorization< gf_polynomial > ddf(const gf_polynomial &f)
{
	factorization< gf_polynomial > ddf_fact;
	my_timer t;
	bool verbose = single_factor< gf_polynomial >::verbose();
	const bigint & q = f.get_field().number_of_elements();

	gf_poly_modulus F(f);
	gf_polynomial h;

	if (verbose) t.start("computing X^q...");
	power_x(h, q, F);
	if (verbose) t.stop();

	ddf(ddf_fact, f, h);
	return ddf_fact;
}



factorization< gf_polynomial > single_factor< gf_polynomial >::ddf() const
{
	return LiDIA::ddf(rep);
}



factorization< gf_polynomial > edf(const gf_polynomial &f, lidia_size_t d)
{
	factorization< gf_polynomial > edf_fact;
	my_timer t;
	bool verbose = single_factor< gf_polynomial >::verbose();
	const bigint & q = f.get_field().number_of_elements();

	gf_poly_modulus F(f);
	gf_polynomial h;

	if (verbose) t.start("computing X^q...");
	power_x(h, q, F);
	if (verbose) t.stop();

	edf(edf_fact, f, h, d);
	return edf_fact;
}



factorization< gf_polynomial > single_factor< gf_polynomial >::edf(lidia_size_t d) const
{
	return LiDIA::edf(rep, d);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
