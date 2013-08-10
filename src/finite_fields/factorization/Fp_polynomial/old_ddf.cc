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
//	Author	:
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/Fp_polynomial.h"
#include	"LiDIA/finite_fields/Fp_polynomial_util.h"
#include	"LiDIA/Fp_poly_modulus.h"

#include	"LiDIA/factorization.h"
#include	"LiDIA/single_factor.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



int OLD_DDF_GCD_BLOCKING_FACTOR = 10;



//
// Task:        appends g, a product of irreducible polynomials of degree m,
//              to u
//              (the exponents of u are the degrees of the irreducible
//               polynomials !)
//

static
void add_factor(factorization< Fp_polynomial > &factors,
		const Fp_polynomial & g, lidia_size_t d)
{
	debug_handler("Fp_polynomial", "add_factor(factorization< Fp_polynomial > &, Fp_polynomial&, lidia_size_t)");

	if (single_factor< Fp_polynomial >::verbose())
		std::cerr << "degree = " << d << ", number = " << g.degree() / d << "\n";
	factors.append(g, d);
}



//
// Task:        splits f by computing gcd(f, buf[i])
//
// Algorithm:   instead of computing gcd(f,buf[i]) for all i, we compute
//              compute t1 = product of tbl[i] for i=0,...,limit-1.
//		If gcd(f,g) = 1, we do not have to compute any further gcd
//		and can return immediately
//

static
void process_table(Fp_polynomial & f, factorization< Fp_polynomial > &factors,
		   const Fp_poly_modulus & F, lidia_size_t limit,
		   const base_vector< Fp_polynomial > &tbl, lidia_size_t d)
{
	debug_handler("Fp_polynomial", "process_table(Fp_polynomial&, factorization< Fp_polynomial > &, Fp_poly_modulus&, lidia_size_t, base_vector< Fp_polynomial > &, lidia_size_t)");

	if (limit == 0)
		return;

	bool verbose = single_factor< Fp_polynomial >::verbose();
	if (verbose) std::cerr << "+";

	Fp_polynomial t1;

	if (limit == 1) {
		gcd(t1, f, tbl[0]);
		if (t1.degree() > 0) {
			add_factor(factors, t1, d);
			divide(f, f, t1);
		}

		return;
	}

	lidia_size_t i;

	t1 = tbl[0];
	for (i = 1; i < limit; i++)
		multiply(t1, t1, tbl[i], F); //Fp_poly_modulus

	gcd(t1, f, t1);

	if (t1.degree() == 0)
		return;

	divide(f, f, t1);

	Fp_polynomial t2;

	i = 0;
	d = d - limit + 1;

	while (2 * d <= t1.degree()) {
		gcd(t2, tbl[i], t1);
		if (t2.degree() > 0) {
			add_factor(factors, t2, d);
			divide(t1, t1, t2);
		}

		i++;
		d++;
	}

	if (t1.degree() > 0)
		add_factor(factors, t1, t1.degree());
}




//*************************************************************************
//
//			the main routine
//
//*************************************************************************

//
// Task:        Performs the distinct degree factorization of f = F.modulus()
//              (i.e. splits f in products of irreducibles of the same
//              degree)
//              when finished, the exponents of 'factors' are the degrees of
//              the irreducible factors of 'f' !!!
//
// Conditions:  f is square-free, monic, deg(f) > 0 and f(0) != 0 (?).
//              h = x^p mod f
//
// Algorithm:   DDF
//

void
old_ddf(factorization< Fp_polynomial > &factors,
	const Fp_polynomial & ff, const Fp_polynomial & hh)
	// Performs distinct-degree factorization.
	// f must be monic, square-free
	// h  = x^p mod f
	// the exponents of 'factors' show the degrees of the irred. polynomials
{
	debug_handler("Fp_polynomial", "old_ddf(factorization< Fp_polynomial > &, Fp_polynomial&, Fp_polynomial&, lidia_size_t)");

	ff.comp_modulus(hh, "old_ddf");

	factors.kill();

	if (ff.degree() == 1) {
		add_factor(factors, ff, 1);
		return;
	}

	Fp_polynomial f(ff);
	Fp_polynomial h(hh);

	lidia_size_t CompTableSize = 2 * square_root(f.degree());
	lidia_size_t gcdTableSize = OLD_DDF_GCD_BLOCKING_FACTOR;

	Fp_poly_modulus F(f);

	poly_argument H;
	H.build(h, F, comparator< lidia_size_t >::min(CompTableSize, f.degree()));

	lidia_size_t i, d, limit, old_n;
	Fp_polynomial g, x;


	base_vector< Fp_polynomial > tbl(gcdTableSize, gcdTableSize);

	x.set_modulus(f);
	x.assign_x();

	i = 0;
	g = h;
	d = 1;
	limit = gcdTableSize;


	while (2 * d <= f.degree()) {
		old_n = f.degree();
		subtract(tbl[i], g, x);
		i++;
		//instead of computing gcd(f, x^(p^i)-x) for all i separately,
		//we call 'process_table' for a block of 'limit' polynomials of
		//the form x^(p^i)-x :
		if (i == limit) {
			process_table(f, factors, F, i, tbl, d);
			i = 0;
		}

		d++;
		if (2 * d <= f.degree()) {
			// we need to go further

			if (f.degree() < old_n) {
				// f has changed

				F.build(f);
				remainder(h, h, F);
				remainder(g, g, F);
				H.build(h, F, comparator< lidia_size_t >::
					min(CompTableSize, f.degree()));
			}

			H.compose(g, g, F);
		}
	}

	//call 'process_table' for the remaining block of polynomials x^(p^i)-x
	process_table(f, factors, F, i, tbl, d - 1);

	if (!f.is_one())
		add_factor(factors, f, f.degree());
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
