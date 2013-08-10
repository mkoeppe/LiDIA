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
#include	"LiDIA/udigit.h"
#include	"LiDIA/factorization.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



static
void split(Fp_polynomial & f1, Fp_polynomial & g1, Fp_polynomial & f2,
	   Fp_polynomial & g2, const Fp_polynomial & f, const Fp_polynomial & g,
	   const base_vector< bigint > &roots, lidia_size_t lo, lidia_size_t mid)
{
	debug_handler("Fp_polynomial", "split(Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, base_vector< bigint > &, lidia_size_t, lidia_size_t)");

	f.comp_modulus(g, "split");

	lidia_size_t i, r = mid - lo + 1;
	base_vector< bigint > lroots(r, r);

	for (i = 0; i < r; i++)
		lroots[i].assign(roots[lo + i]);

	Fp_polynomial a, d;

	if (f.degree() == 0)
		a.assign_one(f.modulus());
	else {
		Fp_polynomial h;
		h.set_modulus(f); //must be set for 'build_from_roots()'
		h.build_from_roots(lroots);

		Fp_poly_modulus F(f);
		compose(a, h, g, F);
	}

	gcd(f1, a, f);
	divide(f2, f, f1);

	remainder(g1, g, f1);
	remainder(g2, g, f2);
}


static
void rec_find_irred_factors(factorization< Fp_polynomial > &factors,
			    const Fp_polynomial & f, const Fp_polynomial & g,
			    const base_vector< bigint > &roots, lidia_size_t lo, lidia_size_t hi)
{
	//only called by find_factors(), therefore no argument checking
	debug_handler("Fp_polynomial", "rec_find_irred_factors(factorization< Fp_polynomial > &, Fp_polynomial&, Fp_polynomial&, base_vector< bigint > &, lidia_size_t, lidia_size_t)");


	lidia_size_t r = hi - lo + 1;

	if (r == 0) return;

	if (r == 1) {
		append_irred_factor(factors, f);
		return;
	}

	Fp_polynomial f1, g1, f2, g2;

	lidia_size_t mid = (lo + hi) / 2;

	split(f1, g1, f2, g2, f, g, roots, lo, mid);

	rec_find_irred_factors(factors, f1, g1, roots, lo, mid);
	rec_find_irred_factors(factors, f2, g2, roots, mid + 1, hi);
}



void
find_irred_factors(factorization< Fp_polynomial > &factors,
		   const Fp_polynomial & f, const Fp_polynomial & g,
		   const base_vector< bigint > &roots)
{
	debug_handler("Fp_polynomial", "find_irred_factors(factorization< Fp_polynomial > &, Fp_polynomial&, Fp_polynomial&, base_vector< bigint > &)");

	f.comp_modulus(g, "find_irred_factors");

	factors.kill();
	lidia_size_t r = roots.size();
	rec_find_irred_factors(factors, f, g, roots, 0, r - 1);

}



#if 0
void
iter_find_factors(factorization< Fp_polynomial > &factors,
		  const Fp_polynomial & f, const Fp_polynomial & g,
		  const base_vector< bigint > &roots)
{
	debug_handler("Fp_polynomial", "iter_find_factors(factorization< Fp_polynomial > &, Fp_polynomial&, Fp_polynomial&, base_vector< bigint > &)");

	f.comp_modulus(g, "iter_find_factors");

	lidia_size_t r = roots.size();
	lidia_size_t i;
	Fp_polynomial h, d;

	factors.kill();
	for (i = 0; i < r; i++) {
		subtract(h, g, roots[i]);
		gcd(d, f, h);
		factors.append(d);
	}
}
#endif



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
