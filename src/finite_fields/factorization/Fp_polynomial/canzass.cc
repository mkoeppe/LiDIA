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



#define LIDIA_DDF ddf
//#define LIDIA_DDF old_ddf
//you can choose between simple 'old_ddf' or faster 'ddf' - routine
//see sf_can_zass_work



//***********************************************************************
//
//			the main routine
//
//***********************************************************************

//
// Literature:	J. von~zur~Gathen and V. Shoup:
//		Computing Frobenius Maps and factoring polynomials,
//		Computational Complexity, 2:187-224, 1992
//
//		V. Shoup
//		A New Polynomial Factorization Algorithm and its Implementation
//		Universitaet des Saarlandes, 1994
//
//		D. G. Cantor and H. Zassenhaus
//		A New Algorithm for Factoring Polynomials over Finite Fields",
//		Math. Comp. 36, 154, 1981, p.587-592
//
//		M. O. Rabin
//		Probabilistic algorithms in finite fields
//		SIAM J. Comput., 1980, 9, p.273-280
//
//		and many more...
//
// Task:        Computes factorization of f = F.modulus()
//
// Conditions:  f is square-free, monic, deg(f) > 0 and f(0) != 0.
//              x_to_the_p = x^p mod f
//
// Algorithm:   DDF: Shoup              if DDF_FLAG = ddf
//		     Cantor/Zassenhaus  if DDF_FLAG = old_ddf
//		EDF: Shoup/von zur Gathen
//

void
sf_can_zass_work(factorization< Fp_polynomial > &factors,
		 const Fp_polynomial &x_to_the_p, const Fp_poly_modulus &F)
{
	debug_handler("Fp_polynomial", "sf_can_zass_work(factorization< Fp_polynomial > &, Fp_polynomial& f, Fp_poly_modulus&, lidia_size_t)");

	const Fp_polynomial &f = F.modulus();

	factors.kill();

	my_timer t;
	bool verbose = single_factor< Fp_polynomial >::verbose();

	factorization< Fp_polynomial > ddf_fact, edf_fact;
	if (verbose) {
		std::cerr << "computing ddf...\n";
		t.start("ddf time: ");
	}


	LIDIA_DDF(ddf_fact, f, x_to_the_p);
	//the exponents of 'ddf_fact' are the degrees of the irred. polynomials !!!


	if (verbose) t.stop();

	Fp_polynomial hh;
	Fp_poly_modulus FF;

	lidia_size_t i;
	for (i = 0; i < ddf_fact.no_of_composite_components(); i++) {
		const Fp_polynomial & g = ddf_fact.composite_base(i).base();
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
				remainder(hh, x_to_the_p, g);
				edf(edf_fact, FF, hh, d);
			}
			multiply(factors, factors, edf_fact);
		}
	}
}



//***********************************************************************
//
//			"interface"
//
//***********************************************************************

void
sf_can_zass(factorization< Fp_polynomial > &factors, const Fp_polynomial & f)
{
	bool verbose = single_factor< Fp_polynomial >::verbose();
	my_timer t;
	Fp_polynomial g;
	Fp_poly_modulus F(f);

	if (verbose) t.start("computing x^p...");
	power_x(g, f.modulus(), F);
	if (verbose) t.stop();

	sf_can_zass_work(factors, g, F);
}



void
can_zass(factorization< Fp_polynomial > &factors, const Fp_polynomial & f)
{
	debug_handler("Fp_polynomial", "can_zass(factorization< Fp_polynomial > &, Fp_polynomial&)");

	factor_generic(factors, f, sf_can_zass);
}



factorization< Fp_polynomial > can_zass(const Fp_polynomial & f)
{
	factorization< Fp_polynomial > F;
	can_zass(F, f);
	return F;
}



factorization< Fp_polynomial >
single_factor< Fp_polynomial >::can_zass() const
{
	factorization< Fp_polynomial > F;
	LiDIA::can_zass(F, rep);
	return F;
}



factorization< Fp_polynomial >
single_factor< Fp_polynomial >::sf_can_zass() const
{
	factorization< Fp_polynomial > F;
	LiDIA::sf_can_zass(F, rep);
	return F;
}



factorization< Fp_polynomial > sf_can_zass(const Fp_polynomial & f)
{
	factorization< Fp_polynomial > F;
	sf_can_zass(F, f);
	return F;
}


#undef LIDIA_OLD_DDF
#undef LIDIA_DDF



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
