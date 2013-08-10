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
#include	"LiDIA/factorization.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



void factor_generic(factorization< Fp_polynomial > &factors,
		    const Fp_polynomial &ff,
		    void (*do_work)(factorization< Fp_polynomial > &, const Fp_polynomial &))
{
	if (ff.is_zero()) {
		lidia_error_handler("Fp_polynomial", "factor(...)::"
				    "input is zero polynomial");
		return; // LC
	}

	bool verbose = single_factor< Fp_polynomial >::verbose();

	Fp_polynomial f(ff);
	factors.kill();

	my_timer t;
	lidia_size_t n = f.degree();
	if (verbose) std::cerr << "factor degree " << n << std::endl;

	if (n <= 1) {
		append_irred_factor(factors, ff); // trivial case
		return;
	}

	Fp_polynomial tmp;
	tmp.set_modulus(f);

	factorization< Fp_polynomial > G;

	if (!f.is_monic()) {
		tmp.assign(f.lead_coeff());
		G.append(tmp); // store leading coeff. in G
		f.make_monic();
	}


	// f is monic

	lidia_size_t i = 0;
	while (i < n && f[i] == 0)		// is f a power of 'x' ?
		i++;
	if (i > 0) {
		tmp.assign_x();
		append_irred_factor(G, tmp, i); // append 'x' to G

		shift_right(f, f, i); // eliminate powers of 'x'
		n -= i;
	}

	if (f.degree() == 2) {
		// factor monic quadr. polyn.
		factor_quadratic_pol(factors, f);
		multiply(factors, factors, G);
		return;
	}


	if (verbose) t.start("square-free decomposition...");
	square_free_decomp(factors, f); // make square-free
	if (verbose) t.stop();

	multiply(factors, factors, G);


	// factor components
	// components are square-free, monic, and have a nonzero constant term

	while (factors.no_of_composite_components() != 0) {
		if (verbose) {
			std::cerr << "factoring multiplicity " << factors.composite_exponent(i);
			std::cerr << ", deg = " << factors.composite_base(i).base().degree()
				  << std::endl;
		}

		do_work(G, factors.composite_base(0).base());

		G.power(factors.composite_exponent(0));
		factors.replace(0, G);
	}
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
