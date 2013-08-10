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
#include	"LiDIA/gf_polynomial.h"
#include	"LiDIA/single_factor.h"
#include	"LiDIA/factorization.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



void append_irred_factor(factorization< gf_polynomial > &F, const gf_polynomial &f, lidia_size_t e)
{
	single_factor< gf_polynomial > tmp(f);
	tmp.set_prime_flag(decomposable_object::prime);
	F.append(tmp, e);
}



void
find_irred_factors(factorization< gf_polynomial > &factors,
		   const gf_polynomial & f, const gf_polynomial & g,
		   const base_vector< gf_element > &roots)
	//assumes f = prod_{i=0}^{r-1} gcd(f,g-roots[i]) is a compl. factorization of f
{
	debug_handler("gf_polynomial", "find_irred_factors(factorization< gf_polynomial > &, gf_polynomial&, gf_polynomial&, base_vector< gf_element > &)");

	lidia_size_t r = roots.size();
	lidia_size_t i;
	gf_polynomial h, d;

	factors.kill();
	for (i = 0; i < r; i++) {
		subtract(h, g, roots[i]);
		gcd(d, f, h);
		append_irred_factor(factors, d, 1);
//may be faster: divide(f,f,d); remainder(g,g,f);
	}
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
