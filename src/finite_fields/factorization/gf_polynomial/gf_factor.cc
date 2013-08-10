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



void factor(factorization< gf_polynomial > &F, const gf_polynomial &f)
{
	const galois_field &K = f.get_field();
	const bigint &q = K.number_of_elements();
	const lidia_size_t k = K.degree();
	lidia_size_t n = f.degree(), i, j;
	if (n < 0)
		lidia_error_handler("gf_polynomial", "factor(...)::input is zero polynomial");

	if (k == 1) {
		Fp_polynomial ff;
		ff.set_modulus(q);
		for (i = n; i >= 0; i--)
			ff[i].assign(f[i].polynomial_rep()[0]);

		gf_polynomial tmp(K);
		factorization< Fp_polynomial > Fact = factor((single_factor< Fp_polynomial > )ff);
		F.kill();
		tmp.set_degree(0);
		tmp[0].set_polynomial_rep(Fact.unit());
		F.append(tmp);

		for (j = 0; j < Fact.no_of_prime_components(); j++) {
			tmp.set_degree(Fact.prime_base(j).base().degree());
			for (i = Fact.prime_base(j).base().degree(); i >= 0; i--) {
				ff.assign(Fact.prime_base(j).base()[i]);
				tmp[i].set_polynomial_rep(ff);
			}
			F.append(tmp, Fact.prime_exponent(j));
		}
		return;
	}
	if (q.bit_length() * n > 100000 || q < 50)
		can_zass(F, f);
	else
		berlekamp(F, f);
}



factorization< gf_polynomial >
single_factor< gf_polynomial >::factor() const
{
	return LiDIA::factor((single_factor< gf_polynomial > )rep);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
