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
#include	"LiDIA/eco_gf2n.h"
#include	"LiDIA/sort_vector.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//-----------------------------------------------------------------
// Assumes, that the prime l is an elkies prime
// and ftau[0], a root of the modular polynomial,
// has been computed. Computes Elkies polynomial
// via an algorithm of Lercier.
//
//
//  The abszisse I of the l-isogeny has the form
//            X P^2(X)
//    I(X) = ----------, where Q(X) is a divisor of the division polynomial
//              Q^2(X)    of degree l-1/2. The Elkies-polynomial.
//
//  This function returns Q(X).

void eco_gf2n::divisor_of_divpol_via_lercier (ff_pol & Q,
				    	      char strategy)
{
	ff_element Eb; // the curve E/C corresponding to l-isogeny

	// j-invariant of curve Eb is root ftau[0] of modular polynomial

	if (!ftau[0].is_zero())
		invert(Eb, ftau[0]);
	else lidia_error_handler("divisor_of_div_poly_via_lercier",
				 "Root of modular polynomial is zero");

	lidia_size_t d = (l-1) / 2; // degree of divisor of division polynomial
	ff_element * p_k;

	p_k = new ff_element[d+1]; // coeff. of P(X) = \sum_{i=0}^{d} p_i^2 X^i

	// compute alpha and beta, sqrt[4] (alpha), sqrt[4](beta) constants

	ff_element sqrt4_alpha, sqrt4_beta;
	ff_element alpha, beta;

	sqrt(alpha, sqrt(A6));
	sqrt(beta, sqrt(Eb));
	sqrt(sqrt4_alpha, sqrt(alpha));
	sqrt(sqrt4_beta, sqrt(beta));

	power_tab< ff_element > alpha_pow;
	alpha_pow.initialize(alpha, 2*d);

	// call Lerciers algorithm that returns p_i, 0<=i<=d

	compute_sqrtP(p_k, alpha_pow, sqrt4_alpha, sqrt4_beta, beta, d, strategy);

	// calculate Coefficients of Q(X) and
	// assign values to the polynomial Q(X) = \sum_i=0^d q[i]^2 X^i

	gf2n temp, temp2, temp_pow;
	lidia_size_t i;
	ff_element qi; // coefficients of fC

	divide(temp, sqrt4_alpha, sqrt4_beta);
	Q.assign(gf2n_polynomial(d));

	for (i = 0; i < d; i++) {
		if(d > 2*i) {
			alpha_pow.get_power(temp_pow, d-2*i);
			temp_pow = sqrt(temp_pow);
		}
		else {
			alpha_pow.get_power(temp_pow, d);
			temp_pow = sqrt(temp_pow);

			alpha_pow.get_power(temp2, i);
			divide(temp_pow, temp_pow, temp2);
		}
		multiply(qi, temp_pow, p_k[d-i]);
		multiply(qi, qi, temp);
		square(qi, qi);
		Q.set_coefficient(qi, i);
	}
	Q.set_coefficient(d);
	delete[] p_k;
	mv_poly::delete_free_list();
}



//-----------------------------------------------------------------------
// In the second version we compute the numerator and the denominator
// of the isogeny. Apart from that the function is the same as the first
// function.
//
//    The abszisse I of the l-isogeny from E(A6_1) to E(A6_2) has the form
//            X P^2(X)
//    I(X) = ----------, where Q(X) is a divisor of the division polynomial
//              Q^2(X)    of degree l-1/2. The Elkies-polynomial.
//
//    This function returns numerator = X P^2(X)
//                      and denominator = Q^2(X)
//

void eco_gf2n::isogeny_via_lercier(ff_pol & numerator, ff_pol & denominator,
			           const ff_element & A6_1,
                                   const ff_element & A6_2,
			           char strategy)
{
	lidia_size_t d = (l-1) / 2; // degree of divisor of division polynomial
	ff_element * p_k;
	p_k = new ff_element[d+1]; // coeff of P(X) = \sum_{i=0}^{d} p_i^2 X^i

	// compute alpha and beta constants
	ff_element alpha, beta, sqrt4_alpha, sqrt4_beta;

	sqrt(alpha, sqrt(A6_1));
	sqrt(beta, sqrt(A6_2));
	sqrt(sqrt4_alpha, sqrt(alpha));
	sqrt(sqrt4_beta, sqrt(beta));

	power_tab< ff_element > alpha_pow;
	alpha_pow.initialize(alpha, 2*d);

	// call Lerciers algorithm that returns p_i, 0<=i<=d

	compute_sqrtP(p_k, alpha_pow, sqrt4_alpha, sqrt4_beta, beta, d, strategy);

	// calculate denominator = Q^2(x)
	// assign values to the polynomial Q^2(X) = \sum_i=0^d q[i]^4 X^{2i}

	ff_element temp, temp2, temp_pow;
	power_tab< ff_element > alpha_inv_pow;

	invert(temp, alpha);
	alpha_inv_pow.initialize(temp, d);
	divide(temp, sqrt4_alpha, sqrt4_beta);

	lidia_size_t i;
	ff_element qi; // coefficients of fC

	denominator.set_size(l);

	for (i = 0; i < d; i++) {
		if(d > 2*i) {
			alpha_pow.get_power(temp_pow, d-2*i);
			temp_pow = sqrt(temp_pow);
		}
		else {
			alpha_pow.get_power(temp_pow, d);
			temp_pow = sqrt(temp_pow);
			alpha_inv_pow.get_power(temp2, i);
			multiply(temp_pow, temp_pow, temp2);
		}
		multiply(qi, temp_pow, p_k[d-i]);
		multiply(qi, qi, temp);
		square(qi, qi);
		square(qi, qi);
		denominator.set_coefficient(qi, 2*i);
	}
	denominator.set_coefficient(2*d);

	// compute numerator = X P^2(X)

	numerator.set_size(l+2);
	for (i = 0; i < d; i++) {
		square(p_k[i], p_k[i]);
		square(p_k[i], p_k[i]);
		numerator.set_coefficient(p_k[i], (2*i) + 1);
	}
	numerator.set_coefficient(l);

	delete[] p_k;
	mv_poly::delete_free_list();
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
