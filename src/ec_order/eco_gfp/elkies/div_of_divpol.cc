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
//	Author	: Frank Lehmann, Markus Maurer, Volker Mueller, Denis Pochuev
//	Changes	: See CVS log
//
//==============================================================================================


// Assumes, that the prime is an elkies prime
// and ftau[0], a root of the modular polynomial,
// has been computed. Computes Elkies polynomial.



#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/eco_prime.h"
#ifdef DEBUG
#include	<cassert>
#endif

#include <iostream>

#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



void eco_prime::compute_divisor_of_division_polynomial (
	ff_pol & div_of_div_pol,
	ff_element & E2a,
	ff_element & E2b)
{
	div_of_div_pol.set_modulus(p);
	ff_element P1;
	base_vector< ff_element > jltau;

	ff_polmod div_of_div_pol_pm; // poly modulus for divisor of
	ff_pol  div_pol; // division polynomial

	bool divisor_found = false;
	lidia_size_t i;

	switch (this->trans_type()) {
	case 1:   // easy case which is deterministic

	  this->tildeEf (E2a, E2b, P1, ftau[0], A, B);
	  this->coefficient_comparison (div_of_div_pol, E2a, E2b, P1);

#ifdef DEBUG
	  div_of_div_pol_pm.build (div_of_div_pol);
	  this->compute_psi (div_pol, static_cast<lidia_size_t>(l), div_of_div_pol_pm);
	  assert(div_pol.is_zero());
#endif
	  break;

	case 2: // compute a vector of possible values of j(l*tau)

		this->guess_jltau (jltau);

		divisor_found = false;

		for (i = 0; i < jltau.get_size() && !divisor_found; i++) {
#ifdef DEBUG
			std::cout << "Testing for j(l*tau): " << i << std::endl;
#endif

			// compute candidate for divisor

			this->tildeEA (E2a, E2b, P1, ftau[0], jltau[i],
				       A, B);

			this->coefficient_comparison (div_of_div_pol, E2a, E2b, P1);

			div_of_div_pol_pm.build (div_of_div_pol);
			this->compute_psi (div_pol, static_cast<lidia_size_t>(l), div_of_div_pol_pm);
			if (div_pol.is_zero())
			  {
			    divisor_found = true;
			  }
		}
		if (!divisor_found)
			lidia_error_handler("eco_prime", "compute_divisor_of_division_polynomial::"
					    "div_of_divpol::no Elkies polynomial found");
		break;

	default:
		lidia_error_handler ("eco_prime", "compute_divisor_of_division_polynomial::"
				     "Invalid transformation type.");
		break;
	}
}



//----------------------------------------------------------------
// determine all candidates for j(ltau) and fill vector
// j(l*tau) is a root of Phi_l(ftau[0], Y) \in Fq[Y]
//


void eco_prime
::guess_jltau (base_vector< ff_element > & jltau)
{
	ff_pol g;
	lidia_size_t i;
	base_vector< bigint > root_vector;

	this->build_poly_in_Y (g, ftau[0]);
	root_vector = find_roots(g);

	jltau.set_capacity(root_vector.get_size());
	for (i = 0; i < root_vector.get_size(); i++) {
		jltau[i] = root_vector[i];
#ifdef DEBUG
		assert(g(jltau[i].mantissa()).is_zero());
#endif
	}
}



//-----------------------------------------------------------
// computes a divisor of the l-th division polynomial knowing
// the isogenous curve (E2a, E2b) and P1
//

void eco_prime
::coefficient_comparison (ff_pol & divisor,
			  const ff_element & E2a,
			  const ff_element & E2b,
			  const ff_element & P1)
{
  lidia_size_t ll = l;
  base_vector< ff_element > weierstrass_series_coefficients;
  base_vector< ff_element > c_k, c_k_tilde;
  dense_power_series< ff_element > weierstrass_series, rhs_series;  // rhs := right hand side
  dense_power_series< ff_element > tmp_weierstrass, weierstrass_to_tmp_degree;

  compute_coefficients_of_weierstrass_p(c_k, ll-3, A, B);
  compute_coefficients_of_weierstrass_p(c_k_tilde, ll-3, E2a, E2b);
  compute_coeffseries(rhs_series, c_k, c_k_tilde, P1, 0);
  compute_weierstrass_p(weierstrass_series, c_k, ll-3);

  lidia_size_t d = (ll-1)/2, tmp_degree;
  ff_element tmp_divisor_coeff;
  dense_power_series< ff_element > *table_ptr;

  precompute_table(table_ptr, weierstrass_series, d);

  for (tmp_degree = d; tmp_degree >= 0; tmp_degree--) {

  

    tmp_divisor_coeff = rhs_series[-2*tmp_degree];
    weierstrass_to_tmp_degree = *(table_ptr + static_cast<int>(tmp_degree));
    multiply(tmp_weierstrass, weierstrass_to_tmp_degree, tmp_divisor_coeff);
    subtract(rhs_series, rhs_series, tmp_weierstrass);
    
    divisor.set_coefficient(tmp_divisor_coeff.mantissa(), tmp_degree);
  }
  delete[] table_ptr;
}



// ----------------------------------------------------------------
// Computes the coefficients c_1, ..., c_kmax of the Weierstrass
// P-function for curve (a, b) over GF(p) in wp_coeff.
//

void eco_prime
::compute_coefficients_of_weierstrass_p (base_vector< ff_element > & wp_coeff,
					 lidia_size_t kmax,
					 const ff_element & a,
					 const ff_element & b)
{
	if (kmax == 0)
		kmax = 1;
	wp_coeff.set_capacity(kmax+1);

	if (kmax == 1) {
		divide(wp_coeff[1], a, ff_element(5));
		wp_coeff[1].negate();
	}
	else
		if (kmax == 2) {
			divide(wp_coeff[1], a, ff_element(5));
			wp_coeff[1].negate();

			divide(wp_coeff[2], b, ff_element(7));
			wp_coeff[2].negate();
		}
		else {
			divide(wp_coeff[1], a, ff_element(5));
			wp_coeff[1].negate();

			divide(wp_coeff[2], b, ff_element(7));
			wp_coeff[2].negate();

			ff_element tmp_sum, fraction, tmp_wp_coeff;
			lidia_size_t denominator, k, h;

			for (k = 3; k < kmax; k++) {
				tmp_sum = 0;
				denominator = (k-2)*(2*k+3);
				divide(fraction, ff_element(3), denominator);

				for (h = 1; h <= k-2; h++) {
					multiply(tmp_wp_coeff, wp_coeff[h], wp_coeff[k-1-h]);
					add(tmp_sum, tmp_sum, tmp_wp_coeff);
				}
				multiply(wp_coeff[k], fraction, tmp_sum);
			}
		}
}



//------------------------------------------------------------------
// Compute the Weierstrass p function wp up to z^M
// using the precomputed coefficients in wp_coeff.
//

void eco_prime
::compute_weierstrass_p (dense_power_series< ff_element > & wp,
			 const base_vector< ff_element > & wp_coeff,
			 lidia_size_t M)
{

	wp.assign_one(M+2);
	wp.multiply_by_xn(-2);
	// wp(z, L) = 1/z^2 + SUM c(k)*z^(2k) {k=1..M}
	// setting coefficient of z^(-2) to 1

	lidia_size_t k;

	if (M >= 2) {
		for (k = 1; k <= M/2; k++)
			wp.set_coeff(wp_coeff[k], 2*k);
	}
}



//--------------------------------------------------------------
// Compute the power-series of the right hand side in r
// up to z^M using the precomputed coefficients in c1 and
// c2.
//

void eco_prime
::compute_coeffseries (dense_power_series< ff_element > & right_hand_side_series,
		       const base_vector< ff_element > & c1,
		       const base_vector< ff_element > & c2,
		       const ff_element & P1, lidia_size_t M)
{

  dense_power_series< ff_element > s_tmp, s_to_the_power_r, s_tmp_sum;
  
  lidia_size_t 	k, tmp_denominator, r, r_fact = 1;
  ff_element	tmp_numerator, tmp_coeff;
  
  // Coefficients used in the serieses:
  // s_tmp, s_to_the_power_r:	[2, ..., M+l-1]

  for (k = 1; k <= static_cast<int>(M+l-3)/2; k++) {
    // s_tmp_sum:			[0, ..., M+l-1]
    tmp_denominator = (2*k+1)*(2*k+2); 	// right_hand_side_series: [-(l-1), .., M]
    multiply(tmp_numerator, ff_element(l), c1[k]);
    subtract(tmp_numerator, tmp_numerator, c2[k]);
    divide(tmp_coeff, tmp_numerator, tmp_denominator);
    s_tmp.set_coeff(tmp_coeff, 2*k+2);
  }
  
  s_tmp.set_coeff(-P1, 2); // series inside the parentheses (s_tmp) is computed

  s_to_the_power_r = s_tmp;
  s_tmp_sum.assign_one(M+l-1); 	//1 is assigned to the sum (case r = 0).
  // The loop starts with r = 1
  for (r = 1; r <= static_cast<int>(M+l-1)/2; r++) {
    r_fact *= r;
    divide(s_to_the_power_r, s_to_the_power_r, ff_element(r));
    add(s_tmp_sum, s_tmp_sum, s_to_the_power_r);
    multiply(s_to_the_power_r, s_to_the_power_r, s_tmp);
  }

  s_tmp_sum.multiply_by_xn(1-l);
  right_hand_side_series = s_tmp_sum;
}



void eco_prime::precompute_table (dense_power_series< ff_element >* &weierstrass_table,
				  const dense_power_series< ff_element > &weierstrass_series,
				  lidia_size_t d)
{
	lidia_size_t i;

	weierstrass_table = new dense_power_series< ff_element >[d+1];
	weierstrass_table[0].assign_one(weierstrass_series.get_last());
	weierstrass_table[1] = weierstrass_series;

	for (i = 2; i <= d; i++) {
		// could be improved.
		multiply(weierstrass_table[i], weierstrass_table[i-1], weierstrass_series);
	}
}



void eco_prime::set_odd_weierstrass (dense_power_series< ff_element >* const table,
				     lidia_size_t odd_number)
{
	dense_power_series< ff_element > tmp_series;
	square(tmp_series, table[odd_number/2] + table[odd_number/2 + 1]);

	tmp_series -= table[odd_number-1];
	tmp_series -= table[odd_number+1];
	divide(table[odd_number], tmp_series, ff_element(2));
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
