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
//	Author	: Volker Mueller
//	Changes	: See CVS log
//
//==============================================================================================



#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/eco_gf2n.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif




//***************************************************************
//  compute_mod_2_power implements a new idea for computing c mod 2^k
//  for some k. The algorithm consists of two steps:
//
//  1. compute a maximal value s such that there is a 2^s torsion point P
//     in E(GF(2^n))
//  2. use x(P) to compute a polynomial g of degree 2^{k-s} that
//     is a divisor of the 2^k division polynomial
//
//  The idea is based on the idea of Menezes, but extended to s > 2.

void eco_gf2n::compute_mod_2_power()
{
	int i;
	bool rc = false;
	lidia_size_t alpha = 0;
	elliptic_curve< gf_element > e = eco_gf2n::generateCurve(A6);
	point< gf_element > P(e);
	point< gf_element > T(e);
	ff_element x, y;

	long l_old, alpha_old;
	char ev_strategy_old;
	ff_pol g;
	ff_pol fpol(2);

	sqrt(x, A6);
	fpol.set_coefficient(x, 0);
	sqrt(x, x);

	if (!solve_quadratic(y, x, x*x*x+A6)) {
		lidia_error_handler("eco_gf2n", "compute_mod_2_power::no 4-TP found");
		return;
	}

	P.assign(eco_gf2n::convertToGFElement(x),
		 eco_gf2n::convertToGFElement(y)); // this is a 4-TP
	l = 4;

	while (sqrt_of_point(T, P) > 0) {
		l <<= 1;
		P = T;
	}
	l_old = 2*l;
	alpha_old = (pn.least_significant_digit() - l + 1) % l_old;

#ifdef DEBUG
	assert((l*P).is_zero());
	assert(!((l/2)*P).is_zero());
#endif

	if (info) {
		std::cout << "\nTwo power torsion points yield directly: c = " << alpha_old;
		std::cout << " mod " << l_old << std::flush;

		if (pn / l_old < 4)  // already close enough, terminate
		{
			c.set_capacity (1);
			l = l_old;
			c[0] = alpha_old;
			if (info)
				std::cout << "\n==> c = " << c[0] << " mod " << l << std::flush;
			return;
		}
	}
	// now P is a l-TP
	g.assign_x();
	g.set_coefficient(eco_gf2n::convertToGF2n(P.get_x()), 0);

	while ((g.degree() < static_cast<int>(degree_max_2_power)) && (pn / l >= 4)) {
		compose_2_torsion_pol(g, g, fpol);
		l <<= 1;
	}

	// Now g is a divisor of the l-th division polynomial of degree
	// approx. max_2_power

#ifdef DEBUG
	ff_polmod fmod;
	ff_pol erg;

	fmod.build(g);
	compute_psi(erg, l, fmod);
	assert(erg.is_zero());
#endif

	lidia_size_t * klist; // compute list of possibilities for ev

	klist = new lidia_size_t [l / l_old + 1];

	klist[0] = l / l_old;
	klist[1] = alpha_old;

	for (i = 2; i <= static_cast<int>(l) / l_old; i++)
		klist[i] = klist[i-1] + l_old;

	ev_strategy_old = ev_strategy; //we always use tables and rat. functions
	ev_strategy = eco_gf2n::EV_RATIONAL_FUNCTION_TABLE;
	rc = schoofpart_rat_function(alpha, g, klist, false);
	ev_strategy = ev_strategy_old;

	if (rc) {
		c.set_capacity (1);
		if (alpha < 0)
			alpha += l;
		c[0] = alpha;
		if (info)
			std::cout << "\n==> c = " << c[0] << " mod " << l << std::flush;
	}
	else
		lidia_error_handler("eco_ff_element", "compute_mod_2_power::no "
				    "eigenvalue found");

	delete [] klist;

}



//**********************************************************************
// The function compose_2_torsion_pol computes
//
//         h = sqrt(f(x^2 + a6/x^2) * x^{2 * deg f})
//           = sqrt(f) (x + sqrt(a6)/x) * x^{deg f}
//
// g must be x^2 + sqrt(a6).

void eco_gf2n::compose_2_torsion_pol(ff_pol & h,
				     const ff_pol & f,
				     const ff_pol & g)
{
	ff_pol res, h1, h2;
	ff_element coeff;
	int j, k;

	res.assign_zero();
	h1.assign_one();
	k = f.degree();

	for (j = 0; j <= k; j++) {
		f.get_coefficient(coeff, j);
		sqrt(coeff, coeff);
		multiply_by_scalar(h2, coeff, h1);
		shift_left(h2, h2, k-j);
		add(res, res, h2);
		if (j < k)
			multiply(h1, h1, g);
	}
	h.assign(res);
}



//***************************************************************
// This function determines for input P all points Q with 2*Q = P and
// returns the number of such points. The result is returned
// in the array 'res' of suitable size.
//

int eco_gf2n::sqrt_of_point(point< gf_element > & res,
			    const point< gf_element > & P)
{
	ff_element root1, root2;
	ff_element c1, c2;

	sqrt(c1, eco_gf2n::convertToGF2n(P.get_x()));
	sqrt(c2, A6);

	if (!solve_quadratic(root1, c1, c2))
		return 0;

	if (root1.relative_degree() > A6.relative_degree())
		return 0;

	add(root2, root1, c1); // root1 and root 2 are the two the roots

	point< gf_element > Q(P.get_curve()), H(P.get_curve());

	if (solve_quadratic(c1, root1, root1*root1*root1 + A6))
		if (c1.relative_degree() <= A6.relative_degree()) {
			Q.assign(eco_gf2n::convertToGFElement(root1),
				 eco_gf2n::convertToGFElement(c1));
			multiply_by_2(H, Q);

			if (H == P) {
				res = Q;
				return 1;
			}

			add(c2, c1, root1);
			Q.assign(eco_gf2n::convertToGFElement(root1),
				 eco_gf2n::convertToGFElement(c2));
			multiply_by_2(H, Q);

			if (H == P) {
				res = Q;
				return 1;
			}
		}

	if (solve_quadratic(c1, root2, root2*root2*root2 + A6))
		if (c1.relative_degree() <= A6.relative_degree()) {
			Q.assign(eco_gf2n::convertToGFElement(root2),
				 eco_gf2n::convertToGFElement(c1));
			multiply_by_2(H, Q);

			if (H == P) {
				res = Q;
				return 1;
			}

			add(c2, c1, root2);
			Q.assign(eco_gf2n::convertToGFElement(root2),
				 eco_gf2n::convertToGFElement(c2));
			multiply_by_2(H, Q);
			if (H == P) {
				res = Q;
				return 1;
			}
		}
	return 0;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
