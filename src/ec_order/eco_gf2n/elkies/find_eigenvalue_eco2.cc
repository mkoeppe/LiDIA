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
//	Author	: Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================



#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/eco_gf2n.h"
#include	"LiDIA/weco2_rat_function.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



// these functions search for the eigenvalue of Frobenius after
// Elkies polynomial is known. Several search algorithms are supported.
// returns true <=> something has been found

//-----------------------------------------------------------------
//  alpha    : eigenvalue to be computed
//  fC       : Elkies polynomial
//  lpower   : compute eigenvalue mod l^lpower
//  prev_ev  : eigenvalue mod l^(lpower-1) (== -1 <==> not yet known)
//
//  test_all : controls the search for alpha in
//             Phi (X, Y) == alpha * (X, Y)
//
//             0 --> stop searching after first match
//             1 --> test all possible values for alpha
//                   or stop after second match
//
//

bool eco_gf2n::find_eigenvalue (lidia_size_t & alpha,
				const ff_pol & fC,
				lidia_size_t prev_ev,
				lidia_size_t lpower,
				bool test_all)
{
	lidia_size_t *klist; // list of possible eigenvalues
	lidia_size_t lpow = l;
	bool  rc =  false;

	//****** computation for l == 3 ******

	if (l == 3 && lpower == 1)   // deg(fC) == (l-1)/2 == 1 --> computation in Fq
		return eigenv_mod3 (alpha, fC);


	//****** computation for l > 3 or lpower > 1 ******

	if (sp_type == ONE_SUBGRP_INV || sp_type == ALL_SUBGRPS_INV) {
		if (lpower == 1) {
			// eigenvalue^2 = p^n = q mod l
			klist = new lidia_size_t [2];
			klist[0] = 1;

			if (jacobi (static_cast<udigit>(eco_gf2n::q), static_cast<udigit>(this->l))) {
				klist[1] = static_cast<lidia_size_t>(sqrt_mod(static_cast<udigit>(eco_gf2n::q),
										     static_cast<udigit>(this->l)));
				if (klist[1] > static_cast<int>(this->l-1)/2)
					klist[1] = this->l - klist[1];

				if (test_x_coordinate()) {
					alpha = klist[1];
					if (test_x_coordinate_uniquely())
						compute_sign(alpha, fC);
					rc = true;
				}
				else {
					char ev = ev_strategy & 0x0f;

					if (ev == eco_gf2n::EV_RATIONAL_FUNCTION ||
					    ev == eco_gf2n::EV_RATIONAL_FUNCTION_TABLE)
						rc = schoofpart_rat_function(alpha, fC, klist, test_all);
					else
						rc = schoofpart(alpha, fC, klist, test_all);

					if (alpha < 0)
						alpha += l;
				}
				delete[] klist;
				return rc;
			}
			else {
				lidia_error_handler ("eco_gf2n::find_eigenvalue",
						     "p^n is not a square mod l");
				return false;
			}
		}
		else {
			lidia_error_handler ("eco_gf2n::find_eigenvalue()",
					     "Sorry, eigenvalue modulo powers of l "
					     "for splitting type (1...1) or (1 l) not yet "
					     "implemented");
			return  false;
		}
	}

	//------------------------------------------------------------
	// now we have Elkies type (11ddd)
	// first we determine list of candidates for ev.


	if (lpower == 1) {
		lidia_size_t * d_list;

		if (sp_degree_known) {
			// compute list of possible d.
			d_list = new lidia_size_t[2];
			d_list[0] = 1;
			d_list[1] = sp_degree;
		}
		else
			compute_d_list (d_list);

		// use this list for possible ev's
		compute_ev_list (d_list, klist);
		delete[] d_list;
	}
	else {
		lpow = static_cast<lidia_size_t>(std::pow(static_cast<double>(l), static_cast<double>(lpower - 1)));

		if (prev_ev == -1) {
			lpow *= l;
			klist = new lidia_size_t [lpow + 1];
			klist[0] = lpow;
			for (int i = 1; i <= lpow; i++)
				klist[i] = i;
		}
		else {
			klist = new lidia_size_t [l + 1];
			klist[0] = l;
			for (unsigned int i = 1; i <= l; i++)
				klist[i] = prev_ev + (i-1) * lpow;
			lpow *= l;
		}
	}

	//---------------------------------------------------------
	// now we search for the ev.

	char ev = ev_strategy & 0x0f;

#ifdef TIMING
	timer t;
	t.set_print_mode();
	t.start_timer();
#endif

	if (ev == eco_gf2n::EV_RATIONAL_FUNCTION ||
	    ev == eco_gf2n::EV_RATIONAL_FUNCTION_TABLE)
		rc = schoofpart_rat_function(alpha, fC, klist, test_all);
	else if (ev == eco_gf2n::EV_DIVISION_POLYNOMIAL)
		rc = schoofpart(alpha, fC, klist, test_all);
	else if (ev == eco_gf2n::EV_BABYSTEP_GIANTSTEP)
		rc = schoofpart_BG(alpha, fC, klist);
	else if (ev == eco_gf2n::EV_FUNNY_BABYSTEP_GIANTSTEP) {
		lidia_size_t old_l = l;
		l = lpow;
		rc = schoofpart_FBG(alpha, fC, klist, false);
		l = old_l;
	}
	else if (ev == eco_gf2n::EV_FUNNY_BABYSTEP_GIANTSTEP_SIGN) {
		lidia_size_t old_l = l;
		l = lpow;
		rc = schoofpart_FBG(alpha, fC, klist, true);
		l = old_l;
	}
	else
		lidia_error_handler("eco_prime",
				    "wrong ev strategy chosen");

#ifdef TIMING
	t.stop_timer();
	if (info)
		std::cout << "\nTotal time needed to find eigenvalue : " << t << std::flush;
#endif

	delete[] klist;

	if (alpha < 0)
		alpha += lpow;

	return rc;
}



//**********************************************************************
//
// this is the special version for finding the eigenvalue of Frobenius
// for l==3. In this case, all computations can be done in GF(2^n), since
// degree of Elkies polynomial is 1.
// fC = (monic) Elkies polynomial for l == 3
//

bool eco_gf2n::eigenv_mod3 (lidia_size_t & alpha, const ff_pol &fC)
{
	ff_element r1, r0, h, x, h1;
	int d = A6.relative_degree(), i;

	r1.assign(fC.const_term()); // compute Y^2
	x.assign(r1);

	square(r0, r1);
	multiply(r0, r0, r1);
	add(r0, r0, A6);
	h.assign(r0);

	for (i = 2; i <= d; i++) {
		// to compute Y^q
		square(h1, r1);
		multiply(r1, h1, x);

		multiply(h1, h1, h);
		square(r0, r0);
		add(r0, r0, h1);
	}

	if (r0.is_zero() && r1.is_one())
		alpha = 1;
	else
		alpha = 2;

	return (true);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
