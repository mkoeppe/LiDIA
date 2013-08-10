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


// 2nd part of Schoofs algorithm for finding the trace of
// Frobenius modulo l, l is assumed to be an Atkin prime,
// several possibilities for c are already determined and stored
// the correct value is returned as c[0]
//
// For convenience, rational functions are used (formulas are
// much easier)



#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/eco_gf2n.h"
#include	"LiDIA/weco2_rat_function.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//-----------------------------------------------------------------
//   returns true iff successfull
//

bool eco_gf2n::schoof_algorithm (int split_degree)
{
	ff_pol xq, xq2, yq0, yq1, yq20, yq21, hh;
	ff_pol div_pol; // the l-th division polynomial
	ff_polmod pdiv_pol; // also as poly_modulus
	int i;

	if (c.size() == 1)
		return true;

	if (info)
		std::cout << "\n\nUsing Schoof's algorithm for l = " << l << " : " << std::flush;

	compute_psi(div_pol, l); // compute the l-th division polynomial
	div_pol.make_monic(); // and it's poly_modulus
	pdiv_pol.build(div_pol);

	Xq(xq, pdiv_pol, A6.relative_degree()); // compute X^q and X^(q^2)

	/* FIX ME:     TRY TO USE SPLIT_DEGREE BETTER */
	if (split_degree != -1)
	  {
	    split_degree = compute_degree(xq, pdiv_pol, div_pol.degree());
	    EDF_one_factor (div_pol, xq, pdiv_pol, split_degree);
	    div_pol.make_monic();           // and it's poly_modulus
	    pdiv_pol.build(div_pol);
	    remainder(xq, xq, div_pol);
	  }

	compose(xq2, xq, xq, pdiv_pol);

	Ytop_f (yq0, yq1, pdiv_pol); // fast exponentiation for Y^q = yq0 + y*yq1
	compose(yq21, yq1, xq, pdiv_pol); // and Y^(q^2)
	multiply(yq20, yq0, yq21, pdiv_pol);
	multiply(yq21, yq21, yq1, pdiv_pol);
	compose(hh, yq0, xq, pdiv_pol);
	add(yq20, hh, yq20);

	weco2_rat_function::initialize(A6, pdiv_pol);
	weco2_rat_function ls, rs, h, xqyq;


	//********** compute ls == phi^2(X, Y) + q*(X, Y) first **********

	ls.assign(xq2, yq20, yq21);
	h.assign_xy();

	multiply(h, q, h);
	add(ls, ls, h);

	if (ls.is_zero()) {
		c.set_size(1);
		c[0] = 0;
		if (info)
			std::cout << "\n\n==> c = 0 mod " << l << std::flush;
		return true;
	}

	//************** now check all the possibilities for c *******

	if (info)
		std::cout << "   " << c[0] << " " << std::flush;

	xqyq.assign(xq, yq0, yq1); // =  phi(X, Y)
	multiply(rs, c[0], xqyq);

	if (ls == rs) {
		c.set_size(1);
		if (info)
			std::cout << "\n\n==> c = " << c[0] << " mod " << l << std::flush;
		return true;
	}
	else {
		negate(h, rs);
		if (ls == h) {
			c[0] = l - c[0];
			c.set_size(1);
			if (info)
				std::cout << "\n\n==> c = " << c[0] << " mod " << l << std::flush;
			return true;
		}
	}

	i = 1;

	while (i < c.size()) {
		if (info)
			std::cout << c[i] << "  " << std::flush;

		multiply(h, c[i] - c[i-1], xqyq);
		add(rs, rs, h);

		if (ls == rs) {
			c[0] = c[i];
			c.set_size(1);
			if (info)
				std::cout << "\n\n==> c = " << c[0] << " mod " << l << std::flush;
			return true;
		}
		else {
			negate(h, rs);
			if (ls == h) {
				c[0] = l - c[i];
				c.set_size(1);
				if (info)
					std::cout << "\n\n==> c = " << c[0] << " mod " << l << std::flush;
				return true;
			}
		}
		i++;
	}
	lidia_error_handler("schoof_algorithm", "no value for trace found");
	return false;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
