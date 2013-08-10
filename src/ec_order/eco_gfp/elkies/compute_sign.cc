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


//  Uses the Dewaghe idea top compute the sign of the eigenvalue



#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/eco_prime.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



void eco_prime::sign_dewaghe (lidia_size_t & alpha,
			      const ff_pol & fC)
{
	if (l % 4 == 1) {
		lidia_error_handler("eco_prime", "sign_dewaghe::Dewaghe's idea only "
				    "usable for l = 3 mod 4 !!");
		return;
	}

	if (alpha < 0)
		alpha = l - alpha;

	Fp_polynomial curve;
	bigmod h;

	curve.set_modulus(A.modulus());
	curve.set_coefficient(3);
	curve.set_coefficient(A.mantissa(), 1);
	curve.set_coefficient(B.mantissa(), 0);

	h = resultant(fC, curve);

	if (jacobi(h.mantissa(), pn) == 1) {
		if (jacobi(static_cast<udigit>(alpha), static_cast<udigit>(l)) == -1)
			alpha = -alpha;
	}
	else
		if (jacobi(static_cast<udigit>(alpha), static_cast<udigit>(l)) == 1)
			alpha = -alpha;

	if (alpha < 0)
		alpha = l + alpha;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
