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
//	Author	: Emre Binisik (EB), Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/single_factor.h"
#include	"LiDIA/factorization.h"
#include	"LiDIA/base/ecm_primes.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



factorization< bigint > single_factor< bigint >::
factor(int size) const
{
	factorization< bigint > f;
	single_factor< bigint > a(*this);

	if (size < 6) {
		lidia_warning_handler("single_factor< bigint >",
				      "factor()::size < 6");
		size = 6;
	}

	if (size > 34) {
		lidia_warning_handler("single_factor< bigint >",
				      "factor()::size > 34");
		size = 34;
	}

	if (size > static_cast<int>(decimal_length(rep))/2)
		size = decimal_length(rep)/2 + 1;

	if (a.rep.is_negative()) {
		f = factorization< bigint > (single_factor< bigint > (-1));
		a.rep.negate();
	}

	if (a.is_prime_factor(1)) {
		f.append(*this);
		return f;
	}
	int D;
	bigint R;

	if ((D = is_power(R, a.rep)) > 1) {
		if (info)
			std::cout << "\nPower Test: factor " << R << " ^ " << D << "\n" << std::flush;
		f.append(single_factor< bigint > (R), D);
		a.rep.assign_one();
		return f;
	}

	D = single_factor::ecm_read_max(size);

	if (D < 1000000)
		D = 1000000;

	ecm_primes prim(1, D+200, 200000);

	if (info)
		std::cout << "\nTrial Division from 2 to 1000000\n" << std::flush;

	a.TrialDiv(f, 2, 1000000, prim);

	if (a.rep.is_one())
		return f;

	size = decimal_length(a.rep) / 2;
	if (size > 10)  size = 10;
	f.concat(f, a.PollardPminus1(size));
	if (a.rep.is_one()) {
		f.refine();
		return f;
	}

	size = decimal_length(a.rep) / 2;
	if (size > 10)  size = 10;
	f.concat(f, a.WilliamsPplus1(size));
	if (a.rep.is_one()) {
		f.refine();
		return f;
	}

	size = decimal_length(a.rep) / 2;

	f.concat(f, a.ECM(size, 6, 3, true));
	if (!a.rep.is_one()) {
		f.append(*this);
		f.refine();
	}

	return f;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
