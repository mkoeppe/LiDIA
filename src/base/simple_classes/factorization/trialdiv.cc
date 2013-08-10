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
//	Author	: Andreas Mueller (AM), Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/rational_factorization.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



void rational_factorization::
trialdiv(lidia_size_t index, unsigned int lower_bound,
	 unsigned int upper_bound, ecm_primes &prim)
{
	bigint N(factors[index].single_base);
	int lexp = factors[index].single_exponent;

	rf_single_factor fact;

	lidia_size_t n = no_of_comp();
	long p, q;
	int count = 0;

	bigint TMP;

	if (info) {
		std::cout << "index " << index << ":   trialdivision up to " << upper_bound << "\n";
		std::cout.flush();
	}

	prim.resetprimes(lower_bound);

	if (lower_bound < 2) {
		while (N.is_even()) {
			count++;
			N.divide_by_2();
		}

		if (count) {
			fact.single_base = bigint(2);
			fact.single_exponent = count * lexp;
			fact.factor_state = rf_single_factor::prime;
			factors[n++] = fact;

			if (info)
				if (count > 1)
					std::cout << "factor: 2 ^ " << count << "\n";
				else
					std::cout << "factor: 2 \n";

			count = 0;
		}
	}

	while ((static_cast<unsigned int>(p = prim.getprimes()) < upper_bound) &&
	       !N.is_one() && p != 1) {
		div_rem(TMP, q, N, p);

		while (!q) {
			N.assign(TMP);
			count++;
			div_rem(TMP, q, N, p);
		}

		if (count) {
			fact.single_base.assign(p);
			fact.single_exponent = count * lexp;
			fact.factor_state = rf_single_factor::prime;

			factors[n++] = fact;

			if (info)
				if (count > 1)
					std::cout << "factor: " << p << " ^ " << count << "\n";
				else
					std::cout << "factor: " << p << " \n";

			count = 0;
		}
	}

	if (N.is_one()) {
		n--;
		factors[index] = factors[n];
		factors.remove_from(n);
	}
	else {
		fact.single_base.assign(N);
		fact.single_exponent = lexp;
		fact.factor_state = rf_single_factor::dont_know;

		factors[index] = fact;
	}
}



void  rational_factorization::
trialdiv_comp(lidia_size_t index, unsigned int upper_bound,
	      unsigned int lower_bound)
{
	if ((index< 0) || (index >= no_of_comp())) {
		lidia_error_handler("rational_factorization", "trialdiv_comp::index out of range");
		return;
	}

	if (rf_single_factor::state(factors[index].factor_state) == rf_single_factor::prime) {
		if (info)
			std::cout << "prime number " << factors[index].single_base << "\n";
		return;
	}

	if (lower_bound > upper_bound) {
		lidia_warning_handler("rational_factorization", "trialdiv_comp::lower_bound > upper_bound");
		return;
	}

	ecm_primes prim(lower_bound, upper_bound+200, 200000);

	trialdiv(index, lower_bound, upper_bound, prim);

	compose();
}



void  rational_factorization::
trialdiv(unsigned int upper_bound, unsigned int lower_bound)
{
	if (lower_bound > upper_bound) {
		lidia_warning_handler("rational_factorization", "trialdiv::lower_bound > upper_bound");
		return;
	}

	ecm_primes prim(lower_bound, upper_bound+200, 200000);

	lidia_size_t index, len = no_of_comp();

	for (index = 0; index < len; index++) {
		if (rf_single_factor::state(factors[index].factor_state) == rf_single_factor::prime) {
			if (info)
				std::cout << "index " << index << " prime number " << factors[index].single_base << "\n";
		}
		else
			trialdiv(index, lower_bound, upper_bound, prim);
	}

	compose();
}



rational_factorization trialdiv(const bigint &N, unsigned int upper_bound,
				unsigned int lower_bound)
{
	rational_factorization f(N);

	f.trialdiv(upper_bound, lower_bound);

	return f;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
