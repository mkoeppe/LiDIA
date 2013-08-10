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
#include	"LiDIA/bigmod.h"
#include	"LiDIA/single_factor.h"
#include	"LiDIA/factorization.h"
#include	"LiDIA/random_generator.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



const int blocking_factor_rho = 100;


void single_factor< bigint >:: PollardRho(factorization< bigint > & f,
					   unsigned int B)
{
	bigint old_modulus(bigmod::modulus());
	bigmod::set_modulus(rep);

	bigmod xs, xc, acc, diff;
	unsigned int range, i;
	random_generator rg;
	int a, counter = 0;

	rg >> a;
	xs.assign(a);
	xc.assign(a);

	range = 10;
	acc.assign_one();

	while (range < B) {
		for (i = 1; i <= range; i++) {
			square(xc, xc); //The algorithm uses the function x^2+1
			inc(xc); //to create the sequence of random integers
		}

		for (i = 1; i <= range; i++) {
			square(xc, xc);
			inc(xc);
			subtract(diff, xc, xs);
			multiply(acc, acc, diff);

			if (counter++ > blocking_factor_rho) {
				counter = 0;
				if (check_gcd(f, acc, rep)) {
					if (! old_modulus.is_zero())
						bigmod::set_modulus(old_modulus);
					return;
				}
			}
		}

		xs.assign(xc);
		range <<= 1;
	}

	check_gcd(f, acc, rep);

	if (! old_modulus.is_zero())
		bigmod::set_modulus(old_modulus);
}



factorization< bigint > single_factor< bigint >:: PollardRho(int size)
{
	factorization< bigint > h;

	if (size < 1) {
		lidia_warning_handler("single_factor< bigint >",
				      "PollardRho::size < 1");
		size = 6;
	}

	if (size > 34) {
		lidia_warning_handler("single_factor< bigint >",
				      "PollardRho::size > 34");
		size = 34;
	}

	if (rep.is_negative()) {
		h.append(single_factor< bigint > (-1));
		rep.negate();
	}

	if (rep.is_one())
		return h;

	if (is_prime_factor(1)) {
		h.append(*this);
		rep.assign_one();
		return h;
	}

	unsigned int B;

	B = static_cast<unsigned int>(std::ceil(1.0308 * std::sqrt(std::pow(10.0, static_cast<double>(size)) / 2.0)));

	if (info)
		std::cout << "\nPollard Rho Algorithm (#iterations = " << B << ")\n";

	PollardRho(h, B);
	return h;
}



factorization< bigint > PollardRho (const bigint & x, int size)
{
	factorization< bigint > f;
	single_factor< bigint > a(x);
	f = a.PollardRho(size);
	f.append(a);
	return f;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
