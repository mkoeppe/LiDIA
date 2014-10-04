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
#include	"LiDIA/random_generator.h"
#include	"LiDIA/power_functions.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



const int blocking_factor_pm1 = 100;
const double scaling_pm1 = 2.0;


bool single_factor< bigint >::
PollardPminus1_VSC(factorization< bigint > & f,
		   const bigmod & Q,
		   unsigned int B1,
		   unsigned int B3,
		   ecm_primes & prim)
{
	const int giant = static_cast<int>(std::ceil(std::sqrt(static_cast<double>(B3)) / 210)) * 210;

	bigmod hQ, Q2, buf, acc;
	bigmod* tab;

	unsigned int p, v;
	int u, counter = 0;
	lidia_size_t i;

	//*********** Babysteps ********************************

	tab = new bigmod[giant];

	hQ.assign(Q);
	square(Q2, Q);

	for (i = 1; i < giant-1; i += 2) {
		if ((i&1)*(i%3)*(i%5)*(i%7) != 0)
			tab[i].assign(hQ);
		// Q^3, Q^5, Q^7, ... , Q^(giant-1)
		multiply(hQ, hQ, Q2);
	}

	tab[i].assign(hQ);
	multiply(hQ, hQ, Q); // hQ = Q^(giant)


	//*********** Giantsteps ************

	v = (B1/giant) + 1;
	power(Q2, hQ, v);
	v *= giant;

	prim.resetprimes(B1);
	p = prim.getprimes();
	acc.assign_one();

	while (p< B3 && p > 1) {
		u = v - p;
		counter ++;

		if (u < 0) {
			v += giant;
			u = v - p;
			multiply(Q2, Q2, hQ);
		}

		subtract(buf, Q2, tab[u]);
		multiply(acc, acc, buf);
		p = prim.getprimes();

		if (counter >= blocking_factor_pm1) {
			counter = 0;
			if (check_gcd(f, acc, rep)) {
				delete[] tab;
				return true;
			}
		}
	}

	delete [] tab;

	if (check_gcd(f, acc, rep))
		return true;

	return false;
}



void single_factor< bigint >:: PollardPminus1(factorization< bigint > & f,
					       unsigned int B1,
					       unsigned int B3,
					       const bigint & a,
					       ecm_primes& prim)
{
	bigint old_modulus(bigmod::modulus());
	bigmod::set_modulus(rep);

	bigmod acc;
	unsigned int s;
	const double log_B2 = static_cast<double>(std::log(static_cast<double>(B1)));
	const unsigned int sqrt_B2 = static_cast<unsigned int>(std::sqrt(static_cast<double>(B1)));
	unsigned long buf;
	lidia_size_t i;
	int exp, count = 0;

	acc.assign(a);

	for (i = 0; i < 32; i++)
		square(acc, acc);

	prim.resetprimes(1);

	while ((s = prim.getprimes()) <= sqrt_B2)      // loop: exponent k > 1
	{
		exp = static_cast<int>(log_B2 / static_cast<double>(std::log(static_cast<double>(s))));
		buf = static_cast<unsigned long>(std::pow(static_cast<double>(s), static_cast<double>(exp)));
		lidia_power_left_to_right(acc, acc, buf);
		count++;
	}

	do                                        // loop: exponent = 1
	{
		lidia_power_left_to_right(acc, acc, static_cast<unsigned long>(s));

		if (count++ > blocking_factor_pm1) {
			count = 0;
			if (check_gcd(f, acc-1, rep)) {
				if (!old_modulus.is_zero())
					bigmod::set_modulus(old_modulus);
				return;
			}
		}
	}
	while ((s = prim.getprimes()) < B1);

	if (check_gcd(f, acc-1, rep)) {
		if (!old_modulus.is_zero())
			bigmod::set_modulus(old_modulus);
		return;
	}


	PollardPminus1_VSC(f, acc, B1, B3, prim);

	if (! old_modulus.is_zero())
		bigmod::set_modulus(old_modulus);
}



factorization< bigint > single_factor< bigint >:: PollardPminus1(int size)
{
	factorization< bigint > h;
	unsigned int B1, B3;
	bigint g;

	if (size < 6) {
		lidia_warning_handler("single_factor< bigint >",
				      "PollardPminus1::size < 6");
		size = 6;
	}

	if (size > 34) {
		lidia_warning_handler("single_factor< bigint >",
				      "PollardPminus1::size > 34");
		size = 34;
	}

	if (rep.is_negative()) {
		h.append(single_factor< bigint > (-1));
		rep.negate();
	}

	if (is_prime_factor()) {
		h.append(*this);
		return h;
	}


	B1 = static_cast<unsigned int>(scaling_pm1 * ecm_params[size-6][0]);
	B3 = static_cast<unsigned int>(scaling_pm1 * ecm_params[size-6][1]);

	if (info) {
		std::cout << "\nPollard p-1 Algorithm with improved Standard Continuation";
		std::cout << "\nBounds:  B1 = " << B1 << ", B3 = " << B3 << "\n" << std::flush;
	}

	ecm_primes prim(1, B3, 200000);

	random_generator rg;
	long m;

	rg >> m;
	bigint s(m);
	g = gcd(s, rep);

	if (!g.is_one()) {
		if (info)
			std::cout << "\nfactor: " << g;

		h.append(g);
		h.append(rep/g);
		rep.assign_one();
		return h;
	}

	PollardPminus1(h, B1, B3, s, prim);
	return h;
}



factorization< bigint >
PollardPminus1(const bigint & x, int size)
{
	factorization< bigint > f;
	single_factor< bigint > a(x);

	f = a.PollardPminus1(size);
	if (!a.is_one())
		f.append(a);
	return f;
}

factorization< bigint >
PollardPminus1(const bigint & x)
{
	factorization< bigint > f;
	single_factor< bigint > a(x);

	f = a.PollardPminus1();
	if (!a.is_one())
		f.append(a);
	return f;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
