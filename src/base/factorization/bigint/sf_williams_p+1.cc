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



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



const int blocking_factor_pp1 = 500;
const double scaling_pp1 = 2.0;


// Description : the following functions determine the Lucas function
// res = V_n(a) mod bigmod::modulus() with fast exponentation variant.

static void lucas(bigmod & res, const bigmod & a, unsigned long n)
{
	bigmod h(a), aa(a);
	res.assign(2);

	while (n) {
		if (n & 1) {
			multiply(res, res, aa);
			subtract(res, res, h);
		}
		else {
			multiply(h, h, aa);
			subtract(h, h, res);
		}

		square(aa, aa);
		subtract(aa, aa, 2);
		n >>= 1;
	}
}



bool single_factor< bigint >::
WilliamsPplus1_VSC(factorization< bigint > & f, const bigmod & Q,
		   unsigned int B1, unsigned int B3, ecm_primes & prim)

{
	const int giant = static_cast<int>(std::ceil(std::sqrt(static_cast<double>(B3)) / 210)) * 210;
	bigmod hQ, Q2, Q3, buf, acc;
	bigmod* tab;

	unsigned int p, v;
	int u, counter = 0;
	lidia_size_t i;

	//*********** Babysteps ********************************

	tab = new bigmod[giant];

	hQ.assign(Q);
	square(Q2, Q); // Q2 = V_2(Q)
	subtract(Q2, Q2, 2); // for i >= 3:  acc = V_(i-2)(Q)
	acc.assign(Q); //              hQ = V_i(Q)

	for (i = 1; i < giant-1; i += 2) {
		if ((i&1)*(i%3)*(i%5)*(i%7) != 0)
			tab[i].assign(hQ);
		// V_3, V_5, V_7, ...
		multiply(buf, hQ, Q2);
		subtract(buf, buf, acc);
		acc.assign(hQ);
		hQ.assign(buf);
	}

	tab[i].assign(hQ);
	lucas(hQ, Q, giant);

	//*********** Giantsteps ************

	v = (B1/giant) + 1;
	lucas(Q2, hQ, v); // Q2 = V_v(hQ), Q3 = V_(v-1)(hQ)
	lucas(Q3, hQ, v-1);
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
			multiply(buf, hQ, Q2);
			subtract(buf, buf, Q3);
			Q3.assign(Q2);
			Q2.assign(buf);
		}

		subtract(buf, Q2, tab[u]);
		multiply(acc, acc, buf);
		p = prim.getprimes();

		if (counter >= blocking_factor_pp1) {
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



void single_factor< bigint >::
WilliamsPplus1(factorization< bigint > & f, unsigned int B1,
	       unsigned int B3, const bigint & a,
	       ecm_primes & prim)
{
	bigint old_modulus(bigmod::modulus());
	bigmod::set_modulus(rep);

	bigmod acc;
	unsigned int s;
	const double log_B2 = static_cast<double>(std::log(static_cast<double>(B1)));
	const double sqrt_B2 = static_cast<int>(std::sqrt(static_cast<double>(B1)));
	long buf, exp, count = 0;

	acc.assign(a);
	lucas(acc, acc, 1073741824);
	prim.resetprimes(1);

	while ((s = prim.getprimes()) <= sqrt_B2)      // loop: exponent k > 1
	{
		exp = static_cast<int>(log_B2 / static_cast<double>(std::log(static_cast<double>(s))));
		buf = static_cast<long>(std::pow(static_cast<double>(s), static_cast<double>(exp)));
		lucas(acc, acc, buf);
		count ++;
	}

	do                                        // loop: exponent = 1
	{
		lucas(acc, acc, s);

		if (count ++ > blocking_factor_pp1) {
			count = 0;
			if (check_gcd(f, acc-2, rep)) {
				if (! old_modulus.is_zero())
					bigmod::set_modulus(old_modulus);
				return;
			}
		}
	}
	while ((s = prim.getprimes()) < B1);

	if (check_gcd(f, acc-2, rep)) {
		if (! old_modulus.is_zero())
			bigmod::set_modulus(old_modulus);
		return;
	}

	WilliamsPplus1_VSC(f, acc, B1, B3, prim);

	if (! old_modulus.is_zero())
		bigmod::set_modulus(old_modulus);
}



factorization< bigint > single_factor< bigint >:: WilliamsPplus1(int size)
{
	factorization< bigint > h;
	unsigned int B1, B3;

	if (size < 6) {
		lidia_warning_handler("single_factor< bigint >",
				      "WilliamsPplus1::size < 6");
		size = 6;
	}

	if (size > 34) {
		lidia_warning_handler("single_factor< bigint >",
				      "WilliamsPplus1::size > 34");
		size = 34;
	}

	if (rep.is_negative()) {
		h = factorization< bigint > (single_factor< bigint > (-1));
		rep.negate();
	}

	if (rep.is_one())
		return h;

	if (is_prime_factor()) {
		h.append(*this);
		rep.assign_one();
		return h;
	}

	B1 = static_cast<unsigned int>(scaling_pp1 * ecm_params[size-6][0]);
	B3 = static_cast<unsigned int>(scaling_pp1 * ecm_params[size-6][1]);

	if (info) {
		std::cout << "\nWilliams' p+1 Method with improved Standard Continuation";
		std::cout << "\nBounds:  B1 = " << B1 << ", B3 = " << B3 << "\n" << std::flush;
	}

	random_generator rg;
	long m;

	rg >> m;

	bigint s(m), g;
	g = gcd(s, rep);

	if (!g.is_one()) {
		if (info)
			std::cout << "\nfactor: " << g;

		h.append(g);
		h.append(rep/g);
		rep.assign_one();
		return h;
	}

	ecm_primes prim(1, B3, 200000);
	WilliamsPplus1(h, B1, B3, s, prim);
	return h;
}



factorization< bigint > WilliamsPplus1(const bigint & x, int size)
{
	single_factor< bigint > a(x);
	factorization< bigint > f;

	f = a.WilliamsPplus1(size);
	if (!a.is_one())
		f.append(a);
	return f;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
