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
# include "config.h"
#endif
#include "LiDIA/bigint.h"
#include "LiDIA/single_factor.h"
#include "LiDIA/factorization.h"
#include "LiDIA/precondition_error.h"


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



void single_factor< bigint >:: Fermat(factorization< bigint > & f,
				       unsigned int iterations)
{
	bigint N(rep);
	bigint m, z, df, root;
	register unsigned int i;

	//*********** Initialization **************************

	sqrt(m, N);
	inc(m);

	if ((N.least_significant_digit() & 3) == 3 && m.is_odd())  // N = 3 mod 4
		inc(m);
	if ((N.least_significant_digit() & 3) == 1 && m.is_even())  // N = 1 mod 4
		inc(m);

	// z = m^2-N

	square(z, m);
	subtract(z, z, N);

	// df = 4*m+4

	multiply(df, m, 4);
	add(df, df, 4);

	for (i = 1; i <= iterations; i++) {
		if (is_square(root, z))
			if (root.is_positive()) {
				subtract(m, df, 4);
				shift_right(m, m, 2);

				if (info)
					std::cout << "\nfactor: " << m-root << std::flush;

				single_factor< bigint > a(m-root);
				if (is_prime(m-root, 8))
					a.set_prime_flag(prime);
				else
					a.set_prime_flag(not_prime);
				f.append(a);

				a = m+root;
				if (is_prime(m+root, 8))
					a.set_prime_flag(prime);
				else
					a.set_prime_flag(not_prime);
				f.append(a);
				rep.assign_one();
				return;
			}
		add(z, z, df);
		add(df, df, 8);
	}
}



factorization< bigint > single_factor< bigint >::
Fermat()
{
	factorization< bigint > f, h;
	single_factor< bigint > a;
	bigint N((*this).rep);

	// Parameter check

	if (N.is_zero()) {
		precondition_error_handler(this->rep, "this->rep",
					   "this->rep != 0",
					   "factorization<bigint> "
					   "single_factor<bigint>::Fermat()",
					   "single_factor< bigint >",
					   "Fermat::input 0 is not allowed");
		return h;
	}

	if (N.is_negative()) {
		h = factorization< bigint > (single_factor< bigint > (-1));
		N.negate();
	}

	if (N.is_one()) {
		rep.assign_one();
		return h;
	}

	if (is_prime_factor()) {
		h.append(*this);
		rep.assign_one();
		return h;
	}

	int exp = 0;

	while (N.is_even()) {
		exp++;
		N.divide_by_2();
	}
	if (exp) {
		a.rep.assign(2);
		a.set_prime_flag(prime);
		h.append(a, exp);
	}

	if (N.is_one()) {
		rep.assign_one();
		return h;
	}

	if (is_prime(N, 8)) {
		a.rep.assign(N);
		a.set_prime_flag(prime);
		h.append(a);
		rep.assign_one();
		return h;
	}

	bigint root;
	exp = is_power(root, N);

	if (exp > 1) {
		a.rep.assign(root);
		if (is_prime(root, 8))
			a.set_prime_flag(prime);
		else
			a.set_prime_flag(not_prime);
		h.append(a, exp);

		if (info)
			std::cout << "factor: " << root << " ^ " << exp;

		(*this).rep.assign_one();
		return h;
	}

	if (info) {
		std::cout << "\nFermat's Method (#iterations = 1000000): \n";
		std::cout.flush();
	}

	a.rep.assign(N);
	a.set_prime_flag(not_prime);

	a.Fermat(f, 1000000);

	if (exp > 1)
		power(f, f, exp);
	multiply(f, f, h);

	*this = a;
	return f;
}



factorization< bigint > Fermat(const bigint & x)
{
	single_factor< bigint > a(x);
	factorization< bigint > f;
	f = a.Fermat();
	if (!a.is_one())
		f.append(a);
	return f;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
