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
# include "config.h"
#endif
#include "LiDIA/single_factor.h"
#include "LiDIA/factorization.h"
#include "LiDIA/base/ecm_primes.h"
#include "LiDIA/precondition_error.h"


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



void single_factor< bigint >::
TrialDiv(factorization< bigint > & f,
	 const unsigned int lower_bound,
	 const unsigned int upper_bound,
	 ecm_primes & prim)
{
	bigint N(rep), TMP;
	single_factor< bigint > d;

	long r;
	unsigned long p;
	int exp = 0;

//--------------------------------------------------------------------------
// Division by 2
//--------------------------------------------------------------------------

	if (lower_bound <= 2)              //The division by 2 is faster then the
	{                                 //division by other numbers,
		while (N.is_even())            //so we do it extra.
		{
			exp++;
			N.divide_by_2();
		}

		if (exp) {
			d.rep = 2;
			d.set_prime_flag(prime);
			f.append(d, exp);
			if (single_factor< bigint >::info)
				if (exp > 1)
					std::cout << "factor: 2 ^ " << exp << "\n";
				else  std::cout << "factor: 2 \n";
			exp = 0;
		}
	}


//--------------------------------------------------------------------------
// Division by all prime numbers between lower_bound and upper_bound
//--------------------------------------------------------------------------

	prim.resetprimes(lower_bound);
	p = prim.getprimes();

	while ((!N.is_one()) && (p < upper_bound) && (p != 1)) {
		div_rem(TMP, r, N, p);
		while (!r)              //As long as the remainder r is zero
		{
			N.assign(TMP);
			div_rem(TMP, r, N, p); //we divide N by the prime number p
			exp++; //and count the divisions in exponent.
		}

		if (exp)              //That means p to the power of exp divides N
		{                    //with an exponent exp > 0.
			d.rep = p; //We have a single_factor< bigint > [p, exp]
			d.set_prime_flag(prime); //and collect it in the Variable f
			f.append(d, exp); //of the type class factorization< bignt >

			if (single_factor< bigint >::info)
				if (exp > 1)
					std::cout << "factor: " << p << " ^ " << exp << "\n";
				else  std::cout << "factor: " << p << " \n";

			exp = 0;
		}

		p = prim.getprimes();
	}

	if (rep != N) {
		rep.assign(N);

		if (is_prime(N, 8)) {
			set_prime_flag(prime);
			f.append(*this);
			rep.assign_one();
		}
		else
			set_prime_flag(not_prime);
	}
}



factorization< bigint > single_factor< bigint >::
TrialDiv(unsigned int upper_bound, unsigned int lower_bound)
{
	factorization< bigint > f;
	bigint N((*this).rep);

	// Parametercheck

	if (lower_bound > upper_bound) {
		precondition_error_handler(upper_bound, "upper_bound",
					   "upper_bound >= lower_bound",
					   lower_bound, "lower_bound",
					   "upper_bound >= lower_bound",
					   "factorization<bigint> "
					   "single_factor< bigint >::"
					   "TrialDiv(unsigned int upper_bound,"
					   " unsigned int lower_bound)",
					   "single_factor< bigint >",
					   "lower_bound > upper_bound");
		return f; // we should never get here!
	}
	if (lower_bound < 1) {
		lidia_warning_handler("single_factor< bigint >",
				      "TrialDiv::lower_bound < 2");
		lower_bound = 2;
	}

	if (N.is_zero()) {
		precondition_error_handler(N, "this->rep", "this-rep != 0",
					   "factorization<bigint> "
					   "single_factor< bigint >::"
					   "TrialDiv(unsigned int upper_bound,"
					   " unsigned int lower_bound)",
					   "single_factor< bigint >",
					   "TrialDiv::input 0 is not allowed");
		return f;
	}

	if (N.is_negative()) {
		f = factorization< bigint > (single_factor< bigint > (-1));
		N.negate();
	}

	if (N.is_one())
		return f;

	if ((*this).is_prime_factor(1)) {
		f.append(*this);
		(*this).rep.assign_one();
		return f;
	}


	single_factor< bigint > a(N);

	ecm_primes prim(lower_bound, upper_bound+200, 200000);

	if (info) {
		std::cout << "\nTrial Division from ";
		std::cout << lower_bound << " to " << upper_bound;
		std::cout << "\n";
	}

	a.TrialDiv(f, lower_bound, upper_bound, prim);
	*this = a;

	return f;
}



factorization< bigint >
TrialDiv(const bigint & x,
	 const unsigned int upper_bound,
	 const unsigned int lower_bound)
{
	factorization< bigint > f;
	single_factor< bigint > a(x);

	f = a.TrialDiv(upper_bound, lower_bound);
	if (!a.rep.is_one())
		f.append(a);
	return(f);
}

factorization< bigint >
TrialDiv(const bigint & x,
	 const unsigned int upper_bound)
{
	factorization< bigint > f;
	single_factor< bigint > a(x);

	f = a.TrialDiv(upper_bound);
	if (!a.rep.is_one())
		f.append(a);
	return(f);
}

factorization< bigint >
TrialDiv(const bigint & x)
{
	factorization< bigint > f;
	single_factor< bigint > a(x);

	f = a.TrialDiv();
	if (!a.rep.is_one())
		f.append(a);
	return(f);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
