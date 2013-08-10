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
//	Author	: Volker Mueller
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/nmbrthry_functions.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// Let n be the integer represented by f.
//
//
// flag = 0; h_divisors computes all positive divisors of n.
// flag = 1; h_divisors computes all positive square-free divisors of n
// flag = 2; h_divisors computes all positive divisors whose square divides n
//

void h_divisors (sort_vector< bigint > & div, rational_factorization & f, int flag)
{
	lidia_size_t e, i, j, k, l, fno, fexp;
	lidia_size_t cap = 1;
	bigint h;
	bigint base;

	if (f.is_prime_factorization() == false) {
		lidia_warning_handler("rational_factorization", "divisors::input no prime factorization");
		do {
			f.factor();
		} while (f.is_prime_factorization() == false);
	}

	fno = f.no_of_comp();

	switch (flag) {
	case 0:
		for (i = 0; i < fno; i++)
			cap *= (1 + f.exponent(i));
		break;

	case 1:
		cap = (1 << fno);
		break;

	case 2:
		for (i = 0; i < fno; i++)
			cap *= (1 + (f.exponent(i) >> 1));
		break;

	default:
		lidia_error_handler ("void h_divisors (sort_vector< bigint > &, rational_factorization &, int)",
				     "Invalid value for parameter int flag");
		break;
	}

	div.set_capacity(cap);
	div.set_size(cap);

	k = 0;
	l = 1;
	div[0] = bigint (1);

	for (i = 0; i < fno; i++) {
		h.assign_one ();
		fexp = f.exponent(i);
		base = f.base(i);

		switch (flag) {
		case 0:  break;
		case 1:  fexp = 1;
			break;
		case 2:  fexp >>= 1;
			break;
		default: lidia_error_handler ("void h_divisors (sort_vector< bigint > &, rational_factorization &, int)",
					      "Invalid value for parameter int flag");
			break;
		}

		for (e = 1; e <= fexp; e++) {
			multiply (h, h, base);

			for (j = 0; j <= k; j++)
				multiply (div[l++], div[j], h);
		}
		k = l - 1;
	}
}



//
// adds the negative divisors
//

void h_all_divisors (sort_vector< bigint > & div, rational_factorization & f, int flag)
{
	lidia_size_t i, div_sz, div_cap;

	h_divisors(div, f, flag);
	div.sort(SORT_VECTOR_DOWN);

	div_sz = div.size();
	div_cap = (static_cast<lidia_size_t>(2))*div_sz;
	div.set_capacity(div_cap);
	div.set_size(div_cap);

	div_cap--;
	for (i = 0; i < div_sz; i++) {
		div[div_cap-i] = div[i];
		div[i].negate();
	}
}



//
// interface functions
//


//
// positive divisors
//

sort_vector< bigint > divisors (rational_factorization & f)
{
	sort_vector< bigint > div;
	h_divisors(div, f);
	div.sort();
	return div;
}



sort_vector< bigint > divisors (const bigint & N)
{
	rational_factorization f(N);
	f.factor();
	return divisors(f);
}



//
// positive and negative divisors
//

sort_vector< bigint > all_divisors (rational_factorization & f)
{
	sort_vector< bigint > div;
	h_all_divisors(div, f);
	return div;
}



sort_vector< bigint > all_divisors (const bigint & N)
{
	rational_factorization f(N);
	f.factor();
	return all_divisors(f);
}



//
// positive square-free divisors
//

sort_vector< bigint > square_free_divisors (rational_factorization & f)
{
	sort_vector< bigint > div;
	h_divisors(div, f, 1);
	div.sort();
	return div;
}



sort_vector< bigint > square_free_divisors (const bigint & N)
{
	rational_factorization f(N);
	f.factor();
	return square_free_divisors(f);
}



//
// positve and negative square-free divisors
//

sort_vector< bigint > all_square_free_divisors (rational_factorization & f)
{
	sort_vector< bigint > div;
	h_all_divisors(div, f, 1);
	div.sort();
	return div;
}



sort_vector< bigint > all_square_free_divisors (const bigint & N)
{
	rational_factorization f(N);
	f.factor();
	return all_square_free_divisors(f);
}



//
// positive divisors whose square divides n
//

sort_vector< bigint > square_divides_n_divisors (rational_factorization & f)
{
	sort_vector< bigint > div;
	h_divisors(div, f, 2);
	div.sort();
	return div;
}



sort_vector< bigint > square_divides_n_divisors (const bigint & N)
{
	rational_factorization f(N);
	f.factor();
	return square_divides_n_divisors(f);
}



//
// positive and negative divisors whose square divides n
//

sort_vector< bigint > all_square_divides_n_divisors (rational_factorization & f)
{
	sort_vector< bigint > div;
	h_all_divisors(div, f, 2);
	div.sort();
	return div;
}



sort_vector< bigint > all_square_divides_n_divisors (const bigint & N)
{
	rational_factorization f(N);
	f.factor();
	return all_square_divides_n_divisors(f);
}



//***********************************************************************
// elementary number theory functions needed for example in
// elliptic curve class.
//***********************************************************************

// return true if f is a divisor of n

bool is_divisor (const bigint & f, const bigint & n)
{
	if (f.is_zero())
		return n.is_zero();
	else
		return (n %f).is_zero();
}



// return maximal exponent e such that f^e | n

int valuation(const bigint & ff, const bigint & nn)
{
	bigint n(abs(nn));
	bigint f(abs(ff));
	bigint h;

	if (n.is_zero() || (f < 2)) {
		lidia_error_handler("file divisors.c::"
				    "int valuation(const bigint&, const bigint&)",
				    "abs(n) == 0 or abs(f) < 2");
		return 0;
	}

	int e = 0;

	div_rem(n, h, n, f);
	while (h.is_zero()) {
		e++;
		div_rem(n, h, n, f);
	}
	return e;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
