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


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigint.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



int jacobi(const bigint &a1, const bigint &b1)
{
	bigint a, b;
	long k = 1, v;

	// Setting up the table for table-lookup

	int table[8] = {0, 1, 0, -1, 0, -1, 0, 1};

	// test trivial cases

	if (b1.is_zero()) {
		a = abs(a1);
		if (a.is_one()) return(1);
		else return(0);
	}

	if (a1.is_even() && b1.is_even())
		return (0);

	a.assign(a1);
	b.assign(b1);

	if (b.is_negative()) {
		b.negate();
		if (a.is_negative())
			k = -1;
	}

	v = 0;
	while (b.is_even()) {
		v++;
		b.divide_by_2();
	}

	if (v&1)
		k *= table[a.least_significant_digit()&7];

	if (a.is_negative()) {
		if (b.least_significant_digit()&2)
			k = -k;
		a.negate();
	}

	while (!a.is_zero()) {
		// main-loop
		v = 0;
		while (a.is_even()) {
			v++;
			a.divide_by_2();
		}
		if (v&1)
			k *= table[b.least_significant_digit()&7];

		if (a.compare(b) < 0) {
			swap(a, b); // swap and correct intermediate result
			if (a.least_significant_digit() & b.least_significant_digit() & 2)
				k = -k;
		}
		subtract(a, a, b);
	}

	if (b.is_one()) return(k);
	else return(0);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
