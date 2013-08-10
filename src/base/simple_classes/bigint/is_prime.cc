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
//	Author	: Andreas M"uller (AM)
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



bool
bigint::is_prime (const int b1) const
{
	static long a[10] = {3, 5, 7, 11, 13, 17, 19, 23, 29, 31};
	long b, j, ok, i, k = 0, sx;

	if (!longify(sx)) {
		// can not be converted to long
		if (sx <= 0)
			return false;

		if (sx == 2)
			return true;

		if (sx <= 31) {
			for (i = 0; i < 10; i++)
				if (sx == a[i])
					return true;

			return false;
		}
	}

	if (is_le_zero() || is_even())
		return false;

	if (b1 <= 0)
		lidia_error_handler("is_prime", "#tests <= 0");

	if (b1 > 9)
		b = 9;
	else
		b = b1;

	bigint erg;
	bigint H(37), Q(*this);
	Q.dec();

	bigint N_minus1(Q);

	while (Q.is_even()) {
		Q.divide_by_2();
		k++;
	}

	for (i = 0; i <= b; i++) {
		power_mod(erg, bigint(a[i]), Q, *this);

		if (!erg.is_one() && erg.compare(N_minus1)) {
			j = k;
			ok = 0;

			while ((j > 0) && !ok) {
				square(erg, erg);
				remainder(erg, erg, *this);

				if (!erg.compare(N_minus1)) ok = 1;

				j--;
			}

			if (!ok) {
				return false;
			}
		}
	}

	for (; i <= b1; i++) {
		if (!compare(H))
			return true;

		power_mod(erg, H, Q, *this);

		if (!erg.is_one() && erg.compare(N_minus1)) {
			j = k;
			ok = 0;

			while ((j > 0) && !ok) {
				square(erg, erg);
				remainder(erg, erg, *this);

				if (!erg.compare(N_minus1)) ok = 1;

				j--;
			}

			if (!ok) {
				return false;
			}
		}

		H = H.next_prime();
	}

	return true;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
