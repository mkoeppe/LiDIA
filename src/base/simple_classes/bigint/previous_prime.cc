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
//	Author	: Volker M"uller (VM)
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



bigint
bigint::previous_prime () const
{
	long zmod3, zmod5, zmod7, zmod11, zmod13, zmod17, zmod19;
	bigint prim;
	long sx;

	if (is_le_zero() || is_one())
		return (bigint(0));

	if (!longify(sx)) {
		// convert x to long
		if (sx <= 23) {
			if (sx == 2)
				return bigint(0);
			if (sx == 3)
				return bigint(2);
			if (sx <= 5)
				return bigint(3);
			if (sx <= 7)
				return bigint(5);
			if (sx <= 11)
				return bigint(7);
			if (sx <= 13)
				return bigint(11);
			if (sx <= 17)
				return bigint(13);
			if (sx <= 19)
				return bigint(17);
			return bigint(19);
		}
	}

	// prim = x-1; make prim odd

	prim.assign(*this);
	prim.dec();
	if (prim.is_even())
		prim.dec();

	// initialize modular counters
	zmod3 = remainder(prim, 3);
	zmod5 = remainder(prim, 5);
	zmod7 = remainder(prim, 7);
	zmod11 = remainder(prim, 11);
	zmod13 = remainder(prim, 13);
	zmod17 = remainder(prim, 17);
	zmod19 = remainder(prim, 19);

	// while not a prime number
	while (!prim.is_prime(4)) {
		do {
			// decrease prim by 2
			subtract(prim, prim, 2UL);

			zmod3 = (zmod3 + 1)%3;
			zmod5 = (zmod5 + 3)%5;
			zmod7 = (zmod7 + 5)%7;
			zmod11 = (zmod11 + 9)%11;
			zmod13 = (zmod13 + 11)%13;
			zmod17 = (zmod17 + 15)%17;
			zmod19 = (zmod19 + 17)%19;
		}
		// until it is not divisible by small primes
		while (!(zmod3 && zmod5 &&
			 zmod7 && zmod11 &&
			 zmod13 && zmod17 &&
			 zmod19));
	}
	return prim;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
