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
//	Author	: Thomas Papanikolaou (TP)
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



#define prime_table_length 100



static int prime_table[prime_table_length] = {
	2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
	67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137,
	139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211,
	223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283,
	293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379,
	383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461,
	463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541};



long
bigint::is_power(bigint & b) const
{
	bigint c;
	int i = 0;

	if (!is_negative()) {
		int l = bit_length();

		for (; i < prime_table_length && prime_table[i] <= l; i++) {
			newton_root(b, *this, prime_table[i]);

			power(c, b, prime_table[i]);
			if (compare(c) == 0)
				return prime_table[i];
		}
		for (i = 545; i <= l; i += 4) {
			// testing 6k-1 and 6k+1 for  k >= 91
			newton_root(b, *this, i);
			power(c, b, i);
			if (compare(c) == 0)
				return i;
			i += 2; // increase from 6k-1 to 6k+1. Do again.
			newton_root(b, *this, i);
			power(c, b, i);
			if (compare(c) == 0)
				return i;
		}
	}
	else {
		i++; // negative numbers can't be even powers
		bigint a2(*this);
		a2.abs();

		int l = a2.bit_length();

		for (; i < prime_table_length && prime_table[i] <= l; i++) {
			newton_root(b, a2, prime_table[i]);
			power(c, b, prime_table[i]);
			if (a2.compare(c) == 0) {
				b.negate();
				return prime_table[i];
			}
		}
		for (i = 545; i <= l; i += 4) {
			// testing 6k-1 and 6k+1 for  k >= 91
			newton_root(b, a2, i);
			power(c, b, i);
			if (a2.compare(c) == 0) {
				b.negate();
				return i;
			}
			i += 2; // increase from 6k-1 to 6k+1. Do again.
			newton_root(b, a2, i);
			power(c, b, i);
			if (a2.compare(c) == 0) {
				b.negate();
				return i;
			}
		}
	}
	return 0;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
