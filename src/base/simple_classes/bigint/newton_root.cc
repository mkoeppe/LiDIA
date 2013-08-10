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



void
newton_root(bigint & b, const bigint & a, int n)
{
	bigint c, d;

	if (n < 0)
		lidia_error_handler("bigint", "newton_root::negative n.");

	if (a.is_negative())
		lidia_error_handler("bigint", "newton_root::negative a.");

	switch (n) {
	case 0:
		b.assign_one();
		break;
	case 1:
		b.assign(a);
		break;
	case 2:
		sqrt(b, a);
		break;
	default:
		b.assign_one();
		shift_left(b, b, static_cast<unsigned int>((a.bit_length() + n - 1) / n));
		do {
			power(c, b, n - 1);
			div_rem(c, d, a, c);
			subtract(c, b, c);
			divide(c, c, static_cast<long>(n));
			subtract(b, b, c);
		} while (c.sign() > 0);
		power(c, b, n);
		if (c.compare(a) > 0)
			dec(b);
		if (b.compare(3UL) == 0) {
			power(c, b, n);
			if (c.compare(a) > 0)
				dec(b);
		}
		break;
	}
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
