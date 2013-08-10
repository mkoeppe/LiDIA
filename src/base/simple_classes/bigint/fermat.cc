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
fermat(const bigint & n)
{
	bigint tmp_a, tmp_n, res;
	register int a = 2;

	if (n < 2)
		return 0;

	if ((n == 2) || (n == 3) || (n == 5) || (n == 7))
		return 1;

	tmp_n.assign(n);
	dec(tmp_n);

	while (a <= 7) {
		if (!remainder(n, a))
			return 0;
		else {
			tmp_a.assign(a);
			power_mod(res, tmp_a, tmp_n, n);
			if (!res.is_one())
				return 0;
			else {
				if (a == 2)
					a += 1;
				else
					a += 2;
			}
		}
	}
	return 1;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
