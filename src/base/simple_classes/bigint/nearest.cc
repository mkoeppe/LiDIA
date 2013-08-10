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
nearest(bigint & z, const bigint & u, const bigint & v)
{
	if (v.is_zero())
		lidia_error_handler("bigint", "nearest()::division by zero.");

	bigint tmp, tmp_u(u), tmp_v(v);

	if (tmp_v.is_lt_zero()) {
		tmp_v.negate();
		tmp_u.negate();
	}

	div_rem(z, tmp, tmp_u, tmp_v);
	tmp.multiply_by_2();
	if (tmp.abs_compare(tmp_v) >= 0) {
		if (tmp_u.is_ge_zero())
			inc(z);
		else
			dec(z);
	}
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
