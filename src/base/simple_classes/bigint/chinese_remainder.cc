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



bigint chinese_remainder(const bigint & a, const bigint & m1,
			 const bigint & b, const bigint & m2)
{
	bigint mm1(abs(m1)), mm2(abs(m2)), t, u, d;

	d = xgcd_left (u, mm1, mm2); // determine lcm(m1, m2)

	if (!d.is_one() && !((a-b) % d).is_zero()) {
		lidia_error_handler("bigint", "chinese_remainder::solution does not exist");
		return 0;
	}
	divide(mm2, mm2, d);
	remainder(t, (b-a)/d*u, mm2);
	remainder(t, a + t*mm1, mm1*mm2);
	if (t.is_negative())
		add(t, t, mm1*mm2);

	return t;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
