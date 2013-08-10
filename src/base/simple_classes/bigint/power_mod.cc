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
power_mod(bigint & res,
          const bigint & a,
          const bigint & n,
          const bigint & m, int err)
{
	if (a.is_one() || n.is_zero()) {
		res.assign_one();
		return;
	}

	bigint tmp_a(a), tmp_n(n), tmp_m(m);

	tmp_m.absolute_value();
	res.assign_one();

	while (tmp_n.sign()) {
		if (tmp_n.is_odd()) {
			multiply(res, res, tmp_a);
			remainder(res, res, tmp_m);
			if (res.is_negative())
				add(res, res, tmp_m);
		}
		tmp_n.divide_by_2();
		if (tmp_n.is_zero())
			break;
		square(tmp_a, tmp_a);
		remainder(tmp_a, tmp_a, tmp_m);
		if (res.is_negative())
			add(res, res, tmp_m);
	}
	if (n.is_negative()) {
		bigint tmp_u, tmp_v;
		tmp_n = xgcd(tmp_u, tmp_v, res, tmp_m);
		if (tmp_n.is_one()) {
			res.assign(tmp_u);
			if (res.is_negative())
				add(res, res, tmp_m);
		}
		else {
			res.assign(tmp_n);
			if (err) {
				warning_handler("bigint", "power_mod::inverse does not exist.");
			}
			else
				lidia_error_handler("bigint", "power_mod::inverse does not exist.");
		}
	}
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
