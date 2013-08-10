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
//	Author	: Volker Mueller (VM), Markus Maurer (MM)
//	Changes	: See CVS log
//
//==============================================================================================



#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/eco_gf2n.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// compute a list d[0, ..., d[0]] of possible values for d
// see page 50, PHD thesis VM
//

void
eco_gf2n::compute_d_list (lidia_size_t * & dlist)
{
	lidia_size_t deg, i, d, h;
	bool iss;

	if (is_elkies())
		deg = l-1;
	else
		deg = l+1;

	dlist = new lidia_size_t [deg+1];
	i=1;

	if (jacobi(pn, bigint(l)) == 1)
		iss = true;
	else iss = false;

	for (d = 2; d <= deg; d++) {
		if (deg % d == 0) {
			h = deg / d;
			if ((iss && !((h*(d-1)) & 1)) || (!iss && ((h*(d-1)) & 1)))
				dlist[i++] = d;
		}
	}
	dlist[0] = i-1;
}



// Compute a list ev[0, ..., ev[0]] of possible values for the eigenvalue.
// The input is list (dlist) of possible splitting degrees d of the
// modular polynomial in increasing order. We know that an eigenvalue alpha
// satiesfies order(alpha^2/q) == d in GF(l).
//

void
eco_gf2n::compute_ev_list (lidia_size_t * dlist, lidia_size_t * & ev)
{
	lidia_size_t i, j, k, n;
	udigit ord;
	ff1_element alpha, qinv;

	n = static_cast<lidia_size_t>(l-1)/2;

	ev = new lidia_size_t[n+1];
	i = 1;
	ff1_element::set_characteristic(l);

	invert(qinv, q);

	for (j = 1; j <= n; j++) {
		// compute order of j^2/p in GF(l)
		//
		alpha = ff1(j*j);
		multiply(alpha, alpha, qinv);
		ord = alpha.multiplicative_order();

		// search for order in list of possible degrees
		//
		k = 1;
		while (ord > static_cast<udigit>(dlist[k]) && k < dlist[0])
			k++;

		// if order is one of the degress, j is a candidate
		// for the eigenvalue
		//
		if (ord == static_cast<udigit>(dlist[k]))
			ev[i++] = j;
	}
	ev[0] = i-1;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
