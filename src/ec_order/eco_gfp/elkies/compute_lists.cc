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
#include	"LiDIA/eco_prime.h"
#include	"LiDIA/Fp_polynomial.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// compute a list d[0, ..., d[0]] of possible values for d
// see page 50, PHD thesis VM
//

void
eco_prime::compute_d_list (lidia_size_t * & dlist)
{
	lidia_size_t deg, i, d, h;
	bool iss;

	if (is_elkies())
		deg = l-1;
	else
		deg = l+1;

	dlist = new lidia_size_t [deg+1];
	i = 1;

	if (jacobi(pn, bigint(l)) == 1)
		iss = true;
	else
		iss = false;

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
// The list ev is ascendingly sorted for absolute values, it has the
// form ev[1] = 1st candidate, ev[2] = - ev[1] if valid, otherwise
// ev[2] > ev[1] and so on.
// The input is list (dlist) of possible splitting degrees d of the
// modular polynomial in increasing order. We know that an eigenvalue alpha
// satiesfies order(alpha^2/q) == d in GF(l).
//

void
eco_prime::compute_ev_list (lidia_size_t * dlist, lidia_size_t * & ev)
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



//-----------------------------------------------------------------------
// Now we add the new idea for primes l == 1 mod 4, see Maurer & Mueller
// "Finding the Eigenvalue in Elkies algorithm".
// For primes l = 1 mod 4, we can filter out essentially (l-1)/2 of the
// l-1 possibilities with a resultant computation.

void eco_prime::refine_ev_list(const ff_pol & fC, lidia_size_t * & ev_list)
{
	if (l % 4 == 3)
		return;

	lidia_size_t i, j, d = (l-1)/2;
	lidia_size_t r, s, a, a2;
	bigint h;
	ff_pol curve;

	curve.set_modulus(A.modulus());
	curve.set_coefficient(3);
	curve.set_coefficient(A.mantissa(), 1);
	curve.set_coefficient(B.mantissa(), 0);

	resultant(h, fC, curve);
	j = jacobi(h, pn);

	for (i = 1; i <= ev_list[0]; i++) {
		a = ev_list[i];
		a2 = (a * a) % l;
		r = 0;

		for (s = 1; s <= d; s++) {
			if ((s*a) % l > static_cast<unsigned int>(d))
				r++;
			if ((s*a2) % l > static_cast<unsigned int>(d))
				r++;
		}
		if (!((j == 1 && r % 2 == 0) || (j == -1 && r % 2 == 1)))
			ev_list[i] = - ev_list[i];
	}

	j = 1;
	for (i = 1; i <= ev_list[0]; i++)
		if (ev_list[i] > 0)
			ev_list[j++] = ev_list[i];
	ev_list[0] = j-1;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
