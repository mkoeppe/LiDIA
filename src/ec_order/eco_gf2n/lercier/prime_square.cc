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
//	Author	: Volker Mueller
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



void eco_gf2n::prime_square(ff_pol & divpol_divisor, const ff_pol & fC,
			    const ff_element & A6_2)
{
	ff_pol numer, denom;
	ff_pol *table;
	ff_pol h;
	ff_element t;
	int d = fC.degree();

	isogeny_via_lercier(numer, denom, A6_2, A6);

#ifdef DEBUG
	ff_polmod fCpm;
	ff_pol hdebug;
	sqrt(hdebug, denom);

	fCpm.build(hdebug);
	eco_gf2n eco2(A6_2);
	eco2.compute_psi (hdebug, static_cast<lidia_size_t>(l), fCpm);
	assert(hdebug.is_zero());
#endif

	table = power_table(denom, d+1);
	h.assign_zero();

	for (int i = d; i > 0; i--) {
		fC.get_coefficient(t, i);
		add(h, h, t * table[d-i]);
		multiply(h, h, numer);
	}
	fC.get_coefficient(t, 0);
	add(divpol_divisor, h, t * table[d]);

	delete[] table;

#ifdef DEBUG
	fCpm.build(divpol_divisor);
	eco2.compute_psi (h, static_cast<lidia_size_t>(l*l), fCpm);
	assert(h.is_zero());
#endif
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
