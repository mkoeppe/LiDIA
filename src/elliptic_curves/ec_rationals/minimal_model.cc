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
//	Author	: Nigel Smart (NS)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/elliptic_curves/ec_arith.h"
#include	"LiDIA/elliptic_curve_bigint.h"
#include	"LiDIA/elliptic_curve.h"
#include	"LiDIA/point_bigint.h"
#include	"LiDIA/curve_isomorphism.h"
#include	"LiDIA/minimal_model.h"
#include	"LiDIA/modular_operations.inl"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



// The following function computes the minimal model (over Z) of an
// elliptic_curve<bigrational>, returning also the isomorphism from
// the original to the minimal model so that points may be mapped
// between them.  For sample usage, see curve_isomorphism_appl.cc

void
minimal_model(elliptic_curve< bigint > & Emin, elliptic_curve< bigrational > & ER,
	      curve_isomorphism< bigrational, bigint > & iso)
{
	bigrational a1r, a2r, a3r, a4r, a6r;
	bigint a1i, a2i, a3i, a4i, a6i;
	bigint a1min, a2min, a3min, a4min, a6min;
	ER.get_ai(a1r, a2r, a3r, a4r, a6r);

	// Step 1: convert to an integral model:

	bigint v = a1r.denominator();
	v = lcm(v, a2r.denominator());
	v = lcm(v, a3r.denominator());
	v = lcm(v, a4r.denominator());
	v = lcm(v, a6r.denominator());
	bigint vi = v;
	a1i = (vi/a1r.denominator()) * a1r.numerator();
	vi *= v;
	a2i = (vi/a2r.denominator()) * a2r.numerator();
	vi *= v;
	a3i = (vi/a3r.denominator()) * a3r.numerator();
	vi *= v;
	a4i = (vi/a4r.denominator()) * a4r.numerator();
	vi *= v;
	a6i = (vi/a6r.denominator()) * a6r.numerator();

	// Step 2: find minimal c4, c6, and scale factor u:

	bigint b2 = a1i*a1i + 4 * a2i;
	bigint b4 = a1i*a3i + 2 * a4i;
	bigint b6 = a3i*a3i + 4 * a6i;
	bigint c4 = b2*b2 - 24 * b4;
	bigint c6 = b2*(36*b4 - b2*b2) - 216*b6;
	bigint delta = (c4*c4*c4 - c6*c6)/1728;

	bigint u;
	sort_vector< bigint > bad_p;

	minimise_c4c6(c4, c6, delta, u, bad_p);

	// Step 3: compute minimal coeffs and create minimal curve

	best_remainder(b2, -c6, 12);
	const bigint& b22 = b2*b2;
	b4 = (b22-c4)/24;
	b6 = (-b2*b22+36*b2*b4-c6)/216;

	a1min = (b2.is_odd() ? 1 : 0);
	a3min = (b6.is_odd() ? 1 : 0);
	a2min = (b2-a1min)/4;
	a4min = (b4-a1min*a3min)/2;
	a6min = (b6-a3min)/4;

	Emin.set_coefficients(a1min, a2min, a3min, a4min, a6min);
	// Step 4: compute the isomorphism parameters

	bigrational uv(u, v); // meaning u/v
	bigrational uv2 = uv*uv;
	bigrational uv3 = uv*uv2;
	bigrational s = (uv*a1min-a1r)/2;
	bigrational r = (uv2*a2min - a2r + s*a1r + s*s)/3;
	bigrational t = (uv3*a3min - a3r - r*a1r)/2;
	iso.init(ER, Emin, uv, r, s, t);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
