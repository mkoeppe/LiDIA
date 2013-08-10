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
//	Author	: Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/elliptic_curves/point_operations_bigint.h"
#include	"LiDIA/elliptic_curve_bigint.h"
#include	"LiDIA/point_bigint.h"
#include	"LiDIA/elliptic_curve_flags.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// constructor / destructor
//

point_operations< bigint >::point_operations()
{
	debug_handler("point_operations< bigint >",
		      "point_operations()");

	_add = NULL;
	_negate = NULL;
	_mult_by_2 = NULL;
}



point_operations< bigint >::~point_operations()
{
	debug_handler("point_operations< bigint >",
		      "~point_operations()");
}



//
// assigment
//

void
point_operations< bigint >::set(elliptic_curve_flags::curve_parametrization cp,
				elliptic_curve_flags::curve_model cm,
				const bigint& characteristic)
{
	switch (cp) {
	case elliptic_curve_flags::SHORT_W:
		_add = add_swnf_projective;
		_negate = negate_swnf_projective;
		_mult_by_2 = mult_by_2_swnf_projective;
		break;

	case  elliptic_curve_flags::LONG_W:
		_add = add_lwnf_projective;
		_negate = negate_lwnf_projective;
		_mult_by_2 = mult_by_2_lwnf_projective;
		break;
	default: lidia_error_handler("point_operations< bigint >::set",
				     "Projective. Wrong type.");
	}
}



point_operations< bigint > &
point_operations< bigint >::operator = (const point_operations< bigint > & p)
{
	if (this != &p) {
		_add = p._add;
		_negate = p._negate;
		_mult_by_2 = p._mult_by_2;
	}
	return *this;
}



//
// Implementation of arithmetic on point<bigint>
// used in point_bigint.c.
//

void
add_lwnf_projective (point< bigint > & erg, const point< bigint > & wp1,
		     const point< bigint > & wp2)
{
	// make sure that both points are on the same curve
	if (*(wp1.ec) != *(wp2.ec)) {
		lidia_error_handler("point< bigint >",
				    "add_lwnf_projective::Can't add points on different curves");
	}
	// NOTE: we don't do any reductions here, but rely on assignment to do it
	// check special cases first
	if (wp1.z.is_zero()) {
		erg.assign(wp2);
		return;
	}
	if (wp2.z.is_zero()) {
		erg.assign(wp1);
		return;
	}

	point< bigint > minusQ;

	negate(minusQ, wp2);
	if (is_equal(wp1, minusQ)) {       // zero
		erg.assign_zero(*(wp1.ec));
		return;
	}
	if (is_equal(wp1, wp2)) {
		mult_by_2_lwnf_projective(erg, wp2);
		return;
	}

	// we now have genuine work to do
	// let's set up some local variables to avoid repeated references
	// coefficients
	bigint a1, a2, a3, a4, a6;
	(wp1.ec)->get_ai(a1, a2, a3, a4, a6);
	// coordinates of P
	const bigint& X1 = wp1.x;
	const bigint& Y1 = wp1.y;
	const bigint& Z1 = wp1.z;
	// coordinates of Q
	const bigint& X2 = wp2.x;
	const bigint& Y2 = wp2.y;
	const bigint& Z2 = wp2.z;
	const bigint& Z12 = Z1 * Z2;

	const bigint& L = - Y2 * Z1 + Y1 * Z2; // lambda
	const bigint& M = - X2 * Z1 + X1 * Z2; // mu   
	const bigint& N = - Y1 * X2 + Y2 * X1; // nu  

	const bigint& Mz = M * M * Z12;

	const bigint& t = L * L * Z12 + M * (a1 * L * Z12 - M * (a2 * Z12
								 + X1 * Z2 + X2 * Z1));
	*(erg.ec) = *(wp1.ec);
	erg.x = M * t;
	erg.y = - (t * (L + a1 * M) +   Mz * (N + a3 * M));
	erg.z = M * Mz;
	erg.height = -1; erg.ord = 0;
	erg.reduce();
}



void
add_swnf_projective (point< bigint > & erg, const point< bigint > & wp1,
		     const point< bigint > & wp2)
{
	// make sure that both points are on the same curve
	if (*(wp1.ec) != *(wp2.ec)) {
		lidia_error_handler("point< bigint >",
				    "add_swnf_projective::Can't add points on different curves");
	}
	// NOTE: we don't do any reductions here, but rely on assignment to do it
	// check special cases first
	if (wp1.z.is_zero()) {
		erg.assign(wp2);
		return;
	}
	if (wp2.z.is_zero()) {
		erg.assign(wp1);
		return;
	}

	point< bigint > minusQ = wp2;

	negate(minusQ, wp2);
	if (is_equal(wp1, minusQ)) {       // zero
		erg.assign_zero(*(wp1.ec));
		return;
	}
	if (is_equal(wp1, wp2)) {
		mult_by_2_swnf_projective(erg, wp2);
		return;
	}

	// we now have genuine work to do
	// let's set up some local variables to avoid repeated references
	// coefficients
	bigint a4 = (wp1.ec)->get_a4();
	bigint a6 = (wp1.ec)->get_a6();
	// coordinates of P
	const bigint& X1 = wp1.x;
	const bigint& Y1 = wp1.y;
	const bigint& Z1 = wp1.z;
	// coordinates of Q
	const bigint& X2 = wp2.x;
	const bigint& Y2 = wp2.y;
	const bigint& Z2 = wp2.z;
	const bigint& Z12 = Z1 * Z2;

	const bigint& L = - Y2 * Z1 + Y1 * Z2; // lambda
	const bigint& M = - X2 * Z1 + X1 * Z2; // mu   
	const bigint& N = - Y1 * X2 + Y2 * X1; // nu  

	const bigint& Mz = M * M * Z12;

	const bigint& t = L * L * Z12 - M * M *(X1 * Z2 + X2 * Z1);
	*(erg.ec) = *(wp1.ec);
	erg.x = M * t;
	erg.y = - (t * L  +   Mz * N);
	erg.z = M * Mz;
	erg.height = -1; erg.ord = 0;
	erg.reduce();
}



void
mult_by_2_lwnf_projective (point< bigint > & erg, const point< bigint > & wp)
{
	// do trivial cases
	point< bigint > minuswp;


	negate(minuswp, wp);
	if (wp.z.is_zero() || is_equal(wp, minuswp)) {
		erg.assign_zero(*(wp.ec));
		return;
	}

	bigint a1, a2, a3, a4, a6;
	(wp.ec)->get_ai(a1, a2, a3, a4, a6);

	bigint x = wp.x, y = wp.y, z = wp.z;

	const bigint& Zsq = z * z;
	const bigint& L = 3*x*x + 2*a2*x*z + a4*Zsq - a1*y*z;
	const bigint& M = 2 * y  +  a1 * x  +  a3 * z;
	const bigint& Mz = M * z;
	const bigint& N = -x*x*x - a3*y*Zsq  + a4*x*Zsq  + 2*a6*z*Zsq;
	const bigint& t = L*L  +  Mz*(a1*L  -  M*(a2*z  + 2*x));

	*(erg.ec) = *(wp.ec);
	erg.x = t * Mz;
	erg.y = - (L * t +  Mz * (a1 * t+ M * (N  + a3 * Mz * z)));
	erg.z = Mz * Mz * Mz;
	if (wp.height.is_ge_zero()) {
		erg.height = wp.height*4;
	}
	else {
		erg.height = -1;
	}
	if ((wp.ord%2) == 0) {
		erg.ord = wp.ord/2;
	}
	else {
		erg.ord = wp.ord;
	}

	erg.reduce();
}



void
mult_by_2_swnf_projective (point< bigint > & erg, const point< bigint > & wp)
{
	// do trivial cases
	point< bigint > minuswp;


	negate(minuswp, wp);
	if (wp.z.is_zero() || is_equal(wp, minuswp)) {
		erg.assign_zero(*(wp.ec));
		return;
	}

	bigint x = wp.x, y = wp.y, z = wp.z;

	bigint a4 = (wp.ec)->get_a4();
	bigint a6 = (wp.ec)->get_a6();

	const bigint& Zsq = z * z;
	const bigint& L = 3*x*x + a4*Zsq;
	const bigint& M = 2 * y;
	const bigint& Mz = M * z;
	const bigint& N = -x*x*x + a4*x*Zsq  + 2*a6*z*Zsq;
	const bigint& t = L*L-Mz*M*2*x;

	*(erg.ec) = *(wp.ec);
	erg.x = t * Mz;
	erg.y = - (L*t +  Mz*M*N);
	erg.z = Mz * Mz * Mz;

	if (wp.height.is_ge_zero()) {
		erg.height = wp.height*4;
	}
	else {
		erg.height = -1;
	}
	if ((wp.ord%2) == 0) {
		erg.ord = wp.ord/2;
	}
	else {
		erg.ord = wp.ord;
	}

	erg.reduce();
}



void
negate_swnf_projective (point< bigint > & erg, const point< bigint > & P)
{
	erg = P;
	erg.y = - erg.y;
	erg.height = P.height;
	erg.ord = P.ord;
}



void
negate_lwnf_projective (point< bigint > & erg, const point< bigint > & P)
{
	erg = P;
	erg.y = -P.y - (P.ec)->get_a1() * P.x - (P.ec)->get_a3();
	erg.height = P.height;
	erg.ord = P.ord;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
