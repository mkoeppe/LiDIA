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
//	Author	: Birgit Henhapl (BHE), Markus Maurer (MM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/gf_element.h"
#include	"LiDIA/elliptic_curves/point_operations.h"
#include	"LiDIA/elliptic_curves/base_point.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// Start of implementation of class point_operations<gf_element>.
//

//
// constructor / destructor
//

point_operations< gf_element >::point_operations ()
{
	debug_handler("point_operations< gf_element >",
		      "point_operations()");

	_add = NULL;
	_negate = NULL;
	_mult_by_2 = NULL;
}



point_operations< gf_element >::~point_operations ()
{
	debug_handler("point_operations< gf_element >",
		      "~point_operations()");
}


//
// assigment
//


void
point_operations< gf_element >::set (elliptic_curve_flags::curve_parametrization cp,
				      elliptic_curve_flags::curve_model cm,
				      const bigint & characteristic)
{

	if (cm == elliptic_curve_flags::AFFINE) {
		switch (cp) {
		case elliptic_curve_flags::SHORT_W:
			_add = add_swnf_affine;
			_negate = negate_swnf_affine;
			_mult_by_2 = mult_by_2_swnf_affine;
			break;

		case  elliptic_curve_flags::LONG_W:
			_add = add_lwnf_affine;
			_negate = negate_lwnf_affine;
			_mult_by_2 = mult_by_2_lwnf_affine;
			break;

			// Parametrization GF2N_F only allowed for characteristic 2.
			//

		case  elliptic_curve_flags::GF2N_F:
			if (characteristic != 2) {
				lidia_error_handler ("point_operations< gf_element >::set",
						     "GF2N_F functions only valid for characteristic 2");
			}
			else {
				_add = add_gf2nf_affine;
				_negate = negate_gf2nf_affine;
				_mult_by_2 = mult_by_2_gf2nf_affine;
			}
			break;
		default: lidia_error_handler("point_operations< gf_element >::set",
					     "Wrong type.");
		}
	}
	else { // projective representation
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

			// Parametrization GF2N_F only allowed for characteristic 2.
			//

		case  elliptic_curve_flags::GF2N_F:
			if (characteristic != 2) {
				lidia_error_handler ("point_operations< gf_element >::set",
						     "GF2N_F functions only valid for characteristic 2");
			}
			else {
				_add = add_gf2nf_projective;
				_negate = negate_gf2nf_projective;
				_mult_by_2 = mult_by_2_gf2nf_projective;
			}
			break;
		default: lidia_error_handler("point_operations< gf_element >::set",
					     "Wrong type.");
		}
	}
}



point_operations< gf_element > &
point_operations< gf_element >::operator = (const point_operations< gf_element > & p)
{
	if (this != &p) {
		_add = p._add;
		_negate = p._negate;
		_mult_by_2 = p._mult_by_2;
	}
	return *this;
}



//*****************************************************************
//
// Implementation of arithmetic for points over finite fields.
// First we describe the operations for GF2N form

// wp1 and wp2 must be different

void
add_gf2nf_affine (base_point< gf_element > & erg,
		  const base_point< gf_element > & wp1,
		  const base_point< gf_element > & wp2)
{
	gf_element lambda, h1, h2;

	add(h1, wp1.x, wp2.x);
	if (h1.is_zero()) {
		erg.assign_zero(wp1.ec);
		return;
	}

	h1.invert();
	add(h2, wp1.y, wp2.y);
	multiply(lambda, h2, h1);
	square(h1, lambda);
	add(h1, h1, lambda);
	add(h1, h1, wp1.x);
	add(h1, h1, wp2.x);
	add(h1, h1, wp2.ec.get_a2());

	add(h2, wp1.x, h1);
	multiply(h2, lambda, h2);
	add(h2, h2, h1);
	add(h2, h2, wp1.y);
#ifdef DEBUG
	erg.assign(h1, h2, wp1.ec);
#else
	erg.assign_no_test(h1, h2, wp1.ec);
#endif
}



//---------------------------------------------------------------

void
mult_by_2_gf2nf_affine (base_point< gf_element > & erg,
			const base_point< gf_element > & wp)
{
	debug_handler("point_operations.c", "mult_by_2_gf2nf_affine");

	gf_element h2, h1, lambda;

	if (wp.x.is_zero()) {
		erg.assign_zero(wp.ec);
		return;
	}

	invert(lambda, wp.x);
	multiply(lambda, lambda, wp.y);
	add(lambda, lambda, wp.x);
	square(h1, lambda);
	add(h1, h1, lambda);
	add(h1, h1, wp.ec.get_a2());

	square(h2, wp.x);
	add(lambda, lambda, 1);
	multiply(lambda, lambda, h1);
	add(h2, lambda, h2);
#ifdef DEBUG
	erg.assign(h1, h2, wp.ec);
#else
	erg.assign_no_test(h1, h2, wp.ec);
#endif
}



//----------------------------------------------------------------

void
negate_gf2nf_affine (base_point< gf_element > & erg,
		     const base_point< gf_element > & P)
{
	debug_handler("point_operations.c", "negate_gf2nf_affine");
	erg = P;
	add(erg.y, P.x, P.y);
}



//-----------------------------------------------------------------
// wp1 and wp2 must be different

void
add_gf2nf_projective (base_point< gf_element > & erg,
		      const base_point< gf_element > & wp1,
		      const base_point< gf_element > & wp2)
{
	debug_handler("point_operations.c", "add_by_2_gf2nf_projective");

	if (wp1.is_zero()) {
		erg.assign(wp2);
		return;
	}

	if (wp2.is_zero()) {
		erg.assign(wp1);
		return;
	}

	gf_element h1 (wp1.x.get_field()), h2(wp1.x.get_field()),
		h3(wp1.x.get_field()), h4(wp1.x.get_field()), h5(wp1.x.get_field()),
		h6(wp1.x.get_field()), h7(wp1.x.get_field());

	if (!wp2.z.is_one()) {
		square(h7, wp2.z);
		multiply(h1, wp1.x, h7); // U_0 = X_0 * Z_1^2
		multiply(h7, wp2.z, h7);
		multiply(h2, wp1.y, h7); // S_0 = Y_0 * Z_1^3
	}
	else {
		h1.assign(wp1.x);
		h2.assign(wp1.y);
	}

	if (!wp1.z.is_one()) {
		square(h7, wp1.z);
		multiply(h4, wp2.x, h7); // U_1 = X_1 * Z_0^2
		multiply(h7, wp1.z, h7);
		multiply(h5, wp2.y, h7); // S_1 = Y_1 * Z_0^3
	}
	else {
		h4.assign(wp2.x);
		h5.assign(wp2.y);
	}

	add(h4, h1, h4); // W = U_0 + U_1
	add(h5, h2, h5); // R = S_0 + S_1

	if (h4.is_zero()) {
		if (h5.is_zero()) {
			mult_by_2_gf2nf_projective(erg, wp1);
		}
		else {
			erg.assign_zero(wp1.ec);
		}
		return;
	}

	multiply(h2, wp1.z, h4); // L = Z_0 * W

	multiply(h1, wp2.x, h5);
	multiply(h6, h2, wp2.y);
	add(h1, h1, h6); // V = R * X_1 + L * Y_1

	multiply(h6, h2, wp2.z); // Z_2 = L * Z_1

	add(h3, h6, h5); // T = R + Z_2

	square(h7, h4);
	multiply(h7, h7, h4);
	multiply(h4, h3, h5);
	add(h4, h4, h7);
	square(h7, h6);
	multiply(h7, h7, wp1.ec.get_a2());
	add(h4, h4, h7); // X_2 = a2 * Z_2^2 + T * R + W^3

	multiply(h7, h3, h4);
	square(h2, h2);
	multiply(h2, h2, h1);
	add(h2, h2, h7); // Y_2 = T * X_2 + V * L^2

#ifdef DEBUG
	erg.assign(h4, h2, h6, wp1.ec);
#else
	erg.assign_no_test(h4, h2, h6, wp1.ec);
#endif
}



//-----------------------------------------------------------------

void
mult_by_2_gf2nf_projective (base_point< gf_element > & erg,
			    const base_point< gf_element > & wp)
{
	debug_handler("point_operations.c", "mult_by_2_gf2nf_projective");

	if (wp.is_zero()) {
		erg.assign_zero(wp.ec);
		return;
	}

	gf_element  h1(wp.x.get_field()), h2(wp.x.get_field()),
		h3(wp.x.get_field()), h4(wp.x.get_field());

	square(h4, wp.z);
	multiply(h3, h4, wp.x); // Z_2 = X_1 * Z_1^2

	multiply(h4, h4, wp.ec.get_sqrt4_a6()); // fourth root of a6!!
	add(h4, h4, wp.x);
	square(h4, h4);
	square(h4, h4); // X_2 = (X_1 + root[4](a6) * Z_1^2)^4

	square(h1, wp.x);
	multiply(h2, wp.y, wp.z);
	add(h2, h2, h1);
	add(h2, h2, h3); // U = Z_2 + X_1^2 + Y_1 * Z_1

	square(h1, h1);
	multiply(h1, h1, h3);
	multiply(h2, h2, h4);
	add(h1, h1, h2); // Y_2 = X_1^4*Z_2 + U * X_2

#ifdef DEBUG
	erg.assign(h4, h1, h3, wp.ec);
#else
	erg.assign_no_test(h4, h1, h3, wp.ec);
#endif
}



//-----------------------------------------------------------------

void
negate_gf2nf_projective (base_point< gf_element > & erg,
			 const base_point< gf_element > & P)
{
	debug_handler("point_operations.c", "negate_gf2nf_projective");
	if (&erg != &P) {
		erg = P;
		multiply(erg.y, P.x, P.z);
		add(erg.y, P.y, erg.y);
	}
	else {
		add(erg.y, erg.y, erg.x * erg.z);
	}
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
