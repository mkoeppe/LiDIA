// -*- C++ -*-
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
//	Author	: Markus Maurer (MM), Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_ELLIPTIC_CURVE_FLAGS_H_GUARD_
#define LIDIA_ELLIPTIC_CURVE_FLAGS_H_GUARD_


#ifndef LIDIA_LIDIA_H_GUARD_
# include	"LiDIA/LiDIA.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class elliptic_curve_flags
{
public:

	enum curve_parametrization {
		SHORT_W = 0, // Y^2 = X^3 + a4*X + a6 (affine)
		// Y^2 = X^3 + a4*X Z^4 + a6 Z^6 (projective)

		LONG_W = 1, // Y^2+a1*X*Y+a3*Y = X^3+a2*X^2+a4*X+a6 (affine)
		// Y^2 +a1*X*Y*Z+a3*Y*Z^3 =
		// X^3+a2*X^2*Z^2+a4*X*Z^4+a6 Z^6 (projective)

		GF2N_F = 2  // Y^2 +X*Y = X^3 + a2*X^2 + a6 (affine)
		// Y^2 + X*Y*Z = X^3 + a2*X^2*Z^2 + a6 Z^6 (projective)
	};

	enum curve_output_mode {
		SHORT = 0, // [a1, a2, a3, a4, a6, X] only, where X is either
		// "A" (for AFFINE) or "P" for (PROJECTIVE), respectively.
		// If ", X" is omitted the model is affine.
		LONG = 1, // as SHORT + b_i, c_j, delta,
		// and more in derived classes;
		PRETTY = 2, // equation form
		TEX = 3     // equation form in TeX format
	};

	enum curve_model {
		AFFINE = 0,
		PROJECTIVE = 1
	};
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_ELLIPTIC_CURVE_FLAGS_H_GUARD_
