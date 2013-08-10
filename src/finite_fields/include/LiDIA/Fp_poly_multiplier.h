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
//	Author	: Victor Shoup (VS), Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_FP_POLY_MULTIPLIER_H_GUARD_
#define LIDIA_FP_POLY_MULTIPLIER_H_GUARD_



#ifndef LIDIA_LIDIA_H_GUARD_
# include	"LiDIA/LiDIA.h"
#endif
#ifndef LIDIA_FP_POLYNOMIAL_H_GUARD_
# include	"LiDIA/Fp_polynomial.h"
#endif
#ifndef LIDIA_FP_POLYNOMIAL_FFT_H_GUARD_
# include	"LiDIA/finite_fields/Fp_polynomial_fft.h"
#endif
#ifndef LIDIA_FP_POLY_MODULUS_H_GUARD_
# include	"LiDIA/Fp_poly_modulus.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



#ifndef HEADBANGER


// If you need to compute a * b % f for a fixed b, but for many a's
// (for example, computing powers of b modulo f), it is
// much more efficient to first build a Fp_poly_multiplier B for b,
// and then use the routine below.

class Fp_poly_multiplier
{
	const Fp_poly_modulus* poly_mod_ptr;

//public:
	Fp_polynomial b;
	bool use_FFT;
	fft_rep B1;
	fft_rep B2;


public:

	Fp_poly_multiplier () : poly_mod_ptr(0)
	{
		// b.assign_zero(); !!!modulus not set!!! */
	}

	Fp_poly_multiplier (const Fp_poly_multiplier& B);

	Fp_poly_multiplier (const Fp_polynomial& x, const Fp_poly_modulus& F)
	{
		build(x, F);
	}

	void build (const Fp_polynomial& x, const Fp_poly_modulus& F);

	Fp_poly_multiplier & operator = (const Fp_poly_multiplier &B);

	const Fp_polynomial& multiplier () const
	{
		return b;
	}

	friend void multiply (Fp_polynomial& x, const Fp_polynomial& a,
			      const Fp_poly_multiplier& B, const Fp_poly_modulus& F);
    	// x = (a * b) % f


	friend void update_map (base_vector< bigint > & x, const base_vector< bigint > & a,
				const Fp_poly_multiplier& B, const Fp_poly_modulus& F);
    	//see compose.cc
    	// computes (a, b), (a, (b*X)%f), ..., (a, (b*X^{n-1})%f),
	// where (,) denotes the vector inner product.
        // This is really a "transposed" MulMod by B.

};

typedef Fp_poly_multiplier poly_multiplier;


#endif



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_FP_POLY_MULTIPLIER_H_GUARD_

