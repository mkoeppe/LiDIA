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
//	Author	: Nigel Smart (NS)
//		  Adaption of John Cremona's code
//                (some of which itself is based on code of
//                Oisin McGuiness and Richard Pinch)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_ELLIPTIC_CURVE_REP_BIGINT_H_GUARD_
#define LIDIA_ELLIPTIC_CURVE_REP_BIGINT_H_GUARD_


#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_BASE_VECTOR_H_GUARD_
# include	"LiDIA/base_vector.h"
#endif
#ifndef LIDIA_SORT_VECTOR_H_GUARD_
# include	"LiDIA/sort_vector.h"
#endif

// To force the use of point < bigint > and point_operations < bigint > !
#ifndef LIDIA_POINT_OPERATIONS_BIGINT_H_GUARD_
# include	"LiDIA/elliptic_curves/point_operations_bigint.h"
#endif
#ifndef LIDIA_POINT_BIGINT_H_GUARD_
# include	"LiDIA/point_bigint.h"
#endif

#ifndef LIDIA_BASE_ELLIPTIC_CURVE_REP_H_GUARD_
# include	"LiDIA/elliptic_curves/base_elliptic_curve_rep.h"
#endif
#ifndef LIDIA_COMPLEX_PERIODS_H_GUARD_
# include	"LiDIA/complex_periods.h"
#endif
#ifndef LIDIA_ELLIPTIC_CURVE_FLAGS_H_GUARD_
# include	"LiDIA/elliptic_curve_flags.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



// Class for holding a rational elliptic curve (minimal model only)

template< class T > class elliptic_curve;
template< class T > class elliptic_curve_rep;
class reduction_type;


template<>
class elliptic_curve_rep< bigint > : public base_elliptic_curve_rep< bigint >
{
#if defined(_MSC_VER)
	friend template elliptic_curve< bigint >;
#else
	friend class elliptic_curve< bigint >;
#endif

	reduction_type* reduct_array; //just a plain c++ array
	bigint N; //the conductor
	sort_vector< bigint > bad_primes; //Sorted vector of bad_primes
	int  conncomp; // number of real components (1 or 2)
	base_vector< point < bigint > > torsion; // Torsion pts (size = 0 if not computed yet)
	complex_periods periods; // norm_code = -1 if nothing computed yet
	double height_constant; // h(P) <= nth(P)+ height_constant
				// = -1 if not computed yet

	bool lock; // for destructor of elliptic_curve.
	//
	// initialization
	//
private:
	void reset();
	void init();
	void init(const bigint& aa1, const bigint& aa2, const bigint& aa3,
		  const bigint& aa4, const bigint& aa6);
	void init(const bigint& cc4, const bigint& cc6);
	void tate();


	//
	// constructors / destructor
	//
public:
	elliptic_curve_rep();
	// construct by Tate's algoritm
        // arg E need not be minimal, but the reduced form will be
	// Will minimalize the curve
	elliptic_curve_rep(const bigint& x1, const bigint& x2, const bigint& x3,
			   const bigint& x4, const bigint& x6,
			   elliptic_curve_flags::curve_model m = elliptic_curve_flags::PROJECTIVE);

	// NB The following constructor constructs a curve from given c4, c6
	// invariants, first checking that they are valid.  This must NOT be
	// confused with the two-parameter base_elliptic_curve constructor which sets
	// a4 and a6, with a1 = a2 = a3 = 0.  This one is more useful for curves/Z
	// where we need to be able to construct them via periods and c4, c6

	elliptic_curve_rep(const bigint& cc4, const bigint& cc6,
			   elliptic_curve_flags::curve_model m = elliptic_curve_flags::PROJECTIVE);

	elliptic_curve_rep(const elliptic_curve_rep< bigint > & E);
	~elliptic_curve_rep();

	//
	// assignment / initialization
	//
	elliptic_curve_rep< bigint > & operator = (const elliptic_curve_rep< bigint > & E);
};



// The following checks whether given c4, c6 are valid elliptic
// curve invariants.   It returns the discriminant via the third
// parameter as it has to compute it anyway

int valid_c4c6(const bigint& c4, const bigint& c6, bigint& disc);

// The following uses Kraus-Laska-Connell; c4 and c6 are input
// and will be overwritten; disc will hold the minimal disc,
// bad_p its prime factors, and u the scaling factor needed.
// It is used in elliptic_curve_bigint.cc and in curve_isomorphism.cc

void minimise_c4c6(bigint& c4, bigint& c6,
		   bigint& disc, bigint& u, sort_vector< bigint > & bad_p);



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_ELLIPTIC_CURVE_REP_BIGINT_H_GUARD_
