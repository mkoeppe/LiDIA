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
//	Author	: Nigel Smart
//		  Adaption of John Cremona's code
//                (some of which itself is based on code of
//                Oisin McGuiness and Richard Pinch)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_ELLIPTIC_CURVE_BIGINT_H_GUARD_
#define LIDIA_ELLIPTIC_CURVE_BIGINT_H_GUARD_


#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_BIGRATIONAL_H_GUARD_
# include	"LiDIA/bigrational.h"
#endif
#ifndef LIDIA_BIGCOMPLEX_H_GUARD_
# include	"LiDIA/bigcomplex.h"
#endif
#ifndef LIDIA_BASE_VECTOR_H_GUARD_
# include	"LiDIA/base_vector.h"
#endif
#ifndef LIDIA_SORT_VECTOR_H_GUARD_
# include	"LiDIA/sort_vector.h"
#endif

// To force the use of point < bigint > !
#ifndef LIDIA_POINT_BIGINT_H_GUARD_
# include	"LiDIA/point_bigint.h"
#endif

// To show the compiler that base_elliptic_curve < bigint >::get_point_operations
// returns the specialization point_operations < bigint > !
#ifndef LIDIA_POINT_OPERATIONS_BIGINT_H_GUARD_
# include	"LiDIA/elliptic_curves/point_operations_bigint.h"
#endif

#ifndef LIDIA_BASE_ELLIPTIC_CURVE_H_GUARD_
# include	"LiDIA/elliptic_curves/base_elliptic_curve.h"
#endif
#ifndef LIDIA_COMPLEX_PERIODS_H_GUARD_
# include	"LiDIA/complex_periods.h"
#endif
#ifndef LIDIA_KODAIRA_CODE_H_GUARD_
# include	"LiDIA/kodaira_code.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



// Class for holding a rational elliptic curve (minimal model only)

template< class T > class elliptic_curve;


template<>
class elliptic_curve< bigint > : public base_elliptic_curve< bigint >
{

private:

	void make_torsion();


public:
	//
	// constructors / destructor
	//

	elliptic_curve();
	// construct by Tate's algoritm
        // arg E need not be minimal, but the reduced form will be
	// Will minimalize the curve
	elliptic_curve(const bigint& x1, const bigint& x2, const bigint& x3,
		       const bigint& x4, const bigint& x6,
		       elliptic_curve_flags::curve_model m = elliptic_curve_flags::PROJECTIVE);

	// NB The following constructor constructs a curve from given c4, c6
	// invariants, first checking that they are valid.  This must NOT be
	// confused with the two-parameter base_elliptic_curve constructor which sets
	// a4 and a6, with a1 = a2 = a3 = 0.  This one is more useful for curves/Z
	// where we need to be able to construct them via periods and c4, c6

	elliptic_curve(const bigint& cc4, const bigint& cc6,
		       elliptic_curve_flags::curve_model m = elliptic_curve_flags::PROJECTIVE);

	elliptic_curve(const elliptic_curve< bigint > & E);
	virtual ~elliptic_curve();


	//
	// assignment / initialization
	//

	elliptic_curve< bigint > & operator = (const elliptic_curve< bigint > & E);
	void assign (const elliptic_curve< bigint > & E);

	void set_coefficients(const bigint & x4, const bigint & x6,
			      elliptic_curve_flags::curve_model m =
			      elliptic_curve_flags::PROJECTIVE);

	void set_coefficients(const bigint & x1, const bigint & x2,
			      const bigint & x3, const bigint & x4, const bigint & x6,
			      elliptic_curve_flags::curve_model m =
			      elliptic_curve_flags::PROJECTIVE);

	//
	// Access Function
	//

	bigrational j_invariant() const;
	bigint get_conductor() const;
	int get_ord_p_discriminant(const bigint& p) const;
	int get_ord_p_N(const bigint& p) const;
	int get_ord_p_j_denom(const bigint& p) const;
	int get_c_p(const bigint& p) const;
	int get_product_cp() const;
	Kodaira_code get_Kodaira_code(const bigint& p) const;

	int get_no_tors();
	base_vector< point < bigint > > get_torsion();

	bigcomplex get_omega_1();
	bigcomplex get_omega_2();

	complex_periods get_periods();
	sort_vector< bigint > get_bad_primes() const;

	// Like Pari's pointell
	void get_xy_coords(bigcomplex& x, bigcomplex& y, const bigcomplex& z);

	double get_height_constant();

	//
	// input / output
	//

	void read(std::istream & in);
	void write(std::ostream & out) const;


	void output_torsion(std::ostream & out); // NOT const
	void output_long(std::ostream & out) const;

};




inline std::istream &
operator >> (std::istream & in, elliptic_curve< bigint > & ec)
{
	debug_handler("elliptic_curve< bigint >", "operator >> (std::istream, ec)");

	ec.read(in);
	return in;
}



inline std::ostream &
operator << (std::ostream & out, const elliptic_curve< bigint > & ec)
{
	ec.write(out);
	return out;
}



// The following computes the real local height given the (real)
// x-coordinate.  Used in point_bigint.cc

bigfloat real_height(const elliptic_curve< bigint > &, const bigfloat&);


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

elliptic_curve< bigint > trans_to_curve(complex_periods& cp); // via makec4c6()



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_ELLIPTIC_CURVE_BIGINT_H_GUARD_
