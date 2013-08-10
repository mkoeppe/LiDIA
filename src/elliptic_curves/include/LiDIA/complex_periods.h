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
//	Author	: Nigel Smart (NS), John Cremona (JC)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_COMPLEX_PERIODS_H_GUARD_
#define LIDIA_COMPLEX_PERIODS_H_GUARD_


#ifndef LIDIA_BIGCOMPLEX_H_GUARD_
# include	"LiDIA/bigcomplex.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class T > class base_elliptic_curve;


inline bigcomplex q(const bigcomplex& z) // q(z) = exp(2*pi*i * z)
{
	bigfloat twopix = 2*Pi()*z.real();
	return exp(-2*Pi()*z.imag()) *
		bigcomplex(cos(twopix), sin(twopix));
}


//
// the complex_periods class:
//


class complex_periods {
	bigcomplex w1, w2;
	bigcomplex tau; // whenever meaningful, = w2/w1
	bigcomplex e1, e2, e3; // 2-division values
	bigcomplex qtau, w1squared, w1cubed, sum3;
	int norm_code; // 0 if not set,
	// 1 if normalized for lattice with Re(tau) = 0.5 "type1"
	//   w1 = 2x real, w2 = x+yi, lattice [2x, x+yi]
	// 2 if "   "   "   "    "  "     "    " = 0   "type2"
	//   w1 = x real, w2 = yi pure imag, lattice [x, yi]
	// 3 if normalized for tau in fundmental region
	// -1 Nothing computed yet

	static void reorder1(bigcomplex& a, bigcomplex& b, bigcomplex& c);
	static void reorder2(bigcomplex& e1, bigcomplex& e2, bigcomplex& e3);
	static void eiperiods(const bigcomplex& e1, const bigcomplex& e2, const bigcomplex& e3, bigcomplex& w1, bigcomplex& w2);
	static void getei(const bigcomplex& a1, const bigcomplex& a2, const bigcomplex& a3,
			  const bigcomplex& a4, const bigcomplex& a6, bigcomplex& e1,
			  bigcomplex& e2, bigcomplex& e3);


public:

	// Constructors
	// nb must give two or no args

	complex_periods() : w1(0), w2(0), tau(0), norm_code(-1) {}
	complex_periods(const bigcomplex& ww1, const  bigcomplex& ww2)
		: w1(ww1), w2(ww2), tau(w2/w1), norm_code(0)
	{}

	// normalize now per code
	complex_periods(const bigcomplex& ww1, const bigcomplex& ww2, const int normize)
		: w1(ww1), w2(ww2), norm_code(0)
	{
		if (normize == 1)
			norm_lattice();
		else
			norm_region();
	}

	// for known lattice types
	complex_periods(const bigfloat& x, const bigfloat& y, const int type_set)
		: norm_code(type_set)
	{
		if (type_set == 1) {
			w1 = 2*x;
			w2 = bigcomplex(x, y);
		}
		else {
			w1 = x;
			w2 = bigcomplex(0, y);
		}
		tau = w2/w1;
	}

	// this will normalize for lattice, type = #comp's
	complex_periods(const base_elliptic_curve< bigint > & E);

	// Copying:
	complex_periods(const complex_periods& cp)
		: w1(cp.w1), w2(cp.w2), tau(cp.tau), norm_code(cp.norm_code)
	{}
	complex_periods & operator = (const complex_periods& cp)
	{
		w1 = cp.w1;
		w2 = cp.w2;
		tau = cp.tau;
		norm_code = cp.norm_code;
		return *this;
	}

	// normalizing:
	void norm_lattice(); // ie for w1 real and w2 pure Im or 2w2-w1 pure Im
	// not yet implemented (and maybe never)
	void norm_region(); // ie for Im(t) > 0, -.5< Re(t) <= .5, and |t| >= 1

	// member access functions
	bigcomplex get_tau()
	{
		return tau;
	}
	// up to caller to check norm_code

	int get_norm_code()
	{
		return norm_code;
	}
	//0 unset, 1 lattice, 2 square lattice, 3 region, -1 not initialised

	int get_omegas(bigcomplex& ww1, bigcomplex& ww2) const
	{
		ww1 = w1;
		ww2 = w2;
		return norm_code;
	}
	bigcomplex get_omega_1() const
	{
		return w1;
	}
	bigcomplex get_omega_2() const
	{
		return w2;
	}

	// reset function
	void reset()
	{
		w1 = 0;
		w2 = 0;
		tau = 0;
		norm_code = -1;
	}

	// I/O
	friend std::ostream& operator << (std::ostream& os, const complex_periods& cp);

	// Weierstrass functions:
	bigcomplex x_coord(const bigcomplex& qz); //qz = q(z), z modulo lattice,
	bigcomplex y_coord(const bigcomplex& qz); //gives coords in Y^2 = 4X^3+... model

	void xy_coords(bigcomplex& x, bigcomplex& y, const bigcomplex& z);

}; // end of class complex_periods def'n



inline std::ostream& operator << (std::ostream& os, const complex_periods& cp)
{
	os << "omega_1: " << cp.w1 << "\nomega_2: " << cp.w2
	   << "\ntau: " << cp.tau << "\t";
	switch(cp.norm_code) {
	case 0:  return os << "(normalization: none)";
	case 1:  return os << "(normalization: lattice type 1)";
	case 2:  return os << "(normalization: lattice type 2)";
	default: return os << "(normalization: fundamental region)";
	}
}



void make_c4_c6(complex_periods& cp, bigcomplex& cc4, bigcomplex& cc6);
// NB this normalizes the arg cp; caller must save old copy if required



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_COMPLEX_PERIODS_H_GUARD_
