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


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/complex_periods.h"
#include	"LiDIA/elliptic_curves/ec_arith.h"
#include	"LiDIA/elliptic_curves/base_elliptic_curve.h"
#include	"LiDIA/elliptic_curve_bigint.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// periods constructed from a curve, normalized to lattice form
//

// uses the following functions:

// reorders 3 bigcomplex nos so real parts are decreasing
void
complex_periods::reorder1(bigcomplex& a,
			  bigcomplex& b,
			  bigcomplex& c)
{
	if (real(a) < real(c))
		swap(a, c);
	if (real(a) < real(b))
		swap(a, b);
	else if (real(b) < real(c))
		swap(b, c);
}



//reorders 3 bigcomplex nos so e1 is real if any (
void
complex_periods::reorder2(bigcomplex& e1,
			  bigcomplex& e2,
			  bigcomplex& e3)
{
	if (e1.is_approx_real())
		return;
	else if (e2.is_approx_real()) {
		swap(e1, e2);
		return;
	}
	else if (e3.is_approx_real()) {
		swap(e1, e3);
		return;
	}
}



//Computes periods of a curve given the 3 2-division points (i.e. the three
//roots of the cubic)
void
complex_periods::eiperiods(const bigcomplex& e1,
			   const bigcomplex& e2,
			   const bigcomplex& e3,
			   bigcomplex& w1,
			   bigcomplex& w2)
{
	bigcomplex a(sqrt(e3-e1));
	bigcomplex b(sqrt(e3-e2));
	bigcomplex c(sqrt(e2-e1));
	bigcomplex agm1 = cagm(a, b);
	bigcomplex agm2 = cagm(a, c);
	bigfloat pi = Pi();
	w1 = bigcomplex(pi, 0)/agm1;
	w2 = bigcomplex(0, pi)/agm2;
}



// Gets the 3 2-division points given the coefficients
void
complex_periods::getei(const bigcomplex& a1,
		       const bigcomplex& a2,
		       const bigcomplex& a3,
		       const bigcomplex& a4,
		       const bigcomplex& a6,
		       bigcomplex& e1,
		       bigcomplex& e2,
		       bigcomplex& e3)
{
	bigcomplex c1 = a2 + a1*a1/4;
	bigcomplex c2 = a4 + a1*a3/2;
	bigcomplex c3 = a6 + a3*a3/4;
	bigcomplex* ei = solve_cubic(c1, c2, c3);
	e1 = ei[0]; e2 = ei[1]; e3 = ei[2];
	delete [] ei;
}



complex_periods::complex_periods(const base_elliptic_curve< bigint > & E)
{
	bigint a1, a2, a3, a4, a6;
	E.get_ai(a1, a2, a3, a4, a6);

	bigcomplex ca1 = bigfloat(a1);
	bigcomplex ca2 = bigfloat(a2);
	bigcomplex ca3 = bigfloat(a3);
	bigcomplex ca4 = bigfloat(a4);
	bigcomplex ca6 = bigfloat(a6);

	complex_periods::getei(ca1, ca2, ca3, ca4, ca6, e1, e2, e3);

	int allrealroots = e1.is_approx_real() && e2.is_approx_real() && e3.is_approx_real();
	if (allrealroots)
		complex_periods::reorder1(e3, e2, e1);
	else
		complex_periods::reorder2(e3, e2, e1);
// this ordering ensures that eiperiods will give w1,w2 per lattice norm:
// if allrealroots then w1 wil be real, w2 pure imag;
// if not, w1 will be real but w2 is not specially determined.

	complex_periods::eiperiods(e1, e2, e3, w1, w2);

	w1 = bigcomplex(w1.real(), 0); // force to be real as it ought to be anyway
	tau = w2/w1;
	if (tau.real() < -0.1) {
		tau += 1;
		w2 += w1;
	}
	norm_code = 1+allrealroots; //ie type 1 or 2 as appropriate
}



void complex_periods::norm_lattice()
{
	lidia_error_handler("complex_periods", "norm_lattice::not implemented yet");
}



// NB In ec_arith.h/cc's normalize and getc4c6 the normalization is such that
//    tau = w1/w2 in fundamental region; here we have tau=w2/w1 (better for
//    lattice normalization) so must call them with w1, w2 reversed.

void complex_periods::norm_region()
{
	if (norm_code != 3) {
		tau = normalize(w2, w1); // NB reverse params
		norm_code = 3;
		qtau = q(tau);
		if (abs(qtau) > 0.99) {
			std::cout << "Warning from complex_periods::norm_region: qtau = "
				  << qtau << " is not small!\n";
		}
		w1squared = w1*w1;
		w1cubed = w1*w1squared;
		bigcomplex term(bigfloat(1.0)), qtm = qtau;
		sum3 = 0;
		for (bigfloat m = 1; ! term.is_approx_zero(); m += 1) {
			term = qtm*m / (1 - qtm);
			qtm *= qtau;
			sum3 += term;
		}
		sum3 = bigfloat(1.0)/12 - 2*sum3;
	}
}



bigcomplex complex_periods::x_coord(const bigcomplex& qz)
{
	bigcomplex sum(sum3), term(1), qtm(1), w;
	while (! term.is_approx_zero()) {
		w = qtm*qz;
		term = w / power((1 - w), 2);
		qtm *= qtau;
		sum += term;
	}
	term = 1; qtm = qtau;
	while (! term.is_approx_zero()) {
		w = qtm / qz;
		qtm *= qtau;
		term = w / power((1 - w), 2);
		sum += term;
	}
	bigcomplex ans = bigcomplex(0, 2*Pi());
	ans = sum*ans*ans;
	return ans;
}



bigcomplex complex_periods::y_coord(const bigcomplex& qz)
{
	bigcomplex sum(0), term(1), qtm(1), w;
	while (!term.is_approx_zero()) {
		w = qtm*qz;
		term = w*(1 + w) / power((1 - w), 3);
		qtm *= qtau;
		sum += term;
	}
	qtm = qtau; term = 1;
	while (!term.is_approx_zero()) {
		w = qtm / qz;
		term = w*(1 + w) / power((w - 1), 3);
		qtm *= qtau;
		sum += term;
	}
	bigcomplex ans = bigcomplex(0, 2*Pi());
	ans = sum * ans*ans*ans;
	return ans;
}



void complex_periods::xy_coords(bigcomplex& x, bigcomplex& y, const bigcomplex& z)
{
	norm_region();
	bigcomplex z1 = z/w1;
	bigcomplex qz = q(z1);
	while (abs(qz) > 0.9)
		qz *= qtau;
	x = x_coord(qz) / w1squared;
	y = y_coord(qz) / w1cubed;
	return;
}



void make_c4_c6(complex_periods& cp, bigcomplex& cc4, bigcomplex& cc6)
{
	cp.norm_region(); // normed for tau in fundamental region
	bigcomplex w1, w2;
	cp.get_omegas(w1, w2);
	getc4c6(w2, w1, cc4, cc6); // from ec_arith.h
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
