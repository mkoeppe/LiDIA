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


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigint.h"
#include	"LiDIA/bigrational.h"
#include	"LiDIA/bigfloat.h"
#include	"LiDIA/bigcomplex.h"

// To force the use of point<bigint>!
#include	"LiDIA/point_bigint.h"

#include	"LiDIA/elliptic_curves/base_elliptic_curve.h"
#include	"LiDIA/elliptic_curves/elliptic_curve_rep_bigint.h"
#include	"LiDIA/elliptic_curve_bigint.h"
#include	"LiDIA/elliptic_curves/ec_arith.h"
#include	"LiDIA/reduction_type.h"
#include	"LiDIA/quartic.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// constructor / destructor
//

elliptic_curve< bigint >::elliptic_curve()
	: base_elliptic_curve< bigint > ()
{
	debug_handler("elliptic_curve< bigint >",
		      "elliptic_curve()");
}



elliptic_curve< bigint >::elliptic_curve(const elliptic_curve< bigint > & E)
	: base_elliptic_curve< bigint > (E.e)
{
	debug_handler("elliptic_curve< bigint >",
		      "elliptic_curve(const elliptic_curve< bigint > &)");
}



// NB The following constructor constructs a curve from given c4, c6
// invariants, first checking that they are valid.  This must NOT be
// confused with the two-parameter base_elliptic_curve constructor which sets
// a4 and a6, with a1=a2=a3=0.  This one is more useful for curves/Z
// where we need to be able to construct them via periods and c4, c6.

elliptic_curve< bigint >::elliptic_curve(const bigint & cc4, const bigint & cc6,
					 elliptic_curve_flags::curve_model m)
	: base_elliptic_curve< bigint > (cc4, cc6, m)
{
	debug_handler("elliptic_curve< bigint >",
		      "elliptic_curve(const bigint&x4, "
		      "const bigint&x6, "
		      "elliptic_curve_flags::curve_model m)");
}



elliptic_curve< bigint >::elliptic_curve(const bigint & x1, const bigint & x2,
					 const bigint & x3, const bigint & x4,
					 const bigint & x6,
					 elliptic_curve_flags::curve_model m)
	: base_elliptic_curve< bigint > (x1, x2, x3, x4, x6, m)
{
	debug_handler("elliptic_curve< bigint >",
		      "elliptic_curve(x1, x2, x3, x4, x6, "
		      "elliptic_curve_flags::curve_model m)");
}



elliptic_curve< bigint >::~elliptic_curve()
{
	debug_handler("elliptic_curve< bigint >",
		      "~elliptic_curve()");

	// Torsion points in (*e) may have references
	// to this curve. Dead lock.
	//
	if (e != NULL) {
		// If lock is set, it is guaranteed, that the ref counter
		// is strictly larger than 1. So, just decrease it, which
		// is done by ~base_elliptic_curve.
		//
		if (e->lock == false) {
			// Otherwise, check whether there are only internal references
			// to the curve. If so remove the curve.
			//
			if (e->torsion.get_size() + 1 == e->get_ref_counter()) {
				// To prevent further calls of this destructor from
				// also checking this condition, set the lock.
				//
				e->lock = true;
				e->torsion.reset();

				// Now, ref counter must be 1.
				//
				if (e->get_ref_counter() != 1)
					lidia_error_handler("elliptic_curve< bigint >::~elliptic_curve()",
							    "Reference counter not 1.");
			}
		}
	}

	// Cleaning done by ~base_elliptic_curve<T>.
}



//
// assignment
//
void
elliptic_curve< bigint >::assign (const elliptic_curve< bigint > & E)
{
	debug_handler ("elliptic_curve< bigint >",
		       "assign(const elliptic_curve< bigint >");

	//std::cout << "elliptic_curve<bigint>::assign(const elliptic_curve<bigint>" << std::endl;

	if (this != &E) {
		if (E.e == NULL)
			lidia_error_handler("elliptic_curve< T >::"
					    "assign(const elliptic_curve< bigint > &)",
					    "e == NULL");
		else if (e == NULL) {
			e = E.e;
			e->inc_ref_counter();
			//std::cout << "elliptic_curve_bigint::assign::inc " << e->get_ref_counter() << std::endl;
		}
		else if (e != E.e) {
			if (e->get_ref_counter() == 1)
				delete e;
			else {
				e->dec_ref_counter();
				//std::cout << "elliptic_curve_bigint::assign::dec " << e->get_ref_counter() << std::endl;
			}

			e = E.e;
			e->inc_ref_counter();
			//std::cout << "elliptic_curve_bigint::assign::inc " << e->get_ref_counter() << std::endl;
		}
	}
}



elliptic_curve< bigint > &
elliptic_curve< bigint >::operator = (const elliptic_curve< bigint > & E)
{
	debug_handler ("elliptic_curve< bigint >",
		       "operator = (const elliptic_curve< bigint >");

	this->assign(E);
	return *this;
}



void
elliptic_curve< bigint >::set_coefficients(const bigint & x4, const bigint & x6,
					   elliptic_curve_flags::curve_model m)
{
	debug_handler ("elliptic_curve< bigint >",
		       "set_coefficients(x4, x6, m)");

	if (e == NULL)
		e = new elliptic_curve_rep< bigint > (x4, x6, m);

	else if (e->get_ref_counter() != 1) {
		lidia_error_handler("elliptic_curve< bigint >::set_coefficients(x4, x6, m)",
				    "Cannot change curve with more than 1 reference.");
	}
	else {
		e->set_coefficients(x4, x6, m);
		e->init(x4, x6);
	}
}



void
elliptic_curve< bigint >::set_coefficients(const bigint & x1, const bigint & x2,
					   const bigint & x3, const bigint & x4, const bigint & x6,
					   elliptic_curve_flags::curve_model m)
{
	debug_handler ("elliptic_curve< bigint >",
		       "set_coefficients(x1, x2, x3, x4, x6, m)");

	if (e == NULL)
		e = new elliptic_curve_rep< bigint > (x1, x2, x3, x4, x6, m);

	else if (e->get_ref_counter() != 1) {
		lidia_error_handler("elliptic_curve< bigint >::set_coefficients(x1, x2, x3, x4, x6, m)",
				    "Cannot change curve with more than 1 reference.");
	}
	else {
		e->set_coefficients(x1, x2, x3, x4, x6, m);
		e->init(x1, x2, x3, x4, x6);
	}
}



//
// Computation of torsion
//

//After Darrin Doud as adapted by JC, Converted to LiDIA by NPS
point< bigint >
make_torsion_point(elliptic_curve< bigint > E,
		   const bigcomplex& z)
{
	bigcomplex cx, cy;
	E.get_xy_coords(cx, cy, z);
	bigint x, y;
	cx.real().bigintify(x);
	cy.real().bigintify(y);
	point< bigint > P(x, y, E);
	return P;
}



void
elliptic_curve< bigint >::make_torsion()
{
	base_vector< point < bigint > > torsion;
	torsion.set_capacity(16); // Upper bound on ntorsion
	torsion.set_size(0);

//
// table[i][] contains a list of possible maximal orders for a point,
// given that the 2-torsion subgroup has order i
// N.B. These are partially ordered by divisiblity
//
// MM, added -1 for AIX xlC.
	static long table[5][5] = {{-1}, {5, 7, 9, 3}, {12, 6, 8, 4, 10}, {-1}, {8, 6, 4}};
	static long nt[5] = {0, 4, 5, 0, 3};

	bigint sa2, sa4, sa6, d, x, y;

	// local copies, MM
	bigint A1(get_a1()), A2(get_a2()), A3(get_a3()), A4(get_a4()), A6(get_a6());

	bigfloat ra1 = bigfloat(A1), ra2 = bigfloat(A2), ra3 = bigfloat(A3);
	long i, j, nroots, ntp, nt2; int scaled_flag;

	base_vector< point < bigint > > two_torsion(4, 0);
	base_vector< point < bigint > > cycle(8, 0);

	// The next function will initialize periods if not done so already
	bigcomplex w1 = get_omega_1(), w2 = get_omega_2();
	bigcomplex z, z2 = w2/2;
	int found;
	e->periods.norm_region();

	ntp = 1; nt2 = 1;
	torsion[0] = point< bigint > (*this); // zero point
	two_torsion[0] = point< bigint > (*this); // zero point

	if (A1.is_zero() && A3.is_zero()) {
		sa2 = A2;
		sa4 = A4;
		sa6 = A6;
		scaled_flag = 0;
	}
	else {
		sa2 = A1*A1 + 4*A2;
		sa4 = 8*A1*A3 + 16*A4;
		sa6 = 16*A3*A3 + 64*A6;
		scaled_flag = 1;
	}
	d = sa2*sa2*(sa4*sa4 - 4*sa2*sa6) + 18*sa2*sa4*sa6
		- 4*sa4*sa4*sa4 - 27*sa6*sa6;
	point< bigint > p(*this), q(*this);


// First test y=0 for points of order 2:

	base_vector< bigint > xlist = int_roots_cubic(sa2, sa4, sa6);
	nroots = xlist.size();
	for (i = 0; i < nroots; i++) {
		x = xlist[i];
		if (scaled_flag)
			p.assign(2*x, - A1*x - 4*A3, bigint(8));
		else
			p.assign(x, 0);

		two_torsion[nt2++] = p;
	}

	for (i = 0, found = 0; (i < nt[nt2]) && (!found); i++) {
		long ni = table[nt2][i]; z = w1/ni;
		p = make_torsion_point(*this, z);

		found = (p.on_curve()) && (order(p, cycle) == ni);
		if (!found && (e->conncomp == 2) && ((ni%2) == 0)) {
			p = make_torsion_point(*this, z+z2);

			found = (p.on_curve() && (order(p, cycle) == ni));
			if (!found && ((ni%4) == 2)) {
				p = make_torsion_point(*this, 2*z+z2);
				found = (p.on_curve()) && (order(p, cycle) == ni);
			}
		}
		if (found) {
			ntp = ni;
			for (j = 1; j < ni; j++) {
				torsion[j] = cycle[j];
			}
		}
	}

	if (ntp == 1) {
		// C2xC2
		e->torsion = two_torsion;
		return;
	}

	if (nt2 == 4) {
		// non-cyclic case, C2xC4, C2xC6 or C2xC8
		// Find a point of order 2 not in the second factor and add it in
		for (i = 1; i < 4; i++) {
			p = two_torsion[i];
			int fl = 0;
			for (j = 0; j < ntp && !fl; j++)
			{ if (torsion[j] == p) { fl = 1; } }
			if (!fl) {
				for (j = 0; j < ntp; j++) {
					torsion[ntp+j] = torsion[j]+p;
				}
				ntp *= 2;
				break;
			}
		}
	}

	e->torsion = torsion;
}



#if 0
// Old Torsion Routine
void
elliptic_curve< bigint >::make_torsion()
{
	bigint sa2, sa4, sa6, d, x, y;
	long i, nroots;
	int scaled_flag;
	torsion.set_capacity(16); // Upper bound on ntorsion
	torsion.set_size(0);
	long ntp = 1;
	torsion[0] = point< bigint > (this); // zero point
	if (a1.is_zero() && a3.is_zero()) {
		sa2 = a2;
		sa4 = a4;
		sa6 = a6;
		scaled_flag = 0;
	}
	else {
		sa2 = a1*a1 + 4*a2;
		sa4 = 8*a1*a3 + 16*a4;
		sa6 = 16*a3*a3 + 64*a6;
		scaled_flag = 1;
	}
	d = sa2*sa2*(sa4*sa4 - 4*sa2*sa6) + 18*sa2*sa4*sa6
		- 4*sa4*sa4*sa4 - 27*sa6*sa6;
	point< bigint > p(this);

	// First test y=0 for points of order 2:

	base_vector< bigint > xlist = int_roots_cubic(sa2, sa4, sa6);
	nroots = xlist.size();
	for (i = 0; i < nroots; i++) {
		x = xlist[i];
		if (scaled_flag)
			p.assign(2*x, - a1*x - 4*a3, 8);
		else
			p.assign(x, 0);
		torsion[ntp++] = p;
	}

	// Now test y such that y^2 divides d:
	base_vector< bigint > possible_y = square_divs(d);
	long j, nj = possible_y.size();
	for (j = 0; j < nj; j++) {
		y = possible_y[j];
		xlist = int_roots_cubic(sa2, sa4, sa6-y*y);
		nroots = xlist.size();
		for (i = 0; i < nroots; i++) {
			x = xlist[i];
			if (scaled_flag)
				p.assign(2*x, y - a1*x - 4*a3, 8);
			else
				p.assign(x, y);
			if (p.get_order() > 0) {
				torsion[ntp++] = p;
				torsion[ntp++] = -p; // N.B. order > 2 here!
			}
		}
	}
}
#endif




//
// Access
//

int
elliptic_curve< bigint >::get_ord_p_discriminant(const bigint& p) const
{
	debug_handler("elliptic_curve< bigint >", "get_ord_p_discriminant(p) const");

	if (e == NULL)
		lidia_error_handler("elliptic_curve< bigint >::"
				    "get_ord_p_discriminant(p) const",
				    "e == NULL");
	int index;
	if (!e->bad_primes.bin_search(p, index)) {
		return 0;
	}
	return (e->reduct_array)[index].get_ord_p_discriminant();
}



int
elliptic_curve< bigint >::get_ord_p_N(const bigint& p) const
{
	debug_handler("elliptic_curve< bigint >", "get_ord_p_N(p) const");

	if (e == NULL)
		lidia_error_handler("elliptic_curve< bigint >::"
				    "get_ord_p_N(p) const",
				    "e == NULL");
	int index;
	if (!e->bad_primes.bin_search(p, index)) {
		return 0;
	}
	return (e->reduct_array)[index].get_ord_p_N();
}



int
elliptic_curve< bigint >::get_ord_p_j_denom(const bigint& p) const
{
	debug_handler("elliptic_curve< bigint >", "get_ord_p_j_denom(p) const");

	if (e == NULL)
		lidia_error_handler("elliptic_curve< bigint >::"
				    "get_ord_p_j_denom(p) const",
				    "e == NULL");
	int index;
	if (!e->bad_primes.bin_search(p, index)) {
		return 0;
	}
	return (e->reduct_array)[index].get_ord_p_j_denom();
}



int
elliptic_curve< bigint >::get_c_p(const bigint& p) const
{
	debug_handler("elliptic_curve< bigint >", "get_c_p(p) const");

	if (e == NULL)
		lidia_error_handler("elliptic_curve< bigint >::"
				    "get_c_p(p) const",
				    "e == NULL");
	int index;
	if (!e->bad_primes.bin_search(p, index)) {
		return 1;
	}
	return (e->reduct_array)[index].get_c_p();
}



int
elliptic_curve< bigint >::get_product_cp() const
{
	debug_handler("elliptic_curve< bigint >", "get_product_cp() const");

	if (e == NULL)
		lidia_error_handler("elliptic_curve< bigint >::"
				    "get_product_cp() const",
				    "e == NULL");

	int i, ans = 1, np = e->bad_primes.get_size();
	for (i = 0; i < np; i++) {
		ans *= (e->reduct_array)[i].get_c_p();
	}
	return ans;
}



Kodaira_code
elliptic_curve< bigint >::get_Kodaira_code(const bigint& p) const
{
	debug_handler("elliptic_curve< bigint >", "get_Kodaira_code(p) const");

	if (e == NULL)
		lidia_error_handler("elliptic_curve< bigint >::"
				    "get_Kodaira_code(p) const",
				    "e == NULL");
	int index;
	if (!e->bad_primes.bin_search(p, index)) {
		return Kodaira_code(0);
	}
	return (e->reduct_array)[index].get_code();
}



bigrational
elliptic_curve< bigint >::j_invariant() const
{
	debug_handler("elliptic_curve< bigint >", "j_invariant() const");

	if (e == NULL)
		lidia_error_handler("elliptic_curve< bigint >::"
				    "j_invariant() const",
				    "e == NULL");

	return bigrational(e->c4*e->c4*e->c4, e->delta);
}



bigint
elliptic_curve< bigint >::get_conductor() const
{
	debug_handler("elliptic_curve< bigint >", "get_conductor() const");

	if (e == NULL)
		lidia_error_handler("elliptic_curve< bigint >::"
				    "get_conductor() const",
				    "e == NULL");
	return e->N;
}



bigcomplex
elliptic_curve< bigint >::get_omega_1()
{
	debug_handler("elliptic_curve< bigint >", "get_omega1() const");

	if (e == NULL)
		lidia_error_handler("elliptic_curve< bigint >::"
				    "get_omega1() const",
				    "e == NULL");

	if (e->periods.get_norm_code() == -1) {
		e->periods = complex_periods(*this);
	}
	return e->periods.get_omega_1();
}



bigcomplex
elliptic_curve< bigint >::get_omega_2()
{
	debug_handler("elliptic_curve< bigint >", "get_omega2() const");

	if (e == NULL)
		lidia_error_handler("elliptic_curve< bigint >::"
				    "get_omega2() const",
				    "e == NULL");

	if (e->periods.get_norm_code() == -1) {
		e->periods = complex_periods(*this);
	}
	return e->periods.get_omega_2();
}



complex_periods
elliptic_curve< bigint >::get_periods()
{
	debug_handler("elliptic_curve< bigint >", "get_periods() const");

	if (e == NULL)
		lidia_error_handler("elliptic_curve< bigint >::"
				    "get_periods() const",
				    "e == NULL");

	if (e->periods.get_norm_code() == -1) {
		e->periods = complex_periods(*this);
	}
	return e->periods;
}



void
elliptic_curve< bigint >::get_xy_coords(bigcomplex& x, bigcomplex& y, const bigcomplex& z)
{
	debug_handler("elliptic_curve< bigint >", "get_xy_coords(x, y, z) const");

	if (e == NULL)
		lidia_error_handler("elliptic_curve< bigint >::"
				    "get_xy_coords(x, y, z) const",
				    "e == NULL");

	if (e->periods.get_norm_code() == -1) {
		e->periods = complex_periods(*this);
	}
	e->periods.xy_coords(x, y, z);
	bigfloat ra1(e->a1), ra2(e->a2), ra3(e->a3);
	x = x-(ra1*ra1+4*ra2)/12;
	y = (y - ra1*x - ra3)/2;
}



int
elliptic_curve< bigint >::get_no_tors()
{
	debug_handler("elliptic_curve< bigint >", "get_no_tors() const");

	if (e == NULL)
		lidia_error_handler("elliptic_curve< bigint >::"
				    "get_no_tors()const",
				    "e == NULL");

	if (e->torsion.size() == 0)
		make_torsion();

	return e->torsion.size();

}



base_vector< point < bigint > >
elliptic_curve< bigint >::get_torsion()
{
	debug_handler("elliptic_curve< bigint >", "get_torsion() const");

	if (e == NULL)
		lidia_error_handler("elliptic_curve< bigint >::"
				    "get_torsion()const",
				    "e == NULL");

	if (e->torsion.size() == 0)
		make_torsion();

	return e->torsion;
}



sort_vector< bigint >
elliptic_curve< bigint >::get_bad_primes() const
{
	debug_handler("elliptic_curve< bigint >", "get_bad_primes() const");

	if (e == NULL)
		lidia_error_handler("elliptic_curve< bigint >::"
				    "get_bad_primes()const",
				    "e == NULL");
	return e->bad_primes;
}



// Functions to compute the height contant
// Silvermans Constant Implemented By John Cremona
// Sikseks    Constant Implemented By Nigel Smart
// Both LiDIA'ified by NS


// ----------------------------------------
// Silvermans Constant Code
// ----------------------------------------

double
logplus(double x)
{
	double ax = std::fabs(x);
	if (ax < 1)
		return 0;
	return std::log(ax);
}



double
hj(const elliptic_curve< bigint > & CD, double& realjay)
{
	bigint c4, c6, njay, djay;
	c4 = CD.get_c4();
	c6 = CD.get_c6();
	power(njay, c4, 3);
	djay = CD.discriminant();
	if (djay.is_zero() || njay.is_zero()) {
		realjay = 0;
		return 0;
	}

	double g = dbl(gcd(njay, djay));
	double xnjay = dbl(njay)/g;
	double xdjay = dbl(djay)/g;

	realjay = xnjay/xdjay;

	double x = std::log(std::fabs(xnjay));
	double y = std::log(std::fabs(xdjay));

	if (x < y)
		return y;
	else
		return x;
}



double
silverman(elliptic_curve< bigint > & E)
{
	double realjay;
	double hjay = hj(E, realjay);
	bigint delta = E.discriminant();
	bigint b2 = E.get_b2();

// NB the constant 1.922 = 2*0.961 below is from Bremner's correction
// to Silverman's paper; Silverman has 0.973 givin 2*0.973 = 1.946.

	double mu = 1.922  + hjay/12
		+ std::log(std::fabs(dbl(delta)))/6
		+ logplus(realjay)/6
		+ logplus(dbl(b2)/12);

	if (!b2.is_zero())
		mu += std::log(2.0);

	return mu;
}



// ----------------------------------------
// Sikseks Constant Code
// ----------------------------------------

int
is_in_int(const bigfloat rr, const bigfloat i1l, const bigfloat i1r,
	  const bigfloat i2l, const bigfloat i2r, const long numint)
{
	if (numint > 0) {
		if (rr <= i1r && rr >= i1l) {
			return 1;
		}
		if (numint == 2 && rr <= i2r && rr >= i2l) {
			return 1;
		}
	}
	return 0;
}



int is_in_int2(const bigfloat rr, bigfloat** intervals, const long numint)
{
	long i;
	for (i = 0; i < numint; i++) {
		if ((rr >= intervals[i][0]) && (rr <= intervals[i][1])) {
			return 1;
		}
	}
	return 0;
}



// This procedure calculates dv' for the infinite prime
bigfloat
calc_dvd_inf(const elliptic_curve< bigint > & EE)
{
	bigfloat b2, b4, b6, b8, rr, rx, rn;

	b2 = bigfloat(EE.get_b2());
	b4 = bigfloat(EE.get_b4());
	b6 = bigfloat(EE.get_b6());
	b8 = bigfloat(EE.get_b8());
	bigcomplex *rt;
	rt = solve_cubic(b2/4.0, b4/2.0, b6/4.0);
	bigint del = EE.discriminant();
	long i, j, numrt;
	bigfloat rrt[4];
	rrt[0] = 0.0;
	if (del < 0) {
		for (i = 0; i < 3; i++) {
			if (rt[i].imag().is_approx_zero()) {
				rr = rt[i].real();
				i = 3;
			}
		}
		if (rr.is_approx_zero()) {
			numrt = 1;
		}
		else {
			rrt[1] = 1.0/rr;
			numrt = 2;
		}
	}
	else {
		numrt = 1;
		for (i = 0; i < 3; i++) {
			if (!rt[i].real().is_approx_zero()) {
				rrt[numrt] = 1.0/rt[i].real();
				numrt++;
			}
		}
	}
	delete[] rt;

	// Sort the real roots
	long fl = 0;
	while (fl == 0) {
		fl = 1;
		for (i = 0; i < numrt-1; i++) {
			if (rrt[i] > rrt[i+1]) {
				rr = rrt[i+1];
				rrt[i+1] = rrt[i];
				rrt[i] = rr;
				fl = 0;
			}
		}
	}

	// Make Intervals
	bigfloat **intervals = new bigfloat*[3];
	for (i = 0; i < 3; i++) {
		intervals[i] = new bigfloat[2];
	}

	bigfloat rhs, lhs;
	long numint = 0;
	if (b6 < 0.0) {
		for (i = 0; i < numrt; i = i+2) {
			intervals[numint][1] = rrt[numrt-i-1];
			intervals[numint][0] = rrt[numrt-i-2];
			numint = numint+1;
		}
	}
	else {
		intervals[numint][1] = rrt[numrt-1]+2;
		intervals[numint][0] = rrt[numrt-1];
		numint = numint+1;
		for (i = 1; i < numrt; i = i+2) {
			intervals[numint][1] = rrt[numrt-i-1];
			if (numrt-i-2 > 0) {
				intervals[numint][0] = rrt[numrt-i-2];
			}
			else {
				intervals[numint][0] = -1.0;
			}
			numint = numint+1;
		}
	}
	for (i = 0; i < numint; i++) {
		if (intervals[i][1] > 1.0) {
			intervals[i][1] = 1.0;
		}
		if (intervals[i][0] < -1.0) {
			intervals[i][0] = -1.0;
		}
	}
	i = 0;
	while (i < numint) {
		fl = 0;
		while (fl == 0 && i < numint) {
			fl = 1;
			if (intervals[i][0] > intervals[i][1]) {
				fl = 0;
				for (j = i; j < numint-1; j++) {
					intervals[j][0] = intervals[j+1][0];
					intervals[j][1] = intervals[j+1][1];
				}
				numint = numint-1;
			}
		}
		i = i+1;
	}

	bigfloat *te = new bigfloat[36];
	long tec = 2*numint;
	for (i = 0; i < numint; i++) {
		te[2*i] = intervals[i][0];
		te[2*i+1] = intervals[i][1];
	}

	if (numint != 0) {
		// Roots of g
		rt = solve_quartic(bigfloat(0.0), -b4, -2.0*b6, -b8);
		for (i = 0; i < 4; i++) {
			if (!rt[i].is_approx_zero()) {
				rt[i] = 1/rt[i];
				if (rt[i].imag().is_approx_zero()) {
					rr = rt[i].real();
					if (is_in_int2(rr, intervals, numint)) {
						te[tec] = rr;
						tec = tec+1;
					}
				}
			}
		}
		delete[] rt;
		// Roots of f
		rr = 0.0;
		if (is_in_int2(rr, intervals, numint)) {
			te[tec] = rr;
			tec = tec+1;
		}
		rt = solve_cubic(b2/4.0, b4/2.0, b6/4.0);
		for (i = 0; i < 3; i++) {
			if (!rt[i].is_approx_zero()) {
				rt[i] = 1/rt[i];
				if (rt[i].imag().is_approx_zero()) {
					rr = rt[i].real();
					if (is_in_int2(rr, intervals, numint)) {
						te[tec] = rr;
						tec = tec+1;
					}
				}
			}
		}
		delete[] rt;
		// Roots of f+g
		rt = solve_quartic(bigfloat(4.0), b2-b4, 2.0*(b4-b6), b6-b8);
		for (i = 0; i < 4; i++) {
			if (!rt[i].is_approx_zero()) {
				rt[i] = 1/rt[i];
				if (rt[i].imag().is_approx_zero()) {
					rr = rt[i].real();
					if (is_in_int2(rr, intervals, numint)) {
						te[tec] = rr;
						tec = tec+1;
					}
				}
			}
		}
		delete[] rt;
		// Roots of f-g
		rt = solve_quartic(bigfloat(-4.0), -b2-b4, 2.0*(-b4-b6), -b6-b8);
		for (i = 0; i < 4; i++) {
			if (!rt[i].is_approx_zero()) {
				rt[i] = 1/rt[i];
				if (rt[i].imag().is_approx_zero()) {
					rr = rt[i].real();
					if (is_in_int2(rr, intervals, numint)) {
						te[tec] = rr;
						tec = tec+1;
					}
				}
			}
		}
		delete[] rt;
		// Roots of f'
		rt = solve_cubic(b2/2.0, 3.0*b4/2.0, b6);
		for (i = 0; i < 3; i++) {
			if (!rt[i].is_approx_zero()) {
				if (rt[i].imag().is_approx_zero()) {
					rr = rt[i].real();
					if (is_in_int2(rr, intervals, numint)) {
						te[tec] = rr;
						tec = tec+1;
					}
				}
			}
		}
		delete[] rt;
		// Roots of g'
		rr = 0.0;
		if (is_in_int2(rr, intervals, numint)) {
			te[tec] = rr;
			tec = tec+1;
		}
		bigfloat dr = 9.0*b6*b6-8.0*b8*b4;
		if (dr >= 0.0) {
			dr = sqrt(dr);
			rr = (3.0*b6+dr)/(-4.0*b8);
			if (is_in_int2(rr, intervals, numint)) {
				te[tec] = rr;
				tec = tec+1;
			}
			rr = (3.0*b6-dr)/(-4.0*b8);
			if (is_in_int2(rr, intervals, numint)) {
				te[tec] = rr;
				tec = tec+1;
			}
		}
	}

	bigfloat dvd = 100000000.0;
	bigfloat f, g;
	for (i = 0; i < tec; i++) {
		rr = te[i];
		f = (((b6*rr+2.0*b4)*rr+b2)*rr+4.0)*rr;
		g = (((-b8*rr-2.0*b6)*rr-b4)*rr)*rr+1.0;
		f = abs(f);
		g = abs(g);
		rn = f;
		if (g > f) {
			rn = g;
		}
		if (i == 0) {
			dvd = rn;
		}
		else if (dvd > rn) {
			dvd = rn;
		}
	}
	delete[] te;
	for (i = 0; i < 3; i++) {
		delete[] intervals[i];
	}
	delete[] intervals;

	return dvd;

}



// This procedure calculates dv for the infinite prime
bigfloat
calc_dv_inf(const elliptic_curve< bigint > & EE)
{
	bigfloat b2, b4, b6, b8, rr, rx, rn;

	b2 = bigfloat(EE.get_b2());
	b4 = bigfloat(EE.get_b4());
	b6 = bigfloat(EE.get_b6());
	b8 = bigfloat(EE.get_b8());
	bigcomplex *rt;
	rt = solve_cubic(b2/4.0, b4/2.0, b6/4.0);
	bigint del = EE.discriminant();
	bigfloat i1l, i1u, i2l, i2u; // Bounds on Intervals
	bigfloat r1l, r1r;
	long i, numint;
	if (del < 0) {
		for (i = 0; i < 3; i++) {
			if (rt[i].imag().is_approx_zero()) {
				rr = rt[i].real();
				i = 3;
			}
		}
		numint = 0;
		if (rr <= 1.0); {
			numint = 1;
			r1l = rr;
			r1r = 1.0;
			if (rr < -1.0) {
				r1l = -1.0;
			}
		}
	}
	else {
		rn = rt[0].real();
		rx = rt[1].real();
		// Want rn<rr<rx
		if (rn > rx) {
			rr = rx;
			rx = rn;
			rn = rr;
		}
		rr = rt[2].real();
		if (rr < rn) {
			rr = rn;
			rn = rt[2].real();
		}
		else if (rx < rr) {
			rr = rx;
			rx = rt[2].real();
		}

		numint = 2;
		i1l = rn;
		i1u = rr;
		i2l = rx;
		i2u = 1.0;
		// Deal With Right Most Interval
		if (i2l > i2u) {
			numint = 1;
		}
		if (i2l < -1.0) {
			numint = 1;
			i1l = -1.0;
			i1u = i2u;
		}
		// Now Deal With Left Most Interval
		if (i1u > 1.0) {
			i1u = 1.0;
		}
		if (i1l < -1.0) {
			i1l = -1.0;
		}
		if (i1l > i1u) {
			numint = numint-1;
			i1l = -1.0;
			i1u = i2u;
		}
	}
	delete[] rt;

	bigfloat* te = new bigfloat[36];
	long tec = 2*numint;
	if (numint > 0) {
		te[0] = i1l;
		te[1] = i1u;
		if (numint == 2) {
			te[2] = i2l;
			te[3] = i2u;
		}
	}

	if (numint != 0) {
		// Roots of g
		rt = solve_quartic(bigfloat(0.0), -b4, -2.0*b6, -b8);
		for (i = 0; i < 4; i++) {
			if (rt[i].imag().is_approx_zero()) {
				rr = rt[i].real();
				if (is_in_int(rr, i1l, i1u, i2l, i2u, numint)) {
					te[tec] = rr;
					tec = tec+1;
				}
			}
		}
		delete[] rt;
		// Roots of f
		rt = solve_cubic(b2/4.0, b4/2.0, b6/4.0);
		for (i = 0; i < 3; i++) {
			if (rt[i].imag().is_approx_zero()) {
				rr = rt[i].real();
				if (is_in_int(rr, i1l, i1u, i2l, i2u, numint)) {
					te[tec] = rr;
					tec = tec+1;
				}
			}
		}
		delete[] rt;
		// Roots of f+g
		rt = solve_quartic(bigfloat(4.0), b2-b4, 2.0*(b4-b6), b6-b8);
		for (i = 0; i < 4; i++) {
			if (rt[i].imag().is_approx_zero()) {
				rr = rt[i].real();
				if (is_in_int(rr, i1l, i1u, i2l, i2u, numint)) {
					te[tec] = rr;
					tec = tec+1;
				}
			}
		}
		delete[] rt;
		// Roots of f-g
		rt = solve_quartic(bigfloat(4.0), b2+b4, 2.0*(b4+b6), b6+b8);
		for (i = 0; i < 4; i++) {
			if (rt[i].imag().is_approx_zero()) {
				rr = rt[i].real();
				if (is_in_int(rr, i1l, i1u, i2l, i2u, numint)) {
					te[tec] = rr;
					tec = tec+1;
				}
			}
		}
		delete[] rt;
		// Roots of f'
		if (-96.0*b4+4.0*b2*b2 >= 0.0) {
			rn = sqrt(-96.0*b4+4*b2*b2);
			rr = (-2*b2+rn)/24.0;
			if (is_in_int(rr, i1l, i1u, i2l, i2u, numint)) {
				te[tec] = rr;
				tec = tec+1;
			}
			rr = (-2*b2-rn)/24.0;
			if (is_in_int(rr, i1l, i1u, i2l, i2u, numint)) {
				te[tec] = rr;
				tec = tec+1;
			}
		}
		// Roots of g'
		rt = solve_cubic(bigfloat(0.0), -b4/2.0, -b6/2.0);
		for (i = 0; i < 3; i++) {
			if (rt[i].imag().is_approx_zero()) {
				rr = rt[i].real();
				if (is_in_int(rr, i1l, i1u, i2l, i2u, numint)) {
					te[tec] = rr;
					tec = tec+1;
				}
			}
		}
		delete[] rt;
	}

	bigfloat dv = 100000000.0;
	for (i = 0; i < tec; i++) {
		rr = te[i]*te[i];
		rn = 4.0*rr*te[i]+b2*rr+2.0*b4*te[i]+b6;
		rn = abs(rn);
		rx = rr*rr-b4*rr-2.0*b6*te[i]-b8;
		rx = abs(rx);
		if (rn > rx) {
			rx = rn;
		}
		if (i == 0) {
			dv = rx;
		}
		else if (dv > rx) {
			dv = rx;
		}
	}
	delete[] te;

	return dv;

}



long
liftu(const elliptic_curve< bigint > & EE, const bigint p,
      const bigint X, const bigint mod, const long i, const long mx)
{
	bigint b2, b4, b6, b8, y, fx, gx, nmod;
	EE.get_bi(b2, b4, b6, b8);

	long nmax;

	nmod = mod*p*p;
	nmax = mx;
	for (y = X; y < nmod; y = y+mod) {
		fx = (((4*y+b2)*y+2*b4)*y+b6) % nmod;
		gx = (((y*y-b4)*y-2*b6)*y-b8) % (mod*p);
		if (fx == 0 && gx == 0) {
			// Square Condition
			quartic q(0, 4, b2, 2*b4, b6);
			if (q.zp_soluble(p, y, 2*i)) {
				if (i > nmax) {
					nmax = i;
				}
				nmax = liftu(EE, p, y, nmod, i+1, nmax);
			}
		}
	}

	return nmax;
}



long
liftv(const elliptic_curve< bigint > & EE, const bigint p, const bigint X, const bigint mod,
      const long i, const long mx)
{
	bigint b2, b4, b6, b8, y, fx, gx, nmod;
	EE.get_bi(b2, b4, b6, b8);

	long nmax;

	nmod = mod*p*p;
	nmax = mx;
	for (y = X; y < nmod; y = y+mod) {
		fx = (((4*y+b2)*y+2*b4)*y+b6) % nmod;
		gx = (((y*y-b4)*y-2*b6)*y-b8) % nmod;
		if (fx == 0 && gx == 0) {
			// Square Condition
			quartic q(0, 4, b2, 2*b4, b6);
			if (q.zp_soluble(p, y, 2*i)) {
				if (i > nmax) {
					nmax = i;
				}
				nmax = liftv(EE, p, y, nmod, i+1, nmax);
			}
		}
	}

	return nmax;
}



// Finite Prime
bigfloat
calc_ln_ev_fin(const elliptic_curve< bigint > & EE, const bigint p)
{
	// Cope with p "good"
	bigint disc = EE.discriminant();
	if (disc%p != 0) {
		return 0.0;
	}

	long ui = liftu(EE, p, 0, 1, 1, 0)+1;
	if (ui == 1) {
		return 0.0;
	}
	bigfloat ee;
	if (p == 2) {
		long vi = liftv(EE, p, 0, 1, 1, 0)+1;
		if (vi < ui) {
			ee = -(2*vi-1);
		}
		else {
			ee = -2*(ui-1);
		}
	}
	else {
		ee = -(ui-1)*2;
	}

	return -ee*log(bigfloat(p));
}



// The curve EE must be minimal
double
siksek(const elliptic_curve< bigint > & EE)
{
	bigfloat dv, dvd, htc, ev;

	// Arch Part
	dv = calc_dv_inf(EE);
	//std::cout << "dv=" << dv << std::endl;
	dvd = calc_dvd_inf(EE);
	//std::cout << "dvd=" << dv << std::endl;
	htc = dvd;
	if (dv < dvd) {
		htc = dv;
	}
	htc = -log(htc)/3;

	sort_vector< bigint > bad_primes = EE.get_bad_primes();

	long cp, np = bad_primes.size();
	for (long i = 0; i < np; i++) {
		bigint p = bad_primes[i];
		cp = EE.get_c_p(p);
		//std::cout << p << "  " << cp;
		if (cp != 1) {
			ev = calc_ln_ev_fin(EE, p);
			//std::cout << "  " <<  ev;
			// Does not use all cases but good enough
			if (cp == 2) {
				htc = htc+ev/4.0;
			}
			else {
				htc = htc+ev/3.0;
			}
		}
		//std::cout << std::endl;
	}
	double ans;
	htc.doublify(ans);

	return ans;
}



// ----------------------------------------
//  Height Constant Code
// ----------------------------------------

double
elliptic_curve< bigint >::get_height_constant()
{
	debug_handler("elliptic_curve< bigint >", "get_height_constant(x, y, z) const");

	if (e == NULL)
		lidia_error_handler("elliptic_curve< bigint >::"
				    "get_height_constant(x, y, z) const",
				    "e == NULL");

	if (e->height_constant >= 0) {
		return e->height_constant;
	}
	e->height_constant = silverman(*this);

// bug in siksek()
//	double te1 = siksek(*this);
//	if (te1 < height_constant) {
//		height_constant = te1;
//	}
	return e->height_constant;
}



//
// input / output
//

void
elliptic_curve< bigint >::read(std::istream& in)
{
	debug_handler("elliptic_curve< bigint >", "input(std::istream, ec)");

	if (e != NULL)
		e->reset();
	base_elliptic_curve< bigint >::read(in);
	e->init();
}



void
elliptic_curve< bigint >::output_long(std::ostream & out) const
{
	debug_handler("elliptic_curve< bigint >", "output_long(out) const");

	if (e == NULL)
		lidia_error_handler("elliptic_curve< bigint >::"
				    "output_long(out) const",
				    "e == NULL");

	output_short(out); out << std::endl;
	out << "b2 = " << e->b2 << "; ";
	out << "b4 = " << e->b4 << "; ";
	out << "b6 = " << e->b6 << "; ";
	out << "b8 = " << e->b8 << std::endl;
	out << "c4 = " << e->c4 << "; ";
	out << "c6 = " << e->c6 << std::endl;
	out << "disc = " << e->delta << "; ";
	out << "bad primes = " << e->bad_primes << "; ";
	out << "# real components = " << e->conncomp << std::endl;
	out << "Conductor = " << e->N << std::endl << std::endl;
	out << "Reduction type at bad primes:\n";
	out << "p\tord(d)\tord(N)\tord(denom(j))\tKodaira\tc_p\n";
	long i, np = e->bad_primes.get_size();
	for (i = 0; i < np; i++)
		out << (e->bad_primes)[i] << "\t" << (e->reduct_array)[i] << std::endl;
	out << "Product of c_p = " << get_product_cp() << std::endl << std::endl;
}



void
elliptic_curve< bigint >::write (std::ostream & out) const
{
	debug_handler("elliptic_curve< bigint >", "operator << (std::ostream, ec)");

	if (e == NULL)
		lidia_error_handler("elliptic_curve< bigint >::"
				    "operator << (std::ostream, ec)",
				    "e == NULL");

	switch(e->get_output_mode()) {
	case elliptic_curve_flags::LONG:
		output_long(out);
		break;
	case elliptic_curve_flags::PRETTY:
		output_pretty(out);
		break;
	case elliptic_curve_flags::TEX:
		output_tex(out);
		break;
	default:
	case elliptic_curve_flags::SHORT:
		output_short(out);
		break;
	}
}



void
elliptic_curve< bigint >::output_torsion(std::ostream & out) // NOT const
{
	if (e->torsion.size() == 0) {
		make_torsion();
	}
	long nt = e->torsion.size(), nt2 = 0, oT, i;
	out << "Torsion order = " << nt << "\n";
	out << "Torsion points:\n";
	for (i = 0; i < nt; i++) {
		out << "\t" << i << ") ";
		point< bigint > T = (e->torsion)[i];
		oT = T.get_order();
		if (oT == 2)
			nt2++;
		out << T << ", \torder " << oT << std::endl;
	}
	if (nt2 > 2)
		out << "Non-cyclic: C2 x C" << (nt/2);
	else
		out << "Cyclic: C" << nt;
	if (nt > 1) {
		if (nt2 > 2)
			out << ", generators " << (e->torsion)[nt/2] << ", " << (e->torsion)[1];
		else
			out << ", generator  " << (e->torsion)[1];
	}
	out << std::endl;
}



//
// real height
//

bigfloat
real_height(const elliptic_curve< bigint > & ec, const bigfloat& xx)
{
	bigint bb2, bb4, bb6, bb8;
	ec.get_bi(bb2, bb4, bb6, bb8);
	bigfloat b2 = bigfloat(bb2), b4 = bigfloat(bb4), b6 = bigfloat(bb6), b8 = bigfloat(bb8);
	bigfloat b2dash = b2 - 12;
	bigfloat b4dash = b4 - b2 + 6;
	bigfloat b6dash = b6 - 2*b4 + b2 - 4;
	bigfloat b8dash = b8 - 3*b6 + 3*b4 - b2 + 3;

	bigfloat t, w, zz, zw;
	bigfloat H = 4.0; //max(4.0, max(abs(b2), max(2*abs(b4), max(2*abs(b6), abs(b8)))));
	t = abs(b2); if (t > H) H = t;
	t = 2*abs(b4); if (t > H) H = t;
	t = 2*abs(b6); if (t > H) H = t;
	t = abs(b8); if (t > H) H = t;
	long precision = bigfloat::get_precision();

	bigfloat tt1 = (5.0/3.0)*precision;
	bigfloat tt2 = log(7.0 + (4.0/3.0)*log(H));
	tt1 = ceil(tt1 + 0.5 + 0.75*tt2);
	long nlim;
	tt1.longify(nlim);

	long beta;
	if (abs(xx) < 0.5) {
		t = 1 / (xx + 1);
		beta = 0;
	}
	else {
		t = 1 / xx;
		beta = 1;
	}

	bigfloat mu = -log(abs(t)), dmu;
	bigfloat f = 1.0;

	for (long n = 0; n <= nlim; n++) {
		f /= 4.0;
		if (beta) {
			w = (((b6*t + 2*b4)*t + b2)*t + 4)*t;
			zz = 1 - t*t*(b4 + t*(2*b6 + t*b8));
			zw = zz + w;
		}
		else {
			w = (((b6dash*t + 2*b4dash)*t + b2dash)*t + 4)*t;
			zz = 1 - t*t*(b4dash + t*(2*b6dash + t*b8dash));
			zw = zz - w;
		}
		if (abs(w) <= 2*abs(zz)) {
			dmu = f*log(abs(zz));
			mu += dmu;
			t = w/zz;
		}
		else {
			dmu = f*log(abs(zw));
			mu += dmu;
			t = w/zw;
			beta = ! beta;
		}
	}
	return mu;
}



//
// transformation to curve
//

elliptic_curve< bigint >
trans_to_curve(complex_periods& cp)
{
	bigint c4, c6;
	bigcomplex cc4, cc6;
	bigfloat rc4, rc6;

	make_c4_c6(cp, cc4, cc6);
	if (abs(cc4.imag()) > 1.0e-10) {
		lidia_error_handler("complex_periods",
				    "trans_to_curve::possible error taking c4 as real");
	}
	if (abs(cc6.imag()) > 1.0e-10) {
		lidia_error_handler("complex_periods",
				    "trans_to_curve::possible error taking c6 as real");
	}

	rc4 = cc4.real();
	rc6 = cc6.real();

	rc4.bigintify(c4);
	rc6.bigintify(c6);

	return elliptic_curve< bigint > (c4, c6);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
