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
//                Adaption of John Cremona's code
//                (some of which itself is based on code of
//                Oisin McGuiness and Richard Pinch)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigint.h"
#include	"LiDIA/base_vector.h"
#include	"LiDIA/sort_vector.h"
#include	"LiDIA/modular_operations.inl"
#include	"LiDIA/nmbrthry_functions.h"
#include	"LiDIA/rational_factorization.h"

#include	"LiDIA/elliptic_curves/point_operations_bigint.h"
#include	"LiDIA/point_bigint.h"

#include	"LiDIA/complex_periods.h"
#include	"LiDIA/reduction_type.h"
#include	"LiDIA/elliptic_curve_flags.h"
#include	"LiDIA/elliptic_curves/base_elliptic_curve_rep.h"
#include	"LiDIA/elliptic_curves/elliptic_curve_rep.h"
#include	"LiDIA/elliptic_curves/elliptic_curve_rep_bigint.h"




#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



// Tates algorithm is now implemented
//
// subroutines -- not general purpose (.cc only)

// the e'th root of aa, mod p
bigint
root (const bigint& aa, int e, const bigint & p)
{
	bigint ans, temp; bool found = false;
	bigint a, i;
	best_remainder(a, aa, p);

	for (i = 1; (! found); i++) {
		ans = i;
		if (e == 2)
			temp = ans*ans - a;
		else
			temp = ans*ans*ans - a;
		found = is_divisor(p, temp);
	}
	return ans;
}



// test if quadratic aX^2 + bX + c = 0 (mod p) has roots

int
roots_exist (const bigint& aa, const bigint& bb, const bigint& cc , const bigint& p)
{
	const bigint& a = aa % p;
	const bigint& b = bb % p;
	const bigint& c = cc % p;
	const bigint& temp = (a*b*c)%p;
	if (p == 2)
		return temp.is_zero();
	if (a.is_zero())
		return 1;
	const bigint& d = b*b - 4*a*c;
	return jacobi(d, p) >= 0; //ie true if jacobi(d, p) = 0, 1
}



// find the number of roots of X^3 + bX^2 + cX + d = 0 (mod p)

int
nrootscubic (const bigint& bb, const bigint& cc, const bigint& dd, const bigint& p)
{
	bigint b, c, d, i;
	best_remainder(b, bb, p);
	best_remainder(c, cc, p);
	best_remainder(d, dd, p);

	bigint root1;
	int found = 0;
	for (i = 0; i < p; i++) {
		bigint temp = (((i+b)*i + c)*i + d)%p;
		if (temp.is_zero()) {
			root1 = i;
			found = 1;
			break;
		}
	}
	if (found) {
		const bigint& e = bb + root1;
		const bigint& f = cc + e*root1;
		if (roots_exist(1, e, f, p))
			return 3;
		else
			return 1;
	}
	else
		return 0;
}



//
// elliptic_curve<bigint> code starts here
//

//
// constructors / destructor
//

elliptic_curve_rep< bigint >::elliptic_curve_rep ()
	: base_elliptic_curve_rep< bigint > (), reduct_array(0)
{
	N = 0; height_constant = -1; conncomp = -1;
	lock = false;
}



// construct by Tate's algoritm
// arg E need not be minimal, but the reduced form will be
// Will minimalize the curve

elliptic_curve_rep< bigint >::elliptic_curve_rep (const bigint& x1,
						  const bigint& x2,
						  const bigint& x3,
						  const bigint& x4,
						  const bigint& x6,
						  elliptic_curve_flags::curve_model m)
	: base_elliptic_curve_rep< bigint > (x1, x2, x3, x4, x6, m), reduct_array(0)
{
	init();

	if (m != elliptic_curve_flags::PROJECTIVE) {
		lidia_error_handler("elliptic_curve_rep< bigint >"
				    "::elliptic_curve_rep(x1, x2, x3, x4, x6, m)",
				    "Only projective model allowed.");
	}

	lock = false;
}



// NB The following constructor constructs a curve from given c4, c6
// invariants, first checking that they are valid.  This must NOT be
// confused with the two-parameter base_elliptic_curve constructor which sets
// a4 and a6, with a1=a2=a3=0.  This one is more useful for curves/Z
// where we need to be able to construct them via periods and c4, c6

elliptic_curve_rep< bigint >::elliptic_curve_rep (const bigint& cc4,
						  const bigint& cc6,
						  elliptic_curve_flags::curve_model m)
	: base_elliptic_curve_rep< bigint > (), reduct_array(0)
{
	init(cc4, cc6);

	if (m != elliptic_curve_flags::PROJECTIVE) {
		lidia_error_handler("elliptic_curve_rep< bigint >"
				    "::elliptic_curve_rep(cc4, cc6, m)",
				    "Only projective model allowed.");
	}
}



elliptic_curve_rep< bigint >::elliptic_curve_rep (const elliptic_curve_rep< bigint > & E)
	: base_elliptic_curve_rep< bigint > (E)
{
	this->operator = (E);
	lock = false;
}



elliptic_curve_rep< bigint >::~elliptic_curve_rep ()
{
	if (reduct_array != 0)
		delete [] reduct_array;
}



//
// Functions for initialization.
//

// The following uses Kraus-Laska-Connell; c4 and c6 are input
// and will be overwritten; disc will hold the minimal disc,
// bad_p its prime factors, and u the scaling factor needed.

void
minimise_c4c6 (bigint& c4, bigint& c6,
	       bigint& disc, bigint& u, sort_vector< bigint > & bad_p)
{
	int u_is_1 = 1;
	const bigint& c62 = c6*c6;
	bigint p, temp;
	long index, np, d, a, b;

	u.assign_one();
	const bigint& g = gcd(c62, disc);
	rational_factorization g_fact = factor(g);
	np = g_fact.no_of_comp();
	if (g > 1)  // NB if g == 1 then g_fact includes p = 1 in its factor base!
	{
		for (index = 0; index < np; index++) {
			p = g_fact.base(index);
			d = static_cast<long>(std::floor(valuation(p, g)/12.0));
			if (p == 2) {
				((c4 >> (4*d)) % 16).longify(a);
				((c6 >> (6*d)) % 32).longify(b);
				if (b < 0)
					b += 32;
				if (((b%4) != 3) && !((a == 0) && ((b == 0) || (b == 8)))) {
					d--;
				}
			}
			//else if (p==3) if (valuation(3,c6)==(6*d + 2)) d--;
			else if (p == 3 && !c6.is_zero()) {
				if (valuation(3, c6) == (6*d + 2))
					d--;
			}

			if (d > 0) {
				power(temp, p, d);
				u *= temp;
				u_is_1 = 0;
			}
		}
		if (!u_is_1) {
			bigint u2, u4, u6, u12;
			square(u2, u);
			square(u4, u2);
			u6 = u4*u2;
			square(u12, u6);
			c4 /= u4;
			c6 /= u6;
			disc /= u12;
		}
	}

	// Now fill the array of bad primes;
	// note that we have already factorised g which is a factor of disc
	// so we use those factors first
	bad_p.set_mode(vector_flags::expand);
	temp = disc;
	if (temp.is_negative())
		temp.negate();
	if (g > 1) {
		for (index = 0; index < np; index++) {
			p = g_fact.base(index);
			if (is_divisor(p, disc)) {
				bad_p.insert(p);
				temp /= p;
				while (is_divisor(p, temp))
					temp /= p;
			}
		}
	}
	// Now factor the rest of disc:
	if (temp > 1) {
		rational_factorization x = factor(temp);
		np = x.no_of_comp();
		for (index = 0; index < np; index++) {
			bad_p.insert(x.base(index));
		}
	}
}



// Reset function: resets everything ready for a new curve
void
elliptic_curve_rep< bigint >::reset()
{
	periods.reset();
	torsion.reset();
	height_constant = -1;
	conncomp = -1;
	bad_primes.reset();
}



// Initialize curve and reduce it

// First init() starts from a1,...,a6

void
elliptic_curve_rep< bigint >::init (const bigint& aa1, const bigint& aa2, const bigint& aa3,
				    const bigint& aa4, const bigint& aa6)
{
	reset();
	a1 = aa1; a2 = aa2; a3 = aa3; a4 = aa4; a6 = aa6;
	compute_invariants();
	init(); //NB init() only needs c4, c6, delta
}



int
valid_c4c6 (const bigint& c4, const bigint& c6, bigint& disc)
{
	disc = c4*c4*c4 - c6*c6;
	if (sign(disc) == 0) return 0; // singular
	if (!is_divisor(1728, disc)) return 0; // need c4^3-c6^2 = 1728D, with D |= 0
	disc /= 1728;
	long x6;

	best_remainder(x6, c6, 27);
	if ((x6 == 9) || (x6 == -9))
		return 0; // need c6 != +-9 (mod 27)
	best_remainder(x6, c6, 4);
	if (x6 == -1)
		return 1; // OK if c6 = -1 (mod 4)
	if (!is_divisor(16, c4))
		return 0; // else need c4 = 0 (mod 16)
	best_remainder(x6, c6, 32); //      and
	return ((x6 == 0) || (x6 == 8)); //      c6 = 0, 8 (mod 32).
}



// Second init() starts from c4,c6, checks for validity

void
elliptic_curve_rep< bigint >::init (const bigint& cc4, const bigint& cc6)
{
	if (!valid_c4c6(cc4, cc6, delta)) {
		lidia_error_handler("elliptic_curve_rep< bigint >",
				    "init::invalid c4, c6 invariants");
		return;
	}
	reset();
	c4 = cc4; c6 = cc6;
	init(); //NB init() only needs c4, c6, delta
}



// init() with no parameters; needs c4, c6, delta to be set already,
// either by base_elliptic_curve constructor or by init(a1,...a6) or init(c4,c6)

void
elliptic_curve_rep< bigint >::init ()
{
	// Laska-Kraus-Connell reduction
	bigint u; // not used here
	minimise_c4c6(c4, c6, delta, u, bad_primes);

	// now compute minimal equation
	best_remainder(b2, -c6, 12);
	const bigint& b22 = b2*b2;
	b4 = (b22-c4)/24;
	b6 = (-b2*b22+36*b2*b4-c6)/216;
	b8 = (b2*b6 - b4*b4) / 4;

	a1 = (b2.is_odd() ? 1 : 0);
	a3 = (b6.is_odd() ? 1 : 0);
	a2 = (b2-a1*a1)/4;
	a4 = (b4-a1*a3)/2;
	a6 = (b6-a3*a3)/4;

	if (a1.is_zero() && a2.is_zero() && a3.is_zero())
		cp = elliptic_curve_flags::SHORT_W;
	else
		cp = elliptic_curve_flags::LONG_W;
	bcp.set(cp, elliptic_curve_flags::PROJECTIVE, 0);

	conncomp = delta.sign() > 0 ? 2 : 1;
	long nbadp = bad_primes.get_size();
	if (reduct_array != NULL)
		delete[] reduct_array;
	reduct_array = new reduction_type[nbadp];

	tate();
}



void
elliptic_curve_rep< bigint >::tate ()
{
// This function (called by init) carrries out Tate's algorithm for
// each bad prime, in order to fill all the entries in reduct_array

	// local variables
	base_elliptic_curve_rep< bigint > C(*this);
	bigint p, halfmodp, temp, r, s, t, b, c, bb, cc, bc, d, w, x, mx, my;
	bigint a2t, a3t, a4t, a6t;
	int ord_p_discr, ord_p_j, c_p = 1, pdiv2, pdiv3, sw, loop, ix, iy;
	long nbadp = bad_primes.get_size();
	long index;

	// main loop - for each of the prime divisors of the discriminant.

	for (index = 0; index < nbadp; index++) {
		p = bad_primes[index];
		ord_p_discr = valuation(p, delta);
		//ord_p_j = ord_p_discr - 3*valuation(p,C.get_c4());
		//if (ord_p_j < 0) ord_p_j = 0;
		if (C.get_c4().is_zero()) {
			ord_p_j = 0;
		}
		else {
			ord_p_j = ord_p_discr - 3*valuation(p, C.get_c4());
			if (ord_p_j < 0) ord_p_j = 0;
		}

		halfmodp = (1 + p) / 2;
		pdiv2 = (p == 2);
		pdiv3 = (p == 3);

		//change coords so that p|C.a3,C.a4,C.a6
		if (pdiv2) {
			if (is_divisor(p, C.get_b2())) {
				r = root(C.get_a4(), 2, p);
				t = root(((r+C.get_a2())*r+C.get_a4())*r+C.get_a6(), 2, p);
			}
			else {
				inv_mod(temp, C.get_a1(), p);
				if (temp.is_negative())
					temp += p;
				r = temp*C.get_a3();
				t = temp*(C.get_a4() + r*r);
			}
		}
		else if (pdiv3) {
			if (is_divisor(p, C.get_b2()))
				r = root(-C.get_b6(), 3, p);
			else {
				inv_mod(r, C.get_b2(), p);
				if (r.is_negative()) r += p;
				r = -r*C.get_b4();
			}
			t = C.get_a1()*r + C.get_a3();
		}
		else {
			if (is_divisor(p, C.get_c4())) {
				inv_mod(r, 12, p);
				if (r.is_negative())
					r += p;
				r = -r*C.get_b2();
			}
			else {
				inv_mod(r, 12*C.get_c4(), p);
				if (r.is_negative()) r += p;
				r = -r*(C.get_c6()+C.get_b2()*C.get_c4());
			}
			t = -halfmodp*(C.get_a1()*r+C.get_a3());
		}
		best_remainder(r, r, p);
		best_remainder(t, t, p);
		C.transform(r, 0, t);

		// test for Types In, II, III, IV
		if (!is_divisor(p, C.get_c4())) {
			temp = -C.get_a2();
			if (roots_exist(1, C.get_a1(), temp, p))
				c_p = ord_p_discr;
			else {
				if ((ord_p_discr%2) != 0)
					c_p = 1;
				else
					c_p = 2;
			}
			reduct_array[index] = reduction_type
				(ord_p_discr, 1, ord_p_j, 10*ord_p_discr, c_p);
			continue;
		}  // Type In (n = ord_p_discr)
		//else if ( valuation(p,C.get_a6()) < 2 )
		else if (!C.get_a6().is_zero() && valuation(p, C.get_a6()) < 2) {
			reduct_array[index] = reduction_type
				(ord_p_discr, ord_p_discr, ord_p_j, 2, 1);
			continue;
		} // Type II
		//else if ( valuation(p,C.get_b8()) < 3 )
		else if (!C.get_b8().is_zero() && valuation(p, C.get_b8()) < 3) {
			reduct_array[index] = reduction_type
				(ord_p_discr, ord_p_discr - 1, ord_p_j, 3, 2);
			continue;
		} // Type III
		//else if ( valuation(p,C.get_b6()) < 3 )
		else if (!C.get_b6().is_zero() && valuation(p, C.get_b6()) < 3) {
			temp = -(C.get_a6()/p)/p;
			bigint temp2 = C.get_a3()/p;
			if (roots_exist(1, temp2, temp, p))
				c_p = 3;
			else
				c_p = 1;
			reduct_array[index] = reduction_type(ord_p_discr, ord_p_discr - 2, ord_p_j, 4, c_p);
			continue;
		} // Type IV

		// else change coords so that
		//     p|C.a1,C.a2, p^2|C.a3,C.a4, p^3|C.a6
		if (pdiv2) {
			s = root(C.get_a2(), 2, p);
			t = p*root((C.get_a6()/p)/p, 2, p);
		}
		else if (pdiv3) {
			s = C.get_a1();
			t = C.get_a3();
		}
		else {
			s = -C.get_a1()*halfmodp;
			t = -C.get_a3()*halfmodp;
		}
		C.transform(0, s, t);

		//                             3     2
		// Analyse roots of the cubic T  + bT  + cT + d = 0, where
		// b=C.a2/p, c=(C.a4/p)/p, d=((C.a6/p)/p)/p
		b = C.get_a2()/p;
		c = (C.get_a4()/p)/p;
		d = ((C.get_a6()/p)/p)/p;
		bb = b*b; cc = c*c; bc = b*c;
		w = 27*d*d - bb*cc + 4*b*bb*d - 18*bc*d + 4*c*cc;
		x = 3*c - bb;

		if (is_divisor(p, w)) {
			if (is_divisor(p, x))
				sw = 3;
			else
				sw = 2;
		}
		else
			sw = 1;

		switch (sw) {
		case 1:
			//Three distinct roots - Type I*0
			reduct_array[index] = reduction_type
				(ord_p_discr, ord_p_discr - 4, ord_p_j, 1, 1+nrootscubic(b, c, d, p));
			break;

		case 2:
			// One double root - Type I*m for some m
			// Change coords so that the double root is T=0 mod p
			if (pdiv2)
				r = root(c, 2, p);
			else if (pdiv3) {
				inv_mod(r, b, p);
				if (r.is_negative())
					r += p;
				r *= c;
			}
			else {
				inv_mod(r, 2*x, p);
				if (r.is_negative())
					r += p;
				r *= (bc - 9*d);
			}
			best_remainder(r, r, p);
			r *= p;
			C.transform(r, 0, 0);

			ix = 3; iy = 3; mx = p*p; my = p*p;
			loop = 1;
			while (loop) {
				a2t = C.get_a2()/p;
				a3t = C.get_a3()/my;
				a4t = (C.get_a4()/p)/mx;
				a6t = (C.get_a6()/mx)/my;
				temp = a3t*a3t + 4*a6t;
				if (is_divisor(p, temp)) {
					if (pdiv2)
						t = my*root(a6t, 2, p);
					else {
						best_remainder(t, -a3t*halfmodp, p);
						t *= my;
					}
					C.transform(0, 0, t);
					my = my*p;
					iy++;
					a2t = C.get_a2()/p;
					a3t = C.get_a3()/my;
					a4t = (C.get_a4()/p)/mx;
					a6t = (C.get_a6()/mx)/my;
					temp = a4t*a4t - 4*a6t*a2t;
					if (is_divisor(p, temp)) {
						if (pdiv2) {
							inv_mod(r, a2t, p);
							if (r.is_negative())
								r += p;
							r = mx*root(a6t*r, 2, p);
						}
						else {
							inv_mod(r, 2*a2t, p);
							if (r.is_negative())
								r += p;
							best_remainder(r, -a4t*r, p);
							r *= mx;
						}
						C.transform(r, 0, 0);
						mx = mx*p;
						ix++; // and stay in loop
					}
					else {
						if (roots_exist(a2t, a4t, a6t, p))
							c_p = 4;
						else
							c_p = 2;
						loop = 0;
					}  // and exit loop
				}
				else {
					temp = -a6t;
					if (roots_exist(1, a3t, temp, p))
						c_p = 4;
					else
						c_p = 2;
					loop = 0;
				}
			}
			reduct_array[index] = reduction_type
				(ord_p_discr, ord_p_discr - ix - iy + 1, ord_p_j,
				 10 * (ix + iy) - 49, c_p);
			break; // Type I*m

		case 3:
			// Triple root
			// change coords so that T=0 mod p
			if (pdiv2)
				r = b;
			else if (pdiv3)
				r = root(-d, 3, p);
			else {
				inv_mod(r, 3, p);
				if (r.is_negative())
					r += p;
				r *= -b;
			}
			best_remainder(r, r, p);
			r *= p;
			C.transform(r, 0, 0);

			a3t = (C.get_a3()/p)/p;
			a6t = (((C.get_a6()/p)/p)/p)/p;

			// test for Type IV*
			temp = a3t*a3t + 4*a6t;
			if (!is_divisor(p, temp)) {
				temp = -a6t;
				if (roots_exist(1, a3t, temp, p))
					c_p = 3;
				else
					c_p = 1;
				reduct_array[index] = reduction_type
					(ord_p_discr, ord_p_discr - 6, ord_p_j, 5, c_p);
				break;
			}
			// change coords so that p^3|C.a3, p^5|C.a6
			if (pdiv2)
				t = -p*p*root(a6t, 2, p);
			else {
				best_remainder(t, -a3t*halfmodp, p);
				t *= p*p;
			}
			C.transform(0, 0, t);

			// test for types III*, II*
			//if ( valuation(p,C.get_a4()) < 4 )
			if (!C.get_a4().is_zero() && valuation(p, C.get_a4()) < 4) {
				reduct_array[index] = reduction_type
					(ord_p_discr, ord_p_discr - 7, ord_p_j, 6, 2);
				break;
			}  // Type III*
			//else if ( valuation(p,C.get_a6()) < 6 )
			else if (!C.get_a6().is_zero() && valuation(p, C.get_a6()) < 6) {
				reduct_array[index] = reduction_type
					(ord_p_discr, ord_p_discr - 8, ord_p_j, 7, 1);
				break;
			} // Type II*

			else {
				lidia_error_handler("elliptic_curve_rep< bigint >", "init: Tate's algorithm reached end of loop");
			}
		}; // end switch
	}     // end primes for-loop

// finally compute conductor from its exponents:

	N.assign_one();
	for (index = 0; index < nbadp; index++) {
		p = bad_primes[index];
		power(temp, p, reduct_array[index].get_ord_p_N());
		N *= temp;
	}
	return;
}



// end of Tate's algorithm



//
// Assignment
//

elliptic_curve_rep< bigint > &
elliptic_curve_rep< bigint >::operator = (const elliptic_curve_rep< bigint > & E)
{
	if (this == &E)
		return *this;

	a1 = E.a1; a2 = E.a2; a3 = E.a3; a4 = E.a4; a6 = E.a6;
	b2 = E.b2; b4 = E.b4; b6 = E.b6; b8 = E.b8;
	c4 = E.c4; c6 = E.c6;
	delta = E.delta;
	bcp = E.bcp; cp = E.cp;

#if 0
	// Copy the torsion across: NB The simpler-looking
	//    assign(torsion,E.torsion);
	// does NOT work properly,  as simply copying the points
	// leaves the copies pointing to the old curve

	long index, nt = E.torsion.get_size();
	torsion.set_capacity(nt);
	torsion.set_size(0);
	for (index = 0; index < nt; index++) {
		torsion[index] = E.torsion[index];
		torsion[index].ec = this;
	}

#endif

	// New structure; does work now. MM
	torsion.assign(E.torsion);

	conncomp = E.conncomp;
	bad_primes = E.bad_primes;
	periods = E.periods;
	height_constant = E.height_constant;

	long nbadp = E.bad_primes.get_size();
	delete[] reduct_array;
	reduct_array = new reduction_type[nbadp];
	long index;
	for (index = 0; index < nbadp; index++)
		reduct_array[index] = E.reduct_array[index];
	N = E.N;

	return *this;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
