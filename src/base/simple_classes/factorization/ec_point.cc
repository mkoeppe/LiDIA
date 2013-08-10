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
//	Author	: Andreas Mueller (AM), Volker Mueller  (VM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/rational_factorization.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



void mecgen16(ec_curve & me, ec_point_M & mp, bigint & factor, long m)
{
	ec_curve we(bigmod(-8L), bigmod(-32L));
	ec_point_W wq(we, bigmod(12L), bigmod(40L));

	bigmod h1, h2;
	bigmod alpha, la, lb, lc, ld;

	multiply(wq, wq, m, factor);

	subtract (h2, wq.X(), 9L);

	factor = h2.invert(1);
	if (!factor.is_one())
		return;

	add(h1, wq.Y(), 25L);
	multiply(h2, h2, h1);
	inc(h2);

	factor = h2.invert(1);
	if (!factor.is_one())
		return;

	alpha.assign(h2);

	square(h2, h2);
	multiply(h2, h2, 8L);
	dec(h2);

	factor = h2.invert(1);
	if (!factor.is_one())
		return;

	multiply(h1, alpha, 4L);
	inc(h1);
	multiply(h1, h1, alpha);
	multiply(h1, h1, 2L);
	multiply(ld, h1, h2);

	h1.assign(ld);
	factor = h1.invert(1);
	if (!factor.is_one())
		return;

	multiply(h2, ld, 2L);
	dec(h2);
	multiply(h1, h1, h2);
	subtract(h2, ld, 1L);
	multiply(lc, h1, h2);

	multiply(lb, lc, ld);

	invert(h2, bigmod(4L));

	multiply(h1, ld, 2L);
	dec(h1);

	multiply(h1, h1, h2);
	h1.negate();

	la = h1;

	subtract(h2, ld, 1L);
	square(h2, h2);

	factor = h2.invert(1);
	if (!factor.is_one())
		return;

	h1 = ld;
	square(h1, h1);

	add(la, la, ld);
	subtract(la, la, h1);
	multiply(la, la, h2);

	mp.assign(la, bigmod(1L));

	subtract(h1, ld, 1L);
	multiply(h1, h1, ld);
	h1.multiply_by_2();

	factor = h1.invert(1);

	h2.assign(ld);
	h2.multiply_by_2();
	dec(h2);
	square(h2, h2);
	multiply(h2, h2, h1);

	square(h2, h2);

	subtract(h2, h2, 2L);

	me.assign(h2, bigmod(1L));

}



void
add(ec_point_M & erg, const ec_point_M & mp1, const ec_point_M & mp2, const ec_point_M & diff)
{
	bigmod t1, t2, h3, h4;

	subtract(t1, mp1.x, mp1.z);
	add(h3, mp2.x, mp2.z);
	multiply(t1, t1, h3);

	add(t2, mp1.x, mp1.z);
	subtract(h3, mp2.x, mp2.z);
	multiply(t2, t2, h3);

	add(h3, t1, t2);
	square(h3, h3);
	multiply(h4, h3, diff.z);

	subtract(h3, t1, t2);
	square(h3, h3);
	multiply(h3, h3, diff.x);

	erg.assign(h4, h3);
}



void
multiply_by_2(ec_point_M & erg, const ec_point_M & mp, const bigmod & a_2)
{
	bigmod h1, h2, h3;

	add(h1, mp.x, mp.z);
	square(h1, h1);

	subtract(h2, mp.x, mp.z);
	square(h2, h2);

	multiply(h3, h1, h2);
	erg.x.assign(h3);

	subtract(h1, h1, h2);
	multiply(h3, h1, a_2);
	add(h3, h3, h2);
	multiply(h3, h1, h3);
	erg.z.assign(h3);
}



void
multiply(ec_point_M & erg, const ec_point_M & mp, const bigmod & a_2, unsigned long p, unsigned long m)
{
	long i;
	unsigned long r;
	unsigned long bit;

	ec_point_M P1 (mp);
	ec_point_M P2, Erg1;

	multiply_by_2(P2, mp, a_2);

	r = integer_log(m) - 1; //  r = erstes gesetztes Bit in m

	bit = 1 << r; //  bit = 2^r

	if (r == 31)
		p = -m; //  p = 2^32 - m = -m
	else {
		p = 1 << (r+1);
		p = p - m; //  p = 2^(r+1) - m
	}

	for (i = r-1; i > 0; i--) {
		add(Erg1, P1, P2, mp); //  Erg1 = P1 + P2 mit P1-P2 = P

		bit >>= 1;
		if (bit & p) {
			//  i.th Bit of p is 1
			multiply_by_2(P1, P1, a_2);
			P2.assign(Erg1);
		}
		else {
			// i.th Bit of p is 0
			multiply_by_2(P2, P2, a_2);
			P1.assign(Erg1);
		}
	}

	add(erg, P1, P2, mp);
}



// Weierstrass parametrization


//
// comparisons
//

bool operator == (const ec_point_W & wp1, const ec_point_W & wp2)
{
	if ((wp1.is_0 == wp2.is_0) &&
	    (wp1.is_0 || (wp1.x == wp2.x && wp1.y == wp2.y)))
		return 1;
	else
		return 0;
}



//
// Procedural versions
//


bigint trans(ec_curve & ell_curve, ec_point_W & wp, ec_point_M & mp)
{
	bigint factor;
	bigmod la, inv_3, inv_t, h1, h2, h3;

	h2 = mp.X();

	square(h2, h2);

	multiply(h1, mp.X(), ell_curve.A());
	add(inv_t, h1, h2);
	inc(inv_t);
	multiply(inv_t, inv_t, mp.X());

	factor = inv_t.invert(1);
	if (!factor.is_one()) return factor;

	wp.assign_one();
	wp.y.assign(inv_t);

	invert(inv_3, bigmod(3));

	multiply(h3, inv_3, ell_curve.A());
	add(h3, h3, mp.X());
	multiply(h3, h3, inv_t);

	wp.x.assign(h3);

	h2 = ell_curve.A();
	square(h2, h2);
	multiply(h1, h2, inv_3);
	h3.assign_one();
	subtract(h3, h3, h1);

	h1 = inv_t;
	square(h1, h1);
	multiply(h3, h3, h1);

	la = h3;

	multiply(h1, h1, inv_t);
	multiply(h1, h1, inv_3);

	multiply(h3, inv_3, inv_3);
	multiply(h2, h2, h3);
	multiply(h2, h2, 2L);
	dec(h2);
	multiply(h2, h2, ell_curve.A());

	multiply(h2, h2, h1);

	ell_curve.assign(la, h2);
	wp.curve = (ec_curve *) &ell_curve;
	return (bigint(1));
}



void
multiply_by_2(ec_point_W & erg, const ec_point_W & wp, bigint & factor)
{
	bigmod h, h1, h2;

	if (wp.is_0) {
		factor.assign_one();
		erg.assign_zero();
		return;
	}

	h2.assign(wp.y);
	h2.multiply_by_2();
	factor = h2.invert(1);
	if (!factor.is_one())  return;

	h1.assign(wp.x);
	square(h1, h1);
	multiply(h1, h1, 3L);
	add(h1, h1, (*(wp.Curve())).A());

	multiply(h, h1, h2);
	h1.assign(h);

	square(h1, h1);
	h2.assign(wp.x);
	h2.multiply_by_2();
	subtract(h2, h1, h2);
	subtract(h1, wp.x, h2);
	multiply(h1, h, h1);
	subtract(h, h1, wp.y);

	erg.assign(h2, h);
}



void
addPQ(ec_point_W & erg, const ec_point_W & wp1, const ec_point_W & wp2, bigint & factor)
{
	bigmod h, h1, h2;

	subtract(h2, wp2.x, wp1.x);

	factor = h2.invert(1);
	if (!factor.is_one()) return;

	subtract(h1, wp2.y, wp1.y);

	multiply(h, h1, h2);
	h1.assign(h);
	square(h1, h1);
	subtract(h1, h1, wp1.x);
	subtract(h1, h1, wp2.x);

	subtract(h2, wp1.x, h1);
	multiply(h2, h2, h);
	subtract(h2, h2, wp1.y);

	erg.assign(h1, h2);
}



void
add(ec_point_W & erg, const ec_point_W & wp1, const ec_point_W & wp2, bigint & factor)
{
	if (wp1.is_0) {
		factor.assign_one();
		erg.assign(wp2);
	}
	else
		if (wp2.is_0) {
			erg.assign(wp1);
			factor.assign_one();
		}
		else
			if (wp1 == wp2) {
				bigmod h, h1, h2;

				h2.assign(wp2.y);
				h2.multiply_by_2();

				factor = h2.invert(1); // h2 = (2y1)^(-1)
				if (!factor.is_one())  return;

				h1.assign(wp2.x);
				square(h1, h1);
				multiply(h1, h1, 3L);
				add(h1, h1, (*(wp2.Curve())).A()); // h1 = 3x1^2 + a

				multiply(h, h1, h2); // h = (3x1^2+a)/(2y1)
				h1.assign(h);
				square(h1, h1); // h1 = h^2
				h2.assign(wp2.x);
				h2.multiply_by_2();
				subtract(h2, h1, h2);

				subtract(h1, wp2.x, h2);
				multiply(h1, h, h1);
				subtract(h, h1, wp2.y);

				erg.assign(h2, h);
			}
			else
				if (wp1.is_negative(wp2)) {
					factor.assign_one();
					erg.assign_zero();
				}
				else {
					bigmod h, h1, h2;

					subtract(h2, wp2.x, wp1.x);

					factor = h2.invert(1);
					if (!factor.is_one()) return;

					subtract(h1, wp2.y, wp1.y);

					multiply(h, h1, h2);
					h1.assign(h);
					square(h1, h1);
					subtract(h1, h1, wp1.x);
					subtract(h1, h1, wp2.x);

					subtract(h2, wp1.x, h1);
					multiply(h2, h2, h);
					subtract(h2, h2, wp1.y);

					erg.assign(h1, h2);
				}
}



void
multiply(ec_point_W & erg, const ec_point_W & wp, long m, bigint & factor)
{
	ec_point_W Z;

	if (m < 5) {
		Z.assign(wp);
		for (; m > 1; m--) {
			add(Z, Z, wp, factor);
			if (!factor.is_one())   return;
		}
	}
	else {
		ec_point_W X(wp);
		Z.assign_zero();

		while (m > 0) {
			while (!(m & 1))              // while m even:
			{                               // doubling point X
				multiply_by_2(X, X, factor);
				if (!factor.is_one())   return;
				m >>= 1;
			}

			add(Z, X, Z, factor);
			if (!factor.is_one())   return;
			m--;
		}
	}

	erg.assign(Z);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
