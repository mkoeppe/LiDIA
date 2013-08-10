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
//	Author	: Nigel Smart, John Cremona
//                Adaption of John Cremona's code; some code came
//                originally from Oisin McGuiness
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigint.h"
#include	"LiDIA/bigfloat.h"
#include	"LiDIA/nmbrthry_functions.h"
#include	"LiDIA/base_vector.h"

#include	"LiDIA/elliptic_curve_bigint.h"
#include	"LiDIA/point_bigint.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// ctors' and d'tor
//

point< bigint >::point()
	: x(0),
	  y(0),
	  z(0),
	  height(0),
	  ord(0)
{
	ec = new elliptic_curve< bigint > ();
}



point< bigint >::point(const bigint & xp,
		       const bigint & yp,
		       const elliptic_curve< bigint > &ecp)

	: x(xp),
	  y(yp),
	  z(1),
	  height(-1),
	  ord(0)

{
	ec = new elliptic_curve< bigint > (ecp);
}



point< bigint >::point(const bigint & xp,
			const bigint & yp,
			const bigint & zp,
			const elliptic_curve< bigint > &ecp)

	: x(xp),
	  y(yp),
	  z(zp),
	  height(-1),
	  ord(0)

{
	ec = new elliptic_curve< bigint > (ecp);
	reduce();
}



point< bigint >::point(const point< bigint > & P)

	: x(P.x),
	  y(P.y),
	  z(P.z),
	  height(P.height),
	  ord(P.ord)

{
	ec = new elliptic_curve< bigint > (*(P.ec));
}



point< bigint >::point(const elliptic_curve< bigint > &e)

	: x(0),
	  y(1),
	  z(0),
	  height(0),
	  ord(1)

{
	ec = new elliptic_curve< bigint > (e);
}



point< bigint >::~point()
{
	delete ec;
}



//
// Access
//

elliptic_curve< bigint >
point< bigint >::get_curve() const
{
	return *ec;
}



//
// Assignment
//

void
point< bigint >::assign_zero(const elliptic_curve< bigint > & e)
{
	ec->assign(e);
	z.assign_zero();
	y.assign_one();
	x.assign_zero();
	height = 0;
	ord = 1;
}



void
point< bigint >::assign(const point< bigint > & P)
{
	if (&P != this) {
		ec->assign(*(P.ec));
		x = P.x;
		y = P.y;
		z = P.z;
		height = P.height;
		ord = P.ord;
	}
}



void
point< bigint >::assign(const bigint & xx, const bigint & yy, const
			elliptic_curve< bigint > & e)
{
	ec->assign(e);
	x = xx;
	y = yy;
	z.assign_one();
	height = -1;
	ord = 0;
}



void
point< bigint >::assign(const bigint & xx, const bigint & yy, const bigint &zz,
			const elliptic_curve< bigint > & e)
{
	ec->assign(e);
	x = xx;
	y = yy;
	z = zz;
	height = -1;
	ord = 0;
}



// assign on last curve (use with care of course)
void
point< bigint >::assign_zero()
{
	z.assign_zero();
	y.assign_one();
	x.assign_zero();
	height = 0;
	ord = 1;
}



void
point< bigint >::assign(const bigint & xx, const bigint & yy)
{
	x = xx;
	y = yy;
	z.assign_one();
	height = -1;
	ord = 0;
}



void
point< bigint >::assign(const bigint & xx, const bigint & yy, const bigint &zz)
{
	x = xx;
	y = yy;
	z = zz;
	height = -1;
	ord = 0;
}



void
point< bigint >::swap(point< bigint > & P)
{
	LiDIA::swap(x, P.x);
	LiDIA::swap(y, P.y);
	LiDIA::swap(z, P.z);
	LiDIA::swap(height, P.height);

	int te = P.ord;
	P.ord = ord;
	ord = te;

	elliptic_curve< bigint > * e = P.ec;
	P.ec = ec;
	ec = e;
}



void
point< bigint >::reduce()
{
	if (z.is_zero()) {
		x = 0;
		y = 1;
		return;
	}
	if (z.is_one()) {
		return;
	}

	bigint d = gcd(x, y);
	if (!d.is_one())
		d = gcd(d, z); // now d = gcd(x, y, z)
	if (d.is_zero()) {
		lidia_error_handler("point< bigint >", "reduce::Invalid point");
	}
	if (!d.is_one()) {
		x /= d;
		y /= d;
		z /= d;
	}
	if (z.is_lt_zero()) {
		x.negate();
		y.negate();
		z.negate();
	}
}



//
// Testing points
//

bool
point< bigint >::is_equal (const point< bigint > & P) const
{
	if (*(P.ec) != *(ec))
		return false;

	bigint temp(P.x*z - P.z*x);

	if (!temp.is_zero()) {
		return false;
	}

	temp = P.y*z - P.z*y;
	if (!temp.is_zero()) {
		return false;
	}

	temp = P.y*x - P.x*y;
	return temp.is_zero();
}



bool
point< bigint >::on_curve() const
{
	if (x.is_zero() && y.is_zero() && z.is_zero()) {
		return 0;
	}
	if (z.is_zero()) {
		return 1;
	}
	// should calculate
	//                       y^2 +a1 x y + a3 y
	//  and
	//                       x^3 + a2 x^2 + a4 x + a6
	// where
	//          x = X/Z, y = Y/Z
	// and verify equality.
	//
	// In homogeneous coordinates:
	//
	// Lhs = Y^2 Z + a1 XYZ + a3 YZ^2 = (YZ) *(Y + a1 X + a3 Z)
	//
	//
	// Rhs = X^3 +a2 X^2 Z + a4 X Z^2 + a6 Z^3
	//
	const bigint& Lhs = y*z*(y + ec->get_a1()*x + ec->get_a3()*z);
	const bigint& Rhs = ec->get_a6()*z*z*z+x*(ec->get_a4()*z*z+x*(ec->get_a2()*z + x));
	return Lhs == Rhs;
}



//
// Arithmetic
//

void
add(point< bigint > & R,
    const point< bigint > & P,
    const point< bigint > & Q)
{
	if (*(P.ec) != *(Q.ec))
		lidia_error_handler("point< bigint >", "add::different elliptic curves");

	if (Q.z.is_zero()) {
		R.assign(P);
		return;
	}
	if (P.z.is_zero()) {
		R.assign(Q);
		return;
	}
	if (P == Q)
		(P.ec)->get_point_operations()._mult_by_2(R, P);
	else
		(P.ec)->get_point_operations()._add(R, P, Q);
}



void
subtract(point< bigint > & R,
	 const point< bigint > & P,
	 const point< bigint > & Q)
{
	point< bigint > H;
	negate(H, Q);
	add(R, P, H);
}



void
multiply_by_2(point< bigint > & R, const point< bigint > & P)
{
	if (P.z.is_zero()) {
		R.assign(P);
		return;
	}
	(P.ec)->get_point_operations()._mult_by_2(R, P);
}



void
multiply (point< bigint > & R, const bigint & mm, const point< bigint > & P)
{
	bigint m(mm);
	point< bigint > Z(P);

	R.assign(P);

	if (m.is_negative()) {
		m.negate();
		negate(Z, Z);
		negate(R, R);
	}

	if (P.z.is_zero() || mm.is_zero()) {
		R.assign(*(P.ec));
		return;
	}

	if (mm.is_one())
		return;

	int i = static_cast<int>(m.bit_length())-2; // note m >= 2

	while (i >= 0) {
		multiply_by_2(R, R);
		if (m.bit(i))
			add(R, R, Z);
		i--;
	}
}



void
negate(point< bigint > & r, const point< bigint > & x)
{
	(r.ec)->assign(*(x.ec));
	(x.ec)->get_point_operations()._negate(r, x);
}



// Advanced Functions

bigfloat
point< bigint >::get_naive_height() const
{
	if (z.is_zero()) {
		return 0;
	}
	bigint g = gcd(x, z);
	bigfloat xx = abs(x/g);
	if (abs(z/g) > abs(x/g)) {
		xx = abs(z/g);
	}
	return log(xx);
}



bigfloat
point< bigint >::get_height()
{
	// WARNING -- no check made of validity of point on curve
	if (z.is_zero())
		return bigfloat(0.0);
	if (height >= 0.0) return height; // already calculated it
	// zero height if torsion
	if (get_order() > 0) {
		height = bigfloat(0.0);
		return height;
	}

	// N.B. So if we ever ask a point its height it will compute its order.
	// otherwise need to calculate it

// Add local heights at finite primes dividing discr(E) OR denom(P).
// The components for primes dividing denom(P) add to log(denom(x(P)));
//   since P=(XZ:Y:Z^3), denom(P)=Z=gcd(XZ,Z^3), called "zroot" here,
//   and so the contribution is log(denom(x(P))) = 2*log(zroot).
//   This avoids factorizing the denominator.

	bigint p;
	bigint zroot(gcd(x, z)); // = cube root of Z
	bigfloat ans(get_realheight() + 2*log(bigfloat(zroot)));

	long i;
	base_vector< bigint > bp = ec->get_bad_primes();

	for (i = 0; i < bp.get_size(); i++) {
		p = bp[i];
		if (!is_divisor(p, zroot))
			ans += get_pheight(p);
	}

	height = ans;
	return ans;
}



bigfloat
point< bigint >::get_pheight(const bigint& pr) const
{
	bigint a1, a2, a3, a4, a6, b2, b4, b6, b8, c4, c6, discr;

	ec->get_ai(a1, a2, a3, a4, a6);
	ec->get_bi(b2, b4, b6, b8);
	ec->get_ci(c4, c6);

	discr = ec->discriminant();
	long n = valuation(pr, discr);

	const bigint& zroot = gcd(x, z); // = cube root of z
	long vpz = 3*valuation(pr, zroot);

	const bigint& x2 = x*x;
	const bigint& z2 = z*z;
	const bigint& xz = x*z;
	const bigint& yz = y*z;
	long a = valuation(pr, 3*x2 + 2*a2*xz + a4*z2 - a1*yz) - 2*vpz;
	long b = valuation(pr, 2*y + a1*x + a3*z) - vpz;
	long c = valuation(pr, 3*x2*x2 + b2*x2*xz + 3*b4*x2*z2 + 3*b6*xz*z2 + b8*z2*z2) -4*vpz;

// some obvious changes enable calculation of lambda as a rational
// some improvements can be made if this is never to be done
// eg in the above, no need to work with projective coords, just use real x/z
	bigfloat halfn = n; halfn /= 2;
	bigfloat lambda;

	if ((a <= 0) || (b <= 0)) {
		lambda = vpz - valuation(pr, x);
		if (lambda < 0)
			lambda = 0;
	}
	else if (!is_divisor(pr, c4)) {
		bigfloat m(b);
		if (halfn < m)
			m = halfn; // m = min(b , halfn);
		lambda = (m*(m-n)) / n;
	}
	else if (c >= (3*b))
		lambda = (-2*b) / bigfloat(3.0);
	else
		lambda = -c / bigfloat(4.0);

	bigfloat ans(lambda * log(bigfloat(pr)));
	return ans;
}



bigfloat
point< bigint >::get_realheight() const
{
	bigfloat xx(bigfloat(x)/bigfloat(z));
	return real_height(*ec, xx);
}



int
point< bigint >::get_order()
{
	if (ord) {
		return ord;
	}
	if (z.is_zero()) {
		ord = 1;
		return 1;
	}
	point< bigint > q = (*this);
	ord = 1;
	bigint eight(8);
	while (!q.is_zero() && (q.z <= eight)) {
		// (worst denom of a torsion point)
		ord++; add(q, q, (*this));
	}
	if (!q.is_zero())
		ord = -1;
	return ord;
}



int
order(point< bigint > & P, base_vector< point < bigint > > & list)
{
	list.set_size(0);
	list.set_capacity(12);
// NB Don't set P.ord until its correct value is known,
// otherwise Q's ord gets set wrongly
	list[0] = point< bigint > (*(P.ec));
	if (P.is_zero()) {
		P.ord = 1;
		return 1;
	}
	long ord = 2; list[1] = P;
	bigint eight(8);
	point< bigint > Q = P;
	while (!(Q.is_zero()) && (Q.z <= eight) && (ord < 13)) {
		add(Q, Q, P);
		if (!(Q.is_zero()))
			list[ord++] = Q;
	}

	if (!(Q.is_zero()))
		ord = -1;
	P.ord = ord;
	return ord;
}



//
// I/O
//

void
point< bigint >::read (std::istream & in)
{
	char c;
	in.flags(in.flags() | std::ios::dec); //force decimal input (bug fix)
	in >> c;
	switch(c) {
	case '(':
		in >> x >> c >> y >> c;
		z.assign_one();
		break;
	case '[':
		in >> x >> c >> y >> c >> z >> c;
		break;
	}
	ord = 0; height = -1.0;
	if (on_curve())
		reduce();
	else
		lidia_error_handler("point< bigint >", "read >>::not on curve");
}



void
point< bigint >::write (std::ostream & out) const
{
	out << "[" << x << " : " << y << " : " << z << "]";
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
