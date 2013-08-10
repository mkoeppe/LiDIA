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


#ifndef LIDIA_CURVE_ISOMORPHISM_CC_GUARD_
#define LIDIA_CURVE_ISOMORPHISM_CC_GUARD_


#ifndef LIDIA_CURVE_ISOMORPHISM_H_GUARD_
# include	"LiDIA/curve_isomorphism.h"
#endif
#ifndef LIDIA_BIGRATIONAL_H_GUARD_
# include	"LiDIA/bigrational.h"
#endif
#ifndef LIDIA_POINT_H_GUARD_
# include	"LiDIA/point.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



template< class S, class T >
bool
find_urst (const elliptic_curve< S > & e1, const elliptic_curve< T > & e2, S& u, S& r, S& s, S& t)
{
	S a1, a2, a3, a4, a6, c4, c6;
	e1.get_ai(a1, a2, a3, a4, a6);
	e1.get_ci(c4, c6);

	T a1d, a2d, a3d, a4d, a6d, c4d, c6d;
	e2.get_ai(a1d, a2d, a3d, a4d, a6d);
	e2.get_ci(c4d, c6d);

	S a1b(a1d), a2b(a2d), a3b(a3d), a4b(a4d), a6b(a6d), c4b(c4d), c6b(c6d);

	// NB We assume (as we may) that both curves are nonsingular
	// The following is then 0 IFF the two j-invariants are equal

	S ee = (c4*c4*c4)*(c6b*c6b) - (c4b*c4b*c4b)*(c6*c6);
	if (!ee.is_zero()) return false;

	// Now the curves are isomorphic over the algebraic closure, but may only
	// be twists over the ground field.  Special treatment is needed if
	// c6 == 0 (quartic twists possible) or
	// c4 == 0 (sextic twists possible);
	// otherwise only quadratic twists are possible.

	// We have u^2 = (c6*c4d)/(c4*c6d) in general,
	//         u^4 = c4/c4d            if c6 = c6d = 0,
	//         u^6 = c6/c6d            if c4 = c4d = 0.


	S u2, u3, u4, u6;

	if (c4b.is_zero()) {
		if (!c4.is_zero()) {
			return false;
		}
		if (c6b.is_zero()) {
			lidia_error_handler("curve_isomorphism", "find_urst::Second curve is singular");
		}
		u6 = c6/c6b;
		if (!square_root(u3, u6)) {
			return false;
		}
		if (!cube_root(u2, u6)) {
			return false;
		}
		{
			u = u3/u2;
		}
	}
	else
		if (c6b.is_zero()) {
			if (!c6.is_zero()) {
				return false;
			}
			if (c4b.is_zero()) {
				lidia_error_handler("curve_isomorphism", "find_urst::Second curve is singular");
			}
			u4 = c4/c4b;
			if (!square_root(u2, u4)) {
				return false;
			}
			if (!square_root(u, u2))
				if (!square_root(u, -u2))
					return false;

		}
		else { // simple quadratic case
			u2 = (c6*c4b)/(c4*c6b);
			if (!square_root(u, u2)) {
				return false;
			}
		}

	// Now they are isomorphic via u and we find r, s, t

	u2 = u*u; u3 = u2*u;
	s = (u*a1b-a1)/2;
	r = (u2*a2b-a2+s*a1+s*s)/3;
	t = (u3*a3b-a3-r*a1)/2;
	return true;
}



template< class S, class T >
curve_isomorphism< S, T >::curve_isomorphism (elliptic_curve< S > & e1, elliptic_curve< T > & e2)
	:ec1(&e1), ec2(&e2)
{
	if (!find_urst(e1, e2, this->u, this->r, this->s, this->t)) {
		lidia_error_handler("curve_isomorphism",
				    "curve_isomorphism::Curves are not isomorphic");
	}
}



template< class S, class T >
bool
are_isomorphic (const elliptic_curve< S > & e1, const elliptic_curve< T > & e2)
{
	S u, r, s, t;
	return find_urst(e1, e2, u, r, s, t);
}



template< class S, class T >
bool
are_isomorphic (const elliptic_curve< S >& e1, const elliptic_curve< T >& e2,
		curve_isomorphism< S, T > & iso)
{
	S u, r, s, t;
	if (!find_urst(e1, e2, u, r, s, t)) return false;
	iso = curve_isomorphism< S, T > (&e1, &e2, u, r, s, t);
	return true;
}



template< class S, class T >
point< S >
curve_isomorphism< S, T >::inverse (const point< T > & p)
{
	if (p.get_curve() != *(this->ec2)) {
		lidia_error_handler("curve_isomorphism< S, T >",
				    "inverse::Point not on correct curve");
	}
	if (p.is_zero()) {
		return point< S > (*this->ec1);
	}

	S x(p.get_x()), y(p.get_y()), newx, newy;

	S u2 = this->u * this->u;
	newx = u2 * x + this->r;
	newy = this->u * u2 * y + u2 * this->s * x + this->t;

	point< S > pnew(newx, newy, *this->ec1);

	return pnew;
}



template< class S, class T >
point< T >
curve_isomorphism< S, T >::map (const point< S > & p)
{
	if (p.get_curve() != *this->ec1) {
		lidia_error_handler("curve_isomorphism< S, T >",
				    "map::Point not on correct curve");
	}
	if (p.is_zero()) {
		return point< T > (*this->ec2);
	}

	S x = p.get_x(), y = p.get_y(), newx, newy;
	S u2 = this->u * this->u;
	newx = (x - this->r) / u2;
	newy = (y - (x - this->r) * this->s - this->t)/ (this->u * u2);

	T nx, ny;
	convert(nx, newx);
	convert(ny, newy);

	point< T > pnew(nx, ny, *this->ec2);

	return pnew;
}


template <>
point< bigrational >
curve_isomorphism< bigrational, bigint >::inverse (const point< bigint > &p)
{
	if (p.get_curve() != *this->ec2) {
		lidia_error_handler("curve_isomorphism< bigrational, bigint >",
				    "inverse::Point not on correct curve");
	}
	if (p.is_zero()) {
		return point< bigrational > (*this->ec1);
	}

	bigint xx, yy, zz;
	p.get_xyz(xx, yy, zz);

	bigrational x(xx, zz), y(yy, zz), newx, newy;
	bigrational u2 = this->u * this->u;
	newx = u2 * x + this->r;
	newy = this->u * u2 * y + u2 * this->s * x + this->t;

	point< bigrational > pnew(newx, newy, *this->ec1);

	return pnew;
}



template <>
point< bigint >
curve_isomorphism< bigrational, bigint >::map (const point< bigrational > & p)
{
	if (p.get_curve() != *this->ec1) {
		lidia_error_handler("curve_isomorphism< bigrational, bigint >",
				    "map::Point not on correct curve");
	}
	if (p.is_zero()) {
		return point< bigint > (*this->ec2);
	}

	bigrational x = p.get_x(), y = p.get_y(), newx, newy;
	newx = (x - this->r) / (this->u * this->u);
	newy = (y - (x - this->r) * this->s - this->t) / (this->u * this->u * this->u);

	bigint nx, ny, nz;
	nz = lcm(newx.denominator(), newy.denominator());
	nx = (nz*newx).numerator();
	ny = (nz*newy).numerator();

	point< bigint > pnew(nx, ny, nz, *this->ec2);

	return pnew;
}



template< class S, class T >
elliptic_curve< S >
make_isomorphic_curve (const elliptic_curve< T > & ee,
		       const S& u, const S& r, const S&s, const S& t)
{
	T t1, t2, t3, t4, t6;
	ee.get_ai(t1, t2, t3, t4, t6);

	S u2 = u * u, u3 = u2 * u, u4 = u2 * u2, u6 = u3 * u3;

	S a1 = (t1 + 2 * s) / u;
	S a2 = (t2 - s * t1 + 3 * r - s * s) / u2;
	S a3 = (t3 + r * t1 + 2 * t) / u3;
	S a4 = (t4 - s * t3 + 2 * r * t2 - (t + r * s) * t1 +
		3 * r * r - 2 * s * t) / u4;
	S a6 = (t6 + r * t4 + r * r * (t2 + r) - t * (t3 + t + r * t1)) / u6;

	elliptic_curve< S > e(a1, a2, a3, a4, a6);

	return e;
}



template< class S, class T >
void
curve_isomorphism< S, T >::write (std::ostream & out) const
{
	out << "[" << this->u << ", " << this->r << ", " << this->s << ", " << this->t << "]";
}



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_CURVE_ISOMORPHISM_CC_GUARD_
