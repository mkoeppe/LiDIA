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
//	Author	: Nigel Smart, John Cremona
//		  Adaption of John Cremona's code; some code came
//                originally from Oisin McGuiness
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_POINT_BIGINT_H_GUARD_
#define LIDIA_POINT_BIGINT_H_GUARD_


#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_BIGFLOAT_H_GUARD_
# include	"LiDIA/bigfloat.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class T > class point;
template< class T > class elliptic_curve;
template< class T > class base_vector;



template<>
class point< bigint >
{
	friend void add_swnf_projective(point< bigint > &, const point< bigint > &,
					const point< bigint > &);
	friend void add_lwnf_projective(point< bigint > &, const point< bigint > &,
					const point< bigint > &);

	friend void negate_swnf_projective(point< bigint > &, const point< bigint > &);
	friend void negate_lwnf_projective(point< bigint > &, const point< bigint > &);

	friend void mult_by_2_swnf_projective(point< bigint > &, const point< bigint > &);
	friend void mult_by_2_lwnf_projective(point< bigint > &, const point< bigint > &);

	bigint x, y, z;
	elliptic_curve< bigint > *ec; // Must be a pointer to avoid include in both directions,
	// point_bigint.h and elliptic_curve < bigint > .
	bigfloat height; // -1.0 if not calculated yet, 0.0 for torsion points
	int ord; // 0 if not calculated yet, -1 if infinite

	void reduce();

public:

	//
	// constructors / destructor
	//
	point();
	point(const bigint & xp, const bigint & yp, const elliptic_curve< bigint > &ecp);
	point(const bigint & xp, const bigint& yp, const bigint& zp, const elliptic_curve< bigint > &ecp);
	point(const point< bigint > & P);
	point(const elliptic_curve< bigint > &e);
	~point();

	//
	// assignment
	//
	point< bigint > & operator = (const point< bigint > & Q);
	void assign_zero(const elliptic_curve< bigint > & e);
	void assign(const point< bigint > & P);
	void assign(const bigint & xx, const bigint & yy, const elliptic_curve< bigint > & e);
	void assign(const bigint & xx, const bigint & yy, const bigint &zz, const elliptic_curve< bigint > & e);

	// assign on last curve (use with care of course)
	void assign_zero();
	void assign(const bigint & xx, const bigint & yy);
	void assign(const bigint & xx, const bigint & yy, const bigint &zz);

	void swap(point< bigint > & P);


	//
	// Testing points
	//
	bool on_curve() const;
	bool is_equal(const point< bigint > &P) const;

	bool is_zero() const;
	bool is_integral() const;
	bool is_negative_of(const point< bigint > & P);

	//
	// Arithmetic
	//
	friend void negate(point< bigint > &, const point< bigint > &);
	friend void add(point< bigint > &, const point< bigint > &, const point< bigint > &);
	friend void subtract(point< bigint > &, const point< bigint > &, const point< bigint > &);
	friend void multiply_by_2(point< bigint > &, const point< bigint > &);
	friend void multiply(point< bigint > &, const bigint&, const point< bigint > &);

	point< bigint > twice() const;


	//
	// I/O
	//
	void write(std::ostream & out) const;
	void read (std::istream & is);


	//
	// Access
	//
	void get_xyz(bigint& xx, bigint& yy, bigint& zz) const
	{
		xx = x;
		yy = y;
		zz = z;
	}
	elliptic_curve< bigint > get_curve() const;

	//
	// Height Functions
	//
	bigfloat get_height();
	bigfloat get_naive_height() const;
	bigfloat get_pheight(const bigint& p) const;
	bigfloat get_realheight() const;

	//
	// Order Functions
	//
	int get_order();
	friend int order(point< bigint > & P, base_vector< point < bigint > > & list);
	// also create and return list of multiples
};



//
// assigners
//

inline point< bigint > &
point< bigint >::operator = (const point< bigint > & P)
{
	assign(P);
	return *this;
}



//
// comparators
//

inline bool
operator == (const point< bigint > & P, const point< bigint > & Q)
{
	return P.is_equal(Q);
}



inline bool
operator != (const point< bigint > & P, const point< bigint > & Q)
{
	return !P.is_equal(Q);
}



inline bool
point< bigint >::is_zero() const
{
	return z.is_zero();
}



inline bool
point< bigint >::is_integral() const
{
	return z.is_one();
}



inline bool
point< bigint >::is_negative_of(const point< bigint > & P)
{
	point< bigint > H;
	negate(H, P);
	return (*this == H);
}



inline bool
is_equal(const point< bigint > & P, const point< bigint > & Q)
{
	return P.is_equal(Q);
}



//
// Operators
//

inline point< bigint > &
operator += (point< bigint > & P, const point< bigint > & Q)
{
	add(P, P, Q);
	return P;
}



inline point< bigint > &
operator -= (point< bigint > & P, const point< bigint > & Q)
{
	subtract(P, P, Q);
	return P;
}



inline point< bigint >
operator + (const point< bigint > & P, const point< bigint > & Q)
{
	point< bigint > R;

	add(R, P, Q);
	return R;
}



inline point< bigint >
operator - (const point< bigint > & P, const point< bigint > & Q)
{
	point< bigint > R;

	subtract(R, P, Q);
	return R;
}



inline point< bigint >
operator - (const point< bigint > & P)
{
	point< bigint > R;

	negate(R, P);
	return R;
}



inline point< bigint >
point< bigint >::twice() const
{
	point< bigint > R;

	multiply_by_2(R, (*this));
	return R;
}



inline point< bigint >
operator * (bigint n, const point< bigint > & P)
{
	point< bigint > R;

	multiply(R, n, P);
	return R;
}



inline void
swap (point< bigint > & P, point< bigint > & Q)
{
	P.swap(Q);
}



inline std::istream &
operator >> (std::istream & in, point< bigint > & P)
{
	P.read(in);
	return in;
}



inline std::ostream &
operator << (std::ostream & out, const point< bigint > & P)
{
	P.write(out);
	return out;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_POINT_BIGINT_H_GUARD_
