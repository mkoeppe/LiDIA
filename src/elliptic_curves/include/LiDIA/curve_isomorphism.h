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
//	Author	: Nigel Smar (NS), John Cremona (JC)
//	Changes	: See CVS log
//
//==============================================================================================

// The idea of this class is so we can keep track of
// curve isomorphisms and check whether too curves are
// isomorphic.
// Note isomorphisms are defined over S, not algebraic closure of S.
// Basic data of an isomorphism: [u, r, s, t] with u non-zero
// Functions: given two curves, test for isomorphism, return [u, r, s, t] if so;
//            given a curve and [u, r, s, t] make the isomorphic curve
//            given a point and [u, r, s, t] make the transformed point


#ifndef LIDIA_CURVE_ISOMORPHISM_H_GUARD_
#define LIDIA_CURVE_ISOMORPHISM_H_GUARD_


#ifndef LIDIA_ELLIPTIC_CURVE_BIGINT_H_GUARD_
# include	"LiDIA/elliptic_curve_bigint.h"
#endif
#ifndef LIDIA_POINT_BIGINT_H_GUARD_
# include	"LiDIA/point_bigint.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



// Here we assume that S has a constructor for entries of type T
//	eg S = bigrational
//	   T = bigint or bigrational
// Class S must also support EXACT division
// and functions should exist for S of the form
//	int square_root(const S& n, S& root);
//	int cube_root(const S& n, S& root);
// which test if n is a square (resp. cube) and put a root in root if so.
//	int is_zero();
// Assumed char S <> 2
// The following function should exist called
//	void convert(T& a, const S& b)
// (Unless you write your own "map" member below)
// This checks whether we can convert b to be of type S, if not error
//	then a becomes the S version of b
// Just to make things easier when S = T we define an obvious one below


template< class S, class T >
bool find_urst(const elliptic_curve< S > & e1, const elliptic_curve< T > & e2, S& u, S& r, S& s, S& t);

template< class T > inline void convert(T& a, const T& b)
{
	a = b;
}



template< class S, class T > class curve_isomorphism
{
protected:

	S u, r, s, t; // The change of variables
	elliptic_curve< S > *ec1;
	elliptic_curve< T > *ec2;

public:

	//
	// constructors / destructor
	//

	curve_isomorphism() :u(0), r(0), s(0), t(0), ec1(0), ec2(0) {}
	// This constructor will give an error if they are not isomorphic

	curve_isomorphism(elliptic_curve< S > &, elliptic_curve< T > &);
	// This constructor is used when the isomorphism is already known
	// (it does not check that is is valid)

	curve_isomorphism(elliptic_curve< S > & e1, elliptic_curve< T > & e2,
			  const S& uu, const S& rr, const S& ss, const S& tt)
		:u(uu), r(rr), s(ss), t(tt), ec1(&e1), ec2(&e2)
	{}

	~curve_isomorphism() {}

	curve_isomorphism(const curve_isomorphism< S, T > & ci)
		:u(ci.u), r(ci.r), s(ci.s), t(ci.t), ec1(ci.ec1), ec2(ci.ec2)
	{}

	curve_isomorphism< S, T > & operator = (const curve_isomorphism< S, T > & ci)
	{
		this->u = ci.u;
		this->r = ci.r;
		this->s = ci.s;
		this->t = ci.t;
		this->ec1 = ci.ec1;
		this->ec2 = ci.ec2;
		return *this;
	}


	//
	// initialization
	//

	void init(elliptic_curve< S > & e1, elliptic_curve< T > & e2,
		  const S& uu, const S& rr, const S& ss, const S& tt)
	{
		this->u = uu;
		this->r = rr;
		this->s = ss;
		this->t = tt;
		this->ec1 = &e1;
		this->ec2 = &e2;
	}

	//
	//
	//

	S get_scale_factor() const {
		return this->u;
	}

	bool is_unimodular() const
	{
		if (this->u == 1)
			return true;
		else
			return false;
	}
	bool is_identity() const
	{
		if ((this->u == 1) && this->r.is_zero() &&
		    this->s.is_zero() && this->t.is_zero())
			return true;
		else
			return false;
	}

	point< S > inverse(const point< T > & p);
	point< T > map(const point< S > & p);


	//
	// output
	//

	void write(std::ostream &out) const;

};



template < class S, class T >
inline std::ostream &
operator << (std::ostream & out, const curve_isomorphism< S, T > & is)
{
	is.write(out);
	return out;
}



// Tests for isomorphism: second one returns the isomorphism via iso parameter

template< class S, class T >
bool are_isomorphic(const elliptic_curve< S > &, const elliptic_curve< T > &);

template< class S, class T >
bool are_isomorphic(const elliptic_curve< S > &, const elliptic_curve< T > &,
		    curve_isomorphism< S, T > & iso);

// Creates an isomorphic curve using given u, r, s, t

template< class S, class T >
elliptic_curve< S > make_isomorphic_curve(const elliptic_curve< T > &,
					  const S& u, const S& r, const S&s, const S& t);





#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/curve_isomorphism.cc"
#endif



#endif	// LIDIA_CURVE_ISOMORPHISM_H_GUARD_
