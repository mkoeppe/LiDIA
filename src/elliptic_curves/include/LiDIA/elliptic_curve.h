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
//	Author	: Markus Maurer (MM), Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_ELLIPTIC_CURVE_H_GUARD_
#define LIDIA_ELLIPTIC_CURVE_H_GUARD_



#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_GF_ELEMENT_H_GUARD_
# include	"LiDIA/gf_element.h"
#endif
#ifndef LIDIA_BASE_ELLIPTIC_CURVE_H_GUARD_
# include	"LiDIA/elliptic_curves/base_elliptic_curve.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class T > class point;



//
// general class
//

template< class T >
class elliptic_curve : public base_elliptic_curve< T >
{
public:

	//
	// constructor / destructor
	//
	elliptic_curve();
	elliptic_curve(const elliptic_curve< T > &);
	elliptic_curve(const T & x4, const T & x6,
		       elliptic_curve_flags::curve_model m =
		       elliptic_curve_flags::AFFINE);

	elliptic_curve(const T & x1, const T & x2,
		       const T & x3, const T & x4, const T & x6,
		       elliptic_curve_flags::curve_model m =
		       elliptic_curve_flags::AFFINE);

	~elliptic_curve();

	//
	// assignment
	//
	elliptic_curve< T > & operator = (const elliptic_curve< T > &);
	void assign(const elliptic_curve< T > &);

	//
	// invariants
	//
	T j_invariant() const;
};



//
// specialization for gf_element
//

template<>
class elliptic_curve< gf_element > : public base_elliptic_curve< gf_element >
{
private:
	rational_factorization order_of_group; // over extension field

public:

	//
	// constructor / destructor
	//
	elliptic_curve();
	elliptic_curve(const elliptic_curve< gf_element > &);
	elliptic_curve(const gf_element & x4, const gf_element & x6,
		       elliptic_curve_flags::curve_model m =
		       elliptic_curve_flags::AFFINE);

	elliptic_curve(const gf_element & x1, const gf_element & x2,
		       const gf_element & x3, const gf_element & x4,
		       const gf_element & x6,
		       elliptic_curve_flags::curve_model m =
		       elliptic_curve_flags::AFFINE);

	~elliptic_curve();

	//
	// assignment
	//
	elliptic_curve< gf_element > & operator = (const elliptic_curve< gf_element > &);
	void assign(const elliptic_curve< gf_element > &);

	//
	// invariants
	//
	gf_element j_invariant() const;

	//
	// properties
	//

	bool is_supersingular();
	unsigned int degree_of_definition() const;
	bool is_cyclic();
	void isomorphism_type(bigint & n1, bigint & n2, bool info = false);
	void isomorphism_type(bigint & n1, bigint & n2,
			      point< gf_element > & P1,
			      point< gf_element > & P2, bool info = false);

	//
	// random points
	//
	point< gf_element > random_point(unsigned int d = 0) const;

	//
	// group order
	//

private:

	bigint small_counting_points_gf2(unsigned int d) const;
	bigint small_counting_points_gfp(unsigned int d) const;
        bigint small_counting_points(unsigned int d) const
          {
	    if (this->get_a6().get_field().characteristic() == 2)
	      return this->small_counting_points_gf2(d);
	    else
	      return this->small_counting_points_gfp(d);
          }


public:

	bigint group_order(int info = 0);
	rational_factorization rf_group_order(bool prime_factorization = false,
					      int info = 0);

	void set_group_order(const rational_factorization & o);
        void set_group_order(const bigint & x)
           {
             set_group_order(rational_factorization(x));
           }

	bool probabilistic_test_of_group_order(const bigint & res,
					       lidia_size_t tests = 5,
					       bool use_internal_info = true);
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/elliptic_curve.cc"
#endif



#endif	// LIDIA_ELLIPTIC_CURVE_H_GUARD_
