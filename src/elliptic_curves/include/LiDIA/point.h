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
//	Author	: Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_POINT_H_GUARD_
#define LIDIA_POINT_H_GUARD_


#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_GF_ELEMENT_H_GUARD_
# include	"LiDIA/gf_element.h"
#endif
#ifndef LIDIA_BASE_POINT_H_GUARD_
# include	"LiDIA/elliptic_curves/base_point.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class T > class elliptic_curve;


template< class T > class point : public base_point< T >
{
	//
	// constructors / destructor
	//
public:
	point();
	point(const base_point< T > &);
	point(const T & xp, const T & yp, const elliptic_curve< T > & e);
	point(const T & xp, const T & yp, const T & zp, const elliptic_curve< T > & e);
	point(const point< T > & P);
	point(const elliptic_curve< T > &e);
	~point();

	//
	// assignments
	//
	point< T > & operator = (const point< T > & P);
};



template<>
class point< gf_element > : public base_point< gf_element >
{
	//
	// constructors / destructor
	//
public:
	point();
	point(const base_point< gf_element > &);
	point(const gf_element & xp, const gf_element & yp, const elliptic_curve< gf_element > & e);
	point(const gf_element & xp, const gf_element & yp, const gf_element & zp, const elliptic_curve< gf_element > & e);
	point(const point< gf_element > & P);
	point(const elliptic_curve< gf_element > &e);
	~point();

	//
	// assignments
	//
	point< gf_element > & operator = (const point< gf_element > & P);
	bool set_x_coordinate(const gf_element & xx);

	//
	// Frobenius map
	//
	friend point< gf_element > frobenius(const point< gf_element > &, unsigned int d);
	void frobenius(unsigned int d = 1);

	//
	// point order
	//
	friend bigint order_point(const point< gf_element > & P);
	friend bigint order_point(const point< gf_element > & P, rational_factorization &);

	friend rational_factorization rf_order_point(const point< gf_element > & P);
	friend rational_factorization rf_order_point(const point< gf_element > & P,
						     rational_factorization &);



	//
	// discrete logarithm
	//
	friend bigint bg_algorithm(const point< gf_element > & P,
				   const point< gf_element > & Q,
				   const bigint & lower, const bigint & upper,
                                   bool info);

	friend bigint bg_algorithm(const point< gf_element > & P,
				   const point< gf_element > & Q,
				   const bigint & lower, const bigint & upper,
				   const bigint & x, const bigint & m,
                                   bool info);

private:
	bigint pollard_rho(point< gf_element > &,
			   const bigint &,
			   bool info = false);
public:

	friend bigint discrete_logarithm (const point< gf_element > &,
					  const point< gf_element > &,
					  bool info);

        friend bigint discrete_logarithm (const point< gf_element > &,
                                          const rational_factorization &,
                                          const point< gf_element > &,
                                          const bigint &,
                                          bool info);
};

// friend functions

//
// Frobenius map
//
point< gf_element > frobenius(const point< gf_element > &, unsigned int d = 1);

//
// point order
//
bigint order_point(const point< gf_element > & P);
bigint order_point(const point< gf_element > & P, rational_factorization &);

rational_factorization rf_order_point(const point< gf_element > & P);
rational_factorization rf_order_point(const point< gf_element > & P,
				      rational_factorization &);

//
// discrete logarithm
//
bigint bg_algorithm(const point< gf_element > & P,
		    const point< gf_element > & Q,
		    const bigint & lower, const bigint & upper,
		    bool info = false);

bigint bg_algorithm(const point< gf_element > & P,
		    const point< gf_element > & Q,
		    const bigint & lower, const bigint & upper,
		    const bigint & x, const bigint & m,
		    bool info = false);

bigint discrete_logarithm (const point< gf_element > &,
			   const point< gf_element > &,
			   bool info = false);

bigint discrete_logarithm (const point< gf_element > &,
			   const rational_factorization &,
			   const point< gf_element > &,
			   const bigint &,
			   bool info = false);



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/point.cc"
#endif



#endif	// LIDIA_POINT_H_GUARD_
