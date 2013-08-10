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


#ifndef LIDIA_POINT_OPERATIONS_H_GUARD_
#define LIDIA_POINT_OPERATIONS_H_GUARD_


#ifndef LIDIA_ELLIPTIC_CURVE_FLAGS_H_GUARD_
# include	"LiDIA/elliptic_curve_flags.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class T > class base_point;
class bigint;
class gf_element;


template< class T >
class point_operations
{
private:

	//
	// types of function pointers
	//
	typedef void (*generic_add_ptr)(base_point< T > &, const base_point< T > &,
					const base_point< T > &);
	typedef void (*generic_negate_ptr)(base_point< T > &, const base_point< T > &);
	typedef void (*generic_mult_by_2_ptr)(base_point< T > &,
					      const base_point< T > &);

	//
	// friends
	//

	friend class base_point< T >;

#if 0
	friend void multiply_by_2< T >(base_point< T > &,
				       const base_point< T > &);

	friend void negate< T >(base_point< T > &,
				const base_point< T > &);

	friend void add< T >(base_point< T > &,
			     const base_point< T > &,
			     const base_point< T > &);

	friend void subtract< T >(base_point< T > &,
				  const base_point< T > &,
				  const base_point< T > &);
#endif

private:

	generic_add_ptr       _add;
	generic_negate_ptr    _negate;
	generic_mult_by_2_ptr _mult_by_2;


	//
	// constructor / destructor
	//
public:

	point_operations ();


	~point_operations ();


	void add(base_point< T > & R,
		 const base_point< T > & P,
		 const base_point< T > & Q) const
	{
		_add(R, P, Q);
	}
	void negate(base_point< T > & R, const base_point< T > & P) const
	{
		_negate(R, P);
	}
	void mult_by_2(base_point< T > & R, const base_point< T > & P) const
	{
		_mult_by_2(R, P);
	}


	//
	// assigment
	//
public:
	void set(elliptic_curve_flags::curve_parametrization cp,
		 elliptic_curve_flags::curve_model cm,
		 const bigint& characteristic);

	point_operations< T > & operator = (const point_operations< T > & p);

};



template<>
class point_operations< gf_element >
{
	//
	// types of function pointers
	//
	typedef void (*generic_add_ptr)(base_point< gf_element > &,
					const base_point< gf_element > &,
					const base_point< gf_element > &);
	typedef void (*generic_negate_ptr)(base_point< gf_element > &,
					   const base_point< gf_element > &);
	typedef void (*generic_mult_by_2_ptr)(base_point< gf_element > &,
					      const base_point< gf_element > &);

	//
	// friends
	//
	friend class base_point< gf_element >;

#if 0
	template< class T > friend void multiply_by_2< T >(base_point< T > &,
				  const base_point< T > &);

	template< class T > friend void negate< T >(base_point< T > &,
			   const base_point< T > &);

	template< class T > friend void add< T >(base_point< T > &,
			const base_point< T > &,
			const base_point< T > &);

	template< class T > friend void subtract< T >(base_point< T > &,
			     const base_point< T > &,
			     const base_point< T > &);
#endif

private:

	generic_add_ptr       _add;
	generic_negate_ptr    _negate;
	generic_mult_by_2_ptr _mult_by_2;



	//
	// constructor / destructor
	//
public:

	point_operations();
	~point_operations();


	void add(base_point< gf_element > & R,
		 const base_point< gf_element > & P,
		 const base_point< gf_element > & Q) const
	{
		_add(R, P, Q);
	}
	void negate(base_point< gf_element > & R, const base_point< gf_element > & P) const
	{
		_negate(R, P);
	}
	void mult_by_2(base_point< gf_element > & R, const base_point< gf_element > & P) const
	{
		_mult_by_2(R, P);
	}

	//
	// assigment
	//
	void set(elliptic_curve_flags::curve_parametrization cp,
		 elliptic_curve_flags::curve_model cm,
		 const bigint & characteristic);


	point_operations< gf_element > & operator = (const point_operations< gf_element > & p);

};


//*********************************************************************
// R = P + Q,
//
// Condition P != +- Q not verified.
//

// Implemented in point_operations.cc
template< class T >
void add_swnf_affine(base_point< T > & R,
		     const base_point< T > & P,
		     const base_point< T > & Q);

// Implemented in point_operations.cc
template< class T >
void add_lwnf_affine(base_point< T > &,
		     const base_point< T > &,
		     const base_point< T > &);

// Implemented in point_operations_gf_element.cc
void add_gf2nf_affine(base_point< gf_element > &,
		      const base_point< gf_element > &,
		      const base_point< gf_element > &);

// Implemented in point_operations.cc
template< class T >
void add_swnf_projective(base_point< T > & R,
			 const base_point< T > & P,
			 const base_point< T > & Q);

// Implemented in point_operations.cc
template< class T >
void add_lwnf_projective(base_point< T > & R,
	       	         const base_point< T > & P,
			 const base_point< T > & Q);

// Implemented in point_operations_gf_element.cc
void add_gf2nf_projective(base_point< gf_element > &,
		          const base_point< gf_element > &,
		          const base_point< gf_element > &);


//***************************************************************
// R = -P
//

// Implemented in point_operations.cc
template< class T >
void negate_swnf_affine(base_point< T > &R,
		        const base_point< T > &P);

// Implemented in point_operations.cc
template< class T >
void negate_lwnf_affine(base_point< T > &,
		        const base_point< T > &);

// Implemented in point_operations.cc
template< class T >
void negate_swnf_projective(base_point< T > &R,
			    const base_point< T > &P);

// Implemented in point_operations.cc
template< class T >
void negate_lwnf_projective(base_point< T > &R,
			    const base_point< T > &P);

// Implemented in point_operations_gf_element.cc
void negate_gf2nf_affine(base_point< gf_element > &,
			 const base_point< gf_element > &);

// Implemented in point_operations_gf_element.cc
void negate_gf2nf_projective(base_point< gf_element > &R,
			     const base_point< gf_element > &P);

//***********************************************************************
// R = 2 * P
//

// Implemented in point_operations.cc
template< class T >
void mult_by_2_swnf_affine(base_point< T > &R,
			   const base_point< T > &P);

// Implemented in point_operations.cc
template< class T >
void mult_by_2_lwnf_affine(base_point< T > &,
			   const base_point< T > &);

// Implemented in point_operations_gf_element.cc
void mult_by_2_gf2nf_affine(base_point< gf_element > &,
			    const base_point< gf_element > &);

// Implemented in point_operations.cc
template< class T >
void mult_by_2_swnf_projective(base_point< T > &R,
			       const base_point< T > &P);

// Implemented in point_operations.cc
template< class T >
void mult_by_2_lwnf_projective(base_point< T > &,
			       const base_point< T > &);

// Implemented in point_operations_gf_element.cc
void mult_by_2_gf2nf_projective(base_point< gf_element > &,
			        const base_point< gf_element > &);



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_POINT_OPERATIONS_H_GUARD_
