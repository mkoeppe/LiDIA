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


#ifndef LIDIA_BASE_POINT_H_GUARD_
#define LIDIA_BASE_POINT_H_GUARD_


#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_ELLIPTIC_CURVE_H_GUARD_
# include	"LiDIA/elliptic_curve.h"
#endif
#ifndef LIDIA_POINT_OPERATIONS_H_GUARD_
# include	"LiDIA/elliptic_curves/point_operations.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class gf_element;

template < class T >
class base_point;

// friend functions of class template base_point

//
// arithmetic
//

template < class T >
void negate(base_point< T > & r, const base_point< T > & x);


template < class T >
void add(base_point< T > & R, const base_point< T > & P1,
	      const base_point< T > & Q);
	

template < class T >
void subtract(base_point< T > & r, const base_point< T > & x,
		       const base_point< T > & y);


template < class T >
void multiply_by_2(base_point< T > & R, const base_point< T > & P);


template < class T >
void multiply(base_point< T > & , const bigint &,
		       const base_point< T > &);

//
// Functions for arithmetic that depend on the curve model
// and that are implemented in point_operations.cc
//


template < class T >
void add_swnf_affine(base_point< T > &, const base_point< T > &,
			      const base_point< T > &);


template < class T >
void add_swnf_projective(base_point< T > &,
				  const base_point< T > &,
				  const base_point< T > &);


template < class T >
void add_lwnf_affine(base_point< T > &, const base_point< T > &,
			      const base_point< T > &);


template < class T >
void add_lwnf_projective(base_point< T > &,
				  const base_point< T > &,
				  const base_point< T > &);


void add_gf2nf_affine(base_point< gf_element > &,
			  const base_point< gf_element > &,
			  const base_point< gf_element > &);


void add_gf2nf_projective(base_point< gf_element > &,
			      const base_point< gf_element > &,
			      const base_point< gf_element > &);


template < class T >
void negate_swnf_affine(base_point< T > &,
				 const base_point< T > &);


template < class T >
void negate_swnf_projective(base_point< T > &,
				     const base_point< T > &);


template < class T >
void negate_lwnf_affine(base_point< T > &, const base_point< T > &);


template < class T >
void negate_lwnf_projective(base_point< T > &,
				     const base_point< T > &);


void negate_gf2nf_affine(base_point< gf_element > &,
			     const base_point< gf_element > &);


void negate_gf2nf_projective(base_point< gf_element > &,
				 const base_point< gf_element > &);


template < class T >
void mult_by_2_swnf_affine(base_point< T > &,
				    const base_point< T > &);


template < class T >
void mult_by_2_swnf_projective(base_point< T > &,
					const base_point< T > &);


template < class T >
void mult_by_2_lwnf_affine(base_point< T > &,
				    const base_point< T > &);


template < class T >
void mult_by_2_lwnf_projective(base_point< T > &,
					const base_point< T > &);


void mult_by_2_gf2nf_affine(base_point< gf_element > &,
				const base_point< gf_element > &);


void mult_by_2_gf2nf_projective(base_point< gf_element > &,
				    const base_point< gf_element > &);



template< class T >
class base_point
{
protected:

	elliptic_curve< T > ec;
	T x, y, z;
	bool is_0;
	int info;



public:

	//
	// constructors / destructor
	//

	base_point();

	base_point(const T & xp, const T & yp,
		   const elliptic_curve< T > & e);


	base_point(const T & xp, const T & yp, const T & zp,
		   const elliptic_curve< T > & e);

	base_point(const base_point< T > & P);

	base_point(const elliptic_curve< T > &e);

	~base_point();


	//
	// assignments
	//
	void set_verbose(int a);

private:

	void assign(const elliptic_curve< T > & e);


public:

	void assign_zero(const elliptic_curve< T > & e);
	void assign(const base_point< T > & P);

	void assign(const T & xx, const T & yy,
		    const elliptic_curve< T > & e);
	void assign(const T & xx, const T & yy,
		    const T & zz, const elliptic_curve< T > & e);

	// for efficiency reasons
	void assign_no_test(const T & xx, const T & yy,
			    const elliptic_curve< T > & e);
	void assign_no_test(const T & xx, const T & yy, const T & zz,
			    const elliptic_curve< T > & e);

	void assign_zero();
	void assign(const T & xx, const T & yy);
	void assign(const T & xx, const T & yy, const T & zz);

	base_point< T > & operator = (const base_point< T > & P);



	//
	// properties
	//
	bool on_curve() const;
	bool is_zero() const;
	bool is_negative_of(const base_point< T > & P) const;

        bool is_affine_point() const;
        bool is_projective_point() const;

	//
	// comparison
	//

	bool is_equal (const base_point< T > & P) const;


	//
	// accessors
	//

	const T & get_x() const;
	const T & get_y() const;
	const T & get_z() const;
	const elliptic_curve< T > & get_curve() const;



	void make_affine_point(bool change_curve = true);
	void make_projective_point(bool change_curve = true);


	//
	// arithmetic
	//
        friend void negate< T >(base_point & r, const base_point & x);

	friend
	void add< T >(base_point & R, const base_point & P1,
			     const base_point & Q);
	
        friend
        void subtract< T >(base_point & r, const base_point & x,
				  const base_point & y);

	friend
        void multiply_by_2< T >(base_point & R, const base_point & P);

	friend
	void multiply< T >(base_point & , const bigint &,
				  const base_point &);

	base_point twice() const;

	void swap(base_point & P);



	//
	// input / output
	//

	// NB It is currently not possible to input the base_point at infinity
	void read(std::istream & in);
	void write(std::ostream & out) const;



	//
	// Functions for arithmetic that depend on the curve model
	// and that are implemented in point_operations.cc
	//

	friend
	void add_swnf_affine< T >(base_point &, const base_point &,
					 const base_point &);

	friend
	void add_swnf_projective< T >(base_point &,
					     const base_point &,
					     const base_point &);

	friend
	void add_lwnf_affine< T >(base_point &, const base_point &,
					 const base_point &);

	friend
	void add_lwnf_projective< T >(base_point &,
					     const base_point &,
					     const base_point &);

	friend
	void add_gf2nf_affine(base_point< gf_element > &,
				     const base_point< gf_element > &,
				     const base_point< gf_element > &);

	friend
	void add_gf2nf_projective(base_point< gf_element > &,
					 const base_point< gf_element > &,
					 const base_point< gf_element > &);

	friend
	void negate_swnf_affine< T >(base_point &,
					    const base_point &);

	friend
	void negate_swnf_projective< T >(base_point &,
						const base_point &);

	friend
	void negate_lwnf_affine< T >(base_point &, const base_point &);

	friend
	void negate_lwnf_projective< T >(base_point &,
						const base_point &);

	friend
	void negate_gf2nf_affine(base_point< gf_element > &,
					const base_point< gf_element > &);

	friend
	void negate_gf2nf_projective(base_point< gf_element > &,
					      const base_point< gf_element > &);

	friend
	void mult_by_2_swnf_affine< T >(base_point &,
					       const base_point &);

	friend
	void mult_by_2_swnf_projective< T >(base_point &,
						   const base_point &);

	friend
	void mult_by_2_lwnf_affine< T >(base_point &,
					       const base_point &);

	friend
	void mult_by_2_lwnf_projective< T >(base_point &,
						   const base_point &);

	friend
	void mult_by_2_gf2nf_affine(base_point< gf_element > &,
					     const base_point< gf_element > &);

	friend
	void mult_by_2_gf2nf_projective(base_point< gf_element > &,
					       const base_point< gf_element > &);

};




//
// c'tors and d'tor
//

template< class T >
inline
base_point< T >::base_point ()
{
	debug_handler("base_point< T >", "base_point()");
	this->info = 0;
}



template< class T >
inline
base_point< T >::base_point (const T & xp, const T & yp,
			     const elliptic_curve< T > & e)
{
	debug_handler("base_point< T >",
		      "base_point(constT&, constT&, const elliptic_curve< T > &)");

	assign(xp, yp, e);
}



template< class T >
inline
base_point< T >::base_point (const T & xp, const T & yp, const T & zp,
			     const elliptic_curve< T > & e)
{
	debug_handler("base_point< T >",
		      "base_point(const T &, const T &, const T &, "
		      "const elliptic_curve< T > &)");

	assign(xp, yp, zp, e);
}



template< class T >
inline
base_point< T >::base_point (const base_point< T > & P)
{
	debug_handler("base_point< T >", "base_point(const base_point< T > &");

	assign(P);
}



template< class T >
inline
base_point< T >::base_point (const elliptic_curve< T > & e)
{
	debug_handler("base_point< T >", "base_point(const elliptic_curve< T > &");
	assign_zero(e);
}



template< class T >
inline
base_point< T >::~base_point ()
{
	debug_handler("base_point< T >", "~base_point()");
}



template< class T >
inline const T &
base_point< T >::get_x () const
{
	if (is_zero())
		lidia_error_handler("base_point< T >::get_x()", "Point at infinity.");
	return this->x;
}



template< class T >
inline const T &
base_point< T >::get_y () const
{
	if (is_zero())
		lidia_error_handler("base_point< T >::get_y()", "Point at infinity.");
	return this->y;
}



template< class T >
inline const T &
base_point< T >::get_z () const
{
	if (is_zero())
		lidia_error_handler("base_point< T >::get_z()", "Point at infinity.");

	if (this->ec.get_model() == elliptic_curve_flags::AFFINE)
		lidia_error_handler("base_point< T >::get_z()", "Affine model.");

	return this->z;
}



template< class T >
inline const elliptic_curve< T > &
base_point< T >::get_curve () const
{
	return this->ec;
}



template< class T >
inline base_point< T > &
base_point< T >::operator = (const base_point< T > & P)
{
	debug_handler("base_point< T >", " operator = const base_point< T > &");
	assign(P);
	return *this;
}



template< class T >
inline bool
base_point< T >::is_zero () const
{
	if (this->ec.get_model() == elliptic_curve_flags::AFFINE)
		return this->is_0;
	else
		return this->z.is_zero();
}

template< class T >
inline bool
base_point< T >::is_affine_point () const
{
	if (this->ec.get_model() == elliptic_curve_flags::AFFINE)
		return true;
	else
		return false;
}

template< class T >
inline bool
base_point< T >::is_projective_point () const
{
	if (this->ec.get_model() == elliptic_curve_flags::PROJECTIVE)
		return true;
	else
		return false;
}



template< class T >
inline bool
base_point< T >::is_negative_of (const base_point< T > & P) const
{
	base_point< T > H;
	negate(H, P);
	return (*this == H);
}



template< class T >
inline bool
operator == (const base_point< T > & P, const base_point< T > & Q)
{
	return P.is_equal(Q);
}



template< class T >
inline bool
operator != (const base_point< T > & P, const base_point< T > & Q)
{
	return !P.is_equal(Q);
}



template< class T >
inline void
swap (base_point< T > & P, base_point< T > & Q)
{
	P.swap(Q);
}



//
// arithmetic via operators
//

template< class T >
inline base_point< T >&
operator += (base_point< T > & P, const base_point< T > & Q)
{
	add(P, P, Q);
	return P;
}



template< class T >
inline base_point< T >&
operator -= (base_point< T > & P, const base_point< T > & Q)
{
	subtract(P, P, Q);
	return P;
}



template< class T >
inline base_point< T >
operator + (const base_point< T > & P, const base_point< T > & Q)
{
	base_point< T > R(P.get_curve());

	add(R, P, Q);
	return R;
}



template< class T >
inline base_point< T >
operator - (const base_point< T > & P, const base_point< T > & Q)
{
	base_point< T > R(P.get_curve());

	subtract(R, P, Q);
	return R;
}



template< class T >
inline base_point< T >
operator - (const base_point< T > & P)
{
	base_point< T > R(P.get_curve());

	negate(R, P);
	return R;
}



template< class T >
inline base_point< T >
operator * (const bigint & n, const base_point< T > & P)
{
	base_point< T > R(P.get_curve());

	multiply(R, n, P);
	return R;
}



template< class T >
inline std::istream &
operator >> (std::istream & in, base_point< T > & P)
{
	P.read(in);
	return in;
}



template< class T >
inline std::ostream &
operator << (std::ostream & out, const base_point< T > & P)
{
	P.write(out);
	return out;
}



#if defined(_MSC_VER) || !defined(LIDIA_NO_FRIEND_PROTOTYPES)
template< class T >
void negate(base_point< T > & r, const base_point< T > & x);

template< class T >
void add(base_point< T > & R, const base_point< T > & P1, const base_point< T > & Q);

template< class T >
void subtract(base_point< T > & r, const base_point< T > & x, const base_point< T > & y);

template< class T >
void multiply_by_2(base_point< T > & R, const base_point< T > & P);

template< class T >
void multiply (base_point< T > & , const bigint &, const base_point< T > &);
#endif



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/base_point.cc"
#endif



#endif	// LIDIA_BASE_POINT_H_GUARD_
