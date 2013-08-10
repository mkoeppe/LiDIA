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


#ifndef LIDIA_BASE_ELLIPTIC_CURVE_H_GUARD_
#define LIDIA_BASE_ELLIPTIC_CURVE_H_GUARD_


#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_ELLIPTIC_CURVE_FLAGS_H_GUARD_
# include	"LiDIA/elliptic_curve_flags.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class gf_element;

template< class T > class elliptic_curve_rep;
template< class T > class point_operations;
template< class T > class base_point;



template< class T >
class base_elliptic_curve
{
protected:
	elliptic_curve_rep< T > * e;

	//
	// constructor / destructor
	//
	base_elliptic_curve();
	base_elliptic_curve(const base_elliptic_curve< T > &);
	base_elliptic_curve(const T & x4, const T & x6,
			    elliptic_curve_flags::curve_model m =
			    elliptic_curve_flags::AFFINE);

	base_elliptic_curve(const T & x1, const T & x2,
			    const T & x3, const T & x4, const T & x6,
			    elliptic_curve_flags::curve_model m =
			    elliptic_curve_flags::AFFINE);

	base_elliptic_curve(elliptic_curve_rep< T > * e1);

	virtual ~base_elliptic_curve();

	//
	// assignment
	//
public:
	base_elliptic_curve< T > & operator = (
		const base_elliptic_curve< T > & e);

	void assign (const base_elliptic_curve< T > & e);

	virtual void set_coefficients(const T & x4, const T & x6,
				      elliptic_curve_flags::curve_model m =
				      elliptic_curve_flags::AFFINE);

	virtual void set_coefficients(const T & x1, const T & x2,
				      const T & x3, const T & x4, const T & x6,
				      elliptic_curve_flags::curve_model m =
				      elliptic_curve_flags::AFFINE);

	void swap (base_elliptic_curve< T > & e);

	//
	// access
	//
	const T & discriminant() const;
	const T & get_a1() const;
	const T & get_a2() const;
	const T & get_a3() const;
	const T & get_a4() const;

	friend void mult_by_2_gf2nf_projective(base_point< gf_element > &,
					       const base_point< gf_element > &);

protected:
	const T & get_sqrt4_a6() const;
public:
	const T & get_a6() const;
	const T & get_b2() const;
	const T & get_b4() const;
	const T & get_b6() const;
	const T & get_b8() const;
	const T & get_c4() const;
	const T & get_c6() const;


	void get_ai(T& t1, T& t2, T& t3, T& t4, T& t6) const;
	void get_bi(T& t2, T& t4, T& t6, T& t8) const;
	void get_ci(T& t4, T& t6) const;

	const point_operations< T > & get_point_operations() const;

	//
	// model setting
	//
	void set_model(elliptic_curve_flags::curve_model m);
	elliptic_curve_flags::curve_model get_model() const;
	elliptic_curve_flags::curve_parametrization
	get_parametrization() const;

	//
	// basic properties
	//
	bool is_null() const;
	bool is_singular() const;

	//
	// comparison
	//
	bool is_equal (const base_elliptic_curve< T > & E) const;


	//
	// Transformation
	//
	void transform(const T& r, const T& s, const T& t);
	// NB  u = 1; the more general case is not implemented here

	//
	// input / output
	//
	void set_verbose (int i = 0);
	int verbose() const;

	void set_output_mode(long i);
	void output_short(std::ostream & out) const;
	void output_long(std::ostream & out) const;
	void output_tex(std::ostream & out) const;
	void output_pretty(std::ostream & out) const;

	void read(std::istream & in);
	void write(std::ostream & out) const;

};



template< class T >
inline base_elliptic_curve< T > &
base_elliptic_curve< T >::operator = (const base_elliptic_curve< T > & E)
{
	debug_handler ("base_elliptic_curve< T >",
		       "operator = (const base_elliptic_curve< T >");

	this->assign(E);
	return *this;
}



template< class T >
inline void
base_elliptic_curve< T >::swap (base_elliptic_curve< T > & ec)
{
	elliptic_curve_rep< T > * E = ec.e;
	ec.e = this->e;
	this->e = E;
}



template< class T >
inline bool
operator == (const base_elliptic_curve< T > & bc1, const base_elliptic_curve< T > & bc2)
{
	return bc1.is_equal(bc2);
}



template< class T >
inline bool
operator != (const base_elliptic_curve< T > & bc1, const base_elliptic_curve< T > & bc2)
{
	return !bc1.is_equal(bc2);
}



template< class T >
inline void
swap (base_elliptic_curve< T > & e1, base_elliptic_curve< T > & e2)
{
	e1.swap(e2);
}



template< class T >
inline std::istream &
operator >> (std::istream & in, base_elliptic_curve< T > & bc)
{
	debug_handler ("base_elliptic_curve< T >",
		       "operator >> (std::istream&, base_elliptic_curve< T > &)");

	bc.read(in);
	return in;
}



template< class T >
inline std::ostream &
operator << (std::ostream & out, const base_elliptic_curve< T > & bc)
{
	debug_handler ("base_elliptic_curve< T >",
		       "operator << (std::ostream&, const base_elliptic_curve< T > &)");

	bc.write(out);
	return out;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/elliptic_curves/base_elliptic_curve.cc"
#endif



#endif	// LIDIA_BASE_ELLIPTIC_CURVE_H_GUARD_
