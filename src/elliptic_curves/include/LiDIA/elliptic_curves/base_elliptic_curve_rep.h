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


#ifndef LIDIA_BASE_ELLIPTIC_CURVE_REP_H_GUARD_
#define LIDIA_BASE_ELLIPTIC_CURVE_REP_H_GUARD_


#ifndef LIDIA_ELLIPTIC_CURVE_FLAGS_H_GUARD_
# include	"LiDIA/elliptic_curve_flags.h"
#endif
#ifndef LIDIA_LIDIA_REFERENCE_COUNTER_H_GUARD_
# include	"LiDIA/lidia_reference_counter.h"
#endif
#ifndef LIDIA_POINT_OPERATIONS_H_GUARD_
# include	"LiDIA/elliptic_curves/point_operations.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class T >
class base_elliptic_curve_rep : public lidia_reference_counter
{
protected:
	T  a1, a2, a3, a4, a6;
	T  b2, b4, b6, b8;
	T  c4, c6;
	T  delta;
	int info;

	elliptic_curve_flags::curve_parametrization cp;
	elliptic_curve_flags::curve_model           cm;
	elliptic_curve_flags::curve_output_mode output_mode;

	point_operations< T > bcp;

	//
	// constructor / destructor
	//

public:
	base_elliptic_curve_rep();

	base_elliptic_curve_rep(const base_elliptic_curve_rep< T > &);

	base_elliptic_curve_rep(const T & x4,
				const T & x6,
				elliptic_curve_flags::curve_model m = elliptic_curve_flags::AFFINE);

	base_elliptic_curve_rep(const T & x1,
				const T & x2,
				const T & x3,
				const T & x4,
				const T & x6,
				elliptic_curve_flags::curve_model m = elliptic_curve_flags::AFFINE);

	~base_elliptic_curve_rep();

	//
	// initialization
	//
protected:
	void compute_invariants();

	//
	// assignment
	//
public:
	base_elliptic_curve_rep< T > & operator = (const base_elliptic_curve_rep< T > & e);
	void assign (const base_elliptic_curve_rep< T > & e);

	void set_coefficients(const T & x4,
			      const T & x6,
			      elliptic_curve_flags::curve_model m = elliptic_curve_flags::AFFINE);

	void set_coefficients(const T & x1,
			      const T & x2,
			      const T & x3,
			      const T & x4,
			      const T & x6,
			      elliptic_curve_flags::curve_model m = elliptic_curve_flags::AFFINE);


	//
	// model setting
	//
	void set_model(elliptic_curve_flags::curve_model m);
	elliptic_curve_flags::curve_model get_model() const;
	elliptic_curve_flags::curve_parametrization
	get_parametrization() const;

	//
	// access
	//
	const T & discriminant() const;
	const T & get_a1() const;
	const T & get_a2() const;
	const T & get_a3() const;
	const T & get_a4() const;
	const T & get_a4_intern() const;
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
	// basic properties
	//
	bool is_null() const;
	bool is_singular() const;

	//
	// Transformation
	//
	void transform(const T& r, const T& s, const T& t);
	// NB  u = 1; the more general case is not implemented here

	//
	// input / output
	//
	void set_verbose (int i);
	int verbose() const;

	void set_output_mode(long i);
	long get_output_mode() const;

	void output_short(std::ostream & out) const;
	void output_long(std::ostream & out) const;
	void output_tex(std::ostream & out) const;
	void output_pretty(std::ostream & out) const;

	void read(std::istream & in);
	void write(std::ostream & out) const;
};



//
// c'tors and d'tor
//

template< class T >
inline
base_elliptic_curve_rep< T >::base_elliptic_curve_rep ()
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "base_elliptic_curve_rep< T > ()");

	this->output_mode = elliptic_curve_flags::SHORT;
	this->info = 0;
}



template< class T >
inline
base_elliptic_curve_rep< T >::base_elliptic_curve_rep (const base_elliptic_curve_rep< T > & e)
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "base_elliptic_curve_rep< T > (const base_elliptic_curve_rep< T > &)");

	this->assign(e);
}



template< class T >
inline
base_elliptic_curve_rep< T >::base_elliptic_curve_rep (const T & x4, const T & x6,
						       elliptic_curve_flags::curve_model m)
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "base_elliptic_curve_rep< T > (const T& x4, const T& x6, "
		       "elliptic_curve_flags::curve_model)");

	this->set_coefficients(x4, x6, m);
}



template< class T >
inline
base_elliptic_curve_rep< T >::base_elliptic_curve_rep (const T & x1, const T & x2, const T & x3,
						       const T & x4, const T & x6,
						       elliptic_curve_flags::curve_model m)
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "base_elliptic_curve_rep< T > (const T&x1, const T&x2, "
		       "const T&x3, const T&x4, const T&x6"
		       "elliptic_curve_flags::curve_model)");

	this->set_coefficients(x1, x2, x3, x4, x6, m);
}



template< class T >
inline
base_elliptic_curve_rep< T >::~base_elliptic_curve_rep ()
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "~base_elliptic_curve_rep< T > ()");
}



//
// accessors
//

template< class T >
inline void
base_elliptic_curve_rep< T >::set_model (elliptic_curve_flags::curve_model m)
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "set_model(elliptic_curve_flags::curve_model)");
	this->cm = m;
	this->bcp.set(this->cp, this->cm, this->a6.characteristic());
}



template< class T >
inline elliptic_curve_flags::curve_model
base_elliptic_curve_rep< T >::get_model () const
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "get_model()");
	return this->cm;
}



template< class T >
inline elliptic_curve_flags::curve_parametrization
base_elliptic_curve_rep< T >::get_parametrization () const
{
	debug_handler ("base_elliptic_curve_rep< T >", "get_parametrization()");
	return this->cp;
}



//
// access
//
template< class T >
inline const T &
base_elliptic_curve_rep< T >::discriminant () const
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "discriminant() const");
	return this->delta;
}



template< class T >
inline const T &
base_elliptic_curve_rep< T >::get_a1 () const
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "get_a1() const");
	return this->a1;
}



template< class T >
inline const T &
base_elliptic_curve_rep< T >::get_a2 () const
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "get_a2() const");
	return this->a2;
}



template< class T >
inline const T &
base_elliptic_curve_rep< T >::get_a3 () const
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "get_a3() const");
	return this->a3;
}



template< class T >
inline const T &
base_elliptic_curve_rep< T >::get_a4_intern () const
{
	return this->a4;
}



template< class T >
inline const T &
base_elliptic_curve_rep< T >::get_a4 () const
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "get_a4() const");
	if (this->cp != elliptic_curve_flags::GF2N_F)
		return this->a4;
	else
		return this->a3;
}



template< class T >
inline const T &
base_elliptic_curve_rep< T >::get_a6 () const
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "get_a6() const");
	return this->a6;
}



template< class T >
inline const T &
base_elliptic_curve_rep< T >::get_b2 () const
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "get_b2() const");
	return this->b2;
}



template< class T >
inline const T &
base_elliptic_curve_rep< T >::get_b4 () const
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "get_b4() const");
	return this->b4;
}



template< class T >
inline const T &
base_elliptic_curve_rep< T >::get_b6 () const
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "get_b6() const");
	return this->b6;
}



template< class T >
inline const T &
base_elliptic_curve_rep< T >::get_b8 () const
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "get_b8() const");
	return this->b8;
}



template< class T >
inline const T &
base_elliptic_curve_rep< T >::get_c4 () const
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "get_c4() const");
	return this->c4;
}



template< class T >
inline const T &
base_elliptic_curve_rep< T >::get_c6 () const
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "get_c6() const");
	return this->c6;
}



template< class T >
inline void
base_elliptic_curve_rep< T >::get_ai (T& t1, T& t2, T& t3, T& t4, T& t6) const
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "get_ai(T&t1, T&t2, T&t3, T&t4, T&t6) const");

	t1 = this->a1;
	t2 = this->a2;
	t3 = this->a3;
	t4 = this->a4;
	t6 = this->a6;
}



template< class T >
inline void
base_elliptic_curve_rep< T >::get_bi (T& t2, T& t4, T& t6, T& t8) const
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "get_bi(T&t2, T&t4, T&t6, T&t8) const");

	t2 = this->b2;
	t4 = this->b4;
	t6 = this->b6;
	t8 = this->b8;
}



template< class T >
inline void
base_elliptic_curve_rep< T >::get_ci (T& t4, T& t6) const
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "get_ci(T&t4, T&t6) const");

	t4 = this->c4;
	t6 = this->c6;
}



template< class T >
inline const point_operations< T > &
base_elliptic_curve_rep< T >::get_point_operations () const
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "get_point_operations()");
	return this->bcp;
}



//
// basic properties
//

template< class T >
inline bool
base_elliptic_curve_rep< T >::is_null () const
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "is_null() const");

	return (this->a1.is_zero() && this->a2.is_zero() && this->a3.is_zero()
		&& this->a4.is_zero() && this->a6.is_zero());
}



template< class T >
inline bool
base_elliptic_curve_rep< T >::is_singular () const
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "is_singular() const");

	return this->delta.is_zero();
}



//
// assigners
//

template< class T >
inline base_elliptic_curve_rep< T > &
base_elliptic_curve_rep< T >::operator = (const base_elliptic_curve_rep< T > & e)
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "operator = (const base_elliptic_curve_rep< T >");

	assign(e);
	return *this;
}



//
// I/O
//

template< class T >
inline std::ostream &
operator << (std::ostream & out, const base_elliptic_curve_rep< T > & bc)
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "operator >> (std::istream&, base_elliptic_curve_rep< T > &)");

	bc.write(out);
	return out;
}



template< class T >
inline std::istream &
operator >> (std::istream & in, base_elliptic_curve_rep< T > & bc)
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "operator >> (std::istream&, base_elliptic_curve_rep< T > &)");

	bc.read(in);
	return in;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/elliptic_curves/base_elliptic_curve_rep.cc"
#endif



#endif	// LIDIA_BASE_ELLIPTIC_CURVE_H_GUARD_
