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


#ifndef LIDIA_BASE_ELLIPTIC_CURVE_CC_GUARD_
#define LIDIA_BASE_ELLIPTIC_CURVE_CC_GUARD_


#ifndef LIDIA_BASE_ELLIPTIC_CURVE_H_GUARD_
# include	"LiDIA/elliptic_curves/base_elliptic_curve.h"
#endif
#ifndef LIDIA_ELLIPTIC_CURVE_REP_H_GUARD_
# include	"LiDIA/elliptic_curves/elliptic_curve_rep.h"
#endif
#ifndef LIDIA_ELLIPTIC_CURVE_REP_BIGINT_H_GUARD_
# include	"LiDIA/elliptic_curves/elliptic_curve_rep_bigint.h"
#endif
#ifndef LIDIA_POINT_OPERATIONS_H_GUARD_
# include	"LiDIA/elliptic_curves/point_operations.h"
#endif
#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



//
// constructor / destructor
//

template< class T >
base_elliptic_curve< T >::base_elliptic_curve()
{
	debug_handler ("base_elliptic_curve< T >",
		       "base_elliptic_curve< T > ()");

	this->e = NULL;
}



template< class T >
base_elliptic_curve< T >::base_elliptic_curve(const base_elliptic_curve< T > & E)
{
	debug_handler ("base_elliptic_curve< T >",
		       "base_elliptic_curve< T > (const base_elliptic_curve< T > &)");

	this->e = E.e;
	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >",
				    "base_elliptic_curve< T > (const base_elliptic_curve< T > &)",
				    "this->e == NULL");
	else {
		this->e->inc_ref_counter();
		//std::cout << "base_elliptic_curve(base_elliptic_curve)::inc " << e->get_ref_counter() << std::endl;
	}
}



template< class T >
base_elliptic_curve< T >::base_elliptic_curve(const T & x4, const T & x6,
					      elliptic_curve_flags::curve_model m)
{
	debug_handler ("base_elliptic_curve< T >",
		       "base_elliptic_curve< T > (const T& x4, const T& x6, "
		       "elliptic_curve_flags::curve_model)");

	this->e = new elliptic_curve_rep< T > (x4, x6, m);
}



template< class T >
base_elliptic_curve< T >::base_elliptic_curve(const T & x1, const T & x2, const T & x3,
					      const T & x4, const T & x6,
					      elliptic_curve_flags::curve_model m)
{
	debug_handler ("base_elliptic_curve< T >",
		       "base_elliptic_curve< T > (x1, x2, x3, x4, x6, "
		       "elliptic_curve_flags::curve_model m)");

	this->e = new elliptic_curve_rep< T > (x1, x2, x3, x4, x6, m);
}



template< class T >
base_elliptic_curve< T >::base_elliptic_curve(elliptic_curve_rep< T > * const e1)
{
	debug_handler ("base_elliptic_curve< T >",
		       "base_elliptic_curve< T > (elliptic_curve_rep< T > *)");

	this->e = e1;
	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "base_elliptic_curve< T > (elliptic_curve_rep< T > *)",
				    "e == NULL");
	else
		this->e->inc_ref_counter();
}



template< class T >
base_elliptic_curve< T >::~base_elliptic_curve()
{
	debug_handler ("base_elliptic_curve< T >",
		       "~base_elliptic_curve< T > ()");

	if (this->e != NULL) {
		if (this->e->get_ref_counter() == 1)
			delete this->e;
		else
			this->e->dec_ref_counter();
	}
}



//
// assignment
//

template< class T >
void
base_elliptic_curve< T >::assign (const base_elliptic_curve< T > & E)
{
	debug_handler ("base_elliptic_curve< T >",
		       "assign(const base_elliptic_curve< T >");

	if (this != &E) {
		if (E.e == NULL)
			lidia_error_handler("base_elliptic_curve< T >::"
					    "assign(const base_elliptic_curve< T > &)",
					    "e == NULL");
		else if (this->e == NULL) {
			this->e = E.e;
			this->e->inc_ref_counter();
		}
		else if (this->e != E.e) {
			if (this->e->get_ref_counter() == 1)
				delete this->e;
			else
				this->e->dec_ref_counter();

			this->e = E.e;
			this->e->inc_ref_counter();
		}
	}
}



template< class T >
void
base_elliptic_curve< T >::set_coefficients(const T & x4, const T & x6,
					   elliptic_curve_flags::curve_model m)

{
	debug_handler ("base_elliptic_curve< T >",
		       "set_coefficients(x4, x6, "
		       "elliptic_curve_flags::curve_model m)");

	if (this->e == NULL)
		this->e = new elliptic_curve_rep< T > (x4, x6, m);
	else
		this->e->set_coefficients(x4, x6, m);
}



template< class T >
void
base_elliptic_curve< T >::set_coefficients(const T & x1, const T & x2, const T & x3,
					   const T & x4, const T & x6,
					   elliptic_curve_flags::curve_model m)
{
	debug_handler ("base_elliptic_curve< T >",
		       "set_coefficients(x1, x2, x3, x4, x6, "
		       "elliptic_curve_flags::curve_model m)");

	if (this->e == NULL)
		this->e = new elliptic_curve_rep< T > (x1, x2, x3, x4, x6, m);
	else
		this->e->set_coefficients(x1, x2, x3, x4, x6, m);
}



//
// access
//

template< class T >
const T &
base_elliptic_curve< T >::discriminant() const
{
	debug_handler ("base_elliptic_curve< T >",
		       "discriminant() const");

	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "discriminant() const",
				    "e == NULL");

	return this->e->discriminant();
}



template< class T >
const T &
base_elliptic_curve< T >::get_a1() const
{
	debug_handler ("base_elliptic_curve< T >",
		       "get_a1() const");
	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "get_a1() const",
				    "e == NULL");
	return this->e->get_a1();
}



template< class T >
const T &
base_elliptic_curve< T >::get_a2() const
{
	debug_handler ("base_elliptic_curve< T >",
		       "get_a2() const");

	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "get_a2() const",
				    "e == NULL");
	return this->e->get_a2();
}



template< class T >
const T &
base_elliptic_curve< T >::get_a3() const
{
	debug_handler ("base_elliptic_curve< T >",
		       "get_a3() const");

	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "get_a3() const",
				    "e == NULL");
	return this->e->get_a3();
}



template< class T >
const T &
base_elliptic_curve< T >::get_sqrt4_a6() const
{
	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "get_a4() const",
				    "e == NULL");
	return this->e->get_a4_intern();
}



template< class T >
const T &
base_elliptic_curve< T >::get_a4() const
{
	debug_handler ("base_elliptic_curve< T >",
		       "get_a4() const");

	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "get_a4() const",
				    "e == NULL");
	return this->e->get_a4();
}



template< class T >
const T &
base_elliptic_curve< T >::get_a6() const
{
	debug_handler ("base_elliptic_curve< T >",
		       "get_a6() const");

	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "get_a6() const",
				    "e == NULL");
	return this->e->get_a6();
}



template< class T >
const T &
base_elliptic_curve< T >::get_b2() const
{
	debug_handler ("base_elliptic_curve< T >",
		       "get_b2() const");
	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "get_b2() const",
				    "e == NULL");
	return this->e->get_b2();
}



template< class T >
const T &
base_elliptic_curve< T >::get_b4() const
{
	debug_handler ("base_elliptic_curve< T >",
		       "get_b4() const");
	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "get_b4() const",
				    "e == NULL");
	return this->e->get_b4();
}



template< class T >
const T &
base_elliptic_curve< T >::get_b6() const
{
	debug_handler ("base_elliptic_curve< T >",
		       "get_b6() const");
	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "get_b6() const",
				    "e == NULL");
	return this->e->get_b6();
}



template< class T >
const T &
base_elliptic_curve< T >::get_b8() const
{
	debug_handler ("base_elliptic_curve< T >",
		       "get_b8() const");
	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "get_b8() const",
				    "e == NULL");
	return this->e->get_b8();
}



template< class T >
const T &
base_elliptic_curve< T >::get_c4() const
{
	debug_handler ("base_elliptic_curve< T >",
		       "get_c4() const");
	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "get_c4() const",
				    "e == NULL");
	return this->e->get_c4();
}



template< class T >
const T &
base_elliptic_curve< T >::get_c6() const
{
	debug_handler ("base_elliptic_curve< T >",
		       "get_c6() const");
	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "get_c6() const",
				    "e == NULL");
	return this->e->get_c6();
}



template< class T >
void
base_elliptic_curve< T >::get_ai(T& t1, T& t2, T& t3, T& t4, T& t6) const
{
	debug_handler ("base_elliptic_curve< T >",
		       "get_ai(T&t1, T&t2, T&t3, T&t4, T&t6) const");

	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "get_ai(t1, t2, t3, t4, t6) const",
				    "e == NULL");

	this->e->get_ai(t1, t2, t3, t4, t6);
}



template< class T >
void
base_elliptic_curve< T >::get_bi(T& t2, T& t4, T& t6, T& t8) const
{
	debug_handler ("base_elliptic_curve< T >",
		       "get_bi(T&t2, T&t4, T&t6, T&t8) const");

	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "get_bi(t2, t4, t6, t8) const",
				    "e == NULL");

	this->e->get_bi(t2, t4, t6, t8);
}



template< class T >
void
base_elliptic_curve< T >::get_ci(T& t4, T& t6) const
{
	debug_handler ("base_elliptic_curve< T >",
		       "get_ci(T&t4, T&t6) const");

	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "get_ci(t4, t6) const",
				    "e == NULL");
	this->e->get_ci(t4, t6);
}



template< class T >
const point_operations< T > &
base_elliptic_curve< T >::get_point_operations() const
{
	debug_handler ("base_elliptic_curve< T >",
		       "get_point_operations()const");

	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "get_point_operations() const",
				    "e == NULL");

	return this->e->get_point_operations();
}



//
// basic properties
//

template< class T >
bool
base_elliptic_curve< T >::is_null() const
{
	debug_handler ("base_elliptic_curve< T >",
		       "is_null() const");

	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "is_null() const",
				    "e == NULL");
	return this->e->is_null();
}



template< class T >
bool
base_elliptic_curve< T >::is_singular() const
{
	debug_handler ("base_elliptic_curve< T >",
		       "is_singular() const");
	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "is_singular() const",
				    "e == NULL");
	return this->e->is_singular();
}



//
// comparison
//

template< class T >
bool base_elliptic_curve< T >::is_equal (const base_elliptic_curve< T > & E) const
{
	debug_handler ("base_elliptic_curve< T >",
		       "operator == (E) const");

	if (this->e == NULL || E.e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "operator == (E) const",
				    "e == NULL");

	return (this->e == E.e);
}



//
// model setting
//

template< class T >
void
base_elliptic_curve< T >::set_model(elliptic_curve_flags::curve_model m)
{
	debug_handler ("base_elliptic_curve< T >",
		       "set_model(elliptic_curve_flags::curve_model)");
	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "set_model(m)",
				    "e == NULL");

	this->e->set_model(m);
}



template< class T >
elliptic_curve_flags::curve_model
base_elliptic_curve< T >::get_model() const
{
	debug_handler ("base_elliptic_curve< T >",
		       "get_model()");
	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "get_model() const",
				    "e == NULL");
	return this->e->get_model();
}



template< class T >
elliptic_curve_flags::curve_parametrization
base_elliptic_curve< T >::get_parametrization() const
{
	debug_handler ("base_elliptic_curve< T >", "get_parametrization()");
	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::get_parametrization() const",
				    "e == NULL");
	return this->e->get_parametrization();
}



//
// Transformation
//

template< class T >
void
base_elliptic_curve< T >::transform(const T& r, const T& s, const T& t)
{
	// NB  u = 1; the more general case is not implemented here

	debug_handler ("base_elliptic_curve< T >",
		       "transform(constT&, constT&, constT&)");

	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "transform(r, s, t)",
				    "e == NULL");
	this->e->transform(r, s, t);
}



//
// input / output
//

template< class T >
void base_elliptic_curve< T >
::set_verbose (int i)
{
	debug_handler ("base_elliptic_curve< T >",
		       "set_verbose(int)");
	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "set_verbose(i)",
				    "e == NULL");
	this->e->set_verbose(i);
}



template< class T >
int
base_elliptic_curve< T >::verbose () const
{
	debug_handler ("base_elliptic_curve< T >",
		       "verbose()");
	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "get_verbose() const",
				    "e == NULL");
	return this->e->verbose();
}



template< class T >
void
base_elliptic_curve< T >::set_output_mode(long i)
{
	debug_handler ("base_elliptic_curve< T >",
		       "set_output_mode(long i)");

	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "set_output_mode(i) const",
				    "e == NULL");

	this->e->set_output_mode(i);
}



template< class T >
void
base_elliptic_curve< T >::output_short(std::ostream & out) const
{
	debug_handler ("base_elliptic_curve< T >",
		       "output_short(std::ostream&)");

	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "output_short(std::ostream&)() const",
				    "e == NULL");

	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "output_short(std::ostream &)",
				    "e == NULL");

	this->e->output_short(out);
}



template< class T >
void
base_elliptic_curve< T >::output_long(std::ostream & out) const
{
	debug_handler ("base_elliptic_curve< T >",
		       "output_long(std::ostream&)");

	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "output_long(std::ostream&)",
				    "e == NULL");

	this->e->output_long(out);
}



template< class T >
void
base_elliptic_curve< T >::output_tex(std::ostream & out) const
{
	debug_handler ("base_elliptic_curve< T >",
		       "output_tex(std::ostream&)");

	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "output_tex(std::ostream&)",
				    "e == NULL");
	this->e->output_tex(out);
}



template< class T >
void
base_elliptic_curve< T >::output_pretty(std::ostream & out) const
{
	debug_handler ("base_elliptic_curve< T >",
		       "output_pretty(std::ostream&)");
	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "output_pretty(std::ostream&)",
				    "e == NULL");

	this->e->output_pretty(out);
}



template< class T >
void
base_elliptic_curve< T >::write (std::ostream & out) const
{
	debug_handler ("base_elliptic_curve< T >",
		       "operator << (std::ostream&, const base_elliptic_curve< T > &)");

	if (this->e == NULL)
		lidia_error_handler("base_elliptic_curve< T >::"
				    "operator << (std::ostream&)",
				    "e == NULL");
	out << *this->e;
}



template< class T >
void
base_elliptic_curve< T >::read(std::istream & in)
{
	debug_handler ("base_elliptic_curve< T >",
		       "input(std::istream&)");

	if (this->e == NULL)
		this->e = new elliptic_curve_rep< T >;

	this->e->read(in);
}



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_BASE_ELLIPTIC_CURVE_CC_GUARD_
