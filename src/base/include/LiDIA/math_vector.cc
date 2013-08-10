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
//	Author	: Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_MATH_VECTOR_CC_GUARD_
#define LIDIA_MATH_VECTOR_CC_GUARD_



#ifndef LIDIA_MATH_VECTOR_H_GUARD_
# include	"LiDIA/math_vector.h"
#endif
#include	"LiDIA/arith.inl"
#ifndef LIDIA_XERROR_H_GUARD_
# include	"LiDIA/xerror.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



//
// debug defines / error defines
//

extern const char *PRT;
extern const char *vector_error_msg[];

#define DV_MV LDBL_VECTOR + 10   // Debug value
#define DM_MV "math_vector"      // Debug message / Error message
#define LMV_ERROR vector_error_msg

//
// debug level
//
//   0 : addtion
//   1 : subtraction
//   2 : multiplication
//   3 : division
//   4 : negation
//   5 : comparisons
//   6 :
//

//
// assignment
//

#if 0
template < class T >
math_vector< T > & math_vector< T >::
operator = (const math_vector< T > & v)
{
	debug_handler_l(DM_BV, "in operator "
			" = (const math_vector< T > &)", DV_BV + 6);

	base_vector< T >::assign(v);
	return *this;
}
#endif



//
// arithmetic procedures
//

//
// addition
//

template< class T >
void math_vector< T >::
add(const math_vector< T > &v, const math_vector< T > &w)
{
	debug_handler_l(DM_MV, "in member - function "
			"add(math_vector< T > &, math_vector< T > &)", DV_MV);

	register lidia_size_t i;
	if (w.length == v.length) {
		this->set_capacity(w.length);
		this->length = w.length;

		for (i = 0; i < w.length; i++)
			LiDIA::add(this->value[i], v.value[i], w.value[i]);
	}
	else
		lidia_error_handler_para(w.length, "w.length", "w.length == v.length",
					 v.length, "v.length", "w.length == v.length",
					 "void math_vector< T >::"
					 "add(const math_vector< T > &v, const math_vector< T > &w)",
					 DM_MV, LMV_ERROR[14]);
}



template< class T >
void math_vector< T >::
add(const math_vector< T > &v, const T &t)
{
	debug_handler_l(DM_MV, "in member - function "
			"add(const math_vector< T > &, const T &)", DV_MV);

	register lidia_size_t i;

	this->set_capacity(v.length);
	this->length = v.length;

	for (i = 0; i < v.length; i++)
		LiDIA::add(this->value[i], v.value[i], t);
}



template< class T >
void math_vector< T >::
add(const T &t, const math_vector< T > &v)
{
	debug_handler_l(DM_MV, "in member - function "
			"add(const T &, const math_vector< T > &)", DV_MV);

	register lidia_size_t i;

	this->set_capacity(v.length);
	this->length = v.length;

	for (i = 0; i < v.length; i++)
		LiDIA::add(this->value[i], t, v.value[i]);
}



//
// subtraction
//

template< class T >
void math_vector< T >::
subtract(const math_vector< T > &v, const math_vector< T > &w)
{
	debug_handler_l(DM_MV, "in member - function "
			"subtract(const math_vector< T > &, const math_vector< T > &)", DV_MV + 1);

	register lidia_size_t i;
	if (w.length == v.length) {
		this->set_capacity(w.length);
		this->length = w.length;

		for (i = 0; i < w.length; i++)
			LiDIA::subtract(this->value[i], v.value[i], w.value[i]);
	}
	else
		lidia_error_handler_para(w.length, "w.length", "w.length == v.length",
					 v.length, "v.length", "w.length == v.length",
					 "void math_vector< T >::"
					 "void math_vector< T >::"
					 "subtract(const math_vector< T > &v, const math_vector< T > &w)",
					 DM_MV, LMV_ERROR[14]);
}



template< class T >
void math_vector< T >::
subtract(const math_vector< T > &v, const T &t)
{
	debug_handler_l(DM_MV, "in member - function "
			"subtract(const math_vector< T > &, const T &)", DV_MV + 1);

	register lidia_size_t i;

	this->set_capacity(v.length);
	this->length = v.length;

	for (i = 0; i < v.length; i++)
		LiDIA::subtract(this->value[i], v.value[i], t);
}



template< class T >
void math_vector< T >::
subtract(const T &t, const math_vector< T > &v)
{
	debug_handler_l(DM_MV, "in member - function "
			"subtract(const T &, const math_vector< T > &)", DV_MV + 1);

	register lidia_size_t i;

	this->set_capacity(v.length);
	this->length = v.length;

	for (i = 0; i < v.length; i++)
		LiDIA::subtract(this->value[i], t, v.value[i]);
}



//
// multiplication
//

template< class T >
void math_vector< T >::
compwise_multiply(const math_vector< T > &v, const math_vector< T > &w)
{
	debug_handler_l(DM_MV, "in member - function "
			"compwise_multiply(const math_vector< T > &, const math_vector< T > &)", DV_MV + 2);

	register lidia_size_t i;
	if (w.length == v.length) {
		this->set_capacity(w.length);
		this->length = w.length;

		for (i = 0; i < w.length; i++)
			LiDIA::multiply(this->value[i], v.value[i], w.value[i]);
	}
	else
		lidia_error_handler_para(w.length, "w.length", "w.length == v.length",
					 v.length, "v.length", "w.length == v.length",
					 "void math_vector< T >::"
					 "compwise_multiply(const math_vector< T > &v, const math_vector< T > &w)",
					 DM_MV, LMV_ERROR[14]);
}



template< class T >
void math_vector< T >::
right_multiply(const math_vector< T > &v, const T &t)
{
	debug_handler_l(DM_MV, "in member - function "
			"multiply(const math_vector< T > &, const T &)" , DV_MV + 2);

	register lidia_size_t i;

	this->set_capacity(v.length);
	this->length = v.length;

	for (i = 0; i < v.length; i++)
		LiDIA::multiply(this->value[i], v.value[i], t);
}



template< class T >
void math_vector< T >::
left_multiply(const T &t, const math_vector< T > &v)
{
	debug_handler_l(DM_MV, "in member - function "
			"multiply(const T &, const math_vector< T > &)", DV_MV + 2);

	register lidia_size_t i;

	this->set_capacity(v.length);
	this->length = v.length;

	for (i = 0; i < v.length; i++)
		LiDIA::multiply(this->value[i], t, v.value[i]);
}



template< class T >
void math_vector< T >::
multiply(T &res, const math_vector< T > &w) const
{
	debug_handler_l(DM_MV, "in member - function "
			"multiply(const math_vector< T > &)", DV_MV + 2);

	register lidia_size_t i;
	T TMP;

	if (this->length == w.length) {
		if (this->length > 0)
			LiDIA::multiply(res, this->value[0], w.value[0]);

		for (i = 1; i < this->length; i++) {
			LiDIA::multiply(TMP, this->value[i], w.value[i]);
			LiDIA::add(res, res, TMP);
		}
	}
	else
		lidia_error_handler_para(this->length, "length", "length == w.length",
					 w.length, "w.length", "length == w.length",
					 "T math_vector< T >::"
					 "multiply(T &, const math_vector< T > &w) const",
					 DM_MV, LMV_ERROR[14]);
}



//
// division
//

template< class T >
void math_vector< T >::
compwise_divide(const math_vector< T > &v, const math_vector< T > &w)
{
	debug_handler_l(DM_MV, "in member - function "
			"compwise_divide(const math_vector< T > &, const math_vector< T > &)", DV_MV + 3);

	register lidia_size_t i;
	if (w.length == v.length) {
		this->set_capacity(w.length);
		this->length = w.length;

		for (i = 0; i < w.length; i++)
			LiDIA::divide(this->value[i], v.value[i], w.value[i]);
	}
	else
		lidia_error_handler_para(w.length, "w.length", "w.length == v.length",
					 v.length, "v.length", "w.length == v.length",
					 "void math_vector< T >::"
					 "compwise_divide(const math_vector< T > &v, const math_vector< T > &w)",
					 DM_MV, LMV_ERROR[14]);
}



template< class T >
void math_vector< T >::
divide(const math_vector< T > &v, const T &t)
{
	debug_handler_l(DM_MV, "in member - function "
			"divide(const math_vector< T > &, const T &)", DV_MV + 4);

	register lidia_size_t i;

	this->set_capacity(v.length);
	this->length = v.length;

	for (i = 0; i < v.length; i++)
		LiDIA::divide(this->value[i], v.value[i], t);
}



template< class T >
void math_vector< T >::
divide(const T &t, const math_vector< T > &v)
{
	debug_handler_l(DM_MV, "in member - function "
			"divide(const T &, const math_vector< T > &)", DV_MV + 4);

	if (0) std::cerr << PRT;
	register lidia_size_t i;

	this->set_capacity(v.length);
	this->length = v.length;

	for (i = 0; i < v.length; i++)
		LiDIA::divide(this->value[i], t, v.value[i]);
}



//
// negation
//

template< class T >
void math_vector< T >::
negate(const math_vector< T > &v)
{
	debug_handler_l(DM_MV, "in member - function "
			"negate(const math_vector< T > &)", DV_MV + 4);

	register lidia_size_t i;
	this->set_capacity(v.length);
	this->length = v.length;

	for (i = 0; i < v.length; i++)
		LiDIA::negate(this->value[i], v.value[i]);
}



//
// sum of squares of two math_vectors
//

template< class T >
T math_vector< T >::
sum_of_squares () const
{
	debug_handler_l(DM_MV, "in member - function "
			"sum_of_squares()", DV_MV + 6);

	register lidia_size_t i;
	T res, TMP;

	if (this->length > 0)
		LiDIA::multiply(res, this->value[0], this->value[0]);
	for (i = 1; i < this->length; i++) {
		LiDIA::multiply(TMP, this->value[i], this->value[i]);
		LiDIA::add(res, res, TMP);
	}
	return res;
}



//
// testing for equality
//

template< class T >
bool math_vector< T >::
equal(const math_vector< T > &v) const
{
	debug_handler_l(DM_MV, "in member - function "
			"equal(const math_vector< T > &)", DV_MV + 5);

	register lidia_size_t i;

	if (this->length != v.length)
		return false;

	for (i = 0; i < this->length; i++)
		if (this->value[i] != v.value[i])
			return false;

	return true;
}



#undef DV_MV
#undef DM_MV
#undef LMV_ERROR



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_MATH_VECTOR_CC_GUARD_
