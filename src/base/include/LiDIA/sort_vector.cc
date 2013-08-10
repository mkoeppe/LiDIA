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


#ifndef LIDIA_SORT_VECTOR_CC_GUARD_
#define LIDIA_SORT_VECTOR_CC_GUARD_



#ifndef LIDIA_SORT_VECTOR_H_GUARD_
# include	"LiDIA/sort_vector.h"
#endif
#ifndef LIDIA_ARRAY_FUNCTIONS_H_GUARD_
# include	"LiDIA/base/array_functions.h"
#endif
#include	"LiDIA/arith.inl"
#ifndef LIDIA_ERROR_H_GUARD_
# include	"LiDIA/error.h"
#endif
#include "LiDIA/precondition_error.h"



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

#define DV_SV LDBL_VECTOR + 10   // Debug value
#define DM_SV "sort_vector"      // Debug message / Error message
#define LSV_ERROR vector_error_msg

//
// debug level
//
//   0 : assignment
//   1 : sort - direction
//   2 : sort - functions
//   3 : insert into a sort vector
//   4 : delete elements
//   5 : lex. compare
//   6 : delete copies
//

//
// assignment
//

template< class T >
sort_vector< T > & sort_vector< T >::
operator = (const sort_vector< T > & v)
{
	debug_handler_l(DM_SV, "in operator "
			" = (const sort_vector< T > &)", DV_SV);

	if (0) std::cerr << PRT;
	if (this != &v) {
		this->set_capacity(v.length);
		this->length = v.length;

		copy_data(this->value, v.value, this->length);

		this->sort_dir = v.sort_dir;
		this->el_cmp = v.el_cmp;
	}
	return *this;
}



//
// sort - direction
//

template< class T >
void sort_vector< T >::
set_sort_direction(unsigned long sd)
{
	debug_handler_l(DM_SV, "in member - function "
			"set_sort_dir(unsigned long)", DV_SV + 1);

	if (sd != vector_flags::sort_vector_up && sd != vector_flags::sort_vector_down)
		precondition_error_handler(sd, "sd", "sd == SORT_VECTOR_UP or sd == SORT_VECTOR_DOWN",
				    "void sort_vector< T >::"
				    "set_sort_direction(unsigned long sd)",
				    DM_SV, LSV_ERROR[10]);
	else {
		this->sort_dir = static_cast<char>(sd);
		this->el_cmp = NULL;
	}
}



template< class T >
void sort_vector< T >::
set_sort_direction(int (*cmp)(const T & a, const T & b))
{
	debug_handler_l(DM_SV, "in member - function "
			"set_sort_dir(int (*cmp)(const T &, const T &)", DV_SV + 1);

	this->sort_dir = vector_flags::sort_vector_cmp;
	this->el_cmp = cmp;
}



//
// sort-function with compare-function as parameter
//

template< class T >
void sort_vector< T >::
sort(int (*cmp)(const T & a, const T & b), lidia_size_t l, lidia_size_t r)
{
	debug_handler_l(DM_SV, "in member - function "
			"sort(int (*cmp)(const T &, const T &), lidia_size_t, lidia_size_t)", DV_SV + 2);

	r = (r == -1) ? this->length - 1 : r;

	register lidia_size_t gap, i, j, len;

	if (this->length == 0 && l == 0 && r == -1)  // empty vector and default - cut
		return;

	if (l< 0 || l >= this->length || r< 0 || r >= this->length || l > r)
		precondition_error_handler(l, "l", "0 <= l < length",
				    r, "r", "0 <= r < length",
				    "void sort_vector< T >::"
				    "sort(int (*cmp)(const T & a, const T & b), lidia_size_t l, lidia_size_t r)",
				    DM_SV, LSV_ERROR[11]);

	T *ptr = &this->value[l];
	len = r - l + 1; // length of selected cut

	for (gap = len / 2; gap > 0; gap /= 2)
		for (i = gap; i < len; i++)
			for (j = i - gap; j >= 0; j -= gap) {
				if (cmp (ptr[j], ptr[j + gap]) < 0)
					break;
				LiDIA::swap(ptr[j], ptr[j + gap]);
			}
}



//
// sort-functions using comparator - functions
//

template< class T >
void sort_vector< T >::
sort(unsigned long sort_direction, lidia_size_t l, lidia_size_t r)
{
	debug_handler_l(DM_SV, "in member - function "
			"sort(unsigned long, lidia_size_t, lidia_size_t)", DV_SV + 2);

	sort_direction = (sort_direction == vector_flags::sort_vector_def) ? this->sort_dir : sort_direction;
	r = (r == -1) ? this->length - 1 : r;

	if (this->length == 0 && l == 0 && r == -1)  // empty vector and default-cut
		return;

	if (l< 0 || l >= this->length || r< 0 || r >= this->length || l > r)
		precondition_error_handler(l, "l", "0 <= l < length",
				    r, "r", "0 <= r < length",
				    "void sort_vector< T >::"
				    "sort(unsigned long sort_direction, lidia_size_t l, lidia_size_t r)",
				    DM_SV, LSV_ERROR[11]);

	// call of corresponding function
	switch(sort_direction) {
	case vector_flags::sort_vector_up:
		this->sort_up(l, r);
		break;

	case vector_flags::sort_vector_down:
		this->sort_down(l, r);
		break;

	case vector_flags::sort_vector_cmp:
		this->sort(this->el_cmp, l, r);
		break;

	default:
		precondition_error_handler(sort_direction, "sort_direction", "sort_direction in "
				    "{SORT_VECTOR_UP, SORT_VECTOR_DOWN, SORT_VECTOR_CMP}",
				    "void sort_vector< T >::"
				    "sort(unsigned long sort_direction, lidia_size_t l, lidia_size_t r)",
				    DM_SV, LSV_ERROR[10]);
	}
}



template< class T >
void sort_vector< T >::
sort_up(lidia_size_t l, lidia_size_t r)
{
	debug_handler_l(DM_SV, "in member - function "
			"sort_up(lidia_size_t, lidia_size_t)", DV_SV + 2);

	lidia_size_t gap, i, j;
	T *ptr = & this->value[l];
	r = r - l + 1; // length of selected cut

	for (gap = r / 2; gap > 0; gap /= 2)
		for (i = gap; i < r; i++)
			for (j = i - gap; j >= 0; j -= gap) {
				if (comparator< T >::less_equal(ptr[j], ptr[j + gap]) > 0)
					break;
				LiDIA::swap(ptr[j], ptr[j + gap]);
			}
}



template< class T >
void sort_vector< T >::
sort_down(lidia_size_t l, lidia_size_t r)
{
	debug_handler_l(DM_SV, "in member - function "
			"sort_down(lidia_size_t, lidia_size_t)", DV_SV + 2);

	register lidia_size_t gap, i, j;
	T *ptr = & this->value[l];
	r = r - l + 1; // length of selected cut

	for (gap = r / 2; gap > 0; gap /= 2)
		for (i = gap; i < r; i++)
			for (j = i - gap; j >= 0; j -= gap) {
				if (comparator< T >::greater_equal(ptr[j], ptr[j + gap]) > 0)
					break;
				LiDIA::swap(ptr[j], ptr[j + gap]);
			}
}



//
// function for linear search in a vector
//

template< class T >
bool sort_vector< T >::
linear_search(const T & x, lidia_size_t & i) const
{
	debug_handler_l(DM_SV, "in member - function "
			"search(const T &, lidia_size_t &)", DV_SV + 2);

	for (i = 0; i< this->length && !comparator < T >::equal(x, this->value[i]); i++);
	return (i < this->length) ? true : false;
}



//******  functions for binary search in sorted vectos  *****

//
//  NOTE : binary search can only be applied to a vector
//         which is sorted in an appropriate order
//

template< class T >
bool sort_vector< T >::
bin_search(const T & x, lidia_size_t & pos, int (*cmp)(const T & a, const T & b), lidia_size_t l, lidia_size_t r) const
{
	debug_handler_l(DM_SV, "in member - function "
			"bin_search(const T &, lidia_size_t, int (*cmp)(const T &, const T &), "
			"lidia_size_t, lidia_size_t)", DV_SV + 2);

	register int c;

	if (r == -1) r = this->length-1;

	if (this->length == 0 && l == 0 && r == -1)  // empty vector and default-cut
	{
		pos = 0;
		return false;
	}

	if (l< 0 || l >= this->length || r< 0 || r >= this->length || l > r)
		precondition_error_handler(l, "l", "0 <= l < length",
				    r, "r", "0 <= r < length",
				    "bool sort_vector< T >::"
				    "bin_search(const T & x, lidia_size_t & pos, int (*cmp)(const T & a, const T & b), "
				    "lidia_size_t l, lidia_size_t r) const",
				    DM_SV, LSV_ERROR[11]);

	// pick center element
	pos = (l + r) / 2;
	while (l <= r) {
		c = cmp(x, this->value[pos]);
		if (c == 0) {
			// search for first occurence of x
			for (; pos > 0 && cmp(x, this->value[pos - 1]) == 0; pos--);
			return true;
		}
		else
			if (c < 0)        // x must occur before value[pos]
				r = pos - 1; // cont. search in value[l], ..., value[pos-1]
			else              // x must occur after value[pos]
				l = pos + 1; // cont. search in value[pos+1], ..., value[r]
		pos = (l + r) / 2;
	}
	pos = r + 1;
	return false;
}



template< class T >
bool sort_vector< T >::
bin_search(const T & x, lidia_size_t & pos, unsigned long sort_direction, lidia_size_t l, lidia_size_t r) const
{
	debug_handler_l(DM_SV, "in member - function "
			"bin_search(const T &, lidia_size_t, unsigned long, lidia_size_t, lidia_size_t)", DV_SV + 2);

	register bool rc = false;

	if (sort_direction == vector_flags::sort_vector_def) sort_direction = this->sort_dir;
	if (r == -1) r = this->length - 1;

	if (this->length == 0 && l == 0 && r == -1)  // empty vector and default-cut
	{
		pos = 0;
		return false;
	}

	if (l< 0 || l >= this->length || r< 0 || r >= this->length || l > r)
		precondition_error_handler (l, "l", "0 <= l < length",
				     r, "r", "0 <= r < length",
				     "bool sort_vector< T >::"
				     "bin_search(const T & x, lidia_size_t & pos, unsigned long sort_direction, lidia_size_t l, "
				     "lidia_size_t r) const",
				     DM_SV, LSV_ERROR[11]);

	switch (sort_direction) {
	case vector_flags::sort_vector_up:             // vector is sorted in ascending order
		rc = this->bin_search_up(x, pos, l, r);
		break;

	case vector_flags::sort_vector_down:           // vector is sorted in descending order
		rc = this->bin_search_down(x, pos, l, r);
		break;

	case vector_flags::sort_vector_cmp:
		rc = this->bin_search(x, pos, this->el_cmp, l, r);
		break;

	default:
		precondition_error_handler(sort_direction, "sort_direction", "sort_direction in "
				    "{SORT_VECTOR_UP, SORT_VECTOR_DOWN, SORT_VECTOR_CMP}",
				    "bool sort_vector< T >::"
				    "bin_search(const T & x, lidia_size_t &pos, unsigned long sort_direction, "
				    "lidia_size_t l, lidia_size_t r) const",
				    DM_SV, LSV_ERROR[10]);
	}
	return rc;
}



template< class T >
bool sort_vector< T >::
bin_search_up(const T & x, lidia_size_t & pos, lidia_size_t l, lidia_size_t r) const
{
	debug_handler_l(DM_SV, "in member - function "
			"bin_search_up(const T &, lidia_size_t &, lidia_size_t, lidia_size_t)", DV_SV + 2);

	// pick center element
	pos = (l + r) / 2;
	while (l <= r) {
		if (comparator< T >::less_equal(x, this->value[pos]))
			if (comparator< T >::equal(this->value[pos], x)) {
				while (pos > 0 && comparator< T >::equal(this->value[pos-1], x))
					pos--;

				return true;
			}
			else                          // x < value[pos]
				r = pos - 1; // cont. search in value[l], ..., value[pos-1]
		else                            // x > value[pos]
			l = pos + 1; // cont. search in value[pos+1], ..., value[r]
		pos = (l + r) / 2;
	}

	pos = r + 1;
	return false;
}



template< class T >
bool sort_vector< T >::
bin_search_down(const T & x, lidia_size_t & pos, lidia_size_t l, lidia_size_t r) const
{
	debug_handler_l(DM_SV, "in member - function "
			"bin_search_down(const T &, lidia_size_t &, lidia_size_t, lidia_size_t)", DV_SV + 2);

	pos = (l + r) / 2;

	while (l <= r) {
		if (comparator< T >::greater_equal(x, this->value[pos]))
			if (comparator< T >::equal(this->value[pos], x)) {
				while (pos > 0 && comparator< T >::equal(this->value[pos - 1], x))
					pos--;

				return true;
			}
			else
				r = pos - 1; // cont. search in value[l], ..., value[pos-1]
		else
			l = pos + 1; // cont. search in value[pos+1], ..., value[r]
		pos = (l + r) / 2;
	}
	pos = r + 1;
	return false;
}



//
// insertion into a sorted vector
//

template< class T >
void sort_vector< T >::
insert(const T & x, int (*cmp)(const T & a, const T & b), lidia_size_t l, lidia_size_t r)
{
	debug_handler_l(DM_SV, "in member - function "
			"insert(const T &, int (*cmp)(const T &, const T &), lidia_size_t, lidia_size_t)",
			DV_SV + 3);

	lidia_size_t pos, i, found;

	r = (r == -1) ? this->length - 1 : r;
	if (this->length == 0 && l == 0 && r == -1)  // empty vector and default-cut
	{
		this->set_size(1);
		this->value[0] = x;
		return;
	}

	if (l< 0 || l >= this->length || r< 0 || r >= this->length || l > r)
		precondition_error_handler(l, "l", "0 <= l < length",
				    r, "r", "0 <= r< this->length and r >= l",
				    "void sort_vector< T >::"
				    "insert(const T & x, int (*cmp)(const T & a, const T & b), lidia_size_t l, "
				    "lidia_size_t r)",
				    DM_SV, LSV_ERROR[11]);

	// examine, where to put element x
	found = this->bin_search(x, pos, cmp, l, r);
	set_size(this->length + 1);

	// if size has been increased, insert at position 'pos'
	for (i = this->length - 1; i > pos; i--)
		this->value[i] = this->value[i - 1];
	this->value[pos] = x;
}



template< class T >
void sort_vector< T >::
insert(const T & x, unsigned long sort_direction, lidia_size_t l, lidia_size_t r)
{
	debug_handler_l(DM_SV, "in member - function "
			"insert(const T &, unsigned long, lidia_size_t, lidia_size_t)", DV_SV + 3);

	register lidia_size_t i;
	register bool found;
	lidia_size_t pos;

	r = (r == -1) ? this->length - 1 : r;
	if (this->length == 0 && l == 0 && r == -1)  // empty vector and default-cut
	{
		this->set_size(1);
		this->value[0] = x;
		return;
	}

	if (l< 0 || l >= this->length || r< 0 || r >= this->length || l > r)
		precondition_error_handler(l, "l", "0 <= l < length",
				    r, "r", "0 <= r< this->length and r >= l",
				    "void sort_vector< T >::"
				    "insert(const T & x, unsigned long sort_direction, lidia_size_t l, lidia_size_t r)",
				    DM_SV , LSV_ERROR[11]);

	// examine, where to put element x
	found = this->bin_search(x, pos, sort_direction, l, r);
	set_size(this->length + 1);

	// if size has been increased, insert at position 'pos'
	for (i = this->length - 1; i > pos; i--)
		this->value[i] = this->value[i - 1];
	this->value[pos] = x;
}



//
// deleting  elements
//

template< class T >
bool sort_vector< T >::
remove(const T & x, int (*cmp)(const T & a, const T & b), lidia_size_t l, lidia_size_t r)
{
	debug_handler_l(DM_SV, "in member - function "
			"remove (T, cmp, lidia_size_t, lidia_size_t)", DV_SV + 4);

	register lidia_size_t i;
	register bool found;
	lidia_size_t pos;

	r = (r == -1) ? this->length - 1 : r;
	if (l< 0 || l >= this->length || r< 0 || r >= this->length || l > r)
		precondition_error_handler(l, "l", "0 <= l < length",
				    r, "r", "0 <= r < this->length and l <= r",
				    "bool sort_vector< T >::"
				    "remove(const T & x, int (*cmp)(const T & a, const T & b), lidia_size_t l, "
				    "lidia_size_t r)",
				    DM_SV, LSV_ERROR[13]);

	// determine first occurence of x in the vector
	found = bin_search(x, pos, cmp, l, r);

	// if x could be found, remove it from position 'pos'
	if (found) {
		for (i = pos; i < this->length-1; i++)
			this->value[i] = this->value[i + 1];
		this->set_size(this->length - 1);
		return true;
	}
	else
		return false;
}



template< class T >
bool sort_vector< T >::
remove(const T &x, unsigned long sort_direction, lidia_size_t l, lidia_size_t r)
{
	debug_handler_l(DM_SV, "in member - function "
			"remove(const T &, unsigned long, lidia_size_t, lidia_size_t)", DV_SV + 4);

	register lidia_size_t i;
	register bool found;
	lidia_size_t pos;

	r = (r == -1) ? this->length - 1 : r;
	if (l< 0 || l >= this->length || r< 0 || r >= this->length || l > r)
		precondition_error_handler(l, "l", "0 <= l < length",
				    r, "r", "0 <= r < length and l <= r",
				    "bool sort_vector< T >::"
				    "remove(const T &x, unsigned long sort_direction, lidia_size_t l, lidia_size_t r)",
				    DM_SV, LSV_ERROR[13]);

	// determine first occurence of x in the vector
	found = this->bin_search(x, pos, sort_direction, l, r);

	// if x could be found, remove it from position 'pos'
	if (found) {
		for (i = pos; i < this->length - 1; i++)
			this->value[i] = this->value[i + 1];
		this->set_size(this->length - 1);
		return true;
	}
	else
		return false;
}



//
// lex. compare
//

template< class T >
int sort_vector< T >::
lex_compare(sort_vector< T > &w) const
{
	debug_handler_l(DM_SV, "in member - function "
			"lex_compare(sort_vector< T > &)", DV_SV + 5);

	register lidia_size_t i, minl = (this->length < w.length) ? this->length : w.length;

	for (i = 0; i< minl && comparator < T >::equal(this->value[i], w.value[i]); i++);
	if (i < minl)
		if (comparator< T >::less_equal(this->value[i], w.value[i]))
			return (comparator< T >::equal(this->value[i], w.value[i])) ? 0 : -1;
		else
			return 1;
	else
		if (i < w.length)
			return -1;
		else
			return (i < this->length) ? 1 : 0;
}



//
// delete_copies:
// removes the copies of elements of this;
// assumes, that the vector is sorted;
//

template< class T >
void sort_vector< T >::
delete_copies()
{
	debug_handler_l(DM_SV, "in member - function "
			"delete_copies()", DV_SV + 6);

	register lidia_size_t i, j;
	int ident;

	for (i = 0, j = 0; i < this->length; j++) {
		this->value[j] = this->value[i];
		i++;
		ident = 1;
		while (i < this->length && ident) {
			ident = comparator< T >::equal(this->value[j], this->value[i]);
			i++;
		}
		if (!ident)
			i--;
	}
	this->length = j;
}



#undef DV_SV
#undef DM_SV
#undef LSV_ERROR



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_SORT_VECTOR_CC_GUARD_
