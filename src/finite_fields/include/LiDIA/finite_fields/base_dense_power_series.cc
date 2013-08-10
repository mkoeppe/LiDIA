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
//	Author	:
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_BASE_DENSE_POWER_SERIES_CC_GUARD_
#define LIDIA_BASE_DENSE_POWER_SERIES_CC_GUARD_


#ifndef LIDIA_BASE_DENSE_POWER_SERIES_H_GUARD_
# include	"LiDIA/finite_fields/base_dense_power_series.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



//
// ***** constructors / destructor
//

template< class T >
base_dense_power_series< T >::base_dense_power_series ()
{
	debug_handler ("base_dense_power_series< T >",
		       "base_dense_power_series()");

	this->coeff = new math_vector< T >;
	this->max_num_coeff = 0;
}



template< class T >
base_dense_power_series< T >::base_dense_power_series (const T & elem, lidia_size_t l)
{
	debug_handler ("base_dense_power_series< T >",
		       "base_dense_power_series(const T&, lidia_size_t)");

	T   zero_elem;
	lidia_size_t i, prec;

	this->max_num_coeff = 0;
	zero_elem.assign_zero ();

	if (l < 0 || elem == zero_elem) {
		this->first = l;
		this->last = l;
		this->coeff = new math_vector< T >;
		this->coeff->set_capacity (1);
		(*this->coeff)[0] = zero_elem;
	}
	else {
		prec = l+1;
		this->first = 0;
		this->last = l;
		this->coeff = new math_vector< T >;
		this->coeff->set_capacity (prec);
		(*this->coeff)[0] = elem;

		for (i = 1; i < prec; i++)
			(*this->coeff)[i] = zero_elem;
	}
}



template< class T >
base_dense_power_series< T >::base_dense_power_series (const base_vector< T > & c, lidia_size_t f)
{
	debug_handler ("base_dense_power_series< T >",
		       "base_dense_power_series(const base_vector< T > &, lidia_size_t)");

	this->coeff = new math_vector< T >;
	this->max_num_coeff = 0;

	lidia_size_t prec, start;
	lidia_size_t i;
	T zero_T;

	zero_T.assign_zero ();
	prec = c.size();
	start = prec;


	// find index of the first non-zero element in c

	for (i = 0; i < prec && start == prec; i ++) {
		if (c[i] != zero_T)
			start = i;
	}


	// all elements are zero

	if (start == prec) {
		this->assign_zero (f + prec - 1);
	}

	// a non-zero element was found

	else {
		this->coeff->set_capacity (prec - start);
		this->first = f + start;
		this->last = f + prec - 1;

		for (i = start; i < prec; i++) {
			(*this->coeff)[i-start] = c[i];
		}
	}
}



template< class T >
base_dense_power_series< T >::base_dense_power_series (const base_dense_power_series< T > & a)
{
	debug_handler ("base_dense_power_series< T >",
		       "base_dense_power_series(const base_dense_power_series< T > &)");

	lidia_size_t non_zero_index;

	this->coeff = new math_vector< T >;
	this->max_num_coeff = a.max_num_coeff;

	if (a.is_zero (non_zero_index)) {
		this->assign_zero (a.last);
	}
	else {
		this->first = a.first + non_zero_index;
		this->last = a.last;

		this->coeff->set_capacity ((a.coeff)->size() - non_zero_index);
		this->coeff->assign (0, *(a.coeff), non_zero_index, (a.coeff)->size()-1);
	}
}



template< class T >
base_dense_power_series< T >::~base_dense_power_series ()
{
	debug_handler ("base_dense_power_series< T >",
		       "~base_dense_power_series()");

	this->coeff->kill();
	delete this->coeff;
}



//
// ***** utility routines *****
//

template< class T >
bool
base_dense_power_series< T >::is_zero (lidia_size_t & non_zero_index) const
{
	debug_handler ("base_dense_power_series< T >", "is_zero()");
	lidia_size_t i, prec;
	bool rc;

	if (this->coeff->size() == 0)
		lidia_error_handler ("base_dense_power_series< T >::is_zero()",
				     "Serie is not initialized.");

	non_zero_index = -1;
	prec = this->coeff->size ();

	for (i = 0; i < prec && non_zero_index == -1; i++) {
		if (! ((*this->coeff)[i].is_zero()))
			non_zero_index = i;
	}

	if (non_zero_index == -1)
		rc = true;
	else
		rc = false;

	return (rc);
}



template< class T >
bool
base_dense_power_series< T >::is_zero () const
{
	debug_handler ("base_dense_power_series< T >", "is_zero()");
	lidia_size_t i, prec;
	lidia_size_t non_zero_index;
	bool rc;

	if (this->coeff->size() == 0)
		lidia_error_handler ("base_dense_power_series< T >::is_zero()",
				     "Serie is not initialized.");

	non_zero_index = -1;
	prec = this->coeff->size ();

	for (i = 0; i < prec && non_zero_index == -1; i++) {
		if (! (*this->coeff)[i].is_zero())
			non_zero_index = i;
	}

	if (non_zero_index == -1)
		rc = true;
	else
		rc = false;

	return (rc);
}



template< class T >
bool
base_dense_power_series< T >::is_one () const
{
	debug_handler ("base_dense_power_series< T >", "is_one()");
	lidia_size_t i, prec;
	lidia_size_t non_zero_index;
	bool rc;

	if (this->coeff->size() == 0)
		lidia_error_handler ("base_dense_power_series< T >::is_one()",
				     "Serie is not initialized.");

	if (this->last < 0)
		rc = false;
	else if (! (*this->coeff)[0].is_one ())
		rc = false;
	else {
		non_zero_index = -1;
		prec = this->coeff->size ();

		for (i = 1; i < prec && non_zero_index == -1; i++) {
			if (! (*this->coeff)[i].is_zero())
				non_zero_index = i;
		}

		if (non_zero_index == -1)
			rc = true;
		else
			rc = false;
	}

	return (rc);
}



// The function allocate() controls the memory-allocation for coeff.  The calling function
// demands cp coefficients to be allocated at least.  If the user has defined an upper bound
// (max_num_coeff) for the number of allocated coefficients and cp exceeds this bound, exactly
// cp elements will be allocated. Otherwise the function tries to allocate n_cp =
// max(2*old-capacity, cp) elements. If this value is greater than max_num_coeff, max_num_coeff
// coefficients will be allocated.  Although this function controls the memory allocation for
// coeff, it is sometimes necessary to allocate memory for a math_vector *c with *c = *coeff;
// this is the reason for the existence of argument c (see function set_coeff(, ,)).

template< class T >
void
base_dense_power_series< T >::allocate (lidia_size_t cp, math_vector< T > * c)
{
	debug_handler ("base_dense_power_series< T >", "allocate(lidia_size_t, math_vector< T > *)");

	lidia_size_t n_cp = c->capacity();

	if (cp > n_cp) {
		n_cp = (n_cp << 1);

		if (cp > n_cp)
			n_cp = cp;

		if (this->max_num_coeff != 0) {
			if (cp > this->max_num_coeff) {
				warning_handler("base_dense_power_series< T >::allocate(lidia_size_t, math_vector< T > *",
						"Number of allocated coefficients exceeds user-defined upper bound.");
				n_cp = cp;
			}
			else {
				if (n_cp > this->max_num_coeff)
					n_cp = this->max_num_coeff;
			}
		}

		c->set_capacity (n_cp);
	}
}



template< class T >
lidia_size_t
base_dense_power_series< T >::get_first () const
{
	debug_handler ("base_dense_power_series< T >", "get_first");

	lidia_size_t non_zero_index;

	if (this->is_zero (non_zero_index))
		return this->last;
	else
		return this->first + non_zero_index;
}



template< class T >
lidia_size_t
base_dense_power_series< T >::get_last () const
{
	debug_handler ("base_dense_power_series< T >", "get_last");

	if (this->coeff->size() == 0)
		lidia_error_handler ("base_dense_power_series< T >::get_last()",
				     "Serie is not initialized.");
	return this->last;
}



template< class T >
void
base_dense_power_series< T >::reduce_last (lidia_size_t l)
{
	debug_handler ("base_dense_power_series< T >", "reduce_last(lidia_size_t)");

	if (this->coeff->size() == 0) {
		warning_handler ("base_dense_power_series< T >::reduce_last(lidia_size_t)",
				 "Serie is not initialized, setting failed.");
	}
	else if (l > this->last) {
		lidia_error_handler ("base_dense_power_series< T >::reduce_last(lidia_size_t)",
				     "New value of this->last is too big, setting failed.");
	}
	else {
		if (l >= this->first) {
			this->last = l;
			this->coeff->set_size ((this->last-first+1));
		}
		else {
			this->assign_zero (l);
		}
	}
}



//
// ***** coefficient handling *****
//

template< class T >
void
base_dense_power_series< T >::get_coeff (T & elem, lidia_size_t e)
{
        debug_handler ("base_dense_power_series< T >", "get_coeff(T&, lidia_size_t)");

        if (this->coeff->size() == 0) {
		lidia_error_handler ("base_dense_power_series< T >::get_coeff(T&, lidia_size_t)",
				     "Serie is not initialized.");
	}
        else if (e > this->last) {
		lidia_error_handler ("base_dense_power_series< T >::get_coeff(T&, lidia_size_t)",
				     "Exponent e is too big.");
	}
        else if (e < this->first) {
		elem.assign_zero ();
	}
        else {
		elem = (*this->coeff)[ e - this->first ];
	}
}



template< class T >
void
base_dense_power_series< T >::get (base_vector< T > & c)
{
        debug_handler ("base_dense_power_series< T >", "get(math_vector< T > &)");

        if (this->coeff->size() == 0)
		lidia_error_handler ("base_dense_power_series< T >::get(math_vector< T > &)",
				     "Serie is not initialized.");

	lidia_size_t non_zero_index;

	if (this->is_zero (non_zero_index)) {
		c.set_capacity (1);
		c[0].assign_zero ();
	}
	else {
		c.set_capacity (this->coeff->size() - non_zero_index);
		c.assign (0, *this->coeff, non_zero_index, this->coeff->size() - 1);
	}
}



template< class T >
void
base_dense_power_series< T >::get (T * & c, lidia_size_t & sz)
{
        debug_handler ("base_dense_power_series< T >", "get(T *&, lidia_szie_t&)");

        if (this->coeff->size() == 0)
		lidia_error_handler ("base_dense_power_series< T >::get(T*&)",
				     "Serie is not initialized.");
	lidia_size_t non_zero_index;
        lidia_size_t i;

        if (c != NULL) {
		delete [] c;
		c = NULL;
	}

	if (this->is_zero (non_zero_index)) {
		c = new T [1];
		c[0].assign_zero();
		sz = 1;
	}
	else {
		sz = this->coeff->size();
		c = new T [ sz - non_zero_index ];

		for (i = non_zero_index; i < sz; i++)
			c[i-non_zero_index] = (*this->coeff)[i];

		sz = sz - non_zero_index;
	}
}



template< class T >
void
base_dense_power_series< T >::set_coeff (const T & elem, lidia_size_t e)
{
	debug_handler ("base_dense_power_series< T >", "set_coeff(const T&, lidia_size_t)");
	set_coeff (elem, e, *this);
}



template< class T >
void
base_dense_power_series< T >::set_coeff (const T & elem, lidia_size_t e, const base_dense_power_series< T > & a)
{
	debug_handler ("base_dense_power_series< T >", "set_coeff(const T&, lidia_size_t, base_dense_power_series&)");

	math_vector< T > *A;
	lidia_size_t aprec;
	int ident;

	A = a.coeff;
	aprec = A->size();

	if (&a == this)
		ident = 1;
	else
		ident = 0;

	if (aprec == 0) {
		this->first = e;
		this->last = e;
		this->allocate (1, this->coeff);
		(*this->coeff)[0] = elem;
	}
	else if (e > a.last) {
		if (a.is_zero ()) {
			this->allocate (1, this->coeff);
			(*this->coeff)[0] = elem;
			this->first = e;
			this->last = e;
		}
		else {
			T zero_elem;
			lidia_size_t i, new_prec;

			zero_elem.assign_zero ();
			new_prec = aprec +  (e - a.last);
			this->allocate (new_prec, this->coeff);

			this->last = e;
			(*this->coeff)[ new_prec - 1 ] = elem;

			new_prec --;

			for (i = aprec; i < new_prec; i++)
				(*this->coeff)[i] = zero_elem;

			if (!ident) {
				this->first = a.first;
				this->coeff->assign (0, *A, 0, aprec-1);
			}
		}
	}
	else if (e < a.first) {
		if (! elem.is_zero ()) {
			if (a.is_zero ()) {
				T zero_elem;
				lidia_size_t i, prec;

				zero_elem.assign_zero ();
				prec = a.first - e + 1;

				this->allocate (prec, this->coeff);

				(*this->coeff)[0] = elem;

				for (i = 1; i < prec; i++)
					(*this->coeff)[i] = zero_elem;

				this->first = e;
				this->last = a.last;
			}
			else {
				math_vector< T > *h;

				if (ident)
					h = new math_vector< T >;
				else
					h = this->coeff;

				T zero_elem;
				lidia_size_t i, new_prec;
				lidia_size_t new_added;

				zero_elem.assign_zero ();
				new_added = (a.first - e);
				new_prec = aprec + new_added;
				this->allocate (new_prec, h);

				(*h)[0] = elem;

				for (i = 1; i < new_added; i++)
					(*h)[i] = zero_elem;

				h->set_size (new_added);
				h->concat (*h, *A);

				this->first = e;

				if (ident) {
					delete this->coeff;
					this->coeff = h;
				}
				else
					this->last = a.last;
			}
		}
		else if (!ident) {
			// assign this = a
			this->first = a.first;
			this->last = a.last;
			*this->coeff = *A;
		}
	}
	else {
		// a.first <= e <= a.last
		T zero_elem;

		zero_elem.assign_zero ();

		if (e == a.first && elem == zero_elem) {
			lidia_size_t delta = 1;
			lidia_size_t sz = aprec;

			while (delta < sz) {
				if ((*A)[delta] == zero_elem)
					delta++;
				else
					sz = delta;
			}

			if (delta == aprec) {
				this->assign_zero (a.last);
			}
			else {
				if (ident) {
					this->coeff->shift_left(delta, delta);
					this->first += delta;
				}
				else {
					this->allocate (aprec-delta, this->coeff);
					this->first = a.first + delta;
					this->last = a.last;
					this->coeff->assign (0, *A, delta, aprec-1);
				}
			}
		}
		else {
			// e != a.first || elem != 0
			if (!ident) {
				this->first = a.first;
				this->last = a.last;
				*this->coeff = *A;
			}

			(*this->coeff)[ e - this->first ] = elem;
		}
	}
}



template< class T >
void
base_dense_power_series< T >::set (base_vector< T > & c, lidia_size_t f)
{
        debug_handler ("base_dense_power_series< T >", "set(const math_vector< T > &, lidia_size_t)");

        lidia_size_t delta;
        lidia_size_t sz = c.size();

        T zero_elem;

	zero_elem.assign_zero ();

        if (sz == 0) {
		warning_handler ("base_dense_power_series< T >::set(const math_vector< T > &, lidia_size_t)",
				 "Argument math_vector< T > not initialized, setting failed!");
	}
        else if (c[0] == zero_elem) {
		delta = 1;

		while (delta < sz) {
			if (c[delta] == zero_elem)
				delta++;
			else
				sz = delta;
		}

		sz = c.size();

		if (delta == sz) {
			this->assign_zero (f + sz - 1);
		}
		else {
			this->coeff->set_capacity(sz - delta);
			this->coeff->assign      (0, c, delta, sz-1);
			this->first = f + delta;
			this->last = f + sz - 1;
		}
	}
        else {
		*this->coeff = c;
		this->first = f;
		this->last = f + sz - 1;
	}
}



template< class T >
void
base_dense_power_series< T >::set (const T *c, lidia_size_t prec, lidia_size_t f)
{
        debug_handler ("base_dense_power_series< T >", "set(const T*, lidia_size_t, lidia_size_t)");

        lidia_size_t i;
        lidia_size_t delta;
        lidia_size_t sz = prec;

        if (sz == 0 || c == NULL) {
		warning_handler ("base_dense_power_series< T >::set(const T*, lidia_size_t, lidia_size_t)",
				 "Argument T* not initialized, setting failed!");
	}
        else if (c[0] == 0) {
		delta = 1;

		while (delta < sz) {
			if (c[delta] == 0)
				delta++;
			else
				sz = delta;
		}

		sz = prec;

		if (delta == sz) {
			this->assign_zero (f + sz - 1);
		}
		else {
			this->coeff->set_capacity(sz - delta);

			this->first = f + delta;
			this->last = f + sz - 1;

			for (i = 0; delta < sz; i++, delta++) {
				(*this->coeff)[i] = c[delta];
			}
		}
	}
        else {
		this->first = f;
		this->last = f + sz - 1;

		this->coeff->set_capacity (sz);

		for (i = 0; i < sz; i++) {
			(*this->coeff)[i] = c[i];
		}
	}
}



template< class T >
void
base_dense_power_series< T >::set (const T & elem, lidia_size_t l)
{
	debug_handler ("base_dense_power_series< T >", "set(T&, lidia_size_t)");

	T    zero_elem;
	lidia_size_t  i, prec;

	zero_elem.assign_zero ();

	if (elem == zero_elem || l < 0)
        {
		this->assign_zero (l);
        }
	else
        {
		prec = l+1;
		this->first = 0;
		this->last = l;
		this->coeff->set_capacity (prec);

		(*this->coeff)[0] = elem;

		for (i = 1; i < prec; i++)
			(*this->coeff)[i] = zero_elem;
        }
}



//
// ***** more utility functions *****
//


template< class T >
void
base_dense_power_series< T >::clear ()
{
        debug_handler ("base_dense_power_series< T >", "clear()");

        this->coeff->set_capacity (0);
        this->max_num_coeff = 0;
}



template< class T >
void
base_dense_power_series< T >::normalize ()
{
        debug_handler ("base_dense_power_series< T >", "normalize()");

        if (this->coeff->size() == 0) {
		warning_handler ("base_dense_power_series< T >::normalize()",
				 "Serie is not initialized, normalization failed.");
	}
        else if (this->first > this->last) {
		warning_handler ("base_dense_power_series< T >::normalize()",
				 "first is greater than this->last, normalization failed.");
	}
        else {
		lidia_size_t non_zero_index;

		if (this->is_zero (non_zero_index))
			this->assign_zero (this->last);
		else {
			if (non_zero_index != 0)
				this->coeff->shift_left (non_zero_index, non_zero_index);

			this->coeff->set_capacity ((this->last - this->first) + 1 - non_zero_index);
			this->first += non_zero_index;
		}
	}
}



template< class T >
void
base_dense_power_series< T >::set_max_num_of_coeff (lidia_size_t sz)
{
	debug_handler ("base_dense_power_series< T >", "set_max_num_of_coeff(lidia_size_t)");

	if (sz <= 0)
        {
		warning_handler ("base_dense_power_series< T >::set_max_num_of_coeff(lidia_size_t)",
				 "Ignoring non-positive size.");
	}
	else
		this->max_num_coeff = sz;
}



template< class T >
lidia_size_t
base_dense_power_series< T >::get_max_num_of_coeff ()
{
	debug_handler ("base_dense_power_series< T >", "get_max_num_of_coeff()");
	return (this->max_num_coeff);
}



template< class T >
void
base_dense_power_series< T >::clear_max_num_of_coeff ()
{
	debug_handler ("base_dense_power_series< T >", "clear_max_num_of_coeff(lidia_size_t)");
	this->max_num_coeff = 0;
}



//
// ***** subscripting *****
//


// Returns the coefficient with exponent e.
// This operator can only be used for
// reading - operations. The exponent e must be
// smaller than or equal to last; otherwise
// the lidia_error_handler() is called.

template< class T >
const T
base_dense_power_series< T >::operator[] (lidia_size_t e) const
{
        debug_handler ("base_dense_power_series< T >", "operator[]");

        if (this->coeff->size() == 0) {
		lidia_error_handler ("base_dense_power_series< T >::operator[]",
				     "Serie is not initialized.");
		return T(0);
	}
        else if (e > this->last) {
		lidia_error_handler ("base_dense_power_series< T >::operator[]",
				     "Exponent e exceeds this->last.");
		return T(0);
	}
        else if (e < this->first) {
		return T(0);
	}
        else {
		return (*this->coeff)[ e - this->first ];
	}
}



// Returns the coefficient with exponent e.
// This operator should only be used for
// writing - operations, because if the exponents e
// exceeds last or falls below first, new coefficients are allocated,
// initialized with zero and last or first is set to zero, respectively.

template< class T >
T &
base_dense_power_series< T >::operator() (lidia_size_t e)
{
        debug_handler ("base_dense_power_series< T >", "operator()");

        if (this->coeff->size() == 0) {
		this->assign_zero (e);
		return (*this->coeff)[0];
	}
        else if (e > this->last) {
		T zero_elem;
		lidia_size_t prec, i;
		lidia_size_t new_prec;

		zero_elem.assign_zero ();
		prec = this->last - this->first + 1;
		new_prec = prec + (e - this->last);

		this->allocate (new_prec, this->coeff);

		for (i = prec; i < new_prec; i++)
			(*this->coeff)[i] = zero_elem;

		this->last = e;

		return (*this->coeff)[ e - this->first ];
	}
        else if (e < this->first) {
		lidia_size_t delta, i;
		math_vector< T > *c;

		c = new math_vector< T >;
		delta = this->first - e;
		c->set_capacity (this->coeff->size() + delta);

		for (i = 0; i < delta; i++)
			(*c)[i].assign_zero ();

		c->assign (delta, *this->coeff, 0, this->coeff->size() -1);

		delete this->coeff;
		this->coeff = c;
		this->first = e;

		return (*this->coeff)[0];
	}
        else {
		return (*this->coeff)[ e - this->first ];
	}
}



// ************************************************
// *********** representation of 0 * x^f **********
// **(internal rep.: prec == 1, coeff[0] == 0, n == f) **
// ************************************************


template< class T >
void
base_dense_power_series< T >::assign_zero (lidia_size_t f)
{
	debug_handler ("base_dense_power_series< T >", "assign_zero(lidia_size_t)");

	this->first = f;
	this->last = f;
	this->coeff->set_capacity (1);
	(*this->coeff)[0].assign_zero ();
}



template< class T >
void
base_dense_power_series< T >::assign_one (lidia_size_t l)
{
	debug_handler ("base_dense_power_series< T >", "assign_one(lidia_size_t)");

	if (l < 0) {
		this->assign_zero (l);
	}
	else {
		lidia_size_t i, prec;

		prec = l+1;
		this->first = 0;
		this->last = l;
		this->coeff->set_capacity (prec);
		(*this->coeff)[0].assign_one ();

		for (i = 1; i < prec; i++)
			(*this->coeff)[i].assign_zero ();
	}
}



template< class T >
void
base_dense_power_series< T >::assign (const base_dense_power_series< T > & a)
{
	debug_handler ("base_dense_power_series< T >", "assign (const base_dense_power_series< T > &)");

	if (this != &a) {
		lidia_size_t non_zero_index;

		if (a.is_zero (non_zero_index)) {
			this->assign_zero (a.last);
		}
		else {
			this->coeff->set_capacity ((a.coeff)->size() - non_zero_index);
			this->coeff->assign (0, *(a.coeff), non_zero_index, (a.coeff)->size() -1);
			this->first = a.first + non_zero_index;
			this->last = a.last;
		}
	}
}



template< class T >
void
base_dense_power_series< T >::assign (const base_sparse_power_series< T > & a)
{
	debug_handler ("base_dense_power_series< T >", "assign (const base_sparse_power_series< T > &)");

	if (a.is_zero ()) {
		this->assign_zero (a.get_last());
	}
	else {
		lidia_size_t  prec, i;

		// initialize variables

		this->first = a.get_first ();
		this->last = a.get_last  ();


		// compute precision of the dense serie

		prec = this->last - this->first + 1;
		this->coeff->set_capacity (prec);


		// copy coefficients

		for (i = this->first; i <= this->last; i++) {
			(*this->coeff)[i-first] = a[i];
		}

		// check number of allocated elements

		if (this->max_num_coeff != 0 && prec > this->max_num_coeff)
			warning_handler ("base_dense_power_series< T >::assign(const base_sparse_power_series&)",
					 "Number of allocated coefficients exceeds user-defined upper bound.");
	}
}



// ************************************************
// *************** cast - operator ****************
// ************************************************

template< class T >
base_dense_power_series< T >::operator base_sparse_power_series< T > ()
{
	debug_handler ("base_dense_power_series< T >", "operator base_sparse_power_series< T > ()");
	return base_sparse_power_series< T > (*this->coeff, this->first);
}



// ************************************************
// **************** comparisons *******************
// ************************************************

template< class T >
bool
base_dense_power_series< T >::is_equal (const base_dense_power_series< T > & a) const
{
	debug_handler ("base_dense_power_series< T >", "is_equal (base_dense_power_series&)");

	if ((a.coeff)->size() == 0 || this->coeff->size() == 0)
		lidia_error_handler ("base_dense_power_series< T >::is_equal (const base_dense_power_series< T > &)",
				     "Serie is not initialized.");

	lidia_size_t non_zero_index_a;
	lidia_size_t non_zero_index_t;
	bool  a_zero, t_zero;
	bool  rc;

	t_zero = this->is_zero (non_zero_index_t);
	a_zero = a.is_zero (non_zero_index_a);

	if (a_zero != t_zero)
		rc = false;
	else if (a_zero && t_zero)
		if (this->last == a.last)
			rc = true;
		else
			rc = false;
	else {
		if (this->last != a.last || (this->first+non_zero_index_t) != (a.first+non_zero_index_a))
			rc = false;
		else {
			lidia_size_t prec_t;
			lidia_size_t i, j;

			prec_t = this->last - this->first + 1;
			rc = 1;

			for (i = non_zero_index_t, j = non_zero_index_a; i < prec_t && rc == true; i++, j++)
				if ((*this->coeff)[i] != (*(a.coeff))[j])
					rc = false;
		}
	}

	return rc;
}



//
// ***** input / output *****
//

template< class T >
int
base_dense_power_series< T >::read (std::istream & in)
{
	debug_handler ("base_dense_power_series< T >",
		       "read(std::istream&)");

	math_vector< T > *tmp;
	char c;
	int  rc;
	T zero_elem;

	in >> std::ws >> c;

	if (c != '[') {
		lidia_error_handler ("base_dense_power_series< T >::read(std::istream&)",
				     "Invalid input format ('[' expected) for base_dense_power_series< T > .");
	}
	else {
		in >> std::ws >> this->first;

		tmp = new math_vector< T >;

		in >> std::ws >> (*tmp);

		if (tmp->size() == 0) {
			delete tmp;
			lidia_error_handler ("base_dense_power_series< T >::read(std::istream&)",
					     "Input does not contain any coefficient.");
		}
		else {
			if ((*tmp)[0].is_zero()) {
				set (*tmp, this->first);
				delete tmp;
			}
			else {
				delete this->coeff;
				this->coeff = tmp;
				this->last = this->first + this->coeff->size() - 1;
			}

			in >> std::ws >> c;

			if (c != ']') {
				lidia_error_handler ("base_dense_power_series< T >::read(std::istream&)",
						     "Invalid input format (']' expected) for base_dense_power_series< T > .");
			}
		}
	}

	for (rc = 0; rc < 0; )
		;
	return (rc);
}



template< class T >
void
base_dense_power_series< T >::write (std::ostream & out) const
{
	debug_handler ("base_dense_power_series", "write(std::ostream&)");
	lidia_size_t non_zero_index;
	lidia_size_t i;

	if (this->is_zero (non_zero_index)) {
		out << "[ " << this->last << " " << "[ " << (*this->coeff)[0] << " ] " << " ]";
	}
	else {
		out << "[ " << this->first+non_zero_index << " " << "[ ";

		for (i = non_zero_index; i< this->coeff->size(); i++)
			out << (*this->coeff)[i] << " ";

		out << "] ]";
	}
}



//
// ***** miscellaneous *****
//

template< class T >
void
base_dense_power_series< T >::multiply_by_xn (lidia_size_t n)
{
	debug_handler ("base_dense_power_series< T >",
                       "multiply_by_xn(lidia_size_t)");
	this->first += n;
	this->last += n;
}



template< class T >
void
base_dense_power_series< T >::compose (lidia_size_t n)
{
	debug_handler ("base_dense_power_series< T >",
                       "compose(lidia_size_t)");
	lidia_size_t prec;
	lidia_size_t i, j;

	T zero_elem;


	prec = this->coeff->size();

	this->coeff->set_capacity (prec*n);

	this->first *= n;
	this->last *= n;

	for (i = prec-1; i >= 0; i--) {
		(*this->coeff)[i*n] = (*this->coeff)[i];
	}

	zero_elem.assign_zero ();

	for (i = 0; i < prec; i++)
		for (j = 1; j < n; j++) {
			(*this->coeff)[j+i*n] = zero_elem;
		}
}



#if 0
template< class T >
void
randomize (base_dense_power_series< T > & a, lidia_size_t f, lidia_size_t l, const T & coeff_bound)
{
	debug_handler ("base_dense_power_series< T >",
                       "randomize(base_dense_power_series< T > &, lidia_size_t, lidia_size_t)");

	lidia_size_t prec = (l - f + 1);
	lidia_size_t i;

	T tmp;

	math_vector< T > * A = a.coeff;

	if (prec <= 0) {
		lidia_error_handler ("base_dense_power_series< T >::randomize(base_dense_power_series&, lidia_size_t, lidia_size_t, T&)",
				     "Non positive number of coefficients!");
        }
	else {
		a.first = f;
		a.last = l;

		A->set_capacity (prec);

		for (i = 0; i < prec; i++) {
			(*A)[i] = randomize (coeff_bound);
		}

		if ((*A)[0].is_zero ()) {
			(*A)[0].assign_one ();
		}
        }
}
#endif



//
// ***** arithmetic via friend functions *****
//

template< class T >
void
base_dense_power_series< T >::add (const base_dense_power_series< T > & a,
				   const base_dense_power_series< T > & b)
{
	debug_handler ("base_dense_power_series< T >",
                       "add(base_dense_power_series< T > &, base_dense_power_series< T > &, base_dense_power_series< T > &");

	lidia_size_t na, nb, nres;
	lidia_size_t pa, pb, pres;
	lidia_size_t l , u , i;

	int  end, ident;

	math_vector< T > *res;
	math_vector< T > *X, *Y;
	T  z;

	na = a.first;
	nb = b.first;
	pa = (a.coeff)->size();
	pb = (b.coeff)->size();

	// check for invalid arguments

	if (pa == 0 || pb == 0)
        {
		lidia_error_handler ("base_dense_power_series< T >::add(base_dense_power_series< T > &, base_dense_power_series< T >, base_dense_power_series< T > &",
				     "Arguments not initialized.");
        }

	// coefficients don't overlap ->the sum is given by the series with
	// the smaller start exponent

	else if (na + pa <= nb)
        {
		*this = a;
        }
	else if (nb + pb <= na)
        {
		*this = b;
        }

	// the coefficients overlap ->compute sum (stored in res)

	else
        {
		ident = 0;

		// check whether input alias output

		if ((this != &a) && (this != &b)) {
			res = this->coeff;
		}
		else {
			ident = 1;
			res = new math_vector< T >;
		}

		// determine the series with the smaller start exponent

		nres = na < nb ? na : nb;

		if (nres == na) {
			l = nb - na;
			pres = (l + ((pa-l < pb) ? pa-l : pb));
			X = a.coeff;
			Y = b.coeff;
		}
		else {
			l = na - nb;
			pres = (l + ((pb-l < pa) ? pb-l : pa));
			X = b.coeff;
			Y = a.coeff;
		}

		// different start exponents

		if (l > 0) {
			// set length of the vector of coefficients of res

			res->set_capacity (pres);

			// The first coefficients of the result are given by
			// the coefficients of the series with the smaller
			// start exponent up to the point the overlapping
			// of the coefficients begins.

			for (i = 0; i < l; i++) {
				res->assign (0, *X, 0, l-1);
			}

			// From the beginning of the overlapping up to the length
			// of the result (pres) the coefficients are added.

			for (i = l; i < pres; i++) {
				LiDIA::add ((*res)[i], (*X)[i], (*Y)[i-l]);
			}
		}

		// same start exponents

		else {
			pres = pa < pb ? pa : pb;
			X = a.coeff;
			Y = b.coeff;

			i = end = 0;

			// ignore leading zero

			LiDIA::add (z, (*X)[0], (*Y)[0]);

			while ((z == 0) && !end) {
				if (i == pres - 1) {
					end = 1;
				}
				else {
					i++;
					LiDIA::add (z, (*X)[i], (*Y)[i]);
				}
			}

			// compute result

			nres = na + i;

			l = i;
			u = pres;
			pres = pres - i;


			// set length of the vector of coefficients of res

			res->set_capacity (pres);


			// add coefficients

			for (i = l; i < u; i++) {
				LiDIA::add ((*res)[i-l], (*X)[i], (*Y)[i]);
			}

		}  // end else - same start exponent

		this->first = nres;
		this->last = nres + pres - 1;

		if (ident) {
			delete this->coeff;
			this->coeff = res;
		}

	} // end else - overlapping
}



template< class T >
void
base_dense_power_series< T >::add (const base_dense_power_series< T > & a,
				   const T & b)
{
	debug_handler ("base_dense_power_series< T >",
                       "add(base_dense_power_series< T > &, base_dense_power_series< T > &, T&)");

	math_vector< T > *A = a.coeff;
	T x;

	// check for invalid arguments

	if (A->size() == 0)
        {
		lidia_error_handler ("base_dense_power_series< T >::add(base_dense_power_series< T > &, base_dense_power_series< T > &, T&)",
				     "Argument not initialized.");
        }


	// add a and b * x^0 + O(x^a.last)

	else if (a.last >= 0)
        {
		if (a.first > 0)
			x = b;
		else
			x = (*A)[-a.first] + b;

		set_coeff (x, 0, a);
        }
	else
		*this = a;
}



template< class T >
void
base_dense_power_series< T >::add (const T & b,
				   const base_dense_power_series< T > & a)
{
	debug_handler ("base_dense_power_series< T >",
                       "add(base_dense_power_series< T > &, T&, base_dense_power_series< T > &)");

	math_vector< T > *A = a.coeff;
	T x;


	// check for invalid arguments

	if (A->size() == 0)
        {
		lidia_error_handler ("base_dense_power_series< T >::add(base_dense_power_series< T > &, T&, base_dense_power_series< T > &)",
				     "Argument not initialized.");
        }


	// add a and b * x^0 + O(x^a.last)

	else if (a.last >= 0)
        {
		if (a.first > 0)
			x = b;
		else
			x = b + (*A)[-a.first];

		set_coeff (x, 0, a);
        }
	else
		*this = a;
}



template< class T >
void
base_dense_power_series< T >::subtract (const base_dense_power_series< T > & a,
					const base_dense_power_series< T > & b)
{
	debug_handler ("base_dense_power_series< T >",
                       "subtract(base_dense_power_series< T > &, base_dense_power_series< T >, base_dense_power_series< T > &");

	lidia_size_t na, nb, nres;
	lidia_size_t pa, pb, pres;
	lidia_size_t l, u, i;
	int end, ident;

	math_vector< T > *res;
	math_vector< T > *X, *Y;
	T  z;

	na = a.first;
	nb = b.first;
	pa = (a.coeff)->size();
	pb = (b.coeff)->size();


	// check for invalid arguments

	if (pa == 0 || pb == 0)
        {
		lidia_error_handler ("base_dense_power_series< T >::subtract(base_dense_power_series< T > &, base_dense_power_series< T >, base_dense_power_series< T > &",
				     "Arguments not initialized.");
        }

	// coefficients don't overlap ->the result is given by the series with
	// the smaller start exponent

	else if (na + pa <= nb)
        {
		*this = a;
        }
	else if (nb + pb <= na)
        {
		*this = b;
        }

	// the coefficients overlap ->compute difference (stored in res)

	else
        {
		ident = 0;

		// check whether input alias output

		if ((this != &a) && (this != &b)) {
			res = this->coeff;
		}
		else {
			ident = 1;
			res = new math_vector< T >;
		}

		// determine the series with the smaller start exponent

		nres = na < nb ? na : nb;

		if (nres == na) {
			l = nb - na;
			pres = (l + ((pa-l < pb) ? pa-l : pb));
			X = a.coeff;
			Y = b.coeff;
		}
		else {
			l = na - nb;
			pres = (l + ((pb-l < pa) ? pb-l : pa));
			X = b.coeff;
			Y = a.coeff;
		}

		// different start exponents

		if (l > 0) {
			// set length of the vector of coefficients of res

			res->set_capacity (pres);

			// The first coefficients of the result are given by
			// the coefficients of the series with the smaller
			// start exponent up to the point the overlapping
			// of the coefficients begins.

			if (X == a.coeff) {

				for (i = 0; i < l; i++) {
					res->assign (0, *X, 0, l-1);
				}
			}
			else {
				for (i = 0; i < l; i++) {
					LiDIA::negate ((*res)[i], (*X)[i]);
				}
			}

			// From the beginning of the overlapping up to the length
			// of the result (pres) the coefficients are subtracted.

			if (X == a.coeff) {
				for (i = l; i < pres; i++) {
					LiDIA::subtract ((*res)[i], (*X)[i], (*Y)[i-l]);
				}
			}
			else {
				for (i = l; i < pres; i++) {
					LiDIA::subtract ((*res)[i], (*Y)[i-l], (*X)[i]);
				}
			}
		}

		// same start exponents

		else {
			pres = pa < pb ? pa : pb;
			X = a.coeff;
			Y = b.coeff;

			i = end = 0;

			// ignore leading zero

			if (X == a.coeff) {
				LiDIA::subtract (z, (*X)[0], (*Y)[0]);

				while ((z == 0) && !end) {
					if (i == pres - 1) {
						end = 1;
					}
					else {
						i++;
						LiDIA::subtract (z, (*X)[i], (*Y)[i]);
					}
				}
			}
			else {
				LiDIA::subtract (z, (*Y)[0], (*X)[0]);

				while ((z == 0) && !end) {
					if (i == pres - 1) {
						end = 1;
					}
					else {
						i++;
						LiDIA::subtract (z, (*Y)[i], (*X)[i]);
					}
				}
			}


			// compute result

			nres = na + i;

			l = i;
			u = pres;
			pres = pres - i;

			// set length of the vector of coefficients of res

			res->set_capacity (pres);

			// subtract coefficients

			if (X == a.coeff) {
				for (i = l; i < u; i++) {
					LiDIA::subtract ((*res)[i-l], (*X)[i], (*Y)[i]);
				}
			}
			else {
				for (i = l; i < u; i++) {
					LiDIA::subtract ((*res)[i-l], (*Y)[i], (*X)[i]);
				}
			}

		}  // end else - same start exponents

		this->first = nres;
		this->last = nres + pres - 1;

		if (ident) {
			delete this->coeff;
			this->coeff = res;
		}

	} // end else - overlapping
}



template< class T >
void
base_dense_power_series< T >::subtract (const base_dense_power_series< T > & a,
					const T & b)
{
	debug_handler ("base_dense_power_series< T >",
                       "subtract(base_dense_power_series< T > &, base_dense_power_series< T > &, T&)");

	math_vector< T > *A = a.coeff;
	T x;


	// check for invalid arguments

	if (A->size() == 0) {
		lidia_error_handler ("base_dense_power_series< T >::subtract(base_dense_power_series< T > &, base_dense_power_series< T > &, T&)",
				     "Argument not initialized.");
        }


	// subtract a and b * x^0 + O(x^a.last)

	else if (a.last >= 0) {
		if (a.first > 0)
			LiDIA::negate (x, b);
		else
			LiDIA::subtract (x, (*A)[-a.first], b);

		set_coeff (x, 0, a);
        }
	else
		*this = a;
}



template< class T >
void
base_dense_power_series< T >::subtract (const T & b,
					const base_dense_power_series< T > & a)
{
	debug_handler ("base_dense_power_series< T >",
                       "subtract(base_dense_power_series< T > &, T&, base_dense_power_series< T > &)");

	math_vector< T > *A = a.coeff;
	T x;


	// check for invalid arguments

	if (A->size() == 0) {
		lidia_error_handler ("base_dense_power_series< T >::subtract(base_dense_power_series< T > &, T&, base_dense_power_series< T > &)",
				     "Argument not initialized.");
        }


	// subtract b * x^0 + O(x^a.last) and a

	else if (a.last >= 0) {
		if (a.first > 0)
			LiDIA::negate (x, b);
		else
			LiDIA::subtract (x, (*A)[-a.first], b);

		set_coeff (x, 0, a);
		this->negate (*this);
        }
	else
		this->negate(a);
}



template< class T >
void
base_dense_power_series< T >::multiply (const base_dense_power_series< T > & a,
					const T & b)
{
	debug_handler ("base_dense_power_series< T >",
                       "multiply(base_dense_power_series< T > &, base_dense_power_series< T > &, T&)");

	math_vector< T > *A = a.coeff;
	math_vector< T > *C = this->coeff;

	lidia_size_t pa = A->size();
	lidia_size_t i;


	// check for invalid arguments

	if (pa == 0) {
		lidia_error_handler ("base_dense_power_series< T >::multiply(base_dense_power_series< T > &, base_dense_power_series< T > &, T&)",
				     "Argument not initialized.");
        }


	// multiply with b

	if (b.is_zero()) {
		C->set_capacity (1);
		(*C)[0].assign_zero ();
		this->first = a.last;
		this->last = a.last;
        }
	else {
		C->set_capacity (pa);

		for (i = 0; i < pa; i++)
			LiDIA::multiply ((*C)[i], (*A)[i], b);

		this->first = a.first;
		this->last = a.last;
        }
}



template< class T >
void
base_dense_power_series< T >::multiply (const T & b,
					const base_dense_power_series< T > & a)
{
	debug_handler ("base_dense_power_series< T >",
                       "multiply(base_dense_power_series< T > &, T&, base_dense_power_series< T > &)");

	math_vector< T > *A = a.coeff;
	math_vector< T > *C = this->coeff;

	lidia_size_t pa = A->size();
	lidia_size_t i;


	// check for invalid arguments

	if (pa == 0) {
		lidia_error_handler ("base_dense_power_series< T >::multiply(base_dense_power_series< T > &, T&, base_dense_power_series< T > &)",
				     "Argument not initialized.");
        }


	// multiply with b

	if (b.is_zero ()) {
		C->set_capacity (1);
		(*C)[0].assign_zero ();
		this->first = a.last;
		this->last = a.last;
        }
	else {
		C->set_capacity (pa);

		for (i = 0; i < pa; i++)
			LiDIA::multiply ((*C)[i], b, (*A)[i]);

		this->first = a.first;
		this->last = a.last;
        }
}



template< class T >
void
base_dense_power_series< T >::divide (const base_dense_power_series< T > & a,
				      const T & b)
{
	debug_handler ("base_dense_power_series< T >",
                       "divide(base_dense_power_series< T > &, base_dense_power_series< T > &, T&)");
	T x;

	invert (x, b);
	this->multiply (a, x);
}



template< class T >
void
base_dense_power_series< T >::negate (const base_dense_power_series< T > & a)
{
	debug_handler ("base_dense_power_series< T >",
                       "negate(base_dense_power_series< T > &, base_dense_power_series< T > &)");

	math_vector< T > *C = this->coeff;
	math_vector< T > *A = a.coeff;

	lidia_size_t pa = A->size();
	lidia_size_t i;


	// check for invalid argument

	if (pa == 0) {
		lidia_error_handler ("base_dense_power_series< T >::negate(base_dense_power_series< T > &, base_dense_power_series< T > &)",
				     "Argument not initialized.");
        }


	// negate coefficients

	C->set_capacity (pa);

	for (i = 0; i < pa; i++)
		LiDIA::negate ((*C)[i], (*A)[i]);

	this->first = a.first;
	this->last = a.last;
}



// ************************************************
// *************** miscellaneous ******************
// ************************************************

template< class T >
void
base_dense_power_series< T >::swap (base_dense_power_series< T > & a)
{
	debug_handler ("base_dense_power_series< T >",
                       "swap(base_dense_power_series< T > &, const base_dense_power_series< T > &)");

	if (this != &a) {
		lidia_size_t      tmp_lidia_size_t;
		math_vector< T > *tmp_base;

		tmp_lidia_size_t = this->first;
		this->first = a.first;
		a.first = tmp_lidia_size_t;

		tmp_lidia_size_t = this->last;
		this->last = a.last;
		a.last = tmp_lidia_size_t;

		tmp_base = this->coeff;
		this->coeff = a.coeff;
		a.coeff = tmp_base;
        }
}



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_BASE_DENSE_POWER_SERIES_CC_GUARD_
