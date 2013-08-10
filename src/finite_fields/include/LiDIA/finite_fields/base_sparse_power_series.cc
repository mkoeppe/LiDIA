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


#ifndef LIDIA_BASE_SPARSE_POWER_SERIES_CC_GUARD_
#define LIDIA_BASE_SPARSE_POWER_SERIES_CC_GUARD_


#ifndef LIDIA_BASE_SPARSE_POWER_SERIES_H_GUARD_
# include	"LiDIA/finite_fields/base_sparse_power_series.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



//
// ****  constructor/destructor functions    ******
//

template< class T >
base_sparse_power_series< T >::base_sparse_power_series ()
{
	debug_handler ("base_sparse_power_series< T >",
		       "base_sparse_power_series()");

	this->coeff = new sort_vector< spc< T > >;
	this->coeff->set_sort_direction (vector_flags::sort_vector_up);

	this->sorted = true;
	this->input_verification = true;
}



template< class T >
base_sparse_power_series< T >::base_sparse_power_series (const T & elem, lidia_size_t l)
{
	debug_handler ("base_sparse_power_series< T >",
		       "base_sparse_power_series(const T&, lidia_size_t)");

	this->coeff = new sort_vector< spc< T > >;
	this->coeff->set_sort_direction (vector_flags::sort_vector_up);

	if (! elem.is_zero ()) {
		this->coeff->set_capacity(1);
		(*this->coeff)[0].coeff = elem;
		(*this->coeff)[0].exp = 0;

		this->first = 0;
		this->last = l;

		this->sorted = true;
	}
	else {
		this->assign_zero (l);
	}
	this->input_verification = true;
}



template< class T >
base_sparse_power_series< T >::base_sparse_power_series (const base_vector< T > & c, lidia_size_t f)
{
	debug_handler ("base_sparse_power_series< T >",
		       "base_sparse_power_series(const base_vector< T > &, lidia_size_t, lidia_size_t)");

	this->coeff = new sort_vector< spc< T > >;
	this->coeff->set_sort_direction (vector_flags::sort_vector_up);
	this->sorted = true;

	lidia_size_t prec, start;
	lidia_size_t num, i;
	T   zero_T;

	zero_T.assign_zero ();
	prec = c.size();
	start = prec;
	num = 1;


	// find first element in c which is not equal to zero
	// and count the number of non-zero elements

	for (i = 0; i < prec && start == prec; i++) {
		if (c[i] != zero_T)
			start = i;
	}

	for (i = start+1; i < prec; i++) {
		if (c[i] != zero_T)
			num ++;
	}


	// all elements are zero

	if (start == prec) {
		this->assign_zero (f + prec - 1);
	}

	// a non-zero element was found

	else {
		this->coeff->set_capacity (num);
		this->first = f + start;
		this->last = f + prec - 1;

		for (i = start, num = 0; i < prec; i++) {
			if (c[i] != zero_T) {
				(*this->coeff)[num].coeff = c[i];
				(*this->coeff)[num].exp = f + i;
				num ++;
			}
		}
	}

	this->input_verification = true;
}



template< class T >
base_sparse_power_series< T >::base_sparse_power_series (const base_sparse_power_series< T > & a)
{
	debug_handler ("base_sparse_power_series< T >",
		       "base_sparse_power_series(const base_sparse_power_series< T > &)");

	a.init_test ("base_sparse_power_series(const base_sparse_power_series< T > &)");

	this->coeff = new sort_vector< spc< T > >;
	this->coeff->set_sort_direction (vector_flags::sort_vector_up);

	this->input_verification = a.input_verification;

	if (a.is_zero()) {
		this->assign_zero (a.last);
	}
	else {
		this->first = a.first;
		this->last = a.last;
		*this->coeff = *a.coeff;

		this->sorted = a.sorted;
	}
}



template< class T >
base_sparse_power_series< T >::~base_sparse_power_series ()
{
	debug_handler ("base_sparse_power_series< T >",
		       "~base_sparse_power_series()");

	delete this->coeff;
}



//
//  *****  initialization functions  *****
//

template< class T >
bool
base_sparse_power_series< T >::is_zero () const
{
	debug_handler ("base_sparse_power_series< T >",
		       "is_zero()");

	init_test ("is_zero()");

	if (this->coeff->size() == 0)
		return true;
	else
		return false;
}



template< class T >
bool
base_sparse_power_series< T >::is_one () const
{
	debug_handler ("base_sparse_power_series< T >",
		       "is_one()");
	bool rc;
	init_test ("is_one()");

	if ((this->coeff->size() != static_cast<lidia_size_t>(1)) ||
	    (this->first != static_cast<lidia_size_t>(0)))
		rc = false;
	else {
		if ((*this->coeff)[0].coeff.is_one())
			rc = true;
		else
			rc = false;
	}

	return rc;
}



//
// assigners
//

template< class T >
void
base_sparse_power_series< T >::assign_zero (lidia_size_t l)
{
	debug_handler ("base_sparse_power_series< T >",
		       "assign_zero(lidia_size_t)");

	this->coeff->set_capacity (1);
	this->coeff->set_size     (0);

	this->first = this->last = l;
	this->sorted = true;
}



template< class T >
void
base_sparse_power_series< T >::assign_one (lidia_size_t l)
{
	debug_handler ("base_sparse_power_series< T >",
		       "assign_one(lidia_size_t)");

	if (l < 0) {
		this->assign_zero (l);
	}
	else {
		this->coeff->set_capacity (1);
		this->coeff->set_size     (1);

		((*this->coeff)[0].coeff).assign_one ();
		(*this->coeff)[0].exp = 0;

		this->first = 0;
		this->last = l;
		this->sorted = true;
	}
}



template< class T >
void
base_sparse_power_series< T >::set (const T & elem, lidia_size_t l)
{
	debug_handler ("base_sparse_power_series< T >",
		       "set (const T &, lidia_size_t)");

	if (! elem.is_zero ()) {
		this->coeff->set_capacity(1);
		(*this->coeff)[0].coeff = elem;
		(*this->coeff)[0].exp = 0;

		this->first = 0;
		this->last = l;

		this->sorted = true;
	}
	else {
		this->assign_zero (l);
	}
}



template< class T >
void
base_sparse_power_series< T >::set (const T * c, const lidia_size_t * e, lidia_size_t num, lidia_size_t l)
{
	debug_handler ("base_sparse_power_series< T >",
		       "set(const T*, const lidia_size_t*, lidia_size_t)");

	lidia_size_t i, j;

	this->coeff->set_capacity (num);

	// --- copying the two vectors into the spc< . > structure

	for (i = 0, j = 0; i < num; i++) {
		if (c[i] != T(0)) {
			(*this->coeff)[j].coeff = c[i];
			(*this->coeff)[j].exp = e[i];
			j ++;
		}

		if (l < e[i])
			lidia_error_handler ("base_sparse_power_series< T >",
					     "set(const T*, ...)::Invalid value for l.");
	}

	if (j == 0) {
		this->assign_zero (l);
	}
	else {
		this->coeff->sort();
		this->first = (*this->coeff)[0].exp;
		this->last = l;
		this->sorted = true;
	}
}



#if 0
template< class T >
void
base_sparse_power_series< T >::set (const base_vector< T > & c,
				    const base_vector< lidia_size_t > & e,
				    lidia_size_t l)
{
	debug_handler ("base_sparse_power_series< T >",
		       "set(const base_vector< T > &, const base_vector< T > &, lidia_size_t)");

	if (c.size() != e.size()) {
		lidia_error_handler ("base_sparse_power_series< T >",
				     "set(2x base_vector< T >, lidia_size_t)::incompatible vector sizes");
	}

	lidia_size_t i, j, num;

	num = c.size();
	this->coeff->set_capacity (num);

	// --- copying the two vectors into the spc< . > structure

	for (i = 0, j = 0; i < num; i++) {
		if (c[i] != T(0)) {
			(*this->coeff)[j].coeff = c[i];
			(*this->coeff)[j].exp = e[i];

			j ++;
		}

		if (l < e[i])
			lidia_error_handler ("base_sparse_power_series< T >",
					     "set(const base_vector< T > &, ...)::Invalid value for l.");
	}

	if (j == 0) {
		this->assign_zero (l);
	}
	else { 
		this->coeff->sort();

		this->first = (*this->coeff)[0].exp;
		this->last = l;
		this->sorted = true;
	}
}
#endif


template< class T >
void
base_sparse_power_series< T >::set_coeff (const T & elem, lidia_size_t exp)
{
	debug_handler ("base_sparse_power_series< T >",
		       "set_coeff (const T &, lidia_size_t)");

	lidia_size_t where, sz;
	bool el_is_zero;
	bool found;
	spc< T > c(elem, exp);

	sz = this->coeff->capacity(); // i.e. not initialized
	el_is_zero = (elem.is_zero ());

	if (sz == 0) {
		// series not yet initialized
		if (el_is_zero) {
			this->assign_zero (exp);
		}
		else {
			this->coeff->set_capacity (1);
			this->first = exp;
			this->last = exp;
			(*this->coeff)[0] = c;

			this->sorted = true;
		}
	}
	else if (this->is_zero ()) {
		if (el_is_zero && exp > this->last) {
			this->assign_zero (exp);
		}
		else {
			this->coeff->set_capacity (1);
			this->first = exp;
			if (exp > this->last)
				this->last = exp;
			(*this->coeff)[0] = c;

			this->sorted = true;
		}
	}
	else {
		// update an already initialized non-zero serie
		if (this->sorted)
			found = this->coeff->bin_search (c, where);
		else
			found = this->coeff->linear_search (c, where);

		if (found) {
			if (el_is_zero) {
				// delete zero-element by overwriting it
				// if ((where < sz-1)) {
				//	  (*this->coeff)[where] = (*coeff)[sz-1];
				//    if ((*coeff)[where+1] < (*coeff)[where])
				//      sorted = false;
				//  }
				// coeff->set_size (sz - 1);
				// first = ???

				this->coeff->remove_from (where);

				if (this->coeff->size() == 0)
					this->assign_zero (this->last);
				else
					this->first = (*this->coeff)[0].exp;
			}
			else {
				// overwrite coeff with new value
				(*this->coeff)[where] = c;
			}
		}
		else if (! el_is_zero) {
			// put element at the end of coeff.-vector or
			// insert it at corresponding position 'where'
			if (this->coeff->capacity() < sz+1)
				this->coeff->set_capacity (alloc_size (sz + 1));

			// (*coeff)[sz] = c;
			// if (sorted && (*coeff)[sz-1].exp > exp)
			//   sorted = false;

			this->coeff->insert_at (c, where);
			this->first = (*this->coeff)[0].exp;
		}

		if (exp > this->last) this->last = exp;
	}
}



template< class T >
void
base_sparse_power_series< T >::reduce_last (lidia_size_t l)
{
	debug_handler ("base_sparse_power_series< T >", "reduce_last (lidia_size_t)");

	init_test ("reduce_last (lidia_size_t)");
	sort_test ();

	lidia_size_t where;
	bool found;
	spc< T > c((*this->coeff)[0].coeff, l);

	if (l > this->last) {
		lidia_error_handler("base_sparse_power_series< T >", "reduce_last(lidia_size_t)::new value exceeded 'last'");
	}
	else if (l < this->first) {
		this->assign_zero (l);
	}
	else {
		found = this->coeff->bin_search (c, where);

		where += (found ? 1 : 0);
		this->coeff->set_size (where);
		this->last = l;
	}

}



template< class T >
void
base_sparse_power_series< T >::clear ()
{
	debug_handler ("base_sparse_power_series< T >", "clear ()");

	this->coeff->kill ();
	this->sorted = false;
}



template< class T >
void
base_sparse_power_series< T >::normalize ()
{
	debug_handler ("base_sparse_power_series< T >", "normalize ()");

	init_test ("normalize()");
	sort_test ();

	if (! this->is_zero ()) this->coeff->set_capacity (this->coeff->size ());
}



//
//  *****  information functions  *****
//

template< class T >
lidia_size_t
base_sparse_power_series< T >::get_first () const
{
	debug_handler ("base_sparse_power_series< T >", "get_first ()");

	init_test ("get_first ()");
	return this->first;
}



template< class T >
lidia_size_t
base_sparse_power_series< T >::get_last () const
{
	debug_handler ("base_sparse_power_series< T >", "get_last ()");

	init_test ("get_last ()");
	return this->last;
}



template< class T >
void
base_sparse_power_series< T >::get (T * & c, lidia_size_t * & e, lidia_size_t & sz) const
{
	debug_handler ("base_sparse_power_series< T >",
		       "get(T* &, lidia_size_t* &, lidia_size_t&)");

	lidia_size_t i;

	init_test ("get(T* &, T* &, lidia_size_t&)");
	sort_test ();
	sz = this->coeff->size();

	if (c != NULL) delete [] c;
	if (e != NULL) delete [] e;


	if (sz != static_cast<lidia_size_t>(0)) {
		c = new T[sz];
		e = new lidia_size_t[sz];
	}
	else {
		c = new T[1];
		e = new lidia_size_t[1];
	}

	memory_handler (c, "base_sparse_power_series< T >",
			"get(T*, lidia_size_t*, lidia_size_t)::out of memory(c)");
	memory_handler (e, "base_sparse_power_series< T >",
			"get(T*, lidia_size_t*, lidia_size_t)::out of memory(e)");

	for (i = 0; i < sz; i++) {
		c[i] = (*this->coeff)[i].coeff;
		e[i] = (*this->coeff)[i].exp;
	}

	if (sz == static_cast<lidia_size_t>(0)) {
		c[0].assign_zero();
		e[0] = this->last;
	}
}



#if 0
template< class T >
void
base_sparse_power_series< T >::get (base_vector< T > & c, base_vector< lidia_size_t > & e) const
{
	debug_handler ("base_sparse_power_series< T >",
		       "get (base_vector< T > &, base_vector< lidia_size_t > &)");

	init_test ("get (base_vector< T > &, base_vector< T > &)");
	sort_test ();

	lidia_size_t i, sz = this->coeff->size();

	if (sz != static_cast<lidia_size_t>(0)) { 
		c.set_capacity (sz);
		c.set_size     (sz);
		e.set_capacity (sz);
		e.set_size     (sz);

		for (i = 0; i < sz; i++) { 
			c[i] = (*this->coeff)[i].coeff;
			e[i] = (*this->coeff)[i].exp;
		}
	}
	else { 
		sz++;
		c.set_capacity (sz);
		c.set_size     (sz);
		e.set_capacity (sz);
		e.set_size     (sz);
		c[0].assign_zero();
		e[0] = this->last;
	}
}
#endif


template< class T >
void
base_sparse_power_series< T >::get_coeff (T & elem, lidia_size_t exp) const
{
	debug_handler ("base_sparse_power_series< T >", "get_coeff (T &, lidia_size_t)");

	init_test ("get_coeff (T &, lidia_size_t &)");
	sort_test ();

	if (exp > this->last)
		lidia_error_handler ("base_sparse_power_series< T >",
				     "get_coeff(T &, lidia_size_t)::exponent exceeded 'last'");

	bool found;
	lidia_size_t where;
	spc< T > c(elem, exp);

	found = this->coeff->bin_search (c, where);

	if (found)
		elem = (*this->coeff)[where].coeff;
	else
		elem.assign_zero ();
}



//
//  *****  subscription operator  *****
//

template< class T >
const T
base_sparse_power_series< T >::operator [] (lidia_size_t exp) const
{
	debug_handler ("base_sparse_power_series< T >", "operator [] (lidia_size_t)");

	init_test ("operator [] (lidia_size_t)");
	sort_test ();

	if (exp > this->last)
		lidia_error_handler ("base_sparse_power_series< T >",
				     "operator [] (lidia_size_t)::exponent exceeded 'last'");

	if (this->is_zero())
		return T(0);

	bool found;
	lidia_size_t where;
	T elem;
	spc< T > c(elem, exp);

	found = this->coeff->bin_search (c, where);

	if (found)
		return (*this->coeff)[where].coeff;
	else
		return T(0);
}



//
//  *****  I/O functions and operators  *****
//

template< class T >
void
base_sparse_power_series< T >
::set_input_verification(bool v)
{
	debug_handler ("base_sparse_power_series< T >", "set_input_verification(bool)");
	this->input_verification = v;
}



template< class T >
int
base_sparse_power_series< T >::read (std::istream & in)
{
	debug_handler ("base_sparse_power_series< T >", "read(std::istream&)");

	int  rc;
	bool found;
	lidia_size_t pos;
	char ch;

	spc< T > zero_elem;

	(zero_elem.coeff).assign_zero ();
	zero_elem.exp = 0;

	in >> ch;
	if (ch != '[')
		lidia_error_handler ("base_sparse_power_series< T >", "read():: character '[' expected");

	in >> *this->coeff;
	in >> this->last;

	in >> ch;
	if (ch != ']')
		lidia_error_handler ("base_sparse_power_series< T >", "read():: character ']' expected");

	if (this->coeff->capacity () == 0) this->assign_zero(last);


	if (this->input_verification) {
		//  ---  delete any zero element in the input  ---

		this->coeff->set_sort_direction (coeff_cmp_zero);
		this->coeff->sort ();

		found = this->coeff->bin_search (zero_elem, pos);

		if (found) {
			this->coeff->set_size (pos);
		}

		this->coeff->set_sort_direction (vector_flags::sort_vector_up);

		rebuild();
		normalize();

		if (this->last< (*this->coeff)[this->coeff->size()-1].exp)
			lidia_error_handler ("base_sparse_power_series",
					     "read::Coefficient with exponent greater than this->last was found.");
	}
	else {
		this->coeff->set_sort_direction (vector_flags::sort_vector_up);
		this->coeff->set_capacity (this->coeff->size ());
		this->sorted = true;

		if (this->coeff->size() == 0)
			this->assign_zero (this->last);
		else
			this->first = (*this->coeff)[0].exp;
	}

	rc = 0;
	return (rc);
}



template< class T >
int
base_sparse_power_series< T >::write (std::ostream & out) const
{
	debug_handler ("base_sparse_power_series< T >", "write(std::ostream&)");

	int rc;
	for (rc = 0; rc < 0; )
		;

	T zero_T(0);

	if (this->is_zero())
		out << "[ " << zero_T << " " << this->last << " ]";
	else
		out << "[ " << (*this->coeff) << " " << this->last << " ]";

	return (rc);
}



template< class T >
void
base_sparse_power_series< T >::info (std::ostream & out) const
{
	debug_handler ("base_sparse_power_series< T >", "info (std::ostream&)");

	out << "{ first = " << this->first << " last = " << this->last;
	out << " size = " << this->coeff->size();
	out << " cap = " << this->coeff->capacity();
	out << (this->sorted ? " TRUE" : " FALSE") << " }";
}



//
//  *****  assignment operator  *****
//

template< class T >
void
base_sparse_power_series< T >::assign (const base_sparse_power_series< T > & a)
{
	debug_handler ("base_sparse_power_series< T >", "operator = (const sp_pow &)");

	a.init_test ("operator = (const sp_pow &)");

	if (this != &a) {
		this->input_verification = a.input_verification;

		if (a.is_zero ())
			this->assign_zero (a.last);
		else {
			this->first = a.first;
			this->last = a.last;
			this->sorted = a.sorted;

			*this->coeff = *a.coeff;
		}
	}
}



//
//  ***** swap - function *****
//

template< class T >
void
base_sparse_power_series< T >::swap (base_sparse_power_series< T > & b)
{
	debug_handler ("base_sparse_power_series< T >", "swap (sp_pow &, sp_pow &)");

	lidia_size_t f, l;
	bool s;
	bool v;
	sort_vector< spc< T > > *c;

	f = this->first;
	l = this->last;
	s = this->sorted;
	c = this->coeff;
	v = this->input_verification;

	this->first = b.first;
	this->last = b.last;
	this->sorted = b.sorted;
	this->coeff = b.coeff;
	this->input_verification = b.input_verification;

	b.first = f;
	b.last = l;
	b.sorted = s;
	b.coeff = c;
	b.input_verification = v;
}



//
//  *****  comparison of series *****
//

template< class T >
bool
base_sparse_power_series< T >::is_equal (const base_sparse_power_series< T > & b) const
{
	debug_handler ("base_sparse_power_series< T >", "is_equal (const sp_pow &)");

	init_test ("is_equal (const sp_pow &)");
	b.init_test ("is_equal (const sp_pow &)");
	sort_test ();
	b.sort_test ();

	return ((this->first == b.first) &&
		(this->last == b.last) &&
		(this->coeff->lex_compare(*b.coeff) == 0));
}



//
//  *****  generating random series  *****
//

#if 0
friend
void
randomize (base_sparse_power_series< T > & a, lidia_size_t f, lidia_size_t l, lidia_size_t all, const T & max)
{
	debug_handler ("base_sparse_power_series< T >", "random (sp_pow &, lidia_size_t, lidia_size_t, lidia_size_t, const T &)");

	lidia_size_t i;
	lidia_size_t e;
	bigint L(l), E;

	T one_T(1);

	if (all > 0) { 
		a.coeff->set_capacity (all);

		if (all > 0) { 
			(*a.coeff)[0].exp = f;
			(*a.coeff)[0].coeff = randomize (max) + one_T;
		}

		for (i = 1; i < all; i++) { 
			do { 
				E = randomize (L);
				e = static_cast<lidia_size_t>(E.least_significant_digit());
			}
			while (e < f);

			(*a.coeff)[i].exp = e;
			(*a.coeff)[i].coeff = randomize (max) + one_T;
		}

		a.last = l;
		a.rebuild ();
	}
	else { 
		a.assign_zero (l);
	}
}
#endif


//
//  *****  some protected member functions  *****
//

template< class T >
void
base_sparse_power_series< T >::init_test (char * mess) const
{
	if (this->coeff->capacity() == 0) {
		char errmess[256];
		strcpy  (errmess, "uninitialized argument in function ");
		strncat (errmess, mess, 200);
		lidia_error_handler ("base_sparse_power_series< T >", errmess);
	}
}



template< class T >
void
base_sparse_power_series< T >::rebuild ()
{
	lidia_size_t i, j;

	this->coeff->sort();
	this->sorted = true;

	// --- only non-zero elements in the vector ---

	for (i = 0; i< this->coeff->size()-1;) {
		j = 1;

		while ((i+j < this->coeff->size()) && ((*this->coeff)[i].exp == (*this->coeff)[i+j].exp)) {
			(*this->coeff)[i].coeff += (*this->coeff)[i+j].coeff;
			j ++;
		}
		if (((*this->coeff)[i].coeff).is_zero ()) {
			this->coeff->remove_from (i, j);
		}
		else if (j > 1) {
			this->coeff->remove_from (i+1, j-1);
			i ++;
		}
		else {
			i ++;
		}
	}

	if (this->coeff->size() == 0) {
		this->assign_zero (this->last);
	}
	else {
		this->first = (*this->coeff)[0].exp;
	}
}



//
//  *****  arithmetical operations  *****
//

template< class T >
void
base_sparse_power_series< T >::multiply_by_xn (lidia_size_t n)
{
	debug_handler ("base_sparse_power_series< T >", "multiply_by_xn (lidia_size_t)");

	init_test ("multiply_by_xn (lidia_size_t)");

	lidia_size_t i, sz;

	if (n != 0) {
		sz = this->coeff->size();
		for (i = 0; i < sz; i++) {
			(*this->coeff)[i].exp += n;
		}

		this->first += n;
		this->last += n;
	}
}



template< class T >
void
base_sparse_power_series< T >::compose (lidia_size_t n)
{
	debug_handler ("base_sparse_power_series< T >", "compose (lidia_size_t)");

	init_test ("compose (lidia_size_t)");

	lidia_size_t i, sz;

	if (n != 0) {
		sz = this->coeff->size();
		for (i = 0; i < sz; i++) {
			(*this->coeff)[i].exp *= n;
		}

		this->first *= n;
		this->last *= n;
	}
}



template< class T >
void
base_sparse_power_series< T >::add (const base_sparse_power_series< T > & a,
				    const base_sparse_power_series< T > & b)
{
	debug_handler ("base_sparse_power_series< T >", "add (sp_pow &, const sp_pow &, const sp_pow &)");

	a.init_test ("add (3 x sp_pow)");
	b.init_test ("add (3 x sp_pow)");

	a.sort_test();
	b.sort_test();

	lidia_size_t i, ap, bp;
	lidia_size_t astop, bstop, cf, cl;

	lidia_size_t  cont, c_all;

	base_sparse_power_series< T > c;
	sort_vector< spc< T > > *ac, *bc, *cc;

	c_all = a.coeff->size() + b.coeff->size(); // max. number of non-zero elts. in result
	cf = comparator< lidia_size_t > ::min (a.first, b.first);
	cl = comparator< lidia_size_t > ::min (a.last, b.last);

	c.coeff->set_capacity (c.alloc_size (c_all));

	i = 0;
	cont = 0;
	ap = 0;
	bp = 0;
	astop = a.coeff->size() - 1;
	bstop = b.coeff->size() - 1;

	ac = a.coeff;
	bc = b.coeff;
	cc = c.coeff;

	// --- determine significant coefficients ---

	while ((astop >= 0) && ((*ac)[astop].exp > cl))
		astop--;
	while ((bstop >= 0) && ((*bc)[bstop].exp > cl))
		bstop--;


	while ((ap <= astop) && (bp <= bstop)) {
		if ((*ac)[ap].exp == (*bc)[bp].exp) {
			// add corresponding coefficients
			LiDIA::add ((*cc)[cont].coeff, (*ac)[ap].coeff, (*bc)[bp].coeff);

			if (! ((*cc)[cont].coeff).is_zero ()) {
				// store non-zero elements only
				(*cc)[cont].exp = (*ac)[ap].exp;
				cont++;
			}

			ap++;
			bp++;
		}
		else if ((*bc)[bp].exp > (*ac)[ap].exp) {
			// copy coeff. that occurs only once; least one first
			// copy coeff of series a

			(*cc)[cont].exp = (*ac)[ap].exp;
			(*cc)[cont].coeff = (*ac)[ap].coeff;

			cont++;
			ap++;
		}
		else if ((*ac)[ap].exp > (*bc)[bp].exp) {
			// copy coeff. that occurs only once; larger one later
			// copy coeff of series b

			(*cc)[cont].exp = (*bc)[bp].exp;
			(*cc)[cont].coeff = (*bc)[bp].coeff;

			cont++;
			bp++;
		}
	}

	while (ap <= astop) {
		// copy remaining coeff of series b
		(*cc)[cont].exp = (*ac)[ap].exp;
		(*cc)[cont].coeff = (*ac)[ap].coeff;

		cont++;
		ap++;
	}
	while (bp <= bstop) {
		// copy remaining coeff of series b
		(*cc)[cont].exp = (*bc)[bp].exp;
		(*cc)[cont].coeff = (*bc)[bp].coeff;

		cont++;
		bp++;
	}

	if (cont == 0) {
		c.assign_zero (cl);
	}
	else {
		c.coeff->set_capacity (cont);
		c.first = (*c.coeff)[0].exp;
		c.last = cl;
	}

	// --- move temporary result to 'res'

	c.swap (*this);
}



template< class T >
void
base_sparse_power_series< T >::subtract (const base_sparse_power_series< T > & a,
					 const base_sparse_power_series< T > & b)
{
	debug_handler ("base_sparse_power_series< T >", "subtract (sp_pow &, const sp_pow &, const sp_pow &)");

	a.init_test ("subtract (3 x sp_pow)");
	b.init_test ("subtract (3 x sp_pow)");

	a.sort_test();
	b.sort_test();


	lidia_size_t i, ap, bp;
	lidia_size_t astop, bstop, cf, cl;

	lidia_size_t  cont, c_all;

	base_sparse_power_series< T > c;
	sort_vector< spc< T > > *ac, *bc, *cc;

	c_all = a.coeff->size() + b.coeff->size(); // max. number of non-zero elts. in result
	cf = comparator< lidia_size_t >::min (a.first, b.first);
	cl = comparator< lidia_size_t >::min (a.last, b.last);

	c.coeff->set_capacity (c.alloc_size (c_all));

	i = 0;
	cont = 0;
	ap = 0;
	bp = 0;
	astop = a.coeff->size() - 1;
	bstop = b.coeff->size() - 1;

	ac = a.coeff;
	bc = b.coeff;
	cc = c.coeff;

	// --- determine significant coefficients ---

	while ((astop >= 0) && ((*ac)[astop].exp > cl))
		astop--;
	while ((bstop >= 0) && ((*bc)[bstop].exp > cl))
		bstop--;


	while ((ap <= astop) && (bp <= bstop)) {
		if ((*ac)[ap].exp == (*bc)[bp].exp) {
			// subtract corresponding coefficients
			LiDIA::subtract ((*cc)[cont].coeff, (*ac)[ap].coeff, (*bc)[bp].coeff);

			if (! ((*cc)[cont].coeff).is_zero ()) {
				// store non-zero elements only
				(*cc)[cont].exp = (*ac)[ap].exp;
				cont++;
			}

			ap++;
			bp++;
		}
		else if ((*bc)[bp].exp > (*ac)[ap].exp) {
			// copy coeff. that occurs only once; least one first
			// copy coeff of series a

			(*cc)[cont].exp = (*ac)[ap].exp;
			(*cc)[cont].coeff = (*ac)[ap].coeff;

			cont++;
			ap++;
		}
		else if ((*ac)[ap].exp > (*bc)[bp].exp) {
			// copy coeff. that occurs only once; larger one later
			// copy negative value of coeff of series b

			(*cc)[cont].exp = (*bc)[bp].exp;
			LiDIA::negate ((*cc)[cont].coeff, (*bc)[bp].coeff);

			cont++;
			bp++;
		}
	}

	while (ap <= astop) {
		// copy remaining coeff of series b
		(*cc)[cont].exp = (*ac)[ap].exp;
		(*cc)[cont].coeff = (*ac)[ap].coeff;

		cont++;
		ap++;
	}
	while (bp <= bstop) {	// copy remaining coeff of series b
		(*cc)[cont].exp = (*bc)[bp].exp;
		LiDIA::negate ((*cc)[cont].coeff, (*bc)[bp].coeff);

		cont++;
		bp++;
	}

	if (cont == 0) {
		c.assign_zero (cl);
	}
	else {
		c.coeff->set_capacity (cont);
		c.first = (*c.coeff)[0].exp;
		c.last = cl;
	}

	// --- move temporary result to 'res'

	c.swap(*this);
}



template< class T >
void
base_sparse_power_series< T >::negate (const base_sparse_power_series< T > & a)
{
	debug_handler ("base_sparse_power_series< T >", "negate (sp_pow &, const sp_pow &)");

	a.init_test ("negate (2 x sp_pow)");

	lidia_size_t i, sa;
	sort_vector< spc< T > > *ac, *rc;

	ac = a.coeff;
	rc = this->coeff;

	sa = ac->size();
	rc->set_capacity (sa);

	for (i = 0; i < sa; i++) {
		LiDIA::negate ((*rc)[i].coeff, (*ac)[i].coeff);
		(*rc)[i].exp = (*ac)[i].exp;
	}

	this->first = a.first;
	this->last = a.last;
	sort();
}



//
//  *****  scalar - operations  *****
//


template< class T >
void
base_sparse_power_series< T >::add (const base_sparse_power_series< T > & a,
				    const T                     & b)
{
	debug_handler ("base_sparse_power_series< T >", "add (sp_pow &, const sp_pow &, const T &)");

	a.init_test ("add (2 x sp_pow, cont T &)");
	a.sort_test();

	lidia_size_t at, where, rsz;
	bool found;

	lidia_size_t e = 0;
	T bb;
	spc< T > c(b, e);

	if ((a.last < e) || (b.is_zero ())) {
		if (this != &a) *this = a;
	}
	else {
		found = a.coeff->bin_search (c, where);

		if (this == &a) {
			if (found) {
				LiDIA::add ((*this->coeff)[where].coeff, (*this->coeff)[where].coeff, b);

				if (((*this->coeff)[where].coeff).is_zero ()) {
					this->coeff->remove_from (where);

					if (this->coeff->size() == 0)
						this->assign_zero (a.last);
					else
						this->first = (*this->coeff)[0].exp;
				}
			}
			else {
				if (this->coeff->capacity()< a.coeff->size() + 1)
					this->coeff->set_capacity (alloc_size (a.coeff->size() + 1));

				this->coeff->insert_at (c, where);
				this->first = (*this->coeff)[0].exp;
			}
		}
		else {
			// res does not alias a
			found = a.coeff->bin_search (c, where);

			rsz = a.coeff->size() + (found ? 0 : 1);
			if (this->coeff->capacity() < rsz)
				this->coeff->set_capacity (alloc_size (rsz));
			this->coeff->set_size (0);

			if (where > 0)
				this->coeff->assign (0, (*a.coeff), 0, where-1);

			if (found) {
				LiDIA::add (bb, (*a.coeff)[where].coeff, b);

				if (! bb.is_zero ()) {
					(*this->coeff)[where].coeff = bb;
					(*this->coeff)[where].exp = e;
					at = where + 1;
				}
				else
					at = where;

				if (where+1< a.coeff->size())
					this->coeff->assign (at, (*a.coeff), where+1, a.coeff->size()-1);

				if (this->coeff->size() == 0)
					this->assign_zero (a.last);
				else
					this->first = (*this->coeff)[0].exp;
			}
			else {
				(*this->coeff)[where] = c;

				if (where< a.coeff->size())
					this->coeff->assign (where+1, (*a.coeff), where, a.coeff->size()-1);

				this->first = (*this->coeff)[0].exp;
			}

			this->last = a.last;
		}
	}
}



template< class T >
void
base_sparse_power_series< T >::add (const T & b,
				    const base_sparse_power_series< T > & a)
{
	debug_handler ("base_sparse_power_series< T >", "add (sp_pow &, const T &, const sp_pow &)");
      	this->add(a, b);
}



template< class T >
void
base_sparse_power_series< T >::subtract (const base_sparse_power_series< T > & a,
					 const T                     & b)
{
	debug_handler ("base_sparse_power_series< T >", "add (sp_pow &, const sp_pow &, const T &)");

	a.init_test ("subtract (2 x sp_pow, const T &)");
	a.sort_test();

	lidia_size_t at, where, rsz;
	bool found;

	T bb;
	spc< T > c;
	lidia_size_t e = 0;

	LiDIA::negate (c.coeff, b);
	c.exp = e;

	if ((a.last < e) || (b.is_zero ())) {
		if (this != &a) *this = a;
	}
	else {
		found = a.coeff->bin_search (c, where);

		if (this == &a) {
			if (found) {
				LiDIA::subtract ((*this->coeff)[where].coeff, (*this->coeff)[where].coeff, b);

				if (((*this->coeff)[where].coeff).is_zero ()) {
					this->coeff->remove_from (where);

					if (this->coeff->size() == 0)
						this->assign_zero (a.last);
					else
						this->first = (*this->coeff)[0].exp;
				}
			}
			else {
				if (this->coeff->capacity()< a.coeff->size() + 1)
					this->coeff->set_capacity (alloc_size (a.coeff->size() + 1));

				this->coeff->insert_at (c, where);
				this->first = (*this->coeff)[0].exp;
			}
		}
		else {
			// res does not alias a
			found = a.coeff->bin_search (c, where);

			rsz = a.coeff->size() + (found ? 0 : 1);
			if (this->coeff->capacity() < rsz)
				this->coeff->set_capacity (alloc_size (rsz));
			this->coeff->set_size (0);

			if (where > 0)
				this->coeff->assign (0, (*a.coeff), 0, where-1);

			if (found) {
				LiDIA::subtract (bb, (*a.coeff)[where].coeff, b);

				if (! bb.is_zero ()) {
					(*this->coeff)[where].coeff = bb;
					(*this->coeff)[where].exp = e;
					at = where + 1;
				}
				else
					at = where;

				if (where+1< a.coeff->size())
					this->coeff->assign (at, (*a.coeff), where+1, a.coeff->size()-1);

				if (this->coeff->size() == 0)
					this->assign_zero (a.last);
				else
					this->first = (*this->coeff)[0].exp;
			}
			else {
				(*this->coeff)[where] = c;

				if (where< a.coeff->size())
					this->coeff->assign (where+1, (*a.coeff), where, a.coeff->size()-1);

				this->first = (*this->coeff)[0].exp;
			}

			this->last = a.last;
		}
	}
}



template< class T >
void
base_sparse_power_series< T >::subtract (const T                     & b,
					 const base_sparse_power_series< T > & a)
{
	debug_handler ("base_sparse_power_series< T >", "add (sp_pow &, const T &, const sp_pow &)");

	a.init_test ("subtract (sp_pow, cont T &, sp_pow)");
	a.sort_test();

	lidia_size_t at, where, rsz;
	bool found;
	lidia_size_t e = 0;
	T bb;
	spc< T > c(b, e);

	if ((a.last < e) || (b.is_zero ())) {
		if (this != &a) *this = a;
	}
	else {
		found = a.coeff->bin_search (c, where);

		if (this == &a) {
			if (found) {
				LiDIA::subtract ((*this->coeff)[where].coeff, b, (*this->coeff)[where].coeff);

				if (((*this->coeff)[where].coeff).is_zero ()) {
					this->coeff->remove_from (where);

					if (this->coeff->size() == 0)
						this->assign_zero (a.last);
					else
						this->first = (*this->coeff)[0].exp;
				}
			}
			else {
				if (this->coeff->capacity()< a.coeff->size() + 1)
					this->coeff->set_capacity (alloc_size (a.coeff->size() + 1));

				this->coeff->insert_at (c, where);
				this->first = (*this->coeff)[0].exp;
			}
		}
		else {
			// res does not alias a
			found = a.coeff->bin_search (c, where);

			rsz = a.coeff->size() + (found ? 0 : 1);
			if (this->coeff->capacity() < rsz)
				this->coeff->set_capacity (alloc_size (rsz));
			this->coeff->set_size (0);

			if (where > 0)
				this->coeff->assign (0, (*a.coeff), 0, where-1);

			if (found) {
				LiDIA::subtract (bb, b, (*a.coeff)[where].coeff);

				if (! bb.is_zero ()) {
					(*this->coeff)[where].coeff = bb;
					(*this->coeff)[where].exp = e;
					at = where + 1;
				}
				else
					at = where;

				if (where+1< a.coeff->size())
					this->coeff->assign (at, (*a.coeff), where+1, a.coeff->size()-1);

				if (this->coeff->size() == 0)
					this->assign_zero (a.last);
				else
					this->first = (*this->coeff)[0].exp;
			}
			else {
				(*this->coeff)[where] = c;

				if (where< a.coeff->size())
					this->coeff->assign (where+1, (*a.coeff), where, a.coeff->size()-1);

				this->first = (*this->coeff)[0].exp;
			}

			this->last = a.last;
		}
	}
}



template< class T >
void
base_sparse_power_series< T >::multiply (const base_sparse_power_series< T > & a,
					 const T                     & b)
{
	debug_handler ("base_sparse_power_series< T >", "multiply (sp_pow &, const sp_pow &, const T &)");

	a.init_test ("multiply (2 x sp_pow, cont T &)");
	a.sort_test ();

	lidia_size_t i, j, sz;
	T  zero_T;

	zero_T.assign_zero ();

	sort_vector< spc< T > > *ac, *rc;

	if (a.is_zero () || b == zero_T) {
		this->assign_zero (a.last);
	}
	else {
		sz = a.coeff->size();
		if ((this != &a) && (this->coeff->capacity() < sz))
			this->coeff->set_capacity(alloc_size (sz));

		ac = a.coeff;
		rc = this->coeff;

		for (i = j = 0; i < sz; i++) {
			LiDIA::multiply ((*rc)[j].coeff, (*ac)[i].coeff, b);

			if ((*rc)[j].coeff != zero_T) {
				(*rc)[j].exp = (*ac)[i].exp;
				j ++;
			}
		}

		if (j == 0)
			this->assign_zero (a.last);
		else {
			this->coeff->set_size (j);
			this->first = (*rc)[0].exp;
			this->last = a.last;

			this->sorted = true;
		}
	}
}



template< class T >
void
base_sparse_power_series< T >::multiply (const T                     & b,
					 const base_sparse_power_series< T > & a)
{
	debug_handler ("base_sparse_power_series< T >", "multiply (sp_pow &, const T &, const sp_pow &)");

	a.init_test ("multiply (sp_pow, cont T &, sp_pow)");
	a.sort_test ();

	lidia_size_t i, j, sz;
	T  zero_T;

	zero_T.assign_zero ();

	sort_vector< spc< T > > *ac, *rc;

	if (a.is_zero () || b == zero_T) {
		this->assign_zero (a.last);
	}
	else {
		sz = a.coeff->size();
		if ((this != &a) && (this->coeff->capacity() < sz))
			this->coeff->set_capacity(alloc_size (sz));

		ac = a.coeff;
		rc = this->coeff;

		for (i = j = 0; i < sz; i++) {
			LiDIA::multiply ((*rc)[j].coeff, b, (*ac)[i].coeff);

			if ((*rc)[j].coeff != zero_T) {
				(*rc)[j].exp = (*ac)[i].exp;
				j ++;
			}
		}

		if (j == 0)
			this->assign_zero (a.last);
		else {
			this->coeff->set_size (j);
			this->first = (*rc)[0].exp;
			this->last = a.last;

			this->sorted = true;
		}
	}
}



template< class T >
void
base_sparse_power_series< T >::divide (const base_sparse_power_series< T > & a,
				       const T                     & b)
{
	debug_handler ("base_sparse_power_series< T >", "divide (sp_pow &, const sp_pow &, const T &)");
	a.init_test ("divide (2 x sp_pow, cont T &)");

	lidia_size_t i, sz;
	sort_vector< spc< T > > *ac, *rc;


	if (b.is_zero ())  // this should be : b is not invertible
		lidia_error_handler ("base_sparse_power_series< T >", "divide (sp_pow &, const sp_pow &, const T &)::non-invertible element");

	sz = a.coeff->size();
	ac = a.coeff;
	rc = this->coeff;

	if (this != &a) {
		if (a.is_zero())
			this->assign_zero (a.last);
		else
			if (rc->capacity() < sz)
				rc->set_capacity (alloc_size (sz));
	}

	for (i = 0; i < sz; i++) {
		LiDIA::divide ((*rc)[i].coeff, (*ac)[i].coeff, b);
		(*rc)[i].exp = (*ac)[i].exp;
	}

	this->first = a.first;
	this->last = a.last;
	this->sorted = a.sorted;
}



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_BASE_SPARSE_POWER_SERIES_CC_GUARD_
