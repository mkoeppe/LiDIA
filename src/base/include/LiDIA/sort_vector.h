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
//	Author	: Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_SORT_VECTOR_H_GUARD_
#define LIDIA_SORT_VECTOR_H_GUARD_



#ifndef LIDIA_COMPARATOR_H_GUARD_
# include	"LiDIA/comparator.h"
#endif
#ifndef LIDIA_BASE_VECTOR_H_GUARD_
# include	"LiDIA/base_vector.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class T >
class sort_vector : public virtual base_vector< T >
{

	//
	// constructors
	//

public:

	sort_vector() :
		base_vector< T > ()
	{
		this->sort_dir = vector_flags::sort_vector_up;
		this->el_cmp = NULL;
	}

	explicit  sort_vector(const vector_flags & md) :
		base_vector< T > (md)
	{
		this->sort_dir = vector_flags::sort_vector_up;
		this->el_cmp = NULL;
	}
	explicit  sort_vector(lidia_size_t all) :
		base_vector< T > (all)
	{
		this->sort_dir = vector_flags::sort_vector_up;
		this->el_cmp = NULL;
	}

	sort_vector(lidia_size_t all, const vector_flags & md) :
		base_vector< T > (all, md)
	{
		this->sort_dir = vector_flags::sort_vector_up;
		this->el_cmp = NULL;
	}
	sort_vector(lidia_size_t all, lidia_size_t len) :
		base_vector< T > (all, len)
	{
		this->sort_dir = vector_flags::sort_vector_up;
		this->el_cmp = NULL;
	}
	sort_vector(lidia_size_t all, lidia_size_t len, const vector_flags & md) :
		base_vector< T > (all, len, md)
	{
		this->sort_dir = vector_flags::sort_vector_up;
		this->el_cmp = NULL;
	}
	sort_vector(const sort_vector< T > &v) :
		base_vector< T > (v)
	{
		this->sort_dir = v.sort_dir;
		this->el_cmp = v.el_cmp;
	}
	sort_vector(const base_vector< T > &v, const vector_flags & md) :
		base_vector< T > (v, md)
	{
		this->sort_dir = vector_flags::sort_vector_up;
		this->el_cmp = NULL;
	}
	sort_vector(const T *v, lidia_size_t len) :
		base_vector< T > (v, len)
	{
		this->sort_dir = vector_flags::sort_vector_up;
		this->el_cmp = NULL;
	}
	sort_vector(const T *v, lidia_size_t len, const vector_flags &md) :
		base_vector< T > (v, len, md)
	{
		this->sort_dir = vector_flags::sort_vector_up;
		this->el_cmp = NULL;
	}

	//
	// destructor
	//

public:

	~sort_vector() {}

	//
	// assign
	//

	sort_vector< T > & operator = (const sort_vector< T > & v);

	//
	// reading & modifying sort-directions
	//

public:

	unsigned long sort_direction() const
	{
		return  this->sort_dir;
	}
	unsigned long get_sort_direction() const
	{
		return  this->sort_dir;
	}

	void set_sort_direction(unsigned long);
	void set_sort_direction(int (*cmp) (const T &, const T &));

	//
	// sort-functions using Quick-Sort routines
	//

public:

	void sort(int (*cmp)(const T &, const T &),
		  lidia_size_t l = 0, lidia_size_t r = -1);
	void sort(unsigned long sort_direction = vector_flags::sort_vector_def,
		  lidia_size_t l = 0, lidia_size_t r = -1);
	void sort(lidia_size_t l, lidia_size_t r)
	{
		sort(vector_flags::sort_vector_def, l, r);
	}

	void sort_down(lidia_size_t l, lidia_size_t r);
	void sort_up(lidia_size_t l, lidia_size_t r);

	//
	// functions for linear search in vectors
	//

public:

	bool linear_search(const T &, lidia_size_t &) const;

	//
	// functions for binary search in sorted vectors
	//

public:

	bool bin_search(const T &, lidia_size_t &,
			int (*cmp)(const T &, const T &),
			lidia_size_t l = 0, lidia_size_t r = -1) const;
	bool bin_search(const T &, lidia_size_t &,
			unsigned long sort_direction = vector_flags::sort_vector_def,
			lidia_size_t l = 0, lidia_size_t r = -1) const;
	bool bin_search(const T & x, lidia_size_t & pos,
			lidia_size_t l, lidia_size_t r) const
	{
		return this->bin_search(x, pos, vector_flags::sort_vector_def,
					l, r);
	}

protected:

	bool bin_search_up(const  T &, lidia_size_t &, lidia_size_t, lidia_size_t) const;
	bool bin_search_down(const  T &, lidia_size_t &, lidia_size_t, lidia_size_t) const;

	//
	// functions to insert new elts. into a vector
	//

public:

	void insert(const T &, int (*cmp)(const T &, const T &),
		    lidia_size_t l = 0, lidia_size_t r = -1);
	void insert(const T &, unsigned long sort_direction = vector_flags::sort_vector_def,
		    lidia_size_t l = 0, lidia_size_t r = -1);
	void insert(const T & x, lidia_size_t l, lidia_size_t r)
	{
		this->insert(x, vector_flags::sort_vector_def, l, r);
	}

	//
	// functions to remove elts. from within a vector
	//

public:

	bool remove(const T &, int (*cmp)(const T &, const T &),
		    lidia_size_t l = 0, lidia_size_t r = -1);
	bool remove(const T &, unsigned long sort_direction = vector_flags::sort_vector_def,
		    lidia_size_t l = 0, lidia_size_t r = -1);
	bool remove(const T & x, lidia_size_t l, lidia_size_t r)
	{
		return this->remove(x, vector_flags::sort_vector_def, l, r);
	}

	int lex_compare(sort_vector< T > &) const;

	//
	// miscellaneous
	//

public:

	void delete_copies();
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/sort_vector.cc"
#endif



#endif	// LIDIA_SORT_VECTOR_H_GUARD_
