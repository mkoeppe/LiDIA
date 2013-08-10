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
//	Author	: Markus Maurer (MM), Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_BASE_POWER_PRODUCT_H_GUARD_
#define LIDIA_BASE_POWER_PRODUCT_H_GUARD_


#ifndef LIDIA_BASE_VECTOR_H_GUARD_
# include	"LiDIA/base_vector.h"
#endif
#ifndef LIDIA_BASE_PPAIR_H_GUARD_
# include	"LiDIA/base_ppair.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class power_product_flags
{
protected:

	enum base_power_product_state {
		normalized = 2
	};
};



template< class T, class exp_type >
class base_power_product : public power_product_flags
{
	//
	//  data elements
	//

protected:

	base_vector< base_ppair< T, exp_type > > component;

	unsigned long attributes; // status information

	//
	//  access functions for 'attributes'
	//

public:

	void clear_attributes()
	{
		attributes &= ~ (normalized);
	}

	//
	//  constructors
	//

public:

	base_power_product();
	base_power_product(const base_power_product< T, exp_type > &);
	base_power_product(const T &);

	//
	// destructor
	//

public:

	~base_power_product();

	//
	//  size control
	//

public:

	void reset();

	void set_exp_ratio(float);
	void set_capacity(lidia_size_t);

	//
	//  assignments
	//

public:

	base_power_product< T, exp_type > & operator = (const base_power_product< T, exp_type > &);
	base_power_product< T, exp_type > & operator = (const T &);

#ifndef HEADBANGER
	void assign(const base_power_product< T, exp_type > &);
	void assign(const T &);
#endif

	void assign_one();

	//
	//  access functions
	//

public:

	const T & get_base(lidia_size_t) const;
	exp_type  get_exponent(lidia_size_t) const;

	lidia_size_t get_no_of_components() const
	{
		return component.size();
	}

	//
	//  swap function
	//

public:

	void swap(base_power_product< T, exp_type > &);

	//
	// input / output
	//

	//
	// I/O - format :
	//
	// [ [base_1, exp_1], [base_2, exp_2], ..., [base_i, exp_i] ]
	//

public:

	void read (std::istream &);
	void write(std::ostream &) const;

	//
	//  functions for computing with base_power_products
	//

public:

#ifndef HEADBANGER
	void concat (const base_power_product< T, exp_type > &a,
		     const base_power_product< T, exp_type > &b);
#endif

	//
	//  high-level functions
	//

public:

	void append(const T &a)
	{
		append(a, 1);
	}

	void append(const T &, exp_type); // appends 'a^exp' to the base_power_product
};



//
// c'tors and d'tor
//

template< class T, class exp_type >
inline
base_power_product< T, exp_type >::base_power_product ()
	: component(0, EXPAND)
{
	// nothing to do
}



template< class T, class exp_type >
inline
base_power_product< T, exp_type >::base_power_product (const T & f)
	: component(0, EXPAND)
{
	assign(f);
}



template< class T, class exp_type >
inline
base_power_product< T, exp_type >::base_power_product (const base_power_product< T, exp_type > &f)
	: component(f.component, EXPAND),
	  attributes(f.attributes)
{
	// nothing to do
}



template< class T, class exp_type >
inline
base_power_product< T, exp_type >::~base_power_product()
{
	// nothing to do
}



//
// size control
//

template< class T, class exp_type >
inline void
base_power_product< T, exp_type >::reset ()
{
	component.reset();
}



template< class T, class exp_type >
inline void
base_power_product< T, exp_type >::set_exp_ratio (float r)
{
	component.set_exp_ratio(r);
}



template< class T, class exp_type >
inline void
base_power_product< T, exp_type >::set_capacity (lidia_size_t c)
{
	component.set_capacity(c);
}



//
//  assignments
//

template< class T, class exp_type >
inline void
base_power_product< T, exp_type >::assign_one ()
{
	component.set_capacity(0);
}



template< class T, class exp_type >
inline void
base_power_product< T, exp_type >::assign(const T & f)
{
	component[0].left() = f;
	component[0].right() = 1;
	component.set_size(1);
}



template< class T, class exp_type >
inline void
base_power_product< T, exp_type >::assign (const base_power_product< T, exp_type > & f)
{
	if (this != &f) {
		component.assign(f.component);
		attributes = f.attributes;
	}
}



template< class T, class exp_type >
inline base_power_product< T, exp_type > &
base_power_product< T, exp_type >::operator = (const base_power_product< T, exp_type > & f)
{
	assign(f);
	return *this;
}



template< class T, class exp_type >
inline base_power_product< T, exp_type > &
base_power_product< T, exp_type >::operator = (const T & f)
{
	assign(f);
	return *this;
}



//
//  access functions
//

template< class T, class exp_type >
inline const T &
base_power_product< T, exp_type >::get_base (lidia_size_t index) const
{
	if ((index< 0) || (index >= get_no_of_components()))
		lidia_error_handler("base_power_product< T, exp_type >",
				    "base::index out of range");

	return component[index].left();
}



template< class T, class exp_type >
inline exp_type
base_power_product< T, exp_type >::get_exponent (lidia_size_t index) const
{
	if ((index< 0) || (index >= get_no_of_components()))
		lidia_error_handler("base_power_product< T, exp_type >",
				    "exponent::index out of range");

	return component[index].right();
}



//
//
//

template< class T, class exp_type >
inline void
base_power_product< T, exp_type >::concat (const base_power_product< T, exp_type > & a,
					   const base_power_product< T, exp_type > & b)
{
	component.concat(a.component, b.component);
	clear_attributes();
}



template< class T, class exp_type >
inline void
base_power_product< T, exp_type >::swap (base_power_product< T, exp_type > & b)
{
	LiDIA::swap(component, b.component);
	LiDIA::swap(attributes, b.attributes);
}



template< class T, class exp_type >
inline void
swap (base_power_product< T, exp_type > & a, base_power_product< T, exp_type > & b)
{
	a.swap(b);
}



template< class T, class exp_type >
inline std::istream &
operator >> (std::istream &is, base_power_product< T, exp_type > & F)
{
	F.read(is);
	return is;
}



template< class T, class exp_type >
inline std::ostream &
operator << (std::ostream &os, const base_power_product< T, exp_type > & F)
{
	F.write(os);
	return os;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#if defined LIDIA_INCLUDE_CC
# include	"LiDIA/base_power_product.cc"
#endif



#endif // LIDIA_BASE_POWER_PRODUCT_H_GUARD_
