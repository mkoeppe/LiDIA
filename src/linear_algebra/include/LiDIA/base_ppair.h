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
//	Author	: Thomas Papanikolaou (TP), Thomas Pfahler (TPf)
//		  Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_BASE_PPAIR_H_GUARD_
#define LIDIA_BASE_PPAIR_H_GUARD_


#ifndef HEADBANGER

#ifndef LIDIA_LIDIA_H_GUARD_
# include	"LiDIA/LiDIA.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class T1, class T2 >
class base_ppair
{
	// MM, for debugging, FIX ME
public:

	T1 *l;
	T2 r;

public:

	//
	// constructors
	//
	base_ppair();
	base_ppair(T1 &pl, T2 &pr);
	base_ppair(const base_ppair< T1, T2 > &p);

	//
	// destructor
	//
	~base_ppair();

	//
	// assignment
	//
	void assign (const base_ppair< T1, T2 > &p);
	base_ppair< T1, T2 > & operator = (const base_ppair< T1, T2 > &p);

	//
	// component access
	//
	const T1 & left() const;
	const T2 & right() const;

	T1 & left();
	T2 & right();

	//
	// INPUT / OUTPUT
	//
	void read(std::istream & in);
	void print(std::ostream & out) const;


	//
	// boolean operators
	//
	bool operator == (const base_ppair< T1, T2 > & p) const;
	bool operator != (const base_ppair< T1, T2 > & p) const;

	//
	// swap function
	//

	void swap(base_ppair< T1, T2 > &);

};



//
// constructor
//

template< class T1, class T2 >
inline
base_ppair< T1, T2 >::base_ppair ()
{
	debug_handler("base_ppair", "base_ppair()");

	l = new T1;
	memory_handler(l, "base_ppair", "base_ppair :: "
		       "Error in memory allocation (l)");
	r = 0;
}



template< class T1, class T2 >
inline
base_ppair< T1, T2 >::base_ppair (T1 & pl, T2 & pr)
{
	l = new T1;
	memory_handler(l, "base_ppair", "base_ppair :: "
		       "Error in memory allocation (l)");
	*l = pl;
	r = pr;
}



template< class T1, class T2 >
inline
base_ppair< T1, T2 >::base_ppair (const base_ppair< T1, T2 > & p)
{
	l = new T1;
	memory_handler(l, "base_ppair", "base_ppair :: "
		       "Error in memory allocation (l)");
	*l = *p.l;
	r = p.r;
}



//
// destructor
//

template< class T1, class T2 >
inline
base_ppair< T1, T2 >::~base_ppair ()
{
	debug_handler("base_ppair", "~base_ppair");
	delete l;
}



//
// component access
//

template< class T1, class T2 >
inline const T1 &
base_ppair< T1, T2 >::left() const
{
	return *l;
}



template< class T1, class T2 >
inline const T2 &
base_ppair< T1, T2 >::right () const
{
	return r;
}



template< class T1, class T2 >
inline T1 &
base_ppair< T1, T2 >::left ()
{
	return *l;
}



template< class T1, class T2 >
inline T2 &
base_ppair< T1, T2 >::right ()
{
	return r;
}



//
// assignment
//

template< class T1, class T2 >
inline void
base_ppair< T1, T2 >::assign (const base_ppair< T1, T2 > & p)
{
	debug_handler("base_ppair", "assign (const base_ppair< T1, T2 > &)");

	*l = *p.l;
	r = p.r;
}



template< class T1, class T2 >
inline base_ppair< T1, T2 > &
base_ppair< T1, T2 >::operator = (const base_ppair< T1, T2 > & p)
{
	debug_handler("base_ppair", "operator = (const base_ppair< T1, T2 > &)");

	assign(p);
	return *this;
}



//
// boolean operator
//

template< class T1, class T2 >
inline bool
base_ppair< T1, T2 >::operator == (const base_ppair< T1, T2 > & p) const
{
	return ((&p == this) || ((*l == *p.l) && (r == p.r)));
}



template< class T1, class T2 >
inline bool
base_ppair< T1, T2 >::operator != (const base_ppair< T1, T2 > &p) const
{
	return !(this->operator == (p));
}



template< class T1, class T2 >
inline void
swap(base_ppair< T1, T2 > &a, base_ppair< T1, T2 > &b)
{
	a.swap(b);
}



template< class T1, class T2 >
inline std::istream &
operator >> (std::istream &is, base_ppair< T1, T2 > &a)
{
	a.read(is);
	return is;
}



template< class T1, class T2 >
inline std::ostream &
operator << (std::ostream &os, const base_ppair< T1, T2 > &a)
{
	a.print(os);
	return os;
}



#endif	// HEADBANGER



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_BASE_PPAIR_H_GUARD_
